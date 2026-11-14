"""Audit saved fog/dust cases and test representative numerical sensitivity.

Run after the explorer has finished. Calculations are sequential. A successful
file audit is distinct from a converged physical prediction.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np

from akutou_sky import (ModelConfig, load_dataset, xyz_from_radiance,
                        validate_sensitivity)

def write_summary(directory, report):
    """Write a readable summary without concealing failed or skipped checks."""
    names = dict(urban='都会', rural='田舎', humid='高湿度', dry='低湿度')
    check_names = dict(angular='角度積分', spectral='波長', vertical='鉛直格子',
                       horizontal='水平格子', iterations='反復停止条件',
                       moments='Mie係数', delta_m='Delta-M次数', single_scatter='一回散乱の独立格子')
    checks = [(n, k, c) for n, s in report['scenarios'].items() for k, c in s['checks'].items()]
    completed = sum(c['status'] == 'completed' for _, _, c in checks)
    passed = sum(c['passes_5percent_Y_and_0_005_xy'] for _, _, c in checks)
    lines = ['# 数値感度の検証（0.3.0 / モデル2.2）', '',
        '伏角−6°〜＋6°、0.5°刻みの25点を保存しています。正の伏角は太陽が地平線の下にあることを表します。', '',
        f'代表条件の検査は {len(checks)} 件中 {completed} 件が完了し、{passed} 件が'
        '「最大輝度差5%以内、かつ色度xy距離0.005以内」を満たしました。', '',
        '格子を変えたときの差であり、実際の空に対する誤差の上限ではありません。'
        '全64条件の収束は保証しません。基準未達や未完了の検査を合格として扱っていません。', '',
        '## 各検査の最大差', '',
        '| 大気 | 検査 | 最大輝度差 | 最大色度差 | 判定 |',
        '|---|---|---:|---:|---|']
    for name, kind, c in checks:
        if c['status'] == 'completed':
            y = f"{100*c['max_relative_Y']:.3f}%"
            xy = f"{c['max_delta_xy']:.6f}"
            status = '基準内' if c['passes_5percent_Y_and_0_005_xy'] else '基準未達'
        else:
            y, xy, status = '—', '—', c['status'] + ': ' + c.get('reason', '').replace('|', '/')
        lines.append(f'| {names[name]} | {check_names[kind]} | {y} | {xy} | {status} |')
    old_path = directory/'history/0.2.0/numerical_sensitivity.json'
    if old_path.exists():
        old = json.loads(old_path.read_text())
        lines += ['', '## 旧版との比較（共通の伏角）', '',
            '旧版の14→26角度点と、新版の110→146角度点の感度を示します。'
            '数値手法と基準格子も変更しているため、実測誤差の改善率ではありません。', '',
            '| 大気 | 伏角 | 旧版の最大輝度差 | 新版の最大輝度差 |',
            '|---|---:|---:|---:|']
        for name in names:
            before = {r['depression_deg']:r['all']['max_relative_Y'] for r in
                      old['scenarios'][name]['checks']['angular'].get('per_depression', [])}
            after = {r['depression_deg']:r['all']['max_relative_Y'] for r in
                     report['scenarios'][name]['checks']['angular'].get('per_depression', [])}
            for d in (0., 6.):
                if d in before and d in after:
                    lines.append(f'| {names[name]} | {d:+g}° | {100*before[d]:.3f}% | {100*after[d]:.3f}% |')
    lines += ['', '## 検査範囲と保存データ', '',
        '- 初期条件4種類について、伏角−6、−3、0、＋3、＋6°、方位0、90、180°、高度0、1、5、12、45、90°を検査しました。地平線の視線も含みます。',
        '- [数値感度の全記録](numerical_sensitivity.json)に、変更設定、伏角ごとの差、地平線付近と高度5°以上の内訳を保存しています。',
        '- [全64条件のデータ検査](artifact_checks.json)で、非負・有限な放射輝度、分光からのXYZ再計算、AODの積分、HTML用データとの一致を確認しています。',
        '- [改善した手法と方程式](numerical-refinement.md)、[全64条件の再現記録](manifest.json)も参照してください。', '',
        '再検証: `python tools/validate_explorer.py --directory results/atmosphere_explorer --cache-dir model_cache`', '']
    (directory/'validation_summary.md').write_text('\n'.join(lines), encoding='utf-8')

def audit(directory):
    manifest = json.loads((directory/'manifest.json').read_text())
    data = json.loads((directory/'data.json').read_text())
    records = manifest['calculations']
    assert manifest['all_64_cases_complete'] and len(records) == 64
    assert len({v['path'] for v in records}) == 64
    for page in [directory/'index.html', *(directory/s['key']/'index.html' for s in data['scenarios'])]:
        embedded = re.search(r'<script id="modelData" type="application/json">(.*?)</script>',
                             page.read_text(), re.S)
        assert embedded is not None, page
        assert json.loads(embedded.group(1)) == data, page
    for scenario in data['scenarios']:
        baseline = load_dataset(directory/scenario['key']/'sky.nc')
        chosen = load_dataset(directory/scenario['cases'][scenario['dust']][scenario['fog']]['path'])
        np.testing.assert_array_equal(baseline.solar_stokes_radiance, chosen.solar_stokes_radiance)
        assert baseline.attrs['config_json'] == chosen.attrs['config_json']
    maxima = []
    clear_reference = None
    max_zenith_relative = 0.
    max_negative_single_fraction = 0.
    for record in records:
        path = directory/record['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256'], path
        run = load_dataset(path)
        cfg = ModelConfig.from_dict(json.loads(run.attrs['config_json']))
        assert cfg.to_dict() == record['config']
        for attr, file in [('transport_source_sha256', 'model.py'),
                           ('microphysics_source_sha256', 'microphysics.py'),
                           ('delta_m_source_sha256', 'delta_m.py')]:
            assert run.attrs[attr] == manifest['source_hashes'][file], path
        radiance = run.solar_stokes_radiance.sel(stokes='I').values
        assert np.isfinite(radiance).all() and (radiance >= 0).all(), path
        if 'solar_single_scatter_radiance' in run:
            single=run.solar_single_scatter_radiance.sel(stokes='I').values
            assert np.isfinite(single).all(),path
            negative=single<0
            if negative.any():
                max_negative_single_fraction=max(max_negative_single_fraction,
                    float(np.max(-single[negative]/np.maximum(radiance[negative],1e-30))))
        xyz = xyz_from_radiance(radiance, run.wavelength.values, run.wavelength_bounds.values)
        np.testing.assert_allclose(xyz, run.total_XYZ, rtol=1e-13, atol=0)
        np.testing.assert_allclose(np.trapezoid(run.aerosol_extinction_550, run.altitude),
                                   run.attrs['aerosol_aod550_ambient'], rtol=1e-13, atol=1e-15)
        np.testing.assert_allclose(np.trapezoid(run.fog_extinction_550, run.altitude),
                                   cfg.fog_aod550, rtol=1e-13, atol=1e-15)
        scenario = next(s for s in data['scenarios'] if s['key'] == path.parent.name)
        i, j = data['dust'].index(cfg.aod550), data['fog'].index(cfg.fog_aod550)
        np.testing.assert_array_equal(scenario['cases'][i][j]['XYZ'], xyz)
        for name, values in [('depression', data['depressions']), ('azimuth', data['azimuths']),
                             ('elevation', data['elevations'])]:
            np.testing.assert_array_equal(run[name], values)
        if cfg.aod550 == 0 and cfg.fog_aod550 == 0:
            if clear_reference is None:
                clear_reference = radiance
            else:
                np.testing.assert_allclose(radiance, clear_reference, rtol=1e-13, atol=0)
        zenith = run.solar_stokes_radiance.sel(stokes='I', elevation=90).values
        relative = np.max(np.abs(zenith[:, 1:] - zenith[:, :1]) / np.maximum(zenith[:, :1], 1e-30))
        max_zenith_relative = max(max_zenith_relative, float(relative))
        maxima.append(dict(path=record['path'], min_Y=float(xyz[..., 1].min()), max_Y=float(xyz[..., 1].max())))
    for i, dust in enumerate(data['dust']):
        if dust > 0:
            for j in range(len(data['fog'])):
                assert data['scenarios'][2]['cases'][i][j]['ambient_aod'] > data['scenarios'][3]['cases'][i][j]['ambient_aod']
    report = dict(all_64_cases_checked=True, finite_nonnegative_radiance=True,
                  XYZ_recomputed_from_spectra=True, vertical_columns_correct=True,
                  embedded_data_matches_NetCDF=True, neutral_atmospheres_identical=True,
                  five_HTML_pages_match_data_JSON=True, baselines_match_selected_cases=True,
                  transport_source_hashes_match=True,
                  max_negative_single_component_fraction_of_total=max_negative_single_fraction,
                  higher_RH_increases_Mie_AOD_at_fixed_particle_count=True,
                  max_zenith_azimuth_relative_difference=max_zenith_relative, Y_ranges=maxima)
    (directory/'artifact_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    return report

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path('results/atmosphere_explorer'))
    parser.add_argument('--cache-dir', type=Path, default=Path('model_cache'))
    parser.add_argument('--max-phase-terms', type=int, default=520_000_000,
                        help='Verification phase-table work cap, not RAM bytes; 0 disables this extra cap')
    args = parser.parse_args()
    audit(args.directory)
    report = dict(scope='Numerical sensitivity of four baseline scenarios on sampled rays; '
                        'not a full convergence certificate or an observational error bound.', scenarios={})
    transport_cache = {}
    for name in ('urban', 'rural', 'dry', 'humid'):
        checks = ('angular',) if name != 'humid' else ('angular', 'spectral', 'vertical', 'horizontal', 'iterations', 'moments', 'delta_m', 'single_scatter')
        print(f'Checking {name}: {", ".join(checks)}', flush=True)
        checked = validate_sensitivity(args.directory/name/'sky.nc', checks=checks,
            output_dir=args.directory/'sensitivity'/name, cache_dir=args.cache_dir, allow_download=False,
            progress=lambda message: print(message, flush=True),
            max_phase_terms=args.max_phase_terms or None, _transport_cache=transport_cache)
        report['scenarios'][name] = checked
        report['all_requested_checks_complete'] = (len(report['scenarios']) == 4 and
            all(v['all_requested_checks_complete'] for v in report['scenarios'].values()))
        report['passes_all_sampled_checks'] = (report['all_requested_checks_complete'] and
            all(v['passes_all_sampled_checks'] for v in report['scenarios'].values()))
        (args.directory/'numerical_sensitivity.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('all_requested_checks_complete', 'passes_all_sampled_checks')}), flush=True)
    write_summary(args.directory, report)
    from akutou_sky.explorer import render_explorer
    render_explorer(args.directory)
    return 0 if report['all_requested_checks_complete'] else 1

if __name__ == '__main__':
    raise SystemExit(main())
