"""Reproducible four-scenario, discrete fog/dust explorer with offline HTML.

Every slider stop has a radiative-transfer result; physically identical
aerosol-free atmospheres share a calculation. Neither spectra nor scenario
colours are interpolated between different atmospheric conditions.
"""
import hashlib
import json
from pathlib import Path

from .config import ModelConfig

DUST_LEVELS = [0., .03, .10, .30]
FOG_LEVELS = [0., .05, .20, .80]
DEPRESSIONS = [i / 2 for i in range(-12, 13)]
ELEVATIONS = [0., .5, 1., 2., 3., 5., 8., 12., 20., 30., 45., 60., 90.]
AZIMUTHS = [0., 90., 180.]
SCENARIOS = [
    dict(key='urban', title='都会', caption='微粒子が多く、光を吸収しやすい条件',
         rh=.60, kappa=.20, index_real=1.50, index_imag=.010, scale_height=1200., dust=3, fog=0),
    dict(key='rural', title='田舎', caption='微粒子が少なく、空気が澄んだ条件',
         rh=.60, kappa=.30, index_real=1.45, index_imag=.003, scale_height=1500., dust=1, fog=0),
    dict(key='humid', title='湿度が高い', caption='水を吸った微粒子と、地表付近の水滴',
         rh=.90, kappa=.30, index_real=1.45, index_imag=.005, scale_height=1500., dust=2, fog=2),
    dict(key='dry', title='湿度が低い', caption='同じ乾燥粒子量で、吸湿成長が小さい条件',
         rh=.30, kappa=.30, index_real=1.45, index_imag=.005, scale_height=1500., dust=2, fog=0),
]

def scenario_config(scenario, dust, fog):
    """Create a full validated config for a scenario and two optical depths."""
    return ModelConfig.from_preset('fog_refined', aod550=dust, fog_aod550=fog,
        relative_humidity=scenario['rh'], hygroscopicity_kappa=scenario['kappa'],
        refractive_index_real=scenario['index_real'], refractive_index_imag=scenario['index_imag'],
        aerosol_scale_height_m=scenario['scale_height'])

def _source_hashes():
    root = Path(__file__).parent
    return {name: hashlib.sha256((root/name).read_bytes()).hexdigest()
            for name in ('model.py', 'microphysics.py', 'delta_m.py', 'config.py', 'colorimetry.py', 'explorer.py')}

def _compatible_signature(actual, expected):
    # Viewer/export code may change without invalidating transport. Scenario
    # configs and all coordinates are independently compared before reuse.
    try:
        left, right = json.loads(actual), json.loads(expected)
        for value in (left, right):
            value['sources'].pop('explorer.py', None)
        return left == right
    except (TypeError, ValueError, KeyError):
        return False

def build_explorer(directory, *, cache_dir=None, allow_download=True, resume=False, progress=None):
    """Compute 64 conditions sequentially, then write an offline viewer.

    With resume=True, cached runs are reused only if configuration, input data
    hashes, source hashes, backend version and coordinate axes match. A failed
    run leaves reusable completed cases, but no partially published HTML.
    """
    import importlib.metadata
    import numpy as np
    from .model import compute_model
    from .colorimetry import add_colorimetry
    from .io import load_dataset, save_dataset, export_results
    from .data import input_paths, resolve_cache_dir
    from ._version import __version__
    directory = Path(directory)
    if directory.exists() and any(directory.iterdir()) and not resume:
        raise FileExistsError(f'{directory} exists; use resume=True / --resume to verify and reuse cases')
    directory.mkdir(parents=True, exist_ok=True)
    sources = _source_hashes()
    cache = resolve_cache_dir(cache_dir)
    inputs = {str(p.relative_to(cache)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in input_paths(cache, allow_download=allow_download)}
    signature = json.dumps(dict(sources=sources, inputs=inputs, backend=importlib.metadata.version('sasktran2'),
                                 library=__version__), sort_keys=True)
    payload = dict(dust=DUST_LEVELS, fog=FOG_LEVELS, depressions=DEPRESSIONS,
                   elevations=ELEVATIONS, azimuths=AZIMUTHS, scenarios=[], library_version=__version__, model_version='2.2')
    records = []
    completed = 0
    transport_cache = {}
    equivalent_runs = {}
    # With no aerosol particles, aerosol-layer RH, composition and scale
    # height do not enter transport (gas is fixed US76; fog is independent).
    inactive_fields = {'relative_humidity', 'hygroscopicity_kappa',
                       'refractive_index_real', 'refractive_index_imag',
                       'aerosol_scale_height_m'}
    for scenario in SCENARIOS:
        cases = []
        for i, dust in enumerate(DUST_LEVELS):
            fog_cases = []
            for j, fog in enumerate(FOG_LEVELS):
                config = scenario_config(scenario, dust, fog)
                relative = Path('calculations')/scenario['key']/f'dust{i}_fog{j}.nc'
                path = directory/relative
                run = None
                if resume and path.exists():
                    candidate = load_dataset(path)
                    if (_compatible_signature(candidate.attrs.get('explorer_signature'), signature)
                        and json.loads(candidate.attrs['config_json']) == config.to_dict()
                        and all(np.array_equal(candidate[name], values) for name, values in
                                [('depression', DEPRESSIONS), ('azimuth', AZIMUTHS), ('elevation', ELEVATIONS)])):
                        stored_inputs = json.loads(candidate.attrs['physical_inputs_json'])
                        if all((cache/v['path']).is_file() and
                               hashlib.sha256((cache/v['path']).read_bytes()).hexdigest() == v['sha256']
                               for v in stored_inputs):
                            run = candidate
                reused = run is not None
                equivalent_key = (json.dumps({k: v for k, v in config.to_dict().items()
                                              if k not in inactive_fields}, sort_keys=True)
                                  if dust == 0 else None)
                if run is None and equivalent_key in equivalent_runs:
                    from .microphysics import aerosol_water_uptake
                    source, original_run = equivalent_runs[equivalent_key]
                    run = original_run.copy(deep=True)
                    growth, radius, _ = aerosol_water_uptake(config)
                    run.attrs.update(config_json=json.dumps(config.to_dict(), sort_keys=True),
                        aerosol_growth_factor=growth, aerosol_wet_median_radius_nm=radius,
                        same_physics_reused_from=source,
                        physical_equivalence_ignored_fields_json=json.dumps(sorted(inactive_fields)))
                    save_dataset(run, path, overwrite=path.exists())
                if run is None:
                    run = add_colorimetry(compute_model(DEPRESSIONS, config=config,
                                   elevations=ELEVATIONS, azimuths=AZIMUTHS,
                                   cache_dir=cache, allow_download=allow_download,
                                   _transport_cache=transport_cache))
                    run.attrs['explorer_signature'] = signature
                    save_dataset(run, path, overwrite=path.exists())
                if equivalent_key is not None:
                    equivalent_runs.setdefault(equivalent_key, (relative.as_posix(), run))
                xyz = run.total_XYZ.values
                if not np.isfinite(xyz).all() or (xyz < 0).any():
                    raise ValueError(f'Invalid XYZ in {relative}; viewer not published')
                fog_cases.append(dict(XYZ=xyz.tolist(), gamut=run.total_gamut_clipped.values.tolist(),
                    ambient_aod=run.attrs['aerosol_aod550_ambient'], growth=run.attrs['aerosol_growth_factor'],
                    wet_radius=run.attrs['aerosol_wet_median_radius_nm'], path=relative.as_posix(), config=config.to_dict()))
                records.append(dict(path=relative.as_posix(), config=config.to_dict(),
                                    sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                    same_physics_reused_from=run.attrs.get('same_physics_reused_from'),
                                    ambient_aod=run.attrs['aerosol_aod550_ambient']))
                completed += 1
                if progress:
                    progress(f'{completed}/64 {scenario["key"]}: dry AOD={dust:g}, fog AOD={fog:g}'
                             + (' (verified cache)' if reused else ' (same aerosol-free atmosphere)'
                                if run.attrs.get('same_physics_reused_from') else ' (computed)'))
                if i == scenario['dust'] and j == scenario['fog']:
                    export_results(run, directory/scenario['key'], overwrite=True)
                    (directory/scenario['key']/'index.html').replace(directory/scenario['key']/'baseline.html')
            cases.append(fog_cases)
        payload['scenarios'].append({**scenario, 'cases': cases})
    manifest = dict(kind='four-scenario-fog-dust-explorer', version=1, library_version=__version__,
        model_version='2.2', source_hashes=sources, atmospheric_inputs=inputs,
        levels=dict(dry_aod550=DUST_LEVELS, fog_aod550=FOG_LEVELS),
        coordinates=dict(depression=DEPRESSIONS, elevation=ELEVATIONS, azimuth=AZIMUTHS),
        scenarios=SCENARIOS, calculations=records, all_64_cases_complete=True,
        distinct_atmospheres=sum(not r['same_physics_reused_from'] for r in records),
        convergence='Check numerical_sensitivity.json for sampled refinement differences; not an observational accuracy certificate.',
        background='Solar only; no night airglow, starlight, moonlight or city light pollution.')
    (directory/'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    (directory/'data.json').write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(',', ':')), encoding='utf-8')
    render_explorer(directory)
    return directory

def render_explorer(directory):
    """Regenerate the five offline pages from completed data without transport."""
    directory = Path(directory)
    payload_text = (directory/'data.json').read_text(encoding='utf-8').replace('</', r'<\/')
    template = (Path(__file__).parent/'templates'/'atmosphere_explorer.html').read_text(encoding='utf-8')
    validation = '<p>数値収束は未確認です。代表条件の数値感度レポートがあれば、下のリンクから確認できます。</p>'
    report_path = directory/'numerical_sensitivity.json'
    bound = False
    if report_path.exists():
        report = json.loads(report_path.read_text(encoding='utf-8'))
        checks = report.get('scenarios', {})
        bound = len(checks) == 4 and all(
            checks.get(s['key'], {}).get('input_sha256') ==
            hashlib.sha256((directory/s['key']/'sky.nc').read_bytes()).hexdigest() for s in SCENARIOS)
        if bound:
            rows = []
            angles = [-6., -3., 0., 3., 6.]
            for scenario in SCENARIOS:
                angular = checks[scenario['key']]['checks'].get('angular', {})
                values = {row['depression_deg']: row['all']['max_relative_Y']
                          for row in angular.get('per_depression', [])}
                cells = ''.join(f'<td>{values[d]*100:,.2f}%</td>' if d in values else '<td>未実施</td>'
                                for d in angles)
                nodes = angular.get('changed_settings', {}).get('nodes', ['?', '?'])
                rows.append(f'<tr><td>{scenario["title"]}</td><td>{nodes[0]} → {nodes[1]}</td>{cells}</tr>')
            all_checks = [c for s in checks.values() for c in s['checks'].values()]
            completed = [c for c in all_checks if c.get('status') == 'completed']
            passed = sum(c['passes_5percent_Y_and_0_005_xy'] for c in completed)
            incomplete = len(all_checks) - len(completed)
            completion_note = (f'計算失敗または未実施が {incomplete} 件あります。' if incomplete else '')
            heading = ''.join(f'<th>{d:+g}°</th>' for d in angles)
            validation = ('<p><strong>新しい表示範囲での数値感度：</strong>角度積分の点数を増やしたときの、'
                          'サンプル視線における最大輝度差です。列は太陽の伏角（正は地平線の下）です。'
                          '実際の空に対する誤差の上限ではありません。</p>'
                          '<div class="tablescroll"><table><thead><tr><th>条件</th><th>角度点数</th>' + heading +
                          '</tr></thead><tbody>' + ''.join(rows) + '</tbody></table></div>'
                          f'<p>完了した検証 {len(completed)} 件中 {passed} 件で、最大輝度差5%以内・色度差0.005以内を満たしました。'
                          + completion_note +
                          '高湿度の初期条件では、波長・鉛直・水平・反復・Mie係数・Delta-M次数・一回散乱の独立格子も検査対象に含めています。'
                          '全64条件・全視線の収束を保証する検査ではありません。'
                          '詳しい差と未達項目は数値感度レポートで確認できます。</p>')
    for number, name in [(-1, None), *enumerate(s['key'] for s in SCENARIOS)]:
        target = directory if name is None else directory/name
        target.mkdir(parents=True, exist_ok=True)
        page = template.replace('__DATA_JSON__', payload_text).replace('__INITIAL_SCENARIO__', str(max(number, 0)))
        page = page.replace('__ROOT_PREFIX__', '' if name is None else '../')
        page = page.replace('__VALIDATION_SUMMARY__', validation)
        (target/'index.html').write_text(page, encoding='utf-8')
    manifest_path = directory/'manifest.json'
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        manifest['presentation'] = dict(template_sha256=hashlib.sha256(template.encode()).hexdigest(),
            renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            validation_bound_to_current_baselines=bound)
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
