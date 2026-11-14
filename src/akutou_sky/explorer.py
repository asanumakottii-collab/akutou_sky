"""Offline comparison page: four skies x dust x stratospheric aerosol, every slider stop computed.

The stratospheric slider controls the purple light of twilight. Each case is
cached as calculations/<scene>/dust<i>_strat<j>.nc and reused only when scene,
numerics and axes are identical. Cases that differ only in the properties of
an absent aerosol layer share one calculation.
"""
import json
from pathlib import Path
import numpy as np

from ._version import __version__
from .color import BANDS, add_color, band_xyz
from .io import load, same_run, save
from .numerics import get_numerics
from .scene import SCENES, Scene
from .solver import MODEL_VERSION, render

DUST_LEVELS = [0., .03, .10, .30]
# Stratospheric sulfate: none, background (Ugolnikov & Maslov 2017), and the 'eruption' scene.
STRATOSPHERE_LEVELS = [
    dict(aod=0., radius_nm=80., width=1.6),
    dict(aod=.005, radius_nm=80., width=1.6),
    dict(aod=.02, radius_nm=150., width=1.5),
]
DEPRESSIONS = [i / 2 for i in range(-12, 13)]
ELEVATIONS = [0., .5, 1., 2., 3., 5., 8., 12., 15., 20., 30., 35., 45., 60., 90.]
AZIMUTHS = [0., 90., 180.]
PAGES = [
    dict(key='urban', title='都会', caption='微粒子が多く、光を吸収しやすい条件'),
    dict(key='rural', title='田舎', caption='微粒子が少なく、空気が澄んだ条件'),
    dict(key='dry', title='乾いた日', caption='吸湿成長が小さい、ふつうの晴天'),
    dict(key='eruption', title='噴火後', caption='成層圏の微粒子が増え、紫光が強く出る条件'),
]
_AEROSOL_FIELDS = ('aerosol_scale_height_m', 'median_radius_nm', 'mode_width', 'refractive_index_real',
                   'refractive_index_imag', 'relative_humidity', 'hygroscopicity_kappa')
_STRATOSPHERE_FIELDS = ('stratospheric_center_m', 'stratospheric_width_m', 'stratospheric_median_radius_nm',
                        'stratospheric_mode_width')


def stratosphere(scene, level):
    return scene.replace(stratospheric_aod550=level['aod'], stratospheric_median_radius_nm=level['radius_nm'],
                         stratospheric_mode_width=level['width'])


def _level_of(scene):
    for j, level in enumerate(STRATOSPHERE_LEVELS):
        if stratosphere(scene, level) == scene:
            return j
    raise ValueError('page scene must use one of STRATOSPHERE_LEVELS')


def _physics_key(scene):
    """Scenes without an aerosol layer are physically identical whatever that layer's properties."""
    data = scene.to_dict()
    defaults = Scene().to_dict()
    if scene.aod550 == 0:
        data.update({k: defaults[k] for k in _AEROSOL_FIELDS})
    if scene.stratospheric_aod550 == 0:
        data.update({k: defaults[k] for k in _STRATOSPHERE_FIELDS})
    return json.dumps(data, sort_keys=True)


def build_explorer(directory, numerics='compare', *, cache_dir=None, allow_download=True, progress=None,
                   accuracy_note='晴天の実測照度とは伏角−6〜＋12°で0.96〜1.58倍で一致しています（docs/精度検証.md）。'
                                 '比較ページは計算時間を抑えた compare 設定（波長20 nm幅）です。'):
    numerics = get_numerics(numerics)
    directory = Path(directory)
    computed = {}
    payload = dict(dust=DUST_LEVELS, strat=STRATOSPHERE_LEVELS, depressions=DEPRESSIONS, elevations=ELEVATIONS,
                   azimuths=AZIMUTHS, bands=[list(b) for b in BANDS], numerics=numerics.to_dict(),
                   version=f'{__version__} / model {MODEL_VERSION}', scenarios=[])
    total, done = len(PAGES) * len(DUST_LEVELS) * len(STRATOSPHERE_LEVELS), 0
    for page in PAGES:
        base = SCENES[page['key']]
        cases = []
        for i, dust in enumerate(DUST_LEVELS):
            row = []
            for j, level in enumerate(STRATOSPHERE_LEVELS):
                scene = stratosphere(base, level).replace(aod550=dust)
                relative = f'calculations/{page["key"]}/dust{i}_strat{j}.nc'
                path = directory / relative
                key = _physics_key(scene)
                run = load(path) if path.exists() else None
                if run is not None and not same_run(run, scene, numerics, DEPRESSIONS, ELEVATIONS, AZIMUTHS):
                    run = None
                status = 'cached'
                if run is None and key in computed:
                    run = computed[key].copy()
                    run.attrs['scene_json'] = json.dumps(scene.to_dict(), sort_keys=True)
                    save(run, path)
                    status = 'shared (identical without the absent layer)'
                if run is None:
                    run = add_color(render(DEPRESSIONS, scene, numerics, ELEVATIONS, AZIMUTHS,
                                           cache_dir=cache_dir, allow_download=allow_download))
                    save(run, path)
                    status = f'computed in {run.attrs["elapsed_seconds"]:.0f} s'
                computed.setdefault(key, run)
                done += 1
                if progress:
                    progress(f'{done}/{total} {page["key"]} dust {dust:g} strat {level["aod"]:g}: {status}')
                row.append(dict(XYZ=np.round(run.XYZ.values, 10).tolist(),
                                bands=np.round(band_xyz(run), 10).tolist(),
                                ambient_aod=run.attrs['aerosol_aod550_ambient'],
                                growth=run.attrs['aerosol_growth_factor'],
                                wet_radius=run.attrs['aerosol_wet_median_radius_nm'],
                                path=relative, scene=scene.to_dict()))
            cases.append(row)
        payload['scenarios'].append(dict(page, rh=base.relative_humidity, kappa=base.hygroscopicity_kappa,
                                         index_real=base.refractive_index_real,
                                         index_imag=base.refractive_index_imag,
                                         dust=DUST_LEVELS.index(base.aod550), strat=_level_of(base),
                                         cases=cases))
    data = json.dumps(payload, ensure_ascii=False, allow_nan=False, separators=(',', ':')).replace('</', r'<\/')
    template = (Path(__file__).parent / 'templates' / 'explorer.html').read_text(encoding='utf-8')
    page = (template.replace('__DATA_JSON__', data).replace('__VERSION__', payload['version'])
            .replace('__ACCURACY_NOTE__', accuracy_note))
    (directory / 'index.html').write_text(page, encoding='utf-8')
    return directory / 'index.html'
