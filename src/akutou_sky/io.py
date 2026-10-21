"""NetCDF persistence and plain CSV export."""
import json
from pathlib import Path
import tempfile
import numpy as np
import xarray as xr
from .color import BANDS, add_color, band_xyz, chromaticity


def save(ds, path):
    """Write compressed NetCDF atomically (replaces an existing file)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f'.{path.stem}-', suffix='.nc', delete=False) as f:
        temporary = Path(f.name)
    try:
        ds.to_netcdf(temporary, engine='netcdf4',
                     encoding={name: {'zlib': True, 'complevel': 4} for name in ds.data_vars})
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def load(path):
    with xr.open_dataset(path) as source:
        return source.load()


def export(ds, directory):
    """sky.nc, colors.csv (every ray), bands.csv (RevC bands), scene.json, numerics.json."""
    directory = Path(directory)
    if 'XYZ' not in ds:
        ds = add_color(ds)
    save(ds, directory / 'sky.nc')
    (directory / 'scene.json').write_text(json.dumps(json.loads(ds.attrs['scene_json']), indent=2) + '\n')
    (directory / 'numerics.json').write_text(json.dumps(json.loads(ds.attrs['numerics_json']), indent=2) + '\n')
    rows = ['depression_deg,azimuth_deg,elevation_deg,X,Y_cd_m2,Z,x,y']
    for i, d in enumerate(ds.depression.values):
        for j, a in enumerate(ds.azimuth.values):
            for k, e in enumerate(ds.elevation.values):
                rows.append(','.join(f'{v:.6g}' for v in (d, a, e, *ds.XYZ.values[i, j, k], *ds.xy.values[i, j, k])))
    (directory / 'colors.csv').write_text('\n'.join(rows) + '\n')
    if ds.elevation.values.min() <= BANDS[0][0] and ds.elevation.values.max() >= BANDS[-1][1]:
        xyz = band_xyz(ds)
        xy = chromaticity(xyz)
        rows = ['depression_deg,azimuth_deg,band,elevation_min_deg,elevation_max_deg,X,Y_cd_m2,Z,x,y']
        for i, d in enumerate(ds.depression.values):
            for j, a in enumerate(ds.azimuth.values):
                for b, (lo, hi) in enumerate(BANDS):
                    rows.append(','.join(f'{v:.6g}' for v in (d, a, b + 1, lo, hi, *xyz[i, j, b], *xy[i, j, b])))
        (directory / 'bands.csv').write_text('\n'.join(rows) + '\n')
    return directory


def same_run(ds, scene, numerics, depressions, elevations, azimuths):
    """True if a saved dataset was computed with exactly these inputs and axes."""
    from .numerics import Numerics
    try:
        saved = Numerics.from_dict(json.loads(ds.attrs['numerics_json']))
        return (json.loads(ds.attrs['scene_json']) == scene.to_dict()
                and saved.resolved(scene) == numerics.resolved(scene)
                and all(np.array_equal(ds[k].values, np.asarray(v, dtype=float)) for k, v in
                        (('depression', depressions), ('elevation', elevations), ('azimuth', azimuths))))
    except (KeyError, ValueError):
        return False
