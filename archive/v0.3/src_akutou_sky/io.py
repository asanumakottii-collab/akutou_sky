"""Explicit persistence: NetCDF alone, or a complete HTML/CSV/NetCDF bundle."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import xarray as xr

EXPORT_FILES = ('sky.nc', 'manifest.json', 'colors.csv', 'spectra.csv', 'bands.csv', 'summary.json', 'index.html')

def load_dataset(path: str | Path) -> xr.Dataset:
    """Load into memory and close the input file before returning."""
    import xarray as xr
    with xr.open_dataset(path) as source:
        return source.load()

def save_dataset(dataset: xr.Dataset, path: str | Path, *, overwrite: bool = False) -> Path:
    """Save compressed NetCDF atomically; refuse overwrite by default."""
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f'{path} exists; pass overwrite=True to replace it')
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f'.{path.stem}-', suffix='.nc', delete=False) as stream:
        temporary = Path(stream.name)
    try:
        encoding = {name: {'zlib': True, 'complevel': 4} for name in dataset.data_vars}
        dataset.to_netcdf(temporary, engine='netcdf4', encoding=encoding)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return path

def export_results(dataset: xr.Dataset, directory: str | Path, *, overwrite: bool = False) -> Path:
    """Write sky.nc, manifest.json, CSVs and a self-contained offline HTML viewer.

    Raw transport datasets receive colour conversion before export. Existing
    result files are refused unless overwrite=True. Validation reports are only
    displayed when their recorded hash matches the exported NetCDF file.
    """
    from .colorimetry import add_colorimetry
    from .export import export
    from .provenance import manifest
    directory = Path(directory)
    occupied = [directory/name for name in EXPORT_FILES if (directory/name).exists()]
    if occupied and not overwrite:
        raise FileExistsError(f'Result files already exist in {directory}; pass overwrite=True')
    if 'total_XYZ' not in dataset:
        dataset = add_colorimetry(dataset)
    directory.mkdir(parents=True, exist_ok=True)
    raw = save_dataset(dataset, directory/'sky.nc', overwrite=overwrite)
    (directory/'manifest.json').write_text(json.dumps(manifest(dataset, raw), ensure_ascii=False, indent=2)+'\n')
    (directory/'requested_config.json').write_text(json.dumps(json.loads(dataset.attrs['config_json']), indent=2)+'\n')
    export(dataset, directory)
    return directory
