"""Result provenance independent of the checkout's directory structure."""
import hashlib
import importlib.metadata
import json
from datetime import datetime, timezone
from pathlib import Path

from ._version import __version__

def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def manifest(dataset, path=None):
    packages = {}
    for name in ['sasktran2', 'numpy', 'scipy', 'colour-science', 'xarray', 'netCDF4']:
        packages[name] = importlib.metadata.version(name)
    package_dir = Path(__file__).parent
    result = {
        'library_version': __version__, 'model_version': dataset.attrs.get('model_version'),
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'config': json.loads(dataset.attrs['config_json']),
        'coordinates': {k: dataset[k].values.tolist() for k in ['depression', 'azimuth', 'elevation']},
        'packages': packages,
        'source_sha256': {p.name: sha256(p) for p in sorted(package_dir.glob('*.py'))},
        'physical_inputs': json.loads(dataset.attrs.get('physical_inputs_json', '[]')),
        'convergence': dataset.attrs.get('convergence', 'not established'),
        'background': None,
    }
    if 'background_path' in dataset.attrs:
        result['background'] = {'path': dataset.attrs['background_path'],
                                'sha256': dataset.attrs['background_sha256']}
    if path is not None:
        result['sky_sha256'] = sha256(path)
    return result
