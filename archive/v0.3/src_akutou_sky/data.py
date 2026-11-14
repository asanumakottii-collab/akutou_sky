"""Explicit SASKTRAN input cache without changing its global environment."""
from __future__ import annotations

import os
from pathlib import Path
import shutil
import sys
import tempfile
from urllib.request import urlopen

INPUT_KEYS = (
    'climatology/fascode/std.atm',
    'cross_sections/o3/dbm.nc',
    'solar/solar_irradiance_hsrs_2022_11_30_extended.nc',
)
DATABASE_URL = 'https://arg.usask.ca/sasktranfiles/sasktran2_db/v_latest/'

def resolve_cache_dir(cache_dir: str | Path | None = None) -> Path:
    """Resolve a path only; never create directories merely to inspect settings."""
    if cache_dir is not None:
        return Path(cache_dir).expanduser().resolve()
    for variable in ('AKUTOU_SKY_CACHE_DIR', 'SASKTRAN2_DATABASE_ROOT'):
        if os.environ.get(variable):
            return Path(os.environ[variable]).expanduser().resolve()
    if os.environ.get('XDG_CACHE_HOME'):
        base = Path(os.environ['XDG_CACHE_HOME']).expanduser()
    elif sys.platform == 'darwin':
        base = Path.home() / 'Library' / 'Caches'
    elif sys.platform == 'win32':
        base = Path(os.environ.get('LOCALAPPDATA', Path.home() / 'AppData' / 'Local'))
    else:
        base = Path.home() / '.cache'
    return (base / 'akutou-sky' / 'model-data').resolve()

def input_paths(cache_dir: Path, *, allow_download: bool = True) -> tuple[Path, ...]:
    paths = tuple(cache_dir / key for key in INPUT_KEYS)
    missing = [path for path in paths if not path.is_file() or path.stat().st_size == 0]
    if missing and not allow_download:
        raise FileNotFoundError('Missing atmospheric data; downloads disabled: ' + ', '.join(map(str, missing)))
    for path in missing:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=path.parent, suffix='.partial', delete=False) as stream:
                temporary = Path(stream.name)
                key = path.relative_to(cache_dir).as_posix()
                with urlopen(DATABASE_URL + key, timeout=120) as response:
                    shutil.copyfileobj(response, stream)
            if temporary.stat().st_size == 0:
                raise OSError(f'Empty atmospheric database download: {path.name}')
            temporary.replace(path)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    return paths
