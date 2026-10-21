"""akutou-sky: sunrise, sunset and twilight sky colours for the akutou planetarium projector.

Importing the package is light; SASKTRAN2 loads on first use of the solver.
"""
from importlib import import_module

from ._version import __version__
from .numerics import Numerics, PRESETS
from .scene import SCENES, Scene

_LAZY = {
    'render': 'solver', 'add_color': 'color', 'band_xyz': 'color', 'BANDS': 'color',
    'led_table': 'led', 'load_primaries': 'led', 'write_header': 'led',
    'save': 'io', 'load': 'io', 'export': 'io', 'build_explorer': 'explorer',
    'convergence': 'check', 'illuminance': 'check',
}
__all__ = ['__version__', 'Scene', 'SCENES', 'Numerics', 'PRESETS', *_LAZY]


def __getattr__(name):
    if name not in _LAZY:
        raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
    value = getattr(import_module(f'.{_LAZY[name]}', __name__), name)
    globals()[name] = value
    return value
