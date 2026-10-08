"""Sunrise, sunset and twilight spectra, photometry and explicit result export.

Importing this package does not load SASKTRAN2, download data, change environment
variables or create output directories. Heavy modules are loaded when needed.
"""
from importlib import import_module

from ._version import __version__
from .config import ModelConfig
from .api import simulate, compute_radiance, ExperimentalPolarizationWarning

_LAZY_EXPORTS = {
    "add_colorimetry": ("colorimetry", "add_colorimetry"),
    "xyz_from_radiance": ("colorimetry", "xyz_from_radiance"),
    "display_values": ("colorimetry", "display_values"),
    "band_mean": ("colorimetry", "band_mean"),
    "save_dataset": ("io", "save_dataset"),
    "load_dataset": ("io", "load_dataset"),
    "export_results": ("io", "export_results"),
    "compare_runs": ("validation", "comparison"),
    "validate_sensitivity": ("validation", "validate_sensitivity"),
    "plot_profiles": ("plotting", "plot_profiles"),
}
__all__ = ["__version__", "ModelConfig", "simulate", "compute_radiance",
           "ExperimentalPolarizationWarning", *_LAZY_EXPORTS]

def __getattr__(name):
    if name not in _LAZY_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module, attribute = _LAZY_EXPORTS[name]
    try:
        value = getattr(import_module(f".{module}", __name__), attribute)
    except ModuleNotFoundError as exc:
        if module == 'plotting' and (exc.name or '').startswith('matplotlib'):
            raise ImportError("PNG plotting requires: pip install 'akutou-sky[plot]'") from exc
        raise
    globals()[name] = value
    return value

def __dir__():
    return sorted(set(globals()) | set(__all__))
