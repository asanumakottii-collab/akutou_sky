"""Stable entry points; no implicit result files or process-wide configuration."""
from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import TYPE_CHECKING
import warnings

from .config import ModelConfig

if TYPE_CHECKING:
    import xarray as xr

class ExperimentalPolarizationWarning(UserWarning):
    """The polarized solver has unresolved coarse-grid symmetry failures."""

def compute_radiance(
    depressions: Sequence[float], *, config: ModelConfig | None = None,
    elevations: Sequence[float] | None = None, azimuths: Sequence[float] | None = None,
    cache_dir: str | Path | None = None, allow_download: bool = True,
    progress: Callable[[str], None] | None = None,
) -> xr.Dataset:
    """Compute solar-only I/Q/U in W m-2 sr-1 nm-1 without colour conversion.

    Depression, elevation and relative azimuth are geometric degrees. Positive
    depression means the solar centre is below the horizon; azimuth 0 is sunward.
    The ordering of the supplied coordinates is preserved.
    """
    configuration = (config or ModelConfig()).validate()
    if configuration.num_stokes == 3:
        warnings.warn(
            "Polarization is experimental: a coarse-grid mirror-symmetry check failed; "
            "quantitative accuracy is not established.",
            ExperimentalPolarizationWarning, stacklevel=2,
        )
    from .model import compute_model
    return compute_model(depressions, configuration, elevations, azimuths, progress,
                         cache_dir=cache_dir, allow_download=allow_download)

def simulate(
    depressions: Sequence[float], *, config: ModelConfig | None = None,
    elevations: Sequence[float] | None = None, azimuths: Sequence[float] | None = None,
    background: str | Path | None = None, cache_dir: str | Path | None = None,
    allow_download: bool = True, progress: Callable[[str], None] | None = None,
) -> xr.Dataset:
    """Compute spectra, absolute XYZ/xy/Y and normalized display sRGB in memory.

    Background is an observed, isotropic, unpolarized at-observer spectrum CSV.
    Use save_dataset or export_results to write results. A single simulation does
    not establish numerical convergence or accuracy against observations.
    """
    configuration = (config or ModelConfig()).validate()
    from .colorimetry import add_colorimetry, background_spectrum
    if background is not None:
        import numpy as np
        from .config import spectral_grid
        _, edges = spectral_grid(configuration.wavelength_step_nm)
        background_spectrum(background, np.c_[edges[:-1], edges[1:]])
    dataset = compute_radiance(depressions, config=configuration, elevations=elevations,
                              azimuths=azimuths, cache_dir=cache_dir,
                              allow_download=allow_download, progress=progress)
    result = add_colorimetry(dataset, background)
    if background is not None:
        from .provenance import sha256
        result.attrs['background_path'] = str(Path(background).resolve())
        result.attrs['background_sha256'] = sha256(background)
    return result
