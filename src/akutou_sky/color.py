"""Photometry and colour: spectral radiance -> CIE 1931 XYZ (Y in cd/m2), xy, bands."""
import numpy as np
import colour
from .optics import integrate_linear

# Elevation bands (deg) of the four RevC LED series, bottom to top.
BANDS = ((0., 5.), (5., 15.), (15., 35.), (35., 60.))
SRGB_FROM_XYZ = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])


def cmf_weights(bounds):
    """CIE 1931 2-degree colour-matching functions integrated over each wavelength bin."""
    cmf = colour.MSDS_CMFS['CIE 1931 2 Degree Standard Observer']
    edges = np.r_[bounds[:, 0], bounds[-1, 1]]
    return np.stack([integrate_linear(cmf.wavelengths, cmf.values[:, i], edges) for i in range(3)], axis=-1)


def xyz_from_radiance(radiance, bounds):
    """X, Y, Z = 683 * sum_bins L * integral(cmf); Y is luminance in cd/m2."""
    radiance = np.asarray(radiance)
    if not np.isfinite(radiance).all() or (radiance < 0).any():
        raise ValueError('Radiance must be finite and nonnegative')
    return 683. * np.einsum('...w,wc->...c', radiance, cmf_weights(np.asarray(bounds)))


def chromaticity(xyz):
    xyz = np.asarray(xyz, dtype=float)
    total = xyz.sum(axis=-1, keepdims=True)
    return np.divide(xyz[..., :2], total, out=np.full_like(xyz[..., :2], np.nan), where=total > 0)


def preview_srgb(xyz):
    """Display-only sRGB with each colour's brightest channel = 1 (no luminance information)."""
    rgb = np.maximum(np.asarray(xyz) @ SRGB_FROM_XYZ.T, 0)
    peak = rgb.max(axis=-1, keepdims=True)
    rgb = np.divide(rgb, peak, out=np.zeros_like(rgb), where=peak > 0)
    return np.clip(np.where(rgb <= .0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - .055), 0, 1)


def add_color(ds):
    """Return ds with XYZ, xy and a display preview for every viewing ray."""
    ds = ds.copy()
    xyz = xyz_from_radiance(ds.radiance.sel(stokes='I').values, ds.wavelength_bounds.values)
    dims = ('depression', 'azimuth', 'elevation')
    ds = ds.assign_coords(tristimulus=['X', 'Y', 'Z'], chromaticity=['x', 'y'], channel=['R', 'G', 'B'])
    ds['XYZ'] = (dims + ('tristimulus',), xyz)
    ds.XYZ.attrs['units'] = 'cd m-2 (Y); X and Z on the same 683 lm/W scale'
    ds['xy'] = (dims + ('chromaticity',), chromaticity(xyz))
    ds['preview_sRGB'] = (dims + ('channel',), preview_srgb(xyz))
    ds.preview_sRGB.attrs['note'] = 'display only; normalized per ray; not LED drive values'
    if ds.sizes['stokes'] == 3:
        i = ds.radiance.sel(stokes='I')
        ds['degree_of_polarization'] = np.hypot(ds.radiance.sel(stokes='Q'), ds.radiance.sel(stokes='U')) / i
    return ds


def band_mean(elevations, xyz, lower, upper):
    """Solid-angle (cos h) weighted mean of XYZ over lower..upper deg of elevation.

    Elevations are the last-but-one axis of xyz; values between samples are
    interpolated linearly in XYZ (light amounts), never in RGB.
    """
    e = np.asarray(elevations, dtype=float)
    order = np.argsort(e)
    e, xyz = e[order], np.asarray(xyz)[..., order, :]
    if not e[0] <= lower < upper <= e[-1]:
        raise ValueError(f'Band {lower}..{upper} deg is outside the computed elevations')
    grid = np.unique(np.r_[lower, e[(e > lower) & (e < upper)], upper])
    flat = xyz.reshape(-1, len(e), 3)
    values = np.stack([[np.interp(grid, e, row[:, c]) for c in range(3)] for row in flat]).transpose(0, 2, 1)
    h = np.deg2rad(grid)
    mean = np.trapezoid(values * np.cos(h)[None, :, None], h, axis=1) / np.trapezoid(np.cos(h), h)
    return mean.reshape(xyz.shape[:-2] + (3,))


def band_xyz(ds, bands=BANDS):
    """XYZ per (depression, azimuth, band) from a coloured dataset."""
    return np.stack([band_mean(ds.elevation.values, ds.XYZ.values, lo, hi) for lo, hi in bands], axis=-2)
