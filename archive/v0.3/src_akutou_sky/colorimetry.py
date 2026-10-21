"""Absolute photometry and explicit display previews, separate from transport."""
import csv
from pathlib import Path
import numpy as np
import colour
import xarray as xr
from .config import integrate_linear

def xyz_from_radiance(radiance, wavelengths, bounds=None):
    rad, w = np.asarray(radiance), np.asarray(wavelengths)
    if rad.shape[-1] != len(w) or not np.isfinite(rad).all() or (rad < 0).any():
        raise ValueError('Radiance must be finite, nonnegative, and match the wavelength grid')
    cmf = colour.MSDS_CMFS['CIE 1931 2 Degree Standard Observer']
    if bounds is None:  # legacy point-sampled data
        cm = np.stack([np.interp(w, cmf.wavelengths, cmf.values[:, i]) for i in range(3)], axis=-1)
        return 683. * np.trapezoid(rad[..., None] * cm, w, axis=-2)
    bounds = np.asarray(bounds)
    if bounds.shape != (len(w), 2) or not np.allclose(bounds[1:, 0], bounds[:-1, 1]):
        raise ValueError('Spectral bounds must form contiguous bins')
    edges = np.r_[bounds[:, 0], bounds[-1, 1]]
    weights = np.stack([integrate_linear(cmf.wavelengths, cmf.values[:, i], edges)
                        for i in range(3)], axis=-1)
    return 683. * np.einsum('...w,wc->...c', rad, weights)

def display_values(xyz):
    xyz = np.asarray(xyz)
    total = xyz.sum(axis=-1, keepdims=True)
    xy = np.divide(xyz[..., :2], total, out=np.full_like(xyz[..., :2], np.nan), where=total > 0)
    rgb = colour.XYZ_to_sRGB(xyz, apply_cctf_encoding=False)
    clipped = (rgb < -1e-12 * np.maximum(1, np.max(np.abs(rgb), axis=-1, keepdims=True))).any(axis=-1)
    positive = np.maximum(rgb, 0)
    peak = positive.max(axis=-1, keepdims=True)
    normalized = np.divide(positive, peak, out=np.zeros_like(positive), where=peak > 0)
    display = np.clip(colour.cctf_encoding(normalized, function='sRGB'), 0, 1)
    return xy, display, clipped

def background_spectrum(path, bounds):
    """At-observer isotropic spectral radiance; never infer a night-sky RGB floor."""
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream))
    try:
        w = np.array([float(r['wavelength_nm']) for r in rows])
        rad = np.array([float(r['radiance_W_m2_sr_nm']) for r in rows])
    except (KeyError, ValueError) as exc:
        raise ValueError('Background CSV requires wavelength_nm,radiance_W_m2_sr_nm') from exc
    if (rad < 0).any():
        raise ValueError('Background radiance must be nonnegative')
    edges = np.r_[bounds[:, 0], bounds[-1, 1]]
    return integrate_linear(w, rad, edges) / np.diff(edges)

def add_colorimetry(ds, background=None):
    ds = ds.copy()
    intensity = ds.solar_stokes_radiance.sel(stokes='I', drop=True)
    dims = ('depression', 'azimuth', 'elevation')
    bounds = ds.wavelength_bounds.values
    bg = np.zeros(ds.sizes['wavelength']) if background is None else background_spectrum(background, bounds)
    ds['background_radiance'] = ('wavelength', bg)
    ds.background_radiance.attrs['units'] = 'W m-2 sr-1 nm-1'
    ds['total_radiance'] = intensity + ds.background_radiance
    ds.total_radiance.attrs['units'] = 'W m-2 sr-1 nm-1'
    ds = ds.assign_coords(tristimulus=['X', 'Y', 'Z'], chromaticity=['x', 'y'], channel=['R', 'G', 'B'])
    for component, radiance in [('solar', intensity), ('total', ds.total_radiance)]:
        xyz = xyz_from_radiance(radiance.values, ds.wavelength.values, bounds)
        xy, rgb, clipped = display_values(xyz)
        ds[f'{component}_XYZ'] = (dims + ('tristimulus',), xyz)
        ds[f'{component}_XYZ'].attrs['units'] = 'cd m-2 (Y); matching 683-scaled X and Z'
        ds[f'{component}_xy'] = (dims + ('chromaticity',), xy)
        ds[f'{component}_preview_sRGB'] = (dims + ('channel',), rgb)
        ds[f'{component}_preview_sRGB'].attrs['normalization'] = 'each ray max linear RGB = 1; no brightness information'
        ds[f'{component}_gamut_clipped'] = (dims, clipped.astype('int8'))
    ds['near_horizon'] = ('elevation', (ds.elevation.values < 5).astype('int8'))
    ds['night_background_important'] = ('depression', (ds.depression.values >= 12).astype('int8'))
    if ds.sizes['stokes'] == 3:
        pol = np.hypot(ds.solar_stokes_radiance.sel(stokes='Q'), ds.solar_stokes_radiance.sel(stokes='U'))
        ds['solar_degree_linear_polarization'] = xr.where(intensity > 0, pol / intensity, np.nan)
    ds.attrs['background'] = 'none' if background is None else 'user supplied, at-observer, isotropic, unpolarized; no extra extinction applied'
    ds.attrs['preview_note'] = 'Display sRGB is not LED PWM, visual adaptation, or predicted dark-adapted appearance.'
    return ds

def band_mean(elevations, xyz, lower, upper):
    """Solid-angle weighted mean, inserting the exact band boundaries."""
    e = np.asarray(elevations)
    order = np.argsort(e)
    e, xyz = e[order], np.asarray(xyz)[..., order, :]
    if not e[0] <= lower < upper <= e[-1]:
        raise ValueError('Band lies outside the computed elevation range')
    grid = np.unique(np.r_[lower, e[(e > lower) & (e < upper)], upper])
    flat = xyz.reshape(-1, len(e), 3)
    interp = np.array([[np.interp(grid, e, row[:, c]) for c in range(3)] for row in flat]).transpose(0, 2, 1)
    radians = np.deg2rad(grid)
    mean = np.trapezoid(interp * np.cos(radians)[None, :, None], radians, axis=1)
    mean /= np.trapezoid(np.cos(radians), radians)
    return mean.reshape(xyz.shape[:-2] + (3,))
