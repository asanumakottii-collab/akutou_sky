"""Spherical twilight radiative transfer with SASKTRAN2 (Geometry2D, successive orders).

    total = single scattering (fine altitude grid, full phase function)
          + diffuse          (coarse grid: successive orders minus their own single scattering)

Scalar runs solve the diffuse field with a consistently truncated delta-M phase
function (Wiscombe 1977) and use TMS single scattering (Nakajima & Tanaka 1988);
this is required for forward-peaked fog droplets. Polarized runs (num_stokes=3,
no fog) solve the full phase function without truncation.
"""
import importlib.metadata
import json
import time
import numpy as np
import sasktran2 as sk
import xarray as xr

from ._version import __version__
from .data import input_paths, resolve_cache_dir
from .numerics import (HORIZONTAL_GROWTH, atmosphere_grid, check_workload, correction_numerics, depression_groups, get_numerics, horizontal_grid,
                       single_scatter_grid, source_grid)
from .optics import Optics
from .scene import get_scene

MODEL_VERSION = '3.0'
EARTH_RADIUS_M = 6371000.
REFRACTION_WAVELENGTH_NM = 550.
DEFAULT_ELEVATIONS = [0., .5, 1., 2., 3., 5., 8., 12., 20., 30., 45., 60., 90.]


def spectral_grid(step):
    edges = np.linspace(380., 780., round(400 / step) + 1)
    return (edges[:-1] + edges[1:]) / 2, edges


def checked_angles(values, name, lower, upper):
    arr = np.atleast_1d(np.asarray(values, dtype=float))
    if arr.ndim != 1 or not len(arr) or not np.isfinite(arr).all():
        raise ValueError(f'{name} must be a nonempty one-dimensional list of numbers')
    if (arr < lower).any() or (arr > upper).any() or len(np.unique(arr)) != len(arr):
        raise ValueError(f'{name} must be unique and within [{lower}, {upper}]')
    return arr


def _config(numerics, scene, *, multiple):
    cfg = sk.Config()
    cfg.num_threads = numerics.threads()
    cfg.num_stokes = numerics.num_stokes
    cfg.num_singlescatter_moments = numerics.moments(scene)
    cfg.solar_refraction = numerics.solar_refraction
    cfg.wavelength_batch_size = 1
    cfg.delta_m_scaling = True  # TMS single scattering on delta-M scaled optical depth
    cfg.num_streams = numerics.delta_m(scene)
    cfg.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    if multiple:
        cfg.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
        cfg.num_successive_orders_incoming = numerics.incoming_nodes
        cfg.num_successive_orders_outgoing = numerics.outgoing_nodes
        cfg.successive_orders_reduced_horizon_quadrature = True
        cfg.successive_orders_altitude_grid_m = source_grid(numerics)
        cfg.num_successive_orders_iterations = 150
        cfg.successive_orders_relative_tolerance = 1e-10
        cfg.successive_orders_absolute_tolerance = 1e-24
    return cfg


def _geometry(height, cfg):
    geo = sk.Geometry2D(cos_sza=0., solar_azimuth=0., earth_radius_m=EARTH_RADIUS_M,
                        altitude_grid_m=height, horizontal_angle_grid_radians=np.deg2rad([-45., 45.]))
    probe = sk.Atmosphere(geo, cfg, wavelengths_nm=np.array([REFRACTION_WAVELENGTH_NM]), calculate_derivatives=False)
    sk.climatology.us76.add_us76_standard_atmosphere(probe)
    geo.refractive_index = sk.optical.refraction.ciddor_index_of_refraction(
        probe.temperature_k, probe.pressure_pa, 0., 420., REFRACTION_WAVELENGTH_NM)
    return geo


def _view(ds, az, es, scene):
    view = sk.ViewingGeometry()
    for d in ds:
        for a in az:
            for e in es:
                view.add_ray(sk.SolarAnglesObserverLocation(
                    -np.sin(np.deg2rad(d)), np.deg2rad(a), np.sin(np.deg2rad(e)), scene.observer_altitude_m))
    return view


def render(depressions, scene='standard', numerics='final', elevations=None, azimuths=(0.,), *,
           cache_dir=None, allow_download=True, progress=None, _native=False):
    """Spectral sky radiance (W m-2 sr-1 nm-1) as an xarray Dataset.

    depressions: geometric depression of the solar centre in degrees (positive
    = below the horizon), -18..24. elevations: viewing elevation 0..90 deg.
    azimuths: relative to the Sun, -180..180 deg (0 = sunward). Order is kept.
    """
    scene, numerics = get_scene(scene), get_numerics(numerics)
    ds = checked_angles(depressions, 'depressions', -18, 24)
    es = checked_angles(DEFAULT_ELEVATIONS if elevations is None else elevations, 'elevations', 0, 90)
    az = checked_angles(azimuths, 'azimuths', -180, 180)
    check_workload(numerics, scene, ds, native=_native)
    split = numerics.num_stokes == 1 and not _native
    wavelengths, edges = spectral_grid(numerics.wavelength_step_nm)
    cache_dir = resolve_cache_dir(cache_dir)
    paths = input_paths(cache_dir, allow_download=allow_download)
    optics = Optics(scene, wavelengths, edges, paths, cache_dir, numerics.moments(scene), numerics.threads())

    coarse, fine = atmosphere_grid(numerics), single_scatter_grid(numerics)
    single = _config(numerics, scene, multiple=False)
    geo_coarse, geo_fine = _geometry(coarse, single), _geometry(fine, single)
    groups = depression_groups(numerics, ds, scene)
    if progress:
        progress(f'{len(ds) * len(az) * len(es)} rays, {len(wavelengths)} wavelengths, {len(groups)} depression '
                 f'group(s); diffuse field {len(coarse)} levels x {len(source_grid(numerics))} sources x '
                 f'{numerics.outgoing_nodes} out / {numerics.incoming_nodes} in directions; single scattering {len(fine)} levels')
    shape = (len(ds), len(az), len(es), len(wavelengths), numerics.num_stokes)
    total, single_part = np.empty(shape), np.empty(shape)
    start = time.monotonic()

    def as_array(radiance, n):
        return radiance.transpose('los', 'wavelength', 'stokes').values.reshape((n,) + shape[1:3] + (1, shape[4]))

    # Observer single scattering: all rays at once on the fine grid (no source field, cheap).
    engine = sk.Engine(single, geo_fine, _view(ds, az, es, scene))
    for i, w in enumerate(wavelengths):
        atmo = sk.Atmosphere(geo_fine, single, wavelengths_nm=np.array([w]), calculate_derivatives=False)
        ozone = optics.add_to(atmo, fine, i, i + 1)
        single_part[..., i:i + 1, :] = as_array(engine.calculate_radiance(atmo).radiance, len(ds))
    total[:] = single_part

    # Diffuse field: coarse grid, one depression group at a time to bound memory.
    for g, group in enumerate(groups):
        view = _view(ds[group], az, es, scene)
        columns = horizontal_grid(numerics, ds[group])
        so = _config(numerics, scene, multiple=True)
        so.num_sza = len(columns)
        so.successive_orders_horizontal_angle_grid_radians = np.deg2rad(columns)
        single_coarse = sk.Engine(single, geo_coarse, view)
        if split:
            from .delta_m import DeltaMTransport
            diffuse_engine = DeltaMTransport(so, geo_coarse, view, numerics.delta_m(scene))
        else:  # polarized: backend-native delta-M for the diffuse field
            diffuse_engine = sk.Engine(so, geo_coarse, view)
        for i, w in enumerate(wavelengths):
            atmo = sk.Atmosphere(geo_coarse, single if split else so, wavelengths_nm=np.array([w]),
                                 calculate_derivatives=False)
            optics.add_to(atmo, coarse, i, i + 1)
            if split:  # the single engine runs first: it fills the delta-M scaled storage the correction uses
                single_coarse.calculate_radiance(atmo)
                diffuse = diffuse_engine.correction(atmo)
            else:
                diffuse = diffuse_engine.calculate_radiance(atmo).radiance - single_coarse.calculate_radiance(atmo).radiance
            total[group, ..., i:i + 1, :] += as_array(diffuse, len(group))
            if progress:
                progress(f'  group {g + 1}/{len(groups)} ({len(columns)} columns): {w:.0f} nm, '
                         f'{time.monotonic() - start:.0f} s')
        del diffuse_engine, single_coarse

    polarization = None
    if numerics.num_stokes == 1 and numerics.polarization_correction:
        # Single scattering of unpolarized sunlight is exact in scalar form; polarization changes the
        # diffuse field. Its effect (vector - scalar, both backend-native delta-M on one coarse grid,
        # so their discretization errors cancel) is added to the fine scalar solution.
        coarse_numerics = correction_numerics(numerics, scene)
        pair = [render(ds, scene, coarse_numerics.replace(num_stokes=n), es, az, cache_dir=cache_dir,
                       allow_download=allow_download, _native=True).radiance.sel(stokes='I').values
                for n in (3, 1)]
        polarization = pair[0] - pair[1]
        total[..., 0] += polarization
        if progress:
            progress(f'  polarization correction: {time.monotonic() - start:.0f} s')

    # Finite phase expansions can leave tiny signed components; never clip, but refuse a negative total.
    if not np.isfinite(total).all() or (total[..., 0] < 0).any():
        raise ValueError('Nonfinite or negative total radiance; refine numerics')
    if numerics.num_stokes == 3 and (np.linalg.norm(total[..., 1:], axis=-1) > total[..., 0] * 1.001 + 1e-20).any():
        raise ValueError('Degree of polarization > 1; refine numerics')

    dims = ('depression', 'azimuth', 'elevation', 'wavelength', 'stokes')
    out = xr.Dataset(
        {'radiance': (dims, total), 'single_scatter_radiance': (dims, single_part),
         **({} if polarization is None else {'polarization_correction': (dims[:-1], polarization)}),
         'wavelength_bounds': (('wavelength', 'bound'), np.c_[edges[:-1], edges[1:]]),
         'solar_irradiance': ('wavelength', optics.irradiance),
         'ozone_number_density': ('altitude', ozone),
         **{f'{k}_extinction_550': ('altitude', v) for k, v in optics.profiles(fine).items()}},
        coords={'depression': ds, 'azimuth': az, 'elevation': es, 'wavelength': wavelengths,
                'stokes': ['I', 'Q', 'U'][:numerics.num_stokes], 'altitude': fine},
        attrs={'library_version': __version__, 'model_version': MODEL_VERSION,
               'sasktran2_version': importlib.metadata.version('sasktran2'),
               'scene_json': json.dumps(scene.to_dict(), sort_keys=True),
               'numerics_json': json.dumps(numerics.to_dict(), sort_keys=True),
               'numerics_resolved_json': json.dumps(numerics.resolved(scene).to_dict(), sort_keys=True),
               'horizontal_growth': HORIZONTAL_GROWTH,
               'aerosol_aod550_ambient': optics.ambient_aod550, 'aerosol_growth_factor': optics.growth,
               'aerosol_wet_median_radius_nm': optics.wet_radius_nm,
               'solver': 'SASKTRAN2 Geometry2D successive orders; ' + (
                   'scalar split delta-M diffuse + TMS single scatter' if split else 'backend-native delta-M')
                   + (' + coarse-grid polarization correction' if polarization is not None else ''),
               'angles': 'geometric; depression of solar centre (+ below horizon); azimuth 0 = sunward',
               'excluded': 'direct solar disk, night-sky background, clouds, terrain, LOS and diffuse refraction'
                           + ('' if numerics.solar_refraction else ', solar refraction'),
               'elapsed_seconds': time.monotonic() - start})
    units = dict(radiance='W m-2 sr-1 nm-1', single_scatter_radiance='W m-2 sr-1 nm-1',
                 **({} if polarization is None else {'polarization_correction': 'W m-2 sr-1 nm-1'}),
                 solar_irradiance='W m-2 nm-1', ozone_number_density='m-3', wavelength='nm',
                 wavelength_bounds='nm', altitude='m', depression='degree', azimuth='degree', elevation='degree')
    for name, unit in units.items():
        out[name].attrs['units'] = unit
    for name in out.data_vars:
        if name.endswith('_extinction_550'):
            out[name].attrs['units'] = 'm-1'
    return out
