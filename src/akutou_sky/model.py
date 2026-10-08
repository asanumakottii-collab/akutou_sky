"""Spherical, polarized twilight radiative transfer; no empirical RGB colouring."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time
import numpy as np
import sasktran2 as sk
import xarray as xr
from scipy.constants import Boltzmann
from .config import ModelConfig, altitude_grid, source_grid, single_scatter_grid, horizontal_grid, spectral_grid, checked_angles, integrate_linear
from .data import input_paths, resolve_cache_dir
from ._version import __version__
from .microphysics import aerosol_water_uptake, WATER_INDEX
MODEL_VERSION = '2.2'
DU_MOLECULES_M2 = 2.687e20

def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def read_fascode(path):
    profiles, key = {}, None
    for line in Path(path).read_text().splitlines():
        line = line.split('!')[0].strip()
        if not line:
            continue
        if line.startswith('*END'):
            break
        if line.startswith('*'):
            key = line[1:].split()[0]
            profiles[key] = []
        elif key is not None:
            profiles[key].extend(float(v) for v in line.replace(',', ' ').split())
    return {key: np.asarray(value) for key, value in profiles.items()}

def normalized_profile(height, shape, column):
    integral = np.trapezoid(shape, height)
    if not np.isfinite(integral) or integral <= 0:
        raise ValueError('Vertical profile has zero or invalid column')
    return shape * column / integral

class BinnedSolar(sk.constituent.base.Constituent):
    def __init__(self, values):
        self.values = values
    def add_to_atmosphere(self, atmo):
        atmo.storage.solar_irradiance[:] = self.values
    def register_derivative(self, atmo, name):
        return {}

def _mie(config, wavelengths, cache_dir, radius=None, width=None, index=None):
    radius = config.median_radius_nm if radius is None else radius
    width = config.mode_width if width is None else width
    index = complex(config.refractive_index_real, -config.refractive_index_imag) if index is None else index
    distribution = sk.mie.distribution.LogNormalDistribution().freeze(median_radius=radius, mode_width=width)
    refractive = sk.mie.refractive.RefractiveIndex(
        lambda w: np.full_like(np.asarray(w, dtype=float), index, dtype=complex),
        f'constant_n{index.real:.12g}_k{-index.imag:.12g}',
    )
    return sk.database.MieDatabase(
        distribution, refractive, wavelengths_nm=np.unique(np.r_[wavelengths, 550.]),
        db_root=cache_dir, max_legendre_moments=config.legendre_moments,
        num_threads=config.num_threads,
    )

def compute_model(depressions, config=None, elevations=None, azimuths=None, progress=None,
                  *, cache_dir=None, allow_download=True, _transport_cache=None, _single_only=False):
    """Compute solar-only spectral I/Q/U as a labeled xarray Dataset.

    Angles are geometric degrees. Positive depression places the Sun below the
    horizon. Cache files are the only implicit writes; result files are explicit.
    The progress callback takes one string. No callback means silent execution.
    """
    config = (config or ModelConfig()).validate()
    ds = checked_angles(depressions, 'depressions', -18, 24)
    es = checked_angles(elevations if elevations is not None else
                        np.r_[np.arange(0, 5, .25), np.arange(5, 91.)], 'elevations', 0, 90)
    az = checked_angles(azimuths if azimuths is not None else [0], 'azimuths', -180, 180)
    work_size = len(ds)*len(az)*len(es)*len(source_grid(config))*config.nsza*config.nodes*config.num_stokes
    if work_size > config.max_view_source_product:
        raise ValueError(f'Viewing/source product {work_size:,} exceeds budget {config.max_view_source_product:,}. '
                         'Use fewer elevations/depressions/azimuths or the overview preset. '
                         'Only raise max_view_source_product with enough RAM; this is a conservative work limit, not a RAM estimate.')
    wavelengths, edges = spectral_grid(config.wavelength_step_nm)
    height = single_scatter_grid(config) if _single_only else altitude_grid(config)
    if _single_only and len(ds)*len(az)*len(es)*len(height)**2 > config.max_view_source_product:
        raise ValueError('Single-scattering ray/altitude work exceeds max_view_source_product; use fewer viewing rays')
    cfg = sk.Config()
    cfg.num_threads = config.num_threads
    cfg.num_stokes = config.num_stokes
    cfg.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
    if config.delta_m_method == 'split' or _single_only:
        cfg.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    cfg.num_successive_orders_incoming = config.nodes
    cfg.num_successive_orders_outgoing = config.nodes
    cfg.successive_orders_reduced_horizon_quadrature = config.horizon_quadrature
    cfg.num_sza = config.nsza
    cfg.successive_orders_altitude_grid_m = source_grid(config)
    cfg.successive_orders_horizontal_angle_grid_radians = np.deg2rad(horizontal_grid(config))
    cfg.num_successive_orders_iterations = config.iterations
    cfg.successive_orders_relative_tolerance = config.relative_tolerance
    cfg.successive_orders_absolute_tolerance = config.absolute_tolerance
    cfg.num_singlescatter_moments = config.legendre_moments
    cfg.delta_m_scaling = config.delta_m_scaling
    if config.delta_m_scaling:
        cfg.num_streams = config.delta_m_order
    cfg.solar_refraction = config.solar_refraction
    cfg.los_refraction = False  # not supported by the Geometry2D solver
    cfg.multiple_scatter_refraction = False  # likewise not supported
    cfg.wavelength_batch_size = 1
    # Explicit, caller-owned batch cache. Only geometry and numerical controls
    # enter the key; each atmosphere and its optical properties are rebuilt.
    from .validation import NUMERICAL_FIELDS
    key_fields = ({'num_threads','legendre_moments','delta_m_scaling','delta_m_order'}
                  if _single_only else NUMERICAL_FIELDS)
    transport_key = (_single_only, tuple(ds), tuple(es), tuple(az), tuple(height),
                     tuple((k, getattr(config, k)) for k in sorted(key_fields)),
                     config.num_stokes, config.observer_altitude_m, config.solar_refraction,
                     config.refraction_wavelength_nm)
    cached_transport = None if _transport_cache is None else _transport_cache.get(transport_key)
    if cached_transport is None and _transport_cache is not None:
        _transport_cache.clear()
    geo = sk.Geometry2D(cos_sza=0., solar_azimuth=0., earth_radius_m=6371000.,
                        altitude_grid_m=height, horizontal_angle_grid_radians=np.deg2rad([-45., 45.]))
    if cached_transport is not None:
        geo = cached_transport[0]
    view = sk.ViewingGeometry()
    for d in ds:
        for a in az:
            for e in es:
                view.add_ray(sk.SolarAnglesObserverLocation(
                    -np.sin(np.deg2rad(d)), np.deg2rad(a), np.sin(np.deg2rad(e)),
                    config.observer_altitude_m))
    cache_dir = resolve_cache_dir(cache_dir)
    fascode_path, ozone_path, solar_path = input_paths(cache_dir, allow_download=allow_download)
    ozone_optical = sk.optical.database.OpticalDatabaseGenericAbsorber(ozone_path)
    std = read_fascode(fascode_path)
    with xr.open_dataset(solar_path) as solar:
        irradiance = integrate_linear(solar.wavelength.values, solar.irradiance.values, edges)
    irradiance /= np.diff(edges) * config.solar_distance_au**2
    optical = None
    growth, wet_radius, wet_index = aerosol_water_uptake(config)
    wet_aod = config.aod550
    mie_paths = []
    if config.aod550 > 0:
        if config.aerosol_model == 'mie':
            dry_optical = _mie(config, wavelengths, cache_dir)
            optical = dry_optical if growth == 1 else _mie(config, wavelengths, cache_dir,
                                                           radius=wet_radius, index=wet_index)
            # Keep the number of dry particles fixed; do not approximate by r^2.
            if growth != 1:
                wet_xs = float(np.asarray(optical.cross_sections([550.], [0.]).extinction).ravel()[0])
                dry_xs = float(np.asarray(dry_optical.cross_sections([550.], [0.]).extinction).ravel()[0])
                if not np.isfinite(wet_xs / dry_xs) or dry_xs <= 0 or wet_xs <= 0:
                    raise ValueError('Invalid hygroscopic Mie extinction cross sections')
                wet_aod = config.aod550 * wet_xs / dry_xs
            mie_paths.extend([Path(dry_optical.path()), Path(optical.path())])
        else:
            hw = np.unique(np.r_[wavelengths, 550.])
            optical = sk.optical.HenyeyGreenstein.from_parameters(
                hw, (hw / 550)**-config.angstrom_exponent,
                np.full_like(hw, config.aerosol_ssa), np.full_like(hw, config.aerosol_g))
    strato = None
    if config.stratospheric_aod550 > 0:
        # Idealized sulfate-like spheres: scenario, not a measured aerosol profile.
        strato = _mie(config, wavelengths, cache_dir, radius=80., width=1.6, index=complex(1.43, -1e-8))
        mie_paths.append(Path(strato.path()))
    fog = None
    if config.fog_aod550 > 0:
        fog = _mie(config, wavelengths, cache_dir, radius=config.fog_median_radius_um * 1000,
                   width=config.fog_mode_width, index=WATER_INDEX)
        mie_paths.append(Path(fog.path()))
    aerosol_ext = normalized_profile(height, np.exp(-height / config.aerosol_scale_height_m), wet_aod)
    fog_ext = normalized_profile(height, np.exp(-height / config.fog_scale_height_m), config.fog_aod550)
    if progress:
        progress(f'Geometry: {len(height)} atmospheric / {len(source_grid(config))} source altitudes, {config.nsza} horizontal columns, '
                 f'{config.nodes} angular nodes, {config.num_stokes} Stokes parameters; '
                 f'{len(ds)*len(az)*len(es)} viewing rays')
    start = time.monotonic()
    result = np.empty((len(ds), len(az), len(es), len(wavelengths), config.num_stokes))
    replace_single = not _single_only and config.single_scatter_step_m > 0
    coarse_single = np.empty_like(result) if replace_single else None
    engine = None if cached_transport is None else cached_transport[1]
    delta_transport = None if cached_transport is None else cached_transport[2]
    for first in range(0, len(wavelengths), config.wavelength_chunk):
        last = min(first + config.wavelength_chunk, len(wavelengths))
        at = sk.Atmosphere(geo, cfg, wavelengths_nm=wavelengths[first:last], calculate_derivatives=False)
        sk.climatology.us76.add_us76_standard_atmosphere(at)
        if config.solar_refraction:
            geo.refractive_index = sk.optical.refraction.ciddor_index_of_refraction(
                at.temperature_k, at.pressure_pa, 0., 420., config.refraction_wavelength_nm)
        at['air'] = sk.constituent.Rayleigh()
        air_density = at.pressure_pa / (Boltzmann * at.temperature_k)
        ozone_density = normalized_profile(height,
            np.interp(height, std['HGT'] * 1000, std['O3'] * 1e-6) * air_density,
            config.ozone_du * DU_MOLECULES_M2)
        at['ozone'] = sk.constituent.VMRAltitudeAbsorber(ozone_optical, height, ozone_density / air_density)
        if optical is not None:
            at['aerosol'] = sk.constituent.ExtinctionScatterer(optical, height, aerosol_ext, 550.)
        if fog is not None:
            at['fog'] = sk.constituent.ExtinctionScatterer(fog, height, fog_ext, 550.)
        if strato is not None:
            ext = normalized_profile(height, np.exp(-.5 * ((height - config.stratospheric_center_m)
                                                           / config.stratospheric_width_m)**2),
                                     config.stratospheric_aod550)
            at['stratospheric_aerosol'] = sk.constituent.ExtinctionScatterer(strato, height, ext, 550.)
        at['solar'] = BinnedSolar(irradiance[first:last])
        at.surface.albedo[:] = config.surface_albedo
        if engine is None:  # initialize after setting refractive index
            engine = sk.Engine(cfg, geo, view)
            if config.delta_m_method == 'split' and not _single_only:
                from .delta_m import DeltaMTransport
                delta_transport = DeltaMTransport(cfg, geo, view, config.delta_m_order)
            if _transport_cache is not None:
                _transport_cache.clear()  # at most one potentially large field
                _transport_cache[transport_key] = (geo, engine, delta_transport)
        output = engine.calculate_radiance(at)
        if replace_single:
            coarse_single[..., first:last, :] = output.radiance.transpose('los', 'wavelength', 'stokes').values.reshape(
                len(ds),len(az),len(es),last-first,config.num_stokes)
        if delta_transport is not None:
            output['radiance'] = output.radiance + delta_transport.correction(at)
        values = output.radiance.transpose('los', 'wavelength', 'stokes').values
        result[..., first:last, :] = values.reshape(len(ds), len(az), len(es), last-first, config.num_stokes)
        if progress:
            progress(f'wavelengths {first+1}..{last}/{len(wavelengths)}, elapsed {time.monotonic()-start:.1f}s')
    fine_single = None
    if replace_single:
        if progress:
            progress(f'Refining observer single scattering on a {config.single_scatter_step_m:g} m upper-atmosphere grid')
        fine_single = compute_model(ds,config,es,az,cache_dir=cache_dir,allow_download=allow_download,_single_only=True)
        result = result - coarse_single + fine_single.solar_stokes_radiance.values
    # Finite Legendre/TMS component approximations can be slightly signed.
    # Never clip them before subtraction/replacement; the published total
    # must still be finite and nonnegative.
    if not np.isfinite(result).all() or (not _single_only and np.any(result[..., 0] < 0)):
        raise ValueError('Nonfinite or negative intensity; refuse to publish this run')
    if config.num_stokes == 3 and np.any(np.linalg.norm(result[..., 1:], axis=-1) > result[..., 0] * 1.001 + 1e-20):
        raise ValueError('Unphysical polarization; increase angular resolution')
    dims = ('depression', 'azimuth', 'elevation', 'wavelength', 'stokes')
    dataset = xr.Dataset(
        {'solar_stokes_radiance': (dims, result),
         'wavelength_bounds': (('wavelength', 'bound'), np.c_[edges[:-1], edges[1:]]),
         'solar_irradiance': ('wavelength', irradiance),
         'pressure': ('altitude', at.pressure_pa), 'temperature': ('altitude', at.temperature_k),
         'ozone_number_density': ('altitude', ozone_density)},
        coords={'depression': ds, 'azimuth': az, 'elevation': es, 'wavelength': wavelengths,
                'stokes': ['I', 'Q', 'U'][:config.num_stokes], 'altitude': height, 'bound': [0, 1]},
        attrs={'library_version': __version__, 'model_version': MODEL_VERSION, 'config_json': json.dumps(config.to_dict(), sort_keys=True),
               'solver': 'SASKTRAN2 Geometry2D SuccessiveOrders', 'engine_version': importlib.metadata.version('sasktran2'),
               'spectral_convention': 'bin-mean solar spectrum times centre-wavelength transport; 380..780 nm',
               'solar_refraction': int(config.solar_refraction), 'los_refraction': 0, 'diffuse_refraction': 0,
               'refraction_note': 'Ciddor dry air, 420 ppm CO2; geometry at configured reference wavelength',
               'solar_disk': 'collimated point source; finite disk and direct disk radiance excluded',
               'background': 'excluded unless an observed at-observer spectrum is supplied at export',
               'angle_convention': 'geometric solar centre depression; geometric elevation; azimuth 0 sunward, 180 antisolar',
               'convergence': 'not established by a single run; use validate_sensitivity',
               'aerosol_aod550_dry': config.aod550, 'aerosol_aod550_ambient': wet_aod,
               'aerosol_growth_factor': growth, 'aerosol_wet_median_radius_nm': wet_radius,
               'humidity_note': 'uniform aerosol-layer RH; kappa growth without Kelvin term; volume-mixed visible constant index; gas US76 unchanged',
               'fog_note': 'independent spherical liquid droplets, exponential surface layer, visible constant water index; no saturation/condensation dynamics',
               'delta_m_scaling': int(config.delta_m_scaling), 'delta_m_order': config.delta_m_order,
               'delta_m_method': config.delta_m_method,
               'radiance_component': 'single_scattering_only' if _single_only else 'total_solar',
               'single_scatter_step_m': config.single_scatter_step_m,
               'delta_m_note': ('Consistently truncated diffuse transport plus full-phase observer single-scattering correction at scaled optical depth'
                                if config.delta_m_method == 'split' else 'Native backend delta-M treatment'),
               'transport_geometry_reused': int(cached_transport is not None),
               'elapsed_seconds': time.monotonic()-start,
               'transport_source_sha256': sha256(Path(__file__)),
               'microphysics_source_sha256': sha256(Path(__file__).with_name('microphysics.py')),
               'delta_m_source_sha256': sha256(Path(__file__).with_name('delta_m.py')),
               'physical_inputs_json': json.dumps([{'path': str(p.relative_to(cache_dir)), 'sha256': sha256(p)}
                                                   for p in dict.fromkeys(map(Path, [fascode_path, ozone_path, solar_path, *mie_paths]))])})
    dataset['aerosol_extinction_550'] = ('altitude', aerosol_ext)
    dataset['fog_extinction_550'] = ('altitude', fog_ext)
    for name in ('aerosol_extinction_550', 'fog_extinction_550'):
        dataset[name].attrs['units'] = 'm-1'
    for name in ('depression', 'azimuth', 'elevation'):
        dataset[name].attrs['units'] = 'degree'
    dataset.wavelength.attrs.update(units='nm', bounds='wavelength_bounds')
    dataset.wavelength_bounds.attrs['units'] = 'nm'
    dataset.solar_stokes_radiance.attrs['units'] = 'W m-2 sr-1 nm-1'
    dataset.solar_irradiance.attrs['units'] = 'W m-2 nm-1'
    dataset.altitude.attrs['units'] = 'm'
    dataset.pressure.attrs['units'] = 'Pa'
    dataset.temperature.attrs['units'] = 'K'
    dataset.ozone_number_density.attrs['units'] = 'm-3'
    dataset = dataset.assign_coords(source_altitude=('source_altitude', source_grid(config)))
    dataset.source_altitude.attrs['units'] = 'm'
    dataset = dataset.assign_coords(source_horizontal_angle=('source_horizontal_angle',horizontal_grid(config)))
    dataset.source_horizontal_angle.attrs['units'] = 'degree'
    if fine_single is not None:
        dataset['solar_single_scatter_radiance'] = (dims, fine_single.solar_stokes_radiance.values)
        dataset.solar_single_scatter_radiance.attrs['units'] = 'W m-2 sr-1 nm-1'
        dataset = dataset.assign_coords(single_scatter_altitude=('single_scatter_altitude',fine_single.altitude.values))
        dataset.single_scatter_altitude.attrs['units'] = 'm'
        dataset.attrs['single_scatter_min_radiance'] = float(fine_single.solar_stokes_radiance.values[...,0].min())
        dataset.attrs['single_scatter_negative_samples'] = int((fine_single.solar_stokes_radiance.values[...,0]<0).sum())
    return dataset


def replace_single_scattering(reference, config, *, cache_dir=None, allow_download=True,
                              coarse_cache=None, fine_cache=None):
    """Reuse a saved split diffuse field and recompute only observer single scatter.

    All settings except the independent single-scatter step must match. This is
    an exact component replacement for the defined split solver, not spectral
    or colour interpolation. A fresh full solve is used as a regression check
    when migrating the explorer. Existing colour variables must be regenerated
    with add_colorimetry by the caller.
    """
    before = ModelConfig.from_dict(json.loads(reference.attrs['config_json']))
    config = config.validate()
    if before.delta_m_method != 'split' or not config.single_scatter_step_m:
        raise ValueError('Single-scatter replacement requires split transport and a separate target grid')
    if any(getattr(before,k) != getattr(config,k) for k in before.to_dict() if k != 'single_scatter_step_m'):
        raise ValueError('Cannot reuse a diffuse field after changing other settings')
    if reference.attrs.get('radiance_component','total_solar') != 'total_solar':
        raise ValueError('Expected total solar radiance, not a diagnostic component')
    args = [reference.depression.values, before, reference.elevation.values, reference.azimuth.values]
    if 'solar_single_scatter_radiance' in reference:
        previous_single = reference.solar_single_scatter_radiance.values
    else:
        previous_single = compute_model(*args,cache_dir=cache_dir,allow_download=allow_download,
            _transport_cache=coarse_cache,_single_only=True).solar_stokes_radiance.values
    args[1] = config
    single = compute_model(*args,cache_dir=cache_dir,allow_download=allow_download,
        _transport_cache=fine_cache,_single_only=True)
    values = reference.solar_stokes_radiance.values - previous_single + single.solar_stokes_radiance.values
    if not np.isfinite(values).all() or (values[...,0] < 0).any():
        raise ValueError('Invalid component replacement; refusing to publish')
    result = reference.copy(deep=True)
    result['solar_stokes_radiance'] = (reference.solar_stokes_radiance.dims,values)
    result.solar_stokes_radiance.attrs['units'] = 'W m-2 sr-1 nm-1'
    result['solar_single_scatter_radiance'] = (reference.solar_stokes_radiance.dims,single.solar_stokes_radiance.values)
    result.solar_single_scatter_radiance.attrs['units'] = 'W m-2 sr-1 nm-1'
    if 'single_scatter_altitude' in result.coords:
        result = result.drop_vars('single_scatter_altitude')
    result = result.assign_coords(single_scatter_altitude=('single_scatter_altitude',single.altitude.values))
    result.single_scatter_altitude.attrs['units'] = 'm'
    result.attrs.update(config_json=json.dumps(config.to_dict(),sort_keys=True),
        single_scatter_step_m=config.single_scatter_step_m, radiance_component='total_solar',
        diffuse_component_reused_from_source_sha256=reference.attrs['transport_source_sha256'],
        single_scatter_replacement='Saved total minus saved/recomputed observer single scatter plus refined single scatter; diffuse transport unchanged',
        transport_source_sha256=sha256(Path(__file__)),
        library_version=__version__,model_version=MODEL_VERSION)
    result.attrs['single_scatter_min_radiance'] = float(single.solar_stokes_radiance.values[...,0].min())
    result.attrs['single_scatter_negative_samples'] = int((single.solar_stokes_radiance.values[...,0]<0).sum())
    return result
