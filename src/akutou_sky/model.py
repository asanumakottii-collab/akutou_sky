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
from .config import ModelConfig, altitude_grid, source_grid, horizontal_grid, spectral_grid, checked_angles, integrate_linear
from .data import input_paths, resolve_cache_dir
from ._version import __version__
MODEL_VERSION = '2.0'
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
                  *, cache_dir=None, allow_download=True):
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
    height = altitude_grid(config)
    cfg = sk.Config()
    cfg.num_threads = config.num_threads
    cfg.num_stokes = config.num_stokes
    cfg.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
    cfg.num_successive_orders_incoming = config.nodes
    cfg.num_successive_orders_outgoing = config.nodes
    cfg.num_sza = config.nsza
    cfg.successive_orders_altitude_grid_m = source_grid(config)
    cfg.successive_orders_horizontal_angle_grid_radians = np.deg2rad(horizontal_grid(config))
    cfg.num_successive_orders_iterations = config.iterations
    cfg.successive_orders_relative_tolerance = config.relative_tolerance
    cfg.successive_orders_absolute_tolerance = config.absolute_tolerance
    cfg.num_singlescatter_moments = config.legendre_moments
    cfg.solar_refraction = config.solar_refraction
    cfg.los_refraction = False  # not supported by the Geometry2D solver
    cfg.multiple_scatter_refraction = False  # likewise not supported
    cfg.wavelength_batch_size = 1
    geo = sk.Geometry2D(cos_sza=0., solar_azimuth=0., earth_radius_m=6371000.,
                        altitude_grid_m=height, horizontal_angle_grid_radians=np.deg2rad([-45., 45.]))
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
    if config.aod550 > 0:
        if config.aerosol_model == 'mie':
            optical = _mie(config, wavelengths, cache_dir)
        else:
            hw = np.unique(np.r_[wavelengths, 550.])
            optical = sk.optical.HenyeyGreenstein.from_parameters(
                hw, (hw / 550)**-config.angstrom_exponent,
                np.full_like(hw, config.aerosol_ssa), np.full_like(hw, config.aerosol_g))
    strato = None
    if config.stratospheric_aod550 > 0:
        # Idealized sulfate-like spheres: scenario, not a measured aerosol profile.
        strato = _mie(config, wavelengths, cache_dir, radius=80., width=1.6, index=complex(1.43, -1e-8))
    if progress:
        progress(f'Geometry: {len(height)} atmospheric / {len(source_grid(config))} source altitudes, {config.nsza} horizontal columns, '
                 f'{config.nodes} angular nodes, {config.num_stokes} Stokes parameters; '
                 f'{len(ds)*len(az)*len(es)} viewing rays')
    start = time.monotonic()
    result = np.empty((len(ds), len(az), len(es), len(wavelengths), config.num_stokes))
    engine = None
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
            ext = normalized_profile(height, np.exp(-height / config.aerosol_scale_height_m), config.aod550)
            at['aerosol'] = sk.constituent.ExtinctionScatterer(optical, height, ext, 550.)
        if strato is not None:
            ext = normalized_profile(height, np.exp(-.5 * ((height - config.stratospheric_center_m)
                                                           / config.stratospheric_width_m)**2),
                                     config.stratospheric_aod550)
            at['stratospheric_aerosol'] = sk.constituent.ExtinctionScatterer(strato, height, ext, 550.)
        at['solar'] = BinnedSolar(irradiance[first:last])
        at.surface.albedo[:] = config.surface_albedo
        if engine is None:  # initialize after setting refractive index
            engine = sk.Engine(cfg, geo, view)
        output = engine.calculate_radiance(at)
        values = output.radiance.transpose('los', 'wavelength', 'stokes').values
        result[..., first:last, :] = values.reshape(len(ds), len(az), len(es), last-first, config.num_stokes)
        if progress:
            progress(f'wavelengths {first+1}..{last}/{len(wavelengths)}, elapsed {time.monotonic()-start:.1f}s')
    if not np.isfinite(result).all() or np.any(result[..., 0] < 0):
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
               'elapsed_seconds': time.monotonic()-start,
               'transport_source_sha256': sha256(Path(__file__)),
               'physical_inputs_json': json.dumps([{'path': str(p.relative_to(cache_dir)), 'sha256': sha256(p)}
                                                   for p in map(Path, [fascode_path, ozone_path, solar_path])])})
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
    return dataset

