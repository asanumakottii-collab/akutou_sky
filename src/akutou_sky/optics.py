"""Scene -> SASKTRAN2 optical constituents (gas, aerosol, fog, sunlight).

Hygroscopic growth: Petters & Kreidenweis (2007), doi:10.5194/acp-7-1961-2007,
without the Kelvin term; one growth factor shifts the lognormal mode while the
particle number fixed by the dry AOD is kept.
"""
from pathlib import Path
import numpy as np
import sasktran2 as sk
import xarray as xr
from scipy.constants import Boltzmann

DU_MOLECULES_M2 = 2.687e20
WATER_INDEX = complex(1.333, -1e-9)          # visible constant-index approximation
SULFATE_INDEX = complex(1.43, -1e-8)         # background stratospheric aerosol
SULFATE_RADIUS_NM, SULFATE_WIDTH = 80., 1.6


def water_uptake(scene):
    """Return growth factor, wet median radius (nm) and volume-mixed index."""
    g3 = 1 + scene.hygroscopicity_kappa * scene.relative_humidity / (1 - scene.relative_humidity)
    dry = complex(scene.refractive_index_real, -scene.refractive_index_imag)
    return g3 ** (1 / 3), scene.median_radius_nm * g3 ** (1 / 3), (dry + WATER_INDEX * (g3 - 1)) / g3


def column_profile(height, shape, column):
    """Scale a vertical shape so that its trapezoidal integral equals column."""
    integral = np.trapezoid(shape, height)
    if not np.isfinite(integral) or integral <= 0:
        raise ValueError('Vertical profile has zero or invalid column')
    return shape * column / integral


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


def integrate_linear(x, y, edges):
    """Exact integrals of a piecewise-linear spectrum over the given bins."""
    x, y, edges = np.asarray(x), np.asarray(y), np.asarray(edges)
    if edges[0] < x[0] or edges[-1] > x[-1] or np.any(np.diff(x) <= 0):
        raise ValueError('Spectrum must be increasing and cover the integration bins')
    cumulative = np.r_[0., np.cumsum(np.diff(x) * (y[1:] + y[:-1]) / 2)]
    i = np.clip(np.searchsorted(x, edges, side='right') - 1, 0, len(x) - 2)
    dx = edges - x[i]
    integral = cumulative[i] + y[i] * dx + np.diff(y)[i] / np.diff(x)[i] * dx**2 / 2
    return np.diff(integral)


class _Solar(sk.constituent.base.Constituent):
    """Bin-mean top-of-atmosphere irradiance."""
    def __init__(self, values):
        self.values = values

    def add_to_atmosphere(self, atmo):
        atmo.storage.solar_irradiance[:] = self.values

    def register_derivative(self, atmo, name):
        return {}


def _mie(wavelengths, cache_dir, radius_nm, width, index, moments, threads):
    distribution = sk.mie.distribution.LogNormalDistribution().freeze(median_radius=radius_nm, mode_width=width)
    refractive = sk.mie.refractive.RefractiveIndex(
        lambda w: np.full_like(np.asarray(w, dtype=float), index, dtype=complex),
        f'constant_n{index.real:.12g}_k{-index.imag:.12g}')
    return sk.database.MieDatabase(distribution, refractive, wavelengths_nm=np.unique(np.r_[wavelengths, 550.]),
                                   db_root=cache_dir, max_legendre_moments=moments, num_threads=threads)


def _extinction_550(database):
    return float(np.asarray(database.cross_sections([550.], [0.]).extinction).ravel()[0])


class Optics:
    """Everything wavelength-dependent that a Scene puts into the atmosphere.

    Built once per render; add_to() then fills any SASKTRAN2 Atmosphere whose
    altitude grid may differ (coarse diffuse grid or fine single-scatter grid).
    """

    def __init__(self, scene, wavelengths, edges, paths, cache_dir, moments, threads):
        fascode, ozone, solar = paths
        self.scene = scene
        self.std = read_fascode(fascode)
        self.ozone_optical = sk.optical.database.OpticalDatabaseGenericAbsorber(ozone)
        with xr.open_dataset(solar) as sun:
            self.irradiance = (integrate_linear(sun.wavelength.values, sun.irradiance.values, edges)
                               / np.diff(edges) / scene.solar_distance_au**2)
        self.growth, self.wet_radius_nm, wet_index = water_uptake(scene)
        self.ambient_aod550 = scene.aod550
        self.aerosol = self.strato = self.fog = None
        mie = lambda r, w, n: _mie(wavelengths, cache_dir, r, w, n, moments, threads)
        if scene.aod550 > 0:
            dry = mie(scene.median_radius_nm, scene.mode_width,
                      complex(scene.refractive_index_real, -scene.refractive_index_imag))
            self.aerosol = dry
            if self.growth != 1:
                self.aerosol = mie(self.wet_radius_nm, scene.mode_width, wet_index)
                self.ambient_aod550 = scene.aod550 * _extinction_550(self.aerosol) / _extinction_550(dry)
        if scene.stratospheric_aod550 > 0:
            self.strato = mie(SULFATE_RADIUS_NM, SULFATE_WIDTH, SULFATE_INDEX)
        if scene.has_fog:
            self.fog = mie(scene.fog_median_radius_um * 1000, scene.fog_mode_width, WATER_INDEX)

    def profiles(self, height):
        """550 nm extinction profiles (m-1) on the given altitude grid."""
        s = self.scene
        out = {}
        if self.aerosol is not None:
            out['aerosol'] = column_profile(height, np.exp(-height / s.aerosol_scale_height_m), self.ambient_aod550)
        if self.strato is not None:
            out['stratospheric_aerosol'] = column_profile(
                height, np.exp(-.5 * ((height - s.stratospheric_center_m) / s.stratospheric_width_m)**2),
                s.stratospheric_aod550)
        if self.fog is not None:
            out['fog'] = column_profile(height, np.exp(-height / s.fog_scale_height_m), s.fog_aod550)
        return out

    def add_to(self, atmo, height, first, last):
        """Populate atmo (wavelength bins first..last) on altitude grid height."""
        sk.climatology.us76.add_us76_standard_atmosphere(atmo)
        atmo['air'] = sk.constituent.Rayleigh()
        air = atmo.pressure_pa / (Boltzmann * atmo.temperature_k)
        ozone = column_profile(height, np.interp(height, self.std['HGT'] * 1000, self.std['O3'] * 1e-6) * air,
                               self.scene.ozone_du * DU_MOLECULES_M2)
        atmo['ozone'] = sk.constituent.VMRAltitudeAbsorber(self.ozone_optical, height, ozone / air)
        databases = dict(aerosol=self.aerosol, stratospheric_aerosol=self.strato, fog=self.fog)
        for name, extinction in self.profiles(height).items():
            atmo[name] = sk.constituent.ExtinctionScatterer(databases[name], height, extinction, 550.)
        atmo['solar'] = _Solar(self.irradiance[first:last])
        atmo.surface.albedo[:] = self.scene.surface_albedo
        return ozone

    def mie_paths(self):
        return [Path(db.path()) for db in (self.aerosol, self.strato, self.fog) if db is not None]
