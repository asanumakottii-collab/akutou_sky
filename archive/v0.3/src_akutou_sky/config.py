"""Validated, serializable physical and numerical inputs for the twilight model."""
from dataclasses import asdict, dataclass, fields
import json
import math
from pathlib import Path
import numpy as np

@dataclass(frozen=True)
class ModelConfig:
    """Immutable physical inputs and numerical controls; lengths are metres.

    Use from_preset for a validated starting point, or dataclasses.replace to
    vary an existing configuration. Presets are workloads, not accuracy grades.
    """
    aod550: float = 0.10
    aerosol_model: str = 'mie'
    aerosol_scale_height_m: float = 1500.
    angstrom_exponent: float = 1.3
    aerosol_ssa: float = .95
    aerosol_g: float = .65
    median_radius_nm: float = 100.
    mode_width: float = 1.6
    refractive_index_real: float = 1.45
    refractive_index_imag: float = .005  # positive absorption; backend uses n - ik
    relative_humidity: float = 0.  # fraction, applied to aerosol, not gas climatology
    hygroscopicity_kappa: float = 0.
    fog_aod550: float = 0.
    fog_scale_height_m: float = 200.
    fog_median_radius_um: float = 5.
    fog_mode_width: float = 1.3
    stratospheric_aod550: float = 0.
    stratospheric_center_m: float = 20000.
    stratospheric_width_m: float = 4000.
    ozone_du: float = 300.
    surface_albedo: float = .1
    observer_altitude_m: float = 2.
    solar_distance_au: float = 1.
    solar_refraction: bool = True
    refraction_wavelength_nm: float = 550.
    num_stokes: int = 1
    wavelength_step_nm: float = 5.
    nodes: int = 38
    horizon_quadrature: bool = True
    nsza: int = 41
    horizontal_sampling: str = 'twilight'
    ground_step_m: float = 125.
    boundary_step_m: float = 0.  # optional refinement below 1 km; 0 disables it
    middle_step_m: float = 500.
    upper_step_m: float = 2000.
    source_ground_step_m: float = 500.
    source_boundary_step_m: float = 0.
    source_surface_m: float = 0.  # 0 keeps legacy midpoints; >0 resolves surface source
    source_middle_step_m: float = 1000.
    source_upper_step_m: float = 4000.
    max_source_samples: int = 120000
    max_view_source_product: int = 350000000
    max_atmosphere_source_product: int = 20000000
    iterations: int = 80
    relative_tolerance: float = 1e-8
    absolute_tolerance: float = 1e-20
    legendre_moments: int = 64
    delta_m_scaling: bool = False
    delta_m_order: int = 16
    delta_m_method: str = 'native'
    single_scatter_step_m: float = 0.  # separate, nested SS grid; 0 uses the transport grid
    wavelength_chunk: int = 1
    num_threads: int = 1

    def validate(self):
        for f in fields(self):
            value = getattr(self, f.name)
            if f.name not in ('aerosol_model', 'solar_refraction', 'horizon_quadrature', 'delta_m_scaling', 'horizontal_sampling', 'delta_m_method'):
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                    raise ValueError(f'{f.name} must be a finite number')
        for name in ('aod550', 'stratospheric_aod550', 'ozone_du', 'refractive_index_imag',
                     'relative_tolerance', 'absolute_tolerance', 'fog_aod550', 'hygroscopicity_kappa'):
            if getattr(self, name) < 0:
                raise ValueError(f'{name} must be nonnegative')
        for name in ('surface_albedo', 'aerosol_ssa'):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f'{name} must be between 0 and 1')
        for name in ('aerosol_scale_height_m', 'stratospheric_width_m', 'solar_distance_au'):
            if getattr(self, name) <= 0:
                raise ValueError(f'{name} must be positive')
        if self.aerosol_model not in ('mie', 'hg'):
            raise ValueError('aerosol_model must be mie or hg')
        if self.horizontal_sampling not in ('twilight','uniform','near_horizon'):
            raise ValueError('horizontal_sampling must be twilight, uniform or near_horizon')
        if self.horizontal_sampling == 'near_horizon' and self.nsza < 25:
            raise ValueError('near_horizon requires nsza >= 25')
        if self.horizontal_sampling == 'twilight' and self.nsza < 11:
            raise ValueError('twilight horizontal sampling requires nsza >= 11')
        if not -.99 < self.aerosol_g < .99:
            raise ValueError('aerosol_g must lie between -0.99 and 0.99')
        if not 1 < self.mode_width <= 2.5 or not 1 <= self.refractive_index_real <= 3:
            raise ValueError('Mie model requires 1 < mode_width <= 2.5 and 1 <= real index <= 3')
        if not 1 <= self.median_radius_nm <= 1000:
            raise ValueError('median_radius_nm must be 1..1000')
        if not 0 <= self.relative_humidity <= .95:
            raise ValueError('relative_humidity must be a fraction in [0, 0.95]; saturation is excluded')
        if self.hygroscopicity_kappa > 1.5:
            raise ValueError('hygroscopicity_kappa must be in [0, 1.5]')
        if self.relative_humidity * self.hygroscopicity_kappa > 0:
            if self.aerosol_model != 'mie' or self.median_radius_nm < 40:
                raise ValueError('Hygroscopic growth requires Mie aerosol with dry median radius >= 40 nm')
        if not 50 <= self.fog_scale_height_m <= 1000:
            raise ValueError('fog_scale_height_m must be 50..1000 m')
        if not .5 <= self.fog_median_radius_um <= 10 or not 1 < self.fog_mode_width <= 1.5:
            raise ValueError('Fog requires radius 0.5..10 um and 1 < mode width <= 1.5')
        for name in ('boundary_step_m', 'source_boundary_step_m'):
            value = getattr(self, name)
            if value != 0 and not 10 <= value <= 1000:
                raise ValueError(f'{name} must be 0 (disabled) or 10..1000 m')
        if self.source_surface_m != 0 and not .1 <= self.source_surface_m <= 10:
            raise ValueError('source_surface_m must be 0 (disabled) or 0.1..10 m')
        if not isinstance(self.delta_m_scaling, bool):
            raise ValueError('delta_m_scaling must be a boolean')
        if self.delta_m_method not in ('native', 'split'):
            raise ValueError('delta_m_method must be native or split')
        if self.delta_m_method == 'split' and (not self.delta_m_scaling or self.num_stokes != 1):
            raise ValueError('Split delta-M requires delta_m_scaling=True and num_stokes=1')
        if self.single_scatter_step_m != 0:
            if not 50 <= self.single_scatter_step_m <= 2000 or self.delta_m_method != 'split':
                raise ValueError('single_scatter_step_m requires split delta-M and a step of 50..2000 m (or 0 disables it)')
        if not isinstance(self.delta_m_order, int) or self.delta_m_order not in (2, 4, 8, 16, 32, 64, 128):
            raise ValueError('delta_m_order must be one of 2, 4, 8, 16, 32, 64, 128')
        if self.delta_m_scaling and self.legendre_moments <= self.delta_m_order:
            raise ValueError('Delta-M requires legendre_moments > delta_m_order')
        if self.fog_aod550 > 0:
            if self.legendre_moments < 256 or not self.delta_m_scaling:
                raise ValueError('Fog requires >= 256 Legendre moments and delta_m_scaling=True; use fog_preview preset')
            if (self.boundary_step_m or self.ground_step_m) > self.fog_scale_height_m / 4:
                raise ValueError('Atmosphere boundary step must be <= fog scale height / 4')
            if (self.source_boundary_step_m or self.source_ground_step_m) > self.fog_scale_height_m / 2:
                raise ValueError('Source boundary step must be <= fog scale height / 2')
        if not 0 <= self.observer_altitude_m < 100000:
            raise ValueError('observer_altitude_m must lie in [0, 100000)')
        if not 0 < self.stratospheric_center_m < 100000:
            raise ValueError('stratospheric_center_m must lie in (0, 100000)')
        if not 380 <= self.refraction_wavelength_nm <= 780:
            raise ValueError('refraction_wavelength_nm must lie in [380, 780]')
        if not .5 <= self.wavelength_step_nm <= 40:
            raise ValueError('wavelength_step_nm must lie in [0.5, 40]')
        if not math.isclose(400 / self.wavelength_step_nm, round(400 / self.wavelength_step_nm)):
            raise ValueError('wavelength_step_nm must divide 400 nm exactly')
        if self.num_stokes not in (1, 3):
            raise ValueError('num_stokes must be 1 or 3')
        if self.nodes not in (6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302):
            raise ValueError('nodes must be a supported Lebedev size (e.g. 26, 50, 74, 110)')
        for name, minimum in (('nsza', 3), ('iterations', 1), ('legendre_moments', 16),
                              ('wavelength_chunk', 1), ('num_threads', 1), ('nodes', 6), ('num_stokes', 1)):
            value = getattr(self, name)
            if not isinstance(value, int) or value < minimum:
                raise ValueError(f'{name} must be an integer >= {minimum}')
        if not isinstance(self.solar_refraction, bool):
            raise ValueError('solar_refraction must be a boolean')
        if not isinstance(self.horizon_quadrature, bool):
            raise ValueError('horizon_quadrature must be a boolean')
        if min(self.ground_step_m, self.middle_step_m, self.upper_step_m) < 10:
            raise ValueError('altitude steps must be >= 10 m')
        if min(self.source_ground_step_m, self.source_middle_step_m, self.source_upper_step_m) < 10:
            raise ValueError('source altitude steps must be >= 10 m')
        if not isinstance(self.max_source_samples, int) or self.max_source_samples < 1:
            raise ValueError('max_source_samples must be a positive integer')
        if not isinstance(self.max_view_source_product, int) or self.max_view_source_product < 1:
            raise ValueError('max_view_source_product must be a positive integer')
        samples = len(source_grid(self)) * self.nsza * self.nodes * self.num_stokes
        if samples > self.max_source_samples:
            raise ValueError(f'Source field has {samples:,} samples, exceeding budget {self.max_source_samples:,}. '
                             'Use a smaller grid or --preset polarized for polarization. '
                             'Only increase --max-source-samples on a machine with sufficient memory; '
                             'this sample budget is not an exact RAM estimate.')
        if not isinstance(self.max_atmosphere_source_product,int) or self.max_atmosphere_source_product<1:
            raise ValueError('max_atmosphere_source_product must be a positive integer')
        transport_size = samples * len(altitude_grid(self))
        if transport_size > self.max_atmosphere_source_product:
            raise ValueError(f'Atmosphere/source product {transport_size:,} exceeds budget '
                             f'{self.max_atmosphere_source_product:,}; reduce vertical or source resolution')
        return self

    def to_dict(self):
        """Return JSON-compatible values for all physical and numerical fields."""
        return asdict(self)

    def to_json(self, path=None):
        """Serialize settings; optionally write to the explicitly supplied path."""
        text = json.dumps(self.to_dict(), indent=2) + '\n'
        if path is not None:
            Path(path).write_text(text, encoding='utf-8')
        return text

    @classmethod
    def from_json(cls, path):
        """Read a settings JSON file and reject unknown or invalid fields."""
        return cls.from_dict(json.loads(Path(path).read_text(encoding='utf-8')))

    @classmethod
    def from_preset(cls, name='standard', **overrides):
        """Create a validated preset, with explicit keyword overrides."""
        if name not in PRESETS:
            raise ValueError(f'Unknown preset {name!r}; choose from {tuple(PRESETS)}')
        return cls.from_dict({**PRESETS[name], **overrides})

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, dict):
            raise ValueError('Model settings must be a JSON object / dictionary')
        unknown = set(data) - {f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f'Unknown model settings: {sorted(unknown)}')
        return cls(**data).validate()

PRESETS = {
    'fog_refined': dict(nodes=110, nsza=41, horizontal_sampling='near_horizon',
                        ground_step_m=250., boundary_step_m=50., middle_step_m=1000.,
                        source_ground_step_m=1000., source_boundary_step_m=100., source_surface_m=1.,
                        source_middle_step_m=2000., source_upper_step_m=4000.,
                        wavelength_step_nm=20., num_stokes=1, iterations=160,
                        relative_tolerance=1e-12, absolute_tolerance=1e-24,
                        legendre_moments=256, delta_m_scaling=True, delta_m_method='split', single_scatter_step_m=500.,
                        max_source_samples=300000, max_atmosphere_source_product=40000000),
    'fog_preview': dict(nodes=14, nsza=21, ground_step_m=250., boundary_step_m=50.,
                        middle_step_m=1000., source_ground_step_m=1000., source_boundary_step_m=100.,
                        source_middle_step_m=2000., source_upper_step_m=4000.,
                        wavelength_step_nm=20., num_stokes=1, iterations=80,
                        legendre_moments=256, delta_m_scaling=True),
    'preview': dict(nodes=14, nsza=21, ground_step_m=250., middle_step_m=1000.,
                    source_ground_step_m=1000., source_middle_step_m=2000., source_upper_step_m=4000.,
                    wavelength_step_nm=10., num_stokes=1, iterations=60),
    'standard': {},
    'overview': dict(nodes=26,nsza=41,source_ground_step_m=500.,source_upper_step_m=2000.,horizontal_sampling='uniform'),
    'polarized': dict(nodes=14, nsza=31, ground_step_m=250., middle_step_m=1000.,
                      source_ground_step_m=500.,source_upper_step_m=2000.,
                      wavelength_step_nm=10., num_stokes=3, wavelength_chunk=1),
}

def altitude_grid(config):
    boundary = np.arange(0, 1000, config.boundary_step_m) if config.boundary_step_m else []
    lower = 1000 if config.boundary_step_m else 0
    return np.unique(np.r_[boundary, np.arange(lower, 10000, config.ground_step_m),
                           np.arange(10000, 30000, config.middle_step_m),
                           np.arange(30000, 100000, config.upper_step_m), 100000.])

def spectral_grid(step):
    edges = np.linspace(380., 780., round(400 / step) + 1)
    return (edges[:-1] + edges[1:]) / 2, edges

def source_grid(config):
    boundary = np.arange(0, 1000, config.source_boundary_step_m) if config.source_boundary_step_m else []
    lower = 1000 if config.source_boundary_step_m else 0
    edges = np.unique(np.r_[boundary, np.arange(lower,10000,config.source_ground_step_m),
                            np.arange(10000,30000,config.source_middle_step_m),
                            np.arange(30000,100000,config.source_upper_step_m),100000.])
    midpoints = (edges[:-1]+edges[1:])/2
    if config.source_surface_m:
        return np.unique(np.r_[config.source_surface_m, midpoints])
    return midpoints

def single_scatter_grid(config):
    """Independent nested grid for the observer single-scattering integral.

    Keep solar-shadow sampling independent of changes to the diffuse grid.
    Halving the single-scatter step also halves the low-altitude steps.
    """
    step = config.single_scatter_step_m
    if not step:
        return altitude_grid(config)
    low = min(step/10, config.fog_scale_height_m/4) if config.fog_aod550 else step/10
    return np.unique(np.r_[np.arange(0,1000,low), np.arange(1000,10000,step/2),
                           np.arange(10000,100000,step),100000.])

def horizontal_grid(config):
    if config.horizontal_sampling == 'uniform':
        return np.linspace(-45.,45.,config.nsza)
    if config.horizontal_sampling == 'near_horizon':
        # Symmetric coverage of sunrise/sunset. 41 -> 73 -> 137 columns gives
        # nested 0.5 -> 0.25 -> 0.125 degree spacing, retaining outer transport.
        return np.r_[[-45., -30., -20., -12.],
                     np.linspace(-8., 8., config.nsza-8), [12., 20., 30., 45.]]
    # At reference SZA=90°, local horizontal angle is minus solar depression.
    # Resolve the steep twilight source field, retaining distant day/night columns.
    return np.r_[[-45.,-35.,-30.,-25.],np.linspace(-20.,5.,config.nsza-7),[15.,30.,45.]]

def checked_angles(values, name, lower, upper):
    arr = np.atleast_1d(values).astype(float)
    if arr.ndim != 1 or not len(arr) or not np.isfinite(arr).all():
        raise ValueError(f'{name} must be a nonempty finite one-dimensional array')
    if (arr < lower).any() or (arr > upper).any() or len(np.unique(arr)) != len(arr):
        raise ValueError(f'{name} must be unique and within [{lower}, {upper}]')
    return arr

def integrate_linear(x, y, edges):
    """Exact integrals of a piecewise-linear spectrum over requested bins."""
    x, y, edges = np.asarray(x), np.asarray(y), np.asarray(edges)
    if x.ndim != 1 or y.shape != x.shape or len(x) < 2 or not np.all(np.diff(x) > 0):
        raise ValueError('Spectrum coordinates must be strictly increasing')
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.isfinite(edges).all():
        raise ValueError('Spectrum must be finite')
    if edges.ndim != 1 or len(edges) < 2 or not np.all(np.diff(edges) > 0) or edges[0] < x[0] or edges[-1] > x[-1]:
        raise ValueError('Integration bins must be ordered and covered by the spectrum')
    cumulative = np.r_[0., np.cumsum(np.diff(x) * (y[1:] + y[:-1]) / 2)]
    indices = np.clip(np.searchsorted(x, edges, side='right') - 1, 0, len(x) - 2)
    dx = edges - x[indices]
    slopes = np.diff(y) / np.diff(x)
    integral = cumulative[indices] + y[indices] * dx + slopes[indices] * dx**2 / 2
    return np.diff(integral)
