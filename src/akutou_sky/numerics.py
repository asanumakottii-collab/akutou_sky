"""Numerical resolution of the solver, kept apart from the physical Scene.

Altitude grids use four layers with their own step:
boundary layer 0-1 km, troposphere 1-10 km, stratosphere 10-30 km, upper 30-100 km.
"""
from dataclasses import asdict, dataclass, fields, replace
import json
import math
import os
from pathlib import Path
import numpy as np

LAYER_EDGES_M = (0., 1000., 10000., 30000., 100000.)
LEBEDEV_SIZES = (6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590)
HORIZONTAL_GROWTH = 2.0  # outer column spacing ratio; 1.5 -> 2.0 changed results <= 1.2%
MAX_SOURCE_WORK = 3.3e7  # incoming directions x sources x columns x levels x Stokes; measured ~10 GB


@dataclass(frozen=True)
class Numerics:
    """Immutable solver resolution. Presets are workloads, not accuracy grades.

    The diffuse (multiple-scattering) field is solved on a coarse grid; the
    observer single scattering, which carries the Earth's shadow edge, is
    integrated separately on a fine grid and replaces the coarse one.
    """
    wavelength_step_nm: float = 10.
    incoming_nodes: int = 302             # Lebedev directions integrating incoming diffuse light
    outgoing_nodes: int = 50              # Lebedev directions sampling the scattered source
    horizontal_step_deg: float = 0.5      # source columns near the observers
    horizontal_margin_deg: float = 3.     # fine columns beyond the observer range
    atmosphere_steps_m: tuple = (50., 500., 2000., 4000.)
    source_steps_m: tuple = (100., 2000., 4000., 8000.)
    single_scatter_steps_m: tuple = (25., 125., 250., 500.)
    legendre_moments: int = 0             # 0: 64 without fog, 256 with fog
    delta_m_order: int = 0                # 0: 16 without fog, 32 with fog (see docs)
    num_stokes: int = 1                   # 3 for I/Q/U (no fog)
    polarization_correction: bool = True  # scalar runs: add (vector - scalar) from a coarse-grid pair
    solar_refraction: bool = False        # see docs: SASKTRAN2 bends the beam but omits refractive dilution
    num_threads: int = 0                  # 0: all CPU cores

    def __post_init__(self):
        for name in ('atmosphere_steps_m', 'source_steps_m', 'single_scatter_steps_m'):
            steps = tuple(float(v) for v in getattr(self, name))
            if len(steps) != 4 or min(steps) < 1 or not all(map(math.isfinite, steps)):
                raise ValueError(f'{name} needs four finite steps >= 1 m (boundary, troposphere, stratosphere, upper)')
            object.__setattr__(self, name, steps)
        if not .5 <= self.wavelength_step_nm <= 40 or not math.isclose(400 / self.wavelength_step_nm,
                                                                         round(400 / self.wavelength_step_nm)):
            raise ValueError('wavelength_step_nm must divide 400 nm and lie in [0.5, 40]')
        for name in ('incoming_nodes', 'outgoing_nodes'):
            if getattr(self, name) not in LEBEDEV_SIZES:
                raise ValueError(f'{name} must be a Lebedev size {LEBEDEV_SIZES}')
        if not 0.05 <= self.horizontal_step_deg <= 5 or not 0 <= self.horizontal_margin_deg <= 20:
            raise ValueError('horizontal_step_deg must be 0.05..5 and horizontal_margin_deg 0..20')
        if self.delta_m_order not in (0, 2, 4, 8, 16, 32, 64, 128):
            raise ValueError('delta_m_order must be 0 (automatic) or a power of two in 2..128')
        if self.legendre_moments and self.delta_m_order and self.legendre_moments <= self.delta_m_order:
            raise ValueError('legendre_moments must exceed delta_m_order')
        if self.num_stokes not in (1, 3):
            raise ValueError('num_stokes must be 1 or 3')
        for name in ('solar_refraction', 'polarization_correction'):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f'{name} must be true or false')
        for name in ('incoming_nodes', 'outgoing_nodes', 'legendre_moments', 'delta_m_order', 'num_stokes', 'num_threads'):
            if isinstance(getattr(self, name), bool) or not isinstance(getattr(self, name), int):
                raise ValueError(f'{name} must be an integer')

    def moments(self, scene):
        return self.legendre_moments or (256 if scene.has_fog else 64)

    def delta_m(self, scene):
        # Fog droplets need a sharper truncated phase function; 32 matches the degree (29)
        # that the 302-direction incoming Lebedev rule integrates exactly.
        return self.delta_m_order or (32 if scene.has_fog else 16)

    def resolved(self, scene):
        """Settings with automatic (0) values replaced by what is used for this scene."""
        return self.replace(legendre_moments=self.moments(scene), delta_m_order=self.delta_m(scene))

    def threads(self):
        return self.num_threads or os.cpu_count() or 1

    def to_dict(self):
        return {k: list(v) if isinstance(v, tuple) else v for k, v in asdict(self).items()}

    def to_json(self, path=None):
        text = json.dumps(self.to_dict(), indent=2) + '\n'
        if path is not None:
            Path(path).write_text(text, encoding='utf-8')
        return text

    @classmethod
    def from_dict(cls, data):
        unknown = set(data) - {f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f'Unknown numerical settings: {sorted(unknown)}')
        return cls(**data)

    def replace(self, **changes):
        return replace(self, **changes)

    def refined(self, name, factor=2):
        """Return settings with one control made finer; used by convergence checks."""
        value = getattr(self, name)
        if isinstance(value, tuple):
            return self.replace(**{name: tuple(v / factor for v in value)})
        if name.endswith('_nodes'):
            larger = [n for n in LEBEDEV_SIZES if n >= value * factor]
            if not larger:
                raise ValueError(f'No Lebedev size above {value} x {factor}')
            return self.replace(**{name: larger[0]})
        if name in ('legendre_moments', 'delta_m_order'):
            if not value:
                raise ValueError(f'{name} is automatic (0); call resolved(scene) first')
            return self.replace(**{name: int(value * factor)})
        return self.replace(**{name: value / factor})


PRESETS = {
    'preview': Numerics(wavelength_step_nm=20., incoming_nodes=110, outgoing_nodes=26, horizontal_step_deg=1.,
                        atmosphere_steps_m=(100., 1000., 4000., 8000.),
                        source_steps_m=(250., 2000., 4000., 8000.),
                        single_scatter_steps_m=(50., 250., 500., 1000.)),
    'compare': Numerics(wavelength_step_nm=20.),
    'final': Numerics(),
}


def correction_numerics(numerics, scene):
    """Coarse settings for the polarization correction (same spectral grid, fog-compatible)."""
    preview = PRESETS['preview'].replace(incoming_nodes=194)
    bl = min(preview.atmosphere_steps_m[0], scene.fog_scale_height_m / 4) if scene.has_fog else preview.atmosphere_steps_m[0]
    sbl = min(preview.source_steps_m[0], scene.fog_scale_height_m / 2) if scene.has_fog else preview.source_steps_m[0]
    return preview.replace(wavelength_step_nm=numerics.wavelength_step_nm, polarization_correction=False,
                           solar_refraction=numerics.solar_refraction, num_threads=numerics.num_threads,
                           atmosphere_steps_m=(bl,) + preview.atmosphere_steps_m[1:],
                           source_steps_m=(sbl,) + preview.source_steps_m[1:])


def get_numerics(name_or_numerics):
    if isinstance(name_or_numerics, Numerics):
        return name_or_numerics
    if name_or_numerics in PRESETS:
        return PRESETS[name_or_numerics]
    path = Path(name_or_numerics)
    if path.suffix == '.json' and path.exists():
        return Numerics.from_dict(json.loads(path.read_text(encoding='utf-8')))
    raise ValueError(f'Unknown numerics {name_or_numerics!r}; choose from {tuple(PRESETS)} or a JSON path')


def layered_grid(steps):
    """Altitude levels 0..100 km with one step per layer (always includes layer edges)."""
    parts = [np.arange(lo, hi, step) for (lo, hi), step in zip(zip(LAYER_EDGES_M, LAYER_EDGES_M[1:]), steps)]
    return np.unique(np.r_[np.concatenate(parts), LAYER_EDGES_M[-1]])


def atmosphere_grid(numerics):
    return layered_grid(numerics.atmosphere_steps_m)


def source_grid(numerics):
    """Diffuse-source altitudes: layer midpoints, plus 1 m to resolve the surface (must be > 0)."""
    edges = layered_grid(numerics.source_steps_m)
    return np.r_[1., (edges[:-1] + edges[1:]) / 2]


def single_scatter_grid(numerics):
    return layered_grid(numerics.single_scatter_steps_m)


def horizontal_grid(numerics, depressions):
    """Source columns (degrees). An observer at depression d sits at angle -d.

    Uniform fine columns cover every observer +- margin; outside, the spacing
    grows by HORIZONTAL_GROWTH per column up to +-45 deg, so the far sunlit and
    night sides are represented without spending memory on them.
    """
    d = np.asarray(depressions, dtype=float)
    step = numerics.horizontal_step_deg
    lo = max(-45., math.floor((-d.max() - numerics.horizontal_margin_deg) / step) * step)
    hi = min(45., math.ceil((-d.min() + numerics.horizontal_margin_deg) / step) * step)
    columns = list(np.arange(lo, hi + step / 2, step))
    for sign, edge in ((-1, lo), (1, hi)):
        spacing = step
        while abs(edge) < 45:
            spacing *= HORIZONTAL_GROWTH
            edge = float(np.clip(edge + sign * spacing, -45, 45))
            if 45 - abs(edge) < spacing / 2:
                edge = 45. * sign
            columns.append(edge)
    return np.unique(np.round(columns, 6))


def workload(numerics, depressions, scene=None):
    """Memory driver of the diffuse field: rays traced from every source point in every incoming direction,
    with phase tables that grow with the delta-M order (measured: order 32 needs about 1.9x order 16)."""
    order = numerics.delta_m(scene) if scene is not None else (numerics.delta_m_order or 16)
    return (numerics.incoming_nodes * len(source_grid(numerics)) * len(horizontal_grid(numerics, depressions))
            * len(atmosphere_grid(numerics)) * numerics.num_stokes * max(1, order / 16) ** 0.9)


def check_workload(numerics, scene, depressions, *, native=False):
    """Validate fog resolution; refuse settings whose diffuse field cannot fit even for one depression."""
    if scene.has_fog:
        if numerics.num_stokes != 1 and not native:
            raise ValueError('Fog needs the scalar split delta-M solver: use num_stokes=1')
        if numerics.atmosphere_steps_m[0] > scene.fog_scale_height_m / 4:
            raise ValueError('Boundary-layer atmosphere step must be <= fog scale height / 4')
        if numerics.source_steps_m[0] > scene.fog_scale_height_m / 2:
            raise ValueError('Boundary-layer source step must be <= fog scale height / 2')
    worst = max(workload(numerics, [d], scene) for d in depressions)
    if worst > MAX_SOURCE_WORK:
        raise ValueError(f'Diffuse field work {worst:,.0f} exceeds {MAX_SOURCE_WORK:,.0f} even for one depression; '
                         'use coarser numerics')
    return worst


def depression_groups(numerics, depressions, scene=None):
    """Split depressions (by value) into groups whose diffuse field fits the memory budget.

    Each group is solved separately; returns lists of indices into depressions.
    """
    order = np.argsort(depressions)
    groups, current = [], []
    for i in order:
        if current and workload(numerics, [depressions[j] for j in current + [i]], scene) > MAX_SOURCE_WORK:
            groups.append(current)
            current = []
        current.append(int(i))
    return groups + [current]
