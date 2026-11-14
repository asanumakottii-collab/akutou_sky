"""Physical atmosphere of one sky: what the air contains, not how it is solved."""
from dataclasses import asdict, dataclass, fields, replace
import json
import math
from pathlib import Path

# Per-parameter limits, each checked for convergence on its own (docs/精度検証.md); combinations at
# several limits at once, and the tropospheric aod550, are not capped. Heavier aerosol loads need more
# scattered-light directions than 16 GB allows to converge: AOD 0.10 / 200 nm and RH 0.95 failed.
MAX_RELATIVE_HUMIDITY = 0.7
MAX_STRATOSPHERIC_AOD = 0.02
STRATOSPHERIC_RADIUS_NM = (40., 150.)
STRATOSPHERIC_WIDTH = (1.2, 1.6)


@dataclass(frozen=True)
class Scene:
    """Immutable physical inputs. Lengths in metres, AOD at 550 nm.

    aod550 is the DRY tropospheric aerosol column; with relative_humidity and
    hygroscopicity_kappa > 0 the particles grow and the ambient AOD increases.
    The stratospheric sulfate layer (Gaussian in altitude, lognormal sizes) is
    on at background level by default and makes the purple light of twilight;
    raising its AOD and particle size, as after a volcanic eruption, makes it
    stronger. Fog droplets are not modelled, and humidity and stratospheric
    loading are capped: beyond these the solution does not converge within
    16 GB of memory (docs/精度検証.md).
    """
    aod550: float = 0.10
    aerosol_scale_height_m: float = 1500.
    median_radius_nm: float = 100.
    mode_width: float = 1.6
    refractive_index_real: float = 1.45
    refractive_index_imag: float = 0.005  # absorption; backend convention is n - ik
    relative_humidity: float = 0.
    hygroscopicity_kappa: float = 0.
    stratospheric_aod550: float = 0.005
    stratospheric_center_m: float = 20000.
    stratospheric_width_m: float = 4000.
    stratospheric_median_radius_nm: float = 80.
    stratospheric_mode_width: float = 1.6
    ozone_du: float = 300.
    surface_albedo: float = 0.1
    observer_altitude_m: float = 2.
    solar_distance_au: float = 1.

    def __post_init__(self):
        for f in fields(self):
            value = getattr(self, f.name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f'{f.name} must be a finite number')
        for name in ('aod550', 'stratospheric_aod550', 'ozone_du',
                     'refractive_index_imag', 'hygroscopicity_kappa'):
            if getattr(self, name) < 0:
                raise ValueError(f'{name} must be nonnegative')
        for name in ('aerosol_scale_height_m', 'stratospheric_width_m', 'solar_distance_au'):
            if getattr(self, name) <= 0:
                raise ValueError(f'{name} must be positive')
        checks = [
            (1 <= self.median_radius_nm <= 1000, 'median_radius_nm must be 1..1000'),
            (1 < self.mode_width <= 2.5, 'mode_width must be in (1, 2.5]'),
            (1 <= self.refractive_index_real <= 3, 'refractive_index_real must be 1..3'),
            (0 <= self.relative_humidity <= MAX_RELATIVE_HUMIDITY,
             f'relative_humidity is a fraction in [0, {MAX_RELATIVE_HUMIDITY}] (converged range)'),
            (self.hygroscopicity_kappa <= 1.5, 'hygroscopicity_kappa must be <= 1.5'),
            (not self.relative_humidity * self.hygroscopicity_kappa or self.median_radius_nm >= 40,
             'hygroscopic growth needs a dry median radius >= 40 nm'),
            (0 < self.stratospheric_center_m < 100000, 'stratospheric_center_m must be in (0, 100 km)'),
            (self.stratospheric_aod550 <= MAX_STRATOSPHERIC_AOD,
             f'stratospheric_aod550 must be <= {MAX_STRATOSPHERIC_AOD} (converged range)'),
            (STRATOSPHERIC_RADIUS_NM[0] <= self.stratospheric_median_radius_nm <= STRATOSPHERIC_RADIUS_NM[1],
             'stratospheric_median_radius_nm must be %g..%g (converged range)' % STRATOSPHERIC_RADIUS_NM),
            (STRATOSPHERIC_WIDTH[0] <= self.stratospheric_mode_width <= STRATOSPHERIC_WIDTH[1],
             'stratospheric_mode_width must be %g..%g (converged range)' % STRATOSPHERIC_WIDTH),
            (0 <= self.surface_albedo <= 1, 'surface_albedo must be 0..1'),
            (0 <= self.observer_altitude_m < 10000, 'observer_altitude_m must be in [0, 10 km)'),
        ]
        for ok, message in checks:
            if not ok:
                raise ValueError(message)

    def to_dict(self):
        return asdict(self)

    def to_json(self, path=None):
        text = json.dumps(self.to_dict(), indent=2) + '\n'
        if path is not None:
            Path(path).write_text(text, encoding='utf-8')
        return text

    @classmethod
    def from_dict(cls, data):
        unknown = set(data) - {f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f'Unknown scene settings: {sorted(unknown)}')
        return cls(**data)

    @classmethod
    def from_json(cls, path):
        return cls.from_dict(json.loads(Path(path).read_text(encoding='utf-8')))

    def replace(self, **changes):
        return replace(self, **changes)


# Named, idealized skies. Comparison examples, not measured regional climatologies.
SCENES = {
    'clear': Scene(aod550=0.05),
    'standard': Scene(),
    'hazy': Scene(aod550=0.30, relative_humidity=0.7, hygroscopicity_kappa=0.3),
    'urban': Scene(aod550=0.30, relative_humidity=0.6, hygroscopicity_kappa=0.2,
                   refractive_index_real=1.50, refractive_index_imag=0.010, aerosol_scale_height_m=1200.),
    'rural': Scene(aod550=0.03, relative_humidity=0.6, hygroscopicity_kappa=0.3,
                   refractive_index_imag=0.003),
    'dry': Scene(aod550=0.10, relative_humidity=0.3, hygroscopicity_kappa=0.3),
    # Strong purple light: standard troposphere, more and larger stratospheric sulfate after an eruption.
    'eruption': Scene(stratospheric_aod550=0.02, stratospheric_median_radius_nm=150., stratospheric_mode_width=1.5),
}


def get_scene(name_or_path):
    """Return a named scene or read a scene JSON file."""
    if isinstance(name_or_path, Scene):
        return name_or_path
    if name_or_path in SCENES:
        return SCENES[name_or_path]
    path = Path(name_or_path)
    if path.suffix == '.json' and path.exists():
        return Scene.from_json(path)
    raise ValueError(f'Unknown scene {name_or_path!r}; choose from {tuple(SCENES)} or a JSON path')
