"""Physical atmosphere of one sky: what the air contains, not how it is solved."""
from dataclasses import asdict, dataclass, fields, replace
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class Scene:
    """Immutable physical inputs. Lengths in metres, AOD at 550 nm.

    aod550 is the DRY tropospheric aerosol column; with relative_humidity and
    hygroscopicity_kappa > 0 the particles grow and the ambient AOD increases.
    Fog is an independent liquid-droplet layer. Background stratospheric
    aerosol is on by default because it shapes the purple light of twilight.
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
    fog_aod550: float = 0.
    fog_scale_height_m: float = 200.
    fog_median_radius_um: float = 5.
    fog_mode_width: float = 1.3
    ozone_du: float = 300.
    surface_albedo: float = 0.1
    observer_altitude_m: float = 2.
    solar_distance_au: float = 1.

    def __post_init__(self):
        for f in fields(self):
            value = getattr(self, f.name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f'{f.name} must be a finite number')
        for name in ('aod550', 'stratospheric_aod550', 'fog_aod550', 'ozone_du',
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
            (0 <= self.relative_humidity <= .95, 'relative_humidity is a fraction in [0, 0.95]'),
            (self.hygroscopicity_kappa <= 1.5, 'hygroscopicity_kappa must be <= 1.5'),
            (not self.relative_humidity * self.hygroscopicity_kappa or self.median_radius_nm >= 40,
             'hygroscopic growth needs a dry median radius >= 40 nm'),
            (0 < self.stratospheric_center_m < 100000, 'stratospheric_center_m must be in (0, 100 km)'),
            (50 <= self.fog_scale_height_m <= 1000, 'fog_scale_height_m must be 50..1000 m'),
            (.5 <= self.fog_median_radius_um <= 10, 'fog_median_radius_um must be 0.5..10'),
            (1 < self.fog_mode_width <= 1.5, 'fog_mode_width must be in (1, 1.5]'),
            (0 <= self.surface_albedo <= 1, 'surface_albedo must be 0..1'),
            (0 <= self.observer_altitude_m < 10000, 'observer_altitude_m must be in [0, 10 km)'),
        ]
        for ok, message in checks:
            if not ok:
                raise ValueError(message)

    @property
    def has_fog(self):
        return self.fog_aod550 > 0

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
    'humid': Scene(aod550=0.10, relative_humidity=0.9, hygroscopicity_kappa=0.3, fog_aod550=0.20),
    'dry': Scene(aod550=0.10, relative_humidity=0.3, hygroscopicity_kappa=0.3),
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
