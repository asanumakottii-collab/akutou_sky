"""Accuracy checks: numerical convergence and comparison with measured twilight illuminance.

Convergence: refine one numerical control at a time and compare what the
projector actually uses, the RevC band-mean XYZ (luminance and xy), plus
individual rays for information.

Illuminance: integrate the computed sky over the hemisphere (plus the direct
Sun when it is above the horizon) and compare with the clear-sky illuminance
model of USNO Circular 171 (Janiczek & DeYoung 1987), a fit to measurements
by Bond & Henderson (1963). It is an empirical clear-sky average, so a factor
of order 1.5 reflects real atmospheric variability rather than model error.
"""
import json
import numpy as np
import sasktran2 as sk

from .color import add_color, band_xyz, chromaticity, cmf_weights
from .numerics import get_numerics
from .scene import get_scene
from .solver import EARTH_RADIUS_M, render

CONVERGENCE_DEPRESSIONS = (-6., -3., 0., 3., 6., 9., 12.)
CONVERGENCE_ELEVATIONS = (0., .5, 1., 2., 3., 5., 8., 12., 20., 30., 45., 60.)
CONVERGENCE_AZIMUTHS = (0., 90., 180.)
REFINEMENTS = {  # control -> refinement factor
    'incoming_nodes': 1.4, 'outgoing_nodes': 1.4, 'horizontal_step_deg': 2, 'horizontal_margin_deg': .5, 'atmosphere_steps_m': 2,
    'source_steps_m': 2, 'single_scatter_steps_m': 2, 'wavelength_step_nm': 2, 'legendre_moments': 2,
    'delta_m_order': 2,
}
# 'polarization' compares the scalar solution plus its coarse-grid correction with a full I/Q/U solve.
CONTROLS = (*REFINEMENTS, 'polarization')
BAND_Y_TOLERANCE, BAND_XY_TOLERANCE = 0.03, 0.003


def _differences(reference, candidate):
    def summary(xyz_a, xyz_b):
        y = np.abs(xyz_b[..., 1] - xyz_a[..., 1]) / np.maximum(xyz_a[..., 1], 1e-30)
        xy = np.linalg.norm(chromaticity(xyz_b) - chromaticity(xyz_a), axis=-1)
        return y, xy
    by, bxy = summary(band_xyz(reference), band_xyz(candidate))
    ry, rxy = summary(reference.XYZ.values, candidate.XYZ.values)
    low = reference.elevation.values < 1
    per_depression = [dict(depression=float(d), band_max_Y=float(by[i].max()), band_max_xy=float(bxy[i].max()),
                           ray_max_Y=float(ry[i][..., ~low].max()), ray_max_xy=float(rxy[i][..., ~low].max()))
                      for i, d in enumerate(reference.depression.values)]
    return dict(band_max_Y=float(by.max()), band_max_xy=float(bxy.max()),
                ray_max_Y_above_1deg=float(ry[..., ~low].max()), ray_max_xy_above_1deg=float(rxy[..., ~low].max()),
                ray_max_Y_below_1deg=float(ry[..., low].max()), ray_max_xy_below_1deg=float(rxy[..., low].max()),
                passes=bool(by.max() <= BAND_Y_TOLERANCE and bxy.max() <= BAND_XY_TOLERANCE),
                per_depression=per_depression)


def convergence(scene='standard', numerics='final', controls=None, *, depressions=CONVERGENCE_DEPRESSIONS,
                cache_dir=None, allow_download=True, progress=None):
    """Refine each control in turn; return a JSON-ready report."""
    scene, numerics = get_scene(scene), get_numerics(numerics)
    controls = [c for c in (controls or CONTROLS)
                if not (c in ('delta_m_order', 'polarization') and numerics.num_stokes == 3)]
    if scene.has_fog and 'polarization' in controls:
        controls.remove('polarization')  # a full polarized solve is not available with fog
    run = lambda n: add_color(render(depressions, scene, n, CONVERGENCE_ELEVATIONS, CONVERGENCE_AZIMUTHS,
                                     cache_dir=cache_dir, allow_download=allow_download))
    base = run(numerics)
    report = dict(scene=scene.to_dict(), numerics=numerics.to_dict(),
                  criterion=f'band-mean luminance <= {BAND_Y_TOLERANCE:.0%} and xy <= {BAND_XY_TOLERANCE}',
                  base_seconds=base.attrs['elapsed_seconds'], checks={})
    for control in controls:
        if control == 'polarization':
            refined = numerics.replace(num_stokes=3, polarization_correction=False)
            entry = dict(base='scalar + coarse-grid correction', refined='full I/Q/U')
        else:
            resolved = numerics.resolved(scene)
            refined = resolved.refined(control, REFINEMENTS[control])
            entry = dict(base=resolved.to_dict()[control], refined=refined.to_dict()[control])
        try:
            entry.update(_differences(base, run(refined)))
        except (ValueError, MemoryError) as exc:
            entry.update(status='not run', reason=str(exc))
        report['checks'][control] = entry
        if progress:
            progress(f'{control}: {entry["base"]} -> {entry["refined"]}: '
                     + (f'band Y {entry["band_max_Y"]:.2%}, xy {entry["band_max_xy"]:.4f}'
                        if 'passes' in entry else entry['reason']))
    report['all_pass'] = all(c.get('passes', False) for c in report['checks'].values())
    return report


def usno_illuminance(solar_altitude_deg):
    """Clear-sky horizontal illuminance (lux), USNO Circular 171 (atmos and refr routines)."""
    h = np.atleast_1d(np.asarray(solar_altitude_deg, dtype=float))
    dr = np.pi / 180
    above = h >= -5 / 6
    h = np.where(above, h + 1 / np.tan((h + 8.6 / (h + 4.42)) * dr) / 60, h)
    x = 753.66156
    u = np.sin(h * dr)
    s = np.arcsin(x * np.cos(h * dr) / (x + 1))
    m = x * (np.cos(s) - u) + np.cos(s)
    return 133775 * (np.exp(-.21 * m) * u + .0289 * np.exp(-.042 * m) * (1 + (h + 90) * u / 57.29577951))


def _direct_illuminance(ds, scene, numerics, cache_dir):
    """Direct-Sun horizontal illuminance for depressions < 0 (point Sun, geometric path)."""
    above = ds.depression.values < 0
    result = np.zeros(len(ds.depression))
    if not above.any():
        return result
    cfg = sk.Config()
    cfg.num_threads = numerics.threads()
    cfg.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    cfg.output_los_optical_depth = True
    cfg.num_singlescatter_moments = numerics.moments(scene)
    height = ds.altitude.values
    geo = sk.Geometry1D(cos_sza=.5, solar_azimuth=0., earth_radius_m=EARTH_RADIUS_M, altitude_grid_m=height,
                        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                        geometry_type=sk.GeometryType.Spherical)
    view = sk.ViewingGeometry()
    for d in ds.depression.values[above]:  # ray from the observer toward the Sun
        view.add_ray(sk.SolarAnglesObserverLocation(.5, 0., np.sin(np.deg2rad(-d)), scene.observer_altitude_m))
    from .optics import Optics
    from .data import input_paths, resolve_cache_dir
    edges = np.r_[ds.wavelength_bounds.values[:, 0], ds.wavelength_bounds.values[-1, 1]]
    cache = resolve_cache_dir(cache_dir)
    optics = Optics(scene, ds.wavelength.values, edges, input_paths(cache, allow_download=False), cache,
                    numerics.moments(scene), numerics.threads())
    atmo = sk.Atmosphere(geo, cfg, wavelengths_nm=ds.wavelength.values, calculate_derivatives=False)
    optics.add_to(atmo, height, 0, len(ds.wavelength))
    tau = sk.Engine(cfg, geo, view).calculate_radiance(atmo).los_optical_depth.transpose('los', 'wavelength').values
    y_bar = cmf_weights(ds.wavelength_bounds.values)[:, 1]
    lux = 683 * (ds.solar_irradiance.values * np.exp(-tau)) @ y_bar
    result[above] = lux * np.sin(np.deg2rad(-ds.depression.values[above]))
    return result


def illuminance(scene='standard', numerics='final', depressions=(-6., -3., 0., 2., 4., 6., 8., 10., 12.), *,
                cache_dir=None, allow_download=True):
    """Model horizontal illuminance vs USNO Circular 171; also zenith chromaticity."""
    scene, numerics = get_scene(scene), get_numerics(numerics)
    elevations = [0., 1., 2., 3., 5., 7., 10., 15., 20., 30., 40., 50., 60., 70., 80., 90.]
    azimuths = [0., 30., 60., 90., 120., 150., 180.]
    ds = add_color(render(depressions, scene, numerics, elevations, azimuths,
                          cache_dir=cache_dir, allow_download=allow_download))
    y = ds.XYZ.sel(tristimulus='Y').values  # (d, a, e)
    e, a = np.deg2rad(ds.elevation.values), np.deg2rad(ds.azimuth.values)
    over_e = np.trapezoid(y * (np.sin(e) * np.cos(e)), e, axis=-1)
    diffuse = 2 * np.trapezoid(over_e, a, axis=-1)  # mirror symmetry about the solar plane
    if np.any(ds.depression.values < 0):
        direct = _direct_illuminance(ds, scene, numerics, cache_dir)
    else:
        direct = np.zeros_like(diffuse)
    reference = usno_illuminance(-ds.depression.values)
    zenith = ds.xy.sel(elevation=90.).isel(azimuth=0).values
    rows = [dict(depression=float(d), diffuse_lux=float(f), direct_lux=float(s), model_lux=float(f + s),
                 usno_lux=float(r), ratio=float((f + s) / r), zenith_x=float(z[0]), zenith_y=float(z[1]),
                 zenith_luminance=float(ds.XYZ.sel(elevation=90., tristimulus='Y').isel(azimuth=0).values[i]))
            for i, (d, f, s, r, z) in enumerate(zip(ds.depression.values, diffuse, direct, reference, zenith))]
    return dict(scene=scene.to_dict(), numerics=numerics.to_dict(),
                reference='USNO Circular 171 (Janiczek & DeYoung 1987) clear-sky illuminance model',
                rows=rows)


def save_report(report, path):
    from pathlib import Path
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
