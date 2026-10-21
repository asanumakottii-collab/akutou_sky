"""Small real-solver checks of source scaling, symmetry, refraction and optional layers."""
from dataclasses import replace
import json
import numpy as np
from pathlib import Path
from akutou_sky import compute_radiance, ModelConfig
ROOT = Path(__file__).resolve().parents[1]

def main():
    c=ModelConfig(nodes=14,nsza=15,ground_step_m=1000.,middle_step_m=2000.,upper_step_m=4000.,
                  source_ground_step_m=1000.,source_middle_step_m=2000.,source_upper_step_m=4000.,
                  wavelength_step_nm=40.,wavelength_chunk=1,num_stokes=3,iterations=60)
    def run(config):
        return compute_radiance([0.,6.,12.],config=config,elevations=[0.,5.,30.,90.],
                                azimuths=[-90.,0.,90.,180.],cache_dir=ROOT/'model_cache',allow_download=False)
    checks = {}
    def compare(name, actual, expected, rtol):
        actual, expected = np.asarray(actual), np.asarray(expected)
        checks[name] = {
            'passed': bool(np.allclose(actual, expected, rtol=rtol, atol=1e-16)),
            'rtol': rtol, 'atol': 1e-16,
            'max_relative_difference': float(np.max(abs(actual-expected)/np.maximum(abs(expected),1e-16))),
            'max_absolute_difference': float(np.max(abs(actual-expected))),
        }
    a=run(c)
    distant=run(replace(c,solar_distance_au=2.))
    compare('inverse_square_irradiance', distant.solar_stokes_radiance, a.solar_stokes_radiance/4, 1e-5)
    intensity=a.solar_stokes_radiance.sel(stokes='I')
    compare('azimuth_mirror_symmetry', intensity.sel(azimuth=-90), intensity.sel(azimuth=90), 1e-8)
    zenith=intensity.sel(elevation=90).values
    compare('zenith_azimuth_invariance', zenith, np.broadcast_to(zenith[:,:1], zenith.shape), 1e-8)
    straight=run(replace(c,solar_refraction=False))
    bend=float(np.max(abs(straight.solar_stokes_radiance.sel(stokes='I').values/intensity.values-1)))
    checks['solar_refraction_changes_solution'] = {'passed': bend>1e-4, 'max_relative_difference': bend}
    haze=run(replace(c,stratospheric_aod550=.01))
    checks['stratospheric_aerosol'] = {'passed': bool(np.isfinite(haze.solar_stokes_radiance).all()
        and (haze.solar_stokes_radiance.sel(stokes='I') >= 0).all()
        and not np.allclose(haze.solar_stokes_radiance.values,a.solar_stokes_radiance.values,rtol=1e-4,atol=1e-16))}
    hg=run(replace(c,aerosol_model='hg'))
    checks['hg_with_polarization'] = {'passed': bool(np.isfinite(hg.solar_stokes_radiance).all()
        and (hg.solar_stokes_radiance.sel(stokes='I') >= 0).all())}
    report={'real_solver_cases':5, 'all_checks_completed': True, 'checks': checks,
            'all_checks_passed': all(item['passed'] for item in checks.values()),
            'config': c.to_dict(),
            'scope':'coarse-grid physical invariants and option smoke checks; not a convergence or observation test',
            'interpretation':'Failed symmetry checks are retained, not corrected by averaging. Polarized off-plane results are experimental; root cause not isolated.'}
    (ROOT/'physics_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    a.to_netcdf(ROOT/'runs/v2_smoke.nc')
    print(json.dumps(report,indent=2),flush=True)
    if not report['all_checks_passed']:
        raise SystemExit(1)

if __name__=='__main__':main()
