"""Conservation, grid nesting and batch-isolation checks for refined transport."""
from dataclasses import replace
from pathlib import Path
import json
import unittest
import numpy as np
from akutou_sky.config import ModelConfig, horizontal_grid, source_grid, single_scatter_grid
from akutou_sky.delta_m import truncated_coefficients
from akutou_sky.explorer import DEPRESSIONS

class RefinementTests(unittest.TestCase):
    def test_single_scattering_grid_is_nested_and_independent_of_diffuse_grid(self):
        cfg=ModelConfig.from_preset('fog_refined',fog_aod550=.2)
        moved=replace(cfg,ground_step_m=cfg.ground_step_m/1.1,
                      middle_step_m=cfg.middle_step_m/1.1,upper_step_m=cfg.upper_step_m/1.1,
                      boundary_step_m=cfg.boundary_step_m/1.1)
        np.testing.assert_array_equal(single_scatter_grid(cfg),single_scatter_grid(moved))
        fine=single_scatter_grid(replace(cfg,single_scatter_step_m=cfg.single_scatter_step_m/2))
        self.assertTrue(np.all(np.isin(single_scatter_grid(cfg),fine)))

    def test_replacement_preserves_the_saved_diffuse_component(self):
        import xarray as xr
        from unittest.mock import patch
        from akutou_sky.model import replace_single_scattering
        cfg=ModelConfig.from_preset('fog_refined',single_scatter_step_m=500.)
        dims=('depression','azimuth','elevation','wavelength','stokes')
        values=np.ones((1,1,1,1,1))
        reference=xr.Dataset({'solar_stokes_radiance':(dims,3*values),
                             'solar_single_scatter_radiance':(dims,values)},
            coords=dict(depression=[6],azimuth=[180],elevation=[45],wavelength=[550],stokes=['I'],altitude=[0,100000]),
            attrs=dict(config_json=json.dumps(cfg.to_dict()),transport_source_sha256='previous'))
        single=reference.copy(deep=True)
        single['solar_stokes_radiance']=(dims,2*values)
        with patch('akutou_sky.model.compute_model',return_value=single) as compute:
            refined=replace_single_scattering(reference,replace(cfg,single_scatter_step_m=250.))
        self.assertTrue(compute.call_args.kwargs['_single_only'])
        np.testing.assert_array_equal(refined.solar_stokes_radiance-refined.solar_single_scatter_radiance,
                                      reference.solar_stokes_radiance-reference.solar_single_scatter_radiance)
        with self.assertRaises(ValueError):
            replace_single_scattering(reference,replace(cfg,single_scatter_step_m=250.,aod550=.2))

    def test_surface_source_brackets_ground_observer(self):
        cfg=ModelConfig.from_preset('fog_refined')
        grid=source_grid(cfg)
        self.assertLess(grid[0],cfg.observer_altitude_m)
        self.assertGreater(grid[1],cfg.observer_altitude_m)
        self.assertTrue(np.all(np.diff(grid)>0))

    def test_saved_batch_isolates_atmospheres(self):
        from akutou_sky.io import load_dataset
        root=Path(__file__).resolve().parents[1]/'results/numerical_refinement'
        if not (root/'cache_isolation.json').exists():
            self.skipTest('Local batch isolation diagnostic is not installed')
        first=load_dataset(root/'final-humid-110.nc')
        again=load_dataset(root/'final-humid-repeated.nc')
        clear=load_dataset(root/'final-clear-reused.nc')
        self.assertEqual(json.loads(first.attrs['config_json']),json.loads(again.attrs['config_json']))
        self.assertTrue(again.attrs['transport_geometry_reused'])
        np.testing.assert_allclose(first.solar_stokes_radiance,again.solar_stokes_radiance,rtol=1e-8,atol=1e-18)
        self.assertFalse(np.allclose(first.solar_stokes_radiance,clear.solar_stokes_radiance))

    def test_split_has_correct_rayleigh_limit(self):
        from akutou_sky.model import compute_model
        cache=Path(__file__).resolve().parents[1]/'model_cache'
        if not (cache/'climatology/fascode/std.atm').exists():
            self.skipTest('Atmospheric input cache not installed')
        cfg=ModelConfig.from_preset('preview',aod550=0.,nodes=14,nsza=11,
            wavelength_step_nm=40.,delta_m_scaling=True,delta_m_order=8,
            relative_tolerance=1e-12,absolute_tolerance=1e-24,iterations=160)
        common=dict(elevations=[0,5,90],azimuths=[0,180],cache_dir=cache,allow_download=False)
        native=compute_model([-6,0,6],cfg,**common)
        split=compute_model([-6,0,6],replace(cfg,delta_m_method='split'),**common)
        np.testing.assert_allclose(native.solar_stokes_radiance,split.solar_stokes_radiance,rtol=1e-8,atol=1e-18)

    def test_delta_m_preserves_normalization_and_low_moments(self):
        order=16; g=.85; f=g**order
        ell=np.arange(33)
        beta=(2*ell+1)*g**ell
        actual=truncated_coefficients((beta/(1-f))[:,None,None], order)[:,0,0]
        expected=(beta[:order]-(2*ell[:order]+1)*f)/(1-f)
        np.testing.assert_allclose(actual,expected,rtol=1e-13,atol=1e-14)
        self.assertEqual(actual[0],1.)
        mu,w=np.polynomial.legendre.leggauss(16)
        self.assertAlmostEqual(np.dot(w,np.polynomial.legendre.legval(mu,actual))/2,1.,places=13)
        omega=.95; extinction=.003
        scaled_ext=extinction*(1-omega*f)
        scaled_omega=omega*(1-f)/(1-omega*f)
        self.assertAlmostEqual(scaled_ext*(1-scaled_omega),extinction*(1-omega),places=16)

    def test_zero_forward_fraction_does_not_change_rayleigh(self):
        beta=np.zeros((32,2,3));beta[0]=1;beta[2]=.5
        np.testing.assert_array_equal(truncated_coefficients(beta,16),beta[:16])

    def test_horizontal_refinement_keeps_existing_nodes(self):
        cfg=ModelConfig.from_preset('fog_preview', horizontal_sampling='near_horizon', nsza=41)
        a=horizontal_grid(cfg);b=horizontal_grid(replace(cfg,nsza=73))
        self.assertEqual(len(a),41);self.assertEqual(len(b),73)
        self.assertTrue(np.all(np.isin(a,b)))
        np.testing.assert_allclose(a,-a[::-1])
        self.assertTrue(np.all(np.isin(DEPRESSIONS,a)))

    def test_depressions_half_degree_both_endpoints(self):
        self.assertEqual(len(DEPRESSIONS),25)
        np.testing.assert_array_equal(DEPRESSIONS,np.arange(-6,6.01,.5))

    def test_split_cannot_silently_disable_scaling_or_accept_vector(self):
        for change in [dict(delta_m_scaling=False),dict(num_stokes=3),dict(delta_m_method='typo')]:
            with self.subTest(change=change),self.assertRaises(ValueError):
                ModelConfig.from_preset('fog_preview',**{**dict(delta_m_method='split'),**change})

if __name__=='__main__':unittest.main()
