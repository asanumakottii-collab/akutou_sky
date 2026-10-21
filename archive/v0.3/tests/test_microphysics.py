"""Water uptake, particle conservation, and surface-layer resolution contracts."""
from dataclasses import replace
import json
from pathlib import Path
import sys
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from akutou_sky import ModelConfig
from akutou_sky.config import altitude_grid, source_grid
from akutou_sky.microphysics import aerosol_water_uptake, WATER_INDEX
from akutou_sky.explorer import build_explorer, _compatible_signature

class WaterTests(unittest.TestCase):
    def test_cache_rejects_physics_changes_but_allows_viewer_changes(self):
        original = dict(sources={'model.py': 'transport', 'explorer.py': 'old-viewer'},
                        inputs={'solar': 'solar-hash'}, backend='1', library='0.2.0')
        viewer = {**original, 'sources': {**original['sources'], 'explorer.py': 'new-viewer'}}
        self.assertTrue(_compatible_signature(json.dumps(original), json.dumps(viewer)))
        for change in [dict(sources={**original['sources'], 'model.py': 'changed-transport'}),
                       dict(inputs={'solar': 'different-input'}), dict(backend='2')]:
            self.assertFalse(_compatible_signature(json.dumps(original), json.dumps({**original, **change})))

    def test_completed_explorer_rebuild_uses_verified_cache(self):
        directory = ROOT/'results/atmosphere_explorer'
        if not (directory/'manifest.json').exists():
            self.skipTest('Full local explorer not available')
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'explorer'
            shutil.copytree(directory/'calculations', out/'calculations')
            with patch('akutou_sky.model.compute_model', side_effect=AssertionError('Unexpected new transport calculation')) as transport:
                build_explorer(out, cache_dir=ROOT/'model_cache', allow_download=False, resume=True)
                transport.assert_not_called()
            self.assertTrue(json.loads((out/'manifest.json').read_text())['all_64_cases_complete'])
            for key in ('urban', 'rural', 'humid', 'dry'):
                self.assertIn('id="modelData"', (out/key/'index.html').read_text())
                self.assertNotIn('id="modelData"', (out/key/'baseline.html').read_text())

    def test_neutral_growth_retains_dry_optics_exactly(self):
        config = ModelConfig()
        self.assertEqual(aerosol_water_uptake(config),
                         (1., 100., complex(config.refractive_index_real, -config.refractive_index_imag)))
        self.assertEqual(aerosol_water_uptake(replace(config, relative_humidity=.90))[0], 1.)

    def test_kappa_growth_and_volume_mixing(self):
        config = ModelConfig(relative_humidity=.90, hygroscopicity_kappa=.30)
        growth, radius, index = aerosol_water_uptake(config)
        self.assertAlmostEqual(growth**3, 3.7, places=14)
        self.assertAlmostEqual(radius, 100*3.7**(1/3))
        dry = complex(1.45, -.005)
        self.assertAlmostEqual(index, (dry+2.7*WATER_INDEX)/3.7)
        low = aerosol_water_uptake(replace(config, relative_humidity=.30))
        self.assertGreater(growth, low[0])
        self.assertLess(abs(index.imag), abs(dry.imag))

    def test_rejects_saturation_and_incompatible_fog_grids(self):
        for values in [dict(relative_humidity=1.), dict(relative_humidity=-.1),
                       dict(hygroscopicity_kappa=2.), dict(hygroscopicity_kappa=-1.),
                       dict(aerosol_model='hg', relative_humidity=.8, hygroscopicity_kappa=.3),
                       dict(fog_aod550=.1), dict(delta_m_scaling='true'),
                       dict(delta_m_scaling=True, delta_m_order=64), dict(boundary_step_m=1.)]:
            with self.subTest(values=values), self.assertRaises(ValueError):
                ModelConfig.from_dict(values)
        valid = ModelConfig.from_preset('fog_preview', fog_aod550=.8)
        for values in [dict(boundary_step_m=100.), dict(source_boundary_step_m=200.),
                       dict(legendre_moments=64), dict(delta_m_scaling=False)]:
            with self.subTest(values=values), self.assertRaises(ValueError):
                replace(valid, **values).validate()

    def test_surface_layer_is_resolved_in_both_grids(self):
        cfg = ModelConfig.from_preset('fog_preview', fog_aod550=.8)
        h = altitude_grid(cfg)
        s = source_grid(cfg)
        np.testing.assert_array_equal(h[:21], np.arange(0, 1001, 50))
        np.testing.assert_array_equal(s[:10], np.arange(50, 1000, 100))
        self.assertEqual(h[-1], 100000)
        self.assertTrue(np.all(np.diff(h) > 0) and np.all(np.diff(s) > 0))

    def test_real_saved_case_preserves_particles_and_columns(self):
        base = ROOT/'results/atmosphere_explorer/calculations/urban'
        if not (base/'dust3_fog3.nc').exists():
            self.skipTest('Explorer has not been computed locally')
        with xr.open_dataset(base/'dust3_fog3.nc') as r:
            cfg = ModelConfig.from_dict(json.loads(r.attrs['config_json']))
            paths = [ROOT/'model_cache'/v['path'] for v in json.loads(r.attrs['physical_inputs_json'])
                     if v['path'].startswith('mie/')]
            with xr.open_dataset(paths[0]) as dry, xr.open_dataset(paths[1]) as wet:
                ratio = float(wet.xs_total.sel(wavelength_nm=550)/dry.xs_total.sel(wavelength_nm=550))
            self.assertAlmostEqual(r.attrs['aerosol_aod550_ambient'], cfg.aod550*ratio)
            self.assertAlmostEqual(np.trapezoid(r.aerosol_extinction_550, r.altitude), cfg.aod550*ratio)
            self.assertAlmostEqual(np.trapezoid(r.fog_extinction_550, r.altitude), cfg.fog_aod550)
            with xr.open_dataset(base/'dust3_fog0.nc') as clear:
                self.assertEqual(r.attrs['aerosol_aod550_ambient'], clear.attrs['aerosol_aod550_ambient'])
                self.assertFalse(np.allclose(r.total_XYZ, clear.total_XYZ))
                # Zenith rays are the same physical direction for all azimuths.
                z = r.solar_stokes_radiance.sel(elevation=90)
                # Forward-peaked fog / deep-twilight ray tracing is sensitive
                # to floating-point differences in the azimuth-degenerate ray.
                np.testing.assert_allclose(z[:,0], z[:,1], rtol=1e-4)
                np.testing.assert_allclose(z[:,0], z[:,2], rtol=1e-4)

if __name__ == '__main__':
    unittest.main()
