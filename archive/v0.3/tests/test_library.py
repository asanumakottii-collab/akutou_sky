"""Library contracts and real-solver regression, distinct from convergence tests."""
from dataclasses import FrozenInstanceError, replace
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import warnings

import numpy as np
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from akutou_sky import (ModelConfig, simulate, save_dataset, load_dataset,
                        export_results, ExperimentalPolarizationWarning, validate_sensitivity)
from akutou_sky.data import input_paths, resolve_cache_dir

FIXTURE = ROOT/'runs/v2_smoke.nc'

class LibraryTests(unittest.TestCase):
    def fixture(self):
        if not FIXTURE.exists():
            self.skipTest('Local solver fixture not available')
        return load_dataset(FIXTURE)

    def test_import_is_light_and_has_no_process_side_effects(self):
        code = '''import os,sys,json
before=dict(os.environ)
import akutou_sky
assert dict(os.environ)==before
assert 'sasktran2' not in sys.modules
assert 'matplotlib' not in sys.modules
print(akutou_sky.__version__)
'''
        with tempfile.TemporaryDirectory() as tmp:
            env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
            result=subprocess.run([sys.executable,'-c',code],cwd=tmp,env=env,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(list(Path(tmp).iterdir()),[])
            from akutou_sky import __version__
            self.assertEqual(result.stdout.strip(),__version__)

    def test_presets_serialization_and_immutability(self):
        config=ModelConfig.from_preset('preview',aod550=.04)
        self.assertEqual(config.nodes,14)
        self.assertEqual(config.aod550,.04)
        with self.assertRaises(FrozenInstanceError):config.aod550=.2
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'config.json'
            config.to_json(path)
            self.assertEqual(ModelConfig.from_json(path),config)
        for data in [[],None,{'unknown':1}]:
            with self.assertRaises(ValueError):ModelConfig.from_dict(data)
        with self.assertRaises(ValueError):ModelConfig.from_preset('unknown')

    def test_explicit_cache_overrides_environment_without_modifying_it(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ,{'AKUTOU_SKY_CACHE_DIR':tmp+'/environment'}):
            self.assertEqual(resolve_cache_dir(),(Path(tmp)/'environment').resolve())
            target=(Path(tmp)/'explicit').resolve()
            self.assertEqual(resolve_cache_dir(target),target)
            self.assertFalse(target.exists())

    def test_offline_missing_data_does_not_create_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'absent'
            with self.assertRaises(FileNotFoundError):input_paths(path,allow_download=False)
            self.assertFalse(path.exists())

    def test_background_error_precedes_transport(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'bad.csv'
            path.write_text('wavelength_nm,radiance_W_m2_sr_nm\n400,1\n700,1\n')
            with patch('akutou_sky.api.compute_radiance') as transport:
                with self.assertRaises(ValueError):simulate([0],background=path)
                transport.assert_not_called()

    def test_netcdf_export_bundle_and_overwrite_refusal(self):
        original=self.fixture().isel(depression=[0],azimuth=[1],elevation=[0,2,3])
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'run'
            self.assertEqual(export_results(original,out),out)
            for name in ['sky.nc','manifest.json','colors.csv','spectra.csv','bands.csv','summary.json','index.html']:
                self.assertTrue((out/name).is_file(),name)
            saved=load_dataset(out/'sky.nc')
            np.testing.assert_array_equal(saved.solar_stokes_radiance,original.solar_stokes_radiance)
            self.assertIn('total_XYZ',saved)
            self.assertNotIn('total_XYZ',original)
            before=(out/'sky.nc').read_bytes()
            with self.assertRaises(FileExistsError):export_results(original,out)
            self.assertEqual((out/'sky.nc').read_bytes(),before)
            with self.assertRaises(FileExistsError):save_dataset(original,out/'sky.nc')
            from akutou_sky import __version__
            self.assertEqual(json.loads((out/'manifest.json').read_text())['library_version'],__version__)

    def test_failed_save_preserves_previous_file_and_cleans_partial(self):
        dataset=self.fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'sky.nc';path.write_bytes(b'previous')
            with patch.object(xr.Dataset,'to_netcdf',side_effect=OSError('simulated write failure')):
                with self.assertRaises(OSError):save_dataset(dataset,path,overwrite=True)
            self.assertEqual(path.read_bytes(),b'previous')
            self.assertEqual(list(Path(tmp).iterdir()),[path])

    def test_real_solver_matches_saved_v2_and_accepts_simple_callback(self):
        reference=self.fixture()
        if any(not (ROOT/'model_cache'/key).exists() for key in
               ['climatology/fascode/std.atm','cross_sections/o3/dbm.nc','solar/solar_irradiance_hsrs_2022_11_30_extended.nc']):
            self.skipTest('Local physical inputs not available')
        config=ModelConfig.from_dict(json.loads(reference.attrs['config_json']))
        messages=[]; before=dict(os.environ)
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter('always',ExperimentalPolarizationWarning)
            actual=simulate(reference.depression.values,config=config,elevations=reference.elevation.values,
                azimuths=reference.azimuth.values,cache_dir=ROOT/'model_cache',allow_download=False,
                progress=messages.append)
        self.assertTrue(any(isinstance(item.message,ExperimentalPolarizationWarning) for item in recorded))
        np.testing.assert_allclose(actual.solar_stokes_radiance,reference.solar_stokes_radiance,rtol=1e-10,atol=1e-16)
        self.assertEqual(dict(os.environ),before)
        self.assertTrue(messages and all(isinstance(message,str) for message in messages))
        self.assertIn('total_XYZ',actual)

    def test_validation_in_memory_returns_report_without_result_files(self):
        reference=self.fixture()
        def candidate(ds,config,elevations,azimuths,progress,**kwargs):
            result=reference.sel(depression=ds,elevation=elevations,azimuth=azimuths).copy()
            result.attrs['config_json']=json.dumps(config.to_dict())
            return result
        with patch('akutou_sky.validation.compute_model',side_effect=candidate),tempfile.TemporaryDirectory() as tmp:
            with patch('akutou_sky.validation.Path.mkdir') as mkdir:
                report=validate_sensitivity(reference,checks=['iterations'],allow_download=False)
                mkdir.assert_not_called()
            self.assertTrue(report['all_requested_checks_complete'])
            self.assertTrue(report['passes_all_sampled_checks'])
            self.assertIsNone(report['input_sha256'])
        with self.assertRaises(ValueError):validate_sensitivity(reference,checks=['unknown'])

    def test_failed_refinement_cannot_be_reported_as_convergence(self):
        reference=self.fixture()
        original=json.loads(reference.attrs['config_json'])
        def candidate(ds,config,elevations,azimuths,progress,**kwargs):
            if config.wavelength_step_nm < original['wavelength_step_nm']:
                raise ValueError('Nonfinite intensity in refinement')
            result=reference.sel(depression=ds,elevation=elevations,azimuth=azimuths).copy()
            result.attrs['config_json']=json.dumps(config.to_dict())
            return result
        with patch('akutou_sky.validation.compute_model',side_effect=candidate):
            report=validate_sensitivity(reference,checks=['spectral','iterations'],allow_download=False)
        self.assertEqual(report['checks']['spectral']['status'],'failed')
        self.assertIn('Nonfinite intensity',report['checks']['spectral']['reason'])
        self.assertEqual(report['checks']['iterations']['status'],'completed')
        self.assertFalse(report['all_requested_checks_complete'])
        self.assertFalse(report['passes_all_sampled_checks'])

    def test_refinement_memory_guard_skips_before_starting_solver(self):
        with patch('akutou_sky.validation.compute_model') as transport:
            report=validate_sensitivity(self.fixture(),checks=['angular'],max_phase_terms=1)
            transport.assert_not_called()
        self.assertEqual(report['checks']['angular']['status'],'skipped')
        self.assertIn('Phase-table work',report['checks']['angular']['reason'])
        self.assertFalse(report['all_requested_checks_complete'])
        self.assertFalse(report['passes_all_sampled_checks'])

    def test_spectral_cache_accepts_verified_mie_tables_but_rejects_changed_inputs(self):
        from types import SimpleNamespace
        from akutou_sky.validation import _compatible_physical_inputs
        from akutou_sky.provenance import sha256
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'mie').mkdir()
            for name in ('gas.dat','mie/coarse.nc','mie/fine.nc'):(root/name).write_text(name)
            item=lambda name:dict(path=name,sha256=sha256(root/name))
            reference=SimpleNamespace(attrs={'physical_inputs_json':json.dumps([item('gas.dat'),item('mie/coarse.nc')])})
            candidate=SimpleNamespace(attrs={'physical_inputs_json':json.dumps([item('gas.dat'),item('mie/fine.nc')])})
            self.assertTrue(_compatible_physical_inputs(candidate,reference,root))
            (root/'mie/fine.nc').write_text('corrupted')
            self.assertFalse(_compatible_physical_inputs(candidate,reference,root))

    def test_plot_does_not_change_global_backend_or_style(self):
        try:
            import matplotlib
            from akutou_sky import plot_profiles
        except ImportError:
            self.skipTest('Optional matplotlib dependency not installed')
        reference=self.fixture()
        backend=matplotlib.get_backend()
        fonts=list(matplotlib.rcParams['font.family'])
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'figure.png'
            self.assertEqual(plot_profiles(reference,path,azimuth=0),path)
            self.assertGreater(path.stat().st_size,1000)
        self.assertEqual(matplotlib.get_backend(),backend)
        self.assertEqual(matplotlib.rcParams['font.family'],fonts)

if __name__=='__main__':unittest.main()
