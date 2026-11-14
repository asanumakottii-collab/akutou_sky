import csv
from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import xarray as xr
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from akutou_sky.config import ModelConfig, altitude_grid, source_grid, spectral_grid, integrate_linear, checked_angles
from akutou_sky.colorimetry import xyz_from_radiance, display_values, band_mean, background_spectrum, add_colorimetry
from akutou_sky.model import normalized_profile, read_fascode
from akutou_sky.validation import comparison

class NumericalTests(unittest.TestCase):
    def test_piecewise_integral_exact_partial_bins(self):
        # y=2x+3, including edges inside source intervals and both endpoints.
        x=np.array([0.,1.,3.,10.]); edges=np.array([0.,.2,2.,10.])
        np.testing.assert_allclose(integrate_linear(x,2*x+3,edges),np.diff(edges**2+3*edges),rtol=1e-14)

    def test_nonuniform_spectrum_integral(self):
        x=np.array([380.,400.,550.,700.,780.]); y=np.array([0.,1.,5.,1.,0.])
        bins=np.array([380.,450.,600.,780.])
        self.assertAlmostEqual(integrate_linear(x,y,bins).sum(),np.trapezoid(y,x))

    def test_no_spectral_extrapolation(self):
        with self.assertRaises(ValueError): integrate_linear([400,700],[1,1],[380,780])

    def test_bins_cover_visible_range_without_double_counting(self):
        w,edges=spectral_grid(5)
        self.assertEqual(len(w),80)
        self.assertEqual(edges[0],380)
        self.assertEqual(edges[-1],780)
        self.assertEqual(np.diff(edges).sum(),400)

    def test_flat_spectrum_photopic_luminance(self):
        w,edges=spectral_grid(5)
        xyz=xyz_from_radiance(np.ones(len(w)),w,np.c_[edges[:-1],edges[1:]])
        # Integral of CIE photopic V(lambda) over visible range is ~106.86 nm.
        self.assertAlmostEqual(xyz[1]/683,106.86,delta=.05)
        self.assertAlmostEqual(xyz[0]/xyz.sum(),1/3,delta=.001)

    def test_flat_spectrum_invariant_to_bin_width(self):
        outputs=[]
        for step in [2.5,5,10,20,40]:
            w,b=spectral_grid(step)
            outputs.append(xyz_from_radiance(np.ones(len(w)),w,np.c_[b[:-1],b[1:]]))
        for xyz in outputs[1:]: np.testing.assert_allclose(xyz,outputs[0],rtol=1e-12)

    def test_black_is_black_undefined_chromaticity(self):
        xy,rgb,gamut=display_values(np.zeros((2,3)))
        self.assertTrue(np.isnan(xy).all())
        self.assertTrue((rgb==0).all())
        self.assertFalse(gamut.any())

    def test_xyz_scaling_preserves_preview(self):
        xyz=np.array([[.4,.3,.2],[.1,.2,.5]])
        a=display_values(xyz); b=display_values(xyz*1e-8)
        np.testing.assert_allclose(a[0],b[0]);np.testing.assert_allclose(a[1],b[1])

    def test_band_constant_field_and_non_sampled_boundaries(self):
        e=np.array([0.,2.,10.,30.,90.]); xyz=np.broadcast_to([2.,3.,4.],(5,3))
        np.testing.assert_allclose(band_mean(e,xyz,5,15),[2,3,4])
        with self.assertRaises(ValueError):band_mean(e,xyz,-1,10)

    def test_band_uses_photon_amount_not_average_rgb(self):
        e=np.array([0.,5.,10.]); xyz=np.array([[100.,0.,0.],[1.,1.,1.],[0.,0.,1.]])
        mean=band_mean(e,xyz,0,10)
        self.assertGreater(mean[0],20*mean[2])

    def test_columns_recover_requested_aod_and_ozone(self):
        h=altitude_grid(ModelConfig())
        profile=normalized_profile(h,np.exp(-h/1500),.123)
        self.assertAlmostEqual(np.trapezoid(profile,h),.123,places=14)
        ozone=normalized_profile(h,np.exp(-.5*((h-24000)/7000)**2),300*2.687e20)
        self.assertAlmostEqual(np.trapezoid(ozone,h)/2.687e20,300,places=10)

    def test_fascode_units_and_profile_lengths(self):
        if not (ROOT/'model_cache/climatology/fascode/std.atm').exists():
            self.skipTest('Local FASCODE input not available')
        data=read_fascode(ROOT/'model_cache/climatology/fascode/std.atm')
        self.assertEqual(len(data['HGT']),len(data['O3']))
        self.assertEqual(data['HGT'][0],0)
        self.assertGreater(data['HGT'][-1],100)
        self.assertTrue((data['O3']>=0).all())

    def test_configuration_rejects_bad_input(self):
        for changes in [{'aod550':-1},{'aod550':float('nan')},{'nsza':1},{'nodes':25},
                        {'num_stokes':2},{'mode_width':1},{'wavelength_step_nm':7},
                        {'solar_refraction':'false'},{'unknown':1},{'surface_albedo':1.1}]:
            with self.subTest(changes=changes),self.assertRaises(ValueError):ModelConfig.from_dict(changes)

    def test_angles_allow_morning_order_reject_duplicates(self):
        np.testing.assert_array_equal(checked_angles([18,6,0,-5],'d',-18,24),[18,6,0,-5])
        for v in [[0,0],[float('nan')],[],[25]]:
            with self.assertRaises(ValueError):checked_angles(v,'d',-18,24)

    def test_large_source_grid_rejected_before_solver(self):
        with self.assertRaisesRegex(ValueError,'exceeding budget'):
            ModelConfig(nodes=50,nsza=91,num_stokes=3).validate()
        self.assertLess(len(source_grid(ModelConfig())),len(altitude_grid(ModelConfig())))

    def test_export_handles_arbitrary_axes_and_background(self):
        path=ROOT/'runs/v2_smoke.nc'
        if not path.exists():self.skipTest('Run solver smoke first')
        from akutou_sky.export import export
        r=xr.load_dataset(path).isel(depression=[1,0],azimuth=[1],elevation=[0,2,3])
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            bg=out/'background.csv'
            bg.write_text('wavelength_nm,radiance_W_m2_sr_nm\n380,0.00001\n780,0.00001\n')
            r=add_colorimetry(r,bg)
            np.testing.assert_allclose(r.total_radiance-r.solar_stokes_radiance.sel(stokes='I'),1e-5,atol=1e-16)
            r.to_netcdf(out/'sky.nc')
            with xr.open_dataset(out/'sky.nc') as saved:
                np.testing.assert_array_equal(saved.total_XYZ,r.total_XYZ)
            export(r,out)
            with (out/'colors.csv').open(encoding='utf-8-sig') as f:
                self.assertEqual(len(list(csv.DictReader(f))),6)
            self.assertIn('application/json',(out/'index.html').read_text())

    def test_background_requires_full_coverage_and_positive_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'bg.csv'
            p.write_text('wavelength_nm,radiance_W_m2_sr_nm\n380,0.1\n780,0.1\n')
            w,b=spectral_grid(10)
            np.testing.assert_allclose(background_spectrum(p,np.c_[b[:-1],b[1:]]),.1)
            p.write_text('wavelength_nm,radiance_W_m2_sr_nm\n400,0.1\n700,0.1\n')
            with self.assertRaises(ValueError):background_spectrum(p,np.c_[b[:-1],b[1:]])

    def test_real_solver_smoke_geometry_and_color_export(self):
        path=ROOT/'runs/v2_smoke.nc'
        if not path.exists():self.skipTest('Run solver smoke first')
        r=xr.load_dataset(path); r=add_colorimetry(r)
        intensity=r.solar_stokes_radiance.sel(stokes='I')
        self.assertTrue(np.isfinite(intensity).all())
        self.assertTrue((intensity>=0).all())
        # Azimuth is degenerate at zenith; all views are the same physical ray.
        zenith=intensity.sel(elevation=90).values
        np.testing.assert_allclose(zenith[:,0,:],zenith[:,1,:],rtol=1e-9)
        np.testing.assert_allclose(zenith[:,0,:],zenith[:,2,:],rtol=1e-9)
        np.testing.assert_allclose(r.solar_XYZ,r.total_XYZ)
        # The transport intensity must be linear in source irradiance.
        rr=r.copy(deep=True);rr.solar_stokes_radiance.values*=.25
        rr=add_colorimetry(rr)
        np.testing.assert_allclose(rr.solar_XYZ,r.solar_XYZ*.25,rtol=1e-12)
        self.assertTrue(comparison(r,r)['passes_5percent_Y_and_0_005_xy'])
        rr.attrs['config_json']=json.dumps(replace(ModelConfig(),aod550=.2).to_dict())
        with self.assertRaises(ValueError):comparison(r,rr)

if __name__=='__main__':unittest.main()
