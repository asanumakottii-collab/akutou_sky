"""Fast tests that need neither SASKTRAN2 runs nor input data."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import colour

from akutou_sky.color import BANDS, band_mean, chromaticity, cmf_weights, xyz_from_radiance
from akutou_sky.led import Primaries, SRGB_PRIMARIES, load_primaries, to_duty, write_header
from akutou_sky.numerics import (MAX_SOURCE_WORK, Numerics, PRESETS, atmosphere_grid, check_workload,
                                 depression_groups, horizontal_grid, layered_grid, single_scatter_grid,
                                 source_grid, workload)
from akutou_sky.optics import column_profile, integrate_linear, water_uptake
from akutou_sky.scene import SCENES, Scene
from akutou_sky.solver import spectral_grid
from akutou_sky.delta_m import truncated_coefficients


class SpectrumAndColour(unittest.TestCase):
    def test_piecewise_linear_integral_is_exact(self):
        x = np.array([0., 1., 3.])
        y = np.array([0., 2., 2.])
        self.assertTrue(np.allclose(integrate_linear(x, y, [0., .5, 2., 3.]), [.25, 2.75, 2.]))

    def test_spectral_bins_cover_visible_range(self):
        w, edges = spectral_grid(10.)
        self.assertEqual(len(w), 40)
        self.assertEqual((edges[0], edges[-1]), (380., 780.))

    def test_flat_spectrum_luminance_matches_cie(self):
        _, edges = spectral_grid(5.)
        bounds = np.c_[edges[:-1], edges[1:]]
        y = xyz_from_radiance(np.ones(len(bounds)), bounds)[1]
        cmf = colour.MSDS_CMFS['CIE 1931 2 Degree Standard Observer']
        mask = (cmf.wavelengths >= 380) & (cmf.wavelengths <= 780)
        expected = 683 * np.trapezoid(cmf.values[mask, 1], cmf.wavelengths[mask])
        self.assertAlmostEqual(y, expected, delta=1e-6 * expected)

    def test_bin_width_does_not_change_flat_spectrum_colour(self):
        values = []
        for step in (5., 20.):
            _, edges = spectral_grid(step)
            values.append(xyz_from_radiance(np.ones(len(edges) - 1), np.c_[edges[:-1], edges[1:]]))
        self.assertTrue(np.allclose(values[0], values[1], rtol=1e-9))

    def test_black_has_undefined_chromaticity(self):
        self.assertTrue(np.isnan(chromaticity(np.zeros(3))).all())

    def test_band_mean_of_constant_field_is_constant(self):
        e = np.array([0., 2., 10., 40., 70.])
        xyz = np.ones((len(e), 3)) * [1., 2., 3.]
        for lo, hi in BANDS:
            self.assertTrue(np.allclose(band_mean(e, xyz, lo, hi), [1., 2., 3.]))

    def test_band_mean_weights_light_not_rgb(self):
        e = np.array([0., 10.])
        xyz = np.array([[0., 0., 0.], [2., 2., 2.]])
        mean = band_mean(e, xyz, 0., 10.)
        self.assertTrue(0.9 < mean[1] < 1.0)  # cos h weights lower elevations slightly more


class SceneAndNumerics(unittest.TestCase):
    def test_scene_validation_and_json_round_trip(self):
        with self.assertRaises(ValueError):
            Scene(relative_humidity=.99)
        with self.assertRaises(ValueError):
            Scene.from_dict({'nonsense': 1})
        scene = SCENES['humid']
        self.assertEqual(Scene.from_dict(json.loads(scene.to_json())), scene)

    def test_numerics_round_trip_and_refinement(self):
        n = Numerics()
        self.assertEqual(Numerics.from_dict(json.loads(n.to_json())), n)
        self.assertEqual(n.refined('incoming_nodes', 1.4).incoming_nodes, 434)
        self.assertEqual(n.refined('atmosphere_steps_m', 2).atmosphere_steps_m[0], n.atmosphere_steps_m[0] / 2)
        self.assertEqual(n.moments(SCENES['standard']), 64)
        self.assertEqual(n.moments(SCENES['humid']), 256)

    def test_grids_contain_layer_edges_and_source_is_inside(self):
        n = Numerics()
        for grid in (atmosphere_grid(n), single_scatter_grid(n), layered_grid(n.source_steps_m)):
            for edge in (0., 1000., 10000., 30000., 100000.):
                self.assertIn(edge, grid)
        s = source_grid(n)
        self.assertTrue((s > 0).all() and (s < 100000).all() and np.all(np.diff(s) > 0))

    def test_horizontal_grid_covers_every_observer(self):
        n = Numerics()
        g = horizontal_grid(n, [-6, 12])
        self.assertTrue(np.all(np.diff(g) > 0) and g[0] == -45 and g[-1] == 45)
        fine = np.diff(g) <= n.horizontal_step_deg + 1e-9
        self.assertTrue(g[fine.argmax()] <= -12 - n.horizontal_margin_deg + n.horizontal_step_deg)

    def test_depression_groups_partition_and_fit(self):
        n = PRESETS['preview']
        d = [x / 2 for x in range(24, -13, -1)]
        groups = depression_groups(n, d)
        self.assertEqual(sorted(i for g in groups for i in g), list(range(len(d))))
        self.assertTrue(all(workload(n, [d[i] for i in g]) <= MAX_SOURCE_WORK for g in groups))

    def test_fog_requires_fine_boundary_layer_and_scalar(self):
        with self.assertRaises(ValueError):
            check_workload(Numerics(num_stokes=3), SCENES['humid'], [0])
        with self.assertRaises(ValueError):
            check_workload(PRESETS['preview'].replace(atmosphere_steps_m=(200., 500., 1000., 2000.)),
                           SCENES['humid'], [0])


class Microphysics(unittest.TestCase):
    def test_no_growth_without_humidity(self):
        g, r, n = water_uptake(Scene())
        self.assertEqual((g, r), (1., 100.))
        self.assertEqual(n, complex(1.45, -0.005))

    def test_kappa_growth(self):
        g, r, n = water_uptake(Scene(relative_humidity=.9, hygroscopicity_kappa=.3))
        self.assertAlmostEqual(g**3, 1 + .3 * .9 / .1)
        self.assertTrue(1.333 < n.real < 1.45)

    def test_column_profile_recovers_column(self):
        h = np.linspace(0, 100000, 2001)
        p = column_profile(h, np.exp(-h / 1500), .1)
        self.assertAlmostEqual(np.trapezoid(p, h), .1)

    def test_delta_m_truncation_normalized(self):
        coefficients = np.array([1.25, 1.5, 1.6, 1.4, 1.2])[:, None]
        t = truncated_coefficients(coefficients, 4)
        self.assertEqual(t[0, 0], 1.)
        self.assertTrue(np.allclose(t[1:, 0], coefficients[1:4, 0] - (2 * np.arange(1, 4) + 1) * .25))


class Led(unittest.TestCase):
    def setUp(self):
        self.p = Primaries.from_xyY(SRGB_PRIMARIES)

    def test_white_is_all_channels(self):
        white = self.p.xyz(np.ones(3))
        self.assertAlmostEqual(white[1], 1.)
        self.assertTrue(np.allclose(chromaticity(white), [.3127, .329], atol=2e-3))

    def test_in_gamut_round_trip(self):
        target = self.p.xyz([.2, .5, .1])
        self.assertTrue(np.allclose(to_duty(target, self.p), [.2, .5, .1]))

    def test_out_of_gamut_keeps_luminance(self):
        x, y, Y = .70, .29, .1  # redder than the sRGB red primary, dim enough to fit
        target = np.array([x / y * Y, Y, (1 - x - y) / y * Y])
        duty = to_duty(target, self.p)
        self.assertTrue((duty >= 0).all())
        self.assertAlmostEqual(self.p.xyz(duty)[1], Y)

    def test_primaries_json_per_band(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'p.json'
            path.write_text(json.dumps({'bands': [SRGB_PRIMARIES] * 4}))
            self.assertEqual(len(load_primaries(path)), 4)

    def test_header_shape(self):
        table = dict(depression=np.array([0., -1.]), azimuth=0., exponent=.5, reference_luminance=100.,
                     duty=np.random.default_rng(0).random((2, 4, 3)))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'a.h'
            write_header(table, path)
            text = path.read_text()
        self.assertIn('AKUTOU_DUTY[AKUTOU_N_DEPRESSION][AKUTOU_N_BAND][3]', text)
        self.assertIn('{-1.00f, 0.00f}', text)  # valid C floats, depression axis sorted ascending


if __name__ == '__main__':
    unittest.main()
