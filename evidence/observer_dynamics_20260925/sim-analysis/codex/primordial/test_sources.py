"""Independent numerical and limiting checks of the primordial source."""
import importlib.util
import math
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.special import gammaincc

MODULE_SPEC = importlib.util.spec_from_file_location("primordial_sources", Path(__file__).with_name("sources.py"))
s = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(s)


class SourceTests(unittest.TestCase):
    def test_finite_age_integral_matches_independent_log_time_quadrature(self):
        # Set 2*kappa=1, so t_min/max equal inverse squared cutoffs.
        beta, kir, kuv, pivot = 1.51697505, 0.003, 0.4, 0.05
        def integral(k):
            return quad(lambda u: math.exp(beta*u-k*k*math.exp(u)),
                        math.log(kuv**-2), math.log(kir**-2), epsabs=1e-10, epsrel=1e-11)[0]
        ks = np.array([1e-5, 0.002, 0.02, 0.05, 0.2, 1.2])
        expected = np.array([(k/pivot)**3*integral(k)/integral(pivot) for k in ks])
        actual = s.history_delta2(ks, As=1, beta=beta, k_ir=kir, k_uv=kuv, k_pivot=pivot)
        np.testing.assert_allclose(actual, expected, rtol=3e-10)

    def test_pivot_and_units(self):
        for row in s.load_spec()["candidates"]:
            f = s.source_callable(row["id"])
            self.assertAlmostEqual(float(f(0.05)), 2.100549e-9, delta=1e-24)
        delta = np.array([2e-9, 3e-9])
        k = np.array([0.01, 0.1])
        power = s.dimensional_curvature_power(k, delta)
        np.testing.assert_allclose(power*k**3/(2*np.pi**2), delta, rtol=1e-15)

    def test_intermediate_slope_and_calibrated_powerlaw(self):
        k = np.geomspace(1e-4, 0.1, 50)
        flat = s.source_callable("history_scale_invariant_wide")(k)
        np.testing.assert_allclose(flat, flat[0], rtol=1e-6)
        history = s.source_callable("history_tilt_calibrated_wide")(k)
        standard = s.source_callable("planck_powerlaw")(k)
        np.testing.assert_allclose(history, standard, rtol=1e-6)
        slope = math.log(history[-1]/history[0])/math.log(k[-1]/k[0])
        self.assertAlmostEqual(slope, 0.9660499-1, delta=2e-7)

    def test_finite_age_extreme_ir_is_white_dimensional_power(self):
        k = np.array([1e-8, 2e-8])
        d = s.history_delta2(k, As=1, beta=1.51697505, k_ir=0.003, k_uv=0.4)
        self.assertAlmostEqual(d[1]/d[0], 8, delta=1e-8)
        power = s.dimensional_curvature_power(k, d)
        self.assertAlmostEqual(power[1]/power[0], 1, delta=1e-9)

    def test_complementary_gamma_prevents_uv_cancellation(self):
        # Both lower regularized Gamma evaluations round to one here.
        actual = float(s.gamma_window(10, 1.5, 0.01, 1))
        expected = float(gammaincc(1.5, 100))
        self.assertGreater(actual, 0)
        self.assertAlmostEqual(actual/expected, 1, places=14)

    def test_floor_is_explicit_and_grid_records_it(self):
        exact = s.source_callable("history_tilt_calibrated_uv")(100)
        floored = s.source_callable("history_tilt_calibrated_uv", numerical_floor=1e-300)(100)
        self.assertEqual(float(exact), 0)
        self.assertEqual(float(floored), 1e-300)
        receipt, grid = s.build_receipts()
        self.assertGreater(receipt["candidates"]["history_tilt_calibrated_uv"]["floored_grid_points"], 0)
        self.assertEqual(len(grid["k_Mpc_inverse"]), receipt["grid_points"])
        self.assertTrue(all(min(v) >= grid["numerical_floor"] for v in grid["delta_R_squared"].values()))

    def test_reject_invalid_parameters(self):
        for k in [0, -1, float("nan"), float("inf")]:
            with self.assertRaises(ValueError):
                s.gamma_window(k, 1.5, 1e-6, 10)
        for beta, kir, kuv in [(0, 1, 2), (1.5, 2, 1), (1.5, 1, 1), (1.5, -1, 2)]:
            with self.assertRaises(ValueError):
                s.gamma_window(0.05, beta, kir, kuv)
        with self.assertRaises(ValueError):
            s.history_delta2(1, As=1, beta=1.5, k_ir=1e-5, k_uv=1e-4, k_pivot=1)


if __name__ == "__main__":
    unittest.main()
