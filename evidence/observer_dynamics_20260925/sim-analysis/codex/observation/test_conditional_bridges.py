import importlib.util
import math
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("bridges", ROOT / "conditional_bridges.py")
BRIDGES = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BRIDGES)


class BridgeTests(unittest.TestCase):
    def test_finite_time_ir_and_stationary_limits_differ(self):
        self.assertEqual(BRIDGES.ou_power(0, 3, noise=2), 12)
        self.assertAlmostEqual(BRIDGES.ou_power(1e-12, 3, noise=2), 12, places=8)
        self.assertAlmostEqual(BRIDGES.ou_power(6, 100, noise=2), 1 / 3)
        self.assertEqual(BRIDGES.ou_power(6, 0), 0)

    def test_variance_evolution_from_independent_time_step(self):
        lam, diffusion, noise, t, dt = 6, 0.3, 2.0, 0.7, 0.4
        old = BRIDGES.ou_power(lam, t, diffusion, noise)
        increment = BRIDGES.ou_power(lam, dt, diffusion, noise)
        expected = math.exp(-2 * diffusion * lam * dt) * old + increment
        self.assertAlmostEqual(BRIDGES.ou_power(lam, t + dt, diffusion, noise), expected)

    def test_conserved_and_nonconserved_stationary_spectra(self):
        for lam in (1, 2, 6, 20, 110):
            self.assertAlmostEqual(BRIDGES.ou_power(lam, 100, noise=2), 2 / lam)
            self.assertAlmostEqual(BRIDGES.ou_power(lam, 100, noise=2, conserved=True), 2)
        self.assertEqual(BRIDGES.ou_power(0, 3, conserved=True), 0)

    def test_history_integral_by_independent_quadrature(self):
        for d in (2, 3):
            for k in (0.25, 1, 2):
                # Substitute t=s^2, then integrate directly by Simpson's rule.
                end = 10 / (math.sqrt(2) * k)
                n, h = 4000, end / 4000
                def integrand(s):
                    return 2 * s**(d-1) * math.exp(-2*k*k*s*s)
                total = integrand(0) + integrand(end)
                total += math.fsum((4 if i % 2 else 2) * integrand(i*h) for i in range(1, n))
                numerical = h * total / 3
                exact = BRIDGES.independent_history_power(k, d)
                self.assertAlmostEqual(numerical / exact, 1, places=9)

    def test_finite_history_has_finite_IR_limit(self):
        self.assertAlmostEqual(BRIDGES.finite_history_zero_mode(0, 4, 3), 16 / 3)
        self.assertAlmostEqual(BRIDGES.finite_history_zero_mode(1, 4, 2), 3)

    def test_distortion_conversion_and_energy_fraction(self):
        data = BRIDGES.result()["conditional_FIRAS_heat_bound"]
        self.assertAlmostEqual(data["mu_per_fractional_energy_release"], 1.4006573255399406)
        self.assertAlmostEqual((55 / 6 - 13 / 2) / (55 / 6), 16 / 55)
        self.assertEqual(data["y_era_fractional_energy_limit_approx"], 6e-5)


if __name__ == "__main__":
    unittest.main()
