"""Independent finite-cloud and analytic normalization controls."""

import cmath
import importlib.util
import math
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("golden_fourier", ROOT / "scripts/codex_golden_fourier.py")
FOURIER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FOURIER)


class FourierTests(unittest.TestCase):
    def test_product_against_explicit_cloud(self):
        for q in (2, 3, 5, 8):
            for mode in ((0, 0, 0), (1, 0, 0), (1, 1, 1), (2, -1, 0), (3, 2, 1)):
                a = FOURIER.amplitude(q, mode)
                b = FOURIER.explicit_cloud_amplitude(q, mode)
                self.assertLess(abs(a - b), 3e-14, (q, mode))

    def test_formula_against_independent_geometric_sum(self):
        phi = (1 + math.sqrt(5)) / 2
        for q in (5, 55, 89, 144):
            for n in (-7, -1, 0, 1, 3, 5):
                values = [cmath.exp(2j * math.pi * n * ((b * phi) % 1)) for b in range(q)]
                direct = complex(math.fsum(z.real for z in values), math.fsum(z.imag for z in values)) / q
                self.assertLess(abs(FOURIER.axis_amplitude(q, n) - direct), 2e-13)

    def test_zero_conjugation_and_axis_permutation(self):
        for q in FOURIER.LEVELS:
            self.assertEqual(FOURIER.amplitude(q, (0, 0, 0)), 1)
            a = FOURIER.amplitude(q, (1, 2, -3))
            b = FOURIER.amplitude(q, (-1, -2, 3))
            self.assertLess(abs(a.conjugate() - b), 1e-18)
            self.assertLess(abs(a - FOURIER.amplitude(q, (-3, 1, 2))), 1e-18)

    def test_equal_norm_modes_are_not_radially_identical(self):
        axis = FOURIER.mode_row(144, (3, 0, 0))
        oblique = FOURIER.mode_row(144, (2, 2, 1))
        self.assertEqual(axis["norm_squared"], oblique["norm_squared"])
        self.assertGreater(axis["normalized_density_power"] / oblique["normalized_density_power"], 1e15)
        # Low-mode suppression is not uniform below the iid benchmark at finite q.
        self.assertGreater(FOURIER.mode_row(144, (5, 0, 0))["structure_factor_N_times_power"], 1)

    def test_fibonacci_identity_and_fixed_mode_asymptotic(self):
        phi = (1 + math.sqrt(5)) / 2
        for q in FOURIER.LEVELS:
            m = FOURIER.fibonacci_index(q)
            p = round(q * phi)
            self.assertAlmostEqual(q * phi - p, (-1)**(m + 1) * phi**(-m), places=12)
            for mode in ((1, 0, 0), (1, 1, 0), (1, 1, 1)):
                r = FOURIER.mode_row(q, mode)
                ratio = r["fixed_mode_fibonacci_asymptotic"]["observed_power_times_q_to_4h_over_constant"]
                delta = (-1)**(m + 1) * phi**(-m)
                # Exact finite-q correction from the Fibonacci identity,
                # rather than an arbitrary closeness-to-limit threshold.
                expected = math.prod(5 * q*q * delta*delta *
                                     (math.sin(math.pi*n*delta)/(math.pi*n*delta))**2
                                     for n in mode if n)
                self.assertAlmostEqual(ratio, expected, places=12)

    def test_input_and_explicit_enumeration_budget(self):
        with self.assertRaises(ValueError):
            FOURIER.axis_amplitude(0, 1)
        with self.assertRaises(ValueError):
            FOURIER.axis_amplitude(4, 1.5)
        with self.assertRaises(ValueError):
            FOURIER.explicit_cloud_amplitude(144, (1, 0, 0))


if __name__ == "__main__":
    unittest.main()
