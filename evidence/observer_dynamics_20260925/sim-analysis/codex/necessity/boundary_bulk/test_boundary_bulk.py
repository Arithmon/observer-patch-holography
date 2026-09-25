import importlib.util
import math
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.special import eval_legendre

SPEC = importlib.util.spec_from_file_location("boundary_bulk_certificate", Path(__file__).with_name("certify.py"))
c = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(c)


class BoundaryBulkTests(unittest.TestCase):
    def test_bessel_gamma_identity_with_tail_control(self):
        for ell in (1,2,6):
            for theta in (-0.5,0,0.5):
                row = c.bessel_check(ell, theta, xmax=2048)
                self.assertGreaterEqual(row["exact_minus_finite"], -1e-11)
                self.assertLessEqual(row["exact_minus_finite"], row["analytical_positive_tail_bound"]+1e-11)
        for ell in (1,2,10,30):
            self.assertAlmostEqual(c.bessel_square_exact(ell), 1/(2*ell*(ell+1)), places=14)

    def test_log_kernel_legendre_spectrum(self):
        for ell in (1,2,3,10):
            observed = 2*np.pi*quad(lambda mu: float(c.log_shell_kernel(mu))*eval_legendre(ell, mu),
                                     -1,1,epsabs=1e-11)[0]
            self.assertAlmostEqual(observed,2*np.pi/(ell*(ell+1)),places=10)

    def test_bump_transform_against_independent_direct_radial_quadrature(self):
        for k in (0,0.01,1,3.99,4.01,12,50):
            direct = 4*np.pi*quad(lambda r: r*r*float(c.bump(r))*np.sinc(k*r/np.pi),
                                  3,4,epsabs=1e-13,epsrel=1e-12,limit=200)[0]
            self.assertAlmostEqual(float(c.bump_fourier(k)), direct, delta=2e-12)

    def test_positive_alternatives_and_exact_shell_null(self):
        k = np.geomspace(1e-6,1e4,1001)
        frac = k**3*c.bump_fourier(k)/(240*np.pi)
        self.assertLess(float(np.max(np.abs(frac))), 0.5)
        for ell in range(8):
            self.assertEqual(c.bump_shell_cl(ell,1),0)
        self.assertNotEqual(c.bump_shell_cl(2,2),0)

    def test_forward_fourier_and_chord_counterexample_agree(self):
        for radius in (1,2):
            for ell in (1,2,4):
                numerical, tail = c.bump_shell_fourier_cl(ell,radius,kmax=512)
                expected = c.bump_shell_cl(ell,radius)
                self.assertLess(abs(numerical-expected),tail+1e-10)

    def test_bump_polynomial_bound_is_correct(self):
        # g=(u+3)u^2(1-u)^2, g''=6-30u+12u^2+20u^3.
        poly = np.polynomial.Polynomial([0,0,3,-5,1,1])
        np.testing.assert_allclose(poly.deriv(2).coef,[6,-30,12,20])
        self.assertEqual(poly(0),0)
        self.assertEqual(poly(1),0)
        self.assertEqual(poly.deriv()(0),0)
        self.assertEqual(poly.deriv()(1),0)
        conservative_l1 = sum(abs(a)/(i+1) for i,a in enumerate(poly.deriv(2).coef))
        self.assertEqual(conservative_l1,30)

    def test_round_covariance_dipole_and_polynomial_ambiguity(self):
        # Any c+d*r^2 on a sphere is in the l=0,1 span, hence invisible
        # after the actual paper's l>=2 projection.
        for radius in (0.5,1,3):
            for ell in (2,3,6):
                integral = quad(lambda mu: (0.7-0.4*radius**2*(2-2*mu))*eval_legendre(ell,mu),
                                -1,1,epsabs=1e-12)[0]
                self.assertAlmostEqual(integral,0,places=12)

    def test_inverse_covariance_preserves_anisotropy(self):
        row = c.anisotropy_fixture()
        self.assertLess(row["baseline_rotation_commutator_norm"],1e-14)
        self.assertGreater(row["perturbed_covariance_commutator_norm"],0.001)
        self.assertLess(row["inverse_commutator_identity_residual"],1e-14)
        self.assertGreater(np.ptp(row["dipole_covariance_eigenvalues"]),0.01)

    def test_reject_divergent_bessel_parameters(self):
        for ell,theta in [(0,0),(1,2),(2,-2),(1.5,0)]:
            with self.assertRaises(ValueError):
                c.bessel_square_exact(ell,theta)


if __name__ == "__main__":
    unittest.main()
