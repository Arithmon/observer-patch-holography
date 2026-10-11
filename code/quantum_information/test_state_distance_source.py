"""Independent source-entry controls for roots and trace distance."""
from fractions import Fraction

import mpmath
import numpy as np
import pytest

from quantum_information.recovery import matrix_sqrt_psd, trace_distance


def test_half_subnormal_hermitian_entries_cannot_be_erased_before_distance():
    tiny = np.nextafter(0., 1.)
    base = np.diag([1., 4 * tiny, 4 * tiny, 4 * tiny, 4 * tiny]).astype(complex)
    perturbed = base.copy()
    perturbed[1, 2] = 1j * tiny
    perturbed[3, 4] = tiny
    # Exact Hermitian averaging gives two blocks with eigenvalues +/-tiny/2.
    # Their half trace norm is tiny, although every averaged entry rounds to 0.
    assert trace_distance(base, perturbed) == tiny
    assert trace_distance(perturbed, base) == tiny
    assert trace_distance(base, base) == 0


@pytest.mark.parametrize("factor, expected", [(1, 0), (3, 2), (5, 2), (7, 4)])
def test_trace_distance_exact_subnormal_midpoints(factor, expected):
    tiny = np.nextafter(0., 1.)
    a = np.diag([1., 8 * tiny, 8 * tiny]).astype(complex)
    b = a.copy()
    b[1, 2] = factor * tiny
    # One block has trace distance factor*tiny/2, an exact rational midpoint.
    if not expected:
        with pytest.raises(ValueError, match="underflow|precision"):
            trace_distance(a, b)
    else:
        assert trace_distance(a, b) == expected * tiny


@pytest.mark.parametrize("exponent", [20, 200, 500, 1000])
def test_complex_distance_is_measured_on_the_difference_scale(exponent):
    mass = 2. ** -exponent
    a = np.diag([1., 4 * mass, 4 * mass]).astype(complex)
    # Both sources have exactly the same trace and are accepted without
    # normalization once 8*mass is below trace tolerance.
    a[0, 0] = 1 - 8 * mass
    b = a.copy()
    b[1, 2], b[2, 1] = complex(mass, mass), complex(mass, -mass)
    ctx = mpmath.mp.clone()
    ctx.dps = 80
    expected = float(ctx.sqrt(2) * ctx.mpf(mass))
    assert trace_distance(a, b) == pytest.approx(expected, rel=3e-15, abs=0)


def test_root_cannot_infer_zero_from_agreeing_low_precision_evaluations():
    # The perturbation is lost at both initial MP precisions if added to the
    # diagonal before solving. sqrt(A)[0,1] = e/(sqrt(.5+e)+sqrt(.5-e)).
    epsilon = 2. ** -1000
    a = np.array([[.5, epsilon], [epsilon, .5]])
    expected = float(Fraction(epsilon)) / np.sqrt(2)
    assert matrix_sqrt_psd(a)[0, 1].real == pytest.approx(expected, rel=3e-15, abs=0)


def test_exact_root_identity_certifies_a_nonstructural_zero():
    # B is PSD and B[0,2]=0, but A=B^2 and A^2 both have nonzero [0,2].
    # Powers alone cannot certify this particular root zero.
    b = np.array([[.5, .125, 0.], [.125, .25, .125], [0., .125, .5]])
    # Add a separate diagonal block so Tr(B^2)=1 without changing this root.
    a = np.zeros((4, 4))
    a[:3, :3] = b @ b
    a[3, 3] = 1 - np.trace(a)
    root = matrix_sqrt_psd(a)
    np.testing.assert_allclose(root[:3, :3], b, atol=0, rtol=3e-15)
    assert root[0, 2] == 0


def test_distance_does_not_hide_an_indefinite_source():
    with pytest.raises(ValueError, match="positive semidefinite"):
        trace_distance(np.array([[1., 2. ** -100], [2. ** -100, 0.]]), np.diag([1., 0.]))


def test_semidefinite_difference_uses_exact_trace_at_an_irrational_spectrum_tie():
    tiny = np.nextafter(0., 1.)
    a = np.diag([1., 4*tiny, 3*tiny])
    a[1, 2] = a[2, 1] = tiny
    b = np.diag([1., tiny, tiny])
    # D's positive eigenvalues are (5 +/- sqrt(5))*tiny/2, so
    # ||D||_1/2=5*tiny/2, exactly halfway between floats 2*tiny and 3*tiny.
    assert trace_distance(a, b) == 2*tiny
    assert trace_distance(b, a) == 2*tiny
