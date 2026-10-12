"""Original-entry controls for source support and the resulting Petz action."""

from fractions import Fraction

import numpy as np
import pytest

from quantum_information.recovery import (
    UnresolvedPetzSupport, matrix_inv_sqrt_psd, matrix_sqrt_psd,
    petz_recovery, petz_recovery_kraus,
)


def _apply(kraus, value):
    return sum(k @ value @ k.conj().T for k in kraus)


def test_complex_projector_root():
    # Every supplied entry is exact dyadic. P^2=P and Tr_C(P)=I/2,
    # hence sqrt(P)=P and its Petz channel is X -> Tr(X) P.
    v = np.array([1, 1j, 1, -1j])
    p = np.outer(v, v.conj()) / 4
    np.testing.assert_array_equal(p @ p, p)
    np.testing.assert_allclose(matrix_sqrt_psd(p), p, atol=3e-15, rtol=0)


def test_entangled_reference_channel_has_correct_action():
    v = np.array([1, 1j, 1, -1j])
    p = np.outer(v, v.conj()) / 4
    kraus = petz_recovery_kraus(p, (2, 2))
    for probe in (np.diag([1., 0.]), np.array([[0., 1j], [-1j, 0.]])):
        np.testing.assert_allclose(_apply(kraus, probe), np.trace(probe)*p,
                                   atol=4e-15, rtol=0)


def test_public_recovery_matches_exact_dyadic_polynomial():
    # Orthogonal vectors give BC^2=BC/2, while B=I/2. Therefore the
    # Petz expression is 4 BC (AB tensor I_C) BC, with no numerical
    # square root, inverse, or support decision in the independent oracle.
    psi = np.array([1, 1, 1, 1j]) / 2
    phi = np.array([1, 1j, -1j, -1j]) / 2
    p, q = np.outer(psi, psi.conj()), np.outer(phi, phi.conj())
    assert np.vdot(psi, phi) == 0
    bc = (p+q)/2
    np.testing.assert_array_equal(bc @ bc, bc/2)
    rho = np.zeros((8, 8), complex)
    rho[:4, :4], rho[4:, 4:] = p/2, q/2
    # Original input entries, indexed directly in A,B,C tensor order.
    ab = np.array([[sum(rho[2*i+c, 2*j+c] for c in range(2))
                    for j in range(4)] for i in range(4)])
    lifted = np.kron(np.eye(2), bc)
    expected = 4*lifted @ np.kron(ab, np.eye(2)) @ lifted
    np.testing.assert_allclose(petz_recovery(rho), expected, atol=5e-15, rtol=0)


@pytest.mark.parametrize("consumer", (matrix_inv_sqrt_psd,
    lambda state: petz_recovery_kraus(state, (3, 1))))
def test_false_float_kernel_cannot_become_a_source_kernel(consumer):
    tau = np.array([[.5, .3, 0.], [.3, .18, 0.], [0., 0., .32]])
    # Sylvester's criterion on the ORIGINAL binary64 entries proves rank 3.
    # Main's complex eigh instead returns a zero whose floating action is
    # exactly zero, silently replacing a full-rank reference by rank 2.
    minor = Fraction(.5)*Fraction(.18)-Fraction(.3)**2
    assert minor > 0
    assert sum(Fraction(x) for x in np.diag(tau)) == 1
    # For C=1 a full-rank Petz channel is the identity. The defective
    # completion changes |0><0| by Frobenius norm 0.706341 despite being TP.
    with pytest.raises(UnresolvedPetzSupport):
        consumer(tau)


def test_true_singular_kernel_completion_is_still_available():
    # This neighboring decimal block really IS rank one over exact entries:
    # .6=2*.3=4*.15. Its nonzero eigenvalue is the original trace s.
    tau = np.array([[.6, .3, 0.], [.3, .15, 0.], [0., 0., .25]])
    assert Fraction(.6)*Fraction(.15) == Fraction(.3)**2
    support = np.array([[.8, .4, 0.], [.4, .2, 0.], [0., 0., 1.]])
    kernel_state = np.eye(3)-support
    # LAPACK can report this exact zero as a tiny positive value (for
    # example 3.8858e-18 on Accelerate). The API deliberately refuses an
    # unresolved floating support; exact source rank must not override
    # that admission guard. Check the refusal instead of demanding that
    # every platform resolve this particular rotated kernel.
    values, vectors = np.linalg.eigh(tau.astype(complex))
    if (np.count_nonzero(values > 0) != 2
            or np.any(tau @ vectors[:, values <= 0] != 0)):
        with pytest.raises(UnresolvedPetzSupport):
            petz_recovery_kraus(tau, (3, 1))
        return
    kraus = petz_recovery_kraus(tau, (3, 1))
    np.testing.assert_allclose(_apply(kraus, kernel_state), tau, atol=3e-15, rtol=0)
    probe = np.diag([0., 0., 1.])
    np.testing.assert_allclose(_apply(kraus, probe), probe, atol=3e-15, rtol=0)


def test_resolved_rotated_kernel_completion_with_source_trace_drift():
    # Equal entries give a directly cancelling numerical kernel while the
    # final diagonal retains exact source trace drift. Completion replaces
    # kernel mass by tau/Tr(tau), and supported action still recovers tau.
    tau = np.array([[.375, .375, 0.], [.375, .375, 0.],
                    [0., 0., np.nextafter(.25, 0.)]])
    trace = sum(Fraction(float(x)) for x in np.diag(tau))
    assert trace != 1
    expected = np.array([[float(Fraction(float(x))/trace) for x in row]
                         for row in tau])
    kernel = np.array([[.5, -.5, 0.], [-.5, .5, 0.], [0., 0., 0.]])
    kraus = petz_recovery_kraus(tau, (3, 1))
    np.testing.assert_allclose(_apply(kraus, kernel), expected, atol=3e-15, rtol=0)
    np.testing.assert_allclose(_apply(kraus, tau), tau, atol=3e-15, rtol=0)


def test_complex_full_rank_reference_with_trivial_c_is_identity():
    tau = np.array([[.5, .125j], [-.125j, .5]])
    kraus = petz_recovery_kraus(tau, (2, 1))
    for probe in (np.diag([1., 0.]), np.array([[0., 1j], [-1j, 0.]])):
        np.testing.assert_allclose(_apply(kraus, probe), probe, atol=3e-15, rtol=0)


@pytest.mark.parametrize("pure_end", ["A", "C"])
def test_exact_complex_product_fixed_point_is_not_made_indefinite(pure_end):
    vector = np.array([1, 1, 1, 1, 1, 1, 1j, 1j])
    rho = np.outer(vector, vector.conj()) / 8
    if pure_end == "A":
        rho = rho.reshape(2, 2, 2, 2, 2, 2).transpose(2, 1, 0, 5, 4, 3).reshape(8, 8)
    # One end is pure and independent. The exact Petz formula therefore
    # recovers this input, regardless of numerical root-factor roundoff.
    np.testing.assert_array_equal(petz_recovery(rho), rho)
