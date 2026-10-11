"""Maintainer-style Petz controls derived from the original source entries."""
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest

from quantum_information.recovery import petz_recovery, petz_recovery_kraus


def _reference(drift, complex_state):
    conditional = (np.array([[.5, .25j], [-.25j, .5]]) if complex_state
                   else np.diag([.75, .25]))
    reference = np.kron(np.diag([1., 0.]), (1 + drift) * conditional)
    # Source entries and trace are dyadic; the prescribed replacement is
    # formed as exact rational quotients before comparison with the API.
    trace = sum(Fraction(float(x.real)) for x in np.diag(reference))
    expected = np.array([[complex(float(Fraction(float(z.real)) / trace),
                                  float(Fraction(float(z.imag)) / trace))
                          for z in row] for row in reference])
    return reference, expected


@pytest.mark.parametrize("drift", (2. ** -42, -2. ** -42))
@pytest.mark.parametrize("complex_state", (False, True))
def test_kernel_replacement_has_unit_trace_for_accepted_source_drift(drift, complex_state):
    reference, expected = _reference(drift, complex_state)
    kraus = petz_recovery_kraus(reference, (2, 2))
    kernel_input = np.diag([0., 1.])
    actual = sum(k @ kernel_input @ k.conj().T for k in kraus)
    np.testing.assert_allclose(actual, expected, rtol=3e-15, atol=0)
    assert np.trace(actual).real == pytest.approx(1., rel=3e-15, abs=0)
    np.testing.assert_allclose(sum(k.conj().T @ k for k in kraus), np.eye(2),
                               rtol=3e-15, atol=3e-15)


@pytest.mark.parametrize("drift", (2. ** -42, -2. ** -42))
@pytest.mark.parametrize("complex_state", (False, True))
def test_supported_petz_action_keeps_original_source_trace(drift, complex_state):
    reference, expected = _reference(drift, complex_state)
    kraus = petz_recovery_kraus(reference, (2, 2))
    marginal = np.diag([1 + drift, 0.])
    actual = sum(k @ marginal @ k.conj().T for k in kraus)
    np.testing.assert_allclose(actual, reference, rtol=3e-15, atol=0)
    # The normalized supported input and a kernel input are both prepared
    # into the same fixed normalized conditional state for this reference.
    supported = sum(k @ np.diag([1., 0.]) @ k.conj().T for k in kraus)
    np.testing.assert_allclose(supported, expected, rtol=3e-15, atol=0)


@pytest.mark.parametrize("drift", (2. ** -42, -2. ** -42))
@pytest.mark.parametrize("complex_state", (False, True))
def test_public_fixed_point_preserves_accepted_trace_drift(drift, complex_state):
    reference, _ = _reference(drift, complex_state)
    np.testing.assert_array_equal(petz_recovery(reference, (1, 2, 2)), reference)


def _conditional_block_oracle(rho, dims):
    """Exact source conditional law for classical A,B and quantum C blocks."""
    a, b, c = dims
    result = np.zeros_like(rho, dtype=complex)
    for i in range(a):
        for j in range(b):
            mass_ab = sum(Fraction(float(rho[(i*b+j)*c+k, (i*b+j)*c+k].real))
                          for k in range(c))
            mass_b = sum(Fraction(float(rho[(h*b+j)*c+k, (h*b+j)*c+k].real))
                         for h in range(a) for k in range(c))
            for k in range(c):
                for l in range(c):
                    real = sum(Fraction(float(rho[(h*b+j)*c+k, (h*b+j)*c+l].real))
                               for h in range(a))
                    imag = sum(Fraction(float(rho[(h*b+j)*c+k, (h*b+j)*c+l].imag))
                               for h in range(a))
                    result[(i*b+j)*c+k, (i*b+j)*c+l] = complex(
                        float(mass_ab * real / mass_b), float(mass_ab * imag / mass_b))
    return result


def test_classical_recovery_rounds_original_conditional_midpoints():
    tiny = np.nextafter(0., 1.)
    rho = np.diag([.25, .25, .25, .25, 0., 5*tiny, 0., 0.])
    expected = _conditional_block_oracle(rho, (2, 2, 2))
    # q100=5s*.25/(.5+5s) is strictly below 2.5s, while q101 is
    # strictly above it. Rounded Kraus roots lose this distinction first.
    assert expected[4, 4] == 2*tiny
    assert expected[5, 5] == 3*tiny
    np.testing.assert_array_equal(petz_recovery(rho), expected)


@pytest.mark.parametrize("phase", (1., 1j))
def test_coherent_recovery_rounds_original_conditional_midpoints(phase):
    tiny = np.nextafter(0., 1.)
    rho = np.zeros((16, 16), complex)
    rho[:4, :4] = np.diag([.25, .125, 0., 0.])
    rho[2:4, 2:4] = [[.0625, .0625*phase], [.0625*np.conj(phase), .0625]]
    rho[4:8, 4:8] = np.eye(4)*.125
    rho[9, 9] = 5*tiny
    # A fixed Hadamard/phase change on C2,C3 creates coherent BC from
    # original dyadic entries; the rare C1 source cell stays unchanged.
    # The source map still obeys the conditional block law, while BC is
    # not diagonal and the general quantum recovery path is required.
    expected = _conditional_block_oracle(rho, (2, 2, 4))
    assert expected[8, 8] == 2*tiny
    assert expected[10, 11] == phase*tiny
    np.testing.assert_array_equal(petz_recovery(rho, (2, 2, 4)), expected)


def test_dense_noncommuting_source_matches_independent_matrix_square_root():
    # Integer Gram entries plus a strictly positive dyadic diagonal give
    # an exactly positive, trace-one dense source before any producer call.
    rng = np.random.default_rng(1729)
    integer = rng.integers(-2, 3, (8, 8)) + 1j*rng.integers(-2, 3, (8, 8))
    gram = integer @ integer.conj().T
    diagonal = (512 - np.trace(gram).real) / 8
    assert diagonal > 0
    rho = (gram + diagonal*np.eye(8)) / 512
    assert sum(Fraction(float(x.real)) for x in np.diag(rho)) == 1

    ctx = mp.mp.clone()
    ctx.dps = 100
    original = ctx.matrix([[ctx.mpc(z.real, z.imag) for z in row] for row in rho])
    # Direct original-entry index sums in A,B,C order, independent of the
    # producer's reductions, rank decisions and eigensystem root evaluator.
    ab = ctx.matrix([[sum(original[2*i+k, 2*j+k] for k in range(2))
                      for j in range(4)] for i in range(4)])
    bc = ctx.matrix([[sum(original[4*k+i, 4*k+j] for k in range(2))
                      for j in range(4)] for i in range(4)])
    middle = ctx.matrix([[sum(bc[2*i+k, 2*j+k] for k in range(2))
                          for j in range(2)] for i in range(2)])

    def tensor(left, right):
        return ctx.matrix([[left[i//right.rows, j//right.cols]
                            * right[i % right.rows, j % right.cols]
                            for j in range(left.cols*right.cols)]
                           for i in range(left.rows*right.rows)])

    supplied = tensor(ab, ctx.eye(2))
    lifted_reference = tensor(ctx.eye(2), bc)
    assert ctx.norm(lifted_reference*supplied - supplied*lifted_reference) > ctx.mpf(".001")
    assert ctx.norm(bc*tensor(middle, ctx.eye(2))
                    - tensor(middle, ctx.eye(2))*bc) > ctx.mpf(".001")
    # mpmath sqrtm uses a matrix square-root algorithm, independently of
    # the producer's exact-support construction and Hermitian eigensolve.
    root_bc, root_b = ctx.sqrtm(bc), ctx.sqrtm(middle)
    assert ctx.norm(root_bc*root_bc - bc) < ctx.mpf("1e-90")
    assert ctx.norm(root_b*root_b - middle) < ctx.mpf("1e-90")
    conditional = root_bc*tensor(root_b**-1, ctx.eye(2))
    lifted = tensor(ctx.eye(2), conditional)
    answer = lifted*supplied*lifted.H
    expected = np.array([[complex(z) for z in row] for row in answer.tolist()])

    actual = petz_recovery(rho)
    np.testing.assert_allclose(actual.real, expected.real, rtol=3e-13, atol=0)
    off_diagonal = ~np.eye(8, dtype=bool)
    assert np.min(np.abs(expected.imag[off_diagonal])) > 1e-7
    np.testing.assert_allclose(actual.imag[off_diagonal], expected.imag[off_diagonal],
                               rtol=3e-13, atol=0)
    np.testing.assert_array_equal(np.diag(actual).imag, np.zeros(8))


def _mixed_complex_pair():
    pair = np.diag([.125, .375, .25, .25]).astype(complex)
    pair[2, 3], pair[3, 2] = .125j, -.125j
    return pair


@pytest.mark.parametrize("product_end", ("A", "C"))
@pytest.mark.parametrize("perturbed", (False, True))
@pytest.mark.parametrize("drift", (0., 2.**-42))
def test_full_rank_marginal_product_has_exact_recovery(product_end, perturbed, drift):
    end = np.array([[.5, .125j], [-.125j, .5]])
    expected = np.kron(_mixed_complex_pair(), end)
    rho = expected.copy()
    if perturbed:
        # The traceless factor on each end leaves both AB and BC unchanged.
        # Petz therefore returns expected, although the supplied ABC differs.
        z = np.diag([1., -1.])
        rho += np.kron(np.kron(z, z), z)/256
    if product_end == "A":
        rho = rho.reshape((2,)*6).transpose(2, 1, 0, 5, 4, 3).reshape(8, 8)
        expected = expected.reshape((2,)*6).transpose(2, 1, 0, 5, 4, 3).reshape(8, 8)
    rho *= 1+drift
    expected *= 1+drift
    if perturbed:
        assert not np.array_equal(rho, expected)
    np.testing.assert_array_equal(petz_recovery(rho), expected)


def _middle_factor_markov():
    right = np.eye(4, dtype=complex)/4
    right[0, 3], right[3, 0] = .125j, -.125j
    # In A,BL,BR,C order, X_ABL tensor Y_BRC is an exact Markov source.
    # B is supplied as one four-dimensional factor, hiding its inner split.
    return np.kron(_mixed_complex_pair(), right)


def test_exact_markov_certificate_resolves_noncommuting_mixed_end_fixed_point():
    rho = _middle_factor_markov()
    ab = np.array([[sum(rho[2*i+k, 2*j+k] for k in range(2))
                    for j in range(8)] for i in range(8)])
    bc = np.array([[sum(rho[8*k+i, 8*k+j] for k in range(2))
                    for j in range(8)] for i in range(8)])
    first, second = np.kron(ab, np.eye(2)), np.kron(np.eye(2), bc)
    assert np.linalg.norm(first@second-second@first) > .01
    np.testing.assert_array_equal(petz_recovery(rho, (2, 4, 2)), rho)


def test_nearby_non_markov_source_cannot_select_exact_fixed_point_fallback():
    fixed = _middle_factor_markov()
    z = np.diag([1., -1.])
    rho = fixed + np.kron(np.kron(z, np.eye(4)), z)/1024
    # This preserves both marginals, so its exact recovered state is fixed,
    # not rho. A candidate certificate proves the whole output even though
    # the perturbed input itself has nonzero conditional information.
    assert not np.array_equal(rho, fixed)
    np.testing.assert_array_equal(petz_recovery(rho, (2, 4, 2)), fixed)


def test_candidate_certificate_requires_exact_marginals_and_markov_identity():
    from quantum_information.recovery_functions import _certified_markov_candidate
    from quantum_information.state_factors import source_state

    fixed = _middle_factor_markov()
    # Independently reduce the original dyadic fixture in A,B,C order.
    ab = np.array([[sum(fixed[2*i+k, 2*j+k] for k in range(2))
                    for j in range(8)] for i in range(8)])
    bc = np.array([[sum(fixed[8*k+i, 8*k+j] for k in range(2))
                    for j in range(8)] for i in range(8)])
    ab_source, bc_source = source_state(ab), source_state(bc)
    assert _certified_markov_candidate(fixed, ab_source, bc_source, (2, 4, 2))
    z = np.diag([1., -1.])
    hostile = fixed + np.kron(np.kron(z, np.eye(4)), z)/1024
    assert not _certified_markov_candidate(hostile, ab_source, bc_source, (2, 4, 2))
    # A different full-rank product is Markov but has the wrong marginals.
    wrong_marginals = np.eye(16)/16
    assert not _certified_markov_candidate(wrong_marginals, ab_source, bc_source, (2, 4, 2))
