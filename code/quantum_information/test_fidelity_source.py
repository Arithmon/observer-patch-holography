"""Source-entry controls for fidelity and positive square roots.

The scalar oracles below precede the producer: pure-state trace identities,
explicit dyadic eigenbases, and a two-dimensional block formula. No oracle
uses the producer's roots, eigenvalues, support choices, or normalization.
"""
from fractions import Fraction

import mpmath
import numpy as np
import pytest

from quantum_information.recovery import (
    matrix_inv_sqrt_psd, matrix_sqrt_psd, petz_recovery_kraus, state_fidelity,
)


def _projectors():
    vectors = [np.array(v, complex) for v in
               ([1, 1j, -1, -1j], [1, 1, 1, 1],
                [1, -1, 1, -1], [1, -1j, -1, 1j])]
    return [np.outer(v, v.conj()) / 4 for v in vectors]


def test_complex_pure_root_preserves_the_original_projector():
    p = _projectors()[0]
    np.testing.assert_allclose(matrix_sqrt_psd(p), p, rtol=3e-15, atol=0)


def _walsh_state(weights, denominator):
    h = np.array([[1, 1, 1, 1], [1, -1, 1, -1],
                  [1, 1, -1, -1], [1, -1, -1, 1]], float) / 2
    state = np.zeros((5, 5))
    state[:4, :4] = (h * np.asarray(weights) / denominator) @ h.T
    state[4, 4] = 1 - sum(weights) / denominator
    return state


@pytest.mark.parametrize("phase", (1., 1j))
def test_irrational_root_has_exact_zero_beyond_source_power_zeros(phase):
    rho = _walsh_state([1, 4, 4, 9], 32).astype(complex)
    rho[0, :] *= phase
    rho[:, 0] *= np.conj(phase)
    # The root eigenvalues are (1,2,2,3)/sqrt(32). Their signed
    # sum for entry(0,3) vanishes although rho(0,3)=1/64 is nonzero.
    root = matrix_sqrt_psd(rho)
    assert root[0, 3] == 0
    assert root[0, 1] == pytest.approx(-phase / (2 * np.sqrt(32)), rel=3e-15, abs=0)
    assert root[0, 0].real == pytest.approx(2 / np.sqrt(32), rel=3e-15, abs=0)


@pytest.mark.parametrize("phase", (1., 1j))
def test_irrational_inverse_root_has_exact_zero_beyond_source_power_zeros(phase):
    rho = _walsh_state([4, 36, 9, 9], 128).astype(complex)
    rho[0, :] *= phase
    rho[:, 0] *= np.conj(phase)
    # Inverse weights sqrt(128)*(1/2,1/6,1/3,1/3) cancel at(0,2).
    inverse = matrix_inv_sqrt_psd(rho)
    assert inverse[0, 2] == 0
    assert inverse[0, 0].real == pytest.approx(np.sqrt(128) / 3, rel=3e-15, abs=0)
    assert inverse[0, 1] == pytest.approx(phase * np.sqrt(128) / 12, rel=3e-15, abs=0)


def test_nearby_nonzero_root_component_is_not_a_structural_zero():
    epsilon = 2. ** -44
    rho = _walsh_state([1, 4, 4, 9], 32)
    rho[:4, :4] += epsilon / 4
    rho[4, 4] -= epsilon
    ctx = mpmath.mp.clone()
    ctx.dps = 90
    expected = float((ctx.sqrt(ctx.mpf(1) / 32 + ctx.mpf(epsilon))
                      - ctx.sqrt(ctx.mpf(1) / 32)) / 4)
    assert matrix_sqrt_psd(rho)[0, 3].real == pytest.approx(expected, rel=3e-15, abs=0)


@pytest.mark.parametrize("rank", (1, 2))
def test_complex_exact_rank_fidelity_against_maximally_mixed(rank):
    p, q, _, _ = _projectors()
    rho = p if rank == 1 else (p + q) / 2
    # A normalized uniform rank-r projector has fidelity r/d against I/d.
    assert state_fidelity(rho, np.eye(4) / 4) == pytest.approx(rank / 4, rel=3e-15, abs=0)


@pytest.mark.parametrize("rank", (1, 2))
def test_complex_orthogonal_support_has_exact_zero_fidelity(rank):
    p, q, r, s = _projectors()
    a, b = (p, q) if rank == 1 else ((p + q) / 2, (r + s) / 2)
    assert state_fidelity(a, b) == 0
    assert state_fidelity(b, a) == 0


@pytest.mark.parametrize("exponent", (20, 40, 50))
def test_full_rank_near_pure_state_retains_each_positive_eigenvalue(exponent):
    epsilon = 2. ** -exponent
    p = _projectors()[0]
    rho = (1 - 3 * epsilon) * p + epsilon * (np.eye(4) - p)
    ctx = mpmath.mp.clone()
    ctx.dps = 90
    e = ctx.mpf(epsilon)
    expected = float((ctx.sqrt(1 - 3 * e) + 3 * ctx.sqrt(e)) ** 2 / 4)
    assert state_fidelity(rho, np.eye(4) / 4) == pytest.approx(expected, rel=3e-15, abs=0)


@pytest.mark.parametrize("exponent", (40, 500, 1000))
def test_tiny_complex_positive_block_is_not_an_eigenvalue_cutoff(exponent):
    mass = 2. ** -exponent
    rho = np.diag([1 - 2 * mass, mass, mass]).astype(complex)
    rho[1, 2], rho[2, 1] = .5j * mass, -.5j * mass
    sigma = np.diag([0., .5, .5])
    ctx = mpmath.mp.clone()
    ctx.dps = 90
    expected = float(ctx.mpf(mass) * (1 + ctx.sqrt(3) / 2))
    assert state_fidelity(rho, sigma) == pytest.approx(expected, rel=3e-15, abs=0)


@pytest.mark.parametrize("mass", (np.nextafter(0., 1.), 2. ** -1000))
def test_rank_two_tiny_fidelity_rounds_on_its_own_scale(mass):
    rho = np.diag([1., mass, mass, 0.])
    sigma = np.diag([0., .5, .5, 0.])
    assert state_fidelity(rho, sigma) == 2 * mass


@pytest.mark.parametrize("multiple,rounded", ((1, 4), (3, 14)))
def test_rank_two_exact_midpoints_round_ties_in_both_directions(multiple, rounded):
    unit = np.nextafter(0., 1.)
    rho = np.diag([1., multiple * unit, 4 * multiple * unit, 0.])
    sigma = np.diag([0., .5, .5, 0.])
    # F=(sqrt(m/2)+sqrt(4m/2))**2=9m/2 exactly. The two
    # controls round 4.5 units down to 4 and 13.5 units up to 14.
    assert state_fidelity(rho, sigma) == rounded * unit


@pytest.mark.parametrize("mass", (np.nextafter(0., 1.), 2. ** -1000))
def test_rank_two_nonsquare_moment_preserves_positive_mass(mass):
    rho = np.diag([1., mass, 2 * mass, 0.])
    sigma = np.diag([0., .5, .5, 0.])
    ctx = mpmath.mp.clone()
    ctx.dps = 90
    expected = float(ctx.mpf(mass) * (3 + 2 * ctx.sqrt(2)) / 2)
    assert state_fidelity(rho, sigma) == expected


def test_nonzero_fidelity_below_float_range_is_refused():
    mass = 2. ** -550
    rho, sigma = np.diag([1., mass, 0.]), np.diag([0., mass, 1.])
    # The sole common diagonal cell gives exact fidelity mass**2 > 0.
    with pytest.raises(ValueError, match="underflow|precision"):
        state_fidelity(rho, sigma)


@pytest.mark.parametrize("exponent", (20, 100, 250, 500))
def test_nearly_orthogonal_complex_pure_states_use_the_overlap_scale(exponent):
    epsilon = 2. ** -exponent
    vector = np.array([1, 1j, -1, -1j], complex)
    a = np.zeros((8, 8), complex)
    a[:4, :4] = np.outer(vector, vector.conj()) / 4
    tilted = np.concatenate([epsilon * vector, vector])
    b = np.outer(tilted, tilted.conj()) / 4
    # Tr(b)=1+epsilon**2 is accepted trace roundoff, not renormalized.
    expected = float(Fraction(epsilon) ** 2)
    assert state_fidelity(a, b) == pytest.approx(expected, rel=3e-15, abs=0)


@pytest.mark.parametrize("drift", (2. ** -42, -2. ** -42))
def test_accepted_trace_roundoff_remains_homogeneous(drift):
    rho = (1 + drift) * _projectors()[0]
    expected = float(Fraction(1 + drift) ** 2)
    assert state_fidelity(rho, rho) == pytest.approx(expected, rel=3e-15, abs=0)


def test_exact_hermitian_average_precedes_subnormal_rounding():
    mass = np.nextafter(0., 1.)
    rho = np.diag([1., mass, mass]).astype(complex)
    rho[1, 2] = 1j * mass
    # Exact symmetrization gives eigenvalues mass/2 and 3*mass/2.
    ctx = mpmath.mp.clone()
    ctx.dps = 90
    low, high = ctx.sqrt(ctx.mpf(mass) / 2), ctx.sqrt(3 * ctx.mpf(mass) / 2)
    expected = float((high - low) / 2)
    root = matrix_sqrt_psd(rho)
    assert root[1, 2].imag == pytest.approx(expected, rel=3e-15, abs=0)
    sigma = np.zeros((3, 3), complex)
    sigma[1:, 1:] = [[.5, .5j], [-.5j, .5]]
    # B is pure, so F=Tr(AB)=3*mass/2, a ties-to-even midpoint.
    assert state_fidelity(rho, sigma) == 2 * mass


@pytest.mark.parametrize("bad", (np.diag([1., -2. ** -44]),
                                  np.array([[1., 2. ** -100], [2. ** -100, 0.]])))
def test_exact_indefinite_roundoff_is_not_a_positive_state(bad):
    with pytest.raises(ValueError, match="positive semidefinite"):
        state_fidelity(bad, np.diag([1., 0.]))
    with pytest.raises(ValueError, match="positive semidefinite"):
        matrix_sqrt_psd(bad)


@pytest.mark.parametrize("bad", ([[True, 0.], [0., 0.]],
                                  np.ma.array(np.eye(2) / 2, mask=[[False, True], [False, False]]),
                                  [[1., 0.], [0.]],
                                  [[float("nan"), 0.], [0., 0.]]))
def test_source_validation_rejects_missing_or_malformed_entries(bad):
    with pytest.raises(ValueError):
        state_fidelity(bad, np.eye(2) / 2)
    with pytest.raises(ValueError):
        matrix_sqrt_psd(bad)


def test_source_conversion_rejects_lost_extended_precision():
    if np.finfo(np.longdouble).nmant <= np.finfo(float).nmant:
        pytest.skip("platform longdouble has no precision beyond binary64")
    rho = np.zeros((2, 2), dtype=np.longdouble)
    rho[0, 0] = np.longdouble(1) + np.longdouble(2) ** -60
    with pytest.raises(ValueError, match="loses information"):
        state_fidelity(rho, np.diag([1., 0.]))
    with pytest.raises(ValueError, match="loses information"):
        matrix_sqrt_psd(rho)


def test_fidelity_rounding_preserves_numpy_error_context():
    mass = np.nextafter(0., 1.)
    rho, sigma = np.diag([1., mass, mass, 0.]), np.diag([0., .5, .5, 0.])
    with np.errstate(all="raise"):
        assert state_fidelity(rho, sigma) == 2 * mass


def test_petz_refuses_unrepresentable_subnormal_reference_coherence():
    mass = np.nextafter(0., 1.)
    reference = np.array([[.5, 1j * mass], [-1j * mass, .5]])
    # With B dimension one this channel must replace input by reference.
    # Rounded binary64 root columns instead double its tiny coherence when
    # multiplied, despite passing a matrix-scale completeness tolerance.
    with pytest.raises(ValueError, match="precision|represent|unresolved"):
        petz_recovery_kraus(reference, (1, 2))


def test_petz_retains_representable_subnormal_replacement_population():
    mass = np.nextafter(0., 1.)
    reference = np.diag([1., mass])
    kraus = petz_recovery_kraus(reference, (1, 2))
    # sqrt(2^-1074)=2^-537 is exactly representable; no blanket
    # subnormal or small-support refusal is justified for this neighbor.
    actual = sum(k @ np.ones((1, 1)) @ k.conj().T for k in kraus)
    assert actual[1, 1] == mass
    assert actual[0, 1] == 0
