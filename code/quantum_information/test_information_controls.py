"""Scalar, matrix-log and exact-product controls independent of the evaluator."""

import itertools
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import logm

from quantum_information import conditional_mutual_information as cmi, mutual_information as mi
from quantum_information import is_markov_exact
from quantum_information.information import weighted_information


def _reduce(a, dims, keep):
    # Independent NumPy contraction, used only for well-conditioned controls.
    dims = list(dims)
    tensor = a.reshape(dims+dims)
    for k in reversed(range(len(dims))):
        if k not in keep:
            tensor = tensor.trace(axis1=k, axis2=k+len(dims))
            dims.pop(k)
    return tensor.reshape(int(np.prod(dims)), int(np.prod(dims)))


def _entropy_logm(a):
    return float(-np.trace(a@logm(a)).real)


@pytest.mark.parametrize("seed", range(5))
def test_quantum_cmi_matches_separate_matrix_log_definition(seed):
    rng = np.random.default_rng(seed+41)
    x = rng.normal(size=(12, 12))+1j*rng.normal(size=(12, 12))
    rho = x@x.conj().T+np.eye(12)
    rho /= np.trace(rho).real
    dims = [2, 3, 2]
    expected = (_entropy_logm(_reduce(rho, dims, [0, 1]))
                + _entropy_logm(_reduce(rho, dims, [1, 2]))
                - _entropy_logm(_reduce(rho, dims, [1]))-_entropy_logm(rho))
    assert cmi(rho, dims, [0], [1], [2]) == pytest.approx(expected, abs=2e-14)


@pytest.mark.parametrize("delta", (2.**-30, 1e-50, 1e-150))
def test_noncommuting_perturbation_of_a_nonuniform_markov_state(delta):
    # AB is correlated; C is independent. The perturbation leaves both AB
    # and BC unchanged, so CMI is exactly the joint entropy loss. Each
    # complement pair gives a separate analytic 2x2 characteristic root.
    diagonal = np.repeat([.375, .125, .0625, .4375], 2)/2
    rho = np.diag(diagonal).astype(complex)
    for i in range(4):
        rho[i, 7-i] = 1j*delta
        rho[7-i, i] = -1j*delta
    ctx = mp.mp.clone()
    ctx.dps = 650
    pieces = []
    for i in range(4):
        a, b, e = ctx.mpf(diagonal[i]), ctx.mpf(diagonal[7-i]), ctx.mpf(delta)
        root = ctx.sqrt((a-b)**2+4*e*e)
        p, q = (a+b+root)/2, (a+b-root)/2
        pieces.append(p*ctx.log(p)+q*ctx.log(q)-a*ctx.log(a)-b*ctx.log(b))
    expected = float(ctx.fsum(pieces))
    assert cmi(rho, [2, 2, 2], [0], [1], [2]) == pytest.approx(expected, rel=3e-13, abs=0)


@pytest.mark.parametrize("delta", (1e-10, 1e-50, 1e-75))
def test_fourth_order_information_is_not_confused_with_numerical_zero(delta):
    x = np.array([[0., 1.], [1., 0.]])
    rho = np.eye(4)/4+delta*(np.kron(x, np.eye(2))+np.kron(np.eye(2), x))
    ctx = mp.mp.clone()
    ctx.dps = 650
    e = ctx.mpf(delta)
    ent = lambda ps: -ctx.fsum(p*ctx.log(p) for p in ps)
    expected = float(2*ent([ctx.mpf('.5')+2*e, ctx.mpf('.5')-2*e])
                     -ent([ctx.mpf('.25')+2*e, ctx.mpf('.25'), ctx.mpf('.25'), ctx.mpf('.25')-2*e]))
    assert expected > 0
    assert mi(rho, [2, 2], [0], [1]) == pytest.approx(expected, rel=3e-13, abs=0)


@pytest.mark.parametrize("seed", range(6))
def test_classical_cmi_matches_direct_probability_enumeration(seed):
    rng = np.random.default_rng(400+seed)
    counts = rng.integers(0, 25, size=(2, 3, 2))
    counts[0, 0, 0] += 1
    values = counts.astype(float)/counts.sum()
    ctx = mp.mp.clone()
    ctx.dps = 100
    p = {key: ctx.mpf(values[key]) for key in itertools.product(range(2), range(3), range(2))}
    ab = {(a, b): sum(p[a, b, c] for c in range(2)) for a in range(2) for b in range(3)}
    bc = {(b, c): sum(p[a, b, c] for a in range(2)) for b in range(3) for c in range(2)}
    middle = {b: sum(ab[a, b] for a in range(2)) for b in range(3)}
    expected = ctx.fsum(mass*ctx.log(mass*middle[b]/(ab[a, b]*bc[b, c]))
                       for (a, b, c), mass in p.items() if mass)
    assert cmi(np.diag(values.ravel()), [2, 3, 2], [0], [1], [2]) == pytest.approx(float(expected), rel=2e-13)


def test_trace_roundoff_uses_one_homogeneous_information_convention():
    # A common factor does not manufacture or remove dependence. For MI
    # the empty-B entropy supplies the otherwise missing T log T term.
    diagonal = np.array([.25+1e-9, .25-1e-9, .25-1e-9, .25+1e-9])
    ctx = mp.mp.clone()
    ctx.dps = 100
    p = [ctx.mpf(x) for x in diagonal]
    total = ctx.fsum(p)
    expected = ctx.fsum(x*ctx.log(4*x/total) for x in p)
    assert mi(np.diag(diagonal), [2, 2], [0], [1]) == pytest.approx(float(expected), rel=2e-13, abs=0)
    assert cmi(np.diag(diagonal), [2, 1, 2], [0], [1], [2]) == mi(np.diag(diagonal), [2, 2], [0], [1])
    assert mi((1+2.**-42)*np.eye(4)/4, [2, 2], [0], [1]) == 0


def test_exact_marginal_accumulation_keeps_small_rare_event_information():
    # Adding the rare event to a unit mass first would round it away.
    mass = 1e-310
    rho = np.diag([1., 0., 0., mass])
    ctx = mp.mp.clone()
    ctx.dps = 400
    e = ctx.mpf(mass)
    expected = float((1+e)*ctx.log1p(e)-e*ctx.log(e))
    assert mi(rho, [2, 2], [0], [1]) == pytest.approx(expected, rel=1e-13, abs=0)


def test_interleaved_parts_and_discarded_environment_preserve_the_measure():
    # Bell(A,C) tensor maximally mixed spectators, in A,B,C,D order.
    bell = np.zeros((4, 4)); bell[0, 0] = bell[3, 3] = bell[0, 3] = bell[3, 0] = .5
    rho = np.kron(bell, np.eye(4)/4).reshape([2]*8)
    rho = rho.transpose([0, 2, 1, 3, 4, 6, 5, 7]).reshape(16, 16)
    assert mi(rho, [2]*4, [2], [0]) == pytest.approx(2*np.log(2), abs=1e-14)
    assert cmi(rho, [2]*4, [2], [3, 1], [0]) == pytest.approx(2*np.log(2), abs=1e-14)
    assert mi(rho, [2]*4, [], [0, 2]) == 0


def test_exact_complex_product_and_singular_markov_states_are_retained():
    left = np.array([[.5, .125j], [-.125j, .5]])
    rho = np.kron(left, np.diag([1., 0.]))
    assert mi(rho, [2, 2], [0], [1]) == 0
    assert cmi(rho, [2, 1, 2], [0], [1], [2]) == 0


def test_exact_markov_zero_survives_an_entangling_rechart_of_the_collar():
    # Rational, nontracial factors have irrational eigenvalues. CNOT on the
    # two middle factors hides their original tensor decomposition; no
    # product across a declared subset of B remains.
    left = np.diag([.375, .125, .25, .25]).astype(complex)
    left[0, 3], left[3, 0] = .0625j, -.0625j
    right = np.diag([.125, .375, .375, .125])
    right[0, 3] = right[3, 0] = .0625
    rho = np.kron(left, right)
    permutation = []
    for a, l, r, c in itertools.product((0, 1), repeat=4):
        permutation.append(8*a+4*l+2*(r^l)+c)
    recharted = rho[np.ix_(permutation, permutation)]
    assert cmi(recharted, [2, 2, 2, 2], [0], [1, 2], [3]) == 0
    # A perturbation preserving the observed marginals must not inherit the
    # exact zero certificate from an almost identical characteristic family.
    x = np.array([[0., 1.], [1., 0.]])
    perturbed = recharted+1e-30*np.kron(np.kron(np.kron(x, x), x), x)
    assert cmi(perturbed, [2, 2, 2, 2], [0], [1, 2], [3]) > 0


def test_exact_multisector_markov_state_with_different_collar_splits():
    left = np.diag([.375, .125, .25, .25]).astype(complex)
    left[0, 3] = left[3, 0] = .0625
    right = np.diag([.125, .375, .25, .25]).astype(complex)
    right[0, 3], right[3, 0] = .03125j, -.03125j
    blocks = [np.kron(left, np.diag([.25, .75])),
              np.kron(np.diag([.75, .25]), right)]
    rho = np.zeros((16, 16), complex)
    for j, block in enumerate(blocks):
        indices = [8*a+2*(2*j+b)+c for a, b, c in itertools.product((0, 1), repeat=3)]
        rho[np.ix_(indices, indices)] = .5*block
    h = np.array([[1, 1, 1, 1], [1, -1, 1, -1],
                  [1, 1, -1, -1], [1, -1, -1, 1]])/2
    u = np.kron(np.kron(np.eye(2), h), np.eye(2))
    rho = u@rho@u.T
    assert cmi(rho, [2, 4, 2], [0], [1], [2]) == 0


def test_information_does_not_modify_global_precision():
    saved = mp.mp.dps
    mi(np.eye(4)/4+1e-50*np.fliplr(np.eye(4)), [2, 2], [0], [1])
    assert mp.mp.dps == saved


def test_hermitian_roundoff_averaging_does_not_erase_subnormal_input():
    rho = np.eye(4)/4
    rho[0, 3] = np.nextafter(0., 1.)
    with pytest.raises(ValueError, match="precision|underflow|unresolved"):
        mi(rho, [2, 2], [0], [1])


def test_exact_moment_test_accepts_product_and_rejects_correlated_state():
    from quantum_information.information import _Reductions, _exact_markov
    x = np.array([[0., 1.], [1., 0.]])
    product_state = np.kron(np.eye(2)/2+.125*x, np.diag([.75, .25]))
    assert _exact_markov(_Reductions(product_state, (2, 2)), (0,), (), (1,))
    assert not _exact_markov(_Reductions(np.eye(4)/4+1e-200*np.kron(x, x), (2, 2)), (0,), (), (1,))
    assert is_markov_exact(product_state, [2, 2], [0], [], [1])
    assert not is_markov_exact(np.eye(4)/4+1e-200*np.kron(x, x), [2, 2], [0], [], [1])


def test_exact_moment_reference_retains_complex_coherences():
    from quantum_information.information import _Reductions, _exact_markov
    x = np.array([[0., 1.], [1., 0.]])
    y = np.array([[0., -1j], [1j, 0.]])
    # A is independent of BC, so CMI vanishes. Transposing the complex BC
    # reference changes its modular moments, whereas its B marginal is real.
    rho = np.kron(np.diag([.75, .25]), np.eye(4)/4+.125*np.kron(y, x))
    assert _exact_markov(_Reductions(rho, (2, 2, 2)), (0,), (1,), (2,))


@pytest.mark.parametrize("bad", (np.eye(4), np.diag([1.1, -.1, 0, 0]),
                               np.full((4, 4), np.nan), [[True, 0], [0, 0]],
                               np.ma.array(np.eye(4)/4, mask=True)))
def test_exact_markov_verifier_rejects_invalid_state_evidence(bad):
    with pytest.raises(ValueError):
        is_markov_exact(bad, [2, 2], [0], [], [1])


def test_exact_moment_verifier_against_all_classical_three_bit_supports():
    from quantum_information.information import _Reductions, _exact_markov
    for mask in range(1, 256):
        counts = np.array([(mask >> k) & 1 for k in range(8)]).reshape(2, 2, 2)
        ab, bc, b = counts.sum(axis=2), counts.sum(axis=0), counts.sum(axis=(0, 2))
        expected = all(counts[a, j, c]*b[j] == ab[a, j]*bc[j, c]
                       for a, j, c in itertools.product((0, 1), repeat=3))
        rho = np.diag(counts.ravel()/counts.sum())
        assert _exact_markov(_Reductions(rho, (2, 2, 2)), (0,), (1,), (2,)) == expected, mask


@pytest.mark.parametrize("weights,values", (([True], [1.]), ([1.], [np.nan]),
                                           ([1.], [-1.]), ([1., 0.], [1.]), ([], [])))
def test_weighted_information_rejects_invalid_evidence(weights, values):
    with pytest.raises(ValueError):
        weighted_information(weights, values)


def test_weighted_information_checks_the_total_without_erasing_small_products():
    with pytest.raises(ValueError, match="underflow|precision"):
        weighted_information([np.nextafter(0., 1.)], [.5])
    tiny = np.nextafter(0., 1.)
    # Individually underflowing products may have a representable sum.
    assert weighted_information([tiny, tiny], [.5, .5]) == tiny
    assert weighted_information([.25, .75], [.5, .125]) == float(Fraction(7, 32))


def test_tolerated_negative_spectrum_cannot_pass_a_product_zero_shortcut():
    # Eigenvalues 1+e and -e; the old state tolerance accepts this matrix.
    e = 2.**-45
    bad = np.array([[.5, .5+e], [.5+e, .5]])
    with pytest.raises(ValueError, match="positive semidefinite"):
        mi(np.kron(bad, np.diag([1., 0.])), [2, 2], [0], [1])
    with pytest.raises(ValueError, match="positive semidefinite"):
        mi(np.kron(bad, np.eye(4)/4), [2, 2, 2], [1], [2])
