"""Independent finite checks of deterministic native scan preparation."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement

import numpy as np
import pytest

from . import sweeps as s


def _basis(full=False):
    pairs = combinations_with_replacement(range(12), 2) if full else combinations(range(12), 2)
    result = []
    for a, b in pairs:
        values = [0]*12
        values[a] += 1
        values[b] += 1
        result.append(tuple(values))
    return tuple(sorted(result))


def _oracle(rho, law, *, full=False, reference=1, order=s.DEFAULT_ORDER):
    """Build endpoint balancing Kraus matrices without native/native images."""
    basis = _basis(full)
    positions = {state: i for i, state in enumerate(basis)}
    value = rho.copy()
    for edge in order:
        a, b = s.SEAMS[edge]
        matrices = defaultdict(lambda: np.zeros((len(basis), len(basis))))
        for i, state in enumerate(basis):
            total = state[a]+state[b]
            low, high = total//2, (total+1)//2
            for coin in (0, 1):
                output = list(state)
                output[a], output[b] = (low, high) if coin == 0 else (high, low)
                record = state[a]-output[a] if law == "C" else state[a]-state[b]
                matrices[coin, record][positions[tuple(output)], i] = 1
        updated = np.zeros_like(value)
        for matrix in matrices.values():
            lifted = np.kron(matrix, np.eye(reference))
            updated += lifted @ value @ lifted.T / 2
        value = updated
    return value


def _reference(rho, dimension, reference):
    return np.einsum("ir is->rs", rho.reshape(dimension, reference, dimension, reference))


def _half_trace(matrix):
    return abs(np.linalg.eigvalsh(matrix)).sum()/2


def test_local_native_and_c0_identities_are_exact():
    report = s.structural_audit()
    assert report["checked_local_matrix_units"] == 30*66*66
    assert report["nonzero_difference_sectors"] == 3102
    assert report["C_surviving_cut_minimum"] == 10
    assert report["C0_weighted_cut_minimum"] == 10
    assert report["C_offdiagonal_per_sweep"] == F(1, 1024)
    assert report["all_offdiagonal_per_sweep"] == F(1, 32)
    # Return values cannot poison the cached structural check.
    report["C_surviving_cut_minimum"] = -1
    assert s.structural_audit()["C_surviving_cut_minimum"] == 10


@pytest.mark.parametrize("order", [s.DEFAULT_ORDER, s.DEFAULT_ORDER[::-1]])
def test_exact_kernel_matches_independent_heat_bath_order(order):
    basis = _basis()
    positions = {state: i for i, state in enumerate(basis)}
    integer = np.eye(66, dtype=object)
    denominator = 1
    for _ in range(2):
        for edge_index in order:
            a, b = s.SEAMS[edge_index]
            permutation = []
            for state in basis:
                swapped = list(state)
                swapped[a], swapped[b] = swapped[b], swapped[a]
                permutation.append(positions[tuple(swapped)])
            integer = integer + integer[permutation, :]
            denominator *= 2
    actual = s.classical_power(2, order)
    assert actual.denominator == denominator
    assert actual.counts == tuple(tuple(int(x) for x in row) for row in integer)
    # Ordered products are generally not reversible.  No eigenvalue gap of
    # the averaged, symmetric random-seam channel can replace this audit.
    assert not np.array_equal(integer, integer.T)
    with pytest.raises(TypeError):
        actual.counts[0][0] = 0


@pytest.mark.parametrize("law", ["C", "E"])
@pytest.mark.parametrize("full", [False, True])
def test_actual_native_execution_matches_independent_kraus_and_keeps_reference(law, full):
    dimension = 78 if full else 66
    rng = np.random.default_rng(314)
    vector = rng.normal(size=dimension*2)+1j*rng.normal(size=dimension*2)
    vector /= np.linalg.norm(vector)
    initial = np.outer(vector, vector.conj())
    actual = s.execute(initial, law, reference_dimension=2, full=full)
    expected = _oracle(initial, law, full=full, reference=2)
    np.testing.assert_allclose(actual, expected, atol=3e-17)
    np.testing.assert_allclose(_reference(actual, dimension, 2), _reference(initial, dimension, 2), atol=2e-16)
    assert abs(np.trace(actual)-np.trace(initial)) < 2e-14
    if full:
        collision_indices = [i*2+r for i, state in enumerate(s.FULL_STATES) if max(state) == 2 for r in range(2)]
        np.testing.assert_array_equal(actual[collision_indices, :], 0)
        np.testing.assert_array_equal(actual[:, collision_indices], 0)


def test_c_offdiagonal_bound_is_attained_and_e_dephases():
    edge_set = {frozenset(edge) for edge in s.SEAMS}
    p, q = next(pair for pair in combinations(range(12), 2) if frozenset(pair) not in edge_set)
    index = {state: i for i, state in enumerate(s.BINARY_STATES)}
    perturbation = np.zeros((66, 66))
    for common in set(range(12))-{p, q}:
        left, right = [0]*12, [0]*12
        left[p] = left[common] = right[q] = right[common] = 1
        i, j = index[tuple(left)], index[tuple(right)]
        perturbation[i, j] = perturbation[j, i] = 1/132
    tau = np.eye(66)/66
    initial = tau+perturbation
    c_output = s.execute(initial, "C")
    e_output = s.execute(initial, "E")
    np.testing.assert_allclose(c_output, tau+perturbation/1024, atol=1e-18)
    np.testing.assert_allclose(e_output, tau, atol=1e-18)
    assert abs(np.linalg.norm(c_output-tau)/np.linalg.norm(perturbation)-1/1024) < 1e-18


def test_all_law_envelope_is_attained_and_survives_coarse_record_phase_feedback():
    """Reconstruct C0 and an admitted adaptive phase law from endpoint repair."""
    basis, positions = _basis(), {state: i for i, state in enumerate(_basis())}
    p, q = 0, 1
    perturbation = np.zeros((66, 66), complex)
    for common in set(range(12))-{p, q}:
        left, right = [0]*12, [0]*12
        left[p] = left[common] = right[q] = right[common] = 1
        i, j = positions[tuple(left)], positions[tuple(right)]
        perturbation[i, j] = perturbation[j, i] = 1/132
    tau = np.eye(66)/66
    initial = tau+perturbation
    envelope = initial.copy()
    histories = {False: initial.copy()}
    for a, b in s.SEAMS:
        matrices = defaultdict(lambda: np.zeros((66, 66), complex))
        for i, state in enumerate(basis):
            total = state[a]+state[b]
            low, high = total//2, (total+1)//2
            probabilities = defaultdict(float)
            for left, right in ((low, high), (high, low)):
                output = list(state)
                output[a], output[b] = left, right
                probabilities[state[a]-left, positions[tuple(output)]] += .5
            for (current, output), probability in probabilities.items():
                matrices[current][output, i] = np.sqrt(probability)
        envelope = sum(matrix @ envelope @ matrix.conj().T for matrix in matrices.values())
        following = {False: np.zeros_like(initial), True: np.zeros_like(initial)}
        for previous_positive, rho in histories.items():
            # The diagonal phase depends only on local total, hence obeys
            # endpoint-swap and local number-phase covariance.  Its choice
            # reads the last coarse current, never hidden loads or a seam.
            phases = np.array([({0: 1, 1: -1j, 2: 1j}[state[a]+state[b]]
                                if previous_positive else 1) for state in basis])
            for current, matrix in matrices.items():
                kraus = matrix*phases[None, :]
                following[current > 0] += kraus @ rho @ kraus.conj().T
        histories = following
    adaptive = sum(histories.values())
    np.testing.assert_allclose(envelope, tau+perturbation/32, atol=3e-16)
    # Absolute coefficients of the feedback history remain bounded by the
    # positive C0 path sum, even though the accumulated phases can cancel.
    assert np.max(abs(adaptive)-envelope.real) <= 3e-16
    np.testing.assert_allclose(np.diag(adaptive), np.diag(envelope), atol=3e-16)
    assert np.linalg.norm(adaptive-np.diag(np.diag(adaptive))) <= np.linalg.norm(perturbation)/32+3e-16


@pytest.mark.parametrize("count", [0, 1, 4, 8, 16])
def test_rational_norm_enclosures_cover_singular_values(count):
    power = s.classical_power(count)
    centered = np.array(power.counts, dtype=float)/power.denominator - np.ones((66, 66))/66
    numerical = np.linalg.svd(centered, compute_uv=False)[0]
    report = s.certificate(count)
    assert numerical <= float(report["diagonal_hilbert_schmidt_upper"])+2e-14
    product = report["classical_row_l1_upper"]*report["classical_column_l1_upper"]
    if product < 1:
        assert report["diagonal_hilbert_schmidt_upper"]**2 >= product


@pytest.mark.parametrize("law", ["C", "E", "all"])
def test_budget_counts_all_initial_absorption_and_mixing_calls(law):
    tolerance = F(1, 1570800)
    report = s.budget(tolerance, law, full=True)
    assert report["sweeps"] == 17
    assert report["binary_mixing_sweeps"] == 16
    assert report["native_attempts"] == 510
    assert report["collision_absorption_sweeps"] == 1
    assert report["trace_distance_upper"] <= tolerance
    assert s.certificate(16, law, full=True)["trace_distance_upper"] > tolerance
    assert s.budget(F(1, 1000), law)["sweeps"] == 8
    assert s.budget(F(1, 1000), law, full=True)["sweeps"] == 9


@pytest.mark.parametrize("law", ["C", "E"])
def test_reference_bound_covers_entangled_input_and_full_collision_input(law):
    rng = np.random.default_rng(2718)
    vector = rng.normal(size=78*3)+1j*rng.normal(size=78*3)
    vector /= np.linalg.norm(vector)
    initial = np.outer(vector, vector.conj())
    reference = _reference(initial, 78, 3)
    tau = np.diag([float(max(state) == 1)/66 for state in s.FULL_STATES])
    target = np.kron(tau, reference)
    output = s.execute(initial, law, sweeps=5, reference_dimension=3, full=True)
    assert _half_trace(output-target) <= float(s.certificate(5, law, full=True)["trace_distance_upper"])


@pytest.mark.parametrize("mass", [0, 2**-40, .375])
def test_subnormalized_branches_are_executed_linearly_without_reset(mass):
    vector = np.zeros(78*2)
    vector[0] = vector[-1] = 1/np.sqrt(2)
    rho = mass*np.outer(vector, vector)
    output = s.execute(rho, "C", sweeps=2, reference_dimension=2, full=True)
    unit = s.execute(np.outer(vector, vector), "C", sweeps=2, reference_dimension=2, full=True)
    np.testing.assert_allclose(output, mass*unit, atol=1e-30)
    assert abs(np.trace(output)-mass) <= 4e-16*mass


@pytest.mark.parametrize("bad", [True, -1, 65, 1.0, None])
def test_bad_sweep_counts_fail_closed(bad):
    with pytest.raises(ValueError):
        s.certificate(bad)


@pytest.mark.parametrize("bad", [(), list(range(29)), list(range(29))+[28], [True]+list(range(1, 30))])
def test_incomplete_or_adaptive_scan_is_not_certified(bad):
    with pytest.raises(ValueError):
        s.classical_power(1, bad)


@pytest.mark.parametrize("bad", [True, 0, -F(1, 2), 1.1, float("nan"), "0.01"])
def test_nonexact_or_invalid_accuracy_is_rejected(bad):
    with pytest.raises(ValueError):
        s.budget(bad)


def test_execution_validation_and_no_implicit_all_law_implementation():
    tau = np.eye(66)/66
    for keyword in ({"law": "all"}, {"full": 1}, {"reference_dimension": True}):
        with pytest.raises(ValueError):
            s.execute(tau, **keyword)
    for bad in (tau*2, np.full((66, 66), np.nan), np.eye(2), tau+1j*np.eye(66),
                np.diag([-1]+[2/65]*65), tau.tolist(), tau.astype(object),
                tau.astype(bool), np.ma.array(tau)):
        with pytest.raises(ValueError):
            s.execute(bad)
