"""Independent finite references and boundary-order checks for certificates."""
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.spatial_interaction import bounds as b


def _cube_reference(hopping, interaction):
    """Build bosons by symmetrizing two distinguishable cube particles."""
    sites = tuple(product(range(2), repeat=3))
    count = len(sites)
    pairs = tuple(combinations_with_replacement(range(count), 2))
    sym = np.zeros((count**2, len(pairs)))
    for column, (x, y) in enumerate(pairs):
        sym[x * count + y, column] = 1 if x == y else 1 / np.sqrt(2)
        if x != y:
            sym[y * count + x, column] = 1 / np.sqrt(2)
    laplacians = [np.zeros((count, count)) for _ in range(6)]
    for x, point in enumerate(sites):
        for axis in range(3):
            neighbor = list(point)
            neighbor[axis] += 1
            if neighbor[axis] >= 2:
                continue
            y = sites.index(tuple(neighbor))
            lap = laplacians[2 * axis + point[axis] % 2]
            lap[x, x] += 1
            lap[y, y] += 1
            lap[x, y] -= 1
            lap[y, x] -= 1
    terms = [float(hopping) * sym.T @ (
        np.kron(lap, np.eye(count)) + np.kron(np.eye(count), lap)) @ sym
        for lap in laplacians]
    contact = np.diag([float(interaction) if x == y else 0 for x, y in pairs])
    terms.append(contact)
    origin = sites.index((0, 0, 0))
    neighbor = sites.index((1, 0, 0))
    doublon = pairs.index((origin, origin))
    split = pairs.index(tuple(sorted((origin, neighbor))))
    initial = np.zeros(len(pairs), complex)
    initial[doublon] = initial[split] = 1 / np.sqrt(2)
    phi = np.zeros_like(initial)
    phi[doublon], phi[split] = 1 / np.sqrt(2), 1j / np.sqrt(2)
    projector = np.outer(phi, phi.conj())
    feedback = np.eye(len(pairs), dtype=complex)
    feedback[split, split] = 1j
    return terms, initial, projector, feedback


def _split(terms, time, steps):
    size = len(terms[0])
    half = [expm(-1j * float(time) / (2 * steps) * term) for term in terms[:-1]]
    word = half + [expm(-1j * float(time) / steps * terms[-1])] + half[::-1]
    step = np.eye(size, dtype=complex)
    for unitary in word:
        step = unitary @ step
    return np.linalg.matrix_power(step, steps)


@pytest.mark.parametrize("hopping,interaction", [(F(1, 8), F(3, 4)), (F(0), F(1))])
def test_strang_certificate_dominates_independent_symmetric_cube_operator(hopping, interaction):
    terms, *_ = _cube_reference(hopping, interaction)
    target = sum(terms)
    norm_sum = sum(np.linalg.norm(term, 2) for term in terms)
    assert norm_sum <= float(b.generator_norm_bound(hopping, interaction)) + 2e-14
    for steps in (1, 2, 4):
        time = F(1, 10)
        actual = _split(terms, time, steps)
        error = np.linalg.norm(actual - expm(-1j * float(time) * target), 2)
        receipt = b.strang_bound(steps, time, hopping, interaction)
        assert error <= float(F(receipt["operator_error_upper_exact"])) + 2e-14
        json.dumps(receipt, allow_nan=False)


def test_whole_history_bound_controls_retained_quantum_outputs_without_branch_factor():
    hopping, interaction = F(1, 8), F(3, 4)
    terms, initial, projector, feedback = _cube_reference(hopping, interaction)
    target = sum(terms)
    segments = [(F(1, 16), 2), (F(1, 16), 3)]
    projectors = (np.eye(len(initial)) - projector, projector)
    protocols = []
    for approximate in (False, True):
        unitaries = [(_split(terms, time, steps) if approximate
                     else expm(-1j * float(time) * target)) for time, steps in segments]
        branches = []
        for first, second in product((0, 1), repeat=2):
            state = projectors[first] @ unitaries[0] @ initial
            if first:
                state = feedback @ state
            branches.append(projectors[second] @ unitaries[1] @ state)
        assert sum(np.vdot(v, v).real for v in branches) == pytest.approx(1, abs=2e-14)
        protocols.append(branches)
    distance = 0.
    total_variation = 0.
    for exact, approximate in zip(*protocols):
        difference = np.outer(exact, exact.conj()) - np.outer(approximate, approximate.conj())
        distance += np.sum(np.abs(np.linalg.eigvalsh(difference))) / 2
        total_variation += abs(np.vdot(exact, exact).real - np.vdot(approximate, approximate).real) / 2
    receipt = b.history_bound(segments, hopping, interaction)
    assert total_variation <= distance + 1e-14
    assert distance <= float(F(receipt["complete_history_tv_upper_exact"])) + 2e-14
    assert receipt["branch_count_multiplier"] == 1
    expected = sum(F(row["operator_error_upper_exact"]) for row in receipt["segments"])
    assert F(receipt["complete_history_tv_upper_exact"]) == expected


def _ordered_initial():
    origin, neighbor = (0, 0, 0), (1, 0, 0)
    return {(origin, origin): 1 / np.sqrt(2),
            (origin, neighbor): .5, (neighbor, origin): .5}


def _neighbors(point, radius):
    for axis in range(3):
        for sign in (-1, 1):
            output = list(point)
            output[axis] += sign
            if max(map(abs, output)) <= radius:
                yield tuple(output)


def _apply_h(vector, radius, hopping=.125, interaction=.75):
    """Sparse ordered-coordinate action, independent of the bosonic code."""
    output = {}

    def add(pair, amplitude):
        output[pair] = output.get(pair, 0) + amplitude

    for (x, y), amplitude in vector.items():
        nx, ny = tuple(_neighbors(x, radius)), tuple(_neighbors(y, radius))
        add((x, y), (hopping * (len(nx) + len(ny)) + interaction * (x == y)) * amplitude)
        for neighbor in nx:
            add((neighbor, y), -hopping * amplitude)
        for neighbor in ny:
            add((x, neighbor), -hopping * amplitude)
    return output


def _ramsey_branch(vector):
    origin, neighbor = (0, 0, 0), (1, 0, 0)
    phi = {(origin, origin): 1 / np.sqrt(2),
           (origin, neighbor): .5j, (neighbor, origin): .5j}
    amplitude = sum(np.conj(value) * vector.get(pair, 0) for pair, value in phi.items())
    return {pair: value * amplitude for pair, value in phi.items()}


def _distance(left, right):
    return math.sqrt(sum(abs(left.get(pair, 0) - right.get(pair, 0))**2
                         for pair in left.keys() | right.keys()))


@pytest.mark.parametrize("radius", [1, 2, 3])
def test_boundary_first_unequal_degree_is_radius_not_radius_plus_one(radius):
    small = large = _ordered_initial()
    for degree in range(radius):
        assert _distance(small, large) < 1e-14
        small, large = _apply_h(small, radius), _apply_h(large, radius + 1)
    assert _distance(small, large) > 1e-6
    receipt = b.boundary_bound(radius, F(1, 8), F(1, 8), F(3, 4))
    assert receipt["first_potentially_unequal_taylor_order"] == radius
    assert receipt["distance_from_apparatus_to_boundary"] == radius - 1


def test_boundary_coefficient_locality_survives_inserted_ramsey_mixing():
    for prefix, suffix in product(range(3), repeat=2):
        if prefix + suffix > 2:
            continue
        states = []
        for radius in (3, 4):
            state = _ordered_initial()
            for _ in range(prefix):
                state = _apply_h(state, radius)
            state = _ramsey_branch(state)
            for _ in range(suffix):
                state = _apply_h(state, radius)
            states.append(state)
        assert _distance(*states) < 1e-14


def test_short_time_interaction_witness_has_correct_readout_sign_and_positive_margin():
    hopping, interaction = F(1, 8), F(3, 4)
    terms, initial, projector, _ = _cube_reference(hopping, interaction)
    receipt = b.interaction_witness_bound(hopping, interaction)
    time = float(F(receipt["model_time_exact"]))
    probabilities = []
    for target in (sum(terms[:-1]), sum(terms)):
        state = expm(-1j * time * target) @ initial
        probabilities.append(np.vdot(state, projector @ state).real)
    lower = F(receipt["interaction_minus_free_probability_lower_exact"])
    assert probabilities[1] - probabilities[0] >= float(lower)
    assert lower == 11 * interaction * F(receipt["model_time_exact"]) / 30
    # The opposite Ramsey orientation has the opposite linear response.
    opposite = projector.conj()
    rho = np.outer(initial, initial.conj())
    derivative = -1j * (terms[-1] @ rho - rho @ terms[-1])
    assert np.trace(projector @ derivative).real == pytest.approx(float(interaction) / 2)
    assert np.trace(opposite @ derivative).real == pytest.approx(-float(interaction) / 2)


def test_exact_receipts_handle_trivial_protocols_and_clip_distances_not_evidence():
    assert b.generator_norm_bound(0, 0) == 0
    assert F(b.strang_bound(1, 0, 1, 1)["operator_error_upper_exact"]) == 0
    assert F(b.history_bound([], 1, 1)["complete_history_tv_upper_exact"]) == 0
    huge = b.boundary_bound(1, 10, 1, 1)
    assert F(huge["unclipped_stinespring_vector_error_upper_exact"]) == 500
    assert F(huge["complete_history_tv_upper_exact"]) == 1
    gap = b.gap_bound(F(3, 5), F(1, 8), F(3, 4))
    assert gap["vacuum_spectral_gap_exact"] == "3/5"
    assert not gap["gap_measured_by_fixed_two_particle_protocol"]
    assert not gap["gap_generated_by_native_repair"]
    json.dumps(gap, allow_nan=False)


@pytest.mark.parametrize("bad", [True, False, .125, float("nan"), float("inf"), -1, "1/8", None])
def test_certificates_reject_inexact_or_invalid_parameters(bad):
    with pytest.raises(ValueError):
        b.generator_norm_bound(bad, 1)
    with pytest.raises(ValueError):
        b.generator_norm_bound(1, bad)
    with pytest.raises(ValueError):
        b.strang_bound(1, bad, 1, 1)


@pytest.mark.parametrize("bad", [True, 0, -1, 1., "1", None])
def test_step_counts_and_radii_are_positive_integers(bad):
    with pytest.raises(ValueError):
        b.strang_bound(bad, 1, 1, 1)
    with pytest.raises(ValueError):
        b.boundary_bound(bad, 1, 1, 1)


@pytest.mark.parametrize("bad", [None, iter([]), [1], [(1,)], [(1, 2, 3)]])
def test_history_requires_explicit_finite_segment_pairs(bad):
    with pytest.raises(ValueError):
        b.history_bound(bad, 1, 1)


def test_strictly_positive_claims_reject_zero_gap_or_zero_interaction():
    with pytest.raises(ValueError):
        b.gap_bound(0, 1, 1)
    with pytest.raises(ValueError):
        b.interaction_witness_bound(1, 0)
