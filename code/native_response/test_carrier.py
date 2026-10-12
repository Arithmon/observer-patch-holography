"""Independent native instrument, operational address and coherence controls."""
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pytest

from native_response.carrier import (
    SeamChannels, antipode, binary_query, cq_state, decode, encode,
    refresh_address, select_address, star_repair, tick,
)
from native_response.source import edges


def random_cq(dimension, seed=1087):
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(12, dimension, dimension)) + 1j * rng.normal(
        size=(12, dimension, dimension))
    state = raw @ raw.conj().transpose(0, 2, 1)
    return state / np.trace(state, axis1=1, axis2=2).sum()


def swap_channels():
    result = {}
    for a, b in edges():
        swap = np.eye(12)
        swap[[a, b]] = swap[[b, a]]
        result[a, b] = np.array([np.eye(12), swap]) / np.sqrt(2)
    return SeamChannels(result)


def native_binary():
    """Reconstruct grouped (coin,current) Kraus maps from integer repair."""
    states = np.array([[int(p in occupied) for p in range(12)]
                       for occupied in combinations(range(12), 2)])
    index = {tuple(state): i for i, state in enumerate(states)}
    families, transitions = {}, {}
    for a, b in edges():
        branches = {}
        kernel = np.zeros((66, 66))
        for i, state in enumerate(states):
            for coin in (0, 1):
                # Use floor/ceiling placement, independently of the response code.
                output = state.copy()
                total = int(state[a] + state[b])
                output[a] = (total + coin) // 2
                output[b] = total - output[a]
                current = int(state[a] - output[a])
                j = index[tuple(output)]
                branches.setdefault((coin, current), np.zeros((66, 66)))[j, i] = 1 / np.sqrt(2)
                kernel[j, i] += .5
        families[a, b] = np.array(list(branches.values()))
        transitions[a, b] = kernel
    return states, SeamChannels(families), transitions


def test_minimal_attachment_retains_native_coherence_and_correlated_inputs():
    state = random_cq(3)
    marginal = state.sum(axis=0)
    np.testing.assert_allclose(decode(state), marginal, atol=1e-15)
    np.testing.assert_allclose(decode(encode(marginal)), marginal, atol=1e-15)
    np.testing.assert_allclose(decode(refresh_address(state)), marginal, atol=1e-15)
    assert abs(marginal[0, 1]) > .01
    assert np.linalg.norm(state - encode(marginal)) > .01
    selected = select_address(state, 7)
    np.testing.assert_allclose(selected[7], marginal, atol=1e-15)
    assert np.count_nonzero(np.trace(selected, axis1=1, axis2=2)) == 1


def test_native_66_global_attempts_match_independent_classical_kernel():
    states, channels, kernels = native_binary()
    rng = np.random.default_rng(78)
    weights = rng.random((12, 66))
    weights /= weights.sum()
    state = np.array([np.diag(row) for row in weights])
    population = weights.sum(axis=0)
    transition = sum(kernels.values()) / 30
    for _ in range(2):
        state = tick(state, channels)
        population = transition @ population
        np.testing.assert_allclose(decode(state), np.diag(population), atol=2e-15)
    # Fresh address is essential: this initial address chooses the local star.
    rho = np.diag((states[:, 0] * states[:, 2]).astype(float))
    local = decode(star_repair(select_address(encode(rho), 0), channels))
    expected = sum(kernels[e] @ np.diag(rho) for e in edges() if 0 in e) / 5
    np.testing.assert_allclose(local, np.diag(expected), atol=2e-15)
    other = decode(star_repair(select_address(encode(rho), antipode(0)), channels))
    assert np.linalg.norm(local - other) > .1


def test_intertwining_multiple_steps_on_quantum_correlated_cq_states():
    channels = swap_channels()
    state = random_cq(12)
    native = state.sum(axis=0)
    for _ in range(3):
        state = tick(state, channels)
        native = sum(channels.apply(edge, native) for edge in edges()) / 30
        np.testing.assert_allclose(decode(state), native, atol=2e-15)
    # Carrying the correlated address forward instead changes the schedule.
    initial = random_cq(12, seed=1088)
    assert np.linalg.norm(decode(star_repair(initial, channels))
                          - decode(tick(initial, channels))) > .001


def test_reference_entanglement_and_conditioned_histories_survive_attachment():
    """The adapter intertwines instruments even with an untouched reference.

    Nonunital, seam-dependent native channels prevent a tracial or uniform
    unitary special case from standing in for the complete channel identity.
    Every branch is retained with its probability; no postselection is used.
    """
    reference = np.eye(2)
    families = {}
    for index, edge in enumerate(edges()):
        probability = .02 + .03 * (index % 5)
        native = (np.diag([1, np.sqrt(1 - probability)]),
                  np.array([[0, np.sqrt(probability)], [0, 0]]))
        families[edge] = [np.kron(kraus, reference) for kraus in native]
    channels = SeamChannels(families)
    plus = np.array([1, 0, 0, 1]) / np.sqrt(2)
    minus = np.array([1, 0, 0, -1]) / np.sqrt(2)
    initial = np.zeros((12, 4, 4), complex)
    initial[0] = .8 * np.outer(plus, plus)
    initial[7] = .2 * np.outer(minus, minus)
    native = initial.sum(axis=0)
    partial_transpose = native.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1)
    assert np.linalg.eigvalsh(partial_transpose.reshape(4, 4))[0] < -.29
    instrument = [np.kron(np.diag(np.sqrt(weights)), reference)
                  for weights in ([.7, .3], [.3, .7])]
    branches = [(1., initial, native)]
    for _ in range(2):
        continued = []
        for weight, state, expected in branches:
            actual = tick(state, channels)
            expected = sum(k @ expected @ k.conj().T
                           for family in families.values() for k in family) / 30
            for kraus in instrument:
                conditional = kraus @ actual @ kraus.conj().T
                probability = np.trace(conditional, axis1=1, axis2=2).sum().real
                conditional /= probability
                expected_branch = kraus @ expected @ kraus.conj().T
                branch_weight = weight * probability
                np.testing.assert_allclose(branch_weight * decode(conditional),
                                           expected_branch, atol=2e-15, rtol=0)
                continued.append((branch_weight, conditional, expected_branch))
        branches = continued
    assert len(branches) == 4
    assert sum(weight for weight, _, _ in branches) == pytest.approx(1)


def test_central_query_records_lueders_branches_and_preserves_within_branch_coherence():
    state = random_cq(3)
    native = decode(state)
    projector = np.diag([1, 1, 0])
    projectors = np.repeat(projector[None, :, :], 12, axis=0)
    port = 4
    result = binary_query(select_address(state, port), projectors)
    complement = np.eye(3) - projector
    np.testing.assert_allclose(result[port], projector @ native @ projector, atol=1e-15)
    np.testing.assert_allclose(result[antipode(port)], complement @ native @ complement,
                               atol=1e-15)
    assert abs(result[port, 0, 1]) > .01
    assert np.trace(result[port]).real == pytest.approx(np.trace(projector @ native).real)
    assert np.linalg.norm(decode(result) - native) > .01


def test_intermediate_read_is_a_different_quantum_experiment():
    angle = np.pi / 8
    unitary = np.cos(angle) * np.eye(2) - 1j * np.sin(angle) * np.array([[0, 1], [1, 0]])
    channels = SeamChannels({edge: unitary[None, :, :] for edge in edges()})
    projectors = np.repeat(np.diag([1, 0])[None, :, :], 12, axis=0)
    initial = encode(np.diag([1, 0]))
    first = tick(initial, channels)
    uninterrupted = decode(tick(first, channels))
    measured = decode(tick(binary_query(select_address(first, 0), projectors), channels))
    assert uninterrupted[1, 1].real == pytest.approx(.5)
    assert measured[1, 1].real == pytest.approx(.25)


def test_retained_endpoint_records_have_the_same_quantum_joint_law():
    angle = np.pi / 8
    unitary = np.cos(angle) * np.eye(2) - 1j * np.sin(angle) * np.array([[0, 1], [1, 0]])
    channels = SeamChannels({edge: unitary[None, :, :] for edge in edges()})
    projectors = np.repeat(np.diag([1, 0])[None, :, :], 12, axis=0)
    query = 3
    records = (antipode(query), query)  # Empty and occupied central outputs.
    initial = binary_query(select_address(encode(np.eye(2) / 2), query), projectors)
    for count in (1, 2):
        joint = np.zeros((2, 2))
        for first, port in enumerate(records):
            weight = np.trace(initial[port]).real
            branch = np.zeros_like(initial)
            branch[port] = initial[port] / weight
            for _ in range(count):
                branch = tick(branch, channels)
            final = binary_query(select_address(branch, query), projectors)
            for last, destination in enumerate(records):
                joint[first, last] = weight * np.trace(final[destination]).real
        stay = np.cos(count * angle)**2
        np.testing.assert_allclose(joint, .5 * np.array([[stay, 1 - stay],
                                                       [1 - stay, stay]]), atol=2e-15)


def test_all_proper_actions_covary_address_repair_and_contextual_query():
    packet = json.loads((Path(__file__).resolve().parents[1]
                        / "source_selection_model/response.json").read_text())
    channels = swap_channels()
    state = random_cq(12)
    projectors = np.array([np.diag(np.eye(12)[p]) for p in range(12)])
    repaired, queried = star_repair(state, channels), binary_query(state, projectors)

    def transport(value, action):
        result = np.zeros_like(value)
        for p in range(12):
            result[action[p]][np.ix_(action, action)] = value[p]
        return result

    for row in packet["closed_paths"]:
        action = row["port_action"]
        assert all(antipode(action[p]) == action[antipode(p)] for p in range(12))
        moved = transport(state, action)
        np.testing.assert_allclose(star_repair(moved, channels), transport(repaired, action),
                                   atol=2e-15)
        np.testing.assert_allclose(binary_query(moved, projectors), transport(queried, action),
                                   atol=2e-15)


@pytest.mark.parametrize("bad", [np.zeros((12, 2, 2)), np.ones((11, 1, 1)) / 11,
                                   np.full((12, 1, 1), np.nan),
                                   np.full((12, 1, 1), True),
                                   np.ma.array(np.ones((12, 1, 1)) / 12)])
def test_malformed_cq_input_rejected(bad):
    with pytest.raises(ValueError):
        cq_state(bad)


def test_invalid_quantum_blocks_maps_and_query_are_rejected():
    bad = np.zeros((12, 2, 2), dtype=complex)
    bad[0] = np.diag([1.1, -.1])
    with pytest.raises(ValueError):
        cq_state(bad)
    bad[0] = np.array([[.5, .1], [0, .5]])
    with pytest.raises(ValueError):
        cq_state(bad)
    valid = encode(np.eye(2) / 2)
    with pytest.raises(ValueError):
        select_address(valid, True)
    with pytest.raises(ValueError):
        SeamChannels({edge: [np.eye(2) / 2] for edge in edges()})
    with pytest.raises(ValueError):
        SeamChannels({(0, 2): [np.eye(2)]})
    malformed_keys = {edge: [np.eye(2)] for edge in edges()}
    malformed_keys[(False, 2)] = malformed_keys.pop((0, 2))
    with pytest.raises(ValueError):
        SeamChannels(malformed_keys)
    with pytest.raises(ValueError):
        binary_query(valid, np.repeat((np.eye(2) / 2)[None, :, :], 12, axis=0))
    with pytest.raises(ValueError):
        star_repair(valid, swap_channels())


def test_supplied_kraus_arrays_are_copied_before_use():
    families = {edge: np.eye(2)[None, :, :] for edge in edges()}
    channels = SeamChannels(families)
    families[(0, 2)][0] = 0
    state = encode(np.array([[.5, .25j], [-.25j, .5]]))
    np.testing.assert_allclose(decode(tick(state, channels)), decode(state), atol=2e-15)


def test_uniform_attachment_refuses_subnormal_population_erasure():
    tiny = np.nextafter(0., 1.)
    with pytest.raises(ValueError, match="underflows"):
        encode(np.diag([1., tiny]))


@pytest.mark.skipif(np.finfo(np.longdouble).tiny >= np.finfo(float).tiny,
                    reason="platform has no extended exponent range")
@pytest.mark.parametrize("imaginary", [False, True])
def test_conversion_refuses_extended_precision_component_erasure(imaginary):
    tiny = np.finfo(np.longdouble).tiny
    state = np.eye(2, dtype=np.clongdouble) / 2
    if imaginary:
        state[0, 1] = np.clongdouble(1j) * tiny
        state[1, 0] = -np.clongdouble(1j) * tiny
    else:
        state[0, 1] = state[1, 0] = tiny
    with pytest.raises(ValueError, match="underflows"):
        encode(state)
