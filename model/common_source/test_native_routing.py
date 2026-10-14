"""Independent all-column routing and complete recorded native-channel checks."""
import math

import numpy as np
import pytest

from ..observer import native
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from . import native_routing as r


def expected_route():
    """Independent physical port swap, with the complement fixed pointwise."""
    expected = np.eye(512, dtype=complex)
    for state in native.states():
        swapped = tuple(state[5] if port == 1 else state[1] if port == 5 else value
                        for port, value in enumerate(state))
        column = circuit.encode(state)
        expected[:, column] = 0
        expected[circuit.encode(swapped), column] = 1
    return expected


def source_columns():
    columns = np.zeros((512, 78), complex)
    columns[circuit.SOURCE_CODES, np.arange(78)] = 1
    return columns


def test_route_is_the_full_phase_exact_source_port_swap_with_padding_retained():
    expected = expected_route()
    forward, backward = r.routing_unitary(), r.routing_unitary(inverse=True)
    for actual in (forward, backward):
        # Absolute matrix comparison includes phases and every padding column.
        assert np.linalg.norm(actual-expected, 2) < 2e-10
        assert np.linalg.norm(actual.conj().T@actual-np.eye(512), 2) < 3e-13
        with pytest.raises(ValueError):
            actual.setflags(write=True)
    assert np.linalg.norm(backward@forward-np.eye(512), 2) < 3e-10
    assert len(r.routing_pairs()) == 11
    assert len(r.routing_edges()) == 93
    assert r.routing_edges(inverse=True) == tuple(reversed(r.routing_edges()))


@pytest.mark.parametrize('candidate', ['C', 'E'])
@pytest.mark.parametrize('seam', [(0, 1), (0, 5)])
def test_complete_routed_repair_maps_match_independent_native_floor_ceiling(candidate, seam):
    branches = r.routed_repair_branches(candidate, seam, source_columns())
    raw_count = 1 << circuit.RECORD_BITS[candidate]
    assert [(b['coin'], b['record_code']) for b in branches] == [
        (coin, code) for coin in (0, 1) for code in range(raw_count)]
    gram = np.zeros((78, 78), complex)
    for branch in branches:
        coin, raw = branch['coin'], branch['record_code']
        expected = np.zeros((512, 78), complex)
        for column, state in enumerate(native.states()):
            a, b = seam
            total = state[a]+state[b]
            low, high = total//2, (total+1)//2
            left, right = (low, high) if coin == 0 else (high, low)
            current = state[a]-left
            label = {0: 0, 1: 1, -1: 2}[current]
            if candidate == 'E':
                label += 4*(total % 2)
            if label == raw:
                repaired = list(state)
                repaired[a], repaired[b] = left, right
                expected[circuit.encode(tuple(repaired)), column] = 1/math.sqrt(2)
        actual = branch['amplitude']
        assert actual.shape == (512, 78)
        assert branch['seam'] == seam
        assert np.linalg.norm(actual-expected, 2) < 3e-10
        gram += actual.conj().T@actual
        if branch['abort']:
            assert np.linalg.norm(actual) < 3e-10
        with pytest.raises(ValueError):
            actual.setflags(write=True)
    assert np.linalg.norm(gram-np.eye(78), 2) < 5e-10


def test_routing_preserves_reference_coherence_and_does_not_replace_the_source():
    rng = np.random.default_rng(1100)
    original = rng.normal(size=(78, 3))+1j*rng.normal(size=(78, 3))
    original /= np.linalg.norm(original)
    encoded = source_columns()@original
    full_maps = r.routed_repair_branches('C', (0, 5), source_columns())
    entangled = r.routed_repair_branches('C', (0, 5), encoded)
    reduced = np.zeros((3, 3), complex)
    for branch, matrix in zip(entangled, full_maps, strict=True):
        np.testing.assert_allclose(branch['amplitude'], matrix['amplitude']@original,
                                   atol=3e-13, rtol=0)
        reduced += branch['amplitude'].conj().T@branch['amplitude']
    np.testing.assert_allclose(reduced, original.conj().T@original, atol=3e-11, rtol=0)
    # Retained unnormalized inputs are never silently normalized by the wrapper.
    scaled = r.routed_repair_branches('C', (0, 5), .25*encoded)
    for actual, expected in zip(scaled, entangled, strict=True):
        np.testing.assert_allclose(actual['amplitude'], .25*expected['amplitude'], atol=1e-14, rtol=0)


def test_controlled_flip_tensor_chart_agrees_with_explicit_source_wire_remapping():
    # A small source-executed word has enough asymmetric wire structure to
    # reject reversed chart conventions without a second large macro replay.
    word = [compiler.pulse('Y', .13, [8]), compiler.pulse('Z', -.27, [1])]
    base = compiler.execute(word, 9)
    for target in (0, 3, 8):
        wires = [bit for bit in range(9) if bit != target]+[target]
        actual = compiler.execute(compiler.remap(word, wires), 9)
        np.testing.assert_allclose(r._chart(base, target), actual, atol=1e-14, rtol=0)


def test_route_resources_charge_both_directions_and_every_cached_macro_use():
    summary = r.routing_word_summary()
    assert summary['canonical_controlled_flip_native_pulses'] == 2044
    assert summary['controlled_flip_invocations_per_direction'] == 93
    costs = summary['resources_per_direction']
    assert costs['native_pulses'] == 93*2044+summary['conjugating_source_x_pulses_per_direction']
    assert costs['proper_code_service_events'] == 2*costs['native_pulses']
    for candidate in ('C', 'E'):
        direct = r.routed_repair_resources(candidate, (0, 1))
        routed = r.routed_repair_resources(candidate, (0, 5))
        for key in r.COUNT_KEYS:
            assert routed[key] == direct[key]+2*costs[key]
        assert routed['route_directions_per_path'] == 2
        assert routed['native_input_preparations'] == 0
        assert routed['retained_classical_bits'] == 1+circuit.RECORD_BITS[candidate]
    assert r.routed_repair_resources('E', (0, 1))['native_pulses'] == (
        r.routed_repair_resources('C', (0, 1))['native_pulses']+3)


@pytest.mark.parametrize('candidate,seam,state', [
    ('bad', (0, 1), np.zeros(512)), ('C', (1, 5), np.zeros(512)),
    ('C', [0, 5], np.zeros(512)), ('C', (False, 1), np.zeros(512)),
    ('C', (0, 5), np.zeros(511)), ('C', (0, 5), np.empty((512, 0))),
    ('C', (0, 5), np.zeros((512, 1, 1))),
    ('C', (0, 5), [1., False]+[0.]*510),
    ('C', (0, 5), [np.array([1.]), np.array([False])]+[np.array([0.])]*510),
    ('C', (0, 5), np.ma.array(np.zeros(512), mask=True)),
    ('C', (0, 5), np.array([math.nan]+[0.]*511)),
    ('C', (0, 5), np.array([math.inf]+[0.]*511)),
])
def test_bad_interfaces_are_rejected_before_source_execution(candidate, seam, state, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('bad inputs must be rejected before routing execution')
    monkeypatch.setattr(r, 'routing_unitary', forbidden)
    with pytest.raises(ValueError):
        r.routed_repair_branches(candidate, seam, state)


@pytest.mark.parametrize('inverse', [1, None, np.bool_(True)])
def test_route_direction_is_an_explicit_boolean(inverse):
    with pytest.raises(ValueError):
        r.routing_unitary(inverse=inverse)
    with pytest.raises(ValueError):
        r.routing_edges(inverse=inverse)
