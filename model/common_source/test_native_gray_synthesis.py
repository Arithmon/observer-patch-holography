"""Exact parity accounting and literal source checks for the shorter macro.

The symbolic oracle reads the emitted word, without reproducing its Gray
construction. Source replay separately checks every complex matrix column
against both the original Walsh compiler and a basis-swap target.
"""
import math

import numpy as np
import pytest

from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from . import native_routing as routing


def test_emitted_phase_word_is_exactly_one_full_register_reflection():
    word = routing._canonical_word()
    assert word[0] == compiler.pulse('Y', math.pi/4, [8])
    assert word[-1] == compiler.pulse('Y', -math.pi/4, [8])
    addresses = np.arange(512)
    phase = np.zeros(512, dtype=np.int64)
    index, phases, cnots = 1, 0, 0
    while index < len(word)-1:
        event = word[index]
        if event['kind'] in ('I', 'Z'):
            coefficient = event['angle']/(-math.pi/512)
            assert coefficient in (-1, 1)
            if event['kind'] == 'I':
                phase += int(coefficient)
            else:
                bit = event['targets'][0]
                phase += int(coefficient)*(1-2*((addresses >> bit) & 1))
            index += 1
            phases += 1
        else:
            # Read the actual phase-exact CNOT's wire pair, then evolve its
            # exact integer basis action. No Gray ordering is assumed here.
            assert len(word[index:index+3]) == 3
            middle = word[index+1]
            assert middle['kind'] == 'J' and middle['angle'] == -math.pi
            control, target = middle['targets']
            assert word[index] == word[index+2] == compiler.pulse('H', math.pi/2, [target])
            addresses ^= ((addresses >> control) & 1) << target
            index += 3
            cnots += 1
    np.testing.assert_array_equal(addresses, np.arange(512))
    expected = np.zeros(512, dtype=np.int64)
    expected[256] = 512
    np.testing.assert_array_equal(phase, expected)
    assert phases == 512 and cnots == 510
    assert len(word) == 2044
    assert routing._canonical_word(True) == compiler.inverse(word)


@pytest.fixture(scope='module')
def literal_words():
    old_word = circuit.controlled_flip_word(9, 8, 0)
    new_word = routing._canonical_word()
    return (compiler.execute(old_word, 9), compiler.execute(new_word, 9),
            compiler.execute(routing._canonical_word(True), 9))


def test_literal_source_replays_preserve_full_phase_padding_and_inverse(literal_words):
    old, new, inverse = literal_words
    expected = np.eye(512, dtype=complex)
    expected[:, [0, 256]] = expected[:, [256, 0]]
    for matrix in (old, new, inverse):
        np.testing.assert_allclose(matrix, expected, atol=2e-12, rtol=0)
    np.testing.assert_allclose(new, old, atol=2e-12, rtol=0)
    np.testing.assert_allclose(inverse@new, np.eye(512), atol=2e-12, rtol=0)
    # The off-code rows and their coherences are included in these arbitrary,
    # unnormalized reference columns; a code compression cannot pass this.
    rng = np.random.default_rng(1102009)
    reference = rng.normal(size=(512, 3))+1j*rng.normal(size=(512, 3))
    np.testing.assert_allclose(new@reference, expected@reference, atol=5e-12, rtol=0)


def test_shorter_word_charges_every_source_pulse_without_new_services():
    old = compiler.resources(circuit.controlled_flip_word(9, 8, 0))
    new = compiler.resources(routing._canonical_word())
    assert old['native_pulses'] == 11782
    assert new['native_pulses'] == 2044
    assert new['proper_code_service_events'] == 2*2044
    assert {bit for event in routing._canonical_word() for bit in event['targets']} == set(range(9))
    assert new == compiler.resources(routing._canonical_word(True))
    for seam in routing.SEAMS:
        summary = routing.routing_word_summary(seam)
        calls = summary['controlled_flip_invocations_per_direction']
        extra = summary['conjugating_source_x_pulses_per_direction']
        counts = summary['resources_per_direction']
        assert counts['native_pulses'] == calls*2044+extra
        assert counts['proper_code_service_events'] == 2*counts['native_pulses']
