"""Source-pulse replay and independent full-carrier circuit tests."""
import math

import numpy as np
import pytest

from model.observer import native
from model.observer.source_scalar import compiler as source
from model.predictions import observer_instrument as prediction
from model.observer.source_discrimination import circuit


def test_embedding_is_exactly_the_native_carrier_and_fixes_preparation_codes():
    assert len(circuit.SOURCE_CODES) == len(set(circuit.SOURCE_CODES)) == 78
    assert tuple(map(len, circuit.SPECTATOR_SHELLS)) == (1, 10, 55)
    assert circuit.encode(prediction.STATES[prediction.X]) == 0
    assert circuit.encode(prediction.STATES[prediction.Z]) == 1
    admitted = []
    for code in range(512):
        if code in circuit.SOURCE_CODES:
            state = circuit.decode(code)
            assert circuit.encode(state) == code
            admitted.append(state)
        else:
            with pytest.raises(ValueError):
                circuit.decode(code)
    assert set(admitted) == set(native.states())


@pytest.mark.parametrize("bits,target,pattern", ((1, 0, 0), (2, 0, 2),
                                                (3, 1, 5), (4, 2, 9)))
def test_controlled_flip_retains_full_phase_on_every_column(bits, target, pattern):
    actual = source.execute(circuit.controlled_flip_word(bits, target, pattern), bits)
    expected = np.eye(1 << bits, dtype=complex)
    a, b = pattern, pattern | (1 << target)
    expected[[a, b]] = expected[[b, a]]
    assert np.linalg.norm(actual - expected, 2) < 1e-12


@pytest.mark.parametrize("candidate,coin", (("C", 0), ("C", 1), ("E", 0), ("E", 1)))
def test_every_local_dilation_column_and_all_78_native_inputs(candidate, coin):
    bits = 3 + circuit.RECORD_BITS[candidate]
    word = circuit.repair_word(candidate, coin)
    actual = source.execute(word, bits)
    current_permutation = list(range(32))
    for a, b in circuit.REPAIR_SWAPS["C", coin]:
        current_permutation[a], current_permutation[b] = b, a
    expected = np.zeros((1 << bits, 1 << bits), complex)
    for column in range(1 << bits):
        output = (column // 32) * 32 + current_permutation[column % 32]
        if candidate == "E" and output >> coin & 1:
            output ^= 32
        expected[output, column] = 1
    # This absolute, phase-sensitive comparison also checks every padding
    # and initially nonblank-record column of the declared unitary.
    assert np.linalg.norm(actual - expected, 2) < 3e-12
    assert len(word) == (5370 if candidate == "C" else 5373)
    record_codes = {value: code for code, value in circuit.record_values(candidate, coin).items()}
    for state in native.states():
        # Independent floor/ceiling characterization, not native.native.
        total = state[0] + state[1]
        low, high = total // 2, (total + 1) // 2
        endpoint = (low, high) if coin == 0 else (high, low)
        repaired = endpoint + state[2:]
        value = state[0] - endpoint[0] if candidate == "C" else state[0] - state[1]
        target = circuit.SEAM_STATES.index(endpoint) + 8 * record_codes[value]
        source_code = circuit.encode(state)
        column = actual[:, source_code % 8]
        assert np.linalg.norm(column - np.eye(1 << bits)[:, target]) < 3e-12
        # The untouched six spectator bits embed the same local identity
        # into the complete carrier, including cross-shell coherences.
        assert circuit.decode((target % 8) + 8 * (source_code // 8)) == repaired


@pytest.mark.parametrize("coin", (0, 1))
def test_contrast_is_exactly_current_plus_one_recorded_parity_bit(coin):
    current_word = circuit.repair_word("C", coin)
    contrast_word = circuit.repair_word("E", coin)
    assert contrast_word == current_word + source.cnot(coin, 5)
    records = circuit.record_values("E", coin)
    assert set(records.values()) == {-2, -1, 0, 1, 2}
    assert len(records) == 5
    current_values = circuit.record_values("C", coin)
    for raw in range(8):
        current = current_values.get(raw % 4)
        if current is None:
            assert raw not in records
            continue
        contrast = 2 * current + (2 * coin - 1) * (raw // 4)
        if abs(contrast) > 2:
            assert raw not in records
        else:
            assert records[raw] == contrast
    for state in native.states():
        total = state[0] + state[1]
        low, high = total // 2, (total + 1) // 2
        endpoint = (low, high) if coin == 0 else (high, low)
        current = state[0] - endpoint[0]
        parity = total % 2
        # This identity makes a single source-qubit CNOT sufficient; no
        # additional many-qubit parity computation is hidden in the word.
        assert (circuit.SEAM_STATES.index(endpoint) >> coin) & 1 == parity
        raw = circuit.CURRENT_CODES[current] + 4 * parity
        assert records[raw] == state[0] - state[1]


@pytest.mark.parametrize("setting", circuit.SETTINGS[:3])
def test_source_preparation_from_blank_has_the_fixed_complex_phase(setting):
    word = circuit.preparation_word(setting)
    actual = source.execute(word, 1, np.array([1., 0.]))
    theta = {"phase_0": 0., "phase_pi": math.pi,
             "phase_pi_over_2": math.pi / 2}[setting]
    expected = np.array([1., np.exp(1j * theta)]) / math.sqrt(2)
    assert np.linalg.norm(actual - expected) < 2e-14
    assert len(word) == 3
    if setting == "phase_0":
        assert sum(event["angle"] == 0 for event in word) == 2


def test_readout_executes_a_lueders_instrument_with_an_untouched_reference():
    rng = np.random.default_rng(7019)
    # The two columns are an arbitrary entangled reference. They include
    # every source code and padding, testing the complement's coherences.
    initial = rng.normal(size=(512, 2)) + 1j * rng.normal(size=(512, 2))
    initial /= np.linalg.norm(initial)
    blank_pointer = np.zeros((1024, 2), complex)
    blank_pointer[:512] = initial
    word = circuit.readout_word()
    actual = source.execute(word, 10, blank_pointer)
    plus = np.zeros_like(initial)
    plus[0] = plus[1] = (initial[0] + initial[1]) / 2
    expected = np.vstack((initial - plus, plus))
    assert len(word) == 26632
    assert np.linalg.norm(actual - expected) < 2e-11
    # The minus branch is a single coherent Q|psi>, not a mixture of
    # separately measured source basis states.
    assert np.linalg.norm(actual[:512] - (initial - plus)) < 2e-11


@pytest.mark.parametrize("function,args", (
    (circuit.encode, ((True,) + (0,) * 11,)),
    (circuit.decode, (True,)),
    (circuit.decode, (512,)),
    (circuit.controlled_flip_word, (3, 1, 2)),
    (circuit.controlled_flip_word, (11, 1, 0)),
    (circuit.basis_swap_word, (2, 2, 3)),
    (circuit.repair_word, ("C", True)),
    (circuit.record_values, ("E", True)),
    (circuit.record_values, ("unknown", 0)),
    (circuit.preparation_word, ("unknown",)),
    (circuit.preparation_word, ("diagonal_control",)),
))
def test_invalid_circuit_requests_fail(function, args):
    with pytest.raises(ValueError):
        function(*args)
