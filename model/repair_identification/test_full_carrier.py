"""Independent emitted-repair controls outside the nine-probe chart."""
import math

import numpy as np
import pytest

from model.observer import native
from model.observer.source_discrimination import circuit
from model.repair_identification import source


def _raw_records(candidate):
    for coin in (0, 1):
        for raw in range(8 if candidate == 'E' else 4):
            if candidate != 'E_collision_half' or raw == 3:
                auxiliaries = ((),)
            elif raw == 0:
                auxiliaries = ((None, 0), (None, 1))
            else:
                auxiliaries = ((0, None), (1, 0), (1, 1))
            for auxiliary in auxiliaries:
                yield coin, raw, auxiliary


def _endpoint_columns(candidate, states, codes):
    """Use endpoint floor/ceiling and recorded parity, not target tables."""
    expected = {record: np.zeros((512, 78), complex)
                for record in _raw_records(candidate)}
    metadata = {}
    for record in expected:
        coin, raw, _ = record
        current = {0: 0, 1: 1, 2: -1}.get(raw % 4)
        contrast = (2*current+(2*coin-1)*(raw//4)
                    if candidate == 'E' and current is not None else None)
        abort = current is None or (contrast is not None and abs(contrast) > 2)
        metadata[record] = ((None, None, True) if abort else (current, contrast, False))
    for column, state in enumerate(states):
        for coin in (0, 1):
            total = state[0] + state[1]
            left = (total + coin) // 2
            endpoint = (left, total-left, *state[2:])
            current = state[0]-left
            target = codes[states.index(endpoint)]
            raw = {0: 0, 1: 1, -1: 2}[current]
            contrast = state[0]-state[1] if candidate == 'E' else None
            if candidate == 'E':
                raw += 4*(abs(contrast) % 2)
            phase = (1j if candidate == 'C_phase_quarter'
                     and current == 0 and total == 1 else 1)
            paths = (((), 1/math.sqrt(2)),)
            if candidate == 'E_collision_half':
                # The actual source CNOT copies output bit selected by coin.
                parity = (target >> coin) & 1
                paths = (((None, parity), 1/math.sqrt(2)),) if current == 0 else (
                    ((0, None), .5), ((1, parity), .5))
            for auxiliary, weight in paths:
                record = coin, raw, auxiliary
                expected[record][target, column] = phase*weight
    return expected, metadata


@pytest.mark.parametrize('candidate', source.CANDIDATES)
def test_emitted_repairs_cover_all_78_columns_and_coherent_references(candidate):
    states = native.states()
    codes = tuple(circuit.encode(state) for state in states)
    assert len(set(codes)) == 78 and max(codes) >= 8
    inputs = np.zeros((512, 78), complex)
    inputs[codes, np.arange(78)] = 1
    expected, metadata = _endpoint_columns(candidate, states, codes)
    actual = source.repair_branches(candidate, inputs)
    assert tuple(branch.record for branch in actual) == tuple(expected)
    for branch in actual:
        assert (branch.current, branch.contrast, branch.abort) == metadata.get(
            branch.record, (None, None, True))
        assert branch.amplitude.shape == (512, 78)
        assert not branch.amplitude.flags.writeable
        # All columns and all 512 output rows, including retained aborts.
        assert np.linalg.norm(branch.amplitude-expected[branch.record], 'fro') < 3e-11
        assert dict(branch.resources)['native_repairs'] == 1
    gram = sum((branch.amplitude.conj().T @ branch.amplitude for branch in actual),
               np.zeros((78, 78), complex))
    assert np.linalg.norm(gram-np.eye(78), 2) < 3e-11

    # Three complex reference columns span every spectator shell. This is
    # neither a diagonal population check nor one of the prepared probe pairs.
    rng = np.random.default_rng(1117)
    coefficients = rng.normal(size=(78, 3)) + 1j*rng.normal(size=(78, 3))
    coefficients /= np.linalg.norm(coefficients)
    direct = source.repair_branches(candidate, inputs @ coefficients)
    for branch, columns in zip(direct, actual, strict=True):
        assert branch.record == columns.record
        assert np.linalg.norm(branch.amplitude-columns.amplitude @ coefficients, 'fro') < 3e-12
        assert np.linalg.norm(branch.amplitude-expected[branch.record] @ coefficients, 'fro') < 3e-11
