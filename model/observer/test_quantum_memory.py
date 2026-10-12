"""Independent floor/ceiling channels and exact memory certificate controls."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement
import json

import pytest

from . import quantum_memory as memory
from . import quantum_mixing as mixing


@lru_cache(maxsize=1)
def _oracle_events():
    states = []
    for a, b in combinations_with_replacement(range(12), 2):
        state = [0] * 12
        state[a] += 1
        state[b] += 1
        states.append(tuple(state))
    states.sort()
    assert tuple(states) == mixing.STATES
    index = {state: i for i, state in enumerate(states)}
    events = []
    for a, b in mixing.seams():
        for coin in (0, 1):
            event = []
            for state in states:
                total = state[a] + state[b]
                output = list(state)
                output[a] = total // 2 if coin == 0 else (total + 1) // 2
                output[b] = total - output[a]
                event.append((index[tuple(output)], state[a] - output[a], state[a] - state[b]))
            events.append(tuple(event))
    return tuple(events)


@lru_cache(maxsize=None)
def _oracle_column(candidate, left, right):
    record = 1 if candidate == "C" else 2
    result = Counter()
    for event in _oracle_events():
        if event[left][record] == event[right][record]:
            result[event[left][0], event[right][0]] += 1
    return tuple(result.items())


def _integer_step(candidate, counts):
    # Keep the common denominator 60**n outside the sparse state, so the
    # n=200 test uses exact integer arithmetic throughout its propagation.
    output = Counter()
    for (left, right), coefficient in counts.items():
        for pair, weight in _oracle_column(candidate, left, right):
            output[pair] += coefficient * weight
    return output


def test_rational_root_enclosure_and_sufficient_stopping_counts():
    lower, upper = memory.contraction_interval()
    assert (60 * lower - 55)**2 < 5 < (60 * upper - 55)**2
    assert 0 < lower < upper < 1
    threshold = F(1, 4355)
    for candidate in ("C", "E"):
        assert memory.choi_purity_interval(candidate, 0) == (F(1), F(1))
        attempts = memory.stopping_bound(candidate)
        assert attempts == 134
        assert memory.choi_purity_interval(candidate, attempts - 1)[0] > threshold
        assert memory.choi_purity_interval(candidate, attempts)[1] < threshold
        assert memory.choi_purity_interval(candidate, attempts + 1)[1] < threshold
    assert lower**177 > threshold > upper**178


@pytest.mark.parametrize("candidate", ["C", "E"])
def test_one_step_normalized_choi_purity_from_independent_matrix_units(candidate):
    binary = [i for i, state in enumerate(mixing.STATES) if max(state) == 1]
    numerator = 0
    for left in binary:
        for right in binary:
            column = _oracle_column(candidate, left, right)
            numerator += sum(count**2 for _, count in column)
    # Choi reshuffling preserves Frobenius norm and the normalized
    # maximally-entangled input supplies exactly the factor 1/D**2.
    purity = F(numerator, 60**2 * len(binary)**2)
    assert F(1, 66**2) < purity < memory.choi_purity_interval(candidate, 1)[1] < 1
    assert purity > F(1, 4355)  # The ball criterion does not yet certify EB.


@pytest.mark.parametrize("candidate, rate", [("C", F(23, 30)), ("E", F(7, 10))])
def test_full_output_partial_transpose_minor_by_200_independent_steps(candidate, rate):
    u, v = mixing.seams()[0]
    pile, binary = [0] * 12, [0] * 12
    pile[u], binary[u], binary[v] = 2, 1, 1
    i, j = mixing.INDEX[tuple(pile)], mixing.INDEX[tuple(binary)]
    pile_population = Counter({(i, i): 1})
    binary_population = Counter({(j, j): 1})
    coherence = Counter({(i, j): 1})
    for attempts in range(201):
        denominator = 60**attempts
        # These three propagated operators, plus the adjoint of coherence,
        # are the complete output of the entangled two-term preparation.
        assert sum(pile_population.values()) == denominator
        assert sum(binary_population.values()) == denominator
        assert all(max(mixing.STATES[a]) == 1 for a, _ in binary_population)
        if attempts in (0, 1, 10, 134, 200):
            top = F(binary_population.get((i, i), 0), 2 * denominator)
            bottom = F(pile_population.get((j, j), 0), 2 * denominator)
            offdiagonal = F(coherence.get((i, j), 0), 2 * denominator)
            assert top == 0
            assert bottom >= 0
            assert offdiagonal == rate**attempts / 2
            determinant = top * bottom - offdiagonal**2
            assert determinant == memory.full_carrier_npt_minor(candidate, attempts) < 0
        if attempts < 200:
            pile_population = _integer_step(candidate, pile_population)
            binary_population = _integer_step(candidate, binary_population)
            coherence = _integer_step(candidate, coherence)


def test_phase_pair_complete_density_outputs_and_readout_are_exact():
    u, v = mixing.seams()[0]
    pairs = []
    initial = Counter()
    for spectator in range(12):
        if spectator in (u, v):
            continue
        left, right = [0] * 12, [0] * 12
        left[u] = left[spectator] = right[v] = right[spectator] = 1
        i, j = mixing.INDEX[tuple(left)], mixing.INDEX[tuple(right)]
        pairs.append((i, j))
        initial.update(((i, i), (j, j), (i, j), (j, i)))
    outputs = {candidate: initial.copy() for candidate in ("C", "E")}
    for attempts in range(7):
        denominator = 20 * 60**attempts
        expected = {}
        coefficient = ((F(5, 6)**attempts - F(7, 10)**attempts) / 20)
        for i, j in pairs:
            expected[i, j] = expected[j, i] = coefficient
        # Compare every matrix entry of the complete output difference;
        # matching a single aggregate contrast would be weaker evidence.
        for pair in set(outputs["C"]) | set(outputs["E"]):
            difference = F(outputs["C"].get(pair, 0) - outputs["E"].get(pair, 0), denominator)
            assert difference == expected.get(pair, F(0))
        for candidate, rate in (("C", F(5, 6)), ("E", F(7, 10))):
            state = outputs[candidate]
            assert sum(count for (i, j), count in state.items() if i == j) == denominator
            expectation = F(sum(state.get((i, j), 0) + state.get((j, i), 0) for i, j in pairs), denominator)
            assert expectation == rate**attempts
        assert sum(expected.values()) / 2 == memory.phase_pair_signal(attempts)
        outputs = {candidate: _integer_step(candidate, state) for candidate, state in outputs.items()}
    assert memory.phase_pair_signal(4) == F(12259, 101250)
    assert memory.phase_pair_signal(3) < memory.phase_pair_signal(4) > memory.phase_pair_signal(5)


def test_receipt_has_checked_rational_evidence_and_explicit_scope():
    result = memory.audit()
    assert result["binary_dimension"] == 66
    assert result["traceless_diagonal_dimensions"] + result["offdiagonal_dimensions"] == 4355
    assert result["normalized_choi_purity_threshold"] == "1/4355"
    assert result["crude_sufficient_attempts"] == 178
    for report in result["binary"].values():
        assert F(report["stopping_spectral_bound_scaled_excess_upper"]) < 1
        assert F(report["previous_attempt_spectral_bound_scaled_excess_lower"]) > 1
        assert report["exact_entanglement_breaking_index_claimed"] is False
    assert result["retained_records_or_private_memory_covered"] is False
    assert result["phase_pair"]["new_source_execution"] is False
    def no_float(value):
        assert type(value) is not float
        if isinstance(value, dict):
            for child in value.values():
                no_float(child)
        elif isinstance(value, (list, tuple)):
            for child in value:
                no_float(child)
    no_float(result)
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("attempts", [-1, True, False, 1.0, "1", None, F(1)])
def test_attempt_parameters_fail_closed(attempts):
    for operation in (lambda: memory.choi_purity_interval("C", attempts),
                      lambda: memory.phase_pair_signal(attempts),
                      lambda: memory.full_carrier_npt_minor("C", attempts)):
        with pytest.raises(ValueError):
            operation()


@pytest.mark.parametrize("candidate", ["c", "unknown", 0, True, None, []])
def test_candidate_parameters_fail_closed(candidate):
    for operation in (lambda: memory.choi_purity_interval(candidate, 1),
                      lambda: memory.stopping_bound(candidate),
                      lambda: memory.full_carrier_witness(candidate)):
        with pytest.raises(ValueError):
            operation()


def test_corrupt_native_witness_row_is_rejected(monkeypatch):
    original = mixing.native_block
    def corrupt(candidate, difference, *, binary=False):
        rows = [list(row) for row in original(candidate, difference, binary=binary)]
        if difference != mixing.ZERO:
            for i in range(len(rows)):
                rows[i][i] = F(0)
        return tuple(map(tuple, rows))
    monkeypatch.setattr(mixing, "native_block", corrupt)
    with pytest.raises(ArithmeticError, match="isolated and positive"):
        memory.full_carrier_witness("C")
