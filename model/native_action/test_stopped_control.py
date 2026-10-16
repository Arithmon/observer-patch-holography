"""Raw endpoint oracle, all matrix units, complete timeout and retained continuation."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import stopped_control as c
from ..observer import native


def _image(mask, edge, coin):
    a, b = edge, edge+1
    values = ((mask >> a) & 1, (mask >> b) & 1)
    target, (_, current) = native.native(values, (0, 1), coin)
    return ((mask & ~((1 << a) | (1 << b))) | target[0] << a | target[1] << b,
            current, values[0]-values[1])


def _oracle(law, left, right, *, drop_remote=False, moving_clock=False, phase=True):
    states = {((), left, right): F(1)}
    finished = defaultdict(F)
    for attempt in range(3):
        following = defaultdict(F)
        for (word, a, b), weight in states.items():
            for edge, coin in product(range(3), (0, 1)):
                x, y = _image(a, edge, coin), _image(b, edge, coin)
                if drop_remote and edge:
                    x, y = (a, 0, 0), (b, 0, 0)
                if x[1] != y[1] or (law == 'E' and x[2] != y[2]):
                    continue
                tick = edge == 0 and (not moving_clock or x[1] != 0)
                newword = (*word, x[1]) if tick else word
                sign = (-1)**((x[0] & 1)+(y[0] & 1)) if phase and tick and not word and x[1] == 0 else 1
                value = weight*sign/6
                if len(newword) == 2 or attempt == 2:
                    status = 'complete' if len(newword) == 2 else 'timeout'
                    finished[status, newword, x[0], y[0]] += value
                else:
                    following[newword, x[0], y[0]] += value
        states = following
    result = {record: {} for record in c.outcomes()}
    for (status, word, a, b), value in finished.items():
        if value:
            result[status, word][a, b] = value
    return result


def _root(column):
    result = {}
    for record, matrix in column.items():
        reduced = defaultdict(F)
        for (a, b), value in matrix.items():
            if a//2 == b//2:
                reduced[a & 1, b & 1] += value
        result[record] = {key: value for key, value in reduced.items() if value}
    return result


@pytest.mark.parametrize('law', ('C', 'E'))
def test_all_matrix_units_match_raw_native_coin_tree_with_feedback_and_timeout(law):
    for a, b in product(range(16), repeat=2):
        actual = c.exact_column(law, a, b)
        assert actual == _oracle(law, a, b)
        assert sum(value for matrix in actual.values() for (i, j), value in matrix.items()
                   if i == j) == int(a == b)


def test_complete_private_schedule_partition_and_observer_alphabet():
    assert len(c.schedules()) == 25
    assert sum(F(1, 3**len(tape)) for tape in c.schedules()) == 1
    assert len(c.outcomes()) == 13
    for cap in range(1, 5):
        success = sum(F(1, 3**len(tape)) for tape in c.schedules(cap) if tape.count(0) == 2)
        assert 1-success == F(2, 3)**cap+cap*F(1, 3)*F(2, 3)**(cap-1)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_native_remote_motion_zero_current_clock_and_feedback_are_essential(law):
    assert _root(c.exact_column(law, 8, 8)) != _root(_oracle(law, 8, 8, drop_remote=True))
    vacuum = c.exact_column(law, 0, 0)
    assert vacuum['complete', (0, 0)][0, 0] == F(7, 27)
    assert vacuum != _oracle(law, 0, 0, moving_clock=True)
    assert c.exact_column(law, 0, 15) != _oracle(law, 0, 15, phase=False)
    root_changed = _root(c.exact_column(law, 8, 9)) != _root(_oracle(law, 8, 9, phase=False))
    assert root_changed == (law == 'C')


@pytest.mark.parametrize('law', ('C', 'E'))
def test_runtime_probe_all_outcomes_exact_error_budget_and_serialized_timeout(law):
    exact, actual = c.exact_probe(law), c.execute(law)
    assert set(actual) == set(c.outcomes())
    distance = 0
    for record, matrix in actual.items():
        target = np.array([[float(exact[record].get((a, b), 0)) for b in range(4)] for a in range(4)])
        distance += np.linalg.norm(matrix-target, ord='nuc')/2
    assert distance < 1e-12
    assert c.continuation(law) < 1e-12
    timeout = sum(value for (status, _), matrix in exact.items() if status == 'timeout'
                  for (a, b), value in matrix.items() if a == b)
    assert timeout == F(20, 27)
    # Normalizing the retained success output loses this full orthogonal flag mass.
    assert sum(np.trace(matrix).real for (status, _), matrix in actual.items()
               if status == 'complete') == pytest.approx(F(7, 27))


def test_receipt_rejects_numerically_normalized_timeout_family(monkeypatch):
    original = c.execute

    def normalized(law, cap=3):
        result = original(law, cap)
        for (status, _), matrix in result.items():
            if status == 'timeout':
                matrix *= 27/20
        return result

    monkeypatch.setattr(c, 'execute', normalized)
    with pytest.raises(ArithmeticError, match='complete probe output'):
        c.run()


@pytest.mark.parametrize('value', (True, 0, 5, 3.0, F(3), np.int64(3)))
def test_original_bounded_cap_required(value):
    with pytest.raises(ValueError):
        c.exact_column('C', 0, 0, value)


@pytest.mark.parametrize('schedule', ([0, 0], (0, False), (0, 0, 1), (0, 0, 0), (3,), ()))
def test_original_stopped_schedule_required(schedule):
    with pytest.raises(ValueError):
        c.program('C', schedule)


@pytest.mark.parametrize('args', (('C', True, 0), ('E', 0, 16), ('X', 0, 0), (True, 0, 0)))
def test_original_candidate_and_matrix_unit_domain(args):
    with pytest.raises(ValueError):
        c.exact_column(*args)
