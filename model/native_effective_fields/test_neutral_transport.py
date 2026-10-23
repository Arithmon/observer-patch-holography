"""Independent native blocks, complete instruments and clock countercontrols."""
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest
import sympy as sp

from . import neutral_transport as n


def test_all_native_phase_columns_attach_to_punctured_graph():
    audit = n.exact_audit()
    assert audit['all_endpoint_pairs'] == 66
    assert audit['source_columns_checked_across_both_laws'] == 1320
    assert audit['every_coherent_movement_uses_one_punctured_source_edge']


@pytest.mark.parametrize('pair', [(0, 1), (0, 2), (0, 3)])
@pytest.mark.parametrize('law', ['C', 'E'])
def test_nonflat_attempt_bands_and_real_stopped_bands(pair, law):
    a = n.adjacency(pair)
    degrees = [sum(a[:, j]) for j in range(10)]
    laplacian = sp.diag(*degrees)-a
    transfer = n.response(pair, law)
    r = 1-sp.Rational(n.killing(pair, law), 30)
    assert transfer*sp.ones(10, 1) == r*sp.ones(10, 1)
    assert transfer != r*sp.eye(10)
    assert (r*sp.eye(10)-transfer)*60 == laplacian
    assert laplacian.rank() == 9
    stopped = n.response(pair, law, clock='movement')
    d = [2*n.killing(pair, law)+x for x in degrees]
    # R=A D^-1 is similar to the real symmetric D^-1/2 A D^-1/2.
    assert stopped*sp.diag(*d) == a
    assert max(sum(stopped[:, j]) for j in range(10)) <= sp.Rational(5, 2*n.killing(pair, law)+5)
    i, j = next((i, j) for i, j in combinations(range(10), 2) if a[i, j])
    assert stopped.extract((i, j), (i, j)).det() == -1/(d[i]*d[j])
    # A negative principal minor of the symmetric similar block excludes
    # a positive transfer and hence a self-adjoint -log(R) spectral claim.
    assert stopped[i, i] == stopped[j, j] == 0


def test_first_movement_uses_actual_edges_and_complete_phase_probabilities():
    pair, common = (0, 1), 2
    col = n.common_ports(pair).index(common)
    for law, divisor in [('C', 15), ('E', 23)]:
        marked = n.stopped.source(law)
        probs = n.exact_probabilities(n.exact_first_current(marked, n.preparation(pair, common)), pair)
        assert sum(probs) == 1 and min(probs) >= 0
        contrasts = tuple(probs[2*k]-probs[2*k+1] for k in range(10))
        assert contrasts == tuple(F(n.adjacency(pair)[k, col], divisor) for k in range(10))
        assert sum(contrasts) == F(5, divisor)
        # Lost contrast remains in legitimate phase/complement outcomes;
        # it is not a herald failure or a probability to renormalize away.
        assert sum(contrasts) < sum(probs)
        negative = n.exact_probabilities(n.exact_first_current(marked, n.preparation(pair, common, -1)), pair)
        assert negative[-1] == probs[-1]
        for k in range(10):
            assert negative[2*k:2*k+2] == tuple(reversed(probs[2*k:2*k+2]))


def test_joint_analyzer_is_a_complete_positive_projective_measurement():
    read = n.readout((0, 1))
    projectors = tuple(group[0] for group in read.kraus)
    assert len(projectors) == 21
    assert np.array_equal(sum(projectors), np.eye(66))
    assert np.trace(projectors[-1]) == 46
    for i, p in enumerate(projectors):
        assert np.array_equal(p, p.conj().T) and np.array_equal(p@p, p)
        for q in projectors[:i]:
            assert not np.any(p@q)


def test_spatial_analyzer_attains_exact_C_E_trace_distance():
    result = n.discrimination({law: n.stopped.source(law) for law in ('C', 'E')})
    assert result['ideal_first_tick_output_trace_distance'] == '4/69'
    assert result['complete_measurement_total_variation'] == '4/69'


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('sign', [1, -1])
def test_executed_complete_capped_history_retains_timeouts(law, sign):
    row = n.execution(n.stopped.source(law), sign=sign)
    assert row['complete_records'] == 168
    assert row['timeout_mass'] == '125/216'
    assert sum(map(F, row['first_movement_joint_outcomes']))+sum(map(F, row['timeout_joint_outcomes'])) == 1
    assert row['terminal_state_and_probability_error_diagnostic'] < 2e-11
    assert row['checkpoint_complete_distance_diagnostic'] < 2e-11


def test_successive_ticks_are_compositions_not_a_renormalized_walk():
    marked, pair, common = n.stopped.source('C'), (0, 1), 2
    state = n.preparation(pair, common)
    column = sp.eye(10)[:, n.common_ports(pair).index(common)]
    r = n.response(pair, 'C', clock='movement')
    for _ in range(3):
        state = n.exact_first_current(marked, state)
        column = r*column
        probs = n.exact_probabilities(state, pair)
        assert sum(probs) == 1
        assert tuple(probs[2*k]-probs[2*k+1] for k in range(10)) == tuple(map(F, column))


def test_identity_idling_cancels_but_native_zero_current_does_not():
    for law in ('C', 'E'):
        t = n.native_response((0, 1), law)
        zero = n.native_response((0, 1), law, 0)
        moving = t-zero
        q = sp.Rational(2, 7)
        assert q*moving*(sp.eye(10)-((1-q)*sp.eye(10)+q*zero)).inv() == n.response((0, 1), law, clock='movement')
        assert zero != sp.eye(10)


@pytest.mark.parametrize('pair', [(True, 1), (1, 1), (1, 0), (0, 12), [0, 1]])
def test_bad_endpoint_domains(pair):
    with pytest.raises(ValueError):
        n.response(pair, 'C')


def test_invalid_clock_preparation_and_cap():
    with pytest.raises(ValueError):
        n.response((0, 1), 'C', clock='seconds')
    with pytest.raises(ValueError):
        n.preparation((0, 1), 0)
    with pytest.raises(ValueError):
        n.preparation((0, 1), 2, True)
    with pytest.raises(ValueError):
        n.program(n.stopped.source('C'), (0, 1), 2, cap=4)
