"""Independent exact controls for native coarse history prediction."""
from fractions import Fraction as F
from itertools import combinations
import json

import pytest

from . import memory as m
from .source import binary_step


def pairing(a, b, size):
    return sum((value*b.get(key, F(0)) for key, value in a.items()), F(0))/size


def exact_native_mean(graph, preparation):
    """Enumerate actual seam/coin branches of a finite configuration mixture."""
    result = {}
    for state, weight in preparation:
        for seam in graph.seams:
            for coin in (0, 1):
                final, _ = binary_step(state, seam, coin)
                for vertex in final:
                    result[vertex] = result.get(vertex, F(0)) + weight/F(2*len(graph.seams))
    return result


def independent_edge_step(graph, values):
    """Apply edge differences, independently of production neighbour rows."""
    result = dict(values)
    denominator = 2*len(graph.seams)
    for a, b in graph.seams:
        flux = (values.get(b, F(0))-values.get(a, F(0)))/denominator
        result[a] = result.get(a, F(0)) + flux
        result[b] = result.get(b, F(0)) - flux
    return {key: value for key, value in result.items() if value}


def test_source_attempt_profile_matches_complete_native_fixed_two_mixture():
    reduction = m.Reduction()
    fibre = reduction.fibres[2]
    mixture = [(state, F(1, 6)) for state in combinations(fibre, 2)]
    expected = exact_native_mean(reduction.graph, mixture)
    initial = {vertex: F(1, 2) for vertex in fibre}
    assert reduction.step(initial) == expected
    assert sum(expected.values()) == 2


def test_exact_block_recurrence_and_memory_convolution_for_arbitrary_signed_profile():
    reduction = m.Reduction()
    initial = {0: F(2, 3), 1: F(-7, 11), 13: F(3, 5), 82: F(7, 4)}
    result = reduction.predict(initial, 5)
    fine = initial
    for n, profile in enumerate(result['profiles']):
        assert profile == reduction.restrict(fine)
        fine = independent_edge_step(reduction.graph, fine)
        if n == 5:
            continue
        # Independently form every D power from the recorded coarse history.
        force = reduction.hidden(initial)
        for _ in range(n):
            force = reduction.d(force)
        terms = [reduction.a(profile), reduction.b(force)]
        for j in range(n):
            history = reduction.c(result['profiles'][j])
            for _ in range(n-1-j):
                history = reduction.d(history)
            terms.append(reduction.b(history))
        assert m._add(*terms) == result['profiles'][n+1]


def test_zero_hidden_initial_profile_still_requires_memory_after_first_attempt():
    reduction = m.Reduction()
    initial = reduction.lift({2: F(1, 2)})
    result = reduction.predict(initial, 3)
    assert reduction.hidden(initial) == {}
    assert all(row['forcing'] == {} for row in result['terms'])
    assert result['terms'][0]['memory'] == {}
    assert result['terms'][1]['memory'] == reduction.b(reduction.c({2: F(1, 2)}))
    assert result['terms'][1]['memory']
    assert result['profiles'][1] == reduction.a(result['profiles'][0])
    assert result['profiles'][2] != reduction.a(reduction.a(result['profiles'][0]))


def test_same_coarse_preparation_requires_initial_hidden_forcing():
    reduction = m.Reduction()
    first, second = {0: F(1), 1: F(1)}, {0: F(1), 13: F(1)}
    assert reduction.restrict(first) == reduction.restrict(second)
    a, b = reduction.predict(first, 1), reduction.predict(second, 1)
    assert a['profiles'][1] != b['profiles'][1]
    assert a['terms'][0]['markov'] == b['terms'][0]['markov']
    difference = m._add(a['profiles'][1], m._negative(b['profiles'][1]))
    assert difference[2] == F(1, 20160)
    assert reduction.fibre_size*difference[2] == F(1, 5040)


def test_uniform_site_adjoint_and_positive_memory_kernel():
    reduction = m.Reduction()
    coarse = {0: F(3, 4), 2: F(-2), 82: F(1, 5)}
    fine = {0: F(2), 1: F(-3), 13: F(7), 82: F(4, 7)}
    hidden = reduction.hidden(fine)
    assert reduction.restrict(hidden) == {}
    assert pairing(reduction.lift(coarse), fine, reduction.graph.ports) == pairing(
        coarse, reduction.restrict(fine), reduction.coarse_size)
    assert pairing(coarse, reduction.b(hidden), reduction.coarse_size) == pairing(
        reduction.c(coarse), hidden, reduction.graph.ports)
    assert m._rms2(reduction.d(hidden), reduction.graph.ports) <= m._rms2(hidden, reduction.graph.ports)
    kernel = reduction.c(coarse)
    for _ in range(4):
        assert pairing(coarse, reduction.b(kernel), reduction.coarse_size) >= 0
        kernel = reduction.d(kernel)


def test_equal_level_restriction_has_zero_memory_for_all_inputs():
    reduction = m.Reduction(0, 0)
    initial = {1: F(1), 27: F(2, 3)}
    result = reduction.predict(initial, 4)
    assert reduction.hidden(initial) == {}
    assert all(row['memory'] == row['forcing'] == {} for row in result['terms'])


@pytest.mark.parametrize('values', ({True: 1}, {0: True}, {0: .5}, {-1: 1}, {960: 1}, [(0, 1)]))
def test_exact_input_boundary(values):
    with pytest.raises(ValueError):
        m.Reduction().step(values)


@pytest.mark.parametrize('steps', (True, -1, 17, 1.0))
def test_horizon_boundary(steps):
    with pytest.raises(ValueError):
        m.Reduction().predict({}, steps)


def test_hidden_domain_and_physical_preparation_guards():
    reduction = m.Reduction()
    for operation in (reduction.b, reduction.d):
        with pytest.raises(ValueError):
            operation({0: 1})
    with pytest.raises(ValueError):
        m.audit_preparation(reduction, {0: F(3, 2)}, 1)
    assert reduction.predict({}, 0) == {'profiles': [{}], 'terms': []}


def test_complete_retained_receipt_is_exact_and_json_portable():
    result = m.run()
    assert json.loads(json.dumps(result)) == result
    block = result['block_constant_fixed_two_preparation']['certificate']
    assert block['initial_hidden_rms_squared']['numerator'] == 0
    assert len(block['finite_horizons']) == 9
    for row in block['finite_horizons']:
        assert row['complete_coarse_profile_residual']['numerator'] == 0
        assert row['total_expected_particles'] == {'numerator': 2, 'denominator': 1}
    assert result['scope']['full_joint_law_or_quantum_instrument_recovered'] is False
