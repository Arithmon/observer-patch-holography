from fractions import Fraction as F
from itertools import combinations
import math

import numpy as np
import pytest

from .construction import (exact_native_channel, graph_certificate,
                           preparation_budget, primitive_execution)


@pytest.mark.parametrize('edges', [((0, 1), (1, 2), (2, 3)),
    ((0, 1), (0, 2), (0, 3)), tuple(combinations(range(4), 2))])
@pytest.mark.parametrize('particles', [1, 2, 3])
@pytest.mark.parametrize('law', ['C', 'E'])
def test_native_all_operator_gap_and_fixed_state(edges, particles, law):
    counts, denominator, basis = exact_native_channel(4, edges, particles, law)
    channel = counts/denominator
    d = len(basis)
    identity = np.eye(d).reshape(-1)
    np.testing.assert_allclose(channel@identity, identity, atol=1e-14)
    np.testing.assert_array_equal(counts, counts.T)
    values = np.linalg.eigvalsh(channel)
    assert values[0] > -1e-14
    assert np.sum(abs(values-1) < 1e-12) == 1
    laplacian = np.zeros((4, 4))
    for a, b in edges:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    expected = np.linalg.eigvalsh(laplacian)[1]/denominator
    assert abs(1-values[-2]-expected) < 2e-14
    lower = F(graph_certificate(4, edges)['quantum_gap_lower_per_attempt'])
    assert expected >= float(lower)-2e-14


@pytest.mark.parametrize('law', ['C', 'E'])
def test_mixing_bound_includes_untouched_entangled_reference(law):
    counts, scale, basis = exact_native_channel(3, ((0, 1), (1, 2)), 1, law)
    d = len(basis)
    initial = np.zeros((d, d, d, d))  # system ket, reference ket, system bra, reference bra
    for i in range(d):
        for j in range(d):
            initial[i, i, j, j] = 1/d
    phi = (counts/scale).reshape(d, d, d, d)
    rho = initial
    target = np.eye(d*d)/(d*d)
    gamma = F(1, 4)
    for n in range(1, 16):
        rho = np.einsum('abij,irjs->arbs', phi, rho)
        dist = np.abs(np.linalg.eigvalsh(rho.reshape(d*d, d*d)-target)).sum()/2
        assert dist <= min(1, d/2*float((1-gamma)**n))+2e-14


@pytest.mark.parametrize('dimension', [2, 6, 66, math.comb(3840, 1920)])
@pytest.mark.parametrize('gap', [F(1), F(1, 2), F(11, 240), F(1, 123456789)])
def test_budget_has_integer_certificate(dimension, gap):
    b = preparation_budget(dimension, gap, 10)
    assert b['half_contraction_block']*gap >= 1
    assert F(b['trace_distance_upper']) <= F(1, 1024)
    assert b['attempts'] == b['blocks']*b['half_contraction_block']
    if b['half_contraction_block'] < 1000:
        assert (1-gap)**b['half_contraction_block'] <= F(1, 2)


def test_sector_input_cannot_be_selected_by_mixing():
    for particles in (0, 4):
        counts, scale, basis = exact_native_channel(4, ((0, 1), (1, 2), (2, 3)), particles, 'C')
        assert counts.tolist() == [[scale]]
        assert len(basis) == 1
        assert sum(basis[0]) == particles
    assert preparation_budget(1, F(1, 10), 10)['attempts'] == 0


@pytest.mark.parametrize('args', [(4, ((0, 1), (2, 3))),
    (3, ((0, 1), (0, 1))), (3, ((1, 0),)), (True, ((0, 1),)),
    (3, ((0, 3),)), (3, ((0, True),))])
def test_invalid_or_disconnected_architecture_rejected(args):
    with pytest.raises(ValueError):
        graph_certificate(*args)


@pytest.mark.parametrize('args', [(0, F(1, 2), 10), (2, .5, 10),
    (2, F(0), 10), (2, F(2), 10), (2, F(1, 2), True)])
def test_uncertified_budget_inputs_rejected(args):
    with pytest.raises(ValueError):
        preparation_budget(*args)


def test_primitive_preparation_and_future_complete_law():
    result = primitive_execution()
    for row in result['executions']:
        assert row['trace_distance_diagnostic'] <= float(F(result['budget']['trace_distance_upper']))
        assert abs(sum(r['prepared_probability'] for r in row['two_current_bit_law'])-1) < 2e-12
        assert row['future_law_TV_diagnostic'] <= row['trace_distance_diagnostic']+2e-12
