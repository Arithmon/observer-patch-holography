"""Actual screen neighborhoods, native sign transport and state domains."""
from collections import defaultdict
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
from pathlib import Path

import pytest

from tools import strict_json
from ..observer.scalar_response.experiment import compare_receipt
from . import construction as local
from . import screen_compatibility as c


@pytest.mark.parametrize('first,second,radius,size', ((1, 3, 3, 45), (2, 4, 6, 114), (3, 5, 9, 222)))
def test_actual_pinned_source_ball_edges_and_addresses(first, second, radius, size):
    a, b = c.root_ball(first, radius), c.root_ball(second, radius)
    mapping = c.identify(a, b)
    assert len(mapping) == size
    assert all(p % 12 == q % 12 for p, q in mapping.items())
    assert len({p//12 for p in mapping}) > 1
    with pytest.raises(ValueError):
        c.identify(c.root_ball(first, radius+1), c.root_ball(second, radius+1))
    rows = list(b.adjacency)
    port, neighbours = rows[0]
    rows[0] = port, neighbours[1:]
    with pytest.raises(ValueError):
        c.identify(a, replace(b, adjacency=tuple(rows)))
    # Identical forged adjacency or truncated address domains cannot certify
    # a radius larger than the actual complete source neighborhood.
    with pytest.raises(ValueError):
        c.identify(replace(a, addresses=a.addresses[:-1]), replace(b, addresses=b.addresses[:-1]))
    with pytest.raises(ValueError):
        c.identify(replace(a, radius=radius+1), replace(b, radius=radius+1))


def test_adjoining_levels_do_not_silently_reuse_fixed_port_frames():
    with pytest.raises(ValueError):
        c.identify(c.root_ball(3, 1), c.root_ball(4, 1))


def _group_column(groups, left, right):
    output = defaultdict(F)
    for weight, matrix in groups:
        a = [(i, row[left]) for i, row in enumerate(matrix) if row[left]]
        b = [(j, row[right]) for j, row in enumerate(matrix) if row[right]]
        for (i, x), (j, y) in product(a, b):
            output[i, j] += weight*x*y
    return dict(output)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_all_source_outcome_matrix_units_with_signed_orientation_transport(law):
    # The actual mapped source paths have reversed global ordering. Forgetting
    # the current-sign transport would therefore fail this all-input check.
    for path in c._probe_paths():
        assert any(a > b for a, b in zip(path, path[1:]))
        for edge in range(3):
            groups = c._groups(law, path, edge)
            for left, right, current in product(range(16), range(16), (-1, 0, 1)):
                assert _group_column(groups[current+1], left, right) == local.quantum.native_unit(
                    law, left, right, (edge, edge+1), current)


def test_coherent_reference_feedback_complete_law_and_serialized_continuation():
    for row in c.probe_execution():
        assert row['complete_current_words'] == 81
        assert row['untouched_reference_dimension'] == 2
        assert row['complete_record_and_terminal_state_distance_diagnostic'] < 1e-12
        assert row['checkpoint_continuation_residual_diagnostic'] < 1e-12


def test_state_compatibility_does_not_replace_fixed_two_by_a_finite_density():
    result = c.initial_state_certificate()
    assert result['endpoint_stationarity_verified_for_both_candidates']
    assert result['branchwise_modular_charge_covariance_checked']
    control = result['uniform_global_total_two_control']
    assert F(control['local_restriction_trace_distance']) > F(2, 100)
    assert all(sum(map(F, row)) == 1 for row in control['local_total_probabilities'])
    assert F(control['vacuum_trace_distance'][1]) < F(control['vacuum_trace_distance'][0])


def test_exact_density_response_common_domain_and_first_discrepancy():
    result = c.density_certificate()
    assert result['common_exact_return_numerators_through_order'] == 22
    assert result['first_difference_order_within_executed_domain'] == 23
    assert F(result['first_difference_witness']['response_coefficient_difference']) != 0
    bound = result['single_initial_root_port_to_single_final_root_port_response_bound']
    assert 0 < F(bound['absolute_probability_difference_bound']) < F(1, 10**17)
    # The first coefficient is independently fixed by the actual seam census.
    first = result['common_exact_return_numerators'][1]
    graph = c.screen(3)
    for i, j in product(range(12), repeat=2):
        p, q = 12*63+i, 12*63+j
        expected = 12-len(graph.graph.neighbour_lists[p]) if i == j else int(q in graph.graph.neighbour_lists[p])
        assert first[i][j] == expected


def test_full_output_bound_has_declared_scheduler_and_state_domain():
    row = c.comparison_certificate(3, 5, 9)
    assert row['unconditional_complete_output_bounds'][1] == {
        'per_seam_Poisson_parameter': '1/4', 'trace_distance_allowance': '2187/5734400'}
    assert c.screen(5).seams == 16*c.screen(3).seams
    # Equal global counts would give unequal expected per-seam exposure.
    assert F(100, c.screen(3).seams) != F(100, c.screen(5).seams)


@pytest.mark.parametrize('value', (True, -1, 6, 3.0, '3'))
def test_invalid_source_level_rejected(value):
    with pytest.raises(ValueError):
        c.screen(value)
    with pytest.raises(ValueError):
        c.probe_program(value, 'C')


def test_retained_screen_receipt_replays():
    compare_receipt(strict_json.load(Path(__file__).with_name('screen_results.json')), c.run())
