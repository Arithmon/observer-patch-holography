"""Retained source custody and scientific scope, separate from numerical replay."""
from fractions import Fraction as F
from pathlib import Path
from copy import deepcopy
import math
import pytest
from tools import strict_json
from .experiment import source_bindings
from .replay import compare_field_receipt, LOCAL_RESIDUAL_REPLAY_ATOL


def test_retained_evidence_binds_every_current_source_and_contract():
    here=Path(__file__).resolve().parent
    receipt=strict_json.load(here/'results.json')
    assert receipt['source_bindings']==source_bindings()
    contract=strict_json.load(here/'campaign.json')
    assert contract['integrated_gate_accepted'] is False
    assert contract['physical_comparison_count_added']==0
    assert contract['base_advance']['head']==receipt['base_advance']['head']
    assert contract['execution']['Gaussian_evolution_of_interacting_state'] is False
    assert contract['execution']['all_requested_routing_pulses_numerically_replayed'] is False
    for row in receipt['levels']:
        assert len(row['complete_two_read_probabilities'])==4
        assert row['evolution_resources']['floating_residual_estimate_outward_certified'] is False
        assert row['apparatus_resources']['floating_residual_estimate_outward_certified'] is False
        assert row['evolution_resources']['routing_FSWAP_calls']>0
        assert row['short_interaction_certificate']['floating_or_hardware_error_certified'] is False
        assert F(row['short_interaction_certificate']['ideal_executed_interaction_minus_free_probability_lower'])>0
        temporal=F(row['temporal_certificate']['complete_history_half_trace_upper'])
        spatial=F(row['spatial_certificate']['complete_history_half_trace_upper'])
        assert F(row['ideal_source_to_infinite_lattice_complete_history_bound'])==min(F(1),temporal+spatial)
    assert receipt['net_node_chirality']==0
    assert len(receipt['naive_Weyl_nodes'])==8
    assert receipt['free_low_momentum_certificate']['interaction_projection_covered'] is False


def retained():
    return strict_json.load(Path(__file__).resolve().parent/'results.json')


def test_replay_allows_only_count_scaled_local_roundoff():
    expected = retained()
    actual = deepcopy(expected)
    resources = actual['levels'][-1]['evolution_resources']
    calls = sum(resources[key] for key in (
        'completed_hopping_blocks', 'completed_contact_blocks', 'routing_FSWAP_calls'))
    value = calls*LOCAL_RESIDUAL_REPLAY_ATOL/2
    assert abs(value-resources['floating_local_residual_sum_with_routing']) > 2e-10
    resources['floating_local_residual_sum_with_routing'] = value
    resources['floating_residual_telescoping_estimate'] = math.expm1(value)
    compare_field_receipt(expected, actual)
    assert expected == retained()  # Comparison never edits either retained input.


@pytest.mark.parametrize('value', [-1e-12, 1., float('inf'), float('nan')])
def test_replay_rejects_invalid_or_excessive_accumulated_residual(value):
    expected = retained()
    actual = deepcopy(expected)
    actual['levels'][-1]['evolution_resources']['floating_local_residual_sum_with_routing'] = value
    with pytest.raises(ValueError, match='per-call cap'):
        compare_field_receipt(expected, actual)


def test_replay_cannot_hide_counts_observables_or_false_certification():
    expected = retained()
    for container, key, value in (
        ('evolution_resources', 'routing_FSWAP_calls', 100_000_000),
        ('evolution_resources', 'floating_residual_telescoping_estimate', 1e-6),
        ('evolution_resources', 'floating_residual_estimate_outward_certified', True),
        (None, 'source_reference_history_TV', 1e-6),
        (None, 'ideal_source_to_infinite_lattice_complete_history_bound', '0'),
    ):
        actual = deepcopy(expected)
        row = actual['levels'][-1]
        (row if container is None else row[container])[key] = value
        with pytest.raises(ValueError):
            compare_field_receipt(expected, actual)


def test_replay_also_rejects_corrupt_retained_diagnostic():
    expected = retained()
    actual = deepcopy(expected)
    expected['levels'][-1]['evolution_resources']['floating_local_residual_sum_with_routing'] = -1e-12
    with pytest.raises(ValueError, match='per-call cap'):
        compare_field_receipt(expected, actual)
