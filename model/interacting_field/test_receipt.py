"""Retained source custody and scientific scope, separate from numerical replay."""
from fractions import Fraction as F
from pathlib import Path
from tools import strict_json
from .experiment import source_bindings


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
