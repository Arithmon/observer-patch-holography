"""Source custody and deliberate physical-premise boundaries."""
from pathlib import Path
from fractions import Fraction as F
from copy import deepcopy
import pytest
from tools import strict_json
from ..interacting_field.replay import compare_field_receipt
from .experiment import source_bindings


def test_retained_receipt_custody_and_scope():
    here = Path(__file__).resolve().parent
    receipt = strict_json.load(here/'results.json')
    campaign = strict_json.load(here/'campaign.json')
    assert receipt['source_bindings'] == source_bindings()
    assert receipt['all_pass'] and all(receipt['checks'].values())
    assert campaign['integrated_gate_accepted'] is False
    assert campaign['physical_comparison_count_added'] == 0
    assert receipt['physical_Lorentzian_or_continuum_attachment_claimed'] is False
    assert receipt['premise_register_dispositions_changed'] is False
    for row in receipt['levels']:
        assert row['regional_dimension'] == 11
        assert len(row['complete_record_probabilities']) == 4
        assert len(row['regional_sector_weights_by_record_N0_N1_N2']) == 4
        for probability, weights in zip(row['complete_record_probabilities'], row['regional_sector_weights_by_record_N0_N1_N2'], strict=True):
            assert len(weights) == 3
            assert abs(sum(weights)-probability) < 1e-12
        assert row['local_sector_postselection'] is False
        assert row['floating_residuals_are_outward_certified'] is False
    pair = receipt['nested_cube_comparison']
    assert pair['isometry_defect_nonzero_entries'] == pair['Hamiltonian_compression_defect_nonzero_entries'] == 0
    assert pair['Hamiltonian_intertwining_defect_Frobenius_norm'] > 0
    assert pair['source_complete_history_cq_distance_after_embedding'] < F(pair['ideal_compiled_complete_history_bound'])
    control = receipt['algebra_controls']
    assert control['multiplicativity_defect_operator_norm'] == 1
    assert control['opposite_phase_states_identical_separated_marginals'] is True
    assert control['opposite_phase_states_joint_half_trace_distance'] == 1


def test_receipt_replay_rejects_dropped_sector_or_changed_physical_claim():
    receipt = strict_json.load(Path(__file__).resolve().parent/'results.json')
    for mutate in (
        lambda r: r['levels'][0]['regional_sector_weights_by_record_N0_N1_N2'][0].pop(),
        lambda r: r.__setitem__('physical_Lorentzian_or_continuum_attachment_claimed', True),
        lambda r: r['nested_cube_comparison'].__setitem__('Hamiltonian_compression_defect_nonzero_entries', 1),
    ):
        changed = deepcopy(receipt); mutate(changed)
        with pytest.raises(ValueError):
            compare_field_receipt(receipt, changed)
