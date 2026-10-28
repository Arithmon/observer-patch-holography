"""Independent native completeness, shared-event and readout controls."""
from fractions import Fraction as F
import hashlib

import numpy as np
import pytest

from tools import strict_json
from model.completion import observer
from model.completion.family import FamilyParameters
from model.clock_defect_energy import continuation, joint_clock
from . import run as shared


def test_four_field_checkpoints_are_actual_native_density_recurrence():
    # Independent exact arc-walk oracle, not the source routine or a fit.
    expected = (F(1, 22), F(19, 110), F(3529, 9900), F(41609, 74250))
    for n, probability in enumerate(expected, 1):
        amplitudes = joint_clock.wave.ramsey_outputs(joint_clock.FIELD, n, -1)
        assert joint_clock.wave.norm_squared(amplitudes[1]) / 330 == probability
        assert joint_clock.wave.norm_squared(amplitudes[0]) / 330 == 1-probability


def test_exact_model_readouts_inhabit_strict_finite_order_but_are_not_affine():
    bounds = {n: F(shared._budget(n)['conditional_joint_read_TV_upper'])
              for n in range(1, 5)}
    row = shared.exact_comparison(bounds)
    clock = [tuple(map(F, row['clock_first_code_probabilities'][str(n)]))
             for n in range(1, 5)]
    field = [F(row['field_Ramsey_minus_probabilities'][str(n)]) for n in range(1, 5)]
    assert all(clock[i][1] < clock[i+1][0] and field[i] < field[i+1]
               for i in range(3))
    assert F(row['positive_affine_two_event_rate'][0]) > 0
    assert F(row['third_event_cross_product_nonzero'][0]) > F(1, 1000)
    assert F(row['source_bound_robust_third_cross_product'][0]) > F(1, 1000)
    assert F(row['source_bound_robust_field_probability_above_affine_prediction'][0]) > F(9, 1000)
    assert row['raw_affine_test_passed'] is False
    s4 = [tuple(map(F, row['parameter_four_control_clock_probabilities'][str(n)]))
          for n in range(1, 5)]
    assert s4[1][0] > s4[2][1]  # inherited parameter-four service aliases


@pytest.mark.parametrize('law', ('C', 'E'))
def test_preferred_native_menu_is_complete_on_work_and_untouched_reference(law):
    params = FamilyParameters(c_weight=F(int(law == 'C')))
    menu = observer.branches((), params, preferred=shared.SEAMS[0])
    columns = np.zeros((66, 2), complex)
    columns[55] = np.array((1, 1j))/2
    columns[65] = np.array((1, 1))/2
    initial = columns.conj().T@columns
    assert abs(initial[0, 1]) > .1
    returned = np.zeros((2, 2), complex)
    supported = 0
    for branch in menu:
        output = np.zeros_like(columns)
        for i, j in branch.pairs:
            if i in continuation.FULL_TO_BINARY and branch.weight:
                output[continuation.FULL_TO_BINARY[j]] += (
                    float(branch.weight)**.5 * columns[continuation.FULL_TO_BINARY[i]])
        returned += output.conj().T@output
        supported += int(np.linalg.norm(output) > 0)
    np.testing.assert_allclose(returned, initial, atol=1e-14, rtol=0)
    assert supported > 1
    assert shared._selected(menu, shared.SEAMS[0], *shared.COINS_LABELS[0]).weight == F(1, 4)
    assert any(branch.record.selector == 'explore' and branch.weight > 0 for branch in menu)


def test_retained_source_campaign_result_and_public_boundaries():
    campaign = shared.contract()
    receipt = strict_json.load(shared.OUTPUT)
    assert receipt['campaign_id'] == campaign['campaign_id']
    assert receipt['campaign_sha256'] == hashlib.sha256(shared.CAMPAIGN.read_bytes()).hexdigest()
    assert receipt['source_sha256'] == shared.source_bindings(receipt['source_sha256'])
    assert receipt['physical_seconds_calibrated'] is False
    assert receipt['preparation_finite_mixed_reference_flags_reused'] is False
    assert receipt['parameter_one_new_clock_words'] == 4
    for law in receipt['source_complete_instruments']:
        assert law['pure_source_probe_is_faithful_reference_preparation'] is False
        events = [row['shared_event_sha256'] for row in law['checkpoints']]
        assert len(events) == len(set(events)) == 4
        for n, row in enumerate(law['checkpoints'], 1):
            assert row['checkpoint'] == n
            assert row['exported_record_count'] == n-1
            assert row['retained_fifo_size'] == 1
            assert row['source_clock_word_pulses'] == 2042
            assert row['terminal_joint_clock_field_read_slots'] == 4
            assert row['clock_export_abort_slots_through_checkpoint'] == n
            assert F(row['source_complete_mass_diagnostic'][0]) <= 1 <= F(row['source_complete_mass_diagnostic'][1])
            selected = tuple(map(F, row['source_selected_path_mass_diagnostic']))
            assert selected[0] < F(1, 4**n) < selected[1]
            assert F(row['bounds']['conditional_joint_read_TV_upper']) < F(1, 1000)
    assert [r['shared_event_sha256'] for r in receipt['source_complete_instruments'][0]['checkpoints']] != [
        r['shared_event_sha256'] for r in receipt['source_complete_instruments'][1]['checkpoints']]


def test_wrong_event_or_time_type_cannot_be_silently_reused():
    row = shared.exact_comparison({n:F(shared._budget(n)['conditional_joint_read_TV_upper'])
                                   for n in range(1, 5)})
    assert row['physical_affine_timelike_atlas_supplied'] is False
    assert row['formal_source_specific_Lean_instance_claimed'] is False
    assert row['calibrated_count_lookup_affine_by_construction'] is True
    with pytest.raises(ValueError):
        shared.continuation._core((0, 0))
    with pytest.raises(ValueError):
        joint_clock.field_probabilities(3)  # old two-attempt readout domain
