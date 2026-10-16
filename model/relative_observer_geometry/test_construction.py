from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest

from . import construction as c


INITIAL = tuple(int(i in (0, 11)) for i in range(12))


def test_every_native_branch_observer_access():
    assert c.branch_audit()['native_branches_checked_exactly'] == 3960


def test_full_runtime_and_actual_checkpoint_continuation():
    result = c.runtime_check()
    assert result['retained_record_words'] == 216
    assert result['positive_record_words'] == 8
    assert result['zero_probability_record_words'] == 208
    assert result['checkpoint_complete_distance_diagnostic'] < 2e-11
    assert abs(result['public_n0_expectation']-result['hidden_n0_audit_expectation']) < 2e-11


def test_instrument_retains_coherence_within_record_sectors():
    seam = (0, 1)
    instrument = c.seam_instrument(seam)
    # A coherent pair is retained within the same committed outcome: the
    # instrumentation does not silently add a configuration measurement.
    rng = np.random.default_rng(1107)
    state = rng.normal(size=66)+1j*rng.normal(size=66)
    state /= np.linalg.norm(state)
    density = np.outer(state, state.conj())
    output = sum(k@density@k.conj().T for group in instrument.kraus for k in group)
    assert abs(np.trace(output)-1) < 1e-12
    assert np.linalg.norm(output-np.diag(np.diag(output))) > .1


def test_boundary_current_is_necessary_and_hidden_initial_state_is_ambiguous():
    controls = c.controls()
    assert controls['omitted_boundary_current'] == {'correct_n0': 0, 'incorrect_n0': 1}
    assert controls['same_record_hidden_preparations']['possible_initial_n0'] == [1, 0]
    law = controls['same_disjoint_regional_marginals_different_global_laws']
    assert law['E_n0_n10_A'] != law['E_n0_n10_B']


def test_finite_complete_history_law_and_checkpoint():
    result = c.history_check()
    assert result['complete_histories'] == 2112
    assert sum(F(row['probability']) for row in result['joint_records']) == 1
    assert result['all_local_checkpoint_continuations_exact']


def test_three_shared_events_distinguish_local_incident_and_shared_seam_clocks():
    result = c.clock_checks()
    assert result['local_incident_tick_pairs'] == [[1, 1, 1], [2, 2, 2], [5, 5, 3]]
    assert result['predesignated_third_event_residual'] == '-2'
    assert result['affine_shared_seam_third_residual'] == '0'


def test_affine_comparison_with_exact_fraction_readings():
    slope, offset, residual = c.affine_comparison(((1, F(1, 3), F(2, 3)),
        (4, F(2, 3), F(4, 3)), (9, F(5, 3), F(10, 3))))
    assert (slope, offset, residual) == (F(2), F(0), (F(0),))


@pytest.mark.parametrize('shared', [(), ((1, 0, 1),), ((1, 0, 1), (1, 1, 2)),
    ((1, 0, 1), (2, 0, 2)), ((1, 0, 1), (2, 1, 1)),
    ((1, 0, 1), (2, True, 2)), ((1, 0., 1), (2, 1, 2))])
def test_clock_domain_rejects_duplicate_events_and_nonclocks(shared):
    with pytest.raises(ValueError):
        c.affine_comparison(shared)


def test_checkpoint_does_not_accept_new_fields_or_duplicate_records():
    _, records = c.execute(INITIAL, c.TAPE, (0, 1, 0, 1, 0))
    view = c.local(INITIAL, records[:2], (0,))
    saved = view.checkpoint()
    changed = deepcopy(saved); changed['hidden_work'] = list(INITIAL)
    with pytest.raises(ValueError):
        c.LocalCounter.restore(changed)
    with pytest.raises(ValueError):
        view.continue_with((records[1],))
    changed = deepcopy(saved); changed['records'][0][4] = True
    with pytest.raises(ValueError):
        c.LocalCounter.restore(changed)
    assert c.LocalCounter.restore(saved) == view


def test_cover_and_overlap_are_required_for_value_gluing():
    a = c.LocalCounter((0, 1), (1, 0))
    b = c.LocalCounter((1, 2), (1, 0))
    with pytest.raises(ValueError, match='overlap'):
        c.glue((a, b), (0, 1, 2))
    with pytest.raises(ValueError, match='cover'):
        c.glue((a,), (0, 1, 2))
    b = c.LocalCounter((1, 2), (0, 1))
    assert c.glue((a, b), (0, 1, 2)) == (1, 0, 1)


@pytest.mark.parametrize('record', [(1, 0, 5, 0, -1), (1, 1, 0, 0, 1),
    (0, 0, 5, 0, 1), (1, 0, 5, True, 1), (1, 0, 5, 0, 2),
    (1, 0, 1, 1, 1)])
def test_false_or_impossible_local_records_rejected(record):
    with pytest.raises(ValueError):
        c.LocalCounter((0,), (1,), (record,))


def test_known_two_endpoints_check_committed_current_against_actual_coin():
    with pytest.raises(ValueError, match='endpoint loads'):
        c.LocalCounter((0, 1), (1, 0), ((1, 0, 1, 0, 0),))
    assert c.LocalCounter((0, 1), (1, 0), ((1, 0, 1, 0, 1),)).values == (0, 1)


def test_record_calibration_graph_is_same_current_interface():
    result = c.geometry_connection()
    assert len(result['recovered_seams']) == 30
    assert result['rank_of_certified_slow_response_Gram'] == 3
    assert not result['calibration_sector_is_same_as_binary_history_preparation']
    assert not result['physical_spatial_dimension_inferred']


@pytest.mark.parametrize('bad', [((False, 1),), ((1, 0),), ([0, 1],)])
def test_wrong_seams_rejected(bad):
    with pytest.raises(ValueError):
        c.execute(INITIAL, bad, (0,))


def test_local_fields_do_not_disclose_the_complement():
    _, records = c.execute(INITIAL, ((0, 1), (0, 5)), (0, 0))
    view = c.local(INITIAL, records, (11,))
    assert view.records == ()
    assert view.initial == (1,)
    assert view.values == (1,)
    assert view.ticks == 0
