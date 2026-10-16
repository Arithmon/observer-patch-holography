from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest

from ..observer_reconstruction import native_geometry
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
    metric = result['exact_finite_spatial_metric']
    assert metric['exact_pair_counts'] == {'edge': 30, 'two_edge': 30, 'antipode': 6}
    assert metric['squared_chord_lengths_qsqrt5'] == {
        'edge': '2-2*sqrt(5)/5', 'two_edge': '2+2*sqrt(5)/5', 'antipode': '4'}
    assert metric['native_seams_are_exactly_shortest_shell']
    assert result['metric_shell_census_covariant_under_port_relabeling']
    assert not result['calibration_sector_is_same_as_binary_history_preparation']
    assert not result['physical_spatial_dimension_inferred']
    refinement = result['source_screen_refinement']
    assert refinement['source_levels'] == [3, 5]
    assert refinement['matched_root_graph_ball'] == {
        'radius_in_native_seam_hops': 9, 'ports': 222, 'induced_seams': 561,
        'canonical_address_and_adjacency_equal': True,
        'root_port_metric_is_exactly_the_same_unitless_chord_metric': True}
    assert refinement['root_graph_hop_ball_port_counts_radii_zero_through_ten'] == [
        12, 15, 30, 45, 54, 84, 114, 132, 186, 222, 243]
    growth = refinement['full_federation_native_hop_growth']
    profiles = growth['profiles']
    assert [row['level'] for row in profiles] == list(range(6))
    assert [row['root_port_eccentricity_in_native_seam_hops'] for row in profiles] == [
        13, 31, 60, 121, 242, 484]
    level5 = profiles[5]
    assert [level5['port_ball_counts_by_radius'][radius]
            for radius in (16, 32, 64, 128, 256)] == [576, 2250, 8805, 37224, 123926]
    assert [row['port_ball_growth_ratio'] for row in growth['exact_doubling_ratios']
            if row['level'] == 5] == ['125/32', '587/150', '12408/2935', '5633/1692']
    assert growth['not_a_calibrated_physical_length_or_asymptotic_dimension']
    assert refinement['first_induced_graph_difference']['radius'] == 10
    assert refinement['existing_federation_refinement_transport'][
        'path_independent_port_frame_or_metric_scaling_recovered'] is False
    for row in refinement['existing_federation_refinement_transport']['level_certificates']:
        assert row['cell_parent_graph_morphism']
        assert not row['port_parent_graph_morphism']
        assert row['uniform_per_edge_laplacian_intertwining_defect_inf'] == '1'
    routing = refinement['existing_intercarrier_routing']
    assert not routing['routing_is_derived_from_native_observer_records']
    assert not routing['record_derived_path_independent_frame_transition_available']
    assert [row['antipodal_consistent_fraction_exact'] for row in routing['level_receipts']] == [
        '113/160', '5951/7680']
    control = refinement['intercarrier_frame_nonuniqueness_control']
    assert control['conditional_equal_shortest_shell_length_used_only_as_control']
    assert control['same_local_gram_and_same_link_length_leave_relative_frame_underdetermined']
    assert control['max_local_gram_residual_diagnostic'] < 2e-12
    assert control['max_cross_carrier_gram_entry_difference_diagnostic'] > 1e-3
    assert control['two_candidate_squared_link_lengths'] == [
        '2-2/sqrt(5)', '2-2/sqrt(5)']
    assert not refinement['physical_global_space_or_continuum_recovered']


def test_exact_spatial_shell_audit_rejects_a_symmetric_wrong_edge_length():
    recovered = native_geometry.recover(native_geometry.execute_calibration())
    gram = [[list(entry) for entry in row] for row in recovered['gram_qsqrt5']]
    seams = tuple((i, j) for i in range(12) for j in range(i+1, 12)
                  if F(recovered['laplacian_exact'][i][j]) == -1)
    # Keep symmetry but move one source edge to the opposite non-edge shell.
    gram[0][1] = gram[1][0] = ['0', '-1/5']
    with pytest.raises(ValueError, match='graph-distance shell'):
        c.spatial_shells(gram, seams)


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
