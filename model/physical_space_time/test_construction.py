from fractions import Fraction as F

import pytest
import sympy as sp

from . import construction as c


@pytest.fixture(scope='module')
def source():
    return c.kernel(0)


def test_actual_native_branches_match_independently_reconstructed_response(source):
    _, laplacian = c.recovered_carrier()
    assert source['native_branch_checks'] == 720
    assert source['attempt'] == sp.eye(12)-laplacian/60
    assert source['tick_probability'] == sp.Rational(1, 6)


def test_remote_word_sum_has_explicit_probability_remainder(source):
    # A finite independently summed family of remote words is a lower bound
    # on the completed law. Its missing column mass is exactly (5/6)^m.
    partial, power = sp.zeros(12), sp.eye(12)
    for _ in range(5):
        partial += source['incident']*power
        power = source['remote']*power
    remainder = source['completed']-partial
    assert all(value >= 0 for value in remainder)
    assert sp.ones(1, 12)*remainder == sp.Rational(5, 6)**5*sp.ones(1, 12)


def test_first_arrival_clock_means_and_all_distance_shells(source):
    attempts = c.mean_arrival(source['attempt'], 0)
    incident = c.mean_arrival(source['completed'], 0)
    distances = c.graph_distances()[0]
    assert incident == attempts/6
    assert {distances[i]: incident[i] for i in range(12)} == {0: 0, 1: 22, 2: 28, 3: 30}
    assert {distances[i]: attempts[i] for i in range(12)} == {0: 0, 1: 132, 2: 168, 3: 180}


def test_first_arrival_capped_records_keep_all_timeout_mass(source):
    law = c.arrival_law(source['completed'], 3, 0, 4)
    assert F(law['first_arrival_probabilities'][0]) == F(1, 222)
    assert sum(map(F, law['first_arrival_probabilities']))+F(law['timeout_probability']) == 1
    assert F(law['timeout_probability']) > 0
    assert not law['postselection_or_timeout_renormalization']
    global_law = c.arrival_law(source['attempt'], 3, 0, 4)
    assert list(map(F, global_law['first_arrival_probabilities'][:2])) == [0, 0]
    assert F(global_law['first_arrival_probabilities'][2]) > 0


def test_state_independent_geometric_tail_is_uniform_over_unhit_states(source):
    assert min(source['completed'][0, i] for i in range(1, 12)) == sp.Rational(1, 222)
    for sender in range(1, 12):
        law = c.arrival_law(source['completed'], sender, 0, 5)
        assert F(law['timeout_probability']) <= F(221, 222)**5


def test_identity_holding_changes_attempt_clock_but_not_local_experiment(source):
    held = c.kernel(0, holding=F(2, 7))
    assert held['completed'] == source['completed']
    assert c.mean_arrival(held['attempt'], 0) == sp.Rational(7, 2)*c.mean_arrival(source['attempt'], 0)
    assert held['tick_probability'] == sp.Rational(1, 21)


def test_nonuniform_exposure_changes_operational_passage_law():
    seams, _ = c.recovered_carrier()
    weights = {edge: F(1+int(edge == seams[0]), 31) for edge in seams}
    first, second = c.kernel(0, seam_weights=weights), c.kernel(1, seam_weights=weights)
    assert c.mean_arrival(first['completed'], 0) != c.mean_arrival(c.kernel(0)['completed'], 0)
    for target, data in ((0, first), (1, second)):
        assert c.mean_arrival(data['completed'], target) == data['tick_probability']*c.mean_arrival(data['attempt'], target)


def test_total_two_witness_invalidates_one_edge_per_incident_tick():
    first, second = (c.native_witness(tuple(int(i in occupied) for i in range(12)), 0,
                    ((2, 10), (0, 10)), (0, 1)) for occupied in ((1, 2), (1, 3)))
    assert first['receiver_field_by_native_attempt'] == [0, 0, 1]
    assert second['receiver_field_by_native_attempt'] == [0, 0, 0]
    assert first['incident_tick_count'] == second['incident_tick_count'] == 1
    assert min(c.graph_distances()[0][i] for i in (2, 3)) == 2
    assert first['accessible_incident_records'][0]['current'] == -1
    assert second['accessible_incident_records'][0]['current'] == 0


def test_reversing_native_signal_word_removes_the_same_clock_signal():
    for occupied in ((1, 2), (1, 3)):
        result = c.native_witness(tuple(int(i in occupied) for i in range(12)), 0,
                                 ((0, 10), (2, 10)), (1, 0))
        assert result['receiver_field_by_native_attempt'] == [0, 0, 0]


@pytest.mark.parametrize('bad', [True, -1, 12, F(0)])
def test_invalid_receiver_is_rejected(bad):
    with pytest.raises(ValueError, match='native integer port'):
        c.kernel(bad)


@pytest.mark.parametrize('bad', [F(0), F(-1), F(2), 1, True])
def test_invalid_holding_service_is_rejected(bad):
    with pytest.raises(ValueError, match='holding acceptance'):
        c.kernel(0, holding=bad)


def test_missing_seam_exposure_and_nonterminating_arrival_are_rejected():
    with pytest.raises(ValueError, match='every native seam'):
        c.kernel(0, seam_weights={(0, 1): F(1)})
    with pytest.raises(ValueError, match='finite mean'):
        c.mean_arrival(sp.eye(12), 0)


def test_receipt_keeps_the_physical_gate_open():
    result = c.run()
    assert result['all_pass']
    assert result['derived_equalities']['centered_passage_gram_rank'] == 11
    assert result['acceptance']['finite_operational_signal_and_local_clock_laws_established']
    assert not result['acceptance']['physical_three_plus_one_dimensional_spacetime_established']
    assert not result['acceptance']['lab_clock_identification_or_empirical_comparison_established']
    assert not result['contract']['federation_required']
