"""Exact source-symmetry selection and finite native scheduling control."""
from fractions import Fraction as F

import pytest

from . import symmetry_schedule as schedule
from ..common_source import slow_geometry


def test_source_automorphisms_make_all_native_seams_one_orbit():
    seams = slow_geometry.source().seams()
    assert len(schedule.edge_orbit(seams)) == 30
    # The first generator alone leaves several edge orbits; the joint
    # certificate, rather than 30 checked uniform weights, does the work.
    assert len(schedule.edge_orbit(seams, schedule.GENERATORS[:1])) == 10
    corrupted = list(schedule.GENERATORS)
    corrupted[0] = tuple(range(11, -1, -1))
    with pytest.raises(ValueError, match='preserve native seams'):
        schedule.edge_orbit(seams, corrupted)


def test_symmetry_and_stationarity_do_not_fix_native_current_history():
    audit = schedule.run()
    iid = audit['iid_uniform']
    persistent = audit['uniform_marginal_persistent']
    assert audit['normalized_symmetry_invariant_one_step_seam_weight'] == '1/30'
    assert audit['same_stationary_uniform_configuration_law']
    assert iid['stationary_moving_current_probability_each_attempt'] == '5/33'
    assert persistent['stationary_moving_current_probability_each_attempt'] == '5/33'
    assert F(iid['two_successive_moving_currents_probability']) == F(23, 990)
    assert F(persistent['two_successive_moving_currents_probability']) == F(49, 990)
    assert not audit['two_attempt_native_current_laws_equal']
    assert not audit['symmetry_or_stationarity_selects_iid_schedule']
    assert audit['repeat_probability_identifiable_from_exact_two_current_law_within_declared_family']
    assert not audit['finite_sample_uncertainty_or_scheduler_family_selected']
    assert not audit['physical_time_or_particle_mass_identified']


@pytest.mark.parametrize('q', [F(0), F(1, 3), F(1, 2), F(3, 4)])
def test_full_two_bit_law_is_positive_with_same_one_bit_marginal(q):
    result = schedule.persistent_seam_current_law(q)
    p = F(result['stationary_moving_current_probability_each_attempt'])
    both = F(result['two_successive_moving_currents_probability'])
    word_law = (1-2*p+both, p-both, p-both, both)
    assert all(x >= 0 for x in word_law)
    assert sum(word_law) == 1
    assert both == (F(23)+52*q)/990
    assert schedule.repeat_probability_from_current_pair(both) == q
    assert result['seam_choices_independent'] == (q == 0)


@pytest.mark.parametrize('q', [0, .5, True, F(-1, 2), F(1), F(3, 2)])
def test_scheduler_countercontrol_requires_exact_ergodic_rational_parameter(q):
    with pytest.raises(ValueError):
        schedule.persistent_seam_current_law(q)


@pytest.mark.parametrize('probability', [0, .04, True, F(1, 100), F(5, 66)])
def test_current_pair_inverse_rejects_incompatible_laws(probability):
    with pytest.raises(ValueError):
        schedule.repeat_probability_from_current_pair(probability)
