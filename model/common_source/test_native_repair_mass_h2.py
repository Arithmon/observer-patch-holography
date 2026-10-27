"""Full-native moment and preservation checks for the fixed two-hop control."""
import sympy as s
import pytest

from . import native_repair_mass_h2 as h


@pytest.mark.parametrize('power,sign', [(1, 1), (1, -1), (2, 1), (2, -1)])
def test_complete_native_rows_reproduce_both_unprojected_moments(power, sign):
    factors = h.complete_native_factors(power, sign)
    f = s.kronecker_product(factors['left_encoder_integer'], factors['right_encoder_integer'])
    f = f*s.SparseMatrix(factors['frame'])/factors['encoder_divisor']
    diagonal = s.SparseMatrix(4356, 4356,
        {(i, i): value for i, value in enumerate(factors['diagonal']) if value})
    assert f.shape == (4356, 6)
    assert f.T*f == h.geometry()['metric']
    m = h.moments(power, sign)
    assert f.T*diagonal*f == m['first']
    assert f.T*diagonal**2*f == m['second']


def test_raw_and_color_complements_have_separate_exact_positive_defects():
    g = h.geometry()
    for sign, full, color in ((1, s.Rational(20353, 23592960), s.Rational(11, 2949120)),
                              (-1, s.Rational(20369, 23592960), s.Rational(13, 2949120))):
        m = h.moments(2, sign)
        assert m['leakage'] == m['syndrome_leakage']+m['color_leakage']
        assert m['leakage'].rank() == m['color_leakage'].rank() == 6
        assert s.trace(g['metric_inverse']*m['leakage']) == full
        assert s.trace(g['metric_inverse']*m['color_leakage']) == color
        assert m['leakage'].is_positive_definite


def test_no_tuning_of_two_local_costs_closes_the_geometric_domain():
    assert h.preserving_coefficient_rank() == 2
    # Scalar compressed H1 is not a closed native evolution either.
    m, g = h.moments(1), h.geometry()
    assert g['metric_inverse']*m['first'] == s.Rational(9, 16)*s.eye(6)
    assert s.trace(g['metric_inverse']*m['leakage']) == s.Rational(79, 768)


def test_report_keeps_native_domain_and_mass_identification_limits():
    r = h.report()
    assert r['complete_binary_native_total_load_four_rows'] == (
        r['diagonal_invariant_binary_load_two_two_rows']+r['retained_zero_binary_rows_in_other_load_sectors'])
    assert not r['scope']['compiled_source_history_executed_here']
    assert not r['scope']['Higgs_family_coupling_or_physical_particle_mass_derived']
    assert not r['scope']['geometric_current_triplets_identified_as_families']


@pytest.mark.parametrize('power,sign', [(True, 1), (3, 1), (2, True), (2, 0)])
def test_original_cost_domain_is_not_silently_changed(power, sign):
    with pytest.raises(ValueError):
        h.complete_native_factors(power, sign)
