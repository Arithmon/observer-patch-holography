"""Independent logarithms, exhaustive gap assignments and attachment controls."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product

import pytest

from . import modular_optical as m


def test_source_gaps_enclose_independent_high_precision_logs():
    with localcontext() as context:
        context.prec = 80
        root = Decimal(3).sqrt()
        plus = -((25+root)/30).ln()
        minus = -((25-root)/30).ln()
        expected = plus, minus, minus-plus
        for actual, (lo, hi) in zip(expected, m.gaps(), strict=True):
            assert Decimal(lo.numerator)/Decimal(lo.denominator) <= actual
            assert actual <= Decimal(hi.numerator)/Decimal(hi.denominator)
            assert hi-lo <= F(1, 10**12)
        for a, b, (lo, hi) in m.assignment_intervals():
            actual = expected[m.LABELS.index(a)]/expected[m.LABELS.index(b)]
            assert Decimal(lo.numerator)/Decimal(lo.denominator) <= actual
            assert actual <= Decimal(hi.numerator)/Decimal(hi.denominator)


def test_every_ordered_pair_and_same_gap_control_is_retained():
    rows = m.assignment_intervals()
    assert {(a, b) for a, b, _ in rows} == set(product(m.LABELS, repeat=2))
    assert len(rows) == 9
    assert all(interval == (F(1), F(1)) for a, b, interval in rows if a == b)
    for _, _, (lo, hi) in rows:
        compatible = m.compare_ratio((lo+hi)/2, error=F(1, 10**9))
        assert not compatible['all_nine_assignments_incompatible']


@pytest.mark.parametrize('isotope,ratio', [(40, '1.05776938758748094'),
                                         (36, '1.05776646273518748')])
def test_published_ratios_have_certified_minimum_discrepancy(isotope, ratio):
    result = m.compare_ratio(F(ratio))
    assert result['all_nine_assignments_incompatible']
    assert F(result['minimum_absolute_ratio_separation']) == F(ratio)-1-m.ERROR_BUDGET
    assert F(result['minimum_absolute_ratio_separation']) > F(57, 1000)
    # A sufficiently permissive attachment-error hypothesis is not rejected.
    assert not m.compare_ratio(F(ratio), error=F(3, 50))['all_nine_assignments_incompatible']


def test_common_units_cancel_but_independent_species_gains_destroy_the_test():
    lo, hi = m.gaps()[0]
    a = (lo+hi)/2
    lo, hi = m.gaps()[2]
    b = (lo+hi)/2
    for scale in (F(1, 17), F(31, 9), F(10**18)):
        assert (scale*a)/(scale*b) == a/b
    measured = F('1.05776938758748094')
    species_gain = measured/(a/b)
    assert species_gain*a/b == measured
    assert species_gain != 1  # This explicitly excluded premise can fit any ratio.


def test_rounding_is_outward_and_exact_on_decimal_grid():
    assert m.outward((F(1, 3), F(2, 3)), 2) == (F(33, 100), F(67, 100))
    assert m.outward((F(3, 10), F(7, 10)), 2) == (F(3, 10), F(7, 10))


@pytest.mark.parametrize('bad', [True, 0.1, -1, 0])
def test_invalid_measured_ratio_is_rejected(bad):
    with pytest.raises(ValueError):
        m.compare_ratio(bad)


def test_noncanonical_fraction_and_invalid_error_or_precision_rejected():
    malformed = F(1, 2)
    malformed._numerator = 2
    for bad in (True, 0.1, -1, malformed, F(2)):
        with pytest.raises(ValueError):
            m.compare_ratio(F(1), error=bad)
    for bad in (0, True, 1.0, 31):
        with pytest.raises(ValueError):
            m.outward((F(1, 3), F(2, 3)), bad)
    for interval in ([F(1), F(2)], (F(2), F(1)), (F(0), F(1))):
        with pytest.raises(ValueError):
            m.outward(interval)


def test_receipt_uses_primary_table_and_makes_only_conditional_claims():
    result = m.run()
    assert result['observation']['primary_url'] == 'https://arxiv.org/pdf/2205.13053v2#page=20'
    assert result['observation']['units'] == 'dimensionless'
    assert result['observation']['source_version_date'] == '2023-09-06'
    assert [row['isotope'] for row in result['comparisons']] == [40, 36]
    assert not result['interpretation']['quoted_standard_uncertainties_treated_as_hard_bounds']
    assert not result['interpretation']['p_value_or_sigma_significance_claimed']
    assert not result['interpretation']['frozen_prospective_prediction']
    assert not result['interpretation']['OPH_or_all_modular_clocks_excluded']
    assert not result['attachment']['atomic_spectrum_or_transition_dictionary_derived']
