"""Collective masks, complete waiting laws, and independent native field replay."""
from fractions import Fraction as F

import pytest
import sympy as sp

from . import regional_clock as r


@pytest.fixture(scope='module')
def receipt():
    return r.run()


def test_renewal_factorization_and_rate_normalization():
    z, a, b, c = sp.symbols('z a b c', positive=True)
    law = (a-c+c*z)/(a+b-c-(b-c)*z)
    factor = (1-c/a+c*z/a)*(a/(a+b-c))/(1-(b-c)*z/(a+b-c))
    assert sp.simplify(law-factor) == 0
    mean = sp.diff(law, z).subs(z, 1)
    variance = sp.diff(law, z, 2).subs(z, 1)+mean-mean**2
    assert sp.simplify(mean-b/a) == 0
    assert sp.simplify(variance-b*(a+b-2*c)/a**2) == 0
    assert sp.simplify(variance*(a/b)**2-(a+b-2*c)/b) == 0


def test_every_nonempty_region_and_ordered_distinct_pair_is_accounted_for(receipt):
    catalogue = r._catalogue()
    assert len(catalogue) == 3433
    all_regions = [region for _, aliases in catalogue for region in aliases]
    assert len(all_regions) == len(set(all_regions)) == 4095
    assert {sum(1 << port for port in region) for region in all_regions} == set(range(1, 4096))
    assert all(r.event_mask(region) == mask for mask, aliases in catalogue for region in aliases)
    classes = r.pair_classes()
    assert len(classes) == 2848
    assert sum(n for _, _, _, n in classes) == 3433*3432 == 11782056
    assert receipt['ordered_distinct_event_set_pairs'] == 11782056
    for a, b, c, _ in classes:
        assert a+b-2*c >= 1 and b <= 30
        assert r.parameters(a, b, c)['normalized_ratio_variance_per_A'] >= F(1, 30)
        if a == b:
            assert a+b-2*c >= 2 and b <= 29
            assert r.parameters(a, b, c)['normalized_ratio_variance_per_A'] >= F(2, 29)


def test_both_fluctuation_floors_are_attained_on_one_source(receipt):
    first = tuple(p for p in range(12) if p not in (0, 1))
    second = tuple(p for p in range(12) if p not in (0, 5))
    complete = tuple(range(12))
    assert r.incidence_counts(first, complete) == (29, 30, 29)
    assert r.parameters(*r.incidence_counts(first, complete))['normalized_ratio_variance_per_A'] == F(1, 30)
    assert r.incidence_counts(first, second) == (29, 29, 28)
    assert receipt['sharp_equal_mean_rate_coefficient'] == '2/29'


def test_clock_aliases_are_pathwise_identical_without_independent_stability_claim(receipt):
    first, second = map(tuple, receipt['alias_witness'])
    assert first != second and r.event_mask(first) == r.event_mask(second) == (1 << 30)-1
    assert receipt['all_event_clock_alias_regions'] == 69
    for count in (1, 3, 9):
        bins, tail = r.count_law(first, second, count, maximum=12)
        assert bins == tuple(F(i == count) for i in range(13)) and tail == 0
        assert r.pgf(first, second, F(3, 7), count) == F(3, 7)**count
    assert r.count_law(first, second, 9, maximum=3) == ((F(0),)*4, F(1))


@pytest.mark.parametrize('first,second', [((0,), (1,)), ((0,), tuple(range(12))),
    (tuple(range(12)), (0,)), ((0, 2, 4), (1, 2, 8))])
def test_seam_by_seam_finite_prefix_is_bounded_by_complete_law(first, second):
    bins, tail = r.count_law(first, second, count=2, maximum=8)
    pending, stopped = {(0, 0): F(1)}, {}
    for _ in range(4):
        next_pending = {}
        for (a_count, b_count), probability in pending.items():
            for seam in r.e._source().seams:
                a = a_count+any(port in first for port in seam)
                b = b_count+any(port in second for port in seam)
                target, key = (stopped, b) if a == 2 else (next_pending, (a, b))
                target[key] = target.get(key, F(0))+probability/30
        pending = next_pending
    remaining = sum(pending.values(), F(0))
    assert sum(stopped.values(), F(0))+remaining == 1
    assert sum(bins)+tail == 1
    for index, probability in enumerate(bins):
        prefix = stopped.get(index, F(0))
        assert prefix <= probability <= prefix+remaining


def test_complete_clock_field_readout_matches_independent_native_branch_oracle(receipt):
    comparisons = receipt['distinct_clock_witness']['joint_field_comparisons']
    assert len(comparisons) == 3
    final = comparisons[-1]
    assert final['A_ticks'] == 2 and final['B_count_PGF_argument'] == '1'
    weights = {(x['initial_A_port_bit'], x['final_B_port_bit']): F(x['weight'])
               for x in final['joint_bit_generating_weights']}
    assert weights[1, 1] == F(101, 5046)
    assert sum(weights.values()) == 1
    assert weights[1, 1]+weights[1, 0] == weights[1, 1]+weights[0, 1] == F(1, 6)
    assert receipt['scope']['collective_union_access_is_same_as_one_primitive_port_tape'] is False


@pytest.mark.parametrize('region', [(), [], (True,), (1, 0), (1, 1), (-1,), (12,), (0.0,)])
def test_region_contract_rejects_ambiguous_or_non_native_masks(region):
    with pytest.raises(ValueError):
        r.event_mask(region)


def test_invalid_counts_and_nonlocal_field_readouts_are_rejected():
    for counts in ((True, 5, 1), (0, 5, 0), (30, 30, 29), (5, 5, -1)):
        with pytest.raises(ValueError):
            r.parameters(*counts)
    with pytest.raises(ValueError):
        r.joint_field((0,), (1,), 2, 1)
    for z in (F(-1), F(2), 0.5, True):
        with pytest.raises(ValueError):
            r.pgf((0,), (1,), z)
