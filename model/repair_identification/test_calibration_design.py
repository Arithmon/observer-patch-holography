"""Fixed sampling design checks and an exact optional-stopping countercontrol."""
from fractions import Fraction as F
from math import isqrt

import pytest

from . import calibration as cal
from . import core


def _evidence(attempted=100):
    return {
        'schema': cal.SCHEMA,
        'settings': [
            {'probe': probe.name, 'attempted': attempted,
             'counts': {label: attempted if label == 'abort' else 0
                        for label in cal.OUTCOMES}}
            for probe in core.PROBES],
    }


def _analyze(design, evidence=None):
    table = {name: 0 for name in ('k', 'g', 'h', 'm', 'n')}
    return cal.analyze_counts(_evidence() if evidence is None else evidence,
                              design, {'interior': table}, 2, (F(0),)*9)


def test_design_consistency_does_not_claim_verified_precommitment():
    design = cal.sampling_design((100,)*9, (F(1, 2),)*9)
    result = _analyze(design)['count_evidence']
    assert result['simultaneous_sampling_failure_bound'] == '9/100'
    assert result['declared_sampling_design'] == {
        'planned_attempts': [100]*9, 'sampling_radii': ['1/2']*9}
    assert result['sampling_design_consistent_with_counts'] is True
    assert result['sampling_design_prior_declaration_verified'] is False
    assert result['premises_established_by_counts'] is False
    assert any('before, and independently of, all analyzed trial outcomes' in premise
               for premise in result['premises'])
    assert result['abort_counts'] == [100]*9


def test_attempts_must_match_the_separate_design():
    design = cal.sampling_design((101,)*9, (F(1, 2),)*9)
    with pytest.raises(ValueError, match='match the declared sampling design'):
        _analyze(design)
    with pytest.raises(ValueError, match='SamplingDesign'):
        _analyze((F(1, 2),)*9)  # The former outcome-derived-radii API is rejected.


@pytest.mark.parametrize('attempts,radii', [
    ([100]*9, (F(1, 2),)*9),
    ((100,)*8, (F(1, 2),)*9),
    ((True,)+(100,)*8, (F(1, 2),)*9),
    ((0,)+(100,)*8, (F(1, 2),)*9),
    ((100,)*9, [F(1, 2)]*9),
    ((100,)*9, (F(1, 2),)*8),
    ((100,)*9, (True,)+(F(1, 2),)*8),
    ((100,)*9, (F(0),)+(F(1, 2),)*8),
    ((100,)*9, (F(2),)+(F(1, 2),)*8),
])
def test_design_factory_rejects_aliases_and_invalid_values(attempts, radii):
    with pytest.raises(ValueError):
        cal.sampling_design(attempts, radii)


@pytest.mark.parametrize('design', [
    cal.SamplingDesign((True,)+(100,)*8, ((1, 2),)*9),
    cal.SamplingDesign((100,)*9, ((True, 2),)+((1, 2),)*8),
    cal.SamplingDesign((100,)*9, ((1, True),)+((1, 2),)*8),
    cal.SamplingDesign((100,)*9, ((2, 4),)+((1, 2),)*8),
    cal.SamplingDesign((100,)*9, ((1, 0),)+((1, 2),)*8),
    cal.SamplingDesign((100,)*9, ((-1, -2),)+((1, 2),)*8),
    cal.SamplingDesign((100,)*9, ((3, 2),)+((1, 2),)*8),
    cal.SamplingDesign((100,)*9, ([1, 2],)+((1, 2),)*8),
    cal.SamplingDesign((100,)*9, ((1, 2),)*8),
])
def test_direct_namedtuple_construction_cannot_bypass_validation(design):
    with pytest.raises(ValueError):
        _analyze(design)


def test_design_detaches_fraction_storage_and_is_immutable():
    radius = F(1, 2)
    design = cal.sampling_design((100,)*9, (radius,)*9)
    object.__setattr__(radius, '_numerator', 2)
    assert design.sampling_radii == ((1, 2),)*9
    assert _analyze(design)['count_evidence']['sampling_radii'] == ['1/2']*9
    assert not hasattr(design, '__dict__')
    with pytest.raises(AttributeError):
        design.planned_attempts = (1,)*9


def test_iid_optional_stopping_can_exceed_the_reported_fixed_design_bound():
    """Exact lower bound, not a stochastic simulation or an asymptotic claim.

    For IID fair +/-1 trials, stop at the first n>=40 with
    |S_n|>ceil(2 sqrt(n)), or at n=1000. At a hit choose radius
    ceil(2 sqrt(n))/(2n); otherwise choose radius one. Every realized
    naive Chebyshev bound is <=1/4, but actual exclusion exceeds 0.3198.

    Integer entries lower-bound surviving probabilities on a 2^-52 grid.
    Combining two parents then flooring their half preserves this order.
    Killed entries are disjoint stopping events; their accumulated mass is
    therefore a rigorous lower bound. Nothing is conditioned or normalized.
    """
    scale = 1 << 52
    alive, origin, hit = [scale], 0, 0
    for n in range(1, 1001):
        new = [0]*(len(alive)+2)
        for index, mass in enumerate(alive):
            new[index] += mass
            new[index+2] += mass
        new = [mass//2 for mass in new]
        origin += 1
        boundary = isqrt(4*n-1)+1
        assert boundary*boundary >= 4*n
        if n >= 40:
            assert F(n, boundary*boundary) <= F(1, 4)
            lo = max(0, origin-boundary)
            hi = min(len(new), origin+boundary+1)
            hit += sum(new[:lo])+sum(new[hi:])
            new = new[lo:hi]
            origin -= lo
        alive = new
        assert hit+sum(alive) <= scale
    exclusion_lower = F(hit, scale)
    assert exclusion_lower == F(360112670136069, 1125899906842624)
    # Eight other fixed settings with m=10^6, radius=1 add only 2e-6.
    naive_nine_setting_bound = F(1, 4)+8*F(1, 4_000_000)
    assert naive_nine_setting_bound == F(125001, 500000)
    assert exclusion_lower > naive_nine_setting_bound
