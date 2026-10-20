"""Complete native moments, coupled clocks and continuous-path controls."""
from fractions import Fraction as F
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp

from . import observation_process as o


def test_each_actual_native_seam_separately_preserves_the_full_binary_sector():
    states, targets, marks, score, count = o.marked_source()
    assert count == 3960 and len(states) == 66
    for seam in range(30):
        assert np.array_equal(np.bincount(targets[:, 2*seam:2*seam+2].ravel(), minlength=66),
                              np.full(66, 2))
    assert marks.sum() == 10
    assert np.all(score.sum(axis=0) == 0) and np.all(score.sum(axis=1) == 0)
    assert not targets.flags.writeable and not marks.flags.writeable and not score.flags.writeable


def test_full_native_two_step_law_matches_every_exact_joint_covariance_entry():
    actual = o.direct_two_steps()
    expected = o.finite_covariances(2)
    assert actual['branches'] == 237600
    assert set(actual['terminal_counts']) == {3600}
    for key, matrix in zip(('attempt_second_moment', 'difference_second_moment',
                            'cross_second_moment'), expected, strict=True):
        assert actual[key] == matrix
    assert actual['difference_lag_one'] == sp.zeros(12)
    assert actual['difference_mean'] == sp.zeros(12, 1)
    assert expected[2] != sp.zeros(12)  # The two clocks are not independent.


def test_one_event_and_joint_clt_covariance_have_distinct_quantifiers():
    a, r, cross = o.finite_covariances(1)
    assert a == o.occupation.CHI*o.occupation.P0
    assert r == 5*a and cross == sp.zeros(12)
    cross_limit, delta = o.joint_clt()
    assert cross_limit != cross and delta == r
    a = o.occupation.field('attempt')['covariance']
    joint = np.array(sp.BlockMatrix([[a, cross_limit], [cross_limit.T, delta]]), dtype=float)
    eigenvalues = np.linalg.eigvalsh(joint)
    assert eigenvalues.min() > -2e-13
    assert (eigenvalues > 1e-10).sum() == 22
    assert a+cross_limit+cross_limit.T+delta == 6*o.occupation.field('incident')['covariance']


def test_all_native_joint_innovations_are_centered_and_have_the_predicted_covariance():
    certificate = o.martingale_certificate()
    assert certificate['conditional_states'] == 66
    assert certificate['native_innovations'] == 3960
    assert certificate['all_conditional_means_zero']
    assert certificate['joint_asymptotic_covariance_checked']


def test_stationary_state_marginals_do_not_replace_independent_seam_scheduling():
    sticky = o.direct_two_steps(sticky=True)
    assert sticky['branches'] == 7920
    assert set(sticky['terminal_counts']) == {120}
    assert sticky['difference_mean'] == sp.zeros(12, 1)
    assert sp.trace(sticky['difference_lag_one']) == sp.Rational(250, 33)
    _, delta = o.joint_clt()
    lag = sticky['difference_lag_one']
    assert sticky['difference_second_moment'] == 2*delta+lag+lag.T
    assert sticky['difference_second_moment'] != 2*delta


def test_nonstationary_native_preparation_invalidates_the_reverse_mean_identity():
    biased = o.direct_two_steps(stationary=False)
    assert biased['branches'] == 3600
    assert biased['difference_mean'] != sp.zeros(12, 1)
    assert sum(biased['difference_mean']) == 0  # Total load remains conserved.
    assert biased['difference_second_moment'] != o.finite_covariances(2)[1]


def test_uncentered_clock_sums_do_not_admit_the_centered_random_time_change():
    # Constant readout1 at a two-attempt prefix: N~Binomial(2,1/6).
    probabilities = (F(25, 36), F(10, 36), F(1, 36))
    noise = tuple(6*k-2 for k in range(3))
    assert sum(p*x for p, x in zip(probabilities, noise, strict=True)) == 0
    assert sum(p*x*x for p, x in zip(probabilities, noise, strict=True)) == 10
    # At a prescribed tick count the uncentered sum is deterministic instead.
    assert len({sum([1]*2) for _ in range(3)}) == 1


@pytest.mark.parametrize('level', (0, 1))
def test_path_clock_constants_match_the_actual_source_seams(level):
    raw, _, _, incident, _, _ = o.continuum.uniform.ancestry.problem(level)
    seams = len(raw.seam_a)
    assert seams == 630*4**level and len(incident) == 33
    assert o.clock_difference(level) == (1-F(len(incident), seams))/(4*len(incident)*4**level)


def test_all_level_hilbert_path_bounds_vanish_with_refinement_and_cutoff():
    for level in (0, 3, 10, 100):
        report = o.path_bounds(level, horizon=F(3, 2), count=3, degree=25)
        d = o.clock_difference(level)
        assert report['finite_prefix_clock_sup_squared_upper'] == 30*d  # ceil(9/2)=5
        assert report['joint_Gaussian_clock_sup_squared_upper'] == 12*d
        assert report['joint_Gaussian_clock_sup_squared_upper'] < F(1, 11*4**level)
        assert report['Gaussian_harmonic_sup_tail_squared_upper'] == F(16389, 26**2)
    assert o.path_bounds(100)['joint_Gaussian_clock_sup_squared_upper'] < F(1, 10**60)
    assert o.path_bounds(0, degree=1000)['Gaussian_harmonic_sup_tail_squared_upper'] < F(11, 1000)


@pytest.mark.parametrize('times', ((F(1),), (F(1, 3), F(2), F(7, 3))))
def test_gaussian_cylinder_action_matches_independent_joint_inverse_and_increments(times):
    result = o.cylinder_action(times)
    covariance = sp.kronecker_product(result['time_covariance'], result['spatial_covariance'])
    inverse = covariance.inv(method='DM')
    assert inverse == sp.kronecker_product(result['time_precision'], result['spatial_precision'])
    difference = sp.eye(len(times))
    for i in range(1, len(times)):
        difference[i, i-1] = -1
    assert difference*result['time_covariance']*difference.T == sp.diag(*result['increments'])
    assert result['dimension'] == 11*len(times)


def test_native_finite_count_observations_are_not_independent_brownian_increments():
    native = o.occupation.field('attempt')['kernel']
    adjacent = o.occupation.CHI*(native*o.occupation.P0)[0, 0]
    assert adjacent == sp.Rational(25, 198) > 0
    assert o.finite_covariances(2)[0] != 2*o.finite_covariances(1)[0]


@pytest.mark.parametrize('bad', (True, -1, 1., '2', F(2)))
def test_integer_hypotheses_fail_closed(bad):
    with pytest.raises(ValueError): o.clock_difference(bad)
    with pytest.raises(ValueError): o.path_bounds(0, count=bad)
    with pytest.raises(ValueError): o.finite_covariances(bad)


@pytest.mark.parametrize('bad', (True, 1., F(), F(-1)))
def test_exact_horizon_hypothesis_fails_closed(bad):
    with pytest.raises(ValueError): o.path_bounds(0, horizon=bad)


def test_noncanonical_fractions_and_time_order_are_rejected():
    bad = F(1, 2); bad._numerator = 2; bad._denominator = 4
    with pytest.raises(ValueError): o.path_bounds(0, horizon=bad)
    for times in ((), (F(1), F(1)), (F(2), F(1)), [F(1)]):
        with pytest.raises(ValueError): o.cylinder_action(times)
    with pytest.raises(ValueError): o.path_bounds(0, count=0)
    with pytest.raises(ValueError): o.direct_two_steps(sticky=1)


def test_uniform_bounds_alone_do_not_force_a_unique_whole_sequence_covariance():
    # A control for the inference, not a substitute native source family.
    sequence = [sp.Rational(5+n % 2) for n in range(8)]
    assert min(sequence) > 4 and max(sequence) < o.continuum.uniform.GREEN_BOUND
    assert abs(sequence[-1]-sequence[-2]) == 1


def test_all_runtime_python_sources_are_pinned_in_a_fresh_process():
    code = '''
from pathlib import Path
import sys
from model.physical_space_time import observation_process as o
o.run()
missing=[]
for module in tuple(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if not path: continue
    try: name=Path(path).resolve().relative_to(o.ROOT).as_posix()
    except ValueError: continue
    if not name.startswith('.') and name.endswith('.py') and name not in o.SOURCE_PATHS: missing.append(name)
assert not missing,sorted(set(missing))
'''
    subprocess.run([sys.executable, '-c', code], cwd=o.ROOT, check=True, timeout=120)
