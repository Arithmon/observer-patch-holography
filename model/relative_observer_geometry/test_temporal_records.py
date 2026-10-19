"""Independent checks of correlation-only temporal reconstruction and scope."""
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest
import sympy as sp

from . import temporal_records as t


@pytest.fixture(scope='module')
def native_correlations():
    return t.source_correlations(20)


@pytest.fixture(scope='module')
def representation(native_correlations):
    return t.recover(tuple(row[0][0]-F(1, 36) for row in native_correlations[:15]))


def test_recovery_consumes_only_record_correlations_and_predicts_unseen_lags(monkeypatch):
    poles, weights = (F(1, 4), F(2, 3)), (F(2, 7), F(3, 11))
    moments = tuple(sum(w*p**k for w, p in zip(weights, poles)) for k in range(5))
    def forbidden():
        raise AssertionError('a hidden native source was accessed by the estimator')
    monkeypatch.setattr(t.core, 'source', forbidden)
    result = t.recover(moments)
    assert result.dimension == 2
    assert result.polynomial == (F(1, 6), F(-11, 12), F(1))
    for lag in (0, 1, 9, 27):
        assert t.predict(result, lag) == sum(w*p**lag for w, p in zip(weights, poles))
    gram, transfer, shifted = map(sp.Matrix, (result.gram, result.transfer, result.shifted_gram))
    assert gram*transfer == shifted == transfer.T*gram


@pytest.mark.parametrize('moments', (
    (F(0), F(0), F(0)), (F(1), F(2), F(1)),
    (F(1), F(1), F(1)), (F(1), F(-1, 2), F(1, 4)),
    (F(1), F(1, 2), F(1, 3)), (F(1), F(1, 2), F(1, 4), F(1, 8), F(1, 17)),
    (F(1), 0.5, F(1, 4)), [F(1), F(1, 2), F(1, 4)]))
def test_invalid_indefinite_noncontractive_or_underdetermined_moments_are_rejected(moments):
    with pytest.raises(ValueError):
        t.recover(moments)


def test_exact_native_law_recovers_seven_visible_modes_and_source_held_out_lags(
        representation, native_correlations):
    assert representation.dimension == 7
    x = sp.Symbol('x')
    actual = sum(sp.Rational(value)*x**i for i, value in enumerate(representation.polynomial))
    expected = (753*x**3-1681*x**2+1241*x-303)*(18147*x**4-54265*x**3+60655*x**2-30039*x+5562)
    assert sp.expand(actual-expected/sp.LC(sp.Poly(expected, x))) == 0
    assert native_correlations[0] == ((F(1, 6), F(1, 66)), (F(1, 66), F(1, 6)))
    assert native_correlations[1][0][0] == F(37, 297)
    for lag in range(15, 21):
        assert t.predict(representation, lag) == native_correlations[lag][0][0]-F(1, 36)


def test_symmetric_and_antisymmetric_records_split_the_visible_modes(native_correlations):
    plus = tuple(sum(row[i][j]-F(1, 36) for i in (0, 1) for j in (0, 1))
                 for row in native_correlations)
    minus = tuple(sum((1 if i == j else -1)*(row[i][j]-F(1, 36))
                     for i in (0, 1) for j in (0, 1)) for row in native_correlations)
    assert t.recover(plus[:9]).dimension == 4
    assert t.recover(minus[:7]).dimension == 3
    for row in native_correlations:
        # <(n0+n1-1/3)(0) (n0-n1)(k)> vanishes; this is observable
        # parity, not an asserted identification with the spatial slow band.
        assert row[0][0]-row[0][1]+row[1][0]-row[1][1] == 0


def test_rational_spectral_intervals_match_independent_native_mean_eigenvalues(representation):
    intervals = t.spectrum(representation)
    roots = np.linalg.eigvals(np.asarray(t._mean_kernel(), dtype=float)).real
    for row in intervals:
        low, high = map(F, row['transfer_interval'])
        candidates = [value for value in roots if float(low)-1e-13 <= value <= float(high)+1e-13]
        assert len(candidates) == 1
        energy_low, energy_high = map(F, row['energy_interval_per_incident_tick'])
        assert float(energy_low)-2e-13 <= -math.log(candidates[0]) <= float(energy_high)+2e-13
    assert len(intervals) == 7


def test_log_bound_is_positive_directed_and_convergent():
    for x in (F(1, 2), F(3, 4), F(9, 10)):
        low, high = t._negative_log_bounds(x, terms=5)
        fine_low, fine_high = t._negative_log_bounds(x, terms=12)
        assert low < fine_low < fine_high < high
        assert float(low) <= -math.log(float(x)) <= float(high)


def test_small_errors_have_a_bound_but_finite_shots_do_not_certify_exact_spectroscopy(representation):
    assert F(5, 10**23) < representation.gram_lower < F(6, 10**23)
    assert t.perturbation_budget(representation, F(0)) == 0
    assert 0 < t.perturbation_budget(representation, F(1, 10**26)) < F(3, 1000)
    assert t.perturbation_budget(representation, F(1, 16)) is None
    assert t.perturbation_budget(representation, representation.gram_lower/7) is None


def test_perturbation_inequality_against_independent_generalized_eigensolve():
    from scipy.linalg import eigvalsh
    moments = tuple(F(1, 3)*F(1, 4)**k+F(2, 3)*F(3, 4)**k for k in range(5))
    result = t.recover(moments)
    gram, shifted = np.array(result.gram, float), np.array(result.shifted_gram, float)
    tolerance = F(1, 10**5)
    bound = float(t.perturbation_budget(result, tolerance))
    for sign in (-1, 1):
        perturb = np.array([[1., -1.], [-1., 1.]])*float(tolerance)*sign
        altered = eigvalsh(shifted+perturb, gram-perturb)
        assert max(abs(altered-np.array([0.25, 0.75]))) <= bound


def test_sampling_bound_has_explicit_conditional_coverage():
    upper = t.sampling_bound(2048, F(1, 16))
    assert 8*math.exp(-16) <= float(upper) < 1/10**6
    assert t.sampling_bound(1, F(1, 16)) == 1
    with pytest.raises(ValueError):
        t.sampling_bound(True, F(1, 16))


def test_accessible_current_counter_executes_declared_finite_lag_predictions(native_correlations):
    rows = t.sampled_records()
    assert [row['joint_ones_count'] for row in rows] == [267, 201, 172, 155]
    assert [row['native_attempts'] for row in rows] == [6956, 13793, 20757, 27571]
    for row in rows:
        exact = native_correlations[row['lag_in_incident_ticks']][0][0]
        assert abs(F(row['joint_ones_count'], row['samples'])-exact) <= F(1, 16)


def test_source_screen_gram_certificate_is_positive_after_all_exported_error_allowances():
    result = t.screen_compression()
    intervals = [tuple(map(F, row)) for row in result['centered_root_field_covariance_intervals']]
    error = max((hi-lo)/2 for lo, hi in intervals)
    centers = [(lo+hi)/2 for lo, hi in intervals]
    assert error == F(result['uniform_covariance_midpoint_error'])
    gram = sp.Matrix(2, 2, lambda i, j: centers[i+j])
    shifted = sp.Matrix(2, 2, lambda i, j: centers[i+j+1])
    for matrix, allowance, field in ((gram, 2*error, 'gram_positive_lower'),
            (shifted, 2*error, 'shifted_gram_positive_lower'),
            (gram-shifted, 4*error, 'gram_minus_shifted_positive_lower')):
        assert t._rational(1/sp.trace(matrix.inv()))-allowance == F(result[field]) > 0
    assert F(result['generalized_eigenvalue_error_upper']) < F(2, 10**5)
    assert result['incident_seams'] == 33
    assert result['compressed_transfer_is_closed_or_equals_isolated_carrier_spectrum'] is False


def test_mutated_reconstruction_cannot_be_used_as_a_prediction_or_spectral_certificate(representation):
    forged = representation._replace(dimension=8)
    with pytest.raises(ValueError):
        t.predict(forged, 20)
    with pytest.raises(ValueError):
        t.spectrum(forged)
    with pytest.raises(ValueError):
        t.perturbation_budget(forged, F(0))


@pytest.mark.parametrize('level', (0, 1, 2))
def test_refined_clock_rates_use_actual_seams_and_ancestral_regions(level):
    census = t.source_clock_incidence(level)
    assert census['port_degree_census'] == {'5': 180*4**level, '6': 60*4**level}
    assert census['one_carrier_incident_seams'] == 33
    assert census['base_face_ports'] == 12*4**level
    assert census['base_face_incident_seams'] == (63*4**level+3*2**level)//2
    assert census['base_face_boundary_seams'] == 3*2**level


def test_uniform_clock_bound_encloses_exact_finite_prefix_failure_law():
    # Enumerate the complete count-state distribution and an absorbing failure
    # bit. This is a separate rational calculation of a maximum-over-time law.
    result = t.incident_clock_bound(0, 105, F(1, 126), F(1, 63))
    a, p, delta = 1260, F(1, 6), F(1, 63)
    pending = {(0, False): F(1)}
    for n in range(1, result['covered_native_attempts']+1):
        updated = {}
        for (count, failed), mass in pending.items():
            for tick, weight in ((0, 1-p), (1, p)):
                after = count+tick
                key = after, failed or abs(after-n*p) > a*p*delta
                updated[key] = updated.get(key, F(0))+mass*weight
        pending = updated
    actual = sum(mass for (_, failed), mass in pending.items() if failed)
    assert sum(pending.values()) == 1
    assert 0 < actual <= F(result['failure_probability_upper']) < 1
    assert result['covered_native_attempts'] == 30


@pytest.mark.parametrize('p', (F(1, 4), F(1, 2), F(3, 4), F(1)))
def test_inverse_clock_inclusion_on_every_admitted_finite_tick_word(p):
    # Check the deterministic event inclusion behind the probability bound.
    # Endpoint/left-limit checks cover the entire piecewise-constant clock.
    a, horizon, delta = 8, F(1, 2), F(1, 2)
    cap = math.ceil(a*(horizon+delta))
    checked = 0
    for word in product((0, 1), repeat=cap):
        counts = [0]
        for tick in word:
            counts.append(counts[-1]+tick)
        if max(abs(counts[n]-n*p) for n in range(cap+1)) > a*p*delta:
            continue
        checked += 1
        targets = math.floor(a*p*horizon)
        assert counts[-1] >= targets
        for k in range(targets+1):
            inverse = F(next(n for n, count in enumerate(counts) if count >= k), a)
            left, right = F(k, a*p), min(horizon, F(k+1, a*p))
            assert abs(inverse-left) <= delta+1/(a*p)
            assert abs(inverse-right) <= delta+1/(a*p)
        for n in range(math.floor(a*horizon)+1):
            displayed = F(counts[n], a*p)
            assert abs(displayed-F(n, a)) <= delta+F(1, a)
            assert abs(displayed-min(horizon, F(n+1, a))) <= delta+F(1, a)
    assert checked > 0


def test_correlated_clocks_share_union_bound_without_independent_replica_assumption():
    source = t.core.source()
    p = F(1, 6)
    for other in (0, 1, 2):
        common = F(sum(0 in seam and other in seam for seam in source.seams), 30)
        covariance = common-p*p
        assert covariance == (F(5, 36) if other == 0 else
                              F(1, 180) if (0, other) in source.seams else F(-1, 36))
    result = t.native_clock_field_dependence()
    assert result['native_branches'] == 3960
    assert F(result['clock_field_covariance']) == -F(25, 2376)
    assert F(result['clock_single_endpoint_covariance']) == 0
    assert result['clock_field_independence_asserted'] is False


@pytest.mark.parametrize('level', (0, 1, 2))
def test_exact_clock_fluctuation_covariance_from_complete_native_seam_distribution(level):
    graph = t.clocks.screen.source.federation(level)
    m = len(graph.seams)
    rows = t.source_clock_fluctuations(level)
    assert [row['incidence_counts'] for row in rows] == [[5, 5, 1], [5, 5, 0], [6, 6, 1]]
    assert [F(row['limiting_correlation_squared_for_fixed_incidence_counts']) for row in rows] == [F(1, 25), F(0), F(1, 36)]
    for row in rows:
        ports = row['ports']
        marks = [tuple(int(p in seam) for p in ports) for seam in graph.seams]
        means = [F(sum(mark[i] for mark in marks), m) for i in (0, 1)]
        direct = [[sum((mark[i]-means[i])*(mark[j]-means[j]) for mark in marks)
                   for j in (0, 1)] for i in (0, 1)]
        reported = [[F(x) for x in values] for values in row['scaled_count_covariance_per_continuum_unit']]
        assert direct == reported
        assert sp.Matrix(reported).is_positive_definite
    # Distinct clocks on disjoint seam sets are negatively correlated at
    # finite source size, despite their independent Gaussian limiting counts.
    assert F(rows[1]['scaled_count_covariance_per_continuum_unit'][0][1]) < 0


def test_finite_count_covariance_retains_aliases_and_deterministic_complete_clock():
    aliases = t.clock_fluctuation_covariance(0, 5, 5, 5)
    assert len({x for row in aliases for x in row}) == 1
    complete = t.clock_fluctuation_covariance(0, 630, 5, 5)
    assert complete[0] == (F(0), F(0)) and complete[1][0] == 0
    assert complete[1][1] > 0
    for args in ((0, 5, 5, True), (0, 630, 5, 4), (0, 5, 5, 6), (0, 5, 5, -1)):
        with pytest.raises(ValueError):
            t.clock_fluctuation_covariance(*args)


def test_every_nonempty_event_set_has_many_ticks_on_density_continuum_scale():
    for level in (0, 1, 2, 8, 100):
        for d in (1, 5, 6, 33, 630*4**level):
            result = t.incident_clock_bound(level, d)
            assert F(result['mean_incident_ticks_per_continuum_unit']) == 2*d*4**level
            if d == 630*4**level:
                assert F(result['failure_probability_upper']) == 0
        # A shrinking rational tolerance supplies an explicit uniform LLN.
        delta = F(1, 2**(level//2))
        bound = t.incident_clock_bound(level, 1, tolerance=delta)
        assert F(bound['failure_probability_upper']) <= min(F(1), F(2, 2**level))


@pytest.mark.parametrize('args', ((True, 5), (101, 5), (0, 0), (0, 631),
    (0, True), (0, 5, 1), (0, 5, F(0)), (0, 5, F(1), F(0))))
def test_clock_bound_rejects_unstated_or_empty_source_exposure(args):
    with pytest.raises(ValueError):
        t.incident_clock_bound(*args)
