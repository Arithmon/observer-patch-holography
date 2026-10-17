"""Independent checks of correlation-only temporal reconstruction and scope."""
from fractions import Fraction as F
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
