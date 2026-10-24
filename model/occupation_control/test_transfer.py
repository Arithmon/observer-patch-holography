"""All-endpoint fixed-apparatus comparisons, provenance and asymmetric dynamics."""
from copy import deepcopy
from fractions import Fraction as F
import gzip
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.anharmonic_observations import reference
from . import transfer as t
from .experiment import run, PARAMETERS, CASES


@pytest.mark.parametrize('reference,approximate', [
    ([F(1, 2), F(1, 2), F(0)], [F(3, 5), F(2, 5), F(0)]),
    ([F(0), F(1), F(0)], [F(1, 5), F(1), F(0)]),
    ([F(3, 10), F(3, 10), F(2, 5)], [F(7, 20), F(7, 10), -F(1, 20)]),
])
def test_trace_aware_bound_covers_signed_unnormalized_paths(reference, approximate):
    difference = [a-r for a, r in zip(approximate, reference, strict=True)]
    epsilon = sum(map(abs, difference))/2
    trace_defect = abs(sum(difference))
    q = approximate[0]+approximate[1]
    actual = abs(approximate[0]/q-reference[0]/(reference[0]+reference[1]))
    assert actual <= t.conditional_error(q, epsilon, trace_defect)


def test_trace_defect_term_is_necessary_and_both_spectral_bounds_are_sharp():
    epsilon = F(1, 10)
    assert t.conditional_error(F(1), epsilon, F(0)) == epsilon
    q = 1+2*epsilon
    # Reference (0,1), approximate (2 epsilon,1). A naive traceless
    # improvement would understate its conditional error by a factor two.
    measured_error = 2*epsilon/q
    assert t.conditional_error(q, epsilon, 2*epsilon) == measured_error
    assert measured_error > epsilon/q
    assert t.conditional_error(epsilon, epsilon, F(0)) is None


def test_probability_radius_does_not_assume_a_positive_approximate_state():
    # A coarse supplied scalar-evaluation error can move the returned number
    # outside the probability interval. The theorem must not clip its radius.
    row = {'mass_drift': 0., 'half_trace_error_upper': .1,
           'acceptance_exact_binary_sum': '1',
           'selected_probability_evaluation_error_upper': 2.,
           'selected_mean_evaluation_error_upper': 0.,
           'estimator_evaluation_error_upper': 0.}
    assert t.trace_aware_readout(row)['selected_probability_error_upper'] > 2


def test_trace_aware_bound_handles_noncommuting_signed_approximation():
    rng = np.random.default_rng(1837)
    qop = np.diag([1., 1., 0., 0.])
    detector = np.array([1., 1j, 0., 0.])/math.sqrt(2)
    effect = np.outer(detector, detector.conj())
    negative_paths = 0
    for _ in range(20):
        vector = rng.normal(size=4)+1j*rng.normal(size=4)
        vector /= np.linalg.norm(vector)
        physical = np.outer(vector, vector.conj())
        raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        error = .002*(raw+raw.conj().T)
        approximate = physical+error
        negative_paths += np.linalg.eigvalsh(approximate).min() < 0
        q = np.trace(qop@approximate).real
        epsilon = np.linalg.norm(error, 'nuc')/2+1e-12
        trace_defect = abs(np.trace(error).real)+1e-12
        radius = t.conditional_error(F(float(q)), F(float(epsilon)), F(float(trace_defect)))
        observed = np.trace(effect@approximate).real/q
        truth = np.trace(effect@physical).real/np.trace(qop@physical).real
        assert radius is not None and abs(observed-truth) <= float(radius)
    assert negative_paths > 0


@pytest.mark.parametrize('arguments', [(True, F(0), F(0)), (F(0), F(0), F(0)),
                                      (F(1), -F(1), F(0)), (F(1), F(0), -F(1)),
                                      (F(1), .1, F(0))])
def test_conditional_bound_requires_exact_valid_inputs(arguments):
    with pytest.raises(ValueError):
        t.conditional_error(*arguments)


@pytest.fixture(scope='module')
def retained():
    return t.read_receipt(t.FOLDER/'transfer_results.json.gz')


def _seal(prediction):
    prediction['seal'] = t.digest({k: v for k, v in prediction.items() if k != 'seal'})


def test_complete_all_case_inventory_and_uncertainty_scope(retained):
    comparison = retained['comparison']
    compared = comparison['runs'][0]
    assert compared['held_out_count'] == 125
    assert {(row['case'], F(row['time_ms'])) for row in compared['held_out_rows']} == {
        (case, F(step, 5)) for case in 'abcde' for step in range(2, 51, 2)}
    assert len(compared['zero_time_diagnostics_not_held_out']) == 5
    assert compared['certificate_endpoints_without_observations'] == []
    assert comparison['experimental_uncertainty_applied'] is None
    assert comparison['full_uncertainty_agreement_established'] is False
    assert comparison['statistical_rejection_claimed'] is False
    assert comparison['accepted_counts_or_raw_trials_authenticated'] is False
    assert retained['prediction']['settings']['new_parameter_fit'] is False
    assert [row['case'] for row in comparison['per_case']] == list('abcde')
    assert all(row['observations'] == 25 for row in comparison['per_case'])
    for result in retained['prediction']['runs']:
        for row in result['records']:
            old, new = row['estimator_error_upper'], row['trace_aware_readout']['estimator_error_upper']
            if old is not None:
                assert new is not None and new <= old
    for row in compared['held_out_rows']:
        if row['status'] == 'certified_interval':
            lower, upper = map(F, row['certified_interval'])
            observed = F(row['observed_decimal'])
            assert F(row['minimum_uncovered_absolute_allowance']) == max(F(0), lower-observed, observed-upper)


def test_every_endpoint_agrees_with_the_independently_executed_finite_source(retained):
    source = json.loads((t.ROOT/'model/anharmonic_transfer/results.json').read_text())
    native = next(r for r in source['prediction']['runs'] if r['name'] == 'occupation_refinement')
    assert native['dimension'] == 64
    rows = {(r['case'], r['steps']): r for r in native['predictions']+native['zero_time_predictions']}
    checked = 0
    for result in retained['prediction']['runs']:
        for row in result['records']:
            target = rows[(result['case'], row['steps'])]
            finite_error = row['exact_finite_path_half_trace_error_upper']
            acceptance = float(F(row['acceptance_exact_binary_sum']))
            assert acceptance > 2*finite_error
            semantic_mean_bound = 4*finite_error/acceptance+row['selected_mean_evaluation_error_upper']
            native_mean_bound = target['ideal_finite_model_error']['selected_readout']['conditional_spin_mean_error_upper']
            delta = abs(row['selected_mean']-target['predicted_spin_mean'])
            # Combine the independently retained finite-law enclosures. The
            # fixed native replay allowance remains only a numerical screen.
            assert delta <= semantic_mean_bound+native_mean_bound+2e-10
            checked += 1
            assert row['drive_elapsed_seconds'] == pytest.approx(target['drive_elapsed_seconds'], abs=1e-17)
    assert checked == 130
    # This is an independent finite implementation diagnostic. It supplies
    # neither a compiler error certificate nor physical parameter accuracy.


def test_prediction_does_not_open_empirical_values(retained, monkeypatch):
    by_case = {r['case']: r for r in retained['prediction']['runs']}
    def forbidden():
        raise AssertionError('empirical values opened during prediction')
    monkeypatch.setattr(t.data, 'load_validation', forbidden)
    monkeypatch.setattr(t.data, 'load_zero_tokens', forbidden)
    calls = []
    def cached_run(**kwargs):
        calls.append(kwargs)
        return deepcopy(by_case[kwargs['case']])
    monkeypatch.setattr(t, 'run', cached_run)
    assert t.predict() == retained['prediction']
    assert [r['case'] for r in calls] == list('abcde')
    assert all(r['record_steps'] == list(t.STEPS) for r in calls)


def test_poisoned_outcome_changes_score_but_no_prediction(retained, monkeypatch):
    observations = t.data.load_validation()
    observations[0]['observed_position_estimator'] = 100.
    observations[0]['observed_position_estimator_decimal'] = '100'
    monkeypatch.setattr(t.data, 'load_validation', lambda: observations)
    altered = t.evaluate(retained['prediction'])
    assert altered['prediction'] == retained['prediction']
    assert altered['comparison'] != retained['comparison']
    assert F(altered['comparison']['per_case'][0]['minimum_common_remaining_allowance']) > 99


@pytest.mark.parametrize('mutation', ['missing_case', 'missing_endpoint', 'duplicate_endpoint',
                                    'case_phase', 'clock', 'radius', 'bool_radius', 'scope', 'trace_radius'])
def test_resealed_scope_or_coverage_forgery_fails(retained, mutation):
    prediction = deepcopy(retained['prediction'])
    result = prediction['runs'][1]
    row = result['records'][1]
    if mutation == 'missing_case':
        prediction['runs'].pop()
    elif mutation == 'missing_endpoint':
        result['records'].pop()
    elif mutation == 'duplicate_endpoint':
        result['records'][1] = deepcopy(result['records'][0])
    elif mutation == 'case_phase':
        result['parameters']['phi'] *= -1
    elif mutation == 'clock':
        row['drive_elapsed_seconds'] = row['published_coordinate_ms']/1000
    elif mutation == 'radius':
        row['estimator_error_upper'] += .1
    elif mutation == 'bool_radius':
        row['estimator_error_upper'] = False
    elif mutation == 'trace_radius':
        row['trace_aware_readout']['estimator_error_upper'] = 0.
    else:
        result['certificate']['apparatus_or_calibration_uncertainty_included'] = True
    _seal(prediction)
    with pytest.raises(ValueError):
        t.evaluate(prediction)


def test_unsealed_prediction_fails_before_measurements(retained, monkeypatch):
    prediction = deepcopy(retained['prediction'])
    prediction['runs'][0]['records'][1]['selected_estimator'] += .1
    def forbidden():
        raise AssertionError('measurements opened before seal validation')
    monkeypatch.setattr(t.data, 'load_validation', forbidden)
    with pytest.raises(ValueError, match='seal'):
        t.evaluate(prediction)


@pytest.mark.parametrize('mutation', ['score', 'summary', 'provenance', 'agreement'])
def test_receipt_cannot_rewrite_derived_empirical_consequences(retained, mutation):
    changed = deepcopy(retained)
    if mutation == 'score':
        changed['comparison']['runs'][0]['held_out_rows'][0]['minimum_uncovered_absolute_allowance'] = '0'
    elif mutation == 'summary':
        changed['comparison']['per_case'][0]['certified'] = 0
    elif mutation == 'provenance':
        changed['data_provenance']['accepted_counts_included'] = True
    else:
        changed['comparison']['full_uncertainty_agreement_established'] = True
    with pytest.raises(ValueError):
        t.replay(retained, changed)


def test_compressed_receipt_roundtrip_is_exclusive_and_rejects_duplicates(retained, tmp_path):
    path = tmp_path/'receipt.json.gz'
    t.write_receipt(path, retained)
    assert t.read_receipt(path) == retained
    with pytest.raises(FileExistsError):
        t.write_receipt(path, retained)
    malformed = tmp_path/'duplicate.json.gz'
    malformed.write_bytes(gzip.compress(b'{"prediction": {}, "prediction": {}}', mtime=0))
    with pytest.raises(ValueError, match='duplicate'):
        t.read_receipt(malformed)


@pytest.mark.parametrize('case', list('bcde'))
def test_asymmetric_finite_carrier_paths_match_independent_unsplit_law(case):
    """Direct dense master equation, including the nonzero carrier dwell.

    This checks signs, preparation sides and both clock conventions using a
    different finite solver. Its numerical agreement does not prove the
    unrestricted bound; the existing full-defect theorem supplies that bound.
    """
    dimension, p = 12, PARAMETERS
    phase, side = CASES[case]
    alpha = {'left': p['alpha_real'], 'right': -p['alpha_real'], 'center': 0.}[side]
    state = reference.preparation(dimension, alpha, p['thermal'])
    _, _, sx = reference.operators(dimension)
    sz = np.diag([1., -1.])
    gap = abs(phase)/p['omega_carrier']
    common = dict(delta=p['delta'], omega=p['omega'], dephasing_rate=p['gamma'])
    for _ in range(2):
        for sign, multiple, angles in ((-1, 1, (-phase, phase)),
                                      (1, 2, (phase, -phase)), (-1, 1, None)):
            state = reference.lindblad_pulse(state, duration=multiple*p['pulse_duration'], sign=sign, **common)
            if angles is not None:
                for index, angle in enumerate(angles):
                    state = reference.lindblad_pulse(state, delta=p['delta'], omega=0.,
                                                      duration=gap, sign=1, dephasing_rate=p['gamma'])
                    unitary = np.kron(expm(-.5j*angle*sx), np.eye(dimension))
                    state = unitary@state@unitary.conj().T
                    if index == 0:
                        unitary = np.kron(expm(-.5j*p['theta']*sz), np.eye(dimension))
                        state = unitary@state@unitary.conj().T
    target = reference.readout(state, [{'u': 0., 'v': p['probe_h'], 'theta': -math.pi/2}])['reads'][0]
    actual = run(dimension=dimension, case=case, steps=2)['records'][-1]
    finite = actual['exact_finite_path_half_trace_error_upper']
    assert np.max(abs(np.array(actual['spin_detector_joint_probabilities'])-
                      np.array(target['spin_detector_joint_probabilities']))) <= 2*finite+2e-11
    assert abs(actual['selected_mean']-target['selected_mean']) < 5e-10
    assert actual['drive_elapsed_seconds'] == 2*4*(p['pulse_duration']+gap)
    assert actual['drive_elapsed_seconds'] > actual['published_coordinate_ms']/1000
