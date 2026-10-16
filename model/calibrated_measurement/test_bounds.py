"""Independent scalar enclosures, complete-interface and hostile controls."""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from copy import deepcopy

import numpy as np
import pytest

from . import bounds as b
from . import experiment as e


ZERO = dict.fromkeys(b.SETTINGS, F(0))
IDEAL = b.InterfaceBounds((F(0),)*6, (F(0),)*2, F(0))


def fixture(zr=F(4, 5), zi=F(0), wr=F(7, 25), wi=F(0), trials=20000):
    offset, gain = F(1, 20), F(17, 20)
    means = dict(zip(b.SETTINGS, (offset+gain, offset-gain, offset+gain*zr,
                    offset+gain*zi, offset+gain*wr, offset+gain*wi), strict=True))
    return e.count_fixture(means, trials)


@pytest.mark.parametrize('value', [F(1), F(2), F(3, 2), F(1200), F(10**9, 7)])
def test_exact_log_encloses_independent_high_precision_evaluation(value):
    lo, hi = b.log_interval(value)
    with localcontext() as context:
        context.prec = 150
        expected = (Decimal(value.numerator)/Decimal(value.denominator)).ln()
        assert Decimal(lo.numerator)/Decimal(lo.denominator) <= expected
        assert expected <= Decimal(hi.numerator)/Decimal(hi.denominator)
        assert hi-lo < F(1, 10**30)


@pytest.mark.parametrize('value', [F(0), F(4, 9), F(2), F(1, 10**100), F(10**40+1)])
def test_root_bounds_are_exact_for_tiny_and_large_rationals(value):
    lo, hi = b.root_interval(value)
    assert lo*lo <= value <= hi*hi
    assert hi-lo <= F(1, b.GRID)


def test_hoeffding_union_bound_uses_outward_radius_and_all_six_allocations():
    for n in (200, 2000, 20000):
        r = b.sampling_radius(n)
        with localcontext() as context:
            context.prec = 100
            exponent = -Decimal(n)*Decimal(r.numerator)**2/Decimal(r.denominator)**2/2
            assert 12*exponent.exp() <= Decimal('.01')


def test_canonical_linear_instrument_completion_bound_for_noncommuting_effects():
    # K0,K1 are individually nonnormal; T is positive and not a scalar.
    k0 = np.array([[.6, .03j], [.02, .55]])
    k1 = np.array([[.78, -.01], [.04j, .82]])
    effects = [k.conj().T@k for k in (k0, k1)]
    total, observable = effects[0]+effects[1], effects[0]-effects[1]
    values, vectors = np.linalg.eigh(total)
    inverse = (vectors/np.sqrt(values))@vectors.conj().T
    corrected = [inverse@effect@inverse for effect in effects]
    np.testing.assert_allclose(sum(corrected), np.eye(2), atol=1e-14)
    c = F(1, 5)
    assert np.linalg.norm(total-np.eye(2), 2) < float(c)
    bound = b.InterfaceBounds((c,)*6, (F(0),)*2, F(0)).completion('h_real')
    assert np.linalg.norm(corrected[0]-corrected[1]-observable, 2) < float(bound)


def test_boundary_state_and_sampling_uncertainty_do_not_cause_false_rejection():
    # U=diag(.8+.6i,.8-.6i), rho=I/2 exactly saturates the moment disk.
    for n in (200, 2000, 20000, 2000000):
        result = b.analyze_counts(fixture(trials=n), IDEAL, service_mean_errors=ZERO)
        assert not result['rejects_joint_premises']
        assert result['status'] == 'not_rejected'
        assert result['quadrature_error_upper']['h_real'] >= abs(result['corrected_means']['h_real']-F(4, 5))


def test_bad_moments_are_rejected_only_when_uncertainty_allows_it():
    counts = fixture(zr=F(9, 10), wr=F(-4, 5))
    resolved = b.analyze_counts(counts, IDEAL, service_mean_errors=ZERO)
    assert resolved['status'] == 'reject'
    assert resolved['rejection_margin_lower'] > 0
    unresolved = b.analyze_counts(fixture(zr=F(9, 10), wr=F(-4, 5), trials=200),
                                  IDEAL, service_mean_errors=ZERO)
    assert not unresolved['rejects_joint_premises']
    uncertain = b.analyze_counts(counts, IDEAL, service_mean_errors=dict.fromkeys(b.SETTINGS, F(1, 4)))
    assert not uncertain['rejects_joint_premises']


def test_prospective_prediction_is_independent_of_later_counts():
    first = fixture()
    prefix = b.analyze_counts({name: first[name] for name in b.FIRST}, IDEAL, service_mean_errors=ZERO)
    changed = dict(first)
    changed['2h_real'], changed['2h_imag'] = [0, 20000, 0, 0], [20000, 0, 0, 0]
    final = b.analyze_counts(changed, IDEAL, service_mean_errors=ZERO)
    assert prefix['status'] == 'prediction_only'
    assert prefix['second_observed_mean_predictions'] == final['second_observed_mean_predictions']
    assert final['rejects_joint_premises']
    predicted = prefix['second_observed_mean_predictions']
    assert predicted['2h_real'][0] <= F(1, 20)+F(17, 20)*F(7, 25) <= predicted['2h_real'][1]
    assert predicted['2h_imag'][0] <= F(1, 20) <= predicted['2h_imag'][1]


def test_effect_sum_and_null_defects_are_included_in_calibration_error():
    c = F(1, 10)
    interface = b.InterfaceBounds((c,)*6, (c,)*2, F(1))
    result = b.analyze_counts(fixture(trials=10**8), interface, service_mean_errors=ZERO)
    # Canonical completion changes (1+c)*.8 to .8. This discrepancy must
    # remain covered even when the sampling radius is very small.
    assert result['quadrature_error_upper']['h_real'] > c*F(4, 5)
    assert result['gain_and_offset_error_upper'] >= c+interface.completion('null_plus')


def test_missing_service_bounds_never_become_zero_and_weak_gain_is_inconclusive():
    counts = fixture(zr=F(9, 10), wr=F(-4, 5))
    services = dict(ZERO); services['2h_imag'] = None
    result = b.analyze_counts(counts, IDEAL, service_mean_errors=services)
    assert result['status'] == 'unbounded_services'
    assert not result['rejects_joint_premises']
    counts['null_minus'] = counts['null_plus']
    result = b.analyze_counts(counts, IDEAL, service_mean_errors=ZERO)
    assert result['status'] == 'inconclusive_gain'
    assert not result['rejects_joint_premises']


def test_both_spin_rows_are_retained_and_only_bits_are_marginalized():
    first = fixture()
    rearranged = {name: [0, 0, row[0], row[1]] for name, row in first.items()}
    assert b.analyze_counts(first, IDEAL, service_mean_errors=ZERO) == b.analyze_counts(rearranged, IDEAL, service_mean_errors=ZERO)


@pytest.fixture(scope='module')
def compiled():
    protocol = e.strict_json.load(e.source.FOLDER/'protocol.json')
    operations, _ = e.source.compile_detectors(4, protocol['settings'], protocol['probe_h'])
    return e.source.complete_raw_operators(4, operations)


def guarantees():
    return {'effect_sum_norm_upper': F(1, 10**12), 'null_observable_norm_upper': F(1, 10**12),
            'moment_defect_upper': F(1, 10**11)}


def test_actual_native_spin_pointer_interface_passes_and_wrong_power_fails(compiled):
    assert e.certify_operators(compiled, guarantees()).moment == F(1, 10**11)
    wrong = dict(compiled); wrong['2h_real'] = compiled['null_minus']
    with pytest.raises(ArithmeticError, match='guarantee'):
        e.certify_operators(wrong, guarantees())


def test_missing_spin_branch_or_concealing_ndarray_cannot_enter_certificate(compiled):
    damaged = dict(compiled); damaged['h_real'] = ((compiled['h_real'][0][0],), compiled['h_real'][1])
    with pytest.raises(ValueError, match='spin'):
        e.certify_operators(damaged, guarantees())
    class Hostile(np.ndarray):
        def __getattribute__(self, name):
            if name in ('shape', 'dtype', 'ndim', 'copy'):
                raise AssertionError('subclass dispatch')
            return super().__getattribute__(name)
    damaged = dict(compiled)
    original = compiled['h_real']
    damaged['h_real'] = ((original[0][0].view(Hostile), original[0][1]), original[1])
    with pytest.raises(ValueError, match='base ndarray'):
        e.certify_operators(damaged, guarantees())


@pytest.mark.parametrize('bad', [[True, 0, 0, 0], [-1, 1, 0, 0], [1., 0, 0, 0], [0, 0, 0, 0], [1, 2]])
def test_invalid_counts_fail_closed(bad):
    counts = fixture(); counts['h_real'] = bad
    with pytest.raises(ValueError):
        b.analyze_counts(counts, IDEAL, service_mean_errors=ZERO)


@pytest.mark.parametrize('n,alpha', [(True,F(1,100)), (0,F(1,100)), (20,True), (20,F(0)), (20,F(1))])
def test_invalid_sampling_inputs_fail_closed(n, alpha):
    with pytest.raises(ValueError):
        b.sampling_radius(n, alpha)


def test_replay_tolerates_metadata_roundoff_but_keeps_hashes_counts_and_bounds_exact():
    from .__main__ import compare_generated
    packet = {'source_expectation_execution_diagnostic': {'probability': .4},
        'source_detector_metadata': {'full_phase_target_residual': 1e-14, 'word_sha256': 'abc', 'pulses': 17},
        'confidence': {'radius': '1/10', 'trials': 20000, 'rejects': False}}
    changed = deepcopy(packet)
    changed['source_detector_metadata']['full_phase_target_residual'] = float(np.nextafter(1e-14, np.inf))
    compare_generated(packet, changed)
    for section, key, value in (('source_detector_metadata','word_sha256','xyz'),
            ('source_detector_metadata','pulses',18), ('confidence','radius','100000000000001/1000000000000000'),
            ('confidence','trials',20001), ('confidence','rejects',True)):
        changed = deepcopy(packet); changed[section][key] = value
        with pytest.raises((ValueError, ArithmeticError)):
            compare_generated(packet, changed)


def test_source_bindings_include_the_strict_verifier_dependency():
    assert 'model/common_source/output_certificate.py' in e.source_bindings()
