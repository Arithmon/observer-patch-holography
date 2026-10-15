"""Independent state, channel and exact-rational enclosure controls."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import bounds as b
from . import experiment as e


def test_moment_inequality_for_arbitrary_non_gaussian_finite_states():
    rng = np.random.default_rng(19271)
    for size in (2, 5, 8):
        raw = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
        rho = raw@raw.conj().T
        rho /= np.trace(rho)
        h = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
        u = expm(1j*(h+h.conj().T))
        first, second = np.trace(rho@u), np.trace(rho@u@u)
        assert not b.moment_witness(first, second, operator_defect=2e-14)['violates_with_allowance']
    assert b.moment_witness(.9, -.8)['violates_with_allowance']


def test_disk_boundary_is_attained_by_a_two_point_phase_distribution():
    for angle in (.1, .9, 1.5):
        z = math.cos(angle)
        w = math.cos(2*angle)
        center, radius = b.moment_disk(z)
        assert abs(w-center) == pytest.approx(radius, abs=1e-15)


def test_nonunitary_compiled_correction_covers_actual_operators():
    rng = np.random.default_rng(141)
    a = np.array([[.9, .07j], [.02, .98j]])
    second = a@a+np.diag([.03, -.05])
    delta = max(np.linalg.norm(a.conj().T@a-np.eye(2), 2), np.linalg.norm(a@a.conj().T-np.eye(2), 2))
    eta = np.linalg.norm(second-a@a, 2)
    for _ in range(10):
        x = rng.normal(size=2)+1j*rng.normal(size=2)
        x /= np.linalg.norm(x)
        assert not b.moment_witness(np.vdot(x, a@x), np.vdot(x, second@x),
                                    operator_defect=delta+eta+1e-14)['violates_with_allowance']


def test_complex_interval_matmul_contains_exact_binary_rational_product():
    a = np.array([[.1+.2j, -.3j], [.7, -.02+.003j]])
    c = np.array([[-.4j, .61], [.009-.8j, -.33j]])
    out = b._multiply(b.iv.complex_exact(a), b.iv.complex_exact(c))
    for i in range(2):
        for j in range(2):
            re = sum((F(a[i,k].real)*F(c[k,j].real)-F(a[i,k].imag)*F(c[k,j].imag)
                      for k in range(2)), F())
            im = sum((F(a[i,k].real)*F(c[k,j].imag)+F(a[i,k].imag)*F(c[k,j].real)
                      for k in range(2)), F())
            assert F(out[0][0][i,j]) <= re <= F(out[0][1][i,j])
            assert F(out[1][0][i,j]) <= im <= F(out[1][1][i,j])


@pytest.fixture(scope='module')
def compiled():
    p = e.strict_json.load(e.FOLDER/'protocol.json')
    return e.compile_detectors(4, p['settings'], p['probe_h'])[0]


def test_native_operator_defects_are_enclosed_and_detect_wrong_doubling(compiled):
    certificate = b.operator_certificates(compiled)
    assert certificate['moment_inequality_defect_upper'] < 1e-9
    altered = dict(compiled)
    altered['2h_real'] = compiled['null_minus']
    mutant = b.operator_certificates(altered)
    assert mutant['B_minus_A_squared_norm_upper'] > 1


def test_certificate_includes_both_spin_record_branches(compiled):
    complete = e.complete_raw_operators(4, compiled)
    certificate = b.operator_certificates(complete)
    assert certificate['moment_inequality_defect_upper'] < 1e-9
    assert all(len(group) == 2 for operation in complete.values() for group in operation)
    damaged = dict(complete)
    groups = complete['h_real']
    damaged['h_real'] = ((groups[0][0], 2*groups[0][1]), groups[1])
    assert b.operator_certificates(damaged)['moment_inequality_defect_upper'] > 1


def test_source_null_controls_do_not_depend_on_preparation_and_correct_readout(compiled):
    noise = e.channel_matrix(.05, .1)
    results = [e.execute(4, a, n, compiled, noise) for a,n in ((-.7,.2),(0.,0.))]
    for result in results:
        records = result['records']
        channel = b.identify_bit_channel(records['null_plus']['mean'], records['null_minus']['mean'])
        assert channel['gain'] == pytest.approx(.85, abs=1e-12)
        assert channel['offset'] == pytest.approx(.05, abs=1e-12)
        assert records['null_imag']['mean'] == pytest.approx(.05, abs=1e-12)
        for row in records.values():
            assert np.sum(row['joint_probabilities']) == pytest.approx(1., abs=1e-12)
            assert row['quantum_output_branches'] == 4
    assert abs(results[0]['records']['h_imag']['mean']-results[1]['records']['h_imag']['mean']) > .1


def test_independent_channel_check_compares_both_spin_rows(compiled, monkeypatch):
    # Put positive mass in both spin records before the native pointer read.
    spin, metadata = e.anharmonic._spin_record(4)
    hadamard = np.kron(np.array([[1., 1.], [1., -1.]])/math.sqrt(2), np.eye(4))
    rotated = tuple(operator@hadamard for operator in spin)
    monkeypatch.setattr(e.anharmonic, '_spin_record', lambda dimension: (rotated, metadata))
    channel = e.channel_matrix(.05, .1)
    executed = e.execute(4, -.7, .2, compiled, channel)
    states, weights, _ = e.native.preparation(4, -.7, .2)
    rho = np.zeros((8, 8), complex)
    rho[:4, :4] = (states*weights)@states.conj().T
    for name, detector in compiled.items():
        expected = e.classical_channel_law(rho, rotated, detector, channel)
        actual = np.asarray(executed['records'][name]['joint_probabilities'])
        assert expected.sum(axis=1) == pytest.approx([.5, .5], abs=1e-12)
        assert e.complete_law_residual(actual, expected) < 1e-12
        damaged = actual.copy()
        damaged[1] = 0.
        assert e.complete_law_residual(damaged, expected) > .2
    with pytest.raises(ValueError, match='2x2'):
        e.complete_law_residual(actual[0], expected[0])


def test_supplied_bit_channel_cannot_be_assumed_when_setting_dependent():
    channel = b.identify_bit_channel(.9, -.8)
    assert channel['offset'] == pytest.approx(.05)
    # Equal sign flips of the two opposite phase means diagnose offset drift;
    # no change in the oscillator state can break the reversal identity.
    forward, reverse = .4+.04, -.3+.04
    assert forward+reverse-2*channel['offset'] == pytest.approx(.08)


def test_uncertainty_allowance_prevents_false_rejection():
    z, w = .8+0j, .28+0j  # Saturates the lower moment-disk boundary.
    assert b.moment_witness(z+.01, w-.01)['violates_with_allowance']
    assert not b.moment_witness(z+.01, w-.01, first_error=.01, second_error=.01)['violates_with_allowance']


def test_sampling_radius_is_for_all_settings_and_calibration_is_charged():
    result = b.sampling_design(20000)
    r = result['simultaneous_raw_mean_radius']
    assert 18*math.exp(-20000*r*r/2) == pytest.approx(.01)
    assert result['corrected_single_quadrature_radius'] == pytest.approx(3*r/.85)
    assert b.sampling_design(1)['positive_true_gain_certifiable'] is False


def test_producer_uses_declared_sampling_confidence_setting_count_and_gain():
    protocol = e.strict_json.load(e.FOLDER/'protocol.json')
    protocol['sampling_design'].update(number_of_settings=12, failure_probability=.025)
    protocol['declared_detector_channel']['flip_zero_to_one'] = .15
    for result in e.protocol_sampling_designs(protocol):
        n, radius = result['fresh_trials_per_setting'], result['simultaneous_raw_mean_radius']
        assert 24*math.exp(-n*radius*radius/2) == pytest.approx(.025)
        assert result['settings'] == 12
        assert result['failure_probability_upper'] == .025
        assert result['corrected_single_quadrature_radius'] == pytest.approx(3*radius/.75)


@pytest.mark.parametrize('p,q', [(-.01,.1),(.1,1.1),(True,.1),(float('nan'),.1)])
def test_invalid_channel_fails_closed(p,q):
    with pytest.raises(ValueError):
        e.channel_matrix(p,q)


@pytest.mark.parametrize('args', [(0,), (True,), (100,0), (100,9,1.), (100,9,.01,0.)])
def test_invalid_sampling_design_fails_closed(args):
    with pytest.raises(ValueError):
        b.sampling_design(*args)
