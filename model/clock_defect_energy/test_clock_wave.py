"""Independent native arithmetic and hostile inputs for the clock/wave join."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
import json
import math

import numpy as np
import pytest
import sympy as sp

from . import clock_wave as w


@pytest.fixture(scope='module')
def result():
    return w.report()


@pytest.fixture(scope='module')
def independent_counts():
    states = tuple(sorted(tuple(int(p in pair) for p in range(12))
                          for pair in combinations(range(12), 2)))
    index = {x: j for j, x in enumerate(states)}
    counts = sp.zeros(66)
    for j, state in enumerate(states):
        for a, b in w.transport.records.seams():
            for coin in (0, 1):
                output = list(state)
                total = output[a]+output[b]
                output[a], output[b] = (total+coin)//2, (total+1-coin)//2
                counts[index[tuple(output)], j] += 1
    return sp.ImmutableMatrix(counts), states


def test_same_carrier_native_counts_and_all_transported_cores(result, independent_counts):
    counts, states = independent_counts
    assert counts == w._source()[0]
    rows = result['transported_clock_vs_wave']
    assert [x['complete_two_attempt_branches'] for x in rows['laws']] == [15, 34]
    for row in rows['core_certificates']:
        weights = tuple((F(2, 225), F(1, 36), F(2, 45))[sum(x[p] for p in row['core'])]
                        for x in states)
        i, j, entry = row['first_matrix_unit_witness']
        assert F(entry) == (weights[i]-weights[j])*int(counts[i, j])/60
        assert F(entry) != 0
        oracle = sum(((weights[i]-weights[j])**2*int(counts[i, j])**2/F(3600)
                      for i in range(66) for j in range(66)), F())
        assert oracle == F(row['squared_frobenius_norm'])


def test_krylov_rank_certificate_exact_without_eigensolver(result, independent_counts):
    counts, _ = independent_counts
    vectors = [sp.eye(66)[:, 0]]
    for _ in range(3):
        vectors.append(counts*vectors[-1])
    matrix = sp.Matrix.hstack(*vectors)
    gram = matrix.T*matrix
    assert gram.det() == 1556992
    assert result['transported_clock_vs_wave']['full_shell_spectrum']['exact_gram_matrix'] == [
        list(map(str, row)) for row in gram.tolist()]
    # A three-level diagonal changed kernel makes the same method fail rank4.
    diagonal = sp.diag(*([0, 1, 2]*22))
    seed = sp.ones(66, 1)
    changed = sp.Matrix.hstack(*(diagonal**j*seed for j in range(4)))
    assert (changed.T*changed).det() == 0


def test_density_subspace_escape_includes_constant_mode(result, independent_counts):
    _, states = independent_counts
    z = sp.Matrix(12, 66, lambda p, j: sp.Rational(states[j][p])-sp.Rational(1, 6))
    p = z.T*z/10+sp.ones(66)/66
    assert p*p == p and sp.trace(p) == 12
    for row in result['transported_clock_vs_wave']['core_certificates']:
        sigma = sp.diag(*[(sp.Rational(2, 225), sp.Rational(1, 36), sp.Rational(2, 45))[
            sum(x[j] for j in row['core'])] for x in states])
        leak = (sp.eye(66)-p)*sigma*z.T
        assert sum(x*x for x in leak) == sp.Rational(2119, 371250)
        assert row['constant_plus_density_reference_leakage_squared_frobenius'] == '2119/371250'


def test_affine_shape_enclosures_against_independent_float_functions(result):
    row = result['three_band_affine_shape']
    omega = [math.acos(1-x/60) for x in (5-math.sqrt(5), 6, 5+math.sqrt(5))]
    shape = (omega[2]-omega[0])/(omega[1]-omega[0])
    values = [shape, math.log(5)/math.log(8/5), math.log(5)/math.log(25/8)]
    keys = ['wave_gap_shape_interval', 'modular_positive_orientation_shape_interval',
            'modular_reversed_orientation_shape_interval']
    for key, value in zip(keys, values):
        lo, hi = map(F, row[key])
        assert float(lo) < value < float(hi)
    assert min(abs(values[0]-v) for v in values[1:]) > float(F(row['positive_separation_lower']))


def test_supplied_two_mode_gibbs_control_and_predual_sign(result, independent_counts):
    counts, states = independent_counts
    mode = sp.Matrix([x[0]-x[1]-x[2]+x[3] for x in states])
    assert counts*mode == 54*mode
    assert sum(mode) == 0
    theta = math.acos(.9)
    beta = math.log(2)/theta
    parameter = -1/beta
    phase = np.exp(-1j*parameter*np.log([2/3, 1/3]))
    assert abs(phase[1]/phase[0]-np.exp(-1j*theta)) < 1e-15
    probability = abs(sum(phase)/2)**2
    assert abs(probability-.95) < 1e-15
    assert result['changed_reference_positive_control']['one_wave_step_balanced_probe_plus_probability_exact'] == '19/20'


def test_tracial_commutation_is_necessary_not_sufficient(result):
    control = result['transported_clock_vs_wave']['tracial_control']
    assert control['commutation_necessary_condition_passes']
    assert control['nonzero_entries'] == 0
    assert control['first_matrix_unit_witness'] is None
    assert result['transported_clock_vs_wave']['full_shell_spectrum']['wave_distinct_generator_levels_lower'] == 4


def test_existing_native_return_bridge_reconstructs_same_kernel(independent_counts, result):
    counts, _ = independent_counts
    kernel = counts/60
    polynomial = sp.eye(66)+kernel+kernel**2/2
    reference = polynomial/sp.trace(polynomial)
    original = w.wave_reference.reference(2)
    assert reference == original['density']
    assert reference*kernel == kernel*reference
    assert reference != original['native_return_state']
    row = result['existing_positive_native_return_bridge']['exact_finite_cap_checks'][0]
    assert row['reconstructed_density_sha256'] == w.fields._digest(reference)
    assert sum(F(row[key]) for key in ('return_probability', 'different_probability', 'cap_probability')) == 1
    assert F(row['reconstructed_minus_actual_return_density_squared_frobenius']) > 0
    assert result['existing_positive_native_return_bridge']['new_theorem_claimed'] is False


def test_readout_calibration_works_without_equal_generators():
    # This exact clock Hamiltonian does not commute with this field Hamiltonian.
    clock = np.array([[0., 1.], [1., 0.]])
    field = np.diag([0., .5])
    assert np.linalg.norm(clock@field-field@clock) > 0
    bound = w.operational_parameter_bound((F(2, 3), F(4, 5)), F(3, 2), F(1, 2))
    middle = float(F(bound['midpoint']))
    epsilon = float(F(bound['unitary_and_cosine_operator_error_upper']))
    for parameter in np.linspace(2/3, 4/5, 31):
        a = np.diag(np.exp(-1j*1.5*parameter*np.diag(field)))
        b = np.diag(np.exp(-1j*1.5*middle*np.diag(field)))
        assert np.linalg.norm(a-b, 2) <= epsilon+1e-15
    assert not bound['clock_and_field_hamiltonians_required_equal']
    backward = w.operational_parameter_bound((F(2, 3), F(4, 5)), F(-3, 2), F(1, 2))
    assert backward['wave_parameter_interval'] == ['-6/5', '-1']
    assert backward['unitary_and_cosine_operator_error_upper'] == bound['unitary_and_cosine_operator_error_upper']


def test_complete_record_inverse_propagates_its_original_joint_failure(result):
    bridge = result['operational_clock_readout_to_field_bound']
    assert bridge['replayed_from_retained_counts_without_programmed_parameter']
    assert not bridge['new_joint_clock_field_instrument_executed']
    for row in bridge['readouts']:
        assert row['joint_clock_confidence_failure_upper'] == '1/128'
        bound = row['conditional_density_band_prediction']
        assert F(bound['radius']) > 0
        assert F(bound['unitary_and_cosine_operator_error_upper']) == F(bound['radius'])/2
        assert not bound['clock_and_field_hamiltonians_required_equal']
    assert result['scope']['earlier_common_clock_test_framing_withdrawn']
    assert not result['scope']['negative_result_excludes_a_common_operational_clock']


@pytest.mark.parametrize('interval,rate,omega', [
    ([F(0), F(1)], F(1), F(1)), ((F(1), F(0)), F(1), F(1)),
    ((F(0), F(1)), F(0), F(1)), ((F(0), F(1)), F(1), F(-1)),
    ((False, F(1)), F(1), F(1)), ((F(0), F(1)), 1., F(1))])
def test_operational_error_bound_rejects_malformed_contracts(interval, rate, omega):
    with pytest.raises(ValueError):
        w.operational_parameter_bound(interval, rate, omega)


@pytest.mark.parametrize('bad', [True, .5, '1/2', complex(1)])
def test_nonexact_interval_arguments_rejected(bad):
    with pytest.raises(ValueError):
        w._sqrt_bounds(bad)
    with pytest.raises(ValueError):
        w._log_bounds(bad)


def test_corrupt_fraction_subclass_and_invalid_domains_rejected():
    class Child(F):
        pass
    bad = F(1, 2)
    bad._numerator = True
    for x in (bad, Child(1, 2)):
        with pytest.raises(ValueError):
            w._sqrt_bounds(x)
    for endpoints in ([F(1), F(2)], (F(2), F(1)), (F(0), F(1)), (F(1), F(13))):
        with pytest.raises(ValueError):
            w._wave_frequency_bounds(endpoints)


def test_changed_kernel_and_reference_rejected(independent_counts):
    counts, _ = independent_counts
    sigma = (F(1, 66),)*66
    with pytest.raises(ValueError):
        w.commutator_certificate(sigma[:-1], counts)
    with pytest.raises(ValueError):
        w.commutator_certificate((F(0), *sigma[1:]), counts)
    changed = sp.MutableDenseMatrix(counts)
    changed[0, 0] += 1
    with pytest.raises(ValueError):
        w.commutator_certificate(sigma, sp.ImmutableMatrix(changed))
    with pytest.raises(ValueError):
        w.commutator_certificate(sigma, sp.ImmutableMatrix([[60]]))


def test_receipt_custody_and_type_sensitive_mutations(result):
    retained = w.strict_json.loads(w.OUTPUT.read_text())
    w.check_receipt(retained, result)
    for key, value in [('physical_c_or_lab_duration_selected', True),
                       ('new_empirical_comparisons', False)]:
        changed = deepcopy(result)
        changed['scope'][key] = value
        with pytest.raises(ValueError):
            w.check_receipt(changed, result)
    changed = deepcopy(result)
    changed['transported_clock_vs_wave']['core_certificates'][0]['squared_frobenius_norm'] = '0'
    with pytest.raises(ValueError):
        w.check_receipt(changed, result)
    assert json.loads(json.dumps(result)) == result


def test_duplicate_json_keys_rejected():
    with pytest.raises(ValueError):
        w.strict_json.loads('{"clock": 1, "clock": 2}')
