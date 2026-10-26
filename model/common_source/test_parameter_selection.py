"""Fixed native-cycle response, correlation compression and scope controls."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import parameter_selection as selection


@pytest.fixture(scope='module')
def evidence():
    return selection.report()


def test_every_port_and_centered_native_family_coordinate_is_scalarized():
    geometry = selection.code_geometry()
    states, v = geometry['states'], geometry['V']
    diagnostic_tolerance = 64*max(v.shape)*np.finfo(float).eps*np.linalg.norm(v, 2)**2
    for port in range(12):
        assert [sum(states[i][port] for i in row) for row in geometry['supports']] == [1]*6
        actual = v.conj().T@np.diag([state[port] for state in states])@v
        np.testing.assert_allclose(actual, np.eye(6)/6, rtol=0, atol=diagnostic_tolerance)
    for axis in range(3):
        values = [sum(geometry['frame'][port, axis]*state[port] for port in range(12)) for state in states]
        actual = v.conj().T@np.diag(list(map(float, values)))@v
        assert all(s.simplify(sum(values[i] for i in row)) == 0 for row in geometry['supports'])
        np.testing.assert_allclose(actual, 0, rtol=0, atol=diagnostic_tolerance*max(1, max(map(abs, map(float, values)))))
    assert geometry['feature'] == s.eye(3)/6


def test_every_function_of_original_core_reference_has_same_compression():
    geometry = selection.code_geometry()
    f0, f1, f2 = s.symbols('f0 f1 f2')
    values = [f0, f1, f2]
    core = selection.preparation.CORE
    compressed = [sum(values[sum(geometry['states'][i][port] for port in core)] for i in row)/6
                  for row in geometry['supports']]
    assert compressed == [(2*f0+f1)/3]*6


def test_arbitrary_native_edge_agreement_cost_is_scalar_not_just_uniform_cost():
    geometry = selection.code_geometry()
    weights = {edge: s.Integer(i+1) for i, edge in enumerate(selection.code.SEAMS)}
    values = [sum(weight*(state[a]-state[b])**2 for (a, b), weight in weights.items())
              for state in geometry['states']]
    means = [sum(values[i] for i in row)/6 for row in geometry['supports']]
    assert means == [sum(weights.values())/3]*6
    assert selection.diagonal_cost(1)['compression'] == s.Rational(9, 32)*s.eye(6)


def test_actual_corrected_channel_on_all36_operators_and_all_raw_records(evidence):
    report, _ = evidence
    assert report['fixed_cycle']['exact_all_code_operator_transfer'] == 'I36'
    assert report['fixed_cycle']['static_seed'] == '0_36x36'
    assert report['fixed_cycle']['positive_source_family_gaps_generated'] == 0
    for law, count in [('C', 240), ('E', 480)]:
        row = report['fixed_cycle']['C_E'][law]
        assert row['native_raw_records_checked'] == count
        assert row['raw_zero_records_retained'] > 0
        assert row['logical_transfer_identity_error'] < 2e-12
        assert sum(F(branch['weight']) for branch in row['common_raw_record_factor']) == 1
    assert report['premises']['terminal_identity_timeout_probability'] == '5/32'
    assert report['fixed_cycle']['C_E']['C']['original_faithful_reference_code_projector_probability_after'] == '241/1500'
    assert report['fixed_cycle']['C_E']['C']['original_reference_half_trace_distance_lower_bound'] == '313/4500'
    assert report['fixed_cycle']['C_E']['E']['original_faithful_reference_code_projector_probability_after'] == '26627/115200'
    assert report['fixed_cycle']['C_E']['E']['original_reference_half_trace_distance_lower_bound'] == '5377/38400'
    assert not report['fixed_cycle']['C_E']['C']['original_same_carrier_faithful_reference_stationary']


def test_wrong_native_map_or_logical_recovery_is_rejected(monkeypatch):
    original = selection.coins.coin_maps
    def omitted_input(law, seam, binary=True):
        rows = original(law, seam, binary)
        if law == 'C' and seam == (0, 1):
            coin, label, mapping = rows[0]
            return ((coin, label, mapping[1:]),)+rows[1:]
        return rows
    monkeypatch.setattr(selection.coins, 'coin_maps', omitted_input)
    with pytest.raises(ArithmeticError, match='floor arithmetic'):
        selection.corrected_cycle()
    monkeypatch.setattr(selection.coins, 'coin_maps', original)
    recovery = selection.code.declared_recovery_matrix66
    def wrong_recovery(law, seam, coin, raw):
        matrix = recovery(law, seam, coin, raw).copy()
        matrix[0] *= -1
        # A native-code row actually used by the first logical matching.
        geometry = selection.code_geometry()
        matrix[geometry['supports'][0][0]] *= -1
        return matrix
    monkeypatch.setattr(selection.code, 'declared_recovery_matrix66', wrong_recovery)
    with pytest.raises(ArithmeticError, match='same-carrier recovery'):
        selection.corrected_cycle()


def test_two_hop_candidate_is_fixed_but_does_not_preserve_code(evidence):
    probe = selection.diagonal_cost(2)
    assert tuple(probe['compression'].diagonal()) == (s.Rational(113, 3072), s.Rational(125, 3072),
        s.Rational(29, 768), s.Rational(113, 3072), s.Rational(37, 1024), s.Rational(55, 1536))
    assert probe['occupation_differences'] == (s.Rational(1, 256), -s.Rational(1, 1024), -s.Rational(1, 3072))
    assert all(value > 0 for value in probe['variance'])
    assert not probe['code_preserving']
    geometry = selection.code_geometry()
    v = geometry['V']
    outside = (np.eye(66)-v@v.conj().T)@np.diag(list(map(float, probe['values'])))@v
    diagnostic_tolerance = 64*max(v.shape)*np.finfo(float).eps*max(map(float, probe['values']))**2*np.linalg.norm(v, 2)**2
    np.testing.assert_allclose(outside.conj().T@outside, np.diag(list(map(float, probe['variance']))), rtol=0, atol=diagnostic_tolerance)
    report, _ = evidence
    for law in ('C', 'E'):
        assert report['fixed_cycle']['C_E'][law]['all_record_first_order_diagonal'] == [str(probe['compression'][j, j]) for j in range(6)]
        assert all(s.Rational(value) > 0 for value in report['fixed_cycle']['C_E'][law]['perturbed_complete_record_syndrome_leakage_quadratic_coefficients'])
    assert not report['first_surviving_correlation']['already_present_as_coherent_term_in_fixed_native_cycle']


def test_graph_invariance_cannot_give_three_positive_occupation_differences():
    result = selection.graph_invariant_control()
    assert result['actual_incidence_automorphisms'] == 120
    assert result['unordered_pair_orbit_sizes'] == [6, 30, 30]
    assert result['occupation_difference_rank_upper_for_fixed_2s_plus_n_factorization'] == 1
    assert result['uniform_two_hop_preserves_code_exactly']
    assert selection.diagonal_cost(2, uniform=True)['occupation_differences'] == (4, 0, 0)


def test_same_operator_different_valid_dictionary_does_not_select_spectrum():
    result = selection.dictionary_control()
    first = result['dictionaries']['existing_2s_plus_n']
    second = result['dictionaries']['alternative_same_operator']
    assert first['occupation_differences'] == ['1/256', '-1/1024', '-1/3072']
    assert second['occupation_differences'] == ['0', '-7/1536', '-1/512']
    assert result['fixed_native_instrument_and_H2_identical']
    assert not result['physical_charge_or_family_dictionary_selected_here']


def test_actual_source_cycles_retain_every_quantum_output_and_padding(evidence):
    report, arrays = evidence
    assert arrays['source_encoded_six_columns'].shape == (512, 6)
    assert arrays['native_C_corrected_outputs'].shape == (8, 512, 6)
    assert arrays['native_E_corrected_outputs'].shape == (16, 512, 6)
    for value in arrays.values():
        assert value.dtype == np.dtype('complex128')
        assert np.isfinite(value).all()
    source = report['source_execution']
    assert source['two_hop_source_operator_error'] < selection.SOURCE_GATE
    assert not source['projection_or_padding_removal_or_probability_normalization_applied']
    assert not source['proper_services_and_physical_errors_authenticated_here']
    assert not report['scope']['calibrated_packets_used']
    assert not report['scope']['Standard_Model_mass_derivation_completed']


def test_receipt_cannot_relabel_baseline_or_hide_leakage(evidence):
    report, _ = evidence
    for path in ('baseline', 'leakage', 'bound'):
        wrong = deepcopy(report)
        if path == 'baseline':
            wrong['fixed_cycle']['positive_source_family_gaps_generated'] = 3
        elif path == 'leakage':
            wrong['first_surviving_correlation']['preserves_code_exactly'] = True
        else:
            wrong['source_execution']['numerical_source_admission'] += 1e-12
        with pytest.raises(ValueError):
            selection._compare(report, wrong)


def test_retained_source_arrays_have_independent_phase_zero_and_padding_guards(evidence):
    _, arrays = evidence
    selection.verify_source_arrays(arrays)
    for kind in ('phase', 'zero', 'padding', 'schema'):
        changed = dict(arrays)
        changed['native_C_corrected_outputs'] = arrays['native_C_corrected_outputs'].copy()
        if kind == 'phase':
            changed['native_C_corrected_outputs'][0] *= -1
        elif kind == 'zero':
            changed['native_C_corrected_outputs'][0] = 0
        elif kind == 'padding':
            changed['native_C_corrected_outputs'][0, -1, 0] = 1e-5
        else:
            changed['extra'] = np.zeros(1)
        with pytest.raises(ValueError):
            selection.verify_source_arrays(changed)


def test_probe_domain_does_not_admit_a_parameter_scan():
    for value in (0, 3, True, 1., '2'):
        with pytest.raises(ValueError):
            selection.diagonal_cost(value)
