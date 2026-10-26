"""Adversarial complete-output and receipt checks for the composed observer."""
from copy import deepcopy
import json

import numpy as np
import pytest

from . import record_field as field, record_observer as observer


@pytest.fixture(scope='module')
def replayed():
    return observer.run()


def test_retained_evidence_replays(replayed):
    retained = json.loads((observer.ROOT/'model/common_source/record_observer_results.json').read_text())
    # Exercise actual producer types, not a fixture that has already crossed JSON.
    observer._compare(json.loads(json.dumps(replayed, allow_nan=False)), replayed)
    observer._compare(retained, replayed)
    assert replayed['optimal_native_code']
    assert len(replayed['complete_field_reference_arrays']) == 17


@pytest.mark.parametrize('mutation', ('nan', 'missing_timeout', 'wrong_shape', 'bool_count',
                                    'false_scope', 'false_bound', 'false_full_carrier_bound', 'source_digest'))
def test_receipt_cannot_relabel_or_drop_evidence(replayed, mutation):
    bad = deepcopy(replayed)
    if mutation == 'nan':
        bad['complete_field_reference_arrays']['initial']['real'][0][0] = float('nan')
    elif mutation == 'missing_timeout':
        del bad['complete_field_reference_arrays']['selector_timeout_second_0_30']
    elif mutation == 'wrong_shape':
        bad['complete_field_reference_arrays']['initial']['shape'] = [128, 1]
    elif mutation == 'bool_count':
        bad['native_source_execution']['executed_code_columns'] = True
    elif mutation == 'false_scope':
        bad['scope']['physical_theory_completed'] = True
    elif mutation == 'false_bound':
        bad['rank15_recovery_obstruction']['E']['entanglement_fidelity_upper'] = '1'
    elif mutation == 'false_full_carrier_bound':
        bad['rank15_recovery_obstruction']['full_total_two']['C']['entanglement_fidelity_upper'] = '43/45'
    else:
        key = next(iter(bad['source_sha256']))
        bad['source_sha256'][key] = '0'*64
    with pytest.raises(ArithmeticError):
        observer._compare(bad, replayed)


def test_same_probabilities_do_not_certify_quantum_outputs(monkeypatch):
    original = field.detector_experiment

    def corrupted():
        summary, arrays = original()
        arrays['initial'][:, 1] *= -1
        # Probability, normalization and reference marginals are unchanged.
        assert np.isclose(np.vdot(arrays['initial'], arrays['initial']), 1)
        return summary, arrays

    monkeypatch.setattr(field, 'detector_experiment', corrupted)
    with pytest.raises(ArithmeticError, match='independent oracle'):
        observer.verify({})  # Its owned replay rejects before receipt comparison.


def test_missing_actual_zero_output_fails_owned_replay(monkeypatch):
    original = field.detector_experiment

    def incomplete():
        summary, arrays = original()
        del arrays['terminal_0_1']
        return summary, arrays

    monkeypatch.setattr(field, 'detector_experiment', incomplete)
    with pytest.raises(ArithmeticError, match='complete oracle'):
        observer.verify({})


def test_false_unitary_walk_cannot_pass_source_execution(monkeypatch):
    original = field.compiler.execute

    def identity_walk(word, qubits, *args, **kwargs):
        return np.eye(64, dtype=complex) if qubits == 6 else original(word, qubits, *args, **kwargs)

    field.source_walk.cache_clear()
    monkeypatch.setattr(field.compiler, 'execute', identity_walk)
    try:
        with pytest.raises(ArithmeticError, match='signed exterior'):
            observer.verify({})
    finally:
        field.source_walk.cache_clear()


def test_source_custody_is_explicit_not_import_order():
    before = observer.SOURCE_PATHS
    import model.common_source.mass_response.model  # An unrelated loaded module.
    assert observer.SOURCE_PATHS == before
    for module in ('record_capacity', 'record_code', 'record_field', 'record_information', 'record_observer'):
        assert f'model/common_source/{module}.py' in before


def test_backend_residuals_use_explicit_gates_not_bit_identity():
    path = ('root.protected_field.complete_source_composition.laws.C.'
            'maximum_complete_two_cycle_operator_error_diagnostic')
    observer._compare(7e-9, 2e-8, path)
    for bad in (1e-3, -1e-12, float('nan'), float('inf')):
        with pytest.raises(ArithmeticError, match='gate'):
            observer._compare(bad, 2e-8, path)
    with pytest.raises(ArithmeticError, match='gate changed'):
        observer._compare(1.1e-8, 1e-8,
                          'root.native_source_execution.maximum_allowed_numerical_execution_error')
    # Actual quantum amplitudes still use the tight numerical comparison.
    with pytest.raises(ArithmeticError, match='numerical replay'):
        observer._compare(0., 2e-8, 'root.complete_field_reference_arrays.initial.real[0][0]')


@pytest.mark.parametrize('path,gate', [
    ('root.native_source_execution.encoder_operator_error', 1e-8),
    ('root.native_source_execution.laws.C.reference_gram_error', 1e-8),
    ('root.native_source_execution.laws.E.maximum_branch_operator_error', 1e-8),
    ('root.native_source_execution.laws.E.complete_retained_output_operator_error', 1e-8),
    ('root.native_source_execution.laws.C.branches[0].operator_error', 1e-8),
    ('root.native_source_execution.forward_full512_source_macro_guard.raw_unitarity_defect', 1e-10),
    ('root.native_source_execution.inverse_full512_source_macro_guard.projection_change', 5.1e-11),
    ('root.native_source_execution.forward_full512_source_macro_guard.restored_unitarity_defect', 2e-13),
    ('root.all_seam_source_composition.forward_full512_macro_polar_change', 5.1e-11),
    ('root.protected_field.dispersion.modes[0].source_residual', 2e-11),
    ('root.protected_field.source.source_walk_max_entry_residual', 2e-11),
    ('root.protected_field.source.walk_resources.maximum_numerical_valid_code_leakage', 2e-11),
    ('root.protected_field.source.feedback_resources.maximum_numerical_valid_code_leakage', 2e-11),
    ('root.protected_field.complete_source_composition.walk_operator_error', 64*2e-11),
    ('root.protected_field.complete_source_composition.feedback_operator_error', 64*2e-11),
    ('root.protected_field.complete_source_composition.proper_code_export_diagnostic.maximum_complement_operator_norm', 2e-11),
    ('root.protected_field.complete_source_composition.proper_code_export_diagnostic.maximum_success_norm_excess', 2e-11),
    ('root.protected_field.dispersion.charged_retarded_anticommutator_response_residual', 3e-11),
    ('root.protected_field.clocks.joint_clock_histories.source_executed_stopped_channel.complete_retained_stopped_operator_error', 5e-11),
    ('root.protected_field.clocks.joint_clock_histories.source_executed_stopped_channel.complete_stopped_gram_error', 1e-10),
])
def test_retained_residuals_obey_the_producer_gate_before_comparison(path, gate):
    observer._compare(gate, gate, path)
    for invalid in (1.01*gate, -gate/100, float('nan'), float('inf')):
        with pytest.raises(ArithmeticError, match='gate'):
            observer._compare(invalid, .99*gate, path)


@pytest.mark.parametrize('name,gate', [
    ('maximum_raw_unitarity_defect', 1e-10),
    ('maximum_projection_change', 5.1e-11),
    ('maximum_restored_unitarity_defect', 2e-13),
])
def test_retained_polar_acceptance_thresholds_are_fixed(name, gate):
    path = 'root.native_source_execution.forward_full512_source_macro_guard.'+name
    observer._compare(gate, gate, path)
    with pytest.raises(ArithmeticError, match='gate changed'):
        observer._compare(gate+1e-15, gate, path)


def _retained():
    return json.loads((observer.ROOT/'model/common_source/record_observer_results.json').read_text())


@pytest.mark.parametrize('mutation', ['one_entry', 'distributed_entries'])
def test_retained_field_arrays_obey_whole_operator_oracle(mutation):
    good = _retained()['complete_field_reference_arrays']
    bad = deepcopy(good)
    if mutation == 'one_entry':
        bad['initial']['real'][0][0] = 1.5e-9
    else:
        # Every entry passes the ordinary 2e-9 comparison, but their operator
        # norm is sqrt(62)*2e-10 > the producer's complete-output gate.
        for index in range(2, 64):
            bad['terminal_0_1']['real'][index][0] = 2e-10
    with pytest.raises(ArithmeticError, match='independent oracle'):
        observer._compare(bad, good, 'root.complete_field_reference_arrays')


@pytest.mark.parametrize('mutation', ['polar_formula', 'macro_sum', 'encoder_exponential',
                                    'native_exponential', 'completion_factor', 'diamond_conversion',
                                    'field_native_join', 'export_failure_join',
                                    'detector_probability', 'total_probability'])
def test_stored_diagnostic_formulas_use_their_own_scale(mutation):
    good = _retained()
    bad = deepcopy(good)
    if mutation == 'polar_formula':
        bad['native_source_execution']['forward_full512_source_macro_guard']['exact_polar_distance_formula_evaluated'] = 0.
    elif mutation == 'macro_sum':
        bad['all_seam_source_composition']['forward_full512_macro_operator_error'] = 0.
    elif mutation == 'encoder_exponential':
        bad['all_seam_source_composition']['encoder_operator_error_bound'] = 0.
    elif mutation == 'native_exponential':
        bad['all_seam_source_composition']['laws']['C']['maximum_corrected_attempt_operator_error'] = 0.
    elif mutation == 'completion_factor':
        bad['all_seam_source_composition']['mixing_angle_synthesis_diagnostics']['6']['uniform_target_unitary_completion_error_upper'] = 0.
    elif mutation == 'diamond_conversion':
        bad['protected_field']['complete_source_composition']['laws']['C']['maximum_complete_cq_diamond_error_diagnostic'] = 0.
    elif mutation == 'field_native_join':
        row = bad['protected_field']['complete_source_composition']['laws']['C']
        epsilon = 2*row['maximum_complete_two_cycle_operator_error_diagnostic']
        row['maximum_complete_two_cycle_operator_error_diagnostic'] = epsilon
        row['maximum_complete_cq_diamond_error_diagnostic'] = epsilon*(2+epsilon)
        row['maximum_complete_cq_half_trace_error_diagnostic'] = epsilon*(1+epsilon/2)
    elif mutation == 'export_failure_join':
        row = bad['protected_field']['complete_source_composition']['laws']['C']
        row['maximum_proper_code_failure_amplitude_diagnostic'] = 0.
        epsilon = row['maximum_successful_block_two_cycle_operator_error_diagnostic']
        row['maximum_complete_two_cycle_operator_error_diagnostic'] = epsilon
        row['maximum_complete_cq_diamond_error_diagnostic'] = epsilon*(2+epsilon)
        row['maximum_complete_cq_half_trace_error_diagnostic'] = epsilon*(1+epsilon/2)
    elif mutation == 'detector_probability':
        bad['protected_field']['two_read_experiment']['complete_detector_outputs'][0]['source_probability_before_selector_weights'] += 1e-9
    else:
        bad['protected_field']['two_read_experiment']['weighted_terminal_source_probability_sum'] += 1e-9
    with pytest.raises(ArithmeticError, match='diagnostic|gate'):
        observer._compare(bad, good)


def test_source_walk_and_returned_field_arrays_cannot_poison_later_outputs():
    walk = field.source_walk()
    with pytest.raises(ValueError):
        walk.setflags(write=True)
    _, first = field.detector_experiment()
    expected = first['terminal_0_0'].copy()
    first['terminal_0_0'][:] = 123
    _, second = field.detector_experiment()
    np.testing.assert_array_equal(second['terminal_0_0'], expected)
