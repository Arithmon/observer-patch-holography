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
    observer._compare(retained, replayed)
    assert replayed['optimal_native_code']
    assert len(replayed['complete_field_reference_arrays']) == 17


@pytest.mark.parametrize('mutation', ('nan', 'missing_timeout', 'wrong_shape', 'bool_count',
                                    'false_scope', 'false_bound', 'source_digest'))
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
