"""Maximal reference compatibility and complete factored source outputs."""
from copy import deepcopy
import math

import numpy as np
import pytest
import sympy as s

from . import native_response_reference as m


@pytest.fixture(scope='module')
def executed():
    return m.execute()


def test_reference_constraint_is_maximal_and_has_the_required_center():
    report = m.exact_certificate()
    assert report['ray_constraint_real_rank'] == 5
    assert report['ray_stabilizer_raw_dimension'] == 13
    assert report['observable_dimension_on_EndM6'] == 12
    assert report['without_class_sum']['raw_dimension'] == 12
    assert report['without_class_sum']['observable_dimension'] == 11
    assert report['coherent_mixed_Bell_reference_control']['raw_dimension'] == 7
    assert report['coherent_mixed_Bell_reference_control']['observable_dimension'] == 6
    assert s.Matrix(report['five_excluded_current_reference_Gram_exact']).det() > 0
    for raw in report['ray_stabilizer_basis_exact']:
        g = s.Matrix(raw)
        B = g[3:, 3:]
        assert B+B.T == 2*s.trace(B)*s.eye(3)/3
    for raw in report['canonical_section_basis_exact']:
        g = s.Matrix(raw)
        assert s.im(g[3:, 3:]) == s.zeros(3)
        assert g[3:, 3:]+g[3:, 3:].T == s.zeros(3)


def test_source_has_nontrivial_compatible_rotation_and_all_five_phase_damage_controls(executed):
    report, arrays = executed
    values = report['execution']['source_reference_survival_probabilities']
    assert min(values[:15]) > 1-1e-8
    assert max(values[15:]) < .95
    assert abs(values[15]-(1+2*math.cos(.5))**2/9) < 1e-8
    assert np.linalg.norm(arrays['controls'][1][3:6, 3:6]-np.eye(3)) > .005
    assert abs(report['execution']['class_sum_cross_triplet_marker_probability']
               -math.cos(math.sqrt(5)/8)**2) < 1e-8
    assert report['error_budget']['complete_flagged_operator_error'] < 1e-5
    assert report['retained_output_contract']['raw_records'] == 16384
    assert report['resources']['declared_services']['native_record_bits'] == 14


def _records():
    result = []
    weight = 1.
    for law in m.LAWS:
        code = m.transport.PAIR_CODES[law]
        record = next((coin, raw) for coin in (0, 1)
            for raw in range(1 << m.circuit.RECORD_BITS[law])
            if code.branch_weight(law, m.transport.SEAMS[law], coin, raw))
        result.append(record)
        weight *= float(code.branch_weight(law, m.transport.SEAMS[law], *record))
    return tuple(result), weight


@pytest.mark.parametrize('control', [0, 1, 13, 15, 19])
def test_actual_raw_source_outputs_match_independent_complete_pointer_law(executed, control):
    report, arrays = executed
    records, weight = _records()
    phi = m._reference_vector().astype(complex)
    U = m._targets()[control]
    moved = np.kron(U, U)@phi
    yes = phi*np.vdot(phi, moved)
    expected = (yes, moved-yes)
    amplitudes = []
    for pointer in (0, 1):
        actual = m.contract_output(arrays, control, records, (0, 0, 0, 0), pointer)
        target = np.zeros(4096, complex)
        for a in range(8):
            for b in range(8):
                target[512*a+8*b] = math.sqrt(weight)*expected[pointer][8*a+b]
        assert actual.shape == (4096,)
        assert np.linalg.norm(actual-target) < 5e-8
        amplitudes.append(actual)
    assert abs(sum(np.vdot(a, a).real for a in amplitudes)-weight) < 1e-8
    high = m.contract_output(arrays, control, records, (1, 0, 0, 0), 0)
    assert np.linalg.norm(high) < 5e-8


def test_arbitrary_reference_columns_and_zero_raw_outputs_retained(executed):
    _, arrays = executed
    records, weight = _records()
    incoming = np.zeros((64, 2), complex)
    incoming[0, 0] = incoming[27, 1] = 1/math.sqrt(2)
    out = m.contract_output(arrays, 1, records, (0, 0, 0, 0), 1, incoming)
    assert out.shape == (4096, 2)
    assert np.linalg.norm(out) > .05*math.sqrt(weight)
    aborted = list(records)
    aborted[0] = (0, 3)
    zero = m.contract_output(arrays, 0, tuple(aborted), (0, 0, 0, 0), 0)
    assert zero.shape == (4096,) and np.linalg.norm(zero) < 1e-8
    complement = m.contract_output(arrays, 0, records, (0, 0, 0, 0), 0, geometric_flag=1)
    assert complement.shape == (4096,) and np.linalg.norm(complement) < 1e-8


def test_tensor_axis_order_with_nonscalar_independent_local_factors():
    arrays = {'controls': np.eye(8, dtype=complex)[None],
        'geometric_preparation': np.eye(64, dtype=complex),
        'geometric_inverse': np.eye(64, dtype=complex),
        'pointer': np.asarray([np.eye(64), np.zeros((64, 64))], complex),
        'geometric_pointer': np.asarray([np.eye(64), np.zeros((64, 64))], complex)}
    for law, start, count in (('C', 1, 8), ('E', 11, 16)):
        arrays['decoded_'+law] = np.zeros((count, 512, 8), complex)
        arrays['decoded_'+law][0, :8] = np.diag(np.arange(start, start+8))
    incoming = np.zeros(64, complex); incoming[10] = 1
    out = m.contract_output(arrays, 0, ((0, 0),)*4, (0,)*4, 0, incoming)
    expected = np.zeros(4096, complex)
    expected[512+16] = 2*11*3*11
    assert np.array_equal(out, expected)


def test_all_tensor_axes_against_independent_complex_einsum_with_reference():
    rng = np.random.default_rng(1129)
    def unitary(n):
        return np.linalg.qr(rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)))[0]
    U, G, V, P, F = unitary(8), unitary(64), unitary(64), unitary(64), unitary(64)
    arrays = {'controls': U[None], 'geometric_preparation': G, 'geometric_inverse': V,
        'pointer': np.asarray([P, P]), 'geometric_pointer': np.asarray([F, F]),
        'decoded_C': np.zeros((8, 512, 8), complex), 'decoded_E': np.zeros((16, 512, 8), complex)}
    factors = [unitary(8) for _ in range(4)]
    for i, (law, factor) in enumerate(zip(m.LAWS, factors, strict=True)):
        arrays['decoded_'+law][0, 8*i:8*(i+1)] = factor
    state = rng.normal(size=(64, 2))+1j*rng.normal(size=(64, 2)); state /= np.linalg.norm(state)
    actual = m.contract_output(arrays, 0, ((0, 0),)*4, (0, 1, 2, 3), 0, state)
    A = G[:, 8*np.arange(8)]
    tensor = np.einsum('ai,ijz,bj->abz', A, (np.kron(U, U)@state).reshape(8, 8, 2), A, optimize=True).reshape(8, 8, 8, 8, 2)
    tensor = np.einsum('pa,qb,rc,sd,abcdz->pqrsz', *factors, tensor, optimize=True)
    tensor = np.einsum('pqab,rscd,abcdz->pqrsz', V.reshape(8, 8, 8, 8), V.reshape(8, 8, 8, 8), tensor, optimize=True)
    tensor = np.einsum('rsbd,abcdz->arcsz', F.reshape(8, 8, 8, 8), tensor, optimize=True)
    expected = np.einsum('prac,abcdz->pbrdz', P.reshape(8, 8, 8, 8), tensor, optimize=True).reshape(4096, 2)
    assert np.linalg.norm(actual-expected) < 3e-14


@pytest.mark.parametrize('kind', ['identity_control', 'wrong_pointer', 'nonfinite', 'labels'])
def test_wrong_but_unitary_or_complete_instruments_fail_independent_target_gates(executed, kind):
    _, arrays = executed
    values = {key: value.copy() for key, value in arrays.items()}
    if kind == 'identity_control':
        values['controls'][15] = np.eye(8)
    elif kind == 'wrong_pointer':
        values['pointer'] = values['pointer'][::-1].copy()
    elif kind == 'nonfinite':
        values['geometric_inverse'][0, 0] = np.nan
    else:
        values['labels_E'][-1] = values['labels_E'][0]
    with pytest.raises(ValueError):
        m.metrics(values)


@pytest.mark.parametrize('control,records,syndromes,pointer', [
    (True, ((0, 0),)*4, (0,)*4, 0), (0, ((False, 0),)*4, (0,)*4, 0),
    (0, ((0, 0),)*4, (-1, 0, 0, 0), 0), (0, ((0, 0),)*4, (0,)*4, True)])
def test_tensor_contraction_rejects_malformed_original_indices(executed, control, records, syndromes, pointer):
    with pytest.raises(ValueError):
        m.contract_output(executed[1], control, records, syndromes, pointer)


def test_fresh_replay_checks_all_arrays_not_survival_only(executed):
    report, arrays = executed
    changed = deepcopy(report)
    changed['scope']['reference_and_class_sum_are_supplied'] = False
    with pytest.raises(ValueError):
        m.compare_execution(changed, arrays, report, arrays)
    m.compare_execution(report, arrays, report, arrays)


@pytest.mark.parametrize('tamper', ['zero_diagnostics', 'zero_budget', 'rescale_budget', 'reinsert_resource_residual'])
def test_coordinated_small_error_tampering_is_rejected(executed, tamper):
    report, arrays = executed
    changed = deepcopy(report)
    if tamper == 'zero_diagnostics':
        for key, value in changed['execution'].items():
            if type(value) is float and ('error' in key or 'defect' in key):
                changed['execution'][key] = 0.
    elif tamper == 'reinsert_resource_residual':
        changed['resources']['pointer']['maximum_numerical_valid_code_leakage'] = 0.
    else:
        scale = 0. if tamper == 'zero_budget' else .5
        for key in ('fixed_admission_successful_block_error', 'raw_macro_adjustment',
                    'successful_history_including_raw_macro_error', 'complete_flagged_operator_error',
                    'complete_flagged_channel_error'):
            changed['error_budget'][key] *= scale
    with pytest.raises(ValueError):
        m.compare_execution(changed, arrays, report, arrays)


def test_proper_code_export_failure_is_charged_and_large_leak_rejected(monkeypatch):
    original = m.compiler.native_pulse
    def leaky(kind, angle):
        value = np.array(original(kind, angle), copy=True)
        value[5, 0] += 1e-6
        return value
    monkeypatch.setattr(m.compiler, 'native_pulse', leaky)
    with pytest.raises(ArithmeticError, match='complement'):
        m._export_diagnostic([[m.compiler.pulse('X', .125, [0])]], 2)
