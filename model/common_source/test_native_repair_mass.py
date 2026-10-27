"""Source law, complete record channel, correlated limit and custody checks."""
from copy import deepcopy

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import native_repair_mass as model


@pytest.fixture(scope='module')
def arrays():
    return model.execute()


def test_actual_event_and_all_raw_outcomes():
    rows = model.native_records()
    assert len(rows) == 128
    assert len({row[:4] for row in rows}) == 128
    assert sum(row[-1] for row in rows) == 1
    assert sum(row[-1] for row in rows if row[-2]) == s.Rational(1, 12)
    assert any(row[-1] == 0 for row in rows)


def test_fixed_native_geometry_and_hierarchy_polynomial():
    d = model.architecture()
    assert d['rotation'] == (0, 5, 6, 3, 2, 11, 8, 1, 9, 4, 7, 10)
    assert d['seed'] == 1
    assert model._zero(d['U']**5-s.eye(3))
    assert s.trace(d['U']) == (1+s.sqrt(5))/2
    assert s.trace(d['Y']) == 1
    assert d['e2'] == s.Rational(1392, 22621)
    assert s.simplify(d['Y'].det()-d['determinant']) == 0
    x = np.linalg.eigvalsh(np.asarray(d['Y'], float))
    assert np.all(x > 0) and np.all(np.diff(x) > 0)


def test_exact_memory_stationary_and_all_transients():
    d = model.architecture(); A = d['transition']; p = d['p']
    assert A*d['stationary'] == d['stationary']
    for old in range(5):
        eta = s.eye(5)[:, old]-d['stationary']
        assert A*eta == p*d['shift']*eta
        assert A**8*eta == p**8*d['shift']**8*eta


def test_feedback_matches_full_exponential():
    X = np.asarray([[0, 1], [1, 0]], complex)
    for Q, U in zip(model.architecture()['projectors'], model.feedback_unitaries(.125), strict=True):
        H = np.kron(X, np.asarray(Q, complex))
        np.testing.assert_allclose(U, expm(-.125j*H), atol=1e-14)
        np.testing.assert_allclose(U.conj().T@U, np.eye(6), atol=1e-14)


def test_literal640_branch_completeness():
    K, labels = model.raw_kraus(.125)
    assert K.shape == (640, 30, 30) and labels.shape == (640, 6)
    assert len(set(map(tuple, labels))) == 640
    gram = np.einsum('aij,aik->jk', K.conj(), K)
    np.testing.assert_allclose(gram, np.eye(30), atol=2e-14)
    for row, (cc, cr, ec, er, old, new) in zip(K, labels, strict=True):
        assert new == ((old+1) % 5 if (cc, cr) == (0, 1) else 0)
        assert not np.any(np.delete(row, np.arange(6*new, 6*new+6), axis=0))


def test_factored_channel_equals_literal_complete_kraus():
    K, _ = model.raw_kraus(.125)
    rng = np.random.default_rng(51)
    state = rng.standard_normal((30, 7))+1j*rng.standard_normal((30, 7))
    full = state@state.conj().T; full /= np.trace(full)
    literal = np.einsum('aij,jk,alk->il', K, full, K.conj(), optimize=True)
    blocks = np.stack([full[6*k:6*k+6, 6*k:6*k+6] for k in range(5)])
    computed = (model.joint_channel(.125)@blocks.reshape(180)).reshape(5, 6, 6)
    expected = np.zeros((30, 30), complex)
    for k in range(5): expected[6*k:6*k+6, 6*k:6*k+6] = computed[k]
    np.testing.assert_allclose(literal, expected, atol=3e-15)


def test_finite_channel_choi_and_correlated_limit(arrays):
    report = model.metrics(arrays)
    rows = report['finite_steps']
    assert all(row['minimum_complete_Choi_eigenvalue_diagnostic'] > -1e-12 for row in rows)
    assert all(row['complete_Choi_Frobenius_diamond_upper_diagnostic']
               < row['analytic_diamond_bound_all_initial_classical_memories'] for row in rows)
    errors = [row['superoperator_spectral_difference_diagnostic'] for row in rows]
    assert all(3.8 < a/b < 4.1 for a, b in zip(errors[:-1], errors[1:], strict=True))
    assert model.averaging_bound(64) == pytest.approx(56/(11*64))
    assert model.averaging_bound(64, initial_memory_stationary=True) == pytest.approx(52/(11*64))


def test_stationary_product_and_initial_memory_are_not_conflated():
    d = model.architecture(); pi = np.asarray(d['stationary'], float).reshape(5)
    blocks = np.zeros((5, 6, 6), complex)
    for k in range(5): blocks[k, k, k] = pi[k]
    marginal = blocks.sum(axis=0)
    assert np.linalg.norm(blocks-pi[:, None, None]*marginal) > .1
    out = model.evolve_cq(blocks, 32)
    assert np.trace(out.sum(axis=0)).real == pytest.approx(1)
    ideal = expm(-1j*np.kron([[0, 1], [1, 0]], np.asarray(d['Y'], complex)))
    delta = out.sum(axis=0)-ideal@marginal@ideal.conj().T
    assert np.linalg.norm(delta, 'nuc') < model.averaging_bound(32)


@pytest.mark.parametrize('steps,time', [(True, 1), (0, 1), (-1, 1), (1.5, 1), (4, -1), (4, float('nan')), (4, True)])
def test_bad_clock_domains_rejected(steps, time):
    with pytest.raises(ValueError): model.averaging_bound(steps, time)


@pytest.mark.parametrize('change', ['missing', 'extra', 'shape', 'dtype', 'nonfinite'])
def test_array_schema_fails_closed(arrays, change):
    values = {k: v.copy() for k, v in arrays.items()}
    if change == 'missing': del values['raw_unitaries']
    if change == 'extra': values['helper'] = np.zeros(1)
    if change == 'shape': values['channels'] = values['channels'][:, :4]
    if change == 'dtype': values['raw_labels'] = values['raw_labels'].astype(float)
    if change == 'nonfinite': values['stationary'][0] = np.nan
    with pytest.raises(ValueError): model._schema(values)


def test_channel_corruption_and_incomplete_records_rejected(arrays):
    bad = {k: v.copy() for k, v in arrays.items()}
    bad['channels'][0, 0] *= .99
    with pytest.raises(ValueError): model.metrics(bad)


def test_diamond_gate_uses_complete_choi_certificate(arrays):
    bad = {k: v.copy() for k, v in arrays.items()}
    pulse = np.diag([np.exp(.001j), 1, 1, 1, 1, 1])
    U = pulse@arrays['averaged_unitary']
    C = np.kron(U, U.conj())
    difference = C-arrays['averaged_channel']
    bound = model.averaging_bound(model.COUNTS[-1])
    assert np.linalg.norm(difference, 2) < bound
    assert 6*np.linalg.norm(difference) > bound
    bad['channels'][-1] = C
    with pytest.raises(ValueError, match='Choi certificate'): model.metrics(bad)


def test_report_types_and_fixed_bounds_are_exact():
    with pytest.raises(ValueError): model._compare({'steps': 5}, {'steps': True})
    with pytest.raises(ValueError): model._compare({'bound': 1.}, {'bound': 1.+1e-12})
    model._compare({'error_diagnostic': 1e-12}, {'error_diagnostic': 2e-12})
    with pytest.raises(ValueError): model._compare({'error_diagnostic': 1e-12}, {'error_diagnostic': float('nan')})


def test_no_empirical_target_import_in_source():
    source = (model.ROOT/'model/common_source/native_repair_mass.py').read_text(encoding='utf-8')
    assert 'sm_calibration' not in source
    assert 'parameter_packet' not in source
    assert 'default_rng' not in source


@pytest.mark.parametrize('helper', ['family', 'dressing', 'feedback', 'cycles', 'bath', 'native_readout'])
@pytest.mark.parametrize('change', ['missing', 'shape', 'dtype', 'nonfinite'])
def test_integrated_helper_array_menu_is_not_optional_when_present(arrays, helper, change):
    from . import native_family_selection as family
    from . import native_mass_dressing as dressing
    from . import native_repair_mass_feedback as source_feedback
    from . import native_mass_record_cycles as cycles
    from . import native_mass_bath as bath
    from . import native_repair_mass_native_readout as native_readout
    schema = {'family': family.ARRAY_SCHEMA, 'dressing': dressing.SCHEMA,
              'feedback': source_feedback.ARRAY_SCHEMA, 'cycles': cycles.SCHEMA,
              'bath': bath.SCHEMA, 'native_readout': native_readout.SCHEMA}[helper]
    values = dict(arrays)
    values.update({key: np.zeros(shape, dtype=dtype) for key, (shape, dtype) in schema.items()})
    model._schema(values)
    key = next(key for key, (_, dtype) in schema.items() if dtype != 'int64')
    if change == 'missing': del values[key]
    if change == 'shape': values[key] = values[key].reshape(-1)
    if change == 'dtype': values[key] = values[key].astype(np.complex64)
    if change == 'nonfinite': values[key].flat[0] = np.nan
    with pytest.raises(ValueError): model._schema(values)


@pytest.mark.parametrize('helper', ['dressing', 'bath'])
def test_retained_pole_summary_and_analytic_bounds_both_bind(monkeypatch, helper):
    from . import native_mass_dressing as dressing
    from . import native_mass_bath as bath
    module = dressing if helper == 'dressing' else bath
    bind = model._bind_dressing if helper == 'dressing' else model._bind_bath
    source = {'charged_channels': [{'light_poles': [.002, .06, .9]}],
              'numerical_operator_admission': 2e-10,
              'scope': {'complete2250_operator_compiled_into_native_source_words': False},
              'exact_no_tuning_shape_bound': {'middle_over_smallest_ratio_strictly_below_32': True}}
    arrays = {key: object() for key in module.SCHEMA}
    seen = []
    def from_arrays(values):
        assert set(values) == set(module.SCHEMA)
        assert all(values[key] is arrays[key] for key in values)
        seen.append(True)
        return deepcopy(source)
    monkeypatch.setattr(module, 'report_from_arrays', from_arrays)
    assert bind(arrays, source) == source and seen
    for field in ('pole', 'bound', 'scope', 'theorem'):
        bad = deepcopy(source)
        if field == 'pole': bad['charged_channels'][0]['light_poles'][0] += 1e-12
        if field == 'bound': bad['numerical_operator_admission'] += 1e-12
        if field == 'scope': bad['scope']['complete2250_operator_compiled_into_native_source_words'] = True
        if field == 'theorem': bad['exact_no_tuning_shape_bound']['middle_over_smallest_ratio_strictly_below_32'] = False
        with pytest.raises(ValueError): bind(arrays, bad)


def test_only_computed_compiler_resource_floats_receive_diagnostic_tolerance():
    for name in ('integrated_processor_norm_upper', 'maximum_numerical_valid_code_leakage'):
        model._compare({'resources': {name: 1.}}, {'resources': {name: 1.+1e-12}})
    with pytest.raises(ValueError):
        model._compare({'resources': {'fixed_channel_allowance': 1.}},
                       {'resources': {'fixed_channel_allowance': 1.+1e-12}})


def test_native_inference_is_bound_to_all_own_source_factors_and_records(monkeypatch):
    from . import native_repair_mass_native_readout as native_readout
    keys = set(native_readout.BASE_SCHEMA) | set(native_readout.SCHEMA)
    arrays = {key: object() for key in keys}
    own = {'all128_pointer_labels_and_native_syndrome_outputs_retained': True,
           'complete_record_l1_reconstruction_error_diagnostic': 0.,
           'inferred_kernel': {'eigenvalue_intervals': [[.001, .003], [.06, .07], [.93, .94]],
                               'matrix_eigenvectors_enter_probes_or_inverse': False}}
    seen = []
    def metrics(values):
        assert set(values) == keys
        assert all(values[key] is arrays[key] for key in keys)
        seen.append(True)
        return deepcopy(own)
    monkeypatch.setattr(native_readout, 'metrics', metrics)
    assert model._bind_native_readout(arrays, own) == own['inferred_kernel'] and seen
    for change in ('interval', 'pointer', 'eigenbasis'):
        bad = deepcopy(own)
        if change == 'interval': bad['inferred_kernel']['eigenvalue_intervals'][0][0] += 1e-6
        if change == 'pointer': bad['all128_pointer_labels_and_native_syndrome_outputs_retained'] = False
        if change == 'eigenbasis': bad['inferred_kernel']['matrix_eigenvectors_enter_probes_or_inverse'] = True
        with pytest.raises(ValueError): model._bind_native_readout(arrays, bad)
