"""Source capacity, reference convention and full finite-history controls."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from . import construction as source
from . import thermal_source as thermal
from .thermal_forms import fixed_code_form


@pytest.fixture(scope='module')
def cold_preparation():
    lattice = source.ProductBand(occupation=8)
    return source.product_gibbs_preparation(lattice, F(4))


def test_original_four_occupation_preparation_keeps_identical_word_and_reference():
    lattice = source.ProductBand()
    raw, target, metadata = source.product_gibbs_preparation(lattice)
    values, basis = source.inherited_source._basis(lattice.hamiltonian.toarray())
    weights = np.exp(-(values-values[0])/4); weights /= weights.sum()
    old_word = source.compiler.compile_real_state(np.sqrt(weights))
    for bit in range(4):
        old_word += source.compiler.cnot(bit, bit+4)
    old_word += source.compiler.compile_real_orthogonal(basis)
    blank = np.zeros(256, complex); blank[0] = 1
    expected = source.compiler.execute(old_word, 8, blank).reshape(16, 16).T
    assert np.array_equal(raw, expected)
    assert np.array_equal(target, basis@np.diag(np.sqrt(weights)))
    assert metadata['source_resources'] == source.compiler.resources(old_word)
    assert 'reference_basis' not in metadata  # Preserve the old receipt schema.


def test_two_occupation_preparation_has_the_correct_copied_register_bits():
    lattice = source.ProductBand(occupation=2)
    raw, target, metadata = source.product_gibbs_preparation(lattice, F(1))
    density = expm(-fixed_code_form(2).orthonormal_matrix().toarray())
    density /= np.trace(density)
    assert raw.shape == target.shape == (4, 4)
    np.testing.assert_allclose(raw@raw.conj().T, density, atol=1e-12, rtol=0)
    assert metadata['source_resources']['maximum_work_qubits'] == 4


def test_cold_eight_occupation_source_uses_canonical_reference(cold_preparation):
    raw, target, metadata = cold_preparation
    hamiltonian = fixed_code_form(8).orthonormal_matrix().toarray()
    square_root = expm(-2*hamiltonian)
    square_root /= np.linalg.norm(square_root)
    assert raw.shape == target.shape == (64, 64)
    np.testing.assert_allclose(target, square_root, atol=2e-14, rtol=0)
    np.testing.assert_allclose(raw, square_root, atol=1e-9, rtol=0)
    assert metadata['reference_basis'] == 'canonical_fixed_frequency_two_Hermite'
    assert metadata['source_resources']['maximum_work_qubits'] == 12
    assert metadata['source_resources']['native_pulses'] == 135171
    assert not metadata['physical_preparation_error_certified']


def test_cold_source_complete_histories_against_independent_dense_oracle(cold_preparation):
    raw, target, _ = cold_preparation
    lattice = source.ProductBand(occupation=8)
    interval, steps = F(1, 16), 64
    split = source.splitting_bound(lattice, interval, steps)
    assert split['two_interval_complete_output_trace_distance_upper'] < F(1, 200)
    for sign in (-1, 1):
        rows, actual = source.history(lattice, sign, initial=raw, interval=interval,
                                      steps=steps, return_states=True)
        expected = thermal.independent_history(8, F(1), target, sign, interval)
        assert [row['bits'] for row in rows] == [[0, 0], [0, 1], [1, 0], [1, 1]]
        assert all(state.shape == (64, 64) for state in actual)
        assert sum(row['probability'] for row in rows) == pytest.approx(1, abs=1e-8)
        assert thermal.pure_record_distance(actual, expected) < 1e-6
        # A detector-only law would miss this deliberately changed retained
        # reference. Our comparison includes that reference in each block.
        changed_reference = [state[:, ::-1] for state in expected]
        assert thermal.pure_record_distance(actual, changed_reference) > .9


def test_complete_distance_keeps_zero_branches_and_unnormalized_mass():
    x = np.array([1., 0.]); y = np.array([0., 1.]); zero = x*0
    assert thermal.pure_record_distance((x/2, zero), (y/2, zero)) == pytest.approx(.25)
    assert thermal.pure_record_distance((x, zero), (zero, x)) == pytest.approx(1)
    with pytest.raises(ValueError):
        thermal.pure_record_distance((x,), (x, zero))


@pytest.mark.parametrize('occupation', (3, 6, 16))
def test_source_preparation_rejects_undeclared_dimensions(occupation):
    with pytest.raises(ValueError, match='occupation two, four or eight'):
        source.product_gibbs_preparation(source.ProductBand(occupation=occupation))


@pytest.mark.parametrize('bad', (np.array([True, False]), np.ma.array([1., 0.]),
                               np.array([np.nan, 0.]), np.array([np.inf, 0.]),
                               np.array([1., 0.], dtype=object), [1., 0.],
                               np.zeros(4097), np.array(1.)))
def test_complete_distance_rejects_hostile_amplitudes(bad):
    with pytest.raises(ValueError):
        thermal.pure_record_distance((bad,), (np.array([1., 0.]),))


def test_complete_distance_rejects_overflow_in_finite_inputs_and_intermediates(monkeypatch):
    huge = np.array([1e308, 0.])
    with pytest.raises(ArithmeticError):
        thermal.pure_record_distance((huge,), (huge,))
    ordinary = np.array([1., 0.])
    monkeypatch.setattr(np.linalg, 'qr', lambda *_args, **_kwargs:
                        (np.eye(2), np.full((2, 2), np.nan)))
    with pytest.raises(ArithmeticError, match='QR'):
        thermal.pure_record_distance((ordinary,), (ordinary,))


def test_complete_distance_rejects_nonfinite_eigenvalue_result(monkeypatch):
    ordinary = np.array([1., 0.])
    monkeypatch.setattr(np.linalg, 'eigvalsh', lambda _matrix: np.array([np.nan]))
    with pytest.raises(ArithmeticError, match='norm'):
        thermal.pure_record_distance((ordinary,), (ordinary,))


@pytest.mark.parametrize('kwargs', ({'occupation': True}, {'occupation': np.int64(8)},
                                  {'occupation': 16}, {'beta': True}, {'beta': 4.},
                                  {'beta': F(17)}, {'beta': F(0)},
                                  {'interval': .0625}, {'interval': F(2)},
                                  {'steps': True}, {'steps': 64.}, {'steps': 4097},
                                  {'include_target': 1}, {'include_target': np.bool_(True)}))
def test_run_rejects_parameters_before_any_source_execution(monkeypatch, kwargs):
    def forbidden(*_args, **_kwargs):
        raise AssertionError('invalid input reached source preparation')
    monkeypatch.setattr(source, 'product_gibbs_preparation', forbidden)
    with pytest.raises(ValueError):
        thermal.run(**kwargs)


@pytest.mark.parametrize('corruption', ('state_nan', 'probability_nan', 'probability_bool',
                                      'probability_wrong_mass', 'record_bool', 'lost_branch'))
def test_run_rejects_corrupted_complete_source_branches(monkeypatch, corruption):
    raw, target, metadata = source.product_gibbs_preparation(source.ProductBand(occupation=2))
    monkeypatch.setattr(source, 'product_gibbs_preparation', lambda *_args: (raw, target, metadata))
    states = [target.copy()/2 for _ in range(4)]
    law = [{'bits': [first, last], 'probability': .25}
           for first in (0, 1) for last in (0, 1)]
    if corruption == 'state_nan':
        states[2][0, 0] = np.nan
    elif corruption == 'probability_nan':
        law[2]['probability'] = float('nan')
    elif corruption == 'probability_bool':
        law[2]['probability'] = True
    elif corruption == 'probability_wrong_mass':
        law[2]['probability'] = .5
    elif corruption == 'record_bool':
        law[2]['bits'] = [True, 0]
    else:
        law.pop(); states.pop()
    monkeypatch.setattr(source, 'history', lambda *_args, **_kwargs: (law, states))
    with pytest.raises((ValueError, ArithmeticError)):
        thermal.run(2)


def test_run_rejects_nonfinite_independent_branch_oracle(monkeypatch):
    raw, target, metadata = source.product_gibbs_preparation(source.ProductBand(occupation=2))
    monkeypatch.setattr(source, 'product_gibbs_preparation', lambda *_args: (raw, target, metadata))
    states = [target.copy()/2 for _ in range(4)]
    law = [{'bits': [first, last], 'probability': .25}
           for first in (0, 1) for last in (0, 1)]
    monkeypatch.setattr(source, 'history', lambda *_args, **_kwargs: (law, states))
    corrupted = [state.copy() for state in states]
    corrupted[0][0, 0] = np.nan
    monkeypatch.setattr(thermal, 'independent_history', lambda *_args: corrupted)
    with pytest.raises(ValueError):
        thermal.run(2)


def test_optional_target_is_the_actual_preparation_target_and_report_stays_serializable():
    import json
    report, target = thermal.run(2, F(1), include_target=True)
    _, expected, _ = source.product_gibbs_preparation(source.ProductBand(occupation=2), F(1))
    assert np.array_equal(target, expected)
    assert target.shape == (4, 4)
    assert 0 < report['target_norm_upper'] < F(1000001, 1000000)
    assert report['target_scaled_splitting_error_upper'] > 0
    json.dumps(source.a.serializable(report), allow_nan=False)


@pytest.mark.parametrize('scale', (1, 2))
def test_source_budget_charges_known_nonnormalized_target_and_density_scaling(scale):
    # A 3-4 vector has exact norm five. A scalar error on every local
    # operation saturates a computable pure-density amplification, exposing
    # a budget that silently replaces the target's norm by one.
    target = scale*np.array([[3., 4.], [0., 0.]])
    report = thermal._source_budget(target, F(0), 3, F(1, 8))
    norm_squared = F(25*scale*scale)
    scalar_gain = (1+F(251, 10**12))**3
    exact_density_error = norm_squared*(scalar_gain**2-1)/2
    assert report['target_norm_upper'] >= 5*scale
    assert report['source_numerical_guard'] >= exact_density_error
    assert report['target_scaled_splitting_error_upper'] >= norm_squared/8
    normalized = thermal._source_budget(np.array([[1., 0.], [0., 0.]]), F(0), 3, F(1, 8))
    assert report['source_numerical_guard'] > 24*scale*scale*normalized['source_numerical_guard']


@pytest.mark.parametrize('bad', (np.full((2, 2), np.inf), np.full((2, 2), np.nan),
                               np.full((2, 2), 1e308), np.ones((2, 2), dtype=bool),
                               np.ma.array(np.eye(2))))
def test_source_budget_rejects_nonfinite_or_overflowing_normalization(bad):
    with pytest.raises((ValueError, ArithmeticError)):
        thermal._source_budget(bad, F(1, 10**9), 772, F(1, 100))
