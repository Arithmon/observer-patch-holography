"""Source-executed replacement, retained unknown inputs and cooling limits."""
from copy import deepcopy
from fractions import Fraction as F
import json

import numpy as np
import pytest
import sympy as s

from . import native_response_cooling as c


@pytest.fixture(scope='module')
def execution():
    return c.execute()


def test_native_clock_has_complete_raw_menu_and_actual_one_twelfth_trigger():
    records = c.raw_records()
    assert [(coin, raw) for coin, raw, _ in records] == [(a, b) for a in (0, 1) for b in range(4)]
    assert sum(row[2] for row in records) == 1
    assert next(w for a, b, w in records if (a, b) == c.trigger_record()) == F(1, 12)
    assert c.clock_data()['records'] == records


def test_executed_swap_keeps_all_unknown_reference_columns_and_old_information(execution):
    report, arrays = execution
    ground = c.ground_vector().astype(complex)
    target = np.kron(np.eye(64), ground[:, None])
    assert np.linalg.norm(arrays['replaced']-target, 2) < c.SOURCE_ADMISSION
    assert np.linalg.norm(arrays['unchanged']-np.kron(ground[:, None], np.eye(64)), 2) < c.SOURCE_ADMISSION
    assert len(c.swap_word()) == 54
    assert report['resources']['SWAP']['native_pulses'] == 54
    assert report['resources']['retained_quantum_bath_and_counter_qubits_per_cycle'] == 15
    assert not report['resources']['old_reference_information_destroyed_or_reset_in_place']


def test_unknown_complex_entanglement_moves_to_bath_without_measurement(execution):
    _, arrays = execution
    columns = np.zeros((64, 3), complex)
    columns[1, 0] = np.sqrt(.2)
    columns[4, 1] = 1j*np.sqrt(.3)
    columns[17, 2] = (1-1j)/2
    output = (arrays['replaced']@columns).reshape(64, 64, 3)
    ground = c.ground_vector()
    expected = np.einsum('br,w->bwr', columns, ground)
    assert np.linalg.norm(output-expected) < 2*c.SOURCE_ADMISSION
    # Tracing the bath transfers correlations there; no input component or
    # imaginary phase was discarded by a normalized selected-branch update.
    reduced = np.einsum('bwr,bvs->wrvs', output, output.conj()).reshape(192, 192)
    ancilla = columns.T@columns.conj()
    assert np.linalg.norm(reduced-np.kron(np.outer(ground, ground), ancilla)) < 3*c.SOURCE_ADMISSION


def test_all_counter_branches_normalize_and_source_blocks_reconstruct(execution):
    report, arrays = execution
    first = np.eye(64, dtype=complex)[:, :2]
    total = np.zeros((64, 64), complex)
    for coin, raw, _ in c.raw_records():
        row = 4*coin+raw
        factor = arrays['replaced'] if (coin, raw) == c.TRIGGER else arrays['unchanged']
        weight = np.vdot(arrays['counter'][row], arrays['counter'][row]).real
        total += weight*(factor.conj().T@factor)
        counter_index = int(np.argmax(abs(arrays['counter'][row])))
        direct = arrays['counter'][row, counter_index]*(factor@first)
        assert np.array_equal(c.output_block(arrays, coin, raw, counter_index, first), direct)
    assert np.linalg.norm(total-np.eye(64), 2) < 6*c.SOURCE_ADMISSION
    assert abs(report['execution']['source_trigger_probability']-1/12) < c.SOURCE_ADMISSION


def test_unconditional_cooling_formula_with_untouched_reference_and_extreme_input():
    r = 2
    state = np.zeros((64, r), complex)
    state[7, 0], state[11, 1] = 1/np.sqrt(2), 1j/np.sqrt(2)
    rho = np.outer(state.ravel(), state.ravel().conj())
    one = c.ideal_marginal(rho, 1, reference_dimension=r)
    two = c.ideal_marginal(one, 1, reference_dimension=r)
    assert np.linalg.norm(two-c.ideal_marginal(rho, 2, reference_dimension=r)) < 3e-16
    target = np.kron(np.outer(c.ground_vector(), c.ground_vector()), np.eye(2)/2)
    q = float(c.remaining_weight(64))
    actual = c.ideal_marginal(rho, 64, reference_dimension=r)
    assert np.linalg.norm(actual-(q*rho+(1-q)*target)) < 3e-16
    # These inputs are orthogonal to the six diagonal ground-support slots,
    # so the trace-distance contraction bound is exactly saturated.
    assert abs(np.linalg.norm(actual-target, ord='nuc')/2-q) < 2e-15


def test_no_go_and_zeno_controls_act_on_the_data_not_only_a_reference():
    report = c.exact_certificate()
    x, z = s.Matrix([[0, 1], [1, 0]]), s.diag(1, -1)
    mixed = s.eye(2)/2
    assert sum((w*u*mixed*u.conjugate().T for w, u in ((s.Rational(1, 12), x),
        (s.Rational(11, 12), z))), s.zeros(2)) == mixed
    p = s.diag(1, 0)
    P = s.kronecker_product(s.eye(2), p)
    H = s.kronecker_product(x, s.eye(2))+s.kronecker_product(s.eye(2), x)
    assert P*H*P == s.kronecker_product(x, p)
    assert report['monitoring_counterexample']['all_data_currents_survive_reference_only_Zeno_limit']
    assert not report['scope']['cooling_alone_selects_data_current_algebra']


@pytest.mark.parametrize('bad', [True, False, -1, 4097, 1., np.int64(1), '1'])
def test_cycle_count_preserves_original_exact_domain(bad):
    with pytest.raises(ValueError):
        c.remaining_weight(bad)


@pytest.mark.parametrize('mutation', ['phase', 'label', 'missing', 'not_finite', 'dishonest_diagnostic', 'type'])
def test_retained_source_or_report_forgery_is_rejected(execution, mutation):
    report, arrays = execution
    bad_report = deepcopy(report)
    bad = {name: value.copy() for name, value in arrays.items()}
    if mutation == 'phase':
        bad['replaced'][:, 7] *= 1j
    elif mutation == 'label':
        bad['labels'][1] = (1, 1)
    elif mutation == 'missing':
        del bad['counter']
    elif mutation == 'not_finite':
        bad['unchanged'][0, 0] = np.nan
    elif mutation == 'dishonest_diagnostic':
        bad_report['execution']['source_trigger_probability'] = .5
    else:
        bad_report['exact']['native_clock']['trigger'][1] = True
    with pytest.raises(ValueError):
        c.compare_execution(bad_report, bad, report, arrays)


def test_original_input_array_subclass_cannot_hide_wrong_reference_columns(execution):
    _, arrays = execution
    class Deceptive(np.ndarray):
        pass
    incoming = np.eye(64, dtype=complex).view(Deceptive)
    with pytest.raises(ValueError, match='floating64-row'):
        c.output_block(arrays, 0, 1, 0, incoming)


def test_complete_report_roundtrips_and_keeps_added_architecture_explicit(execution):
    report, arrays = execution
    copy = json.loads(json.dumps(report, allow_nan=False))
    c.compare_execution(copy, arrays, report, arrays)
    assert copy['exact']['scope']['new_bath_coupling_is_a_supplied_architecture']
    assert not copy['exact']['scope']['passive_C_E_already_cools_or_selects_currents']
    assert copy['error_budget']['numerical_not_outward_rounded_or_physical']
