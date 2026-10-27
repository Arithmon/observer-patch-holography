"""Complete clock/control histories, independently derived outputs and mutations."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import native_response_controller as m


@pytest.fixture(scope='module')
def executed():
    ref, ra = m.reference.execute()
    cool, ca = m.cooling.execute()
    report, arrays = m.execute(ra, ca, reference_report=ref, cooling_report=cool)
    return report, arrays, ra, ca


def test_all_raw_clock_histories_counted_without_conditioning():
    rows = m.cooling.raw_records()
    expected = [F(0)]*3
    for history in product(rows, repeat=3):
        records = tuple(row[:2] for row in history)
        expected[m.terminal_count(records, 2)] += math.prod(row[2] for row in history)
    assert expected == m.distribution(3, 2)
    assert sum(expected) == 1
    assert expected[0] == F(11, 12)**3
    assert m.terminal_count([(0, 1)]*19) == 16
    certificate = m.exact_certificate()
    assert 4.28e-7 < float(F(certificate['timeout_probability_exact'])) < 4.29e-7
    assert not certificate['finite_exact_current_dimension_is_twelve']
    assert certificate['finite_feedback_remains_full_rank']


@pytest.mark.parametrize('bad', [True, -1, 1.0, 4097])
def test_invalid_clock_domain_fails_closed(bad):
    with pytest.raises(ValueError):
        m.distribution(bad)


@pytest.mark.parametrize('bad', [[(True, 1)], [(0, 4)], [[0, 1]], [(0., 1)]])
def test_invalid_raw_record_not_reinterpreted(bad):
    with pytest.raises(ValueError):
        m.terminal_count(bad)


def test_feedback_is_full_Gram_recurrence_without_target_projection():
    gram = m._source_data()[1]
    rng = np.random.default_rng(824)
    vector = rng.normal(size=18)
    result = vector.copy()
    for _ in range(m.CAP):
        result = m.update(result)
    exact = np.asarray(m.dynamics.feedback_matrix(m._source_data()[0], m.CAP).evalf(30), float)
    assert np.linalg.norm(result-exact@vector) < 3e-14
    assert np.linalg.norm(m.update(vector)-(vector-gram@vector/2)) < 1e-15
    assert np.linalg.matrix_rank(exact, tol=1e-12) == 18
    assert np.linalg.norm(exact[:, 12]) < 2.4e-8
    # With an unbounded input budget the attenuated direction is still usable.
    amplified = np.eye(18)[:, 12]/exact[12, 12]
    assert np.linalg.norm(exact@amplified-np.eye(18)[:, 12]) < 1e-12


def test_complete_source_controls_and_whole_history_budget(executed):
    report, arrays, ra, ca = executed
    assert arrays['controls'].shape == (17, 18, 8, 8)
    assert max(np.ravel(report['execution']['complete_control_operator_errors'])) < 1e-9
    assert report['execution']['all_clock_records_and_timeouts_retained']
    assert report['error_budget']['ideal_complete_history_channel_error'] < 1e-6
    assert report['error_budget']['complete_source_history_vs_limiting_control_channel_error'] < .01
    assert report['resources']['raw_clock_record_bits'] == 1536
    assert report['resources']['retained_bath_and_counter_qubits'] == 7680
    assert np.linalg.norm(arrays['controls'][0, 12]-np.eye(8)) > .1
    assert np.linalg.norm(arrays['controls'][16, 12]-np.eye(8)) < 4e-9
    assert np.linalg.norm(arrays['controls'][16, 3]-np.eye(8)) > .1
    m.compare_execution(report, arrays, report, arrays, reference_arrays=ra,
        cooling_arrays=ca, fresh_reference_arrays=ra, fresh_cooling_arrays=ca)


def test_clock_probabilities_come_from_actual_parent_factors(executed):
    report, arrays, ra, ca = executed
    changed = {key: value.copy() for key, value in ca.items()}
    changed['counter'][1] *= 1.001
    with pytest.raises(ValueError):
        m.metrics(arrays, ra, changed)


def test_unselected_parent_bath_column_is_validated(executed):
    _, arrays, ra, ca = executed
    changed = {key: value.copy() for key, value in ca.items()}
    changed['replaced'][:, 7] *= 1j
    with pytest.raises(ValueError):
        m.metrics(arrays, ra, changed)


@pytest.mark.parametrize('kind', ['phase', 'identity', 'coefficient', 'dropped_padding', 'nan'])
def test_controller_trash_or_unitary_wrong_answers_rejected(executed, kind):
    report, arrays, ra, ca = executed
    bad = {key: value.copy() for key, value in arrays.items()}
    if kind == 'phase':
        bad['controls'][16, 3] *= np.exp(.02j)
    elif kind == 'identity':
        bad['controls'][:] = np.eye(8)
    elif kind == 'coefficient':
        bad['coefficients'][4, 3, 3] *= -1
    elif kind == 'dropped_padding':
        bad['controls'] = bad['controls'][:, :, :6, :6]
    else:
        bad['controls'][0, 0, 0, 0] = np.nan
    with pytest.raises(ValueError):
        m.metrics(bad, ra, ca)


def test_non_scalar_tensor_oracle_catches_wire_reversal(executed):
    _, arrays, ra, _ = executed
    rng = np.random.default_rng(321)
    changed = {key: value.copy() for key, value in ra.items()}
    A = rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8))
    B = rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8))
    changed['decoded_C'][1, 8:16] = A
    changed['decoded_E'][10, 16:24] = B
    incoming = rng.normal(size=(8, 3))+1j*rng.normal(size=(8, 3))
    actual = m.work_output(arrays, changed, 7, 3, ((0, 1), (1, 2)), (1, 2), incoming)
    expected = changed['geometric_inverse']@np.kron(A, B)@changed['geometric_preparation'][:, 8*np.arange(8)]@arrays['controls'][7, 3]@incoming
    assert np.linalg.norm(actual-expected) < 2e-12
    reversed_output = changed['geometric_inverse']@np.kron(B, A)@changed['geometric_preparation'][:, 8*np.arange(8)]@arrays['controls'][7, 3]@incoming
    assert np.linalg.norm(actual-reversed_output) > 1


def test_actual_data_outputs_include_abort_and_high_syndromes(executed):
    _, arrays, ra, _ = executed
    incoming = np.eye(8, 6, dtype=complex)
    for records, syndromes in ((((0, 3), (0, 1)), (0, 0)), (((0, 1), (0, 1)), (1, 0))):
        out = m.work_output(arrays, ra, 16, 3, records, syndromes, incoming)
        assert out.shape == (64, 6)
        assert np.linalg.norm(out) < 1e-8


def test_sequential_baths_and_data_share_the_same_records_and_entangled_input(executed):
    _, arrays, ra, ca = executed
    rng = np.random.default_rng(129)
    changed = {key: value.copy() for key, value in ca.items()}
    # Nonscalar complex blocks make swaps of reference/bath axes visible.
    A = (rng.normal(size=(64, 64))+1j*rng.normal(size=(64, 64)))/64
    B = (rng.normal(size=(64, 64))+1j*rng.normal(size=(64, 64)))/64
    changed['replaced'].reshape(64, 64, 64)[4] = A
    changed['unchanged'].reshape(64, 64, 64)[5] = B
    changed['counter'][1, 7] = .3+.2j
    changed['counter'][0, 9] = .4-.1j
    incoming = rng.normal(size=(64, 8, 2))+1j*rng.normal(size=(64, 8, 2))
    work_records = tuple(next((c, r) for c in (0, 1) for r in range(4 if law == 'C' else 8)
        if m.transport.PAIR_CODES[law].branch_weight(law, m.transport.SEAMS[law], c, r)) for law in ('C', 'E'))
    actual = m.history_output_block(arrays, ra, changed,
        clock_records=((0, 1), (0, 0)), counter_rows=(7, 9), bath_rows=(4, 5),
        primitive=12, work_records=work_records, work_syndromes=(0, 0), incoming=incoming)
    data = m.work_output(arrays, ra, 1, 12, work_records, (0, 0), np.eye(8, dtype=complex))
    expected = ((.3+.2j)*(.4-.1j)*np.kron(B@A, data)@incoming.reshape(512, 2)).reshape(64, 64, 2)
    assert np.linalg.norm(expected) > .01
    assert np.linalg.norm(actual-expected) < 2e-14
    wrong_data = m.work_output(arrays, ra, 0, 12, work_records, (0, 0), np.eye(8, dtype=complex))
    wrong = ((.3+.2j)*(.4-.1j)*np.kron(B@A, wrong_data)@incoming.reshape(512, 2)).reshape(64, 64, 2)
    assert np.linalg.norm(actual-wrong) > 1e-4


def test_bad_claim_and_deleted_tiny_budget_rejected(executed):
    report, arrays, ra, ca = executed
    for name in ('claim', 'budget'):
        bad = deepcopy(report)
        if name == 'claim':
            bad['exact']['finite_exact_current_dimension_is_twelve'] = True
        else:
            bad['error_budget']['ideal_complete_history_channel_error'] = 0.
        with pytest.raises(ValueError):
            m.compare_execution(bad, arrays, report, arrays, reference_arrays=ra,
                cooling_arrays=ca, fresh_reference_arrays=ra, fresh_cooling_arrays=ca)
