"""Complete native labels, dynamic work codes and signed-kick controls."""
from dataclasses import replace
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from model.interacting_band_extension import construction as source
from model.interacting_band_extension import thermal_controller as controller


def target(occupation):
    # The controller's all-input comparison is independent of preparation.
    # Production integration passes the existing executed raw purification;
    # this test uses an independent canonical state without replaying it.
    hamiltonian = source.ProductBand(occupation=occupation).hamiltonian.toarray()
    value = expm(-2*hamiltonian)
    return value/np.linalg.norm(value)


def test_cold_controller_retains_all_columns_native_labels_and_timeout_prefixes():
    initial = target(8)
    row = controller.audit(initial, initial)
    assert row['retained_work_input_columns'] == 64
    assert row['retained_preparation_reference_columns'] == 64
    assert row['complete_native_tapes'] == 3600
    assert row['complete_native_field_outputs'] == 4900
    assert row['maximum_active_instructions'] == 775
    assert row['all_input_source_Kraus_comparison_residual_numerical'] < 2e-12
    assert row['all_input_channel_completeness_residual_numerical'] < 2e-11
    assert row['complete_prepared_work_reference_prefix_distance_numerical'] < 1e-7
    prefixes = row['complete_clock_field_prefix_law']
    assert sum(x['probability'] for x in prefixes) == pytest.approx(1, abs=1e-9)
    for tick, weight in enumerate((F(25, 36), F(5, 18), F(1, 36))):
        selected = [x for x in prefixes if x['incident_ticks'] == tick]
        assert len(selected) == 2**tick
        assert all(x['native_tape_weight'] == weight for x in selected)
    assert row['timeout_probability'] == F(35, 36)
    assert not row['source_preparation_executed_again']
    assert not row['new_d8_quantum_access_algebra_rank_claimed']
    assert row['valid_code_service_failure_probability_exact'] == 0
    assert not row['nonzero_service_failure_on_actual_preparation_executed']
    assert not row['arbitrary_unfinished_mid_service_output_bound_claimed']


@pytest.mark.parametrize('sign', (-1, 1))
@pytest.mark.parametrize('occupation', (2, 8))
def test_source_prefix_words_match_all_controller_columns_for_both_signs(occupation, sign):
    program = source.ProductControllerInstrument(occupation=occupation, sign=sign, steps=1)
    identity = np.eye(occupation**2)
    direct = source._product_clock_prefixes(program.lattice, identity, source=True,
                                            interval=program.interval, steps=1, sign=sign)
    (bits, kick, _), = program.service_closure('kick', ())
    assert bits == ()
    previous = {(): kick}
    np.testing.assert_array_equal(kick, direct[0][()])
    for tick in (1, 2):
        following = {}
        for prefix, operation in previous.items():
            for bits, local, _ in program.service_closure('evolve', prefix):
                following[bits] = local@operation
                np.testing.assert_allclose(following[bits], direct[tick][bits], atol=2e-12, rtol=0)
        previous = following
    q = source.normalized(source.a.polynomial(occupation, 1))
    np.testing.assert_allclose(program.kick, expm(1j*sign*q/4), atol=2e-10, rtol=0)
    np.testing.assert_allclose(program.feedback, expm(1j*q/4), atol=2e-10, rtol=0)
    if sign == -1:
        assert np.linalg.norm(program.kick-program.feedback) > .1


@pytest.mark.parametrize('kwargs', ({'occupation': True}, {'occupation': np.int64(8)},
                                  {'occupation': 16}, {'sign': True}, {'sign': np.int64(1)},
                                  {'sign': 0}, {'sign': '-1'}))
def test_controller_rejects_nonoriginal_settings(kwargs):
    with pytest.raises(ValueError):
        source.ProductControllerInstrument(**kwargs)


@pytest.mark.parametrize('corruption', ('duplicate', 'drop'))
def test_complete_audit_rejects_duplicate_or_dropped_native_field_labels(monkeypatch, corruption):
    completed = source.ProductControllerInstrument.completed
    def corrupted(self, *args, **kwargs):
        rows = list(completed(self, *args, **kwargs))
        if corruption == 'duplicate':
            rows[1] = rows[0]  # Same length, but one distinct output is lost.
        else:
            rows.pop()
        return tuple(rows)
    monkeypatch.setattr(source.ProductControllerInstrument, 'completed', corrupted)
    initial = target(2)
    with pytest.raises(ArithmeticError, match='label|all-input'):
        controller.audit(initial, initial, occupation=2, steps=1)


def test_larger_controller_rejects_erasing_a_committed_field_prefix(monkeypatch):
    transition = source.product_controller_transition
    def erase(state, outcome=None, *, steps=64):
        following = transition(state, outcome, steps=steps)
        if state.instruction == 'feedback' and outcome is None:
            return replace(following, bits=(0,))
        return following
    monkeypatch.setattr(source, 'product_controller_transition', erase)
    program = source.ProductControllerInstrument(occupation=8, steps=1)
    with pytest.raises(ArithmeticError, match='committed prefix'):
        program.service_closure('evolve', ())


def test_controller_rejects_undercharged_instruction_counts(monkeypatch):
    completed = source.ProductControllerInstrument.completed
    def undercharged(self, *args, **kwargs):
        return tuple((state, operation, weight, count-1)
                     for state, operation, weight, count in completed(self, *args, **kwargs))
    monkeypatch.setattr(source.ProductControllerInstrument, 'completed', undercharged)
    initial = target(2)
    with pytest.raises(ArithmeticError, match='instruction budget'):
        controller.audit(initial, initial, occupation=2, steps=1)
