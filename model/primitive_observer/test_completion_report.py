"""Joined execution must preserve the complete theorem-attached observation."""
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import completion_report as report
from . import runtime


def test_joined_report_reproduces_exact_laws_and_complete_reference_history():
    result = report.run()
    assert result['all_pass'] is True
    assert result['exact_attachment']['native_matrix_units_checked'] == 4356
    assert [row['ideal_C_power_probabilities'] for row in result['runs']] == [
        ['5/8', '3/8'], ['1/2', '1/2'], ['13/32', '19/32'], ['193/512', '319/512']]
    assert [row['exact_total_variation_and_half_diamond'] for row in result['runs']] == [
        '1/4', '1/8', '1/32', '1/512']
    for row in result['runs']:
        assert row['complete_records'] == [[0], [1]]
        assert row['diagnostics']['compressed_complete_cq_distance'] < 2e-14
        assert row['diagnostics']['independent_terminal_density_residual'] < 2e-14
        assert row['diagnostics']['iterated_plus_probability'] == pytest.approx(
            float(F(row['ideal_C_power_probabilities'][0])), abs=2e-14)
    reference = result['entangled_reference_and_checkpoint']
    assert reference['reference_dimension'] == 2
    assert reference['quantum_output_dimension_per_record'] == 132
    assert reference['compressed_complete_cq_distance_diagnostic'] < 2e-14
    assert reference['checkpoint_complete_cq_distance_diagnostic'] == 0
    assert reference['checkpoint_json_independently_replayed']
    assert not result['scope']['native_pulse_words_executed']
    assert not result['scope']['physical_event_duration_reduced']
    assert not result['scope']['empirical_data_received_or_fitted']
    json.dumps(result, allow_nan=False)


def test_complete_terminal_minus_space_is_compared_not_just_pointer_probability():
    events = 2
    history = runtime.execute(report._program(events, 'iterated'))
    rho = report._analytic_terminal(events)
    # Corrupt the opposite-orientation spectator population without changing
    # the plus effect. The complete terminal comparison must reject it.
    a, u, v = report.completion.witness_indices()
    other = next(i for i in range(66) if i not in (a, u, v))
    wrong = rho.copy()
    wrong[v, v] -= .01
    wrong[other, other] += .01
    assert np.trace(wrong) == 1
    assert np.trace(report._preparation()@wrong) == np.trace(report._preparation()@rho)
    with pytest.raises(ArithmeticError, match='terminal density'):
        report._terminal_error(history, wrong)


def test_collapsing_every_finite_step_to_same_completion_is_rejected(monkeypatch):
    original = report.completion.compressed_kraus
    monkeypatch.setattr(report.completion, 'compressed_kraus',
                        lambda events: report.completion.binary_kraus('E'))
    with pytest.raises(ArithmeticError, match='probability'):
        report.run()
    monkeypatch.setattr(report.completion, 'compressed_kraus', original)


@pytest.mark.parametrize('bad', [float('nan'), float('inf'), -1., 1e-5])
def test_nonfinite_or_failed_diagnostics_cannot_be_exported(bad):
    with pytest.raises(ArithmeticError):
        report._bounded(bad, 'hostile diagnostic')


def test_analytic_terminal_reference_is_a_normalized_positive_mixed_state():
    # It is not the prepared projector: repair has mixed the endpoint orientation.
    rho = report._analytic_terminal(4)
    assert np.trace(rho) == 1
    assert np.linalg.eigvalsh(rho)[0] >= 0
    assert np.trace(rho@rho).real < 1
    a, u, v = report.completion.witness_indices()
    assert rho[a, a] == .5 and rho[u, u] == rho[v, v] == .25
    assert rho[a, u] == rho[u, a] == 1/32
