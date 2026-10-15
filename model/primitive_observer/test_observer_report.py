"""The retained observer report is bound to actual complete continuing states."""
import json
from dataclasses import replace

import numpy as np
import pytest

from . import algebra, common, observer_report as report, runtime
from ..common_source import protocol


@pytest.fixture(scope='module')
def result():
    return report.run()


@pytest.fixture(scope='module')
def work():
    return common.execute(protocol.continuation_program(report.DIMENSION, report.STEPS),
                          report.canonical_environment())


def test_report_retains_both_native_laws_active_passive_controls_and_complete_raw_outcomes(result):
    assert result['all_pass'] and all(result['checks'].values())
    rows = result['native_and_reference_reports']
    assert {(row['law'], row['mode']) for row in rows} == {
        (law, mode) for law in ('C', 'E') for mode in ('protected_qubit', 'passive_qubit')}
    for row in rows:
        reference = row['joint_reference_output']
        assert reference['reference_dimension'] == 3
        assert reference['joint_probability_mass'] == pytest.approx(1., abs=2e-9)
        assert reference['reference_marginal_residual'] < 2e-9
        assert reference['complete_raw_outcomes'] == (770 if row['law'] == 'C' else 2566)
        assert row['native_channel_audit']['first_abort_outcome_count'] > 0
        assert row['native_channel_audit']['final_abort_outcome_count'] > 0
        oracle = row['independent_original_wrapper_comparison']
        assert oracle['complete_reference_output_trace_distance'] < 2e-9
        assert oracle['maximum_work_block_trace_distance'] < 2e-9
        assert oracle['every_native_map_and_raw_record_agrees']
        assert oracle['source_resources_agree']
    json.dumps(result, allow_nan=False)


def test_report_attaches_views_and_forecast_to_same_native_work_program(result, work):
    attached = result['same_native_work_history_access']
    assert attached['work_history_canonical_declaration_digest'] == protocol.continuation_program(
        report.DIMENSION, report.STEPS).digest
    assert attached['visible_record_count'] == 16
    assert len(attached['record_law']) == 16
    actual = {tuple(row['record']): row['probability'] for row in attached['record_law']}
    assert actual == pytest.approx({branch.record: branch.probability for branch in work.branches}, abs=2e-12)
    assert attached['regional_GNS_dimension'] == 64
    assert attached['regional_GNS_isometry_residual'] < 2e-9
    assert attached['regional_GNS_lift_residual'] < 2e-9
    assert attached['full_work_forecast_probability_residual'] < 2e-9
    assert attached['regional_forecast_probability_residual'] < 2e-9
    assert attached['checkpoint_complete_output_distance'] < 2e-9
    assert attached['old_record_marginal_residual'] < 2e-9
    assert attached['conditional_expectation_multiplicativity_failure'] == pytest.approx(1.)


def test_modular_motion_is_nontrivial_but_supplied_time_and_physical_gaps_stay_explicit(result):
    attached = result['same_native_work_history_access']
    assert attached['smallest_normalized_regional_block_eigenvalue_diagnostic'] > 1e-12
    assert attached['maximum_modular_probe_motion'] > 1e-4
    assert attached['modular_state_pairing_residual'] < 2e-9
    assert attached['modular_public_record_residual'] < 2e-9
    assert not attached['modular_parameter_identified_with_physical_time']
    assert not attached['regional_autonomous_interacting_future_inferred']
    assert not result['empirical_data_consumed']
    assert not result['physical_service_error_certified']
    assert not result['spatial_metric_or_physical_clock_recovered']
    assert not result['full_consensus_tower_or_OperationalObserver_instantiated']
    assert result['noise_deletion_complete_output_distance'] > 1e-5


def test_wrong_full_probability_forecast_cannot_pass_same_history_report(monkeypatch, work):
    original = algebra.forecast
    def wrong(*args, **kwargs):
        values = original(*args, **kwargs)
        first = next(iter(values))
        values[first] += .1
        return values
    monkeypatch.setattr(algebra, 'forecast', wrong)
    with pytest.raises(ArithmeticError, match='attachment failed'):
        report.same_history_report(work)


def test_nonfinite_pairing_cannot_hide_behind_numerical_diagnostics(monkeypatch, work):
    monkeypatch.setattr(algebra, 'expectation', lambda *args: complex(float('nan'), 0))
    with pytest.raises(ArithmeticError, match='finite nonnegative'):
        report.same_history_report(work)


def test_report_refuses_incomplete_prefix_and_wrong_work_algebra(work):
    prefix = runtime.execute(work.program, stop_after=1)
    with pytest.raises(ValueError, match='complete canonical'):
        report.same_history_report(prefix)
    simple = runtime.execute(runtime.Program('wrong_work', 2, np.eye(2)/2, ()))
    with pytest.raises(ValueError, match='two-scalar'):
        report.same_history_report(simple)


def test_canonical_declaration_cannot_relabel_a_different_executed_program(work):
    # This altered program preserves its canonical metadata and complete
    # record menu, but adds an actual intervention before the first read.
    intervention = runtime.Event('extra_intervention', instruments=(runtime.Instrument(
        kraus=((np.diag([-1., *([1.]*15)]),),)),))
    altered = replace(work.program, events=(intervention, *work.program.events))
    history = runtime.execute(altered)
    with pytest.raises(ValueError, match='executed canonical program'):
        report.same_history_report(history)


def test_canonical_report_refuses_a_supplied_preparation_even_if_numerically_equal(work):
    history = runtime.execute(work.program, initial=work.program.initial_density)
    with pytest.raises(ValueError, match='declared preparation'):
        report.same_history_report(history)


def test_declared_reference_input_has_nontrivial_complex_entanglement():
    initial = report.reference_input()
    assert initial.shape == (2, 3)
    assert np.linalg.norm(initial) == pytest.approx(1.)
    assert np.linalg.matrix_rank(initial) == 2
    reference = initial.T@initial.conj()
    assert np.max(np.abs(reference.imag)) > .01
