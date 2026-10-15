"""Independent finite composition witnesses and receipt/verifier attacks."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from ..common_source import protocol, scalar_coupling
from ..primitive_observer import runtime
from . import clock, construction as c, experiment as e


@pytest.fixture(scope='module')
def small():
    work = c.execute_work(1)
    native = clock.execute_native('C', 1)
    bounds = e.work_certificates(work)
    return work, native, bounds


def _reseal(history, branches):
    changed = replace(history, branches=tuple(branches), seal='')
    return replace(changed, seal=runtime._seal(changed))


def _nuclear(matrix):
    return float(np.linalg.svd(matrix, compute_uv=False).sum())


def test_complete_work_composition_against_independent_continuous_source(small):
    work, native, bounds = small
    initial = np.zeros(16, complex); initial[0] = 1
    waiting = expm(-1j*float(work.delta)*scalar_coupling.ideal_hamiltonian(2, .25))
    expected_gauge = 0.
    native_probabilities = clock.statistics('C', 1)['probabilities']
    for n, current, history in work.successes:
        target = protocol.execute(protocol.continuation_program(2, work.steps, native_current=current),
                                  backend='continuous', initial=waiting@initial)
        distance = 0.
        for branch, ideal in zip(history.branches, target.branches, strict=True):
            assert branch.record == ideal.record
            density = ideal.state@ideal.state.conj().T
            distance += _nuclear(branch.density-density)/2
            if branch.record[3]:
                expected_gauge += float(native_probabilities[n, current])*float(np.trace(density).real)
        assert 0 < distance <= float(bounds['errors'][n, current])+2e-13
    wanted_timeout = waiting@initial
    timeout_error = _nuclear(work.prefixes[1].branches[0].density-
                            np.outer(wanted_timeout, wanted_timeout.conj()))/2
    assert timeout_error <= float(bounds['errors'][1, 0])+2e-13
    proof = e.joint_certificate(work, bounds, 'C', native, cap=1)
    observed = e.exact_observations(work, native, cap=1)
    assert abs(float(observed['final_gauge_one_and_success'])-expected_gauge) <= 2*float(proof['finite_error'])
    assert proof['eventual_error'] == proof['finite_error']+native_probabilities[1, 0]
    assert 0 < proof['finite_error'] <= e.ERROR_GATE


def test_retained_state_reset_negative_control_is_resolved_without_long_horizon(small):
    work, native, bounds = small
    baseline = e.joint_certificate(work, bounds, 'C', native, cap=1)
    reset = e.joint_certificate(work, bounds, 'C', native, cap=1, reset=True)
    observed = e.exact_observations(work, native, cap=1)
    reset_observed = e.exact_observations(work, native, cap=1, reset=bounds['reset'])
    difference = abs(observed['first_slot_gauge_one']-reset_observed['first_slot_gauge_one'])
    assert difference > 2*(baseline['finite_error']+reset['finite_error'])


def test_actual_density_certificate_rejects_reset_outputs_with_valid_source_headers(small):
    work, _, bounds = small
    wrong = []
    for n, current, history in work.successes:
        reset = {branch.record: branch.density for branch in bounds['reset'][current].branches}
        branches = [replace(branch, density=reset[branch.record]) for branch in history.branches]
        wrong.append((n, current, _reseal(history, branches)))
    fabricated = replace(work, successes=tuple(wrong))
    # Source headers, resources, positivity, total mass and recomputed seals
    # all survive. The exact arithmetic comparison must still reject reset.
    c.validate_work(fabricated)
    with pytest.raises((ValueError, ArithmeticError)):
        e.work_certificates(fabricated)


def test_native_postselection_and_wrong_trace_preserving_outputs_cannot_pass(small):
    work, native, bounds = small
    permutation = np.roll(np.eye(66), 1, axis=0)
    wrong = {record: permutation@density@permutation.T for record, density in native.items()}
    assert abs(sum(np.trace(value) for value in wrong.values())-1) < 1e-12
    with pytest.raises(ArithmeticError):
        e.joint_certificate(work, bounds, 'C', wrong, cap=1)
    mass = sum(np.trace(value).real for (_, current), value in native.items() if current)
    selected = {record: value/mass if record[1] else np.zeros_like(value) for record, value in native.items()}
    assert abs(sum(np.trace(value) for value in selected.values())-1) < 1e-12
    with pytest.raises(ArithmeticError):
        e.joint_certificate(work, bounds, 'C', selected, cap=1)


def test_resource_upper_bound_matches_every_actual_source_path(small):
    work, _, _ = small
    result = e.resources(work, 1)
    for history in (*[h for _, _, h in work.successes], work.prefixes[1]):
        for branch in history.branches:
            for name, count in branch.resources:
                assert count <= result['maximum_work_source_resources'][name]
    assert result['maximum_work_source_resources']['work_preparations'] == 1
    assert result['maximum_work_source_resources']['fresh_pointer_blanks'] == 4
    assert result['joint_complete_record_alphabet_size'] == 33
    assert result['encoded_clock_record_bits'] == 2
    assert result['native_channel_source_pulse_compilation_claimed'] is False
    stats = clock.statistics('C', 1)
    step = dict(work.prefixes[0].program.events[0].instruments[0].resources)['native_pulses']
    prep = dict(work.prefixes[0].program.preparation_resources)['native_pulses']
    assert F(result['expected_work_source_pulses_upper']) >= prep+step*stats['expected_controller_slots']


def test_resource_report_rejects_resealed_native_pulse_undercounts(small):
    work, _, _ = small
    wrong = []
    for n, current, history in work.successes:
        branches = [replace(branch, resources=tuple((key, 0 if key == 'native_pulses' else value)
                    for key, value in branch.resources)) for branch in history.branches]
        wrong.append((n, current, _reseal(history, branches)))
    fabricated = replace(work, successes=tuple(wrong))
    with pytest.raises((ValueError, ArithmeticError)):
        e.resources(fabricated, 1)


def test_work_certificate_requires_complete_suffix_outcome_menu(small):
    work, _, _ = small
    first = work.successes[0]
    missing = _reseal(first[2], first[2].branches[1:])
    fabricated = replace(work, successes=((first[0], first[1], missing), *work.successes[1:]))
    with pytest.raises((ValueError, ArithmeticError)):
        e.work_certificates(fabricated)


def test_joint_certificate_cannot_reuse_bounds_for_another_valid_source_history(small):
    _, native, bounds = small
    different = c.execute_work(1, delta=F(1, 2))
    c.validate_work(different)
    with pytest.raises((ValueError, ArithmeticError)):
        e.joint_certificate(different, bounds, 'C', native, cap=1)


def test_observations_reject_aliased_native_record_types_and_wrong_original_shapes(small):
    work, native, _ = small
    boolean_records = {(True, current): array for (_, current), array in native.items()}
    with pytest.raises((ValueError, ArithmeticError)):
        e.exact_observations(work, boolean_records, cap=1)
    # These 1x1 matrices preserve exactly the original traces. An observable
    # or mass-only check cannot notice the missing complete native state.
    collapsed = {record: np.array([[np.trace(array)]]) for record, array in native.items()}
    with pytest.raises((ValueError, ArithmeticError)):
        e.exact_observations(work, collapsed, cap=1)
    masked = {record: np.ma.array(array, mask=False) for record, array in native.items()}
    with pytest.raises((ValueError, ArithmeticError)):
        e.exact_observations(work, masked, cap=1)


def test_reset_observation_override_requires_current_specific_reset_source_history(small):
    work, native, bounds = small
    # Valid complete source histories which *include waiting* cannot stand in
    # for the reset histories to which the separate reset certificate applies.
    wrong_programs = {current: history for _, current, history in work.successes}
    with pytest.raises((ValueError, ArithmeticError)):
        e.exact_observations(work, native, cap=1, reset=wrong_programs)
    aliased = {-1: bounds['reset'][-1], True: bounds['reset'][1]}
    with pytest.raises((ValueError, ArithmeticError)):
        e.exact_observations(work, native, cap=1, reset=aliased)


def test_original_complex_traces_and_observable_widening_do_not_discard_defects():
    native = np.zeros((66, 66), complex); native[0, 0] = 2+.25j
    work = np.zeros((16, 16), complex); work[0, 0] = .5-.125j
    assert e.probability(native, work) == (F(33, 32), F(-1, 8))
    # A general non-Hermitian report has no trace-zero difference guarantee.
    # The real effect expectation can move by the full nuclear norm, 2*error.
    assert tuple(map(F, e.interval(F(1, 2), F(1, 10)))) == (F(3, 10), F(7, 10))
    assert tuple(map(F, e.interval(F(1, 10), F(1, 10)))) == (F(0), F(3, 10))


@pytest.mark.parametrize('bad', (True, False, 0, .25, F(-1)))
def test_intermediate_bound_rounding_rejects_ambiguous_or_negative_values(bad):
    with pytest.raises(ValueError):
        e.up(bad)


def test_receipt_comparison_rejects_truthy_scalars_missing_fields_and_changed_guarantees(small):
    work, native, bounds = small
    proof = e.joint_certificate(work, bounds, 'C', native, cap=1)
    actual = {'all_pass': True, 'error': e.report.outward(proof['finite_error']),
              'resources': e.resources(work, 1), 'observations': [e.interval(F(1, 2), proof['finite_error'])]}
    e.report.compare_receipt(deepcopy(actual), actual)
    changes = []
    changed = deepcopy(actual); changed['all_pass'] = 1; changes.append(changed)
    changed = deepcopy(actual); changed['error'] = '0'; changes.append(changed)
    changed = deepcopy(actual); changed['resources']['native_preparations'] = False; changes.append(changed)
    changed = deepcopy(actual); changed['observations'].clear(); changes.append(changed)
    changed = deepcopy(actual); changed['extra'] = True; changes.append(changed)
    for changed in changes:
        with pytest.raises(ArithmeticError):
            e.report.compare_receipt(changed, actual)
