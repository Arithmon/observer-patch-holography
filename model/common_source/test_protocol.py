"""Executed composition: legacy agreement, retained states, and prefix laws."""
from dataclasses import replace
from fractions import Fraction as F
from unittest.mock import patch

import numpy as np
import pytest

from . import experiment, scalar_coupling, smooth_protocol
from . import protocol as p


def vector(row):
    value = np.asarray(row['unnormalized_terminal_state'])
    return (value[:, 0]+1j*value[:, 1])[:, None]


@pytest.mark.parametrize('dimension,current,feedback', [
    (2, -1, True), (2, 0, False), (2, 1, True), (4, 0, True)])
def test_program_replays_actual_existing_smooth_instrument(dimension, current, feedback):
    program = p.smooth_program(dimension, 12, native_current=current, feedback=feedback)
    result = p.execute(program)
    retained = smooth_protocol.history(dimension, 12, native_current=current, feedback=feedback)
    assert result.records == ((0, 0), (0, 1), (1, 0), (1, 1))
    for actual, old in zip(result.branches, retained['branches']):
        assert actual.record == tuple(old['record'])
        np.testing.assert_allclose(actual.state, vector(old), atol=2e-12, rtol=0)
        assert abs(actual.probability-old['probability']) < 2e-12
    assert dict(result.branches[0].resources)['work_preparations'] == 1


@pytest.mark.parametrize('dimension,current', [(2, 0), (4, 1)])
def test_explicit_sharp_compatibility_keeps_its_original_record_meaning(dimension, current):
    result = p.execute(p.sharp_program(dimension, 12, native_current=current))
    retained = experiment.history(dimension, 12, current=current)
    for actual, old in zip(result.branches, retained['branches']):
        np.testing.assert_allclose(actual.state, vector(old), atol=2e-12, rtol=0)
    meanings = result.program.specification()['bit_meanings']
    assert meanings['smooth_read'].startswith('0: (I+sin Q0)')
    assert meanings['sharp_read'].startswith('1: Q0>0')


@pytest.fixture(scope='module')
def continuation():
    program = p.continuation_program(2, 24, native_current=1)
    return program, p.execute(program)


def test_continuation_preserves_every_prefix_effect_on_the_whole_work_algebra(continuation):
    program, final = continuation
    for cursor in (2, 5, 7):
        prefix = p.execute(program, stop_after=cursor)
        length = len(prefix.branches[0].record)
        effects = p.prefix_effects(final, length)
        for branch in prefix.branches:
            np.testing.assert_allclose(effects[branch.record], branch.kraus.conj().T@branch.kraus,
                                       atol=2e-12, rtol=0)
    assert len(final.branches) == 16
    assert final.completeness_error < 2e-12
    assert abs(final.probability_mass-1) < 2e-12


def test_complete_source_branch_maps_match_independent_product_evolution(continuation):
    program, source = continuation
    product = p.execute(program, backend='product')
    comparison = p.compare_histories(source, product)
    assert comparison['maximum_kraus_operator_error'] < 3e-12
    assert comparison['complete_flagged_vector_error'] < 3e-12
    assert comparison['classical_total_variation'] < 3e-12
    assert len(source.stages) == len(program.operations)


def test_supplied_reference_columns_are_retained_without_normalizing_branches(continuation):
    program, _ = continuation
    rng = np.random.default_rng(721)
    initial = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    initial /= np.linalg.norm(initial)
    result = p.execute(program, initial=initial)
    product = p.execute(program, initial=initial, backend='product')
    for branch in result.branches:
        assert branch.state.shape == (16, 3)
        np.testing.assert_allclose(branch.state, branch.kraus@initial, atol=2e-13, rtol=0)
    # Tracing out work after the complete local instrument fixes the reference.
    after = sum(branch.state.conj().T@branch.state for branch in result.branches)
    np.testing.assert_allclose(after, initial.conj().T@initial, atol=2e-12, rtol=0)
    assert p.compare_histories(result, product)['complete_flagged_vector_error'] < 3e-12
    assert result.preparation_mode == 'supplied_input'
    assert dict(result.preparation_resources)['work_preparations'] == 0


def test_pause_resume_at_every_operation_has_identical_states_maps_and_costs(continuation):
    program, whole = continuation
    for cursor in range(len(program.operations)+1):
        prefix = p.execute(program, stop_after=cursor)
        resumed = p.execute(program, checkpoint=prefix)
        assert resumed.seal == whole.seal
        for actual, expected in zip(resumed.branches, whole.branches):
            np.testing.assert_array_equal(actual.state, expected.state)
            np.testing.assert_array_equal(actual.kraus, expected.kraus)
            assert actual.resources == expected.resources


def test_serialized_roundtrip_reexecutes_prefix_then_continues(continuation):
    program, whole = continuation
    prefix = p.execute(program, stop_after=7)
    text = p.checkpoint_to_json(prefix)
    restored = p.checkpoint_from_json(program, text)
    assert p.checkpoint_to_json(restored) == text
    assert p.execute(program, checkpoint=restored).seal == whole.seal
    with pytest.raises(ValueError):
        restored.branches[0].state.setflags(write=True)
    with pytest.raises(ValueError):
        restored.branches[0].kraus[0, 0] = 0


def test_resource_accounting_distinguishes_one_realized_path_from_expanded_tree(continuation):
    _, result = continuation
    for branch in result.branches:
        counts = dict(branch.resources)
        assert counts['work_preparations'] == 1
        assert counts['fresh_pointer_blanks'] == 4
        assert counts['committed_classical_bits'] == 4
        assert counts['native_pulses'] > 0
        assert counts['proper_code_service_events'] == 2*counts['native_pulses']
    tree = dict(result.expanded_tree_resources)
    assert tree['work_preparations'] == 1
    assert tree['fresh_pointer_blanks'] == 1+2+4+8
    assert tree['native_pulses'] > dict(result.branches[0].resources)['native_pulses']


def test_split_error_telescope_covers_the_complete_continued_instrument():
    coarse = p.continuation_program(2, 12)
    fine = replace(coarse, steps=24)
    errors = []
    for program in (coarse, fine):
        split = p.execute(program, backend='product')
        continuous = p.execute(program, backend='continuous')
        # Stack full branch maps: this measures the complete flagged dilation,
        # including phases and every input column, rather than selected counts.
        difference = np.vstack([a.kraus-b.kraus
                                for a, b in zip(split.branches, continuous.branches)])
        error = float(np.linalg.norm(difference, 2))
        bound = sum(F(scalar_coupling.time_error_bound(
            program.dimension, program.steps, operation.time, program.coupling
        )['ideal_split_operator_error_upper']) for operation in program.operations
                    if operation.kind == 'evolve')
        assert error <= float(bound)
        errors.append(error)
    assert 3.8 < errors[0]/errors[1] < 4.2


def test_new_feedback_really_uses_second_smooth_record(continuation):
    program, on = continuation
    off_ops = tuple(replace(op, enabled=False) if index == 7 else op
                    for index, op in enumerate(program.operations))
    off = p.execute(replace(program, operations=off_ops))
    # No earlier outcome is changed, but the later gauge clock record responds.
    for key, value in p.prefix_effects(on, 3).items():
        np.testing.assert_allclose(value, p.prefix_effects(off, 3)[key], atol=2e-12, rtol=0)
    on_last = sum(b.probability for b in on.branches if b.record[-1] == 1)
    off_last = sum(b.probability for b in off.branches if b.record[-1] == 1)
    assert abs(on_last-off_last) > 1e-4


def test_execution_budget_is_checked_before_source_preparation(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('preparation must not occur before budget validation')
    monkeypatch.setattr(scalar_coupling, 'preparation', forbidden)
    with pytest.raises(ValueError, match='budget'):
        p.execute(p.continuation_program(), max_branches=8)


def test_both_legacy_receipt_producers_use_the_canonical_execution_engine():
    with patch.object(p, 'execution_prefixes', wraps=p.execution_prefixes) as execute:
        experiment.history(2, 8)
        smooth_protocol.history(2, 8)
    assert execute.call_count == 2
    assert [call.args[0].operations[-1].kind for call in execute.call_args_list] == [
        'sharp_read', 'smooth_read']
