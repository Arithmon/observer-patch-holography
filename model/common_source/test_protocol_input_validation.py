"""Canonical program admission and detachment at public execution boundaries."""
from dataclasses import replace
from fractions import Fraction as F

import numpy as np
import pytest

from . import protocol as p
from . import smooth_history


@pytest.mark.parametrize('backend', ['source', 'product', 'continuous'])
def test_execute_revalidates_native_current_before_using_a_mutated_program(backend):
    program = p.Program()
    program.__dict__['native_current'] = 2
    with pytest.raises(ValueError):
        p.execute(program, backend=backend)


def test_continuous_execution_rejects_mutated_negative_time():
    operation = p.Operation('evolve', F(1, 8))
    program = p.Program(operations=(operation,))
    operation.__dict__['time'] = F(-1, 8)
    with pytest.raises(ValueError):
        p.execute(program, backend='continuous')


def test_prefix_execution_revalidates_future_feedback_even_before_first_read():
    feedback = p.Operation('feedback', record=0)
    program = p.Program(operations=(p.Operation('gauge_read'), feedback))
    feedback.__dict__['record'] = 1
    with pytest.raises(ValueError):
        p.execute(program, backend='continuous', stop_after=0)


@pytest.mark.parametrize('component', ['program', 'operation', 'time', 'coupling'])
def test_history_does_not_share_program_inputs_with_its_caller(component):
    operation = p.Operation('evolve', F(1, 17))
    program = p.Program(steps=3, coupling=F(-2, 7), native_current=-1,
                        operations=(operation, p.Operation('gauge_read')))
    history = p.execute(program, backend='continuous')
    original = history.program.specification()
    expected_program = replace(program, operations=tuple(replace(op) for op in program.operations))
    if component == 'program':
        program.__dict__['native_current'] = 1
    elif component == 'operation':
        operation.__dict__['time'] = F(2, 17)
    elif component == 'time':
        object.__setattr__(operation.time, '_numerator', 2)
    else:
        object.__setattr__(program.coupling, '_numerator', -1)
    assert history.program.specification() == original
    p.validate_checkpoint(expected_program, history, 'continuous')


@pytest.mark.parametrize('field,value', [
    ('dimension', True), ('dimension', 3), ('steps', False), ('steps', 1.),
    ('coupling', .25), ('coupling', True), ('native_current', np.int64(0)),
])
def test_program_contracts_are_checked_before_any_preparation(monkeypatch, field, value):
    program = p.Program()
    program.__dict__[field] = value

    def forbidden(*args, **kwargs):
        raise AssertionError('invalid program reached source preparation')

    monkeypatch.setattr(p.scalar, 'preparation', forbidden)
    with pytest.raises(ValueError):
        p.execute(program)


@pytest.mark.parametrize('field,value', [
    ('kind', 'unknown'), ('time', .5), ('time', True), ('record', 0),
    ('enabled', False), ('enabled', 1),
])
def test_nested_contracts_are_checked_before_any_preparation(monkeypatch, field, value):
    operation = p.Operation('gauge_read')
    program = p.Program(operations=(operation,))
    operation.__dict__[field] = value

    def forbidden(*args, **kwargs):
        raise AssertionError('invalid operation reached source preparation')

    monkeypatch.setattr(p.scalar, 'preparation', forbidden)
    with pytest.raises(ValueError):
        p.execute(program, stop_after=0)


def test_program_container_and_size_guards_precede_entry_traversal(monkeypatch):
    class UntrustedTuple(tuple):
        def __iter__(self):
            raise AssertionError('noncanonical container was traversed')

    def forbidden(*args, **kwargs):
        raise AssertionError('operation was reconstructed before the container gate')

    monkeypatch.setattr(p, 'replace', forbidden)
    for entries in (UntrustedTuple(), [p.Operation('evolve')],
                    (p.Operation('evolve'),)*(p.MAX_OPERATIONS+1)):
        program = p.Program()
        program.__dict__['operations'] = entries
        with pytest.raises(ValueError, match='bounded immutable tuple'):
            p.validate_program(program)


def test_canonical_budget_endpoints_and_nondefault_rationals_remain_admitted():
    for operations in ((p.Operation('evolve'),)*p.MAX_OPERATIONS,
                       (p.Operation('gauge_read'),)*p.MAX_READS):
        original = p.Program(dimension=4, steps=3, coupling=F(-2, 7),
                             native_current=-1, operations=operations)
        detached = p.validate_program(original)
        assert detached == original and detached is not original
        assert all(a is not b for a, b in zip(detached.operations, original.operations))
        assert detached.coupling is not original.coupling
    program = p.Program()
    program.__dict__['operations'] = (p.Operation('gauge_read'),)*(p.MAX_READS+1)
    with pytest.raises(ValueError, match='eight-read'):
        p.validate_program(program)


@pytest.mark.parametrize('boundary', ['validate', 'resume', 'serialize'])
def test_checkpoint_boundaries_revalidate_the_saved_program_even_with_a_consistent_seal(boundary):
    program = p.Program()
    history = p.execute(program, backend='continuous')
    history.program.__dict__['native_current'] = False  # Equal to 0, but inadmissible.
    history = replace(history, seal=p._seal(history))
    with pytest.raises(ValueError):
        if boundary == 'validate':
            p.validate_checkpoint(program, history, 'continuous')
        elif boundary == 'resume':
            p.execute(program, backend='continuous', checkpoint=history)
        else:
            p.checkpoint_to_json(history)


@pytest.mark.parametrize('boundary', ['validate', 'resume', 'deserialize'])
def test_checkpoint_boundaries_return_detached_programs(boundary):
    program = p.Program(native_current=1, operations=(p.Operation('gauge_read'),))
    prefix = p.execute(program, backend='continuous', stop_after=0)
    if boundary == 'validate':
        result = p.validate_checkpoint(program, prefix, 'continuous')
    elif boundary == 'resume':
        result = p.execute(program, backend='continuous', stop_after=0, checkpoint=prefix)
    else:
        result = p.checkpoint_from_json(program, p.checkpoint_to_json(prefix), backend='continuous')
    expected = result.program.specification()
    program.__dict__['native_current'] = -1
    prefix.program.operations[0].__dict__['kind'] = 'smooth_read'
    assert result.program.specification() == expected
    p.validate_checkpoint(result.program, result, 'continuous')


def test_checkpoint_import_rejects_a_mutated_requested_program_before_json_parsing(monkeypatch):
    program = p.Program()
    program.__dict__['steps'] = 0

    def forbidden(*args, **kwargs):
        raise AssertionError('invalid program reached the JSON parser')

    monkeypatch.setattr(p.json, 'loads', forbidden)
    with pytest.raises(ValueError):
        p.checkpoint_from_json(program, '{}')


@pytest.mark.parametrize('backend', ['source', 'product', 'continuous'])
@pytest.mark.parametrize('dimension', [2, 4])
@pytest.mark.parametrize('reference_columns', [1, 3])
def test_nondefault_program_preserves_original_work_reference_preparation_semantics(
        backend, dimension, reference_columns):
    program = p.Program(dimension=dimension, steps=3, coupling=F(-2, 7),
                        native_current=-1, operations=(p.Operation('gauge_read'),
                        p.Operation('feedback', record=0, enabled=False),
                        p.Operation('evolve', F(1, 17))))
    initial = np.arange(1, program.work_dimension*reference_columns+1, dtype=float)
    initial = initial + 1j*initial[::-1]
    initial /= np.linalg.norm(initial)
    initial = initial.reshape(program.work_dimension, reference_columns)
    supplied = initial[:, 0] if reference_columns == 1 else initial
    history = p.execute(program, backend=backend, initial=supplied)
    assert history.program == program and history.program is not program
    assert history.preparation_mode == 'supplied_input'
    assert dict(history.preparation_resources)['work_preparations'] == 0
    assert history.initial_state.shape == initial.shape
    np.testing.assert_array_equal(history.initial_state, initial)
    for branch in history.branches:
        np.testing.assert_allclose(branch.state, branch.kraus@initial, atol=2e-12, rtol=0)
    # Reference columns still denote a pure joint amplitude, not a density matrix.
    with pytest.raises(ValueError, match='unit norm'):
        p.execute(program, backend=backend, initial=np.eye(program.work_dimension)/program.work_dimension)


def test_smooth_certificate_and_execution_use_the_same_program_validator(monkeypatch):
    original = p.validate_program
    calls = []

    def checked(program):
        calls.append(program)
        return original(program)

    monkeypatch.setattr(p, 'validate_program', checked)
    program = p.Program(operations=(p.Operation('gauge_read'),))
    resolution = smooth_history.Resolution(cutoff=2, preparation=0, degree=1,
                                            steps_per_unit=2, bits=32)
    enclosure = smooth_history.enclose(program, dimension=2, resolution=resolution)
    history = p.execute(program, backend='continuous')
    assert calls == [program, program]
    assert enclosure.program == history.program == program
    assert enclosure.program is not history.program and history.program is not program
