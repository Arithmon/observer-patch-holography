"""Independent full-history controls for the prepared-output certificate.

The numerical oracle starts from the declared Fock operators and original
preparation parameters. It does not reuse producer polynomial actions, norms,
support cuts, pointer maps, or runtime histories.
"""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
from functools import lru_cache
import math
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from . import protocol, smooth_history as s


FAST = s.Resolution(cutoff=8, preparation=6, degree=12,
                    steps_per_unit=8, bits=64)


@lru_cache(maxsize=2)
def _original_operators(d):
    annihilation = np.diag(np.sqrt(np.arange(1, d)), 1)
    q = (annihilation+annihilation.T)/math.sqrt(2)
    p = (annihilation-annihilation.T)/(1j*math.sqrt(2))
    boundary = np.zeros((d, d))
    boundary[-1, -1] = d/2
    onsite = ((p@p).real+boundary+(1+9/32)*(q@q+boundary))/2
    x = np.array([[0., 1.], [1., 0.]])
    gauge_x = np.kron(np.eye(2), x)
    gauge = gauge_x+np.kron(x, np.eye(2))+np.diag([0., 1., 1., 0.])
    hamiltonian = np.kron(np.kron(np.eye(d), onsite)
                         +np.kron(onsite, np.eye(d))
                         -(9/32)*np.kron(q, q), np.eye(4))
    hamiltonian += np.kron(np.eye(d*d), gauge)
    hamiltonian += (math.sqrt(2)/4)*np.kron(np.kron(np.eye(d), q),
                                         np.diag([1., -1., 1., -1.]))
    coordinate = np.kron(np.kron(np.eye(d), q), np.eye(4))
    minus, plus = expm(-.5j*coordinate), expm(.5j*coordinate)
    pointer = ((minus-1j*plus)/2, (minus+1j*plus)/2)
    vacuum = np.zeros((d, d))
    for pair in range(d//2):
        occupation = 2*pair
        amplitude = (-1)**pair*math.sqrt(math.comb(occupation, pair))/(2**pair*9**pair)
        for n0 in range(occupation+1):
            n1 = occupation-n0
            vacuum[n1, n0] = amplitude*(-1)**n1*math.sqrt(math.comb(occupation, n0)/2**occupation)
    vacuum /= np.linalg.norm(vacuum)
    prepared = np.kron(vacuum.ravel(), [1., 0., 0., 0.]).astype(complex)
    return hamiltonian, pointer, np.kron(np.eye(d*d), gauge_x), prepared


@lru_cache(maxsize=6)
def _original_history(d, current):
    h, pointer, x, prepared = _original_operators(d)
    kick = np.eye(len(prepared)) if current == 0 else .6*np.eye(len(prepared))-.8j*current*x
    feedback = .6*np.eye(len(prepared))-.8j*x
    states = {(): kick@prepared}
    # Literal canonical schedule: it must not follow a damaged producer stage
    # dispatch, feedback address, or record decoder.
    for kind, value in (('evolve', 1), ('gauge', None), ('feedback', 0),
                        ('evolve', .5), ('smooth', None), ('evolve', .5),
                        ('smooth', None), ('feedback', 2), ('evolve', .5),
                        ('gauge', None)):
        if kind == 'evolve':
            unitary = expm(-1j*value*h)
            states = {key: unitary@state for key, state in states.items()}
        elif kind == 'feedback':
            states = {key: feedback@state if key[value] else state
                      for key, state in states.items()}
        else:
            following = {}
            for key, state in states.items():
                for bit in (0, 1):
                    output = pointer[bit]@state if kind == 'smooth' else state.copy()
                    if kind == 'gauge':
                        output[1-bit::2] = 0
                    following[(*key, bit)] = output
            states = following
    return states


def _polynomial_vector(enclosure, polynomial):
    d = enclosure.dimension
    vector = np.zeros((d, d, 4), complex)
    divisor = polynomial.denominator*math.sqrt(float(enclosure.normalization))
    for (n0, n1, gauge), (real, imag) in polynomial.coefficients.items():
        vector[n1, n0, gauge] = complex(real, imag)*math.sqrt(
            math.factorial(n0)*math.factorial(n1))/divisor
    return vector.ravel()


def _complete_distance(actual, expected):
    """Half trace norm, stably evaluated on each rank-at-most-two block."""
    assert actual.keys() == expected.keys()
    total = 0.
    for key, vector in actual.items():
        reference = expected[key]
        difference = vector-reference
        norm_difference = float(np.vdot(difference, difference).real)
        overlap_difference = np.vdot(reference, difference)
        mass = float(np.vdot(reference, reference).real)
        discriminant = ((2*overlap_difference.real+norm_difference)**2
                        +4*(mass*norm_difference-abs(overlap_difference)**2))
        total += .5*math.sqrt(max(0., discriminant))
    return total


def _assert_original_output(enclosure):
    actual = {key: _polynomial_vector(enclosure, vector)
              for key, vector in enclosure.branches.items()}
    expected = _original_history(enclosure.dimension, enclosure.program.native_current)
    assert len(actual) == 16
    assert abs(sum(float(np.vdot(v, v).real) for v in expected.values())-1) < 2e-13
    assert _complete_distance(actual, expected) <= float(enclosure.output_error)+2e-13
    return actual, expected


@pytest.mark.parametrize('dimension', [2, 4])
@pytest.mark.parametrize('current', [-1, 0, 1])
def test_every_complete_output_matches_original_input_oracle(dimension, current):
    enclosure = s.enclose(protocol.continuation_program(dimension, native_current=current),
                          dimension=dimension, resolution=FAST)
    actual, expected = _assert_original_output(enclosure)
    for bit in (0, 1):
        measured = sum(float(np.vdot(v, v).real) for key, v in actual.items() if key[3] == bit)
        target = sum(float(np.vdot(v, v).real) for key, v in expected.items() if key[3] == bit)
        assert abs(measured-target) < 2*float(enclosure.output_error)+2e-13


@pytest.mark.parametrize('mutation', ['record_swap', 'feedback_disabled', 'erased_excited_support'])
def test_full_output_oracle_rejects_semantic_mutations(monkeypatch, mutation):
    if mutation == 'record_swap':
        original = s._pointer

        def changed(*args, **kwargs):
            outputs, errors, support = original(*args, **kwargs)
            return outputs[::-1], errors, support

        monkeypatch.setattr(s, '_pointer', changed)
    elif mutation == 'feedback_disabled':
        # Current zero leaves the first kick intact; both real recorded
        # feedback operations are removed without changing the outcome menu.
        monkeypatch.setattr(s, '_kick', lambda vector, current: vector)
    else:
        original = s._evolve

        def changed(vector, time, normalization, resolution, dimension):
            # Erasing excited support is an undeclared projection, even if
            # later evolution and the declared error sum remain unchanged.
            if time == F(1, 2):
                vector = deepcopy(vector)
                vector.coefficients = {
                    (0, 0, 0): vector.coefficients.get((0, 0, 0), (0, 0))}
            return original(vector, time, normalization, resolution, dimension)

        monkeypatch.setattr(s, '_evolve', changed)
    enclosure = s.enclose(protocol.continuation_program(2), dimension=2, resolution=FAST)
    with pytest.raises(AssertionError):
        _assert_original_output(enclosure)


def test_record_probabilities_do_not_certify_quantum_output():
    expected = _original_history(2, 0)
    # The final gauge read makes this phase constant on a branch; use instead
    # a relative scalar phase, leaving every record probability exactly fixed.
    phase = np.repeat([1., -1., 1., -1.], 4)
    wrong = {key: phase*vector for key, vector in expected.items()}
    assert all(abs(np.vdot(wrong[k], wrong[k])-np.vdot(v, v)) < 1e-15
               for k, v in expected.items())
    assert _complete_distance(wrong, expected) > .1
    harmless = {key: (1j if key[0] else -1)*vector for key, vector in expected.items()}
    assert _complete_distance(harmless, expected) < 2e-8


@pytest.fixture
def finite_enclosure():
    return s.enclose(protocol.continuation_program(2), dimension=2, resolution=FAST)


@pytest.mark.parametrize('field', ['normalization', 'vector_error', 'preparation_error',
                                  'maximum_support', 'dimension', 'program', 'resolution',
                                  'branches', 'stages'])
def test_mutable_certificate_metadata_cannot_redefine_verified_target(finite_enclosure, field):
    changed = deepcopy(finite_enclosure)
    if field == 'normalization':
        value = changed.normalization*2
    elif field == 'vector_error':
        value = F(0)
    elif field == 'preparation_error':
        value = F(1, 10)
    elif field == 'maximum_support':
        value = changed.maximum_support+1
    elif field == 'dimension':
        value = None
    elif field == 'program':
        value = replace(changed.program, native_current=1)
    elif field == 'resolution':
        value = replace(changed.resolution, bits=65)
    elif field == 'branches':
        value = dict(changed.branches)
        value.pop(next(iter(value)))
    else:
        value = deepcopy(changed.stages)
        value[0]['kind'] = 'feedback'
    changed = replace(changed, **{field: value})
    for consumer in (s.summary, lambda x: s.distance_interval(x, finite_enclosure)):
        with pytest.raises((ValueError, ArithmeticError)):
            consumer(changed)


def test_matching_omissions_on_both_sides_are_not_a_complete_channel(finite_enclosure):
    changed = deepcopy(finite_enclosure)
    changed.branches.pop((0, 0, 0, 0))
    with pytest.raises((ValueError, ArithmeticError)):
        s.distance_interval(changed, deepcopy(changed))


def test_continuous_target_accepts_distinct_source_step_counts_and_cutoffs():
    a = s.enclose(protocol.continuation_program(2, steps=64), dimension=2, resolution=FAST)
    b = s.enclose(protocol.continuation_program(4, steps=512), dimension=4, resolution=FAST)
    low, high = s.distance_interval(a, b)
    expected = _original_history(2, 0)
    larger = _original_history(4, 0)
    embedded = {}
    for key, vector in expected.items():
        full = np.zeros((4, 4, 4), complex)
        full[:2, :2] = vector.reshape(2, 2, 4)
        embedded[key] = full.ravel()
    diagnostic = _complete_distance(embedded, larger)
    assert float(low)-2e-13 <= diagnostic <= float(high)+2e-13


@pytest.mark.parametrize('operations', [
    (protocol.Operation('evolve', F(4001, 1000)),),
    (protocol.Operation('smooth_read'),)*3,
    (protocol.Operation('gauge_read'),)*5,
])
def test_certificate_work_budget_rejects_before_polynomial_execution(monkeypatch, operations):
    def forbidden(*args, **kwargs):
        raise AssertionError('rejected workload reached polynomial preparation')
    monkeypatch.setattr(s.ih, '_initial', forbidden)
    with pytest.raises(ValueError):
        s.enclose(protocol.Program(operations=operations), dimension=2, resolution=FAST)


@pytest.mark.parametrize('bad', [True, 2., '2', None, float('nan')])
def test_resolution_original_types_are_checked(bad):
    with pytest.raises(ValueError):
        s.Resolution(cutoff=bad)


@pytest.mark.parametrize('broken_backend', ['source', 'product'])
def test_execution_join_rejects_nonfinite_branch_states(finite_enclosure, monkeypatch, broken_backend):
    from . import composition, continuation_bounds
    program = finite_enclosure.program
    histories = {backend: composition._executed(program, backend)
                 for backend in ('source', 'product')}
    temporal = F(continuation_bounds.program_bound(program)['ideal_temporal_instrument_error_upper'])
    healthy = s._execution_check(finite_enclosure, temporal)
    assert healthy['source_product_operator_and_branch_checks_pass']
    broken = histories[broken_backend]
    state = broken.branches[0].state.copy()
    state[0, 0] = complex(float('nan'), 0)
    histories[broken_backend] = replace(
        broken, branches=(replace(broken.branches[0], state=state), *broken.branches[1:]))
    monkeypatch.setattr(composition, '_executed',
                        lambda supplied, backend='source': histories[backend])
    with pytest.raises((ValueError, ArithmeticError)):
        s._execution_check(finite_enclosure, temporal)


def test_temporal_slack_cannot_hide_a_changed_preparation(finite_enclosure, monkeypatch):
    from . import composition, continuation_bounds
    program = finite_enclosure.program
    histories = {backend: composition._executed(program, backend)
                 for backend in ('source', 'product')}
    temporal = F(continuation_bounds.program_bound(program)['ideal_temporal_instrument_error_upper'])
    assert s._execution_check(finite_enclosure, temporal)['product_vs_exact_continuous_polynomial_check_pass']
    initial = histories['product'].initial_state.copy()
    initial[1, 0] = 1e-6
    initial /= np.linalg.norm(initial)
    defect = float(np.linalg.norm(initial-histories['product'].initial_state))
    assert 1e-8 < defect < float(temporal)
    for backend, history in histories.items():
        changed = tuple(replace(branch, state=branch.kraus@initial)
                        for branch in history.branches)
        assert all(a.kraus is b.kraus for a, b in zip(changed, history.branches))
        changed_history = replace(history, initial_state=initial, branches=changed)
        changed_history = replace(changed_history, seal=protocol._seal(changed_history))
        # Simulate a producer consistently using the wrong normalized input.
        # This is a valid checkpoint, so only the independent preparation
        # comparison can reject it; an obsolete integrity seal must not mask
        # that scientific assertion.
        protocol.validate_checkpoint(program, changed_history, backend)
        histories[backend] = changed_history
    monkeypatch.setattr(composition, '_executed',
                        lambda supplied, backend='source': histories[backend])
    with pytest.raises((ValueError, ArithmeticError)):
        s._execution_check(finite_enclosure, temporal)


@pytest.mark.parametrize('backend', ['source', 'product'])
def test_preparation_override_cannot_inherit_declared_preparation_service(finite_enclosure, monkeypatch, backend):
    from . import composition, continuation_bounds
    program = finite_enclosure.program
    histories = {kind: composition._executed(program, kind) for kind in ('source', 'product')}
    temporal = F(continuation_bounds.program_bound(program)['ideal_temporal_instrument_error_upper'])
    assert s._execution_check(finite_enclosure, temporal)['source_product_operator_and_branch_checks_pass']
    histories[backend] = protocol.execute(program, backend=backend,
                                         initial=histories[backend].initial_state)
    assert histories[backend].preparation_mode == 'supplied_input'
    assert dict(histories[backend].preparation_resources)['work_preparations'] == 0
    monkeypatch.setattr(composition, '_executed',
                        lambda supplied, backend='source': histories[backend])
    with pytest.raises((ValueError, ArithmeticError)):
        s._execution_check(finite_enclosure, temporal)


def _receipt(enclosure):
    return {'schema': 'oph.common_source.smooth_history_enclosure.v1',
            'all_pass': True, 'certificate': s.summary(enclosure),
            'complete_native_composition': [{'dimension': 2,
                'maximum_current_ideal_product_output_error_upper': '3/10',
                'physical_output_error_upper': None}],
            'scope': {'retained_records_per_work_history': 16,
                      'complete_quantum_outputs_compared': True,
                      'arbitrary_initial_native_work_correlations_included': False}}


def test_cli_replay_checks_every_leaf_and_complete_record_structure(finite_enclosure, tmp_path, monkeypatch):
    actual = _receipt(finite_enclosure)
    path = tmp_path/'receipt.json'
    monkeypatch.setattr(s, 'run', lambda: actual)
    monkeypatch.setattr(sys, 'argv', ['smooth_history', '--check', str(path)])
    path.write_text(s.ih.canonical(actual), encoding='utf-8')
    s.main()

    def leaves(value, address=()):
        if isinstance(value, dict):
            for key, child in value.items():
                yield from leaves(child, (*address, key))
        elif isinstance(value, list):
            for index, child in enumerate(value):
                yield from leaves(child, (*address, index))
        else:
            yield address, value

    checked = 0
    for address, value in leaves(actual):
        wrong = deepcopy(actual)
        parent = wrong
        for component in address[:-1]:
            parent = parent[component]
        parent[address[-1]] = (not value if type(value) is bool else value+1
                              if type(value) is int else 'changed')
        path.write_text(s.ih.canonical(wrong), encoding='utf-8')
        with pytest.raises(ValueError):
            s.main()
        checked += 1
    assert checked > 100
    for mutation in ('missing_branch', 'duplicate_branch', 'extra_key'):
        wrong = deepcopy(actual)
        if mutation == 'missing_branch':
            wrong['certificate']['branches'].pop()
        elif mutation == 'duplicate_branch':
            wrong['certificate']['branches'].append(wrong['certificate']['branches'][0])
        else:
            wrong['unexpected'] = True
        path.write_text(s.ih.canonical(wrong), encoding='utf-8')
        with pytest.raises(ValueError):
            s.main()


@pytest.mark.parametrize('payload', [
    '{"all_pass": false, "all_pass": true}',
    '{"certificate":{"mass":"0","mass":"1"}}',
    '{"all_pass": NaN}', '{"all_pass": Infinity}',
    '{"all_pass": -Infinity}', '{"all_pass": 1.0}',
])
def test_cli_replay_rejects_ambiguous_or_inexact_json(payload, finite_enclosure, tmp_path, monkeypatch):
    path = tmp_path/'bad.json'
    path.write_text(payload, encoding='utf-8')
    monkeypatch.setattr(s, 'run', lambda: _receipt(finite_enclosure))
    monkeypatch.setattr(sys, 'argv', ['smooth_history', '--check', str(path)])
    with pytest.raises(ValueError):
        s.main()


@pytest.mark.parametrize('component', ['resolution', 'nested_stages'])
def test_exported_summary_is_detached_from_its_enclosure(finite_enclosure, component):
    exported = s.summary(finite_enclosure)
    original = deepcopy(exported)
    resolution_snapshot = deepcopy(finite_enclosure.resolution.__dict__)
    try:
        if component == 'resolution':
            exported['resolution']['bits'] = 0
        else:
            exported['stages'][0]['component_l2_upper']['Taylor'] = 'edited for display'
        assert s.summary(finite_enclosure) == original
    finally:
        finite_enclosure.resolution.__dict__.clear()
        finite_enclosure.resolution.__dict__.update(resolution_snapshot)


def test_editing_summary_cannot_poison_later_default_calls():
    # Preserve the pre-fix shared default even when this regression fails,
    # so its intended failure cannot contaminate unrelated tests.
    previous_default = s.enclose.__kwdefaults__['resolution']
    saved = deepcopy(previous_default.__dict__) if previous_default is not None else None
    program = protocol.Program(operations=(protocol.Operation('evolve', F(0)),))
    try:
        first = s.enclose(program, dimension=2)
        exported = s.summary(first)
        exported['resolution']['bits'] = -1
        later = s.enclose(program, dimension=2)
        assert later.resolution.bits == 112
        assert later.vector_error == 0
    finally:
        if saved is not None:
            previous_default.__dict__.clear()
            previous_default.__dict__.update(saved)


@pytest.mark.parametrize(('field', 'value'), [
    ('cutoff', 0), ('preparation', 3), ('degree', 0),
    ('steps_per_unit', 0), ('bits', -1),
])
def test_explicit_resolution_is_revalidated_at_use(field, value):
    resolution = s.Resolution()
    resolution.__dict__[field] = value
    with pytest.raises(ValueError):
        s.enclose(protocol.Program(), dimension=2, resolution=resolution)


def test_caller_owned_resolution_cannot_change_completed_enclosure():
    resolution = s.Resolution()
    enclosure = s.enclose(protocol.Program(), dimension=2, resolution=resolution)
    original = s.summary(enclosure)
    resolution.__dict__['bits'] = 0
    assert s.summary(enclosure) == original


@pytest.mark.parametrize('dimension', [None, 4])
@pytest.mark.parametrize('operations', [(), (protocol.Operation('gauge_read'),),
    (protocol.Operation('gauge_read'), protocol.Operation('feedback', record=0))])
def test_reported_support_includes_the_initial_polynomial(dimension, operations):
    program = protocol.Program(dimension=4, native_current=1, operations=operations)
    enclosure = s.enclose(program, dimension=dimension)
    expected = 12 if dimension is None else 2
    actual = max(n0+n1 for state in enclosure.branches.values()
                 for n0, n1, gauge in state.coefficients)
    assert actual == expected
    assert s.summary(enclosure)['maximum_expanded_total_occupation'] == expected


@pytest.mark.parametrize('current', [-1, 0, 1])
def test_reported_support_includes_finite_pointer_intermediate_power(current):
    # B_2 sends |0> to |1>, while the next power B_2**2 returns to |0>.
    # A remainder-only support check therefore misses an actually executed
    # polynomial coefficient and understates the declared work support.
    program = protocol.Program(native_current=current,
                               operations=(protocol.Operation('smooth_read'),))
    resolution = s.Resolution(cutoff=2, preparation=0, degree=1,
                              steps_per_unit=2, bits=32)
    enclosure = s.enclose(program, dimension=2, resolution=resolution)
    occupied = {n0+n1 for state in enclosure.branches.values()
                for n0, n1, gauge in state.coefficients}
    assert occupied == {0, 1}
    assert s.summary(enclosure)['maximum_expanded_total_occupation'] == 1


@pytest.mark.parametrize('expected_backend', ['source', 'product'])
def test_execution_join_requires_each_actual_backend_role(finite_enclosure, monkeypatch, expected_backend):
    from . import composition, continuation_bounds
    program = finite_enclosure.program
    histories = {backend: composition._executed(program, backend)
                 for backend in ('source', 'product')}
    temporal = F(continuation_bounds.program_bound(program)['ideal_temporal_instrument_error_upper'])
    assert s._execution_check(finite_enclosure, temporal)['source_product_operator_and_branch_checks_pass']
    substituted_backend = 'product' if expected_backend == 'source' else 'source'
    replacement = histories[substituted_backend]
    # This is an intact, healthy executed history, including its original
    # valid seal. Only the source/independent execution role is incorrect.
    protocol.validate_checkpoint(program, replacement, substituted_backend)
    histories[expected_backend] = replacement
    monkeypatch.setattr(composition, '_executed',
                        lambda supplied, backend='source': histories[backend])
    with pytest.raises((ValueError, ArithmeticError)):
        s._execution_check(finite_enclosure, temporal)


@pytest.mark.parametrize(('field', 'value'), [
    ('native_current', 2), ('native_current', True),
    ('dimension', 3), ('steps', 0),
])
def test_program_fields_are_revalidated_before_preparation(monkeypatch, field, value):
    program = protocol.Program()
    program.__dict__[field] = value

    def forbidden(*args, **kwargs):
        raise AssertionError('invalid program reached polynomial preparation')

    monkeypatch.setattr(s.ih, '_initial', forbidden)
    with pytest.raises(ValueError):
        s.enclose(program, resolution=FAST)


@pytest.mark.parametrize('operations', [
    [], [protocol.Operation('gauge_read')], None, 'gauge_read',
    ({'kind': 'gauge_read'},), (object(),),
    (protocol.Operation('evolve'),)*(protocol.MAX_OPERATIONS+1),
])
def test_program_requires_its_original_immutable_typed_sequence(monkeypatch, operations):
    program = protocol.Program()
    program.__dict__['operations'] = operations

    def forbidden(*args, **kwargs):
        raise AssertionError('invalid operation sequence reached polynomial preparation')

    monkeypatch.setattr(s.ih, '_initial', forbidden)
    with pytest.raises(ValueError):
        s.enclose(program, dimension=2, resolution=FAST)


@pytest.mark.parametrize(('kind', 'field', 'value'), [
    ('evolve', 'time', F(-1)), ('evolve', 'time', .5),
    ('evolve', 'time', True), ('evolve', 'kind', 'unknown'),
    ('gauge_read', 'time', F(1)), ('gauge_read', 'enabled', False),
    ('feedback', 'record', -1), ('feedback', 'record', 1),
    ('feedback', 'record', True), ('feedback', 'record', '0'),
    ('feedback', 'enabled', 1),
])
def test_nested_operation_fields_are_revalidated_before_preparation(monkeypatch, kind, field, value):
    operation = (protocol.Operation(kind, record=0) if kind == 'feedback'
                 else protocol.Operation(kind))
    prefix = (protocol.Operation('gauge_read'),) if kind == 'feedback' else ()
    program = protocol.Program(operations=(*prefix, operation))
    operation.__dict__[field] = value

    def forbidden(*args, **kwargs):
        raise AssertionError('invalid operation reached polynomial preparation')

    monkeypatch.setattr(s.ih, '_initial', forbidden)
    with pytest.raises(ValueError):
        s.enclose(program, dimension=2, resolution=FAST)


@pytest.mark.parametrize('component', ['program', 'operation'])
def test_caller_owned_program_cannot_change_completed_enclosure(component):
    operation = protocol.Operation('gauge_read')
    program = protocol.Program(operations=(operation,))
    enclosure = s.enclose(program, dimension=2, resolution=FAST)
    original = s.summary(enclosure)
    specification = deepcopy(enclosure.program.specification())
    if component == 'program':
        program.__dict__['native_current'] = 1
    else:
        operation.__dict__['kind'] = 'smooth_read'
    assert enclosure.program.specification() == specification
    assert s.summary(enclosure) == original
