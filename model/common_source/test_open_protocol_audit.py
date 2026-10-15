"""Independent mixed-history transport and public-boundary audit controls."""
from dataclasses import replace
from fractions import Fraction as F
import json

import numpy as np
import pytest
from scipy.linalg import expm

from . import open_environment as env
from . import open_protocol as mixed
from . import protocol, scalar_coupling


@pytest.mark.parametrize('field', [
    'steps', 'native_current', 'feedback_record', 'feedback_enabled',
    'environment_after', 'environment_site', 'environment_order',
])
def test_serialized_identity_preserves_original_scalar_types(field):
    program = protocol.Program(steps=1, operations=(
        protocol.Operation('gauge_read'), protocol.Operation('feedback', record=0)))
    schedule = env.Environment((env.PhaseNoise(0, 0, F(0), 1),))
    history = mixed.execute(program, schedule)
    payload = mixed.checkpoint_to_json(history)
    healthy = mixed.checkpoint_from_json(program, schedule, payload)
    assert mixed.complete_distance(history, healthy) < 1e-14
    value = json.loads(payload)
    if field.startswith('environment_'):
        name = field.removeprefix('environment_')
        value['environment']['events'][0][name] = name == 'order'
    elif field.startswith('feedback_'):
        name = field.removeprefix('feedback_')
        value['program']['operations'][1][name] = 1 if name == 'enabled' else False
    else:
        value['program'][field] = field == 'steps'
    with pytest.raises(ValueError):
        mixed.checkpoint_from_json(program, schedule, json.dumps(value))


def _input(d, reference=3):
    """Two mixed, complex, non-product vectors with a non-qubit reference."""
    size = 4*d*d
    first, second = np.zeros((2, size*reference), complex)
    first[0], first[(size-1)*reference+1] = 1/np.sqrt(2), 1j/np.sqrt(2)
    second[5*reference+2], second[3*reference] = np.sqrt(.3), np.sqrt(.7)*np.exp(.4j)
    return .3*np.outer(first, first.conj())+.7*np.outer(second, second.conj())


def _program(d=2):
    return protocol.Program(dimension=d, native_current=-1, coupling=F(-2, 7), steps=3,
        operations=(protocol.Operation('gauge_read'),
                    protocol.Operation('feedback', record=0),
                    protocol.Operation('evolve', F(1, 7)),
                    protocol.Operation('smooth_read'),
                    protocol.Operation('feedback', record=1, enabled=False),
                    protocol.Operation('sharp_read'),
                    protocol.Operation('feedback', record=2),
                    protocol.Operation('evolve', F(1, 11))))


def _finite_reference(program, schedule, initial):
    """Full tensor products and explicit two-point environment mixtures.

    No runtime conjugation, noise Schur matrices, operation maps or output
    probabilities enter the oracle. Order two has nodes -1,+1, weights 1/2.
    Product evolution uses the independent matrix exponential implementation.
    """
    d, size = program.dimension, program.work_dimension
    reference = len(initial)//size
    x = np.zeros((size, size))
    x[np.arange(size)^1, np.arange(size)] = 1
    kick = expm(-1j*program.native_current*np.arctan2(4, 3)*x)
    lower = np.diag(np.sqrt(np.arange(1, d)), 1)
    coordinate = np.kron(np.kron(np.eye(d), (lower+lower.T)/np.sqrt(2)), np.eye(4))
    values, vectors = np.linalg.eigh(coordinate)
    sharp = tuple((vectors*((values > 0) == bit))@vectors.conj().T for bit in (0, 1))
    minus, plus = expm(-.5j*coordinate), expm(.5j*coordinate)
    smooth = ((minus-1j*plus)/2, (minus+1j*plus)/2)
    gauge = tuple(np.diag((np.arange(size)%2 == bit).astype(float)) for bit in (0, 1))

    def apply(matrix, rho):
        full = np.kron(matrix, np.eye(reference))
        return full@rho@full.conj().T

    def after(cursor, rows):
        for event in schedule.events:
            if event.after != cursor:
                continue
            assert event.order == 2
            occupation = (np.arange(size)//(4*d**event.site)) % d
            maps = tuple(np.diag(np.exp(sign*1j*np.sqrt(float(event.variance))*occupation))
                         for sign in (-1, 1))
            rows = {record: sum(apply(matrix, rho) for matrix in maps)/2
                    for record, rho in rows.items()}
        return rows

    rows = after(0, {(): apply(kick, initial)})
    for cursor, operation in enumerate(program.operations, 1):
        following = {}
        for record, rho in rows.items():
            if operation.kind == 'evolve':
                maps = (scalar_coupling.independent_split_unitary(
                    d, program.steps, float(operation.time), float(program.coupling)),)
            elif operation.kind == 'feedback':
                angle = np.arctan2(4, 3) if operation.enabled and record[operation.record] else 0.
                maps = (expm(-1j*angle*x),)
            else:
                maps = {'gauge_read': gauge, 'smooth_read': smooth, 'sharp_read': sharp}[operation.kind]
            for bit, matrix in enumerate(maps):
                following[(*record, bit) if operation.kind in protocol.READS else record] = apply(matrix, rho)
        rows = after(cursor, following)
    return rows


@pytest.mark.parametrize('d', [2, 4])
def test_all_source_operations_with_mixed_three_level_reference(d):
    program = _program(d)
    schedule = env.Environment((env.PhaseNoise(0, 0, F(1, 5), 2),
                                env.PhaseNoise(2, 1, F(2, 7), 2),
                                env.PhaseNoise(4, 0, F(1, 5), 2),
                                env.PhaseNoise(8, 1, F(1, 9), 2)))
    initial = _input(d)
    expected = _finite_reference(program, schedule, initial)
    history = mixed.execute(program, schedule, initial=initial)
    assert history.reference_dimension == 3
    assert history.records == tuple(expected)
    for branch in history.branches:
        np.testing.assert_allclose(branch.density, expected[branch.record], atol=4e-13, rtol=4e-12)
    reduced = lambda rho: np.trace(rho.reshape(program.work_dimension, 3,
                                             program.work_dimension, 3), axis1=0, axis2=2)
    np.testing.assert_allclose(reduced(sum(b.density for b in history.branches)),
                               reduced(initial), atol=8e-13)


def test_complete_flagged_choi_output_matches_independent_finite_instrument():
    program = _program()
    schedule = env.Environment((env.PhaseNoise(3, 0, F(1, 3), 2),))
    # Reference dimension equal to work dimension tests the complete finite
    # instrument, including off-diagonal input matrix units, in one input.
    vector = np.eye(program.work_dimension, dtype=complex).ravel()/np.sqrt(program.work_dimension)
    initial = np.outer(vector, vector.conj())
    expected = _finite_reference(program, schedule, initial)
    actual = mixed.execute(program, schedule, initial=initial)
    assert actual.records == tuple(expected)
    assert actual.reference_dimension == program.work_dimension
    for branch in actual.branches:
        np.testing.assert_allclose(branch.density, expected[branch.record], atol=4e-13, rtol=4e-12)


def test_repeated_source_draws_are_independent_and_not_amplitude_averaged():
    # Two identical phase services separated by an identity observer action.
    # The correct coherence is cos(sqrt(v))**2/2; one reused draw gives
    # cos(2*sqrt(v))/2. Averaging amplitudes would also change diagonal mass.
    program = protocol.Program(operations=(protocol.Operation('evolve', F(0)),))
    variance = F(1, 2)
    schedule = env.Environment((env.PhaseNoise(0, 0, variance, 2),
                                env.PhaseNoise(1, 0, variance, 2)))
    vector = np.zeros(16, complex)
    vector[0], vector[4] = 1/np.sqrt(2), 1j/np.sqrt(2)
    history = mixed.execute(program, schedule, initial=np.outer(vector, vector.conj()))
    actual = history.branches[0].density
    expected = -.5j*np.cos(np.sqrt(.5))**2
    assert actual[0, 4] == pytest.approx(expected, abs=3e-14)
    assert abs(actual[0, 4]-(-.5j*np.cos(2*np.sqrt(.5)))) > .1
    np.testing.assert_allclose(actual.diagonal()[[0, 4]], [.5, .5], atol=2e-14)
    assert dict(history.branches[0].resources)['independent_environment_draws'] == 2


@pytest.mark.parametrize('probability', [0., 1e-250, 1e-15, .25])
def test_zero_and_tiny_visible_blocks_retain_their_original_scale(probability):
    program = protocol.Program(operations=(protocol.Operation('gauge_read'),
                                          protocol.Operation('gauge_read')))
    rho = np.zeros((16, 16))
    rho[0, 0], rho[1, 1] = 1-probability, probability
    history = mixed.execute(program, backend='product', initial=rho)
    assert history.records == ((0, 0), (0, 1), (1, 0), (1, 1))
    assert history.branches[1].probability == history.branches[2].probability == 0
    assert history.branches[3].probability == pytest.approx(probability, rel=1e-14, abs=0)
    assert history.branches[3].density[1, 1] == pytest.approx(probability, rel=1e-14, abs=0)
    restored = mixed.checkpoint_from_json(program, env.Environment(), mixed.checkpoint_to_json(history),
                                         backend='product')
    assert restored.branches[3].probability == history.branches[3].probability


@pytest.mark.parametrize('defect', ['negative', 'nonhermitian', 'nonfinite'])
def test_tiny_bad_blocks_cannot_hide_under_total_mass_tolerance(defect):
    program = protocol.Program(operations=(protocol.Operation('gauge_read'),))
    initial = np.zeros((16, 16))
    initial[0, 0] = 1
    history = mixed.execute(program, initial=initial, backend='product')
    bad = np.zeros((16, 16), complex)
    if defect == 'negative':
        bad[0, 0], bad[1, 1] = 1e-250, -1e-251
    elif defect == 'nonhermitian':
        bad[0, 1] = 1e-250j
    else:
        bad[0, 0] = np.nan
    branch = replace(history.branches[1], density=mixed._immutable(bad))
    history = replace(history, branches=(history.branches[0], branch))
    history = replace(history, seal=mixed._seal(history))
    with pytest.raises((ValueError, ArithmeticError)):
        mixed.execute(program, checkpoint=history, backend='product')


def test_one_preparation_and_exact_node_resource_count_through_resume(monkeypatch):
    calls = []
    prepare = scalar_coupling.preparation
    def counted(d):
        calls.append(d)
        return prepare(d)
    monkeypatch.setattr(scalar_coupling, 'preparation', counted)
    program = protocol.Program(steps=2, operations=(
        protocol.Operation('gauge_read'), protocol.Operation('feedback', record=0),
        protocol.Operation('smooth_read'), protocol.Operation('sharp_read')))
    schedule = env.Environment((env.PhaseNoise(0, 0, F(0), 2),
                                env.PhaseNoise(1, 0, F(1, 9), 2),
                                env.PhaseNoise(1, 1, F(0), 2),
                                env.PhaseNoise(3, 1, F(1, 9), 2),
                                env.PhaseNoise(4, 0, F(0), 2)))
    history = mixed.execute(program, schedule, stop_after=0)
    for cursor in range(1, len(program.operations)+1):
        history = mixed.execute(program, schedule, checkpoint=history, stop_after=cursor)
    assert calls == [2]
    assert {dict(b.resources)['work_preparations'] for b in history.branches} == {1}
    assert dict(history.expanded_tree_resources)['work_preparations'] == 1
    assert {dict(b.resources)['independent_environment_draws'] for b in history.branches} == {5}
    # Events after 0,1,1,3,4 see respectively 1,2,2,4,8 visible nodes.
    assert dict(history.expanded_tree_resources)['independent_environment_draws'] == 17
    # Read nodes occur before branching: 1 gauge + 2 smooth + 4 sharp.
    assert dict(history.expanded_tree_resources)['fresh_pointer_blanks'] == 7
    assert dict(history.expanded_tree_resources)['committed_classical_bits'] == 7
    assert {dict(b.resources)['committed_classical_bits'] for b in history.branches} == {3}
    quiet = mixed.execute(program)
    assert dict(history.expanded_tree_resources)['native_pulses']-dict(quiet.expanded_tree_resources)['native_pulses'] == 34
    assert {dict(noisy.resources)['native_pulses']-dict(baseline.resources)['native_pulses']
            for noisy, baseline in zip(history.branches, quiet.branches, strict=True)} == {10}


@pytest.mark.parametrize('layer', ['program', 'operation', 'environment', 'event', 'branch', 'stage'])
def test_no_op_resume_detaches_every_caller_owned_envelope(layer):
    program = protocol.Program(operations=(protocol.Operation('evolve', F(1, 13)),))
    schedule = env.Environment((env.PhaseNoise(0, 0, F(1, 7), 2),))
    original = mixed.execute(program, schedule)
    retained = mixed.execute(program, schedule, checkpoint=original)
    before = mixed.checkpoint_to_json(retained)
    if layer == 'program':
        original.program.__dict__['native_current'] = 99
    elif layer == 'operation':
        original.program.operations[0].__dict__['time'] = F(-1)
    elif layer == 'environment':
        original.environment.__dict__['events'] = ()
    elif layer == 'event':
        original.environment.events[0].__dict__['site'] = 99
    elif layer == 'branch':
        original.branches[0].__dict__['density'] = np.zeros((16, 16))
    else:
        original.stages[0].__dict__['probability_mass'] = 0
    assert mixed.checkpoint_to_json(retained) == before


def test_complete_reference_budget_is_checked_before_noise_execution(monkeypatch):
    program = protocol.Program(operations=tuple(protocol.Operation('gauge_read') for _ in range(3)))
    schedule = env.Environment((env.PhaseNoise(0, 0, F(1, 8), 2),))
    # A small test cap reproduces the same arithmetic as an oversized actual
    # work/reference register, without allocating that register or its tree.
    monkeypatch.setattr(mixed, 'MAX_DENSITY_ENTRIES', 1024)
    initial = np.eye(16)/16
    assert initial.size <= mixed.MAX_DENSITY_ENTRIES < 8*initial.size
    def forbidden(*args, **kwargs):
        pytest.fail('source environment executed before rejecting the declared complete allocation')
    monkeypatch.setattr(env, 'apply', forbidden)
    with pytest.raises(ValueError, match='allocation'):
        mixed.execute(program, schedule, initial=initial)


@pytest.mark.parametrize('corruption', ['reference_order', 'selected_noise', 'feedback_bit'])
def test_independent_oracle_detects_wrong_but_complete_channels(monkeypatch, corruption):
    program = _program()
    schedule = env.Environment((env.PhaseNoise(2, 0, F(1, 2), 2),))
    initial = _input(2)
    expected = _finite_reference(program, schedule, initial)
    healthy = mixed.execute(program, schedule, initial=initial)
    for branch in healthy.branches:
        np.testing.assert_allclose(branch.density, expected[branch.record], atol=4e-13, rtol=4e-12)
    if corruption == 'reference_order':
        def wrong_reference(matrix, rho, reference_dimension):
            joint = np.kron(np.eye(reference_dimension), matrix)
            return joint@rho@joint.conj().T
        monkeypatch.setattr(mixed, 'conjugate', wrong_reference)
    elif corruption == 'selected_noise':
        def select_draw(rho, program, event, backend='source'):
            occupation = (np.arange(program.work_dimension)//(4*program.dimension**event.site)) % program.dimension
            phases = np.repeat(np.exp(-1j*np.sqrt(float(event.variance))*occupation),
                               len(rho)//program.work_dimension)
            return phases[:, None]*rho*phases.conj()[None, :]
        monkeypatch.setattr(env, 'apply', select_draw)
    else:
        maps = protocol._operation_maps
        def flip_feedback(program, operation, selected, backend):
            if operation.kind == 'feedback':
                selected = 1-selected
            return maps(program, operation, selected, backend)
        monkeypatch.setattr(protocol, '_operation_maps', flip_feedback)
    corrupted = mixed.execute(program, schedule, initial=initial)
    # Each substitution is itself trace preserving. Normalization and the
    # record menu cannot replace the independent full-output comparison.
    assert corrupted.probability_mass == pytest.approx(1, abs=2e-12)
    assert corrupted.records == healthy.records
    maximum = max(np.linalg.norm(branch.density-expected[branch.record], 2)
                  for branch in corrupted.branches)
    assert maximum > 1e-3
