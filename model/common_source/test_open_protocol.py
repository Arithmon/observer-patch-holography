"""Independent complete-channel and hostile continuation controls."""
from dataclasses import replace
from fractions import Fraction as F
import json

import numpy as np
import pytest
from scipy.linalg import expm

from . import open_environment as env
from . import open_protocol as mixed
from . import protocol, scalar_coupling, smooth_protocol
from ..observer.source_scalar import field


def short_program(d=2):
    return protocol.Program(dimension=d, steps=8, native_current=1, operations=(
        protocol.Operation('evolve', F(1, 8)), protocol.Operation('gauge_read'),
        protocol.Operation('feedback', record=0), protocol.Operation('smooth_read'),
        protocol.Operation('evolve', F(1, 4)), protocol.Operation('gauge_read')))


def noise():
    return env.Environment((env.PhaseNoise(0, 1, F(1, 16), 4),
                            env.PhaseNoise(2, 0, F(1, 8), 4),
                            env.PhaseNoise(4, 1, F(1, 16), 4)))


def entangled(d=2):
    size = 4*d*d
    vector = np.zeros((size, 2), complex)
    vector[0, 0], vector[-1, 1], vector[5, 1] = 1/np.sqrt(2), 1j/2, 1/2
    flat = vector.ravel()
    return np.outer(flat, flat.conj())


def direct_gaussian(program, schedule, initial):
    """Full tensor matrices, direct Hamiltonian and analytic local channel.

    This oracle never calls the mixed runtime, its conjugation, the source
    environment or the canonical operation-map executor.
    """
    d, size = program.dimension, program.work_dimension
    r = len(initial)//size
    gauge_x = np.kron(np.eye(d*d), np.kron(np.eye(2), np.array([[0., 1.], [1., 0.]])))
    kick = expm(-1j*program.native_current*np.arctan2(4, 3)*gauge_x)
    kick = np.kron(kick, np.eye(r))
    result = {(): kick@initial@kick.conj().T}
    q = np.kron(np.kron(np.eye(d), field.quadratures(d)[0]), np.eye(4))
    v, u = np.linalg.eigh(q)
    smooth = tuple((u*(np.cos(v/2)+sign*np.sin(v/2))/np.sqrt(2))@u.conj().T
                   for sign in (1, -1))
    gauge = tuple(np.diag([float(i % 2 == bit) for i in range(size)]) for bit in (0, 1))

    def environment_at(cursor, rows):
        for event in schedule.events:
            if event.after == cursor:
                n = np.repeat((np.arange(size)//(4*d**event.site)) % d, r)
                matrix = np.exp(-float(event.variance)*(n[:, None]-n[None, :])**2/2)
                rows = {key: matrix*rho for key, rho in rows.items()}
        return rows

    result = environment_at(0, result)
    for cursor, operation in enumerate(program.operations, 1):
        output = {}
        for record, rho in result.items():
            if operation.kind == 'evolve':
                maps = (expm(-1j*float(operation.time)*scalar_coupling.ideal_hamiltonian(d)),)
            elif operation.kind == 'feedback':
                angle = np.arctan2(4, 3) if operation.enabled and record[operation.record] else 0.
                maps = (expm(-1j*angle*gauge_x),)
            else:
                maps = gauge if operation.kind == 'gauge_read' else smooth
            for bit, matrix in enumerate(maps):
                key = (*record, bit) if operation.kind in protocol.READS else record
                joint = np.kron(matrix, np.eye(r))
                output[key] = joint@rho@joint.conj().T
        result = environment_at(cursor, output)
    return result


@pytest.mark.parametrize('d', [2, 4])
def test_complete_reference_output_matches_independent_channel(d):
    program, schedule, initial = short_program(d), noise(), entangled(d)
    actual = mixed.execute(program, schedule, backend='continuous', noise_backend='gaussian', initial=initial)
    expected = direct_gaussian(program, schedule, initial)
    assert actual.records == tuple(expected)
    for branch in actual.branches:
        np.testing.assert_allclose(branch.density, expected[branch.record], atol=2e-13, rtol=2e-12)
    reference_before = np.trace(initial.reshape(program.work_dimension, 2, program.work_dimension, 2), axis1=0, axis2=2)
    final = sum(b.density for b in actual.branches)
    reference_after = np.trace(final.reshape(program.work_dimension, 2, program.work_dimension, 2), axis1=0, axis2=2)
    np.testing.assert_allclose(reference_after, reference_before, atol=3e-13)


@pytest.mark.parametrize('backend', ['source', 'product', 'continuous'])
def test_empty_environment_recovers_every_pure_branch(backend):
    program = short_program()
    pure = protocol.execute(program, backend=backend)
    actual = mixed.execute(program, backend=backend)
    assert actual.records == pure.records
    for a, b in zip(actual.branches, pure.branches, strict=True):
        flat = b.state.ravel()
        np.testing.assert_allclose(a.density, np.outer(flat, flat.conj()), atol=3e-13)
        assert dict(a.resources)['independent_environment_draws'] == 0
        assert {key: dict(a.resources)[key] for key in protocol.RESOURCE_KEYS} == dict(b.resources)


def test_source_environment_mixes_state_and_all_branches_survive():
    program = protocol.Program()
    vector = np.zeros(16, complex)
    vector[0] = vector[4] = 1/np.sqrt(2)
    initial = np.outer(vector, vector.conj())
    schedule = env.Environment((env.PhaseNoise(0, 0, F(1), 8),))
    result = mixed.execute(program, schedule, initial=initial)
    rho = result.branches[0].density
    assert result.probability_mass == pytest.approx(1)
    assert np.trace(rho@rho).real < .7
    assert rho[0, 4].real == pytest.approx(np.exp(-.5)/2, abs=1e-9)
    assert dict(result.branches[0].resources)['independent_environment_draws'] == 1
    # An amplitude average loses mass and cannot stand in for a mixed output.
    assert np.linalg.matrix_rank(rho, tol=1e-10) == 2


@pytest.mark.parametrize('stop', range(7))
def test_resume_preserves_state_reference_prefix_and_resources(stop):
    program, schedule, initial = short_program(), noise(), entangled()
    full = mixed.execute(program, schedule, initial=initial)
    prefix = mixed.execute(program, schedule, initial=initial, stop_after=stop)
    continued = mixed.execute(program, schedule, checkpoint=prefix)
    assert mixed.complete_distance(full, continued) < 1e-14
    assert full.expanded_tree_resources == continued.expanded_tree_resources
    assert [b.resources for b in full.branches] == [b.resources for b in continued.branches]
    prefix_length = len(prefix.records[0])
    assert mixed.prefix_probabilities(prefix, prefix_length) == pytest.approx(
        mixed.prefix_probabilities(continued, prefix_length), abs=1e-12)


def test_checkpoint_roundtrip_replays_before_continuation():
    program, schedule = short_program(), noise()
    history = mixed.execute(program, schedule, stop_after=4)
    encoded = mixed.checkpoint_to_json(history)
    recovered = mixed.checkpoint_from_json(program, schedule, encoded)
    assert mixed.checkpoint_to_json(recovered) == encoded
    assert mixed.complete_distance(mixed.execute(program, schedule, checkpoint=recovered),
                                   mixed.execute(program, schedule)) < 1e-14


@pytest.mark.parametrize('change', ['state', 'cost', 'initial', 'stages'])
def test_resealed_forgery_rejected_by_independent_import_replay(change):
    program, schedule = short_program(), noise()
    history = mixed.execute(program, schedule, stop_after=4)
    if change == 'state':
        # Keep positivity, branch probabilities and total mass exactly; alter
        # a retained quantum state by a unitary, then recompute the public seal.
        matrix = np.diag([1., -1.]*8)
        branches = tuple(replace(b, density=mixed._immutable(matrix@b.density@matrix)) for b in history.branches)
        history = replace(history, branches=branches)
    elif change == 'cost':
        cost = dict(history.branches[0].resources)
        cost['native_pulses'] += 1
        branch = replace(history.branches[0], resources=tuple(cost.items()))
        history = replace(history, branches=(branch, *history.branches[1:]))
    elif change == 'initial':
        history = replace(history, initial_density=mixed._immutable(np.eye(16)/16))
    else:
        stage = replace(history.stages[0], probability_mass=1+1e-9)
        history = replace(history, stages=(stage, *history.stages[1:]))
        # A discrepancy inside numerical replay tolerance is not a forgery
        # claim. Change resources instead to make a valid strict rejection.
        history = replace(history, preparation_resources=mixed._cost(work_preparations=2))
    history = replace(history, seal=mixed._seal(history))
    payload = mixed.checkpoint_to_json(history)
    with pytest.raises(ValueError, match='reproduce'):
        mixed.checkpoint_from_json(program, schedule, payload)


@pytest.mark.parametrize('value', [np.eye(16), np.eye(16)*np.nan, np.eye(16, dtype=bool),
                                  np.zeros((16, 15)), np.zeros((3, 3)),
                                  np.ma.array(np.eye(16)/16), [[True]]])
def test_invalid_initial_states_rejected(value):
    with pytest.raises(ValueError):
        mixed.execute(protocol.Program(), initial=value)


def test_nonpositive_and_nonhermitian_inputs_rejected():
    diagonal = np.zeros(16)
    diagonal[:2] = [1.1, -.1]
    with pytest.raises(ValueError, match='positive'):
        mixed.execute(protocol.Program(), initial=np.diag(diagonal))
    matrix = np.eye(16, dtype=complex)/16
    matrix[0, 1] = 1j/8
    with pytest.raises(ValueError, match='Hermitian'):
        mixed.execute(protocol.Program(), initial=matrix)


@pytest.mark.parametrize('change', ['environment', 'backend', 'rewind', 'input', 'branch_budget'])
def test_invalid_resume_and_allocation_rejected(change):
    program, schedule = short_program(), noise()
    prefix = mixed.execute(program, schedule, stop_after=4)
    kwargs = {'checkpoint': prefix}
    if change == 'environment':
        schedule = env.Environment(())
    elif change == 'backend':
        kwargs['noise_backend'] = 'gaussian'
    elif change == 'rewind':
        kwargs['stop_after'] = 1
    elif change == 'input':
        kwargs['initial'] = np.eye(16)/16
    else:
        kwargs['max_branches'] = 4
    with pytest.raises(ValueError):
        mixed.execute(program, schedule, **kwargs)


def test_checkpoint_bad_records_duplicate_keys_and_nonfinite_payloads_rejected():
    program, schedule = short_program(), noise()
    result = mixed.execute(program, schedule, stop_after=2)
    payload = mixed.checkpoint_to_json(result)
    value = json.loads(payload)
    value['branches'][1]['record'] = [0]
    with pytest.raises(ValueError):
        mixed.checkpoint_from_json(program, schedule, json.dumps(value))
    with pytest.raises(ValueError):
        mixed.checkpoint_from_json(program, schedule, payload.replace('{', '{"cursor":2,', 1))
    value = json.loads(payload)
    value['branches'][0]['density']['values'][0][0] = float('nan')
    with pytest.raises(ValueError):
        mixed.checkpoint_from_json(program, schedule, json.dumps(value))


def test_outputs_do_not_alias_caller_or_allow_numpy_writes():
    program, schedule, initial = short_program(), noise(), entangled()
    history = mixed.execute(program, schedule, initial=initial)
    original = mixed.checkpoint_to_json(history)
    initial[:] = 0
    program.operations[0].__dict__['time'] = F(99)
    schedule.events[0].__dict__['variance'] = F(99)
    assert mixed.checkpoint_to_json(history) == original
    with pytest.raises(ValueError):
        history.branches[0].density.flags.writeable = True


def test_resource_draws_are_paths_not_visible_or_hidden_tree_size():
    program = protocol.Program(operations=(protocol.Operation('gauge_read'), protocol.Operation('gauge_read')))
    schedule = env.Environment((env.PhaseNoise(1, 0, F(0), 4), env.PhaseNoise(2, 0, F(0), 4)))
    result = mixed.execute(program, schedule)
    assert {dict(b.resources)['independent_environment_draws'] for b in result.branches} == {2}
    assert dict(result.expanded_tree_resources)['independent_environment_draws'] == 6
    baseline = mixed.execute(program)
    for noisy, quiet in zip(result.branches, baseline.branches, strict=True):
        assert dict(noisy.resources)['native_pulses']-dict(quiet.resources)['native_pulses'] == 4
    assert mixed.complete_distance(result, baseline) < 1e-12
