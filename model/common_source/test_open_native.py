"""Independent endpoint, mixed-output and reference controls for the attachment."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from ..anharmonic_transfer import environment as phase_source
from ..observer.source_discrimination import circuit
from . import encoded_memory as memory
from . import encoded_observer
from . import open_environment as noise
from . import open_native as joint
from . import protocol, quantum_native as pure


SUPPORTS = (((0, 2), (1, 3)), ((0, 3), (1, 2)))


def configuration(pair):
    return tuple(int(port in pair) for port in range(12))


def endpoint(state, seam, coin, law):
    a, b = seam
    low, parity = divmod(state[a]+state[b], 2)
    final = list(state)
    final[a], final[b] = low+coin*parity, low+(1-coin)*parity
    current = state[a]-final[a]
    raw = {0: 0, 1: 1, -1: 2}[current]+(4*parity if law == 'E' else 0)
    return tuple(final), raw


def code():
    result = np.zeros((512, 2), complex)
    for logical, pairs in enumerate(SUPPORTS):
        for pair in pairs:
            result[circuit.encode(configuration(pair)), logical] = 1/math.sqrt(2)
    return result


def weight(law, seam, coin, raw):
    counts = [sum(endpoint(configuration(pair), seam, coin, law)[1] == raw
                  for pair in pairs) for pairs in SUPPORTS]
    assert counts[0] == counts[1]
    return F(counts[0], 4)


def expected_native(record, law, protected):
    c1, raw1, bits, c2, raw2 = record
    result = np.zeros((512, 2), complex)
    if not bits:
        return result
    seam = (2, 3) if bits[2] == 0 else (0, 1)
    if protected:
        return math.sqrt(float(weight(law, (0, 1), c1, raw1)
                               *weight(law, seam, c2, raw2)))*code()
    # Direct two-event endpoint action, with both actual fair-coin factors.
    for logical, pairs in enumerate(SUPPORTS):
        for pair in pairs:
            middle, actual1 = endpoint(configuration(pair), (0, 1), c1, law)
            final, actual2 = endpoint(middle, seam, c2, law)
            if actual1 == raw1 and actual2 == raw2:
                result[circuit.encode(final), logical] += 1/(2*math.sqrt(2))
    return result


def current(law, coin, raw):
    j = {0: 0, 1: 1, 2: -1}.get(raw % 4)
    return None if (j is None or law == 'E' and abs(2*j+(2*coin-1)*(raw//4)) > 2) else j


def complete_menu(law):
    result = {}
    for c1, raw1 in product((0, 1), range(4 if law == 'C' else 8)):
        if current(law, c1, raw1) is None:
            result[c1, raw1, (), None, None] = True
        else:
            for bits in product((0, 1), repeat=4):
                for c2, raw2 in product((0, 1), range(4 if law == 'C' else 8)):
                    result[c1, raw1, bits, c2, raw2] = current(law, c2, raw2) is None
    return result


@pytest.fixture(scope='module')
def environment():
    return noise.Environment((noise.PhaseNoise(4, 0, F(1, 8), 4),
                              noise.PhaseNoise(6, 1, F(1, 8), 4)))


@pytest.fixture(scope='module', params=('C', 'E'))
def histories(request, environment):
    return tuple(joint.execute(environment, request.param, 2, 4, native_mode=mode)
                 for mode in ('protected_qubit', 'passive_qubit'))


def test_complete_raw_native_menu_and_independent_endpoint_maps(histories):
    for h in histories:
        assert {b.record: b.abort for b in h.branches} == complete_menu(h.law)
        assert len(h.branches) == (770 if h.law == 'C' else 2566)
        checked = set()
        for b in h.branches:
            key = (*b.record[:2], b.record[2][2] if b.record[2] else None, *b.record[3:])
            if key not in checked:
                np.testing.assert_allclose(b.native_map, expected_native(
                    b.record, h.law, h.native_mode == 'protected_qubit'), atol=2e-11, rtol=0)
                checked.add(key)
            assert b.native_map.shape == (512, 2)
            assert b.work_density.shape == (16, 16)
        report = joint.audit(h)
        assert report['prepared_work_native_input_completeness_error_diagnostic'] < 2e-10
        assert report['floating_or_physical_error_certified'] is False


def test_native_reference_output_preserves_complete_entanglement(histories):
    h = histories[0]
    # Complex, nonmaximal entanglement and a spectator reference column.
    initial = np.array([[1, 2j, 0], [1-1j, 0, 3]], complex)
    initial /= np.linalg.norm(initial)
    encoded = code()@initial
    # Only four code rows enter the overlap with the independent target;
    # the full 512 rows remain in the actual streamed output and mass check.
    target_rows = np.flatnonzero(np.linalg.norm(encoded, axis=1))
    reduced = np.zeros((len(target_rows)*3,)*2, complex)
    mass = 0.
    for _, amplitude, density in joint.output(h, initial):
        p = np.trace(density).real
        vector = amplitude[target_rows].ravel()
        reduced += p*np.outer(vector, vector.conj())
        mass += p*np.vdot(amplitude, amplitude).real
        assert amplitude.shape == (512, 3)
        assert not amplitude.flags.writeable and not density.flags.writeable
    target = encoded[target_rows].ravel()
    np.testing.assert_allclose(reduced, np.outer(target, target.conj()), atol=2e-10, rtol=0)
    assert abs(mass-1) < 2e-10
    assert abs(joint.logical_entanglement_fidelity(h)-1) < 2e-10


def test_passive_loss_and_active_memory_use_actual_noisy_work_weights(histories):
    active, passive = histories
    probabilities = {h.program.native_current: math.fsum(np.trace(b.density).real
                     for b in h.branches if b.record[2] == 0) for h in active.work_histories}
    qbar = (probabilities[-1]+2*probabilities[0]+probabilities[1])/4
    assert .1 < qbar < .9
    witness = joint.memory_witness(active, passive)
    assert abs(witness['second_smooth_plus_probability_diagnostic']-qbar) < 2e-12
    assert abs(joint.phase_distance(passive)-(1-qbar)) < 2e-10
    assert abs(joint.phase_distance(active)-1) < 2e-10
    assert joint.logical_entanglement_fidelity(passive) < .9
    assert witness['native_storage_or_record_fault_tolerance_certified'] is False


def test_actual_density_is_mixed_and_never_renormalized(histories):
    h = histories[0]
    ranks, masses = [], []
    for worker in h.work_histories:
        for b in worker.branches:
            values = np.linalg.eigvalsh(b.density)
            ranks.append(sum(values > 1e-8))
            masses.append(np.trace(b.density).real)
    assert max(ranks) > 1
    assert max(masses) < .5
    assert abs(math.fsum(masses)-3) < 2e-10
    # Every raw native continuation shares the corresponding retained density.
    workers = {w.program.native_current: {b.record: b.density for b in w.branches}
               for w in h.work_histories}
    for b in h.branches:
        expected = (h.work_histories[1].initial_density if not b.record[2]
                    else workers[current(h.law, *b.record[:2])][b.record[2]])
        assert b.work_density is expected


def test_terminal_noise_is_cp_sum_not_average_amplitude():
    event = noise.PhaseNoise(10, 0, F(1, 4), 4)
    h = joint.execute(noise.Environment((event,)), 'C', 2, 4)
    worker = h.work_histories[2]
    quiet = protocol.execute(worker.program)
    channel = phase_source.channel(2, 1., .25, 4)
    occupation = np.arange(16)//4 % 2
    wrong_mass, difference = 0., 0.
    for actual, old in zip(worker.branches, quiet.branches, strict=True):
        assert actual.record == old.record
        alternatives = [np.diag(u)[occupation]*old.state.ravel() for u in channel.unitaries]
        correct = sum(p*np.outer(a, a.conj()) for p, a in zip(channel.weights, alternatives, strict=True))
        average = sum(p*a for p, a in zip(channel.weights, alternatives, strict=True))
        wrong = np.outer(average, average.conj())
        np.testing.assert_allclose(actual.density, correct, atol=2e-11, rtol=0)
        wrong_mass += np.trace(wrong).real
        difference += np.linalg.norm(correct-wrong, 'fro')
    assert wrong_mass < .98 and difference > .01


def test_empty_noise_matches_prior_complete_pure_history():
    h = joint.execute(law='C', dimension=2, steps=4)
    old = pure.execute('C', 2, 4, native_mode='protected_qubit')
    for new, prior in zip(h.branches, old.branches, strict=True):
        assert new.record == prior.record and new.abort == prior.abort
        np.testing.assert_allclose(new.native_map, prior.native_map, atol=0, rtol=0)
        np.testing.assert_allclose(new.work_density, prior.work_state@prior.work_state.conj().T,
                                   atol=2e-11, rtol=0)
    # Unlike the old diagnostic receipt, this stops before input-specific
    # three-pulse preparation, terminal source decoding and destructive read.
    old_cost = encoded_observer.resources('C', 2, 4)['maximum_native_pulses_per_attempt']
    new_cost = joint.resources(h)['maximum_realized_path_resources']['native_pulses']
    assert old_cost-new_cost == 3+memory.decoding_resources()['native_pulses']+1


def test_full_carrier_alternative_still_retains_78_input_columns(environment):
    h = joint.execute(environment, 'C', 2, 4, native_mode='full_carrier')
    report = joint.audit(h)
    assert report['prepared_work_native_input_completeness_error_diagnostic'] < 2e-10
    initial = np.zeros((78, 2), complex)
    initial[0, 0], initial[-1, 1] = 1/math.sqrt(2), 1j/math.sqrt(2)
    mass = math.fsum(float(np.vdot(a, a).real*np.trace(rho).real)
                     for _, a, rho in joint.output(h, initial))
    assert abs(mass-1) < 2e-10
    with pytest.raises(ValueError):
        joint.logical_entanglement_fidelity(h)


def test_resources_pay_every_noise_invocation_but_no_diagnostic_read(histories):
    h = histories[0]
    resources = joint.resources(h)
    counts = resources['maximum_realized_path_resources']
    assert counts['work_preparations'] == 1
    assert counts['fresh_pointer_blanks'] == counts['committed_classical_bits'] == 4
    assert counts['independent_environment_draws'] == len(h.environment.events)
    assert counts['proper_code_service_events'] == 2*counts['native_pulses']
    assert resources['native_encoder_count'] == 1
    assert resources['native_logical_storage_through_work_required'] is True
    assert resources['diagnostic_input_preparation_and_terminal_read_pulses'] == 0
    assert resources['retained_classical_bits_including_two_status_bits'] == (12 if h.law == 'C' else 14)
    assert resources['environment_draws_are_observer_records'] is False


@pytest.mark.parametrize('mutation', ('omit_abort', 'duplicate', 'relabel', 'phase',
                                      'wrong_mixed_branch', 'normalized_branch', 'masked_density'))
def test_integrity_checks_reject_changed_complete_output(histories, mutation):
    h = histories[0]
    rows = list(h.branches)
    index = next(i for i, b in enumerate(rows) if b.record[2]
                 and np.linalg.norm(b.native_map) > .1 and np.trace(b.work_density).real > .01)
    b = rows[index]
    if mutation == 'omit_abort':
        rows = [b for b in rows if b.record[2]]
    elif mutation == 'duplicate':
        rows.append(rows[-1])
    elif mutation == 'relabel':
        rows[index] = replace(b, record=(bool(b.record[0]), *b.record[1:]))
    elif mutation == 'phase':
        rows[index] = replace(b, native_map=1j*b.native_map)
    elif mutation == 'wrong_mixed_branch':
        other = next(r.work_density for r in rows if r.record[2] and r.record[2] != b.record[2])
        rows[index] = replace(b, work_density=other)
    elif mutation == 'normalized_branch':
        rows[index] = replace(b, work_density=b.work_density/np.trace(b.work_density))
    else:
        rows[index] = replace(b, work_density=np.ma.array(b.work_density, mask=False))
    with pytest.raises(ValueError):
        joint.output(replace(h, branches=tuple(rows)), [1, 0])


@pytest.mark.parametrize('bad', ([True, 0], [2, 0], [float('nan'), 0], np.ones((2, 2, 2)),
                                np.ma.array([1., 0.], mask=[False, True])))
def test_native_input_contract_rejects_invalid_or_joint_work_payload(histories, bad):
    with pytest.raises(ValueError):
        joint.output(histories[0], bad)


def test_frozen_environment_mutation_is_rejected_before_execution(monkeypatch):
    env = noise.Environment((noise.PhaseNoise(0, 0, F(1, 8), 4),))
    object.__setattr__(env.events[0], 'variance', F(-1))
    monkeypatch.setattr(joint.mixed, 'execute', lambda *a, **k: pytest.fail('invalid input executed'))
    with pytest.raises(ValueError):
        joint.execute(env, steps=4)


@pytest.mark.parametrize('mutable_arrays', (False, True))
def test_stream_snapshots_frozen_envelope_and_only_detaches_mutable_arrays(histories, mutable_arrays):
    h = histories[0]
    original = next(b for b in h.branches if b.record[2] and np.linalg.norm(b.native_map) > .1)
    branch = replace(original, native_map=original.native_map.copy() if mutable_arrays else original.native_map.view(),
                     work_density=original.work_density.copy() if mutable_arrays else original.work_density.view())
    changed = replace(h, branches=tuple(branch if b is original else b for b in h.branches))
    # Equal-content replacement retains the existing digest. Snapshot before
    # mutation, then consume lazily; neither metadata nor storage may leak in.
    rows = joint.output(changed, [1, 0])
    if mutable_arrays:
        branch.native_map[:] = 0
        branch.work_density[:] = 0
    else:
        branch.native_map.shape = (-1,)
        branch.work_density.shape = (-1,)
    object.__setattr__(branch, 'record', ('changed',))
    object.__setattr__(branch, 'native_map', np.zeros_like(original.native_map))
    actual = next((a, rho) for record, a, rho in rows if record == original.record)
    np.testing.assert_allclose(actual[0], original.native_map[:, :1], atol=0, rtol=0)
    np.testing.assert_allclose(actual[1], original.work_density, atol=0, rtol=0)
    if not mutable_arrays:
        assert np.shares_memory(actual[1], original.work_density)
    assert not actual[1].flags.writeable


def test_first_abort_retains_pre_kick_pre_noise_work():
    env = noise.Environment((noise.PhaseNoise(0, 0, F(1, 2), 4),))
    h = joint.execute(env, steps=4)
    expected = h.work_histories[1].initial_density
    for b in h.branches:
        if not b.record[2]:
            assert b.work_density is expected
            assert b.record[3:] == (None, None) and b.abort
    assert np.linalg.norm(h.work_histories[1].branches[0].density-expected) > .01
