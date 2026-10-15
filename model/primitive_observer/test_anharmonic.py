"""Independent compatibility and continuation checks for the apparatus adapter."""
from dataclasses import replace
import math

import numpy as np
import pytest

from ..anharmonic_observations import source as primitive
from ..anharmonic_transfer import source, environment
from . import anharmonic as a, runtime


@pytest.mark.parametrize('phi', [0., math.pi/20])
@pytest.mark.parametrize('steps', [0, 1, 3])
def test_complete_laboratory_law_agrees_with_existing_executor(phi, steps):
    apparatus = a.Apparatus(dimension=4, subdivisions=2, phi=phi)
    history = runtime.execute(a.program(apparatus, -.7, .13, steps))
    actual = a.readout(history)
    expected = source.run(4, -.7, .13, [steps], **apparatus.legacy_arguments())
    row = expected['records'][0]
    np.testing.assert_allclose(actual['spin_detector_joint_probabilities'],
                               row['reads'][0]['spin_detector_joint_probabilities'],
                               rtol=0, atol=2e-11)
    assert abs(actual['acceptance_probability']-row['spin_probabilities'][0]) < 2e-11
    assert tuple(branch.record for branch in history.branches) == ((0,), (1,), (2,), (3,))
    assert abs(sum(branch.probability for branch in history.branches)-1) < 2e-10
    # The spin commitment is implemented here with an additional three-pulse
    # native CNOT. Both complete spin branches remain available to the observer.
    spin, metadata = a._spin_record(4)
    assert metadata['resources']['native_pulses'] == 3
    np.testing.assert_allclose(sum(k.conj().T@k for k in spin), np.eye(8), atol=1e-12)


def test_reference_entangled_history_matches_independent_dense_kraus_execution():
    apparatus = a.Apparatus(dimension=2, subdivisions=1, phi=math.pi/20)
    program = a.program(apparatus, -.7, .13, 1)
    # Work-major/reference-minor, entangled in occupation and reference.
    ket = np.zeros(8, complex)
    ket[0], ket[3] = 1/math.sqrt(2), 1j/math.sqrt(2)
    initial = np.outer(ket, ket.conj())
    history = runtime.execute(program, initial=initial)
    rho = initial.copy()
    duration = apparatus.pulse_duration
    pos, _ = primitive.pulse_unitary(2, apparatus.delta, apparatus.omega, duration, 1)
    neg, _ = primitive.pulse_unitary(2, apparatus.delta, apparatus.omega, duration, -1)
    plus, _ = source.spin_gap(2, apparatus.phi, apparatus.delta, apparatus.omega_carrier)
    minus, _ = source.spin_gap(2, -apparatus.phi, apparatus.delta, apparatus.omega_carrier)
    z, _ = primitive.spin_rotation('Z', apparatus.theta)
    z = np.kron(z, np.eye(2))
    noise = environment.channel(2, apparatus.dephasing_rate, duration/2)
    gap_noise = environment.channel(2, apparatus.dephasing_rate, apparatus.gap_duration)

    def unitary(value, matrix):
        lifted = np.kron(matrix, np.eye(2))
        return lifted@value@lifted.conj().T

    def diffuse(value, channel):
        return sum(weight*unitary(value, np.kron(np.eye(2), phase))
                   for weight, phase in zip(channel.weights, channel.unitaries, strict=True))

    for pulse, count, gaps in ((neg, 1, (minus, plus)), (pos, 2, (plus, minus)),
                                (neg, 1, None)):
        for _ in range(count):
            rho = diffuse(unitary(diffuse(rho, noise), pulse), noise)
        if gaps:
            rho = diffuse(unitary(rho, gaps[0]), gap_noise)
            rho = unitary(rho, z)
            rho = diffuse(unitary(rho, gaps[1]), gap_noise)
    detector, _ = primitive.detector(2, **source.position_setting())
    for branch in history.branches:
        spin, bit = divmod(branch.record[0], 2)
        projection = np.diag([int(index == spin) for index in range(2)])
        expected = unitary(rho, np.kron(projection, detector[bit]))
        np.testing.assert_allclose(branch.density, expected, atol=2e-11, rtol=0)
    # The untouched reference marginal survives all complete operations.
    reference = sum(branch.density for branch in history.branches).reshape(4, 2, 4, 2)
    np.testing.assert_allclose(np.einsum('iaib->ab', reference), np.eye(2)/2, atol=2e-11)


def test_checkpoint_after_measurement_resumes_disturbed_work_without_reset():
    apparatus = a.Apparatus(dimension=4, subdivisions=2, h=1.2)
    program = a.program(apparatus, -.7, .13, 2, read_steps=(1, 2))
    prefix = runtime.execute(program, stop_after=6)
    assert len(prefix.branches) == 4
    payload = runtime.checkpoint_to_json(prefix)
    restored = runtime.checkpoint_from_json(program, payload)
    resumed = runtime.execute(program, checkpoint=restored)
    complete = runtime.execute(program)
    for left, right in zip(resumed.branches, complete.branches, strict=True):
        assert left.record == right.record
        np.testing.assert_allclose(left.density, right.density, atol=2e-12, rtol=0)
    assert len(complete.branches) == 16
    for label in range(4):
        before = prefix.branches[label].probability
        after = sum(branch.probability for branch in complete.branches if branch.record[0] == label)
        assert abs(before-after) < 2e-11
    no_earlier_read = runtime.execute(a.program(apparatus, -.7, .13, 2))
    disturbed = np.array(a.readout(complete, h=1.2)['spin_detector_joint_probabilities'])
    original = np.array(a.readout(no_earlier_read, h=1.2)['spin_detector_joint_probabilities'])
    assert np.sum(np.abs(disturbed-original))/2 > 1e-4


@pytest.mark.parametrize('changes', [{'dimension': True}, {'subdivisions': True},
                                    {'phi': float('nan')}, {'h': 0.},
                                    {'dephasing_rate': -1.}, {'finite_spin_duration': 1}])
def test_invalid_apparatus_rejected_before_cached_source_calls(changes):
    with pytest.raises(ValueError):
        a.Apparatus(**changes)


def test_time_coordinates_remain_distinct_and_reads_are_explicit():
    apparatus = a.Apparatus(dimension=2, subdivisions=1, phi=math.pi/20)
    assert apparatus.drive_step_duration > apparatus.published_step_duration
    assert a.Apparatus(phi=0.).drive_step_duration == apparatus.published_step_duration
    for read_steps in ((1, 1), (2, 1), (-1,), (True,), (3,), []):
        with pytest.raises(ValueError):
            a.program(apparatus, -.7, .13, 2, read_steps=read_steps)
    with pytest.raises(ValueError):
        a.program(replace(apparatus), -.7, .13, True)


def test_report_cannot_relabel_or_rescale_an_executed_detector():
    apparatus = a.Apparatus(dimension=2, h=1.2)
    history = runtime.execute(a.program(apparatus, -.7, .13, 0))
    report = a.readout(history)
    assert report['position_estimator'] == report['selected_mean']/(math.sqrt(2)*1.2)
    with pytest.raises(ValueError, match='differs'):
        a.readout(history, h=.4)
    generic = runtime.Program('other', 2, np.diag([1., 0.]), (
        runtime.Event('other_read', instruments=(runtime.Instrument(
            kraus=((np.diag([1., 0.]),), (np.diag([0., 1.]),)), visible=True),)),))
    with pytest.raises(ValueError, match='anharmonic'):
        a.readout(runtime.execute(generic))


def test_snapshot_reuse_preserves_separate_experiments_and_preparation_cost():
    apparatus = a.Apparatus(dimension=4, subdivisions=2, phi=math.pi/20)
    histories = tuple(a.snapshots(apparatus, -.7, .13, (0, 1, 3)))
    for step, reused in histories:
        fresh = runtime.execute(a.program(apparatus, -.7, .13, step))
        assert runtime.complete_distance(reused, fresh) < 2e-11
        assert len(reused.branches) == 4
        assert dict(reused.program.preparation_resources)['native_pulses'] > 0
        restored = runtime.checkpoint_from_json(reused.program, runtime.checkpoint_to_json(reused))
        assert runtime.complete_distance(restored, fresh) < 2e-11
