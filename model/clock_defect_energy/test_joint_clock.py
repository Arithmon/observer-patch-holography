"""Joint correlations, full quantum instrument, and changed-premise controls."""
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import joint_clock as joint


@pytest.fixture(scope='module', params=('C', 'E'))
def law_frames(request):
    return request.param, joint.clock_frames(request.param, raw=False)


def test_native_arc_probabilities_are_existing_exact_continuation_not_fitted():
    assert joint.field_probabilities(1) == (F(21, 22), F(1, 22))
    assert joint.field_probabilities(2) == (F(91, 110), F(19, 110))
    assert joint.field_probabilities(2, 1) != joint.field_probabilities(2)
    assert joint.wave.norm_squared(joint.wave.embed(joint.FIELD)) == 330


def test_both_laws_retained_joint_law_and_identical_marginal_countercontrol(law_frames):
    law, _ = law_frames
    rows = joint.exact_clock_rows(law)
    grouped = {n: joint.ZERO for n in (1, 2)}
    joint_dist = {}
    for key, coefficients in rows.items():
        n = key[-1]
        grouped[n] = joint._plus(grouped[n], coefficients)
        for sign, p in zip((1, -1), joint.field_probabilities(n), strict=True):
            joint_dist[key, sign] = joint._scale(coefficients, p)
    assert grouped == {1: (F(307, 360), F(), F()), 2: (F(53, 360), F(), F())}
    fplus = sum((grouped[n][0]*joint.field_probabilities(n)[0] for n in (1, 2)), F())
    control = {(key, sign): joint._scale(coefficients, fplus if sign == 1 else 1-fplus)
               for key, coefficients in rows.items() for sign in (1, -1)}
    for key, value in rows.items():
        assert joint._plus(control[key, 1], control[key, -1]) == value
        assert joint._plus(joint_dist[key, 1], joint_dist[key, -1]) == value
    for sign in (1, -1):
        sums = []
        for table in (joint_dist, control):
            total = joint.ZERO
            for (_, actual_sign), value in table.items():
                if actual_sign == sign:
                    total = joint._plus(total, value)
            sums.append(total)
        assert sums[0] == sums[1]
    assert joint_dist != control
    report = joint.exact_joint(law)
    assert F(report['independent_field_resampling_control']['joint_total_variation_exact']) > F(3, 100)
    assert F(report['clock_negative_current_field_plus_covariance_interval'][1]) < -F(6, 1000)
    assert F(report['one_phase_word_control_complete_clock_TV_interval'][0]) > F(1, 10000)


def test_all_clock_operator_effects_and_every_ideal_abort(law_frames):
    _, frames = law_frames
    effect = sum((f.matrix.conj().T@f.matrix for f in frames), np.zeros((66, 66), complex))
    np.testing.assert_allclose(effect, np.eye(66), atol=2e-14, rtol=0)
    assert all(np.count_nonzero(f.matrix) == 0 for f in frames if f.abort)
    assert {f.ticks for f in frames if not f.abort} == {1, 2}


def test_actual_native_branch_ratios_reuse_transport_and_available_s4_words(law_frames):
    law, _ = law_frames
    admitted = {joint.clock.ratios()}
    for candidate in ('C', 'E'):
        admitted.update(joint.clock.ratios(b, law=candidate) for b in joint.clock.branches(candidate))
    for first in joint.clock.transport.word_instrument(law, (joint.clock.FIRST_EDGE,)):
        assert joint.clock.ratios(first, law=law) in admitted
    sigma = joint.clock.transport.reference()
    for branch in joint.clock.branches(law):
        extension = joint.clock.transport.extended_reference(law, branch, sigma)
        assert all(extension[j] == sigma[i] for i, j in zip(branch.inputs, branch.outputs, strict=True))


def test_every_raw_native_register_word_has_its_zero_or_record_block(law_frames):
    law, frames = law_frames
    certificate = joint.clock.raw_alphabet_certificate(law, {f.key: f.matrix for f in frames})
    assert certificate['attained_semantic_record_slots'] == sum(not f.abort for f in frames)
    assert certificate['zero_blocks_are_defined_on_work_and_untouched_reference']


def test_interleaved_transport_equals_count_dependent_input_clock_on_all_columns(law_frames):
    law, frames = law_frames
    codes = np.array(joint.clock.full.support()[0])
    native_part = lambda array: array[np.ix_(codes, codes)]
    probe = native_part(joint.clock.rotation(joint.clock.probe_pair(), 1, False))
    branches = {b.history: b for b in joint.clock.branches(law)}
    reads = {(r.coin, r.label, r.current): r
             for r in joint.clock.interface.fixed_maps(law, joint.clock.READ_EDGE)}
    for frame in frames:
        if frame.abort:
            continue
        _, history, coin, label, current = frame.key
        branch = branches.get(history)
        expected = np.zeros((66, 66), complex)
        if branch is not None:
            k = np.zeros((66, 66), complex)
            for i, j in zip(branch.inputs, branch.outputs, strict=True):
                k[j, i] = .5
            initial_clock = native_part(joint.clock.clock(F(4*frame.ticks), joint.clock.ratios(), False))
            expected = k@initial_clock@probe
            if all(i in branch.inputs for i in joint.clock.probe_pair()):
                pair = tuple(branch.outputs[branch.inputs.index(i)] for i in joint.clock.probe_pair())
                expected = native_part(joint.clock.rotation(pair, -1, False))@expected
            expected = joint._matrix(reads[coin, label, current], False)@expected
        np.testing.assert_allclose(frame.matrix, expected, atol=2e-15, rtol=0)


def test_exact_law_matches_independent_executed_native_instrument(law_frames):
    law, frames = law_frames
    rho = np.diag(list(map(float, joint.clock.transport.reference()))).astype(complex)
    expected = joint.exact_clock_rows(law)
    for frame in frames:
        if frame.abort:
            continue
        value = expected.get((*frame.key[1:], frame.ticks), joint.ZERO)
        lo, hi = joint.evaluate_interval(value)
        actual = np.trace(frame.matrix@rho@frame.matrix.conj().T).real
        assert float(lo)-1e-13 <= actual <= float(hi)+1e-13
    probabilities = joint.factor_probabilities(frames, rho)
    assert abs(sum(probabilities.values())-1) < 2e-14


def test_joint_entangled_field_source_and_offdiagonal_untouched_reference(law_frames):
    _, frames = law_frames
    vertices = (tuple([1]*66), joint.FIELD)
    columns = np.zeros((66, 2, 2), complex)
    a, b = joint.clock.probe_pair()
    columns[a, :, 0] = (1/math.sqrt(264), 1j/math.sqrt(264))
    columns[b, :, 1] = (1/math.sqrt(1320), 1/math.sqrt(1320))
    identity = joint.Frame(('identity',), 0, True, np.eye(66, dtype=complex))
    expected = joint.schmidt_reference(identity, columns, vertices)
    assert abs(expected[0, 1]) > .1
    assert abs(np.trace(expected)-1) < 1e-14
    total = np.zeros((2, 2), complex)
    for frame in frames:
        signs = (None,) if frame.abort else (1, -1)
        for sign in signs:
            output = joint.schmidt_reference(frame, columns, vertices, sign)
            assert np.linalg.eigvalsh(output).min() > -1e-14
            total += output
    np.testing.assert_allclose(total, expected, atol=2e-14, rtol=0)


def test_nonzero_export_abort_keeps_untouched_reference_and_field(monkeypatch):
    original = joint.clock.rotation
    def leaky(pair, sign, raw):
        target = original(pair, sign, False).copy()
        if sign == 1:
            a = joint.clock.full.support()[0][pair[0]]
            outside = next(i for i in range(512) if i not in joint.clock.full.support()[0])
            plane = np.eye(512, dtype=complex)
            plane[np.ix_((a, outside), (a, outside))] = ((.8, .6), (-.6, .8))
            target = plane@target
        return target
    monkeypatch.setattr(joint.clock, 'rotation', leaky)
    frames = joint.clock_frames('C', raw=False)
    vertices = (tuple([1]*66),)
    columns = np.zeros((66, 2, 1), complex)
    columns[joint.clock.probe_pair()[0], :, 0] = (1/math.sqrt(132), 1j/math.sqrt(132))
    abort = next(f for f in frames if f.key == ('abort', 'probe'))
    output = joint.schmidt_reference(abort, columns, vertices)
    assert np.trace(output).real > .17
    assert abs(output[0, 1]) > .08
    total = sum((joint.schmidt_reference(f, columns, vertices, sign)
                 for f in frames for sign in ((None,) if f.abort else (1, -1))), np.zeros((2, 2), complex))
    np.testing.assert_allclose(total, [[.5, -.5j], [.5j, .5]], atol=2e-14, rtol=0)


@pytest.mark.parametrize('value', (True, 0, 3, F(1)))
def test_reject_undeclared_incident_count(value):
    with pytest.raises(ValueError):
        joint.field_probabilities(value)
    with pytest.raises(ValueError):
        joint.field_probabilities(value, -1)


def test_reject_hidden_event_policy_or_nonpositive_density(law_frames):
    _, frames = law_frames
    with pytest.raises(ValueError):
        joint.tick_count('C', (((0, 5, 0, 0)), ((0, 7, 0, 0))))
    bad = np.eye(66, dtype=complex)
    bad[0, 0] = -1
    with pytest.raises(ValueError):
        joint.factor_probabilities(frames, bad)


def test_cosine_enclosure_independent_complex_exponential():
    for n, (lo, hi) in enumerate(joint.cosine_intervals(), 1):
        value = math.cos(4*n*math.log(8/5))
        assert float(lo)-2e-15 <= value <= float(hi)+2e-15
