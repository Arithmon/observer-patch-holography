"""Hostile semantic controls for the complete coherent native/field join."""
from dataclasses import replace

import numpy as np
import pytest

from ..observer.source_discrimination import circuit
from . import native_routing as routing
from . import quantum_native as q


@pytest.fixture(scope='module')
def baseline():
    return q.execute('C', 2, 8)


def test_wrong_raw_current_decoder_cannot_redefine_the_independent_field_target(baseline, monkeypatch):
    original = circuit.record_values

    def reversed_current(law, coin):
        values = original(law, coin)
        return {raw: -current for raw, current in values.items()} if law == 'C' else values

    q._native_tree.cache_clear()
    monkeypatch.setattr(circuit, 'record_values', reversed_current)
    try:
        changed = q.execute('C', 2, 8)
        # This mutation preserves all raw record labels and native maps. Its
        # sole change is which physical field preparation the current chooses.
        assert all(a.record == b.record for a, b in zip(baseline.branches, changed.branches, strict=True))
        assert max(np.linalg.norm(a.work_map-b.work_map) for a, b in
                   zip(baseline.branches, changed.branches, strict=True)) > .1
        with pytest.raises(ArithmeticError, match='independent targets'):
            q.audit(changed)
        # Paired verifier control: the former checker consulted the same
        # corrupted decoder as execution and consequently accepted this run.
        def obsolete_decoder(law, coin, raw):
            value = circuit.record_values(law, coin).get(raw)
            if value is None or law == 'C':
                return value
            return (value-(2*coin-1)*(abs(value) % 2))//2

        with monkeypatch.context() as weakened:
            weakened.setattr(q, '_record_current', obsolete_decoder)
            report = q.audit(changed)
            assert report['full_outcome_count_including_abort'] == 770
    finally:
        q._native_tree.cache_clear()


def test_constant_routing_cannot_pass_the_full_channel_audit(baseline, monkeypatch):
    q._native_tree.cache_clear()
    monkeypatch.setattr(routing, 'routing_unitary', lambda *, inverse=False: np.eye(512))
    try:
        changed = q.execute('C', 2, 8)
        with pytest.raises(ArithmeticError, match='independent targets'):
            q.audit(changed)
    finally:
        q._native_tree.cache_clear()


def test_hidden_native_phase_is_rejected_even_when_basis_probabilities_agree(baseline, monkeypatch):
    original = routing.routing_unitary
    phase = np.ones(512, complex)
    phase[circuit.SOURCE_CODES[0]] = 1j
    wrong = original()*phase[None, :]
    q._native_tree.cache_clear()
    monkeypatch.setattr(routing, 'routing_unitary',
                        lambda *, inverse=False: original(inverse=True) if inverse else wrong)
    try:
        changed = q.execute('C', 2, 8)
        # A phase fault cannot be rejected from basis transition probabilities.
        for a, b in zip(baseline.branches, changed.branches, strict=True):
            np.testing.assert_allclose(abs(a.native_map)**2, abs(b.native_map)**2,
                                       atol=2e-10, rtol=0)
        with pytest.raises(ArithmeticError, match='independent targets'):
            q.audit(changed)
    finally:
        q._native_tree.cache_clear()


def test_readout_refuses_nonfinite_aborted_probability_from_finite_forged_factors(baseline, monkeypatch):
    index = next(i for i, b in enumerate(baseline.branches) if not b.record[2])
    branch = baseline.branches[index]
    # Every stored entry is finite, but multiplying independently overflowing
    # factor norms can produce NaN. A total NaN must never pass normalization.
    changed = replace(branch, native_map=np.zeros_like(branch.native_map),
                      work_state=np.full_like(branch.work_state, 1e200))
    branches = list(baseline.branches)
    branches[index] = changed
    forged = replace(baseline, branches=tuple(branches))
    with np.errstate(over='ignore', invalid='ignore'):
        with pytest.raises((ArithmeticError, ValueError)):
            q.readout(forged, q.prepared_native('phase_0'))
        # Removing the finite guards reproduces the former false acceptance;
        # comparison with one alone does not reject NaN probability mass.
        with monkeypatch.context() as weakened:
            weakened.setattr(q.math, 'isfinite', lambda value: True)
            result = q.readout(forged, q.prepared_native('phase_0'))
            assert np.isnan(result['probability_mass'])
            assert np.isnan(result['abort_probability'])
            assert result['complete_outcome_count'] == 295106


@pytest.mark.parametrize('joint', (False, True))
@pytest.mark.parametrize('amplitude', (1e200, 1e200+1e200j))
def test_public_inputs_reject_overflowing_norms_of_finite_amplitudes(baseline, joint, amplitude):
    shape = (78, 16, 2) if joint else (78, 2)
    initial = np.zeros(shape, complex)
    initial.flat[0] = amplitude
    assert np.isfinite(initial).all()
    consumer = q.output_on_joint_input if joint else q.output
    with np.errstate(over='ignore', invalid='ignore'):
        with pytest.raises(ValueError, match='normalized'):
            consumer(baseline, initial)


@pytest.mark.parametrize('joint', (False, True))
def test_nan_norm_reduction_guard_has_a_deterministic_mutation_control(baseline, joint, monkeypatch):
    shape = (78, 16, 2) if joint else (78, 2)
    initial = np.zeros(shape, complex)
    initial.flat[0] = 1.
    consumer = q.output_on_joint_input if joint else q.output
    # BLAS implementations may return either inf or NaN as the real part of
    # an overflowing complex inner product. Force its NaN reduction here so
    # the guard mutation is independent of that implementation detail.
    monkeypatch.setattr(q.np, 'vdot', lambda left, right: complex(float('nan'), 0))
    with pytest.raises(ValueError, match='normalized'):
        consumer(baseline, initial)
    with monkeypatch.context() as weakened:
        weakened.setattr(q.math, 'isfinite', lambda value: True)
        assert consumer(baseline, initial) is not None


@pytest.mark.parametrize('joint', (False, True))
def test_public_inputs_preserve_valid_complex_reference_amplitudes(baseline, joint):
    shape = (78, 16, 2) if joint else (78, 2)
    initial = np.zeros(shape, complex)
    initial.flat[0], initial.flat[-1] = .5+.5j, .5-.5j
    if joint:
        rows = q.output_on_joint_input(baseline, initial)
        mass = 0.
        for branch, (record, value) in zip(baseline.branches, rows, strict=True):
            mass += float(np.vdot(value, value).real)
            expected = np.zeros((512, 16, 2), complex)
            expected[:, :, 0] = (.5+.5j)*np.outer(branch.native_map[:, 0], branch.work_map[:, 0])
            expected[:, :, 1] = (.5-.5j)*np.outer(branch.native_map[:, -1], branch.work_map[:, -1])
            assert record == branch.record
            np.testing.assert_allclose(value, expected, atol=2e-12, rtol=0)
    else:
        rows = q.output(baseline, initial)
        mass = sum(float(np.vdot(a, a).real)*float(np.vdot(v, v).real)
                   for a, v in rows.values())
        for branch in baseline.branches:
            expected = np.column_stack(((.5+.5j)*branch.native_map[:, 0],
                                        (.5-.5j)*branch.native_map[:, -1]))
            np.testing.assert_allclose(rows[branch.record][0], expected, atol=2e-12, rtol=0)
    assert abs(mass-1) < 2e-10
