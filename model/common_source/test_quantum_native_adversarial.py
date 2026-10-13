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
