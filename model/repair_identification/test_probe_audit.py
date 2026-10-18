"""Independent final-audit controls for complete one-event producer output."""
import math

import numpy as np
import pytest

from . import bridge, core, experiment, source
from . import test_identification as oracle


@pytest.fixture(scope='module')
def h_probe():
    return source.execute_probe('C', 'h_re')


def _descriptor():
    return next(probe for probe in core.PROBES if probe.name == 'h_re')


def _probabilities(execution):
    result = np.zeros(7)
    for branch in execution.branches:
        label = (None, None) if branch.abort else (branch.current, branch.record[-1])
        result[bridge.MENU.index(label)] += np.vdot(branch.amplitude, branch.amplitude).real
    return result


def _replace(execution, replacements):
    branches = tuple(branch._replace(amplitude=replacements[i]) if i in replacements else branch
                     for i, branch in enumerate(execution.branches))
    return execution._replace(branches=branches)


def _produce(monkeypatch, execution):
    monkeypatch.setattr(source, 'execute_probe', lambda candidate, probe: execution)
    return experiment._probe('C', _descriptor(), core.controls()['C'])


@pytest.mark.parametrize('candidate', source.CANDIDATES)
def test_new_single_event_oracle_agrees_with_original_native_endpoints(candidate):
    # Reuse the retained independent native.native oracle, not the production
    # Gram-table map or production read helper. Include every zero outcome.
    table = experiment.candidate_tables()[candidate]
    for descriptor in core.PROBES:
        original = oracle._native_read(candidate, descriptor,
                                      oracle._joint_probe_input(descriptor, 1), 1)
        expected = {(i,): oracle._encoded_density(original[label], 1)
                    for i, label in enumerate(bridge.MENU[:-1])}
        expected[6,] = np.zeros((8, 8), complex)
        actual = bridge.exact_probe(table, descriptor.name)
        assert tuple(actual) == tuple(expected)
        for record in expected:
            assert np.array_equal(actual[record], expected[record])


def test_single_event_coarsening_retains_all_seven_quantum_blocks(h_probe):
    blocks, diagnostics = bridge.coarse_probe(h_probe)
    assert tuple(blocks) == tuple((i,) for i in range(7))
    expected = {key: np.zeros((8, 8), complex) for key in blocks}
    full_mass = 0.
    for branch in h_probe.branches:
        index = 6 if branch.abort else bridge.MENU.index((branch.current, branch.record[-1]))
        expected[index,] += np.outer(branch.amplitude[:8], branch.amplitude[:8].conj())
        full_mass += np.vdot(branch.amplitude, branch.amplitude).real
    for key in blocks:
        assert np.array_equal(blocks[key], expected[key])
    assert diagnostics['full_probability_mass'] == pytest.approx(full_mass, abs=1e-15)


def test_producer_rejects_relative_phase_error_with_every_probability_unchanged(h_probe, monkeypatch):
    index = max((i for i, branch in enumerate(h_probe.branches)
                 if np.count_nonzero(np.abs(branch.amplitude) > 1e-4) >= 2),
                key=lambda i: np.linalg.norm(h_probe.branches[i].amplitude))
    original = h_probe.branches[index].amplitude
    changed = original.copy()
    changed[int(np.argmax(np.abs(changed)))] *= -1
    forged = _replace(h_probe, {index: changed})
    assert np.array_equal(_probabilities(forged), _probabilities(h_probe))
    difference = np.outer(original[:8], original[:8].conj())-np.outer(changed[:8], changed[:8].conj())
    assert np.linalg.svd(difference, compute_uv=False).sum()/2 > 0.49
    with pytest.raises(ArithmeticError):
        _produce(monkeypatch, forged)


def test_producer_accepts_an_unobservable_whole_raw_branch_phase(h_probe, monkeypatch):
    index = max(range(len(h_probe.branches)), key=lambda i: np.linalg.norm(h_probe.branches[i].amplitude))
    unchanged = _replace(h_probe, {index: -h_probe.branches[index].amplitude})
    result = _produce(monkeypatch, unchanged)
    assert result['complete_source_to_exact_distance_bound_numeric'] < experiment.TOLERANCE


def test_producer_rejects_mass_moved_to_abort_without_changing_selected_statistic(h_probe, monkeypatch):
    selected = bridge.MENU.index((_descriptor().current, 1))
    donor = max((i for i, branch in enumerate(h_probe.branches)
                 if not branch.abort and (branch.current, branch.record[-1]) != bridge.MENU[selected]),
                key=lambda i: np.linalg.norm(h_probe.branches[i].amplitude))
    failure = next(i for i, branch in enumerate(h_probe.branches) if branch.abort)
    amplitude = h_probe.branches[donor].amplitude/math.sqrt(2)
    forged = _replace(h_probe, {donor: amplitude, failure: amplitude.copy()})
    before, after = _probabilities(h_probe), _probabilities(forged)
    assert after[selected] == before[selected]
    assert abs(sum(before)-sum(after)) < 1e-12
    assert after[6]-before[6] > 0.01
    with pytest.raises(ArithmeticError):
        _produce(monkeypatch, forged)


def test_producer_charges_full_source_coherences_even_when_tail_mass_is_tiny(h_probe, monkeypatch):
    index = max(range(len(h_probe.branches)), key=lambda i: np.linalg.norm(h_probe.branches[i].amplitude))
    changed = np.zeros(512, complex)
    changed[:8] = h_probe.branches[index].amplitude[:8]
    changed[511] = 1e-6
    forged = _replace(h_probe, {index: changed})
    blocks, diagnostics = bridge.coarse_probe(forged)
    original, _ = bridge.coarse_probe(h_probe)
    assert all(np.array_equal(blocks[key], original[key]) for key in blocks)
    assert diagnostics['tail_mass'] < experiment.TOLERANCE
    assert diagnostics['full_output_projection_bound_numeric'] > 100*experiment.TOLERANCE
    with pytest.raises(ArithmeticError):
        _produce(monkeypatch, forged)
