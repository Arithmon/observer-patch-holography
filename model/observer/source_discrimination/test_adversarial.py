"""Failure-retention and phase-sensitive guards beyond the fixed counts."""
import math

import numpy as np
import pytest

from model.observer.source_scalar import compiler as source
from . import circuit, experiment


@pytest.mark.parametrize("candidate", ("C", "E"))
def test_diagonal_resource_bound_covers_the_nonzero_phase_branch(candidate):
    # Both randomizer alternatives have equal pulse counts, but phase_pi
    # has two nonzero phase requests. A branch maximum must charge these.
    diagonal = experiment.resources(candidate, "diagonal_control")
    phase_pi = experiment.resources(candidate, "phase_pi")
    required = phase_pi["integrated_processor_norm_upper"] + math.pi / 4
    assert diagonal["integrated_processor_norm_upper"] >= required - 1e-9


def test_terminal_density_guard_rejects_a_probability_invisible_phase_error():
    original = experiment.events("C", "phase_0")
    phase_word = [source.pulse("I", math.pi / 4, [0]),
                  source.pulse("Z", -math.pi / 4, [0])]
    changed = [{**event, "amplitude": source.execute(phase_word, 9, event["amplitude"])}
               for event in original]
    # A final unitary preserves every record/outcome probability, so even
    # complete classical counts cannot detect this wrong terminal channel.
    for before, after in zip(original, changed, strict=True):
        assert abs(experiment.probability(before) - experiment.probability(after)) < 1e-14
    with pytest.raises(ArithmeticError, match="branch densities"):
        experiment._check_events("C", "phase_0", changed)


def test_terminal_guard_rejects_nonfinite_amplitudes():
    changed = experiment.events("C", "phase_0")
    changed[0]["amplitude"][0] = complex(float("nan"), 0.)
    with pytest.raises(ArithmeticError, match="finite"):
        experiment._check_events("C", "phase_0", changed)


def test_terminal_channel_guard_accepts_an_unobservable_global_phase():
    original = experiment.events("C", "phase_0")
    phase_word = [source.pulse("I", math.pi / 4, [0])]
    changed = [{**event, "amplitude": source.execute(phase_word, 9, event["amplitude"])}
               for event in original]
    result = experiment._check_events("C", "phase_0", changed)
    assert result["maximum_branch_density_frobenius_error"] < experiment.TOLERANCE


@pytest.mark.parametrize("candidate", ("C", "E"))
def test_record_adapter_retains_a_unit_probability_reserved_record(monkeypatch, candidate):
    # Inject a local unitary reaching a reserved code to test the branch
    # adapter itself. This is not a physical-error model or calibration.
    bits = 3 + circuit.RECORD_BITS[candidate]
    reserved = (1 << circuit.RECORD_BITS[candidate]) - 1
    unitary = np.eye(1 << bits, dtype=complex)
    unitary[[0, 8 * reserved]] = unitary[[8 * reserved, 0]]
    monkeypatch.setattr(experiment, "_repair", lambda _candidate, _coin: unitary)
    incoming = np.zeros(512, complex)
    incoming[0] = 1
    branches = experiment.repair_branches(candidate, incoming)
    assert len(branches) == 2 * (1 << circuit.RECORD_BITS[candidate])
    abort = sum(experiment.probability(branch) for branch in branches if branch["abort"])
    success = sum(experiment.probability(branch) for branch in branches if not branch["abort"])
    assert abs(abort - 1) < 1e-14
    assert success == 0
