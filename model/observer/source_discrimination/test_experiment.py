"""Independent complete-channel and frozen-protocol checks of source execution."""
from collections import defaultdict
from fractions import Fraction
import math
import random

import numpy as np
import pytest

from model.predictions import observer_instrument as exact
from model.predictions import protocol
from . import circuit, experiment


DIMENSION = 1 << circuit.SOURCE_QUBITS
TOLERANCE = 2e-10


@pytest.fixture(scope="module")
def executed_events():
    # Reuse completed executions between assertions; expected distributions
    # and density operators below come from the separate exact implementation.
    return {(candidate, setting): experiment.events(candidate, setting)
            for candidate in exact.CANDIDATES for setting in exact.SETTINGS}


def _label(candidate, event):
    return (event["coin"], event["current" if candidate == "C" else "contrast"],
            event["outcome"])


def _exact_density(event):
    density = np.zeros((DIMENSION, DIMENSION), complex)
    for entry in event["source_density"]:
        row = circuit.encode(exact.STATES[entry["row"]])
        column = circuit.encode(exact.STATES[entry["column"]])
        density[row, column] = complex(float(Fraction(entry["real"])),
                                      float(Fraction(entry["imaginary"])))
    return density


def _densities(candidate, events):
    output = defaultdict(lambda: np.zeros((DIMENSION, DIMENSION), complex))
    for event in events:
        if not event["abort"]:
            vector = event["amplitude"]
            output[_label(candidate, event)] += np.outer(vector, vector.conj())
    return dict(output)


@pytest.mark.parametrize("candidate", exact.CANDIDATES)
@pytest.mark.parametrize("setting", exact.SETTINGS)
def test_every_record_and_terminal_density_matches_exact_native_instrument(
        executed_events, candidate, setting):
    events = executed_events[candidate, setting]
    randomizers = (0, 1) if setting == "diagonal_control" else (None,)
    record_codes = range(1 << circuit.RECORD_BITS[candidate])
    # Empty branches and unused record codes remain explicit. This also
    # detects dropped minus outcomes or a success-only denominator.
    expected_keys = {(r, coin, record, outcome) for r in randomizers
                     for coin in (0, 1) for record in record_codes
                     for outcome in (-1, 1)}
    actual_keys = [(event["preparation_record"], event["coin"],
                    event["record_code"], event["outcome"]) for event in events]
    assert len(actual_keys) == len(set(actual_keys))
    assert set(actual_keys) == expected_keys
    total = 0.
    abort = 0.
    for event in events:
        vector = event["amplitude"]
        assert vector.shape == (DIMENSION,)
        assert np.all(np.isfinite(vector))
        probability = float(np.vdot(vector, vector).real)
        total += probability
        if event["abort"]:
            abort += probability
        else:
            assert (event["coin"], event["current"], event["outcome"]) in protocol.BINS
    assert abs(total - 1) < TOLERANCE
    assert any(event["abort"] for event in events)
    assert abort < TOLERANCE
    expected = {}
    for event in exact.event_outputs(candidate, setting):
        record = event["record"]
        label = (record["coin"], record["current" if candidate == "C" else "contrast"],
                 event["outcome"])
        expected[label] = _exact_density(event)
    actual = _densities(candidate, events)
    assert set(actual) == set(expected)
    for label, density in expected.items():
        # Frobenius norm checks every complex matrix entry, including source
        # coherences and padding; no probability or branch normalization.
        assert np.linalg.norm(actual[label] - density) < TOLERANCE, label


@pytest.mark.parametrize("candidate", exact.CANDIDATES)
def test_repair_preserves_all_spectator_and_reference_coherences(candidate):
    rng = np.random.default_rng(88271)
    incoming = np.zeros((DIMENSION, 3), complex)
    coefficients = rng.normal(size=(78, 3)) + 1j * rng.normal(size=(78, 3))
    coefficients /= np.linalg.norm(coefficients)
    for state, value in zip(exact.STATES, coefficients, strict=True):
        incoming[circuit.encode(state)] = value
    expected = defaultdict(lambda: np.zeros_like(incoming))
    for state in exact.STATES:
        total = state[0] + state[1]
        low, high = total // 2, (total + 1) // 2
        for coin in (0, 1):
            endpoints = (low, high) if coin == 0 else (high, low)
            repaired = endpoints + state[2:]
            value = (state[0] - endpoints[0] if candidate == "C"
                     else state[0] - state[1])
            expected[coin, value][circuit.encode(repaired)] += (
                incoming[circuit.encode(state)] / math.sqrt(2))
    actual = experiment.repair_branches(candidate, incoming)
    observed = set()
    for branch in actual:
        amplitude = branch["amplitude"]
        assert amplitude.shape == incoming.shape
        if branch["abort"]:
            assert np.linalg.norm(amplitude) < TOLERANCE
            continue
        label = (branch["coin"], branch["current" if candidate == "C" else "contrast"])
        assert label not in observed
        observed.add(label)
        # Phase-sensitive vector equality is stronger than reduced density
        # equality and protects cross-shell spectator/reference interference.
        assert np.linalg.norm(amplitude - expected[label]) < TOLERANCE
    assert observed == set(expected)
    assert abs(sum(np.linalg.norm(branch["amplitude"])**2 for branch in actual) - 1) < TOLERANCE


def test_readout_wrapper_preserves_coherent_minus_subspace_and_reference():
    incoming = np.zeros((DIMENSION, 2), complex)
    x, z = circuit.encode(exact.STATES[exact.X]), circuit.encode(exact.STATES[exact.Z])
    outside = [code for code in circuit.SOURCE_CODES if code not in (x, z)][:2]
    # One reference column contains the phase-sensitive two-state witness;
    # the other contains a coherent state entirely in the minus subspace.
    incoming[x, 0], incoming[z, 0] = 1 / 2, 1j / 2
    incoming[outside, 1] = (1 / 2, -1j / 2)
    plus = np.zeros_like(incoming)
    plus[x] = plus[z] = (incoming[x] + incoming[z]) / 2
    branches = experiment.readout_branches(incoming)
    assert set(branches) == {-1, 1}
    assert np.linalg.norm(branches[1] - plus) < TOLERANCE
    assert np.linalg.norm(branches[-1] - (incoming - plus)) < TOLERANCE
    minus_density = branches[-1] @ branches[-1].conj().T
    a, b = outside
    assert abs(minus_density[a, b] - 1j / 4) < TOLERANCE
    # Measuring the full source basis would instead remove this entry.
    overmeasured = np.diag(np.diag(minus_density))
    assert np.linalg.norm(minus_density - overmeasured) > .3
    assert abs(sum(np.linalg.norm(vector)**2 for vector in branches.values()) - 1) < TOLERANCE


@pytest.mark.parametrize("candidate", exact.CANDIDATES)
def test_diagonal_preparation_retains_both_phase_randomizer_outcomes(executed_events, candidate):
    diagonal = executed_events[candidate, "diagonal_control"]
    for preparation_record, setting in ((0, "phase_0"), (1, "phase_pi")):
        retained = [event for event in diagonal
                    if event["preparation_record"] == preparation_record]
        actual = _densities(candidate, retained)
        expected = _densities(candidate, executed_events[candidate, setting])
        assert set(actual) == set(expected)
        for label in actual:
            assert np.linalg.norm(actual[label] - expected[label] / 2) < TOLERANCE
        assert abs(sum(np.linalg.norm(event["amplitude"])**2 for event in retained) - .5) < TOLERANCE
    plus = {setting: sum(np.linalg.norm(event["amplitude"])**2
                        for event in executed_events[candidate, setting]
                        if event["outcome"] == 1 and not event["abort"])
            for setting in exact.SETTINGS}
    assert abs(plus["phase_pi_over_2"] - plus["diagonal_control"]) < TOLERANCE
    assert abs(plus["phase_0"] - plus["phase_pi"] - (.5 if candidate == "C" else 0.)) < TOLERANCE


@pytest.mark.parametrize("candidate", exact.CANDIDATES)
def test_frozen_protocol_discriminates_counts_sampled_from_executed_source(executed_events, candidate):
    plan = protocol.checked_protocol()
    original_digest = protocol.digest(plan)
    counts = protocol.data_template(plan)
    counts["data_kind"] = "synthetic"
    counts["provenance"] = "Independent seeded draws from executed native source-pulse branch norms."
    rng = random.Random(92301 if candidate == "C" else 92302)
    for row in counts["settings"]:
        events = executed_events[candidate, row["setting"]]
        cumulative = []
        total = 0.
        for event in events:
            total += float(np.vdot(event["amplitude"], event["amplitude"]).real)
            cumulative.append(total)
        assert abs(total - 1) < TOLERANCE
        measured = dict.fromkeys(protocol.BINS, 0)
        aborts = 0
        for _ in range(protocol.TRIALS):
            draw = rng.random()
            # Do not call synthetic_counts, use the ideal probabilities, or
            # renormalize completed outcomes. Source numerical mass is checked
            # above; a draw in any missing tail would fail this test explicitly.
            selected = next((event for event, boundary in zip(events, cumulative, strict=True)
                             if draw < boundary), None)
            assert selected is not None
            if selected["abort"]:
                aborts += 1
            else:
                measured[selected["coin"], selected["current"], selected["outcome"]] += 1
        row["aborts"] = aborts
        for event in row["events"]:
            event["count"] = measured[event["coin"], event["current"], event["outcome"]]
    verdict = protocol.evaluate(counts, plan)
    other = "E" if candidate == "C" else "C"
    assert verdict["status"] == "comparison_complete"
    assert verdict["candidate_verdicts"][candidate]["verdict"] == "not_rejected"
    assert verdict["candidate_verdicts"][other]["verdict"] == "rejected_under_contract"
    assert verdict["data_kind"] == "synthetic"
    assert verdict["empirical_confirmation"] is False
    assert verdict["provenance_verified"] is False
    assert protocol.digest(plan) == original_digest


@pytest.mark.parametrize("candidate", exact.CANDIDATES)
@pytest.mark.parametrize("setting", exact.SETTINGS)
def test_resource_bound_counts_the_complete_attempted_trial(candidate, setting):
    # A classical branch executes one repair, never the sum of both coin
    # branches. Its maximum must still include both separate recorded coin
    # preparations when the diagonal control uses a retained phase randomizer.
    preparation_settings = ("phase_0", "phase_pi") if setting == "diagonal_control" else (setting,)
    expected = (1 + max(len(circuit.preparation_word(name)) for name in preparation_settings)
                + max(len(circuit.repair_word(candidate, coin)) for coin in (0, 1))
                + len(circuit.readout_word()) + int(setting == "diagonal_control"))
    resources = experiment.resources(candidate, setting)
    assert resources["native_pulses"] == expected
    assert resources["proper_code_load_events"] == expected
    assert resources["proper_code_export_events"] == expected
    assert resources["proper_code_service_events"] == 2 * expected
    assert resources["service_failure_outcomes_retained"] is True
    assert resources["physical_loading_and_service_errors_certified"] is False
    # Source, coin, repair record and pointer each need a blank register;
    # diagonal_control additionally initializes its independent randomizer.
    randomized = int(setting == "diagonal_control")
    record_bits = circuit.RECORD_BITS[candidate]
    assert resources["blank_resets"] == 4 + randomized
    assert resources["record_readout_channels"] == 3 + randomized
    assert resources["retained_classical_bits"] == 2 + record_bits + randomized
    assert resources["maximum_live_code_qubits"] == max(circuit.SOURCE_QUBITS + record_bits,
                                                        circuit.SOURCE_QUBITS + 1)
    assert resources["maximum_work_qubits"] == resources["maximum_live_code_qubits"]
    assert resources["processor_dimension"] == 6
    # Equal pulse counts do not imply equal integrated control: the phase-pi
    # branch contributes two additional nonzero rotations in the mixture.
    norm = lambda word: math.fsum(abs(event["angle"]) for event in word)
    maximum_integrated_norm = (
        math.pi / 4 + max(norm(circuit.preparation_word(name)) for name in preparation_settings)
        + max(norm(circuit.repair_word(candidate, coin)) for coin in (0, 1))
        + norm(circuit.readout_word()) + randomized * math.pi / 4)
    assert abs(resources["integrated_processor_norm_upper"] - maximum_integrated_norm) < 2e-9
    # The receipt names the phase-pi branch for a diagonal trial. Its zero
    # count is branch-specific, while every zero-angle request still incurs
    # one pulse plus load/export. Both fair-coin preparations have nonzero angle.
    representative = "phase_pi" if randomized else setting
    resource_word = (circuit.preparation_word(representative)
                     + circuit.repair_word(candidate, 0) + circuit.readout_word())
    assert resources["zero_angle_requests"] == sum(event["angle"] == 0 for event in resource_word)
