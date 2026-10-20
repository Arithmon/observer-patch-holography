"""Independent controls for the prepared state's continuing native clock."""
from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest

from ..observer import quantum_mixing
from . import prepared_clock as clock


def _endpoint(state, seam, coin):
    """Enumerate the classical endpoint operation without the native executor."""
    a, b = seam
    total = state[a]+state[b]
    result = list(state)
    result[a] = total//2 + coin*(total % 2)
    result[b] = total-result[a]
    return tuple(result), state[a]-result[a]


def _independent_tapes():
    groups = defaultdict(list)
    seams = quantum_mixing.seams()
    for edge1, seam1 in enumerate(seams):
        for coin1 in (0, 1):
            middle, current1 = _endpoint(clock.INITIAL, seam1, coin1)
            first = edge1, coin1, current1, int(10 in seam1)
            # Derive the controller solely from the publicly committed current.
            change10 = current1*(int(seam1[1] == 10)-int(seam1[0] == 10))
            feedback = int(change10 < 0)
            for edge2, seam2 in enumerate(seams):
                for coin2 in (0, 1):
                    final, current2 = _endpoint(middle, seam2, coin2)
                    second = edge2, coin2, current2, int(10 in seam2)
                    groups[first[3], feedback, second[3]].append((first, second, final))
    return tuple((key, tuple(rows)) for key, rows in sorted(groups.items()))


def _laplacian():
    result = np.zeros((12, 12), dtype=int)
    for a, b in quantum_mixing.seams():
        result[a, a] += 1
        result[b, b] += 1
        result[a, b] -= 1
        result[b, a] -= 1
    return result


def _mass(vector):
    return float(np.vdot(vector, vector).real)


def test_native_tapes_reconstruct_both_repairs_on_the_same_carrier():
    expected = _independent_tapes()
    assert clock.native_histories() == expected
    assert sum(len(rows) for _, rows in expected) == 3600
    assert len({(a[:2], b[:2]) for _, rows in expected for a, b, _ in rows}) == 3600
    reset_disagrees = 0
    for _, rows in expected:
        for first, second, final in rows:
            assert sum(final) == 2 and set(final) <= {0, 1}
            reset, current = _endpoint(clock.INITIAL,
                                       quantum_mixing.seams()[second[0]], second[1])
            reset_disagrees += (reset, current) != (final, second[2])
    assert reset_disagrees > 0, "resetting before the second repair must be detected"


def test_single_port_clocks_and_feedback_have_independent_exact_census():
    counts = Counter()
    for (tick1, feedback, tick2), rows in clock.native_histories():
        counts[tick1, feedback, tick2] += len(rows)
    assert set(counts) == {(a, b, c) for a, b in ((0, 0), (1, 0), (1, 1))
                           for c in (0, 1)}
    assert sum(F(n, 3600)*a for (a, _, _), n in counts.items()) == F(1, 6)
    assert sum(F(n, 3600)*c for (_, _, c), n in counts.items()) == F(1, 6)
    assert sum(F(n, 3600)*b for (_, b, _), n in counts.items()) == F(1, 15)
    assert sum(F(n, 3600)*a*c for (a, _, c), n in counts.items()) == F(1, 36)
    assert sum(F(n, 3600)*b*c for (_, b, c), n in counts.items()) == F(1, 90)
    for coin in (0, 1):
        assert _endpoint(clock.INITIAL, (0, 1), coin) == (clock.INITIAL, 0)


@pytest.fixture(scope="module")
def prepared_pair():
    return clock.prepare(), clock.prepare(ideal=True)


@pytest.fixture(scope="module")
def histories():
    return (clock.execute(axis="Y", intervention=1),
            clock.execute(axis="Y", intervention=1, ideal=True))


def test_ideal_preparation_keeps_every_unnormalized_first_failure(prepared_pair):
    _, prepared = prepared_pair
    lazy = np.eye(12)-_laplacian()/120
    previous = np.eye(12)/math.sqrt(12)
    for n, actual in enumerate(prepared.successes):
        assert np.allclose(actual[:12], previous, rtol=0, atol=2e-14)
        assert np.count_nonzero(actual[12:]) == 0
        if n < len(prepared.failures):
            following = lazy@previous
            assert _mass(prepared.failures[n]) == pytest.approx(
                _mass(previous)-_mass(following), abs=4e-14)
            assert prepared.failures[n].shape == (32756, 12)
            previous = following
    complete = float(prepared.rejection)
    for n, weight in enumerate(prepared.weights):
        complete += float(weight)*(
            sum(_mass(a) for a in prepared.failures[:n])+_mass(prepared.successes[n]))
    assert complete == pytest.approx(1, abs=6e-14)
    assert prepared.rejection > F(1, 10)
    assert len(prepared.weights) == 5 and len(prepared.failures) == 4


def test_requested_mixture_uses_success_weights_without_component_normalization(prepared_pair):
    _, prepared = prepared_pair
    summed = sum(float(w)*(a@a.conj().T)
                 for w, a in zip(prepared.weights, prepared.successes, strict=True))
    normalized = summed/np.trace(summed)
    wrong = sum(float(w)*(a@a.conj().T)/_mass(a)
                for w, a in zip(prepared.weights, prepared.successes, strict=True))
    wrong /= float(sum(prepared.weights))
    assert np.linalg.norm(normalized-wrong, "fro") > 1e-4
    exact_success = F(994135427220848655906122603138321,
                      1239337859764309689882378240000000)
    assert float(np.trace(summed).real) == pytest.approx(float(exact_success), abs=5e-14)


def test_raw_preparation_and_complete_history_stay_close_without_erasing_failures(histories):
    actual, ideal = histories
    assert clock.complete_distance(actual, ideal) < 1e-6
    observed = clock.probabilities(actual)
    expected = clock.probabilities(ideal)
    assert observed.keys() == expected.keys()
    assert sum(observed.values()) == pytest.approx(1, abs=1e-6)
    assert sum(expected.values()) == pytest.approx(1, abs=1e-12)
    assert len([key for key in observed if key[0] == "herald_failure"]) == 10
    assert len([key for key in observed if key[0] == "read"]) == 5*6*16
    assert observed["sampler_reject",] > .17
    assert sum(v for k, v in observed.items() if k[0] == "herald_failure") > .02


@pytest.mark.parametrize("axis", ("X", "Y"))
@pytest.mark.parametrize("intervention", (0, 1))
def test_pair_response_uses_the_same_prepared_gap_and_second_incident_tick(axis, intervention):
    history = clock.execute(axis=axis, intervention=intervention, ideal=True)
    weights = history.preparation.weights
    m0 = float(sum(weights))
    m6 = sum(float(w)*(19/20)**(2*n) for n, w in enumerate(weights))
    gap = math.log(m0/m6)
    for (n, tick1, feedback, tick2), amplitude in history.successful.items():
        # Address2 is one chosen lambda=6 vector, rather than all five modes.
        contrast = _mass(amplitude[0])-_mass(amplitude[2])
        asymmetry = (1-(19/20)**(2*n))/12
        wanted = 0.
        if intervention and feedback:
            wanted = asymmetry*(math.sin(gap*tick2) if axis == "X"
                                else math.cos(gap*tick2))
        assert contrast == pytest.approx(wanted, abs=4e-13)
    if axis == "X" and intervention:
        probabilities = clock.probabilities(history)
        contrast = sum(value*(1 if key[-1] == 0 else -1)
                       for key, value in probabilities.items()
                       if key[0] == "read" and key[-1] in (0, 2))
        assert contrast > 1e-7


def test_reference_phase_corruption_changes_quantum_output_without_any_probability_change(histories):
    actual, ideal = histories
    phases = np.array([1, -1]*6, dtype=complex)
    changed = replace(actual, successful={key: value*phases
                                         for key, value in actual.successful.items()})
    assert clock.probabilities(changed) == clock.probabilities(actual)
    assert clock.complete_distance(changed, actual) > 0.1
    with pytest.raises(ArithmeticError, match="quantum outputs"):
        clock.audit(changed, ideal)


def test_failure_reference_phase_and_missing_failure_are_detected(histories):
    actual, ideal = histories
    phases = np.array([1, -1]*6, dtype=complex)
    changed_preparation = replace(actual.preparation, failures=tuple(
        value*phases for value in actual.preparation.failures))
    changed = replace(actual, preparation=changed_preparation)
    assert clock.probabilities(changed) == clock.probabilities(actual)
    assert clock.complete_distance(changed, actual) > 1e-4
    with pytest.raises(ArithmeticError, match="quantum outputs"):
        clock.audit(changed, ideal)
    omitted = replace(actual, preparation=replace(actual.preparation,
                                                  failures=actual.preparation.failures[:-1]))
    with pytest.raises(ValueError, match="preparation"):
        clock.probabilities(omitted)


def test_global_phase_of_each_complete_quantum_block_is_unobservable(histories):
    actual, _ = histories
    changed = replace(actual, successful={key: -value for key, value in actual.successful.items()},
                      preparation=replace(actual.preparation,
                          failures=tuple(-value for value in actual.preparation.failures)))
    assert clock.complete_distance(changed, actual) < 1e-15
    assert clock.probabilities(changed) == clock.probabilities(actual)


def test_complete_comparison_retains_raw_mass_defect_instead_of_normalizing(histories):
    actual, _ = histories
    changed = replace(actual, successful={key: 1.1*value
                                         for key, value in actual.successful.items()})
    assert sum(clock.probabilities(changed).values()) > 1.1
    assert clock.complete_distance(changed, actual) > .05


@pytest.mark.parametrize("mutation", ("missing", "nonfinite", "real_dtype", "boolean_key"))
def test_complete_output_validator_rejects_missing_or_aliased_fields(histories, mutation):
    actual, _ = histories
    outputs = dict(actual.successful)
    key = next(iter(outputs))
    if mutation == "missing":
        del outputs[key]
    elif mutation == "nonfinite":
        outputs[key] = outputs[key].copy()
        outputs[key][0, 0] = np.nan
    elif mutation == "real_dtype":
        outputs[key] = outputs[key].real.copy()
    else:
        value = outputs.pop(key)
        outputs[(bool(key[0]), *key[1:])] = value
    with pytest.raises(ValueError):
        clock.probabilities(replace(actual, successful=outputs))


def test_pure_block_upper_bound_covers_independent_rank_two_trace_norm():
    random = np.random.default_rng(10102026)
    for _ in range(12):
        a = random.normal(size=12)+1j*random.normal(size=12)
        b = random.normal(size=12)+1j*random.normal(size=12)
        exact = math.sqrt(max(0., (_mass(a)+_mass(b))**2
                             - 4*abs(np.vdot(a, b))**2))/2
        assert exact <= clock._pure_distance(a, b)+1e-12
    a = np.array([1, 1j], dtype=complex)
    assert clock._pure_distance(a, -a) == 0


@pytest.fixture(scope="module")
def saved():
    return clock.checkpoint(intervention=1)


def test_checkpoint_resumes_the_serialized_live_state_with_completed_failures(histories, saved):
    actual, _ = histories
    resumed = clock.execute(axis="Y", intervention=1, saved=saved)
    assert resumed.resumed and not actual.resumed
    assert clock.complete_distance(resumed, actual) < 1e-15
    assert clock.probabilities(resumed) == clock.probabilities(actual)
    for key in actual.successful:
        assert np.array_equal(resumed.successful[key], actual.successful[key])
    assert all(np.array_equal(a, b) for a, b in
               zip(resumed.preparation.failures, actual.preparation.failures, strict=True))


@pytest.mark.parametrize("mutation", ("source", "control", "amplitude", "failure_log", "missing"))
def test_checkpoint_rejects_changed_source_controls_amplitudes_or_completed_failures(saved, mutation):
    payload = json.loads(saved)
    if mutation == "source":
        first = next(iter(payload["source_sha256"]))
        payload["source_sha256"][first] = "0"*64
    elif mutation == "control":
        payload["intervention"] = 0
    elif mutation == "amplitude":
        entry = payload["live_arrays_complex128_little_endian"][0]
        raw = bytearray.fromhex(entry[1])
        raw[0] ^= 1
        entry[1] = raw.hex()
    elif mutation == "failure_log":
        payload["completed_preparation_sha256"] = "0"*64
    else:
        payload["live_arrays_complex128_little_endian"].pop()
    corrupted = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    with pytest.raises(ValueError, match="checkpoint"):
        clock.execute(saved=corrupted)


@pytest.mark.parametrize("axis,intervention", (("Z", 1), ("Y", True), ("Y", 2), (False, 1)))
def test_protocol_choices_reject_implicit_or_unimplemented_controls(axis, intervention):
    with pytest.raises(ValueError):
        clock.execute(axis=axis, intervention=intervention)


def test_resumption_marker_cannot_alias_a_boolean(histories):
    actual, _ = histories
    with pytest.raises(ValueError, match="Boolean"):
        clock.probabilities(replace(actual, resumed=1))


def test_actual_producer_roundtrip_and_stored_negative_norm_refusal():
    # Execute the actual producer, not a serialized fixture substituted for it.
    result = clock.run()
    recorded = json.loads(json.dumps(result, allow_nan=False))
    clock.check_receipt(recorded, result)
    recorded["error_budget"]["filter_complete_Frobenius_residual"] = -1e-10
    with pytest.raises(ValueError, match="scientific guard"):
        clock.check_receipt(recorded, result)
    recorded = json.loads(json.dumps(result, allow_nan=False))
    recorded["executions"][0]["complete_quantum_half_trace_upper_diagnostic"] = -1e-10
    with pytest.raises(ValueError, match="scientific guard"):
        clock.check_receipt(recorded, result)


def test_sampler_intervals_cover_all_coin_words_and_preserve_rejections():
    rows = clock.sampler_contract()
    assert rows[0]["first_word"] == 0 and rows[-1]["past_last_word"] == 1 << 64
    assert rows[-1]["request"] is None
    assert sum((r["probability"] for r in rows), F()) == 1
    for index, row in enumerate(rows):
        assert F(row["past_last_word"]-row["first_word"], 1 << 64) == row["probability"]
        if index:
            assert rows[index-1]["past_last_word"] == row["first_word"]
