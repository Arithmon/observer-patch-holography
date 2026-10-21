"""Independent controls for the prepared state's continuing native clock."""
from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm, logm

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


@pytest.fixture(scope="module")
def independent_terminal_outputs():
    # Derive preparation and its clock directly from the native Laplacian;
    # neither producer prefix assembly nor its ideal maps enter this oracle.
    lazy = np.eye(12)-_laplacian()/120
    weights = clock.preparation.subordinated_filter_weights()
    mixture = sum(float(w)*np.linalg.matrix_power(lazy, 2*n)
                  for n, w in enumerate(weights))
    hamiltonian = -logm(mixture/float(sum(weights)))
    tick = np.eye(16, dtype=complex)
    tick[:12, :12] = expm(-1j*hamiltonian)
    # Reuse only the declared eigenvectors that define the chosen pair and
    # address readout. The clock above does not use this spectral basis.
    basis = clock.source.native_spectral_data()["basis"]
    pair = 2  # The declared single lambda=6 address.
    zero, selected = basis[:, 0], basis[:, pair]
    pair_x = np.outer(zero, selected)+np.outer(selected, zero)
    kick = expm(1j*np.pi*pair_x/4)
    detectors = {axis: basis.T.astype(complex) for axis in ("X", "Y")}
    detectors["X"][0] = (zero+selected)/math.sqrt(2)
    detectors["X"][pair] = (-zero+selected)/math.sqrt(2)
    detectors["Y"][0] = (zero-1j*selected)/math.sqrt(2)
    detectors["Y"][pair] = (-1j*zero+selected)/math.sqrt(2)
    controls = [control for control, _ in _independent_tapes()]
    outputs = {}
    for axis, detector in detectors.items():
        for intervention in (0, 1):
            expected = {}
            for n in range(len(weights)):
                prepared = np.zeros((16, 12), dtype=complex)
                prepared[:12] = np.linalg.matrix_power(lazy, n)/math.sqrt(12)
                for first, feedback, second in controls:
                    expected[n, first, feedback, second] = (
                        detector@np.linalg.matrix_power(tick, second)
                        @np.linalg.matrix_power(kick, intervention*feedback)
                        @np.linalg.matrix_power(tick, first)@prepared)
            outputs[axis, intervention] = expected
    return outputs


def _reference_density_rows(amplitude):
    return amplitude[:, :, None]*amplitude[:, None, :].conj()


def _assert_terminal_reference_outputs(actual, expected):
    assert actual.successful.keys() == expected.keys()
    for key, amplitude in actual.successful.items():
        np.testing.assert_allclose(
            _reference_density_rows(amplitude), _reference_density_rows(expected[key]),
            atol=2e-11, rtol=0,
            err_msg=f"terminal reference density for {actual.axis}, {actual.intervention}, {key}")


@pytest.mark.parametrize("axis", ("X", "Y"))
@pytest.mark.parametrize("intervention", (0, 1))
def test_complete_terminal_reference_matches_independent_composition(
        independent_terminal_outputs, axis, intervention):
    actual = clock.execute(axis=axis, intervention=intervention)
    _assert_terminal_reference_outputs(actual, independent_terminal_outputs[axis, intervention])


def test_independent_reference_oracle_detects_omitted_first_clock_despite_equal_probabilities(
        histories, independent_terminal_outputs, monkeypatch):
    actual, _ = histories
    expected = independent_terminal_outputs["Y", 1]
    _assert_terminal_reference_outputs(actual, expected)

    def omit_first_clock(amplitude, control, intervention, maps):
        return (maps[f"kick_plus_{clock.BAND}"]@amplitude
                if control[1] and intervention else amplitude)

    monkeypatch.setattr(clock, "_prefix", omit_first_clock)
    changed = clock.execute(axis="Y", intervention=1)
    actual_law, changed_law = clock.probabilities(actual), clock.probabilities(changed)
    assert actual_law.keys() == changed_law.keys()
    np.testing.assert_allclose(list(actual_law.values()), list(changed_law.values()),
                               atol=2e-11, rtol=0)
    difference = max(np.linalg.norm(
        _reference_density_rows(changed.successful[key])-_reference_density_rows(amplitude))
        for key, amplitude in expected.items())
    assert difference > 1e-4
    with pytest.raises(AssertionError, match="terminal reference density"):
        _assert_terminal_reference_outputs(changed, expected)


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


@pytest.mark.parametrize("request_index", (1, 4))
def test_overflowing_retained_preparation_cannot_bypass_audit(histories, request_index):
    actual, ideal = histories
    successes = list(actual.preparation.successes)
    successes[request_index] = np.full((16, 12), 1e308, dtype=complex)
    changed = replace(actual, preparation=replace(actual.preparation,
                                                  successes=tuple(successes)))
    # The saved intermediate preparation is exposed data but not reused when
    # reading these already-computed terminal outputs. Its audit must stand
    # independently of the unchanged, well-behaved terminal distribution.
    assert clock.probabilities(changed) == clock.probabilities(actual)
    with pytest.raises(ArithmeticError, match="norm calculation is not finite"):
        clock.audit(changed, ideal)


def test_distance_rejects_nonfinite_arithmetic_on_finite_vectors():
    # Each separate squared norm fits; the norm of their difference does not.
    a = np.array([1.2e154, 0], dtype=complex)
    b = np.array([0, 1.2e154], dtype=complex)
    assert np.isfinite(np.linalg.norm(a)) and np.isfinite(np.linalg.norm(b))
    with pytest.raises(ArithmeticError, match="distance calculation is not finite"):
        clock._pure_distance(a, b)


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
    recorded = json.loads(json.dumps(result, allow_nan=False))
    recorded["shared_state_field"]["executions"][0]["complete_quantum_half_trace_upper_diagnostic"] = -1e-10
    with pytest.raises(ValueError, match="scientific guard"):
        clock.check_receipt(recorded, result)


def test_pointer_power_is_recomputed_with_its_propagated_bound_not_fixed_tolerance():
    # Use the actual retained receipt, allowing this regression to remain
    # cheap without replacing any producer execution in the roundtrip test.
    result = json.loads((clock.ROOT/"model/predictions/prepared_clock_results.json").read_text())
    result["shared_state_field"].pop("raw_pointer_estimator_power", None)
    recorded = json.loads(json.dumps(result))
    field = recorded["shared_state_field"]
    old = complex(*field["raw_pointer_estimator"])
    field["raw_pointer_estimator"][0] += 1e-12
    field["raw_pointer_estimator"][1] -= 1e-12
    new = complex(*field["raw_pointer_estimator"])
    assert abs(new**field["angle_copy_parameter"]-old**field["angle_copy_parameter"]) > 2e-10
    clock.check_receipt(recorded, result)
    field["raw_pointer_estimator"][0] += 1e-6
    with pytest.raises(ValueError):
        clock.check_receipt(recorded, result)
    # Matching stored/current malformed metadata must not bypass the domain
    # guard merely because the ordinary structural comparator sees equality.
    with pytest.raises(ValueError, match="propagated scientific guard"):
        clock.check_receipt(recorded, recorded)


def test_sampler_intervals_cover_all_coin_words_and_preserve_rejections():
    rows = clock.sampler_contract()
    assert rows[0]["first_word"] == 0 and rows[-1]["past_last_word"] == 1 << 64
    assert rows[-1]["request"] is None
    assert sum((r["probability"] for r in rows), F()) == 1
    for index, row in enumerate(rows):
        assert F(row["past_last_word"]-row["first_word"], 1 << 64) == row["probability"]
        if index:
            assert rows[index-1]["past_last_word"] == row["first_word"]


def test_resource_counts_include_external_selection_and_scope_storage():
    result = clock.resources()
    success = F(clock.source.prepared_reference()["normalization"])/12
    # Two IID selections among 30 seams cannot be generated exactly from a
    # fixed finite fair-bit word: 30 does not divide any power of two.
    assert result["uniform_thirty_seam_selection_service_calls_maximum"] == 2
    assert result["uniform_thirty_seam_selection_service_calls_ideal_mean"] == 2*success
    assert result["exact_uniform_seam_selector_is_a_supplied_service"]
    assert not result["finite_worst_case_fair_bit_budget_for_seam_selection_claimed"]
    assert result["native_carrier_placement_coin_events_maximum"] == 64+2
    assert result["native_carrier_placement_coin_events_ideal_mean"] == 64+2*success
    thresholds = [r["past_last_word"] for r in clock.sampler_contract()
                  if r["request"] is not None]
    assert len(thresholds) == 5
    assert all(int.from_bytes(t.to_bytes(8, "little"), "little") == t for t in thresholds)
    assert result["sampler_cumulative_threshold_storage_bits"] == 320
    size = result["source_word_packet_utf8_bytes"]
    bits = result["source_word_packet_byte_cursor_bits"]
    assert 1 << (bits-1) < size <= 1 << bits
    assert result["expanded_source_pulse_ordinal_bits"] == 20 < bits
    assert result["observer_setting_storage_bits"] == 2
    assert not result["controller_interpreter_working_memory_included"]


def test_one_prepared_state_fixes_all_eleven_fluctuation_pairs_and_clock_rates():
    data = clock.shared_fluctuation_modes()
    reference = clock.source.prepared_reference()
    rho = data["density"]
    h = np.diag([*reference["modular_gaps"], 0., 0., 0., 0.])
    tick = expm(-1j*.73*h)
    for mode in data["modes"]:
        q, p, g = mode["Q"], mode["P"], mode["gap"]
        assert abs(np.trace(rho@(q@p-p@q))-1j) < 2e-14
        assert np.trace(rho@q) == np.trace(rho@p) == 0
        np.testing.assert_allclose(tick.conj().T@q@tick,
            math.cos(.73*g)*q+math.sin(.73*g)*p, atol=1e-13, rtol=0)
        # A finite matrix commutator is not the identity operator. Its
        # product-copy centered L2 error falls as 1/N, with this coefficient.
        centered = q@p-p@q-1j*np.eye(16)
        variance = np.trace(rho@centered.conj().T@centered).real
        assert variance == pytest.approx(mode["equal_mode_commutator_L2_variance_times_copy_count"])
        assert variance > 1
    np.testing.assert_allclose(data["symplectic"],
        np.kron(np.eye(11), [[0., 1.], [-1., 0.]]), atol=2e-14, rtol=0)
    assert np.count_nonzero(abs(data["covariance"]-np.diag(np.diag(data["covariance"]))) > 1e-14) == 0


@pytest.mark.parametrize("copies", (1, 4, 64, 4096))
def test_ordered_weyl_finite_copy_bound_against_independent_matrix_exponentials(copies):
    data = clock.shared_fluctuation_modes()
    law = clock.ordered_weyl_law(copies=copies)
    ordered = np.eye(16, dtype=complex)
    parity = np.diag([1., *([-1.]*11), *([1.]*4)])
    for matrix in law["matrices"]:
        np.testing.assert_allclose(parity@matrix@parity, -matrix, atol=0, rtol=0)
        ordered = expm(1j*matrix/math.sqrt(copies))@ordered
    single = np.trace(data["density"]@ordered)
    assert abs(single-law["one_copy_characteristic"]) < 3e-16
    assert abs(single**copies-law["gaussian_characteristic"]) <= law["finite_copy_error_upper"]+1e-12
    variances = {row["address"]: row["variance"] for row in data["modes"]}
    independent_c = .5*(.02*variances[4]+.0004*variances[1]+.000625*variances[2])-.005j
    assert abs(independent_c-law["gaussian_exponent"]) < 3e-16
    reversed_law = clock.ordered_weyl_law(tuple(reversed(clock.WEYL_FACTORS)), copies)
    assert abs(reversed_law["gaussian_characteristic"]-law["gaussian_characteristic"].conjugate()) < 1e-15


def test_two_copy_factorization_is_a_tensor_law_not_independent_mode_preparation():
    rho = clock.shared_fluctuation_modes()["density"][:12, :12]
    law = clock.ordered_weyl_law(copies=2)
    # Explicitly construct the two-copy observable; the factors within each
    # copy share the constant vector and do not form independent subsystems.
    unitary = law["one_copy_unitary_spectral"][:12, :12]
    actual = np.trace(np.kron(rho, rho)@np.kron(unitary, unitary))
    assert abs(actual-law["independent_copy_characteristic"]) < 3e-16
    q1 = clock.shared_fluctuation_modes()["modes"][0]["Q"]
    q2 = clock.shared_fluctuation_modes()["modes"][1]["Q"]
    assert np.linalg.norm(q1@q2-q2@q1) > 1


@pytest.fixture(scope="module")
def pointer_histories():
    return {axis: (clock.execute_weyl(axis=axis), clock.execute_weyl(axis=axis, ideal=True))
            for axis in ("X", "Y")}


@pytest.mark.parametrize("axis", ("X", "Y"))
def test_pointer_complete_reference_oracle_and_all_failure_mass(pointer_histories, axis):
    actual, target = pointer_histories[axis]
    basis = clock.source.native_spectral_data()["basis"]
    law = clock.ordered_weyl_law()
    work = basis@law["one_copy_unitary_spectral"]@basis.T
    # Construct the Hadamard-test branches directly, independently of the
    # controlled-Pauli compiler and its target constructor.
    outputs = []
    lazy = np.eye(12)-_laplacian()/120
    for n in range(5):
        amplitude = np.zeros((16, 12), complex)
        amplitude[:12] = np.linalg.matrix_power(lazy, n)/math.sqrt(12)
        if axis == "X":
            outputs.append(np.vstack((amplitude+work@amplitude, -amplitude+work@amplitude))/2)
        else:
            outputs.append(np.vstack((amplitude-1j*work@amplitude, -1j*amplitude+work@amplitude))/2)
    for observed, expected in zip(actual.successful, outputs, strict=True):
        np.testing.assert_allclose(_reference_density_rows(observed),
                                   _reference_density_rows(expected), atol=3e-11, rtol=0)
    probabilities = clock.pointer_probabilities(actual)
    assert len(probabilities) == 1+10+5*32
    assert abs(sum(probabilities.values())-1) < 2e-10
    audit = clock.pointer_audit(actual, target)
    expected = law["one_copy_characteristic"].real if axis == "X" else law["one_copy_characteristic"].imag
    success = float(clock.source.prepared_reference()["normalization"])/12
    assert abs(audit["unconditional_pointer_contrast"]-success*expected) < 2e-10
    # A complete terminal probability law does not certify the references.
    changed = replace(actual, successful=tuple(value*np.array([1, -1]*6)
                                               for value in actual.successful))
    assert clock.pointer_probabilities(changed) == probabilities
    with pytest.raises(ArithmeticError, match="pointer quantum outputs"):
        clock.pointer_audit(changed, target)


@pytest.mark.parametrize("factors,copies", (
    (((True, "X", F(1)),), 1), (((1, "Z", F(1)),), 1),
    (((1, "X", 0.1),), 1), (clock.WEYL_FACTORS, True),
    (clock.WEYL_FACTORS, 0), (clock.WEYL_FACTORS, (1 << 40)+1)))
def test_fluctuation_menu_rejects_implicit_or_unimplemented_choices(factors, copies):
    with pytest.raises(ValueError):
        clock.ordered_weyl_law(factors, copies)


def test_prepared_gap_enclosures_use_rational_positive_log_series():
    import mpmath as mp
    with mp.workdps(90):
        weights = [mp.mpf(w.numerator)/w.denominator
                   for w in clock.preparation.subordinated_filter_weights()]
        eigenvalues = (5-mp.sqrt(5), mp.mpf(6), 5+mp.sqrt(5))
        for lam, (lo, hi) in zip(eigenvalues, clock.prepared_gap_intervals(), strict=True):
            mass = sum(w*(1-lam/120)**(2*n) for n, w in enumerate(weights))
            expected = mp.log(sum(weights)/mass)
            assert mp.mpf(lo.numerator)/lo.denominator <= expected
            assert expected <= mp.mpf(hi.numerator)/hi.denominator
            assert hi-lo < F(1, 10**23)


def test_prepared_optical_comparison_covers_every_fixed_gap_assignment_without_refit():
    result = clock.prepared_optical_comparison()
    assert len(result["assignments"]) == 9
    assert {(row["numerator_band"], row["denominator_band"]) for row in result["assignments"]} == {
        (i, j) for i in (1, 2, 3) for j in (1, 2, 3)}
    assert all(row["all_nine_direct_assignments_incompatible"] for row in result["comparisons"])
    assert min(row["minimum_absolute_ratio_separation_lower"] for row in result["comparisons"]) > F(57, 1000)
    assert result["declared_absolute_per_gap_attachment_error"] > 10*result["numerically_checked_clock_operator_guard"]
    assert not result["prospective_prediction_or_OPH_exclusion"]
    assert not result["standard_uncertainties_or_numerical_guards_are_physical_hard_bounds"]


def test_rank_supported_fourth_derivative_and_physical_field_norm():
    from itertools import product
    data = clock.shared_fluctuation_modes()
    rho = data["density"][:12, :12]
    random = np.random.default_rng(21102026)
    matrices, radii = [], []
    for _ in range(3):
        f, h = random.normal(size=(2, 11))/100
        a, radius = np.zeros((12, 12), complex), 0.
        expected_norm_squared = 0.
        for j, mode in enumerate(data["modes"]):
            omega = mode["gap"]
            a += (f[j]*mode["Q"][:12, :12]/math.sqrt(omega)
                  +h[j]*mode["P"][:12, :12]*math.sqrt(omega))
            term = f[j]**2/omega+h[j]**2*omega
            radius += term
            expected_norm_squared += term/(2*mode["delta"])
        assert np.linalg.norm(a, 2)**2 == pytest.approx(expected_norm_squared)
        matrices.append(a)
        radii.append(math.sqrt(radius))
    span = np.column_stack((np.eye(12)[:, 0], *[a[:, 0] for a in matrices]))
    orthogonal = np.linalg.qr(span)[0]
    projection = orthogonal@orthogonal.conj().T
    for a in matrices:
        np.testing.assert_allclose(projection@a@projection, a, atol=2e-14, rtol=0)
    weight = np.trace(rho@projection).real
    assert weight <= 4*rho[0, 0].real+1e-15
    derivative = np.zeros((12, 12), complex)
    for powers in product(range(5), repeat=3):
        if sum(powers) != 4:
            continue
        coefficient = math.factorial(4)/math.prod(math.factorial(k) for k in powers)
        term = np.eye(12, dtype=complex)
        for a, k in zip(matrices, powers, strict=True):
            term = np.linalg.matrix_power(1j*a, k)@expm(.31j*a)@term
        derivative += coefficient*term
    np.testing.assert_allclose(projection@derivative@projection, derivative, atol=2e-12, rtol=0)
    norm_sum = sum(np.linalg.norm(a, 2) for a in matrices)
    assert abs(np.trace(rho@derivative)) <= weight*norm_sum**4
    copies = 4096
    unitary = np.eye(12, dtype=complex)
    for a in matrices:
        unitary = expm(1j*a/math.sqrt(copies))@unitary
    c = sum(np.trace(rho@a@a)/2 for a in matrices)
    c += sum(np.trace(rho@matrices[j]@matrices[i]) for j in range(3) for i in range(j))
    bound = clock.growing_fluctuation_bound(12, copies,
        tuple(F(math.ceil(r*10**10), 10**10) for r in radii),
        gap_floor=F(1, 100), inverse_temperature=F(1))
    assert abs(np.trace(rho@unitary)**copies-np.exp(-c)) <= float(bound["finite_copy_gaussian_error_upper"])


def test_growing_dimension_requires_copy_growth_and_correct_field_normalization():
    result = clock.growing_field_controls()
    rows = result["finite_family"]
    # |cos sqrt(x)-cos sqrt(y)|<=|x-y|/2 and the log1p
    # remainder give an elementary 2/d bound for d>=16.
    assert abs(rows[-1]["copy_controls"][0]["characteristic"]-result["N_equals_dimension_limit"]) < 2/rows[-1]["dimension"]
    assert abs(result["N_equals_dimension_limit"]-result["gaussian_characteristic"]) > .01
    linear_copy_errors = [row["copy_controls"][0]["gaussian_error"] for row in rows]
    assert linear_copy_errors[-1] > .01
    assert all(b > a > 0 for a, b in zip(linear_copy_errors, linear_copy_errors[1:]))
    errors = [row["copy_controls"][1]["gaussian_error"] for row in rows]
    assert all(a > b for a, b in zip(errors, errors[1:]))
    for row in result["finite_prepared_clock_normalization"]:
        raw = complex(*row["raw_operator_norm_scaled_two_point"])*row["operator_norm_squared"]
        target = complex(*row["ideal_ordered_two_point"])
        assert abs(raw-target) <= row["source_numerical_error_upper"]
        assert row["normalized_retarded_derivative_at_zero"] == 1
        assert row["unscaled_Q_retarded_derivative_at_zero_control"] < .05
    totals = [row["total_ordered_characteristic_error_upper"]
              for row in result["sufficient_full_port_resource_schedule"]]
    assert all(a > b > 0 for a, b in zip(totals, totals[1:]))
    assert not result["finite_clock_frequencies_identified_with_sqrt_native_L"]


def test_conditioning_and_partition_bounds_charge_growing_resources():
    kwargs = dict(gap_floor=F(1, 2), inverse_temperature=F(1),
                  complete_preparation_error=F(1, 10**9), success_lower=F(1, 1000))
    low = clock.growing_fluctuation_bound(1000, 1000, (F(1),), **kwargs)
    assert low["conditional_state_trace_norm_error_upper"] == F(1, 250000)
    assert low["N_copy_raw_source_error_upper"] == F(1, 249)
    assert low["expected_preparation_attempts_upper"] == 1000000
    traced = clock.growing_fluctuation_bound(1000, 1000, (F(1),), partition_upper=F(2), **kwargs)
    assert traced["finite_copy_gaussian_error_upper"] < low["finite_copy_gaussian_error_upper"]
    assert traced["N_copy_raw_source_error_upper"] == low["N_copy_raw_source_error_upper"]
    exhausted = clock.growing_fluctuation_bound(1000, 10**6, (F(1),), **kwargs)
    assert exhausted["N_copy_raw_source_error_upper"] is None
    assert exhausted["total_ordered_characteristic_error_upper"] is None
    no_success = clock.growing_fluctuation_bound(12, 10, (F(1),), gap_floor=F(1),
        inverse_temperature=F(1), complete_preparation_error=F(1, 4), success_lower=F(1, 2))
    assert no_success["conditional_state_trace_norm_error_upper"] is None


def test_gap_and_time_error_bound_controls_the_whole_ordered_exponent():
    # Independent phase-space calculation; the change in prepared gaps also
    # changes thermal noise, not only the rotation frequency.
    random = np.random.default_rng(2102026)
    frequencies = np.array([.3, .7, 1.2])
    gap_error, beta, duration, gamma = F(1, 100), F(2), F(3), F(1, 2)
    gaps = 2*frequencies+np.array([.01, -.01, .005])
    vectors = random.normal(size=(3, 3, 2))/10
    times = [-2., .5, 3.]
    radii = tuple(F(math.ceil(np.linalg.norm(v)*10**10), 10**10) for v in vectors)
    def exponent(source_gaps):
        rotated = []
        for v, time in zip(vectors, times, strict=True):
            angles = source_gaps*time/float(beta)
            c, s = np.cos(angles), np.sin(angles)
            rotated.append(np.column_stack((c*v[:, 0]-s*v[:, 1], s*v[:, 0]+c*v[:, 1])))
        noise = .5/np.tanh(source_gaps/2)
        def two(a, b):
            return np.sum(noise*np.sum(a*b, axis=1)
                          +.5j*(a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]))
        return sum(two(v, v)/2 for v in rotated)+sum(two(rotated[j], rotated[i])
                            for j in range(3) for i in range(j))
    bound = clock.growing_fluctuation_bound(4, 1000, radii, gap_floor=gamma,
        inverse_temperature=beta, clock_gap_error=gap_error, maximum_time=duration)
    assert abs(exponent(gaps)-exponent(2*frequencies)) < float(bound["clock_and_common_spectral_error_upper"])


def test_vanishing_single_copy_state_error_does_not_ensure_joint_field_limit():
    # A positive coherent preparation error O(N^-1/2) vanishes on each
    # register but leaves a finite collective mean phase. Parity is broken.
    p0, p1, strength = 2/3, 1/3, .1
    qnorm, gaussian = math.sqrt(1.5), math.exp(-.75)
    for copies in (100, 1000, 10000):
        perturbation = strength/math.sqrt(copies)
        assert p0*p1 > perturbation**2
        angle = qnorm/math.sqrt(copies)
        single = complex(math.cos(angle), 2*perturbation*math.sin(angle))
        assert abs(single**copies-gaussian) > .1
    expected = gaussian*np.exp(2j*strength*qnorm)
    assert abs(single**copies-expected) < 1e-4


@pytest.mark.parametrize("changes", (
    {"dimension": True}, {"copies": 0}, {"factor_radii": [F(1)]},
    {"factor_radii": (0.1,)}, {"gap_floor": F(0)},
    {"inverse_temperature": F(0)}, {"success_lower": F(0)},
    {"partition_upper": F(13)}, {"clock_gap_error": F(-1)}))
def test_growing_field_contract_rejects_unstated_or_nonpositive_domains(changes):
    values = dict(dimension=12, copies=100, factor_radii=(F(1),),
                  gap_floor=F(1), inverse_temperature=F(1))
    values.update(changes)
    with pytest.raises(ValueError):
        clock.growing_fluctuation_bound(**values)
