"""Independent closed-form/CAR/frame and hierarchy controls for transport."""
from fractions import Fraction as F
from copy import deepcopy
import math
import json
from unittest.mock import patch

import numpy as np
import pytest
from scipy.linalg import expm

from . import calibrated_transport as t


@pytest.mark.parametrize("mass,e", [(F(7, 10), F(3, 5)), (F(7, 10), F(-3, 5)),
    (F(1, 10**6), F(2)), (F(170), F(-2)), (F(0), F(1)), (F(0), F(0))])
def test_exact_pair_all_four_columns_and_odd_global_phase(mass, e):
    time = F(3, 7)
    block = t.moving_pair(mass, e, time)
    h = np.array([[0, 0, 0, float(mass)], [0, float(e), 0, 0],
                  [0, 0, float(e), 0], [float(mass), 0, 0, 2*float(e)]], complex)
    np.testing.assert_allclose(block["unitary"], expm(-1j*float(time)*h), atol=4e-13, rtol=0)
    for odd in (1, 2):
        assert block["unitary"][odd, odd] == pytest.approx(np.exp(-1j*float(e*time)), abs=2e-14)
    assert block["exact_binary_operator_error_upper"] < F(1, 10**12)
    assert block["resources"]["native_pulses"] == 79
    assert not block["global_phase_removed_or_fitted"]
    with pytest.raises(ValueError):
        block["unitary"].flags.writeable = True


def test_large_winding_needs_no_trotter_or_mass_dependent_word_length():
    for mass in (F(1, 10**6), F(1), F(170)):
        block = t.moving_pair(mass, F(2), F(2**35, 512))
        assert block["resources"]["native_pulses"] == 79
        assert block["exact_binary_operator_error_upper"] < t.ADMISSION
        assert block["exact_binary_effective_read_operator_error_upper"] < t.ADMISSION
    # Dropping the odd-sector scalar phase would be an all-input error.
    block = t.moving_pair(F(1), F(2), F(1, 5))
    phase_removed = block["unitary"]*np.exp(2j/5)
    assert t.observer._operator_error(phase_removed, t._target_boxes(F(1), F(2), F(1, 5))) > F(1, 10)


@pytest.mark.parametrize("axis", ("x", "y", "z"))
def test_actual_spin_frame_has_correct_axis_and_fermionic_two_particle_sector(axis):
    block = t.helicity_frame(axis)
    S, U = block["spinor_unitary"], block["Fock_unitary"]
    sigma = t.kinetic.pauli_matrices()[("x", "y", "z").index(axis)]
    np.testing.assert_allclose(S.conj().T@sigma@S, np.diag([1, -1]), atol=4e-14, rtol=0)
    np.testing.assert_allclose(S.T@t.kinetic.spin_epsilon()@S, t.kinetic.spin_epsilon(), atol=4e-14, rtol=0)
    assert abs(U[3, 3]-np.linalg.det(S)) < 1e-14


def test_all_180_calibrated_Nambu_columns_share_spin_family_and_current_frame():
    rows = t.symbol_intertwining()
    assert len(rows) == 6
    assert all(row["complete_Nambu_columns"] == 180 for row in rows)
    assert max(row["source_frame_symbol_difference_diagnostic"] for row in rows) < 2e-10
    assert len(t.momentum_menu()) == 13
    assert t.momentum_menu().count((F(), F(), F())) == 1


@pytest.mark.parametrize("axis", ("x", "y", "z"))
def test_actual_256_column_opposite_momentum_source_matches_independent_CAR(axis):
    result = t.charged_momentum_pair(F(7, 10), F(1), axis, F(1, 8))
    np.testing.assert_allclose(result["unitary"], result["target"], atol=2e-10, rtol=0)
    assert result["physical_CAR_modes"] == 8 and result["all_Fock_columns"] == 256
    assert not result["global_phase_removed_or_fitted"]


@pytest.mark.parametrize("bad", (True, .7, float("nan"), -1))
def test_original_mass_domain_refuses_floating_or_negative(bad):
    with pytest.raises(ValueError):
        t.moving_pair(bad, F(1), F(1, 8))


@pytest.mark.parametrize("mass,p", [(F(7, 10), F(1)), (F(1, 10**6), F(2)), (F(170), F(1)), (F(), F(1))])
def test_fixed_probe_and_three_pointer_quadratures_match_independent_Bloch_solution(mass, p):
    time = F(3, 10)
    block = t.moving_pair(mass, -p, time)
    prep = t.fixed_probe()["unitary"]
    np.testing.assert_allclose(prep[:, 0], np.array([1, 0, 0, 1j])/math.sqrt(2), atol=2e-14, rtol=0)
    E = math.hypot(float(mass), float(p))
    expected = {"X": -float(p)/E*math.sin(2*E*float(time)),
                "Y": math.cos(2*E*float(time)), "Z": float(mass)/E*math.sin(2*E*float(time))}
    for axis in t.QUADRATURES:
        outputs, _ = t.fixed_read(block["unitary"]@prep, axis)
        actual = float(t.observer._norm(outputs[0][:, 0])-t.observer._norm(outputs[1][:, 0]))
        assert actual == pytest.approx(expected[axis], abs=1e-12)


def test_signed_norm_through_zero_and_conflicting_sign_refusal():
    assert t.signed_sine_interval((F(3, 10),)*2, (F(4, 10),)*2) == (F(1, 2),)*2
    assert t.signed_sine_interval((F(-3, 10),)*2, (F(-4, 10),)*2) == (F(-1, 2),)*2
    lo, hi = t.signed_sine_interval((F(-1, 1000), F(1, 1000)), (F(-1, 1000), F(1, 1000)))
    assert lo < 0 < hi and hi < F(1, 500)
    with pytest.raises(ValueError, match="opposite certified signs"):
        t.signed_sine_interval((F(1, 10), F(2, 10)), (F(-2, 10), F(-1, 10)))


def test_classical_rotation_preserves_large_winding_phase_without_quantum_controls():
    phase = F(10**10)+F(1, 3)
    _, residual, _ = t.observer._reduced((phase, phase))
    cosine, sine = t._trig(residual, True), t._trig(residual)
    c, q = t.rotate_intervals(cosine, sine, phase)
    assert c[0] <= 1 <= c[1] and q[0] <= 0 <= q[1]


def test_all_modes_disjoint_and_literal_global_occupation_wiring_refuses_warm_corruption():
    layout = t.source_layout()
    assert sorted(layout["canonical_to_original"]) == list(range(1170))
    assert len(layout["pairs"]) == 546 and len(layout["neutral"]) == 78
    word = t.moving_pair(F(7, 10), F(-1), F(1, 4))["word"]
    literal = t.bind_source_word(word, (91, 92))
    assert {target for event in literal for target in event["targets"]} <= {9*91, 9*92}
    def wrong(w, sites):
        return t.COMPILER.remap(w, tuple(9*site+1 for site in sites))
    with patch.object(t.protected, "decoded_source_word", side_effect=wrong):
        with pytest.raises(ArithmeticError, match="syndrome wire"):
            t.bind_source_word(word, (91, 92))


def test_source_route_metadata_preserves_exact_json_types():
    route = t._route_plan((91, 200, 93))
    t.observer._same_report(route, json.loads(json.dumps(route)))


@pytest.fixture(scope="module")
def independent_records():
    """Correlatednine-lane analytic counts, not action/source outputs."""
    levels = 14
    counts = np.zeros((3, levels, 3, 32, 512, 4), dtype=np.int64)
    rng = np.random.default_rng(2468)
    for protocol, p in enumerate((0, 1, 2)):
        E = math.hypot(1, p)
        for level in range(levels):
            time = 2**level/512
            theta = 2*E*time
            means = (-p/E*math.sin(theta), math.cos(theta), 1/E*math.sin(theta))
            clock_means = (.75*math.cos(time), .75*math.sin(time))
            # Clockbit0meansS, highbitmeansC; scorebit0=-1,bit1=+1.
            clocks = np.array([math.prod((1+(1 if (record >> bit)&1 else -1)*clock_means[1-bit])/2
                                          for bit in range(2)) for record in range(4)])
            for setting, mean in enumerate(means):
                law = np.r_[(15/16)*np.outer([(1+mean)/2, (1-mean)/2], clocks).ravel(), 1/16]
                sample = rng.multinomial(t.observer.TRIALS_PER_SETTING, law)
                counts[protocol, level, setting, 0, 0] = sample[:4]
                counts[protocol, level, setting, 0, 511] = sample[4:8]
                counts[protocol, level, setting, 31, 0, 0] = sample[8]
    return counts


def test_count_only_inverse_with_unknown_shared_amplitudes_and_no_packet_source_access(independent_records):
    with patch.object(t.matter, "parameters", side_effect=AssertionError("packet read")), \
         patch.object(t.observer, "catalogue", side_effect=AssertionError("mass catalogue read")), \
         patch.object(t.family_source, "source_frame", side_effect=AssertionError("source read")):
        decoded = t.decode_records(independent_records)
    for protocol, p in enumerate((0, 1, 2)):
        for box in decoded[protocol]["energy_GeV_intervals"].values():
            assert float(box[0]) < math.hypot(1, p) < float(box[1])


def test_bad_original_counts_and_absent_visibility_refuse(independent_records):
    for bad in (independent_records.astype(float), np.ma.array(independent_records), independent_records.tolist(), independent_records[:, :, :, :31]):
        with pytest.raises(ValueError):
            t.decode_records(bad)
    bad = np.zeros((3, 1, 3, 32, 512, 4), np.int64)
    bad[:, :, :, 31, 0, 0] = t.observer.TRIALS_PER_SETTING
    with pytest.raises(ValueError, match="unresolved"):
        t.decode_records(bad)


def test_checker_rejects_complex_phase_zero_extra_arrays_and_fixed_admission_forgery():
    report = {"error_contract": {"maximum_binary_pair_error_upper": str(t.ADMISSION), "maximum_binary_complete_score_bias_upper": str(t.observer.NUMERIC_SCORE_ADMISSION),
                "backend_pair_error_diagnostic": "0", "backend_complete_score_bias_diagnostic": "0",
                "binary_all_column_pair_error_admission": str(t.ADMISSION)},
              "frame_source_word_grammar": {"backend_FSWAP_error_diagnostic": "0", "exact_binary_FSWAP_operator_error_upper": "1/10000000000000", "four_frame_score_error_upper_conditional": "1/10000000"}}
    arrays = {"complex": np.array([[1, 1j]], complex), "counts": np.array([1, 2], np.int64)}
    for kind in ("phase", "zero", "extra"):
        forged = {key: value.copy() for key, value in arrays.items()}
        if kind == "extra":
            forged["extra"] = np.zeros(1)
        else:
            forged["complex"][0, 1] = -1j if kind == "phase" else 0
        with pytest.raises(ValueError):
            t.compare_execution(report, forged, report, arrays)
    forged = deepcopy(report)
    forged["error_contract"]["binary_all_column_pair_error_admission"] = "1"
    with pytest.raises(ValueError):
        t.compare_execution(forged, arrays, report, arrays)
    forged = deepcopy(report)
    forged["error_contract"]["maximum_binary_pair_error_upper"] = "0"
    with pytest.raises(ValueError):
        t.compare_execution(forged, arrays, report, arrays)
    forged = deepcopy(report)
    forged["frame_source_word_grammar"]["four_frame_score_error_upper_conditional"] = "0"
    with pytest.raises(ValueError):
        t.compare_execution(forged, arrays, report, arrays)


def test_checker_rejects_duplicate_nonfinite_and_wrong_custody_before_source_execution(tmp_path):
    path = tmp_path/"result.json"
    path.write_text('{"source_bindings":{},"source_bindings":{}}')
    with pytest.raises(ValueError, match="duplicate JSON"):
        t.check(path)
    path.write_text('{"bad":NaN}')
    with pytest.raises(ValueError, match="nonfinite"):
        t.check(path)
    path.write_text('{"source_bindings":{}}')
    with patch.object(t, "source_bindings", return_value={"own": "sha"}), patch.object(t, "execute", side_effect=AssertionError("must not execute")):
        with pytest.raises(ValueError, match="source bindings"):
            t.check(path)
