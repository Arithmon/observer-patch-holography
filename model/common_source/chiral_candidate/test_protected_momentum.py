"""Independent complete-output, protected-wire and frequency controls."""
from copy import deepcopy
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import action, protected_momentum as p


@pytest.fixture(scope="module")
def executed():
    return p.execute(steps=4)


def test_one_complete_momentum_family_retains_every_CAR_component(executed):
    report, arrays = executed
    assert report["field"]["total_CAR_modes"] == 630
    assert [row["CAR_modes"] for row in report["layouts"]] == [90, 180, 180, 180]
    assert all(not row["zero_momentum_duplicated"] and row["modes_discarded"] == 0 for row in report["layouts"])
    assert len(report["complete_complex_tensor_factors"]) == 192
    assert all(array.dtype.kind == "c" and np.isfinite(array).all() for array in arrays.values())
    assert report["same_field_current_attachment"]["same_retained_internal_generators_per_momentum"] == 12
    assert report["same_field_current_attachment"]["all_Fock_column_source_model_electric_charge_residual"] < 1e-10
    # Every global carrier appears exactly once per interval, and its
    # occupation wire has the actual 9-bit native-register stride.
    for interval in (0, 1):
        factors = [v for k, v in report["complete_complex_tensor_factors"].items() if k.startswith(f"t{interval}_")]
        assert sorted(j for row in factors for j in row["global_native_CAR_indices"]) == list(range(630))
        assert all(row["literal_decoded_occupation_wires"] == [9*j for j in row["global_native_CAR_indices"]] for row in factors)


def test_actual_source_words_have_one_time_step_and_three_axis_dispersion(executed):
    report, _ = executed
    assert report["source_interval_relation"]["long_steps"] == 2*report["source_interval_relation"]["short_steps"]
    assert report["source_interval_relation"]["all_Fock_column_doubled_word_operator_error"] < 1e-9
    for key in ("1", "2", "3"):
        observed = report["dispersion_readout"][key]
        assert abs(observed["measured_energy_squared_minus_rest_squared_over_clock_squared"]-1/36) < 2e-6
        assert report["independent_comparisons"][key]["source_frequency_over_clock_difference"] < 5e-6
        assert .28 < report["decoded"][key]["measured_visibility"] < .32
    assert not report["scope"]["operator_admission_certifies_frequency_uncertainty_or_a_physical_pole"]


def test_all_timeout_records_and_clock_reference_factors_are_explicit(executed):
    report, arrays = executed
    for key, packet in report["protocols"].items():
        rows = packet["terminal"]
        assert len(rows) == 40
        assert {row["flag"] for row in rows} == {"completed", "first_selector_timeout", "second_selector_timeout"}
        assert abs(math.fsum(row["probability"] for row in rows)-1) < 1e-7
        assert packet["work_preparations"] == 1 and packet["work_resets"] == 0
        assert arrays[key+"_probe_vectors"].shape[0] == 40
        assert arrays[key+"_clock_reference_columns"].shape == (40, 4, 4)
        # No reference purifier is measured or reconstructed from a number
        # expectation: every terminal retains its original complex columns.
        for flag in ("first_selector_timeout", "second_selector_timeout"):
            assert {row["selector_raw"] for row in rows if row["flag"] == flag} == {30, 31}


def test_decoder_does_not_consult_action_momenta_or_simulator_time(executed, monkeypatch):
    report, _ = executed
    def forbidden(*args, **kwargs):
        raise AssertionError("frequency decoder consulted hidden Hamiltonian data")
    monkeypatch.setattr(action, "mass_kernel", forbidden)
    monkeypatch.setattr(p.kinetic, "coordinate_speed", forbidden)
    monkeypatch.setattr(p.model, "parameters", forbidden)
    actual = p.infer_frequency(report["protocols"]["p1_t0"], report["protocols"]["p1_t1"],
                               report["protocols"]["p1_t0"]["clock_calibration"], report["phase_contract"])
    assert actual == report["decoded"]["1"]
    assert actual["local_inverse_sensitivity_diagnostics"]["absolute_derivative_wrt_long_probability"] > 1000
    assert not actual["recovered_frequency_interval_or_physical_pole_certified"]


def test_omitted_timeout_or_incompatible_clock_is_refused(executed):
    report, _ = executed
    short, long = (deepcopy(report["protocols"][key]) for key in ("p1_t0", "p1_t1"))
    short["terminal"] = [row for row in short["terminal"] if row["flag"] != "first_selector_timeout"]
    with pytest.raises(ValueError, match="complete unnormalized"):
        p.infer_frequency(short, long, short["clock_calibration"], report["phase_contract"])
    with pytest.raises(ValueError, match="alias window"):
        p.infer_frequency(long, long, long["clock_calibration"], {})


@pytest.mark.parametrize("fault", ["syndrome", "reordered", "extra_phase"])
def test_current_literal_wiring_is_checked_after_warming_all_caches(monkeypatch, fault):
    p.field_block((0, 1/3, 0), steps=1)
    original = p.decoded_source_word
    def wrong(word, sites):
        literal = original(word, sites)
        if fault == "syndrome":
            literal[0]["targets"][0] += 1
        elif fault == "reordered" and len(sites) == 2:
            a, b = (9*j for j in sites)
            literal = [{**event, "targets": [b if j == a else a if j == b else j for j in event["targets"]]} for event in literal]
        elif fault == "extra_phase":
            literal += [p.compiler.pulse("I", math.pi/2, [9*sites[0]])]
        return literal
    monkeypatch.setattr(p, "decoded_source_word", wrong)
    with pytest.raises(ArithmeticError, match="protected source word|factor order"):
        p.field_block((0, 1/3, 0), steps=1)


def test_mutable_metadata_cannot_poison_source_caches():
    block = p.field_block((0, 0, 0), steps=1)
    resources = dict(block.resources)
    factor_resources = dict(block.factors[0].resources)
    block.resources["native_pulses"] = -99
    block.factors[0].resources["native_pulses"] = -98
    block.layout["CAR_modes"] = -97
    fresh = p.field_block((0, 0, 0), steps=1)
    assert fresh.resources == resources and fresh.factors[0].resources == factor_resources
    assert fresh.layout["CAR_modes"] == 90
    with pytest.raises(ValueError):
        fresh.factors[0].unitary.flags.writeable = True


def test_global_native_resource_and_export_budget_is_scaled_to630_modes(executed):
    report, _ = executed
    for budget in report["protection_and_service_budgets"].values():
        assert budget["native_carriers"] == 630 and budget["quantum_work_code_bits"] == 5670
        assert not budget["floating_diagnostics_are_outward_certificates"]
        for law in ("C", "E"):
            row = budget["laws"][law]
            assert row["native_attempts_per_two_cycle_history"] == 1260
            assert row["source_resources"]["native_pulses"] > 3_000_000_000
            assert 0 < row["first_export_failure_amplitude_diagnostic"] < 1e-6
            assert 0 < row["complete_flagged_operator_error_diagnostic"] < row["fixed_complete_operator_admission"] == .01
            assert row["maximum_archive_and_export_ordinal_bits"] > 3800
            assert not row["physical_service_error_values_supplied"]
            alphabet = 8 if law == "C" else 16
            factors = budget["native_record_factors"][law]
            assert factors["native_carrier_factors_per_accepted_seam"] == 630
            assert factors["complete_native_record_addresses_per_accepted_cycle"] == 30*alphabet**630


def test_phase_corruption_with_identical_probabilities_fails_complete_output_oracle():
    blocks = tuple(p.field_block(momentum, steps=2) for momentum in p.momentum_menu())
    actual = p.protocol(blocks, 1, F(1, 4))
    damaged = {**actual, "terminal": [{**row, "probe_vector": 1j*row["probe_vector"]} for row in actual["terminal"]]}
    assert [row["probability"] for row in damaged["terminal"]] == [row["probability"] for row in actual["terminal"]]
    with pytest.raises(ArithmeticError, match="complete complex"):
        p._conditional_probe_comparison(damaged, actual)


def _receipt(tmp_path, executed):
    import hashlib, json
    report, arrays = executed
    document = p.canonical_report(report)
    document["source_bindings"] = p.source_bindings()
    artifact = tmp_path/"protected.npz"
    np.savez_compressed(artifact, **arrays)
    document["artifact"] = {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
    path = tmp_path/"protected.json"
    path.write_text(json.dumps(document, allow_nan=False))
    return path, document, arrays


def test_retained_original_custody_and_full_replay_comparison(tmp_path, executed):
    path, _, _ = _receipt(tmp_path, executed)
    retained, stored = p.verify_retained(path)
    p.compare_execution(retained, stored, *executed)
    # A backend's tiny last-bit array change need not match its byte hash.
    altered = dict(executed[1])
    key = next(iter(altered))
    altered[key] = altered[key]+1e-12j
    p.compare_execution(retained, stored, executed[0], altered)


def test_fresh_replay_rejects_coordinated_phase_tamper_and_rescaled_gate(tmp_path, executed):
    import hashlib, json
    path, document, original = _receipt(tmp_path, executed)
    changed = dict(original)
    key = "t0_p1_f1"
    changed[key] = original[key].copy()
    changed[key][2] *= 1j
    document["complete_complex_tensor_factors"][key]["complex_factor_sha256"] = hashlib.sha256(changed[key].tobytes()).hexdigest()
    artifact = path.with_suffix(".npz")
    np.savez_compressed(artifact, **changed)
    document["artifact"]["sha256"] = hashlib.sha256(artifact.read_bytes()).hexdigest()
    path.write_text(json.dumps(document))
    retained, stored = p.verify_retained(path)
    with pytest.raises(ValueError, match="complete complex"):
        p.compare_execution(retained, stored, *executed)
    document["Gaussian_and_full_Fock_factorization_admission"] = 0.
    with pytest.raises(ValueError, match="admission"):
        p.compare_execution(document, original, *executed)


@pytest.mark.parametrize("invalid", ["duplicate", "nan", "state_mass"])
def test_strict_retained_loader_refuses_ambiguous_or_damaged_evidence(tmp_path, executed, invalid):
    import json
    path, document, _ = _receipt(tmp_path, executed)
    if invalid == "duplicate":
        path.write_text('{"schema":"bad",'+path.read_text()[1:])
    elif invalid == "nan":
        document["decoded"]["0"]["energy_over_clock"] = float("nan")
        path.write_text(json.dumps(document))
    else:
        document["protocols"]["p1_t0"]["terminal"][0]["probability"] = .9
        path.write_text(json.dumps(document))
    with pytest.raises(ValueError):
        p.verify_retained(path)
