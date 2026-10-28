"""Response identifiability, charge-compatible apparatus and continued-state controls."""
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path, PureWindowsPath
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from . import record_two_point as r


def test_source_binding_keys_are_portable_under_windows_paths(monkeypatch):
    inherited = "model/common_source/chiral_candidate/protected_momentum.py"
    monkeypatch.setattr(r.protected, "source_bindings", lambda: {inherited: "unused"})
    relative_to = Path.relative_to
    with monkeypatch.context() as patch:
        patch.setattr(Path, "relative_to", lambda self, *args, **kwargs:
                      PureWindowsPath(relative_to(self, *args, **kwargs)))
        actual = r.source_bindings()
    names = (inherited, "model/common_source/chiral_candidate/record_two_point.py",
             "model/common_source/chiral_candidate/test_record_two_point.py",
             "model/common_source/hierarchy_diagnosis/audit.py",
             "model/common_source/replay.py", "tools/strict_json.py")
    assert actual == {name: hashlib.sha256((r.ROOT/name).read_bytes()).hexdigest() for name in names}


@pytest.fixture(scope="module")
def executed():
    return r.execute()


def test_all_detector_branches_carry_their_unchanged_quantum_outputs(executed):
    report, arrays = executed
    for key, packet in report["records"].items():
        work = {name[len(key)+1:]: value for name, value in arrays.items() if name.startswith(key+"_")}
        r.validate_outputs(packet, work)
        assert len(packet["contexts"]) == 16
        for context in packet["contexts"]:
            assert len(context["terminal"]) == 16
            assert abs(math.fsum(row["probability"] for row in context["terminal"])-1) < 1e-10
        contract = packet["instrument_contract"]
        assert contract["physical_mode_charges"] == [-1, 1, -1]
        assert contract["work_preparations_per_context"] == 1
        assert contract["work_resets_per_context"] == 0
        assert not contract["charged_external_phase_reference_required"]
        assert not contract["native_attempt_abort_and_export_failure_tree_included_here"]
    calibration = {name[len("calibration_"):]: value for name, value in arrays.items() if name.startswith("calibration_")}
    r.validate_outputs(report["equal_time_calibration"], calibration)


def test_decoder_is_blind_to_source_coefficient_and_programmed_duration(executed, monkeypatch):
    report, _ = executed
    def forbidden(*args, **kwargs):
        raise AssertionError("inverse consulted the source or apparatus")
    for name in ("source_pair", "source_unitary", "apparatus", "measure"):
        monkeypatch.setattr(r, name, forbidden)
    monkeypatch.setattr(r.protected, "field_layout", forbidden)
    monkeypatch.setattr(r.history, "_clock_preparation", forbidden)
    for key, packet in report["records"].items():
        assert r.decode(packet) == report["decoded"][key]


@pytest.mark.parametrize("fault", ["omitted", "duplicate", "renormalized", "clock", "phase_prior"])
def test_missing_records_and_incompatible_inverse_contract_are_rejected(executed, fault):
    packet = deepcopy(executed[0]["records"]["0"])
    if fault == "omitted":
        packet["contexts"][0]["terminal"].pop()
    elif fault == "duplicate":
        packet["contexts"][0]["terminal"][-1] = packet["contexts"][0]["terminal"][0]
    elif fault == "renormalized":
        packet["contexts"][0]["terminal"][0]["probability"] += .1
    elif fault == "clock":
        packet["clock_calibration"]["contrast"] /= 2
    else:
        packet["branch_contract"]["generator_eigenphases_strictly_inside_principal_interval"] = False
    with pytest.raises(ValueError):
        r.decode(packet)


def test_zero_time_calibrates_full_identity_car_response(executed):
    calibration = r.decode(executed[0]["equal_time_calibration"], zero_time=True)
    np.testing.assert_allclose(r._unpack(calibration["transfer"]), np.eye(2), atol=1e-11)
    for row in executed[0]["decoded"].values():
        residues = list(map(r._unpack, row["canonical_residues"]))
        np.testing.assert_allclose(sum(residues), np.eye(2), atol=1e-12)
        for residue in residues:
            assert np.linalg.eigvalsh(residue).min() > -1e-12
        assert row["generator_identified_only_within_supplied_alias_window"]
        assert not row["interacting_wavefunction_renormalization_recovered"]


def test_normalized_records_with_nonunitary_transfer_are_refused(executed):
    packet = deepcopy(executed[0]["records"]["0"])
    rows = packet["contexts"][0]["terminal"]
    occupied = sum(row["probability"] for row in rows if row["first"])
    for row in rows:
        row["probability"] *= .9 if row["first"] else (1-.9*occupied)/(1-occupied)
    assert abs(sum(row["probability"] for row in rows)-1) < 1e-10
    with pytest.raises(ValueError, match="closed canonical"):
        r.decode(packet)


def test_sign_indistinguishable_in_rabi_marginals_is_resolved_by_neutral_reference():
    h, delta, _ = r.source_pair(.25, 1.)
    plus = r.protected._small_evolution(h, delta, float(r.INTERVAL), 32).unitary
    minus = r.protected._small_evolution(h, -delta, float(r.INTERVAL), 32).unitary
    # A pair initially empty has exactly the same occupation probabilities.
    np.testing.assert_allclose(abs(plus[:, 0])**2, abs(minus[:, 0])**2, atol=2e-12)
    decoded = []
    for small in (plus, minus):
        whole = np.kron(np.eye(2), small)
        packet, arrays = r.measure(whole)
        r.validate_outputs(packet, arrays)
        decoded.append(r.decode(packet))
    generators = [r._unpack(row["generator_over_clock"]) for row in decoded]
    assert generators[0][0, 1].real > .1
    assert generators[1][0, 1].real < -.1
    np.testing.assert_allclose(decoded[0]["poles_over_clock"], decoded[1]["poles_over_clock"], atol=2e-11)


def test_continuation_state_error_is_detected_even_when_every_record_weight_agrees(executed):
    report, stored = executed
    packet = report["records"]["0"]
    arrays = {name[2:]: value for name, value in stored.items() if name.startswith("0_")}
    damaged = {name: value.copy() for name, value in arrays.items()}
    for name, array in damaged.items():
        if name.endswith("_work"):
            array[:, 6] *= -1
    # The fault is neutral, preserves all norms and leaves the tomography
    # marginals unchanged. A phase-sensitive subsequent probe could see it.
    r.validate_outputs(packet, damaged)
    assert r.complete_output_distance(packet, arrays, packet, damaged) > .1
    global_phase = {name: (1j*value if name.endswith("_work") else value)
                    for name, value in arrays.items()}
    assert r.complete_output_distance(packet, arrays, packet, global_phase) < 1e-12
    damaged["i0_j0_X+_work"][0] *= 2
    with pytest.raises(ValueError, match="probability/state"):
        r.validate_outputs(packet, damaged)


def test_joint_response_predicts_withheld_momenta_higgs_controls_and_fixed_step_time(executed):
    report, _ = executed
    fit = report["joint_response"]
    np.testing.assert_allclose(r._unpack(fit["momentum_derivative"]), np.diag([-.5, .5]), atol=1e-7)
    assert fit["maximum_training_operator_residual"] < 2e-8
    for row in report["held_out_validation"]:
        assert row["operator_residual"] < 7e-8
        assert row["second_duration_operator_residual"] < 7e-8
        assert row["complete_continued_record_quantum_output_distance"] < 4e-8
        assert row["second_duration_complete_continued_output_distance"] < 5e-8
        assert row["same_elementary_word_three_short_vs_two_long_operator_error"] < 2e-10
    assert report["source_interval_relation"]["long_steps"] == 3*report["source_steps"]//2


def test_two_hop_kinetics_pass_training_but_fail_withheld_shell(executed):
    report, _ = executed
    for key, setting in enumerate(r.TRAINING):
        alternative, _ = r.source_unitary(*setting, law="two_hop_control")
        alternative_packet, _ = r.measure(alternative)
        np.testing.assert_allclose(r._unpack(r.decode(alternative_packet)["transfer"]),
                                   r._unpack(report["decoded"][str(key)]["transfer"]), atol=2e-11)
    alternative, _ = r.source_unitary(*r.HELD_OUT[0], law="two_hop_control")
    packet, _ = r.measure(alternative)
    assert np.linalg.norm(r._unpack(r.decode(packet)["generator_over_clock"])
                          -r._unpack(report["decoded"]["5"]["generator_over_clock"]), 2) > .1


def test_two_commensurate_times_do_not_remove_unbounded_generator_aliases(executed):
    row = executed[0]["decoded"]["0"]
    generator = r._unpack(row["generator_over_clock"])
    time = row["clock_interval"]
    residue = r._unpack(row["canonical_residues"][0])
    aliased = generator+4*math.pi/time*residue
    assert np.linalg.norm(aliased-generator, 2) > 30
    for duration in (time, 1.5*time):
        np.testing.assert_allclose(expm(-1j*duration*aliased), expm(-1j*duration*generator), atol=2e-10)
    # These observations cannot verify the independently supplied norm prior.
    assert executed[0]["source_interval_relation"]["two_commensurate_times_do_not_remove_unbounded_generator_aliases"]


def test_apparatus_cache_cannot_be_poisoned_by_a_consumer():
    basis, charge, rotations, resources = r.apparatus()
    with pytest.raises(ValueError):
        basis.flags.writeable = True
    resources[0, "X+"]["native_pulses"] = -1
    rotations.clear()
    _, _, fresh, costs = r.apparatus()
    assert len(fresh) == 8 and costs[0, "X+"]["native_pulses"] > 0
    for forward, reverse in fresh.values():
        np.testing.assert_allclose(charge@forward, forward@charge, atol=1e-12)
        np.testing.assert_allclose(reverse@forward, np.eye(8), atol=1e-12)


def test_rank_deficient_joint_interventions_are_rejected(executed):
    row = executed[0]["decoded"]["0"]
    with pytest.raises(ValueError, match="independent momentum and Higgs"):
        r.fit_joint_response([(0, 1)]*3, [row]*3)


def test_unresolved_poles_are_not_given_arbitrary_rank_one_physical_labels():
    packet, _ = r.measure(np.eye(8))
    with pytest.raises(ValueError, match="whole-subspace residues"):
        r.decode(packet)


def test_rehashed_quantum_output_tampering_still_fails_fresh_replay(executed, tmp_path, monkeypatch):
    report, arrays = executed
    receipt = deepcopy(report)
    receipt["source_bindings"] = r.source_bindings()
    changed = {name: value.copy() for name, value in arrays.items()}
    changed["0_i0_j0_X+_work"][:, 6] *= -1
    artifact = tmp_path/"records.npz"
    np.savez_compressed(artifact, **changed)
    receipt["artifact"] = {"filename": artifact.name, "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest()}
    path = tmp_path/"records.json"
    path.write_text(json.dumps(receipt, allow_nan=False))
    monkeypatch.setattr(sys, "argv", ["record_two_point", "--check", str(path)])
    with pytest.raises(ValueError, match="quantum array differs"):
        r.main()
