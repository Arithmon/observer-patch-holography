"""Fresh all66 execution, independent phase controls and complete failures."""
from fractions import Fraction as F
import json

import numpy as np
import pytest

from model.native_modular_reference import full_flow, local_preparation
from . import clock_services as s


@pytest.fixture(scope="module")
def report():
    # One full path executes AND certifies every original source array.
    return s.run()


def reference_marginal(density):
    r = len(density)//66
    return np.einsum("aras->rs", density.reshape(66, r, 66, r))


def test_all66_concrete_source_words_and_their_resource_counts(report):
    assert report["inventory"]["source_words_executed_and_certified"] == 66
    assert report["inventory"]["additional_words"] == 62
    assert len({row["word_sha256"] for row in report["words"]}) == 66
    for row in report["words"]:
        assert row["original_array_entries_checked"] == 512**2
        assert row["native_pulses"] == 2042
        assert row["pulse_kind_counts"] == {"H": 1020, "I": 1, "J": 510, "Z": 511}
        assert row["proper_code_load_events"]+row["proper_code_export_events"] == 4084
        assert F(row["operator_error_certified_below"]) == F(1, 10**9)
    assert report["resources"]["source_pulses_for_24_clock_requests"] == 49008
    assert report["resources"]["load_export_service_calls_for_24_clock_requests"] == 98016
    assert report["resources"]["source_program_utf8_bytes_total"] == sum(row["source_program_utf8_bytes"] for row in report["words"])


@pytest.mark.parametrize("core", [(1, 7), (1, 3), (10, 11)])
def test_new_core_array_against_independent_occupation_phase_oracle(core):
    # Independently use occupation counts, rather than the producer ratios().
    target = np.ones(512, complex)
    weights = (F(2, 225), F(1, 36), F(2, 45))
    from model.observer import native
    for address, state in zip(full_flow.support()[0], native.states(True), strict=True):
        weight = weights[state[core[0]]+state[core[1]]]
        target[address] = np.exp(-4j*np.log(float(60*weight)))
    actual = s.matrix(core)
    assert np.linalg.norm(actual-np.diag(target), "fro") < 1e-10
    assert s.certify(core, actual) < s.OPERATOR_GATE


def test_wrong_core_and_changed_source_phase_fail_the_certificate():
    with pytest.raises(ArithmeticError, match="operator gate"):
        s.certify((1, 7), s.matrix((0, 1)))
    tape = s.word((1, 7))
    changed = [dict(event) for event in tape]
    next(event for event in changed if event["kind"] == "Z")["angle"] += .125
    wrong = full_flow.compiler.execute(changed, 9)
    with pytest.raises(ArithmeticError, match="operator gate"):
        s.certify((1, 7), wrong)


def test_complete_service_on_complex_entangled_reference():
    rng = np.random.default_rng(1136)
    vector = rng.normal(size=198)+1j*rng.normal(size=198)
    vector /= np.linalg.norm(vector)
    rho = np.outer(vector, vector.conj())
    outputs = s.service((1, 7), rho)
    assert set(outputs) == {"success", "export_abort"}
    assert abs(np.trace(outputs["success"])+np.trace(outputs["export_abort"])-1) < 1e-10
    assert np.linalg.norm(reference_marginal(outputs["success"])+outputs["export_abort"]-reference_marginal(rho)) < 1e-10
    assert np.trace(outputs["export_abort"]).real <= float(s.OPERATOR_GATE**2)


def test_injected_padding_leakage_is_an_explicit_reference_bearing_abort():
    codes = full_flow.support()[0]
    padding = next(i for i in range(512) if i not in codes)
    raw = np.eye(512, dtype=complex)
    a, b = codes[0], padding
    raw[np.ix_((a, b), (a, b))] = [[np.sqrt(3)/2, -1/2], [1/2, np.sqrt(3)/2]]
    vector = np.zeros(132, complex)
    vector[0] = vector[3] = 1/np.sqrt(2)
    rho = np.outer(vector, vector.conj())
    outputs = s.export_outputs(raw, rho)
    assert np.allclose(outputs["export_abort"], np.diag([1/8, 0]), atol=1e-15)
    assert np.trace(outputs["success"]).real == pytest.approx(7/8)
    assert np.allclose(reference_marginal(outputs["success"])+outputs["export_abort"], np.eye(2)/2)
    with pytest.raises(ArithmeticError, match="operator gate"):
        s.certify((1, 7), raw)


def test_preparation_failure_blocks_bypass_clock_unchanged():
    vector = s._probe().ravel()
    rho = np.outer(vector, vector.conj())
    blocks = {flag: rho/4 for flag in local_preparation.FLAGS}
    outputs = s.after_preparation(blocks, (1, 7))
    for flag in local_preparation.FLAGS:
        if flag != "success":
            assert np.array_equal(outputs["preparation", flag], blocks[flag])
    assert set(outputs) == {("clock", "success"), ("clock", "export_abort"),
                            *(('preparation', flag) for flag in local_preparation.FLAGS if flag != "success")}
    assert abs(sum(np.trace(value) for value in outputs.values())-1) < 1e-10
    with pytest.raises(ValueError, match="all predecessor"):
        s.after_preparation({"success": rho}, (1, 7))
    wrong_reference = dict(blocks)
    wrong_reference["herald_failure"] = np.eye(66, dtype=complex)
    with pytest.raises(ValueError, match="share one H66 tensor reference"):
        s.after_preparation(wrong_reference, (1, 7))


def test_complete_record_covariance_and_exported_history_have_new_services(report):
    for law in ("C", "E"):
        attempt = report["complete_one_attempt"][law]
        assert attempt["semantic_native_slots_retained"] == 497
        assert [1, 7] in [list(core) for core in attempt["new_cores_exercised"]]
        assert attempt["clock_after_export_abort_retained_per_native_outcome"]
        history = report["exported_history"][law]
        assert history["final_core"] == [1, 3]
        assert history["exported_records"] == 3 and history["retained_records"] == 1
        assert history["history_corner_rank"] >= 2
        assert F(history["branch_weight"]) == F(1, 120**4)
        assert not history["reference_reprepared"]


def test_arbitrary_reference_error_composes_without_branch_count_factor():
    one, many = s.error_budget(1), s.error_budget(24)
    gain = F(one["per_word_channel_norm_upper"])
    delta = F(one["per_word_half_diamond_upper"])
    assert F(many["complete_interleaved_half_diamond_upper"]) == delta*sum(gain**j for j in range(24))
    assert F(many["complete_interleaved_half_diamond_upper"]) < F(25, 10**9)
    with pytest.raises(ValueError):
        s.error_budget(True)


def test_source_array_copy_cannot_poison_a_later_service():
    first = s.matrix((1, 7))
    old = first[0, 0]
    first[0, 0] = 100
    assert s.matrix((1, 7))[0, 0] == old


def test_receipt_replay_and_scoped_finite_theorem_binding(report):
    s.check_receipt(s.OUTPUT, report)
    assert all(path in report["source_sha256"] for path in s.THEOREMS)
    assert report["scope"]["service_inventory_extended_under_successor_contract"]
    assert not report["scope"]["new_native_primitive"]
    assert not report["scope"]["integrated_cp01_accepted"]


def test_retained_receipt_does_not_accept_type_changed_service_count(tmp_path):
    expected = json.loads(s.OUTPUT.read_text())
    changed = json.loads(s.OUTPUT.read_text())
    changed["inventory"]["additional_words"] = True
    path = tmp_path/"changed.json"
    path.write_text(json.dumps(changed))
    with pytest.raises(ValueError):
        s.check_receipt(path, expected)
