"""Independent complete-output and source-word tests for mobile recovery."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import mobile, mobile_service as service


def test_independent_primitive_addresses_and_compiled_layout():
    result = service.compiled_layout_certificate()
    assert result["core_encoders"] == 60
    assert result["encoder_slots"] == 360
    assert result["native_routes_full78_and_padding"] == 30
    assert result["positive_branch_inverse_schedules"] > 10000
    assert result["logical_holonomy_on_every_closed_corrected_core_path"] == "I2"
    assert result["primitive_address_contract"]["binary_states"] == 66


def test_all_target_encoders_and_recoveries_are_full_unitary_extensions():
    # Integer schedule verification is independent of floating source replay.
    # On each positive branch the two decoded uniform columns land on the
    # same two logical slots, with a normalized target-code encoder afterwards.
    resets = motions = 0
    for args in service._addresses():
        rows = mobile.branch_images(*args)
        if not rows[0]:
            continue
        status = mobile.transition(*args)
        target = tuple(status["target_core"])
        assert target in mobile.legal_cores()
        assert len(rows[0]) == len(rows[1])
        decoded = service._integer_permutation(service.branch_permutation(*args), True)
        assert [[decoded[x] for x in row] for row in rows] == [
            [logical+8*term for term in range(len(rows[0]))] for logical in (0, 1)]
        resets += status["reset_performed"]
        motions += target != args[0] and not status["reset_performed"]
    assert resets > 0 and motions > resets


@pytest.fixture(scope="module")
def replay():
    return service.source_replay_report()


def test_native_executed_motion_reset_and_closed_exchange_keep_all_outputs(replay):
    assert {tuple(row["core"]) for row in replay["source_executed_prototypes"]} == {(0, 1), (0, 2)}
    histories = replay["complete_stopped_histories"]
    moving = histories["motion_then_forbidden_reset"]
    assert moving["steps"][0]["core_before"] != moving["steps"][0]["target_core"]
    assert not moving["steps"][0]["reset_performed"]
    assert moving["steps"][1]["reset_performed"]
    assert moving["steps"][1]["target_core"] == moving["steps"][1]["seam"]
    loop = histories["closed_core_exchange"]
    assert loop["steps"][0]["core_before"] == loop["steps"][-1]["target_core"]
    assert not any(row["reset_performed"] for row in loop["steps"])
    for row in histories.values():
        assert row["terminal_exact_probability_sum"] == "1"
        assert row["terminal_outputs"] == sum(step["retained_raw_slots"] for step in row["steps"])-len(row["steps"])+1
        assert row["maximum_terminal_source_operator_error"] < 1e-8
        assert row["selected_terminal_decoder_error"] < 1e-8


def test_unknown_entangled_reference_and_zero_abort_extensions():
    reference = np.asarray([[1, 1j, .25], [1j, -1, .75]], complex)
    for core in service.PROTOTYPES:
        output = service.decoding(core, service.encoding(core, reference))
        assert np.linalg.norm(output[:2]-reference) < 1e-8
        assert np.linalg.norm(output[2:]) < 1e-8
    core = service.PROTOTYPES[0]
    state = np.zeros((512, 2), complex)
    state[211, 0] = 1
    state[501, 1] = 1j
    args = next(args for args in service._addresses() if args[0] == core
        and mobile.transition(*args)["admitted"] and not mobile.transition(*args)["positive"])
    assert np.array_equal(service.apply_recovery(*args, state), state)
    with pytest.raises(ValueError):
        service.apply_recovery(core, "C", (0, 1), 0, 3, state)
    with pytest.raises(ValueError):
        service.encoding(core, [True, False])


def test_every_mobile_word_has_counted_source_error_and_resources():
    bound = service.source_composition_report()
    for law in ("C", "E"):
        row = bound["laws"][law]
        assert row["encoded_t_attempt_operator_error_upper"]["100"] < 1e-5
        assert row["all_sixty_cores_thirty_seams_and_raw_slots_counted"]
    costs = service.resource_report()
    assert len(costs["encoders"]) == 60
    assert costs["additional_recovery_scratch_qubits"] == 0
    for law in ("C", "E"):
        row = costs["laws"][law]
        assert F(row["stationary_mean_native_plus_recovery_pulses"]) > 0
        assert row["maximum_one_attempt_native_plus_recovery_pulses"] > F(row["stationary_mean_native_plus_recovery_pulses"])
        assert row["native_tape_bits_per_attempt"] == (8 if law == "C" else 9)
        assert row["native_blank_preparations_per_attempt"] == {
            "coin_qubits": 1, "record_qubits": 2 if law == "C" else 3}
        counted = row["stationary_mean_native_plus_recovery_resources"]
        assert counted["proper_code_load_events"] == counted["native_pulses"]
        assert counted["proper_code_export_events"] == counted["native_pulses"]
        assert row["stationary_core_reset_rate"] == "1/150"


def test_wrong_encoder_permutation_is_rejected(monkeypatch):
    original = service.encoder_permutation
    def wrong(core):
        result = list(original(core))
        result[0], result[1] = result[1], result[0]
        return tuple(result)
    monkeypatch.setattr(service, "encoder_permutation", wrong)
    with pytest.raises(ArithmeticError, match="logical support"):
        service.compiled_layout_certificate()


def test_wrong_recovery_permutation_is_rejected(monkeypatch):
    original = service.branch_permutation
    def wrong(*args):
        result = list(original(*args))
        result[0], result[1] = result[1], result[0]
        return tuple(result)
    monkeypatch.setattr(service, "branch_permutation", wrong)
    with pytest.raises(ArithmeticError, match="branch logical coordinate"):
        service.compiled_layout_certificate()


def test_mutated_emitted_schedule_does_not_reuse_a_warm_target_admission(monkeypatch):
    images = service.encoder_permutation((0, 1))
    service._integer_permutation(images)
    monkeypatch.setattr(service.routing, "permutation_edges", lambda *args, **kwargs: ())
    with pytest.raises(ArithmeticError, match="declared map"):
        service._integer_permutation(images)


def test_mutated_native_address_cannot_relabel_the_oracle(monkeypatch):
    service.encoder_permutation((0, 1))  # Deliberately warm the public cache.
    original = service.circuit.encode
    monkeypatch.setattr(service.circuit, "encode", lambda state: original(state)^1)
    with pytest.raises(ArithmeticError, match="independent register layout"):
        service.primitive_address_certificate()
    with pytest.raises(ArithmeticError, match="independent register layout"):
        service.encoder_permutation((0, 1))
    with pytest.raises(ArithmeticError, match="independent register layout"):
        service.code_matrix((0, 1))


def test_wrong_recovery_source_output_fails_complete_replay(monkeypatch, replay):
    original = service.apply_recovery
    def wrong(*args):
        value = original(*args).copy()
        value[[0, 1]] = value[[1, 0]]
        # A phase on one logical column is visible against an untouched
        # reference even if its classical branch probability is unchanged.
        value[:, 0] *= -1
        return value
    monkeypatch.setattr(service, "apply_recovery", wrong)
    with pytest.raises(ArithmeticError, match="complete source branch"):
        service.execute_stopped_history(service.witness_schedules()["motion_then_forbidden_reset"])


def test_nonfinite_source_output_is_not_hidden_in_a_maximum(monkeypatch, replay):
    monkeypatch.setattr(service, "apply_recovery", lambda *args: np.full_like(args[-1], np.nan))
    with pytest.raises((ArithmeticError, np.linalg.LinAlgError)):
        service.execute_stopped_history(service.witness_schedules()["motion_then_forbidden_reset"])


@pytest.mark.parametrize("warm", (False, True))
@pytest.mark.parametrize("bad_core", ((False, 1), (0, 1.), [0, 1], (0, 3)))
def test_public_core_apis_reject_invalid_scalar_aliases_before_caches(warm, bad_core):
    service._encoder_permutation.cache_clear()
    service._code_matrix.cache_clear()
    if warm:
        service.encoder_permutation((0, 1))
        service.code_matrix((0, 1))
    logical = np.asarray([1., 0.])
    source = np.eye(512, 2, dtype=complex)
    calls = (
        lambda: service.encoder_permutation(bad_core),
        lambda: service.code_matrix(bad_core),
        lambda: service.ideal_code_matrix(bad_core),
        lambda: service.encoding(bad_core, logical),
        lambda: service.decoding(bad_core, source),
        lambda: service.encoding_resources(bad_core),
        lambda: service.repair_branches(bad_core, "C", (0, 1), source),
    )
    for invoke in calls:
        with pytest.raises(ValueError):
            invoke()


@pytest.mark.parametrize("warm", (False, True))
@pytest.mark.parametrize("position,bad", (
    (0, (False, 1)), (0, (0, 1.)), (0, (0, 3)),
    (1, True), (2, (False, 1)), (2, (0, 1.)),
    (3, False), (3, 0.), (4, False), (4, 0.),
))
def test_public_branch_apis_reject_invalid_addresses_before_caches(warm, position, bad):
    args = [(0, 1), "C", (0, 1), 0, 0]
    service._branch_permutation.cache_clear()
    if warm:
        service.branch_permutation(*args)
    args[position] = bad
    source = np.eye(512, 2, dtype=complex)
    for invoke in (
            lambda: service.branch_permutation(*args),
            lambda: service.apply_recovery(*args, source),
            lambda: service.recovery_resources(*args)):
        with pytest.raises(ValueError):
            invoke()


@pytest.mark.parametrize("damage", ("drop_abort", "duplicate_abort", "bad_current", "boolean_coin"))
def test_complete_native_grid_rejects_zero_omissions_and_metadata_mutation(monkeypatch, replay, damage):
    original = service.recorded.repair_branches
    def wrong(law, state):
        rows = [dict(row) for row in original(law, state)]
        abort = next(i for i, row in enumerate(rows) if row["coin"] == 0 and row["record_code"] == 3)
        if damage == "drop_abort":
            rows.pop(abort)
        elif damage == "duplicate_abort":
            rows[0] = dict(rows[abort])
        elif damage == "bad_current":
            rows[0]["current"] = 1
        else:
            rows[0]["coin"] = False
        return rows
    monkeypatch.setattr(service.recorded, "repair_branches", wrong)
    with pytest.raises(ArithmeticError, match="native raw output"):
        service.execute_stopped_history(service.witness_schedules()["motion_then_forbidden_reset"][:1])


@pytest.mark.parametrize("schedule", ((), (("C", (0, 1), False, 0),),
                                      (("C", (0, 1), 0, False),),
                                      (("C", (False, 1), 0, 0),)))
def test_stopped_history_requires_a_strict_complete_continuation_address(schedule):
    with pytest.raises(ValueError):
        service.execute_stopped_history(schedule)
