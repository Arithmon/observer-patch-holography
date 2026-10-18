"""Independent interval, joint-origin, native-law and provenance controls."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path

import pytest

from tools import strict_json
from .attachment import Endpoint, analyse, fit_rate, offset_bounds, predict_seconds, validate
from .experiment import deadline_survival, run
from ..native_clock import clock
from ..primitive_observer import transfer


def point(name, n, t, error=0):
    return Endpoint(name, n, n, F(t)-F(error), F(t)+F(error))


@pytest.fixture(scope="module")
def result():
    return run()


@pytest.fixture(scope="module")
def reduced():
    return clock.lumping(transfer.native_counts())


def test_exact_origin_and_rate_are_recovered_without_validation():
    cal = (point("a", 0, 7), point("b", 20, 17))
    fitted = fit_rate(cal)
    assert fitted.lower == fitted.upper == F(1, 2)
    assert offset_bounds(cal, F(1, 2)) == (F(7), F(7))
    assert predict_seconds(cal, 100, 100) == (F(57), F(57))
    rejected = validate(cal, (point("v", 100, 58),))
    assert rejected["joint_validation"] == "incompatible_with_every_calibration"
    assert rejected["calibration"] == fitted.report()


def test_timestamp_prediction_keeps_offset_rate_correlation():
    cal = (point("a", 0, 0, 1), point("b", 10, 10, 1))
    fitted = fit_rate(cal)
    assert (fitted.lower, fitted.upper) == (F(4, 5), F(6, 5))
    assert predict_seconds(cal, 20, 20) == (F(17), F(23))
    # Treating origin in [-1,1] and rate in [.8,1.2] independently would give
    # [-1+20*.8, 1+20*1.2] = [15,25], which is unnecessarily loose.
    assert predict_seconds(cal, 5, 5) == (F(4), F(6))
    assert predict_seconds(cal, 19, 21) == (F(81, 5), F(121, 5))


def test_separately_compatible_validation_points_can_fail_one_common_clock(result):
    rows = result["joint_validation_control"]
    assert all(p["compatible_with_calibration"] for p in rows["per_endpoint"])
    assert rows["joint_validation"] == "incompatible_with_every_calibration"
    assert rows["calibration"]["seconds_per_attempt_lower"] == "4/5"


def test_count_detection_errors_are_not_silently_zero(result):
    control = result["missed_event_control"]
    assert control["no_missing_events"]["joint_validation"] == "incompatible_with_every_calibration"
    assert control["up_to_five_missing_events"]["joint_validation"] == "compatible_with_some_calibration"
    example = strict_json.load(Path(__file__).with_name("endpoints_example.json"))
    example["endpoints"][2]["count_error"] = None
    checked = analyse(example)
    assert checked["numerical_test"] is None
    assert checked["unresolved_error_endpoints"] == ["v0"]


def test_uncertain_integer_count_relaxation_is_only_a_necessary_condition():
    # Exact endpoints force tau=1 and b=0. A third timestamp 5/2 lies in the
    # real relaxation [2,3] but corresponds to no integer event count.
    cal = (point("a", 0, 0), point("b", 1, 1))
    uncertain = Endpoint("v", 2, 3, F(5, 2), F(5, 2))
    assert fit_rate(cal+(uncertain,)).feasible
    assert not any(fit_rate(cal+(point("v", n, F(5, 2)),)).feasible for n in (2, 3))
    report = validate(cal, (uncertain,))
    assert report["count_bounds_relaxed_to_real_intervals"]
    assert not report["integer_count_realization_certified"]
    assert not report["true_count_order_enforced_with_uncertain_counts"]
    assert "outer feasibility only" in report["compatibility_scope"]


def test_no_clock_scale_from_one_event_and_zero_or_negative_rates_fail():
    assert fit_rate((point("a", 0, 4),)).feasible
    assert not fit_rate((point("a", 0, 4),)).identified
    assert not fit_rate((point("a", 0, 4), point("b", 10, 4))).feasible
    assert not fit_rate((point("a", 0, 4), point("b", 10, 3))).feasible
    assert not fit_rate((point("a", 0, 4), point("b", 0, 5))).feasible
    with pytest.raises(ValueError, match="bound its positive rate"):
        predict_seconds((point("a", 0, 4),), 1, 1)


def test_absolute_time_origin_change_does_not_change_rate():
    source = (point("a", 0, 2, F(1, 1000)), point("b", 1000, 3, F(1, 1000)))
    translated = tuple(Endpoint(p.name, p.count_lo, p.count_hi,
                               p.seconds_lo+100, p.seconds_hi+100) for p in source)
    assert fit_rate(source) == fit_rate(translated)
    assert predict_seconds(translated, 2000, 2000) == tuple(
        x+100 for x in predict_seconds(source, 2000, 2000))


def test_native_deadline_uses_complete_waiting_tail_not_mean_interval(reduced):
    exact = deadline_survival(reduced, (F(1, 1000), F(1, 1000)), F(1, 100))
    _, tail = clock.waiting_law(reduced, 10)
    assert exact["attempts_completed_bounds"] == [10, 10]
    assert exact["probability_first_current_after_deadline"] == [str(tail)]*2
    assert 0 < tail < 1  # replacing the mean wait by a deterministic time loses this tail
    uncertain = deadline_survival(reduced, (F(9, 10000), F(11, 10000)), F(1, 100))
    lo, hi = map(F, uncertain["probability_first_current_after_deadline"])
    assert uncertain["attempts_completed_bounds"] == [9, 11]
    assert lo <= tail <= hi
    for start in ("stationary", "palm"):
        zero = deadline_survival(reduced, (F(1), F(2)), F(0), start=start)
        assert zero["probability_first_current_after_deadline"] == ["1", "1"]


def test_calibration_control_cannot_refit_and_host_accuracy_stays_unresolved(result):
    baseline = result["synthetic_endpoint_fixture"]
    changed = result["synthetic_held_out_rate_change"]
    assert baseline["numerical_test"]["calibration"] == changed["numerical_test"]["calibration"]
    assert baseline["numerical_test"]["joint_validation"] == "compatible_with_some_calibration"
    assert changed["numerical_test"]["joint_validation"] == "incompatible_with_every_calibration"
    actual = result["retained_host_attachment"]
    assert actual["all_native_records_replayed"] and actual["same_record_tape"]
    assert actual["native_attempts_each_run"] == [264532, 264532]
    host = actual["actual_host_measurement_with_unsupplied_accuracy"]
    assert host["numerical_test"] is None
    assert len(host["unresolved_error_endpoints"]) == 5
    assert not host["laboratory_correspondence_established"]
    assert not actual["new_timestamps_acquired"]
    assert actual["zero_error_counterfactual"]["baseline"]["joint_validation"] == "incompatible_with_every_calibration"


def test_provenance_strings_do_not_authenticate_physical_measurement():
    example = strict_json.load(Path(__file__).with_name("endpoints_example.json"))
    example["data_kind"] = "laboratory_measurement"
    example["provenance"] = {key: "unverified external assertion" for key in example["provenance"]}
    report = analyse(example)
    assert report["unresolved_provenance"] == []
    assert not report["provenance_authenticated"]
    assert not report["laboratory_correspondence_established"]


@pytest.mark.parametrize("field,value", [("attempt_count", False), ("count_error", True),
    ("seconds", 1.5), ("seconds_error", "-1"), ("role", "training")])
def test_malformed_counter_and_timestamp_contracts_fail(field, value):
    example = strict_json.load(Path(__file__).with_name("endpoints_example.json"))
    example["endpoints"][1][field] = value
    with pytest.raises(ValueError):
        analyse(example)


def test_disjoint_calibration_and_validation_required():
    cal = (point("a", 0, 0), point("b", 1, 1))
    with pytest.raises(ValueError, match="splits disjoint"):
        validate(cal, (point("a", 2, 2),))


def test_receipt_replays_exactly(result):
    expected = strict_json.load(Path(__file__).with_name("results.json"))
    assert result == expected
