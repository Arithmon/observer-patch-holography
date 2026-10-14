"""Fresh end-to-end audit of record access, prediction order and retained laws."""
from copy import deepcopy
from itertools import product
import json

import numpy as np
import pytest

from . import experiment, recovery, source


@pytest.fixture(scope="module")
def executed():
    """Run once, observing actual calls without replacing their outputs."""
    events, field_calls, inverse_inputs = [], [], []
    original_field = source.field_readout
    original_inverse = recovery.recover
    original_prediction = recovery.predict_chi
    original_target = source.supplied_target

    def forbidden_backend(*args, **kwargs):
        raise AssertionError("inverse or prediction attempted backend/target access")

    def field(dimension, duration, x, y, gamma, eta, **kwargs):
        events.append(("field", dimension, duration, x, y, gamma, eta,
                       kwargs.get("anharmonic", 0.)))
        probabilities, resources = original_field(dimension, duration, x, y, gamma, eta, **kwargs)
        field_calls.append({"dimension": dimension, "duration": duration,
                            "x": x, "y": y, "gamma": gamma, "eta": eta,
                            "anharmonic": kwargs.get("anharmonic", 0.),
                            "raw_probabilities": probabilities.copy(),
                            "resources": deepcopy(resources)})
        return probabilities, resources

    def infer(*args, **kwargs):
        events.append(("inverse",))
        saved = tuple(np.asarray(value).copy() for value in args)
        with pytest.MonkeyPatch.context() as block:
            block.setattr(source, "supplied_target", forbidden_backend)
            block.setattr(source, "modes", forbidden_backend)
            block.setattr(source.c, "execute", forbidden_backend)
            result = original_inverse(*args, **kwargs)
        assert all(np.array_equal(before, after) for before, after in zip(saved, args, strict=True))
        inverse_inputs.append(saved)
        return result

    def predict(stiffness, time, x, y, **kwargs):
        events.append(("prediction", time, x, y, kwargs["gamma"], kwargs["eta"]))
        with pytest.MonkeyPatch.context() as block:
            block.setattr(source, "supplied_target", forbidden_backend)
            block.setattr(source, "modes", forbidden_backend)
            block.setattr(source.c, "execute", forbidden_backend)
            return original_prediction(stiffness, time, x, y, **kwargs)

    def target():
        events.append(("audit_target",))
        return original_target()

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(source, "field_readout", field)
        patch.setattr(recovery, "recover", infer)
        patch.setattr(recovery, "predict_chi", predict)
        patch.setattr(source, "supplied_target", target)
        result = experiment.run()
    return result, events, field_calls, inverse_inputs


def test_executed_records_are_the_only_inverse_inputs_and_are_not_rewritten(executed):
    result, _, _, inverse_inputs = executed
    assert len(inverse_inputs) == len(result["levels"]) == 2
    for level, arguments in zip(result["levels"], inverse_inputs, strict=True):
        first, second, clock_first, clock_second = arguments
        for argument, table in zip((first, second), level["calibration_records"], strict=True):
            assert np.array_equal(argument, np.asarray(table["probabilities"])[..., 0])
        for argument, table in zip((clock_first, clock_second),
                                   level["calibration_clock_probabilities"], strict=True):
            assert np.array_equal(argument, np.asarray(table)[:, 0])


def test_each_withheld_prediction_precedes_its_source_read_and_target_audit(executed):
    result, events, _, _ = executed
    cursor = 0
    for level in result["levels"]:
        dimension = level["occupation_per_normal_mode"]
        calibration = events[cursor:cursor+256]
        assert len(calibration) == 256
        assert all(event[0] == "field" and event[1] == dimension for event in calibration)
        assert {event[2] for event in calibration} == set(experiment.CALIBRATION_DURATIONS)
        cursor += 256
        assert events[cursor] == ("inverse",)
        cursor += 1
        for setting, retained in zip(experiment.WITHHELD, level["withheld_predictions"], strict=True):
            duration, x, y, gamma, eta = setting
            predicted, observed = events[cursor:cursor+2]
            assert predicted[0] == "prediction"
            assert predicted[1] == pytest.approx(duration, abs=2e-12)
            assert predicted[2:] == (x, y, gamma, eta)
            assert observed == ("field", dimension, duration, x, y, gamma, eta, 0.)
            assert retained["clock_reading"] == predicted[1]
            assert duration not in experiment.CALIBRATION_DURATIONS
            cursor += 2
        assert events[cursor] == ("audit_target",)
        cursor += 1
    assert all(event[0] != "audit_target" for event in events[cursor:])


def test_raw_complete_field_laws_and_per_trial_resource_counts_survive_receipt(executed):
    result, _, calls, _ = executed
    position = 0
    for level in result["levels"]:
        for retained in level["calibration_records"]:
            block = calls[position:position+128]
            assert len(block) == 128
            raw = np.array([call["raw_probabilities"] for call in block]).reshape(2, 8, 8, 2, 2)
            assert np.array_equal(raw, retained["probabilities"])
            assert raw.min() >= -2e-12
            np.testing.assert_allclose(raw.sum(axis=-1), 1., atol=2e-10, rtol=0)
            assert retained["fresh_phase_trials"] == 256
            for key in source.COUNT_KEYS:
                expected = sum(call["resources"]["resources_for_two_fresh_phase_trials"][key]
                               for call in block)
                assert retained["resources"][key] == expected
            position += 128
        for retained in level["withheld_predictions"]:
            call = calls[position]
            assert np.array_equal(call["raw_probabilities"], retained["executed_probabilities"])
            assert retained["execution"] == call["resources"]
            assert retained["execution"]["complete_pointer_outcomes"] == 2
            resources = retained["execution"]["resources_for_two_fresh_phase_trials"]
            assert resources["native_pulses"] > 0
            assert resources["proper_code_service_events"] == 2*resources["native_pulses"]
            np.testing.assert_allclose(call["raw_probabilities"].sum(axis=-1), 1., atol=2e-10, rtol=0)
            position += 1
    assert position == 544
    assert len(calls) == position+3  # separately labelled changed-action controls


def test_recovered_generator_graph_and_heldout_probabilities_pass_independent_audit(executed):
    result, _, _, _ = executed
    positions = np.asarray(list(product((0, 1), repeat=3)))[list(source.PERMUTATION)]
    squared_distances = np.sum((positions[:, None]-positions[None, :])**2, axis=-1)
    adjacency = (squared_distances == 1).astype(int)
    expected = source.MASS**2*np.eye(8)+source.KAPPA*(3*np.eye(8)-adjacency)
    for level in result["levels"]:
        recovered = level["recovered"]
        assert all(level["checks"].values()) and level["all_pass"]
        np.testing.assert_allclose(recovered["stiffness"], expected, atol=2e-8, rtol=0)
        np.testing.assert_array_equal(recovered["adjacency"], adjacency)
        np.testing.assert_array_equal(recovered["graph_distances"], squared_distances)
        np.testing.assert_allclose(recovered["gram"], (positions-.5)@(positions-.5).T, atol=1e-13, rtol=0)
        assert recovered["mass"] == pytest.approx(source.MASS, abs=2e-8)
        assert recovered["coupling"] == pytest.approx(source.KAPPA, abs=2e-8)
        assert recovered["spatial_dimension"] == 3
        assert level["maximum_withheld_probability_error"] < 2e-9
        assert len(level["withheld_predictions"]) == len(experiment.WITHHELD) == 16
        for row in level["withheld_predictions"]:
            difference = np.max(abs(np.asarray(row["predicted_probabilities"])
                                    - np.asarray(row["executed_probabilities"])))
            assert row["maximum_probability_error"] == difference


def test_controls_and_native_geometry_preserve_the_claim_boundaries(executed):
    result, _, calls, _ = executed
    controls = result["controls"]
    assert all(controls["checks"].values())
    periodic = controls["periodic_clock"]
    assert periodic["unwrapped_window_is_a_protocol_premise"]
    np.testing.assert_allclose(periodic["first_probabilities"],
                               periodic["one_cycle_later_probabilities"], atol=2e-12, rtol=0)
    loop = controls["weyl_loop"]
    forward, reverse = np.asarray(loop["forward_probabilities"]), np.asarray(loop["reverse_probabilities"])
    assert np.max(abs(forward-reverse)) > .1
    assert loop["commuting_control_probability_margin"] > .05
    assert not loop["untrusted_device_quantumness_certified"]
    changed = controls["changed_action"]
    assert not changed["is_local_site_phi4"]
    assert not changed["free_occupation_certificate_applies"]
    for row, call in zip(changed["rows"], calls[-3:], strict=True):
        assert call["anharmonic"] == .25
        assert np.array_equal(row["executed_anharmonic_probabilities"], call["raw_probabilities"])
    assert max(row["maximum_probability_disagreement"] for row in changed["rows"]) > 1e-3
    native = result["native_geometry"]
    assert native["all_pass"] and native["recovery"]["rank_three_unit_gram_recovered"]
    assert native["source_execution"]["native_branch_calls_main"] == 720
    assert not result["scope"]["native_carrier_identified_with_field_grid"]
    assert not result["scope"]["native_repair_selects_field_action"]
    assert not result["scope"]["interacting_relativistic_QFT_recovered"]
    assert not result["scope"]["physical_completion_or_new_empirical_prediction"]
    assert result["all_pass"]
    json.dumps(result, allow_nan=False)
