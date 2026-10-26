"""Retained scientific outputs cannot pass CI on custody checks alone."""
import copy
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from .replay import compare_arrays, compare_report


def _mass_report(word_resources):
    word_resources = {"integrated_processor_norm_upper": 0.0, **word_resources}
    return {"environment": {"python": "3.12.6", "numpy": "2.3.5", "scipy": "1.16.3"},
            "source_bindings": {"producer.py": "same-source-hash"},
            "cases": {"interacting_n32": {
                "elapsed_host_seconds": 1.0, "record_probability": 0.25,
                "resources_per_realized_nonabort_history": {
                    "preparation_source_word": word_resources}}}}


def test_roundoff_zero_classification_keeps_every_pulse_and_quantum_output():
    from .mass_response import history
    from ..observer.source_scalar import compiler
    word = [compiler.pulse("Z", 0., [0]), compiler.pulse("H", 0.4, [0])]
    perturbed = [compiler.pulse("Z", 1e-15, [0]), word[1]]
    blank = np.array([[1.], [0.]], complex)
    before, retained_cost = history.execute_word(word, 1, blank)
    after, current_cost = history.execute_word(perturbed, 1, blank)
    assert retained_cost["zero_angle_requests"] == 1
    assert current_cost["zero_angle_requests"] == 0
    assert retained_cost["native_pulses"] == current_cost["native_pulses"] == 2
    assert retained_cost["proper_code_service_events"] == current_cost["proper_code_service_events"] == 4
    compare_arrays({"complete_output": before}, {"complete_output": after})
    retained, current = _mass_report(retained_cost), _mass_report(current_cost)
    compare_report(retained, current, host_timing=True)
    # The input receipts retain their actual host diagnostic; only temporary
    # comparison copies are normalized.
    assert retained_cost["zero_angle_requests"] == 1
    assert current_cost["zero_angle_requests"] == 0
    compare_report(retained, current)


@pytest.mark.parametrize("field,new_value", [
    ("native_pulses", 3), ("proper_code_service_events", 6),
    ("integrated_processor_norm_upper", 7.0),
    ("service_failure_outcomes_retained", False),
])
def test_zero_diagnostic_exception_preserves_meaningful_resource_checks(field, new_value):
    cost = {"zero_angle_requests": 1, "native_pulses": 2,
            "proper_code_service_events": 4, "integrated_processor_norm_upper": 0.4,
            "service_failure_outcomes_retained": True}
    retained = _mass_report(cost)
    current = copy.deepcopy(retained)
    word = current["cases"]["interacting_n32"]["resources_per_realized_nonabort_history"]["preparation_source_word"]
    word["zero_angle_requests"] = 0
    word[field] = new_value
    with pytest.raises(ValueError, match="replay diagnostic|exact value changed|integrated processor norm"):
        compare_report(retained, current, host_timing=True)


def test_preparation_integrated_norm_is_checked_as_a_host_envelope():
    cost = {"zero_angle_requests": 1, "native_pulses": 2,
            "integrated_processor_norm_upper": 0.4,
            "proper_code_service_events": 4}
    retained = _mass_report(cost)
    current = copy.deepcopy(retained)
    word = current["cases"]["interacting_n32"]["resources_per_realized_nonabort_history"]["preparation_source_word"]
    word["integrated_processor_norm_upper"] = 5.0
    compare_report(retained, current, host_timing=True)
    compare_report(retained, current)
    for bad in (-1.0, 7.0, True, "5"):
        word["integrated_processor_norm_upper"] = bad
        with pytest.raises(ValueError, match="integrated processor norm"):
            compare_report(retained, current, host_timing=True)
    word["integrated_processor_norm_upper"] = float("nan")
    with pytest.raises(ValueError, match="not JSON compliant"):
        compare_report(retained, current, host_timing=True)


def test_nested_continuing_field_preparation_uses_same_validated_diagnostics():
    word = {"native_pulses": 2, "zero_angle_requests": 1,
            "integrated_processor_norm_upper": 0.4,
            "proper_code_service_events": 4}
    retained = {"continuing_field_read": {
        "resources_per_realized_nonabort_history": {
            "preparation_source_word": word}, "physical_response": 0.25}}
    current = copy.deepcopy(retained)
    actual = current["continuing_field_read"]["resources_per_realized_nonabort_history"]["preparation_source_word"]
    actual["zero_angle_requests"] = 0
    actual["integrated_processor_norm_upper"] = 5.0
    compare_report(retained, current)
    actual["proper_code_service_events"] = 6
    with pytest.raises(ValueError, match="exact value changed"):
        compare_report(retained, current)
    actual["proper_code_service_events"] = 4
    current["continuing_field_read"]["physical_response"] = 0.5
    with pytest.raises(ValueError, match="floating replay diagnostic"):
        compare_report(retained, current)


@pytest.mark.parametrize("count", [-1, 3, True, 0.5, "0", None])
def test_invalid_zero_angle_diagnostic_is_not_ignored(count):
    retained = _mass_report({"zero_angle_requests": 1, "native_pulses": 2})
    current = _mass_report({"zero_angle_requests": count, "native_pulses": 2})
    with pytest.raises(ValueError, match="zero-angle diagnostic"):
        compare_report(retained, current, host_timing=True)


def test_missing_zero_count_and_changed_source_binding_still_fail():
    retained = _mass_report({"zero_angle_requests": 1, "native_pulses": 2})
    current = _mass_report({"native_pulses": 2})
    with pytest.raises(ValueError, match="zero-angle diagnostic"):
        compare_report(retained, current, host_timing=True)
    current = _mass_report({"zero_angle_requests": 0, "native_pulses": 2})
    current["source_bindings"]["producer.py"] = "different-source-hash"
    with pytest.raises(ValueError, match="exact value changed"):
        compare_report(retained, current, host_timing=True)


def test_host_runtime_may_change_but_record_physics_may_not():
    original = {"environment": {"python": "3.12.6", "numpy": "2.5.3", "scipy": "1.18.1"},
                "base_release_commit": "same-declared-base",
                "cases": {"interacting": {
        "elapsed_host_seconds": 1.0, "record_probability": 0.25}}}
    replayed = copy.deepcopy(original)
    replayed["cases"]["interacting"]["elapsed_host_seconds"] = 20.0
    replayed["environment"] = {"python": "3.12.12", "numpy": "2.3.5", "scipy": "1.16.3"}
    compare_report(original, replayed, host_timing=True)
    replayed["base_release_commit"] = "different-base"
    with pytest.raises(ValueError, match="exact value changed"):
        compare_report(original, replayed, host_timing=True)
    replayed["base_release_commit"] = original["base_release_commit"]
    replayed["cases"]["interacting"]["record_probability"] = 0.5
    with pytest.raises(ValueError, match="floating replay diagnostic"):
        compare_report(original, replayed, host_timing=True)


def test_missing_environment_provenance_cannot_be_ignored():
    bad = {"environment": {"python": "3.12"}, "cases": {}}
    with pytest.raises(ValueError, match="environment metadata"):
        compare_report(bad, bad, host_timing=True)


def test_fresh_mass_producer_preserves_its_exact_computational_source_closure():
    root = Path(__file__).resolve().parents[2]
    program = ("import json; "
               "from model.common_source.mass_response.run import source_bindings; "
               "print(json.dumps(source_bindings()))")
    child = subprocess.run([sys.executable, "-c", program], cwd=root,
                           capture_output=True, text=True, check=True)
    retained = json.loads((root/"model/common_source/mass_response/results.json").read_text())
    assert json.loads(child.stdout) == retained["source_bindings"]


def test_record_labels_and_small_quantum_amplitudes_are_compared():
    retained = {"records": np.array(["00", "01"]),
                "outputs": np.array([1.0, 0.0], dtype=complex)}
    compare_arrays(retained, retained)
    changed = {**retained, "records": np.array(["00", "10"])}
    with pytest.raises(ValueError, match="records"):
        compare_arrays(retained, changed)
    changed = {**retained, "outputs": np.array([1.0, 1e-5j])}
    with pytest.raises(ValueError, match="outputs"):
        compare_arrays(retained, changed)


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_nonfinite_quantum_output_is_never_equal(value):
    arrays = {"output": np.array([value])}
    with pytest.raises(ValueError, match="output"):
        compare_arrays(arrays, arrays)
