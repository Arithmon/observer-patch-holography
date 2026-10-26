"""Retained scientific outputs cannot pass CI on custody checks alone."""
import copy

import numpy as np
import pytest

from .replay import compare_arrays, compare_report


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
