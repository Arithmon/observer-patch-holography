"""Check the finite witnesses and the fixed-floor counterexample separately."""
from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "kernel_admissibility_diagnostic",
    HERE / "derive_family_transport_kernel_admissibility_freedom_theorem.py",
)
diagnostic = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(diagnostic)


def test_positive_parameter_algebra_does_not_imply_floor_acceptance():
    floor = Fraction.from_float(diagnostic.GAP_FLOOR)
    # Exact arithmetic exposes the obstruction before matrix roundoff enters.
    for ratio in (Fraction(1, 10), Fraction(1), Fraction(10)):
        for span in (floor / 2, floor, 2 * floor):
            lower = span * ratio / (1 + ratio)
            upper = span / (1 + ratio)
            assert lower > 0 and upper > 0
            assert lower / upper == ratio
            assert lower + upper == span
            assert min(lower, upper) <= floor
    # Equality also fails because the battery uses a strict inequality.
    assert not (floor > floor)


def test_small_span_is_not_reported_as_passing():
    witness = diagnostic.build_witness(1.0, 1e-12)
    assert not witness["certificates"]["simple_centered_spectrum"]
    assert not witness["all_certificates_pass"]


def test_shrink_can_pass_defect_while_failing_overlap_floor():
    witness = diagnostic.build_witness(1.0, 1e-8)
    assert witness["certificates"]["simple_centered_spectrum"]
    assert witness["certificates"]["riesz_margin_defect_below_half_gap"]
    assert not witness["certificates"]["edge_amplitudes_above_floor"]
    assert not witness["all_certificates_pass"]


@pytest.mark.parametrize("field,value", [("GAP_FLOOR", 100.),
                                          ("AMPLITUDE_FLOOR", .99)])
def test_mutated_floor_rejects_the_grid(monkeypatch, field, value):
    monkeypatch.setattr(diagnostic, field, value)
    with pytest.raises(AssertionError, match="failed the certificate battery"):
        diagnostic.build()


@pytest.mark.parametrize("ratio,span", [(0., 1.), (-1., 1.), (1., 0.),
                                       (1., -1.), (float("nan"), 1.),
                                       (1., float("inf"))])
def test_invalid_parameters_rejected(ratio, span):
    with pytest.raises(ValueError, match="positive finite"):
        diagnostic.build_witness(ratio, span)


def test_committed_grid_remains_a_finite_numerical_diagnostic():
    report = diagnostic.build()
    recorded = json.loads(diagnostic.DEFAULT_OUT.read_text(encoding="utf-8"))
    for data in (report, recorded):
        assert data["proof_status"] == "finite_grid_numerical_witnesses_only"
        assert data["row_class"] == "numerical_diagnostic"
        assert data["proof_kind"] == "finite_binary64_witness_evaluation"
        assert "theorem_statement" not in data
        assert not data["guards"]["universal_tolerance_battery_surjectivity"]
        assert not data["guards"]["physical_no_go"]
        assert not data["guards"]["grid_is_source_selected"]
        assert data["guards"]["kernel_template_consumed_for_context"]
        assert not data["witness_grid"]["includes_on_disk_operating_point"]
        assert len(data["witnesses"]) == 12
        assert data["numeric_battery"]["strict_gap_floor"] == 1e-12
        assert data["numeric_battery"]["strict_squared_overlap_floor"] == 1e-12
        obstruction = data["numeric_battery"]["exact_span_obstruction"]
        assert Fraction(obstruction["gap_floor"]) == Fraction.from_float(
            data["numeric_battery"]["strict_gap_floor"])
        assert Fraction(obstruction["excluded_span_upper_inclusive"]) == (
            2 * Fraction(obstruction["gap_floor"]))
    for fresh, old in zip(report["witnesses"], recorded["witnesses"], strict=True):
        assert fresh["target"] == old["target"]
        assert fresh["certificates"] == old["certificates"]
        assert fresh["all_certificates_pass"]
        for key in ("r", "s", "rho_ord", "x2"):
            assert fresh["emitted"][key] == pytest.approx(
                old["emitted"][key], abs=1e-9, rel=1e-12)


def test_context_is_recomputed_and_not_mistaken_for_an_exact_grid_point():
    kernel = json.loads(diagnostic.KERNEL_PATH.read_text(encoding="utf-8"))
    level = max(kernel["refinements"], key=lambda row: row["level"])
    matrix = (np.asarray(level["hermitian_descendant"]["real"])
              + 1j * np.asarray(level["hermitian_descendant"]["imag"]))
    eigenvalues = np.linalg.eigvalsh(matrix)
    ratio = np.diff(eigenvalues)[0] / np.diff(eigenvalues)[1]
    span = eigenvalues[-1] - eigenvalues[0]
    context = diagnostic.build()["on_disk_kernel_context"]
    assert context["raw_gap_ratio_r"] == pytest.approx(ratio, abs=1e-13)
    assert context["spectral_span"] == pytest.approx(span, abs=1e-13)
    assert context["nearest_r_grid_distance"] > 1e-6
    assert context["nearest_s_grid_distance"] > 1e-8
