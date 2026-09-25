"""Fast contract checks, with no CAMB evolution or network access."""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent


def load_module(filename, name):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


forward = load_module("run.py", "codex_boltzmann_forward_contract")
audit = load_module("interpolation_audit.py", "codex_boltzmann_interpolation_contract")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_comparison_column_order_and_finite_te_at_zero():
    # Independent hand calculation, including a TE zero and negative TE.
    reference = np.array([[30, 4, 0, 9], [31, 16, -6, 25]], dtype=float)
    predicted = np.array([[30, 4.4, .6, 9.9, 999, 999, 999],
                          [31, 14.4, -8, 20, 888, 888, 888]])
    result = forward.comparison(predicted, reference)
    assert result["TT_fractional_rms"] == pytest.approx(.1)
    assert result["TE_auto_normalized_rms"] == pytest.approx(.1)
    assert result["EE_fractional_rms"] == pytest.approx(np.sqrt(.025))
    assert result["TE_auto_normalized_max_abs"] == pytest.approx(.1)
    assert result["EE_fractional_max_abs"] == pytest.approx(.2)
    assert result["multipoles"] == 2


@pytest.mark.parametrize("drop", [1, 2])
def test_comparison_rejects_missing_interior_or_final_multipole(drop):
    reference = np.array([[30, 4, 0, 9], [31, 4, 0, 9], [32, 4, 0, 9]], dtype=float)
    with pytest.raises(ValueError, match="all reference multipoles"):
        forward.comparison(np.delete(reference, drop, axis=0), reference)


def test_comparison_selects_declared_range_and_rejects_bad_auto_power():
    reference = np.array([[29, 4, 0, 9], [30, 4, 0, 9], [31, 4, 0, 9]], dtype=float)
    predicted = reference.copy()
    predicted[0, 1] = 4e8
    result = forward.comparison(predicted, reference, ell_min=30, ell_max=30)
    assert result["TT_fractional_rms"] == 0
    assert result["multipoles"] == 1
    reference[1, 3] = 0
    with pytest.raises(ValueError, match="auto-power"):
        forward.comparison(predicted, reference)


def test_retained_native_fixture_uses_dell_and_expected_column_mapping():
    spec = json.loads((HERE / "spec.json").read_text())
    reference = np.loadtxt(HERE / spec["theory"])[:, :4]
    predicted = np.loadtxt(HERE / "results/native_planck.txt")
    expected = json.loads((HERE / "results/native_planck_receipt.json").read_text())
    actual = forward.comparison(predicted, reference)
    for key, value in expected["raw_official_theory_comparison"].items():
        assert actual[key] == pytest.approx(value, rel=1e-7, abs=1e-12)
    assert actual["TT_fractional_rms"] < .005
    assert actual["EE_fractional_rms"] < .005
    assert actual["TE_auto_normalized_rms"] < .005
    # A producer accidentally exporting C_l instead of D_l must be detected.
    wrong = predicted.copy()
    ell = wrong[:, 0]
    wrong[:, 1:] *= (2 * np.pi / (ell * (ell + 1)))[:, None]
    assert forward.comparison(wrong, reference)["TT_fractional_rms"] > .99
    swapped = predicted.copy()
    swapped[:, [2, 3]] = swapped[:, [3, 2]]
    assert forward.comparison(swapped, reference)["EE_fractional_rms"] > .1


def test_retained_splined_powerlaw_agrees_with_native_forward_fixture():
    native = np.loadtxt(HERE / "results/native_planck.txt")
    splined = np.loadtxt(HERE / "results/planck_powerlaw.txt")
    result = forward.comparison(splined, native[:, :4])
    assert result["TT_fractional_rms"] < 1e-5
    assert result["EE_fractional_rms"] < 1e-5
    assert result["TE_auto_normalized_rms"] < 1e-5


@pytest.mark.parametrize("floor", [1e-50, 1e-40])
def test_actual_camb_regularized_uv_interpolation_is_positive(floor):
    camb = pytest.importorskip("camb")
    expected = json.loads((HERE / "spec.json").read_text())["camb_version"]
    if camb.__version__ != expected:
        pytest.skip(f"Pinned interpolation audit uses CAMB {expected}")
    result = audit.probe_case(7200, floor, dense_samples=10001)
    assert result["negative_count"] == 0
    assert result["min_spline"] > 0
    assert result["core_relative_interpolation_error_max"] < 1e-6


def test_interpolation_receipt_retains_failure_and_regularization_controls():
    receipt = json.loads((HERE / "interpolation_audit.json").read_text())
    for name, expected in receipt["inputs_sha256"].items():
        assert sha(HERE / name) == expected
    cases = {(r["N_min"], r["floor"], r["regularization"]): r for r in receipt["cases"]}
    assert cases[1800, 1e-300, "max"]["negative_count"] > 0
    assert cases[7200, 1e-300, "max"]["negative_count"] > 0
    for floor in (1e-50, 1e-40):
        assert cases[7200, floor, "max"]["negative_count"] == 0
        assert cases[7200, floor, "max"]["min_spline"] > 0


def test_retained_forward_receipt_inputs_and_artifacts_match_hashes():
    spec = json.loads((HERE / "spec.json").read_text())
    receipts = sorted((HERE / "results").glob("*_receipt.json"))
    assert len(receipts) >= 5
    inputs = {
        "run.py": [HERE / "run.py", HERE / "initial_attempt/run.py"],
        "spec.json": [HERE / "spec.json"],
        "primordial/sources.py": [HERE.parent / "primordial/sources.py"],
        "primordial/spec.json": [HERE.parent / "primordial/spec.json"],
        "planck_2018.ini": [HERE / spec["ini"]],
        "official_theory": [HERE / spec["theory"]],
    }
    for path in receipts:
        receipt = json.loads(path.read_text())
        for key, expected in receipt["inputs_sha256"].items():
            assert any(p.exists() and sha(p) == expected for p in inputs[key]), (path, key)
        case = receipt["case"]
        assert sha(path.parent / f"{case}.txt") == receipt["spectrum_sha256"]
        assert sha(path.parent / f"{case}_params.txt") == receipt["parameters_sha256"]


def test_control_receipts_preserve_inputs_lensing_mode_and_covariance():
    receipts = sorted((HERE / "control_results").glob("*_receipt.json"))
    assert receipts, "The declared controls must retain their receipts"
    for path in receipts:
        receipt = json.loads(path.read_text())
        for name, expected in receipt["inputs_sha256"].items():
            candidates = [HERE / name]
            if name == "run.py":
                candidates.append(HERE / "initial_attempt/run.py")
            assert any(p.exists() and sha(p) == expected for p in candidates), (path, name)
        assert receipt["interpolant_positivity_probe"]["negative_points"] == 0
        assert receipt["interpolant_positivity_probe"]["min_power"] > 0
        if receipt["status"] == "CAMB_failure":
            assert receipt["error"]
            continue
        assert receipt["status"] == "computed"
        name = receipt["name"]
        spectrum = path.parent / f"{name}.txt"
        params = path.parent / f"{name}_params.txt"
        assert sha(spectrum) == receipt["spectrum_sha256"]
        assert sha(params) == receipt["parameters_sha256"]
        parameter_text = params.read_text()
        assert " DoLensing = True\n" in parameter_text
        assert " WantTensors = False\n" in parameter_text
        if not receipt["nonlinear_matter_lensing"]:
            assert " NonLinear = NonLinear_none\n" in parameter_text
        table = np.loadtxt(spectrum)
        assert np.isfinite(table).all()
        assert (table[:, [1, 3, 4, 6]] > 0).all()
        for tt, te, ee in [(1, 2, 3), (4, 5, 6)]:
            assert (table[:, te]**2 <= table[:, tt] * table[:, ee] * (1 + 1e-8)).all()
