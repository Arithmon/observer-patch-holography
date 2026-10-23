"""Cross-scale diagnostic: exact transport, isolation, and retained-data replay."""
from decimal import localcontext
from fractions import Fraction as F
import json
from pathlib import Path

import pytest

from model.predictions import shared_parameter as m


@pytest.fixture(scope="module")
def result():
    return m.build()


def test_exact_affine_image_includes_both_endpoints_and_reverses_order():
    lo, hi = F(3, 2), F(7, 4)
    lower, upper = m.tilt_interval(lo, hi)
    assert lower == F(185, 192) and upper == F(31, 32)
    assert upper - lower == (hi - lo) / 48
    for point in (lo, (lo+hi)/2, hi):
        assert lower <= 1-point/48 <= upper


@pytest.mark.parametrize("lo,hi", [(F(2), F(1)), (F(0), F(1)), (F(1), F(48)), (1.0, 2.0)])
def test_invalid_or_float_source_interval_is_refused(lo, hi):
    with pytest.raises(ValueError):
        m.tilt_interval(lo, hi)


def test_exact_tail_enclosure_keeps_endpoint_equality_and_multiplicity():
    rows = [["2", "0.9"], ["3", "0.91"], ["5", "0.92"], ["7", "0.93"]]
    found = m.weighted_positions(rows, F("0.91"), F("0.92"), F("0.915"))
    assert found == {"multiplicity": 17, "forward_tail_weight_interval": [5, 10], "inherited_tail_weight": 5}
    for invalid in ([["1.5", "0.96"]], [["0", "0.96"]], [[1, "0.96"]], []):
        with pytest.raises(ValueError):
            m.weighted_positions(invalid, F("0.91"), F("0.92"), F("0.915"))


def test_predictions_need_no_observations_and_cannot_read_act_or_target(monkeypatch):
    original = Path.read_bytes

    def guarded(path):
        assert path.relative_to(m.ROOT).as_posix() in m.endpoint.SOURCE_PINS
        return original(path)

    monkeypatch.setattr(Path, "read_bytes", guarded)
    monkeypatch.setattr(m.endpoint, "target_fixture", lambda: pytest.fail("target loaded during prediction"))
    prediction = m.frozen_forward_prediction()
    with localcontext() as context:
        context.prec = 4
        assert m.frozen_forward_prediction() == prediction


def test_one_forward_parameter_differs_from_inherited_target_calibration(result):
    prediction = result["prediction"]
    lo, hi = map(F, prediction["P_interval"])
    lower, upper = map(F, prediction["scalar_tilt_interval"])
    assert m.INHERITED_P < lo <= hi
    assert (lower, upper) == (1-hi/48, 1-lo/48)
    assert result["prediction_sha256_before_comparison"] == m.digest(prediction)
    assert upper < F(result["cosmology"]["inherited_scalar_tilt"])
    assert result["coverage"]["one_forward_P_used_in_both_maps"]
    assert result["coverage"]["joint_likelihood"] is None
    assert not result["coverage"]["full_common_physical_construction_established"]


def test_actual_chain_replay_and_retained_receipt(result):
    assert json.loads(m.OUTPUT.read_text()) == result
    cmb = result["cosmology"]
    existing = json.loads((m.ROOT / "code/cosmology/tilt_posterior/results.json").read_text())
    assert cmb["inherited_P"] == existing["prediction"]["P_exact"]
    assert cmb["inherited_scalar_tilt"] == existing["prediction"]["edge_center_ns_exact"]
    assert cmb["multiplicity"] == int(existing["pooled"]["total_multiplicity"])
    old = existing["pooled"]["positions"]["edge_center"]
    assert cmb["inherited_tail_weight"] == int(old["weight_strictly_below"]) + int(old["weight_equal"])
    assert cmb["forward_tail_weight_interval"] == [3114, 3114]
    assert cmb["retained_rows"] == 379175
    assert cmb["forward_and_inherited_have_same_tail_weight"]
    assert result["coupling"]["zero_missing_correction_verdict"] == "inconsistent_with_declared_comparison_band"


def test_mutated_projection_is_refused_before_it_can_improve_tail(tmp_path):
    target = tmp_path / m.PROJECTION
    target.parent.mkdir(parents=True)
    target.write_bytes(b"substituted observations")
    with pytest.raises(ValueError, match="projection hash"):
        m.compare_projection(m.frozen_forward_prediction(), root=tmp_path)


def test_source_pin_failure_precedes_observation_access(tmp_path):
    for path in m.endpoint.SOURCE_PINS:
        copied = tmp_path / path
        copied.parent.mkdir(parents=True, exist_ok=True)
        copied.write_bytes((m.ROOT / path).read_bytes())
    source = tmp_path / m.endpoint.RECEIPT
    source.write_bytes(source.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="Source pin mismatch"):
        m.build(root=tmp_path)
