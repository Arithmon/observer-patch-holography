"""Independent expanded-sample and original-decimal posterior controls."""
from fractions import Fraction as F
import importlib.util
import json
import math
from pathlib import Path
import statistics

import pytest

PATH = Path(__file__).with_name("audit.py")
SPEC = importlib.util.spec_from_file_location("oph_tilt_audit", PATH)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_weighted_results_match_independently_expanded_sample():
    rows = [(F(3), F(91, 100)), (F(1), F(99, 100)), (F(4), F(101, 100)),
            (F(2), F(97, 100))]
    expanded = sorted(x for w, x in rows for _ in range(int(w)))
    target = F(99, 100)
    result = mod.summary(rows, {"control": target})
    assert result["ns_mean"] == float(statistics.mean(expanded))
    assert result["ns_posterior_sd"] == math.sqrt(float(statistics.pvariance(expanded)))
    for q in mod.QUANTILES:
        expected = expanded[max(0, math.ceil(q * len(expanded)) - 1)]
        assert result["step_cdf_quantiles"][f"{float(q):.3f}"] == float(expected)
    observed = result["positions"]["control"]
    assert F(observed["posterior_mass_le_target_exact"]) == F(
        sum(x <= target for x in expanded), len(expanded))
    assert observed["weight_equal"] == "1"


def test_burn_in_is_per_chain_stored_rows_not_multiplicity():
    rows = mod.parse_rows([["1000", "0.9"], ["1", "0.95"],
                           ["1", "1.0"], ["1", "1.05"], ["1", "1.1"]])
    kept, dropped = mod.discard_burn_in(rows)
    assert dropped == 2
    assert kept == rows[2:]
    assert mod.summary(kept, {})["total_multiplicity"] == "3"


def test_original_decimal_tail_boundary_survives_float_rounding():
    rows = mod.parse_rows([["2", "0.99999999999999999999"],
                           ["3", "1"], ["5", "1.00000000000000000001"]])
    assert len({float(x) for _, x in rows}) == 1
    result = mod.summary(rows, {"boundary": F(1)})["positions"]["boundary"]
    assert result["weight_strictly_below"] == "2"
    assert result["weight_equal"] == "3"
    assert result["weight_strictly_above"] == "5"
    assert result["posterior_mass_le_target_exact"] == "1/2"


@pytest.mark.parametrize("bad", [[], [["0", "0.96"]], [["-1", "0.96"]],
                                 [["0.5", "0.96"]], [["1", "NaN"]],
                                 [["1", "inf"]], [["1", "3"]],
                                 [["1/2", "0.96"]]])
def test_invalid_or_non_cobaya_samples_fail(bad):
    with pytest.raises(ValueError):
        mod.parse_rows(bad)


def test_prediction_is_frozen_rational_not_refitted():
    assert mod.EDGE_TILT == F(463690317905960407, 480000000000000000)
    assert mod.predictions()["edge_center"] == mod.EDGE_TILT


def test_ordering_keeps_excursions_separate_from_multiplicity():
    rows = [(F(100), F(9, 10)), (F(2), F(9, 10)), (F(1), F(11, 10)),
            (F(4), F(9, 10)), (F(1), F(11, 10))]
    result = mod.ordering_diagnostics(rows, {"boundary": F(1)})
    assert result["ns_changes_between_stored_rows"] == 3
    assert result["contiguous_tail_excursions"]["boundary"] == 2
    assert result["effective_sample_size_estimator"] is False


def test_unequal_chain_pooling_and_omission_match_expanded_observations():
    inputs = [
        [["30", "0.7"], ["2", "0.95"], ["1", "1.02"]],
        [["10", "0.8"], ["1", "1.01"], ["7", "0.94"], ["2", "1.02"]],
        [["9", "0.8"], ["1", "0.9"], ["1", "1.01"], ["2", "0.96"], ["4", "1.02"]],
        [["20", "0.7"], ["3", "0.93"], ["1", "1.03"]],
    ]
    payload = {"schema": "oph.act_dr6.ns_projection.v1",
               "archive_sha256": mod.SOURCE_ARCHIVE_SHA256,
               "burn_in_fraction": "0.5",
               "chains": [{"member": f"synthetic-{i}", "source_sha256": "0" * 64,
                           "rows": rows} for i, rows in enumerate(inputs)]}
    result = mod.analyze(payload)
    expanded = [[F(x) for w, x in rows[len(rows)//2:] for _ in range(int(w))]
                for rows in inputs]
    for omitted, observed in enumerate(result["leave_one_chain_out"]):
        samples = [x for i, values in enumerate(expanded) if i != omitted for x in values]
        for key, target in mod.predictions().items():
            expected = F(sum(x <= target for x in samples), len(samples))
            assert observed["posterior_mass_le_target"][key] == float(expected)
    all_samples = [x for samples in expanded for x in samples]
    assert result["pooled"]["ns_mean"] == float(statistics.mean(all_samples))


def test_retained_actual_data_replay():
    result = mod.build()
    assert result == json.loads(mod.RECEIPT.read_text())
    pooled = result["pooled"]
    assert pooled["stored_rows"] == 379175
    assert pooled["total_multiplicity"] == "508985"
    # Separate archive-based float parser produced these integer event totals.
    edge = pooled["positions"]["edge_center"]
    assert edge["weight_strictly_below"] == "3114"
    assert edge["weight_equal"] == "0"
    assert F(edge["posterior_mass_le_target_exact"]) == F(3114, 508985)
    assert result["total_tail_excursions_across_chains"]["edge_center"] == 85
    assert result["analysis"]["tail_mass_is_p_value"] is False
    assert result["analysis"]["FZ15_prospective_verdict"] is None


def test_rewritten_projection_is_rejected_before_analysis(tmp_path):
    altered = tmp_path / "altered.gz"
    altered.write_bytes(b"arbitrary replacement archive")
    with pytest.raises(ValueError, match="projection hash differs"):
        mod.build(altered)
