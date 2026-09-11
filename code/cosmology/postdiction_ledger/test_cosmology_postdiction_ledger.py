"""Tests for the cosmology postdiction ledger producer, receipt, and verifier."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))

import build_cosmology_postdiction_ledger as ledger  # noqa: E402

PRODUCER = HERE / "build_cosmology_postdiction_ledger.py"
VERIFIER = HERE / "verify_cosmology_postdiction_ledger_independent.py"
RECEIPT = HERE / "runtime" / "cosmology_postdiction_ledger.json"
MARKDOWN = HERE / "COSMOLOGY_POSTDICTION_LEDGER.md"
PUBLIC = HERE / "public_inputs.json"
CORPUS = HERE / "corpus_inputs.json"


def run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *args], capture_output=True, text=True, cwd=str(ROOT), check=False
    )


def committed_receipt() -> dict:
    return json.loads(RECEIPT.read_text(encoding="utf-8"))


def rows(receipt: dict) -> dict[str, dict]:
    return {row["row_id"]: row for section in receipt["sections"] for row in section["rows"]}


def test_rebuild_into_temp_dir_matches_committed_bytes(tmp_path: Path) -> None:
    result = run(str(PRODUCER), "--output-dir", str(tmp_path))
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "runtime" / "cosmology_postdiction_ledger.json").read_bytes() == RECEIPT.read_bytes()
    assert (tmp_path / "COSMOLOGY_POSTDICTION_LEDGER.md").read_bytes() == MARKDOWN.read_bytes()
    receipt_bytes, markdown = ledger.build_outputs()
    assert receipt_bytes == RECEIPT.read_bytes()
    assert markdown == MARKDOWN.read_bytes()
    check = run(str(PRODUCER), "--check")
    assert check.returncode == 0, check.stdout + check.stderr


def test_independent_verifier_passes_and_fails_closed(tmp_path: Path) -> None:
    result = run(str(VERIFIER))
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.startswith("PASS:")

    tampered = committed_receipt()
    row = rows(tampered)["ns_edge_center_P_C_planck2018_vi_ns"]
    row["sigma_distance"] = -row["sigma_distance"]
    tampered_path = tmp_path / "tampered.json"
    tampered_path.write_bytes(ledger.render_bytes(tampered))
    failed = run(str(VERIFIER), "--receipt", str(tampered_path))
    assert failed.returncode == 1
    assert failed.stdout.startswith("FAIL:")

    reclassified = committed_receipt()
    rows(reclassified)["tensor_zero_bk18_r_0p05"]["class"] = "shared_baseline"
    reclassified_path = tmp_path / "reclassified.json"
    reclassified_path.write_bytes(ledger.render_bytes(reclassified))
    assert run(str(VERIFIER), "--receipt", str(reclassified_path)).returncode == 1


def test_public_value_mutation_changes_digest_and_sigma() -> None:
    public = json.loads(PUBLIC.read_text(encoding="utf-8"))
    public["items"]["planck2018_vi_ns"]["value"] = 0.9749
    mutated = ledger.build_receipt(
        json.dumps(public, indent=2).encode("utf-8"), CORPUS.read_bytes()
    )
    committed = committed_receipt()
    assert mutated["rows_sha256"] != committed["rows_sha256"]
    original = rows(committed)["ns_edge_center_P_C_planck2018_vi_ns"]
    changed = rows(mutated)["ns_edge_center_P_C_planck2018_vi_ns"]
    assert changed["sigma_distance"] != original["sigma_distance"]
    ns = Fraction(original["oph_value_exact_rational"])
    assert changed["sigma_distance"] == float((ns - Fraction("0.9749")) / Fraction("0.0042"))
    assert changed["verdict"] == "tension"
    assert mutated["pins"][f"{ledger.PACKAGE}/public_inputs.json"] != committed["pins"][f"{ledger.PACKAGE}/public_inputs.json"]


def test_sigma_sign_convention_is_theory_minus_measurement() -> None:
    assert ledger.signed_sigma(Fraction(1), Fraction(0), Fraction(1)) == 1
    assert ledger.signed_sigma(Fraction(0), Fraction(1), Fraction(2)) == Fraction(-1, 2)
    table = rows(committed_receipt())
    planck = table["ns_edge_center_P_C_planck2018_vi_ns"]
    assert planck["oph_value"] > planck["measurement"]["value"] and planck["sigma_distance"] > 0
    act = table["ns_edge_center_P_C_act_dr6_p_act_lb_ns"]
    assert act["oph_value"] < act["measurement"]["value"] and act["sigma_distance"] < 0
    for row in table.values():
        if row["sigma_distance_kind"] == "signed_one_dimensional":
            difference = row["oph_value"] - row["measurement"]["value"]
            assert (difference > 0) == (row["sigma_distance"] > 0)


def test_upper_bound_rows_are_consistent_only_below_the_bound() -> None:
    assert ledger.verdict_from_bound(Fraction(0), Fraction(36, 1000)) == "consistent"
    assert ledger.verdict_from_bound(Fraction(36, 1000), Fraction(36, 1000)) == "exceeds_stated_confidence_bound"
    assert ledger.verdict_from_bound(Fraction(1, 10), Fraction(36, 1000)) == "exceeds_stated_confidence_bound"
    bound_rows = [
        row
        for row in rows(committed_receipt()).values()
        if row["measurement"]["kind"] == "upper_bound" and row["class"] != "not_evaluable"
    ]
    assert len(bound_rows) == 7
    for row in bound_rows:
        assert row["sigma_distance"] is None
        below = row["oph_value"] < row["measurement"]["upper_bound"]
        assert row["verdict"] == ("consistent" if below else "exceeds_stated_confidence_bound")


def test_verdict_rule_thresholds() -> None:
    assert ledger.verdict_from_sigma(Fraction(1999, 1000)) == "consistent"
    assert ledger.verdict_from_sigma(Fraction(-1999, 1000)) == "consistent"
    assert ledger.verdict_from_sigma(Fraction(2)) == "tension"
    assert ledger.verdict_from_sigma(-2.5) == "tension"
    assert ledger.verdict_from_sigma(Fraction(2999, 1000)) == "tension"
    assert ledger.verdict_from_sigma(3.0) == "exceeds_three_sigma_diagnostic"
    for row in rows(committed_receipt()).values():
        if row["sigma_distance"] is not None:
            assert row["verdict"] == ledger.verdict_from_sigma(row["sigma_distance"])


def test_each_row_has_exactly_one_allowed_class_and_the_required_fields() -> None:
    receipt = committed_receipt()
    allowed = set(ledger.ALLOWED_CLASSES)
    assert set(receipt["allowed_classes"]) == allowed
    for row in rows(receipt).values():
        assert set(ledger.ROW_FIELDS) <= set(row)
        assert set(row) <= set(ledger.ROW_FIELDS) | {"extras"}
        class_keys = [key for key in row if key == "class"]
        assert class_keys == ["class"]
        assert isinstance(row["class"], str) and row["class"] in allowed
        assert row["verdict"] in ledger.ALLOWED_VERDICTS
        assert isinstance(row["discriminates_from_baseline"], bool)
        assert row["prospective_data"]


def test_receipt_contract_and_expected_verdicts() -> None:
    receipt = committed_receipt()
    assert receipt["schema"] == "oph.cosmology.postdiction_ledger.v1"
    assert receipt["promotes_nothing"] is True
    assert receipt["frozen_prediction_rows"] == 0
    assert receipt["classification"].startswith("retrospective postdiction ledger")
    text = RECEIPT.read_text(encoding="utf-8")
    assert "\r" not in text and text.endswith("\n") and text.isascii()
    assert not re.search(r"/Users/|/home/|[A-Za-z]:\\\\Users", text)
    assert not re.search(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}", text)
    table = rows(receipt)
    assert receipt["row_count"] == len(table) == 53
    expected = {
        "ns_edge_center_P_C_planck2018_vi_ns": "consistent",
        "ns_edge_center_P_C_act_dr6_p_act_lb_ns": "tension",
        "ns_edge_center_P_fwd_act_dr6_p_act_lb_ns": "tension",
        "ns_clock_branch_P_C_planck2018_vi_ns": "consistent",
        "running_zero_planck2018_running": "consistent",
        "tensor_zero_bk18_r_0p05": "consistent",
        "isocurvature_zero_planck2018_x_cdi_general_100beta_klow": "consistent",
        "curvature_zero_desi_dr2_cmb_omegak": "tension",
        "w0wa_fixed_capacity_desi_dr2_bao_cmb": "tension",
        "w0wa_fixed_capacity_desi_dr2_bao_cmb_pantheonplus": "tension",
        "w0wa_fixed_capacity_desi_dr2_bao_cmb_union3": "exceeds_three_sigma_diagnostic",
        "w0wa_fixed_capacity_desi_dr2_bao_cmb_desy5": "exceeds_three_sigma_diagnostic",
        "lambda_capacity_candidate_a_planck2018_central": "not_evaluable",
        "lambda_lp2_capacity_candidate_a_desi_dr2_display": "tension",
        "lambda_lp2_capacity_candidate_b_desi_dr2_display": "tension",
        "a0_sparc_full_rar_upsilon_disk_0p5": "not_evaluable",
        "btfr_exponent_four_sparc_ols": "tension",
        "ringdown_integer_k_comb_fz14": "not_evaluable",
        "sum_mnu_desi_dr2_not_evaluable": "not_evaluable",
    }
    for row_id, verdict in expected.items():
        assert table[row_id]["verdict"] == verdict, row_id
    discriminating = {row_id for row_id, row in table.items() if row["discriminates_from_baseline"]}
    assert all(
        row_id.startswith(("ns_", "running_zero_", "tensor_zero_", "w0wa_")) for row_id in discriminating
    )
    assert table["a0_sparc_full_rar_upsilon_disk_0p5"]["extras"]["xi_a0_over_a_dS"] == pytest.approx(0.2144, abs=5e-4)
    assert receipt["derived_constants"]["a_dS_m_s2"] == pytest.approx(5.415e-10, rel=1e-3)
    assert receipt["derived_constants"]["ns_edge_center_P_C"] == pytest.approx(0.96602, abs=5e-6)


def test_prose_avoids_banned_vocabulary() -> None:
    banned = re.compile(
        r"\b(now|already|still|currently|recently|previously|no longer|not yet|remains open|so far|future work|next step)\b|—|\bhonest|\bnot\b[^.!?\n]{0,120}\bbut\b",
        re.IGNORECASE,
    )
    for path in (RECEIPT, MARKDOWN, PUBLIC, CORPUS):
        text = path.read_text(encoding="utf-8")
        assert not banned.search(text), f"{path.name}: {banned.search(text).group(0)!r}"
