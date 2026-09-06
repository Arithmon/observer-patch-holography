"""Independent permutation/path controls for the finite K1 survey."""

from itertools import permutations, product
from pathlib import Path
from types import SimpleNamespace
import json
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_k1_carrier_readout_survey as survey


ELEMENTS = tuple(permutations(range(3)))
INDEX = {p: i for i, p in enumerate(ELEMENTS)}


def inverse(p):
    return tuple(p.index(i) for i in range(3))


def conjugacy_class(p):
    fixed = sum(p[i] == i for i in range(3))
    return {3: 0, 1: 1, 0: 2}[fixed]


CLASSES = np.array([conjugacy_class(p) for p in ELEMENTS])
INVERSES = np.array([INDEX[inverse(p)] for p in ELEMENTS])
MULTIPLICATION = np.array([
    [INDEX[tuple(p[q[i]] for i in range(3))] for q in ELEMENTS]
    for p in ELEMENTS
])


def oracle_counts(left, right, labels):
    """Enumerate distinct vertex triples and compose permutations directly."""
    arrows = {}
    for i, j, label in zip(left, right, labels):
        p = ELEMENTS[label]
        arrows[int(i), int(j)] = p
        arrows[int(j), int(i)] = inverse(p)
    counts = np.zeros((3, 3), dtype=np.int64)
    vertices = set(left) | set(right)
    for i, j, k in permutations(vertices, 3):
        if (i, j) in arrows and (j, k) in arrows:
            p, q = arrows[i, j], arrows[j, k]
            composed = tuple(p[q[n]] for n in range(3))
            counts[conjugacy_class(p), conjugacy_class(composed)] += 1
    return counts


def counted(left, right, labels):
    return survey.fusion_counts(
        np.asarray(left), np.asarray(right), np.asarray(labels),
        CLASSES, MULTIPLICATION, INVERSES)


@pytest.mark.parametrize("first,second", product(range(6), repeat=2))
def test_all_two_edge_label_pairs_match_permutation_oracle(first, second):
    left, right, labels = [0, 2], [1, 1], [first, second]
    np.testing.assert_array_equal(
        counted(left, right, labels), oracle_counts(left, right, labels))


def test_three_cycle_reverse_edge_cancels_to_identity():
    cycle = INDEX[(1, 2, 0)]
    counts = counted([0, 2], [1, 1], [cycle, cycle])
    expected = np.zeros((3, 3))
    expected[2, 0] = 2
    np.testing.assert_array_equal(counts, expected)
    # Copying the forward label into reverse adjacency counts a^2 here.
    assert counts[2, 2] == 0


LEFT = np.array([0, 0, 1, 2, 3])
RIGHT = np.array([1, 2, 2, 3, 0])
LABELS = np.array([INDEX[(1, 2, 0)], INDEX[(1, 0, 2)],
                   INDEX[(2, 0, 1)], INDEX[(0, 2, 1)], 0])


@pytest.mark.parametrize("mask", range(1 << len(LEFT)))
def test_every_edge_storage_reversal_preserves_exact_counts(mask):
    left, right, labels = LEFT.copy(), RIGHT.copy(), LABELS.copy()
    flip = np.array([bool(mask & (1 << i)) for i in range(len(left))])
    left[flip], right[flip] = right[flip], left[flip]
    labels[flip] = INVERSES[labels[flip]]
    np.testing.assert_array_equal(
        counted(left, right, labels), oracle_counts(LEFT, RIGHT, LABELS))


def test_edge_order_and_backtracking_exclusion():
    expected = oracle_counts(LEFT, RIGHT, LABELS)
    for order in (np.arange(len(LEFT))[::-1], np.array([2, 4, 0, 3, 1])):
        np.testing.assert_array_equal(
            counted(LEFT[order], RIGHT[order], LABELS[order]), expected)
    degrees = np.bincount(np.concatenate([LEFT, RIGHT]))
    assert expected.sum() == sum(d * (d - 1) for d in degrees)
    np.testing.assert_array_equal(counted([0], [1], [3]), np.zeros((3, 3)))


def oracle_readout(left, right, labels):
    counts = oracle_counts(left, right, labels)
    row_sums = counts.sum(axis=1, keepdims=True)
    transition = counts / np.where(row_sums > 0, row_sums, 1)
    gram = transition @ transition.T
    eigenvalues = np.linalg.eigvalsh(gram - np.trace(gram) / 3 * np.eye(3))
    gaps = np.diff(eigenvalues)
    return eigenvalues, 3 * gaps[1] / (2 * gaps[1] + gaps[0])


def test_family_b_physical_and_shuffle_outputs_use_oriented_counts(tmp_path, monkeypatch):
    run = tmp_path / "runs" / "fixture"
    run.mkdir(parents=True)
    np.savez(run / "s3_gauge_state.npz", left=LEFT, right=RIGHT, gauge=LABELS)
    monkeypatch.setattr(survey, "SIM_ROOT", tmp_path)
    monkeypatch.setattr(survey, "GAUGE_RUNS", ("fixture",))
    monkeypatch.setattr(survey, "_sim_module", lambda *args: SimpleNamespace(
        S3_CLASS=CLASSES, S3_MUL=MULTIPLICATION, S3_INV=INVERSES))
    cell, = survey.family_b_cells()
    expected = [(cell["physical"], LABELS)]
    for seed in survey.NULL_SEEDS:
        expected.append((cell["nulls"][f"label_shuffle_seed{seed}"],
                         np.random.default_rng(seed).permutation(LABELS)))
    for row, labels in expected:
        eigenvalues, rho = oracle_readout(LEFT, RIGHT, labels)
        np.testing.assert_allclose(row["eigenvalues_centered"], eigenvalues,
                                   rtol=0, atol=1e-14)
        assert row["rho_ord"] == pytest.approx(rho, abs=1e-14)


@pytest.fixture
def declared_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(survey, "SIM_ROOT", tmp_path)
    monkeypatch.setattr(survey, "SCAN_PATH", tmp_path / "scan.json")
    monkeypatch.setattr(survey, "RESPONSE_RUNS", ("response",))
    monkeypatch.setattr(survey, "GAUGE_RUNS", ("gauge",))
    monkeypatch.setattr(survey, "TOWER_RUNS", ("tower",))
    monkeypatch.setattr(survey, "TIMELINE_RUNS", ("timeline",))
    relative = (
        "scan.json",
        "runs/response/modular_response_kernel_cache.json",
        "runs/response/modular_response_kernel_payload.npz",
        "runs/gauge/s3_gauge_state.npz",
        "runs/tower/s3_gauge_state.npz",
        "runs/tower/freezeout_fields.npz",
        "runs/timeline/defect_timeline_report.json",
    )
    paths = tuple(tmp_path / p for p in relative)
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"preflight presence fixture")
    return paths


def test_all_declared_inputs_pass_presence_preflight(declared_inputs):
    survey.require_declared_inputs()


@pytest.mark.parametrize("index", range(7))
def test_each_missing_declared_input_is_rejected(declared_inputs, index):
    missing = declared_inputs[index]
    missing.unlink()
    with pytest.raises(FileNotFoundError, match="refusing a partial survey") as exc:
        survey.require_declared_inputs()
    assert missing.as_posix() in str(exc.value)


def test_missing_inputs_cannot_replace_existing_report(tmp_path, monkeypatch):
    monkeypatch.setattr(survey, "SIM_ROOT", tmp_path / "missing-sim")
    monkeypatch.setattr(survey, "SCAN_PATH", tmp_path / "missing-scan.json")
    report = tmp_path / "historical.json"
    original = b'{"historical": "do not replace with an empty survey"}\n'
    report.write_bytes(original)
    monkeypatch.setattr(sys, "argv", ["survey", "--output", str(report)])
    with pytest.raises(FileNotFoundError, match="refusing a partial survey"):
        survey.main()
    assert report.read_bytes() == original


def test_historical_affected_rows_remain_explicitly_unreplayed():
    report = json.loads(survey.OUT_PATH.read_text(encoding="utf-8"))
    correction = report["orientation_correction"]
    assert report["generated_utc"] == "2026-07-12T07:23:35Z"
    assert len(report["cells"]) == 19
    assert report["shortlist"] == []
    assert report["guards"]["public_promotion_allowed"] is False
    assert correction["corrected_orientation_replay_performed"] is False
    assert correction["historical_shortlist_valid_for_corrected_estimator"] is False
    affected = [c for c in report["cells"] if c["family"] == "B_class_fusion_law"]
    assert [c["cell_id"] for c in affected] == correction["affected_cell_ids"]
    assert [c["physical"]["rho_ord"] for c in affected] == [
        1.499999520686712, 1.4999885814885956, 1.4999982221156978]
    for row in affected:
        assert row["current_corrected_orientation_evidence_valid"] is False
        assert row["orientation_evidence_status"] == "historical_uncorrected_reverse_edge_labels"
