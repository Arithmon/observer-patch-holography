#!/usr/bin/env python3
"""Verify the archived particle-simulation receipts without simulator imports.

Checks, in order: the manifest schema and the exact file inventory; the
byte count and SHA-256 of every archived receipt; strict JSON (no duplicate
keys, no non-finite constants); the calibration-null audit (both 65,536-patch
runs completed every cycle, every upstream pin is a SHA-256 digest, neither
run lands in the legacy band, and the receipt permits no quark prediction);
and the permutation-transport assay, recomputed from its stored singular-value
ratios: each state's distance to both POFT targets and its Haar rank-one
distance, each spectral-match and Haar flag against the stored thresholds,
and the aggregate receipts and verdict from those flags.

The upstream run arrays are not archived here, so this replays the recorded
report values and verdict logic, not the simulator runs themselves.

    python3 evidence/particle_simulation_receipts/verify_particle_simulation_receipts.py
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "manifest.json"
MANIFEST_SCHEMA = "oph.particle-simulation-receipt-archive.v1"
CONTROL_FILES = {"README.md", "manifest.json", "verify_particle_simulation_receipts.py"}
K1_AUDIT = "k1_64k_completion_audit.json"
POFT_ASSAY = "poft_transport_emission_targeted_20260712.json"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
UPSTREAM_FILES = {"config.yml", "defect_timeline_report.json", "freezeout_fields.npz", "s3_gauge_state.npz"}
POFT_VERDICT = "CURRENT_S3_EDGE_CARRIER_HAAR_EQUILIBRATED_NOT_POFT"


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def _reject_constant(value: str) -> None:
    raise VerificationError(f"non-finite JSON constant: {value}")


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_strict(path: Path) -> Any:
    return json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_no_duplicates,
        parse_constant=_reject_constant,
    )


def verify_manifest() -> None:
    manifest = load_strict(MANIFEST)
    require(manifest.get("schema") == MANIFEST_SCHEMA, "manifest schema drift")
    listed = [entry["path"] for entry in manifest["files"]]
    require(sorted(listed) == sorted({K1_AUDIT, POFT_ASSAY}), "manifest file list drift")
    present = {path.name for path in ROOT.iterdir() if path.is_file()}
    require(present == set(listed) | CONTROL_FILES, f"inventory drift: {sorted(present)}")
    for entry in manifest["files"]:
        raw = (ROOT / entry["path"]).read_bytes()
        require(len(raw) == entry["bytes"], f"{entry['path']}: byte count drift")
        require(hashlib.sha256(raw).hexdigest() == entry["sha256"], f"{entry['path']}: sha256 drift")


def verify_k1_audit() -> None:
    audit = load_strict(ROOT / K1_AUDIT)
    require(audit["artifact"] == "oph_k1_64k_completion_audit", "k1 artifact name drift")
    require(audit["public_quark_prediction_allowed"] is False, "k1 audit permits a quark prediction")
    require(audit["row_class"] == "compare_only_calibration_diagnostic", "k1 row class drift")
    for name in ("dense_run", "fusion_run"):
        run = audit[name]
        progress = run["progress"]
        require("64k" in run["run_id"], f"{name}: run id does not name a 65,536-patch run")
        require(progress["completed_cycles"] == progress["cycles"] > 0, f"{name}: incomplete run")
        require(progress["stage"] == "base_repair_loop_complete", f"{name}: stage drift")
        require(set(run["hashes_sha256"]) == UPSTREAM_FILES, f"{name}: upstream pin set drift")
        require(all(SHA256.match(value) for value in run["hashes_sha256"].values()), f"{name}: malformed pin")
        require(run["inside_legacy_band"] is False, f"{name}: lands in the legacy band")
    verdict = audit["verdict"]
    require(verdict["tower_result"] == "outside_legacy_band_and_at_calibration_null", "k1 tower verdict drift")
    require("cannot validate a Yukawa prediction" in verdict["physical_reading"], "k1 physical boundary drift")


def verify_poft_assay() -> None:
    assay = load_strict(ROOT / POFT_ASSAY)
    require(assay["artifact"] == "oph_poft_transport_emission_targeted_assay_v1", "assay artifact name drift")
    require(assay["comparison_uses_quark_masses"] is False, "assay compares against quark masses")
    thresholds = assay["thresholds"]
    haar_max = thresholds["haar_rank_one_max_nontrivial_ratio"]
    match_max = thresholds["poft_singular_ratio_max_abs"]
    targets = {key: assay["poft_targets"][key]["singular_value_ratios"] for key in ("T0", "T1")}
    states = assay["states"]
    require(sorted(state["node_count_seen"] for state in states) == [4096, 4096, 65536, 65536], "state census drift")
    seen = {"T0": False, "T1": False}
    for state in states:
        label = state["label"]
        ratios = state["singular_value_ratios"]
        require(len(ratios) == 3 and ratios[0] == 1.0, f"{label}: singular ratios are not normalized")
        require(SHA256.match(state["source_sha256"]) is not None, f"{label}: malformed source pin")
        haar = max(ratios[1:])
        require(haar == state["haar_rank_one_distance"], f"{label}: Haar distance does not replay")
        require(state["haar_rank_one_compatible"] is (haar <= haar_max), f"{label}: Haar flag does not replay")
        for key, target in targets.items():
            distance = max(abs(a - b) for a, b in zip(ratios, target))
            require(distance == state["max_abs_singular_ratio_distance"][key], f"{label}: {key} distance does not replay")
            match = distance <= match_max
            require(state[f"poft_{key}_necessary_spectral_match"] is match, f"{label}: {key} match flag does not replay")
            seen[key] = seen[key] or match
    receipts = assay["receipts"]
    require(receipts["direct_T0_necessary_spectral_shape_seen"] is seen["T0"], "aggregate T0 shape flag does not replay")
    require(receipts["direct_T1_necessary_spectral_shape_seen"] is seen["T1"], "aggregate T1 shape flag does not replay")
    require(not any(receipts.values()), "a direct POFT emission receipt is true")
    all_haar = all(state["haar_rank_one_compatible"] for state in states)
    require(all_haar and not (seen["T0"] or seen["T1"]), "verdict premises do not replay")
    require(assay["verdict"] == POFT_VERDICT, "assay verdict drift")


def main() -> int:
    try:
        verify_manifest()
        verify_k1_audit()
        verify_poft_assay()
    except (VerificationError, KeyError, TypeError) as error:
        print(f"FAIL: {error}")
        return 1
    print("particle simulation receipts: manifest, calibration-null audit, and transport assay replay OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
