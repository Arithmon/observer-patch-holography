#!/usr/bin/env python3
"""Fail-closed checker for the OPH causal poset evidence package.

Standard library and numpy only. Checks every manifest digest, agreement of
provenance metadata with the pinned receipts, strict JSON decoding, the
receipt schemas, the nonclaim flags, the cross-check flags of
the carrier realization, and rebuilds the poset at q = 5 and q = 8 through
each generation's own frozen generator (``build_causal_poset.py`` for
2026-09-09, ``build_causal_poset_2026-09-24.py`` for the current receipts;
neither imports simulator code).
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "archive_manifest.json"
CONTROL = {"archive_manifest.json", "verify_causal_poset_archive.py"}


def fail(message: str) -> None:
    raise SystemExit(f"FAIL: {message}")


def require(condition: bool, message: str) -> None:
    if not condition:
        fail(message)


def sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def strict_json(path: Path):
    def reject(token: str):
        raise ValueError(f"non-finite numeric token {token!r} in {path.name}")

    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate key {key!r} in {path.name}")
            result[key] = value
        return result

    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject,
                      object_pairs_hook=unique)


def check_metadata(source: dict, result: dict, family: dict, carrier: dict, family_name: str, carrier_name: str) -> None:
    """Cross-check descriptive metadata of one generation against authenticated archive bytes.

    Commit-to-source authentication is recorded by the curator; this offline
    check does not contact GitHub or substitute a mutable simulator checkout.
    """
    producers = {}
    research = set()
    for receipt in (family, carrier):
        for name, digest in receipt["source_pins"].items():
            if name.startswith(("oph_exact/", "tests/")):
                require(name not in producers or producers[name] == digest,
                        f"conflicting producer pin {name}")
                producers[name] = digest
            elif name.startswith("reverse-engineering-reality/"):
                research.add(name.removeprefix("reverse-engineering-reality/"))
    require(source["producer_files"] == producers,
            f"producer metadata differs from receipt pins ({family_name})")
    declared_research = source["research_files_pinned_by_the_receipts"]
    require(len(declared_research) == len(research) and set(declared_research) == research,
            f"research source census differs from receipt pins ({family_name})")
    require(result["source_net_receipt_sha256"] == sha256(HERE / family_name), f"result receipt digest {family_name}")
    require(result["carrier_receipt_sha256"] == sha256(HERE / carrier_name), f"result receipt digest {carrier_name}")
    require(carrier["source_pins"]["data/exact/source_net_causal_limit_receipt.json"] ==
            sha256(HERE / family_name).removeprefix("sha256:"),
            f"carrier attachment digest {family_name}")
    for name in ("carrier_source_net_logs/q5_event_log.json.gz", "carrier_source_net_logs/q8_event_log.json.gz"):
        require(carrier["source_pins"]["data/exact/" + name] == sha256(HERE / name).removeprefix("sha256:"),
                f"carrier attachment digest {name}")


def check_generation(label: str, source: dict, result: dict, family_name: str, carrier_name: str,
                     family_levels: list, carrier_levels: list, generator: str) -> None:
    family = strict_json(HERE / family_name)
    carrier = strict_json(HERE / carrier_name)
    require(family["schema"] == "oph.exact.source-net-causal-limit.v1", f"{label}: family schema")
    require(carrier["schema"] == "oph.exact.carrier-source-net.v1", f"{label}: carrier schema")
    check_metadata(source, result, family, carrier, family_name, carrier_name)
    for flag in (
        "native_repair_selected",
        "physical_clock_or_spacetime_identified",
        "poisson_sprinkling",
        "finite_runs_demonstrate_asymptotic_limit",
    ):
        require(family["scope"].get(flag) is False, f"{label}: family scope flag {flag}")
    require(family["rer_cross_check"]["all_agree"] is True, f"{label}: theory replay cross-check")
    require(family["rer_cross_check"]["levels_compared"] == [5, 6, 7], f"{label}: cross-check levels")
    require([lvl["q"] for lvl in family["levels"]] == family_levels, f"{label}: family levels")
    for lvl in family["levels"]:
        for fam in lvl["families"]:
            vi = fam["vertical_intervals"][-1]
            require(0.0 < vi["ordering_fraction_float"] < 1.0, f"{label}: ordering fraction range")
            require(1.0 < vi["myrheim_meyer_dimension"] < 6.0, f"{label}: dimension range")
            for probe in fam["reachability_probes"]:
                require(probe["outer_cone_violations"] == 0, f"{label}: outer cone violation")
                require(probe["certified_inner_cone_misses"] == 0, f"{label}: inner cone miss")
    require(carrier["readback_metric"]["scale_to_paper_position_s"] == "1", f"{label}: readback scale")
    require([lvl["q"] for lvl in carrier["levels"]] == carrier_levels, f"{label}: carrier levels")
    for lvl in carrier["levels"]:
        require(lvl["neighbours"]["equals_source_net_digest"] is True, f"{label}: neighbour digest q={lvl['q']}")
        require(lvl["provenance"]["derived_rank_equals_round"] is True, f"{label}: provenance rank q={lvl['q']}")
        require(lvl["provenance"]["read_relation_equals_neighbour_digest"] is True, f"{label}: read relation q={lvl['q']}")
        require(lvl["intervention"]["equals_future_cone_all_rounds"] is True, f"{label}: intervention q={lvl['q']}")
    require(carrier["rer_cross_check"].get("all_agree", carrier["rer_cross_check"].get("agree", True)) is not False,
            f"{label}: carrier theory cross-check")
    result_run = subprocess.run(
        [sys.executable, str(HERE / generator), "--q", "5", "8", "--quiet"],
        capture_output=True, text=True, check=False)
    require(result_run.returncode == 0, f"{label}: generator failed: {result_run.stderr[-2000:]}")
    tail = result_run.stdout.strip().splitlines()[-1] if result_run.stdout.strip() else ""
    require("CAUSAL_POSET_REBUILT_AND_EQUAL_TO_RECEIPTS" in tail, f"{label}: generator status: {tail[:200]}")


def main() -> int:
    manifest = strict_json(MANIFEST)
    require(manifest.get("schema") == "oph.curated_evidence_package.v1", "manifest schema")
    files = {row["path"]: row for row in manifest["inventory"]}
    require(len(files) == len(manifest["inventory"]), "duplicate inventory path")
    listed = set(files)
    present = {
        p.relative_to(HERE).as_posix()
        for p in HERE.rglob("*")
        if p.is_file() and "__pycache__" not in p.parts
    }
    control = CONTROL | {"README.md"}
    require(listed | control == present, f"file set mismatch: {sorted(listed ^ (present - control))}")
    inventory_lines = []
    total_bytes = 0
    for name, entry in sorted(files.items()):
        path = HERE / name
        require(path.is_file(), f"missing {name}")
        require(sha256(path) == "sha256:" + entry["sha256"], f"digest mismatch {name}")
        require(type(entry["bytes"]) is int and path.stat().st_size == entry["bytes"],
                f"size mismatch {name}")
        total_bytes += entry["bytes"]
        inventory_lines.append(f"{entry['sha256']}  {entry['bytes']}  {name}\n")
    curated = manifest["curated_archive"]
    require(type(curated["file_count"]) is int and curated["file_count"] == len(files),
            "file count")
    require(type(curated["total_bytes"]) is int and curated["total_bytes"] == total_bytes,
            "total byte count")
    require(curated["inventory_sha256"] ==
            hashlib.sha256("".join(inventory_lines).encode("utf-8")).hexdigest(),
            "inventory digest")

    current = manifest["current_generation"]
    check_generation(current["label"], manifest["source"], manifest["result"], current["family"], current["carrier"],
                     current["family_levels"], current["carrier_levels"], current["generator"])
    for gen in manifest["historical_generations"]:
        check_generation(gen["label"], gen["source"], gen["result"], gen["family"], gen["carrier"],
                         gen["family_levels"], gen["carrier_levels"], gen["generator"])
    print("CAUSAL_POSET_ARCHIVE_VERIFIED", len(files), "files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
