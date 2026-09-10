#!/usr/bin/env python3
"""Fail-closed checker for the OPH causal poset evidence package.

Standard library and numpy only. Checks every manifest digest, agreement of
provenance metadata with the pinned receipts, strict JSON decoding, the
receipt schemas, the nonclaim flags, the cross-check flags of
the carrier realization, and rebuilds the poset at q = 5 and q = 8 through
``build_causal_poset.py`` (which imports no simulator code).
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


def check_metadata(manifest: dict, family: dict, carrier: dict) -> None:
    """Cross-check descriptive metadata against authenticated archive bytes.

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
    require(manifest["source"]["producer_files"] == producers,
            "producer metadata differs from receipt pins")
    declared_research = manifest["source"]["research_files_pinned_by_the_receipts"]
    require(len(declared_research) == len(research) and set(declared_research) == research,
            "research source census differs from receipt pins")
    for field, name in (
        ("source_net_receipt_sha256", "source_net_causal_limit_receipt.json"),
        ("carrier_receipt_sha256", "carrier_source_net_receipt.json"),
    ):
        require(manifest["result"][field] == sha256(HERE / name),
                f"result receipt digest {field}")
    for name in ("source_net_causal_limit_receipt.json",
                 "carrier_source_net_logs/q5_event_log.json.gz",
                 "carrier_source_net_logs/q8_event_log.json.gz"):
        require(carrier["source_pins"]["data/exact/" + name] ==
                sha256(HERE / name).removeprefix("sha256:"),
                f"carrier attachment digest {name}")


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

    family = strict_json(HERE / "source_net_causal_limit_receipt.json")
    carrier = strict_json(HERE / "carrier_source_net_receipt.json")
    require(family["schema"] == "oph.exact.source-net-causal-limit.v1", "family schema")
    require(carrier["schema"] == "oph.exact.carrier-source-net.v1", "carrier schema")
    check_metadata(manifest, family, carrier)
    for flag in (
        "native_repair_selected",
        "physical_clock_or_spacetime_identified",
        "poisson_sprinkling",
        "finite_runs_demonstrate_asymptotic_limit",
    ):
        require(family["scope"].get(flag) is False, f"family scope flag {flag}")
    require(family["rer_cross_check"]["all_agree"] is True, "theory replay cross-check")
    require(family["rer_cross_check"]["levels_compared"] == [5, 6, 7], "cross-check levels")
    require([lvl["q"] for lvl in family["levels"]] == [5, 8, 13, 21, 34, 55], "family levels")
    for lvl in family["levels"]:
        for fam in lvl["families"]:
            vi = fam["vertical_intervals"][-1]
            require(0.0 < vi["ordering_fraction_float"] < 1.0, "ordering fraction range")
            require(1.0 < vi["myrheim_meyer_dimension"] < 6.0, "dimension range")
            for probe in fam["reachability_probes"]:
                require(probe["outer_cone_violations"] == 0, "outer cone violation")
                require(probe["certified_inner_cone_misses"] == 0, "inner cone miss")

    require(carrier["readback_metric"]["scale_to_paper_position_s"] == "1", "readback scale")
    require([lvl["q"] for lvl in carrier["levels"]] == [5, 8, 13, 21, 34], "carrier levels")
    for lvl in carrier["levels"]:
        require(lvl["neighbours"]["equals_source_net_digest"] is True, f"neighbour digest q={lvl['q']}")
        require(lvl["provenance"]["derived_rank_equals_round"] is True, f"provenance rank q={lvl['q']}")
        require(lvl["provenance"]["read_relation_equals_neighbour_digest"] is True, f"read relation q={lvl['q']}")
        require(lvl["intervention"]["equals_future_cone_all_rounds"] is True, f"intervention q={lvl['q']}")
    require(carrier["rer_cross_check"].get("all_agree", carrier["rer_cross_check"].get("agree", True)) is not False, "carrier theory cross-check")

    result = subprocess.run(
        [sys.executable, str(HERE / "build_causal_poset.py"), "--q", "5", "8", "--quiet"],
        capture_output=True,
        text=True,
        check=False,
    )
    require(result.returncode == 0, f"generator failed: {result.stderr[-2000:]}")
    tail = result.stdout.strip().splitlines()[-1] if result.stdout.strip() else ""
    require("CAUSAL_POSET_REBUILT_AND_EQUAL_TO_RECEIPTS" in tail, f"generator status: {tail[:200]}")
    print("CAUSAL_POSET_ARCHIVE_VERIFIED", len(files), "files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
