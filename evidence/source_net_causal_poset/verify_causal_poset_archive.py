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



def _profiles(K: int):
    import math
    h = 1.0 / (2.0 * K)
    g = 1.0 / K
    m = (math.sqrt(2.0) - 1.0) / K
    return {
        "constant": (lambda j: 1.0, lambda k: float(k)),
        "de_sitter": (lambda j: 1.0 / (1.0 - h * j), lambda k: -math.log(1.0 - h * k) / h),
        "radiation": (lambda j: 1.0 + g * j, lambda k: k + g * k * k / 2.0),
        "matter": (lambda j: (1.0 + m * j) ** 2, lambda k: ((1.0 + m * k) ** 3 - 1.0) / (3.0 * m)),
    }


def _sig12(x: float) -> float:
    return 0.0 if x == 0 else float(f"{float(x):.12g}")


def _near(a: float, b: float, rel: float = 1e-9) -> bool:
    return abs(a - b) <= rel * max(abs(a), abs(b), 1e-300) + 1e-12


def check_flrw_readout(name: str, family_name: str, generator: str) -> None:
    """Recompute the FLRW record-density readout from the family receipt's per-layer counts."""
    readout = strict_json(HERE / name)
    require(readout["schema"] == "oph.exact.flrw-record-density-readout.v1", "FLRW readout schema")
    require(readout["source_receipt"]["sha256"] == sha256(HERE / family_name).removeprefix("sha256:"), "FLRW readout pins the family receipt")
    require((HERE / generator).is_file(), "FLRW readout generator present")
    family = strict_json(HERE / family_name)
    require(len(readout["levels"]) == len(family["levels"]), "FLRW readout level count")
    for row, lv in zip(readout["levels"], family["levels"]):
        fam = lv["families"][0]
        require(fam["dimension"] == 3 and row["q"] == lv["q"], "FLRW readout level order")
        K = int(fam["layer_steps"])
        ladder = {int(e["layers"]): [int(c) for c in e["counts_by_layer"]] for e in fam["vertical_intervals"]}
        kI, kJ = int(fam["count_clock"]["interval_layers"]), int(fam["count_clock"]["reference_layers"])
        require(row["K"] == K and row["interval_layers"] == kI and row["reference_layers"] == kJ, "FLRW readout layers")
        nI, nJ = ladder[kI], ladder[kJ]
        require(row["counts_by_layer_interval"] == nI and row["counts_by_layer_reference"] == nJ, "FLRW readout counts")
        flat = (sum(nI) / sum(nJ)) ** 0.25
        require(_near(row["flat_count_clock"], _sig12(flat)), "FLRW flat clock")
        for pname, (sigma, tau) in _profiles(K).items():
            p = row["profiles"][pname]
            mass = lambda counts, off=0: sum(sigma(off + j) ** 4 * n for j, n in enumerate(counts))
            for entry in p["ladder"]:
                k = int(entry["layers"])
                require(_near(entry["mass"], _sig12(mass(ladder[k]))), f"FLRW mass q={lv['q']} {pname} k={k}")
            MI, MJ = mass(nI), mass(nJ)
            physical = (MI / MJ) ** 0.25
            cc = p["expanding_count_clock"]
            require(_near(cc["physical_clock"], _sig12(physical)), f"FLRW physical clock q={lv['q']} {pname}")
            require(_near(cc["proper_time_ratio"], _sig12(tau(kI) / tau(kJ))), f"FLRW proper-time ratio q={lv['q']} {pname}")
            sI = [sigma(j) for j in range(kI + 1)]
            sJ = [sigma(j) for j in range(kJ + 1)]
            lo, hi = (min(sI) / max(sJ)) * flat, (max(sI) / min(sJ)) * flat
            require(lo * (1 - 1e-12) <= physical <= hi * (1 + 1e-12) and cc["enclosure_holds"] is True, f"FLRW enclosure q={lv['q']} {pname}")
            n_e, n_0 = mass(nJ, 0), mass(nJ, K - kJ)
            reading = (n_0 / n_e) ** 0.25
            require(_near(p["redshift"]["one_plus_z_reading"], _sig12(reading)), f"FLRW redshift q={lv['q']} {pname}")
            require(p["sandwich_holds_all"] is True and p["redshift"]["enclosure_holds"] is True, f"FLRW flags q={lv['q']} {pname}")
    s = readout["summary"]
    require(s["sandwich_holds_all"] and s["clock_enclosure_holds_all"] and s["redshift_enclosure_holds_all"], "FLRW summary flags")
    require(len(readout["nonclaims"]) >= 4, "FLRW nonclaims")


def check_sampled_extension(name: str, generator: str, q: int, fibonacci_index: int) -> None:
    """A standalone sampled level beyond the exact family: schema, level, counting modes, clock arithmetic."""
    receipt = strict_json(HERE / name)
    require(receipt["schema"] == "oph.exact.source-net-causal-limit.v1", "sampled extension schema")
    require((HERE / generator).is_file(), "sampled extension generator present")
    require("extension_note" in receipt, "sampled extension note")
    require([lv["q"] for lv in receipt["levels"]] == [q], "sampled extension level")
    lv = receipt["levels"][0]
    require(lv["fibonacci_index"] == fibonacci_index, "sampled extension Fibonacci index")
    for fam in lv["families"]:
        top = fam["vertical_intervals"][-1]
        require(0.0 < top["ordering_fraction_float"] < 1.0 and 1.0 < top["myrheim_meyer_dimension"] < 6.0, "sampled extension ranges")
        for probe in fam["reachability_probes"]:
            require(probe["outer_cone_violations"] == 0 and probe["certified_inner_cone_misses"] == 0, "sampled extension cones")
        cc = fam["count_clock"]
        num, den = cc["clock_exponent"].split("/")
        require(int(den) == fam["dimension"] + 1 and int(num) == 1, "sampled extension clock exponent is the (dim + 1)-th root")
        require(_near(cc["count_clock"], (cc["interval_count"] / cc["reference_count"]) ** (int(num) / int(den)), 1e-9), "sampled extension clock arithmetic")
        if fam["dimension"] == 3:
            require(fam["vertical_pair_counting"] == "stratified_sample" and top["ordering_fraction_standard_error"] > 0, "sampled extension counting mode")
            require(fam["moving_tip_interval"]["pair_counting"] == "stratified_sample", "sampled extension moving counting")
            fractions = [e["ordering_fraction_float"] for e in fam["vertical_intervals"] if e["layers"] % 2 == 0 and e["layers"] >= 4]
            require(all(a < b for a, b in zip(fractions, fractions[1:])), "sampled extension even-K ladder monotone")
        else:
            require(fam["vertical_pair_counting"] == "exact_all_pairs", "sampled extension controls exact")
    for flag in ("native_repair_selected", "physical_clock_or_spacetime_identified", "poisson_sprinkling", "finite_runs_demonstrate_asymptotic_limit"):
        require(receipt["scope"].get(flag) is False, f"sampled extension scope flag {flag}")


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
    for ext in manifest["extensions"]:
        kind = ext.get("kind")
        if kind == "flrw_record_density_readout":
            check_flrw_readout(ext["file"], ext["source_receipt"], ext["generator"])
        elif kind == "sampled_level_extension":
            check_sampled_extension(ext["file"], ext["generator"], int(ext["q"]), int(ext["fibonacci_index"]))
    print("CAUSAL_POSET_ARCHIVE_VERIFIED", len(files), "files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
