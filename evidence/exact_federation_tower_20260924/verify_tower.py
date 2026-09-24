"""Standalone verifier for the exact-federation tower receipts (levels eight to ten).

Imports no simulator code.  For every level receipt it recomputes from the declared rules alone
what a receipt can be held to without the geometry cache and the terminal arrays, which stay
on the producing host because of their size:

* the carrier, port and seam counts of the tower level (``20 * 4^L`` carriers, twelve ports,
  thirty intra seams per carrier, three glued ports per carrier);
* the loads ``default_rng(20260909 + level).integers(0, 6, size=ports)`` and, from them and the
  declared single component, the balanced-class minimum, the expected component-multiset hash
  and the mean-law minimum;
* for every integer schedule: the seed rule, termination, ``V_terminal`` equal to the minimum,
  the quotient hash equal to the expected hash, the ledger monotone from ``V_initial`` to the
  minimum with one entry per sweep, the attempt index inside the last sweep, the digest chains
  of the declared length, and the summary fields recomputed from the entries;
* for the mean schedule: the withheld-hash rule and the flags;
* for the response kernels: symmetry, trace twelve, the slow-band share recomputed from the
  matrix with the carrier's slow-band projector rebuilt here from the icosahedral seam graph,
  and the summaries;
* the per-schedule ``schedule.json`` files against the receipt entries, and the verification
  records ``verify_*.json`` written by the full replay on the producing host (which checked the
  terminal arrays, every draw digest and a layered replay from the last checkpoint), with
  their verdicts;
* the manifest digests of every file in the package.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PORTS = 12
INTRA = 30
LOAD_MAX = 5
LOAD_SEED_BASE = 20260909
SCHEDULE_SEED_BASE = 909000


class Failure(AssertionError):
    pass


def require(c: bool, m: str) -> None:
    if not c:
        raise Failure(m)


def canonical(x) -> bytes:
    return (json.dumps(x, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False) + "\n").encode("ascii")


def sha256_of(x) -> str:
    return hashlib.sha256(canonical(x)).hexdigest()


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _sig(x: float, digits: int) -> float:
    return 0.0 if x == 0 else float(f"{float(x):.{digits}g}")


# --- the carrier's slow-band projector from the icosahedral seam graph --------------------


CARRIER_SEAMS = [(0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (1, 5), (1, 7), (1, 8), (1, 9), (2, 3), (2, 4), (2, 6), (2, 10), (2, 11), (3, 4), (3, 6), (3, 8), (3, 9), (4, 5), (4, 9), (4, 11), (5, 9), (5, 11), (6, 7), (6, 8), (6, 10), (7, 8), (7, 10), (8, 9), (10, 11)]  # carrier.seams() of the simulator, ports 0..11 in the carrier's labelling


def slow_band_projector() -> np.ndarray:
    """Spectral projector of the carrier seam Laplacian onto its 5 - sqrt5 band (rank three)."""

    lap = np.zeros((PORTS, PORTS))
    for a, b in CARRIER_SEAMS:
        lap[a, a] += 1
        lap[b, b] += 1
        lap[a, b] -= 1
        lap[b, a] -= 1
    require(len(CARRIER_SEAMS) == INTRA and np.allclose(np.diag(lap), 5.0), "the carrier seam graph is not five-regular on twelve ports")
    w, u = np.linalg.eigh(lap)
    band = np.isclose(w, 5 - 5**0.5)
    require(int(band.sum()) == 3, "the 5 - sqrt5 band of the carrier seam Laplacian is not three-dimensional")
    b = u[:, band]
    return b @ b.T


# --- expectations from the loads alone ---------------------------------------------------


def expected(level: int, ports: int) -> dict:
    loads = np.random.default_rng(LOAD_SEED_BASE + level).integers(0, LOAD_MAX + 1, size=ports, dtype=np.int64)
    total = int(loads.sum())
    q, r = divmod(total, ports)
    vc = []
    if ports - r > 0:
        vc.append([q, ports - r])
    if r > 0:
        vc.append([q + 1, r])
    return {
        "V_initial": int(np.dot(loads, loads)),
        "balanced_minimum": (ports - r) * q * q + r * (q + 1) ** 2,
        "expected_hash": sha256_of({"canonicalizer": "component_multiset", "components": [[ports, vc]]}),
        "mean_minimum": Fraction(total * total, ports),
        "loads_int8_sha256": hashlib.sha256(loads.astype("<i1").tobytes()).hexdigest(),
        "total": total,
    }


def check_level(directory: Path, p_slow: np.ndarray) -> dict:
    receipt = json.loads((directory / "receipt.json").read_text())
    require(receipt["schema"] == "oph.exact.federation-huge.v1", "schema")
    level = int(receipt["level"])
    carriers = 20 * 4**level
    ports = PORTS * carriers
    inter = 3 * carriers // 2
    seams = INTRA * carriers + inter
    require(receipt["carriers"] == carriers and receipt["ports"] == ports and receipt["inter_seams"] == inter and receipt["seams"] == seams, f"L{level}: counts")
    require(receipt["components"] == 1, f"L{level}: the tower gluing is connected; a different component count needs the geometry")
    exp = expected(level, ports)
    require(receipt["expected"]["expected_integer_quotient_hash"] == exp["expected_hash"], f"L{level}: expected hash")
    require(receipt["expected"]["balanced_minimum"] == exp["balanced_minimum"], f"L{level}: balanced minimum")
    require(abs(receipt["expected"]["mean_minimum"] - float(exp["mean_minimum"])) <= 1e-9 * float(exp["mean_minimum"]), f"L{level}: mean minimum")
    entries = receipt["integer_law"]["entries"]
    seeds = [SCHEDULE_SEED_BASE + 1000 * level + 100 + k for k in range(len(entries))]
    require([e["seed"] for e in entries] == seeds == receipt["seeds"]["schedules"]["integer"], f"L{level}: seeds")
    for e in entries:
        require(e["terminated"] and e["V_terminal"] == e["V_minimum"] == exp["balanced_minimum"], f"L{level} seed {e['seed']}: V at termination")
        require(e["V_initial"] == exp["V_initial"], f"L{level} seed {e['seed']}: V_initial")
        require(e["quotient_hash"] == exp["expected_hash"] and e["quotient_hash_equals_expected"], f"L{level} seed {e['seed']}: quotient hash")
        ledger = e["V_ledger"]
        require(ledger[0] == e["V_initial"] and ledger[-1] == e["V_terminal"] and len(ledger) == e["sweeps"] + 1, f"L{level} seed {e['seed']}: ledger ends")
        require(all(a >= b for a, b in zip(ledger, ledger[1:])), f"L{level} seed {e['seed']}: ledger monotone")
        require((e["sweeps"] - 1) * seams < e["attempts_to_balanced_class"] <= e["sweeps"] * seams == e["attempts"], f"L{level} seed {e['seed']}: attempts")
        require(len(e["draw_sha256_per_sweep"]) == e["sweeps"] and all(len(d) == 2 for d in e["draw_sha256_per_sweep"]), f"L{level} seed {e['seed']}: draw digests")
        require(len(e["state_sha256_per_sweep"]) == e["sweeps"] + 1 and e["state_sha256_per_sweep"][0] == exp["loads_int8_sha256"], f"L{level} seed {e['seed']}: state digest chain")
        require(e["state_sha256_per_sweep"][-1] == e["terminal_state"]["sha256"], f"L{level} seed {e['seed']}: terminal digest")
        require(e["unit_transfer_decrement_identity_violations"] == 0 and e["conservation_exact"] and e["descent_ledger_exact"], f"L{level} seed {e['seed']}: flags")
        sched_path = directory / f"integer_{e['seed']}.json"
        if sched_path.is_file():
            sched = json.loads(sched_path.read_text())
            for key in ("seed", "sweeps", "quotient_hash", "attempts_to_balanced_class", "descents", "swaps", "waits", "unit_transfers", "V_ledger", "draw_sha256_per_sweep", "state_sha256_per_sweep"):
                require(sched[key] == e[key], f"L{level} seed {e['seed']}: schedule file differs on {key}")
            require(len(sched["per_sweep"]) == e["sweeps"] and all(p["V"] == v for p, v in zip(sched["per_sweep"], e["V_ledger"][1:])), f"L{level} seed {e['seed']}: per-sweep ledger")
            require(sum(p["descents"] for p in sched["per_sweep"]) == e["descents"] and sum(p["unit_transfers"] for p in sched["per_sweep"]) == e["unit_transfers"], f"L{level} seed {e['seed']}: per-sweep counts")
    s = receipt["integer_law"]
    require(s["schedules"] == len(entries) and s["all_terminated"] and s["quotient_hash_equals_expected_all"] and s["unique_quotient_hash_count"] == 1, f"L{level}: integer summary")
    require(s["sweeps"] == {"min": min(e["sweeps"] for e in entries), "max": max(e["sweeps"] for e in entries)}, f"L{level}: sweep summary")
    require(s["violations_total"] == 0 and s["conservation_exact_all"], f"L{level}: violations")
    for e in receipt["mean_law_float"]["entries"]:
        require((e["terminal_quotient_hash"] is not None) == e["lattice_snap_unambiguous"], f"L{level}: withheld-hash rule")
        require(e["strict_descent_violations"] == 0 and e["descent_ledger_relative_error_below_1e-9"] == (e["descent_ledger_relative_error"] < 1e-9), f"L{level}: mean flags")
        require(len(e["V_ledger"]) == e["sweeps"] + 1 and all(a >= b - 1e-6 * abs(a) for a, b in zip(e["V_ledger"], e["V_ledger"][1:])), f"L{level}: mean ledger")
    m = receipt["mean_law_float"]
    require(m["ambiguous_terminal_count"] == sum(1 for e in m["entries"] if e["terminal_quotient_hash"] is None), f"L{level}: ambiguous count")
    rk = receipt["response_kernels"]
    worst = 0.0
    for entry in rk["cells"]:
        for n, kk in entry["kernels"].items():
            k = np.asarray(kk)
            require(k.shape == (PORTS, PORTS) and np.allclose(k, k.T, atol=1e-12) and abs(np.trace(k) - PORTS) < 1e-9, f"L{level}: kernel shape/symmetry/trace cell {entry['cell']} n={n}")
            share = float(np.trace(p_slow @ k @ p_slow) / np.trace(k))
            worst = max(worst, abs(share - entry["slow_band_share"][n]))
    require(worst < 1e-9, f"L{level}: slow-band shares")
    for n, summary in rk["glued_slow_band_share"].items():
        shares = [e["slow_band_share"][n] for e in rk["cells"]]
        require(abs(summary["min"] - _sig(min(shares), 6)) < 1e-9 and abs(summary["max"] - _sig(max(shares), 6)) < 1e-9, f"L{level}: share summary n={n}")
    records = sorted(directory.glob("verify_*.json"))
    verdicts = [json.loads(p.read_text()) for p in records]
    require(all(v["verdict"] == "PASS" and v["level"] == level for v in verdicts), f"L{level}: host verification records")
    covered = sorted({r["seed"] for v in verdicts for r in v["integer"]})
    require(covered == seeds, f"L{level}: host verification does not cover every schedule")
    require(all(r["draw_sweeps_checked"] == r["sweeps"] and r["replayed_sweeps"] >= 1 for v in verdicts for r in v["integer"]), f"L{level}: host verification checked fewer draws or replayed nothing")
    return {"level": level, "carriers": carriers, "seams": seams, "schedules": len(entries), "sweeps": s["sweeps"],
            "host_verification_records": len(records), "replayed_sweeps_min": min(r["replayed_sweeps"] for v in verdicts for r in v["integer"])}


def main() -> int:
    manifest = json.loads((HERE / "manifest.json").read_text())
    for rel, digest in manifest["files"].items():
        require(file_sha256(HERE / rel) == digest, f"manifest digest differs: {rel}")
    p_slow = slow_band_projector()
    results = []
    for level_dir in sorted(HERE.glob("L*")):
        if (level_dir / "receipt.json").is_file():
            results.append(check_level(level_dir, p_slow))
            print(f"L{results[-1]['level']}: PASS ({results[-1]['carriers']} carriers, {results[-1]['schedules']} schedules, sweeps {results[-1]['sweeps']}, {results[-1]['host_verification_records']} host records)")
    require(len(results) == len(manifest["levels"]), "manifest levels differ from the package")
    # reserve-generator receipts (issue 985): verified by their own standalone verifier against the settlement receipts and textures
    reserve_dir = HERE / "reserve"
    for name in manifest.get("reserve_receipts", []):
        receipt = json.loads((reserve_dir / name).read_text())
        level = int(receipt["level"])
        engine = HERE / f"L{level}" / "receipt.json"
        textures = sorted(str(p) for p in HERE.glob("L*/texture.json"))
        cmd = [sys.executable, str(reserve_dir / "verify_federation_reserve_independent.py"), str(reserve_dir / name), "--textures", *textures]
        if engine.is_file():
            cmd += ["--engine", str(engine)]
        run = subprocess.run(cmd, capture_output=True, text=True)
        require(run.returncode == 0 and run.stdout.strip().endswith("PASS"), f"reserve receipt {name}: {run.stdout[-300:]} {run.stderr[-300:]}")
        print(f"reserve L{level}: PASS ({len(receipt['schedules'])} schedules)")
    print("TOWER_PACKAGE_PASS")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Failure as error:
        print(f"FAIL: {error}")
        sys.exit(1)
