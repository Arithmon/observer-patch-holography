#!/usr/bin/env python3
"""Verify the archived closure-loop receipt and event logs without simulator imports.

Checks, in order: every manifest digest and the inventory; strict and
canonical JSON; the receipt's self digest and schema; the log pins of the
firewall block; the invariant-vector equality flags of every loop; the
negative-control flags; and a recomputation from the canonical source log
alone (the seam set, the degree sequence, the distance-three antipode, the
rule class of every event, the seam-count chi-square, the conserved total,
the strict descent of the centered squared norm, the seam-form increases,
and the 120-digit probe readback as the exact power of the recovered
expectation operator).  The same recovery of ports, seams, pairing and rule
class is run on every stored log and compared with the loop invariants.

    python3 evidence/closure_loop/verify_closure_loop_archive.py
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "archive_manifest.json"
RECEIPT = ROOT / "closure_loop_receipt.json"
LOG_DIR = ROOT / "closure_loop_logs"
CONTROL_FILES = {"README.md", "archive_manifest.json", "verify_closure_loop_archive.py"}
MANIFEST_SCHEMA = "oph.curated_evidence_package.v1"
RECEIPT_SCHEMA = "oph.exact.closure-loop.v1"
LOG_SCHEMA = "oph.exact.closure-loop.event-log.v1"
CANONICAL_LOG = "icosahedron_seam_mean__inhabited_structure.json"
ICOSAHEDRON_ANTIPODE = [3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8]
INVARIANT_KEYS = (
    "ports", "seams", "degree_sequence", "pairing_type", "farthest_pairing_distance",
    "automorphism_order", "rotation_order", "rule_class", "nonconservative_events_present",
    "tie_law", "schedule_law", "gram_eigenvalues", "gram_rank", "gram_step_index",
    "probe_prediction_within_tolerance", "lie_split", "lie_assignment_count",
    "trivial_isotypic_dimension", "strict_descent_violation_count",
    "terminal_equals_initial_mean", "terminal_schedule_independent",
    "quotient_multiset_schedule_independent", "terminal_invariant_type",
    "components", "components_isomorphic", "global_schedule_law",
)
EXPECTED_NOT_CLAIMED = (
    "universe-level closure",
    "existence or uniqueness of the cosmic fixed point",
    "physical identification of any recovered invariant",
    "glued federation (work in progress under lane L1)",
    "records as the only observer access; the probe is a direct readback",
    "refinement across levels",
)


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def _reject_constant(value: str) -> None:
    raise VerificationError(f"non-finite JSON constant: {value}")


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in pairs:
        if key in out:
            raise VerificationError(f"duplicate JSON key: {key}")
        out[key] = value
    return out


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="ascii"), object_pairs_hook=_no_duplicates,
                      parse_constant=_reject_constant)


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                       allow_nan=False) + "\n").encode("ascii")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------
# Manifest and inventory
# --------------------------------------------------------------------------


def verify_manifest() -> dict:
    manifest = load_json(MANIFEST)
    require(manifest.get("schema") == MANIFEST_SCHEMA, "manifest schema changed")
    require(manifest.get("archive_id") == "closure_loop", "archive id changed")
    expected = {row["path"]: row for row in manifest["inventory"]}
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*")
              if p.is_file() and p.name not in CONTROL_FILES and "__pycache__" not in p.parts}
    require(actual == set(expected), f"inventory names differ: {sorted(actual ^ set(expected))}")
    lines = []
    total = 0
    for name in sorted(expected):
        row = expected[name]
        path = ROOT / name
        size = path.stat().st_size
        digest = sha256_file(path)
        require(size == row["bytes"], f"byte count mismatch: {name}")
        require(digest == row["sha256"], f"sha256 mismatch: {name}")
        total += size
        lines.append(f"{digest}  {size}  {name}\n")
    curated = manifest["curated_archive"]
    require(curated["file_count"] == len(expected), "file count mismatch")
    require(curated["total_bytes"] == total, "total byte count mismatch")
    require(hashlib.sha256("".join(lines).encode()).hexdigest() == curated["inventory_sha256"],
            "inventory digest mismatch")
    for key, value in manifest["known_boundaries"].items():
        require(value is False, f"known boundary {key} is not recorded as open")
    require(manifest["source"]["commit"] == "uncommitted-working-tree", "source commit field changed")
    return manifest


# --------------------------------------------------------------------------
# Receipt structure
# --------------------------------------------------------------------------


def verify_receipt_structure(manifest: dict) -> dict:
    raw = RECEIPT.read_bytes()
    receipt = load_json(RECEIPT)
    require(canonical(receipt) == raw, "receipt is not canonical JSON")
    require(receipt.get("schema") == RECEIPT_SCHEMA, "receipt schema changed")
    payload = dict(receipt)
    digest = payload.pop("receipt_sha256")
    require(digest == "sha256:" + hashlib.sha256(canonical(payload)).hexdigest(), "receipt self digest differs")
    for key in ("firewall", "loops", "canonical_loop", "negative_controls", "federation", "scope",
                "claim_boundary", "pins", "status", "canonical_recovered_specification", "parameters"):
        require(key in receipt, f"missing block {key}")
    for pin in receipt["pins"].values():
        producer = manifest["source"]["producer_files"].get(pin["path"])
        require(producer is not None and "sha256:" + producer == pin["sha256"],
                f"manifest producer pin differs from the receipt pin: {pin['path']}")
    require(list(receipt["scope"]["not_claimed"]) == list(EXPECTED_NOT_CLAIMED), "not-claimed list changed")
    require(receipt["federation"]["glued_federation"].startswith("work in progress"),
            "glued federation is not recorded as work in progress")
    return receipt


def verify_log_pins(receipt: dict) -> dict[str, dict]:
    logs = {}
    for name, row in receipt["firewall"]["logs"].items():
        path = LOG_DIR / name
        if not row["stored"]:
            require(row["path"] is None and not path.exists(), f"unstored log present: {name}")
            continue
        require(row["path"] == f"data/exact/closure_loop_logs/{name}", f"log path changed: {name}")
        require(path.is_file(), f"missing log {name}")
        require("sha256:" + sha256_file(path) == row["sha256"], f"log digest mismatch: {name}")
        log = load_json(path)
        require(canonical(log) == path.read_bytes(), f"log is not canonical JSON: {name}")
        require(log["schema"] == LOG_SCHEMA, f"log schema changed: {name}")
        require(len(log["events"]) == row["events"], f"log event count changed: {name}")
        logs[name] = log
    require(len(logs) == 15, "fifteen stored logs expected")
    return logs


# --------------------------------------------------------------------------
# Loops, invariant vectors, negative controls
# --------------------------------------------------------------------------


def invariant_vector(inv: dict) -> tuple:
    return tuple(canonical(inv[k]) for k in INVARIANT_KEYS)


def verify_loops(receipt: dict, logs: dict[str, dict]) -> None:
    loops = receipt["loops"]
    for name, loop in loops.items():
        stages = loop["stages"]
        require(loop["source_spec_name"] == name, f"loop name mismatch: {name}")
        require(len(stages) == loop["iterations"] + 1, f"stage count mismatch: {name}")
        vectors = [invariant_vector(s["invariants"]) for s in stages]
        agree = all(v == vectors[0] for v in vectors)
        require(agree == loop["invariant_vectors_agree"], f"invariant-vector agreement flag wrong: {name}")
        terminal = stages[0]["invariants"]["terminal_invariant_type"]
        fixed = agree and terminal in ("state", "quotient_multiset")
        require(fixed == loop["fixed_point"], f"fixed-point flag wrong: {name}")
        require(loop["fixed_point_level"] == (terminal if fixed else "none"), f"fixed-point level wrong: {name}")
        require(loop["receipt"] == ("CLOSURE_FIXED_POINT_AT_CARRIER_SCALE" if fixed else "NO_FIXED_POINT"),
                f"loop receipt string wrong: {name}")
        for stage in stages:
            log_name = f"{name}__{stage['role']}.json"
            row = receipt["firewall"]["logs"][log_name]
            require(row["sha256"] == stage["log_sha256"], f"stage log digest differs from the firewall pin: {log_name}")
            require(row["events"] == stage["events"], f"stage event count differs: {log_name}")
    require(receipt["canonical_loop"] == loops["icosahedron_seam_mean"], "canonical loop differs from its entry")
    nc = receipt["negative_controls"]
    ov = loops["icosahedron_overwrite"]
    require(nc["overwrite"]["fixed_point"] is False and ov["fixed_point"] is False, "overwrite must have no fixed point")
    require(nc["overwrite"]["invariant_vectors_agree"] is False, "overwrite invariant vectors must differ")
    require(nc["overwrite"]["nonconservative_events"] > 0, "overwrite must show nonconservative events")
    require(nc["overwrite"]["terminal_invariant_type"] == "none", "overwrite terminal invariant must be none")
    for carrier, ports, seams, auto, rot in (("tetrahedron", 4, 6, 24, 12), ("octahedron", 6, 12, 48, 24)):
        row = nc[carrier]
        require(row["ports"] == ports and row["seams"] == seams, f"{carrier} port or seam count changed")
        require(row["automorphism_order"] == auto and row["rotation_order"] == rot, f"{carrier} group orders changed")
        require(row["pairing_type"] == "none" and row["gram_rank"] == 3 and row["fixed_point"] is True,
                f"{carrier} control row changed")
    require(nc["tetrahedron"]["farthest_pairing_distance"] is None, "tetrahedron has no farthest pairing")
    require(nc["octahedron"]["farthest_pairing_distance"] == 2, "octahedron farthest pairing is distance two")
    require(nc["lie_split_by_carrier"] == {"icosahedron": "1+3+8", "octahedron": "none", "tetrahedron": "1+3"},
            "Lie split table changed")
    require(nc["gram_rank_alone_does_not_select_icosahedron"] is True, "Gram-rank selection flag changed")
    require(nc["pairing_plus_rotation_order_selects_icosahedron"] is True, "pairing-plus-rotation flag changed")
    integer = nc["integer_nearest_agreement"]
    require(integer["fixed_point_level"] == "quotient_multiset" and integer["tie_law"] == "uniform",
            "integer control row changed")
    require(receipt["federation"]["components"] == 20 and receipt["federation"]["fixed_point"] is True,
            "isolated federation block changed")
    require(receipt["status"].startswith("CLOSURE_FIXED_POINT_AT_CARRIER_SCALE_ATTAINED"), "status changed")


# --------------------------------------------------------------------------
# Recomputation from the logs
# --------------------------------------------------------------------------


def events_table(log: dict) -> list[tuple[int, int, Fraction, Fraction, Fraction, Fraction]]:
    ports = int(log["ports"])
    table = []
    for ev in log["events"]:
        changed = ev["changed"]
        require(len(changed) == 2, "an event touches a number of ports other than two")
        (a, va), (b, vb) = sorted(((int(k), v) for k, v in changed.items()), key=lambda kv: kv[0])
        require(0 <= a < b < ports, "event port outside the port range")
        table.append((a, b, Fraction(va[0]), Fraction(va[1]), Fraction(vb[0]), Fraction(vb[1])))
    return table


def graph_distances(ports: int, seams) -> list[list[int]]:
    adj = [set() for _ in range(ports)]
    for i, j in seams:
        adj[i].add(j)
        adj[j].add(i)
    dist = [[-1] * ports for _ in range(ports)]
    for s in range(ports):
        dist[s][s] = 0
        frontier = [s]
        d = 0
        while frontier:
            d += 1
            nxt = []
            for u in frontier:
                for v in adj[u]:
                    if dist[s][v] < 0:
                        dist[s][v] = d
                        nxt.append(v)
            frontier = nxt
    return dist


def involution(candidates: list[list[int]]) -> list[int] | None:
    if any(len(c) != 1 for c in candidates):
        return None
    pairing = [c[0] for c in candidates]
    if any(pairing[p] == p or pairing[pairing[p]] != p for p in range(len(pairing))):
        return None
    return pairing


def recover(log: dict) -> dict:
    """Ports, seams, degrees, pairing and rule class from the event log alone."""
    ports = int(log["ports"])
    table = events_table(log)
    seams = sorted({(a, b) for a, b, *_ in table})
    degree = [0] * ports
    for a, b in seams:
        degree[a] += 1
        degree[b] += 1
    dist = graph_distances(ports, seams)
    three = involution([[q for q in range(ports) if dist[p][q] == 3] for p in range(ports)])
    ecc = [max(row) for row in dist]
    farthest = None
    if len(set(ecc)) == 1 and ecc[0] > 0:
        pairing = involution([[q for q in range(ports) if dist[p][q] == ecc[0]] for p in range(ports)])
        if pairing is not None:
            farthest = ecc[0]
    conservative = symmetric = mean_form = overwrite_form = nearest = 0
    for a, b, ba, aa, bb, ab in table:
        if aa + ab == ba + bb:
            conservative += 1
        if aa == ab:
            symmetric += 1
        if aa == ab == (ba + bb) / 2:
            mean_form += 1
        if (aa == ba and ab == ba) or (ab == bb and aa == bb):
            overwrite_form += 1
        s = ba + bb
        if s.denominator == 1 and {aa, ab} == {Fraction(int(s) // 2), Fraction(-(-int(s) // 2))} and aa + ab == s:
            nearest += 1
    n = len(table)
    if mean_form == n:
        rule = "seam_mean"
    elif nearest == n and conservative == n:
        rule = "integer_nearest_agreement"
    elif overwrite_form == n:
        rule = "overwrite"
    else:
        rule = "unclassified"
    return {"ports": ports, "seams": seams, "seam_count": len(seams), "degree_sequence": sorted(degree),
            "pairing_type": "distance_three_involution" if three is not None else "none",
            "distance_three_pairing": three, "farthest_pairing_distance": farthest,
            "rule_class": rule, "nonconservative_events_present": conservative < n, "event_count": n,
            "table": table}


def readback_string(numerator: int, denominator: int, digits: int) -> str:
    """The receipt's rounding: ``digits`` significant decimals, half to even, as ``<mantissa>E<exponent>``."""
    if numerator == 0:
        return "0"
    sign = "-" if numerator < 0 else ""
    num = abs(numerator)
    den = denominator
    estimate = int((num.bit_length() - den.bit_length()) * 0.30102999566398)
    shift = digits - 1 - estimate
    lower, upper = 10 ** (digits - 1), 10 ** digits
    while True:
        if shift >= 0:
            divisor = den
            q, r = divmod(num * 10 ** shift, divisor)
        else:
            divisor = den * 10 ** (-shift)
            q, r = divmod(num, divisor)
        if q >= upper:
            shift -= 1
            continue
        if q < lower:
            shift += 1
            continue
        break
    twice = 2 * r
    if twice > divisor or (twice == divisor and q % 2 == 1):
        q += 1
        if q == upper:
            q //= 10
            shift -= 1
    return f"{sign}{q}E{-shift}"


def matrix_power(m: list[list[int]], e: int) -> list[list[int]]:
    n = len(m)
    result = [[int(i == j) for j in range(n)] for i in range(n)]
    base = [row[:] for row in m]

    def mul(x, y):
        return [[sum(x[i][k] * y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

    while e:
        if e & 1:
            result = mul(result, base)
        e >>= 1
        if e:
            base = mul(base, base)
    return result


def verify_canonical_log(receipt: dict, log: dict) -> dict:
    rec = recover(log)
    comp = receipt["canonical_recovered_specification"]["component"]
    require(rec["ports"] == 12 and rec["seam_count"] == 30, "canonical carrier is not (12, 30)")
    require(rec["degree_sequence"] == [5] * 12, "canonical carrier is not five-regular")
    require([list(s) for s in rec["seams"]] == comp["seams"], "recovered seam set differs from the receipt")
    require(rec["distance_three_pairing"] == comp["pairing"]["distance_three_pairing"] == ICOSAHEDRON_ANTIPODE,
            "recovered antipode differs from the receipt or the carrier antipode")
    require(rec["rule_class"] == "seam_mean" and comp["rule"]["class"] == "seam_mean", "rule class is not seam_mean")
    require(rec["event_count"] == 360 == comp["event_count"], "canonical event count changed")
    # Schedule: seam event counts against the uniform law.
    counts = Counter((a, b) for a, b, *_ in rec["table"])
    per_seam = 360 // 30
    chi = sum(Fraction((counts[s] - per_seam) ** 2, per_seam) for s in rec["seams"])
    require(str(chi) == comp["schedule"]["chi_square"], "seam-count chi-square differs")
    require(min(counts.values()) == comp["schedule"]["min_seam_count"]
            and max(counts.values()) == comp["schedule"]["max_seam_count"], "seam-count extremes differ")
    # Conservation, descent of the centered squared norm, seam-form increases, terminal deviation.
    state = [Fraction(v) for v in log["initial"]]
    total = sum(state)
    seams = rec["seams"]

    def centered(x):
        return sum(v * v for v in x) - sum(x) ** 2 / len(x)

    def seam_form(x):
        return sum((x[i] - x[j]) ** 2 for i, j in seams)

    violations = 0
    increases = 0
    stutters = 0
    for a, b, ba, aa, bb, ab in rec["table"]:
        require(state[a] == ba and state[b] == bb, "event before-values do not follow the replayed state")
        before_c, before_s = centered(state), seam_form(state)
        state[a], state[b] = aa, ab
        after_c, after_s = centered(state), seam_form(state)
        if ba != bb:
            require(before_c - after_c == (ba - bb) ** 2 / 2, "centered-norm decrement is not d^2/2")
            if after_c >= before_c:
                violations += 1
        else:
            stutters += 1
        if after_s > before_s:
            increases += 1
    require([str(v) for v in state] == log["final"], "replayed final state differs from the log")
    require(sum(state) == total, "the total is not conserved")
    require(violations == comp["descent"]["strict_descent_violations"] == 0, "strict-descent violations differ")
    require(stutters == comp["descent"]["stutter_events"], "stutter-event count differs")
    require(increases == comp["descent"]["laplacian_form_increase_events"], "seam-form increase count differs")
    mean = total / 12
    max_dev = max(abs(v - mean) for v in state)
    require(abs(float(max_dev) - comp["terminal"]["log_final_max_deviation_from_mean_float"]) < 1e-9,
            "terminal deviation from the mean differs")
    # The probe: T^(2n) e_p with T = I - L/60, at 120 significant digits.
    lap = [[0] * 12 for _ in range(12)]
    for i, j in seams:
        lap[i][i] += 1
        lap[j][j] += 1
        lap[i][j] -= 1
        lap[j][i] -= 1
    m = [[(60 if r == c else 0) - lap[r][c] for c in range(12)] for r in range(12)]
    probe = log["probe"]
    require(probe["steps"] == receipt["parameters"]["probe_steps"] and probe["readback_digits"] == 120,
            "probe steps or resolution changed")
    for n in probe["steps"]:
        power = matrix_power(m, 2 * n)
        den = 60 ** (2 * n)
        for p in range(12):
            expected = {str(r): readback_string(power[r][p], den, 120) for r in range(12) if power[r][p] != 0}
            require(probe["readback"][str(n)][p] == expected, f"probe readback differs at step {n}, port {p}")
    return {"chi_square": str(chi), "seam_form_increases": increases, "probe_steps_verified": probe["steps"]}


def verify_all_logs(receipt: dict, logs: dict[str, dict]) -> int:
    checked = 0
    for name, loop in receipt["loops"].items():
        for stage in loop["stages"]:
            log_name = f"{name}__{stage['role']}.json"
            if log_name not in logs:
                continue
            rec = recover(logs[log_name])
            inv = stage["invariants"]
            for key in ("ports", "pairing_type", "farthest_pairing_distance", "rule_class",
                        "nonconservative_events_present"):
                require(rec[key] == inv[key], f"{log_name}: recovered {key} differs from the loop invariant")
            require(rec["seam_count"] == inv["seams"], f"{log_name}: recovered seam count differs")
            require(rec["degree_sequence"] == inv["degree_sequence"], f"{log_name}: degree sequence differs")
            checked += 1
    require(checked == 15, "fifteen stored logs expected in the loops")
    return checked


def verify() -> dict:
    manifest = verify_manifest()
    receipt = verify_receipt_structure(manifest)
    logs = verify_log_pins(receipt)
    verify_loops(receipt, logs)
    canonical_result = verify_canonical_log(receipt, logs[CANONICAL_LOG])
    checked = verify_all_logs(receipt, logs)
    return {"status": "VERIFIED_CLOSURE_LOOP_ARCHIVE",
            "manifest_sha256": sha256_file(MANIFEST),
            "archived_files_verified": manifest["curated_archive"]["file_count"],
            "loops": {name: loop["receipt"] for name, loop in receipt["loops"].items()},
            "canonical_fixed_point_level": receipt["canonical_loop"]["fixed_point_level"],
            "canonical_log_recomputation": canonical_result,
            "logs_recovered_and_compared": checked,
            "glued_federation": receipt["federation"]["glued_federation"]}


def main() -> int:
    try:
        result = verify()
    except (OSError, UnicodeError, json.JSONDecodeError, KeyError, ValueError, VerificationError) as exc:
        print(json.dumps({"status": "FAIL", "reason": f"{type(exc).__name__}: {exc}"}, sort_keys=True))
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
