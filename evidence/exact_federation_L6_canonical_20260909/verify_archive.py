#!/usr/bin/env python3
"""Fail-closed verifier for an exact-federation evidence archive.

Imports numpy and the standard library only.  Checks every manifest digest,
rebuilds the port-graph components from the seam arrays, recomputes the
expected terminal multisets from the initial loads, checks every archived
terminal vector, replays schedule 0 of the integer law from its seed with
numpy PCG64 and the declared move rule, checks the V ledgers, the float
terminal state, the lambda_2 eigenvector and the isolated kernel limit.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
CONTROL_FILES = {"README.md", "archive_manifest.json", "verify_archive.py"}
SCHEMA = "oph.exact_federation_archive.v1"
PORTS = 12
SEAMS_PER_CARRIER = 30
KERNEL_STEPS = (1, 5, 30, 100, 300)
INTEGER_JSON = {"port_pair": "integer_law_port_pair.json", "isolated": "integer_law_isolated.json"}
INTEGER_NPZ = {"port_pair": "integer_terminal_port_pair.npz", "isolated": "integer_terminal_isolated.npz"}


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False) + "\n"


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("ascii")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def array_sha256(values: np.ndarray, dtype: str) -> str:
    arr = np.ascontiguousarray(np.asarray(values), dtype=np.dtype(dtype).newbyteorder("<"))
    return hashlib.sha256(arr.tobytes()).hexdigest()


def load_json(name: str) -> Any:
    with (ROOT / name).open("r", encoding="ascii") as handle:
        return json.load(handle)


def load_npz(name: str) -> dict[str, np.ndarray]:
    with np.load(ROOT / name, allow_pickle=False) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def close(a: float, b: float, rel: float, absolute: float = 0.0) -> bool:
    return abs(float(a) - float(b)) <= rel * max(abs(float(a)), abs(float(b))) + absolute


# --------------------------------------------------------------------------
# Inventory
# --------------------------------------------------------------------------


def verify_inventory(manifest: dict[str, Any]) -> None:
    expected = {row["path"]: row for row in manifest["inventory"]}
    actual = {path.name for path in ROOT.iterdir() if path.is_file() and path.name not in CONTROL_FILES}
    require(actual == set(expected), "archive inventory names differ")
    lines = []
    total = 0
    for name in sorted(expected):
        row = expected[name]
        path = ROOT / name
        size = path.stat().st_size
        digest = sha256_file(path)
        require(size == row["bytes"], f"byte count mismatch: {name}")
        require(digest == row["sha256"], f"SHA-256 mismatch: {name}")
        total += size
        lines.append(f"{digest}  {size}  {name}\n")
    curated = manifest["curated_archive"]
    require(len(expected) == curated["file_count"], "file count mismatch")
    require(total == curated["total_bytes"], "total byte count mismatch")
    require(hashlib.sha256("".join(lines).encode("utf-8")).hexdigest() == curated["inventory_sha256"], "inventory digest mismatch")
    require(sha256_file(ROOT / "verify_archive.py") == manifest["exact_verification"]["verifier_sha256"], "verifier digest mismatch")


# --------------------------------------------------------------------------
# Geometry, components and expected multisets
# --------------------------------------------------------------------------


def multiset_hash(size: np.ndarray, q: np.ndarray, r: np.ndarray) -> str:
    entries = []
    for m, qq, rr in zip(size.tolist(), q.tolist(), r.tolist()):
        multiset = []
        if m - rr:
            multiset.append([qq, m - rr])
        if rr:
            multiset.append([qq + 1, rr])
        entries.append([m, multiset])
    return sha256_json({"canonicalizer": "component_multiset", "components": entries})


def components_by_union_find(carriers: int, inter: np.ndarray) -> np.ndarray:
    parent = list(range(carriers))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for u, v in zip(inter[:, 0].tolist(), inter[:, 2].tolist()):
        ru, rv = find(u), find(v)
        if ru != rv:
            if ru < rv:
                parent[rv] = ru
            else:
                parent[ru] = rv
    roots = np.asarray([find(i) for i in range(carriers)], dtype=np.int64)
    distinct = np.unique(roots)  # sorted, so components are ordered by their lowest carrier
    relabel = np.empty(carriers, dtype=np.int64)
    relabel[distinct] = np.arange(distinct.size)
    return relabel[roots]


def verify_geometry(config: dict[str, Any], prim: dict[str, np.ndarray]) -> dict[str, Any]:
    template = prim["carrier_seam_template"].astype(np.int64)
    require(template.shape == (SEAMS_PER_CARRIER, 2), "carrier seam template shape")
    require(bool(np.all(template[:, 0] < template[:, 1])) and template.min() >= 0 and template.max() < PORTS, "template ports")
    require(len({(int(i), int(j)) for i, j in template}) == SEAMS_PER_CARRIER, "duplicate template seams")
    require(template.tolist() == config["carrier_seam_template"], "template differs from config")
    require(bool(np.all(np.bincount(template.ravel(), minlength=PORTS) == 5)), "template degree")
    adjacency = np.zeros((PORTS, PORTS), dtype=bool)
    adjacency[template[:, 0], template[:, 1]] = True
    adjacency |= adjacency.T
    reach = np.zeros(PORTS, dtype=bool)
    reach[0] = True
    for _ in range(PORTS):
        reach |= adjacency[reach].any(axis=0)
    require(bool(reach.all()), "template is not connected")
    carriers = int(config["carriers"])
    ports = carriers * PORTS
    require(prim["initial_loads"].shape == (carriers, PORTS), "initial loads shape")
    require(config["ports"] == ports, "port count")
    require(np.array_equal(prim["intra_seam_carrier"], np.repeat(np.arange(carriers), SEAMS_PER_CARRIER)), "intra seam carriers")
    require(np.array_equal(prim["intra_seam_port_a"], np.tile(template[:, 0], carriers)), "intra seam ports a")
    require(np.array_equal(prim["intra_seam_port_b"], np.tile(template[:, 1], carriers)), "intra seam ports b")
    inter = prim["inter_seam"].astype(np.int64)
    m = int(inter.shape[0])
    require(inter.shape == (m, 4), "inter seam shape")
    require(config["seams"]["intra"] == carriers * SEAMS_PER_CARRIER and config["seams"]["inter"] == m, "seam counts")
    require(config["seams"]["total"] == carriers * SEAMS_PER_CARRIER + m, "total seam count")
    require(inter[:, [0, 2]].min() >= 0 and inter[:, [0, 2]].max() < carriers, "inter seam carriers")
    require(inter[:, [1, 3]].min() >= 0 and inter[:, [1, 3]].max() < PORTS, "inter seam ports")
    require(bool(np.all(inter[:, 0] != inter[:, 2])), "inter seam inside one carrier")
    slots = np.concatenate([inter[:, 0] * PORTS + inter[:, 1], inter[:, 2] * PORTS + inter[:, 3]])
    require(np.unique(slots).size == slots.size, "a carrier port is glued more than once")
    glued = np.bincount(np.concatenate([inter[:, 0], inter[:, 2]]), minlength=carriers)
    production = config["gluing"]["production"]
    require(bool(np.all(glued == production["glued_ports_per_carrier"])), "glued ports per carrier")
    usage = np.bincount(np.concatenate([inter[:, 1], inter[:, 3]]), minlength=PORTS)
    require(usage.tolist() == production["port_usage_histogram"], "port usage histogram")
    require(sha256_json(inter.tolist()) == production["port_pairs_sha256"], "port pair digest")
    loads = prim["initial_loads"].astype(np.int64).ravel()
    require(loads.min() >= 0 and loads.max() <= config["seeds"]["loads"]["load_max"], "load range")
    rng = np.random.default_rng(config["seeds"]["loads"]["seed"])
    require(np.array_equal(rng.integers(0, config["seeds"]["loads"]["load_max"] + 1, size=ports), loads), "initial loads differ from the declared seed stream")
    carrier_labels = components_by_union_find(carriers, inter)
    labels = np.repeat(carrier_labels, PORTS)
    require(np.array_equal(labels, prim["component_of_port"].astype(np.int64)), "component labels")
    count = int(labels.max()) + 1
    require(config["components"]["port_pair"] == count, "component count")
    size = np.bincount(labels, minlength=count).astype(np.int64)
    total = np.rint(np.bincount(labels, weights=loads.astype(float), minlength=count)).astype(np.int64)
    q = total // size
    r = total - q * size
    for name, arr in (("component_size", size), ("component_total", total), ("component_q", q), ("component_r", r)):
        require(np.array_equal(prim[name].astype(np.int64), arr), f"{name} differs")
    iso_labels = np.arange(ports, dtype=np.int64) // PORTS
    iso_size = np.full(carriers, PORTS, dtype=np.int64)
    iso_total = loads.reshape(carriers, PORTS).sum(axis=1)
    iso_q = iso_total // PORTS
    iso_r = iso_total - iso_q * PORTS
    for name, arr in (("isolated_component_total", iso_total), ("isolated_component_q", iso_q), ("isolated_component_r", iso_r)):
        require(np.array_equal(prim[name].astype(np.int64), arr), f"{name} differs")
    require(config["components"]["isolated"] == carriers, "isolated component count")
    intra_a = prim["intra_seam_carrier"].astype(np.int64) * PORTS + prim["intra_seam_port_a"].astype(np.int64)
    intra_b = prim["intra_seam_carrier"].astype(np.int64) * PORTS + prim["intra_seam_port_b"].astype(np.int64)
    seam_a = np.concatenate([intra_a, inter[:, 0] * PORTS + inter[:, 1]])
    seam_b = np.concatenate([intra_b, inter[:, 2] * PORTS + inter[:, 3]])
    return {
        "carriers": carriers, "ports": ports, "loads": loads, "template": template,
        "port_pair": {"labels": labels, "size": size, "total": total, "q": q, "r": r, "seam_a": seam_a, "seam_b": seam_b},
        "isolated": {"labels": iso_labels, "size": iso_size, "total": iso_total, "q": iso_q, "r": iso_r, "seam_a": intra_a, "seam_b": intra_b},
    }


# --------------------------------------------------------------------------
# Integer law: terminal vectors, ledgers and the replay of schedule 0
# --------------------------------------------------------------------------


def unit_transfer_decrements(max_mismatch: int) -> np.ndarray:
    """``decrement[d]`` = sum of 2(d_k - 1) over the unit transfers of a descent with mismatch d."""

    table = np.zeros(max_mismatch + 1, dtype=np.int64)
    for d in range(2, max_mismatch + 1):
        total = 0
        dk = d
        while dk >= 2:
            total += 2 * (dk - 1)
            dk -= 2
        require(total == (d * d - (d % 2)) // 2, "closed form of the unit-transfer ledger")
        table[d] = total
    return table


def apply_sweep(x: np.ndarray, seam_a: np.ndarray, seam_b: np.ndarray, seq: np.ndarray, coin: np.ndarray, v: int, v_min: int,
                decrement: np.ndarray, chunk: int = 65536) -> tuple[tuple[int, int, int, int], int, int]:
    """Sequential semantics of one sweep, applied in dependency layers of pairwise disjoint seams."""

    ports = x.size
    descents = swaps = waits = transfers = 0
    first = -1
    for start in range(0, seq.size, chunk):
        ids = seq[start : start + chunk]
        ea = seam_a[ids]
        eb = seam_b[ids]
        ec = coin[start : start + chunk]
        r = ids.size
        delta = np.zeros(r, dtype=np.int64)
        remaining = np.arange(r, dtype=np.int64)
        while remaining.size:
            ra = ea[remaining]
            rb = eb[remaining]
            m = remaining.size
            cat = np.empty(2 * m, dtype=np.int64)
            cat[0::2] = ra
            cat[1::2] = rb
            order = np.argsort(cat, kind="stable")
            sorted_ports = cat[order]
            starts = np.r_[True, sorted_ports[1:] != sorted_ports[:-1]]
            earliest = np.full(ports, m, dtype=np.int64)
            earliest[sorted_ports[starts]] = (order // 2)[starts]
            pos = np.arange(m, dtype=np.int64)
            ready = (earliest[ra] == pos) & (earliest[rb] == pos)
            la = ra[ready]
            lb = rb[ready]
            lc = ec[remaining[ready]]
            xa = x[la]
            xb = x[lb]
            tot = xa + xb
            lo = tot >> 1
            hi = tot - lo
            na = np.where(lc == 1, hi, lo)
            nb = tot - na
            d = xa - xb
            ad = np.abs(d)
            unit = ad == 1
            wait = (d == 0) | (unit & (na == xa))
            swap = unit & ~wait
            desc = ad >= 2
            change = na * na + nb * nb - xa * xa - xb * xb
            require(bool(np.all(change[desc] == -decrement[ad[desc]])), "a descent move violates the unit-transfer decrement identity")
            require(bool(np.all(change[~desc] == 0)), "a wait or swap changed V")
            waits += int(np.count_nonzero(wait))
            swaps += int(np.count_nonzero(swap))
            descents += int(np.count_nonzero(desc))
            transfers += int(np.sum((ad >> 1)[desc]))
            x[la] = na
            x[lb] = nb
            delta[remaining[ready]] = change
            remaining = remaining[~ready]
        running = v + np.cumsum(delta)
        if first < 0:
            hit = np.flatnonzero(running == v_min)
            if hit.size:
                first = start + int(hit[0])
        v = int(running[-1])
    return (descents, swaps, waits, transfers), v, first


def replay_schedule_zero(config: dict[str, Any], geometry: dict[str, Any], gluing: str, schedule: dict[str, Any], archived: np.ndarray) -> None:
    g = geometry[gluing]
    seam_a, seam_b = g["seam_a"], g["seam_b"]
    seams = int(seam_a.size)
    loads = geometry["loads"]
    x = loads.astype(np.int64).copy()
    v = int(np.dot(x, x))
    v_min = int(schedule["V_minimum"])
    decrement = unit_transfer_decrements(int(loads.max() - loads.min()))
    rng = np.random.default_rng(int(schedule["seed"]))
    draw_hashes = schedule.get("draw_sha256_per_sweep")
    require(draw_hashes is not None and len(draw_hashes) == schedule["sweeps"], "schedule 0 carries one draw digest per sweep")
    first_attempt = 0 if v == v_min else -1
    sweep = 0
    while first_attempt < 0 and sweep < schedule["sweeps"]:
        seq = rng.integers(0, seams, size=seams, dtype=np.int64)
        coin = rng.integers(0, 2, size=seams, dtype=np.int64)
        require(
            [array_sha256(seq, "int64"), array_sha256(coin, "int64")] == draw_hashes[sweep],
            f"RNG stream differs at sweep {sweep} (numpy {np.__version__} against archived numpy {config['rng']['numpy_version']})",
        )
        counts, v, first = apply_sweep(x, seam_a, seam_b, seq, coin, v, v_min, decrement)
        if first >= 0:
            first_attempt = sweep * seams + first + 1
        sweep += 1
        row = schedule["per_sweep"][sweep - 1]
        require(v == schedule["V_ledger"][sweep], f"replayed V differs at sweep {sweep}")
        require(v == int(np.dot(x, x)), f"replayed ledger differs from the state at sweep {sweep}")
        require(counts == (row["descents"], row["swaps"], row["waits"], row["unit_transfers"]), f"replayed move counts differ at sweep {sweep}")
    require(sweep == schedule["sweeps"] and first_attempt == schedule["attempts_to_balanced_class"], "replayed termination differs")
    require(np.array_equal(x, archived.astype(np.int64)), "replayed terminal vector differs from the archived vector")
    require(array_sha256(x.astype(np.int8), "int8") == schedule["terminal_sha256"], "replayed terminal digest differs")


def verify_integer_law(config: dict[str, Any], geometry: dict[str, Any], gluing: str, *, replay: bool) -> dict[str, Any]:
    block = load_json(INTEGER_JSON[gluing])
    arrays = load_npz(INTEGER_NPZ[gluing])
    g = geometry[gluing]
    labels, size, q, r = g["labels"], g["size"], g["q"], g["r"]
    seam_a, seam_b = g["seam_a"], g["seam_b"]
    seams = int(seam_a.size)
    loads = geometry["loads"]
    schedules = block["schedules"]
    require(block["gluing"] == gluing and block["seams"] == seams, "integer block identity")
    require(len(schedules) == config["schedule_count"], "integer schedule count")
    require([s["seed"] for s in schedules] == config["seeds"]["schedules"][gluing], "integer schedule seeds")
    expected = multiset_hash(size, q, r)
    require(block["expected_quotient_hash"] == expected, "expected multiset hash")
    v0 = int(np.dot(loads, loads))
    v_min = int(np.sum((size - r) * q * q + r * (q + 1) * (q + 1)))
    require(block["V_initial"] == v0 and block["V_minimum"] == v_min, "V endpoints of the integer block")
    dmax = int(loads.max() - loads.min())
    decrement = unit_transfer_decrements(dmax)
    q_of = q[labels]
    for index, s in enumerate(schedules):
        tag = f"{gluing} schedule {index}"
        require(s["terminated"] is True, f"{tag} did not terminate")
        require(s["unit_transfer_decrement_identity_violations"] == 0 and s["strict_descent_violations"] == 0, f"{tag} violations")
        require(s["quotient_hash"] == expected and s["quotient_hash_equals_expected"] is True, f"{tag} quotient hash")
        require(s["V_initial"] == v0 and s["V_minimum"] == v_min and s["V_terminal"] == v_min, f"{tag} V endpoints")
        ledger = s["V_ledger"]
        rows = s["per_sweep"]
        require(len(ledger) == s["sweeps"] + 1 and len(rows) == s["sweeps"], f"{tag} ledger length")
        require(ledger[0] == v0 and ledger[-1] == v_min, f"{tag} ledger endpoints")
        require(s["attempts"] == s["sweeps"] * seams, f"{tag} attempts")
        require(sum(row["descents"] for row in rows) == s["descents"], f"{tag} descent total")
        require(sum(row["swaps"] for row in rows) == s["swaps"], f"{tag} swap total")
        require(sum(row["waits"] for row in rows) == s["waits"], f"{tag} wait total")
        require(sum(row["unit_transfers"] for row in rows) == s["unit_transfers"], f"{tag} unit transfer total")
        require(sum(row["violations"] for row in rows) == 0, f"{tag} per-sweep violations")
        for k, row in enumerate(rows):
            require(row["descents"] + row["swaps"] + row["waits"] == seams, f"{tag} sweep {k + 1} attempt count")
            drop = ledger[k] - ledger[k + 1]
            require(drop >= 0, f"{tag} sweep {k + 1} raises V")
            require((2 * row["unit_transfers"] <= drop <= 2 * (dmax - 1) * row["unit_transfers"]) if dmax >= 2 else drop == 0, f"{tag} sweep {k + 1} ledger bounds")
            require(row["unit_transfers"] >= row["descents"], f"{tag} sweep {k + 1} transfers below descents")
        first = s["attempts_to_balanced_class"]
        require((first == 0 and s["sweeps"] == 0) or ((s["sweeps"] - 1) * seams < first <= s["sweeps"] * seams), f"{tag} terminating attempt lies outside the last sweep")
        x8 = arrays[s["terminal_array"]]
        require(x8.dtype == np.int8 and x8.shape == (geometry["ports"],), f"{tag} terminal array shape")
        require(array_sha256(x8, "int8") == s["terminal_sha256"], f"{tag} terminal digest")
        x = x8.astype(np.int64)
        require(x.min() >= loads.min() and x.max() <= loads.max(), f"{tag} terminal range")
        require(int(np.dot(x, x)) == v_min, f"{tag} terminal V")
        require(np.array_equal(np.rint(np.bincount(labels, weights=x.astype(float), minlength=size.size)).astype(np.int64), g["total"]), f"{tag} conservation")
        require(bool(np.all((x == q_of) | (x == q_of + 1))), f"{tag} terminal reading outside the balanced class")
        high = np.bincount(labels, weights=(x == q_of + 1).astype(float), minlength=size.size)
        require(np.array_equal(np.rint(high).astype(np.int64), r), f"{tag} balanced split")
        low = size - np.rint(high).astype(np.int64)
        require(multiset_hash(size, q, np.rint(high).astype(np.int64)) == expected and bool(np.all(low + np.rint(high).astype(np.int64) == size)), f"{tag} terminal multiset hash")
        d = np.abs(x[seam_a] - x[seam_b])
        require(int(d.max()) <= 1 and int(d.max()) == s["max_seam_difference_at_termination"], f"{tag} seam difference")
        require(int(np.count_nonzero(d == 1)) == s["odd_tie_seams_at_termination"], f"{tag} odd-tie seam count")
    require(block["unique_quotient_hash_count"] == 1 and block["quotient_hash_equals_expected_all"] is True, "integer block verdict")
    if replay:
        started = time.perf_counter()
        replay_schedule_zero(config, geometry, gluing, schedules[0], arrays[schedules[0]["terminal_array"]])
        block["_replay_seconds"] = time.perf_counter() - started
    return block


# --------------------------------------------------------------------------
# Mean law: the budgeted float schedule and the isolated control
# --------------------------------------------------------------------------


def verify_float_law(config: dict[str, Any], geometry: dict[str, Any]) -> dict[str, Any] | None:
    if not (ROOT / "float_law_port_pair.json").exists():
        return None
    block = load_json("float_law_port_pair.json")
    state = load_npz("float_terminal_port_pair.npz")[block["terminal_state_array"]]
    g = geometry["port_pair"]
    seam_a, seam_b, labels = g["seam_a"], g["seam_b"], g["labels"]
    seams = int(seam_a.size)
    loads = geometry["loads"]
    require(block["gluing"] == "port_pair" and block["seams"] == seams, "float block identity")
    require(block["seed"] == config["seeds"]["schedules"]["float_port_pair"][0], "float seed")
    require(block["sweeps"] == config["budgets"]["float_port_pair_sweeps"] and block["attempts"] == block["sweeps"] * seams, "float budget")
    require(state.dtype == np.float64 and state.shape == (geometry["ports"],), "float terminal state shape")
    require(array_sha256(state, "float64") == block["terminal_state_sha256"], "float terminal digest")
    require(bool(np.all(np.isfinite(state))), "float terminal state is not finite")
    d0 = loads[seam_a] - loads[seam_b]
    phi0 = int(np.dot(d0, d0))
    v0 = int(np.dot(loads, loads))
    require(close(block["Phi_initial"], phi0, 1e-9) and close(block["V_initial"], v0, 1e-9), "float initial Phi and V")
    d = state[seam_a] - state[seam_b]
    phi = float(np.dot(d, d))
    v = float(np.dot(state, state))
    require(close(block["Phi_terminal"], phi, 1e-6), "float terminal Phi")
    require(close(block["V_terminal"], v, 1e-9), "float terminal V")
    require(close(block["total_terminal"], float(np.sum(state)), 1e-9) and close(float(np.sum(state)), float(np.sum(loads)), 1e-9), "float conservation")
    count = int(labels.max()) + 1
    totals = np.bincount(labels, weights=loads.astype(float), minlength=count)
    sizes = np.bincount(labels, minlength=count)
    v_min = float(np.sum(totals * totals / sizes))
    require(close(block["V_minimum"], v_min, 1e-9), "float V minimum")
    means = np.bincount(labels, weights=state, minlength=count) / sizes
    deviation = float(np.max(np.abs(state - means[labels])))
    require(close(block["max_abs_deviation_from_component_mean"], deviation, 1e-5, 1e-9), "float deviation from the component mean")
    rows = block["per_sweep"]
    require(len(rows) == block["sweeps"], "float per-sweep rows")
    previous = float(v0)
    ledger = 0.0
    waits = violations = raising = 0
    for k, row in enumerate(rows):
        require(row["sweep"] == k + 1, "float sweep index")
        require(row["V"] <= previous * (1 + 1e-12), f"float V rises at sweep {k + 1}")
        require(abs((previous - row["V"]) - row["ledger_increment"]) <= 1e-9 * v0 + 1e-6 * abs(row["ledger_increment"]), f"float ledger arithmetic at sweep {k + 1}")
        require(row["violations"] == 0, f"float strict-descent violation at sweep {k + 1}")
        require(0 <= row["waits"] <= seams, f"float waits at sweep {k + 1}")
        previous = row["V"]
        ledger += row["ledger_increment"]
        waits += row["waits"]
        violations += row["violations"]
        raising += row["phi_raising_moves"]
    require(close(previous, v, 1e-9), "float ledger end differs from the terminal state")
    require(abs((v0 - v) - ledger) <= 1e-8 * v0, "float total ledger")
    require(block["waits"] == waits and block["strict_descent_violations"] == violations == 0 and block["phi_raising_moves"] == raising, "float totals")
    require(block["non_wait_moves"] == block["attempts"] - waits, "float non-wait moves")
    require(block["terminated"] is (phi < block["threshold"]), "float termination flag")
    return block


def verify_float_isolated(config: dict[str, Any], geometry: dict[str, Any]) -> dict[str, Any]:
    block = load_json("float_law_isolated.json")
    g = geometry["isolated"]
    carriers = geometry["carriers"]
    schedules = block["schedules"]
    require(block["gluing"] == "isolated" and block["seams"] == int(g["seam_a"].size), "isolated float identity")
    require(len(schedules) == config["schedule_count"], "isolated float schedule count")
    require([s["seed"] for s in schedules] == config["seeds"]["schedules"]["float_isolated"], "isolated float seeds")
    payload = {
        "canonicalizer": "component_lattice_snap",
        "component_sizes": g["size"].tolist(),
        "component_of_port": g["labels"].tolist(),
        "q": g["total"][g["labels"]].tolist(),
    }
    expected = sha256_json(payload)
    require(block["expected_terminal_quotient_hash"] == expected, "isolated float expected hash")
    for index, s in enumerate(schedules):
        tag = f"isolated float schedule {index}"
        require(s["terminated"] is True and s["phi_terminal"] < block["threshold"], f"{tag} termination")
        require(s["terminal_quotient_hash"] == expected, f"{tag} terminal hash")
        require(s["strict_descent_violations"] == 0, f"{tag} violations")
        require(s["descent_ledger_relative_error_below_1e-9"] is True, f"{tag} ledger")
        require(s["component_mean_within_1e-9"] is True and s["max_abs_deviation_from_component_mean"] < 1e-9, f"{tag} deviation")
        require(s["lattice_residual_max"] < 0.5, f"{tag} lattice snap margin")
        require(s["attempts"] == s["sweeps"] * block["seams"], f"{tag} attempts")
    require(block["unique_terminal_hash_count"] == 1 and block["terminal_hash_equals_expected_all"] is True and block["all_terminated"] is True, "isolated float verdict")
    require(block["components"] == carriers, "isolated float components")
    return block


# --------------------------------------------------------------------------
# Spectrum
# --------------------------------------------------------------------------


def carrier_laplacian(template: np.ndarray) -> np.ndarray:
    lap = np.zeros((PORTS, PORTS))
    for i, j in template.tolist():
        lap[i, j] -= 1.0
        lap[j, i] -= 1.0
        lap[i, i] += 1.0
        lap[j, j] += 1.0
    return lap


def verify_spectrum(config: dict[str, Any], geometry: dict[str, Any]) -> dict[str, Any] | None:
    template = geometry["template"]
    iso_values = np.linalg.eigvalsh(carrier_laplacian(template))
    require(abs(iso_values[1] - (5.0 - 5.0**0.5)) < 1e-12, "isolated gap differs from 5 - sqrt(5)")
    if not (ROOT / "spectrum_port_pair.json").exists():
        return None
    block = load_json("spectrum_port_pair.json")
    vec = load_npz("fiedler_port_pair.npz")[block["eigenvector_array"]]
    g = geometry["port_pair"]
    seam_a, seam_b, labels = g["seam_a"], g["seam_b"], g["labels"]
    ports = geometry["ports"]
    require(vec.dtype == np.float64 and vec.shape == (ports,), "eigenvector shape")
    require(array_sha256(vec, "float64") == block["eigenvector_sha256"], "eigenvector digest")
    norm = float(np.linalg.norm(vec))
    require(abs(norm - 1.0) < 1e-9, "eigenvector norm")
    degree = np.bincount(seam_a, minlength=ports) + np.bincount(seam_b, minlength=ports)
    lv = degree * vec - (np.bincount(seam_a, weights=vec[seam_b], minlength=ports) + np.bincount(seam_b, weights=vec[seam_a], minlength=ports))
    lam = float(block["laplacian_lambda_2"])
    residual = float(np.linalg.norm(lv - lam * vec))
    require(residual < 1e-8, f"eigenvector residual {residual}")
    dv = vec[seam_a] - vec[seam_b]
    rayleigh = float(np.dot(dv, dv))
    require(close(rayleigh, lam, 1e-8, 1e-14), "Rayleigh quotient differs from lambda_2")
    require(close(block["eigenvector_rayleigh_quotient"], rayleigh, 1e-8, 1e-14), "archived Rayleigh quotient")
    count = int(labels.max()) + 1
    require(float(np.max(np.abs(np.bincount(labels, weights=vec, minlength=count)))) < 1e-8, "eigenvector is not orthogonal to the component indicators")
    denominator = 2 * seam_a.size / geometry["carriers"]
    require(close(block["D"], denominator, 1e-12), "D")
    require(close(block["t_fed_second_eigenvalue"], 1.0 - lam / denominator, 1e-9), "second eigenvalue of T_fed")
    require(lam > 0 and block["laplacian_lambda_max"] >= lam, "spectral ordering")
    require(block["laplacian_lambda_max"] <= 2 * float(degree.max()) + 1e-9, "lambda_max above the degree bound")
    require(close(block["isolated_lambda_2"], iso_values[1], 1e-9), "isolated lambda_2")
    block["_residual"] = residual
    return block


# --------------------------------------------------------------------------
# Response kernels
# --------------------------------------------------------------------------


def slow_band_projector(template: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(carrier_laplacian(template))
    band = np.abs(values - (5.0 - 5.0**0.5)) < 1e-9
    require(int(np.count_nonzero(band)) == 3, "slow band multiplicity")
    v = vectors[:, band]
    return v @ v.T


def isolated_kernels_by_propagation(prim: dict[str, np.ndarray], cell: int) -> dict[int, np.ndarray]:
    """Direct matrix powers with rescaling on the carrier's own block, from the archived intra seams."""

    rows = slice(cell * SEAMS_PER_CARRIER, (cell + 1) * SEAMS_PER_CARRIER)
    require(bool(np.all(prim["intra_seam_carrier"][rows] == cell)), "intra seam rows of the sample carrier")
    template = np.stack([prim["intra_seam_port_a"][rows], prim["intra_seam_port_b"][rows]], axis=1).astype(np.int64)
    lap = carrier_laplacian(template)
    q = np.eye(PORTS) - np.ones((PORTS, PORTS)) / PORTS
    y = q.copy()
    out = {}
    target = {2 * n: n for n in KERNEL_STEPS}
    for step in range(1, max(target) + 1):
        y = y - (lap @ y) / 60.0
        scale = np.max(np.abs(y))
        if scale > 0:
            y /= scale
        if step in target:
            c = q @ y
            c = 0.5 * (c + c.T)
            out[target[step]] = PORTS * c / np.trace(c)
    return out


def verify_kernels(config: dict[str, Any], geometry: dict[str, Any], prim: dict[str, np.ndarray]) -> dict[str, Any]:
    block = load_json("kernel_readout.json")
    mats = load_npz("kernel_matrices.npz")
    cells = [int(c) for c in block["cells"]]
    require(cells == config["seeds"]["kernel_sample"]["cells"] and mats["cells"].tolist() == cells, "kernel sample cells")
    require(mats["steps"].tolist() == list(KERNEL_STEPS) and block["steps"] == list(KERNEL_STEPS), "kernel steps")
    for name in ("cells", "steps", "glued", "isolated"):
        entry = block["matrices"][name]
        require(mats[name].shape == tuple(entry["shape"]) and str(mats[name].dtype) == entry["dtype"], f"kernel array {name} shape")
        require(array_sha256(mats[name], entry["dtype"]) == entry["sha256"], f"kernel array {name} digest")
    require(mats["glued"].shape == (len(cells), len(KERNEL_STEPS), PORTS, PORTS), "glued kernel shape")
    p_slow = slow_band_projector(geometry["template"])
    gram = 4.0 * p_slow
    pentagonal = prim["pentagonal_cell"]
    vertex_cells = prim["vertex_cells"].astype(np.int64)
    require(vertex_cells.shape == (12, 5), "vertex cell table")
    covered = sorted({int(v) for v in range(12) for c in cells if c in vertex_cells[v].tolist()})
    require(covered == block["vertices_covered"], "vertices covered by the sample")
    require(block["pentagonal_cells_in_sample"] == int(sum(bool(pentagonal[c]) for c in cells)), "pentagonal cells in the sample")
    iso_dev = 0.0
    gram_dev = 0.0
    for i, entry in enumerate(block["per_cell"]):
        cell = int(entry["cell"])
        require(cell == cells[i] and entry["pentagonal"] == bool(pentagonal[cell]), f"per-cell identity {i}")
        for kind in ("glued", "isolated"):
            for j, (n, row) in enumerate(zip(KERNEL_STEPS, entry[kind])):
                k = mats[kind][i, j]
                require(row["n"] == n, "kernel step order")
                require(float(np.max(np.abs(k - k.T))) < 1e-12, f"{kind} kernel not symmetric (cell {cell}, n {n})")
                require(abs(float(np.trace(k)) - PORTS) < 1e-9, f"{kind} kernel trace (cell {cell}, n {n})")
                eig = np.sort(np.linalg.eigvalsh(k))[::-1]
                require(float(eig[-1]) > -1e-9, f"{kind} kernel not positive semidefinite (cell {cell}, n {n})")
                require(np.allclose(eig, np.asarray(row["eigenvalues"]), atol=1e-8, rtol=0), f"{kind} eigenvalues (cell {cell}, n {n})")
                share = float(np.trace(p_slow @ k @ p_slow) / np.trace(k))
                require(abs(share - row["slow_band_share"]) < 1e-8 and -1e-9 <= share <= 1 + 1e-9, f"{kind} slow-band share (cell {cell}, n {n})")
                require(int(np.linalg.matrix_rank(k, tol=1e-9)) == row["rank_1e-9"], f"{kind} rank (cell {cell}, n {n})")
                if kind == "isolated" and n == 300:
                    gram_dev = max(gram_dev, float(np.max(np.abs(k - gram))))
        require(float(np.max(np.abs(mats["isolated"][i] - mats["isolated"][0]))) < 1e-12, f"isolated kernels differ between carriers (cell {cell})")
    require(gram_dev < 1e-12, f"isolated K_300 differs from 4 P_slow by {gram_dev}")
    samples = [cells[0], cells[-1]] if len(cells) > 1 else [cells[0]]
    for cell in samples:
        i = cells.index(cell)
        recomputed = isolated_kernels_by_propagation(prim, cell)
        for j, n in enumerate(KERNEL_STEPS):
            dev = float(np.max(np.abs(recomputed[n] - mats["isolated"][i, j])))
            iso_dev = max(iso_dev, dev)
            require(dev < 1e-9, f"recomputed isolated kernel differs (cell {cell}, n {n}): {dev}")
            eig = np.sort(np.linalg.eigvalsh(recomputed[n]))[::-1]
            require(np.allclose(eig, np.asarray(block["per_cell"][i]["isolated"][j]["eigenvalues"]), atol=1e-8, rtol=0), f"recomputed isolated eigenvalues (cell {cell}, n {n})")
        limit = np.sort(np.linalg.eigvalsh(recomputed[300]))[::-1]
        require(np.allclose(limit[:3], 4.0, atol=1e-12) and np.allclose(limit[3:], 0.0, atol=1e-12), f"isolated limit eigenvalues (cell {cell})")
    block["_samples"] = samples
    block["_recompute_deviation"] = iso_dev
    block["_gram_deviation"] = gram_dev
    return block


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def main() -> int:
    started = time.perf_counter()
    try:
        manifest = load_json("archive_manifest.json")
        require(manifest["schema"] == SCHEMA, "archive schema")
        config = load_json("config.json")
        require(config["schema"] == SCHEMA and config["archive_id"] == manifest["archive_id"], "config identity")
        verify_inventory(manifest)
        prim = load_npz("primitives.npz")
        for name, entry in config["primitive_arrays"]["arrays"].items():
            require(name in prim and list(prim[name].shape) == entry["shape"] and str(prim[name].dtype) == entry["dtype"], f"primitive array {name}")
            require(array_sha256(prim[name], entry["dtype"]) == entry["sha256"], f"primitive array digest {name}")
        geometry = verify_geometry(config, prim)
        snapshot = manifest["run_snapshot"]
        require(snapshot["level"] == config["level"] and snapshot["carriers"] == geometry["carriers"], "snapshot identity")
        integer = verify_integer_law(config, geometry, "port_pair", replay=True)
        isolated = verify_integer_law(config, geometry, "isolated", replay=False)
        require(snapshot["integer_law_port_pair"]["expected_quotient_hash"] == integer["expected_quotient_hash"], "snapshot integer hash")
        require(snapshot["integer_law_isolated"]["expected_quotient_hash"] == isolated["expected_quotient_hash"], "snapshot isolated hash")
        flt = verify_float_law(config, geometry)
        flt_iso = verify_float_isolated(config, geometry)
        spectrum = verify_spectrum(config, geometry)
        kernels = verify_kernels(config, geometry, prim)
        if spectrum is not None:
            require(close(snapshot["spectrum_port_pair"]["laplacian_lambda_2"], spectrum["laplacian_lambda_2"], 1e-12), "snapshot lambda_2")
        require(manifest["known_boundaries"]["gluing_is_declared_convention"] is True and manifest["known_boundaries"]["physical_identification"] is False, "boundary flags")
    except (OSError, KeyError, TypeError, ValueError, IndexError, VerificationError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    seconds = time.perf_counter() - started
    print(
        f"PASS: {geometry['carriers']} carriers, {int(geometry['port_pair']['seam_a'].size)} seams, {int(geometry['port_pair']['size'].size)} component(s); "
        f"inventory and primitive digests verified; integer law port_pair {len(integer['schedules'])} schedules in "
        f"{integer['sweeps']['min']}..{integer['sweeps']['max']} sweeps with one multiset hash {integer['expected_quotient_hash'][:16]} "
        f"(schedule 0 replayed from seed {integer['schedules'][0]['seed']} in {integer['_replay_seconds']:.1f} s, terminal vector identical); "
        f"integer law isolated {len(isolated['schedules'])} schedules verified; "
        + (f"float port_pair {flt['sweeps']} sweeps, Phi {flt['Phi_initial']:.6g} -> {flt['Phi_terminal']:.6g}, V ledger exact to 1e-9, 0 violations; " if flt else "float port_pair absent; ")
        + f"float isolated {len(flt_iso['schedules'])} schedules terminated at the component mean; "
        + (f"lambda_2 {spectrum['laplacian_lambda_2']:.6g} with eigenvector residual {spectrum['_residual']:.2e}; " if spectrum else "lambda_2 absent; ")
        + f"kernels: {len(kernels['cells'])} cells, isolated limit recomputed for carriers {kernels['_samples']} within {kernels['_recompute_deviation']:.1e}, "
        f"K_300 = 4 P_slow within {kernels['_gram_deviation']:.1e}; {seconds:.1f} s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
