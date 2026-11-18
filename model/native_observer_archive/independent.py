"""Sparse reconstruction from the written OPH-NOA-1 contract.

Run this file directly: no model initializer or production module is needed.
All amplitudes in basis columns are the common exact symbol 1/sqrt(12).
"""

from collections import Counter, defaultdict
from contextlib import contextmanager
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib.abc
import json
from math import comb, isqrt, sqrt
from pathlib import Path
import random
import sys
from typing import NamedTuple


HERE = Path(__file__).resolve().parent
EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
ADDRESSES = tuple(p for a, b in EDGES for p in ((a, b), (b, a)))
BACKGROUNDS = ((0, 0, 0, 0), (1, 0, 0, 0))
VARIANTS = ("native", "uphill", "no_communication", "no_feedback", "no_publish")


def require(condition, predicate):
    if not condition:
        raise ValueError(predicate)


def integer(value, upper):
    require(type(value) is int and 0 <= value < upper, "integer domain")


def admit():
    contract = json.loads((HERE / "observer_independent_contract.json").read_text())
    raw = (HERE / "protocol.json").read_bytes()
    require(hashlib.sha256(raw).hexdigest() == contract["audited_sha256"][
        "model/native_observer_archive/protocol.json"], "protocol custody")
    p = json.loads(raw)
    for key, expected in {
        "rail": {"sites": 4096, "offset": 2048, "order": 64},
        "window": [256, 512], "read_times": [0, 256, 384, 512, 5120],
        "packet_directions": [1, -1], "sources": [0, 1, 85],
        "backgrounds": [list(b) for b in BACKGROUNDS], "variants": list(VARIANTS),
        "second_address_preparation": 0,
        "literal_check": {"sites": 48, "offset": 8, "order": 4, "time": 24},
    }.items():
        require(p[key] == expected, "declared protocol " + key)
    return contract, p


@lru_cache(maxsize=8192)
def loads(word, background):
    q = list(background)
    for e, (left, right) in enumerate(EDGES):
        flux = ((word >> (2 * e + 1)) & 1) - ((word >> (2 * e)) & 1)
        q[left] += flux
        q[right] -= flux
    return tuple(q)


def event(word, address, background, law="native"):
    integer(word, 4096)
    integer(address, 12)
    require(background in BACKGROUNDS and law in ("native", "uphill"), "event domain")
    left, right = ADDRESSES[address]
    q = loads(word, background)
    difference = q[left] - q[right]
    if law == "uphill":
        difference = -difference
    old = (word >> address) & 1
    current = int(not old and difference > 0) - int(old == 1 and difference < 0)
    return word ^ ((1 << address) if current else 0), old, current


def reversible(word, address, memory, background, law="native", inverse=False):
    integer(memory, 2)
    if not inverse:
        target, old, _ = event(word, address, background, law)
        return target ^ (memory << address), old
    # The output memory is the old source bit, irrespective of its successor.
    old_word = (word & ~(1 << address)) | (memory << address)
    native, _, _ = event(old_word, address, background, law)
    return old_word, ((word ^ native) >> address) & 1


class Registers(NamedTuple):
    sources: tuple
    addresses: tuple
    couriers: tuple
    pa: int
    pb: int
    memories: tuple
    action: int


def prepared(word, address):
    integer(word, 4096)
    integer(address, 12)
    return Registers((word, 0, 0, 0, 0), (address, 0, 0, 0), (0, 0, 0),
                     0, 0, (0, 0), 0)


def unpack(data):
    result = dict(zip(("s0", "s1", "s2", "s3", "s4"), data.sources))
    result.update(zip(("a0", "a1", "a2", "a3"), data.addresses))
    result.update(zip(("c0", "c1", "c2"), data.couriers))
    result.update(pa=data.pa, pb=data.pb, m0=data.memories[0],
                  m1=data.memories[1], action=data.action)
    return result


def pack(r):
    return Registers(tuple(r[f"s{i}"] for i in range(5)),
                     tuple(r[f"a{i}"] for i in range(4)),
                     tuple(r[f"c{i}"] for i in range(3)), r["pa"], r["pb"],
                     (r["m0"], r["m1"]), r["action"])


def instructions(stage, variant):
    require(variant in VARIANTS, "program variant")
    integer(stage, 4)
    move = lambda a, b: ("swap", a, b)
    if stage == 0:
        return [("event", "s0", "a0", "m0"), move("s0", "s1"), move("a0", "a1")]
    if stage == 1:
        return ([] if variant == "no_publish" else [("publish",)]) + (
            [] if variant == "no_communication" else [("xor", "c0", "pa")]) + [
                move("s1", "s2"), move("a1", "a2"), move("c0", "c1")]
    if stage == 2:
        return [move("s2", "s3"), move("a2", "a3"), move("c1", "c2")]
    return [("xor", "pb", "c2")] + (
        [] if variant == "no_feedback" else [("add",)]) + [
            ("event", "s3", "action", "m1"), move("s3", "s4")]


def gate(data, stage, background, variant="native", inverse=False):
    r = unpack(data)
    ops = instructions(stage, variant)
    for op in reversed(ops) if inverse else ops:
        kind, *args = op
        if kind == "swap":
            a, b = args
            r[a], r[b] = r[b], r[a]
        elif kind == "xor":
            a, b = args
            r[a] ^= r[b]
        elif kind == "publish":
            r["pa"] ^= 2 + ((r["s1"] >> r["a1"]) & 1)
        elif kind == "add":
            amount = 6 * (r["pb"] & 1) * ((r["pb"] >> 1) & 1)
            r["action"] = (r["action"] + (-amount if inverse else amount)) % 12
        elif kind == "event":
            s, a, m = args
            r[s], r[m] = reversible(r[s], r[a], r[m], background,
                                    "uphill" if variant == "uphill" else "native", inverse)
        else:
            raise ValueError("unknown instruction")
    return pack(r)


def prefixes(word, address, background, variant="native"):
    out = [prepared(word, address)]
    for stage in range(4):
        out.append(gate(out[-1], stage, background, variant))
    return tuple(out)


def identity(data):
    """Injective mixed-radix encoding of all registers, including blanks."""
    value = 0
    for entries, radix in ((data.sources, 4096), (data.addresses, 12),
                           (data.couriers, 4), ((data.pa, data.pb), 4),
                           (data.memories, 2), ((data.action,), 12)):
        for entry in entries:
            integer(entry, radix)
            value = value * radix + entry
    return value


def pointer(data):
    return data.pa, data.pb


def dirty_states():
    rng = random.Random(12041203)
    for _ in range(256):
        yield Registers(tuple(rng.randrange(4096) for _ in range(5)),
                        tuple(rng.randrange(12) for _ in range(4)),
                        tuple(rng.randrange(4) for _ in range(3)),
                        rng.randrange(4), rng.randrange(4),
                        (rng.randrange(2), rng.randrange(2)), rng.randrange(12))


def tv(left, right):
    require(sum(left.values()) == sum(right.values()) == 12, "complete source law")
    return F(sum(abs(left[x] - right[x]) for x in left.keys() | right.keys()), 24)


def clock_bound(protocol):
    n, a, m = (protocol["rail"][k] for k in ("sites", "offset", "order"))
    lo, hi = protocol["window"]
    coefficients = [comb(m, j) for j in range(m + 1)]
    norm = sum(x * x for x in coefficients)
    require(norm == comb(2 * m, m), "binomial normalization")
    mean = sum(F(j * b * b, norm) for j, b in enumerate(coefficients))
    variance = sum(F(b * b, norm) * (j - mean) ** 2 for j, b in enumerate(coefficients))
    overlaps = [F(sum(coefficients[j] * coefficients[j + d]
                      for j in range(m + 1 - d)), norm) for d in (1, 2)]
    velocity = overlaps[0]
    vv = (1 + overlaps[1]) / 2 - velocity ** 2
    require(mean == F(m, 2) and variance == F(m * m, 4 * (2 * m - 1)), "position moments")
    require(overlaps == [F(m, m + 1), F(m * (m - 1), (m + 1) * (m + 2))], "overlap moments")
    covariance = sum(F(coefficients[j] * coefficients[j + 1], norm) *
                     (F(2 * j + 1, 2) - mean) for j in range(m))
    require(covariance == 0, "symmetrized covariance")
    gap0 = m + 3 - mean
    require(gap0 > 0 and velocity * lo > gap0 and vv > 0, "window monotonicity hypotheses")
    # d/dt [(variance+vv*t^2)/(velocity*t-gap0)^2]
    # = -2*(vv*gap0*t+velocity*variance)/(velocity*t-gap0)^3 < 0.
    spread = variance + vv * lo ** 2
    q = spread / (spread + (velocity * lo - gap0) ** 2)
    distance = min(a + 1, n - a - m)
    ratio = F(3 * hi, distance)
    require(ratio < 1 and hi < distance + 1, "boundary tail hypotheses")
    tail = 2 * ratio ** distance / (1 - F(hi, distance + 1))
    power = tail.denominator.bit_length() - tail.numerator.bit_length() - 1
    dyadic = F(1, 2 ** power)
    require(tail <= dyadic, "boundary enclosure")
    # Integer square-root interval, independent of production's error formula.
    scale = 2 ** 200
    lower = F(isqrt(norm * scale * scale), scale)
    upper = lower + F(1, scale)
    require(lower ** 2 <= norm < upper ** 2, "root enclosure")
    binary = [F.from_float(b / sqrt(norm)) for b in coefficients]
    epsilon = sum(max(abs(x - F(b, 1) / lower), abs(x - F(b, 1) / upper))
                  for b, x in zip(coefficients, binary))
    require(q + dyadic + 2 * epsilon < F(1, 1000), "combined readiness bound")
    return {"normalization": str(norm), "mean": str(mean), "variance": str(variance),
            "overlaps": list(map(str, overlaps)), "velocity_variance": str(vv),
            "cantelli": str(q), "boundary_distance": distance,
            "boundary_ratio": str(ratio), "boundary_dyadic_power": power,
            "binary_vector_error_upper": float(epsilon),
            "combined_failure_upper": float(q + dyadic + 2 * epsilon)}


def coherence_witness():
    """Declared exploratory entangled witness; unnormalized public blocks."""
    # determinant 1-i is nonzero: the source and reference are entangled.
    coefficients = {(0, 0): .5, (0, 1): .5j, (3, 0): .5, (3, 1): .5}
    blocks = defaultdict(dict)
    for (source, reference), coefficient in coefficients.items():
        for a in range(12):
            data = prefixes(source, a, BACKGROUNDS[0])[-1]
            blocks[pointer(data)][data, reference] = coefficient / sqrt(12)
    trace = 0.0
    reference_density = defaultdict(complex)
    maximum = (0.0, None)
    removed_squared = 0.0
    for outcome, vector in blocks.items():
        trace += sum(abs(v) ** 2 for v in vector.values())
        for (left, r), x in vector.items():
            for (right, s), y in vector.items():
                entry = x * y.conjugate()
                if left == right:
                    reference_density[r, s] += entry
                # Private archive dephasing, preserving pointers and source.
                archive_l = (left.addresses, left.memories, left.action, left.couriers)
                archive_r = (right.addresses, right.memories, right.action, right.couriers)
                if archive_l != archive_r:
                    removed_squared += abs(entry) ** 2
                    if abs(entry) > maximum[0]:
                        maximum = (abs(entry), {"outcome": outcome,
                            "left": [left._asdict(), r], "right": [right._asdict(), s],
                            "coherence": [entry.real, entry.imag]})
    expected_reference = {(0, 0): .5, (1, 1): .5, (0, 1): .25 - .25j, (1, 0): .25 + .25j}
    require(abs(trace - 1) < 2e-14, "public instrument normalization")
    require(all(abs(reference_density[k] - v) < 2e-14 for k, v in expected_reference.items()),
            "untouched reference")
    require(maximum[0] > 0 and removed_squared > 0, "private coherence sensitivity")
    return {"trace": trace, "dephasing_public_probability_change": 0,
            "block_density_difference_frobenius": sqrt(removed_squared),
            "nonzero_entry": maximum[1]}


def exact_audit():
    contract, protocol = admit()
    require(ADDRESSES[0] == (0, 1) and ADDRESSES[1] == (1, 0)
            and ADDRESSES[6] == (1, 2), "address orientation")
    require(loads(1, BACKGROUNDS[0]) == (-1, 1, 0, 0), "bit order and incidence")
    native_dirty = 0
    for background in BACKGROUNDS:
        for law in ("native", "uphill"):
            for word in range(4096):
                for a in range(12):
                    for memory in range(2):
                        for reverse in (False, True):
                            out = reversible(word, a, memory, background, law, reverse)
                            require(reversible(out[0], a, out[1], background, law, not reverse)
                                    == (word, memory), "dirty native inverse")
                        native_dirty += 1
    summaries = []
    source_laws = {}
    for background in BACKGROUNDS:
        for variant in VARIANTS:
            seen = [set() for _ in range(5)]
            public = [Counter() for _ in range(5)]
            noops = 0
            for word in range(4096):
                for address in range(12):
                    history = prefixes(word, address, background, variant)
                    first = event(word, address, background, "uphill" if variant == "uphill" else "native")
                    noops += first[2] == 0
                    if variant == "native":
                        y = first[0]
                        r = (y >> address) & 1
                        z, old, _ = event(y, 6 * r, background)
                        require(history[-1] == Registers((0, 0, 0, 0, z), (0, 0, 0, address),
                                (0, 0, 2 + r), 2 + r, 2 + r, (first[1], old), 6 * r),
                                "derived completed formula")
                    for stage, output in enumerate(history):
                        code = identity(output)
                        require(code not in seen[stage], "isometry disjoint column support")
                        seen[stage].add(code)
                        public[stage][pointer(output)] += 1
                        valid = stage == 4 and variant not in ("no_communication", "no_publish")
                        require(bool(output.pb & 2) == valid, "actual receiver validity")
                        if variant == "no_publish":
                            require(output.pa == output.pb == 0, "invalid agreement excluded")
                        if stage:
                            require(gate(output, stage - 1, background, variant, True) == history[stage - 1],
                                    "prepared prefix inverse")
                    if word == 85:
                        source_laws.setdefault((background, variant), Counter())[history[-1].sources[4]] += 1
            for data in dirty_states():
                for stage in range(4):
                    for reverse in (False, True):
                        output = gate(data, stage, background, variant, reverse)
                        require(gate(output, stage, background, variant, not reverse) == data,
                                "bounded dirty program inverse")
            require(all(len(s) == 4096 * 12 for s in seen), "complete column count")
            summaries.append({"background": background, "variant": variant,
                              "columns": 4096, "entries_per_prefix": [len(s) for s in seen],
                              "first_noops": noops,
                              "public_partitions": [{str(k): v for k, v in c.items()} for c in public]})
    feedback = tv(source_laws[BACKGROUNDS[0], "native"], source_laws[BACKGROUNDS[0], "no_feedback"])
    uphill = tv(source_laws[BACKGROUNDS[0], "native"], source_laws[BACKGROUNDS[0], "uphill"])
    require(feedback >= F(protocol["gates"]["feedback_source_tv_minimum_on_source85_background0"]),
            "completed-source feedback threshold")
    require(uphill >= F(protocol["gates"]["native_uphill_source_tv_minimum_on_source85_background0"]),
            "completed-source uphill threshold")
    return {"target": contract["upstream_revision"], "status": "INDEPENDENT_CONSTRUCTION_CHECKED",
            "native_dirty_states_both_inverse_orders": native_dirty,
            "program_dirty_states_per_configuration": 256, "configurations": summaries,
            "source85_background0_completed_source_laws": {v: dict(source_laws[BACKGROUNDS[0], v]) for v in VARIANTS},
            "completed_source_tv": {"feedback": str(feedback), "uphill": str(uphill)},
            "clock": clock_bound(protocol), "coherence": coherence_witness()}


@contextmanager
def isolated():
    """Deny production imports and every repository read except two inputs.

    An audit hook also covers io.open/os.open. No tracked file is moved.
    The hook stays installed but is inactive after the context exits.
    """
    class Block(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            if fullname == "model" or fullname.startswith("model."):
                raise ImportError("production imports blocked: " + fullname)

    block = Block()
    allowed = {(HERE / n).resolve() for n in ("protocol.json", "observer_independent_contract.json")}
    root = HERE.parents[1]
    active = [True]

    def guard(event_name, arguments):
        if active[0] and event_name == "open" and isinstance(arguments[0], (str, bytes)):
            path = Path(arguments[0]).resolve()
            if path.is_relative_to(root) and path not in allowed:
                raise PermissionError("undeclared repository input: " + str(path))

    sys.addaudithook(guard)
    sys.meta_path.insert(0, block)
    try:
        for name in ("model.native_memory_clock.native", "model.native_observer_archive.program"):
            try:
                __import__(name)
            except ImportError:
                pass
            else:
                raise ValueError("production import isolation failed")
        for name in ("results.json.gz", "extended_results.json.gz"):
            try:
                (HERE / name).read_bytes()
            except PermissionError:
                pass
            else:
                raise ValueError("receipt isolation failed")
        yield
    finally:
        active[0] = False
        sys.meta_path.remove(block)


if __name__ == "__main__":
    with isolated():
        result = exact_audit()
    result["isolation"] = "production imports blocked; repository reads limited to contract and protocol"
    print(json.dumps(result, indent=2))
