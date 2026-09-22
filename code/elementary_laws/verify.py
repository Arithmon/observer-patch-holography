"""Independent replay of the elementary-law receipt; imports no producer code."""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def fail(msg):
    raise ValueError(msg)


def check_zeno(block):
    for row in block["zeno_sequence"]:
        n = row["n"]
        c = Q(n * n - 1, n * n + 1)
        if Q(row["cos_step"]) != c:
            fail("zeno cos")
        if Q(row["survival_with_reads"]) != c ** (2 * n):
            fail("zeno survival")
        if Q(row["one_minus_survival"]) > Q(4, n):
            fail("zeno bound")
        # Chebyshev recursion for cos(n theta) from cos theta.
        t0, t1 = Q(1), c
        for _ in range(n - 1):
            t0, t1 = t1, 2 * c * t1 - t0
        if Q(row["survival_without_reads"]) != t1 * t1:
            fail("zeno no-read survival")
    if Q(block["reread_weight"]) != 1 or not block["idempotent"]:
        fail("lueders")


def check_celestial(block):
    for boost in block["boosts"]:
        k = Q(boost["k"])
        v, g = Q(boost["velocity"]), Q(boost["gamma"])
        if v != (k * k - 1) / (k * k + 1) or g != (k * k + 1) / (2 * k):
            fail("boost parameters")
        for row in boost["directions"]:
            ct = Q(row["cos_theta"])
            if Q(row["scale"]) != g * (1 - v * ct):
                fail("doppler")
            if Q(row["cos_theta_boosted"]) != (ct - v) / (1 - v * ct):
                fail("aberration")
    if any(Q(x) != 1 for x in block["rotation_scales"]):
        fail("rotation scale")
    va = block["velocity_addition"]
    v1, v2, v12 = (Q(va[key]) for key in ("v1", "v2", "composed"))
    if v12 != (v1 + v2) / (1 + v1 * v2):
        fail("velocity addition")


def check_occupation(rows):
    for row in rows:
        x = Q(row["boltzmann_factor"])
        if Q(row["fermi_dirac"]) != x / (1 + x) or Q(row["bose_einstein"]) != x / (1 - x):
            fail("occupation")
        if Q(row["truncated_mean_occupation"]["1"]) != Q(row["fermi_dirac"]):
            fail("truncation at one is fermi-dirac")


def check_coarse(block):
    p = [Q(x) for x in block["records"]]
    q = [Q(x) for x in block["repaired"]]
    if sum(p) != 1 or sum(q) != 1 or len(q) >= len(p):
        fail("coarse graining shape")
    if Q(block["collision_sum_after"]) < Q(block["collision_sum_before"]):
        fail("collision monotonicity")
    if Q(block["kernel_chi_square_after"]) > Q(block["kernel_chi_square_before"]):
        fail("kernel monotonicity")


def check_ports(block):
    for row in block["per_port"]:
        if Q(row["gram_determinant"]) != 0 or row["max_distinguishable"] != row["dimension"]:
            fail("port bound")


def verify(receipt):
    if receipt["schema"] != "oph.elementary_laws.receipt.v1":
        fail("schema")
    digest = hashlib.sha256(
        (ROOT / "paper/tex_fragments/ELEMENTARY_LAW_COROLLARIES.tex").read_bytes()).hexdigest()
    if receipt["fragment_sha256"] != digest:
        fail("fragment pin")
    check_zeno(receipt["lueders_zeno"])
    check_ports(receipt["port_bound"])
    check_celestial(receipt["celestial"])
    check_occupation(receipt["occupation"])
    check_coarse(receipt["coarse_graining"])
    return {"ok": True, "zeno_rows": len(receipt["lueders_zeno"]["zeno_sequence"])}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", default=str(HERE / "receipt.json"))
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt))))


if __name__ == "__main__":
    main()
