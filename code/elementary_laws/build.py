"""Exact finite controls for the elementary corollaries of the record and screen structure.

Every number is a rational computed with ``fractions.Fraction``. The receipt
records the controls named in the paper fragment
``paper/tex_fragments/ELEMENTARY_LAW_COROLLARIES.tex``: the Lueders reread and
the discrete Zeno sequence, the port distinguishability bound, the Moebius
aberration, Doppler and rapidity-addition identities on the celestial sphere,
the single-mode occupation laws under declared exchange sectors, and
coarse-graining monotonicity of the record distribution. It selects no
physical clock, energy, frequency, or statistics law.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = "oph.elementary_laws.receipt.v1"


def s(x: F) -> str:
    return str(F(x))


def mat_mul(a, b):
    n, m, k = len(a), len(b[0]), len(b)
    return [[sum(a[i][t] * b[t][j] for t in range(k)) for j in range(m)] for i in range(n)]


def mat_pow(a, p):
    out = [[F(int(i == j)) for j in range(len(a))] for i in range(len(a))]
    for _ in range(p):
        out = mat_mul(out, a)
    return out


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def lueders_and_zeno():
    # A real symmetric rational density matrix on three records and the
    # projector onto the first two.
    rho = [[F(1, 2), F(1, 8), F(1, 8)],
           [F(1, 8), F(1, 4), F(0)],
           [F(1, 8), F(0), F(1, 4)]]
    assert trace(rho) == 1
    proj = [[F(1), F(0), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(0)]]
    weight = trace(mat_mul(rho, proj))
    prp = mat_mul(mat_mul(proj, rho), proj)
    updated = [[x / weight for x in row] for row in prp]
    reread = trace(mat_mul(updated, proj))
    twice = mat_mul(mat_mul(proj, updated), proj)
    idempotent = twice == updated
    # Discrete Zeno control: a rotation by the rational angle with
    # tan(theta_n / 2) = 1 / n, read through the rank-one projector on e1
    # after every step. cos(theta_n) = (n^2 - 1) / (n^2 + 1).
    zeno = []
    for n in range(2, 13):
        c = F(n * n - 1, n * n + 1)
        sn = F(2 * n, n * n + 1)
        assert c * c + sn * sn == 1
        with_reads = c ** (2 * n)
        rot = [[c, -sn], [sn, c]]
        cos_total = mat_pow(rot, n)[0][0]
        without_reads = cos_total * cos_total
        bound = F(4, n)
        assert 1 - with_reads <= bound
        zeno.append({
            "n": n,
            "cos_step": s(c),
            "survival_with_reads": s(with_reads),
            "survival_without_reads": s(without_reads),
            "one_minus_survival": s(1 - with_reads),
            "bound_four_over_n": s(bound),
        })
    return {
        "density": [[s(x) for x in row] for row in rho],
        "event_weight": s(weight),
        "reread_weight": s(reread),
        "idempotent": idempotent,
        "zeno_sequence": zeno,
    }


def port_bound():
    # In dimension d the standard basis plus the all-ones vector is linearly
    # dependent (zero Gram determinant), so no d + 1 nonzero vectors are
    # pairwise orthogonal: at most d perfectly distinguishable states per port.
    rows = []
    for d in range(2, 6):
        vectors = [[F(int(i == j)) for j in range(d)] for i in range(d)]
        vectors.append([F(1)] * d)
        gram = [[sum(u[t] * v[t] for t in range(d)) for v in vectors] for u in vectors]
        rows.append({"dimension": d, "gram_determinant": s(det(gram)), "max_distinguishable": d})
    return {"per_port": rows, "twelve_ports_dimension_two_bits_per_read": 12}


def det(m):
    m = [row[:] for row in m]
    n = len(m)
    result = F(1)
    for c in range(n):
        pivot = next((r for r in range(c, n) if m[r][c] != 0), None)
        if pivot is None:
            return F(0)
        if pivot != c:
            m[c], m[pivot] = m[pivot], m[c]
            result = -result
        result *= m[c][c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            m[r] = [x - f * y for x, y in zip(m[r], m[c])]
    return result


def cabs2(z):
    return z[0] * z[0] + z[1] * z[1]


def cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def cadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def moebius_scale(a, b, c, d, zeta):
    # Null-vector scale factor (|a + b z|^2 + |c + d z|^2) / (1 + |z|^2) for
    # the spinor (1, z), the convention of the paper fragment.
    num = cabs2(cadd(a, cmul(b, zeta))) + cabs2(cadd(c, cmul(d, zeta)))
    return num / (1 + cabs2(zeta))


def celestial():
    zetas = [(F(0), F(0)), (F(1, 3), F(0)), (F(1), F(0)), (F(2), F(0)), (F(1, 2), F(1, 3))]
    boosts = []
    for k in (F(4), F(9)):
        # diag(sqrt k, 1 / sqrt k) with sqrt k rational for k = 4, 9.
        root = {F(4): F(2), F(9): F(3)}[k]
        aa, dd = (1 / root, F(0)), (root, F(0))  # diag(k^-1/2, k^1/2)
        zero = (F(0), F(0))
        gamma = (k + 1 / k) / 2
        v = (k - 1 / k) / (k + 1 / k)
        rows = []
        for z in zetas:
            z2 = cabs2(z)
            cos_theta = (1 - z2) / (1 + z2)
            zp = (k * z[0], k * z[1])  # (c + d z) / (a + b z) with b = c = 0
            zp2 = cabs2(zp)
            cos_theta_p = (1 - zp2) / (1 + zp2)
            scale = moebius_scale(aa, zero, zero, dd, z)
            assert scale == gamma * (1 - v * cos_theta)
            assert cos_theta_p == (cos_theta - v) / (1 - v * cos_theta)
            rows.append({"zeta": [s(z[0]), s(z[1])], "cos_theta": s(cos_theta),
                         "cos_theta_boosted": s(cos_theta_p), "scale": s(scale)})
        boosts.append({"k": s(k), "velocity": s(v), "gamma": s(gamma), "directions": rows})
    # A rotation in SU(2) with rational entries leaves the scale factor at one.
    alpha, beta = (F(3, 5), F(0)), (F(4, 5), F(0))
    rot_scales = [s(moebius_scale(alpha, (-beta[0], beta[1]), beta, alpha, z)) for z in zetas]
    assert all(F(x) == 1 for x in rot_scales)
    # Rapidity addition: diag boosts multiply, k = k1 k2.
    k1, k2 = F(4), F(9)
    v1 = (k1 - 1 / k1) / (k1 + 1 / k1)
    v2 = (k2 - 1 / k2) / (k2 + 1 / k2)
    k12 = k1 * k2
    v12 = (k12 - 1 / k12) / (k12 + 1 / k12)
    assert v12 == (v1 + v2) / (1 + v1 * v2)
    return {"boosts": boosts, "rotation_scales": rot_scales,
            "velocity_addition": {"v1": s(v1), "v2": s(v2), "composed": s(v12)}}


def occupation():
    rows = []
    for x in (F(1, 3), F(1, 5), F(2, 3)):
        fd = x / (1 + x)
        be = x / (1 - x)
        truncated = {}
        for m in (1, 2, 3, 10):
            num = sum(n * x ** n for n in range(m + 1))
            den = sum(x ** n for n in range(m + 1))
            truncated[str(m)] = s(num / den)
        assert F(truncated["1"]) == fd
        assert all(F(truncated[str(m)]) <= be for m in (1, 2, 3, 10))
        rows.append({"boltzmann_factor": s(x), "fermi_dirac": s(fd), "bose_einstein": s(be),
                     "truncated_mean_occupation": truncated})
    return rows


def coarse_graining():
    p = [F(1, 2), F(1, 4), F(1, 8), F(1, 8)]
    merge = [0, 1, 2, 2]  # records 2 and 3 repair to one public record
    q = [F(0)] * 3
    for i, j in enumerate(merge):
        q[j] += p[i]
    collision_before = sum(x * x for x in p)
    collision_after = sum(x * x for x in q)
    assert collision_after >= collision_before
    # A stochastic repair kernel with the uniform law stationary: the
    # chi-square divergence from the stationary law does not increase.
    kernel = [[F(1, 2), F(1, 4), F(1, 4)], [F(1, 4), F(1, 2), F(1, 4)], [F(1, 4), F(1, 4), F(1, 2)]]
    pi = [F(1, 3)] * 3
    assert [sum(pi[i] * kernel[i][j] for i in range(3)) for j in range(3)] == pi
    start = [F(1), F(0), F(0)]
    after = [sum(start[i] * kernel[i][j] for i in range(3)) for j in range(3)]
    chi = lambda r: sum((r[i] - pi[i]) ** 2 / pi[i] for i in range(3))
    assert chi(after) <= chi(start)
    return {"records": [s(x) for x in p], "repaired": [s(x) for x in q],
            "collision_sum_before": s(collision_before), "collision_sum_after": s(collision_after),
            "kernel_chi_square_before": s(chi(start)), "kernel_chi_square_after": s(chi(after))}


def build():
    receipt = {
        "schema": SCHEMA,
        "scope": ("Exact finite controls for elementary corollaries. No physical clock, energy, "
                  "frequency, statistics law, or source selection is asserted."),
        "lueders_zeno": lueders_and_zeno(),
        "port_bound": port_bound(),
        "celestial": celestial(),
        "occupation": occupation(),
        "coarse_graining": coarse_graining(),
    }
    receipt["fragment_sha256"] = hashlib.sha256(
        (ROOT / "paper/tex_fragments/ELEMENTARY_LAW_COROLLARIES.tex").read_bytes()).hexdigest()
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = build()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.write:
        (HERE / "receipt.json").write_text(text, encoding="utf-8")
    else:
        print(text)


if __name__ == "__main__":
    main()
