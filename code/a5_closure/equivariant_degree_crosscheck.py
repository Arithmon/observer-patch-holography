"""Exact 2I characters and degree arithmetic; not a continuous-map proof gate."""

from fractions import Fraction
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "receipts/equivariant_degree_crosscheck.json"
ZERO, ONE = (0, 0), (1, 0)


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[0], -x[1]


def mul(x, y):
    a, b = x
    c, d = y
    return a*c+b*d, a*d+b*c+b*d


def sigma(x):
    return x[0]+x[1], -x[1]


def need(ok, why):
    if not ok:
        raise ValueError(why)


def sum_field(xs):
    out = ZERO
    for x in xs:
        out = add(out, x)
    return out


def inner(x, y, sizes):
    a, b = sum_field(mul((w, 0), mul(u, v)) for w, u, v in zip(sizes, x, y))
    need(b == 0 and a % 120 == 0, "integral character inner product")
    return a // 120


def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def remainder(a, b):
    a, b = trim(a), trim(b)
    while len(a) >= len(b):
        shift, scale = len(a)-len(b), a[-1]/b[-1]
        for j, v in enumerate(b):
            a[j+shift] -= scale*v
        a = trim(a)
    return a


def hms_control():
    # Houghton--Manton--Sutcliffe (3.18), coefficient order low to high.
    p = [-1, 0, -7, 0, 0, -7, 0, 1]
    q = [1, 0, -7, 0, 0, 7, 0, 1]
    a, b = list(map(Fraction, p)), list(map(Fraction, q))
    degrees = [7, 7]
    while b:
        a, b = b, remainder(a, b)
        degrees.append(len(b)-1)
    need(len(a) == 1 and a[0] != 0, "HMS common finite zero")
    need((p[-1], q[-1]) != (0, 0), "HMS common infinity zero")
    return dict(equation="HMS (3.18)", numerator=p, denominator=q,
                euclidean_remainder_degrees=degrees, gcd_degree=0,
                homogeneous_degree=7, common_infinity_zero=False,
                source_typed_equivariance="Not identified with OPH port coordinates by this polynomial gcd check")


def build():
    # Q trace classes of the binary lift. Opposite signs are distinct lifts.
    traces = [(2, 0), (-2, 0), (0, 0), (1, 0), (-1, 0),
              (0, 1), (0, -1), (-1, 1), (1, -1)]
    sizes = [1, 1, 30, 20, 20, 12, 12, 12, 12]
    sym = [[ONE]*9, traces]
    for n in range(2, 95):
        sym.append([add(mul(t, a), neg(b)) for t, a, b in zip(traces, sym[-1], sym[-2])])
    qprime = list(map(sigma, traces))
    threeprime = [add(mul(t, t), neg(ONE)) for t in qprime]
    four = [add(a, neg(b)) for a, b in zip(sym[6], threeprime)]
    chars = [*sym[:6], four, qprime, threeprime]
    labels = ["1", "Q", "Sym2Q", "Sym3Q", "Sym4Q", "Sym5Q", "4_even", "Q_prime", "Sym2Q_prime"]
    gram = [[inner(a, b, sizes) for b in chars] for a in chars]
    need(gram == [[int(i == j) for j in range(9)] for i in range(9)], "complete irreducible orthogonality")
    dimensions = [c[0][0] for c in chars]
    need(sum(d*d for d in dimensions) == 120, "complete irreducible dimension squares")
    adjacency = [[inner([mul(t, x) for t, x in zip(traces, a)], b, sizes)
                  for b in chars] for a in chars]
    multiplicities = [[inner(c, q, sizes) for c in sym] for q in (traces, qprime)]
    for got, exponents in zip(multiplicities, ([1, 11, 19, 29], [7, 13, 17, 23])):
        expected = [sum(n == e+12*a+20*b for e in exponents for a in range(8) for b in range(5)) for n in range(95)]
        need(got == expected, "Hilbert series through 94")
    residues = {name: [d for d in range(30) if d % 2 == 1 and d % 3 and d % 5 in allowed]
                for name, allowed in (("untwisted", (1, 4)), ("twin", (2, 3)))}
    minima = {name: min(abs(d) for d in range(-30, 31) if d % 30 in rs) for name, rs in residues.items()}
    implementation = {p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in (
        "equivariant_degree_crosscheck.py", "verify_equivariant_degree_crosscheck.py")}
    return dict(schema="oph.equivariant_degree_crosscheck.v1", implementation_sha256=implementation,
                coefficient_field="Q(phi), phi^2=phi+1",
                class_Q_traces=traces, class_sizes=sizes, labels=labels, dimensions=dimensions,
                characters=chars, character_gram=gram, mckay_adjacency=adjacency,
                maximum_checked_degree=94, multiplicities=dict(Q=multiplicities[0], Q_prime=multiplicities[1]),
                hilbert_numerator_exponents=dict(Q=[1, 11, 19, 29], Q_prime=[7, 13, 17, 23]),
                denominator_degrees=[12, 20], minimum_homogeneous_degrees=dict(Q=1, Q_prime=7),
                necessary_residues_mod_30=residues, minimum_absolute_degree_from_congruences=minima,
                hms=hms_control(),
                continuous_theorem="FROZEN_PAIR_AND_ALGEBRAIC_CROSSCHECK_ONLY__INDEPENDENT_PROOF_REVIEW_PENDING",
                bernhard_original_proof_and_scripts="NOT_LOCATED__NO_REPLAY_CLAIM",
                attachment="ATTACHMENT_REMAINS_OPEN", physical_selection=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.write:
        RECEIPT.write_text(payload)
    elif args.check:
        need(RECEIPT.read_text() == payload, "algebraic receipt replay")
        print("EQUIVARIANT_DEGREE_CROSSCHECK_REPLAY_VALID")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
