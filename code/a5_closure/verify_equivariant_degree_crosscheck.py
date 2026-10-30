"""Independent affine-E8 recurrence and modular polynomial cross-control."""

import json
import hashlib
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
DEFAULT = HERE / "receipts/equivariant_degree_crosscheck.json"
LABELS = ["1", "Q", "Sym2Q", "Sym3Q", "Sym4Q", "Sym5Q", "4_even", "Q_prime", "Sym2Q_prime"]
# The two-dimensional standard node is adjacent to the trivial affine node.
# The other two-dimensional node is seven edges from the affine node.
EDGES = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (5, 8)]


def require(ok, why):
    if not ok:
        raise ValueError(why)


def mod_gcd(a, b, prime=101):
    def clean(p):
        p = [x % prime for x in p]
        while p and p[-1] == 0:
            p.pop()
        return p
    a, b = clean(a), clean(b)
    while b:
        while len(a) >= len(b):
            q, offset = a[-1]*pow(b[-1], -1, prime) % prime, len(a)-len(b)
            for j, x in enumerate(b):
                a[j+offset] = (a[j+offset]-q*x) % prime
            a = clean(a)
        a, b = b, a
    return len(a)-1


def verify(packet=None):
    p = json.loads(DEFAULT.read_text()) if packet is None else packet
    require(p["schema"] == "oph.equivariant_degree_crosscheck.v1", "algebraic schema")
    require(p["implementation_sha256"] == {path: hashlib.sha256((HERE/path).read_bytes()).hexdigest()
            for path in ("equivariant_degree_crosscheck.py", "verify_equivariant_degree_crosscheck.py")}, "implementation custody")
    require(p["labels"] == LABELS and p["maximum_checked_degree"] == 94, "labels/coverage")
    dims = [1, 2, 3, 4, 5, 6, 4, 2, 3]
    require(p["dimensions"] == dims and sum(d*d for d in dims) == 120, "dimensions")
    adjacency = [[0]*9 for _ in range(9)]
    for i, j in EDGES:
        adjacency[i][j] = adjacency[j][i] = 1
    require(p["mckay_adjacency"] == adjacency, "affine E8 independent graph")
    require(all(sum(a*d for a, d in zip(row, dims)) == 2*dims[i]
                for i, row in enumerate(adjacency)), "McKay dimension equation")
    previous, current = [0]*9, [1]+[0]*8
    vectors = [current]
    for degree in range(1, 95):
        following = [sum(x*y for x, y in zip(row, current))-previous[i]
                     for i, row in enumerate(adjacency)]
        require(all(x >= 0 for x in following) and sum(x*d for x, d in zip(following, dims)) == degree+1,
                "symmetric-power dimensions and positivity")
        previous, current = current, following
        vectors.append(current)
    for name, node, exponents in (("Q", 1, [1, 11, 19, 29]), ("Q_prime", 7, [7, 13, 17, 23])):
        counts = [v[node] for v in vectors]
        require(p["multiplicities"][name] == counts, "independent multiplicities " + name)
        # Multiply by the denominator instead of enumerating monomials.
        numerator = [counts[n]-(counts[n-12] if n >= 12 else 0)
                     -(counts[n-20] if n >= 20 else 0)+(counts[n-32] if n >= 32 else 0)
                     for n in range(95)]
        require(numerator == [int(n in exponents) for n in range(95)], "Hilbert numerator")
        require(p["hilbert_numerator_exponents"][name] == exponents, "serialized numerator")
        require(p["minimum_homogeneous_degrees"][name] == next(i for i, c in enumerate(counts) if c), "minimum degree")
    # Validate the producer's whole character table independently against the
    # E8 eigenvector equation, including the marked Q and Q-prime traces.
    traces = [[2, 0], [-2, 0], [0, 0], [1, 0], [-1, 0], [0, 1], [0, -1], [-1, 1], [1, -1]]
    weights = [1, 1, 30, 20, 20, 12, 12, 12, 12]
    require(p["class_Q_traces"] == traces and p["class_sizes"] == weights, "binary lift classes")
    require(p["coefficient_field"] == "Q(phi), phi^2=phi+1", "field convention")
    chars = p["characters"]
    require(len(chars) == 9 and all(len(c) == 9 for c in chars), "character shape")
    require(chars[0] == [[1, 0]]*9 and chars[1] == traces, "trivial/standard labels")
    require(chars[7] == [[a+b, -b] for a, b in traces], "Galois twin label")
    for c, t in enumerate(traces):
        for node in range(9):
            a, b = chars[node][c]
            u, v = t
            require([u*a+v*b, u*b+v*a+v*b] ==
                    [sum(adjacency[node][j]*chars[j][c][part] for j in range(9)) for part in range(2)],
                    "character McKay equation")
    # Full character orthogonality in the phi basis.
    for i in range(9):
        for j in range(9):
            a = sum(w*(x[0]*y[0]+x[1]*y[1]) for w, x, y in zip(weights, chars[i], chars[j]))
            b = sum(w*(x[0]*y[1]+x[1]*y[0]+x[1]*y[1]) for w, x, y in zip(weights, chars[i], chars[j]))
            require((a, b) == (120*int(i == j), 0), "orthogonality recomputation")
    require(p["character_gram"] == [[int(i == j) for j in range(9)] for i in range(9)], "Gram table")
    require(p["denominator_degrees"] == [12, 20], "denominator")
    residues = {"untwisted": [1, 11, 19, 29], "twin": [7, 13, 17, 23]}
    require(p["necessary_residues_mod_30"] == residues, "residues")
    for name, rs in residues.items():
        sign_choices = {1, -1} if name == "untwisted" else {2, -2}
        recomputed = sorted({d % 30 for d in range(-60, 61)
                             if d % 2 and d % 3 and any((d-s) % 5 == 0 for s in sign_choices)})
        require(rs == recomputed, "CRT recomputation")
    require(p["minimum_absolute_degree_from_congruences"] == {"untwisted": 1, "twin": 7}, "minimum absolute degrees")
    h = p["hms"]
    require(h["numerator"] == [-1, 0, -7, 0, 0, -7, 0, 1] and
            h["denominator"] == [1, 0, -7, 0, 0, 7, 0, 1], "HMS equation 3.18 coefficients")
    require(mod_gcd(h["numerator"], h["denominator"]) == 0, "coprime modulo 101, hence over Q")
    require(h["gcd_degree"] == 0 and h["homogeneous_degree"] == 7 and h["common_infinity_zero"] is False, "basepoint/degree")
    require(h["equation"] == "HMS (3.18)" and h["euclidean_remainder_degrees"] == [7, 7, 5, 2, 0, -1], "Euclidean control metadata")
    require(h["source_typed_equivariance"] == "Not identified with OPH port coordinates by this polynomial gcd check", "no source identification")
    require(p["continuous_theorem"] == "FROZEN_PAIR_AND_ALGEBRAIC_CROSSCHECK_ONLY__INDEPENDENT_PROOF_REVIEW_PENDING", "no numerical universal proof")
    require(p["bernhard_original_proof_and_scripts"] == "NOT_LOCATED__NO_REPLAY_CLAIM", "provenance gap")
    require(p["attachment"] == "ATTACHMENT_REMAINS_OPEN" and p["physical_selection"] is False, "physical boundary")
    require(set(p) == {"schema", "implementation_sha256", "coefficient_field", "class_Q_traces", "class_sizes", "labels", "dimensions", "characters", "character_gram", "mckay_adjacency", "maximum_checked_degree", "multiplicities", "hilbert_numerator_exponents", "denominator_degrees", "minimum_homogeneous_degrees", "necessary_residues_mod_30", "minimum_absolute_degree_from_congruences", "hms", "continuous_theorem", "bernhard_original_proof_and_scripts", "attachment", "physical_selection"}, "closed schema")
    return dict(checked_through=94, minimum_homogeneous_degrees=p["minimum_homogeneous_degrees"],
                residues=residues, continuous_theorem=p["continuous_theorem"])


if __name__ == "__main__":
    require(len(sys.argv) == 1, "no CLI arguments")
    print(json.dumps(verify(), sort_keys=True))
