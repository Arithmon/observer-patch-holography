"""Construct response histories using the existing exact carrier arithmetic.

This is an explicitly defined mathematical source, not a capture of the
registered amplitude driver. The checker reconstructs D by a different formula.
"""
from fractions import Fraction
from types import SimpleNamespace
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "a5_closure"))
import port_current_inner_certificate as algebra
import response_grammar_completeness_certificate as incidence


def encode_real(x):
    return [str(x.a), str(x.b)]


def encode_complex(x):
    return encode_real(x.re) + encode_real(x.im)


def encode_blocks(blocks):
    return [[[encode_complex(x) for x in row] for row in block] for block in blocks]


def identity():
    return [[algebra.C5(algebra.ONE if i == j else algebra.ZERO)
             for j in range(3)] for i in range(3)]


def inverse(matrix):
    # The Cayley factors here are real. Solve all columns over Q(sqrt(5)).
    n = len(matrix)
    rows = [[x.re for x in row] + [algebra.ONE if i == j else algebra.ZERO
                                 for j in range(n)]
            for i, row in enumerate(matrix)]
    reduced, pivots = algebra.rref(rows)
    if pivots[:n] != list(range(n)):
        raise ValueError("singular Cayley denominator")
    return [[algebra.C5(reduced[i][n+j]) for j in range(n)] for i in range(n)]


def multiply(left, right):
    return tuple(algebra.cmul(a, b) for a, b in zip(left, right, strict=True))


def geometric_response(field, frame):
    """A primitive law on six internal coordinates; no Lie-type input."""
    even, odd = frame.even_odd(field)
    mean = sum(even, algebra.ZERO) / algebra.F5(6)
    sym = algebra.rzeros(3, 3)
    for weight, vector in zip(even, frame.axis_vectors):
        for i in range(3):
            for j in range(3):
                sym[i][j] = sym[i][j] + (weight-mean)*vector[i]*vector[j]
    for i in range(3):
        sym[i][i] = sym[i][i] + mean
    return (algebra.to_cmat(algebra.hat(frame.frame_map(odd)), sym),
            algebra.to_cmat(algebra.hat(frame.galois_frame_map(odd))))


def ordered_jet(first, second):
    """Execute four response factors in the ring Q[s,t]/(s^2,t^2).

    A Cayley response with control s*D/2 is I+sD+O(s^2).
    The mixed coefficient is extracted from the whole ordered word.
    """
    unit = (identity(), identity())
    zero = (algebra.czeros(3), algebra.czeros(3))
    word = [((1, 0), first), ((0, 1), second),
            ((1, 0), tuple([[ -x for x in row] for row in b] for b in first)),
            ((0, 1), tuple([[ -x for x in row] for row in b] for b in second))]
    result = {(0, 0): unit}
    for exponent, coefficient in word:
        following = dict(result)
        for (a, b), value in result.items():
            power = (a+exponent[0], b+exponent[1])
            if max(power) <= 1:
                term = multiply(value, coefficient)
                following[power] = tuple(algebra.cadd(x, y)
                                         for x, y in zip(following.get(power, zero), term))
        result = following
    return result[(1, 1)]


def source_packet():
    verts, adjacency, _, antipode = incidence.port_model()
    group = incidence.rotation_permutations(verts, adjacency)
    carrier = SimpleNamespace(antipode=antipode)
    frame = algebra.FrameRealization(carrier, tuple(range(12)), verts)
    generators = [geometric_response(f, frame) for f in algebra.BASIS_FIELDS]
    flat = [algebra.flatten(g) for g in generators]
    tomography = []
    for g in generators:
        full = [[g[i//3][i%3][j%3] if i//3 == j//3 else algebra.C5()
                 for j in range(6)] for i in range(6)]
        observed = []
        for i in range(6):
            for j in range(6):
                if i != j:
                    probe = algebra.czeros(6)
                    probe[j][j] = algebra.C5(algebra.ONE)
                    observed.append(algebra.commutator(full, probe)[i][j])
        for i in range(5):
            probe = algebra.czeros(6)
            probe[i][5] = algebra.C5(algebra.ONE)
            observed.append(algebra.commutator(full, probe)[i][5])
        tomography.append([encode_complex(x) for x in observed])
    jets = []
    for p in range(12):
        for q in range(p+1, 12):
            coefficient = ordered_jet(generators[p], generators[q])
            solution = algebra.solve_in_span(flat, algebra.flatten(coefficient))
            jets.append({"ports": [p, q], "mixed_coefficient": encode_blocks(coefficient),
                         "in_port_basis": [encode_real(x) for x in solution]})
    unit = tuple(range(12))
    thirds = [g for g in group if g != unit and
              tuple(g[g[g[p]]] for p in range(12)) == unit]
    words = {unit: []}
    for g in thirds:
        words[g] = [g]
    for g in thirds:
        for h in thirds:
            words.setdefault(tuple(h[g[p]] for p in range(12)), [g, h])
    if set(words) != set(group):
        raise ValueError("order-three factors do not cover the proper group")
    factors = []
    executed = {}
    for g in thirds:
        rotation = frame.rotation_of(g)
        target = (algebra.to_cmat(rotation),
                  algebra.to_cmat([[x.conj() for x in row] for row in rotation]))
        k = tuple(algebra.cmul(algebra.csub(r, identity()),
                              inverse(algebra.cadd(r, identity()))) for r in target)
        coefficients = algebra.solve_in_span(flat, algebra.flatten(k))
        # Execute the source response through its Cayley law, not the target.
        actual = tuple(algebra.cmul(algebra.cadd(identity(), b),
                                   inverse(algebra.csub(identity(), b))) for b in k)
        executed[g] = actual
        factors.append({"port_action": list(g),
                        "control": [encode_real(x) for x in coefficients],
                        "observed_response": encode_blocks(actual)})
    by_action = {tuple(f["port_action"]): f for f in factors}
    paths = []
    for g in group:
        actual = (identity(), identity())
        for move in words[g]:
            actual = multiply(executed[move], actual)
        paths.append({"port_action": list(g),
                      "factor_indices": [factors.index(by_action[h]) for h in words[g]],
                      "observed_response": encode_blocks(actual)})
    return {"schema": "oph.source-selection-response.v1",
            "generators": [encode_blocks(g) for g in generators],
            "matrix_unit_response_probes": tomography,
            "ordered_mixed_responses": jets,
            "cayley_factors": factors, "closed_paths": paths}


if __name__ == "__main__":
    import json
    destination = Path(__file__).with_name("response.json")
    destination.write_text(json.dumps(source_packet(), separators=(",", ":")) + "\n",
                           encoding="utf-8", newline="\n")
    print(destination.name, destination.stat().st_size)
