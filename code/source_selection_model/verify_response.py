"""Independent formula, matrix identities and complete response-path checks.

Only the existing exact number-field arithmetic is shared with the producer.
No producer model, frame constructor, group generator or solver is called.
"""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "a5_closure"))
from port_current_inner_certificate import F5, C5

Z, O = F5(), F5(1)
CZ, CO = C5(), C5(O)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def keys(value, expected):
    require(type(value) is dict and set(value) == set(expected), "object keys")


def rat(value):
    require(type(value) is str and len(value) < 128, "rational type/size")
    x = Q(value)
    require(str(x) == value, "noncanonical rational")
    return x


def real(value):
    require(type(value) is list and len(value) == 2, "field element")
    return F5(rat(value[0]), rat(value[1]))


def block(value):
    require(type(value) is list and len(value) == 3, "matrix rows")
    out = []
    for row in value:
        require(type(row) is list and len(row) == 3, "matrix columns")
        line = []
        for x in row:
            require(type(x) is list and len(x) == 4, "complex field element")
            line.append(C5(real(x[:2]), real(x[2:])))
        out.append(line)
    return out


def blocks(value):
    require(type(value) is list and len(value) == 2, "response block count")
    return tuple(block(x) for x in value)


def eye():
    return [[CO if i == j else CZ for j in range(3)] for i in range(3)]


def add(a, b):
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]


def neg(a):
    return [[-x for x in row] for row in a]


def mm(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(3)), CZ)
             for j in range(3)] for i in range(3)]


def star(a):
    return [[a[j][i].conj() for j in range(3)] for i in range(3)]


def bmul(a, b):
    return tuple(mm(x, y) for x, y in zip(a, b, strict=True))


def combine(gs, coefficients):
    require(len(coefficients) == 12, "control dimension")
    return tuple([[sum((C5(coefficients[p])*gs[p][k][i][j]
                        for p in range(12)), CZ) for j in range(3)]
                  for i in range(3)] for k in range(2))


def rank(rows):
    rows = [list(r) for r in rows]
    lead = 0
    for col in range(len(rows[0])):
        candidate = next((i for i in range(lead, len(rows)) if rows[i][col] != Z), None)
        if candidate is None:
            continue
        rows[lead], rows[candidate] = rows[candidate], rows[lead]
        scale = rows[lead][col]
        rows[lead] = [x/scale for x in rows[lead]]
        for i in range(len(rows)):
            if i != lead:
                scale = rows[i][col]
                rows[i] = [x-scale*y for x, y in zip(rows[i], rows[lead])]
        lead += 1
        if lead == len(rows):
            break
    return lead


def vertices():
    phi = F5(Q(1, 2), Q(1, 2))
    base = [(Z, F5(s), F5(t)*phi) for s in (1, -1) for t in (1, -1)]
    return base + [(y, z, x) for x, y, z in base] + [(z, x, y) for x, y, z in base]


def dot(v, w):
    return sum((a*b for a, b in zip(v, w)), Z)


def determinant(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def reconstructed_generators():
    vs = vertices()
    norm = dot(vs[0], vs[0])
    out = []
    for v in vs:
        both = []
        for conjugate in (False, True):
            w = [x.conj()/F5(2) if conjugate else x/F5(2) for x in v]
            skew = [[Z, -w[2], w[1]], [w[2], Z, -w[0]], [-w[1], w[0], Z]]
            both.append([[C5(skew[i][j], Z if conjugate else
                             v[i]*v[j]/F5(2) +
                             (F5(Q(1, 12))-norm/F5(6) if i == j else Z))
                          for j in range(3)] for i in range(3)])
        out.append(tuple(both))
    return out


def permutation(value):
    require(type(value) is list and len(value) == 12 and
            all(type(x) is int for x in value) and sorted(value) == list(range(12)),
            "port permutation")
    return tuple(value)


def check_port_action(u, g):
    vs = vertices()
    for b in u:
        require(all(x.im == Z for row in b for x in row), "proper rotation must be real")
        require(mm(star(b), b) == eye(), "rotation orthogonality")
        require(determinant([[x.re for x in row] for row in b]) == O, "orientation")
    require([[x.re.conj() for x in row] for row in u[0]] ==
            [[x.re for x in row] for row in u[1]], "Galois response pairing")
    for p, v in enumerate(vs):
        actual = tuple(sum((u[0][i][j].re*v[j] for j in range(3)), Z) for i in range(3))
        require(actual == vs[g[p]], "executed response port action")


def verify(packet):
    keys(packet, ("schema", "generators", "matrix_unit_response_probes",
                  "ordered_mixed_responses", "cayley_factors", "closed_paths"))
    require(packet["schema"] == "oph.source-selection-response.v1", "schema")
    gs = [blocks(x) for x in packet["generators"]]
    require(gs == reconstructed_generators(), "source generator formula")
    require(rank([[component for b in g for row in b for x in row
                   for component in (x.re, x.im)] for g in gs]) == 12, "faithful response")
    for g in gs:
        require(all(star(b) == neg(b) for b in g), "skew adjointness")
    tomography = packet["matrix_unit_response_probes"]
    require(type(tomography) is list and len(tomography) == 12, "tomography port coverage")
    response_rows, restricted_rows = [], []
    for g, row in zip(gs, tomography):
        require(type(row) is list and len(row) == 35, "tomography observable coverage")
        observed = []
        for value in row:
            require(type(value) is list and len(value) == 4, "tomography complex value")
            observed.append(C5(real(value[:2]), real(value[2:])))
        full = [[g[i//3][i%3][j%3] if i//3 == j//3 else CZ for j in range(6)] for i in range(6)]
        expected = [full[i][j] for i in range(6) for j in range(6) if i != j]
        expected += [full[i][i]-full[5][5] for i in range(5)]
        require(observed == expected, "response tomography")
        response_rows.append([part for x in observed for part in (x.re, x.im)])
        restricted = [b[i][j] for b in g for i in range(3) for j in range(3) if i != j]
        restricted += [b[i][i]-b[2][2] for b in g for i in range(2)]
        restricted_rows.append([part for x in restricted for part in (x.re, x.im)])
    require(rank(response_rows) == 12, "operational response faithfulness")
    require(rank(restricted_rows) == 11, "cross-block coherence deletion control")
    expected_pairs = list(combinations(range(12), 2))
    rows = packet["ordered_mixed_responses"]
    require(type(rows) is list and len(rows) == 66, "complete ordered probes")
    for row, (p, q) in zip(rows, expected_pairs):
        keys(row, ("ports", "mixed_coefficient", "in_port_basis"))
        require(row["ports"] == [p, q] and all(type(i) is int for i in row["ports"]), "probe order")
        expected = tuple(add(mm(a, b), neg(mm(b, a))) for a, b in zip(gs[p], gs[q]))
        require(blocks(row["mixed_coefficient"]) == expected, "ordered response jet")
        require(combine(gs, [real(x) for x in row["in_port_basis"]]) == expected, "bracket closure")
    # All independent pairs of real skew blocks occur in the odd image.
    vs = vertices()
    anti = [vs.index(tuple(-x for x in v)) for v in vs]
    odd = []
    for p in range(12):
        if p < anti[p]:
            odd.append([gs[p][b][i][j].re-gs[anti[p]][b][i][j].re
                        for b in range(2) for i, j in ((0, 1), (0, 2), (1, 2))])
    require(rank(odd) == 6, "same-response rotation logarithms")
    factors = packet["cayley_factors"]
    require(type(factors) is list and len(factors) == 20, "complete order-three factors")
    executed, factor_actions = [], []
    for f in factors:
        keys(f, ("port_action", "control", "observed_response"))
        g = permutation(f["port_action"])
        require(g != tuple(range(12)) and tuple(g[g[g[p]]] for p in range(12)) == tuple(range(12)),
                "factor order")
        k = combine(gs, [real(x) for x in f["control"]])
        u = blocks(f["observed_response"])
        for kb, ub in zip(k, u):
            require(star(kb) == neg(kb), "Cayley generator skew")
            # I-K is invertible for skew K; multiplication certifies its unique Cayley output.
            require(mm(ub, add(eye(), neg(kb))) == add(eye(), kb), "Cayley execution")
        check_port_action(u, g)
        executed.append(u)
        factor_actions.append(g)
    require(len(set(factor_actions)) == 20, "duplicate primitive response factor")
    paths = packet["closed_paths"]
    require(type(paths) is list and len(paths) == 60, "holonomy coverage")
    table = {}
    directed_images = set()
    phi = F5(Q(1, 2), Q(1, 2))
    edge = next((p, q) for p, q in combinations(range(12), 2) if dot(vs[p], vs[q]) == phi)
    for path in paths:
        keys(path, ("port_action", "factor_indices", "observed_response"))
        g = permutation(path["port_action"])
        word = path["factor_indices"]
        require(type(word) is list and len(word) <= 2 and
                all(type(i) is int and 0 <= i < 20 for i in word), "closed path word")
        u, action = (eye(), eye()), tuple(range(12))
        for i in word:
            u = bmul(executed[i], u)
            action = tuple(factor_actions[i][action[p]] for p in range(12))
        require(action == g and u == blocks(path["observed_response"]), "closed path replay")
        check_port_action(u, g)
        require(g not in table, "duplicate closed path action")
        table[g] = u
        directed_images.add((g[edge[0]], g[edge[1]]))
        for p in range(12):
            conjugated = bmul(bmul(u, gs[p]), tuple(star(b) for b in u))
            require(conjugated == gs[g[p]], "response naturality")
    all_edges = {(p, q) for p in range(12) for q in range(12) if dot(vs[p], vs[q]) == phi}
    require(directed_images == all_edges and len(all_edges) == 60, "all proper incidence actions")
    for g, ug in table.items():
        for h, uh in table.items():
            composite = tuple(g[h[p]] for p in range(12))
            require(composite in table and bmul(ug, uh) == table[composite], "path functor composition")
    return {"response_rank": 12, "observable_response_rank": 12,
            "block_diagonal_observable_control_rank": 11,
            "ordered_brackets": 66, "cayley_factors": 20,
            "closed_path_actions": 60, "response_naturality_pairs": 720,
            "path_compositions": 3600}


def strict_load(path):
    def unique(pairs):
        out = {}
        for k, v in pairs:
            require(k not in out, "duplicate JSON key")
            out[k] = v
        return out
    require(path.stat().st_size < 200_000, "packet size")
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique,
                      parse_constant=lambda _: (_ for _ in ()).throw(ValueError("nonfinite JSON")))


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) == 2 else Path(__file__).with_name("response.json")
    print(json.dumps(verify(strict_load(path)), sort_keys=True))
