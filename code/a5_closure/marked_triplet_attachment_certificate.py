"""Exact fixed-action attachment obstruction from live native incidence.

No external golden target is consumed. Output is a finite representation
certificate; it asserts no nonselection theorem for arbitrary source laws.
"""
import ast
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

from galois_port_frame_certificate import F, Z, O, S, PHI, encode, mm, transpose, dot, require, coordinates

ROOT = Path(__file__).resolve().parents[2]
INCIDENCE = 'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'
RESPONSE = 'code/source_selection_model/response.json'
INPUTS = [INCIDENCE, RESPONSE, 'model/common_source/slow_geometry.py',
          'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
          'code/source_selection_model/verify_response.py',
          'code/a5_closure/manifests/galois_port_frame_reference.json']
BASELINE = '879ed37ad2c1933ecb6dd8f54271362ee3f4d41c'
INTEGRATION_BASELINE = '7b0997b4a8a425679d9e2946f82c86cd92289059'


def rank(rows):
    rows = [list(row) for row in rows]
    pivot = 0
    for col in range(len(rows[0])):
        k = next((k for k in range(pivot, len(rows)) if rows[k][col] != Z), None)
        if k is None:
            continue
        rows[pivot], rows[k] = rows[k], rows[pivot]
        scale = rows[pivot][col]
        rows[pivot] = [x/scale for x in rows[pivot]]
        for i in range(pivot+1, len(rows)):
            scale = rows[i][col]
            if scale != Z:
                rows[i] = [x-scale*y for x, y in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def native_incidence():
    tree = ast.parse((ROOT/INCIDENCE).read_text())
    faces = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                 and any(isinstance(t, ast.Name) and t.id == 'ORIENTED_BASE_FACES' for t in n.targets))
    edges = {tuple(sorted((a, b))) for face in faces for a in face for b in face if a != b}
    require(len(faces) == 20 and len(edges) == 30, 'native incidence census')
    ds = []
    for i in range(12):
        row = [0 if i == j else 99 for j in range(12)]
        for _ in range(3):
            row = [min(row[j], min(row[k]+1 for k in range(12) if tuple(sorted((j,k))) in edges))
                   for j in range(12)]
        ds.append(row)
    require(all(sorted(row) == [0]+[1]*5+[2]*5+[3] for row in ds), 'distance census')
    return ds


def chart_maps(ds, swap=False):
    """All graph automorphisms or all edge/distance-two exchanges."""
    target = [[{0:0, 1:2, 2:1, 3:3}[d] if swap else d for d in row] for row in ds]
    out = []
    def extend(p):
        if len(p) == 12:
            out.append(tuple(p[i] for i in range(12)))
            return
        i = max((i for i in range(12) if i not in p),
                key=lambda i: (sum(ds[i][j] == 1 for j in p), -i))
        for a in range(12):
            if a not in p.values() and all(ds[i][j] == target[a][b] for j,b in p.items()):
                extend({**p, i:a})
    extend({})
    return sorted(out)


def raw_frame():
    """Bind the coordinate witness to the live source's restricted AST."""
    carrier = ast.parse((ROOT/INPUTS[3]).read_text())
    golden = next(n.value for n in carrier.body if isinstance(n,ast.AnnAssign)
                  and isinstance(n.target,ast.Name) and n.target.id == 'GOLDEN')
    require(isinstance(golden,ast.Tuple) and len(golden.elts) == 2, 'native golden pair')
    def rational(node):
        require(isinstance(node,ast.Call) and isinstance(node.func,ast.Name)
                and node.func.id == 'Fraction' and not node.keywords
                and all(isinstance(x,ast.Constant) and type(x.value) is int for x in node.args),
                'native golden rational')
        return Q(*(x.value for x in node.args))
    p = F(*(rational(x) for x in golden.elts))
    require(p == PHI, 'native positive golden value')
    tree = ast.parse((ROOT/INPUTS[2]).read_text())
    fun = next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name == 'raw_frame')
    assignments = [n for n in fun.body if isinstance(n,ast.Assign)]
    expected = ast.parse('c = source()\nzero, one, minus, phi = c.q5(0), c.q5(1), c.q5(-1), c.GOLDEN\nnegative_phi = c.q5_sub(zero, phi)').body
    require([ast.dump(n) for n in assignments] == [ast.dump(n) for n in expected], 'native frame local bindings')
    returned = next(n.value for n in fun.body if isinstance(n,ast.Return))
    names = {'zero':Z,'one':O,'minus':-O,'phi':p,'negative_phi':-p}
    require(isinstance(returned,ast.Tuple) and len(returned.elts) == 12
            and all(isinstance(row,ast.Tuple) and len(row.elts) == 3
                    and all(isinstance(x,ast.Name) and x.id in names for x in row.elts)
                    for row in returned.elts), 'native raw frame shape')
    return [[names[x.id] for x in row.elts] for row in returned.elts]


def response_vertices():
    base = [[Z,F(s),F(t)*PHI] for s in (1,-1) for t in (1,-1)]
    return base + [[y,z,x] for x,y,z in base] + [[z,x,y] for x,y,z in base]


def source_marked_action(frame, vertices, relabel):
    """Transport the #1088 frozen order-five generator into the live chart."""
    packet = json.loads((ROOT/INPUTS[5]).read_text())
    old_to_response = [vertices.index(v) for v in coordinates()]
    response_to_old = [old_to_response.index(i) for i in range(12)]
    old = packet['generator']
    require(type(old) is list and len(old) == 12 and
            all(type(i) is int for i in old) and sorted(old) == list(range(12)),
            'frozen #1088 generator permutation')
    response_action = [old_to_response[old[response_to_old[i]]] for i in range(12)]
    response_to_native = [relabel.index(i) for i in range(12)]
    native_action = tuple(response_to_native[response_action[relabel[i]]]
                          for i in range(12))
    power = tuple(range(12))
    for n in range(1, 6):
        power = tuple(native_action[power[i]] for i in range(12))
        require((power == tuple(range(12))) == (n == 5), '#1088 transported order five')
    return native_action


def hom_rank(actions, left, right):
    equations = []
    for a, b in zip(left, right, strict=True):
        for i in range(3):
            for j in range(3):
                equations.append([(b[i][k] if l == j else Z) - (a[l][j] if k == i else Z)
                                  for k in range(3) for l in range(3)])
    return rank(equations)


def build(response=None):
    ds = native_incidence()
    # Both projectors are derived from the graph distances. The strict slow
    # ordering follows 5-sqrt(5) < 5+sqrt(5) in the declared real field.
    projectors = []
    for sign in (1,-1):
        table = [O, F(0,Q(sign,5)), F(0,Q(-sign,5)), -O]
        p = [[table[d]/4 for d in row] for row in ds]
        require(mm(p,p) == p and sum((p[i][i] for i in range(12)),Z) == F(3), 'rank-three projector')
        lap = [[F(5 if i == j else -int(ds[i][j] == 1)) for j in range(12)] for i in range(12)]
        require(mm(lap,p) == [[(F(5)-sign*S)*x for x in row] for row in p], 'native eigenband')
        projectors.append(p)
    frame, vertices = raw_frame(), response_vertices()
    relabel = [vertices.index(v) for v in frame]
    norm = dot(frame[0],frame[0])
    require([[dot(x,y)/(4*norm) for y in frame] for x in frame] == projectors[0], 'native slow-frame witness')
    packet = response if response is not None else json.loads((ROOT/RESPONSE).read_text())
    rows = packet['closed_paths']
    require(len(rows) == 60, 'all live response paths')
    actions, upper, lower, chars = [], [], [], []
    for row in rows:
        action = row['port_action']
        require(sorted(action) == list(range(12)), 'response permutation')
        p = tuple(relabel.index(action[relabel[i]]) for i in range(12))
        require(all(ds[p[i]][p[j]] == ds[i][j] for i in range(12) for j in range(12)), 'response incidence binding')
        blocks = []
        for block in row['observed_response']:
            require(all(x[2:] == ['0','0'] for r in block for x in r), 'real response rotations')
            blocks.append([[F(Q(x[0]), Q(x[1])) for x in r] for r in block])
        a,b = blocks
        require(b == [[x.sigma() for x in r] for r in a], 'Galois paired internal actions')
        require(mm(a,transpose(a)) == [[O if i == j else Z for j in range(3)] for i in range(3)], 'orthogonal action')
        require(all([dot(r,v) for r in a] == vertices[action[i]] for i,v in enumerate(vertices)), 'marked upper-frame intertwiner')
        cs = [sum((projector[i][p[i]] for i in range(12)),Z) for projector in projectors]
        require(cs == [sum((m[i][i] for i in range(3)),Z) for m in blocks], 'spatial/internal character attachment')
        actions.append(p); upper.append(a); lower.append(b); chars.append(cs)
    require(len(set(actions)) == 60, 'distinct proper actions')
    for i,p in enumerate(actions):
        for j,q in enumerate(actions):
            pq = tuple(p[q[k]] for k in range(12))
            require(pq in actions and mm(upper[i],upper[j]) == upper[actions.index(pq)], 'representation composition')
    require(hom_rank(actions,upper,lower) == 9, 'zero fixed-action Hom')
    require(hom_rank(actions,upper,upper) == 8, 'same-triplet Hom control')
    charts, twists = chart_maps(ds), chart_maps(ds,True)
    require(len(charts) == len(twists) == 120, 'full chart census')
    for p in charts:
        require(all(projectors[0][p[i]][p[j]] == projectors[0][i][j] for i in range(12) for j in range(12)), 'chart preserves slow band')
    twist = twists[0]
    require(all(projectors[0][twist[i]][twist[j]] == projectors[1][i][j] for i in range(12) for j in range(12)), 'twist exchanges bands')
    inverse = [twist.index(i) for i in range(12)]
    require(all(chars[actions.index(tuple(twist[p[inverse[i]]] for i in range(12)))][0] == chars[k][1]
                for k,p in enumerate(actions)), 'outer action twist character')
    generator = source_marked_action(frame, vertices, relabel)
    require(generator in actions, '#1088 action belongs to live response action')
    marked = actions.index(generator)
    return encode({'schema':'oph.marked-triplet-attachment.v1', 'baseline':BASELINE,
        'integration_baseline':INTEGRATION_BASELINE,
        'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in INPUTS},
        'native_to_response_ports':relabel, 'projectors':projectors,
        'actions':[list(p) for p in actions], 'characters':chars,
        'fixed_action_hom_dimensions':{'slow_to_upper':1, 'slow_to_lower':0, 'fast_to_lower':1},
        'source_marked_order_five_action':list(generator),
        'marking_reference':'code/a5_closure/manifests/galois_port_frame_reference.json:generator',
        'marked_action_index':marked, 'marked_traces':chars[marked],
        'incidence_recharting_count':len(charts), 'incidence_rechartings_exchanging_bands':0,
        'outer_twist_example':list(twist), 'edge_distance_two_exchange_count':len(twists),
        'scalar_clock_rescaling_changes_representation':False,
        'scope':'fixed marked native incidence action; scalar spectral/clock reparameterizations and incidence rechartings',
        'missing_datum':'source-derived action twist or explicit representation attachment across spatial and internal domains',
        'verdict':'SPATIAL_SLOW_AND_INTERNAL_LOWER_TRIPLETS_NOT_SAME_MARKED_OBJECT',
        'physical_golden_embedding_selected':False, 'public_data_selection_input':False})


if __name__ == '__main__':
    print(json.dumps(build(), indent=2, sort_keys=True))
