"""Independent fixed-action audit in Q[phi]/(phi^2-phi-1).

The producer uses distance-table projectors, Hom equation ranks and graph
backtracking. This replay uses Laplacian interpolation, character inner
products and the unique distance signatures of an ordered triangle. It
imports only the separate Galois-port verifier's exact field arithmetic.
"""
import ast
from collections import deque
from fractions import Fraction as Q
import hashlib
from itertools import permutations
import json
from pathlib import Path
import sys

from verify_galois_port_frame_independent import (
    ZERO, ONE, PHI, add, sub, mul, div, sigma, scalar, total, read, matrix,
    dot, determinant, need,
)

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT/'code/a5_closure/manifests/marked_triplet_attachment_reference.json'
INCIDENCE = 'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'
RESPONSE = 'code/source_selection_model/response.json'
INPUTS = [INCIDENCE, RESPONSE, 'model/common_source/slow_geometry.py',
          'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
          'code/source_selection_model/verify_response.py',
          'code/a5_closure/manifests/galois_port_frame_reference.json']
BASELINE = '879ed37ad2c1933ecb6dd8f54271362ee3f4d41c'
INTEGRATION_BASELINE = '7b0997b4a8a425679d9e2946f82c86cd92289059'
S = sub(mul(scalar(2), PHI), ONE)
KEYS = {'schema', 'baseline', 'integration_baseline', 'source_sha256', 'native_to_response_ports',
        'projectors', 'actions', 'characters', 'fixed_action_hom_dimensions',
        'source_marked_order_five_action', 'marking_reference', 'marked_action_index', 'marked_traces', 'incidence_recharting_count',
        'incidence_rechartings_exchanging_bands', 'outer_twist_example',
        'edge_distance_two_exchange_count', 'scalar_clock_rescaling_changes_representation',
        'scope', 'missing_datum', 'verdict', 'physical_golden_embedding_selected',
        'public_data_selection_input'}


def unique_pairs(rows):
    """Reject duplicate keys before any receipt or source interpretation."""
    result = {}
    for key, value in rows:
        if key in result:
            raise ValueError('duplicate JSON key: '+key)
        result[key] = value
    return result


def load(path):
    def reject_constant(value):
        raise ValueError('nonfinite JSON constant: '+value)
    return json.loads(Path(path).read_text(), object_pairs_hook=unique_pairs,
                      parse_constant=reject_constant)


def mm(a, b):
    return [[dot(row, col) for col in zip(*b)] for row in a]


def eye(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def diagonal_shift(a, x):
    return [[sub(t, x) if i == j else t for j, t in enumerate(row)]
            for i, row in enumerate(a)]


def incidence(root):
    tree = ast.parse((root/INCIDENCE).read_text())
    fs = next(ast.literal_eval(node.value) for node in tree.body
              if isinstance(node, ast.Assign) and
              any(isinstance(t, ast.Name) and t.id == 'ORIENTED_BASE_FACES'
                  for t in node.targets))
    need(len(fs) == 20 and len(set(fs)) == 20, 'twenty distinct live faces')
    edges = {frozenset((f[i], f[(i+1) % 3])) for f in fs for i in range(3)}
    need(len(edges) == 30, 'thirty live edges')
    adjacent = [{j for j in range(12) if frozenset((i,j)) in edges} for i in range(12)]
    ds = []
    for i in range(12):
        row = [None]*12
        row[i] = 0
        queue = deque([i])
        while queue:
            a = queue.popleft()
            for b in adjacent[a]:
                if row[b] is None:
                    row[b] = row[a]+1
                    queue.append(b)
        need(sorted(row) == [0]+[1]*5+[2]*5+[3], 'native distance census')
        ds.append(row)
    return fs, ds


def interpolate(lap, eigenvalue):
    """Lagrange spectral idempotent, using the other three exact bands."""
    result = eye(12)
    for other in (ZERO, sigma(eigenvalue), scalar(6)):
        factor = [[div(x, sub(eigenvalue, other)) for x in row]
                  for row in diagonal_shift(lap, other)]
        result = mm(result, factor)
    need(mm(result, result) == result, 'interpolated idempotent')
    need(total(result[i][i] for i in range(12)) == scalar(3), 'rank-three band')
    need(mm(lap, result) == [[mul(eigenvalue,x) for x in row] for row in result],
         'interpolated eigenspace')
    return result


def triangle_maps(ds, swap=False):
    """An ordered face's twelve distance signatures uniquely name all ports."""
    source = next(t for t in permutations(range(12), 3)
                  if all(ds[t[i]][t[j]] == 1 for i in range(3) for j in range(i)))
    signatures = [tuple(ds[i][a] for a in source) for i in range(12)]
    need(len(set(signatures)) == 12, 'triangle resolves all native vertices')
    change = {0:0, 1:2, 2:1, 3:3} if swap else {i:i for i in range(4)}
    result = []
    for target in permutations(range(12), 3):
        if not all(ds[target[i]][target[j]] == change[1]
                   for i in range(3) for j in range(i)):
            continue
        lookup = {tuple(change[ds[b][a]] for a in target):b for b in range(12)}
        need(len(lookup) == 12 and set(lookup) == set(signatures), 'target triangle resolves ports')
        p = tuple(lookup[s] for s in signatures)
        need(all(ds[p[i]][p[j]] == change[ds[i][j]] for i in range(12) for j in range(12)),
             'complete triangle extension')
        result.append(p)
    need(len(result) == len(set(result)) == 120, 'triangle map census')
    return sorted(result)


def source_frame(root):
    """Reconstruct the frame from the live GOLDEN constant and local bindings."""
    carrier = ast.parse((root/INPUTS[3]).read_text())
    golden = next(n.value for n in carrier.body if isinstance(n, ast.AnnAssign)
                  and isinstance(n.target, ast.Name) and n.target.id == 'GOLDEN')
    need(isinstance(golden, ast.Tuple) and len(golden.elts) == 2, 'live GOLDEN pair')
    def rational(node):
        need(isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
             and node.func.id == 'Fraction' and not node.keywords
             and all(isinstance(x, ast.Constant) and type(x.value) is int for x in node.args),
             'live GOLDEN rational')
        return Q(*(x.value for x in node.args))
    golden_pair = tuple(rational(x) for x in golden.elts)
    # carrier.QS5 uses a + b*sqrt(5); this verifier uses a + b*phi.
    native_phi = (golden_pair[0]-golden_pair[1], 2*golden_pair[1])
    need(native_phi == PHI, 'live GOLDEN equals frozen field generator')
    tree = ast.parse((root/'model/common_source/slow_geometry.py').read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'raw_frame')
    assignments = [n for n in function.body if isinstance(n, ast.Assign)]
    expected = ast.parse('c = source()\nzero, one, minus, phi = c.q5(0), c.q5(1), c.q5(-1), c.GOLDEN\nnegative_phi = c.q5_sub(zero, phi)').body
    need([ast.dump(n) for n in assignments] == [ast.dump(n) for n in expected],
         'live raw-frame bindings')
    returned = next(n.value for n in function.body if isinstance(n, ast.Return))
    names = {'zero':ZERO, 'one':ONE, 'minus':scalar(-1), 'phi':native_phi,
             'negative_phi':sub(ZERO,native_phi)}
    need(isinstance(returned, ast.Tuple) and len(returned.elts) == 12, 'live native frame shape')
    need(all(isinstance(row, ast.Tuple) and len(row.elts) == 3 and
             all(isinstance(x, ast.Name) and x.id in names for x in row.elts)
             for row in returned.elts), 'live native coordinate symbols')
    return [[names[x.id] for x in row.elts] for row in returned.elts]


def response_vertices():
    base = [[ZERO, scalar(s), mul(scalar(t),PHI)] for s in (1,-1) for t in (1,-1)]
    return base + [[y,z,x] for x,y,z in base] + [[z,x,y] for x,y,z in base]


def source_marked_action(root, vertices, relabel):
    """Independently transport the committed #1088 generator into this chart."""
    t, nt = PHI, sub(ZERO,PHI)
    z, o, no = ZERO, ONE, scalar(-1)
    old_vertices = [[z,o,t],[o,t,z],[t,z,o],[no,t,z],[z,no,t],[t,z,no],
                    [nt,z,o],[z,o,nt],[o,nt,z],[nt,z,no],[no,nt,z],[z,no,nt]]
    marker = load(root/INPUTS[5])
    old_to_response = [vertices.index(v) for v in old_vertices]
    response_to_old = [old_to_response.index(i) for i in range(12)]
    old = marker.get('generator')
    need(type(old) is list and len(old) == 12 and all(type(i) is int for i in old)
         and sorted(old) == list(range(12)), 'committed #1088 generator type')
    response_action = [old_to_response[old[response_to_old[i]]] for i in range(12)]
    response_to_native = [relabel.index(i) for i in range(12)]
    native_action = tuple(response_to_native[response_action[relabel[i]]]
                          for i in range(12))
    power = tuple(range(12))
    for n in range(1,6):
        power = tuple(native_action[power[i]] for i in range(12))
        need((power == tuple(range(12))) == (n == 5), 'transported #1088 order five')
    return native_action


def permutation(p):
    need(type(p) is list and len(p) == 12 and all(type(i) is int for i in p)
         and sorted(p) == list(range(12)), 'permutation shape')
    return tuple(p)


def verify(packet, root=ROOT, response=None):
    need(type(packet) is dict and set(packet) == KEYS, 'closed packet schema')
    need(all(type(packet[k]) is str for k in
             ('schema','baseline','scope','missing_datum','verdict','marking_reference')),
         'string field types')
    need(packet['marking_reference'] ==
         'code/a5_closure/manifests/galois_port_frame_reference.json:generator',
         'source marker provenance')
    need(packet['schema'] == 'oph.marked-triplet-attachment.v1' and
         packet['baseline'] == BASELINE and
         packet['integration_baseline'] == INTEGRATION_BASELINE, 'source and integration baselines')
    need(type(packet['source_sha256']) is dict and packet['source_sha256'] == {
        name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in INPUTS},
        'live-source byte pins')
    need(all(type(v) is str for v in packet['source_sha256'].values()), 'source hash types')
    fs, ds = incidence(root)
    lap = [[scalar(5 if i == j else -int(ds[i][j] == 1)) for j in range(12)] for i in range(12)]
    slow = interpolate(lap, sub(scalar(5),S))
    fast = interpolate(lap, add(scalar(5),S))
    ps = [slow, fast]
    need(type(packet['projectors']) is list and len(packet['projectors']) == 2 and
         [matrix(p,12,12) for p in packet['projectors']] == ps, 'interpolated projector replay')
    frame, vertices = source_frame(root), response_vertices()
    relabel = [vertices.index(v) for v in frame]
    need(type(packet['native_to_response_ports']) is list and
         all(type(x) is int for x in packet['native_to_response_ports']) and
         packet['native_to_response_ports'] == relabel, 'native response port attachment')
    norm = dot(frame[0],frame[0])
    need([[div(dot(x,y),mul(scalar(4),norm)) for y in frame] for x in frame] == slow,
         'source frame realizes the native slow eigenspace')
    response = load(root/RESPONSE) if response is None else response
    rows = response['closed_paths']
    need(type(rows) is list and len(rows) == 60, 'sixty live response paths')
    actions, internal, spatial = [], [], []
    blocks_by_action = {}
    for row in rows:
        action = permutation(row['port_action'])
        p = tuple(relabel.index(action[relabel[i]]) for i in range(12))
        encoded = row['observed_response']
        need(len(encoded) == 2 and all(len(b) == 3 and all(len(r) == 3 for r in b) for b in encoded),
             'two internal triplets')
        need(all(type(x) is list and len(x) == 4 and x[2:] == ['0','0']
                 for b in encoded for r in b for x in r), 'real response elements')
        a,b = [matrix([[x[:2] for x in r] for r in block],3,3) for block in encoded]
        need(b == [[sigma(x) for x in r] for r in a], 'Galois-paired internal triplets')
        need(determinant(a) == ONE and mm(a,list(map(list,zip(*a)))) == eye(3),
             'proper orthogonal internal action')
        need(all([dot(r,v) for r in a] == vertices[action[i]] for i,v in enumerate(vertices)),
             'independent upper action on all response vertices')
        need(all(ds[p[i]][p[j]] == ds[i][j] for i in range(12) for j in range(12)),
             'fixed native incidence action')
        actions.append(p)
        blocks_by_action[p] = (a,b)
        internal.append([total(m[i][i] for i in range(3)) for m in (a,b)])
        spatial.append([total(projector[i][p[i]] for i in range(12)) for projector in ps])
    need(len(set(actions)) == 60 and spatial == internal, 'all marked character attachments')
    oriented = {tuple(face[k:]+face[:k]) for face in map(list,fs) for k in range(3)}
    need(all(all(tuple(p[i] for i in face) in oriented for face in fs) for p in actions),
         'sixty orientation-preserving incidence actions')
    for p in actions:
        for q in actions:
            pq = tuple(p[q[i]] for i in range(12))
            need(pq in blocks_by_action, 'action group closure')
            need(all(mm(x,y) == z for x,y,z in zip(blocks_by_action[p],blocks_by_action[q],
                                                   blocks_by_action[pq],strict=True)),
                 'both representations compose')
    need(type(packet['actions']) is list and
         [permutation(p) for p in packet['actions']] == actions, 'marked action list')
    need(type(packet['characters']) is list and len(packet['characters']) == 60,
         'character table JSON shape')
    need([[read(x) for x in row] for row in packet['characters']] == spatial,
         'character table replay')
    hom = [[div(total(mul(s[i],t[j]) for s,t in zip(spatial,internal,strict=True)),scalar(60))
            for j in range(2)] for i in range(2)]
    need(hom == [[ONE,ZERO],[ZERO,ONE]], 'character inner-product Hom dimensions')
    need(type(packet['fixed_action_hom_dimensions']) is dict and
         all(type(value) is int for value in packet['fixed_action_hom_dimensions'].values()) and
         set(packet['fixed_action_hom_dimensions']) ==
         {'slow_to_upper','slow_to_lower','fast_to_lower'} and
         packet['fixed_action_hom_dimensions'] ==
         {'slow_to_upper':1,'slow_to_lower':0,'fast_to_lower':1}, 'Hom dimensions')
    charts, twists = triangle_maps(ds), triangle_maps(ds,True)
    need(all(all(slow[p[i]][p[j]] == slow[i][j] for i in range(12) for j in range(12))
             for p in charts), 'all incidence charts preserve slow band')
    need(all(all(slow[p[i]][p[j]] == fast[i][j] for i in range(12) for j in range(12))
             for p in twists), 'all edge exchanges swap bands')
    twist = permutation(packet['outer_twist_example'])
    need(twist == twists[0] and twist not in charts, 'twist is not an incidence recharting')
    inverse = [twist.index(i) for i in range(12)]
    need(all(spatial[actions.index(tuple(twist[p[inverse[i]]] for i in range(12)))][0]
             == spatial[k][1] for k,p in enumerate(actions)), 'outer-twisted character attachment')
    source_generator = source_marked_action(root, vertices, relabel)
    need(source_generator in actions, '#1088 action is present in the live action group')
    marked = actions.index(source_generator)
    need(permutation(packet['source_marked_order_five_action']) == source_generator,
         'source-frozen order-five permutation')
    need(type(packet['marked_action_index']) is int and packet['marked_action_index'] == marked,
         'incidence-defined marked action index')
    need([read(x) for x in packet['marked_traces']] == spatial[marked], 'computed marked traces')
    need(type(packet['incidence_recharting_count']) is int and
         type(packet['incidence_rechartings_exchanging_bands']) is int and
         type(packet['edge_distance_two_exchange_count']) is int and
         packet['incidence_recharting_count'] == 120 and
         packet['incidence_rechartings_exchanging_bands'] == 0 and
         packet['edge_distance_two_exchange_count'] == 120, 'exhaustive finite counts')
    need(type(packet['scalar_clock_rescaling_changes_representation']) is bool and
         type(packet['physical_golden_embedding_selected']) is bool and
         type(packet['public_data_selection_input']) is bool, 'boolean field types')
    need(packet['scope'] == 'fixed marked native incidence action; scalar spectral/clock reparameterizations and incidence rechartings',
         'bounded nonselection scope')
    need(packet['missing_datum'] == 'source-derived action twist or explicit representation attachment across spatial and internal domains',
         'exact missing attachment')
    need(packet['verdict'] == 'SPATIAL_SLOW_AND_INTERNAL_LOWER_TRIPLETS_NOT_SAME_MARKED_OBJECT', 'typed verdict')
    need(packet['scalar_clock_rescaling_changes_representation'] is False and
         packet['physical_golden_embedding_selected'] is False and
         packet['public_data_selection_input'] is False, 'no physical promotion or target input')
    return {'verified':True, 'proper_actions':60, 'incidence_rechartings':120,
            'edge_distance_two_exchanges':120, 'same_action_Hom_dimension':0,
            'method':'Laplacian interpolation and character inner products over Q(phi)'}


if __name__ == '__main__':
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else REFERENCE
    print(json.dumps(verify(load(path)),indent=2,sort_keys=True))
