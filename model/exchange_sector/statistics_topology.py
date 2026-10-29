"""Exact topology of the native two-particle path complex.

The existing 66/300/495 complex X is the two-skeleton of D(K)/swap,
where D(K) consists of products of vertex-disjoint simplices of the native
triangulated sphere K.  This module checks that identification, including
integer cell orientations.  The extra cells are a proof device, not extra
native operations or extra relations on paths.

Here is the continuous argument, independent of floating rank estimates.
For distinct points with global barycentric coordinates l and r, put
m_i=min(l_i,r_i), M=sum(m_i)<1.  The pair

    ((l-t*m)/(1-t*M), (r-t*m)/(1-t*M)),  0 <= t <= 1,

stays in the original two simplices, remains distinct, and ends with
disjoint supports.  It fixes D(K) and commutes with swap.  Thus it is an
equivariant strong deformation retraction of F_2(|K|) onto D(K).
The exact connected, closed, oriented surface/link certificate below and
Euler characteristic two identify |K| with S^2 by surface classification.
For x,y in S^2, let m=(x+y)/2, d=(x-y)/2.  They are orthogonal, d != 0,
and ((1-t)*m +/- d)/sqrt((1-t)^2*||m||^2+||d||^2) gives an equivariant
strong deformation retraction onto antipodal pairs.  The quotient of that
last space is RP^2.  Consequently pi_1(D(K)/swap)=Z/2.  Attaching cells of
dimension at least three does not change pi_1, so pi_1(X)=Z/2 as well.

Only the full deleted product has the homotopy type of RP^2: X itself has
Euler characteristic 261.  The result is a fundamental-group computation,
not an assertion that the existing two-complex has no higher topology.

Every finite-dimensional flat unitary transport on X, up to vertexwise
unitary gauge, is therefore determined by one U with U^2=I.  Unitarity
makes U self-adjoint, so it is a direct sum of +1 and -1 characters.
The ordered lift detects the nontrivial generator.  No higher-dimensional
irreducible exchange representation or arbitrary U(1) exchange angle is
available for this particular two-particle deformation domain.  The
topology admits both characters; it does not choose a quantum instrument.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product

from ..clock_defect_energy import defect_history as history
from ..observer import native_loop_readout


def _connected(vertices, neighbors):
    reached = {next(iter(vertices))}
    queue = list(reached)
    for vertex in queue:
        for other in neighbors[vertex]:
            if other not in reached:
                reached.add(other)
                queue.append(other)
    return reached == set(vertices)


def sphere_triangulation_certificate(vertices, edges, faces):
    """Check the exact finite hypotheses of the closed-surface theorem.

    This validator also accepts altered incidence data for negative controls.
    Face input orientations are immaterial; coherent orientations are found
    from incidence.  The conclusion uses classification of closed surfaces,
    rather than treating an Euler count alone as proof of a sphere.
    """
    if (type(vertices) is not tuple or not vertices
            or any(type(v) is not int for v in vertices)
            or len(set(vertices)) != len(vertices)):
        raise ValueError("distinct integer vertices required")
    vertex_set = set(vertices)
    if type(edges) is not tuple or type(faces) is not tuple:
        raise ValueError("tuple edge and triangle inventories required")
    if any(type(e) is not tuple or len(e) != 2
           or any(type(v) is not int or v not in vertex_set for v in e)
           or e[0] >= e[1] for e in edges):
        raise ValueError("sorted distinct-vertex edges required")
    if len(set(edges)) != len(edges):
        raise ValueError("duplicate edge")
    if any(type(f) is not tuple or len(f) != 3
           or any(type(v) is not int or v not in vertex_set for v in f)
           or len(set(f)) != 3 for f in faces):
        raise ValueError("three distinct vertices per face required")
    triangles = tuple(sorted(tuple(sorted(f)) for f in faces))
    if len(set(triangles)) != len(triangles):
        raise ValueError("duplicate triangle")
    edge_faces = {e: [] for e in edges}
    neighbors = {v: set() for v in vertices}
    for a, b in edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
    if not _connected(vertices, neighbors):
        raise ValueError("surface one-skeleton is disconnected")
    for i, triangle in enumerate(triangles):
        for j in range(3):
            edge = triangle[:j] + triangle[j+1:]
            if edge not in edge_faces:
                raise ValueError("face boundary edge absent")
            edge_faces[edge].append((i, (-1)**j))
    if any(len(incident) != 2 for incident in edge_faces.values()):
        raise ValueError("each edge must have two incident triangles")
    link_lengths = []
    for vertex in vertices:
        link_edges = [tuple(v for v in f if v != vertex)
                      for f in triangles if vertex in f]
        link_vertices = set().union(*(set(e) for e in link_edges))
        links = {v: set() for v in link_vertices}
        for a, b in link_edges:
            links[a].add(b)
            links[b].add(a)
        if (len(link_vertices) < 3
                or any(len(adjacent) != 2 for adjacent in links.values())
                or not _connected(link_vertices, links)):
            raise ValueError("every vertex link must be a connected circle")
        link_lengths.append(len(link_vertices))
    face_neighbors = defaultdict(list)
    for (a, sa), (b, sb) in edge_faces.values():
        face_neighbors[a].append((b, -sa*sb))
        face_neighbors[b].append((a, -sa*sb))
    orientations = {0: 1}
    queue = [0]
    for face in queue:
        for other, change in face_neighbors[face]:
            sign = orientations[face]*change
            if other in orientations and orientations[other] != sign:
                raise ValueError("surface has no coherent face orientation")
            if other not in orientations:
                orientations[other] = sign
                queue.append(other)
    if len(orientations) != len(triangles):
        raise ValueError("triangle adjacency graph is disconnected")
    euler = len(vertices)-len(edges)+len(triangles)
    if euler != 2:
        raise ValueError("closed oriented surface must have Euler characteristic two")
    return {"vertices": len(vertices), "edges": len(edges), "faces": len(faces),
            "euler_characteristic": euler,
            "vertex_link_lengths": tuple(link_lengths),
            "coherent_face_orientations": tuple(orientations[i] for i in range(len(triangles))),
            "connected_closed_oriented_surface": True,
            "surface_classification_conclusion": "S^2"}


@lru_cache(maxsize=1)
def simplices():
    edges, faces = native_loop_readout.carrier()
    sphere_triangulation_certificate(tuple(range(12)), edges, faces)
    return (tuple((v,) for v in range(12)), edges,
            tuple(sorted(tuple(sorted(face)) for face in faces)))


def _simplex_key(simplex):
    return len(simplex), simplex


def _cell(value):
    if type(value) is not tuple or len(value) != 2:
        raise ValueError("a pair of disjoint native simplices required")
    inventory = simplices()
    for simplex in value:
        if (type(simplex) is not tuple or not 1 <= len(simplex) <= 3
                or any(type(v) is not int for v in simplex)
                or simplex not in inventory[len(simplex)-1]):
            raise ValueError("a sorted nonempty native simplex required")
    if set(value[0]) & set(value[1]):
        raise ValueError("cell simplices must be vertex-disjoint")
    return value


def quotient_cell(cell):
    """Project an oriented ordered product cell, retaining its Koszul sign."""
    left, right = _cell(cell)
    if _simplex_key(left) < _simplex_key(right):
        return (left, right), 1
    return (right, left), (-1)**((len(left)-1)*(len(right)-1))


@lru_cache(maxsize=2)
def _deleted_cells(ordered):
    result = [[] for _ in range(5)]
    flat = tuple(simplex for layer in simplices() for simplex in layer)
    for left, right in product(flat, repeat=2):
        if set(left).isdisjoint(right) and (ordered or _simplex_key(left) < _simplex_key(right)):
            result[len(left)+len(right)-2].append((left, right))
    return tuple(tuple(sorted(layer)) for layer in result)


def deleted_cells(*, ordered=False):
    """The actual product cells, indexed by dimension zero through four."""
    if type(ordered) is not bool:
        raise ValueError("ordered must be a boolean")
    return _deleted_cells(ordered)


@lru_cache(maxsize=None)
def _boundary(cell, ordered):
    left, right = cell
    terms = Counter()
    for side, simplex in enumerate(cell):
        if len(simplex) == 1:
            continue
        for i in range(len(simplex)):
            face = simplex[:i]+simplex[i+1:]
            child = (face, right) if side == 0 else (left, face)
            coefficient = (-1)**(i+(len(left)-1 if side == 1 else 0))
            if not ordered:
                child, sign = quotient_cell(child)
                coefficient *= sign
            terms[child] += coefficient
    return tuple(sorted((child, value) for child, value in terms.items() if value))


def cell_boundary(cell, *, ordered=False):
    """Sparse integral boundary in the chosen oriented cell basis."""
    if type(ordered) is not bool:
        raise ValueError("ordered must be a boolean")
    cell = _cell(cell)
    if not ordered and quotient_cell(cell)[0] != cell:
        raise ValueError("unordered cell must use its canonical representative")
    return _boundary(cell, ordered)


def _chain_certificate(ordered):
    layers = deleted_cells(ordered=ordered)
    for dimension, layer in enumerate(layers):
        for cell in layer:
            twice = Counter()
            for face, coefficient in cell_boundary(cell, ordered=ordered):
                if dimension == 0 or face not in layers[dimension-1]:
                    raise ArithmeticError("boundary left its cell inventory")
                for edge, second in cell_boundary(face, ordered=ordered):
                    twice[edge] += coefficient*second
            if any(twice.values()):
                raise ArithmeticError("integer product boundary does not square to zero")
    return tuple(map(len, layers))


def _skeleton_certificate():
    pairs, edges, cells, _, _, _ = history.geometry()
    index = {pair: i for i, pair in enumerate(pairs)}
    layers = deleted_cells()

    def vertex(cell):
        return index[tuple(sorted((cell[0][0], cell[1][0])))]

    edge_projection = {}
    for cell in layers[1]:
        endpoints = {vertex(face): coefficient for face, coefficient in cell_boundary(cell)}
        a, b = sorted(endpoints)
        if endpoints[a] != -endpoints[b] or abs(endpoints[b]) != 1:
            raise ArithmeticError("invalid product-edge endpoints")
        edge_projection[cell] = (edges.index((a, b)), endpoints[b])
    if set(index) != {tuple(sorted((a[0], b[0]))) for a, b in layers[0]}:
        raise ArithmeticError("deleted-product vertices differ from native configurations")
    if {i for i, _ in edge_projection.values()} != set(range(len(edges))):
        raise ArithmeticError("deleted-product edges differ from native movements")
    native_cells = {frozenset(vertices): history.boundary(vertices, edges) for _, vertices in cells}
    if len(native_cells) != len(cells):
        raise ArithmeticError("native polygon vertices do not distinguish cells")
    seen = set()
    for cell in layers[2]:
        vertices = frozenset(index[tuple(sorted(pair))] for pair in product(*cell))
        chain = [0]*len(edges)
        for edge, coefficient in cell_boundary(cell):
            i, sign = edge_projection[edge]
            chain[i] += coefficient*sign
        expected = native_cells.get(vertices)
        if expected is None or tuple(chain) not in (expected, tuple(-v for v in expected)):
            raise ArithmeticError("oriented two-cell differs from an existing native polygon")
        seen.add(vertices)
    if seen != set(native_cells):
        raise ArithmeticError("deleted product omitted a native local cell")
    return {"vertices": len(layers[0]), "edges": len(layers[1]), "two_cells": len(layers[2]),
            "all_oriented_boundaries_match_up_to_cell_orientation": True}


def _cover_certificate():
    layers = deleted_cells(ordered=True)
    quotient = deleted_cells()
    for dimension, layer in enumerate(layers):
        multiplicities = Counter(quotient_cell(cell)[0] for cell in layer)
        if set(multiplicities) != set(quotient[dimension]) or set(multiplicities.values()) != {2}:
            raise ArithmeticError("ordered deleted product is not a free two-sheeted cover")
        for cell in layer:
            base, sign = quotient_cell(cell)
            projected = Counter()
            for child, coefficient in cell_boundary(cell, ordered=True):
                target, orientation = quotient_cell(child)
                projected[target] += coefficient*orientation
            expected = dict((child, sign*coefficient) for child, coefficient in cell_boundary(base))
            if {child: value for child, value in projected.items() if value} != expected:
                raise ArithmeticError("oriented cover projection is not a chain map")
    neighbors = {cell: set() for cell in layers[0]}
    pairs = history.configurations()
    for cell in layers[1]:
        endpoints = tuple(child for child, _ in cell_boundary(cell, ordered=True))
        a, b = endpoints
        neighbors[a].add(b)
        neighbors[b].add(a)
        charts = [int(left[0] > right[0]) for left, right in endpoints]
        base = [pairs.index(tuple(sorted((left[0], right[0])))) for left, right in endpoints]
        if charts[0] ^ charts[1] != history.edge_parity(*base):
            raise ArithmeticError("cover lift disagrees with the native ordered-particle parity")
    if not _connected(layers[0], neighbors):
        raise ArithmeticError("ordered configuration cover is disconnected")
    return {"ordered_cell_counts": tuple(map(len, layers)), "sheets": 2,
            "connected": True, "projection_is_integral_chain_map": True,
            "native_lift_matches_every_ordered_edge": True}


def topology_certificate():
    """Compact exact finite premises plus the continuous theorem conclusions."""
    edges, faces = native_loop_readout.carrier()
    surface = sphere_triangulation_certificate(tuple(range(12)), edges, faces)
    unordered = _chain_certificate(False)
    ordered = _chain_certificate(True)
    if unordered != (66, 300, 495, 360, 100) or ordered != tuple(2*n for n in unordered):
        raise ArithmeticError("native deleted-product cell census changed")
    return {"surface": surface, "unordered_cell_counts": unordered,
            "integer_boundary_squared_zero_ordered_and_unordered": True,
            "native_two_skeleton": _skeleton_certificate(), "ordered_cover": _cover_certificate(),
            "full_deleted_product_euler_characteristic": sum((-1)**i*n for i, n in enumerate(unordered)),
            "native_two_skeleton_euler_characteristic": sum((-1)**i*n for i, n in enumerate(unordered[:3])),
            "native_fundamental_group": "Z/2",
            "ordered_cover_fundamental_group": "trivial",
            "full_deleted_product_homotopy_type": "RP^2",
            "unitary_irreducible_characters": (1, -1),
            "native_quantum_character_selected": False,
            "negative_control_remove_triangle_relations": graph_only_certificate()}


def graph_only_certificate():
    """Changed-domain control: omit the declared one-particle face fillings.

    We compute H^1 over F_2, not the graph-only fundamental group.  An
    explicit integer edge cochain additionally exhibits a continuous U(1)
    family: its sum vanishes around every disjoint square, but equals one
    around a removed spectator triangle.  Exponentiating by any angle is
    flat on the changed domain and violates the original triangle relation
    unless that angle is an integer multiple of 2*pi.
    """
    pairs, edges, cells, _, chords, _ = history.geometry()
    squares = tuple(history.boundary(vertices, edges)
                    for kind, vertices in cells if kind == "disjoint_square")
    vertices = next(vertices for kind, vertices in cells if kind == "one_particle_face")
    before, after = pairs[vertices[0]], pairs[vertices[1]]
    origin, destination = next(iter(set(before)-set(after))), next(iter(set(after)-set(before)))
    cochain = []
    for i, j in edges:
        removed = next(iter(set(pairs[i])-set(pairs[j])))
        added = next(iter(set(pairs[j])-set(pairs[i])))
        cochain.append(int((removed, added) == (origin, destination))
                       - int((removed, added) == (destination, origin)))
    if any(sum(a*b for a, b in zip(cochain, chain)) for chain in squares):
        raise ArithmeticError("one-particle edge phase is not flat on disjoint squares")
    triangle = history.boundary(vertices, edges)
    exponent = sum(a*b for a, b in zip(cochain, triangle))
    rank = history.gf2_rank(squares)
    if exponent != 1 or rank != 215:
        raise ArithmeticError("triangle-deletion negative control changed")
    return {"retained_square_cells": len(squares), "omitted_triangle_cells": 180,
            "cycle_rank": len(chords), "square_boundary_rank_mod2": rank,
            "first_cohomology_mod2_dimension": len(chords)-rank,
            "graph_only_fundamental_group_computed": False,
            "continuous_character_witness": {
                "directed_one_particle_edge": (origin, destination),
                "all_square_exponents": 0,
                "removed_triangle_occupied_pairs": tuple(pairs[v] for v in vertices+vertices[:1]),
                "removed_triangle_exponent": exponent,
                "holonomy": "exp(i*theta)",
                "flat_on_original_domain_only_when": "theta is an integer multiple of 2*pi"},
            "changes_declared_deformation_contract": True}


def certificate():
    """Deterministic public entry point for the statistics packet."""
    return topology_certificate()


def _rational(value):
    if type(value) is not int and type(value) is not F:
        raise ValueError("exact integer or Fraction coordinates required")
    return F(value)


def _time(value):
    value = _rational(value)
    if not 0 <= value <= 1:
        raise ValueError("deformation time must lie in [0,1]")
    return value


def _barycentric(point):
    if type(point) is not tuple or len(point) != 12:
        raise ValueError("twelve barycentric coordinates required")
    point = tuple(map(_rational, point))
    if any(x < 0 for x in point) or sum(point) != 1:
        raise ValueError("nonnegative barycentric coordinates summing to one required")
    support = tuple(i for i, x in enumerate(point) if x)
    if len(support) > 3 or support not in simplices()[len(support)-1]:
        raise ValueError("barycentric support must be a native simplex")
    return point


def deleted_product_deformation(left, right, time):
    """Evaluate the exact rational restriction of the continuous PL retraction."""
    left, right, time = _barycentric(left), _barycentric(right), _time(time)
    if left == right:
        raise ValueError("configuration points must be distinct")
    common = tuple(min(a, b) for a, b in zip(left, right))
    denominator = 1-time*sum(common)
    return (tuple((a-time*m)/denominator for a, m in zip(left, common)),
            tuple((b-time*m)/denominator for b, m in zip(right, common)))


def antipodal_deformation(left, right, time):
    """Exact sphere homotopy on rational unit inputs, without rounded radicals.

    Return (left_numerator, right_numerator, common_norm_squared).  Actual
    unit vectors are both numerators divided by sqrt(common_norm_squared).
    The docstring proof applies to all real sphere inputs; this evaluator
    supplies exact checks on rational points, not a sampled topological proof.
    """
    def unit(point):
        if type(point) is not tuple or len(point) != 3:
            raise ValueError("three unit-sphere coordinates required")
        point = tuple(map(_rational, point))
        if sum(x*x for x in point) != 1:
            raise ValueError("unit-sphere point required")
        return point
    left, right, time = unit(left), unit(right), _time(time)
    if left == right:
        raise ValueError("sphere configuration points must be distinct")
    middle = tuple((a+b)/2 for a, b in zip(left, right))
    difference = tuple((a-b)/2 for a, b in zip(left, right))
    a = tuple((1-time)*m+d for m, d in zip(middle, difference))
    b = tuple((1-time)*m-d for m, d in zip(middle, difference))
    norm = sum(x*x for x in a)
    if norm <= 0 or norm != sum(x*x for x in b):
        raise ArithmeticError("antipodal homotopy normalization failed")
    return a, b, norm
