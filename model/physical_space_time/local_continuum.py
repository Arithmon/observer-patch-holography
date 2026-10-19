"""Compact, strongly local native spatial limits with spectral dimension two.

The source topology and its ancestry dictionary are the existing ones.  The
result upgrades their common weak Green subsequence to operator-norm
convergence, with an injective compact limit and a positive inverse generator.
Full-port energy lifts supply its conservative strongly local Dirichlet form;
an explicit source interpolation identifies its regular H1 form core, and
triangle tests bound both sides of its spectral growth. A selected intrinsic
metric and physical spatial identification are not asserted.
The all-level proof is in OCCUPATION_FIELD.md.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path
import sympy as sp

from . import continuum_field as continuum
from . import source_green as source
from . import uniform_green as uniform
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*continuum.SOURCE_PATHS,
    'model/physical_space_time/local_continuum.py'))))
GREEN_BOUND = uniform.GREEN_BOUND
PATCH_POINCARE_CONSTANT = F(26)
COARSE_ENERGY_TAIL_CONSTANT = F(624)


def _level(value):
    if type(value) is not int or value < 0:
        raise ValueError('an original nonnegative integer is required')
    return value


def compactness_tail_squared(coarse_level):
    """Squared norm bound for (I-P_k)G_l, at every l>=k."""
    return COARSE_ENERGY_TAIL_CONSTANT*GREEN_BOUND/4**_level(coarse_level)


def covariance_eigenvalue_tail(coarse_level):
    """Bound on positive eigenvalue number 20*4^k and every later one."""
    return COARSE_ENERGY_TAIL_CONSTANT/4**_level(coarse_level)


def continuum_eigenvalue_enclosure(index):
    """All positive eigenvalues, counted from one, for every selected limit."""
    if type(index) is not int or index < 1:
        raise ValueError('an original positive integer eigenvalue index is required')
    return {'lower': F(index, 49920), 'upper': F(3*index, 2)}


def tent_eigenvalue_bound(coarse_level):
    """Min-max bound from mean-zero combinations of ancestral triangle tents."""
    scale = 4**_level(coarse_level)
    return {'index': 20*scale-1, 'upper': F(3*scale, 2)}


def harmonic_lift_tail(level):
    """Squared internal fluctuation per unit scaled full-port energy."""
    return F(1, 2*4**_level(level))


def cycle_variance_control(valence):
    """Exact rational PSD certificate for variance <= edge energy on C5/C6."""
    if type(valence) is not int or valence not in (5, 6):
        raise ValueError('a source vertex has original integer valence five or six')
    laplacian = sp.zeros(valence)
    for i in range(valence):
        j = (i+1) % valence
        laplacian[i, i] += 1; laplacian[j, j] += 1
        laplacian[i, j] -= 1; laplacian[j, i] -= 1
    residual = laplacian-sp.eye(valence)+sp.ones(valence)/valence
    pivots = []
    for _ in range(valence):
        pivot = residual[0, 0]
        if pivot < 0 or (pivot == 0 and any(residual[0, j] != 0
                                            for j in range(residual.cols))):
            raise ArithmeticError('the exact vertex-cycle variance certificate failed')
        pivots.append(pivot)
        tail = residual[1:, 1:]
        if pivot:
            tail -= residual[1:, :1]*residual[:1, 1:]/pivot
        residual = tail
    return {'valence': valence, 'rational_nonnegative_pivots': [str(x) for x in pivots],
        'variance_is_bounded_by_cycle_edge_energy': True}


def interpolation_geometry_control():
    """Exact six-triangle energy and mass identities for the source H1 lift."""
    vertex, midpoint, centroid = sp.symbols('vertex midpoint centroid', real=True)
    root = sp.sqrt(3)
    # One subtriangle has coordinates (0,0), (1/2,0), (1/2,sqrt(3)/6).
    gradient = sp.Matrix([2*(midpoint-vertex), 2*root*(centroid-midpoint)])
    integrated = sp.expand(root/24*(gradient.T*gradient)[0])
    target = root/6*(vertex-midpoint)**2+root/2*(centroid-midpoint)**2
    if sp.simplify(integrated-target) != 0:
        raise ArithmeticError('the exact six-triangle gradient identity changed')
    x, y = sp.symbols('x y', nonnegative=True)
    error = (1-x-y)*(vertex-centroid)+x*(midpoint-centroid)
    # The reference simplex has area1/2; division gives its normalized mean.
    mean_square = 2*sp.integrate(sp.integrate(error**2, (y, 0, 1-x)), (x, 0, 1))
    expected = ((vertex-centroid)**2+(midpoint-centroid)**2
                +(vertex-centroid)*(midpoint-centroid))/6
    if sp.simplify(mean_square-expected) != 0:
        raise ArithmeticError('the exact affine interpolation mass identity changed')
    return {'unit_side_subtriangle_energy': 'sqrt(3)/6*(v-m)^2+sqrt(3)/2*(c-m)^2',
        'subtriangle_mean_squared_error': '((v-c)^2+(m-c)^2+(v-c)*(m-c))/6',
        'source_vertex_cycle_certificates': [cycle_variance_control(n) for n in (5, 6)],
        'normalized_gradient_times_side_squared_per_dual_energy_upper': '7/30',
        'normalized_mass_error_times_cell_count_per_dual_energy_upper': '5/18'}


def source_interpolation_control(level):
    """Bind the interpolation, vertex cycles and all edge sums to native joins."""
    if _level(level) > 3:
        raise ValueError('source interpolation controls use levels zero through three')
    vertices, faces = uniform.affine_prefix(level)
    edge_cells, vertex_cells = {}, {v: set() for v in range(len(vertices))}
    for cell, face in enumerate(faces):
        for v in face:
            vertex_cells[v].add(cell)
        for a, b in zip(face, (*face[1:], face[0]), strict=True):
            edge_cells.setdefault(tuple(sorted((a, b))), []).append(cell)
    if any(len(cells) != 2 for cells in edge_cells.values()):
        raise ArithmeticError('source interpolation requires two cells at each primal edge')
    edges = set(uniform.dual_edges(faces))
    if edges != {tuple(sorted(cells)) for cells in edge_cells.values()}:
        raise ArithmeticError('interpolation and native dual-edge dictionaries differ')
    valences = Counter(map(len, vertex_cells.values()))
    if set(valences)-{5, 6} or valences[5] != 12:
        raise ArithmeticError('source interpolation vertex valences changed')
    adjacency = {cell: set() for cell in range(len(faces))}
    for a, b in edges:
        adjacency[a].add(b); adjacency[b].add(a)
    for star in vertex_cells.values():
        if any(len(adjacency[cell] & star) != 2 for cell in star):
            raise ArithmeticError('a source vertex star is not a cycle')
        reached, pending = set(), {next(iter(star))}
        while pending:
            cell = pending.pop()
            reached.add(cell)
            pending.update((adjacency[cell] & star)-reached)
        if reached != star:
            raise ArithmeticError('a source vertex cycle is disconnected')
    values = [F(17*cell % 23-11, 5) for cell in range(len(faces))]
    vertex_means = {v: sum(values[c] for c in cells)/len(cells)
                    for v, cells in vertex_cells.items()}
    dual_energy = sum((values[a]-values[b])**2 for a, b in edges)
    vm, cm, vc, mass = F(0), F(0), F(0), F(0)
    for edge, cells in edge_cells.items():
        midpoint = sum(values[c] for c in cells)/2
        for cell in cells:
            for vertex in edge:
                dv, dm = vertex_means[vertex]-values[cell], midpoint-values[cell]
                vm += (vertex_means[vertex]-midpoint)**2
                cm += dm*dm
                vc += dv*dv
                mass += (dv*dv+dm*dm+dv*dm)/(36*len(faces))
    gradient = vm/30+cm/10  # Normalized integral times original side squared.
    if (cm != dual_energy or vm > 4*dual_energy or vc > 4*dual_energy
            or gradient > F(7, 30)*dual_energy
            or mass > F(5, 18*len(faces))*dual_energy):
        raise ArithmeticError('the native interpolation energy or mass bound failed')
    return {'level': level, 'cells': len(faces),
        'vertex_valence_counts': {str(n): count for n, count in sorted(valences.items())},
        'all_native_vertex_stars_are_the_certified_cycles': True,
        'exact_dual_energy': str(dual_energy),
        'normalized_gradient_times_side_squared': str(gradient),
        'normalized_gradient_times_side_squared_upper': str(F(7, 30)*dual_energy),
        'exact_normalized_interpolation_mass_error': str(mass),
        'normalized_interpolation_mass_error_upper': str(F(5, 18*len(faces))*dual_energy)}


def tent_geometry_control():
    """Exact local honeycomb tensor and continuous barycentric tent moments."""
    root = sp.sqrt(3)
    basis = sp.Matrix([[1, sp.Rational(1, 2)], [0, root/2]])
    offsets = [basis*sp.Matrix(v)/3 for v in ((1, 1), (-2, 1), (1, -2))]
    tensor = sum((v*v.T for v in offsets), sp.zeros(2))
    if sp.simplify(tensor-sp.eye(2)/2) != sp.zeros(2):
        raise ArithmeticError('the exact honeycomb centroid energy tensor changed')
    # One of the three tent regions has corners (0,0),(1,0),(1/3,1/3).
    # The affine simplex coordinates have determinant1/3 and tent value v.
    u, v = sp.symbols('u v', nonnegative=True)
    moments = [sp.integrate(sp.integrate(v**power/3, (v, 0, 1-u)), (u, 0, 1))*6
               for power in (1, 2)]
    if moments != [sp.Rational(1, 3), sp.Rational(1, 6)]:
        raise ArithmeticError('the integrated barycentric tent moments changed')
    gradient_squared = sp.simplify((3/(root/2))**2)
    if gradient_squared != 12:
        raise ArithmeticError('the equilateral tent gradient changed')
    quotient = gradient_squared/(48*moments[1])
    if quotient != sp.Rational(3, 2):
        raise ArithmeticError('the continuum tent energy quotient changed')
    return {'unit_side_centroid_tensor': 'I_2/2',
        'normalized_tent_mean': str(moments[0]),
        'normalized_tent_second_moment': str(moments[1]),
        'unit_side_tent_gradient_squared': str(gradient_squared),
        'native_limit_tent_energy_quotient_at_coarse_level_zero': str(quotient)}


def tent_refinement_control(level):
    """Exact source-mesh energies for all twenty continuous ancestral tents."""
    level = _level(level)
    if level > 4:
        raise ValueError('finite tent refinement controls use levels zero through four')
    m = 2**level
    vertices = ((0, 0, 0, 0, 0, 0), (2, 0, 0, 0, 0, 0), (0, 2, 0, 0, 0, 0))
    faces = ((0, 1, 2),)
    for _ in range(level):
        vertices, faces = uniform.refine(vertices, faces)
    numerators = []
    for face in faces:
        x = sum(vertices[p][0]//2 for p in face)
        y = sum(vertices[p][1]//2 for p in face)
        numerators.append(min(x, y, 3*m-x-y))
    if any(not 0 < value <= m for value in numerators):
        raise ArithmeticError('a centroid is outside its ancestral tent')
    _, actual_faces = uniform.affine_prefix(level)
    edges = uniform.dual_edges(actual_faces)
    count = 20*m*m
    # Every score is at most m, and there are30*m² edges. This bound
    # certifies all integer products and partial sums before they execute.
    if 30*m**4 >= 2**40:
        raise ArithmeticError('the finite exact tent control exceeds its arithmetic budget')
    scores = np.zeros((count, 20), dtype=np.int64)
    for base in range(20):
        scores[base*m*m:(base+1)*m*m, base] = numerators
    energy = np.zeros((20, 20), dtype=np.int64)
    for a, b in edges:
        difference = scores[a]-scores[b]
        energy += np.outer(difference, difference)
    gram = scores.T@scores
    if (not np.array_equal(gram, np.eye(20, dtype=np.int64)*sum(x*x for x in numerators))
            or np.any(energy.sum(axis=1) < 0)):
        raise ArithmeticError('continuous tent norm or source energy control failed')
    diagonal = tuple(F(int(energy[i, i]), 20*m*m) for i in range(20))
    off_diagonal = max(abs(F(int(energy[i, j]), 20*m*m))
                       for i in range(20) for j in range(20) if i != j)
    return {'level': level, 'source_cells': count,
        'exact_normalized_patch_tent_mean': str(F(sum(numerators), m**3)),
        'exact_normalized_patch_tent_second_moment': str(F(sum(x*x for x in numerators), m**4)),
        'max_dual_energy_diagonal_error_from_three_twentieths': str(max(abs(x-F(3, 20)) for x in diagonal)),
        'max_dual_energy_off_diagonal': str(off_diagonal),
        'integer_source_adjacency_and_all_twenty_tent_energies_checked': True}


def full_port_form_control():
    """Exact harmonic lift, normal contraction and cutoff on the native source."""
    graph, _, r, integer_f, denominator = source.matrices(0)
    _, carrier = source.passage.recovered_carrier()
    local_green = (carrier+sp.ones(12)/12).inv()-sp.ones(12)/12
    ri, fi = sp.Matrix(r), sp.Matrix(integer_f)
    values = sp.Matrix([sp.Rational(i % 7-3, 3) for i in range(graph.cells)])
    edge_response = fi.inv()*ri*values*denominator
    drive = sp.zeros(graph.ports, 1)
    for j, (a, b) in enumerate(graph.inter_seams):
        drive[a] += edge_response[j]; drive[b] -= edge_response[j]
    internal = sp.zeros(graph.ports, 1)
    for cell in range(graph.cells):
        internal[12*cell:12*(cell+1), 0] = -local_green*drive[12*cell:12*(cell+1), 0]
    lifted = sp.Matrix([values[p//12]+internal[p] for p in range(graph.ports)])
    def energy(vector):
        return sum((vector[a]-vector[b])**2 for a, b in graph.seams)
    raw_energy = energy(lifted)
    raw_internal = sum(x*x for x in internal)
    if (any(sum(internal[12*c:12*(c+1), 0]) != 0 for c in range(graph.cells))
            or raw_energy != (values.T*ri.T*edge_response)[0]
            or 2*raw_internal > raw_energy):
        raise ArithmeticError('the full-port harmonic energy lift failed')
    contracted = sp.Matrix([max(sp.S.Zero, min(sp.S.One, x)) for x in lifted])
    coarse_contracted = sp.Matrix([max(sp.S.Zero, min(sp.S.One, values[p//12]))
                                   for p in range(graph.ports)])
    contracted_energy = energy(contracted)
    contraction_error = sum((x-y)**2 for x, y in zip(contracted, coarse_contracted, strict=True))
    if contracted_energy > raw_energy or contraction_error > raw_internal:
        raise ArithmeticError('the native full-port normal contraction bound failed')
    # A separate exact field checks the product inequality without assuming
    # that a cell-mean Schur matrix is a Markov generator.
    test = sp.Matrix([sp.Rational(7*p % 19-9, 5) for p in range(graph.ports)])
    cutoff = [sp.Rational((p//12) % 5, 4) for p in range(graph.ports)]
    altered = sp.Matrix([x*y for x, y in zip(test, cutoff, strict=True)])
    transition = {p for a, b in graph.inter_seams if cutoff[a] != cutoff[b] for p in (a, b)}
    omega = max((cutoff[a]-cutoff[b])**2 for a, b in graph.inter_seams)
    epsilon = sp.Rational(1, 3)
    allowance = ((1+epsilon)*energy(test)
                 +(1+1/epsilon)*omega*sum(test[p]**2 for p in transition))
    if energy(altered) > allowance:
        raise ArithmeticError('the exact native cutoff product bound failed')
    return {'level': 0, 'ports': graph.ports,
        'full_port_harmonic_energy_identity_checked': True,
        'internal_norm_squared_over_raw_energy': str(raw_internal/raw_energy),
        'contracted_energy_over_raw_energy': str(contracted_energy/raw_energy),
        'contraction_error_squared_over_internal_norm_squared': str(contraction_error/raw_internal),
        'cutoff_energy': str(energy(altered)), 'cutoff_energy_upper': str(allowance),
        'constant_zero_mode_adjoined_before_normal_contraction': True}


def separated_action_bound(level, cell_distance):
    """Rescaled Schur bilinear bound for unit vectors at cell distance d>=1."""
    level = _level(level)
    if type(cell_distance) is not int or cell_distance < 1:
        raise ValueError('a positive integer cell distance is required')
    return F(4**level, 2**cell_distance)


def macroscopic_action_exponent(level, separation):
    """Base-two exponent of the separated-support action upper bound.

    Separation is intrinsic distance on the auxiliary unit-circumradius flat
    polyhedron.  Adjacent centroids are less than 2^-l apart; each fine cell
    has diameter less than 2*2^-l.  This is a comparison length, not metres.
    An exponent avoids allocating a doubly exponential denominator.
    """
    level = _level(level)
    if type(separation) is not F or separation <= 0:
        raise ValueError('a positive exact Fraction separation is required')
    lower = separation*2**level-4
    distance = max(1, -((-lower.numerator)//lower.denominator))
    if lower <= 0:
        raise ValueError('the cells are not yet certified disjoint at this level')
    return 2*level-distance


def triangle_patch(depth):
    """Exact upright/downward cells and dual edges in one descendant triangle."""
    depth = _level(depth)
    if depth > 5:
        raise ValueError('the finite stencil control is limited to depth five')
    m = 2**depth
    up = tuple((i, j) for i in range(m) for j in range(m-i))
    down = tuple((i, j) for i in range(m-1) for j in range(m-1-i))
    vertices = tuple(('u', *v) for v in up)+tuple(('d', *v) for v in down)
    index = {v: i for i, v in enumerate(vertices)}
    edges = tuple(sorted(tuple(sorted((index[('d', i, j)], index[('u', x, y)])))
        for i, j in down for x, y in ((i, j), (i+1, j), (i, j+1))))
    return m, vertices, edges


def reflected_square_control(depth):
    """The square has exactly twice each upright-triangle edge energy."""
    m, vertices, edges = triangle_patch(depth)
    def reflect(i, j):
        return (i, j) if i+j < m else (m-1-j, m-1-i)
    upright_edges = Counter()
    for i in range(m):
        for j in range(m-i):
            for x, y in ((i+1, j), (i, j+1)):
                if x+y < m:
                    upright_edges[tuple(sorted(((i, j), (x, y))))] = 2
    reflected_edges = Counter()
    for i in range(m):
        for j in range(m):
            for x, y in ((i+1, j), (i, j+1)):
                if x < m and y < m:
                    reflected_edges[tuple(sorted((reflect(i, j), reflect(x, y))))] += 1
    if reflected_edges != upright_edges:
        raise ArithmeticError('square reflection changed the exact edge multiplicities')
    multiplicity = Counter(reflect(i, j) for i in range(m) for j in range(m))
    if (set(multiplicity) != {(i, j) for kind, i, j in vertices if kind == 'u'}
            or not all(n in (1, 2) for n in multiplicity.values())):
        raise ArithmeticError('upright cells are not covered once or twice')
    return {'depth': depth, 'side_subdivisions': m, 'descendant_cells': len(vertices),
        'internal_dual_edges': len(edges), 'square_edges': sum(reflected_edges.values()),
        'square_edge_energy_is_twice_upright_energy': True,
        'all_upright_cells_have_multiplicity_one_or_two': True}


def source_patch_control(depth):
    """Bind every base-face descendant adjacency to the independent stencil."""
    if _level(depth) > 3:
        raise ValueError('native ancestry binding controls use depths zero through three')
    m, vertices, edges = triangle_patch(depth)
    # The first two slots here are independent combinatorial triangle
    # coordinates; they are not the physical x/y slots of the carrier chart.
    simple = ((0, 0, 0, 0, 0, 0), (2, 0, 0, 0, 0, 0), (0, 2, 0, 0, 0, 0))
    faces = ((0, 1, 2),)
    for _ in range(depth):
        simple, faces = uniform.refine(simple, faces)
    lookup = {v: i for i, v in enumerate(vertices)}
    labels = []
    for face in faces:
        corners = tuple((simple[v][0]//2, simple[v][1]//2) for v in face)
        i, j = min(x for x, _ in corners), min(y for _, y in corners)
        kind = 'u' if (i, j) in corners else 'd'
        labels.append(lookup[(kind, i, j)])
    _, actual_faces = uniform.affine_prefix(depth)
    actual_edges = uniform.dual_edges(actual_faces)
    expected = set(edges)
    for base in range(20):
        lower, upper = base*m*m, (base+1)*m*m
        internal = {tuple(sorted((labels[a-lower], labels[b-lower])))
                    for a, b in actual_edges if lower <= a < upper and lower <= b < upper}
        if internal != expected:
            raise ArithmeticError('native descendants differ from the triangular patch stencil')
    return {'depth': depth, 'base_patches_checked': 20,
        'cell_ancestry_order_and_all_internal_edges_checked': True}


def exact_patch_energy_control(depth):
    """Independent rational LDL certificate for the full patch variance bound."""
    m, vertices, edges = triangle_patch(depth)
    if depth > 3:
        raise ValueError('exact matrix energy controls are limited to depth three')
    n = len(vertices)
    laplacian = sp.zeros(n)
    for a, b in edges:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    constant = sp.ones(n)/n
    certificate = 26*m*m*laplacian-sp.eye(n)+2*constant
    # Constant direction has eigenvalue one; the other directions certify
    # 26*m² L_patch >= I-P_constants without a floating eigenvalue test.
    _, pivots = certificate.LDLdecomposition(hermitian=False)
    if any(pivots[i, i] <= 0 for i in range(n)):
        raise ArithmeticError('the exact triangular patch Poincare bound failed')
    return {'depth': depth, 'cells': n,
        'variance_energy_multiplier': 26*m*m,
        'rational_LDL_positive_pivots': n}


def native_operator_control(level):
    """Native support, compactness and separated-action finite controls."""
    if _level(level) > 2:
        raise ValueError('native operator controls use levels zero through two')
    graph, _, r, integer_f, denominator = source.matrices(level)
    cell_l = r.T@r
    adjacency = -(cell_l-np.diag(np.diag(cell_l)))
    distances = shortest_path(adjacency, directed=False, unweighted=True)
    numerator_t = 2*denominator*np.eye(len(integer_f), dtype=np.int64)-integer_f
    term = r.copy()
    support_rows = []
    row_norm = max(sum(abs(int(value)) for value in row) for row in numerator_t)
    incidence_column_norm = max(sum(abs(int(value)) for value in column) for column in r.T)
    integer_limit = int(np.iinfo(np.int64).max)
    for power in range(4):
        term_norm = max(abs(int(value)) for value in term.flat)
        if incidence_column_norm*term_norm > integer_limit:
            raise ArithmeticError('exact polynomial coefficient could overflow int64')
        coefficient = r.T@term
        if np.any(coefficient[distances > power+1]):
            raise ArithmeticError('Schur polynomial has a forbidden long-range coefficient')
        support_rows.append({'power': power, 'maximum_cell_range': power+1,
            'integer_polynomial_support_checked': True})
        if power < 3:
            if row_norm*term_norm > integer_limit:
                raise ArithmeticError('exact polynomial power could overflow int64')
            term = numerator_t@term
    schur = r.T@np.linalg.solve(integer_f/denominator, r)/12
    green = source._green(schur)/4**level
    compact_rows = []
    for coarse in range(level+1):
        descendants = 4**(level-coarse)
        projected = np.repeat(green.reshape(20*4**coarse, descendants, graph.cells).mean(axis=1),
                              descendants, axis=0)
        squared = float(np.linalg.norm(green-projected, 2)**2)
        if squared > float(compactness_tail_squared(coarse))+1e-9:
            raise ArithmeticError('native Green violates the coarse compactness bound')
        compact_rows.append({'coarse_level': coarse,
            'tail_norm_squared_diagnostic': squared,
            'tail_norm_squared_upper': str(compactness_tail_squared(coarse))})
    distance = int(distances.max())
    a, b = map(int, np.argwhere(distances == distance)[0])
    coupling = float(abs(4**level*schur[a, b]))
    bound = separated_action_bound(level, distance)
    if coupling > float(bound)+1e-12:
        raise ArithmeticError('native distant action exceeded its exact bound')
    return {'level': level, 'cells': graph.cells,
        'exact_polynomial_support': support_rows, 'compactness_controls': compact_rows,
        'distant_native_pair': [a, b], 'cell_distance': distance,
        'absolute_rescaled_coupling_diagnostic': coupling,
        'absolute_rescaled_coupling_upper': str(bound)}


def non_markov_control():
    """Averages are not a vertex Kron reduction: a positive off-diagonal survives."""
    _, _, r, integer_f, denominator = source.matrices(0)
    ri, fi = sp.Matrix(r), sp.Matrix(integer_f)
    entry = (ri[:, 5].T*fi.inv()*ri[:, 9])[0]*denominator/12
    if entry <= 0:
        raise ArithmeticError('the actual source mean-elimination countercontrol changed')
    return {'level': 0, 'cells': [5, 9], 'exact_positive_Schur_off_diagonal': str(entry),
        'finite_cell_mean_generator_is_a_Markov_generator': False,
        'continuum_locality_or_nonlocality_decided_by_this_control': False}


def run():
    return {'schema': 'oph.native_compact_spatial_generator.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'contract': {'source': 'existing spherical shared-midpoint carrier joins with distinct ports',
            'space': 'mean-zero L2 of normalized auxiliary flat-polyhedron area',
            'observation_dictionary': 'normalized ancestral cell means',
            'limit': 'the same common spatial subsequence as the native Green field',
            'physical_spatial_metric_or_radial_coordinate_supplied_by_this_result': False},
        'exact_patch_controls': [reflected_square_control(k) | source_patch_control(k)
                                 | exact_patch_energy_control(k) for k in range(3)],
        'native_operator_controls': [native_operator_control(k) for k in range(3)],
        'mean_elimination_control': non_markov_control(),
        'full_port_form_control': full_port_form_control(),
        'exact_interpolation_geometry': interpolation_geometry_control(),
        'exact_source_interpolation_controls': [source_interpolation_control(k) for k in range(4)],
        'exact_tent_geometry': tent_geometry_control(),
        'exact_tent_refinement_controls': [tent_refinement_control(k) for k in range(5)],
        'all_level_bounds': {
            'patch_variance_energy_multiplier': '26*4^(l-k)',
            'coarse_projection_energy_multiplier': '624*4^(-k)',
            'Green_norm_upper': str(GREEN_BOUND),
            'Green_tail_squared': {str(k): str(compactness_tail_squared(k)) for k in (0, 4, 20)},
            'covariance_eigenvalue_number_20_times_4_power_k_upper': {
                str(k): str(covariance_eigenvalue_tail(k)) for k in (0, 4, 20)},
            'separated_cell_support_bilinear_action': '4^l*2^(-d)*||u||*||v||, d>=1',
            'macroscopic_separation_one_quarter_action_base_two_exponents': {
                str(k): macroscopic_action_exponent(k, F(1, 4)) for k in (5, 8, 12)},
            'full_port_internal_tail_per_unit_energy': {
                str(k): str(harmonic_lift_tail(k)) for k in (0, 4, 20)},
            'interpolation_gradient_squared_upper': '112/a^2*E_l(f_l)',
            'interpolation_mass_error_upper': '20/(3*4^l)*E_l(f_l)',
            'continuum_form_comparison_on_H1': 'a^2/112*integral|grad f|^2 <= E(f) <= a^2/48*integral|grad f|^2',
            'continuum_positive_eigenvalue_enclosures': {
                str(n): {key: str(value) for key, value in continuum_eigenvalue_enclosure(n).items()}
                for n in (1, 19, 20, 79, 80, 1000000)},
            'mean_zero_heat_trace_lower': '1/(exp(3*t/2)-1), t>0',
            'mean_zero_heat_trace_upper': '1/(exp(t/49920)-1), t>0'},
        'theorem': {'common_Green_subsequence_converges_in_operator_norm': True,
            'limit_is_compact_positive_and_injective_on_mean_zero_L2': True,
            'inverse_is_positive_self_adjoint_spatial_generator': True,
            'finite_reduced_forms_Mosco_converge_on_full_L2_with_constant_zero_mode': True,
            'limit_is_conservative_strongly_local_Dirichlet_form': True,
            'form_domain_is_H1_of_the_auxiliary_flat_surface': True,
            'continuous_piecewise_affine_functions_are_a_regular_form_core': True,
            'heat_trace_spectral_dimension': 2,
            'same_dictionary_scaled_compact_smooth_second_order_three_dimensional_elliptic_identification_excluded': True,
            'finite_rescaled_action_couplings_at_separated_continuum_supports_vanish': True,
            'unique_whole_sequence_generator_selected': False,
            'exact_intrinsic_metric_identified': False,
            'limit_identified_as_round_sphere_Laplacian': False,
            'three_dimensional_space_or_Lorentzian_metric_derived': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native compact spatial-generator bounds and source controls verified.')
    else:
        print(canonical(result))


if __name__ == '__main__':
    main()
