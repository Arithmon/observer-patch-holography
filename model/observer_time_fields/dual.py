"""Exact directed backward evolution of squarefree native field monomials.

For a binary source configuration, Q_e n_X = (n_X+n_{swap_e X})/2.
The one- and two-label exclusion lifts are therefore exact observable maps,
not an independent-particle approximation. Private sparse integer matrices
implement the stopped backward map sum_j A**j B, with downward division.
"""
from functools import lru_cache
from collections import Counter
from itertools import combinations, product
from typing import NamedTuple

import numpy as np
from scipy import sparse
from ..observer import native

SOURCE_PATHS = ('model/observer/native.py', 'model/observer_time_fields/dual.py')
MAX_PORTS = 256
MAX_EDGES = 4096
MAX_BITS = 52


class Operator(NamedTuple):
    ports: tuple
    edges: tuple
    incident: tuple
    degree: int


def build(ports, edges, incident, degree=1):
    """Declare a finite graph using original, canonical immutable integers.

    Only this immutable declaration escapes. Sparse arrays and mutable cache
    envelopes stay private; callers cannot poison subsequent evolutions.
    """
    if (type(ports) is not tuple or not 2 <= len(ports) <= MAX_PORTS
            or any(type(p) is not int or p < 0 for p in ports)):
        raise ValueError('two through 256 original nonnegative integer ports required')
    if ports != tuple(sorted(set(ports))):
        raise ValueError('distinct increasing port labels required')
    if type(degree) is not int or degree not in (1, 2):
        raise ValueError('an exact monomial degree one or two required')
    domain = set(ports)
    for menu in (edges, incident):
        if type(menu) is not tuple or not 1 <= len(menu) <= MAX_EDGES:
            raise ValueError('a nonempty bounded original seam tuple required')
        for edge in menu:
            if (type(edge) is not tuple or len(edge) != 2
                    or any(type(p) is not int for p in edge)
                    or edge[0] >= edge[1] or not set(edge) <= domain):
                raise ValueError('original ordered endpoint pairs inside the graph required')
        if menu != tuple(sorted(set(menu))):
            raise ValueError('distinct increasing seams required')
    if not set(incident) <= set(edges):
        raise ValueError('every clock seam must belong to the complete graph')
    return Operator(ports, edges, incident, degree)


def source_carrier_certificate(ports, edges, incident):
    """Check the actual root carrier giving the strict incident bound 15/M.

    On three binary ports the original native count C=2*sum_e Q_e has
    spectrum in {3,6}, hence each triangle contributes at least 3/2 I.
    The twenty source triangles cover each of thirty internal seams twice,
    giving sum_internal Q_e >= 15 I. Tensoring with arbitrary other binary
    ports preserves the inequality; further incident Q_e are PSD.
    """
    operator = build(ports, edges, incident)
    root = set(range(12))
    internal = {edge for edge in operator.edges if set(edge) <= root}
    faces = tuple(face for face in combinations(range(12), 3)
                  if set(combinations(face, 2)) <= internal)
    coverage = Counter(edge for face in faces for edge in combinations(face, 2))
    if (not root <= set(operator.ports) or len(internal) != 30 or len(faces) != 20
            or not internal <= set(operator.incident)):
        raise ValueError('the complete incident root carrier needs thirty seams and twenty faces')
    if set(coverage) != internal or set(coverage.values()) != {2}:
        raise ValueError('every root seam must occur in exactly two source faces')
    binary = tuple(product((0, 1), repeat=3))
    index = {state: i for i, state in enumerate(binary)}
    counts = np.zeros((8, 8), dtype=np.int64)
    for column, state in enumerate(binary):
        for seam in combinations(range(3), 2):
            outputs = []
            for coin in (0, 1):
                final, _ = native.native(state, seam, coin)
                outputs.append(final)
                counts[index[final], column] += 1
            swapped = list(state)
            a, b = seam
            swapped[a], swapped[b] = swapped[b], swapped[a]
            if sorted(outputs) != sorted((state, tuple(swapped))):
                raise ArithmeticError('the original binary native map is not identity plus swap')
    identity = np.eye(8, dtype=np.int64)
    if (not np.array_equal(counts, counts.T) or np.any(counts.sum(axis=0) != 6)
            or np.any((counts-3*identity)@(counts-6*identity))):
        raise ArithmeticError('the original binary triangle lost its exact spectral certificate')
    return {'root_ports': 12, 'internal_seams': 30, 'source_triangles': 20,
            'triangle_edge_multiplicity': 2, 'binary_triangle_states': 8,
            'triangle_count_minimum_eigenvalue': 3,
            'triangle_count_maximum_eigenvalue': 6,
            'unnormalized_incident_transfer_lower': 15}


def _checked(operator):
    if type(operator) is not Operator:
        raise ValueError('a typed immutable native dual operator required')
    return build(*operator)


@lru_cache(maxsize=16)
def _states(operator):
    return tuple(combinations(operator.ports, operator.degree))


def states(operator):
    """Canonical monomial labels, with no repeated occupied port."""
    return _states(_checked(operator))


@lru_cache(maxsize=8)
def _matrices(operator):
    labels = _states(operator)
    index = {label: i for i, label in enumerate(labels)}
    incident = set(operator.incident)
    neighbours = {p: [] for p in operator.ports}
    for a, b in operator.edges:
        flag = int((a, b) in incident)
        neighbours[a].append((b, flag)); neighbours[b].append((a, flag))
    rows, columns = [[], []], [[], []]
    degrees = [np.zeros(len(labels), dtype=np.int64) for _ in range(2)]
    for column, occupied in enumerate(labels):
        for old in occupied:
            for new, flag in neighbours[old]:
                if new in occupied:
                    continue
                target = tuple(sorted(new if p == old else p for p in occupied))
                rows[flag].append(index[target]); columns[flag].append(column)
                degrees[flag][column] += 1
    counts = (len(operator.edges)-len(incident), len(incident))
    matrices = []
    for flag, count in enumerate(counts):
        values = [1]*len(rows[flag])
        diagonal = 2*count-degrees[flag]
        for i, value in enumerate(diagonal):
            if value:
                rows[flag].append(i); columns[flag].append(i); values.append(int(value))
        matrix = sparse.coo_matrix((values, (rows[flag], columns[flag])),
            shape=(len(labels), len(labels)), dtype=np.int64).tocsr()
        if (np.any(matrix.data < 0) or (matrix-matrix.T).nnz
                or np.any(np.asarray(matrix.sum(axis=0)).ravel() != 2*count)
                or np.any(np.asarray(matrix.sum(axis=1)).ravel() != 2*count)):
            raise ArithmeticError('the exact exclusion lift lost its symmetric constant mass')
        for array in (matrix.data, matrix.indices, matrix.indptr):
            array.flags.writeable = False
        matrices.append(matrix)
    return tuple(matrices)


def evolve(operator, initial_integer_vector, ticks, cap=512, bits=52):
    """Return lower coefficients of the completed backward field polynomial.

    Input and output use the common scale 2**bits. The nonnegative input may
    have partial mass. Each interval retains cap terms A**j B, j=0..cap-1;
    no unresolved mass is normalized back into the answer. If p is the clock
    seam fraction, the ideal capped coefficient mass is input_mass times
    (1-(1-p)**cap)**ticks. Componentwise domination makes the difference
    between that mass and the returned mass the exact finite-cap L1 error.

    A row's integer coefficient sum is at most 2M. Thus 2M*2**bits<2**63
    bounds every multiplication and accumulation; all partial completed
    vectors have mass at most the input mass. No floating operation occurs.
    """
    operator = _checked(operator)
    if (type(ticks) is not int or not 0 <= ticks <= 16
            or type(cap) is not int or not 1 <= cap <= 2048
            or type(bits) is not int or not 8 <= bits <= MAX_BITS):
        raise ValueError('bounded exact ticks, interval cap and grid precision required')
    denominator, grid = 2*len(operator.edges), 1 << bits
    if denominator*grid >= 1 << 63:
        raise ValueError('the declared grid exceeds the proved int64 product bound')
    if (type(initial_integer_vector) is not tuple
            or len(initial_integer_vector) != len(_states(operator))
            or any(type(value) is not int or value < 0 for value in initial_integer_vector)
            or sum(initial_integer_vector) > grid):
        raise ValueError('an original nonnegative integer tuple of mass at most the grid required')
    if ticks == 0:
        return initial_integer_vector
    remote, incident = _matrices(operator)
    vector = np.array(initial_integer_vector, dtype=np.int64)
    for _ in range(ticks):
        maximum_mass = int(vector.sum())
        waiting = (incident@vector)//denominator
        completed = waiting.copy()
        for _ in range(1, cap):
            waiting = (remote@waiting)//denominator
            completed += waiting
        if np.any(completed < 0) or int(completed.sum()) > maximum_mass:
            raise ArithmeticError('directed backward evolution exceeded its integer mass budget')
        vector = completed
    return tuple(map(int, vector))
