"""Exact native quantum instruments on a source incident-attempt clock.

The one-particle carrier retains populations and every complex coherence.
For cap h, the symmetric successful maps are
sum_{l+r<h} c_l c_r A**l B_m A**r, c_n=binom(2n,n)/4**n.
Together with A**h they form a complete rational CP instrument. The native
chronology instead uses sum_{n<h} B_m A**n, with the same timeout.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import math
from typing import NamedTuple

import numpy as np

from ..observer import native
from ..observer_time_fields import core
from ..record_laws import family

REGIONS = ((0, 1), (0, 5))
UNION = (0, 1, 5)
PAIR = (0, 1)
MARKS = ((0, 1, 0), (1, 0, 0), (1, 1, -1), (1, 1, 0), (1, 1, 1))
TIMEOUT = 'timeout'
SOURCE_PATHS = tuple(sorted(set((*core.SOURCE_PATHS,
    'model/observer_correspondence/quantum_clock.py',
    'model/record_laws/family.py', 'model/predictions/observer_instrument.py'))))


class Channel(NamedTuple):
    population: tuple
    coherence: tuple


class Source(NamedTuple):
    candidate: str
    root: tuple
    seams: tuple
    remote: Channel
    incident: tuple
    probability: F
    incident_hs_lower: F


def _parameters(law, root, cap=None, effective=False):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('the original C or E native instrument is required')
    if (type(root) is not tuple or any(type(p) is not int for p in root)
            or root not in (UNION, PAIR)):
        raise ValueError('the canonical union or pair clock region is required')
    if cap is not None and (type(cap) is not int or not 1 <= cap <= 32):
        raise ValueError('an original integer interval cap from one through thirty-two is required')
    if type(effective) is not bool:
        raise ValueError('an original Boolean representation selector is required')


def _integers(matrix):
    return tuple(tuple(int(value) for value in row) for row in matrix)


def _fresh(population, coherence, denominator):
    # Fractions have writable protected slots. Cached certificate state holds
    # integers only; even those slots in a public result cannot poison it.
    return Channel(*(tuple(tuple(F(value, denominator) for value in row) for row in matrix)
                     for matrix in (population, coherence)))


def _mark(seam, current):
    return (*(int(bool(set(seam).intersection(region))) for region in REGIONS),
            current*(int(seam[1] == 0)-int(seam[0] == 0)))


@lru_cache(maxsize=4)
def _source_counts(law, root):
    seams = core.source().seams
    states = tuple(tuple(int(i == p) for i in range(12)) for p in range(12))
    channels = {key: (np.zeros((12, 12), dtype=object), np.zeros((12, 12), dtype=object))
                for key in (None, *MARKS)}
    for seam in seams:
        incident = bool(set(seam).intersection(root))
        for coin in (0, 1):
            images = []
            for state in states:
                output, record = native.native(state, seam, coin)
                if (record != (coin, state[seam[0]]-output[seam[0]])
                        or sum(output) != 1 or any(x not in (0, 1) for x in output)):
                    raise ArithmeticError('the original binary native branch changed')
                images.append((output.index(1), record[1],
                    record[1] if law == 'C' else state[seam[0]]-state[seam[1]]))
            for left, right in product(range(12), repeat=2):
                a, b = images[left], images[right]
                if a[2] != b[2]:
                    continue
                if a[1] != b[1]:
                    raise ArithmeticError('a candidate branch does not have a definite current')
                key = _mark(seam, a[1]) if incident else None
                population, coherence = channels[key]
                if left == right:
                    population[a[0], left] += 1
                else:
                    if (a[0], b[0]) != (left, right):
                        raise ArithmeticError('one-particle off-diagonal action is not scalar')
                    coherence[left, right] += 1
    count = sum(bool(set(edge).intersection(root)) for edge in seams)
    unread = sum((channels[key][0] for key in MARKS), np.zeros((12, 12), dtype=object))
    remote = channels[None][0]
    if (len(seams) != 30 or count != (12 if root == UNION else 9)
            or not np.array_equal(unread, unread.T) or not np.array_equal(remote, remote.T)
            or set(unread.sum(axis=0)) != {2*count}
            or set(remote.sum(axis=0)) != {60-2*count}):
        raise ArithmeticError('the canonical source clock or symmetric transfer changed')
    if (min(remote[i, i]-sum(remote[i, j] for j in range(12) if i != j)
            for i in range(12)) < 0
            or any(not 0 <= channels[None][1][i, j] <= 60-2*count
                   for i, j in product(range(12), repeat=2))
            or any(not 0 <= sum(channels[key][1][i, j] for key in (None, *MARKS)) <= 60
                   for i, j in product(range(12), repeat=2))):
        raise ArithmeticError('the remote CP transfer lost its Hilbert-Schmidt contraction bounds')
    lower = min(int(unread[i, i]-sum(unread[i, j] for j in range(12) if i != j))
                for i in range(12))
    lower = min(lower, *(sum(channels[key][1][i, j] for key in MARKS)
                        for i, j in product(range(12), repeat=2) if i != j))
    expected = {('C', UNION): 14, ('E', UNION): 6, ('C', PAIR): 8, ('E', PAIR): 0}
    if lower != expected[law, root]:
        raise ArithmeticError('the exact one-particle incident spectral floor changed')
    return seams, count, int(lower), tuple((key, _integers(p), _integers(c))
                                          for key, (p, c) in channels.items())


def source(law, root=UNION):
    """Original native source, with detached exact matrices and no cache aliases."""
    _parameters(law, root)
    seams, count, lower, rows = _source_counts(law, root)
    maps = {key: _fresh(p, c, 60) for key, p, c in rows}
    return Source(law, root, seams, maps[None], tuple((key, maps[key]) for key in MARKS),
                  F(count, 30), F(lower, 60))


@lru_cache(maxsize=16)
def _capped_counts(law, root, cap, effective):
    _, _, _, rows = _source_counts(law, root)
    counts = {key: (np.array(p, dtype=object), np.array(c, dtype=object)) for key, p, c in rows}
    a, ac = counts[None]
    powers = [np.eye(12, dtype=object)]
    for _ in range(cap):
        powers.append(a@powers[-1])
    series = sum((60**(cap-1-n)*powers[n] for n in range(cap)),
                 np.zeros((12, 12), dtype=object))
    central = tuple(math.comb(2*n, n) for n in range(cap))
    if any(sum(central[l]*central[n-l] for l in range(n+1)) != 4**n for n in range(cap)):
        raise ArithmeticError('the exact half-wait convolution identity failed')
    denominator = 60*240**(cap-1) if effective else 60**cap
    result = []
    for mark in MARKS:
        b, bc = counts[mark]
        if effective:
            population = np.zeros((12, 12), dtype=object)
            for left in range(cap):
                product_left = powers[left]@b
                for right in range(cap-left):
                    population += (central[left]*central[right]*240**(cap-1-left-right))*(
                        product_left@powers[right])
        else:
            population = b@series
        # All coherence blocks are scalars, so the same central-binomial
        # convolution reduces both chronologies to their geometric sum.
        coherence = np.zeros((12, 12), dtype=object)
        for i, j in product(range(12), repeat=2):
            if i != j:
                coherence[i, j] = bc[i, j]*sum(
                    int(ac[i, j])**n*60**(cap-1-n) for n in range(cap))
        result.append((mark, _integers(population), _integers(coherence), denominator, 60**cap))
    timeout_coherence = np.zeros((12, 12), dtype=object)
    for i, j in product(range(12), repeat=2):
        if i != j:
            timeout_coherence[i, j] = int(ac[i, j])**cap
    result.append((TIMEOUT, _integers(powers[cap]), _integers(timeout_coherence), 60**cap, 60**cap))
    return tuple(result)


def capped_maps(law, cap, effective=False, root=UNION):
    """Every finite CP outcome, including zero maps and the unrepaired timeout.

    Native and symmetric maps have the same state-independent timeout q**h.
    No successful output is conditioned or renormalized. One interval uses
    at most h native operations, including the incident operation if any.
    """
    if type(cap) is not int:
        raise ValueError('an original integer interval cap is required')
    _parameters(law, root, cap, effective)
    return {mark: Channel(
        tuple(tuple(F(value, pd) for value in row) for row in population),
        tuple(tuple(F(value, cd) for value in row) for row in coherence))
        for mark, population, coherence, pd, cd in _capped_counts(law, root, cap, effective)}


def _exact_matrix(value, size):
    if type(value) is np.ndarray:
        if value.shape != (size, size) or value.dtype != object:
            raise ValueError('an original square object array of exact scalars is required')
        rows = value.tolist()
    elif type(value) is tuple and len(value) == size:
        if any(type(row) is not tuple or len(row) != size for row in value):
            raise ValueError('original square exact matrix rows are required')
        rows = value
    else:
        raise ValueError('an original exact matrix is required')
    if any(type(x) not in (int, F) for row in rows for x in row):
        raise ValueError('original integer or rational matrix entries are required')
    return np.array([[F(x.numerator, x.denominator) if type(x) is F else F(x)
                      for x in row] for row in rows], dtype=object)


def apply(channel, density, reference_dimension=1):
    """Exact real-linear action; apply separately to real/imaginary operator parts.

    The reference is the minor tensor factor. This is an operator helper,
    not a positivity test or a certificate that a supplied map is native.
    """
    if type(channel) is not Channel or type(reference_dimension) is not int or not 1 <= reference_dimension <= 4:
        raise ValueError('a typed channel and reference dimension one through four are required')
    population = _exact_matrix(channel.population, 12)
    coherence = _exact_matrix(channel.coherence, 12)
    if any(coherence[i, i] for i in range(12)):
        raise ValueError('coherence coefficients have a zero diagonal by convention')
    reference = reference_dimension
    original = _exact_matrix(density, 12*reference).reshape(12, reference, 12, reference)
    output = np.empty_like(original)
    output.fill(F(0))
    for i, j in product(range(12), repeat=2):
        if i == j:
            for k in range(12):
                output[i, :, i, :] += population[i, k]*original[k, :, k, :]
        else:
            output[i, :, j, :] = coherence[i, j]*original[i, :, j, :]
    return output.reshape(12*reference, 12*reference)


@lru_cache(maxsize=1)
def _total_two_data():
    states, seams = core.source().states, core.source().seams
    index = {state: i for i, state in enumerate(states)}
    sectors = {}
    for i, j in product(range(66), repeat=2):
        sectors.setdefault(tuple(a-b for a, b in zip(states[i], states[j])), []).append((i, j))
    images = []
    for seam in seams:
        if not set(seam).intersection(UNION):
            continue
        for coin in (0, 1):
            images.append(tuple((index[native.native(state, seam, coin)[0]],
                native.native(state, seam, coin)[1][1], state[seam[0]]-state[seam[1]])
                for state in states))
    checked, exceptional, smallest = 0, 0, None
    for law, label, shift in (('C', 1, 8), ('E', 2, 1)):
        for pairs in sectors.values():
            positions, size = {pair: i for i, pair in enumerate(pairs)}, len(pairs)
            matrix = [[0]*size for _ in range(size)]
            for column, (i, j) in enumerate(pairs):
                for image in images:
                    left, right = image[i], image[j]
                    if left[label] == right[label]:
                        matrix[positions[left[0], right[0]]][column] += 1
            if any(matrix[i][j] != matrix[j][i] for i, j in product(range(size), repeat=2)):
                raise ArithmeticError('incident operator block is not self-adjoint')
            floor = min(row[i]-sum(value for j, value in enumerate(row) if j != i)
                        for i, row in enumerate(matrix))
            if floor < shift:
                if law != 'E' or size != 10:
                    raise ArithmeticError('an unexpected block requires a new positivity proof')
                m = [[F(value-int(i == j)) for j, value in enumerate(row)] for i, row in enumerate(matrix)]
                for k in range(size):
                    pivot = m[k][k]
                    if pivot <= 0:
                        raise ArithmeticError('exact incident spectral-floor elimination failed')
                    smallest = pivot if smallest is None else min(smallest, pivot)
                    for i in range(k+1, size):
                        for j in range(i, size):
                            m[i][j] -= m[i][k]*m[k][j]/pivot
                            m[j][i] = m[i][j]
                exceptional += 1
            checked += size
    dimensions = tuple(sorted(Counter(map(len, sectors.values())).items()))
    if checked != 2*66**2 or dimensions != ((1, 2970), (10, 132), (66, 1)):
        raise ArithmeticError('the complete total-two operator partition changed')
    return checked, dimensions, exceptional, (smallest.numerator, smallest.denominator)


def total_two_certificate():
    checked, dimensions, exceptional, pivot = _total_two_data()
    return {'all_matrix_units_per_candidate': checked//2, 'block_dimensions': dimensions,
        'incident_seams': 12, 'incident_hs_lower': {'C': F(2, 15), 'E': F(1, 60)},
        'exceptional_exact_eliminations': exceptional, 'minimum_elimination_pivot': F(*pivot)}


def triangle_certificate():
    """Exact full-binary operator certificate, including inter-number coherences.

    Count matrices sum the two original coin branches without the factor
    one half. Their symmetric annihilators certify their spectra without
    a numerical eigensolver. The single-edge check permits adding every
    remaining incident edge; spectator ports contribute identity factors.
    """
    ports, edges = (0, 1, 5), ((0, 1), (0, 5), (1, 5))
    seams = core.source().seams
    if (len(seams) != 30 or any(edge not in seams for edge in edges)
            or any(not set(edge).intersection(UNION) for edge in edges)):
        raise ArithmeticError('the full-binary incident triangle is not in the actual source')
    bits = tuple(product((0, 1), repeat=3))
    states = tuple(tuple(row[ports.index(i)] if i in ports else 0 for i in range(12))
                   for row in bits)
    index = {state: i for i, state in enumerate(states)}
    identity = np.eye(64, dtype=object)
    matrices = {}
    spectra = {}
    for law in ('C', 'E'):
        edge_counts = []
        for edge in edges:
            matrix = np.zeros((64, 64), dtype=object)
            for coin in (0, 1):
                images = []
                for state in states:
                    out, record = native.native(state, edge, coin)
                    if (out not in index or sum(out) != sum(state)
                            or record != (coin, state[edge[0]]-out[edge[0]])):
                        raise ArithmeticError('the original binary triangle branch changed')
                    label = record[1] if law == 'C' else state[edge[0]]-state[edge[1]]
                    images.append((index[out], label))
                for i, j in product(range(8), repeat=2):
                    if images[i][1] == images[j][1]:
                        matrix[8*images[i][0]+images[j][0], 8*i+j] += 1
            edge_counts.append(matrix)
        count = sum(edge_counts, np.zeros((64, 64), dtype=object))
        roots = (2, 3, 4, 6) if law == 'C' else (0, 2, 3, 6)
        single_roots = (0, 1, 2) if law == 'C' else (0, 2)
        for matrix, allowed in [(m, single_roots) for m in edge_counts]+[(count, roots)]:
            annihilator = identity.copy()
            for root in allowed:
                annihilator = annihilator@(matrix-root*identity)
            if not np.array_equal(matrix, matrix.T) or np.any(annihilator):
                raise ArithmeticError('the exact binary quantum annihilator failed')
        # Trace preservation is checked on all matrix units, including
        # off-diagonal units connecting different particle-number sectors.
        trace_row = np.array([int(i == j) for i, j in product(range(8), repeat=2)], dtype=object)
        if not np.array_equal(trace_row@count, 6*trace_row):
            raise ArithmeticError('the triangle omitted native coin or record branches')
        multiplicities = []
        for root in roots:
            projector = identity.copy()
            denominator = 1
            for other in roots:
                if other != root:
                    projector = projector@(count-other*identity)
                    denominator *= root-other
            multiplicity = F(int(np.trace(projector)), denominator)
            if multiplicity.denominator != 1 or multiplicity <= 0:
                raise ArithmeticError('the exact triangle spectral support changed')
            multiplicities.append((root, multiplicity.numerator))
        if sum(size for _, size in multiplicities) != 64:
            raise ArithmeticError('the triangle spectral certificate is incomplete')
        matrices[law], spectra[law] = count, tuple(multiplicities)
    left, right = bits.index((1, 0, 0)), bits.index((0, 1, 0))
    witness = identity[:, 8*left+right]
    if (np.any(matrices['E']@witness)
            or not np.array_equal(matrices['C']@witness, 2*witness)):
        raise ArithmeticError('the C floor or E zero-mode control changed')
    return {'triangle_ports': ports, 'triangle_seams': edges,
        'all_binary_states': 8, 'all_matrix_units_per_candidate': 64,
        'count_spectra_with_multiplicity': spectra,
        'single_seam_count_annihilator_roots': {'C': (0, 1, 2), 'E': (0, 2)},
        'triangle_count_annihilator_roots': {'C': (2, 3, 4, 6), 'E': (0, 2, 3, 6)},
        'scheduled_count_denominator': 2*len(seams),
        'full_binary_C_incident_hs_lower': F(1, len(seams)),
        'witness_triangle_bits': (bits[left], bits[right]),
        'C_witness_count_eigenvalue': 2, 'E_witness_count_eigenvalue': 0}


def record_family_obstruction():
    """An exact CP current instrument can still have negative clock spectrum."""
    values = {}
    for law in ('C', 'E', 'C_signed'):
        instrument = family.INSTRUMENTS[law]()
        if not family.completeness(instrument) or not family.definite_currents(instrument):
            raise ArithmeticError('the classified current instrument lost its CP normalization')
        table = family.extract_table(instrument)
        g = table['g']
        if g.b or g.c or g.d:
            raise ArithmeticError('the retained clock witness requires a rational multiplier')
        # Eight pair-incident seams meet one of0,1; seam(0,1) kills this
        # coherence. All21 remote seams fix it. The same-source table g is
        # the empty-versus-waiting multiplier at each of those eight seams.
        pairs = source('C', PAIR)
        if sum(bool(set(e).intersection(PAIR)) for e in pairs.seams) != 9 or (0, 1) not in pairs.seams:
            raise ArithmeticError('the canonical pair witness graph changed')
        values[law] = 8*g.a/9
    if values != {'C': F(4, 9), 'E': F(0), 'C_signed': F(-4, 9)}:
        raise ArithmeticError('the exact zero/negative effective-clock controls changed')
    return {'pair_clock_coherence_eigenvalues': values,
        'all_three_are_complete_CP_current_instruments': True,
        'E_pair_has_no_finite_full_operator_log': True,
        'C_signed_pair_is_not_positive_on_Hilbert_Schmidt_space': True}
