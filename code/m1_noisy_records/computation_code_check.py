"""Independent truth-table intersection construction and GF(2) checks."""

import hashlib
import itertools
import math
from .independent import keys, need


def monomials(m, r):
    all_points = (1 << (1 << m))-1
    coordinates = []
    for bit in range(m):
        stripe = (1 << (1 << bit))-1
        coordinates.append(sum(stripe << start
                               for start in range(1 << bit, 1 << m, 2 << bit)))
    result = []
    for subset in range(1 << m):
        if subset.bit_count() > r:
            continue
        value = all_points
        for bit in range(m):
            if subset >> bit & 1:
                value &= coordinates[bit]
        result.append(value)
    return result


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            leading = row.bit_length()-1
            if leading not in pivots:
                pivots[leading] = row
                break
            row ^= pivots[leading]
    return len(pivots)


def distance_bound(m, r):
    """Induction f=p+x_m q: q=0 doubles; q!=0 bounds by weight(q)."""
    if r == 0:
        return 1 << m
    if r == m:
        return 1
    return min(2*distance_bound(m-1, r), distance_bound(m-1, r-1))


def small_complete_distances():
    result = []
    for m, r in ((3, 1), (4, 1), (5, 2)):
        code = [0]
        for row in monomials(m, r):
            code += [word ^ row for word in code]
        actual = min(word.bit_count() for word in code if word)
        need(actual == distance_bound(m, r), 'complete small RM distance')
        result.append((m, r, len(code), actual))
    return result


def verify(row):
    keys(row, 'variables degree extended_length extended_dimension extended_distance '
         'quantum_length quantum_dimension quantum_distance corrects procedure_spread_bound '
         'correction_then_gate_spread_bound tolerated_faults_per_rectangle generator_sha256')
    need(all(type(v) is int for k, v in row.items() if k != 'generator_sha256'), 'integer code parameters')
    need(row['variables'] == 11 and row['degree'] == 5, 'fixed computation code')
    m, r = row['variables'], row['degree']
    generators = monomials(m, r)
    length, dimension = 1 << m, sum(math.comb(m, j) for j in range(r+1))
    need(rank(generators) == dimension == length//2, 'self-dual dimension')
    need(all(g.bit_count() % 4 == 0 for g in generators), 'doubly even generators')
    need(all((a & b).bit_count() % 2 == 0 for a in generators for b in generators),
         'self-orthogonality')
    distance = distance_bound(m, r)
    # The degree-r monomial is a distance witness and is 1 at the puncture.
    witness = monomials(m, r)[(1 << r)-1]
    need(witness.bit_count() == distance and witness >> (length-1) == 1, 'distance witness')
    mask = (1 << (length-1))-1
    punctured = [g & mask for g in generators]
    # Every generator is 1 at the removed all-ones point. Pairing each with
    # the constant generator gives the shortened stabilizer code basis.
    shortened = [g ^ punctured[0] for g in punctured[1:]]
    need(rank(punctured) == dimension and rank(shortened) == dimension-1, 'puncture/shorten rank')
    need(all(g.bit_count() % 4 == 0 for g in shortened), 'CSS phase class')
    need(all((a & b).bit_count() % 2 == 0 for a in shortened for b in punctured),
         'shortened code equals punctured dual')
    corrects = (distance-2)//2
    for key, expected in dict(extended_length=length, extended_dimension=dimension,
                             extended_distance=distance, quantum_length=length-1,
                             quantum_dimension=(length-1)-2*(dimension-1),
                             quantum_distance=distance-1, corrects=corrects,
                             procedure_spread_bound=4, correction_then_gate_spread_bound=4*4,
                             tolerated_faults_per_rectangle=corrects//16).items():
        need(row[key] == expected, 'derived computation-code parameter '+key)
    raw = b''.join(g.to_bytes(length//8, 'little') for g in generators)
    need(row['generator_sha256'] == hashlib.sha256(raw).hexdigest(), 'generator custody')
    small_complete_distances()
