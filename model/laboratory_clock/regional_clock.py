"""Collective incident clocks as unions of visible records on one carrier.

Overlapping port regions can expose the same event set. Such aliases are one
clock, even when their port labels differ. Distinct event sets have a derived
renewal law and a positive normalized fluctuation floor; no noise parameter
or independent replica is added.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json

import sympy as sp

from . import event_geometry as e

SOURCE_PATHS = tuple(sorted(set((*e.SOURCE_PATHS,
    'model/laboratory_clock/regional_clock.py'))))


def _region(region):
    if (type(region) is not tuple or not region or len(region) > 12
            or any(type(port) is not int or not 0 <= port < 12 for port in region)
            or tuple(sorted(set(region))) != region):
        raise ValueError('a nonempty increasing tuple of distinct native ports is required')
    return region


def event_mask(region):
    """A seam ticks once when either endpoint is in the region."""
    region = _region(region)
    return sum(1 << index for index, seam in enumerate(e._source().seams)
               if any(port in region for port in seam))


def incidence_counts(first, second):
    left, right = event_mask(first), event_mask(second)
    return left.bit_count(), right.bit_count(), (left & right).bit_count()


def parameters(a, b, common):
    """Exact event-set law from integer incidence counts, including aliases."""
    if (any(type(x) is not int for x in (a, b, common))
            or not 1 <= a <= 30 or not 1 <= b <= 30
            or not max(0, a+b-30) <= common <= min(a, b)):
        raise ValueError('valid exact counts of two nonempty thirty-seam event sets required')
    return {'mean_B_per_A': F(b, a),
            'variance_B_per_A': F(b*(a+b-2*common), a*a),
            'normalized_ratio_variance_per_A': F(a+b-2*common, b),
            'terminal_B_probability': F(common, a),
            'geometric_success_probability': F(a, a+b-common),
            'identical_event_sets': a == b == common}


def pgf(first, second, z=F(1), count=1):
    """B-count generating function at count A ticks; all waiting tails summed."""
    z = e._argument(z); e._integer(count, 1, 64)
    a, b, common = incidence_counts(first, second)
    return ((a-common+common*z)/(a+b-common-(b-common)*z))**count


def count_law(first, second, count=1, maximum=12):
    """Complete finite bins and an explicit overflow event for B's count."""
    e._integer(count, 1, 64); e._integer(maximum, 0, 128)
    a, b, common = incidence_counts(first, second)
    p = parameters(a, b, common)
    terminal, success = p['terminal_B_probability'], p['geometric_success_probability']
    q = 1-success
    single = [(1-terminal)*success]
    single.extend((1-terminal)*success*q**k+terminal*success*q**(k-1)
                  for k in range(1, maximum+1))
    result = [F(1)]+[F(0)]*maximum
    for _ in range(count):
        result = [sum((result[j]*single[k-j] for j in range(k+1)), F(0))
                  for k in range(maximum+1)]
    return tuple(result), 1-sum(result, F(0))


@lru_cache(maxsize=1)
def _catalogue():
    """Integer-only cache: every region, every distinct event set and alias."""
    masks = {}
    for word in range(1, 1 << 12):
        region = tuple(port for port in range(12) if word >> port & 1)
        masks.setdefault(event_mask(region), []).append(region)
    return tuple((mask, tuple(sorted(regions, key=lambda x: (len(x), x))))
                 for mask, regions in sorted(masks.items()))


@lru_cache(maxsize=1)
def pair_classes():
    """All ordered distinct event-set pairs, grouped by (a,b,common)."""
    catalogue = _catalogue()
    classes = Counter()
    for index, (left, _) in enumerate(catalogue):
        a = left.bit_count()
        for right, _ in catalogue[index+1:]:
            b, common = right.bit_count(), (left & right).bit_count()
            classes[a, b, common] += 1
            classes[b, a, common] += 1
    return tuple((*key, amount) for key, amount in sorted(classes.items()))


def joint_field(first, second, initial_port, final_port, z=F(1), count=2, *, oracle=False):
    """Selected local-bit PGF on the same collective clocks and native source.

    Initial and final ports belong to their respective visible regions. The
    preparation is uniform binary total-two. The full-state oracle enumerates
    actual native branches, checks own-current readback, then sums all waits.
    The reduced prediction derives each seam's density map independently.
    """
    first, second = _region(first), _region(second)
    e._ports(initial_port, final_port)
    if initial_port not in first or final_port not in second or type(oracle) is not bool:
        raise ValueError('local field ports must be admitted by their declared regions')
    z = e._argument(z); e._integer(count, 1, 8)
    size = 66 if oracle else 12
    remote, incident = sp.zeros(size), sp.zeros(size)
    if oracle:
        states = tuple(e.native.states(binary=True)); index = {state: i for i, state in enumerate(states)}
        for column, before in enumerate(states):
            for seam in e._source().seams:
                tick_a, tick_b = any(p in first for p in seam), any(p in second for p in seam)
                target = incident if tick_a else remote
                for coin in (0, 1):
                    after, (_, current) = e.native.native(before, seam, coin)
                    for port in (initial_port, final_port):
                        decoded = before[port]+current*(int(port == seam[1])-int(port == seam[0]))
                        if decoded != after[port]:
                            raise ArithmeticError('collective local current readback failed')
                    target[index[after], column] += sp.Rational(z)**int(tick_b)/60
        state = sp.Matrix(66, 2, lambda i, bit: sp.Rational(int(states[i][initial_port] == bit), 66))
        inverse = (sp.eye(66)-remote).inv(method='DM')
        for _ in range(count):
            state = incident*inverse*state
        result = {(x, y): F(sum(state[i, x] for i in range(66) if states[i][final_port] == y))
                  for x in (0, 1) for y in (0, 1)}
    else:
        for seam in e._source().seams:
            a, b = seam
            vector = sp.zeros(12, 1); vector[a] = 1; vector[b] = -1
            step = sp.Rational(z)**int(any(p in second for p in seam))*(sp.eye(12)-vector*vector.T/2)/30
            if any(p in first for p in seam):
                incident += step
            else:
                remote += step
        kernel = incident*(sp.eye(12)-remote).inv()
        mass = pgf(first, second, z, count)
        both = mass/36+F(e.CHI*(kernel**count*e.P0)[final_port, initial_port])
        result = {(1, 1): both, (1, 0): mass/6-both,
                  (0, 1): mass/6-both, (0, 0): 2*mass/3+both}
    if any(x < 0 for x in result.values()) or sum(result.values()) != pgf(first, second, z, count):
        raise ArithmeticError('complete local-bit generating weights lost mass')
    return result


def run():
    catalogue, classes = _catalogue(), pair_classes()
    if sum(len(regions) for _, regions in catalogue) != 4095:
        raise ArithmeticError('some nonempty region was omitted')
    unequal = min(parameters(a, b, c)['normalized_ratio_variance_per_A'] for a, b, c, _ in classes)
    equal = min(parameters(a, b, c)['normalized_ratio_variance_per_A'] for a, b, c, _ in classes if a == b)
    if unequal != F(1, 30) or equal != F(2, 29):
        raise ArithmeticError('sharp collective clock fluctuation bounds changed')
    first = tuple(p for p in range(12) if p not in (0, 1))
    second = tuple(p for p in range(12) if p not in (0, 5))
    if incidence_counts(first, second) != (29, 29, 28):
        raise ArithmeticError('the maximally overlapping distinct-clock witness changed')
    comparisons = []
    for count, z in ((1, F(0)), (1, F(1, 2)), (2, F(1))):
        law = joint_field(first, second, 5, 1, z, count)
        if law != joint_field(first, second, 5, 1, z, count, oracle=True):
            raise ArithmeticError('collective clock/field prediction failed native replay')
        comparisons.append({'A_ticks': count, 'B_count_PGF_argument': str(z),
            'joint_bit_generating_weights': [{'initial_A_port_bit': a, 'final_B_port_bit': b,
                'weight': str(p)} for (a, b), p in sorted(law.items())]})
    aliases = next(regions for mask, regions in catalogue if mask == (1 << 30)-1)
    alias_first, alias_second = aliases[:2]
    if any(pgf(alias_first, alias_second, z, n) != z**n
           for n in (1, 2, 5) for z in (F(0), F(1, 2), F(1))):
        raise ArithmeticError('identical event sets no longer define one clock')
    finite, tail = count_law(first, second, count=2, maximum=8)
    return {'definition': 'a region clock ticks once per incident source event, including zero current',
        'renewal_PGF': '(a-c+c*z)/(a+b-c-(b-c)*z)',
        'mean_B_at_N_A_ticks': 'N*b/a',
        'normalized_ratio': 'a*B(N)/(b*N)',
        'normalized_ratio_variance': '(a+b-2*c)/(b*N)',
        'regions_examined': 4095, 'distinct_event_sets': len(catalogue),
        'ordered_distinct_event_set_pairs': sum(row[3] for row in classes),
        'distinct_incidence_classes': len(classes),
        'class_inventory_sha256': hashlib.sha256(json.dumps(classes, separators=(',', ':')).encode()).hexdigest(),
        'sharp_normalized_variance_coefficient': str(unequal),
        'sharp_equal_mean_rate_coefficient': str(equal),
        'unequal_rate_floor_witness': {'A': list(first), 'B': list(range(12)),
            'incidence_counts': list(incidence_counts(first, tuple(range(12))))},
        'all_event_clock_alias_regions': len(aliases),
        'alias_witness': [list(alias_first), list(alias_second)],
        'alias_count_identity': 'B(N)=N pathwise; different region names are not independent clocks',
        'distinct_clock_witness': {'A': list(first), 'B': list(second),
            'incidence_counts': [29, 29, 28], 'initial_field_port': 5, 'final_field_port': 1,
            'count_two_A_ticks_bins_zero_through_eight': list(map(str, finite)),
            'count_overflow_probability': str(tail), 'joint_field_comparisons': comparisons},
        'scope': {'same_twelve_port_source': True, 'independent_replica_averaging_added': False,
            'uniform_independent_seam_scheduler': True, 'ordered_regional_tape_union_and_coincidences_admitted': True,
            'collective_union_access_is_same_as_one_primitive_port_tape': False,
            'field_preparation_uniform_binary_total_two': True,
            'selected_classical_joint_field_only': True,
            'coherent_C_E_selection_or_atomic_clock_identification': False}}
