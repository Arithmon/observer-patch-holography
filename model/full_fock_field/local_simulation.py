"""Local observer schedules for the supplied all-number Weyl/contact field.

All error arithmetic is exact rational arithmetic.  The guarantees concern
ideal source words at fixed lattice spacing and Hamiltonian model time.
They neither prepare the input restriction nor bound floating/hardware error.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from math import factorial, isqrt, prod

from .bounds import _parameters


def _time(value):
    if type(value) not in (int, F):
        raise ValueError('exact signed Hamiltonian model duration required')
    return F(value)


def _steps(value):
    if type(value) is not int or value < 1:
        raise ValueError('positive integer split count required')
    return value


def factorial_tail_upper(value, start):
    """Exact enclosure of sum_(n>=start) value^n/n!, with a geometric tail."""
    if (type(value) not in (int, F) or not 0 <= value <= 32
            or type(start) is not int or not 1 <= start <= 4097):
        raise ValueError('exact exponent in 0..32 and integer tail order in 1..4097 required')
    value = F(value)
    twice = 2*value
    stop = max(start, -(-twice.numerator//twice.denominator))
    term, total = value**start/factorial(start), F()
    for n in range(start, stop):
        total += term
        term *= value/(n+1)
    return total+term/(1-value/(stop+1))


def local_tail_boundary_bound(radius, support_radius, support_sites, duration, speed=F(1, 8)):
    """Integrate the actual commutator-chain tail, retaining its factorial decay.

    Bound ||tau_ambient(A)-tau_radius(A)|| for ||A||<=1, also at every
    untouched-reference matrix level.  The supplied action is unchanged.
    """
    if (type(radius) is not int or type(support_radius) is not int
            or not 0 <= support_radius < radius <= 4096
            or type(support_sites) is not int
            or not 1 <= support_sites <= (2*support_radius+1)**3
            or type(speed) not in (int, F) or speed < 0):
        raise ValueError('strictly nested bounded cubes and exact nonnegative speed required')
    time, speed = abs(_time(duration)), F(speed)
    tail = factorial_tail_upper(24*speed*time, radius-support_radius+1)
    return min(F(2), F(support_sites*(2*radius+1)**2, 4)*tail)


def extensive_split_coefficient(sides, speed=F(1, 8), interaction=F(1)):
    """K with ||S(t/n)^n-exp(-itH)|| <= K |t|^3/n^2.

    The degree bound uses only overlapping even CAR supports and is linear
    in spatial edge count.  This arithmetic family has no Fock-array or
    FieldGeometry allocation; those executors have separate finite domains.
    """
    _, edges, c, u = _parameters(sides, speed, interaction)
    incident = 3*c+u
    return edges*c*incident*(c/6+2*incident)


@lru_cache(maxsize=8)
def _overlap_counts(sides):
    # The catalogue is the actual field schedule, not a surrogate lattice.
    # Norms are c/2 for hopping terms and U for the commuting onsite terms.
    from .field import FieldGeometry
    geometry = FieldGeometry(sides, speed=1., interaction=1.)
    terms = [(frozenset((p, q)), 0) for group in geometry.groups for p, q, _ in group]
    hopping_count = len(terms)
    terms.extend((frozenset((2*x, 2*x+1)), 1) for x in range(geometry.volume))
    incident = [set() for _ in range(geometry.modes)]
    for i, (support, _) in enumerate(terms):
        for mode in support:
            incident[mode].add(i)
    repeated, suffix = [0]*3, [0]*3
    for i in range(hopping_count):
        support_i = terms[i][0]
        neighbors = set().union(*(incident[p] for p in support_i))
        for j in neighbors:
            if j <= i:
                continue
            support_j, contact_j = terms[j]
            repeated[contact_j] += 1
            outside = set().union(*(incident[p] for p in support_i | support_j))
            for k in outside:
                if k > i:
                    suffix[contact_j+terms[k][1]] += 1
    # Contacts follow every hop, and their suffix consists of disjoint
    # commuting contacts.  Thus they contribute no outer Strang term.
    return tuple(repeated), tuple(suffix), hopping_count, geometry.volume


def overlap_split_certificate(sides, speed=F(1, 8), interaction=F(1)):
    """Exact ordered support count on the existing bounded FieldGeometry.

    The count is an upper bound, not an assertion that every overlapping
    commutator is nonzero.  No global operator matrix is constructed.
    """
    _, _, c, u = _parameters(sides, speed, interaction)
    repeated, suffix, hops, contacts = _overlap_counts(sides)
    weights = tuple((c/2)**(3-q)*u**q for q in range(3))
    coefficient = sum((F(a, 6)+F(b, 3))*w
                      for a, b, w in zip(repeated, suffix, weights, strict=True))
    coarse = extensive_split_coefficient(sides, c, u)
    if coefficient > coarse:
        raise ArithmeticError('ordered support count exceeds its degree bound')
    return {'coefficient': coefficient, 'degree_coefficient': coarse,
            'repeated_outer_counts_by_contact_degree': repeated,
            'suffix_outer_counts_by_contact_degree': suffix,
            'hopping_terms': hops if c else 0, 'contact_terms': contacts,
            'all_number_sectors': True, 'global_operator_allocated': False}


def split_bound(sides, duration, steps, speed=F(1, 8), interaction=F(1), *, refined=False):
    """Uncapped ideal operator error; it also bounds half diamond distance."""
    duration, steps = _time(duration), _steps(steps)
    if type(refined) is not bool:
        raise ValueError('explicit Boolean support-refinement choice required')
    coefficient = (overlap_split_certificate(sides, speed, interaction)['coefficient']
                   if refined else extensive_split_coefficient(sides, speed, interaction))
    return coefficient*abs(duration)**3/steps**2


def history_split_bound(sides, segments, speed=F(1, 8), interaction=F(1), *, refined=False):
    """Complete finite-box history, with all instrument outcomes retained."""
    if type(segments) not in (tuple, list) or not segments:
        raise ValueError('nonempty finite sequence of duration/count pairs required')
    total = F()
    for segment in segments:
        if type(segment) not in (tuple, list) or len(segment) != 2:
            raise ValueError('duration/count pair required')
        total += split_bound(sides, *segment, speed, interaction, refined=refined)
    return min(F(1), total)


def source_schedule_counts(sides, steps, speed=F(1, 8)):
    """Exact counts of the existing routed templates, without a global array.

    A completed hopping or FSWAP block has (22 single,4 pair) pulses;
    contact has (8 single,2 pair).  Every request, including zero angles,
    incurs its existing load/export pair.  Preparation, region restriction,
    observer interventions and readout remain separately supplied services.
    """
    volume, edges, c, _ = _parameters(sides, speed, F())
    steps = _steps(steps)
    hopping = 4*steps*edges if c else 0
    route_one_direction = 0
    if c:
        for axis in range(3):
            axis_edges = (sides[axis]-1)*prod(sides[j] for j in range(3) if j != axis)
            stride = prod(sides[axis+1:])
            route_one_direction += axis_edges*(8*stride-4)
    swaps, contacts = 2*steps*route_one_direction, 2*steps*volume
    single, pair = 22*(hopping+swaps)+8*contacts, 4*(hopping+swaps)+2*contacts
    pulses = single+pair
    return {'native_pulses': pulses, 'single_code_pulses': single,
            'pair_code_pulses': pair, 'proper_code_load_events': pulses,
            'proper_code_export_events': pulses, 'proper_code_service_events': 2*pulses,
            'completed_hopping_blocks': hopping, 'completed_contact_blocks': contacts,
            'routing_FSWAP_calls': swaps, 'source_qubits': 2*volume,
            'geometry_shape_admitted': max(sides) <= 128 and volume <= 4096,
            'explicit_state_geometry_admitted': 2*volume <= 16,
            'floating_parameter_admission_checked': False,
            'global_Fock_array_allocated': False,
            'preparation_restriction_intervention_readout_costs_included': False,
            'floating_or_hardware_error_included': False}


def _ceil_sqrt(value):
    result = isqrt(value.numerator//value.denominator)
    return result+int(result*result*value.denominator < value.numerator)


def local_observer_schedule(durations, epsilon, *, support_radius=0, support_sites=1,
                            speed=F(1, 8), interaction=F(1), max_radius=4096,
                            refined=False):
    """A complete local history with successively smaller simulation cubes.

    Durations are chronological.  Even local interventions act on the fixed
    observer support and retained memory; every outcome remains present.
    Working backwards, each duration gets a buffer around the full support
    of the later simulated effect.  This is essential for history control.
    The returned half-diamond bound is uniform in ambient box, number sector
    and untouched reference.  It assumes the correct initial restriction.
    """
    if type(durations) not in (tuple, list) or not durations:
        raise ValueError('nonempty finite chronological duration sequence required')
    times = tuple(_time(t) for t in durations)
    if (type(epsilon) not in (int, F) or not 0 < epsilon <= 1
            or type(support_radius) is not int or type(support_sites) is not int
            or not 0 <= support_radius < 4096
            or not 1 <= support_sites <= (2*support_radius+1)**3
            or type(max_radius) is not int or not support_radius <= max_radius <= 4096
            or type(refined) is not bool):
        raise ValueError('exact error budget, bounded nested cubes and explicit refinement required')
    _, _, c, u = _parameters((1, 1, 1), speed, interaction)
    epsilon, per_stage = F(epsilon), F(epsilon)/len(times)
    if any(48*c*abs(t) > 64 for t in times):
        raise ValueError('duration exceeds the exact locality enclosure domain; subdivide it')
    outer, sites, backwards = support_radius, support_sites, []
    for duration in reversed(times):
        if not c or not duration:
            radius, boundary = outer, F()
        else:
            for radius in range(outer+1, max_radius+1):
                boundary = local_tail_boundary_bound(radius, outer, sites, duration, c)
                if boundary <= per_stage:
                    break
            else:
                raise ValueError('no padding radius meets the budget in the declared arithmetic domain')
        sides = (2*radius+1,)*3
        admitted = max(sides) <= 128 and prod(sides) <= 4096
        use_refined = refined and admitted
        coefficient = (overlap_split_certificate(sides, c, u)['coefficient'] if use_refined
                       else extensive_split_coefficient(sides, c, u))
        steps = max(1, _ceil_sqrt(2*coefficient*abs(duration)**3/per_stage))
        split = coefficient*abs(duration)**3/steps**2
        bound = boundary/2+split
        if bound > per_stage:
            raise ArithmeticError('exact local history allocation exceeded its budget')
        backwards.append({'duration': duration, 'input_radius': radius, 'output_radius': outer,
                          'tested_output_support_sites': sites, 'steps': steps,
                          'boundary_operator_bound': boundary, 'split_operator_bound': split,
                          'half_diamond_bound': bound, 'split_coefficient': coefficient,
                          'coefficient_method': 'ordered overlaps' if use_refined else 'local degree',
                          'resources': source_schedule_counts(sides, steps, c)})
        outer, sites = radius, prod(sides)
    stages = tuple(reversed(backwards))
    return {'stages': stages, 'initial_restriction_radius': outer,
            'terminal_support_radius': support_radius,
            'half_diamond_bound': sum((stage['half_diamond_bound'] for stage in stages), F()),
            'requested_epsilon': epsilon,
            'native_pulses': sum(stage['resources']['native_pulses'] for stage in stages),
            'ambient_volume_parameter_required': False,
            'initial_restricted_state_prepared': False,
            'all_number_sectors_and_references_retained': True,
            'outcome_conditioning_or_postselection': False,
            'spatial_continuum_or_physical_clock_claimed': False}
