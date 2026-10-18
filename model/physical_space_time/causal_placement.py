"""A source-record placement in the registered rank-three response carrier.

The spatial readout is a declared seam midpoint. Time is dependency height,
not seconds or a primitive observer's incident clock. Forward cone inclusion
is derived; the converse is tested on actual native branches and fails.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

import sympy as sp

from ..common_source import slow_geometry
from ..observer import native
from . import construction as passage

SOURCE_PATHS = tuple(sorted(set((*passage.SOURCE_PATHS,
    'model/physical_space_time/causal_placement.py',
    'Lean/Screen/PortFrameGram.lean'))))
SPEED_SQUARED = (1+sp.sqrt(5)/5)/2


def _seam(edge):
    if (type(edge) not in (tuple, list) or len(edge) != 2
            or any(type(i) is not int for i in edge)
            or tuple(edge) not in passage.recovered_carrier()[0]):
        raise ValueError('an increasing native seam is required')
    return tuple(edge)


@lru_cache(maxsize=1)
def gram():
    """Existing record-recovered Gram, also bound to the registered Lean chart."""
    slow_geometry.registered_binding()
    c = slow_geometry.source()
    value = sp.ImmutableMatrix([[sp.Rational(a)+sp.Rational(b)*sp.sqrt(5)
                                for a, b in row] for row in c.intrinsic_gram_exact()])
    _, laplacian = passage.recovered_carrier()
    if ((laplacian*value-(5-sp.sqrt(5))*value).applyfunc(sp.expand) != sp.zeros(12)
            or (value*value-4*value).applyfunc(sp.expand) != sp.zeros(12)
            or set(value.diagonal()) != {1}):
        raise ArithmeticError('the source Gram does not match the recovered slow carrier')
    return value


def distance_squared(first, second):
    """Squared distance between seam midpoints, using only the response Gram."""
    first, second = _seam(first), _seam(second)
    coefficients = [sp.Rational(first.count(i)-second.count(i), 2) for i in range(12)]
    support = [i for i, value in enumerate(coefficients) if value]
    return sp.expand(sum(coefficients[i]*gram()[i, j]*coefficients[j]
                         for i in support for j in support))


def dependencies(tape):
    """Latest endpoint-use parents from a complete ordered seam tape.

    The tape includes zero-current events. Restricting it to one observer's
    incident records does not preserve these full-source parent relations.
    This is the endpoint DAG for an exogenous scheduler. Adaptive controller
    reads or preparation dependencies require additional edges and bounds.
    """
    if type(tape) not in (tuple, list):
        raise ValueError('a complete finite seam tape is required')
    last, rows = {}, []
    for number, original in enumerate(tape):
        edge = _seam(original)
        parents = tuple(sorted({last[port] for port in edge if port in last}))
        height = 0 if not parents else 1+max(rows[p]['height'] for p in parents)
        ancestors = set(parents)
        for parent in parents:
            ancestors.update(rows[parent]['ancestors'])
        rows.append({'seam': edge, 'parents': parents, 'height': height,
                     'ancestors': tuple(sorted(ancestors))})
        for port in edge:
            last[port] = number
    return rows


def interval_squared(first, second):
    """The (+---) interval with dimensionless time sqrt(c²)*height."""
    for row in (first, second):
        if (type(row) is not dict or 'seam' not in row or 'height' not in row
                or type(row['height']) is not int or row['height'] < 0):
            raise ValueError('a native seam and original nonnegative integer height are required')
    return sp.expand(SPEED_SQUARED*(second['height']-first['height'])**2
                     -distance_squared(first['seam'], second['seam']))


def finite_readout():
    """All actual fixed-total states, plus every native current telescoping step."""
    c = slow_geometry.source()
    states = native.states()
    readouts = {state: slow_geometry.load_readout(state) for state in states}
    checked = 0
    for before in states:
        for seam in passage.recovered_carrier()[0]:
            for coin in (0, 1):
                after, (_, current) = native.native(before, seam, coin)
                difference = tuple(c.q5_sub(b, a) for a, b in
                                   zip(readouts[before], readouts[after], strict=True))
                if difference != slow_geometry.current_displacement(seam, current):
                    raise ArithmeticError('actual current readout does not telescope')
                checked += 1
    return {'total_two_states': len(states), 'distinct_rank_three_images': len(set(readouts.values())),
            'binary_states': len(native.states(binary=True)),
            'distinct_binary_images': len({readouts[s] for s in native.states(binary=True)}),
            'native_branch_telescoping_checks': checked,
            'unit_centroid_normalization': 'raw_load_readout / (2*sqrt((5+sqrt(5))/2))',
            'normalized_load_centroid_in_convex_hull_of_unit_ports': True,
            'waiting_on_fixed_total_carrier_realizes_full_signed_control_completion': False}


def run():
    seams, _ = passage.recovered_carrier()
    sharing = [(a, b) for a in seams for b in seams if set(a) & set(b)]
    distances = [distance_squared(a, b) for a, b in sharing]
    if any(sp.simplify(SPEED_SQUARED-d).is_nonnegative is not True for d in distances):
        raise ArithmeticError('a native endpoint dependency exceeds the placement bound')
    sharp = next((pair for pair, d in zip(sharing, distances, strict=True)
                  if sp.expand(d-SPEED_SQUARED) == 0), None)
    if sharp is None:
        raise ArithmeticError('the stated edge bound is not sharp')
    tape = ((0, 1), (2, 3), (0, 1), (0, 1))
    rows = dependencies(tape)
    if (rows[3]['ancestors'] != (0, 2) or rows[3]['height']-rows[1]['height'] != 2
            or interval_squared(rows[1], rows[3]) != 0):
        raise ArithmeticError('the actual-seam cone-converse witness changed')
    for later in rows:
        for previous in later['ancestors']:
            if interval_squared(rows[previous], later).is_nonnegative is not True:
                raise ArithmeticError('a generated ancestor left the forward cone')
    before = tuple(int(i in (0, 2)) for i in range(12))
    branches = []
    for coins in product((0, 1), repeat=len(tape)):
        state, records = before, []
        for seam, coin in zip(tape, coins, strict=True):
            state, record = native.native(state, seam, coin)
            records.append(list(record))
        branches.append({'coins': list(coins), 'records': records,
                         'terminal_occupied_ports': [i for i, value in enumerate(state) if value]})
    return {'contract': {
                'spatial_readout': 'declared seam midpoint in the native response Gram',
                'time_parameter': 'canonical endpoint-dependency height times sqrt(c_squared)',
                'complete_source_seam_tape_including_zero_current_events_required': True,
                'predeclared_exogenous_schedule_required_for_endpoint_only_dependency_graph': True,
                'arbitrary_controller_or_preparation_dependency_edges_covered': False,
                'endpoint_parents_are_structural_not_proven_signal_influences': True,
                'physical_seconds_or_unique_spatial_readout_inferred': False},
            'forward_cone': {'c_squared': str(SPEED_SQUARED),
                'ordered_sharing_seam_pairs_checked': len(sharing),
                'sharp_pair': [list(s) for s in sharp],
                'argument': 'Each parent shares an endpoint and raises height by at least one; '
                            'the Gram norm triangle inequality extends the edge bound to every ancestor.',
                'every_finite_complete_native_tape_satisfies_forward_placement': True},
            'converse_control': {'seam_tape': [list(s) for s in tape],
                'parents': [list(r['parents']) for r in rows],
                'heights': [r['height'] for r in rows],
                'tested_pair': [1, 3], 'future_interval_squared': '0',
                'earlier_is_ancestor': False, 'all_coin_branches': branches,
                'zero_current_events_kept': True},
            'finite_current_readout': finite_readout(),
            'acceptance': {'forward_placement_geometric_hypotheses_checked': True,
                'authenticated_semantic_log_or_Lean_term_emitted': False,
                'cone_implies_ancestry': False, 'physical_Lorentz_symmetry_or_clock_derived': False,
                'continuum_spacetime_derived': False}, 'all_pass': True}
