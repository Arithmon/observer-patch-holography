"""A source-record placement in the registered rank-three response carrier.

The spatial readout is a declared seam midpoint. Time is dependency height,
not seconds or a primitive observer's incident clock. Forward cone inclusion
is derived; the converse is tested on actual native branches and fails in
general, while explicit four-event fork-and-join diamonds on the recovered
twelve-port carrier realize two-way causal order equivalence and matching
metric causal intervals across all sixteen coin branches.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product
import ast
from pathlib import Path
import re

import sympy as sp

from ..common_source import slow_geometry
from ..observer import native
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import construction as passage
from tools import strict_json

SOURCE_PATHS = tuple(sorted(set((*passage.SOURCE_PATHS,
    'model/physical_space_time/causal_placement.py',
    'Lean/Geometry/NativeRepairSemanticPlacement.lean',
    'Lean/Geometry/SourceDerivedSpacetimeCarrier.lean',
    'Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean',
    'Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean',
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


def _loads(values):
    if (type(values) not in (tuple, list) or len(values) != 12
            or any(type(v) is not int or v < 0 for v in values)
            or sum(values) != 2):
        raise ValueError('a native nonnegative total-two twelve-port state is required')
    return tuple(values)


def authenticated_parent_ports(parent, child):
    """The four clauses of Lean ReadAfterWriteAt, without trusting citations.

    This relation alone does not validate a log's execution. In particular,
    an arbitrary writer label does not establish either write membership or
    continuity of the value actually read.
    """
    return tuple(p for p in child['causal_support']
                 if p in parent['write_set']
                 and child['before_writers'][p] == parent['event_id']
                 and child['before_values'][p] == parent['after_values'][p])


def verify_semantic_execution(initial, log, *, tape=None, coins=None):
    """Independently check the complete native CommitSpec execution contract.

    The producer calls native.native; this verifier reconstructs its payload
    by floor/ceiling averaging. Endpoint version stamps describe retained
    instrument invocations, including zero current. They are full execution
    instrumentation, not an intrinsic primitive observer's read interface.

    Without an external tape the check certifies the supplied finite execution,
    not that a longer experiment was fully disclosed. Supplying tape and coins
    also checks coverage of every declared invocation.
    """
    state, writers = _loads(initial), (None,)*12
    if type(log) not in (tuple, list):
        raise ValueError('a complete finite semantic log is required')
    if (tape is None) != (coins is None):
        raise ValueError('external coverage requires both the seam and coin tapes')
    if tape is not None:
        if (type(tape) not in (tuple, list) or type(coins) not in (tuple, list)
                or len(tape) != len(log) or len(coins) != len(log)):
            raise ValueError('the semantic log does not cover the complete declared experiment')
        if (any(type(coin) is not int or coin not in (0, 1) for coin in coins)
                or any(type(entry) is not dict or entry.get('seam') != list(_seam(seam))
                       or entry.get('coin') != coin
                       for entry, seam, coin in zip(log, tape, coins, strict=True))):
            raise ValueError('the semantic log does not match the declared seam/coin sequence')
    fields = {'event_id', 'seam', 'coin', 'record', 'before_values', 'after_values',
              'before_writers', 'after_writers', 'read_set', 'write_set', 'causal_support'}
    rows = []
    for event_id, entry in enumerate(log):
        if type(entry) is not dict or set(entry) != fields:
            raise ValueError('the complete semantic commit fields are required')
        if type(entry['event_id']) is not int or entry['event_id'] != event_id:
            raise ValueError('event identifiers must be fresh consecutive execution positions')
        seam = _seam(entry['seam'])
        coin = entry['coin']
        if type(coin) is not int or coin not in (0, 1):
            raise ValueError('a binary exogenous coin is required')
        for name in ('read_set', 'write_set', 'causal_support'):
            if (type(entry[name]) not in (list, tuple)
                    or any(type(p) is not int for p in entry[name]) or tuple(entry[name]) != seam):
                raise ValueError('the admitted native read, write and support sets are both endpoints')
        if _loads(entry['before_values']) != state:
            raise ValueError('the before-value snapshot does not thread the actual execution')
        for name in ('before_writers', 'after_writers'):
            if (type(entry[name]) not in (list, tuple) or len(entry[name]) != 12
                    or any(v is not None and (type(v) is not int or v < 0)
                           for v in entry[name])):
                raise ValueError('writer versions must be twelve optional event identifiers')
        if tuple(entry['before_writers']) != writers:
            raise ValueError('the before-writer snapshot does not thread the actual execution')
        a, b = seam
        total = state[a]+state[b]
        low, high = total//2, (total+1)//2
        expected = list(state)
        expected[a], expected[b] = (low, high) if coin == 0 else (high, low)
        record = (coin, state[a]-expected[a])
        if (type(entry['record']) not in (list, tuple) or len(entry['record']) != 2
                or any(type(v) is not int for v in entry['record'])
                or tuple(entry['record']) != record or _loads(entry['after_values']) != tuple(expected)):
            raise ValueError('the payload or current label is not the native endpoint repair')
        stamped = list(writers)
        for p in seam:
            stamped[p] = event_id
        if tuple(entry['after_writers']) != tuple(stamped):
            raise ValueError('the endpoint stamp or outside-write frame law failed')
        parents = tuple(p for p, parent in enumerate(log[:event_id])
                        if authenticated_parent_ports(parent, entry))
        # A fresh threaded execution leaves no unauthenticated earlier writer.
        cited = tuple(sorted({writers[p] for p in seam if writers[p] is not None}))
        if parents != cited:
            raise ValueError('a current writer citation lacks an authenticated value/write witness')
        height = 0 if not parents else 1+max(rows[p]['height'] for p in parents)
        ancestors = set(parents)
        for parent in parents:
            ancestors.update(rows[parent]['ancestors'])
        rows.append({'seam': seam, 'parents': parents, 'height': height,
                     'ancestors': tuple(sorted(ancestors))})
        state, writers = tuple(expected), tuple(stamped)
    return rows


def semantic_execution(initial, tape, coins):
    """Emit native branch payloads with authenticated endpoint provenance.

    The seam and coin sequence is exogenous. The declared log retains every
    invocation and its record; it does not add feedback or claim that the
    local observer can inspect hidden source values or all writer versions.
    """
    state = _loads(initial)
    if type(tape) not in (tuple, list) or type(coins) not in (tuple, list) or len(tape) != len(coins):
        raise ValueError('equal-length complete seam and coin tapes are required')
    writers, log = (None,)*12, []
    for event_id, (original, coin) in enumerate(zip(tape, coins, strict=True)):
        seam = _seam(original)
        if type(coin) is not int or coin not in (0, 1):
            raise ValueError('a binary exogenous coin is required')
        after, record = native.native(state, seam, coin)
        stamped = list(writers)
        for p in seam:
            stamped[p] = event_id
        log.append({'event_id': event_id, 'seam': list(seam), 'coin': coin,
                    'record': list(record), 'before_values': list(state), 'after_values': list(after),
                    'before_writers': list(writers), 'after_writers': stamped,
                    'read_set': list(seam), 'write_set': list(seam), 'causal_support': list(seam)})
        state, writers = after, tuple(stamped)
    if verify_semantic_execution(initial, log, tape=tape, coins=coins) != dependencies(tape):
        raise ArithmeticError('authenticated native parents differ from the endpoint DAG')
    return log


def semantic_support():
    """Finite support certificate relative to the admitted two-register read.

    Equal local reads give equal local payload/current. Each endpoint is
    necessary within this read interface, witnessed on actual total-two
    states. This does not assert a unique minimal support among global
    read interfaces related by the fixed-total constraint.
    """
    checks, necessary = 0, 0
    for seam in passage.recovered_carrier()[0]:
        for coin in (0, 1):
            local = {}
            for before in native.states():
                log = semantic_execution(before, [seam], [coin])
                after, record = log[0]['after_values'], log[0]['record']
                key = tuple(before[p] for p in seam)
                value = (tuple(after[p] for p in seam), tuple(record))
                if key in local and local[key] != value:
                    raise ArithmeticError('native accepted payload depends on a hidden outside read')
                local[key] = value
                checks += 1
            for port in (0, 1):
                if not any(x[1-port] == y[1-port] and x[port] != y[port] and local[x] != local[y]
                           for x in local for y in local):
                    raise ArithmeticError('a declared endpoint has no local support witness')
                necessary += 1
    return {'actual_native_branch_checks': checks, 'necessary_endpoint_witnesses': necessary,
            'support_sound_and_irredundant_within_two_endpoint_read_interface': True,
            'unique_minimal_support_among_global_read_interfaces_claimed': False}


def lean_branch_binding(branches):
    """Compare the complete native branch table proved in the Lean witness.

    Python checks the literal executable-data boundary. Lean CI checks
    branch_table_exact and the placement theorem; this parser is not a
    substitute for kernel proof checking.
    """
    source = (Path(__file__).resolve().parents[2]/
              'Lean/Geometry/NativeRepairSemanticPlacement.lean').read_text(encoding='utf-8')
    match = re.search(r'def branchTable\s*:.*?:=\s*(\[.*?\])\n\n', source, re.S)
    if match is None:
        raise ValueError('the registered Lean literal branch table is absent')
    table = ast.literal_eval(match.group(1))
    actual = [(branch['coins'], [r[1] for r in branch['records']],
               branch['semantic_log'][-1]['after_values']) for branch in branches]
    if table != actual:
        raise ArithmeticError('the Lean native execution table and actual source branches differ')
    declared = ' '.join(source.split())
    if ('quotientEquivVec3.symm ![0, Real.goldenRatio / Real.sqrt (2+Real.goldenRatio), 0]' not in declared
            or 'def spatial (e : Event) : FrameQuotient := if e = 1 then -midpoint else midpoint' not in declared
            or 'def timeScale : ℝ := Real.goldenRatio / Real.sqrt (2+Real.goldenRatio)' not in declared):
        raise ValueError('the registered Lean midpoint and scale readout has changed')
    phi = (1+sp.sqrt(5))/2
    raw = slow_geometry.raw_frame()
    q5 = lambda x: sp.Rational(x[0])+sp.Rational(x[1])*sp.sqrt(5)
    mid = [sum(q5(raw[p][d]) for p in (0, 1))/2 for d in range(3)]
    opposite = [sum(q5(raw[p][d]) for p in (2, 3))/2 for d in range(3)]
    if (any(sp.simplify(x-y) != 0 for x, y in zip(mid, [0, phi, 0], strict=True))
            or any(sp.simplify(x+y) != 0 for x, y in zip(mid, opposite, strict=True))
            or sp.simplify(phi**2/(2+phi)-SPEED_SQUARED) != 0):
        raise ArithmeticError('the formal midpoint chart differs from the actual source Gram placement')
    return {'branches_compared': len(actual),
            'exact_table_theorem': 'OPH.NativeRepairSemanticPlacement.branch_table_exact',
            'placement_theorem': 'OPH.NativeRepairSemanticPlacement.native_generatedBeforeEq_sourceCausalLE',
            'python_table_comparison_replaces_Lean_kernel_check': False,
            'source_midpoint_chart_and_squared_time_scale_checked': True,
            'formal_instance_scope': 'the declared four-event tape, every one of its sixteen coin branches'}


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


def _verify_diamond_witness(tape, initial, *, require_null_edges):
    """Verify a four-event fork-and-join diamond and its 16-branch causal order."""
    rows = dependencies(tape)
    expected_parents = [(), (0,), (0,), (1, 2)]
    expected_heights = [0, 1, 1, 2]
    expected_ancestors = [(), (0,), (0,), (0, 1, 2)]
    if ([r['parents'] for r in rows] != expected_parents
            or [r['height'] for r in rows] != expected_heights
            or [r['ancestors'] for r in rows] != expected_ancestors):
        raise ArithmeticError('the four-event tape is not a fork-and-join causal diamond')

    diamond_edges = {(0, 1), (0, 2), (1, 3), (2, 3)}
    for a, b in sorted(diamond_edges):
        edge_interval = sp.simplify(interval_squared(rows[a], rows[b]))
        if require_null_edges and edge_interval != 0:
            raise ArithmeticError('a parent edge of the null diamond is not lightlike')
        if edge_interval.is_nonnegative is not True:
            raise ArithmeticError('a parent edge of the diamond left the forward cone')

    timelike_03 = sp.expand(interval_squared(rows[0], rows[3]))
    spacelike_12 = sp.expand(interval_squared(rows[1], rows[2]))
    if timelike_03.is_positive is not True or spacelike_12.is_negative is not True:
        raise ArithmeticError('diamond apex (0,3) or antichain (1,2) lost its causal sign')

    def generated_before_eq(a, b):
        return a == b or a in rows[b]['ancestors']

    def source_causal_le(a, b):
        return bool(
            sp.simplify(interval_squared(rows[a], rows[b])).is_nonnegative is True
            and rows[b]['height'] >= rows[a]['height'])

    order_pairs = []
    intervals = {}
    for a in range(4):
        for b in range(4):
            gen_le = generated_before_eq(a, b)
            metric_le = source_causal_le(a, b)
            if gen_le != metric_le:
                raise ArithmeticError('generated causal order and source metric cone order differ')
            gen_interval = [
                c for c in range(4)
                if generated_before_eq(a, c) and generated_before_eq(c, b)]
            metric_interval = [
                c for c in range(4)
                if source_causal_le(a, c) and source_causal_le(c, b)]
            if gen_interval != metric_interval:
                raise ArithmeticError('generated causal interval and metric diamond differ')
            intervals[f'{a}->{b}'] = gen_interval
            order_pairs.append({
                'pair': [a, b],
                'generated_before_eq': gen_le,
                'source_causal_le': metric_le,
                'interval_squared': str(sp.expand(interval_squared(rows[a], rows[b]))),
            })

    if (intervals['0->3'] != [0, 1, 2, 3]
            or intervals['1->2'] != [] or intervals['2->1'] != []
            or intervals['0->1'] != [0, 1] or intervals['0->2'] != [0, 2]
            or intervals['1->3'] != [1, 3] or intervals['2->3'] != [2, 3]):
        raise ArithmeticError('diamond causal intervals do not match fork-and-join structure')

    branch_summaries = []
    zero_current_events = 0
    for coins in product((0, 1), repeat=4):
        log = semantic_execution(initial, tape, coins)
        verified = verify_semantic_execution(initial, log, tape=tape, coins=coins)
        if verified != rows:
            raise ArithmeticError('authenticated semantic log differs from diamond dependencies')
        for a in range(4):
            for b in range(a+1, 4):
                ports = authenticated_parent_ports(log[a], log[b])
                if bool(ports) != ((a, b) in diamond_edges):
                    raise ArithmeticError('authenticated parent ports differ from diamond edges')
        for entry in log:
            if entry['record'][1] == 0:
                zero_current_events += 1
        branch_summaries.append({
            'coins': list(coins),
            'records': [entry['record'] for entry in log],
            'terminal_occupied_ports': [
                i for i, value in enumerate(log[-1]['after_values']) if value],
        })

    return {
        'seam_tape': [list(s) for s in tape],
        'parents': [list(r['parents']) for r in rows],
        'heights': [r['height'] for r in rows],
        'ancestors': [list(r['ancestors']) for r in rows],
        'parent_edge_intervals_squared': {
            f'{a}->{b}': str(sp.expand(interval_squared(rows[a], rows[b])))
            for a, b in sorted(diamond_edges)
        },
        'apex_0_to_3_interval_squared': str(timelike_03),
        'antichain_1_to_2_interval_squared': str(spacelike_12),
        'causal_intervals': intervals,
        'order_equivalence_pairs_checked': len(order_pairs),
        'coin_branches_authenticated': len(branch_summaries),
        'zero_current_commits_retained': zero_current_events,
        'all_coin_branch_summaries': branch_summaries,
    }


def faithful_diamond_placement():
    """Verify two-way causal order and diamond intervals on the 12-port carrier."""
    seams, _ = passage.recovered_carrier()
    initial = tuple(int(i in (0, 2)) for i in range(12))

    # Find the lexicographically first 4-event fork-and-join diamond where all
    # four parent edges are exact null vectors (distance_squared == SPEED_SQUARED).
    null_diamond = None
    for e0 in seams:
        for e1 in seams:
            if not (set(e0) & set(e1)) or sp.expand(distance_squared(e0, e1)-SPEED_SQUARED) != 0:
                continue
            for e2 in seams:
                if (e2 <= e1 or not (set(e0) & set(e2))
                        or sp.expand(distance_squared(e0, e2)-SPEED_SQUARED) != 0):
                    continue
                for e3 in seams:
                    if (sp.expand(distance_squared(e1, e3)-SPEED_SQUARED) != 0
                            or sp.expand(distance_squared(e2, e3)-SPEED_SQUARED) != 0):
                        continue
                    candidate = (e0, e1, e2, e3)
                    rows = dependencies(candidate)
                    if ([r['parents'] for r in rows] == [(), (0,), (0,), (1, 2)]
                            and [r['height'] for r in rows] == [0, 1, 1, 2]):
                        null_diamond = candidate
                        break
                if null_diamond is not None:
                    break
            if null_diamond is not None:
                break
        if null_diamond is not None:
            break

    if null_diamond != ((0, 1), (0, 10), (1, 8), (0, 1)):
        raise ArithmeticError('the lexicographically first null causal diamond changed')

    # Also find the lexicographically first 4-distinct-seam 4-cycle diamond
    # with disjoint apex seams e0 and e3.
    cycle_diamond = None
    for e0 in seams:
        u, v = e0
        for e1 in seams:
            if e1 == e0 or u not in e1 or v in e1:
                continue
            w1 = e1[0] if e1[1] == u else e1[1]
            for e2 in seams:
                if e2 in (e0, e1) or v not in e2 or u in e2 or w1 in e2:
                    continue
                w2 = e2[0] if e2[1] == v else e2[1]
                e3 = (min(w1, w2), max(w1, w2))
                if e3 in seams and set(e0).isdisjoint(e3):
                    cycle_diamond = (e0, e1, e2, e3)
                    break
            if cycle_diamond is not None:
                break
        if cycle_diamond is not None:
            break

    if cycle_diamond != ((0, 1), (0, 5), (1, 9), (5, 9)):
        raise ArithmeticError('the lexicographically first four-seam cycle diamond changed')

    primary = _verify_diamond_witness(null_diamond, initial, require_null_edges=True)
    distinct_cycle = _verify_diamond_witness(cycle_diamond, initial, require_null_edges=False)

    return {
        **primary,
        'all_parent_edges_are_null': True,
        'distinct_seam_cycle_diamond': {
            'seam_tape': distinct_cycle['seam_tape'],
            'parents': distinct_cycle['parents'],
            'heights': distinct_cycle['heights'],
            'ancestors': distinct_cycle['ancestors'],
            'parent_edge_intervals_squared': distinct_cycle['parent_edge_intervals_squared'],
            'apex_0_to_3_interval_squared': distinct_cycle['apex_0_to_3_interval_squared'],
            'antichain_1_to_2_interval_squared': distinct_cycle['antichain_1_to_2_interval_squared'],
            'order_equivalence_pairs_checked': distinct_cycle['order_equivalence_pairs_checked'],
            'coin_branches_authenticated': distinct_cycle['coin_branches_authenticated'],
        },
        'two_way_causal_order_equivalence_on_witness_checked': True,
        'all_sixteen_causal_intervals_match_metric_diamonds': True,
        'all_sixteen_coin_branches_authenticated': True,
        'arbitrary_native_tape_satisfies_converse_cone_implication': False,
    }


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
        log = semantic_execution(before, tape, coins)
        branches.append({'coins': list(coins), 'records': [entry['record'] for entry in log],
                         'terminal_occupied_ports': [i for i, value in enumerate(log[-1]['after_values']) if value],
                         'semantic_log': log})
    return {'contract': {
                'spatial_readout': 'declared seam midpoint in the native response Gram',
                'time_parameter': 'canonical endpoint-dependency height times sqrt(c_squared)',
                'complete_source_seam_tape_including_zero_current_events_required': True,
                'predeclared_exogenous_schedule_required_for_endpoint_only_dependency_graph': True,
                'arbitrary_controller_or_preparation_dependency_edges_covered': False,
                'endpoint_parents_are_structural_not_proven_signal_influences': True,
                'complete_value_version_log_is_execution_instrumentation': True,
                'intrinsic_primitive_observer_receives_complete_value_version_log': False,
                'classical_branch_logs_identify_arbitrary_coherent_quantum_histories': False,
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
            'faithful_diamond_control': faithful_diamond_placement(),
            'finite_current_readout': finite_readout(),
            'semantic_execution_support': semantic_support(),
            'Lean_native_branch_binding': lean_branch_binding(branches),
            'acceptance': {'forward_placement_geometric_hypotheses_checked': True,
                'authenticated_native_semantic_log_emitted': True,
                'kernel_checked_Lean_native_placement_instance': True,
                'faithful_diamond_two_way_causal_order_witness_checked': True,
                'cone_implies_ancestry': False, 'physical_Lorentz_symmetry_or_clock_derived': False,
                'continuum_spacetime_derived': False}, 'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--output', type=Path)
    modes.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native causal placement, diamond order witness and Lean branch controls verified.')
    elif args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
