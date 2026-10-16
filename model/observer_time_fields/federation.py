"""Complete field histories at local source opportunities on the federation.

The carrier is all 240 ports of the actual level-zero federation, with two
particles globally. Only the readout region comprises two adjacent carriers.
Remote native operations continue while the observer waits for its next
incident attempt. Initial values of the two read ports are disclosed.

The executor acts on diagonal native states. Integer downward rounding gives
componentwise lower bounds for every complete terminal density, including
timeouts. No floating probability, hidden normalization, dense resolvent or
independent-particle replacement enters these bounds.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
import math

import numpy as np
from scipy import sparse

from ..federated_observer import source as native_source, scaling


SOURCE_PATHS = tuple(sorted(set((*native_source.SOURCE_PATHS, *scaling.SOURCE_PATHS,
    'model/observer_time_fields/federation.py'))))
MAX_CAP = 512
MAX_TICKS = 4
DEFAULT_CAP = 256
DEFAULT_BITS = 51
INITIAL_LABELS = tuple(product((0, 1), repeat=2))


@dataclass(frozen=True)
class Source:
    ports: int
    cells: int
    seams: tuple
    region: tuple
    incident: tuple
    remote: tuple
    read_ports: tuple
    states: tuple
    denominator: int


@lru_cache(maxsize=1)
def _source():
    graph = native_source.federation(0)
    seam = graph.inter_seams[0]
    cells = {port//12 for port in seam}
    region = tuple(port for port in range(graph.ports) if port//12 in cells)
    incident = tuple(e for e in graph.seams if any(port in region for port in e))
    remote = tuple(e for e in graph.seams if e not in incident)
    if (graph.ports != 240 or graph.cells != 20 or len(graph.seams) != 630
            or len(region) != 24 or len(incident) != 65 or len(remote) != 565
            or any(neighbour not in region for port in seam for neighbour in graph.neighbours[port])):
        raise ArithmeticError('the canonical complete federation/readout interface changed')
    return Source(graph.ports, graph.cells, tuple(graph.seams), region, incident, remote,
                  tuple(seam), tuple(combinations(range(graph.ports), 2)), 2*len(graph.seams))


def source():
    """Fresh envelope; all nested source data are immutable integer tuples."""
    return replace(_source())


def _index(a, b, ports=240):
    if a > b:
        a, b = b, a
    return a*(2*ports-a-1)//2+b-a-1


@lru_cache(maxsize=1)
def _operators():
    data = _source()
    neighbours = [[] for _ in range(data.ports)]
    for a, b in data.seams:
        neighbours[a].append(b)
        neighbours[b].append(a)
    incident = set(data.incident)
    rows, cols, entries = [[], []], [[], []], [[], []]
    # The established actual endpoint rule has fair average (I+Swap)/2.
    # Verify that premise from the source, before assembling its sector lift.
    for loads in product((0, 1), repeat=2):
        actual = [native_source.native.native(loads, (0, 1), coin)[0] for coin in (0, 1)]
        if sorted(actual) != sorted((loads, loads[::-1])):
            raise ArithmeticError('the source fair endpoint rule ceased to be identity/swap')
    for column, (a, b) in enumerate(data.states):
        degrees = [0, 0]
        for occupied, other in ((a, b), (b, a)):
            for target in neighbours[occupied]:
                if target == other:
                    continue
                flag = int(tuple(sorted((occupied, target))) in incident)
                rows[flag].append(_index(target, other, data.ports))
                cols[flag].append(column)
                entries[flag].append(1)
                degrees[flag] += 1
        for flag, edges in ((0, data.remote), (1, data.incident)):
            rows[flag].append(column)
            cols[flag].append(column)
            entries[flag].append(2*len(edges)-degrees[flag])
    size = len(data.states)
    maps = [sparse.coo_matrix((entries[z], (rows[z], cols[z])), shape=(size, size),
                             dtype=np.int64).tocsr() for z in (0, 1)]
    repair_rows, repair_cols = [], []
    for column, occupied in enumerate(data.states):
        for coin in (0, 1):
            target, _ = native_source.binary_step(occupied, data.read_ports, coin)
            repair_rows.append(_index(*target, data.ports))
            repair_cols.append(column)
    maps.append(sparse.coo_matrix((np.ones(len(repair_rows), dtype=np.int64),
                                  (repair_rows, repair_cols)), shape=(size, size)).tocsr())
    for matrix, weight in zip(maps, (2*len(data.remote), 2*len(data.incident), 2), strict=True):
        if (matrix.dtype != np.int64 or np.any(matrix.data < 0)
                or np.any(np.asarray(matrix.sum(axis=0)).ravel() != weight)
                or np.any(np.asarray(matrix.sum(axis=1)).ravel() != weight)
                or (matrix-matrix.T).nnz):
            raise ArithmeticError('the complete integer native sector operator changed')
        matrix.data.flags.writeable = False
        matrix.indices.flags.writeable = False
        matrix.indptr.flags.writeable = False
    return tuple(maps)


def operators():
    """Detached integer CSR matrices (remote, incident, feedback repair).

    Columns act on populations. The first two divide by 1260 and have column
    sums 1130 and 130; the actual fair feedback repair divides by two.
    """
    return tuple(matrix.copy() for matrix in _operators())


def _parameters(cap, ticks, grid_bits, feedback, initial):
    if type(cap) is not int or not 1 <= cap <= MAX_CAP:
        raise ValueError('an integer per-interval cap from one through 512 is required')
    if type(ticks) is not int or not 1 <= ticks <= MAX_TICKS:
        raise ValueError('one through four complete scheduled local tick reads are required')
    if type(grid_bits) is not int or not 16 <= grid_bits <= DEFAULT_BITS:
        raise ValueError('an integer directed grid precision from 16 through 51 is required')
    if type(feedback) is not bool:
        raise ValueError('an explicit record-controlled repair Boolean is required')
    data = _source()
    if (initial is not None and (type(initial) is not tuple or len(initial) != 2
            or any(type(port) is not int for port in initial)
            or not 0 <= initial[0] < initial[1] < data.ports)):
        raise ValueError('a canonical occupied pair or the stationary fixed-two preparation is required')
    scale = 1 << grid_bits
    # Every product row has sum at most 2M, and every nonnegative input entry
    # is at most scale. Thus every int64 multiply/accumulate is below 2^62.
    if data.denominator*scale > 1 << 62:
        raise ValueError('the declared integer grid would exceed its checked arithmetic bound')
    return scale


def records(ticks=MAX_TICKS):
    if type(ticks) is not int or not 1 <= ticks <= MAX_TICKS:
        raise ValueError('one through four complete scheduled local tick reads are required')
    timeouts = tuple((label, word, 'timeout') for length in range(ticks)
                     for label in INITIAL_LABELS for word in product((0, 1), repeat=length))
    complete = tuple((label, word, 'complete') for label in INITIAL_LABELS
                     for word in product((0, 1), repeat=ticks))
    return timeouts+complete


def _floor_apply(matrix, masses, denominator):
    """Internal source-validated nonnegative integer action with downward division.

    _operators validates all matrix coefficients and row/column sums;
    _parameters bounds multiplication. Mass conservation and downward rounding
    preserve the total-mass bound scale at every call, including feedback.
    """
    return (matrix@masses)//denominator


def _immutable_int(array):
    return np.frombuffer(array.tobytes(), dtype=np.int64).reshape(array.shape)


@dataclass(frozen=True)
class Execution:
    cap: int
    ticks: int
    grid_bits: int
    feedback: bool
    initial: tuple | None
    records: tuple
    counts: np.ndarray
    seal: str


def _seal(result):
    header = (result.cap, result.ticks, result.grid_bits, result.feedback,
              result.initial, result.records)
    return hashlib.sha256(json.dumps(header, separators=(',', ':')).encode()
                          +result.counts.tobytes()).hexdigest()


@lru_cache(maxsize=4)
def _compute(*, cap, ticks, grid_bits, feedback, initial):
    """Execute every initial-labelled field word and unresolved timeout block.

    A read mask is the proven current-counter value, not an additional
    occupation measurement. After the first read, bit one requests one actual
    fair repair of the canonical inter-carrier seam; bit zero does nothing.
    That extra source operation is charged and also increments the incident
    opportunity clock. Subsequent uniform-scheduler waits retain remote moves.
    """
    scale = _parameters(cap, ticks, grid_bits, feedback, initial)
    data = _source()
    remote, incident, repair = _operators()
    states = np.asarray(data.states)
    occupied = tuple((states[:, 0] == port) | (states[:, 1] == port) for port in data.read_ports)
    prior = np.zeros((len(states), 4), dtype=np.int64)
    label_indices = 2*occupied[0].astype(np.int64)+occupied[1].astype(np.int64)
    if initial is None:
        prior[np.arange(len(states)), label_indices] = scale//len(states)
    else:
        index = _index(*initial, data.ports)
        prior[index, label_indices[index]] = scale
    labels = tuple((label, ()) for label in INITIAL_LABELS)
    blocks, menu = [], []
    for tick in range(ticks):
        waiting, stopped = prior, np.zeros_like(prior)
        for _ in range(cap):
            stopped += _floor_apply(incident, waiting, data.denominator)
            waiting = _floor_apply(remote, waiting, data.denominator)
        blocks.append(waiting)
        menu.extend((label, word, 'timeout') for label, word in labels)
        following = np.empty((len(states), 2*len(labels)), dtype=np.int64)
        mask = occupied[tick % 2][:, None]
        following[:, 0::2] = stopped*(~mask)
        following[:, 1::2] = stopped*mask
        labels = tuple((label, (*word, bit)) for label, word in labels for bit in (0, 1))
        if feedback and tick == 0:
            selected = [i for i, (_, word) in enumerate(labels) if word[0] == 1]
            following[:, selected] = _floor_apply(repair, following[:, selected], 2)
        prior = following
    blocks.append(prior)
    menu.extend((label, word, 'complete') for label, word in labels)
    if tuple(menu) != records(ticks):
        raise ArithmeticError('the complete initial-labelled field/timeout menu changed')
    counts = _immutable_int(np.concatenate(blocks, axis=1))
    if np.any(counts < 0) or sum(map(int, counts.ravel())) > scale:
        raise ArithmeticError('directed native propagation exceeded its nonnegative mass budget')
    result = Execution(cap, ticks, grid_bits, feedback, initial, tuple(menu), counts, '')
    return replace(result, seal=_seal(result))


def execute(*, cap=DEFAULT_CAP, ticks=MAX_TICKS, grid_bits=DEFAULT_BITS,
            feedback=False, initial=None):
    """Execute the complete source family and detach its original integer output.

    Exact parameter validation precedes caching. Cached envelopes and ndarray
    objects never escape: even immutable array data have mutable shape metadata.
    """
    _parameters(cap, ticks, grid_bits, feedback, initial)
    result = _compute(cap=cap, ticks=ticks, grid_bits=grid_bits, feedback=feedback, initial=initial)
    return replace(result, counts=_immutable_int(result.counts))


def validate_execution(result):
    """Detach and compare a complete output with the actual source replay.

    An integrity hash alone is insufficient. The expected digest comes only
    from _compute for the exact validated parameters; a cache miss recomputes
    the source execution. Caller-authored or resealed wrong laws are rejected.
    """
    if type(result) is not Execution:
        raise ValueError('a typed computed federation field execution is required')
    if type(result.counts) is not np.ndarray or result.counts.dtype != np.int64:
        raise ValueError('an original integer array is required')
    # A caller may supply a writable base-array alias. Detach before checking
    # its values or seal, then inspect exactly the snapshot we will return.
    # Otherwise a mutation between hashing and copying could evade replay.
    result = replace(result, counts=_immutable_int(result.counts))
    scale = _parameters(result.cap, result.ticks, result.grid_bits, result.feedback, result.initial)
    if (type(result.records) is not tuple or result.records != records(result.ticks)
            or type(result.counts) is not np.ndarray or result.counts.dtype != np.int64
            or result.counts.shape != (len(_source().states), len(result.records))
            or np.any(result.counts < 0) or sum(map(int, result.counts.ravel())) > scale
            or type(result.seal) is not str or result.seal != _seal(result)):
        raise ValueError('the complete original integer output or its execution binding changed')
    expected = _compute(cap=result.cap, ticks=result.ticks, grid_bits=result.grid_bits,
                        feedback=result.feedback, initial=result.initial)
    if result.seal != expected.seal:
        raise ValueError('the claimed complete integer law differs from source replay')
    return result


def budget(result):
    """Exact total-output error from componentwise domination and complete mass.

    The ideal finite output has mass one. Every stored diagonal entry is a
    lower bound for its own complete block, so the total nuclear discrepancy
    is exactly the missing stored mass. The eventual comparison additionally
    allows the exact finite timeout population, without postselection.
    """
    result = validate_execution(result)
    scale = 1 << result.grid_bits
    masses = tuple(F(int(value), scale) for value in result.counts.sum(axis=0))
    defect = 1-sum(masses, F(0))
    data = _source()
    survival = F(len(data.remote), len(data.seams))**result.cap
    timeout = 1-(1-survival)**result.ticks
    success = sum((mass for record, mass in zip(result.records, masses, strict=True)
                   if record[2] == 'complete'), F(0))
    return {'per_record_mass_lower': masses, 'complete_l1_error': defect,
            'finite_half_nuclear_error': defect/2,
            'exact_timeout_probability': timeout,
            'eventual_half_nuclear_error': defect/2+timeout,
            'eventual_field_missing_mass': 1-success}


def _partitions(items):
    if not items:
        return ((),)
    first, *rest = items
    result = []
    for partition in _partitions(tuple(rest)):
        result.append(((first,), *partition))
        for index, block in enumerate(partition):
            result.append((*partition[:index], (first, *block), *partition[index+1:]))
    return tuple(result)


def field_statistics(result):
    """Directed eventual moments and cumulant from all complete field words."""
    proof = budget(result)
    probabilities = {word: F(0) for word in product((0, 1), repeat=result.ticks)}
    for (_, word, status), mass in zip(result.records, proof['per_record_mass_lower'], strict=True):
        if status == 'complete':
            probabilities[word] += mass
    missing = proof['eventual_field_missing_mass']
    moments = {}
    for size in range(1, result.ticks+1):
        for indices in combinations(range(result.ticks), size):
            lower = sum((mass for word, mass in probabilities.items()
                         if all(word[i] for i in indices)), F(0))
            moments[indices] = (lower, min(F(1), lower+missing))
    cumulant = None
    if result.ticks == 4:
        lower = upper = F(0)
        for partition in _partitions(tuple(range(4))):
            coefficient = (-1)**(len(partition)-1)*math.factorial(len(partition)-1)
            lo = math.prod(moments[tuple(sorted(block))][0] for block in partition)
            hi = math.prod(moments[tuple(sorted(block))][1] for block in partition)
            lower += coefficient*(lo if coefficient > 0 else hi)
            upper += coefficient*(hi if coefficient > 0 else lo)
        cumulant = (lower, upper)
    return {'field_word_mass_lower': probabilities, 'moments': moments,
            'eventual_field_missing_mass': missing, 'fourth_cumulant_interval': cumulant}


def read_clock_values(word, *, feedback=False):
    if (type(word) is not tuple or not word or len(word) > MAX_TICKS
            or any(type(bit) is not int or bit not in (0, 1) for bit in word)
            or type(feedback) is not bool):
        raise ValueError('a nonempty bounded exact binary field word and Boolean policy required')
    return tuple(index+1+int(feedback and index > 0 and word[0]) for index in range(len(word)))


def resources(result):
    result = validate_execution(result)
    data = _source()
    q = F(len(data.incident), len(data.seams))
    hit = 1-(1-q)**result.cap
    expected_uniform = sum((hit**n/q for n in range(1, result.ticks+1)), F(0))
    if result.initial is None:
        first_one = F(2, data.ports)
    else:
        port = data.read_ports[0]
        neighbours = [b if a == port else a for a, b in data.incident if port in (a, b)]
        first_one = F(int(port in result.initial))
        first_one -= F(len(neighbours)*int(port in result.initial)
                       -sum(other in result.initial for other in neighbours), 2*len(data.incident))
    # The first read port and every neighbour lie inside R, hence no remote
    # waiting operation changes the loads entering this first conditional mean.
    expected_repair = hit*first_one if result.feedback else F(0)
    return {'native_preparations': 1, 'disclosed_initial_field_bits': 2,
            'maximum_uniform_scheduler_attempts': result.cap*result.ticks,
            'maximum_feedback_native_repairs': int(result.feedback),
            'maximum_total_native_attempts': result.cap*result.ticks+int(result.feedback),
            'expected_uniform_scheduler_attempts_exact': str(expected_uniform),
            'expected_feedback_native_repairs_exact': str(expected_repair),
            'expected_total_native_attempts_exact': str(expected_uniform+expected_repair),
            'maximum_scheduled_incident_reads': result.ticks,
            'maximum_total_incident_opportunities': result.ticks+int(result.feedback),
            'complete_output_alphabet_size': len(result.records),
            'fixed_terminal_record_bits': (len(result.records)-1).bit_length(),
            'per_interval_attempt_counter_bits': result.cap.bit_length(),
            'retained_native_diagonal_entries': len(data.states)*len(result.records),
            'retained_terminal_integer_bytes': 8*len(data.states)*len(result.records),
            'source_sparse_integer_nonzeros': sum(int(matrix.nnz) for matrix in _operators()),
            'feedback_repair_also_increments_incident_clock': True,
            'source_pulse_compilation_or_physical_clock_calibration': False}


def readback(initial_label, history):
    """Two visible field counters receive only disclosed values and incident records."""
    data = _source()
    if (type(initial_label) is not tuple or len(initial_label) != 2
            or any(type(bit) is not int or bit not in (0, 1) for bit in initial_label)
            or type(history) is not tuple or len(history) > MAX_CAP*MAX_TICKS+1):
        raise ValueError('two disclosed binary fields and a bounded immutable incident record stream required')
    values, previous = list(initial_label), 0
    for record in history:
        if (type(record) is not tuple or len(record) != 5
                or any(type(value) is not int for value in record)):
            raise ValueError('typed event/seam/coin/current records are required')
        event, a, b, coin, current = record
        if (event <= previous or (a, b) not in data.incident
                or (coin, current) not in ((0, 0), (0, 1), (1, -1), (1, 0))):
            raise ValueError('an authentic ordered incident native record is required')
        if (a, b) == data.read_ports:
            _, (_, expected) = native_source.native.native(tuple(values), (0, 1), coin)
            if current != expected:
                raise ValueError('the disclosed endpoint fields contradict this native record')
        for index, port in enumerate(data.read_ports):
            values[index] += current*(int(port == b)-int(port == a))
        if any(value not in (0, 1) for value in values):
            raise ValueError('the current stream contradicts the disclosed binary fields')
        previous = event
    return tuple(values)


def readback_audit():
    """Lift actual source endpoint records to every state/seam/coin address.

    Source evaluation on four endpoint patterns supplies the exact local
    action; vectorized sector lifting checks both tracked fields on all
    36,136,800 global addresses. This is the same locality decomposition used
    by native_source.binary_step, not a hidden-state input to the controller.
    """
    data = _source()
    states = np.asarray(data.states)
    before = tuple((states[:, 0] == p).astype(np.int64)+(states[:, 1] == p)
                   for p in data.read_ports)
    table = {}
    for loads in product((0, 1), repeat=2):
        for coin in (0, 1):
            table[*loads, coin] = native_source.native.native(loads, (0, 1), coin)
    for seam in data.seams:
        a, b = seam
        at_a = ((states[:, 0] == a) | (states[:, 1] == a)).astype(np.int64)
        at_b = ((states[:, 0] == b) | (states[:, 1] == b)).astype(np.int64)
        for coin in (0, 1):
            current = np.array([table[x, y, coin][1][1] for x, y in product((0, 1), repeat=2)],
                               dtype=np.int64)[2*at_a+at_b]
            for index, port in enumerate(data.read_ports):
                if port in seam:
                    target = np.array([table[x, y, coin][0][int(port == b)]
                                       for x, y in product((0, 1), repeat=2)],
                                      dtype=np.int64)[2*at_a+at_b]
                else:
                    target = before[index]
                reconstructed = before[index]+current*(int(port == b)-int(port == a))
                if not np.array_equal(reconstructed, target):
                    raise ArithmeticError('current-only field integration disagrees with the native source')
    return {'source_endpoint_branches': 8,
            'all_global_state_seam_coin_addresses': len(states)*len(data.seams)*2,
            'feedback_addresses_already_included': len(states)*2,
            'tracked_fields': list(data.read_ports),
            'current_readback_identity_holds': True,
            'initial_disclosure_and_authentic_complete_incident_stream_required': True}


def transfer_certificate():
    """Verify the actual triangle-cover premises of the full-sector floor."""
    data = _source()
    scaling.triangle_gram()
    seam_set, incident = set(data.seams), set(data.incident)
    cells = tuple(sorted({port//12 for port in data.region}))
    certificates = []
    for cell in cells:
        ports = tuple(range(12*cell, 12*(cell+1)))
        internal = {edge for edge in seam_set if edge[0] in ports and edge[1] in ports}
        faces = tuple(face for face in combinations(ports, 3)
                      if set(combinations(face, 2)) <= internal)
        multiplicity = Counter(edge for face in faces for edge in combinations(face, 2))
        if (len(internal) != 30 or len(faces) != 20 or set(multiplicity) != internal
                or set(multiplicity.values()) != {2} or not internal <= incident):
            raise ArithmeticError('an actual incident carrier lacks the required triangle cover')
        certificates.append({'cell': cell, 'internal_seams': len(internal),
                             'triangles': len(faces), 'edge_triangle_multiplicity': 2})
    lower = F(15*len(cells), len(data.seams))
    if lower != F(1, 21):
        raise ArithmeticError('the full-sector incident transfer floor changed')
    return {'actual_carrier_triangle_covers': certificates,
            'exact_eight_state_triangle_gram_identity_verified': True,
            'B_lower': str(lower), 'symmetric_T_lower': str(lower),
            'symmetric_T_upper': '1', 'all_binary_number_sectors': True,
            'local_field_multipliers_commute_with_remote_A': True}


def _summary_execution(result):
    proof = budget(result)
    fields = field_statistics(result)
    data = _source()
    return {'ports': data.ports, 'elementary_carriers': data.cells,
            'global_particle_number': 2, 'complete_native_state_dimension': len(data.states),
            'region_ports': list(data.region), 'read_ports': list(data.read_ports),
            'read_port_word': [data.read_ports[i % 2] for i in range(result.ticks)],
            'incident_seams': len(data.incident), 'remote_seams': len(data.remote),
            'initial_preparation': 'uniform fixed-two with disclosed two-field label'
                if result.initial is None else {'occupied': list(result.initial)},
            'cap_per_scheduled_interval': result.cap, 'scheduled_reads': result.ticks,
            'feedback': result.feedback, 'grid_bits': result.grid_bits,
            'int64_multiply_accumulate_upper': str(data.denominator*(1 << result.grid_bits)),
            'complete_output_count': len(result.records),
            'complete_l1_rounding_error': str(proof['complete_l1_error']),
            'finite_half_nuclear_error': str(proof['finite_half_nuclear_error']),
            'exact_timeout_probability': str(proof['exact_timeout_probability']),
            'eventual_half_nuclear_error': str(proof['eventual_half_nuclear_error']),
            'eventual_field_missing_mass': str(proof['eventual_field_missing_mass']),
            'fourth_cumulant_interval': None if fields['fourth_cumulant_interval'] is None
                else list(map(str, fields['fourth_cumulant_interval'])),
            'field_means': [list(map(str, fields['moments'][i,])) for i in range(result.ticks)],
            'complete_outcome_mass_lower': [
                {'initial_label': list(label), 'field_word': list(word), 'status': status,
                 'mass_lower': str(mass)}
                for (label, word, status), mass in zip(result.records, proof['per_record_mass_lower'], strict=True)],
            'resources': resources(result),
            'scope': {'diagonal_native_sector': True, 'all_remote_native_moves_retained': True,
                      'independent_particle_approximation': False,
                      'initial_field_label_retained': True,
                      'native_source_channels_not_compiled_physical_pulses': True,
                      'field_masks_use_proved_current_readback': True,
                      'calibrated_physical_or_modular_clock_execution': False}}


def summary(result=None):
    """Portable summaries, or the complete baseline/feedback production packet."""
    if result is not None:
        return _summary_execution(result)
    source_proof = transfer_certificate()
    readback_proof = readback_audit()
    baseline = execute()
    intervention = execute(feedback=True)
    baseline_fields, feedback_fields = field_statistics(baseline), field_statistics(intervention)
    baseline_proof, feedback_proof = budget(baseline), budget(intervention)
    if max(baseline_proof['complete_l1_error'], feedback_proof['complete_l1_error']) > F(1, 10**6):
        raise ArithmeticError('the complete directed federation output budget is unresolved')
    cumulant = baseline_fields['fourth_cumulant_interval']
    if cumulant is None or cumulant[1] >= 0:
        raise ArithmeticError('the connected four-field witness is unresolved')
    before, after = baseline_fields['moments'][1,], feedback_fields['moments'][1,]
    difference = (after[0]-before[1], after[1]-before[0])
    if difference[0] <= 0:
        raise ArithmeticError('the record-controlled cross-carrier repair witness is unresolved')
    return {'transfer_certificate': source_proof, 'current_readback_certificate': readback_proof,
            'baseline': _summary_execution(baseline),
            'record_controlled_repair': _summary_execution(intervention),
            'second_field_feedback_minus_baseline_interval': list(map(str, difference)),
            'nonlinear_cumulant_bound': 'exact interval arithmetic over all 15 set partitions; '
                'a coarse independent bound is 75 times the complete eventual field missing mass',
            'every_initial_label_and_timeout_retained': True,
            'all_pass': True}
