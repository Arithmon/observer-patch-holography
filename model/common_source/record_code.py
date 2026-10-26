"""Sharp six-dimensional memory for the complete recorded native seam menu.

Every logical word is the uniform superposition of one perfect matching of
the nonedge graph.  The six matchings partition all 36 nonnative pairs.
The recorded C/E branch Gram matrices are scalar, although the unread work
channel does not preserve this code.  Recovery uses the existing source
compiler and its full-register Gray permutation macro.  All 512 rows,
raw abort records and arbitrary reference columns are retained.

This is active encoded memory under declared coherent control services,
not a physical field identification or a passive native invariant sector.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import math

import numpy as np

from ..observer import native
from ..observer.source_discrimination import circuit, experiment as recorded
from ..observer.source_scalar import compiler
from ..repair_scheduler import coins
from . import native_routing as routing, slow_geometry


MATCHINGS = (
    ((0, 2), (1, 3), (4, 6), (5, 7), (8, 10), (9, 11)),
    ((0, 3), (1, 2), (4, 7), (5, 6), (8, 11), (9, 10)),
    ((0, 4), (1, 6), (2, 8), (3, 11), (5, 10), (7, 9)),
    ((0, 6), (1, 4), (2, 9), (3, 10), (5, 8), (7, 11)),
    ((0, 8), (1, 11), (2, 5), (3, 7), (4, 10), (6, 9)),
    ((0, 9), (1, 10), (2, 7), (3, 5), (4, 8), (6, 11)),
)
SEAMS = tuple(slow_geometry.source().seams())
SUPPORTS = tuple(tuple(tuple(int(p in pair) for p in range(12))
                       for pair in matching) for matching in MATCHINGS)
DIMENSION = 6
SIZE = 512


def _seam(seam):
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(p) is not int for p in seam) or seam not in SEAMS):
        raise ValueError('an ordered native seam from the thirty-seam menu required')
    return seam


def _label(law, seam, coin, raw):
    law, seam = routing._candidate(law), _seam(seam)
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError('an integer native coin zero or one required')
    if type(raw) is not int or not 0 <= raw < 1 << circuit.RECORD_BITS[law]:
        raise ValueError('a raw integer record in the complete native menu required')
    return law, seam, coin, raw


def branch_images(law, seam, coin, raw):
    """Independent endpoint enumeration; all six rows have the same length."""
    law, seam, coin, raw = _label(law, seam, coin, raw)
    images = []
    a, b = seam
    for support in SUPPORTS:
        row = []
        for state in support:
            low, parity = divmod(state[a]+state[b], 2)
            left = low+coin*parity
            current = state[a]-left
            code = {0: 0, 1: 1, -1: 2}[current]
            if law == 'E':
                code += 4*parity
            if code == raw:
                final = list(state)
                final[a], final[b] = left, low+(1-coin)*parity
                row.append(circuit.encode(tuple(final)))
        images.append(tuple(row))
    if (len({len(row) for row in images}) != 1
            or len({i for row in images for i in row}) != sum(map(len, images))):
        raise ArithmeticError('native matching code lost its exact scalar Gram identity')
    return tuple(images)


def branch_weight(law, seam, coin, raw):
    """Exact outcome probability including fair coin and six code supports."""
    return F(len(branch_images(law, seam, coin, raw)[0]), 12)


def _encoder_permutation():
    return routing.complete_permutation(tuple(
        (logical+8*term, circuit.encode(state))
        for logical, support in enumerate(SUPPORTS)
        for term, state in enumerate(support)))


def _branch_permutation(law, seam, coin, raw):
    return routing.complete_permutation(tuple(
        (logical+8*term, address)
        for logical, row in enumerate(branch_images(law, seam, coin, raw))
        for term, address in enumerate(row)))


def _mixing_word(count, inverse=False):
    if type(count) is not int or count not in (1, 4, 5, 6):
        raise ValueError('a declared matching-support size required')
    if type(inverse) is not bool:
        raise ValueError('inverse must be an explicit Boolean')
    vector = [1/math.sqrt(count)]*count+[0.]*(8-count)
    word = compiler.remap(compiler.compile_real_state(vector), (3, 4, 5))
    return compiler.inverse(word) if inverse else word


@lru_cache(maxsize=2)
def _macro_charts(inverse):
    # These are executed source matrices, not ideal controlled-flip targets.
    macro, _ = routing._executed_macro(inverse)
    return tuple(routing._immutable(routing._chart(macro, target)) for target in range(9))


def _apply_permutation(images, state, inverse=False):
    """Apply the complete emitted source schedule, retaining numerical leakage.

    Reuse only the source-executed macro and its literal wire charts.  Each
    invocation is charged by permutation_resources even when numerically
    cached.  No ideal permutation matrix replaces executed amplitudes.
    """
    result = state
    charts = _macro_charts(inverse)
    for target, pattern in routing.permutation_edges(images, inverse=inverse):
        before = routing._conjugation_word(pattern)
        result = compiler.execute(before, 9, result)
        result = charts[target]@result
        result = compiler.execute(compiler.inverse(before), 9, result)
    return result


def _encode(value):
    return _apply_permutation(_encoder_permutation(),
                              compiler.execute(_mixing_word(6), 9, value))


def _finite(value):
    if not np.isfinite(value).all():
        raise ArithmeticError('matching-code execution produced nonfinite amplitudes')
    return routing._immutable(value)


@lru_cache(maxsize=1)
def code_matrix():
    """Execute the encoder on six unknown logical input columns."""
    return _finite(_encode(np.eye(SIZE, DIMENSION, dtype=complex)))


def ideal_code_matrix():
    """Independent declaration used only as a comparison, never preparation."""
    matrix = np.zeros((SIZE, DIMENSION), complex)
    for logical, support in enumerate(SUPPORTS):
        for state in support:
            matrix[circuit.encode(state), logical] = 1/math.sqrt(6)
    return routing._immutable(matrix)


def declared_recovery_matrix66(law, seam, coin, raw):
    """Algebraic recovery on exactly66 states, separate from source execution.

    First invert native movement by a complete endpoint swap.  Each of the
    six disjoint matching blocks then uses a Householder reflection taking
    its normalized surviving subset to the full uniform matching.  All30
    native-edge pair rows remain in the full matrix.  This provides a
    same-carrier extension; source execution uses the independently compiled
    full512 permutation/preparation extension in apply_recovery instead.
    """
    law, seam, coin, raw = _label(law, seam, coin, raw)
    if raw not in circuit.record_values(law, coin):
        raise ValueError('invalid raw native records halt instead of selecting recovery')
    rows = branch_images(law, seam, coin, raw)
    count = len(rows[0])
    if not count:
        return routing._immutable(np.eye(66))
    states = native.states(binary=True)
    index = {state: i for i, state in enumerate(states)}
    source_index = {circuit.encode(state): i for i, state in enumerate(states)}
    value = circuit.record_values(law, coin)[raw]
    current = value if law == 'C' else (value-(2*coin-1)*(abs(value) % 2))//2
    permutation = list(range(66))
    if current:
        a, b = seam
        for i, state in enumerate(states):
            previous = list(state)
            previous[a], previous[b] = previous[b], previous[a]
            permutation[i] = index[tuple(previous)]
    inverse_movement = np.zeros((66, 66), complex)
    inverse_movement[permutation, np.arange(66)] = 1
    reflection = np.eye(66, dtype=complex)
    for logical, outputs in enumerate(rows):
        difference = np.zeros(66, complex)
        difference[[permutation[source_index[x]] for x in outputs]] = 1/math.sqrt(count)
        difference[[index[state] for state in SUPPORTS[logical]]] -= 1/math.sqrt(6)
        norm_squared = float(np.vdot(difference, difference).real)
        if norm_squared:
            reflection -= 2*np.outer(difference, difference.conj())/norm_squared
    return routing._immutable(reflection@inverse_movement)


def encoding(state):
    """Supply a six-level logical state and six blank code bits."""
    routing._original_scalars(state)
    raw = np.asarray(state)
    if (raw.ndim not in (1, 2) or raw.shape[0] != DIMENSION
            or (raw.ndim == 2 and raw.shape[1] < 1)
            or raw.dtype.kind not in 'iufc' or not np.isfinite(raw).all()):
        raise ValueError('finite six-component logical amplitudes or reference columns required')
    value = np.array(raw, dtype=complex, copy=True)
    if not np.isfinite(value).all():
        raise ValueError('logical amplitudes must remain finite in execution precision')
    return _finite(code_matrix()@value)


def decoding(state):
    """Execute the inverse on every source row; syndrome is not selected."""
    value = routing._state(state)
    value = _apply_permutation(_encoder_permutation(), value, True)
    return _finite(compiler.execute(_mixing_word(6, True), 9, value))


def apply_recovery(law, seam, coin, raw, state):
    """Complete record-selected source unitary; invalid raw records halt."""
    law, seam, coin, raw = _label(law, seam, coin, raw)
    value = routing._state(state)
    if raw not in circuit.record_values(law, coin):
        raise ValueError('invalid raw native records halt instead of selecting recovery')
    count = len(branch_images(law, seam, coin, raw)[0])
    if not count:
        return _finite(value)
    value = _apply_permutation(_branch_permutation(law, seam, coin, raw), value, True)
    value = compiler.execute(_mixing_word(count, True), 9, value)
    return _finite(_encode(value))


def recovery_unitary(law, seam, coin, raw):
    """Source-execute all 512 input columns when a complete matrix is needed."""
    _label(law, seam, coin, raw)
    return apply_recovery(law, seam, coin, raw, np.eye(SIZE, dtype=complex))


def _route_permutation(seam):
    """Route the requested physical ports to source seam (0,1)."""
    a, b = _seam(seam)
    images = list(range(SIZE))
    for state in native.states():
        final = list(state)
        final[0], final[a] = final[a], final[0]
        final[1], final[b] = final[b], final[1]
        images[circuit.encode(state)] = circuit.encode(tuple(final))
    return routing._permutation(tuple(images))


def repair_branches(law, seam, state, *, recover=True):
    """Execute all30 native seams with every raw source record retained."""
    law, seam, value = routing._candidate(law), _seam(seam), routing._state(state)
    if type(recover) is not bool:
        raise ValueError('active recovery must be an explicit Boolean')
    route = _route_permutation(seam)
    value = _apply_permutation(route, value)
    branches = recorded.repair_branches(law, value)
    for branch in branches:
        coin, raw = branch['coin'], branch['record_code']
        amplitude = _apply_permutation(route, branch['amplitude'], True)
        applied = recover and not branch['abort']
        if applied:
            amplitude = apply_recovery(law, seam, coin, raw, amplitude)
        branch.update(amplitude=_finite(amplitude), seam=seam,
                      recovery_applied=applied,
                      code_branch_weight=branch_weight(law, seam, coin, raw),
                      recovery_resources=recovery_resources(law, seam, coin, raw)
                      if applied else dict.fromkeys(routing.COUNT_KEYS, 0))
    return branches


def _counts(*rows):
    return {key: sum(row[key] for row in rows) for key in routing.COUNT_KEYS}


@lru_cache(maxsize=1)
def _macro_resource_counts():
    # Cache only immutable per-invocation counts, never a caller-owned report.
    macro = compiler.resources(routing._canonical_word())
    one_x = compiler.resources([compiler.pulse('X', math.pi/2, [0])])
    return (tuple(macro[key] for key in routing.COUNT_KEYS),
            tuple(one_x[key] for key in routing.COUNT_KEYS))


def _permutation_resources(images):
    edges = routing.permutation_edges(images)
    conjugations = 2*sum(pattern.bit_count() for _, pattern in edges)
    macro, one_x = _macro_resource_counts()
    counts = {key: len(edges)*m+conjugations*x
              for key, m, x in zip(routing.COUNT_KEYS, macro, one_x, strict=True)}
    return {**counts, 'controlled_flip_invocations': len(edges),
            'conjugating_source_x_pulses': conjugations, 'whole_register_dimension': SIZE}


@lru_cache(maxsize=8)
def _mixing_resource_counts(count, inverse=False):
    resources = compiler.resources(_mixing_word(count, inverse))
    return tuple(resources[key] for key in routing.COUNT_KEYS)


def _mixing_resources(count, inverse=False):
    return dict(zip(routing.COUNT_KEYS, _mixing_resource_counts(count, inverse), strict=True))


def encoding_resources():
    counts = _counts(_permutation_resources(_encoder_permutation()),
                     _mixing_resources(6))
    return {**counts, 'supplied_logical_levels': 6, 'supplied_register_qubits': 3,
            'blank_code_qubits': 6, 'new_classical_records': 0,
            'full_register_dimension': SIZE}


def decoding_resources():
    return {**encoding_resources(), 'blank_code_qubits': 0,
            'output_rows_retained': SIZE, 'syndrome_postselection_applied': False}


def recovery_resources(law, seam, coin, raw):
    law, seam, coin, raw = _label(law, seam, coin, raw)
    count = len(branch_images(law, seam, coin, raw)[0])
    admitted = raw in circuit.record_values(law, coin)
    counts = dict.fromkeys(routing.COUNT_KEYS, 0)
    if admitted and count:
        counts = _counts(encoding_resources(),
                         _permutation_resources(_branch_permutation(law, seam, coin, raw)),
                         _mixing_resources(count, True))
    return {**counts, 'code_branch_weight_exact': str(F(count, 12)),
            'raw_abort': not admitted,
            'identity_on_valid_zero_weight_branch': admitted and not count,
            'blank_code_qubits': 0, 'new_classical_records': 0,
            'archive_record_rewrites': 0, 'input_dependent_normalization': False}


@lru_cache(maxsize=2)
def _native_repair_counts(law):
    resources = routing.routed_repair_resources(law, (0, 1))
    return tuple(resources[key] for key in routing.COUNT_KEYS)


def repair_resources(law, seam):
    law, seam = routing._candidate(law), _seam(seam)
    route = _permutation_resources(_route_permutation(seam))
    repair = dict(zip(routing.COUNT_KEYS, _native_repair_counts(law), strict=True))
    return {**_counts(route, route, repair), 'retained_classical_bits': 1+circuit.RECORD_BITS[law],
            'full_register_dimension': SIZE,
            'recovery_cost_is_record_dependent': True}


def exact_certificate():
    """Exact combinatorial witness and the arbitrary-code sharp capacity proof.

    With the record retained the KL conditions are V*D_r V=p_r I, since
    different output labels are orthogonal.  For each edge ab, C and E both
    contain D_ab=n_a(1-n_b) and D_ba.  Thus compressed n_a-n_b is scalar.
    Connectivity and sum_a n_a=2I force every compressed n_a scalar.  Then
    compressed n_a n_b is scalar and rank at most one, hence zero if k>1.
    The code avoids all30 native-edge basis rows.  Each remaining n_a has
    rank6, and sum_a n_a=2I forces some nonzero scalar compression: k<=6.
    This argument concerns arbitrary complex subspaces and arbitrary CPTP
    recovery with full classical records and unrestricted reference systems.
    """
    all_pairs = set(combinations(range(12), 2))
    flattened = tuple(pair for row in MATCHINGS for pair in row)
    nonedges = all_pairs-set(SEAMS)
    if len(set(flattened)) != 36 or set(flattened) != nonedges:
        raise ArithmeticError('six matchings no longer partition the source nonedges')
    if any(sorted(p for pair in row for p in pair) != list(range(12)) for row in MATCHINGS):
        raise ArithmeticError('a declared logical support is not a perfect matching')
    neighbors = {v: {b if a == v else a for a, b in SEAMS if v in (a, b)} for v in range(12)}
    reached = {0}
    for _ in range(12):
        reached |= set().union(*(neighbors[v] for v in reached))
    if len(reached) != 12 or {len(neighbors[v]) for v in neighbors} != {5}:
        raise ArithmeticError('full native seam graph lost connectivity or degree5')
    basis = native.states(binary=True)
    index = {state: i for i, state in enumerate(basis)}
    code_indices = tuple({index[state] for state in row} for row in SUPPORTS)
    law_reports = {}
    for law in ('C', 'E'):
        counts = Counter()
        for seam in SEAMS:
            total = F(0)
            for coin, label, mapping in coins.coin_maps(law, seam):
                ins = {i for i, _ in mapping}
                overlaps = tuple(len(ins & row) for row in code_indices)
                if len(set(overlaps)) != 1:
                    raise ArithmeticError('source-enumerated code Gram is nonscalar')
                state = basis[next(iter(ins))]
                _, emitted = native.native(state, seam, coin)
                raw = {0: 0, 1: 1, -1: 2}[emitted[1]]
                if law == 'E':
                    raw += 4*((state[seam[0]]+state[seam[1]]) % 2)
                independent = branch_images(law, seam, coin, raw)
                actual = tuple(tuple(circuit.encode(basis[o]) for i, o in mapping if i in row)
                               for row in code_indices)
                if any(set(a) != set(b) for a, b in zip(actual, independent, strict=True)):
                    raise ArithmeticError('independent endpoint oracle differs from native maps')
                weight = F(overlaps[0], 12)
                counts[str(weight)] += 1
                total += weight
            if total != 1:
                raise ArithmeticError('source branches do not preserve complete probability')
        law_reports[law] = {'nonzero_record_branches': sum(counts.values()),
                            'probability_counts_across_thirty_seams': dict(sorted(counts.items()))}
    return {'logical_dimension': 6, 'native_binary_carrier_dimension': 66,
            'source_execution_register_dimension': SIZE,
            'matchings': [[list(pair) for pair in row] for row in MATCHINGS],
            'seams_checked': 30, 'native_graph_connected': True,
            'native_graph_degree': 5, 'nonedge_graph_degree': 6,
            'all_exact_branch_gram_identities': True,
            'maximum_exact_record_correctable_code_dimension': 6,
            'rank_fifteen_exact_record_correctable': False,
            'arbitrary_complex_codes_and_references_in_upper_bound': True,
            'record_conditioned_arbitrary_CPTP_recovery_in_upper_bound': True,
            'complete_unitary_recovery_exists_on_binary66_itself': True,
            'same66_recovery': 'Inverse endpoint swap followed by six disjoint matching-block Householder reflections; no additional quantum carrier in this algebraic extension.',
            'source_compilation_extension': 'The independently compiled permutation/preparation recovery uses the existing full512 source register and may use padding during intermediate controls; no fresh quantum ancilla is added per recovery.',
            'upper_bound_scope': 'One fixed input code precedes state-independent selection with positive probability for every native seam; the seam, coin and C/E value are retained. Each branch is corrected before the next attempt. Seam-dependent re-encoding before the native measurement is a different contract.',
            'lower_bound_scope': 'Six explicit matching superpositions with full source-word active recovery; coherent controls, code load/export and initialization remain specified services.',
            'laws': law_reports}


def execution_report():
    """Numerical source replay on two declared seams, with all code columns.

    The exact certificate covers every seam.  This separate floating audit
    executes the existing source words for C at (0,1) and E at (2,10),
    including general port routing, every raw record and all512 output rows.
    It does not call numerical replay a physical accuracy certificate.
    """
    expected = ideal_code_matrix()
    executed = code_matrix()
    rows = {}
    for law, seam in (('C', (0, 1)), ('E', (2, 10))):
        branches = repair_branches(law, seam, executed)
        gram = sum((b['amplitude'].conj().T@b['amplitude'] for b in branches),
                   np.zeros((DIMENSION, DIMENSION), complex))
        branch_rows = []
        for branch in branches:
            probability = branch['code_branch_weight']
            target = math.sqrt(float(probability))*expected
            branch_rows.append({
                'coin': branch['coin'], 'raw': branch['record_code'],
                'abort': branch['abort'], 'probability_exact': str(probability),
                'operator_error': float(np.linalg.norm(branch['amplitude']-target, 2)),
                'source_recovery_resources': branch['recovery_resources'],
            })
        rows[law] = {
            'seam': list(seam), 'raw_branches_retained': len(branches),
            'reference_gram_error': float(np.linalg.norm(gram-np.eye(DIMENSION), 2)),
            'maximum_branch_operator_error': max(b['operator_error'] for b in branch_rows),
            'complete_retained_output_operator_error': float(np.linalg.norm(
                np.concatenate([b['amplitude']-math.sqrt(float(b['code_branch_weight']))*expected
                                for b in branches], axis=0), 2)),
            'source_repair_resources': repair_resources(law, seam),
            'branches': branch_rows,
        }
    _, forward = routing._executed_macro(False)
    _, inverse = routing._executed_macro(True)
    encoder_error = float(np.linalg.norm(executed-expected, 2))
    diagnostics = [encoder_error, *(row[key] for row in rows.values() for key in (
        'reference_gram_error', 'maximum_branch_operator_error',
        'complete_retained_output_operator_error'))]
    if any(not math.isfinite(value) or value > 1e-8 for value in diagnostics):
        raise ArithmeticError('complete matching-code source replay exceeds its numerical gate')
    return {'encoder_operator_error': encoder_error,
            'maximum_allowed_numerical_execution_error': 1e-8,
            'encoder_resources': encoding_resources(),
            'executed_code_columns': DIMENSION, 'output_rows_retained': SIZE,
            'all_thirty_seam_numerical_source_replay': False,
            'all_thirty_seam_exact_native_certificate': True,
            'forward_full512_source_macro_guard': dict(forward),
            'inverse_full512_source_macro_guard': dict(inverse),
            'code_or_padding_projection_applied': False,
            'physical_service_or_roundoff_certificate': False,
            'laws': rows}


def _ideal_primitive(kind, angle):
    """Analytical comparison only; never supplies an executed gate."""
    if kind == 'J':
        return np.diag([np.exp(1j*angle)]*3+[1.])
    identity = np.eye(2, dtype=complex)
    generators = {
        'X': np.asarray([[0, 1j], [1j, 0]]),
        'Y': np.asarray([[0, 1], [-1, 0]], complex),
        'Z': np.diag([1j, -1j]),
        'I': 1j*identity,
    }
    generators['H'] = (generators['X']+generators['Z'])/math.sqrt(2)
    return math.cos(angle)*identity+math.sin(angle)*generators[kind]


@lru_cache(maxsize=2048)
def _primitive_error(kind, angle):
    dimension = 4 if kind == 'J' else 2
    actual = compiler.native_pulse(kind, angle)[:dimension, :dimension]
    return float(np.linalg.norm(actual-_ideal_primitive(kind, angle), 2))


def _word_log_error(word):
    return math.fsum(math.log1p(_primitive_error(e['kind'], e['angle'])) for e in word)


def _mixing_synthesis_report(count):
    """Separate angle synthesis from source realization of requested angles.

    An exactly unitary requested word U can be changed to a unitary T with
    T|0> equal to the declared uniform vector by a rotation on their span.
    Its operator change is at most twice the first-column discrepancy.  The
    same bound applies to U* versus T*.  This comparison supplies no gate
    to the execution path.  Values are numerical diagnostics, not intervals.
    """
    target = np.asarray([1/math.sqrt(count)]*count+[0.]*(8-count), complex)
    word = compiler.compile_real_state(target.real)
    requested = np.eye(8, dtype=complex)
    for event in word:
        requested = compiler._apply(requested, _ideal_primitive(event['kind'], event['angle']),
                                    event['targets'], 3)
    column_error = float(np.linalg.norm(requested[:, 0]-target))
    correction = 2*column_error
    residual = float(np.linalg.norm(requested.conj().T@requested-np.eye(8), 2))
    if any(not math.isfinite(value) or value > 1e-10 for value in (correction, residual)):
        raise ArithmeticError('uniform-support angle synthesis exceeds its numerical gate')
    return {'support_count': count,
            'requested_word_uniform_column_error': column_error,
            'uniform_target_unitary_completion_error_upper': correction,
            'requested_word_unitarity_residual': residual,
            'forward_and_inverse_use_the_same_target_completion_bound': True}


def _native_dilation_error(law, coin):
    """Compare every local source column to its independent integer dilation."""
    dimension = 1 << (3+circuit.RECORD_BITS[law])
    lower = list(range(32))
    for left, right in circuit.REPAIR_SWAPS['C', coin]:
        lower[left], lower[right] = lower[right], lower[left]
    images = []
    for address in range(dimension):
        final = (address & ~31)+lower[address & 31]
        if law == 'E' and (final >> coin) & 1:
            final ^= 32
        images.append(final)
    expected = np.zeros((dimension, dimension), complex)
    expected[images, np.arange(dimension)] = 1
    return float(np.linalg.norm(recorded._repair(law, coin)-expected, 2))


@lru_cache(maxsize=1)
def _source_composition_report():
    """Uniform numerical error budget from every all30 source schedule.

    If A_j differs from an exactly unitary U_j by delta_j, telescoping
    bounds the full product error by product_j(1+delta_j)-1.  Every emitted
    macro and conjugation is counted, with no ideal target used in execution.
    Each macro error includes both its stabilized target residual and the
    full-register polar change, so every invocation covers the raw executed
    source word as well as the stabilized numerical replay.
    The maximum over orthogonal classical control blocks costs a maximum,
    not a sum over their probabilities.  The calculated residuals and sums
    use floating arithmetic: this is a numerical composition diagnostic,
    not a directed-rounding, hardware or physical-service certificate.
    """
    exact_flip = np.eye(SIZE, dtype=complex)
    exact_flip[[0, 256]] = exact_flip[[256, 0]]
    macro_errors, stabilized_errors, projection_changes = {}, {}, {}
    for inverse in (False, True):
        actual, guard = routing._executed_macro(inverse)
        stabilized_errors[inverse] = float(np.linalg.norm(actual-exact_flip, 2))
        # compiler.stabilize_executed_unitary evaluates this change with the
        # spectral norm on all512 rows and columns.  The triangle inequality
        # charges it once per macro, not once per entire composed history.
        projection_changes[inverse] = guard['projection_change']
        macro_errors[inverse] = stabilized_errors[inverse]+projection_changes[inverse]
        if (guard['dimension'] != SIZE or not math.isfinite(macro_errors[inverse])
                or not 0 <= projection_changes[inverse] <= guard['maximum_projection_change']):
            raise ArithmeticError('full-register macro stabilization guard is invalid')
    exact_ix = np.asarray([[0, 1j], [1j, 0]])
    xerror = max(float(np.linalg.norm(compiler.native_pulse('X', sign*math.pi/2)[:2, :2]
                                      -sign*exact_ix, 2)) for sign in (1, -1))
    synthesis = {count: _mixing_synthesis_report(count) for count in (1, 4, 5, 6)}

    mixing_logs = {(count, inverse): (_word_log_error(_mixing_word(count, inverse))
                   +math.log1p(synthesis[count]['uniform_target_unitary_completion_error_upper']))
                   for count in (1, 4, 5, 6) for inverse in (False, True)}

    def mixing_log(count, inverse=False):
        return mixing_logs[count, inverse]

    def permutation_log(images, inverse=False):
        edges = routing.permutation_edges(images, inverse=inverse)
        xcount = 2*sum(pattern.bit_count() for _, pattern in edges)
        return len(edges)*math.log1p(macro_errors[inverse])+xcount*math.log1p(xerror)

    encoder_log = permutation_log(_encoder_permutation())+mixing_log(6)
    decoder_log = permutation_log(_encoder_permutation(), True)+mixing_log(6, True)
    encoder_resources = encoding_resources()
    exact_coin = np.asarray([[1., -1.], [1., 1.]])/math.sqrt(2)
    coin_error = float(np.linalg.norm(compiler.native_pulse('Y', -math.pi/4)[:2, :2]-exact_coin, 2))
    coin_log = math.log1p(coin_error)
    rows = {}
    for law in ('C', 'E'):
        native_errors = [_native_dilation_error(law, coin) for coin in (0, 1)]
        native_log = math.log1p(max(native_errors))
        worst_log = 0.
        worst_resources = dict.fromkeys(routing.COUNT_KEYS, 0)
        worst_recovery_pulses = 0
        for seam in SEAMS:
            route = _route_permutation(seam)
            route_log = permutation_log(route)+permutation_log(route, True)
            largest_recovery_log = 0.
            largest_recovery_resources = dict.fromkeys(routing.COUNT_KEYS, 0)
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    count = len(branch_images(law, seam, coin, raw)[0])
                    if not count:
                        continue
                    current_log = (permutation_log(_branch_permutation(law, seam, coin, raw), True)
                                   +mixing_log(count, True)+encoder_log)
                    largest_recovery_log = max(largest_recovery_log, current_log)
                    resources = recovery_resources(law, seam, coin, raw)
                    for key in routing.COUNT_KEYS:
                        largest_recovery_resources[key] = max(largest_recovery_resources[key], resources[key])
            worst_log = max(worst_log, route_log+coin_log+native_log+largest_recovery_log)
            costs = _counts(repair_resources(law, seam), largest_recovery_resources)
            for key in routing.COUNT_KEYS:
                worst_resources[key] = max(worst_resources[key], costs[key])
            worst_recovery_pulses = max(worst_recovery_pulses, largest_recovery_resources['native_pulses'])
        rows[law] = {
            'maximum_corrected_attempt_operator_error': math.expm1(worst_log),
            'maximum_encoded_corrected_attempt_operator_error': math.expm1(encoder_log+worst_log),
            'maximum_corrected_attempt_log_error': worst_log,
            'maximum_corrected_attempt_resources': worst_resources,
            'maximum_recovery_native_pulses': worst_recovery_pulses,
            'complete_local_native_dilation_operator_errors_by_coin': native_errors,
            'all_thirty_seam_schedules_and_raw_recoveries_counted': True,
        }
    if any(not math.isfinite(row['maximum_encoded_corrected_attempt_operator_error'])
           or row['maximum_encoded_corrected_attempt_operator_error'] > 1e-6
           for row in rows.values()):
        raise ArithmeticError('all-seam source composition exceeds its numerical gate')
    return {
        'method': 'product(1+factor_operator_error)-1; maximum over retained classical control blocks',
        'forward_full512_macro_operator_error': macro_errors[False],
        'inverse_full512_macro_operator_error': macro_errors[True],
        'forward_full512_stabilized_macro_operator_error': stabilized_errors[False],
        'inverse_full512_stabilized_macro_operator_error': stabilized_errors[True],
        'forward_full512_macro_polar_change': projection_changes[False],
        'inverse_full512_macro_polar_change': projection_changes[True],
        'raw_and_stabilized_macro_conventions_both_covered': True,
        'macro_error_convention': 'Stabilized-to-exact operator residual plus raw-to-stabilized spectral norm change, charged at every emitted macro invocation.',
        'source_X_conjugation_operator_error': xerror,
        'fair_coin_preparation_operator_error': coin_error,
        'mixing_angle_synthesis_diagnostics': {str(k): v for k, v in synthesis.items()},
        'encoder_operator_error_bound': math.expm1(encoder_log),
        'encoder_log_error': encoder_log,
        'decoder_operator_error_bound': math.expm1(decoder_log),
        'decoder_log_error': decoder_log,
        'encoder_resources': encoder_resources,
        'numerical_source_composition_bound_not_interval_or_hardware_certificate': True,
        'unrestricted_reference_extension_has_same_operator_bound': True,
        'complete_raw_classical_control_blocks_retained': True,
        'laws': rows,
    }


def source_composition_report():
    """Return a fresh report; caller mutations cannot poison cached bounds."""
    return deepcopy(_source_composition_report())


source_execution_report = execution_report
