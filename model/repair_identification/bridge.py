"""Current-only adaptation to the canonical continuing observer runtime.

Raw source execution remains on 512 amplitudes. The runtime uses the eight
active columns, including both padding states. Its numerical reduction is
measured explicitly; neither amplitudes nor probabilities are normalized.
"""
from fractions import Fraction as F
from itertools import product
import math

import numpy as np

from . import core, source
from ..record_laws import family
from ..primitive_observer import runtime
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler


MENU = tuple(product((-1, 0, 1), (-1, 1)))+((None, None),)
TOLERANCE = 2e-9


def outcome(current, pointer):
    return MENU.index((current, pointer))


def maximum_cost(rows):
    rows = [dict(row) for row in rows]
    return tuple((key, max(row.get(key, 0) for row in rows))
                 for key in sorted(set().union(*(row.keys() for row in rows))))


def _norm(matrix):
    return float(np.linalg.norm(matrix, 2))


def reduction_diagnostics(matrices):
    """Measured full-column reduction error, not a roundoff certificate.

    If V stacks the complete raw Kraus maps and W their embedded active
    rows, ||V-W|| (||V||+||W||)/2 bounds the exact-arithmetic half trace
    distance of the represented maps, including any untouched reference.
    The reported norms themselves are evaluated in binary64.
    """
    if not matrices or any(type(k) is not np.ndarray or k.shape != (512, 8)
                           or not np.isfinite(k).all() for k in matrices):
        raise ValueError('every complete finite emitted 512 by 8 Kraus map is required')
    full = sum((k.conj().T@k for k in matrices), np.zeros((8, 8), complex))
    active = sum((k[:8].conj().T@k[:8] for k in matrices), np.zeros((8, 8), complex))
    tail = sum((k[8:].conj().T@k[8:] for k in matrices), np.zeros((8, 8), complex))
    leakage = math.sqrt(_norm(tail))
    bound = leakage*(math.sqrt(_norm(full))+math.sqrt(_norm(active)))/2
    result = {'full_tp_residual': _norm(full-np.eye(8)),
              'active_tp_residual': _norm(active-np.eye(8)),
              'stacked_leakage_norm': leakage, 'projection_change_bound_numeric': bound}
    if max(result.values()) > TOLERANCE:
        raise ArithmeticError('source reduction exceeds the declared numerical guard')
    return result


def program(candidate, probe, reference=False):
    """Two complete compound events; only public current/pointer drives Z."""
    source._candidate(candidate); source.settings(probe); source._flags(reference)
    raw = source.event_matrices(candidate, probe)
    diagnostics = reduction_diagnostics([branch.amplitude for branch in raw])
    groups = [[] for _ in MENU]
    for branch in raw:
        index = outcome(None, None) if branch.abort else outcome(branch.current, branch.record[-1])
        groups[index].append(branch.amplitude[:8])
    if any(not group for group in groups):
        raise ArithmeticError('a complete raw outcome group disappeared')
    costs = maximum_cost(branch.resources for branch in raw)
    binding = {'candidate': candidate, 'probe': probe, 'words': dict(source.word_hashes(probe, reference)),
               'visible_menu': [list(item) for item in MENU],
               'resource_semantics': 'componentwise upper bound over private raw outcomes'}
    normal = runtime.Instrument(tuple(tuple(group) for group in groups), True,
                                resources=costs, binding=binding)
    feedback = compiler.execute(source.feedback_word(), 3)
    controlled = runtime.Instrument(tuple(tuple(k@feedback for k in group) for group in groups), True,
        resources=runtime.add_resources(costs, source._resources(source.feedback_word(),
                                         feedback_phase_controls=1)), binding={**binding, 'feedback': 'Z bit 0'})
    idle = runtime.Instrument(tuple((np.eye(8) if i == 6 else np.zeros((8, 8)),)
                                    for i in range(7)), True, binding={'after_abort': 'retain and halt'})
    prepared = source.preparation(probe, reference)
    active = prepared[:8].reshape(8, -1)
    density = active@active.conj().T
    declared = runtime.Program('complete repair identification', 8, density,
        (runtime.Event('repair and read 1', (normal,)),
         runtime.Event('repair and read 2', tuple(idle if i == 6 else
             controlled if label == (0, 1) else normal for i, label in enumerate(MENU)), control=0)),
        source.preparation_resources(probe, reference), binding)
    initial = np.outer(active.ravel(), active.ravel().conj()) if reference else None
    diagnostics['preparation_tail_mass'] = float(np.vdot(prepared[8:], prepared[8:]).real)
    return declared, initial, diagnostics


def run_history(candidate, probe, reference=False):
    declared, initial, diagnostics = program(candidate, probe, reference)
    whole = runtime.execute(declared, initial=initial)
    prefix = runtime.execute(declared, initial=initial, stop_after=1)
    payload = runtime.checkpoint_to_json(prefix)
    restored = runtime.checkpoint_from_json(declared, payload)
    resumed = runtime.execute(declared, checkpoint=restored)
    delta = distance({branch.record: branch.density for branch in whole.branches},
                     {branch.record: branch.density for branch in resumed.branches})
    if (delta != 0 or whole.resources != resumed.resources or
            any(a.record != b.record or a.resources != b.resources or not np.array_equal(a.density, b.density)
                for a, b in zip(whole.branches, resumed.branches, strict=True))):
        raise ArithmeticError('canonical checkpoint continuation changed the complete history')
    return whole, {**diagnostics, 'checkpoint_replay_distance': delta,
                   'checkpoint_bytes': len(payload.encode('utf-8'))}


def _complex(value):
    return complex(float(value.a)+float(value.b)*math.sqrt(2),
                   float(value.c)+float(value.d)*math.sqrt(2))


def _read(operator, pair, pointer):
    def column(index):
        if index not in pair:
            return {index: F(1)} if pointer == -1 else {}
        return {j: F(1, 2) if pointer == 1 or j == index else -F(1, 2) for j in pair}
    result = {}
    for ((x, r), (y, s)), value in operator.items():
        for left, a in column(x).items():
            for right, b in column(y).items():
                key = ((left, r), (right, s))
                result[key] = result.get(key, family.S())+value*a*b
    return {key: value for key, value in result.items() if value}


def exact_history(table, probe, reference=False):
    """Exact 78-state sparse oracle, retaining a minor reference factor."""
    table, descriptor = core.validate_table(table), core._probe(probe)
    source._flags(reference)
    left, right = descriptor.input_pair
    phase = family.S(1) if descriptor.phase == 0 else family.S(0, 0, 1)
    amplitudes = {(left, 0): family.S(1), (right, int(reference)): phase}
    state = {(x, y): a*b.conj()*F(1, 2) for x, a in amplitudes.items() for y, b in amplitudes.items()}
    branches = {(): state}
    for step in range(2):
        output = {}
        for record, rho in branches.items():
            if step and MENU[record[0]] == (0, 1):
                rho = {(x, y): value*((-1)**((circuit.SOURCE_CODES[x[0]]&1)+(circuit.SOURCE_CODES[y[0]]&1)))
                       for (x, y), value in rho.items()}
            for current in (-1, 0, 1):
                moved = {}
                for r, s in product(range(2 if reference else 1), repeat=2):
                    work = {(x, y): value for ((x, a), (y, b)), value in rho.items() if (a, b) == (r, s)}
                    moved.update({((x, r), (y, s)): value
                                  for (x, y), value in core.apply_current(table, current, work).items()})
                for pointer in (-1, 1):
                    output[(*record, outcome(current, pointer))] = _read(moved, descriptor.output_pair, pointer)
            output[(*record, 6)] = {}
        branches = output
    ref = 2 if reference else 1
    result = {}
    for record, rho in branches.items():
        matrix = np.zeros((8*ref, 8*ref), complex)
        for ((x, r), (y, s)), value in rho.items():
            if max(circuit.SOURCE_CODES[x], circuit.SOURCE_CODES[y]) >= 8:
                raise ArithmeticError('the exact continuing oracle escaped its active chart')
            matrix[circuit.SOURCE_CODES[x]*ref+r, circuit.SOURCE_CODES[y]*ref+s] += _complex(value)
        result[record] = matrix
    return result


def coarse_raw(history):
    """Keep all raw mass; bound discarded rows and their coherent cross terms.

    For each unnormalized raw |v>=(a,b), the half trace distance to
    |a><a| embedded in the full space is <=||a||||b||+||b||^2/2.
    Summing before any record coarsening gives a conservative complete-output
    bound without allocating 49 dense 1024 by 1024 matrices.
    """
    validate_raw_history(history)
    ref = history.reference_dimension
    result = {record: np.zeros((8*ref, 8*ref), complex) for record in product(range(7), repeat=2)}
    reduction, tail, mass, abort = [], [], [], []
    for branch in history.branches:
        labels = tuple(6 if raw is None or current is None else outcome(current, raw[-1])
                       for current, raw in zip(branch.currents, branch.record, strict=True))
        value = branch.amplitude.reshape(512, ref)
        active = value[:8].ravel()
        small, rest = float(np.vdot(active, active).real), float(np.vdot(value[8:], value[8:]).real)
        result[labels] += np.outer(active, active.conj())
        tail.append(rest); mass.append(small+rest)
        abort.append(small+rest if branch.abort else 0.)
        reduction.append(math.sqrt(small)*math.sqrt(rest)+rest/2)
    return result, {'full_probability_mass': math.fsum(mass), 'tail_mass': math.fsum(tail),
                    'full_output_projection_bound_numeric': math.fsum(reduction),
                    'abort_probability': math.fsum(abort), 'raw_history_count': len(history.branches)}


def distance(first, second):
    for blocks in (first, second):
        if (type(blocks) is not dict or not blocks or
                any(type(key) is not tuple or any(type(i) is not int for i in key) for key in blocks)):
            raise ValueError('original nonempty record-to-matrix dictionaries are required')
        for matrix in blocks.values():
            if (type(matrix) is not np.ndarray or matrix.ndim != 2 or not matrix.shape[0]
                    or matrix.shape[0] != matrix.shape[1] or matrix.size > runtime.MAX_DENSITY_ENTRIES
                    or matrix.dtype not in (np.dtype('float64'), np.dtype('complex128'))
                    or not np.isfinite(matrix).all()):
                raise ValueError('original finite bounded numeric quantum matrices are required')
    if first.keys() != second.keys():
        raise ValueError('matching complete visible record domains required')
    if any(first[key].shape != second[key].shape for key in first):
        raise ValueError('matching quantum output dimensions are required')
    return math.fsum(float(np.linalg.svd(a-b, compute_uv=False).sum())/2
                     for key, a in first.items() for b in (second[key],))


def _original_equal(first, second):
    if type(first) is not type(second):
        return False
    if type(first) is tuple:
        return len(first) == len(second) and all(_original_equal(a, b) for a, b in zip(first, second))
    return first == second


def validate_raw_history(history):
    """Require every raw outcome, including zero aborts, and exact service costs."""
    if type(history) is not source.HistoryExecution:
        raise ValueError('an original complete source history is required')
    source._candidate(history.candidate); source.settings(history.probe)
    source._flags(feedback=history.feedback)
    if type(history.reference_dimension) is not int or history.reference_dimension not in (1, 2):
        raise ValueError('the declared untouched reference has dimension one or two')
    reference = history.reference_dimension == 2
    if not _original_equal(history.word_hashes, source.word_hashes(history.probe, reference)):
        raise ValueError('the continuing history must retain its emitted word bindings')
    raw = source.event_matrices(history.candidate, history.probe)
    expected = []
    for first in raw:
        costs = runtime.add_resources(source.preparation_resources(history.probe, reference), first.resources)
        if first.abort:
            expected.append(((first.record, None), (first.current, None), True, costs))
            continue
        if history.feedback and (first.current, first.record[-1]) == (0, 1):
            costs = runtime.add_resources(costs, source._resources(source.feedback_word(), feedback_phase_controls=1))
        expected.extend(((first.record, final.record), (first.current, final.current), final.abort,
                         runtime.add_resources(costs, final.resources)) for final in raw)
    if type(history.branches) is not tuple or len(history.branches) != len(expected):
        raise ValueError('all raw continuing outcomes, including zero failures, are required')
    shape = (512, 2) if reference else (512,)
    for branch, header in zip(history.branches, expected, strict=True):
        if (type(branch) is not source.HistoryBranch
                or not _original_equal((branch.record, branch.currents, branch.abort, branch.resources), header)
                or type(branch.amplitude) is not np.ndarray or branch.amplitude.shape != shape
                or branch.amplitude.dtype != np.dtype(complex) or not np.isfinite(branch.amplitude).all()):
            raise ValueError('raw branch records, costs and full finite amplitudes must match the declared protocol')


def validate_probe_execution(execution, candidate, probe):
    """Validate the full single-probe grammar before accumulating any statistic."""
    if (type(execution) is not source.ProbeExecution or
            not _original_equal((execution.candidate, execution.probe, execution.reference_dimension),
                                (candidate, probe, 1)) or
            not _original_equal(execution.word_hashes, source.word_hashes(probe))):
        raise ValueError('the executed probe must bind the requested candidate, preparation and words')
    expected = source.event_matrices(candidate, probe)
    if type(execution.branches) is not tuple or len(execution.branches) != len(expected):
        raise ValueError('all single-probe raw outcomes, including zero failures, are required')
    for branch, declared in zip(execution.branches, expected, strict=True):
        resources = runtime.add_resources(source.preparation_resources(probe), declared.resources)
        header = (declared.record, declared.current, declared.contrast, declared.abort, resources)
        if (type(branch) is not source.Branch or
                not _original_equal((branch.record, branch.current, branch.contrast, branch.abort, branch.resources), header)
                or type(branch.amplitude) is not np.ndarray or branch.amplitude.shape != (512,)
                or branch.amplitude.dtype != np.dtype(complex) or not np.isfinite(branch.amplitude).all()):
            raise ValueError('probe branch records, costs and finite full amplitudes must match the protocol')
