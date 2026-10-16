"""Actual signed native stopping branches and their exact finite targets.

The native C/E instrument is reused from native_relational_readout, with its
moving outcome split by the actual signed current of each partial injection.
Exact targets use its integer zero-current Schur counts and moving maps.
Numerical execution applies the original Kraus arrays to the original phase
preparation; exact targets are never substituted for those executed arrays.

An outcome (n, -1) or (n, +1) means the first current at controller slot n.
The single (N, 0) outcome is the unresolved tail at the declared cap N. These
records encode every possible stopped word 0**(n-1),j or the all-zero word;
no successful outcome is conditionally normalized. The counter is not a
physical duration. The coherent preparation and uniform scheduler are inputs.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from ..common_source import exact_output as exact
from ..native_relational_readout import construction as relational
from ..native_quantum_history import construction as quantum
from ..primitive_observer import runtime


CURRENTS = (-1, 0, 1)
MAX_HORIZON = 128
SOURCE_PATHS = tuple(dict.fromkeys((
    'model/clocked_common_observer/clock.py',
    'model/common_source/exact_output.py', *relational.SOURCE_PATHS,
)))


def _horizon(value):
    if type(value) is not int or not 1 <= value <= MAX_HORIZON:
        raise ValueError('an integer native-clock cap from one through 128 is required')
    return value


def _phase(value):
    if type(value) is not int or value not in (-1, 1):
        raise ValueError('an integer preparation phase sign minus or plus one is required')
    return value


def _activity(value):
    if type(value) is not F or not 0 < value <= 1:
        raise ValueError('an exact Fraction activity in (0,1] is required')
    return value


def _contract(law, count, sign, activity):
    return quantum.candidate(law), _horizon(count), _phase(sign), _activity(activity)


def _array(value, shape):
    """Snapshot original base arrays before any subclass method dispatch.

    A subclass can override copy, nonzero or tobytes and thereby conceal
    original entries from an exact comparison. Reject it before the shared
    dtype/domain validator or any array operation inspects supplied data.
    """
    if type(value) is not np.ndarray:
        raise ValueError('an original base numeric ndarray is required')
    return exact._array(value, shape)


@dataclass(frozen=True)
class SignedSource:
    law: str
    no_current_counts: tuple
    # Each item is (signed current, tuple((input, output), ...)).
    moving_maps: tuple
    instrument: runtime.Instrument


def _current(mapping):
    """Recover the actual current and seam from native input/output loads."""
    seams = {frozenset(seam): seam for seam in quantum.mixing.seams()}
    address = None
    for initial, final in mapping:
        before, after = quantum.STATES[initial], quantum.STATES[final]
        changed = frozenset(i for i in range(12) if before[i] != after[i])
        if len(changed) != 2 or changed not in seams:
            raise ArithmeticError('a moving native partial injection lacks one source seam')
        a, b = seams[changed]
        current = before[a]-after[a]
        if current not in (-1, 1) or before[b]-after[b] != -current:
            raise ArithmeticError('a moving native partial injection lost signed current conservation')
        candidate = (a, b, current)
        if address is not None and candidate != address:
            raise ArithmeticError('a committed native Kraus map mixes signed current outcomes')
        address = candidate
    if address is None:
        raise ArithmeticError('a moving native branch has no input support')
    return address[2]


@lru_cache(maxsize=2)
def _source(law):
    original = relational.source(law)
    groups = {-1: [], 0: list(original.instrument.kraus[0]), 1: []}
    maps = []
    for mapping, matrix in zip(original.moving_maps, original.instrument.kraus[1], strict=True):
        current = _current(mapping)
        groups[current].append(matrix)
        maps.append((current, mapping))
    instrument = runtime.Instrument(kraus=tuple(tuple(groups[j]) for j in CURRENTS), visible=True,
        resources=(('native_attempts', 1), ('native_current_readouts', 1)),
        binding={'candidate': law, 'current_values_in_outcome_order': list(CURRENTS),
                 'scheduler': 'uniform thirty source seams and fair placement coin',
                 'forgotten_labels': ['seam', 'coin', 'committed candidate repair label'],
                 'native_channel_semantics': True, 'native_pulse_words_emitted': False})
    # Completeness and the uniform tail contraction hold on arbitrary native
    # inputs, including a reference: every Kraus map is a partial injection,
    # so its effect is diagonal. Check every effect entry exactly.
    for i in range(66):
        moving = sum(any(left == i for left, _ in mapping) for _, mapping in maps)
        if moving+original.no_current_counts[i][i] != 60:
            raise ArithmeticError('the signed native instrument lost a complete basis-state effect')
    if max(original.no_current_counts[i][i] for i in range(66)) != 52:
        raise ArithmeticError('the native all-input tail contraction changed')
    return SignedSource(law, original.no_current_counts, tuple(maps), instrument)


def source(law):
    """Detach actual signed Kraus groups after validating the original label.

    Frozen dataclass envelopes still expose a mutable __dict__, and a
    read-only ndarray still has mutable shape metadata. Reconstruct the
    instrument through the existing validator, which takes fresh immutable
    byte-backed matrix snapshots. Callers never receive a cached envelope
    or ndarray object; the exact integer tuples can safely remain shared.
    """
    cached = _source(quantum.candidate(law))
    instrument = runtime.Instrument(kraus=cached.instrument.kraus,
        visible=cached.instrument.visible, resources=cached.instrument.resources,
        binding=cached.instrument.binding)
    return SignedSource(cached.law, cached.no_current_counts, cached.moving_maps, instrument)


def records(count):
    """Complete first-event/time-out menu, in native execution order."""
    count = _horizon(count)
    return tuple((n, j) for n in range(1, count+1) for j in (-1, 1))+((count, 0),)


def _trace(state):
    return sum((value for (i, j), value in state.items() if i == j), F(0))


def _moving(state, maps, activity):
    output = {-1: {}, 1: {}}
    for current, entries in maps:
        mapping = dict(entries)
        target = output[current]
        for (i, j), value in state.items():
            if i in mapping and j in mapping:
                key = mapping[i], mapping[j]
                target[key] = target.get(key, F(0))+value*activity/60
    return {current: {key: value for key, value in state.items() if value}
            for current, state in output.items()}


def exact_outputs(law, count, *, sign=1, activity=F(1)):
    """Fresh sparse rational unnormalized outputs for the declared preparation.

    Optional activity q inserts an independent identity idle on the native
    register: T0'=(1-q)Id+q*T0 and Tj'=q*Tj. It does not assert that work also
    idles. In particular, work advancing every slot changes the joint stopped
    observation under this modification, despite native eventual invariance.
    """
    law, count, sign, activity = _contract(law, count, sign, activity)
    marked, state = source(law), relational.exact_preparation(sign)
    outputs = {}
    for n in range(1, count+1):
        for current, stopped in _moving(state, marked.moving_maps, activity).items():
            outputs[n, current] = stopped
        state = {key: value*(1-activity+activity*F(marked.no_current_counts[key[0]][key[1]], 60))
                 for key, value in state.items()}
    outputs[count, 0] = state
    if tuple(outputs) != records(count) or sum(map(_trace, outputs.values()), F(0)) != 1:
        raise ArithmeticError('the complete stopped native instrument lost an outcome or mass')
    return outputs


def exact_eventual_outputs(law, *, sign=1, activity=F(1)):
    """Exact convergent-series comparison, not a finite-budget execution.

    The native-only stopped instrument is invariant under inserted identity
    idling: q/[1-((1-q)+q*w)] = 1/(1-w). This cancellation does not identify
    the joint clock/work instrument when work evolves during native idles.
    """
    law, _, sign, activity = _contract(law, 1, sign, activity)
    marked = source(law)
    resolvent = {key: value/(activity*(1-F(marked.no_current_counts[key[0]][key[1]], 60)))
                 for key, value in relational.exact_preparation(sign).items()}
    result = _moving(resolvent, marked.moving_maps, activity)
    if sum(map(_trace, result.values()), F(0)) != 1:
        raise ArithmeticError('the exact eventual native instrument lost mass')
    return result


def _numerical_steps(law, count, sign, activity):
    """One sequential execution, retaining each unresolved numerical state."""
    marked = source(law)
    state = quantum.phase_apparatus(sign)[0]
    q = float(activity)
    for n in range(1, count+1):
        negative, waiting, positive = runtime._apply(marked.instrument, state, 1)
        stopped = { (n, -1): runtime._immutable(q*negative),
                    (n, 1): runtime._immutable(q*positive) }
        state = runtime._immutable(waiting if activity == 1 else (1-q)*state+q*waiting)
        yield n, stopped, state


def execute_native(law, count, *, sign=1, activity=F(1)):
    """Execute original numerical native maps, returning every stopped block.

    Each iteration applies the complete existing Kraus instrument to the
    retained unresolved density. Reported stopped arrays are never replaced
    by exact_targets, normalized, projected or symmetrized. The prescribed
    preparation has no external reference; the map/tail argument in source
    applies to arbitrary native/reference inputs independently of this run.
    """
    law, count, sign, activity = _contract(law, count, sign, activity)
    outputs = {}
    for _, stopped, state in _numerical_steps(law, count, sign, activity):
        outputs.update(stopped)
    outputs[count, 0] = state
    return outputs


def execute_prefixes(law, count, *, sign=1, activity=F(1)):
    """Complete output dicts at every cap, from one continuing execution.

    Return ``{cap: {(n,current): original_density, ...}}`` for cap=1..count.
    Earlier stopped outputs share immutable arrays. Each cap's own timeout
    is the actual unresolved state at that slot, retained before continuation;
    neither later terminal states nor exact targets stand in for this state.
    """
    law, count, sign, activity = _contract(law, count, sign, activity)
    prefixes, stopped_so_far = {}, {}
    for n, stopped, state in _numerical_steps(law, count, sign, activity):
        stopped_so_far.update(stopped)
        prefixes[n] = {**stopped_so_far, (n, 0): state}
    return prefixes


def compare_outputs(law, count, arrays, *, sign=1, activity=F(1)):
    """Directed full-output error and norm bounds for original native arrays.

    For each reported 66-by-66 matrix B and exact positive target R, compute
    ||B-R||F squared from exact original binary components. Rank <=66 gives
    half ||B-R||1 <= sqrt(66*||B-R||F**2)/2. Consequently ||B||1 is at most
    Tr(R)+2*error, including any non-Hermitian floating defect or norm gain.
    No positivity assumption or repair is imposed on the reported matrix.
    These bounds cover the prescribed phase-prepared output. The separate
    uniform native-tail theorem, not this comparison, covers arbitrary inputs.
    """
    law, count, sign, activity = _contract(law, count, sign, activity)
    menu = records(count)
    if (type(arrays) is not dict or any(type(record) is not tuple or len(record) != 2
            or any(type(value) is not int for value in record) for record in arrays)
            or set(arrays) != set(menu)):
        raise ValueError('every complete native first-event and timeout original array is required')
    targets, results = exact_outputs(law, count, sign=sign, activity=activity), {}
    for record in menu:
        array = _array(arrays[record], (66, 66))
        target, squared = targets[record], F(0)
        # Structural zeros contribute exactly zero, not a numerical cutoff.
        # The union includes every original nonzero, including leakage into
        # entries absent from the ideal support and arbitrarily tiny defects.
        support = set(target) | {(int(i), int(j)) for i, j in zip(*array.nonzero(), strict=True)}
        for i, j in support:
            real, imag = exact._component(array[i, j])
            squared += (real-target.get((i, j), F(0)))**2+imag*imag
        error = exact._sqrt_interval(66*squared)[1]/2
        mass = _trace(target)
        results[record] = {'distance_upper': error, 'ideal_mass': mass,
                           'reported_nuclear_norm_upper': mass+2*error,
                           'frobenius_error_squared': squared}
    return {'distance_upper': sum((row['distance_upper'] for row in results.values()), F(0)),
            'reported_nuclear_norm_upper': sum((row['reported_nuclear_norm_upper']
                                                for row in results.values()), F(0)),
            'per_record': results}


def statistics(law, count, *, sign=1, activity=F(1)):
    """Exact output distribution, all-input tail bound and resource moments.

    No-current effects obey T0*(I)<=13/15 I. With activity q the coefficient
    is 1-2q/15, so the tail bound survives arbitrary reference entanglement
    and any trace-preserving work processing which does not alter native
    scheduling. Expected active attempts equal q times expected active
    controller slots: an independent activity draw occurs at each such slot.
    """
    law, count, sign, activity = _contract(law, count, sign, activity)
    outputs = exact_outputs(law, count, sign=sign, activity=activity)
    probabilities = {record: _trace(state) for record, state in outputs.items()}
    slots = sum((n*probability for (n, _), probability in probabilities.items()), F(0))
    a, b = 1-F(2, 15)*activity, 1-F(1, 6)*activity
    expected = (3*(1-a**count)+F(18, 5)*(1-b**count))/activity
    tail = probabilities[count, 0]
    if slots != expected or tail != F(2, 5)*a**count+F(3, 5)*b**count:
        raise ArithmeticError('signed native outputs differ from the independent phase-preparation clock law')
    return {'probabilities': probabilities, 'timeout_probability': tail,
            'uniform_tail_bound': a**count,
            'expected_controller_slots': slots,
            'expected_native_attempts': activity*slots,
            'expected_native_current_readouts': activity*slots,
            'maximum_controller_slots': count, 'maximum_native_attempts': count,
            'outcome_count': 2*count+1,
            'uniform_bound_all_finite_references': True,
            'physical_duration_or_native_pulse_cost_derived': False}


def outcome_resources(count, record, *, activity=F(1)):
    """Native primitive resource counts; hidden idles require upper bounds.

    The compressed (n,j) record preserves the complete possible stopped
    signed word. If q<1 its zeros do not reveal whether a native attempt
    occurred, so its actual native-attempt cost is not determined. Return
    explicit upper bounds in that case, rather than charging a false count.
    """
    count, activity = _horizon(count), _activity(activity)
    if (type(record) is not tuple or len(record) != 2
            or any(type(value) is not int for value in record) or record not in records(count)):
        raise ValueError('a complete native first-event or unresolved-tail record is required')
    n, _ = record
    suffix = '' if activity == 1 else '_upper'
    return (('controller_slots', n), ('native_attempts'+suffix, n),
            ('native_current_readouts'+suffix, n), ('retained_clock_outcomes', 1))
