"""Local detector laws while the complete supplied scalar field evolves.

Every calibration trial executes the same fixed symmetric onsite/bond
schedule on all eight four-level sites.  No pair or onsite generator is
isolated for inference.  The nine preparations are contracted as independent
state columns to save host work; resource counts retain all 88 fresh trials
at each duration.  There is no coherent superposition of those preparations.

Only opaque labels and complete binary detector laws reach the inverse.
The action family, code services, number preparations, phase-sensitive
instruments, scheduling and independent clock remain supplied apparatus.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
from itertools import combinations
import json
import math
from numbers import Integral, Real

import numpy as np

from ..interacting_reconstruction import source as preparation_source
from ..observer.source_scalar import compiler as c
from ..relativistic_scalar import source as scalar_source
from ..relativistic_scalar.experiment import evolve
from ..relativistic_scalar.lattice import ScalarLattice, ScalarParameters


PERMUTATION = preparation_source.PERMUTATION
PAIRS = tuple(combinations(range(8), 2))
PHASE_SIGNS = (-1, 1)
DIMENSIONS = (4,)*8
COUNT_KEYS = scalar_source.COUNT_KEYS


def _integer(value, name, lower, upper):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or not lower <= value <= upper):
        raise ValueError(f'{name} must be an integer between {lower} and {upper}')
    return int(value)


def _duration(value):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not math.isfinite(value) or value < 0):
        raise ValueError('model-coordinate duration must be finite and nonnegative')
    return float(value)


def _permutation(value):
    try:
        result = tuple(_integer(x, 'site permutation entry', 0, 7) for x in value)
    except TypeError as exc:
        raise ValueError('site permutation must contain each label exactly once') from exc
    if len(result) != 8 or set(result) != set(range(8)):
        raise ValueError('site permutation must contain each label exactly once')
    return result


def _parameters(value):
    if value is None:
        return ScalarParameters()
    if not isinstance(value, ScalarParameters):
        raise ValueError('parameters must be a ScalarParameters instance')
    return value


def _execute(stages, qubits, target):
    matrix, metadata = scalar_source._execute(stages, qubits, target)
    word = [event for _, stage in stages for event in stage]
    metadata['word_sha256'] = hashlib.sha256(json.dumps(
        word, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return matrix, metadata


@lru_cache(maxsize=6)
def _instrument(d, initial, other, sign):
    # Both a local four-level factor and a two-factor sixteen-level register
    # fit the compiler's complete real-basis limit.  Adding the binary pointer
    # needs at most five qubits, within its complete-unitary execution limit.
    basis = np.eye(d)
    basis[:, [0, initial]] = basis[:, [initial, 0]]
    other_column = int(np.argmax(basis[other]))
    first, second = basis[:, 0].copy(), basis[:, other_column].copy()
    basis[:, 0] = (first+second)/math.sqrt(2)
    basis[:, other_column] = (-first+second)/math.sqrt(2)
    phases = np.zeros(d)
    phases[other] = sign*math.pi/2
    change = c.compile_real_orthogonal(basis)+c.compile_diagonal(phases)
    pointer = d.bit_length()-1
    pointer_phases = np.zeros(2*d)
    pointer_phases[d] = math.pi
    hadamard = [c.pulse('H', math.pi/2, [pointer])]
    controlled_flip = hadamard+c.compile_diagonal(pointer_phases)+c.inverse(hadamard)
    vector = np.zeros(d, dtype=complex)
    vector[initial], vector[other] = 1/math.sqrt(2), sign*1j/math.sqrt(2)
    yes = np.outer(vector, vector.conj())
    no = np.eye(d)-yes
    matrix, metadata = _execute([
        ('pointer_inverse_probe_basis', c.inverse(change)),
        ('pointer_controlled_flip', controlled_flip),
        ('pointer_restore_probe_basis', change),
    ], d.bit_length(), np.block([[no, yes], [yes, no]]))
    kraus = (matrix[:d, :d], matrix[d:, :d])
    metadata['completeness_residual'] = float(np.linalg.norm(
        sum(k.conj().T@k for k in kraus)-np.eye(d), 2))
    return kraus, metadata


def instrument(d, initial, other, sign):
    """Retain both branches of (|initial>+sign*i|other>)/sqrt(2).

    For two four-level factors, little-endian index 5 denotes |1,1>.
    The complement includes all other occupations; none is postselected out.
    """
    if (isinstance(d, (bool, np.bool_)) or not isinstance(d, Integral)
            or d not in (4, 16)):
        raise ValueError('local probe register must have dimension four or sixteen')
    d = int(d)
    initial = _integer(initial, 'initial occupation', 0, d-1)
    other = _integer(other, 'other occupation', 0, d-1)
    sign = _integer(sign, 'projector sign', -1, 1)
    if initial == other or sign == 0:
        raise ValueError('distinct occupations and projector sign -1 or +1 required')
    kraus, metadata = _instrument(d, initial, other, sign)
    return tuple(k.copy() for k in kraus), deepcopy(metadata)


def reduced_density(state, dimensions, sites):
    """Exact spectator trace of a pure source state, without renormalization."""
    if (not isinstance(dimensions, (tuple, list)) or not dimensions
            or any(isinstance(d, (bool, np.bool_)) or not isinstance(d, Integral)
                   or d < 1 for d in dimensions)):
        raise ValueError('positive tensor-factor dimensions required')
    if (not isinstance(sites, (tuple, list)) or not sites
            or any(isinstance(i, (bool, np.bool_)) or not isinstance(i, Integral)
                   or not 0 <= i < len(dimensions) for i in sites)
            or len(set(sites)) != len(sites)):
        raise ValueError('distinct accessible tensor factors required')
    if np.ma.isMaskedArray(state):
        raise ValueError('state cannot be masked')
    state = np.asarray(state)
    if (state.ndim != 1 or state.shape[0] != math.prod(dimensions)
            or state.dtype.kind not in 'fiuc' or not np.isfinite(state).all()):
        raise ValueError('one finite amplitude per tensor occupation required')
    tensor = state.reshape(dimensions, order='F')
    local = np.moveaxis(tensor, sites, tuple(range(len(sites))))
    block = local.reshape((math.prod(dimensions[i] for i in sites), -1), order='F')
    return block@block.conj().T


def _preparations(lattice, permutation):
    gates = {n: preparation_source.preparation(4, n) for n in (0, 1)}
    columns, metadata = [], []
    for excited in (None, *range(8)):
        state = lattice.vacuum_state()
        requests = []
        for label, physical in enumerate(permutation):
            matrix, info = gates[int(label == excited)]
            state = lattice.apply_local(state, [physical], matrix)
            requests.append(info)
        columns.append(state)
        metadata.append({
            'resources': {key: sum(item['resources'][key] for item in requests) for key in COUNT_KEYS},
            'maximum_local_target_residual': max(item['full_phase_target_residual'] for item in requests),
        })
    return np.column_stack(columns), metadata


def _evolution_metadata(lattice, duration, steps):
    dt = duration/steps
    requests = [(2*steps, scalar_source.onsite(d, int(lattice.degree[i]), dt/2,
                                             lattice.parameters)[1])
                for i, d in enumerate(lattice.dimensions)]
    requests.extend((steps, scalar_source.bond(lattice.dimensions[i], lattice.dimensions[j],
                                               dt, lattice.parameters)[1]) for i, j in lattice.edges)
    return {
        'resources': {key: sum(count*info['resources'][key] for count, info in requests) for key in COUNT_KEYS},
        'maximum_local_target_residual': max(info['full_phase_target_residual'] for _, info in requests),
        'local_generators_per_trial': sum(count for count, _ in requests),
    }


def calibration_records(duration, parameters=None, *, permutation=PERMUTATION, steps=1):
    """Execute 88 fresh binary trials with the complete field coupled.

    Record axes match the isolated-generator predecessor: onsite axes are
    opaque site, transition (0->2,1->3), phase (-i,+i), outcome (0,1);
    pair axes are listed opaque pair, phase, outcome.  Unlike that predecessor,
    all onsite and bond generators participate in every trial.  The inverse
    must use a bound for the complete fixed schedule at both clock durations.
    """
    duration, parameters = _duration(duration), _parameters(parameters)
    permutation = _permutation(permutation)
    if isinstance(steps, (bool, np.bool_)) or not isinstance(steps, Integral) or steps < 1:
        raise ValueError('positive integer fixed schedule step count required')
    steps = int(steps)
    lattice = ScalarLattice(DIMENSIONS, parameters)
    prepared, preparations = _preparations(lattice, permutation)
    states = evolve(lattice, prepared, duration, steps)
    evolution = _evolution_metadata(lattice, duration, steps)
    onsite, bonds = np.zeros((8, 2, 2, 2)), np.zeros((len(PAIRS), 2, 2))
    resources = {key: 0 for key in COUNT_KEYS}
    maximum_residual = max(evolution['maximum_local_target_residual'],
                           *(item['maximum_local_target_residual'] for item in preparations))
    word_hashes = set()

    def read(column, sites, initial, other):
        nonlocal maximum_residual
        density = reduced_density(states[:, column], DIMENSIONS, sites)
        laws = []
        for sign in PHASE_SIGNS:
            kraus, metadata = instrument(len(density), initial, other, sign)
            laws.append([float(np.trace(k@density@k.conj().T).real) for k in kraus])
            for key in COUNT_KEYS:
                resources[key] += (preparations[column]['resources'][key]
                                   +evolution['resources'][key]+metadata['resources'][key])
            maximum_residual = max(maximum_residual, metadata['full_phase_target_residual'])
            word_hashes.add(metadata['word_sha256'])
        return laws

    for label, physical in enumerate(permutation):
        onsite[label, 0] = read(0, [physical], 0, 2)
        onsite[label, 1] = read(label+1, [physical], 1, 3)
    for index, (left, right) in enumerate(PAIRS):
        bonds[index] = read(0, [permutation[left], permutation[right]], 0, 5)
    for laws in (onsite, bonds):
        if (not np.isfinite(laws).all() or laws.min() < -2e-12
                or np.max(abs(laws.sum(axis=-1)-1)) > 2e-10):
            raise ArithmeticError('complete source detector laws lost normalization')
    return {
        'onsite_probabilities': onsite, 'bond_probabilities': bonds,
        'pairs': PAIRS, 'phase_signs': PHASE_SIGNS,
        'resources': resources,
        'evolution_resources_per_trial': evolution['resources'],
        'preparation_resources_per_trial': preparations[0]['resources'],
        'maximum_local_target_residual': maximum_residual,
        'pointer_word_sha256': sorted(word_hashes),
        'fresh_trials': 88, 'independent_preparation_columns_contracted': 9,
        'complete_pointer_outcomes': 2,
        'occupation_dimensions': DIMENSIONS,
        'fixed_symmetric_steps': steps,
        'local_generators_per_trial': evolution['local_generators_per_trial'],
        'complete_field_schedule_in_every_trial': True,
        'independently_isolated_action_services': False,
        'physical_service_error_certified': False,
    }
