"""Source-executed calibration probes of a hidden-label local scalar action.

The experimenter may isolate a named onsite service or a named pair service.
This is a supplied control capability, stronger than passive observation of
an uncontrolled field.  All preparations, evolution and binary instruments
below are contractions of emitted twelve-port source words.  Independent
finite targets check those words; they never replace executed matrices.

Only complete detector laws and opaque labels cross to the inverse problem.
The action, oscillator basis, code services and model-coordinate scheduling
are supplied.  The separate ``supplied_target`` accessor is for later audit.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
from itertools import combinations
import hashlib
import json
import math
from numbers import Integral, Real

import numpy as np

from ..observer.source_scalar import compiler as c
from ..relativistic_scalar import source as scalar_source
from ..relativistic_scalar.lattice import CUBE_EDGES, ScalarParameters


PERMUTATION = (5, 2, 7, 0, 6, 1, 4, 3)
PAIRS = tuple(combinations(range(8), 2))
PHASE_SIGNS = (-1, 1)
COUNT_KEYS = scalar_source.COUNT_KEYS


def _integer(value, name, lower, upper):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or not lower <= value <= upper):
        raise ValueError(f'{name} must be an integer between {lower} and {upper}')
    return int(value)


def _dimension(d):
    if isinstance(d, (bool, np.bool_)) or not isinstance(d, Integral) or d not in (2, 4, 8):
        raise ValueError('source occupation dimension must be 2, 4 or 8')
    return int(d)


def _duration(value):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not math.isfinite(value) or value < 0):
        raise ValueError('model-coordinate duration must be finite and nonnegative')
    return float(value)


def _parameters(value):
    if value is None:
        return ScalarParameters()
    if not isinstance(value, ScalarParameters):
        raise ValueError('parameters must be a ScalarParameters instance')
    return value


def _permutation(value):
    try:
        result = tuple(_integer(x, 'site permutation entry', 0, 7) for x in value)
    except TypeError as exc:
        raise ValueError('site permutation must contain each label exactly once') from exc
    if len(result) != 8 or set(result) != set(range(8)):
        raise ValueError('site permutation must contain each label exactly once')
    return result


def _edges(permutation):
    inverse = {original: label for label, original in enumerate(permutation)}
    return {tuple(sorted((inverse[a], inverse[b]))) for a, b in CUBE_EDGES}


def supplied_target(parameters=None, *, permutation=PERMUTATION):
    """Independent action/graph audit; never used to produce detector laws."""
    parameters, permutation = _parameters(parameters), _permutation(permutation)
    adjacency = np.zeros((8, 8), dtype=int)
    for a, b in _edges(permutation):
        adjacency[a, b] = adjacency[b, a] = 1
    diagonal = parameters.mass**2+parameters.kappa*adjacency.sum(axis=1)
    return {'adjacency': adjacency,
            'stiffness': np.diag(diagonal)-parameters.kappa*adjacency,
            'onsite_quadratic': diagonal,
            'mass_squared': parameters.mass**2,
            'quartic': parameters.g,
            'bond_coupling': parameters.kappa}


def _execute(stages, qubits, target):
    matrix, metadata = scalar_source._execute(stages, qubits, target)
    word = [event for _, stage in stages for event in stage]
    metadata['word_sha256'] = hashlib.sha256(json.dumps(
        word, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return matrix, metadata


@lru_cache(maxsize=14)
def _preparation(d, occupation):
    basis = np.eye(d)
    basis[:, [0, occupation]] = basis[:, [occupation, 0]]
    return _execute([('number_state_preparation', c.compile_real_orthogonal(basis))],
                    d.bit_length()-1, basis)


def preparation(d, occupation=0):
    """Execute a complete basis permutation taking a blank to |occupation>."""
    d = _dimension(d)
    occupation = _integer(occupation, 'occupation', 0, d-1)
    matrix, metadata = _preparation(d, occupation)
    return matrix.copy(), deepcopy(metadata)


def _pointer_word(d):
    pointer = d.bit_length()-1
    phases = np.zeros(2*d)
    phases[d] = math.pi
    hadamard = [c.pulse('H', math.pi/2, [pointer])]
    return hadamard+c.compile_diagonal(phases)+c.inverse(hadamard)


@lru_cache(maxsize=32)
def _instrument(d, initial, other, sign):
    # The first basis column is the requested real superposition.  A later
    # diagonal phase makes its second component imaginary.  The whole basis,
    # including every complement column, is compiled and executed.
    basis = np.eye(d)
    basis[:, [0, initial]] = basis[:, [initial, 0]]
    other_column = 0 if other == initial else int(np.argmax(basis[other]))
    cosine = math.sqrt(.5)
    first, second = basis[:, 0].copy(), basis[:, other_column].copy()
    basis[:, 0] = cosine*(first+second)
    basis[:, other_column] = cosine*(-first+second)
    phases = np.zeros(d)
    phases[other] = sign*math.pi/2
    change = c.compile_real_orthogonal(basis)+c.compile_diagonal(phases)
    vector = np.zeros(d, dtype=complex)
    vector[initial], vector[other] = cosine, sign*1j*cosine
    yes = np.outer(vector, vector.conj())
    no = np.eye(d)-yes
    target = np.block([[no, yes], [yes, no]])
    matrix, metadata = _execute([
        ('pointer_inverse_probe_basis', c.inverse(change)),
        ('pointer_controlled_flip', _pointer_word(d)),
        ('pointer_restore_probe_basis', change),
    ], d.bit_length(), target)
    kraus = (matrix[:d, :d], matrix[d:, :d])
    metadata['completeness_residual'] = float(np.linalg.norm(
        sum(k.conj().T@k for k in kraus)-np.eye(d), 2))
    return kraus, metadata


def instrument(d, initial, other, sign):
    """Both unnormalized branches; bit 1 reads (|initial>+sign*i|other>)/sqrt2."""
    d = _dimension(d)
    initial = _integer(initial, 'initial occupation', 0, d-1)
    other = _integer(other, 'other occupation', 0, d-1)
    sign = _integer(sign, 'projector sign', -1, 1)
    if initial == other or sign == 0:
        raise ValueError('two distinct occupations and projector sign -1 or +1 required')
    kraus, metadata = _instrument(d, initial, other, sign)
    return tuple(k.copy() for k in kraus), deepcopy(metadata)


def _summary(gates, trials):
    return {
        'resources': {key: sum(gate['resources'][key] for gate in gates) for key in COUNT_KEYS},
        'maximum_local_target_residual': max(gate['full_phase_target_residual'] for gate in gates),
        'fresh_trials': trials,
        'complete_pointer_outcomes': 2,
        'physical_service_error_certified': False,
    }


def _probe(unitary, evolution_metadata, initial, other):
    prepared, preparation_metadata = preparation(len(unitary), initial)
    state = unitary@prepared[:, 0]
    probabilities, gates = [], []
    for sign in PHASE_SIGNS:
        kraus, pointer_metadata = instrument(len(unitary), initial, other, sign)
        probabilities.append([float(np.vdot(k@state, k@state).real) for k in kraus])
        gates.extend([preparation_metadata, evolution_metadata, pointer_metadata])
    probabilities = np.array(probabilities)
    if (not np.isfinite(probabilities).all() or probabilities.min() < 0
            or np.max(abs(probabilities.sum(axis=1)-1)) > 2e-10):
        raise ArithmeticError('source detector law lost normalization')
    return probabilities, _summary(gates, 2)


def onsite_probe(site, transition, duration, parameters=None, *, permutation=PERMUTATION):
    """Two fresh phase reads after isolating one four-level onsite service.

    Transition 0 prepares |0> and reads its coherence with |2>; transition 1
    prepares |1> and reads |3>.  Phase order is (-i,+i).  Therefore
    p1(-i)-p1(+i) = 2*t*h[other,initial] + O(t^3) for real symmetric h.
    The full onsite diagonal from the gradient action remains during isolation.
    """
    site = _integer(site, 'site label', 0, 7)
    transition = _integer(transition, 'onsite transition', 0, 1)
    duration, parameters = _duration(duration), _parameters(parameters)
    edges = _edges(_permutation(permutation))
    degree = sum(site in edge for edge in edges)
    unitary, metadata = scalar_source.onsite(4, degree, duration, parameters)
    return _probe(unitary, metadata, transition, transition+2)


def bond_probe(pair, duration, parameters=None, *, permutation=PERMUTATION):
    """Isolate an opaque pair, prepare |00>, and read |00> +/- i|11>.

    Each site has two occupation levels.  Missing pairs execute the same bond
    compilation schedule with zero duration.  Edge presence is disclosed by
    the detector law, not by a list of successful service requests.
    """
    try:
        pair = tuple(_integer(x, 'pair site label', 0, 7) for x in pair)
    except TypeError as exc:
        raise ValueError('pair must name two distinct site labels') from exc
    if len(pair) != 2 or pair[0] == pair[1]:
        raise ValueError('pair must name two distinct site labels')
    duration, parameters = _duration(duration), _parameters(parameters)
    edges = _edges(_permutation(permutation))
    unitary, metadata = scalar_source.bond(2, 2,
        duration if tuple(sorted(pair)) in edges else 0., parameters)
    return _probe(unitary, metadata, 0, 3)


def calibration_records(duration, parameters=None, *, permutation=PERMUTATION):
    """Execute all 88 fresh trials at one supplied model-coordinate duration.

    ``onsite_probabilities`` axes are site, transition (0->2,1->3), phase
    (-i,+i), outcome (0,1).  ``bond_probabilities`` axes are the 28 listed
    opaque pairs, phase, outcome.  The duration, action coefficients and
    graph are deliberately absent from returned records.  A separate clock
    service is required by the inverse problem.
    """
    duration, parameters = _duration(duration), _parameters(parameters)
    permutation = _permutation(permutation)
    onsite = np.zeros((8, 2, 2, 2))
    bonds = np.zeros((len(PAIRS), 2, 2))
    summaries = []
    for site in range(8):
        for transition in range(2):
            onsite[site, transition], metadata = onsite_probe(
                site, transition, duration, parameters, permutation=permutation)
            summaries.append(metadata)
    for index, pair in enumerate(PAIRS):
        bonds[index], metadata = bond_probe(pair, duration, parameters, permutation=permutation)
        summaries.append(metadata)
    return {
        'onsite_probabilities': onsite,
        'bond_probabilities': bonds,
        'pairs': PAIRS,
        'phase_signs': PHASE_SIGNS,
        'resources': {key: sum(item['resources'][key] for item in summaries) for key in COUNT_KEYS},
        'maximum_local_target_residual': max(item['maximum_local_target_residual'] for item in summaries),
        'fresh_trials': sum(item['fresh_trials'] for item in summaries),
        'complete_pointer_outcomes': 2,
        'independently_isolated_action_services': True,
        'physical_service_error_certified': False,
    }
