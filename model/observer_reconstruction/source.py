"""Source-executed probes of a hidden-label finite canonical free field.

The backend supplies the action, normal-mode vacuum, coherent controls and
clock. Only probabilities cross its interface to the inverse problem.
Normal-mode occupation truncation is distinct from the existing local-site
regulator. Partial tracing a sequential common-pointer circuit contracts its
full local source matrices; no analytic Gaussian target supplies a probability.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
import json
import math

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar import compiler as c
from ..relativistic_scalar.lattice import CUBE_EDGES, local_operators
from ..scalar_observations.experimental_readout import canonical_source as canonical


PERMUTATION = (5, 2, 7, 0, 6, 1, 4, 3)
MASS = 1.25
KAPPA = .75
COUNT_KEYS = ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
              'identity_phase_pulses', 'proper_code_load_events',
              'proper_code_export_events', 'proper_code_service_events')


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(f'{name} must be a finite real number')
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f'{name} must be a finite real number')
    return result


def _dimension(d):
    if type(d) is not int or d not in (8, 16):
        raise ValueError('the retained normal-mode source codes have dimension 8 or 16')
    return d


def _site(site):
    if type(site) is not int or not 0 <= site < 8:
        raise ValueError('site label must be an integer from zero through seven')
    return site


def supplied_target():
    """Independent audit target; never passed to the reconstruction routine."""
    laplacian = np.zeros((8, 8))
    for a, b in CUBE_EDGES:
        laplacian[a, a] += 1
        laplacian[b, b] += 1
        laplacian[a, b] -= 1
        laplacian[b, a] -= 1
    return (MASS*MASS*np.eye(8)+KAPPA*laplacian)[np.ix_(PERMUTATION, PERMUTATION)]


def modes():
    # Analytic Walsh basis fixes degenerate eigenspaces and column signs.
    basis = np.array([[(-1)**((site & mode).bit_count())/math.sqrt(8)
                       for mode in range(8)] for site in PERMUTATION])
    frequencies = np.array([math.sqrt(MASS*MASS+2*KAPPA*j.bit_count()) for j in range(8)])
    return frequencies, basis


def _basis(matrix):
    values, basis = np.linalg.eigh(matrix)
    for column in basis.T:
        pivot = int(np.flatnonzero(abs(column) > 1e-12)[0])
        if column[pivot] < 0:
            column *= -1
    return values, basis


def _execute(word, bits, target):
    raw = c.execute(word, bits)
    residual = float(np.linalg.norm(raw-target, 2))
    if not math.isfinite(residual) or residual > 2e-10:
        raise ArithmeticError('source word failed its independent full-phase target check')
    matrix, stability = c.stabilize_executed_unitary(raw)
    matrix.flags.writeable = False
    resources = c.resources(word)
    return matrix, {
        'word_sha256': hashlib.sha256(json.dumps(word, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
        'resources': {key: resources[key] for key in COUNT_KEYS},
        'full_phase_target_residual': residual,
        'numerical_stabilization': stability,
    }


@lru_cache(maxsize=2)
def _coordinate_basis_word(d):
    return canonical.coordinate_word(d)


@lru_cache(maxsize=256)
def _coordinate(d, omega, angle, controlled):
    q = local_operators(d, omega)['q']
    # Reuse the portable fixed-arithmetic Gray/Walsh angle producer. This
    # selects a reproducible requested word, not its executed matrix.
    context = canonical._context()
    values, _ = canonical._coordinate_data(d)
    factor = (canonical._parameter(context, angle, 'angle')
              *context.sqrt(2/canonical._parameter(context, omega, 'frequency')))
    change = _coordinate_basis_word(d)
    phases = [factor*value for value in values]
    if controlled:
        phases = [context.zero]*d+phases
    word = c.inverse(change)+canonical._diagonal_word(context, phases)+change
    target = expm(1j*angle*q)
    if controlled:
        target = np.block([[np.eye(d), np.zeros((d, d))], [np.zeros((d, d)), target]])
    return _execute(word, d.bit_length()-1+int(controlled), target)


def coordinate(d, omega, angle, *, controlled=False):
    d = _dimension(d)
    omega, angle = _real(omega, 'basis frequency'), _real(angle, 'coordinate angle')
    if omega <= 0 or type(controlled) is not bool:
        raise ValueError('positive basis frequency and boolean control required')
    matrix, info = _coordinate(d, omega, angle, controlled)
    return matrix.copy(), deepcopy(info)


@lru_cache(maxsize=256)
def _evolution(d, omega, duration, anharmonic):
    h = np.diag(omega*np.arange(d))
    if anharmonic:
        # Deliberate model-change control: quartic in NORMAL MODE coordinates.
        # This is not a local site phi^4 interaction or a continuum claim.
        h = h+anharmonic*local_operators(d, omega)['q4']
        values, basis = _basis(h)
        change = c.compile_real_orthogonal(basis)
        word = c.inverse(change)+c.compile_diagonal(-duration*values)+change
    else:
        word = c.compile_diagonal(-duration*np.diag(h))
    return _execute(word, d.bit_length()-1, expm(-1j*duration*h))


def evolution(d, omega, duration, *, anharmonic=0.):
    d = _dimension(d)
    omega, duration = _real(omega, 'basis frequency'), _real(duration, 'duration')
    anharmonic = _real(anharmonic, 'anharmonic control')
    if omega <= 0 or duration < 0 or anharmonic < 0:
        raise ValueError('positive frequency and nonnegative duration/control required')
    matrix, info = _evolution(d, omega, duration, anharmonic)
    return matrix.copy(), deepcopy(info)


@lru_cache(maxsize=16)
def _pointer_basis(phase):
    word = c.compile_diagonal([0., phase])+c.inverse([c.pulse('H', math.pi/2, [0])])
    h = np.array([[1., 1.], [1., -1.]])/math.sqrt(2)
    return _execute(word, 1, -1j*h@np.diag([1., np.exp(1j*phase)]))


@lru_cache(maxsize=1)
def _pointer_preparation():
    word = [c.pulse('H', math.pi/2, [0])]
    h = np.array([[1., 1.], [1., -1.]])/math.sqrt(2)
    return _execute(word, 1, 1j*h)


def _transfer(rho, state, controlled):
    """Trace one mode after its actual common-pointer controlled source word.

    All four pointer blocks are used, including floating off-block residuals.
    This is the exact tensor contraction of a star circuit, not a replacement
    of the source matrix with its analytic target or diagonal blocks.
    """
    d = len(state)
    vectors = np.einsum('aibj,j->abi', controlled.reshape(2, d, 2, d), state)
    return np.einsum('bc,abi,dci->ad', rho, vectors, vectors.conj())


def pointer_probabilities(states, controls, phase):
    if len(states) != len(controls) or not states:
        raise ValueError('one controlled source matrix per prepared mode required')
    prep, _ = _pointer_preparation()
    rho = np.outer(prep[:, 0], prep[:, 0].conj())
    for state, control in zip(states, controls, strict=True):
        rho = _transfer(rho, state, control)
    read, _ = _pointer_basis(_real(phase, 'pointer phase'))
    rho = read@rho@read.conj().T
    probabilities = np.diag(rho).real
    if (not np.isfinite(probabilities).all() or probabilities.min() < -2e-12
            or abs(probabilities.sum()-1) > 2e-10):
        raise ArithmeticError('complete pointer law lost normalization')
    return probabilities


def field_readout(d, duration, kick_site, read_site, gamma, eta, *, anharmonic=0.):
    """Two Ramsey phases, both outcomes, on fresh copies of the same protocol."""
    _dimension(d)
    kick_site, read_site = _site(kick_site), _site(read_site)
    gamma, eta = _real(gamma, 'kick'), _real(eta, 'probe')
    frequencies, basis = modes()
    states, controls, gates = [], [], []
    for mode, omega in enumerate(frequencies):
        kick, ki = coordinate(d, omega, gamma*basis[kick_site, mode])
        free, fi = evolution(d, omega, duration, anharmonic=anharmonic)
        control, ci = coordinate(d, omega, eta*basis[read_site, mode], controlled=True)
        states.append(free@kick[:, 0])
        controls.append(control)
        gates.extend([ki, fi, ci])
    probabilities = np.array([pointer_probabilities(states, controls, phase)
                              for phase in (0., math.pi/2)])
    prep = _pointer_preparation()[1]
    reads = [_pointer_basis(phase)[1] for phase in (0., math.pi/2)]
    resources = {key: 2*(sum(g['resources'][key] for g in gates)+prep['resources'][key])
                 + sum(g['resources'][key] for g in reads) for key in COUNT_KEYS}
    return probabilities, {'resources_for_two_fresh_phase_trials': resources,
        'maximum_local_target_residual': max(g['full_phase_target_residual'] for g in [*gates, prep, *reads]),
        'working_qubits_per_trial': 8*(d.bit_length()-1)+1,
        'complete_pointer_outcomes': 2,
        'source_words': sorted({g['word_sha256'] for g in [*gates, prep, *reads]})}


def response_records(d, duration, gamma=.3, eta=.4):
    records = np.empty((2, 8, 8, 2, 2))
    resources = {key: 0 for key in COUNT_KEYS}
    residual = 0.
    words = set()
    for sign_index, sign in enumerate((1, -1)):
        for x in range(8):
            for y in range(8):
                row, meta = field_readout(d, duration, x, y, sign*gamma, eta)
                records[sign_index, x, y] = row
                for key in COUNT_KEYS:
                    resources[key] += meta['resources_for_two_fresh_phase_trials'][key]
                residual = max(residual, meta['maximum_local_target_residual'])
                words.update(meta['source_words'])
    return {'probabilities': records, 'resources': resources,
            'maximum_local_target_residual': residual, 'source_words': sorted(words),
            'fresh_phase_trials': 256}


def clock_readout(duration):
    """Independent two-level phase clock; its declared positive gap sets units.

    We choose the orientation exp(+i t N). p0(phi)=(1+cos(t+phi))/2.
    The inverse receives probabilities, never this host-supplied duration.
    """
    duration = _real(duration, 'clock model duration')
    if duration < 0:
        raise ValueError('nonnegative clock duration required')
    probabilities, metadata = [], []
    for phase in (0., math.pi/2):
        word = ([c.pulse('H', math.pi/2, [0])]
                +c.compile_diagonal([0., duration+phase])
                +c.inverse([c.pulse('H', math.pi/2, [0])]))
        h = np.array([[1., 1.], [1., -1.]])/math.sqrt(2)
        matrix, info = _execute(word, 1, h@np.diag([1., np.exp(1j*(duration+phase))])@h)
        probabilities.append(abs(matrix[:, 0])**2)
        metadata.append(info)
    return np.array(probabilities), metadata


def weyl_loop(d=16, *, reverse=False, a=.3, b=.4):
    """Ordered Q/P loop on the same oscillator source; both pointer outcomes."""
    _dimension(d)
    if type(reverse) is not bool:
        raise ValueError('boolean loop orientation required')
    a, b = _real(a, 'Q displacement'), _real(b, 'P displacement')
    omega = MASS
    qplus, _ = coordinate(d, omega, a, controlled=True)
    qminus, _ = coordinate(d, omega, -a, controlled=True)
    pplus_q, _ = coordinate(d, omega, b*omega, controlled=True)
    pminus_q, _ = coordinate(d, omega, -b*omega, controlled=True)
    rotation, _ = evolution(d, omega, math.pi/(2*omega))
    phases = math.pi*np.arange(d)/2
    inverse_rotation, _ = _execute(c.compile_diagonal(phases), d.bit_length()-1,
                                   np.diag(np.exp(1j*phases)))
    # R* Q R=P/omega, with R=exp(-i pi N/2).
    lift = np.kron(np.eye(2), rotation)
    inverse_lift = np.kron(np.eye(2), inverse_rotation)
    pplus = inverse_lift@pplus_q@lift
    pminus = inverse_lift@pminus_q@lift
    factors = [qplus, pplus, qminus, pminus]
    if reverse:
        factors = [pplus, qplus, pminus, qminus]
    loop = np.eye(2*d, dtype=complex)
    for matrix in factors:
        loop = matrix@loop
    blank = np.eye(d)[:, 0]
    return np.array([pointer_probabilities([blank], [loop], phase)
                     for phase in (0., math.pi/2)])
