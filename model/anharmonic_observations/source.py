"""Native source execution of the published continuously detuned pulse schedule.

The supplied apparatus uses a truncated oscillator and a two-level spin.  Its
Schrodinger-frame pulse Hamiltonian is delta*N - sign*Omega*X*(a+a†)/2.
All unitary matrices below are contractions of twelve-port source words;
independent SciPy matrices only check those contractions.  Phase diffusion is
an explicitly supplied channel, split symmetrically around source unitaries.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import math
from numbers import Integral

import numpy as np
from scipy.linalg import expm

from ..observer.source_scalar import compiler as c
from ..scalar_observations.source import _basis
from ..scalar_observations.experimental_readout import canonical_source as canonical
from ..scalar_observations.experimental_readout import executable_source as native


X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])


def _integer(value, name, minimum=0):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f'{name} must be an integer at least {minimum}')
    return int(value)


def _counts(stages):
    return {key: sum(stage['resources'][key] for stage in stages) for key in native.COUNT_KEYS}


@lru_cache(maxsize=16)
def _forced_basis(dimension, ratio):
    """Canonical real eigenbasis of N-ratio*(a+a†)/2, executed both ways."""
    context = canonical._context()
    ratio_mp = canonical._parameter(context, ratio, 'force/frequency ratio')
    ladder = canonical._ladder(context, dimension)
    hamiltonian = context.diag(list(range(dimension))) - ratio_mp*(ladder+ladder.T)/2
    values_mp, vectors = context.eigsy(hamiltonian)
    for column in range(dimension):
        pivot = next(row for row in range(dimension) if abs(vectors[row, column]) > context.mpf('1e-12'))
        if vectors[pivot, column] < 0:
            for row in range(dimension):
                vectors[row, column] = -vectors[row, column]
    word = canonical._orthogonal_word(context, vectors)
    a = native._ladder(dimension)
    _, target = _basis(np.diag(np.arange(dimension))-ratio*(a+a.T)/2)
    bits = dimension.bit_length()-1
    basis, metadata = native._contract(word, bits, target, name='forced_oscillator_basis')
    inverse, inverse_metadata = native._contract(c.inverse(word), bits, target.T,
                                               name='forced_oscillator_inverse_basis')
    return tuple(context.nstr(value, canonical.DECIMAL_DIGITS) for value in values_mp), basis, inverse, metadata, inverse_metadata


@lru_cache(maxsize=128)
def spin_rotation(axis, angle):
    """Execute exp(-i*angle*Pauli/2) from a native two-level source pulse."""
    if axis not in ('X', 'Y', 'Z'):
        raise ValueError('spin axis must be X, Y or Z')
    angle = native._real(angle, 'spin rotation angle')
    word = [c.pulse(axis, -angle/2, [0])]
    return native._contract(word, 1, expm(-.5j*angle*dict(X=X, Y=Y, Z=Z)[axis]),
                            name=f'spin_rotation_{axis}')


@lru_cache(maxsize=8)
def _controlled_parity(dimension):
    # Occupation parity is the least significant oscillator bit. Retain the
    # full scalar phase: the native J(-pi) realizes minus controlled Z.
    dimension = native._dimension(dimension)
    spin_bit = dimension.bit_length()-1
    phases = np.r_[np.zeros(dimension), math.pi*(np.arange(dimension) % 2)]
    word = [c.pulse('J', -math.pi, [spin_bit, 0]),
            c.pulse('I', math.pi, [0])]
    return native._contract(word, dimension.bit_length(), np.diag(np.exp(1j*phases)),
                            name='spin_controlled_oscillator_parity')


@lru_cache(maxsize=128)
def pulse_unitary(dimension, delta, omega, duration, sign=1):
    """Execute the exact finite, continuously detuned Hamiltonian pulse.

    The canonical eigenbasis is shared across signs and durations.  Reusing a
    contracted matrix is host arithmetic reuse; its full word cost is counted
    on every physical application.  No displacement/CCR identity at the cutoff
    is assumed.
    """
    dimension = native._dimension(dimension)
    delta, omega, duration = (native._real(value, name) for value, name in
                              ((delta, 'delta'), (omega, 'omega'), (duration, 'duration')))
    if delta <= 0 or omega < 0 or duration < 0 or type(sign) is not int or sign not in (-1, 1):
        raise ValueError('positive delta, nonnegative omega/duration and sign +/-1 required')
    values, basis, inverse, bmeta, imeta = _forced_basis(dimension, omega/delta)
    context = canonical._context()
    delta_mp = canonical._parameter(context, delta, 'delta')
    duration_mp = canonical._parameter(context, duration, 'duration')
    phases_mp = [-delta_mp*duration_mp*context.mpf(value) for value in values]
    word = canonical._diagonal_word(context, phases_mp)
    phases = np.array([float(value) for value in phases_mp])
    diagonal, dmeta = native._contract(word, dimension.bit_length()-1,
                                       np.diag(np.exp(1j*phases)), name='forced_eigenphases')
    axis, ameta = spin_rotation('Y', math.pi/2)
    inverse_axis, aimeta = spin_rotation('Y', -math.pi/2)
    parity, pmeta = _controlled_parity(dimension)
    # In the spin-X eigenbasis the negative-force sector is a parity conjugate.
    single = basis @ diagonal @ inverse
    central = parity @ np.kron(np.eye(2), single) @ parity
    left_axis = np.kron(axis, np.eye(dimension))
    right_axis = np.kron(inverse_axis, np.eye(dimension))
    raw = left_axis @ central @ right_axis
    # Evaluate the same ideal inverse-word conjugations as identity plus an
    # increment. Independently contracted B and B^-1 otherwise contribute a
    # tiny, duration-independent B B^-1-I at every short subpulse. Both actual
    # contractions remain in this expression and in the raw mandatory check.
    # Parity is applied twice, not followed by its inverse: retain P^2-I,
    # including its literal binary-angle phase and contraction roundoff.
    identity = np.eye(2*dimension)
    increment = basis @ (diagonal-np.eye(dimension)) @ inverse
    increment = (parity @ parity-identity
                 + parity @ np.kron(np.eye(2), increment) @ parity)
    increment = left_axis @ increment @ right_axis
    stages = [aimeta, pmeta, imeta, dmeta, bmeta, pmeta, ameta]
    if sign == -1:
        flip, fmeta = spin_rotation('Z', math.pi)
        inverse_flip, fimeta = spin_rotation('Z', -math.pi)
        left_flip = np.kron(flip, np.eye(dimension))
        right_flip = np.kron(inverse_flip, np.eye(dimension))
        raw = left_flip @ raw @ right_flip
        increment = left_flip @ increment @ right_flip
        stages = [fimeta] + stages + [fmeta]
    matrix = identity+increment
    a = native._ladder(dimension)
    hamiltonian = delta*np.kron(np.eye(2), np.diag(np.arange(dimension))) - sign*omega/2*np.kron(X, a+a.T)
    target = expm(-1j*duration*hamiltonian)
    if not all(np.all(np.isfinite(value)) for value in (raw, matrix, target)):
        raise ArithmeticError('source pulse has nonfinite contraction data')
    raw_residual = float(np.linalg.norm(raw-target, 2))
    if not math.isfinite(raw_residual) or raw_residual > 3e-10:
        raise ArithmeticError('source pulse raw contraction failed its complete finite Hamiltonian check')
    change = float(np.linalg.norm(matrix-raw, 2))
    # The same fixed numerical scale as native word contraction, not a limit
    # fitted to the observed trajectory or to this pulse's target residual.
    change_limit = 2e-10
    if not math.isfinite(change) or change > change_limit:
        raise ArithmeticError('source pulse inverse-word representation change exceeds its fixed guard')
    increment_residual = float(np.linalg.norm(matrix-target, 2))
    if not math.isfinite(increment_residual) or increment_residual > 3e-10:
        raise ArithmeticError('source pulse increment failed its complete finite Hamiltonian check')
    # Identity-plus-increment arithmetic removes duration-independent inverse
    # drift, but independently contracted inverse matrices need not be exact
    # adjoints. Restore the ideal pulse's unitary representation with the same
    # fixed polar guard used for every native word. The raw source contraction
    # and cancellation change have already passed their independent guards;
    # neither a target matrix nor an evolved state enters this correction.
    matrix, stabilization = c.stabilize_executed_unitary(matrix)
    if not np.all(np.isfinite(matrix)):
        raise ArithmeticError('source pulse has nonfinite stabilized contraction data')
    complete_change = float(np.linalg.norm(matrix-raw, 2))
    if not math.isfinite(complete_change) or complete_change > change_limit:
        raise ArithmeticError('source pulse complete representation change exceeds its fixed guard')
    residual = float(np.linalg.norm(matrix-target, 2))
    if not math.isfinite(residual) or residual > 3e-10:
        raise ArithmeticError('source pulse failed its complete finite Hamiltonian check')
    defect = float(np.linalg.norm(matrix.conj().T @ matrix-identity, 2))
    if not math.isfinite(defect):
        raise ArithmeticError('source pulse has nonfinite unitarity diagnostic')
    matrix.flags.writeable = False
    return matrix, {'name': 'continuous_detuned_pulse', 'dimension': dimension,
                    'delta': delta, 'omega': omega, 'duration_seconds': duration, 'sign': sign,
                    'resources': _counts(stages), 'stages': deepcopy(stages),
                    'full_phase_target_residual': residual,
                    'inverse_word_cancellation': {
                        'method': 'identity plus executed conjugation increments, then guarded polar representation',
                        'raw_full_phase_target_residual': raw_residual,
                        'increment_full_phase_target_residual': increment_residual,
                        'representation_change': change,
                        'representation_change_limit': change_limit,
                        'complete_representation_change': complete_change,
                        'post_cancellation_unitary_stabilization': stabilization,
                        'represented_unitarity_defect': defect,
                        'literal_parity_square_retained': True,
                        'all_original_source_stages_and_resources_retained': True,
                        'physical_inverse_service_errors_cancelled': False,
                        'outward_roundoff_or_physical_error_certificate': False,
                    },
                    'matrix_origin': 'contracted_native_source_pulses',
                    'cutoff_CCR_identity_used': False}


def phase_diffusion(rho, dimension, rate, duration):
    """Supplied exp(duration*D[sqrt(rate)*N]), including spin coherences."""
    rate, duration = native._real(rate, 'dephasing rate'), native._real(duration, 'dephasing duration')
    if rate < 0 or duration < 0:
        raise ValueError('nonnegative dephasing rate and duration required')
    numbers = np.tile(np.arange(dimension), 2)
    return rho*np.exp(-rate*duration*(numbers[:, None]-numbers[None, :])**2/2)


@lru_cache(maxsize=512)
def detector(dimension, u, v, theta):
    """Return both source-executed oscillator Kraus operators, without selection."""
    dimension = native._dimension(dimension)
    setting = native._setting({'u': u, 'v': v, 'theta': theta})
    u, v, theta = (setting[key] for key in ('u', 'v', 'theta'))
    values, basis, inverse, bmeta, imeta = native._coordinate_basis(dimension)
    _, _, rotation_word, conditional_word = canonical.quadrature_words(dimension, u, v, theta)
    phi = 0. if math.hypot(u, v) == 0 else math.atan2(-u, v)
    rotation, rmeta = native._contract(rotation_word, dimension.bit_length()-1,
                                      np.diag(np.exp(-1j*phi*np.arange(dimension))), name='readout_quadrature_rotation')
    angles = math.hypot(u, v)*values+theta/2
    conditional, cmeta = native._contract(conditional_word, dimension.bit_length(),
                                         np.diag(np.exp(1j*np.r_[-angles, angles])), name='readout_conditional_phase')
    axis, ameta = spin_rotation('X', -math.pi/2)
    inverse_axis, aimeta = spin_rotation('X', math.pi/2)
    pointer = np.kron(axis, np.eye(dimension)) @ conditional @ np.kron(inverse_axis, np.eye(dimension))
    kraus = tuple(basis @ pointer[index*dimension:(index+1)*dimension, :dimension] @ inverse @ rotation
                  for index in (0, 1))
    for value in kraus:
        value.flags.writeable = False
    stages = [rmeta, imeta, aimeta, cmeta, ameta, bmeta]
    return kraus, {'setting': setting, 'resources': _counts(stages), 'stages': deepcopy(stages)}


def run(dimension, alpha_real, thermal, step_counts, settings, *,
        delta=2*math.pi*500, omega=math.pi/150e-6, pulse_duration=50e-6,
        theta=.8, phi=0., dephasing_rate=18., subdivisions=1):
    """Execute preparation, full spin-motion evolution and selected readout.

    No observed values are accepted.  Step k denotes simulated elapsed time
    4*k*pulse_duration.  The field is represented in the simulated Schrodinger
    frame.  A supplied phase-diffusion channel brackets each executed subpulse;
    varying subdivisions checks this channel/unitary splitting numerically.
    The spin is retained coherently through every step, then both final spin
    probabilities are kept.  Conditional detector probabilities explicitly
    divide the retained spin-zero joint branch by its acceptance probability.
    """
    dimension = native._dimension(dimension)
    if not isinstance(step_counts, (list, tuple)) or not step_counts:
        raise ValueError('nonempty list of step counts required')
    steps = [_integer(value, 'step count') for value in step_counts]
    if steps != sorted(set(steps)):
        raise ValueError('step counts must be sorted and unique')
    subdivisions = _integer(subdivisions, 'subdivisions', 1)
    if not isinstance(settings, (list, tuple)) or not settings:
        raise ValueError('nonempty list of detector settings required')
    settings = [native._setting(setting) for setting in settings]
    theta, phi, dephasing_rate = (native._real(value, name) for value, name in
                                ((theta, 'theta'), (phi, 'phi'), (dephasing_rate, 'dephasing rate')))
    if dephasing_rate < 0:
        raise ValueError('nonnegative dephasing rate required')
    states, weights, prepmeta = native.preparation(dimension, alpha_real, thermal)
    oscillator = (states*weights) @ states.conj().T
    rho = np.zeros((2*dimension, 2*dimension), complex)
    rho[:dimension, :dimension] = oscillator
    duration = native._real(pulse_duration, 'pulse duration')/subdivisions
    if duration <= 0:
        raise ValueError('positive pulse duration required')
    positive, pmeta = pulse_unitary(dimension, delta, omega, duration, 1)
    negative, nmeta = pulse_unitary(dimension, delta, omega, duration, -1)
    zrotation, zmeta = spin_rotation('Z', theta)
    xrotation, xmeta = spin_rotation('X', phi)
    xinverse, ximeta = spin_rotation('X', -phi)
    # The two axes agree with R_{-phi}(theta) then R_{phi}(theta).
    first_spin = xrotation @ zrotation @ xinverse
    second_spin = xinverse @ zrotation @ xrotation
    first = np.kron(first_spin, np.eye(dimension))
    second = np.kron(second_spin, np.eye(dimension))
    records = []
    detector_operations = [detector(dimension, **setting) for setting in settings]
    a = native._ladder(dimension)
    q = (a+a.T)/math.sqrt(2)
    for step in range(steps[-1]+1):
        if step in steps:
            blocks = [rho[index*dimension:(index+1)*dimension, index*dimension:(index+1)*dimension]
                      for index in (0, 1)]
            spin_probabilities = [float(np.trace(block).real) for block in blocks]
            if spin_probabilities[0] <= 0:
                raise ArithmeticError('zero selected spin probability')
            reads = []
            for setting, (kraus, metadata) in zip(settings, detector_operations, strict=True):
                joint = [[float(np.trace(operator @ block @ operator.conj().T).real)
                          for operator in kraus] for block in blocks]
                conditional = np.array(joint[0])/spin_probabilities[0]
                reads.append({**setting, 'spin_detector_joint_probabilities': joint,
                              'selected_detector_probabilities': conditional.tolist(),
                              'selected_mean': float(conditional[0]-conditional[1]),
                              'resources': deepcopy(metadata['resources'])})
            records.append({'steps': step, 'simulated_duration_seconds': 4*step*pulse_duration,
                            'spin_probabilities': spin_probabilities,
                            'total_probability_mass': float(np.trace(rho).real),
                            'selected_position_mean': float(np.trace(q @ blocks[0]).real/spin_probabilities[0]),
                            'reads': reads})
        if step == steps[-1]:
            break
        for unitary, count, spin_after in ((negative, subdivisions, first),
                                            (positive, 2*subdivisions, second),
                                            (negative, subdivisions, None)):
            for _ in range(count):
                rho = phase_diffusion(rho, dimension, dephasing_rate, duration/2)
                rho = unitary @ rho @ unitary.conj().T
                rho = phase_diffusion(rho, dimension, dephasing_rate, duration/2)
            if spin_after is not None:
                rho = spin_after @ rho @ spin_after.conj().T
    perstep = _counts([nmeta]* (2*subdivisions)+[pmeta]*(2*subdivisions)+[zmeta]*2+[xmeta, ximeta]*2)
    return {'schema': 'source_executed_anharmonic_pulses_v1', 'dimension': dimension,
            'preparation': prepmeta, 'records': records,
            'operations': {'positive_subpulse': deepcopy(pmeta), 'negative_subpulse': deepcopy(nmeta),
                           'spin_z': deepcopy(zmeta), 'spin_x': deepcopy(xmeta), 'spin_x_inverse': deepcopy(ximeta),
                           'detectors': [deepcopy(metadata) for _, metadata in detector_operations]},
            'resources_per_evolution_step': perstep,
            'settings': {'delta': delta, 'omega': omega, 'pulse_duration_seconds': pulse_duration,
                         'theta': theta, 'phi': phi, 'dephasing_rate_per_second': dephasing_rate,
                         'subdivisions_per_pulse': subdivisions},
            'time_parameter': 'supplied simulated Schrodinger-frame duration in seconds',
            'final_spin_postselection': True, 'both_final_spin_branches_retained': True,
            'intermediate_spin_projection_or_reset': False,
            'dephasing_channel_supplied_not_source_derived': True,
            'dephasing_splitting': 'half diffusion, executed Hamiltonian subpulse, half diffusion',
            'cutoff_error_certified_here': False, 'physical_service_error_certified': False,
            'coordinate_phase_convention': 'phi_m=0; x cosine; down=Z+; virtual Rz(+theta)',
            'prediction_origin': 'source_executed_unitaries_and_explicit_supplied_phase_diffusion',
            'source_summary': c.source_summary(), 'canonical_compilation': canonical.policy()}
