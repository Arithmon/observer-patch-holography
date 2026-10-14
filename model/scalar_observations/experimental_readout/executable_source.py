"""Preparation and readout by contracted twelve-port observer source words.

The classical thermal draw, blank register and coherent proper-code services
are supplied apparatus. Every branch's Fock preparation, displacement,
quadrature rotation and two-outcome detector is nevertheless executed from
native source pulses. No fitted density matrix or analytic characteristic
function is inserted in the prediction path.

The d=64 basis has 1,162,054 native pulses per direction. A complete cold
run can take several minutes; repeated settings share already contracted
operation matrices. The original pulse counts and full-phase guards remain
unchanged. Single-threaded BLAS avoids excessive overhead on these small
matrices (for example OPENBLAS_NUM_THREADS=1).
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
import json
import math
from numbers import Integral, Real

import numpy as np
from scipy.linalg import expm

from ...observer.source_scalar import compiler as c
from ..source import _basis


MAXIMUM_DIMENSION = 64
COUNT_KEYS = ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
              'identity_phase_pulses', 'proper_code_load_events',
              'proper_code_export_events', 'proper_code_service_events')


def _dimension(value):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or value < 2 or value > MAXIMUM_DIMENSION or value & (value - 1)):
        raise ValueError('a power-of-two occupation dimension from 2 through 64 is required')
    return int(value)


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f'{name} must be a finite real number')
    return float(value)


def _word_hash(word):
    payload = json.dumps(word, sort_keys=True, separators=(',', ':'), allow_nan=False)
    return hashlib.sha256(payload.encode()).hexdigest()


def _large_real_orthogonal(matrix):
    """The existing Gray-Givens algorithm extended from four to six qubits.

    Each Gray-neighbour rotation is its full commuting Walsh-Y expansion in
    the existing compiler's tensor-Pauli words. No new source primitive is
    introduced. At dimensions <=16 we call the established implementation.
    """
    value = np.asarray(matrix, dtype=float)
    dimension = _dimension(len(value))
    if value.shape != (dimension, dimension) or not np.isfinite(value).all():
        raise ValueError('a finite real square basis is required')
    if np.linalg.norm(value.T @ value - np.eye(dimension), 2) > 2e-11:
        raise ValueError('an orthogonal basis is required')
    if dimension <= 16:
        return c.compile_real_orthogonal(value)
    return _extended_gray_word(value)


def _extended_gray_word(value):
    """Same complete Gray elimination, with dimension validation by caller."""
    dimension = len(value)
    bits = dimension.bit_length() - 1
    gray = [j ^ (j >> 1) for j in range(dimension)]
    residual = value[np.ix_(gray, gray)].copy()
    rotations = []
    for column in range(dimension - 1):
        for row in range(dimension - 1, column, -1):
            x, y = residual[row - 1, column], residual[row, column]
            theta = 0. if x == 0 and y == 0 else math.atan2(y, x)
            cosine, sine = math.cos(theta), math.sin(theta)
            old = residual[[row - 1, row], :].copy()
            residual[row - 1], residual[row] = cosine * old[0] + sine * old[1], -sine * old[0] + cosine * old[1]
            rotations.append((gray[row - 1], gray[row], theta))
    signs = np.where(np.diag(residual) >= 0, 1., -1.)
    if np.linalg.norm(residual - np.diag(signs), 2) > 3e-11:
        raise ArithmeticError('Gray elimination did not retain the complete basis')
    phases = np.zeros(dimension)
    for j, sign in enumerate(signs):
        phases[gray[j]] = 0 if sign == 1 else math.pi
    word = c.compile_diagonal(phases)
    for a, b, theta in reversed(rotations):
        difference = a ^ b
        target = difference.bit_length() - 1
        orientation = (-1) ** ((a >> target) & 1)
        for mask in range(1 << bits):
            if mask & difference:
                continue
            label = ['Y' if j == target else 'Z' if mask >> j & 1 else 'I' for j in range(bits)]
            coefficient = -orientation * theta * (-1) ** ((a & mask).bit_count()) / (1 << (bits - 1))
            word.extend(c.pauli_rotation(label, coefficient))
    return word


def _contract(word, qubits, target, *, name):
    """Only the native contraction supplies the returned operation matrix."""
    raw = c.execute(word, qubits)
    error = float(np.linalg.norm(raw - target, 2))
    if not math.isfinite(error) or error > 2e-10:
        raise ArithmeticError(f'{name}: native contraction differs from its independent target')
    matrix, stabilization = c.stabilize_executed_unitary(raw)
    matrix.flags.writeable = False
    return matrix, {
        'name': name, 'word_sha256': _word_hash(word),
        'resources': c.resources(word), 'register_qubits': qubits,
        'full_phase_target_residual': error,
        'numerical_stabilization': stabilization,
        'matrix_origin': 'contracted_native_source_pulses',
        'physical_service_error_certified': False,
    }


def _ladder(dimension):
    return np.diag(np.sqrt(np.arange(1, dimension)), 1)


@lru_cache(maxsize=32)
def _displacement(dimension, alpha_real):
    ladder = _ladder(dimension)
    target = expm(alpha_real * (ladder.T - ladder))
    word = _large_real_orthogonal(target)
    return _contract(word, dimension.bit_length() - 1, target, name='real_displacement')


@lru_cache(maxsize=8)
def _coordinate_basis(dimension):
    ladder = _ladder(dimension)
    coordinate = (ladder + ladder.T) / 2
    values, basis = _basis(coordinate)
    word = _large_real_orthogonal(basis)
    bits = dimension.bit_length() - 1
    forward, metadata = _contract(word, bits, basis, name='coordinate_basis')
    backward, inverse_metadata = _contract(c.inverse(word), bits, basis.T, name='coordinate_inverse_basis')
    return values, forward, backward, metadata, inverse_metadata


@lru_cache(maxsize=8)
def _fock_preparations(dimension):
    bits = dimension.bit_length() - 1
    columns = []
    records = []
    blank = np.zeros(dimension, complex)
    blank[0] = 1
    for occupation in range(dimension):
        # exp(-i*pi/2*I) exp(i*pi/2*X)=X, with the scalar phase retained.
        word = []
        for bit in range(bits):
            if occupation >> bit & 1:
                word.extend([c.pulse('X', math.pi / 2, [bit]),
                             c.pulse('I', -math.pi / 2, [bit])])
        state = c.execute(word, bits, blank)
        target = np.zeros(dimension, complex)
        target[occupation] = 1
        residual = float(np.linalg.norm(state - target))
        if residual > 2e-12:
            raise ArithmeticError('native Fock preparation failed')
        columns.append(state)
        records.append({'occupation': occupation, 'word_sha256': _word_hash(word),
                        'resources': c.resources(word), 'full_phase_target_residual': residual})
    matrix = np.column_stack(columns)
    matrix.flags.writeable = False
    return matrix, records


def thermal_weights(dimension, thermal):
    """A supplied capped geometric draw, with no rejected preparation branch.

    Occupations >= dimension-1 are routed to the last code state. The trace
    distance to the unbounded thermal mixture before displacement is r**d.
    The finite classical weights are defined directly, not by conditioning
    on an accepted draw or renormalizing a quantum state.
    """
    dimension = _dimension(dimension)
    thermal = _real(thermal, 'thermal occupation')
    if thermal < 0:
        raise ValueError('thermal occupation must be nonnegative')
    ratio = thermal / (1 + thermal)
    weights = [(1 - ratio) * ratio ** n for n in range(dimension - 1)]
    weights.append(ratio ** (dimension - 1))
    return np.array(weights), {
        'rule': 'capped_geometric_classical_preparation',
        'mean_thermal_occupation_parameter': thermal,
        'ratio': ratio,
        'classical_branch_probabilities': weights,
        'initial_thermal_trace_distance': ratio ** dimension,
        'thermal_draw_supplied': True,
        'no_preparation_branch_discarded': True,
        'state_or_outcome_renormalization': False,
    }


def preparation(dimension, alpha_real, thermal):
    """Return actually prepared state columns and supplied classical weights."""
    dimension = _dimension(dimension)
    alpha_real = _real(alpha_real, 'real displacement')
    weights, thermal_metadata = thermal_weights(dimension, thermal)
    fock, fock_metadata = _fock_preparations(dimension)
    displacement, displacement_metadata = _displacement(dimension, alpha_real)
    states = displacement @ fock
    metadata = {
        'dimension': dimension, 'alpha_real': alpha_real,
        'thermal': thermal_metadata,
        'fock_branches': deepcopy(fock_metadata),
        'displacement': deepcopy(displacement_metadata),
        'weighted_prepared_mass': float(np.dot(weights, np.sum(abs(states) ** 2, axis=0))),
        'blank_register_and_classical_draw_supplied': True,
        'initialization_storage_readout_physical_costs_included': False,
        'postselection': False,
    }
    return states, weights, metadata


def _setting(setting):
    if not isinstance(setting, dict) or set(setting) != {'u', 'v', 'theta'}:
        raise ValueError('each setting must contain exactly u, v and theta')
    return {key: _real(setting[key], key) for key in ('u', 'v', 'theta')}


def source_run(dimension, alpha_real, thermal, settings):
    """Emit binary probabilities from source-executed preparation and detector.

    ``u+i*v`` denotes D(u+i*v), and the bit difference is
    Re(exp(i*theta)*chi(u+i*v)); theta=0 and -pi/2 yield real and imaginary
    components. The quadrature rotation angle is a setting, not model time.
    This function accepts no observed means or withheld response values.
    """
    dimension = _dimension(dimension)
    if not isinstance(settings, (tuple, list)) or not settings:
        raise ValueError('a nonempty list of detector settings is required')
    settings = [_setting(setting) for setting in settings]
    states, weights, preparation_metadata = preparation(dimension, alpha_real, thermal)
    values, basis, inverse_basis, basis_metadata, inverse_metadata = _coordinate_basis(dimension)
    bits = dimension.bit_length() - 1
    pointer = bits
    inverse_axis_word = [c.pulse('X', -math.pi / 4, [pointer])]
    axis_word = [c.pulse('X', math.pi / 4, [pointer])]
    inverse_axis, inverse_axis_metadata = _contract(
        inverse_axis_word, bits + 1,
        np.kron(expm(-1j * math.pi / 4 * np.array([[0, 1], [1, 0]])), np.eye(dimension)),
        name='pointer_inverse_axis')
    axis, axis_metadata = _contract(
        axis_word, bits + 1,
        np.kron(expm(1j * math.pi / 4 * np.array([[0, 1], [1, 0]])), np.eye(dimension)),
        name='pointer_axis')
    records = []
    for setting in settings:
        u, v, theta = (setting[key] for key in ('u', 'v', 'theta'))
        radius = math.hypot(u, v)
        phi = 0. if radius == 0 else math.atan2(-u, v)
        phases = -phi * np.arange(dimension)
        rotation_word = c.compile_diagonal(phases)
        rotation, rotation_metadata = _contract(rotation_word, bits, np.diag(np.exp(1j * phases)),
                                                name='quadrature_rotation')
        rotated = inverse_basis @ rotation @ states
        joint = np.vstack([rotated, np.zeros_like(rotated)])
        angles = radius * values + theta / 2
        pointer_phases = np.r_[-angles, angles]
        conditional_word = c.compile_diagonal(pointer_phases)
        conditional, conditional_metadata = _contract(
            conditional_word, bits + 1, np.diag(np.exp(1j * pointer_phases)),
            name='conditional_detector_phase')
        joint = axis @ conditional @ inverse_axis @ joint
        branch_zero, branch_one = basis @ joint[:dimension], basis @ joint[dimension:]
        branch_probabilities = np.array([np.sum(abs(branch_zero) ** 2, axis=0),
                                        np.sum(abs(branch_one) ** 2, axis=0)])
        probabilities = branch_probabilities @ weights
        stages = [preparation_metadata['displacement'], rotation_metadata,
                  inverse_metadata, inverse_axis_metadata, conditional_metadata,
                  axis_metadata, basis_metadata]
        common_counts = {key: sum(stage['resources'][key] for stage in stages)
                         for key in COUNT_KEYS}
        branch_counts = [{key: common_counts[key] + branch['resources'][key]
                          for key in COUNT_KEYS}
                         for branch in preparation_metadata['fock_branches']]
        resources = {key: {
            'maximum_per_shot': max(branch[key] for branch in branch_counts),
            'classical_weighted_mean_per_shot': float(sum(weight * branch[key]
                       for weight, branch in zip(weights, branch_counts, strict=True))),
        } for key in COUNT_KEYS}
        records.append({
            **setting, 'quadrature_angle': phi, 'detector_strength_eta': 2 * radius,
            'probabilities': probabilities.tolist(),
            'mean': float(probabilities[0] - probabilities[1]),
            'probability_mass': float(sum(probabilities)),
            'probabilities_by_initial_occupation': branch_probabilities.T.tolist(),
            'rotation': rotation_metadata, 'conditional_detector_phase': conditional_metadata,
            'chronological_word_hashes_after_fock_preparation': [stage['word_sha256'] for stage in stages],
            'resources': resources,
            'native_pulses_by_initial_occupation': [branch['native_pulses'] for branch in branch_counts],
            'both_detector_outcomes_retained': True,
            'postselection': False,
        })
    return {
        'schema': 'source_executed_oscillator_readout_v1',
        'dimension': dimension,
        'preparation': preparation_metadata,
        'shared_detector_operations': {
            'coordinate_basis': deepcopy(basis_metadata),
            'coordinate_inverse_basis': deepcopy(inverse_metadata),
            'pointer_inverse_axis': inverse_axis_metadata, 'pointer_axis': axis_metadata,
        },
        'records': records,
        'prediction_origin': 'norms_of_source_executed_pointer_branches',
        'analytic_characteristic_used_for_prediction': False,
        'source_summary': c.source_summary(),
        'supplied_apparatus': ['blank register', 'classical thermal draw',
                              'coherent proper-code load/export', 'pointer readout'],
        'service_failure_outcomes_retained': True,
        'valid_code_failure_probability_in_exact_arithmetic': 0,
        'physical_service_error_certified': False,
        'occupation_error_certified_here': False,
        'time_parameter': 'none; quadrature rotation is an apparatus setting',
    }
