"""Native asymmetric pulses with finite spin duration and finite environment.

The oscillator/spin primitives are imported from the frozen symmetric packet.
Every dephasing application now contracts a positive, source-executed random
phase instrument. Neither its classical weights nor the apparatus are derived.
"""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
import math

import numpy as np

from ..anharmonic_observations import source as primitive
from ..observer.source_scalar import compiler as c
from ..scalar_observations.experimental_readout import canonical_source as canonical
from ..scalar_observations.experimental_readout import executable_source as native
from . import bounds, environment


def position_setting(h=.4):
    """Literal beta=i*h; division by sqrt(2)*h happens after the detector."""
    h = native._real(h, 'finite-difference coordinate')
    if h <= 0:
        raise ValueError('positive finite-difference coordinate required')
    return {'u': 0., 'v': h, 'theta': -math.pi/2}


def _sum_resources(items):
    return {key: sum(item[key] for item in items) for key in native.COUNT_KEYS}


@lru_cache(maxsize=64)
def spin_gap(dimension, angle, delta, omega_carrier, finite_spin_duration=True):
    """Execute exp[-i*(delta*N+sign(angle)*Omega0*X/2)*|angle|/Omega0].

    The commuting number and spin words are contracted separately. Setting
    finite_spin_duration=False retains the same spin rotation at zero gap.
    """
    dimension = native._dimension(dimension)
    angle, delta, omega_carrier = (native._real(value, name) for value, name in
                                  ((angle, 'spin angle'), (delta, 'delta'),
                                   (omega_carrier, 'carrier frequency')))
    if delta <= 0 or omega_carrier <= 0 or type(finite_spin_duration) is not bool:
        raise ValueError('positive frequencies and boolean duration flag required')
    duration = abs(angle)/omega_carrier if finite_spin_duration else 0.
    spin, spin_meta = primitive.spin_rotation('X', angle)
    context = canonical._context()
    phase = (canonical._parameter(context, delta, 'delta')
             * canonical._parameter(context, duration, 'gap duration'))
    # N=sum_b 2^b*(I-Z_b)/2; the scalar pulse preserves the full phase.
    # This exact identity avoids a generic Walsh expansion of a linear phase.
    word = [c.pulse('I', float(-phase*(dimension-1)/2), [0])]
    word.extend(c.pulse('Z', float(phase*2**bit/2), [bit])
                for bit in range(dimension.bit_length()-1))
    target = np.diag(np.exp(-1j*delta*duration*np.arange(dimension)))
    free, free_meta = native._contract(word, dimension.bit_length()-1, target,
                                      name='spin_gap_free_number_phase')
    unitary = np.kron(spin, free)
    unitary.flags.writeable = False
    return unitary, {
        'name': 'finite_spin_rotation_and_commuting_number_drift',
        'angle': angle, 'duration_seconds': duration,
        'omega_carrier_per_second': omega_carrier,
        'finite_spin_duration': finite_spin_duration,
        'resources': _sum_resources([spin_meta['resources'], free_meta['resources']]),
        'stages': [deepcopy(spin_meta), deepcopy(free_meta)],
        'matrix_origin': 'contracted_native_source_pulses',
        'noise_commutes_with_unitary': True,
    }


def run(dimension, alpha_real, thermal, step_counts, settings=None, *,
        delta=2*math.pi*500, omega=math.pi/150e-6, pulse_duration=50e-6,
        theta=.8, phi=0., dephasing_rate=18., subdivisions=128,
        omega_carrier=math.pi/35e-6, finite_spin_duration=True,
        quadrature_order=environment.DEFAULT_QUADRATURE_ORDER):
    """Execute the declared apparatus, without receiving observed outcomes.

    Chronology: H-(tau), Rx(-phi), Rz(theta), Rx(phi), H+(2*tau),
    Rx(phi), Rz(theta), Rx(-phi), H-(tau). The potential phase in the
    small-angle zero-drift limit is Phi=-phi. Published coordinate duration
    is 4*k*tau; the continuously advancing drive clock also includes Rx gaps.
    """
    dimension = native._dimension(dimension)
    if not isinstance(step_counts, (list, tuple)) or not step_counts:
        raise ValueError('nonempty list of step counts required')
    steps = [primitive._integer(value, 'step count') for value in step_counts]
    if steps != sorted(set(steps)):
        raise ValueError('step counts must be sorted and unique')
    subdivisions = primitive._integer(subdivisions, 'subdivisions', 1)
    if settings is None:
        settings = [position_setting()]
    if not isinstance(settings, (list, tuple)) or not settings:
        raise ValueError('nonempty list of detector settings required')
    settings = [native._setting(setting) for setting in settings]
    delta, omega, pulse_duration, theta, phi, dephasing_rate, omega_carrier = (
        native._real(value, name) for value, name in
        ((delta, 'delta'), (omega, 'omega'), (pulse_duration, 'pulse duration'),
         (theta, 'theta'), (phi, 'phi'), (dephasing_rate, 'dephasing rate'),
         (omega_carrier, 'carrier frequency')))
    if (delta <= 0 or omega < 0 or pulse_duration <= 0 or dephasing_rate < 0
            or omega_carrier <= 0 or type(finite_spin_duration) is not bool):
        raise ValueError('valid frequencies, durations, noise rate and boolean gap flag required')
    duration = pulse_duration/subdivisions
    states, weights, prepmeta = native.preparation(dimension, alpha_real, thermal)
    oscillator = (states*weights) @ states.conj().T
    rho = np.zeros((2*dimension, 2*dimension), complex)
    rho[:dimension, :dimension] = oscillator
    positive, pmeta = primitive.pulse_unitary(dimension, delta, omega, duration, 1)
    negative, nmeta = primitive.pulse_unitary(dimension, delta, omega, duration, -1)
    rotation_z, zmeta = primitive.spin_rotation('Z', theta)
    zjoint = np.kron(rotation_z, np.eye(dimension))
    plus_gap, gpmeta = spin_gap(dimension, phi, delta, omega_carrier, finite_spin_duration)
    minus_gap, gmmeta = spin_gap(dimension, -phi, delta, omega_carrier, finite_spin_duration)
    gap_duration = gpmeta['duration_seconds']
    drive_step_duration = 4*(pulse_duration+gap_duration)
    diffusion = environment.channel(dimension, dephasing_rate, duration/2, quadrature_order)
    gap_diffusion = environment.channel(dimension, dephasing_rate, gap_duration, quadrature_order)
    environment_per_step = 8*subdivisions+4
    branch_resources = _sum_resources(
        [diffusion.metadata['resources_per_realization']]*(8*subdivisions)
        + [gap_diffusion.metadata['resources_per_realization']]*4)
    environment_all_branches = _sum_resources(
        [diffusion.metadata['resources_for_all_branch_contractions']]*(8*subdivisions)
        + [gap_diffusion.metadata['resources_for_all_branch_contractions']]*4)
    unitary_resources = _sum_resources(
        [pmeta['resources'], nmeta['resources']]*(2*subdivisions)
        + [zmeta['resources']]*2 + [gpmeta['resources'], gmmeta['resources']]*2)
    total_resources = _sum_resources([unitary_resources, branch_resources])
    operations = [primitive.detector(dimension, **setting) for setting in settings]
    a = native._ladder(dimension)
    q = (a+a.T)/math.sqrt(2)
    records = []
    step_set = set(steps)
    error_half = Fraction(diffusion.metadata['ideal_quadrature_error']['exact_half_diamond_error_upper'])
    error_gap = Fraction(gap_diffusion.metadata['ideal_quadrature_error']['exact_half_diamond_error_upper'])
    for step in range(steps[-1]+1):
        if step in step_set:
            blocks = [rho[index*dimension:(index+1)*dimension,
                          index*dimension:(index+1)*dimension] for index in (0, 1)]
            spin_probabilities = [float(np.trace(block).real) for block in blocks]
            acceptance = spin_probabilities[0]
            if not 0 < acceptance <= 1+3e-9:
                raise ArithmeticError('invalid selected spin probability')
            reads = []
            for setting, (kraus, metadata) in zip(settings, operations, strict=True):
                joint = [[float(np.trace(operator @ block @ operator.conj().T).real)
                          for operator in kraus] for block in blocks]
                conditional = np.asarray(joint[0])/acceptance
                reads.append({**setting, 'spin_detector_joint_probabilities': joint,
                              'selected_detector_probabilities': conditional.tolist(),
                              'selected_mean': float(conditional[0]-conditional[1]),
                              'resources': deepcopy(metadata['resources'])})
            budget = bounds.execution_error_budget(
                dimension, step, subdivisions, min(1., acceptance),
                environment_applications=[(8*step*subdivisions, error_half), (4*step, error_gap)],
                delta=delta, omega=omega, pulse_duration=pulse_duration,
                dephasing_rate=dephasing_rate)
            records.append({
                'steps': step, 'simulated_duration_seconds': 4*step*pulse_duration,
                'drive_elapsed_duration_seconds': step*drive_step_duration,
                'spin_probabilities': spin_probabilities,
                'total_probability_mass': float(np.trace(rho).real),
                'selected_position_mean': float(np.trace(q @ blocks[0]).real/acceptance),
                'reads': reads, 'execution_error_budget': budget,
                'independent_environment_draws': step*environment_per_step,
                'evolution_resources_per_realization': {
                    key: step*value for key, value in total_resources.items()},
            })
        if step == steps[-1]:
            break
        for unitary, count, spin_gaps in (
                (negative, subdivisions, (minus_gap, plus_gap)),
                (positive, 2*subdivisions, (plus_gap, minus_gap)),
                (negative, subdivisions, None)):
            for _ in range(count):
                rho = diffusion.apply(rho)
                rho = unitary @ rho @ unitary.conj().T
                rho = diffusion.apply(rho)
            if spin_gaps is not None:
                first, second = spin_gaps
                rho = first @ rho @ first.conj().T
                rho = gap_diffusion.apply(rho)
                rho = zjoint @ rho @ zjoint.conj().T
                rho = second @ rho @ second.conj().T
                rho = gap_diffusion.apply(rho)
    return {
        'schema': 'source_executed_asymmetric_transfer_v1', 'dimension': dimension,
        'preparation': prepmeta, 'records': records,
        'operations': {'positive_subpulse': deepcopy(pmeta), 'negative_subpulse': deepcopy(nmeta),
                       'spin_z': deepcopy(zmeta), 'positive_spin_gap': deepcopy(gpmeta),
                       'negative_spin_gap': deepcopy(gmmeta),
                       'sideband_half_environment': deepcopy(diffusion.metadata),
                       'spin_gap_environment': deepcopy(gap_diffusion.metadata),
                       'detectors': [deepcopy(metadata) for _, metadata in operations]},
        'resources_per_evolution_step': total_resources,
        'unitary_resources_per_evolution_step': unitary_resources,
        'environment_resources_per_realization_per_step': branch_resources,
        'environment_resources_for_all_branch_contractions_per_step': environment_all_branches,
        'independent_environment_draws_per_step': environment_per_step,
        'settings': {'delta': delta, 'omega': omega, 'pulse_duration_seconds': pulse_duration,
                     'theta': theta, 'phi': phi, 'small_angle_potential_phase': -phi,
                     'dephasing_rate_per_second': dephasing_rate,
                     'subdivisions_per_pulse': subdivisions, 'omega_carrier': omega_carrier,
                     'finite_spin_duration': finite_spin_duration,
                     'spin_gap_duration_seconds': gap_duration,
                     'drive_elapsed_duration_per_step_seconds': drive_step_duration,
                     'published_coordinate_duration_per_step_seconds': 4*pulse_duration,
                     'environment_quadrature_order': quadrature_order},
        'time_parameter': 'published coordinate duration and continuously advancing drive elapsed duration retained separately',
        'drive_phase_continues_during_spin_gaps': finite_spin_duration,
        'gap_phase_compensation_authenticated_from_experiment': False,
        'final_spin_postselection': True, 'both_final_spin_branches_retained': True,
        'intermediate_spin_projection_or_reset': False,
        'all_environment_phase_unitaries_source_executed': True,
        'environment_distribution_and_rate_supplied': True,
        'analytic_diffusion_channel_inserted': False,
        'cutoff_error_certified_here': False, 'physical_service_error_certified': False,
        'prediction_origin': 'source_executed_unitaries_and_finite_source_phase_environment',
        'source_summary': c.source_summary(), 'canonical_compilation': canonical.policy(),
    }
