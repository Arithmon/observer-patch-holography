"""The existing anharmonic apparatus as a complete primitive-observer program.

The native pulse compiler, preparation and finite environment are unchanged.
Each published duration names a separate destructively read experiment. A
multi-read program is a different, explicitly requested intervention protocol.
No empirical outcome or fitted continuation enters this module.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import math

import numpy as np

from tools import strict_json
from ..anharmonic_observations import source as primitive
from ..anharmonic_transfer import source, environment
from ..observer.source_scalar import compiler
from ..scalar_observations.experimental_readout import executable_source as native
from . import runtime


@dataclass(frozen=True)
class Apparatus:
    dimension: int = 32
    subdivisions: int = 128
    delta: float = 2*math.pi*500
    omega: float = math.pi/150e-6
    pulse_duration: float = 50e-6
    theta: float = .8
    phi: float = 0.
    dephasing_rate: float = 18.
    omega_carrier: float = math.pi/35e-6
    finite_spin_duration: bool = True
    quadrature_order: int = 16
    h: float = .4

    def __post_init__(self):
        object.__setattr__(self, 'dimension', native._dimension(self.dimension))
        object.__setattr__(self, 'subdivisions',
                           primitive._integer(self.subdivisions, 'subdivisions', 1))
        if self.subdivisions > 4096:
            raise ValueError('subdivision allocation budget exceeded')
        for name in ('delta', 'omega', 'pulse_duration', 'theta', 'phi',
                     'dephasing_rate', 'omega_carrier', 'h'):
            object.__setattr__(self, name, native._real(getattr(self, name), name))
        if (self.delta <= 0 or self.omega < 0 or self.pulse_duration <= 0
                or self.dephasing_rate < 0 or self.omega_carrier <= 0 or self.h <= 0
                or type(self.finite_spin_duration) is not bool):
            raise ValueError('valid frequencies, durations, readout and duration flag required')
        environment._order(self.quadrature_order)

    @property
    def gap_duration(self):
        return abs(self.phi)/self.omega_carrier if self.finite_spin_duration else 0.

    @property
    def published_step_duration(self):
        return 4*self.pulse_duration

    @property
    def drive_step_duration(self):
        return 4*(self.pulse_duration+self.gap_duration)

    def legacy_arguments(self):
        return {key: getattr(self, key) for key in self.__dataclass_fields__
                if key not in ('dimension', 'h')}


def _resources(metadata, **extra):
    counts = metadata.get('resources', {})
    return tuple(sorted(({key: counts[key] for key in native.COUNT_KEYS if key in counts}
                         | extra).items()))


def _unitary(matrix, metadata):
    return runtime.Instrument(kraus=((matrix,),), visible=False,
                              resources=_resources(metadata), binding=metadata)


def _diffusion(channel):
    # Work is spin x occupation; reference is the LAST tensor factor in the
    # runtime. The runtime lifts this work-basis Schur matrix over that factor.
    return runtime.Instrument(
        schur=np.tile(channel.schur, (2, 2)), visible=False,
        resources=tuple(sorted((channel.metadata['resources_per_realization'] |
                                {'independent_environment_draws': 1}).items())),
        binding=channel.metadata)


def _spin_record(dimension):
    """Contract a fresh native CNOT pointer; retain both projective outputs."""
    work = 2*dimension
    spin_bit, pointer_bit = dimension.bit_length()-1, dimension.bit_length()
    word = compiler.cnot(spin_bit, pointer_bit)
    target = np.zeros((2*work, 2*work), complex)
    for index in range(2*work):
        target[index ^ (((index >> spin_bit) & 1) << pointer_bit), index] = 1
    unitary, metadata = native._contract(word, pointer_bit+1, target,
                                        name='primitive_observer_spin_record')
    return tuple(unitary[s*work:(s+1)*work, :work] for s in (0, 1)), metadata


def read_instrument(apparatus):
    """Complete four-outcome instrument: label=2*spin+detector_bit."""
    apparatus = replace(apparatus)
    spin, spin_meta = _spin_record(apparatus.dimension)
    detector, detector_meta = primitive.detector(
        apparatus.dimension, **source.position_setting(apparatus.h))
    groups = tuple((np.kron(np.eye(2), detector[bit]) @ spin[choice],)
                   for choice in (0, 1) for bit in (0, 1))
    resources = source._sum_resources([spin_meta['resources'], detector_meta['resources']])
    return runtime.Instrument(
        kraus=groups, visible=True,
        resources=tuple(sorted((resources | {'fresh_pointer_blanks': 2,
                                              'committed_classical_bits': 2}).items())),
        binding={'schema': 'oph.primitive_observer.spin_pointer_read.v1',
                 'h': apparatus.h, 'spin_pointer': spin_meta, 'oscillator_pointer': detector_meta,
                 'outcome_meaning': 'label=2*spin+detector; spin=0 accepted',
                 'pointer_blank_and_commitment_services_supplied': True})


def compiled_cycle(apparatus):
    """Five fixed events, preserving the original pulse/noise chronology."""
    apparatus = replace(apparatus)
    a, n = apparatus, apparatus.subdivisions
    duration = a.pulse_duration/n
    positive, pm = primitive.pulse_unitary(a.dimension, a.delta, a.omega, duration, 1)
    negative, nm = primitive.pulse_unitary(a.dimension, a.delta, a.omega, duration, -1)
    z, zm = primitive.spin_rotation('Z', a.theta)
    plus_gap, gpm = source.spin_gap(a.dimension, a.phi, a.delta, a.omega_carrier,
                                   a.finite_spin_duration)
    minus_gap, gmm = source.spin_gap(a.dimension, -a.phi, a.delta, a.omega_carrier,
                                    a.finite_spin_duration)
    noise = _diffusion(environment.channel(a.dimension, a.dephasing_rate,
                                           duration/2, a.quadrature_order))
    gap_noise = _diffusion(environment.channel(a.dimension, a.dephasing_rate,
                                               a.gap_duration, a.quadrature_order))
    pos, neg = _unitary(positive, pm), _unitary(negative, nm)
    plus, minus = _unitary(plus_gap, gpm), _unitary(minus_gap, gmm)
    spin_z = _unitary(np.kron(z, np.eye(a.dimension)), zm)
    return (
        runtime.Event('negative_sideband_first', sequence=(noise, neg, noise), repeat=n),
        runtime.Event('first_spin_rotation', sequence=(minus, gap_noise, spin_z, plus, gap_noise)),
        runtime.Event('positive_sideband', sequence=(noise, pos, noise), repeat=2*n),
        runtime.Event('second_spin_rotation', sequence=(plus, gap_noise, spin_z, minus, gap_noise)),
        runtime.Event('negative_sideband_last', sequence=(noise, neg, noise), repeat=n))


def program(apparatus, alpha_real, thermal, steps, *, read_steps=None):
    """Build one fixed experiment; repeated destructive reads require opt-in."""
    if type(apparatus) is not Apparatus:
        raise ValueError('a typed apparatus required')
    apparatus = replace(apparatus)
    steps = primitive._integer(steps, 'steps')
    if steps > 100:
        raise ValueError('finite experiment allocation budget exceeded')
    if read_steps is None:
        read_steps = (steps,)
    if (type(read_steps) is not tuple
            or any(type(value) is not int or not 0 <= value <= steps for value in read_steps)
            or tuple(sorted(set(read_steps))) != read_steps):
        raise ValueError('distinct increasing read coordinates within the experiment required')
    if len(read_steps) > 4:
        raise ValueError('at most four complete four-outcome reads supported')
    states, weights, metadata = native.preparation(apparatus.dimension, alpha_real, thermal)
    # The capped thermal draw is hidden. Its possible Fock words have different
    # lengths, so a single integer per-path ledger uses an explicit upper bound.
    prep_counts = {
        key: metadata['displacement']['resources'][key]
             + max(branch['resources'][key] for branch in metadata['fock_branches'])
        for key in native.COUNT_KEYS}
    oscillator = (states*weights) @ states.conj().T
    density = np.zeros((2*apparatus.dimension, 2*apparatus.dimension), complex)
    density[:apparatus.dimension, :apparatus.dimension] = oscillator
    cycle = compiled_cycle(apparatus) if steps else ()
    detector = read_instrument(apparatus)
    events = []
    for step in range(steps+1):
        if step in read_steps:
            events.append(runtime.Event(f'spin_and_pointer_read_at_{step}', instruments=(detector,)))
        if step < steps:
            events.extend(cycle)
    binding = {
        'schema': 'oph.primitive_observer.anharmonic.v1',
        'apparatus': apparatus.__dict__, 'steps': steps, 'read_steps': read_steps,
        'preparation': metadata,
        'preparation_resource_semantics': 'componentwise upper bound per hidden thermal realization',
        'published_coordinate_seconds': steps*apparatus.published_step_duration,
        'drive_elapsed_seconds': steps*apparatus.drive_step_duration,
        'clock_status': 'two supplied apparatus coordinates; no recovered physical clock',
        'record_dictionary': '2*final_spin+pointer_bit; 0,1 accepted; 2,3 rejected',
        'destructive_read_changes_continuation': True,
        'empirical_outcomes_received': False,
        'physical_noise_preparation_and_detector_parameters_supplied': True,
        'matrix_origin': 'existing contracted native source primitives',
    }
    return runtime.Program(
        'anharmonic_native_observer', 2*apparatus.dimension, density, tuple(events),
        preparation_resources=tuple(sorted((prep_counts | {'work_preparations': 1}).items())),
        binding=binding)


def snapshots(apparatus, alpha_real, thermal, step_counts):
    """Yield separate complete read experiments, reusing unread host prefixes.

    The unread trunk is never updated with a measured state. Each returned
    history contains the whole prefix and its own destructive final read,
    and can be replayed independently from its original preparation.
    """
    if (type(step_counts) is not tuple or not step_counts
            or any(type(step) is not int or step < 0 for step in step_counts)
            or tuple(sorted(set(step_counts))) != step_counts):
        raise ValueError('distinct increasing integer snapshot coordinates required')
    trunk = program(apparatus, alpha_real, thermal, step_counts[-1], read_steps=())
    read = runtime.Event('spin_and_pointer_read', instruments=(read_instrument(apparatus),))
    history = runtime.execute(trunk, stop_after=0)
    for step in step_counts:
        history = runtime.execute(trunk, checkpoint=history, stop_after=5*step)
        measured = runtime.branch_experiment(
            history, (read,), name=f'anharmonic_independent_endpoint_{step}',
            binding={'requested_steps': step,
                     'published_coordinate_seconds': step*apparatus.published_step_duration,
                     'drive_elapsed_seconds': step*apparatus.drive_step_duration,
                     'earlier_snapshot_read_is_not_part_of_this_experiment': True})
        yield step, measured


def readout(history, *, record_index=-1, h=None):
    """Report conditioning without dropping any branch from the history."""
    history = runtime.validate_history(history)
    if type(record_index) is not int:
        raise ValueError('an integer record index required')
    reads = [event for event in history.program.events[:history.cursor] if event.visible]
    if not reads or not -len(reads) <= record_index < len(reads):
        raise ValueError('the requested committed read must exist')
    meanings = [strict_json.loads(instrument.binding)
                for instrument in reads[record_index].instruments]
    if (any(item.get('schema') != 'oph.primitive_observer.spin_pointer_read.v1'
            for item in meanings) or any(item['h'] != meanings[0]['h'] for item in meanings)):
        raise ValueError('an executed anharmonic spin/pointer read with one calibrated coordinate required')
    executed_h = native._real(meanings[0]['h'], 'executed readout coordinate')
    if h is not None and native._real(h, 'readout coordinate') != executed_h:
        raise ValueError('reported coordinate differs from the executed detector')
    h = executed_h
    if h <= 0:
        raise ValueError('positive executed coordinate required')
    probabilities = np.zeros(4)
    for branch in history.branches:
        if not branch.record or not -len(branch.record) <= record_index < len(branch.record):
            raise ValueError('the requested committed record must exist')
        label = branch.record[record_index]
        if type(label) is not int or label not in range(4):
            raise ValueError('a four-outcome spin/pointer record required')
        probabilities[label] += branch.probability
    acceptance = float(probabilities[0]+probabilities[1])
    if acceptance <= 0:
        raise ValueError('cannot condition on zero acceptance')
    mean = float((probabilities[0]-probabilities[1])/acceptance)
    return {'spin_detector_joint_probabilities': probabilities.reshape(2, 2).tolist(),
            'acceptance_probability': acceptance,
            'rejection_probability': float(probabilities[2]+probabilities[3]),
            'selected_mean': mean, 'position_estimator': mean/(math.sqrt(2)*h),
            'total_probability_mass': float(probabilities.sum())}
