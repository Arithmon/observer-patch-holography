"""A declared rate-transfer service on the same complete interacting field.

Compatibility is fixed before execution: H=dGamma(h)+U D and the mathematical
reference is exp(-beta H)/Z with beta*U=log(25/8). The native count decoder
supplies an estimate of this single rate. It does not select H, prepare its
Gibbs state, or identify the field with the native carrier. Existing local
source controls execute H on every N=2 state, not a fitted band Hamiltonian.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product
import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from ..bound_pair.protection import PeriodicPair, canonical_band, packet_coefficients
from ..bound_pair.fiber import PairFiber
from ..bound_pair.spectrum import complete_spectrum
from ..interacting_field import source, bounds
from ..interacting_field.lattice import SIGMA
from ..interacting_field.experiment import pure_cq_distance
from ..particle_clock.clock import infer_clock
from ..local_repair_family import inference


class ClockField(PeriodicPair):
    """Same periodic action, with source-executable disjoint bond groups.

Even periodic side makes each spatial matching a disjoint mode matching.
Wrap bonds are ordered by mode index and their coefficient is conjugated;
this orientation matters for the Weyl symbol and fermionic routing signs.
"""
    def __init__(self, side=4, speed=.125, interaction=1.):
        if type(side) is not int or side not in (4, 6):
            raise ValueError('an even periodic side4 or6 required')
        super().__init__(side, speed, interaction)
        groups = [[] for _ in range(6)]
        for i, x in enumerate(self.sites):
            for axis in range(3):
                y = list(x); y[axis] = (y[axis]+1) % side
                j = self.site_index[tuple(y)]
                for a, b in product(range(2), repeat=2):
                    z = complex(1j*self.speed*SIGMA[axis][a,b]/2)
                    if not z:
                        continue
                    p, q = 2*i+a, 2*j+b
                    if p > q:
                        p, q, z = q, p, z.conjugate()
                    groups[2*axis+x[axis] % 2].append((p,q,z))
        self.groups = tuple(tuple(g) for g in groups)
        rows, cols, values = [], [], []
        for group in self.groups:
            modes = [m for p,q,_ in group for m in (p,q)]
            if len(set(modes)) != len(modes):
                raise ArithmeticError('periodic matching has overlapping modes')
            for p,q,z in group:
                rows.extend((p,q)); cols.extend((q,p)); values.extend((z,z.conjugate()))
        reconstructed = sparse.coo_matrix((values,(rows,cols)),shape=self.single_h.shape).tocsr()
        if np.max(abs((reconstructed-self.single_h).data), initial=0) > 1e-14:
            raise ArithmeticError('source bond groups changed the independent periodic action')


def rate_decision(counts, parameter=F(-1)):
    """Compile time from complete native endpoint counts alone.

The fixed rate window bounds all admitted source-word durations. Stops keep
the training tape and input state unchanged in both comparison protocols.
There is no retry or conditioning on an accepted estimate.
"""
    if type(parameter) not in (int, F) or abs(parameter) > 1:
        raise ValueError('exact rational modular parameter with |s|<=1 required')
    result = infer_clock(counts)
    if result['status'] != 'accept':
        return result
    ratio = F(result['energy_ratio_estimate'])
    if not F(3) <= ratio <= F(10,3):
        return {'status':'stop', 'reason':'outside_predeclared_rate_window'}
    # log(10/3)<5/4: the first four exponential terms already exceed10/3.
    beta = math.log(ratio.numerator)-math.log(ratio.denominator)
    return {'status':'accept', 'energy_ratio_estimate':str(ratio),
            'beta_estimate':beta, 'modular_parameter':str(parameter),
            'Hamiltonian_model_duration':-float(parameter)*beta,
            'absolute_model_duration_upper':str(F(5,4)*abs(parameter)),
            'beta_error_on_confidence_event':result['beta_times_U_error_on_confidence_event'],
            'statistical_failure_upper':result['statistical_failure_upper'],
            'rate_conversion_U': '1',
            'physical_seconds_inferred':False}


def error_budget(parameter=F(-1), steps=32):
    """Exact ideal-word bound on complete outputs, uniform over H_N2 tensor R.

Each hopping matching has N=2 norm<=c and D has norm1. Thus S=6c+U.
The symmetric product error is S^3 |t|^3/(3 n^2). Centering H by U/2
gives ||H-U/2||<=U/2+2sqrt(3)c<15/16. A two-edge native count path
supplies |beta_hat-beta|<=3/80 on the confidence event. Contractivity
covers all subsequent common instruments and arbitrary untouched R.
"""
    if type(parameter) not in (int, F) or abs(parameter)>1:
        raise ValueError('exact rational modular parameter with |s|<=1 required')
    if type(steps) is not int or steps<1:
        raise ValueError('positive split count required')
    s=abs(F(parameter)); t=F(5,4)*s
    split=F(bounds.history(((t,steps),))['complete_history_half_trace_upper'])
    clock=F(15,16)*F(3,80)*s
    failure=inference.fixed_budget()['statistical_failure_upper']
    return {'split_steps':steps, 'absolute_model_duration_upper':str(t),
            'split_generator_norm_sum_upper':'7/4',
            'centered_full_H_operator_norm_upper':'15/16',
            'ideal_split_halfdiamond_upper':str(split),
            'rate_error_halfdiamond_upper_on_confidence_event':str(clock),
            'statistical_failure_upper':str(failure),
            'complete_average_halfdiamond_upper':str(min(F(1),split+clock+failure)),
            'all_N2_states_and_arbitrary_untouched_reference':True,
            'complete_instruments_and_stop_outputs_retained':True,
            'conditioned_on_acceptance_or_pointer_outcome':False,
            'floating_or_physical_source_error_included':False}


def reference_spectrum(side=4, beta=math.log(25/8)):
    """One faithful Gibbs reference on the full N=2 algebra, including coherences.

Momentum fibers are only a computational decomposition. A SINGLE partition
function supplies the cross-fiber weights; no fiber is omitted or normalized
into an independent state. No measured input is replaced by this reference.
"""
    if side != 4 or type(side) is not int:
        raise ValueError('complete retained reference uses side4')
    if type(beta) not in (int,float) or not math.isfinite(beta) or beta<=0:
        raise ValueError('positive finite dimensionless beta required')
    energies, offsets, partitions = [], [0], []
    for index in product(range(side), repeat=3):
        f=PairFiber(list(2*np.pi*np.array(index)/side),side=side)
        spectrum=complete_spectrum(f)
        e=spectrum['energies']
        energies.extend(e.tolist()); offsets.append(len(energies))
        partitions.append(float(np.sum(np.exp(-beta*e))))
    energy=np.array(energies)
    if len(energy)!=math.comb(2*side**3,2):
        raise ArithmeticError('some physical field states were omitted')
    log_z=math.log(float(np.sum(np.exp(-beta*energy))))
    log_p=-beta*energy-log_z
    probability=np.exp(log_p)
    if not np.all(probability>0) or abs(float(sum(probability))-1)>2e-13:
        raise ArithmeticError('full field reference lost faithfulness or normalization')
    return {'energy':energy,'log_probability':log_p,'probability':probability,
            'fiber_offsets':offsets,'fiber_partitions':np.array(partitions),
            'log_partition':log_z,'beta':beta}


def reference_report():
    ref=reference_spectrum()
    energy=ref['energy']; logs=ref['log_probability']; beta=ref['beta']; s=-1.
    # One scalar removes normalization from the full unitary, including every
    # off-diagonal momentum block. Equality only within fibers is insufficient.
    error=float(np.max(abs(np.exp(1j*s*(logs+ref['log_partition']))-
                           np.exp(-1j*s*beta*energy))))
    z=ref['fiber_partitions']; offset=ref['fiber_offsets']
    lo,hi=int(np.argmin(z)),int(np.argmax(z))
    i,j=offset[lo],offset[hi]
    correct=np.exp(1j*s*(logs[i]-logs[j]))
    independently_normalized=np.exp(1j*s*(-beta*(energy[i]-energy[j])-math.log(z[lo]/z[hi])))
    return {'dimension':len(energy),'momentum_fibers':64,
            'fiber_dimensions':[b-a for a,b in zip(offset,offset[1:])],
            'minimum_probability':float(min(ref['probability'])),
            'maximum_probability':float(max(ref['probability'])),
            'full_trace':float(sum(ref['probability'])),
            'global_partition':math.exp(ref['log_partition']),
            'full_algebra_modular_Hamiltonian_unitary_residual':error,
            'independently_normalized_fibers_cross_coherence_error':float(abs(correct-independently_normalized)),
            'selected_cross_fiber_indices':[lo,hi],
            'Gibbs_reference_prepared_by_native_repair':False,
            'Gibbs_reference_dynamically_selected':False,
            'reference_input_is_the_observed_state':False,
            'common_faithful_reference_is_a_supplied_compatibility_premise':True}


def packet_with_reference():
    """Two translated coherent band packets entangled with an untouched qubit.

Preparation is supplied; the action, reference test, and detector never fit
or discard momentum amplitudes. All64 physical band fibers occur.
"""
    _,band,energy,_,_=canonical_band()
    a=packet_coefficients()
    kx=2*np.pi*np.array(list(product(range(4),repeat=3)))[:,0]/4
    coefficients=np.column_stack((a,a*np.exp(-1j*kx)))/math.sqrt(2)
    return band@coefficients,band,energy,coefficients


def pointers(space,state):
    """Complete two-bit occupation instrument on the origin's two modes."""
    return tuple(second for first in source.read_number(space,state,0)
                 for second in source.read_number(space,first,1))


def execute(counts,state,parameter=F(-1),steps=32,space=None):
    """Actual local source words; every state column is an untouched R basis."""
    space=ClockField() if space is None else space
    if not isinstance(space,ClockField) or space.speed != .125 or space.interaction != 1.:
        raise ValueError('canonical supplied clock field required')
    state=space.state(state)
    decision=rate_decision(counts,parameter)
    if decision['status']!='accept':
        return decision,{('stop',decision['reason']):state.copy()}
    error_budget(parameter,steps)
    out=source.evolve(space,state,decision['Hamiltonian_model_duration'],steps)
    return decision,{('read',a,b):v for (a,b),v in zip(product((0,1),repeat=2),pointers(space,out),strict=True)}


def exact_state(space,state,duration):
    return expm_multiply(-1j*duration*space.hamiltonian,state)


def cq_distance(first,second):
    return pure_cq_distance([v.ravel() for v in first],[v.ravel() for v in second])
