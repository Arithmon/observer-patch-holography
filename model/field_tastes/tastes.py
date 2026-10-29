"""The contact interaction couples every retained Weyl momentum pocket.

These are supplied momentum projectors in the same periodic two-fermion
action. They are not source-selected particle species or compiled detectors.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import product
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from ..bound_pair.fiber import PairFiber


def pocket_labels(side=6):
    """Exact cos-sign cells; side=2 mod4 excludes ambiguous pi/2 boundaries."""
    if type(side) is not int or not 6 <= side <= 18 or side%4 != 2:
        raise ValueError('side 6..18 congruent to two modulo four required')
    indices=np.array(list(product(range(side),repeat=3)),dtype=int)
    bits=((4*indices>side)&(4*indices<3*side)).astype(int)
    return 4*bits[:,0]+2*bits[:,1]+bits[:,2]


def setup(side=6):
    labels=pocket_labels(side)
    f=PairFiber([0.,0.,0.],side=side)
    return f,labels


def momentum_operator(fiber,interaction=1.):
    """Complete four-spin fiber; U=0 is an explicit action-deletion control."""
    if interaction not in (0.,1.):
        raise ValueError('retained interaction is zero or one')
    return sparse.block_diag(fiber.blocks,format='csr')+interaction*sparse.csr_matrix(
        np.outer(fiber.d,fiber.d.conj()))


def mask_for(labels,pockets):
    if (not isinstance(pockets,(tuple,list)) or not pockets
            or any(type(p) is not int or p not in range(8) for p in pockets)
            or len(set(pockets))!=len(pockets)):
        raise ValueError('distinct pocket indices 0..7 required')
    return np.repeat(np.isin(labels,pockets),4)


def initial(fiber,labels,pockets=(0,)):
    mask=mask_for(labels,pockets)
    weight=F(int(np.count_nonzero(mask)),fiber.dimension)
    vector=mask*fiber.d/math.sqrt(float(weight))
    if fiber.exchange_residual(vector)>2e-13:
        raise ArithmeticError('pocket failed the physical fermion exchange constraint')
    return vector,mask,weight


def commutator_certificate(side=6,pockets=(0,)):
    f,labels=setup(side)
    a,mask,w=initial(f,labels,pockets)
    qd=mask*f.d
    # Evaluate [H,Q]a directly, independently of its two-dimensional proof.
    acted=f.apply(mask*a)-mask*f.apply(a)
    coefficient=w*(1-w)
    residual=abs(float(np.vdot(acted,acted).real)-float(coefficient))
    free_comm=f.free_apply(mask*a)-mask*f.free_apply(a)
    # On span{Qd,(1-Q)d}, eigenvalues of [P,Q]*[P,Q] are w(1-w).
    if w<1:
        b=(~mask)*f.d/math.sqrt(float(1-w))
        basis=np.column_stack((a,b))
        comm=np.column_stack([f.apply(mask*v)-mask*f.apply(v) for v in basis.T])
        singular=np.linalg.svd(basis.conj().T@comm,compute_uv=False)
        restricted_residual=float(np.max(abs(singular**2-float(coefficient))))
    else:
        restricted_residual=0.
    return {'side':side,'selected_pockets':list(pockets),'contact_weight':str(w),
            'commutator_operator_norm_squared_exact':str(coefficient),
            'short_time_outside_probability_t_squared_coefficient_exact':str(coefficient),
            'free_commutator_norm':float(np.linalg.norm(free_comm)),
            'direct_action_squared_norm_residual':residual,
            'two_dimensional_singular_value_residual':restricted_residual,
            'physical_exchange_residual':f.exchange_residual(a),
            'interacting_pocket_reduces_full_action':w==1}


def complete_probabilities(vector,labels):
    v=np.asarray(vector,complex)
    if v.shape!=(4*len(labels),) or not np.isfinite(v).all():
        raise ValueError('one complete finite momentum-fiber vector required')
    density=np.sum(abs(v.reshape(-1,4))**2,axis=1)
    return [float(np.sum(density[labels==i])) for i in range(8)]


def leakage_bound(duration=F(1,32)):
    """Rational short-time bound using ||H|| <= 3/2 and integral remainder."""
    if not isinstance(duration,F) or duration<=0:
        raise ValueError('positive exact rational model duration required')
    # ||(1-Q)H a||=sqrt(7)/8 >=5/16, with upper bound1/3.
    # ||e^-itH-I+itH|| <= t²||H||²/2 for self-adjoint H.
    b=F(3,2)
    remainder=duration**2*b*b/2
    lower=max(F(0),duration*F(5,16)-remainder)**2
    upper=min(F(1),(duration*F(1,3)+remainder)**2)
    return {'model_time':str(duration),'Hamiltonian_norm_upper':str(b),
            'amplitude_remainder_upper':str(remainder),
            'unconditioned_outside_probability_lower':str(lower),
            'unconditioned_outside_probability_upper':str(upper),
            'rational_bound_covers_ideal_Hamiltonian_evolution':True,
            'floating_or_hardware_error_certified':False}


def evolution_report(side=6):
    f,labels=setup(side)
    start,mask,w=initial(f,labels)
    h=momentum_operator(f)
    free=momentum_operator(f,0.)
    rows=[]
    for duration in (F(1,32),F(1,2),F(2)):
        changed=expm_multiply(-1j*float(duration)*h,start)
        deleted=expm_multiply(-1j*float(duration)*free,start)
        probabilities=complete_probabilities(changed,labels)
        free_probabilities=complete_probabilities(deleted,labels)
        rows.append({'model_time':str(duration),'all_eight_pocket_probabilities':probabilities,
                     'free_control_all_eight_probabilities':free_probabilities,
                     'normalization_residual':abs(sum(probabilities)-1),
                     'physical_exchange_residual':f.exchange_residual(changed),
                     'unconditioned_outside_probability':sum(probabilities[1:]),
                     'omitting_other_pockets_loses_probability':1-probabilities[0],
                     'interaction_deletion_total_variation':sum(abs(a-b) for a,b in zip(
                         probabilities,free_probabilities,strict=True))/2})
    nodes=[]
    for i,bits in enumerate(product((0,1),repeat=3)):
        nodes.append({'pocket':i,'node_pi_coordinates':list(bits),
                      'Weyl_jacobian_chirality':-(-1)**sum(bits),
                      'relative_momentum_points':int(np.count_nonzero(labels==i)),
                      'contact_weight_exact':'1/8'})
    return {'total_momentum_pi_coordinates':[0,0,0],'torus_side':side,
            'full_four_spin_fiber_dimension':f.dimension,
            'initial_pocket':0,'initial_state_supplied':True,
            'field_momentum_detector_is_supplied_not_source_compiled':True,
            'all_nodes':nodes,'net_Weyl_node_chirality':sum(n['Weyl_jacobian_chirality'] for n in nodes),
            'commutator':commutator_certificate(side),
            'two_pocket_union_commutator':commutator_certificate(side,(0,1)),
            'all_pockets_control':commutator_certificate(side,tuple(range(8))),
            'short_time_bound':leakage_bound(),'complete_evolutions':rows}
