"""Exact Gaussian-integer nonflat-link Ward and pair-transfer certificate."""
from itertools import product
import json

import numpy as np
from scipy import sparse

from ..field_stress.energy import GaussianInteger as GI, commutator
from ..interacting_field.lattice import SIGMA
from .links import GaugedField


ROOTS = (1,1j,-1,-1j)


def exact_certificate():
    """The rational field is assembled directly, without rounding trig output.

    16h has Gaussian-integer coefficients; 256 PTQ QTP does too. A complete
    radius-one box retains 1,431 states, with nonflat fourth-root link phases.
    """
    field = GaugedField()
    modes = field.modes
    h = sparse.lil_matrix((modes,modes),dtype=complex)
    phase = {edge:ROOTS[(3*x+y+x*y)%4] for edge in field.axes for x,y in [edge]}
    # An explicit single-link twist guarantees a nonflat plaquette control.
    first_edge = next(iter(field.axes))
    phase[first_edge] *= 1j
    for (x,y),axis in field.axes.items():
        block = phase[x,y]*1j*SIGMA[axis]
        h[2*x:2*x+2,2*y:2*y+2] = block
        h[2*y:2*y+2,2*x:2*x+2] = block.conj().T
    h = h.tocsr()
    t = field.lift_one_body(h)
    contact = np.zeros(field.dimension,dtype=int);contact[field.contact_indices] = 16
    total = GI.from_exact(t+sparse.diags(contact,dtype=np.int64))
    currents = {}
    divergence = [GI.zero(field.dimension) for _ in field.sites]
    for x,y in field.axes:
        one = sparse.lil_matrix(h.shape,dtype=complex)
        one[2*x:2*x+2,2*y:2*y+2] = h[2*x:2*x+2,2*y:2*y+2]
        forward = GI.from_exact(field.lift_one_body(one.tocsr()))
        current = (forward-forward.adjoint()).times_i()
        currents[x,y] = current
        divergence[x] = divergence[x]+current
        divergence[y] = divergence[y]-current
    for nx,div in zip(field.local_numbers(),divergence,strict=True):
        if not (commutator(total,GI.from_exact(nx)).times_i()+div).is_zero():
            raise ArithmeticError('nonflat exact local Ward identity failed')

    alpha = np.array([ROOTS[(x*x+2*x+1)%4] for x in range(field.volume)])
    source_phase = alpha[field.first//2]*alpha[field.second//2]
    transformed_h = sparse.diags(np.repeat(alpha,2))@h@sparse.diags(np.repeat(alpha.conj(),2))
    transformed = GI.from_exact(field.lift_one_body(transformed_h)+sparse.diags(contact,dtype=np.int64))
    g = GI.from_exact(sparse.diags(source_phase))
    if not (transformed-g@total@g.adjoint()).is_zero():
        raise ArithmeticError('exact full-sector gauge covariance failed')

    p = field.contact_indices
    q = np.setdiff1d(np.arange(field.dimension),p)
    b = t[p][:,q]
    product_block = GI.from_exact(b@b.conj().T)
    independent = np.zeros((field.volume,field.volume),complex)
    for x,y in field.axes:
        # Every Weyl bond has two unit-modulus integer entries and determinant
        # equal to its link phase squared. Two virtual CAR paths contribute.
        independent[x,x] += 2; independent[y,y] += 2
        independent[x,y] = 2*phase[x,y]**2
        independent[y,x] = independent[x,y].conjugate()
    if not (product_block-GI.from_exact(independent)).is_zero():
        raise ArithmeticError('exact determinant pair-transfer identity failed')
    pair_g = GI.from_exact(sparse.diags(alpha**2))
    transformed_t = field.lift_one_body(transformed_h)
    bt = transformed_t[p][:,q]
    if not (GI.from_exact(bt@bt.conj().T)-pair_g@product_block@pair_g.adjoint()).is_zero():
        raise ArithmeticError('pair link does not carry twice the site phase')

    def oriented(x,y):
        return phase[x,y] if x<y else phase[y,x].conjugate()
    holonomies = []
    for start in field.sites:
        for a,baxis in ((0,1),(0,2),(1,2)):
            xa=list(start);xa[a]+=1
            xb=list(start);xb[baxis]+=1
            xab=xa.copy();xab[baxis]+=1
            if all(tuple(z) in field.site_index for z in (xa,xb,xab)):
                cycle=[field.site_index[tuple(z)] for z in (start,xa,xab,xb)]
                holonomy=np.prod([oriented(cycle[j],cycle[(j+1)%4]) for j in range(4)])
                holonomies.append(holonomy)
    if not any(z != 1 for z in holonomies):
        raise ArithmeticError('background accidentally has no nonflat control')
    return {
        'complete_N2_dimension':field.dimension,'sites':field.volume,
        'nonzero_bond_currents':len(currents),'exact_local_Ward_identities':field.volume,
        'full_sector_gauge_covariance_exact':True,
        'second_order_pair_transfer_2det_h_exact':True,
        'pair_site_character':2,'pair_effective_link_character':2,
        'square_holonomies_checked':len(holonomies),
        'nontrivial_square_holonomies':sum(int(z != 1) for z in holonomies),
        'source_selected_gauge_background':False,
        'physical_electromagnetic_charge_assigned':False,
        'full_band_replaced_by_leading_pair_action':False,
        'time_parameter':'Hamiltonian model time; supplied hbar=1',
    }


if __name__ == '__main__':
    print(json.dumps(exact_certificate(),indent=2))
