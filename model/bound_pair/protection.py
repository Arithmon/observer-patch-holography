"""Closed-action spectral protection and complete periodic pair transport.

All spatial two-fermion amplitudes are retained.  The periodic box and the
bound-band packet preparation are supplied; no native encoder is claimed.
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import product
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from ..interacting_field.lattice import SIGMA, TwoFermions, real
from .fiber import PairFiber, parameters


class PeriodicPair(TwoFermions):
    """Independent position-space CAR lift, with every periodic bond once."""
    def __init__(self, side=4, speed=.125, interaction=1.):
        if type(side) is not int or not 4 <= side <= 6:
            raise ValueError('retained periodic side must be 4..6')
        self.side = side
        self.speed = real(speed, 'coordinate speed', nonnegative=True)
        self.interaction = real(interaction, 'repulsion', nonnegative=True)
        self.sites = tuple(product(range(side), repeat=3))
        self.site_index = {x:i for i,x in enumerate(self.sites)}
        self.volume = side**3
        super().__init__(2*self.volume)
        rows,cols,vals = [],[],[]
        for i,x in enumerate(self.sites):
            for axis in range(3):
                y=list(x); y[axis]=(y[axis]+1)%side
                j=self.site_index[tuple(y)]
                for a,b in product(range(2),repeat=2):
                    z=1j*self.speed*SIGMA[axis][a,b]/2
                    if z:
                        p,q=2*i+a,2*j+b
                        rows.extend((p,q)); cols.extend((q,p))
                        vals.extend((z,z.conjugate()))
        self.single_h=sparse.coo_matrix((vals,(rows,cols)),shape=(self.modes,self.modes)).tocsr()
        self.free_h=self.lift_one_body(self.single_h)
        self.contact_indices=np.array([self.index(2*i,2*i+1) for i in range(self.volume)])
        diagonal=np.zeros(self.dimension)
        diagonal[self.contact_indices]=self.interaction
        self.hamiltonian=self.free_h+sparse.diags(diagonal)
        xs=np.array(self.sites)[self.first//2]
        ys=np.array(self.sites)[self.second//2]
        delta=np.abs(xs-ys)
        self.relative_distance=np.sum(np.minimum(delta,side-delta),axis=1)

    def fiber_state(self, momentum_index, vector):
        """Unitary Fourier inclusion into the full antisymmetric N=2 space."""
        index=np.asarray(momentum_index)
        if index.shape!=(3,) or index.dtype.kind not in 'iu' or np.any(index<0) or np.any(index>=self.side):
            raise ValueError('one reciprocal-grid momentum index required')
        v=np.asarray(vector,complex)
        if v.shape!=(4*self.volume,) or not np.isfinite(v).all():
            raise ValueError('one complete finite fiber vector required')
        k=np.array(self.sites,dtype=int)
        q=(index-k)%self.side
        tensor=np.zeros((self.side,)*3+(2,)+(self.side,)*3+(2,),complex)
        tensor[k[:,0],k[:,1],k[:,2],:,q[:,0],q[:,1],q[:,2],:]=v.reshape(-1,2,2)
        matrix=np.fft.ifftn(tensor,axes=(0,1,2,4,5,6),norm='ortho').reshape(self.modes,self.modes)
        if np.linalg.norm(matrix+matrix.T)>2e-12:
            raise ValueError('fiber vector violates fermion exchange')
        return math.sqrt(2)*matrix[self.first,self.second]

    def observations(self,state):
        state=self.state(state)
        if state.ndim!=1:
            raise ValueError('one spatial pair state required')
        probability=np.abs(state)**2
        number=np.bincount(self.first//2,weights=probability,minlength=self.volume)
        number+=np.bincount(self.second//2,weights=probability,minlength=self.volume)
        contact=probability[self.contact_indices]
        return {'trace':float(np.sum(probability)),
                'site_particle_number':number.tolist(),
                'site_contact_probability':contact.tolist(),
                'separated_probability':float(np.sum(probability[self.relative_distance>0]))}


def canonical_band():
    """All 64 isolated eigenvectors of the retained L=4 periodic action."""
    space=PeriodicPair()
    vectors,energies,residues,exchanges=[],[],[],[]
    for km in space.sites:
        f=PairFiber(list(2*np.pi*np.array(km)/space.side),side=space.side)
        pole=f.pole(); v=f.eigenstate()
        vectors.append(space.fiber_state(km,v))
        energies.append(pole['energy']); residues.append(pole['contact_residue'])
        exchanges.append(f.exchange_residual(v))
    band=np.column_stack(vectors)
    return space,band,np.array(energies),np.array(residues),max(exchanges)


def packet_coefficients():
    """A predeclared broad packet directed along the x coordinate."""
    k=2*np.pi*np.array(list(product(range(4),repeat=3)))/4
    coefficient=np.exp(.5*(np.cos(k[:,0]-np.pi/2)+np.cos(k[:,1]-np.pi)+np.cos(k[:,2]-np.pi)))
    return coefficient/np.linalg.norm(coefficient)


def tail_bound(radius):
    """Exact canonical probability bound for relative distance >= radius.

    ||T||=sqrt(3)/4 < 7/16; this rational relaxation is outward valid.
    """
    if type(radius) is not int or radius<0:
        raise ValueError('nonnegative integer relative radius required')
    b=F(7,16)
    return min(F(1),(1/(1-b))**2*b**(2*radius))


def report():
    space,band,energy,residue,exchange=canonical_band()
    a=packet_coefficients()
    states=[band@(a*np.exp(-1j*energy*t)) for t in (0.,8.,16.)]
    records=[]
    for t,state in zip((0.,8.,16.),states,strict=True):
        obs=space.observations(state)
        obs['model_duration']=t
        obs['relative_tails']=[{
            'distance_at_least':r,
            'probability':float(np.sum(np.abs(state[space.relative_distance>=r])**2)),
            'rational_upper':str(tail_bound(r))} for r in range(1,5)]
        records.append(obs)
    exact=expm_multiply(-16j*space.hamiltonian,states[0])
    initial=np.array(records[0]['site_contact_probability'])
    final=np.array(records[-1]['site_contact_probability'])
    # Site contact outcomes and their one complementary separated outcome
    # form a complete distribution, with no success conditioning.
    complete_tv=(float(np.sum(abs(final-initial)))+
                 abs(records[-1]['separated_probability']-records[0]['separated_probability']))/2
    plus=space.site_index[(1,0,0)]; minus=space.site_index[(3,0,0)]
    flat=band@(a*np.exp(-16j*np.ones_like(energy)))
    return {
        'scope':'SUPPLIED_PERIODIC_ACTION_COMPLETE_N2_CLOSED_EVOLUTION',
        'side':4,'spatial_sites':64,'fermion_modes':128,'full_two_fermion_dimension':space.dimension,
        'bound_band_dimension':len(energy),'coordinate_speed':'1/8','repulsion':'1',
        'separating_gap_formula':'1-sqrt(3)/4',
        'separating_gap_rational_lower':'9/16',
        'band_energy_min':float(energy.min()),'band_energy_max':float(energy.max()),
        'contact_residue_min':float(residue.min()),'contact_residue_max':float(residue.max()),
        'complete_spatial_exchange_residual':float(exchange),
        'band_orthonormality_residual':float(np.linalg.norm(band.conj().T@band-np.eye(64))),
        'position_action_eigenstate_residual':float(np.linalg.norm(space.hamiltonian@band-band*energy)),
        'position_action_evolution_residual':float(np.linalg.norm(exact-states[-1])),
        'complete_contact_record_TV':complete_tv,
        'initial_x_neighbor_asymmetry':float(initial[plus]-initial[minus]),
        'final_x_neighbor_asymmetry':float(final[plus]-final[minus]),
        'flattened_band_density_change':float(np.max(abs(np.abs(flat)**2-np.abs(states[0])**2))),
        'records':records,
        'topology':{
            'band_complex_line_bundle':'TRIVIAL_BY_GLOBAL_POSITIVE_CONTACT_SECTION',
            'homotopy':'H_lambda=lambda*T+U*P, 0<=lambda<=1; gap>=U-2sqrt(3)c',
            'real_space_topological_defect_claimed':False,
            'Weyl_node_charge_identified_with_pair_charge':False},
        'native_all_branch_protection_proved':False,
        'packet_preparation_source_compiled':False,
        'physical_time_or_mass_calibrated':False,
        'boundary':'Finite closed N=2 periodic lattice, with supplied action and packet; not arbitrary source noise or a continuum-vacuum particle.'}
