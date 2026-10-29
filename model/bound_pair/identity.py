"""Actual charge, current and spin maps of the supplied Weyl-contact pair.

The conserved U(1) is CAR particle number. Neither its name nor its integer
normalization identifies electromagnetic charge, Standard Model matter or
physical Lorentz spin. All derivatives use the declared Hamiltonian model time.
"""
from __future__ import annotations
from fractions import Fraction as F
import numpy as np
from scipy import sparse
from ..interacting_field.lattice import WeylLattice, SIGMA
from ..field_stress.energy import GaussianInteger, commutator, sum_exact
from .fiber import PairFiber, SINGLET

H_SCALE = 16


def full_fock_annihilators(modes=4):
    """Independent occupation-bit CAR matrices; no N=2 projection involved."""
    if type(modes) is not int or not 1 <= modes <= 8:
        raise ValueError('1..8 full Fock modes required')
    operators=[]
    for p in range(modes):
        matrix=np.zeros((2**modes,2**modes),dtype=np.int64)
        for mask in range(2**modes):
            if mask & (1 << p):
                sign=(-1)**((mask & ((1 << p)-1)).bit_count())
                matrix[mask^(1 << p),mask]=sign
        operators.append(matrix)
    return tuple(operators)


def comm(first,second):
    return first@second-second@first


def car_pair_certificate():
    a=full_fock_annihilators()
    ident=np.eye(16,dtype=np.int64)
    n=[op.T@op for op in a]
    number=sum(n)
    parity=np.diag([(-1)**mask.bit_count() for mask in range(16)])
    bdag=a[0].T@a[1].T
    b=bdag.T
    other=a[2].T@a[3].T
    equations={
        'charge_two':comm(number,bdag)-2*bdag,
        'even_parity':comm(parity,bdag),
        'hard_core':bdag@bdag,
        'hard_core_commutator':comm(b,bdag)-(ident-n[0]-n[1]),
        'disjoint_pair_creation_commutes':comm(bdag,other),
        'disjoint_pair_annihilation_commutes':comm(bdag,other.T),
    }
    for p in range(4):
        for q in range(4):
            equations[f'CAR_annihilator_{p}_{q}']=a[p]@a[q]+a[q]@a[p]
            equations[f'CAR_mixed_{p}_{q}']=a[p]@a[q].T+a[q].T@a[p]-(ident if p==q else 0)
    if not all(np.count_nonzero(value)==0 for value in equations.values()):
        raise ArithmeticError('full-Fock CAR/pair identity failed')
    drive=bdag+b
    odd_drive=a[0].T+a[0]
    return {'full_Fock_dimension':16,'exact_integer_identities_checked':len(equations),
            'particle_number_charge_of_pair_creation':2,
            'pair_has_even_fermion_parity':True,
            'disjoint_pair_operators_commute':True,
            'pair_creation_is_hard_core':True,
            'canonical_boson_CCR_claimed':False,
            'pair_drive_preserves_parity':not bool(np.count_nonzero(comm(parity,drive))),
            'pair_drive_number_commutator_HS_squared':int(np.sum(comm(number,drive)**2)),
            'single_fermion_drive_parity_commutator_HS_squared':int(np.sum(comm(parity,odd_drive)**2)),
            'identities_exact':True}


def local_numbers(lattice):
    """Integer N_x on the full retained two-fermion sector."""
    return [GaussianInteger.from_exact(sparse.diags(
        (lattice.first//2==site).astype(int)+(lattice.second//2==site).astype(int)))
        for site in range(lattice.volume)]


def oriented_bond(lattice,x,y):
    """A_xy=sum_(a,b) c_(x,a)^* h_(xa,yb) c_(y,b); direction x<-y."""
    if type(x) is not int or type(y) is not int or not 0 <= x < y < lattice.volume:
        raise ValueError('ordered distinct lattice sites required')
    single=lattice.single_h.tocoo()
    take=(single.row//2==x)&(single.col//2==y)
    h=sparse.coo_matrix((single.data[take],(single.row[take],single.col[take])),shape=single.shape)
    return lattice.lift_one_body(h)


def number_currents(lattice=None):
    """Return independently assembled J_x->y=i(A_xy-A_xy*) and exact certificate."""
    lattice=WeylLattice(1) if lattice is None else lattice
    if lattice.speed!=1/8 or lattice.interaction not in (0.,1.):
        raise ValueError('exact certificate covers c=1/8 and U=0 or 1')
    local=local_numbers(lattice)
    total_number=sum_exact(local,lattice.dimension)
    twice=GaussianInteger.from_exact(2*sparse.eye(lattice.dimension))
    if not (total_number-twice).is_zero():
        raise ArithmeticError('local numbers do not sum to N=2')
    h=GaussianInteger.from_exact(H_SCALE*lattice.hamiltonian)
    hopping=GaussianInteger.zero(lattice.dimension)
    flux=[GaussianInteger.zero(lattice.dimension) for _ in local]
    currents={}
    for x in range(lattice.volume):
        for y in range(x+1,lattice.volume):
            forward=GaussianInteger.from_exact(H_SCALE*oriented_bond(lattice,x,y))
            hopping=hopping+forward+forward.adjoint()
            current=(forward-forward.adjoint()).times_i()
            if current.is_zero():
                continue
            if not (current-current.adjoint()).is_zero():
                raise ArithmeticError('non-Hermitian number current')
            if sum(abs(a-b) for a,b in zip(lattice.sites[x],lattice.sites[y],strict=True))!=1:
                raise ArithmeticError('number current is not a nearest-neighbor bond')
            currents[x,y]=current
            flux[x]=flux[x]+current
            flux[y]=flux[y]-current
    interaction=h-hopping
    expected_contact=np.zeros(lattice.dimension)
    expected_contact[lattice.contact_indices]=H_SCALE*lattice.interaction
    if not (interaction-GaussianInteger.from_exact(sparse.diags(expected_contact))).is_zero():
        raise ArithmeticError('action differs from declared hopping/contact terms')
    for nx,divergence in zip(local,flux,strict=True):
        if not commutator(interaction,nx).is_zero():
            raise ArithmeticError('contact violates local number conservation')
        if not (commutator(h,nx).times_i()+divergence).is_zero():
            raise ArithmeticError('exact number continuity failed')
    if not sum_exact(flux,lattice.dimension).is_zero():
        raise ArithmeticError('number boundary flux failed to cancel')
    region={i for i,x in enumerate(lattice.sites) if x[0]<=0}
    regional=sum_exact([local[x] for x in region],lattice.dimension)
    boundary=GaussianInteger.zero(lattice.dimension)
    crossing=0
    for (x,y),value in currents.items():
        if (x in region)!=(y in region):
            boundary=boundary+value if x in region else boundary-value
            crossing+=1
    if not (commutator(h,regional).times_i()+boundary).is_zero():
        raise ArithmeticError('regional number continuity failed')
    certificate={'sites':lattice.volume,'complete_two_fermion_dimension':lattice.dimension,
                 'all_unordered_site_pairs_checked':lattice.volume*(lattice.volume-1)//2,
                 'local_continuity_equations_checked':lattice.volume,
                 'nonzero_bond_currents':len(currents),'maximum_current_range_lattice_steps':1,
                 'current_denominator':H_SCALE,'total_number_exact':2,
                 'contact_commutes_with_each_local_number':True,
                 'region_sites':len(region),'crossing_number_currents':crossing,
                 'identities_exact':True,'time_parameter':'Hamiltonian model time, supplied hbar=1'}
    return certificate,local,currents


def hs_squared_exact(matrix,denominator):
    return str(F(int(np.sum(matrix.real.data**2)+np.sum(matrix.imag.data**2)),denominator**2))


def spin_certificate(lattice=None):
    lattice=WeylLattice(1) if lattice is None else lattice
    h=GaussianInteger.from_exact(H_SCALE*lattice.hamiltonian)
    rows=[]
    for axis,sigma in enumerate(SIGMA):
        spin=lattice.lift_one_body(sparse.kron(sparse.eye(lattice.volume),sigma))
        twice_spin=GaussianInteger.from_exact(spin)
        changed=commutator(h,twice_spin)
        # Independent full-Fock two-site spin computation uses sigma without
        # the factor1/2, so every identity still has Gaussian-integer entries.
        a=full_fock_annihilators()
        local_spin=sum(sigma[p,q]*(a[p].T@a[q]) for p in range(2) for q in range(2))
        pair=a[0].T@a[1].T
        if np.count_nonzero(comm(local_spin,pair)):
            raise ArithmeticError('onsite pair is not an internal spin singlet')
        total_spin=np.kron(sigma,np.eye(2))+np.kron(np.eye(2),sigma)
        if np.linalg.norm(total_spin@SINGLET)>1e-14:
            raise ArithmeticError('fiber singlet dictionary failed')
        if changed.is_zero():
            raise ArithmeticError('declared spin-orbit action unexpectedly conserves standalone spin')
        rows.append({'axis':axis,'commutator_H_S_HS_squared_exact':hs_squared_exact(changed,2*H_SCALE)})
    bands=[]
    for label,momentum in [('eta',[np.pi]*3),('away_eta',[0,np.pi,np.pi])]:
        fiber=PairFiber(momentum,side=4)
        vector=fiber.eigenstate().reshape(-1,4)
        singlet=(vector@SINGLET.conj())[:,None]*SINGLET
        triplet=vector-singlet
        bands.append({'name':label,'momentum_over_pi':[float(k/np.pi) for k in momentum],
                      'periodic_side':4,'triplet_weight':float(np.vdot(triplet,triplet).real),
                      'norm_residual':abs(float(np.vdot(vector,vector).real)-1),
                      'fermion_exchange_residual':fiber.exchange_residual(vector.reshape(-1)),
                      'pole_energy':fiber.pole()['energy']})
    return {'onsite_pair_internal_spin_singlet_exact':True,
            'standalone_internal_spin_conserved':False,
            'full_action_spin_commutators':rows,'bound_band_spin_content':bands,
            'physical_Lorentz_spin_assigned':False}


def report():
    lattice=WeylLattice(1)
    current,_,_=number_currents(lattice)
    return {'scope':'same supplied Weyl-contact action; CAR identity and conserved number current',
            'car_pair':car_pair_certificate(),'number_current':current,
            'spin':spin_certificate(lattice),
            'assignment':{'CAR_particle_number':2,'fermion_parity':'even',
                          'composite_exchange_algebra':'disjoint onsite pair operators commute; hard-core onsite',
                          'electromagnetic_charge':None,'Standard_Model_representation':None,
                          'family_generation':None,'physical_Lorentz_spin':None,
                          'real_space_topological_defect_charge':None},
            'electromagnetic_identification_supplied':False,
            'source_selected_CAR_statistics':False,'physical_comparison_count_added':0}
