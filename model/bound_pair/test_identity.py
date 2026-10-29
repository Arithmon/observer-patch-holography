"""Independent complete-sector current and full-Fock identity checks."""
import numpy as np
import pytest
from scipy import sparse
from ..interacting_field.lattice import WeylLattice
from ..field_stress.energy import GaussianInteger,commutator,sum_exact
from .identity import (H_SCALE,full_fock_annihilators,car_pair_certificate,
                       number_currents,oriented_bond,spin_certificate,report)


@pytest.fixture(scope='module')
def packet():
    lattice=WeylLattice(1)
    return lattice,number_currents(lattice)


def test_complete_fock_pair_statistics_and_number_changing_countercontrols():
    cert=car_pair_certificate()
    assert cert['identities_exact']
    assert cert['exact_integer_identities_checked']==38
    assert cert['particle_number_charge_of_pair_creation']==2
    assert cert['pair_drive_preserves_parity']
    assert cert['pair_drive_number_commutator_HS_squared']>0
    assert cert['single_fermion_drive_parity_commutator_HS_squared']>0
    # Hard-core commutator is not the canonical boson CCR, on all16 inputs.
    a=full_fock_annihilators()
    bd=a[0].T@a[1].T
    assert not np.array_equal(bd.T@bd-bd@bd.T,np.eye(16))


def test_exact_number_currents_and_region_balance_on_full_action(packet):
    lattice,(cert,local,currents)=packet
    assert cert['complete_two_fermion_dimension']==1431
    assert cert['local_continuity_equations_checked']==27
    assert cert['all_unordered_site_pairs_checked']==351
    assert cert['nonzero_bond_currents']==54
    assert cert['crossing_number_currents']==9
    assert (sum_exact(local,lattice.dimension)-GaussianInteger.from_exact(
        2*sparse.eye(lattice.dimension))).is_zero()
    assert all(not value.is_zero() for value in currents.values())


def test_contact_deletion_leaves_number_current_operators_identical(packet):
    _,(_,_,currents)=packet
    _,_,free=number_currents(WeylLattice(1,interaction=0))
    assert currents.keys()==free.keys()
    assert all((value-free[key]).is_zero() for key,value in currents.items())


def test_number_flux_has_independent_evolution_and_probe_orientation(packet):
    lattice,(_,local,currents)=packet
    x=lattice.site_index[(0,0,0)]
    y=lattice.site_index[(1,0,0)]
    forward=oriented_bond(lattice,x,y)
    current=currents[x,y].matrix(H_SCALE)
    # A Peierls link phase couples to this current, not its negative.
    dt=1e-5
    plus=np.exp(1j*dt)*forward+np.exp(-1j*dt)*forward.getH()
    minus=np.exp(-1j*dt)*forward+np.exp(1j*dt)*forward.getH()
    derivative=(plus-minus)/(2*dt)
    assert np.max(abs((derivative-current).data),initial=0)<2e-12
    # One known spin-changing hopping transition has a nonzero number flux.
    psi=lattice.basis_state(lattice.origin,lattice.origin+1)
    psi=psi+1j*lattice.basis_state(lattice.origin,lattice.neighbor)
    psi/=np.linalg.norm(psi)
    # Select phase independently so this imaginary hopping bond carries flux.
    if abs(np.vdot(psi,current@psi))<1e-4:
        psi=lattice.basis_state(lattice.origin,lattice.origin+1)
        psi=psi+lattice.basis_state(lattice.origin,lattice.neighbor)
        psi/=np.linalg.norm(psi)
    assert abs(np.vdot(psi,current@psi))>.01
    flux=GaussianInteger.zero(lattice.dimension)
    for (p,q),value in currents.items():
        if p==x: flux=flux+value
        elif q==x: flux=flux-value
    before,after=lattice.exact_evolve(psi,-dt),lattice.exact_evolve(psi,dt)
    nx=local[x].matrix()
    actual=(np.vdot(after,nx@after)-np.vdot(before,nx@before)).real/(2*dt)
    predicted=-np.vdot(psi,flux.matrix(H_SCALE)@psi).real
    assert abs(actual-predicted)<2e-10


def test_site_phase_covariance_retains_same_local_contact_action(packet):
    lattice,(_,local,_)=packet
    angles=np.linspace(-.3,.8,lattice.volume)
    state_phases=np.exp(1j*(angles[lattice.first//2]+angles[lattice.second//2]))
    one=lattice.single_h.tocoo()
    changed=sparse.coo_matrix((one.data*np.exp(1j*(angles[one.row//2]-angles[one.col//2])),
                              (one.row,one.col)),shape=one.shape)
    gauged=lattice.lift_one_body(changed)
    diagonal=np.zeros(lattice.dimension); diagonal[lattice.contact_indices]=1
    gauged+=sparse.diags(diagonal)
    original=lattice.hamiltonian.tocoo()
    conjugated=sparse.coo_matrix((original.data*state_phases[original.row]*state_phases[original.col].conj(),
                                 (original.row,original.col)),shape=original.shape).tocsr()
    assert np.max(abs((gauged-conjugated).data),initial=0)<4e-15


def test_changed_action_is_rejected_instead_of_inheriting_conservation():
    lattice=WeylLattice(1)
    lattice.hamiltonian+=sparse.csr_matrix(([1/16],([0],[1])),shape=lattice.hamiltonian.shape)
    with pytest.raises(ArithmeticError,match='action differs'):
        number_currents(lattice)


def test_spin_singlet_at_eta_is_not_a_conserved_spin_species(packet):
    lattice,_=packet
    cert=spin_certificate(lattice)
    assert cert['onsite_pair_internal_spin_singlet_exact']
    assert not cert['standalone_internal_spin_conserved']
    eta,away=cert['bound_band_spin_content']
    assert eta['triplet_weight']<1e-28
    assert away['triplet_weight']>.01
    assert away['fermion_exchange_residual']<1e-13
    assert not cert['physical_Lorentz_spin_assigned']


def test_identity_report_does_not_assign_unproved_physical_species():
    cert=report()
    for field in ('electromagnetic_charge','Standard_Model_representation',
                  'family_generation','physical_Lorentz_spin','real_space_topological_defect_charge'):
        assert cert['assignment'][field] is None
    assert not cert['source_selected_CAR_statistics']
    assert cert['physical_comparison_count_added']==0
