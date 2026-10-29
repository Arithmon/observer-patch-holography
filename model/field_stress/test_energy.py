"""Hostile checks of the exact local-energy map and source evolution."""
from fractions import Fraction
import numpy as np
import pytest
from scipy import sparse
from ..interacting_field.lattice import WeylLattice
from .energy import GaussianInteger, SCALE, certify, commutator, sum_exact
from .experiment import dynamics


@pytest.fixture(scope='module')
def packet():
    lattice = WeylLattice(1)
    receipt,local,currents = certify(lattice)
    return lattice,receipt,local,currents


def test_complete_exact_operator_partition_and_continuity(packet):
    lattice,receipt,local,currents = packet
    assert receipt['operator_identities_exact']
    assert receipt['all_unordered_site_pairs_checked'] == 351
    assert receipt['local_continuity_equations_checked'] == 27
    assert receipt['nonzero_unordered_currents'] == 81
    assert receipt['maximum_current_range_in_lattice_steps'] == 2
    assert receipt['boundary_current_pairs'] == 18
    total = sum_exact(local,lattice.dimension)
    assert (total-GaussianInteger.from_exact(SCALE*lattice.hamiltonian)).is_zero()
    assert not currents[min(currents)].is_zero()


def test_dropping_range_two_boundary_currents_breaks_regional_balance(packet):
    lattice,_,local,currents = packet
    region = {i for i,x in enumerate(lattice.sites) if x[0] <= 0}
    nearest_only = GaussianInteger.zero(lattice.dimension)
    for (x,y),current in currents.items():
        distance = sum(abs(a-b) for a,b in zip(lattice.sites[x],lattice.sites[y],strict=True))
        if distance == 1 and ((x in region) != (y in region)):
            nearest_only = nearest_only+current if x in region else nearest_only-current
    total = sum_exact(local,lattice.dimension)
    regional = sum_exact([local[x] for x in region],lattice.dimension)
    assert not (commutator(total,regional).times_i()+nearest_only).is_zero()


def test_exact_map_rejects_a_changed_global_action():
    lattice = WeylLattice(1)
    perturbation = sparse.csr_matrix(([1/32],([0],[0])),shape=lattice.hamiltonian.shape)
    lattice.hamiltonian = lattice.hamiltonian+perturbation
    with pytest.raises(ArithmeticError,match='independently assembled action'):
        certify(lattice)


def test_current_sign_matches_independent_time_derivative(packet):
    lattice,_,local,currents = packet
    state = lattice.basis_state(lattice.origin,lattice.origin+1)
    state += lattice.basis_state(lattice.origin,lattice.neighbor)
    state /= np.sqrt(2)
    site = lattice.site_index[(0,0,0)]
    flux = sum_exact([commutator(local[site],other).times_i() for other in local],lattice.dimension)
    operator = local[site].matrix(SCALE)
    dt = 1e-5
    before,after = lattice.exact_evolve(state,-dt),lattice.exact_evolve(state,dt)
    derivative = (np.vdot(after,operator@after)-np.vdot(before,operator@before)).real/(2*dt)
    predicted = -np.vdot(state,flux.matrix(SCALE**2)@state).real
    assert abs(predicted) > .01
    assert abs(derivative-predicted) < 1e-9


def test_interaction_deletion_changes_actual_energy_transfer(packet):
    lattice,_,local,currents = packet
    _,_,free = certify(WeylLattice(1,interaction=0))
    changed = [key for key,current in currents.items() if key not in free or not (current-free[key]).is_zero()]
    assert len(changed) == 54
    for x,y in changed:
        assert sum(abs(a-b) for a,b in zip(lattice.sites[x],lattice.sites[y],strict=True)) == 1


@pytest.mark.parametrize('value',[.1,1j/3,float('inf')])
def test_integer_certificate_does_not_round_noninteger_inputs(value):
    with pytest.raises(ValueError):
        GaussianInteger.from_exact(sparse.csr_matrix([[value]]))


def test_integer_multiplication_bound_prevents_overflow():
    big = GaussianInteger.from_exact(sparse.eye(4)*2**30)
    with pytest.raises(OverflowError):
        _ = big@big


def test_source_evolution_retains_finite_energy_drift_and_resources(packet):
    lattice,_,local,currents = packet
    receipt = dynamics(lattice,local,currents)
    assert abs(receipt['initial_origin_outward_energy_flux']) > .01
    assert receipt['sparse_exact_energy_drift'] < 2e-13
    assert receipt['preparation_resources']['routing_FSWAP_calls'] > 0
    for row in receipt['source_runs']:
        assert row['source_energy_drift'] < float(Fraction(row['ideal_energy_drift_upper']))
        assert row['source_state_to_sparse_exact_norm'] < float(Fraction(row['ideal_unitary_error_upper']))
        assert row['requested_evolution_resources']['routing_FSWAP_calls'] > 0
        assert not row['outward_certified_floating_or_hardware_error']
