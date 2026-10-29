"""Independent full-sector, full-local-word and effective-current controls."""
from itertools import product
import json

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import eigh

from ..bound_pair.identity import full_fock_annihilators
from ..common_source.chiral_candidate import gaussian_source
from ..interacting_field import source
from ..interacting_field.lattice import WeylLattice
from .certificate import exact_certificate
from .links import (GaugedField, contact_blocks, determinant_pair_action,
                    second_order_pair, feshbach, feshbach_remainder_bound,
                    pair_current_response)


SQUARE = tuple(product(range(2),range(2),(0,)))


def maximum(matrix):
    return float(np.max(abs(matrix.data),initial=0)) if sparse.issparse(matrix) else float(np.max(abs(matrix),initial=0))


def test_zero_links_are_the_existing_full_action_without_changed_couplings():
    actual,old=GaugedField(),WeylLattice(1)
    assert actual.sites == old.sites and actual.groups == old.groups
    assert maximum(actual.hamiltonian-old.hamiltonian) == 0


def test_exact_nonflat_link_ward_and_determinant_certificate():
    result=exact_certificate()
    assert json.loads(json.dumps(result)) == result
    assert result['complete_N2_dimension'] == 1431
    assert result['exact_local_Ward_identities'] == 27
    assert result['nonzero_bond_currents'] == 54
    assert result['nontrivial_square_holonomies'] > 0
    assert not result['physical_electromagnetic_charge_assigned']


def test_local_car_two_virtual_paths_fix_the_determinant_sign():
    # Unrelated generic integer bond, including spin mixing, on all16 Fock
    # columns. This oracle does not use TwoFermions' CAR-routing formula.
    a=full_fock_annihilators(4)
    h=np.array([[1+2j,3-1j],[-2+1j,2-3j]])
    forward=sum(h[p,q]*(a[p].T@a[q+2]) for p,q in product(range(2),repeat=2))
    hopping=forward+forward.conj().T
    # c0*c1*|vac> is mask3, c2*c3*|vac> is mask12 in ascending convention.
    assert (hopping@hopping)[3,12] == 2*(h[0,0]*h[1,1]-h[0,1]*h[1,0])
    assert abs((hopping@hopping)[3,3]-np.vdot(h,h)) < 1e-12


def test_nonflat_arbitrary_links_preserve_full_covariance_and_local_Ward():
    field=GaugedField(SQUARE,angles={(0,1):.31,(0,2):-.47,(1,3):.28})
    changed,phases=field.gauge_transform([.2,-.8,.7,1.1])
    g=sparse.diags(phases)
    assert maximum(changed.hamiltonian-g@field.hamiltonian@g.conj().T) < 1e-15
    n=field.local_numbers()
    flux=[sparse.csr_matrix(field.hamiltonian.shape,dtype=complex) for _ in n]
    for edge in field.axes:
        j=field.current(edge);x,y=edge
        flux[x]+=j;flux[y]-=j
        assert maximum(changed.current(edge)-g@j@g.conj().T) < 1e-15
    for nx,div in zip(n,flux,strict=True):
        assert maximum(1j*(field.hamiltonian@nx-nx@field.hamiltonian)+div) < 1e-15
    assert maximum(sum(flux)) < 1e-15


def test_link_derivative_is_the_same_oriented_complete_current():
    edge=(0,1);a=.37;eps=1e-5
    field=GaugedField(SQUARE,angles={edge:a})
    plus=GaugedField(SQUARE,angles={edge:a+eps})
    minus=GaugedField(SQUARE,angles={edge:a-eps})
    assert maximum((plus.hamiltonian-minus.hamiltonian)/(2*eps)-field.current(edge)) < 2e-12


def test_completed_local_source_word_is_gauge_covariant_on_all_Fock_inputs():
    alpha,beta=.37,-.81
    z=.07+.09j;duration=.43
    initial=gaussian_source.local_gate('hopping',z,duration)
    changed=gaussian_source.local_gate('hopping',z*np.exp(1j*(alpha-beta)),duration)
    g=np.diag(np.exp(1j*np.array([0,alpha,beta,alpha+beta])))
    assert maximum(changed.fock_unitary-g@initial.fock_unitary@g.conj().T) < 2e-12
    # Link phase changes a coefficient of an admitted word, not its primitive set.
    for key in source.COUNT_KEYS:
        assert initial.resources[key] == changed.resources[key]


def test_complete_source_history_covariance_keeps_all_pointer_and_reference_outputs():
    field=GaugedField(SQUARE,angles={(0,1):.4,(0,2):-.2})
    changed,phase=field.gauge_transform([.1,-.7,.4,1.3])
    rng=np.random.default_rng(851)
    state=rng.normal(size=(field.dimension,3))+1j*rng.normal(size=(field.dimension,3))
    state/=np.linalg.norm(state)
    output=source.evolve(field,state,.31,3)
    gauged=source.evolve(changed,phase[:,None]*state,.31,3)
    assert np.linalg.norm(gauged-phase[:,None]*output)<2e-11
    branches=source.read_number(field,output,0)
    changed_branches=source.read_number(changed,gauged,0)
    for left,right in zip(branches,changed_branches,strict=True):
        assert np.linalg.norm(right-phase[:,None]*left)<2e-11
    assert np.linalg.norm(sum(b.conj().T@b for b in branches)-state.conj().T@state)<2e-11
    assert abs(sum(np.vdot(b,b).real for b in branches)-1)<2e-11
    first,second=source.resources(field,.31,3),source.resources(changed,.31,3)
    assert all(first[k]==second[k] for k in source.COUNT_KEYS)


def test_all_virtual_scattering_states_give_charge_two_pair_links():
    field=GaugedField(SQUARE,angles={(0,1):.37,(0,2):-.21})
    actual=second_order_pair(field).toarray()
    assert maximum(actual-determinant_pair_action(field))<1e-15
    for x,y in field.axes:
        expected=field.speed**2/(2*field.interaction)*np.exp(2j*field.angles[x,y])
        assert abs(actual[x,y]-expected)<1e-15
    changed,_=field.gauge_transform([.7,-.4,.1,.2])
    g=np.diag(np.exp(2j*np.array([.7,-.4,.1,.2])))
    assert maximum(second_order_pair(changed).toarray()-g@actual@g.conj().T)<1e-15


def test_exact_contact_resolvent_and_uniform_second_order_remainder():
    field=GaugedField(SQUARE,angles={(0,1):.37})
    p,q,b,d=contact_blocks(field)
    energy=1.4
    h=field.hamiltonian.toarray()
    projected=np.linalg.inv(energy*np.eye(field.dimension)-h)[np.ix_(p,p)]
    effective=feshbach(field,energy)
    assert np.linalg.norm(projected-np.linalg.inv(energy*np.eye(field.volume)-effective))<2e-14
    alpha=np.array([.2,.6,-.4,.9])
    changed,_=field.gauge_transform(alpha)
    pair_g=np.diag(np.exp(2j*alpha))
    assert np.linalg.norm(feshbach(changed,energy)-pair_g@effective@pair_g.conj().T)<2e-14
    error=np.linalg.norm(effective-second_order_pair(field).toarray(),2)
    assert 0 < error <= feshbach_remainder_bound(field,energy)
    eigs,vectors=eigh(h)
    assert np.count_nonzero(eigs>=field.interaction-1e-12)==field.volume
    assert eigs[-field.volume]>=field.interaction-1e-12
    assert eigs[-field.volume-1]<=field.kinetic_norm_upper()+1e-12
    for energy,psi in zip(eigs[-field.volume:],vectors[:,-field.volume:].T,strict=True):
        assert np.linalg.norm(feshbach(field,energy)@psi[p]-energy*psi[p])<2e-14
        assert min(abs(np.linalg.eigvalsh(second_order_pair(field).toarray())-energy)) <= feshbach_remainder_bound(field,energy)


def test_full_pair_current_requires_the_resolvent_norm_metric():
    edge=(0,1);angle=.37
    field=GaugedField(SQUARE,angles={edge:angle})
    energies,vectors=eigh(field.hamiltonian.toarray())
    energy,psi=energies[-1],vectors[:,-1]
    p,q,_,_=contact_blocks(field)
    u=psi[p]/np.linalg.norm(psi[p])
    derivative,metric,x=pair_current_response(field,energy,edge)
    norm=np.vdot(u,metric@u).real
    exact=np.vdot(psi,field.current(edge)@psi).real
    response=np.vdot(u,derivative@u).real/norm
    assert norm>1.001 and abs(exact)>1e-5
    assert abs(response-exact)<2e-15
    assert abs(np.linalg.norm(psi[p])**2-1/norm)<2e-14
    # Dropping the separated-state norm gives a detectably different current.
    assert abs(np.vdot(u,derivative@u).real-exact)>1e-6
    eps=1e-5
    plus=np.linalg.eigvalsh(GaugedField(SQUARE,angles={edge:angle+eps}).hamiltonian.toarray())[-1]
    minus=np.linalg.eigvalsh(GaugedField(SQUARE,angles={edge:angle-eps}).hamiltonian.toarray())[-1]
    assert abs((plus-minus)/(2*eps)-exact)<5e-10


def test_half_flux_periodicity_of_leading_pair_action_is_not_full_field_periodicity():
    a=GaugedField(SQUARE,angles={(0,1):.37})
    b=GaugedField(SQUARE,angles={(0,1):.37+np.pi})
    assert maximum(second_order_pair(a)-second_order_pair(b))<1e-15
    assert np.linalg.norm(np.linalg.eigvalsh(a.hamiltonian.toarray())-
                          np.linalg.eigvalsh(b.hamiltonian.toarray()))>1e-3


@pytest.mark.parametrize('options',[
    {'sites':[(0,0,0),(0,0,0)]}, {'sites':[(False,0,0),(1,0,0)]},
    {'sites':SQUARE,'angles':{(0,3):.1}}, {'sites':SQUARE,'angles':{(False,1):.1}},
    {'sites':SQUARE,'angles':{(0,1):float('nan')}},
    {'sites':SQUARE,'angles':{(0,1):True}},
])
def test_invalid_geometry_and_links_fail_before_evolution(options):
    with pytest.raises(ValueError):GaugedField(**options)


def test_outside_resolvent_and_zero_interaction_do_not_receive_pair_certificate():
    field=GaugedField(SQUARE)
    for energy in (field.kinetic_norm_upper(),-.1,float('nan'),True):
        with pytest.raises(ValueError):feshbach(field,energy)
    free=GaugedField(SQUARE,interaction=0)
    with pytest.raises(ValueError):second_order_pair(free)
    with pytest.raises(ValueError):feshbach_remainder_bound(free,2.)


def test_resolvent_host_limit_does_not_exclude_the_algebraic_field():
    large=GaugedField(tuple((2*x,0,0) for x in range(65)))
    assert large.dimension == 8385
    assert second_order_pair(large).shape == (65,65)
    with pytest.raises(ValueError,match='host resolvent'):
        feshbach(large,2.)
    with pytest.raises(ValueError,match='host resolvent'):
        pair_current_response(large,2.,(0,1))
