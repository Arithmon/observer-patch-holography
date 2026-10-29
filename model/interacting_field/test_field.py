"""Independent Fock-space, CAR, source-word and full-outcome controls."""
from fractions import Fraction as F
import itertools
import math
import numpy as np
import pytest
from scipy.linalg import expm
from .lattice import TwoFermions, WeylLattice, node_table, symbol
from . import source, bounds
from .experiment import slater_witness, protocol, level, clock_obstruction
from ..common_source.chiral_candidate import gaussian_source as gaussian
from ..observer.source_scalar import compiler


def annihilators(modes):
    result=[]
    for p in range(modes):
        a=np.zeros((2**modes,2**modes),complex)
        for j in range(2**modes):
            if j & (1<<p):
                a[j^(1<<p),j]=(-1)**((j & ((1<<p)-1)).bit_count())
        result.append(a)
    return result


def embedding(space):
    out=np.zeros((2**space.modes,space.dimension),complex)
    for k,(p,q) in enumerate(zip(space.first,space.second,strict=True)):
        out[(1<<p)+(1<<q),k]=1
    return out


@pytest.mark.parametrize('p,q',[(0,1),(0,5),(2,5),(1,4)])
def test_lifted_hopping_against_independent_full_fock_and_reference(p,q):
    s=TwoFermions(6); a=annihilators(6); j=embedding(s)
    z=.3+.2j;t=.37
    h=z*a[p].conj().T@a[q]+z.conjugate()*a[q].conj().T@a[p]
    rng=np.random.default_rng(3);v=rng.normal(size=(s.dimension,3))+1j*rng.normal(size=(s.dimension,3))
    v/=np.linalg.norm(v)
    actual=source.hop(s,v,p,q,z,t)
    np.testing.assert_allclose(j@actual,expm(-1j*t*h)@j@v,atol=2e-13)
    # Wrong bosonic routing loses the occupied-spectator sign.
    if q>p+1:
        bad=np.zeros_like(h);bad=h.copy()
        for x in range(64):
            for y in range(64):
                if bad[x,y]!=0:bad[x,y]=abs(bad[x,y])
        assert np.linalg.norm(expm(-1j*t*bad)@j@v-j@actual)>1e-3


def test_second_quantization_against_dense_CAR_oracle():
    s=TwoFermions(5);a=annihilators(5);j=embedding(s)
    rng=np.random.default_rng(32);z=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5));h=z+z.conj().T
    exact=sum(h[p,q]*a[p].conj().T@a[q] for p in range(5) for q in range(5))
    np.testing.assert_allclose(j.conj().T@exact@j,s.lift_one_body(h).toarray(),atol=1e-14)


def test_actual_FSWAP_route_agrees_with_direct_CAR_lift():
    s=TwoFermions(6);rng=np.random.default_rng(11)
    state=rng.normal(size=s.dimension)+1j*rng.normal(size=s.dimension);state/=np.linalg.norm(state)
    p,q=0,5;sw=gaussian.local_gate('fswap').fock_unitary
    gate=gaussian.local_gate('hopping',.7j,.3).fock_unitary
    routed=state.copy();route=[(j-1,j) for j in range(q,p+1,-1)]
    for pair in route:routed=s.apply_pair(routed,*pair,sw)
    routed=s.apply_pair(routed,p,p+1,gate)
    for pair in reversed(route):routed=s.apply_pair(routed,*pair,sw)
    np.testing.assert_allclose(routed,s.apply_pair(state,p,q,gate),atol=5e-14)
    # The double-occupancy FSWAP minus sign is physically material.
    bad=sw.copy();bad[3,3]*=-1
    wrong=state.copy()
    for pair in route:wrong=s.apply_pair(wrong,*pair,bad)
    wrong=s.apply_pair(wrong,p,p+1,gate)
    for pair in reversed(route):wrong=s.apply_pair(wrong,*pair,bad)
    assert np.linalg.norm(wrong-routed)>.01


def test_contact_word_has_genuine_nonquadratic_phase_on_all_inputs():
    gate,info=source.contact(.25)
    np.testing.assert_allclose(gate,np.diag([1,1,1,np.exp(-.25j)]),atol=3e-14)
    # Number-preserving two-mode Gaussian gate obeys det(U_1)=U_0 U_2.
    assert abs(np.linalg.det(gate[1:3,1:3])-gate[0,0]*gate[3,3])>.2
    assert info['contracted_Fock_columns']==4


def test_non_gaussian_four_point_is_exact_and_free_control_vanishes():
    result=slater_witness()
    assert abs(result['initial_slater_defect'])<1e-12
    assert abs(result['source_final_slater_defect']-math.sin(.25)**2)<1e-12
    assert abs(result['connected_four_point']+math.sin(.25)**2/4)<1e-12
    assert result['source_final_slater_defect']>float(F(result['ideal_slater_defect_lower_exact']))
    s=TwoFermions(4);v=s.basis_state(0,1)
    for p,q in ((0,2),(1,3)):v=source.hop(s,v,p,q,-1j,math.pi/4)
    null,_=source.contact(.25,0.)
    for p,q in ((0,1),(2,3)):v=s.apply_pair(v,p,q,null)
    g=s.one_body(v)
    assert abs(2-np.trace(g@g).real)<1e-12


def test_complete_ramsey_instrument_with_complex_untouched_reference():
    s=TwoFermions(6);a=annihilators(6);j=embedding(s)
    rng=np.random.default_rng(10);v=rng.normal(size=(15,2))+1j*rng.normal(size=(15,2));v/=np.linalg.norm(v)
    p,q=1,4
    b=(a[p]-1j*a[q])/math.sqrt(2);n=b.conj().T@b
    branches=source.read_ramsey(s,v,p,q)
    for bit,w in enumerate(branches):
        np.testing.assert_allclose(j@w,(n if bit else np.eye(64)-n)@j@v,atol=4e-13)
    assert abs(sum(np.linalg.norm(w)**2 for w in branches)-1)<1e-12
    before=v.conj().T@v;after=sum(w.conj().T@w for w in branches)
    np.testing.assert_allclose(after,before,atol=1e-12)


def test_open_lattice_has_declared_bulk_Weyl_symbol_and_all_matchings():
    l=WeylLattice(1);k=(.3,-.2,.4);spin=np.array([1.,.3j]);plane=np.array([
        np.exp(1j*np.dot(k,x))*spin for x in l.sites]).reshape(-1)
    observed=(l.single_h@plane)[l.origin:l.origin+2]
    np.testing.assert_allclose(observed,symbol(k)@spin,atol=1e-14)
    for group in l.groups:
        modes=[m for p,q,_ in group for m in (p,q)]
        assert len(modes)==len(set(modes))
    assert np.linalg.norm((l.hamiltonian-l.hamiltonian.conj().T).data)==0


def test_split_and_sparse_evolution_with_entangled_reference():
    l=WeylLattice(1);v=np.zeros((l.dimension,2),complex)
    v[l.index(l.origin,l.origin+1),0]=1/math.sqrt(2)
    v[l.index(l.origin,l.neighbor),1]=1j/math.sqrt(2)
    actual=source.evolve(l,v,.0625,8);reference=l.exact_evolve(v,.0625)
    certified=float(F(bounds.history(((F(1,16),8),))['complete_history_half_trace_upper']))
    assert np.linalg.norm(actual-reference)<certified
    first=source.read_ramsey(l,actual,l.origin+1,l.neighbor+1)
    assert abs(sum(np.linalg.norm(x)**2 for x in first)-1)<1e-10


def test_complete_history_interaction_and_feedback_controls():
    row=level(1,8)
    assert len(row['complete_two_read_probabilities'])==4
    assert row['history_normalization_residual']<1e-10
    assert row['source_reference_history_TV']<2e-8
    assert row['interaction_deletion_history_TV']>.03
    assert row['feedback_deletion_history_TV']>.2
    assert row['source_short_interaction_probability_contrast']>float(F(row['short_interaction_certificate'][
        'ideal_executed_interaction_minus_free_probability_lower']))
    assert row['evolution_resources']['routing_FSWAP_calls']>0


def test_all_eight_nodes_retained_with_opposite_jacobian_chiralities():
    rows=node_table();assert len(rows)==8
    assert sum(r['jacobian_chirality'] for r in rows)==0
    for row in rows:
        k=[math.pi*x for x in row['pi_coordinates']]
        assert np.linalg.norm(symbol(k))<1e-15
        derivatives=[]
        pauli=(np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))
        for axis in range(3):
            plus=k.copy();minus=k.copy();plus[axis]+=1e-5;minus[axis]-=1e-5
            d=(symbol(plus)-symbol(minus))/2e-5
            derivatives.append([np.trace(s@d).real/2 for s in pauli])
        assert np.sign(np.linalg.det(derivatives))==row['jacobian_chirality']


def test_low_momentum_bound_does_not_claim_interacting_patch_closure():
    cert=bounds.low_momentum(F(1,4),F(1,8))
    assert cert['interaction_projection_covered'] is False
    for k in itertools.product((-.25,-.1,.1,.25),repeat=3):
        from .lattice import SIGMA
        continuum=-.125*sum(x*s for x,s in zip(k,SIGMA,strict=True))
        assert np.linalg.norm(symbol(k)-continuum,2)<float(F(cert['one_body_operator_difference_upper']))


def test_record_derived_clock_cannot_be_silently_identified_with_field():
    row=clock_obstruction()
    assert sorted(row['reference_probabilities_and_multiplicities'].values())==[1,20,45]
    assert row['stationary_target_modular_energy_levels']==3
    assert row['unrestricted_unitary_algebra_identification_possible'] is False
    assert row['embedding_or_compression_excluded'] is False


@pytest.mark.parametrize('bad',[True,float('nan'),float('inf'),'1'])
def test_invalid_source_parameters_rejected_even_after_cache_warmup(bad):
    source.contact(1.,1.)
    with pytest.raises(ValueError):source.contact(bad,1.)
    with pytest.raises(ValueError):source.contact(1.,bad)


def test_reject_number_leaking_gate_and_bad_domains():
    s=TwoFermions(4);v=s.basis_state(0,1);bad=np.eye(4);bad[0,3]=.1
    with pytest.raises(ValueError):s.apply_pair(v,0,1,bad)
    with pytest.raises(ValueError):WeylLattice(True)
    with pytest.raises(ValueError):bounds.history(((F(1,16),True),))
    with pytest.raises(ValueError):bounds.interaction_witness(interaction=F(0))
    with pytest.raises(ValueError):s.state(np.full(6,np.nan))


def test_record_inferred_modular_flow_remains_gaussian_under_native_occupation_dictionary():
    row=clock_obstruction()['occupation_preserving_CAR_dictionary']
    assert row['all_positive_port_weights_covered'] is True
    assert row['contact_interaction_generated_by_this_clock'] is False
    for control in row['controls']:
        assert abs(control['source_executed_modular_Slater_defect'])<1e-11


def test_cached_local_source_outputs_cannot_be_mutated():
    gate,info=source.contact(.25)
    with pytest.raises(ValueError):gate.setflags(write=True)
    info['resources']['native_pulses']=-1
    info['word'].clear()
    assert source.contact(.25)[1]['resources']['native_pulses']>0
    assert source.contact(.25)[1]['word']
    kraus,info=source.occupation()
    with pytest.raises(ValueError):kraus[0].setflags(write=True)
    info['resources'].clear()
    assert source.occupation()[1]['resources']['native_pulses']>0


def test_wrong_source_contact_word_cannot_reuse_a_target_certificate(monkeypatch):
    source._contact.cache_clear()
    monkeypatch.setattr(compiler,'execute',lambda word,qubits:np.eye(2**qubits,dtype=complex))
    with pytest.raises(ArithmeticError):source.contact(.37,1.)
    source._contact.cache_clear()


def test_zero_kinetic_control_remains_a_well_defined_contact_process():
    l=WeylLattice(1,speed=0.)
    assert sum(map(len,l.groups))==0
    initial=l.basis_state(l.origin,l.origin+1)
    np.testing.assert_allclose(l.exact_evolve(initial,.2),np.exp(-.2j)*initial,atol=1e-14)
    np.testing.assert_allclose(source.evolve(l,initial,.2,3),np.exp(-.2j)*initial,atol=1e-12)


def test_complete_quantum_distance_stable_for_nearly_equal_outputs():
    from .experiment import pure_cq_distance
    a=np.array([1.,0.,0.],complex)
    b=np.array([math.sqrt(1-1e-18),1e-9,0.],complex)
    assert abs(pure_cq_distance([a],[b])-1e-9)<1e-20
    rng=np.random.default_rng(34)
    a=(rng.normal(size=7)+1j*rng.normal(size=7))/5
    b=(rng.normal(size=7)+1j*rng.normal(size=7))/5
    exact=np.linalg.svd(np.outer(a,a.conj())-np.outer(b,b.conj()),compute_uv=False).sum()/2
    assert abs(pure_cq_distance([a],[b])-exact)<1e-14
