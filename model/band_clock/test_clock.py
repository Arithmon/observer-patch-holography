"""Independent full-action, source, reference and stop controls."""
from fractions import Fraction as F
import math
import numpy as np
import pytest
from scipy import sparse
from . import clock
from .experiment import retained_counts, run
from ..interacting_field import source


def test_periodic_wrap_signs_and_all_matching_norms():
    space=clock.ClockField()
    assert [len(group) for group in space.groups]==[64]*6
    rebuilt=sparse.csr_matrix(space.single_h.shape,dtype=complex)
    for group in space.groups:
        rows=[];cols=[];values=[]
        for p,q,z in group:
            assert p<q
            rows.extend((p,q));cols.extend((q,p));values.extend((z,z.conjugate()))
        h=sparse.coo_matrix((values,(rows,cols)),shape=space.single_h.shape).tocsr()
        assert np.linalg.norm(h.toarray(),2)<=1/16+1e-14
        rebuilt+=h
    assert np.max(abs((rebuilt-space.single_h).data),initial=0)==0
    assert np.max(abs((space.hamiltonian-space.free_h-sparse.diags(
        np.isin(np.arange(space.dimension),space.contact_indices).astype(float))).data),initial=0)==0


def test_common_full_reference_includes_cross_momentum_observables():
    result=clock.reference_report()
    assert result['dimension']==8128
    assert sum(result['fiber_dimensions'])==8128
    assert result['minimum_probability']>0
    assert abs(result['full_trace']-1)<1e-12
    assert result['full_algebra_modular_Hamiltonian_unitary_residual']<2e-14
    assert result['independently_normalized_fibers_cross_coherence_error']>.05
    assert not result['Gibbs_reference_prepared_by_native_repair']


def test_source_whole_field_packet_and_untouched_reference():
    report=run()
    output=report['source_execution']
    assert output['source_vs_same_rate_exact_vector_norm']<3e-5
    assert output['independent_position_vs_momentum_evolution_norm']<2e-12
    assert output['complete_CQ_half_trace_vs_true_rate']<.001
    assert output['untouched_reference_marginal_residual']<2e-8
    assert abs(sum(output['complete_pointer_probabilities'])-1)<2e-8
    assert output['complete_contact_site_plus_separated_motion_TV']>1e-4
    assert output['zero_clock_complete_CQ_contrast']>.001
    assert output['deleted_contact_complete_CQ_contrast']>.1
    assert output['requested_resources']['routing_FSWAP_calls']>0
    assert not report['unfulfilled_joins']['new_full_Gibbs_reference_prepared_or_selected']


def test_local_source_gate_mutation_changes_full_history(monkeypatch):
    counts=retained_counts();space=clock.ClockField()
    state,_,_,_=clock.packet_with_reference()
    _,before=clock.execute(counts,state,steps=4,space=space)
    original=source.contact
    # A source-word replacement that deletes ONLY the contact operation must
    # change the same complete history. No inferred rate or packet is refitted.
    monkeypatch.setattr(source,'contact',lambda duration,interaction:original(duration,0.))
    _,after=clock.execute(counts,state,steps=4,space=space)
    assert clock.cq_distance(before.values(),after.values())>.1


def test_uniform_bound_is_reference_and_volume_independent():
    for steps in (1,4,32):
        result=clock.error_budget(F(-1),steps)
        assert F(result['rate_error_halfdiamond_upper_on_confidence_event'])==F(9,256)
        expected=min(F(1),F(7,4)**3*F(5,4)**3/(3*steps**2))
        assert F(result['ideal_split_halfdiamond_upper'])==expected
        assert result['all_N2_states_and_arbitrary_untouched_reference']
    zero=clock.error_budget(0)
    assert F(zero['ideal_split_halfdiamond_upper'])==0
    assert F(zero['rate_error_halfdiamond_upper_on_confidence_event'])==0
    assert sum(F(5,4)**j/math.factorial(j) for j in range(4))>F(10,3)


def test_outside_rate_window_stop_preserves_quantum_input_and_reference():
    counts={edge:(n,n//2) for edge,(n,_) in retained_counts().items()}
    space=clock.ClockField()
    state=np.zeros((space.dimension,2),complex)
    state[0,0]=1/math.sqrt(2);state[-1,1]=1j/math.sqrt(2)
    decision,outputs=clock.execute(counts,state,space=space)
    assert decision=={'status':'stop','reason':'outside_predeclared_rate_window'}
    assert len(outputs)==1
    assert np.array_equal(next(iter(outputs.values())),state)
    assert next(iter(outputs.values())) is not state


def test_source_budget_bound_also_covers_scattering_input():
    counts=retained_counts();space=clock.ClockField()
    state=np.zeros((space.dimension,2),complex)
    state[space.index(0,7),0]=1/math.sqrt(2)
    state[space.index(1,126),1]=1j/math.sqrt(2)
    decision,outputs=clock.execute(counts,state,F(-1,2),steps=8,space=space)
    ideal=clock.pointers(space,clock.exact_state(space,state,math.log(25/8)/2))
    assert clock.cq_distance(outputs.values(),ideal)<float(F(
        clock.error_budget(F(-1,2),8)['complete_average_halfdiamond_upper']))
    assert decision['status']=='accept'


@pytest.mark.parametrize('value',[True,1.1,F(3,2),float('nan')])
def test_parameter_validation(value):
    with pytest.raises(ValueError):clock.rate_decision(retained_counts(),value)
    with pytest.raises(ValueError):clock.error_budget(value)


def test_odd_periodic_matching_and_wrong_reference_rejected():
    with pytest.raises(ValueError):clock.ClockField(5)
    with pytest.raises(ValueError):clock.reference_spectrum(beta=0.)
    with pytest.raises(ValueError):clock.reference_spectrum(side=6)
