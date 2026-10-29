from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import clock as m


def counts():
    n = m.inference.TRIALS_PER_SEAM
    w = (25,25)+(8,)*10
    return {(a,b):(n,round(F(n*w[a],w[a]+w[b]))) for a,b in m.native.SEAMS}


def weights():
    return tuple(F(v) for v in (25,25)+(8,)*10)


@pytest.fixture(scope='module')
def receipt():
    return m.report()


def test_record_inference_uses_no_truth_table_or_reference(monkeypatch):
    observed = counts()
    def forbidden(*args, **kwargs):
        raise AssertionError('private source weight/reference access')
    monkeypatch.setattr(m.native_clock, 'parameters', forbidden)
    monkeypatch.setattr(m.native, 'reference', forbidden)
    result = m.infer_clock(observed)
    assert result['status'] == 'accept'
    assert result['tree_distance'] == 2
    assert F(result['beta_times_U_error_on_confidence_event']) == F(3,80)
    assert F(result['clock_halfdiamond_error_per_unit_parameter']) == F(3,160)
    lo, hi = map(F, result['energy_ratio_interval'])
    assert lo <= F(25,8) <= hi
    assert abs(result['beta_times_U_estimate_diagnostic']-math.log(25/8)) < .001


def test_exact_Gibbs_identity_on_the_declared_M2_and_not_full_field():
    rho = np.diag([25/33,8/33])
    h = np.diag([0.,1.]); beta = math.log(25/8)
    assert np.linalg.norm(np.diag(-np.log(np.diag(rho)))-beta*h-math.log(33/25)*np.eye(2)) < 1e-15
    for s in (-1.,.3,1.):
        modular = np.diag(np.diag(rho)**(-1j*s))
        model = np.diag(np.exp(1j*beta*s*np.diag(h)))
        assert np.linalg.norm(modular/modular[0,0]-model) < 5e-16
    # Restoring any unused field direction makes the encoded density nonfaithful.
    assert np.linalg.matrix_rank(np.diag([25/33,8/33,0])) == 2


def test_clock_isometry_uses_the_actual_interacting_two_fermion_fiber():
    from ..bound_pair.fiber import PairFiber
    for side in (4,6):
        attachment=m.clock_subspace(side=side)
        f=PairFiber(attachment['momentum'],side=side)
        r,d=attachment['basis'].T
        assert np.linalg.norm(f.apply(r))<1e-13
        assert np.linalg.norm(f.apply(d)-d)<1e-13
        rho=(25*np.outer(r,r.conj())+8*np.outer(d,d.conj()))/33
        assert np.linalg.norm(attachment['basis'].conj().T@rho@attachment['basis']-np.diag([25/33,8/33]))<1e-13
        assert max(attachment['residuals'].values())<1e-13


def test_pointer_is_complete_on_all_addresses_and_preserves_corner_coherence():
    (outside, inside), info = m.pointer()
    mask = np.zeros(512); mask[list(m._corner_codes())] = 1
    assert np.linalg.norm(inside-np.diag(mask),'fro') < float(m.POINTER_GATE)
    assert np.linalg.norm(outside-np.diag(1-mask),'fro') < float(m.POINTER_GATE)
    assert np.linalg.norm(inside.conj().T@inside+outside.conj().T@outside-np.eye(512),'fro') < 3e-11
    assert info['full_input_columns'] == 512
    assert info['source_pulses_per_herald'] == 2068
    # Complex entanglement inside the corner survives without resolving energy.
    v = np.zeros((512,2),complex); a,b=m._corner_codes()
    v[a,0]=1/math.sqrt(2); v[b,1]=1j/math.sqrt(2)
    assert np.linalg.norm(inside@v-v) < 1e-11
    assert np.linalg.norm(outside@v) < 1e-11
    inside[:] = 0
    assert np.linalg.norm(m.pointer()[0][1]) > 1


def test_warm_pointer_cache_cannot_hide_a_changed_literal_program(monkeypatch):
    m.pointer()
    original = m.native_routing._canonical_word
    def wrong(*args, **kwargs):
        word = original(*args, **kwargs)
        word[0] = dict(word[0], angle=word[0]['angle']+.2)
        return word
    monkeypatch.setattr(m.native_routing, '_canonical_word', wrong)
    with pytest.raises(ArithmeticError, match='macro'):
        m.pointer()


@pytest.mark.parametrize('law',['C','E'])
def test_complete_map_and_readout_from_independent_raw_basis_evolution(law):
    rho = np.zeros((66,66),complex)
    a,b=m.corner_indices(); rho[a,a]=25/33;rho[b,b]=8/33
    out=m.execute(rho,weights(),law,raw=False)
    alpha=F(1519,1980) if law=='C' else F(37,60)
    expected=85/198+float(F(17,66)*alpha)*math.cos(math.log(25/8))
    assert m.record_probability(out,law)==pytest.approx(expected,abs=1e-13)
    assert sum(np.trace(v).real for v in out.values())==pytest.approx(1,abs=1e-13)
    assert sum(np.trace(v).real for k,v in out.items() if k[0]=='outside')==pytest.approx(14/99,abs=1e-13)
    assert m.record_probability(out,law)+m.record_probability(out,law,1)==pytest.approx(85/99,abs=1e-13)
    zero=m.raw_output(out,law,m.native.SEAMS[0],0,3,0,3)
    assert not np.any(zero)
    live=next(k for k in out if k[0]=='read')
    removed=dict(out); del removed[live]
    with pytest.raises(ValueError,match='missing live'):
        m.raw_output(removed,law,*live[1:])


def test_complete_instrument_preserves_untouched_complex_reference_without_conditioning():
    a,b=m.corner_indices()
    vector=np.zeros((66,2),complex)
    vector[a,0]=.5;vector[b,1]=.5j;vector[0,0]=.5j;vector[0,1]=.5
    density=np.outer(vector.ravel(),vector.ravel().conj())
    out=m.execute(density,weights(),'E',raw=False)
    marginal=np.zeros((2,2),complex)
    for block in out.values():
        if len(block)==2:
            marginal+=block
        else:
            marginal+=np.einsum('irjs->rs',block.reshape(66,2,66,2)*np.eye(66)[:,None,:,None])
    expected=vector.T@vector.conj()
    assert np.linalg.norm(marginal-expected)<1e-13
    # Rejected source states and their nontrivial reference correlations survive.
    outside=out['initial_outside',]
    assert np.trace(outside).real==pytest.approx(.5)
    assert abs(outside[0,1])>.2


def test_stopped_training_retains_carrier_and_does_not_prepare(monkeypatch):
    observed=counts(); edge=next(iter(observed));observed[edge]=(1<<20,0)
    def forbidden(*args,**kwargs):
        raise AssertionError('stopped protocol executed a later service')
    monkeypatch.setattr(m.record_clock,'continued_source_density',forbidden)
    monkeypatch.setattr(m,'source_budget',forbidden)
    decision,outputs=m.continue_from_records(observed,17,'C')
    assert decision['status']=='stop'
    state=next(iter(outputs.values()))
    assert state[17,17]==1 and np.trace(state)==1


def test_compiler_failure_is_a_retained_stop(monkeypatch):
    def rejected(*args,**kwargs):
        raise ArithmeticError('wrong source unitary')
    monkeypatch.setattr(m,'source_budget',rejected)
    decision,outputs=m.continue_from_records(counts(),17,'E')
    assert decision=={'status':'stop','reason':'source_certification_failed'}
    assert next(iter(outputs.values()))[17,17]==1


def test_full_retained_native_tape_and_source_history(receipt):
    assert receipt['training_attempts']==31457280
    assert F(receipt['complete_average_half_trace_upper'])<F(19,1000)
    assert not receipt['error_bound_conditioned_on_herald_or_training_acceptance']
    for law,row in receipt['complete_readouts'].items():
        assert F(row['complete_output_mass_error_upper'])<F(1,1000000)
        assert F(row['initial_outside_probability_interval'][0])==F(24,25)
        assert F(row['energy_r_record_probability_interval'][0])>0
        assert row['exact_ideal_response_and_complete_quantum_outputs_checked']
    assert receipt['source_sha256']==m.bindings()
    assert not receipt['scope']['physical_temperature_time_energy_normalization']
