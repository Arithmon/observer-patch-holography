"""Complete stopped phase/work laws, independent targets and adverse controls."""
from dataclasses import replace
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from ..clocked_common_observer import construction as work, experiment as old_budget
from ..common_source import protocol, scalar_coupling
from ..native_relational_readout import construction as relational
from ..primitive_observer import runtime
from . import construction as c, experiment as e


@pytest.fixture(scope='module')
def small():
    w = work.execute_work(2)
    return w, old_budget.work_certificates(w)


@pytest.mark.parametrize('law', ('C', 'E'))
@pytest.mark.parametrize('cap', (1, 4, 16))
def test_exact_phase_law_reuses_previous_first_current_result_without_postselection(law, cap):
    outputs = c.exact_outputs(law, cap)
    assert tuple(outputs) == c.records(cap)
    assert sum(map(c.clock._trace, outputs.values())) == 1
    assert c.clock._trace(outputs[0,0,0,-1]) == F(28,33)
    stats = e.exact_native_statistics(law,cap)
    old = relational.exact_stopped_state(relational.source(law),through=cap)
    assert stats['phase_match_and_success'] == c.YIELD*relational.phase_probabilities(old)[0]
    assert stats['timeout'] == c.YIELD*c.clock.statistics(law,cap)['timeout_probability']
    assert stats['phase_match_and_success'] < c.YIELD
    for n in range(1,cap+1):
        for j in (-1,1):
            for b in (0,1):
                assert c.clock._trace(outputs[1,n,j,b]) == c.clock._trace(outputs[2,n,j,b])


def test_phase_discrimination_is_absent_at_first_slot_then_positive_with_identical_work_control_law():
    for cap in (1,4,16):
        left,right=[e.exact_native_statistics(law,cap) for law in ('C','E')]
        lp,rp=left['probabilities'],right['probabilities']
        tv=sum(abs(lp[r]-rp[r]) for r in lp)/2
        assert tv == left['phase_match_and_success']-right['phase_match_and_success']
        assert (tv == 0) if cap == 1 else (tv > 0)
        for h in (1,2):
            for n in range(1,cap+1):
                for j in (-1,1):
                    assert sum(lp[h,n,j,b] for b in (0,1)) == sum(rp[h,n,j,b] for b in (0,1))


@pytest.mark.parametrize('law', ('C','E'))
def test_original_numerical_densities_match_exact_complete_native_instrument(law):
    native=c.execute_native(law,4)
    certificate=c.certify_native(native)
    assert 0 < certificate['error'] < F(1,10**12)
    assert len(native.outputs)==35
    assert all(not a.flags.writeable for _,a in native.outputs)
    assert sum(c.clock._trace(s) for s in c.exact_outputs(law,4).values())==1
    # The initial execution really reads the fixed-seam current on I/66.
    prep=native.preparation
    assert prep.records==((0,),(1,),(2,))
    assert [round(b.probability,12) for b in prep.branches] == [round(float(x),12) for x in (F(28,33),F(5,66),F(5,66))]
    resumed=runtime.execute(prep.program,checkpoint=runtime.checkpoint_from_json(
        prep.program,runtime.checkpoint_to_json(runtime.execute(prep.program,stop_after=1))))
    assert runtime.complete_distance(prep,resumed)==0


def test_actual_native_execution_never_substitutes_exact_target_or_normalized_phase_state(monkeypatch):
    def fail(*args,**kwargs):
        raise AssertionError('an exact target or supplied coherent preparation was consumed')
    monkeypatch.setattr(c,'exact_outputs',fail)
    monkeypatch.setattr(relational,'exact_preparation',fail)
    monkeypatch.setattr(c.clock.quantum,'phase_apparatus',fail)
    native=c.execute_native('C',2)
    assert abs(sum(np.trace(a) for _,a in native.outputs)-1)<1e-13


def test_cap_extension_keeps_stopped_outputs_but_preserves_unmeasured_timeout():
    small=dict(c.execute_native('C',2).outputs)
    large=dict(c.execute_native('C',4).outputs)
    for r,a in small.items():
        if r[2] or r[0]==0:
            np.testing.assert_array_equal(a,large[r])
    for h in (1,2):
        assert not np.array_equal(small[h,2,0,-1],large[h,4,0,-1])
        assert (h,2,0,0) not in small
        # An unauthorized phase analyzer would disturb the retained tail.
        b=c.interface.mixer()
        assert np.linalg.norm(b@small[h,2,0,-1]@b.conj().T-small[h,2,0,-1])>1e-3


def test_joint_complete_outputs_include_failure_work_vacuum_and_all_stopped_records(small):
    w,bounds=small
    native=c.execute_native('C',2)
    rows=list(c.output_rows(native,w))
    assert len(rows)==259
    assert rows[0][0]==(0,0,0,(),-1)
    np.testing.assert_array_equal(rows[0][2],w.prefixes[0].branches[0].density)
    assert rows[-1][0]==(2,2,0,(),-1)
    assert abs(sum(np.trace(nr)*np.trace(wr) for _,nr,wr in rows)-1)<1e-12
    obs=e.observations(native,w)
    assert 0<float(obs['phase_match_and_final_gauge_one'])<float(obs['phase_match_and_success'])
    assert float(obs['phase_match_and_success'])<float(c.YIELD)
    proof=e.joint_bound(native,w,bounds)
    assert 0<proof['finite_error']<F(1,10**5)


def test_joint_measurement_matches_independent_continuous_work_evolution(small):
    w,bounds=small
    initial=np.zeros(16,complex); initial[0]=1
    unitary=expm(-.25j*scalar_coupling.ideal_hamiltonian(2,.25))
    predictions=[]
    for law in ('C','E'):
        exact=c.exact_outputs(law,2)
        wanted=0.
        for n in (1,2):
            for j in (-1,1):
                history=protocol.execute(protocol.continuation_program(2,128,native_current=j),
                    backend='continuous',initial=np.linalg.matrix_power(unitary,n)@initial)
                gauge=sum(float(np.linalg.norm(branch.state)**2) for branch in history.branches if branch.record[3])
                wanted+=gauge*float(sum(c.clock._trace(exact[h,n,j,0]) for h in (1,2)))
        native=c.execute_native(law,2)
        observed=e.observations(native,w)
        bound=e.joint_bound(native,w,bounds)
        assert abs(float(observed['phase_match_and_final_gauge_one'])-wanted)<=2*float(bound['finite_error'])
        predictions.append(observed)
    assert abs(float(predictions[0]['final_gauge_one_and_success']-predictions[1]['final_gauge_one_and_success']))<1e-13
    assert predictions[0]['phase_match_and_final_gauge_one']>predictions[1]['phase_match_and_final_gauge_one']


def test_resource_counts_charge_failed_herald_and_one_preparation(small):
    w,_=small
    resources=e.resources(w,2)
    native=c.clock.statistics('C',2)
    assert resources['complete_joint_records']==259
    assert resources['maximum_native_attempts_including_herald_and_pointer']==4
    assert resources['failed_herald_cost']['work_preparations']==1
    assert resources['failed_herald_cost']['balanced_mixer_calls']==0
    assert resources['timeout_after_n_slots']['phase_reads']==0
    assert F(resources['expected_native_attempts_exact'])==1+c.YIELD*(native['expected_native_attempts']+1-native['timeout_probability'])


def test_bad_work_binding_wrong_states_and_postselection_fail_closed(small):
    w,bounds=small
    native=c.execute_native('C',2)
    different=work.execute_work(2,delta=F(1,2))
    with pytest.raises(ValueError):
        e.joint_bound(native,different,bounds)
    outputs=list(native.outputs)
    record,a=outputs[1]
    changed=a.copy(); changed[0,1]+=1j/16
    outputs[1]=(record,changed)
    with pytest.raises(ArithmeticError):
        e.joint_bound(replace(native,outputs=tuple(outputs)),w,bounds)
    selected=tuple((r,np.zeros_like(a) if r[0]==0 else a/float(c.YIELD)) for r,a in native.outputs)
    with pytest.raises(ArithmeticError):
        e.joint_bound(replace(native,outputs=selected),w,bounds)


@pytest.mark.parametrize('change',('missing','duplicate','boolean','masked','shape'))
def test_complete_original_native_array_and_record_contract(change):
    native=c.execute_native('C',1)
    outputs=list(native.outputs)
    r,a=outputs[0]
    if change=='missing': outputs.pop()
    elif change=='duplicate': outputs[-1]=outputs[0]
    elif change=='boolean': outputs[0]=((False,*r[1:]),a)
    elif change=='masked': outputs[0]=(r,np.ma.array(a))
    else: outputs[0]=(r,np.zeros((1,1)))
    with pytest.raises(ValueError):
        c.certify_native(replace(native,outputs=tuple(outputs)))


def test_original_uint64_and_subnormal_components_are_not_coerced_away():
    native=c.execute_native('C',1)
    outputs=list(native.outputs)
    r,_=outputs[0]
    integer=np.zeros((66,66),dtype=np.uint64); integer[0,0]=2**63+1
    outputs[0]=(r,integer)
    detached=c.validate_native(replace(native,outputs=tuple(outputs)))
    assert detached.outputs[0][1].dtype==np.uint64
    assert int(detached.outputs[0][1][0,0])==2**63+1
    with pytest.raises(ValueError):
        c._balanced_exact({(c.interface.phase_pairs()[0][0],65):F(1)})


def test_full_source_implementation_bound_retains_padding_and_norm_gain():
    from ..compiled_native_interface import construction as compiled
    error,metadata=e.compiled_service()
    certificate=compiled.certify()
    assert error==certificate['half_diamond_upper']*(1+certificate['channel_norm_upper'])
    assert error>2*certificate['half_diamond_upper']
    assert error<F(1,10**11)
    assert metadata['full_source_carrier_dimension']==512
    assert metadata['source_matrix_compressed_to_native_sector'] is False
    assert metadata['floating_full_512D_history_executed_or_certified'] is False
