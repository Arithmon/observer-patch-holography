"""Independent source, stopping, coherent and continuing-history controls."""
from fractions import Fraction as F
import json

import numpy as np
import pytest
import sympy as sp

from . import construction as c
from . import kernel as k
from ..observer import native
from ..primitive_observer import runtime


def independent_counts(bias=False):
    result = np.zeros((11, 66, 66), dtype=np.int64)
    for i, state in enumerate(k.quantum.STATES):
        for edge, seam in enumerate(k.quantum.mixing.seams()):
            for coin in (0, 1):
                output, _ = native.native(state, seam, coin)
                delta = output[k.CLOCK_PORT]-state[k.CLOCK_PORT]
                neighbour = next((v for v in seam if v != k.CLOCK_PORT), -1)
                record = (neighbour, delta) if delta else k.WAIT
                result[k.marks().index(record), k.quantum.INDEX[output], i] += 2 if bias and edge == 0 else 1
    return result


def test_actual_source_records_and_remote_non_tick_motion():
    counts = independent_counts()
    np.testing.assert_array_equal(counts, k.block_counts('C'))
    np.testing.assert_array_equal(counts, k.block_counts('E'))
    assert np.sum(counts[0]-np.diag(np.diag(counts[0]))) > 0
    np.testing.assert_array_equal(counts.sum(axis=(0,1)), np.full(66, 60))
    assert k.CLOCK_PORT == 2
    assert k.intervention_seam() == (0, 10)
    assert k.CLOCK_PORT not in k.intervention_seam()


def test_stopping_resolvent_against_independent_solve_and_wrong_order():
    counts = independent_counts().astype(float)/60
    inverse = np.linalg.inv(np.eye(66)-counts[0])
    moving = counts[1:].sum(axis=0)
    target = moving@inverse
    np.testing.assert_allclose(np.array(k.event_block('C'), float), target, atol=2e-14, rtol=0)
    assert np.max(np.abs(inverse@moving-target)) > .01
    np.testing.assert_allclose(target.sum(axis=0), 1, atol=1e-14, rtol=0)


def test_palm_is_actual_fixed_attempt_herald_not_first_hit_or_full_faithful_state():
    failed, success, palm = k.source_palm()
    assert k.trace(failed) == F(193, 198)
    assert k.trace(success) == F(5, 198)
    assert k.trace(palm) == 1 and len(palm) == 51
    assert k.event_apply('C', palm) == palm
    first = k.event_apply('C', k.diagonal((F(1, 66),)*66))
    assert sum(abs(first.get((i,i),F())-palm.get((i,i),F())) for i in range(66))/2 == F(1, 3)


def test_zero_immediate_hazard_does_not_mean_a_dead_clock():
    counts = independent_counts()
    hazard = counts[1:].sum(axis=(0,1))
    zero = np.flatnonzero(hazard == 0)
    assert len(zero) == 15
    matrix = k.event_block('C')
    for i in zero:
        assert sum(row[int(i)] for row in matrix) == 1
    tail = k.tail_certificate((1,2))
    assert F(tail[0]['maximum_timeout_probability']) == 1
    assert F(tail[1]['maximum_timeout_probability']) == F(1799, 1800)
    no = counts[0]
    assert max((no.T@no.T@np.ones(66, dtype=np.int64)).tolist()) == 3598


@pytest.mark.parametrize('law', ('C','E'))
def test_full_quantum_clock_channel_preserves_trace_with_reference(law):
    source = k.source(law).instrument
    pairs = k.quantum.phase_apparatus()[3]
    a,b = pairs[0]
    psi = np.zeros(132, complex); psi[2*a] = 1/np.sqrt(2); psi[2*b+1] = 1j/np.sqrt(2)
    rho = np.outer(psi,psi.conj())
    actual = runtime._apply(source, rho, 2)
    assert len(actual) == 11
    assert sum(np.trace(x).real for x in actual) == pytest.approx(1,abs=2e-14)
    assert all(np.linalg.eigvalsh(x).min() >= -2e-14 for x in actual)
    assert all(np.count_nonzero(x) <= 100 for x in actual)


@pytest.mark.parametrize('law', ('C','E'))
def test_exact_selected_coherence_blocks_match_direct_kraus_action(law):
    source = k.source(law).instrument
    phase = tuple(1 if i==0 else -1 if i==1 else 0 for i in range(12))
    for delta in (phase, tuple(-x for x in phase)):
        pairs = k.sectors()[delta]
        blocks = k.block_counts(law,delta)
        for column,(i,j) in enumerate(pairs):
            density=np.zeros((66,66),complex); density[i,j]=1
            actual=runtime._apply(source,density,1)
            for mark in range(11):
                target=np.zeros((66,66),complex)
                for row,key in enumerate(pairs): target[key]=blocks[mark][row][column]/60
                np.testing.assert_allclose(actual[mark],target,atol=2e-15,rtol=0)


def test_two_event_joint_laws_against_independent_dense_solve():
    counts=independent_counts().astype(float)/60
    inverse=np.linalg.inv(np.eye(66)-counts[0])
    event=[m@inverse for m in counts[1:]]
    uniform=np.full(66,1/66)
    success=counts[1:].sum(axis=0)@uniform
    palm=success/success.sum()
    actual=c.joint_states()
    for first in range(10):
        for second in range(10):
            expected=event[second]@event[first]@palm
            state=actual[first+1,-1,second+1]
            np.testing.assert_allclose([float(state.get((i,i),F())) for i in range(66)],expected,atol=2e-15,rtol=0)
    # Memory survives even after enforcing the unavoidable binary sign alternation.
    same=sum(k.trace(state) for (a,_,b),state in actual.items() if k.marks()[a][0]==k.marks()[b][0])
    assert same == F(5428819044974590805373,19780089058627344141475)
    assert same > F(1,5)


def test_event_triggered_intervention_is_native_and_changes_future_records():
    baseline=c.joint_states()
    changed=c.joint_states(intervention=True)
    assert len(baseline)==100 and len(changed)==300
    assert sum(map(k.trace,changed.values()),F())==1
    pooled={(a,b):sum(k.trace(changed[a,mark,b]) for mark in range(3))
            for a in range(1,11) for b in range(1,11)}
    distance=sum(abs(k.trace(state)-pooled[a,b]) for (a,_,b),state in baseline.items())/2
    assert distance>0
    for state in c.coherent_initial('C'):
        actual=runtime._apply(k.fixed_intervention('C'),k.dense(state),1)
        expected=k.exact_intervention('C',state)
        for a,e in zip(actual,expected,strict=True):
            np.testing.assert_allclose(a,k.dense(e),atol=2e-15,rtol=0)


def test_global_identity_holds_cancel_but_finite_cap_and_relative_weights_do_not():
    counts=independent_counts().astype(float)/60
    no,moving=counts[0],counts[1:].sum(axis=0)
    original=moving@np.linalg.inv(np.eye(66)-no)
    for activity in (.5,2/3):
        held=(1-activity)*np.eye(66)+activity*no
        np.testing.assert_allclose(activity*moving@np.linalg.inv(np.eye(66)-held),original,atol=4e-14,rtol=0)
    state=c.coherent_initial('C')[1]
    base=k.exact_step('C',state)
    lazy=k.exact_step('C',state,activity=F(1,2))
    assert sum(map(k.trace,lazy[1:])) == sum(map(k.trace,base[1:]))/2
    assert k.trace(lazy[0]) > k.trace(base[0])
    biased=c.joint_states(bias=True)
    assert any(k.trace(biased[key])!=k.trace(state) for key,state in c.joint_states().items())


@pytest.mark.parametrize('law,first,second', (
    ('C',F(3804737,7050540),F(3588793,7050540)),
    ('E',F(289535633,552465980),F(278113337,552465980))))
def test_phase_sensitive_prediction_on_actual_local_event_clock(law,first,second):
    for flag in (1,2):
        state=c.coherent_initial(law)[flag]
        for probability in (first,second):
            state=k.event_apply(law,state)
            actual=c.prepared._pointer_exact(c.prepared._balanced_exact(state),law)
            assert k.trace(state)==F(5,66)
            assert k.trace(actual[0 if flag==1 else 1])/k.trace(state)==probability


@pytest.mark.parametrize('intervention', (False,True))
def test_complete_finite_joint_execution_and_all_timeouts(intervention):
    report=c.finite_joint(intervention,cap=1)
    rows=report['complete_records']
    assert len(rows)==2+10*(3 if intervention else 1)*11
    assert sum(F(row['probability']) for row in rows)==1
    assert any(F(row['probability'])==0 for row in rows)
    assert rows[0]['record']==[0,0,0,-1,0,0]
    assert F(rows[0]['probability'])==F(193,198)
    assert rows[-1]['record']==[1,2,0,-1,0,0]
    assert F(rows[-1]['probability'])>0
    assert report['full_state_maximum_frobenius_residual_diagnostic']<2e-14
    assert report['reset_between_local_clock_events'] is False


@pytest.mark.parametrize('law', ('C','E'))
def test_complete_phase_execution_includes_native_herald_failure_and_tail(law):
    result=c.coherent_execution(law,cap=1)
    assert len(result['complete_records'])==43
    assert sum(F(row['probability']) for row in result['complete_records'])==1
    assert F(result['complete_records'][0]['probability'])==F(28,33)
    assert result['maximum_frobenius_residual_diagnostic']<2e-14


def test_serialized_continuation_replays_same_source_clock_history():
    instrument=k.source('C').instrument
    program=runtime.Program('two_local_clock_reads',66,np.eye(66)/66,
        tuple(runtime.Event(f'local_source_{i}',instruments=(instrument,)) for i in range(2)))
    whole=runtime.execute(program)
    prefix=runtime.execute(program,stop_after=1)
    checkpoint=runtime.checkpoint_from_json(program,runtime.checkpoint_to_json(prefix))
    resumed=runtime.execute(program,checkpoint=checkpoint)
    assert runtime.complete_distance(whole,resumed)==0
    assert len(whole.branches)==121


def test_public_exact_results_are_json_and_source_bindings_exist():
    assert all((c.ROOT/path).is_file() for path in c.SOURCE_PATHS)
    report=c.palm_report()
    assert json.loads(json.dumps(report))==report
    assert report['Palm_is_faithful_on_full_M66'] is False
    with pytest.raises(TypeError): k.sectors()[k.ZERO]=()


@pytest.mark.parametrize('law', ('bad',1,True,None))
def test_invalid_candidate_rejected(law):
    with pytest.raises(ValueError): k.source(law)


@pytest.mark.parametrize('activity', (0,F(0),F(-1),F(2),.5,True))
def test_invalid_activity_rejected(activity):
    with pytest.raises(ValueError): k.exact_step('C',{},activity=activity)


@pytest.mark.parametrize('value', (0,-1,True,1.0,17))
def test_invalid_coherent_caps_rejected(value):
    with pytest.raises(ValueError): c.coherent_execution('C',value)


def test_invalid_state_sector_mark_and_protocol_rejected():
    with pytest.raises(ValueError): k.sparse({(0,0):.5})
    with pytest.raises(ValueError): k.sparse({(-1,0):F(1)})
    with pytest.raises(ValueError): k.block_counts('C',(1,))
    with pytest.raises(ValueError): k.event_block('C',mark=0)
    with pytest.raises(ValueError): k.event_block('C',mark=True)
    with pytest.raises(ValueError): c.joint_states(intervention=1)
    with pytest.raises(ValueError): c.finite_joint(False,cap=True)


def test_faithful_local_palm_corner_and_modular_covariance_obstruction():
    result=c.local_modular_corner()
    assert result['corner_dimension']==51
    assert result['faithful_corner_reference_weight_census']=={
        '1/100':30,'1/50':10,'1/25':5,'1/20':6}
    assert result['completed_event_channel_maps_entire_native_algebra_into_corner']
    assert result['no_tick_map_can_leave_corner']
    assert result['stopped_channel_commutes_with_corner_modular_flow'] is False
    for row,coefficient in zip(result['covariance_witnesses'],(F(1,14),F(1,24)),strict=True):
        assert F(row['eventual_channel_coefficient'])==coefficient
        assert F(row['input_modular_weight_ratio'])==F(5,4)
        assert F(row['output_modular_weight_ratio'])==F(1,2)
        unit=np.zeros((66,66),complex)
        unit[tuple(row['input_matrix_unit'])]=1
        direct=runtime._apply(k.source(row['candidate']).instrument,unit,1)
        output=direct[k.marks().index((4,-1))][tuple(row['output_matrix_unit'])]
        assert output==pytest.approx(1/60,abs=2e-15)


def test_cached_exact_apis_reject_bool_aliases_and_empty_state_bad_controls():
    k.block_counts('C'); k.event_block('C',mark=1)
    for function in (k.block_counts,k.resolvent,k.event_block):
        with pytest.raises(ValueError): function('C',bias=0)
    with pytest.raises(ValueError): k.event_block('C',mark=True)
    with pytest.raises(ValueError): k.event_apply('C',{},mark=True)
    with pytest.raises(ValueError): k.event_apply('C',{},bias=0)
