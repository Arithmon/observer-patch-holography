"""Adversarial tests of the nonlinear detector consumer."""
from copy import deepcopy
from fractions import Fraction
import pytest
import readout
import verify_readout as verifier
from rational import I, phi


@pytest.fixture(scope='module')
def packet():
    return verifier.load()


def test_reproduction(packet):
    assert packet == readout.bounds()
    assert verifier.verify(packet)['mathematical_replay'] is True


@pytest.mark.parametrize('field', ['quantum_dynamics_or_born_outcomes',
    'source_action_population_clock_selected','whole_state_numerical_trajectory_enclosed',
    'first_kick_is_exact_time_sample','full_fermion_generation_executed','laboratory_or_empirical_result'])
def test_false_scope_rejected(packet,field):
    p=deepcopy(packet); p['scope'][field]=True
    with pytest.raises(ValueError,match='scope'):
        verifier.verify(p,replay_parent=False)


@pytest.mark.parametrize('mutation', ['drop','duplicate','value','error','writer_event','time','zero_majorant','signal','window'])
def test_false_green_mutations(packet,mutation):
    p=deepcopy(packet)
    if mutation=='drop': p['checkpoints'].pop()
    if mutation=='duplicate': p['checkpoints'][-1]=p['checkpoints'][0]
    if mutation=='value': p['checkpoints'][-1]['decoded_intensity']='1/25'
    if mutation=='error': p['checkpoints'][-1]['error_upper']='0'
    if mutation=='writer_event': p['checkpoints'][-1]['event']=526
    if mutation=='time': p['checkpoints'][-1]['time']='1/50'
    if mutation=='zero_majorant': p['majorants']['scalar_acceleration']='0'
    if mutation=='signal': p['neighbor_response']['magnitude_lower_coefficient']='1'
    if mutation=='window': p['window']=['0','1']
    with pytest.raises(ValueError):
        verifier.verify(p,replay_parent=False)


def test_intermediate_split_is_excluded(packet):
    assert len(packet['checkpoints'])==9
    assert all(r['event'] not in [525,526] for r in packet['checkpoints'])


def test_pasted_frame_cannot_replace_readback():
    parent=readout.load_parent(); run=deepcopy(parent['runs'][1]); cp=run['checkpoints'][-1]
    expected=verifier.readback(run,cp)
    cp['values']['psi/37']=[99.,99.]
    assert verifier.readback(run,cp)==expected
    event=run['events'][cp['event']]
    read=next(r for r in event['reads'] if r['port']=='psi/37')
    read['writer']=0
    with pytest.raises(ValueError,match='writer'):
        verifier.readback(run,cp)


def test_exact_interval_arithmetic_against_rational_samples():
    for a,b,c,d in [(-3,2,-5,-1),(-2,3,-4,5),(1,4,2,7),(-4,-2,1,3)]:
        x,y=I(a,b),I(c,d)
        for u in [Fraction(a),Fraction(a+b,2),Fraction(b)]:
            for v in [Fraction(c),Fraction(c+d,2),Fraction(d)]:
                assert (x*y).lo<=u*v<=(x*y).hi
                if not c<=0<=d:
                    assert (x/y).lo<=u/v<=(x/y).hi
            assert (x**2).lo<=u*u<=(x**2).hi
    f=phi()
    assert f.lo*f.lo-f.lo-1<0<f.hi*f.hi-f.hi-1
    with pytest.raises(ValueError): I(1)/I(-1,1)


def test_current_and_detector_are_resolved(packet):
    assert Fraction(packet['current']['lower_on_window'])>Fraction(15,1000)
    assert Fraction(packet['electric_feedback']['lower_at_end'])>Fraction(15,100000)
    assert Fraction(packet['neighbor_response']['magnitude_lower_at_end'])==Fraction(45,1000000)
    assert max(Fraction(r['error_upper']) for r in packet['checkpoints'])<Fraction(6,1000000)
