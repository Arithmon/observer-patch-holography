"""Exact fail-closed controls and independent pointer/Kraus matrix checks."""
import cmath
from copy import deepcopy
from fractions import Fraction as F
import importlib.util
import itertools
import json
import math
from pathlib import Path

import pytest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('sequential_instrument_verifier',HERE/'verify_sequential_instrument.py')
verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)


@pytest.fixture(scope='module')
def packet(): return verifier.load()


def test_complete_independent_mathematical_replay(packet):
    result=verifier.verify(packet)
    assert result['mathematical_replay'] is True and result['parent_events_replayed']==5888
    assert result['uniformly_resolved_steps']==[5,6,7,8,9,10,11,13,14,15,16,17,18,20,21]
    assert result['readout_times']==21 and result['prior_pair_brackets']==210


def test_internal_arithmetic_does_not_claim_parent_replay(packet):
    result=verifier.verify_arithmetic(packet)
    assert result['mathematical_replay'] is False and result['parent_events_replayed'] is None


def test_positive_sequential_signal_and_measured_reference(packet):
    row=packet['rows'][15]
    assert F(row['sequential_split_response_abs_lower'])>F('0.00950581279')
    assert F(row['sequential_probability_error_upper'])<F('0.00083418122')
    assert F(row['continuous_sequential_response_abs_lower'])>F('0.00867163157')
    assert 1-F(packet['rows'][-1]['discrete_attenuation_interval'][0])<F('0.000044641')
    assert packet['scope']['branch_conditioned_clock_reconstructed'] is False
    assert 'same centered instrument' in packet['scope']['reference']


def test_full_action_first_two_commutators(packet):
    v=verifier.current_module('code/source_scalar_execution/verify_source_scalar_execution.py','_test_instrument_scalar')
    mass,action,_,g,_=v.model(); R=v.R
    g0=sum((m*x*x for m,x in zip(mass,g)),R())
    ag=[sum((a*g[j] for j,a in row),R()) for row in action]
    g1=sum((m*x*y for m,x,y in zip(mass,g,ag)),R())
    assert packet['cross_time_brackets'][0]['commutator_over_i_Qphi']==(v.TAU*g0).encode()
    assert packet['cross_time_brackets'][1]['commutator_over_i_Qphi']==(2*v.TAU*g0-v.TAU**3*g1).encode()
    assert packet['energy_injected_per_readout_Qphi']==(g0/8).encode()


def test_continuous_bound_uses_worst_opposite_clock_endpoints(packet):
    g=F(packet['detector_mass_norm_squared_upper'])
    for j,row in enumerate(packet['rows']):
        b=F(row['recovered_elapsed_time_interval'][1])
        worst=sum((g*g*(b-F(old['recovered_elapsed_time_interval'][0]))**2/8 for old in packet['rows'][:j]),F(0))
        loss=1-F(row['continuous_attenuation_interval'][0])
        assert worst<=loss<worst+F(1,10**15)
        assert F(row['attenuation_difference_upper'])==max(loss,1-F(row['discrete_attenuation_interval'][0]))


@pytest.mark.parametrize('mutation',[
    'commutator_sign','missing_half_factor','discard_time','discard_pair','zero_backaction',
    'continuous_equals_discrete','drop_clock_error','flip_sine_sign','zero_energy_injection',
    'fake_observed_outcome','branch_clock','independent_trials','no_born_input',
    'wrong_pointer_basis','wrong_local_coefficient','nonedge_cnot','duplicate_cnot',
    'false_operation_count','boolean_step','extra_top_field',
])
def test_mutation_rejected(packet,mutation):
    bad=deepcopy(packet)
    if mutation=='commutator_sign':
        bad['cross_time_brackets'][0]['commutator_over_i_Qphi']=[str(-F(x)) for x in bad['cross_time_brackets'][0]['commutator_over_i_Qphi']]
    elif mutation=='missing_half_factor': bad['cross_time_brackets'][0]['square_upper']=str(4*F(bad['cross_time_brackets'][0]['square_upper']))
    elif mutation=='discard_time': bad['rows'].pop()
    elif mutation=='discard_pair': bad['cross_time_brackets'].pop()
    elif mutation=='zero_backaction': bad['rows'][15]['discrete_attenuation_interval']=['1','1']
    elif mutation=='continuous_equals_discrete': bad['rows'][15]['continuous_attenuation_interval']=bad['rows'][15]['discrete_attenuation_interval']
    elif mutation=='drop_clock_error': bad['rows'][15]['sequential_probability_error_upper']='0'
    elif mutation=='flip_sine_sign': bad['instrument']['unitary']='exp(i Z tensor Phi(g)/2)'
    elif mutation=='zero_energy_injection': bad['energy_injected_per_readout_Qphi']=['0','0']
    elif mutation=='fake_observed_outcome': bad['scope']['physical_clock_or_observed_outcomes']=True
    elif mutation=='branch_clock': bad['scope']['branch_conditioned_clock_reconstructed']=True
    elif mutation=='independent_trials': bad['instrument']['all_record_marginals_independent']=True
    elif mutation=='no_born_input': bad['scope']['born_instrument_probability_law_supplied']=False
    elif mutation=='wrong_pointer_basis': bad['regional_gadget']['local_pointer_readout_bases'][0][1]='X'
    elif mutation=='wrong_local_coefficient': bad['regional_gadget']['local_controlled_field_coefficients_Qphi'][0][1]=['0','0']
    elif mutation=='nonedge_cnot': bad['regional_gadget']['ghz_cnot_tree'][0]=[32,63]
    elif mutation=='duplicate_cnot': bad['regional_gadget']['ghz_cnot_tree'][1]=bad['regional_gadget']['ghz_cnot_tree'][0]
    elif mutation=='false_operation_count': bad['regional_gadget']['total_ideal_readout_operations']=21
    elif mutation=='boolean_step': bad['rows'][0]['step']=True
    elif mutation=='extra_top_field': bad['sampled_history']=[0]*21
    with pytest.raises(ValueError): verifier.verify_arithmetic(bad)


@pytest.mark.parametrize('data',['{"a":1,"a":2}','{"a":0.1}','{"a":NaN}'])
def test_noncanonical_json_rejected(data,tmp_path):
    path=tmp_path/'bad.json';path.write_text(data)
    with pytest.raises(ValueError): verifier.load(path)


@pytest.mark.parametrize('n',[1,2,3,5])
def test_direct_ghz_fine_record_amplitudes_and_sine_sign(n):
    # Contract the two nonzero GHZ amplitudes with each actual pointer
    # measurement bra. This shares no algebra with the certificate verifier.
    phi=[(i+1)/7 for i in range(n)]; total=sum(phi); probabilities=[0.,0.]
    for bits in itertools.product((0,1),repeat=n):
        a=cmath.exp(-0.5j*total)/math.sqrt(2)
        b=cmath.exp(0.5j*total)/math.sqrt(2)
        for i,bit in enumerate(bits):
            a/=math.sqrt(2)
            b*=((-1j if i==0 else 1)*(-1)**bit)/math.sqrt(2)
        parity=sum(bits)%2
        k=(cmath.exp(-0.5j*total)+(-1j if parity==0 else 1j)*cmath.exp(0.5j*total))/2
        assert abs((a+b)-k/2**((n-1)/2))<2e-15
        probabilities[parity]+=abs(a+b)**2
    assert abs(probabilities[0]-(1+math.sin(total))/2)<2e-15
    # The explicitly chosen U sign and +Y label must give certainty at pi/2.
    k0=(cmath.exp(-0.25j*math.pi)-1j*cmath.exp(0.25j*math.pi))/2
    assert abs(abs(k0)**2-1)<2e-15


def mm(a,b):
    n=len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def adj(a):return [list(map(complex.conjugate,row)) for row in zip(*a)]
def scale(a,c):return [[c*x for x in row] for row in a]
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def tr(a):return sum(a[i][i] for i in range(len(a)))
def sandwich(u,r):return mm(mm(u,r),adj(u))


def test_sequential_kraus_tree_matches_characteristic_product():
    # Finite clock/shift Weyl matrices: independent non-Gaussian mixed-state
    # and joint-record test of the centered channel's phase multiplier.
    n=16; root=cmath.exp(2j*math.pi/n)
    def unitary(a,b):return [[root**(b*i) if i==(j+a)%n else 0j for j in range(n)] for i in range(n)]
    ops=[unitary(0,1),unitary(1,0),unitary(1,1),unitary(2,1)]
    psi=[complex(i%3+1,(i+1)%4) for i in range(n)]; norm=math.sqrt(sum(abs(x)**2 for x in psi));psi=[x/norm for x in psi]
    initial=[[x*y.conjugate() for y in psi] for x in psi]
    branches=[('',initial)]; previous=[]
    for u in ops:
        w=mm(u,u); characteristic=tr(mm(initial,w)); attenuation=1.
        for old in previous:
            conjugated=mm(mm(adj(old),w),old)
            phase=tr(mm(conjugated,adj(w)))/n
            assert max(abs(conjugated[i][j]-phase*w[i][j]) for i in range(n) for j in range(n))<5e-14
            attenuation*=phase.real
        k=[scale(add(adj(u),scale(u,-1j)),.5),scale(add(adj(u),scale(u,1j)),.5)]
        branches=[(history+str(bit),sandwich(k[bit],rho)) for history,rho in branches for bit in (0,1)]
        actual=sum(tr(rho).real for history,rho in branches if history[-1]=='0')
        assert abs(actual-(1+attenuation*characteristic.imag)/2)<2e-13
        assert all(tr(rho).real>=-2e-14 for _,rho in branches)
        assert abs(sum(tr(rho).real for _,rho in branches)-1)<2e-13
        previous.append(u)


def test_full_verifier_rejects_missing_parent_replay(packet,monkeypatch):
    original=verifier.current_module
    class BadParent:
        load=staticmethod(lambda path:{})
        verify=staticmethod(lambda *args,**kwargs:{'full_clock_replay':False})
    monkeypatch.setattr(verifier,'current_module',lambda path,name:BadParent if path==verifier.CLOCK_VERIFIER else original(path,name))
    with pytest.raises(ValueError,match='full parent clock and quantum replay'): verifier.verify(packet)
