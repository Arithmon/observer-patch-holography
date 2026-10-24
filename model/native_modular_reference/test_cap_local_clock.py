"""Bare-cap algebra, spectator coherence, full source and complete history."""
from copy import deepcopy
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import cap_local_clock as c


def _random_density(seed=20261024, ref=2):
    rng = np.random.default_rng(seed)
    vector = rng.normal(size=66*ref)+1j*rng.normal(size=66*ref)
    vector /= np.linalg.norm(vector)
    return np.outer(vector, vector.conj())


def _outside(density, cap):
    ref = len(density)//66
    result = {}
    for sector in c.common.cut.cut_sectors(cap):
        output = np.zeros((len(sector.complement_states)*ref,)*2, complex)
        for row in sector.indices:
            ix = [i*ref+r for i in row for r in range(ref)]
            output += density[np.ix_(ix, ix)]
        result[sector.charge] = output
    return result


def _native_branches(law, seam):
    states = c.common.regional.source()[0].states
    positions = {x:i for i,x in enumerate(states)}
    groups = {}
    for coin in (0, 1):
        for i, x in enumerate(states):
            out, (_, current) = c.local.native.native(x, seam, coin)
            label = current if law == 'C' else x[seam[0]]-x[seam[1]]
            groups.setdefault((coin,label,current), np.zeros((66,66),complex))[positions[out],i] = 1/math.sqrt(2)
    return groups


def _native_rotation(index, sign):
    # Read occupations directly, without production pairs or encoded targets.
    states = c.local.native.states(binary=True)
    positions = {x:i for i,x in enumerate(states)}
    a,b = c.SEAMS[index]
    result = np.eye(66,dtype=complex)
    for i,x in enumerate(states):
        if x[a] == 1 and x[b] == 0:
            y=list(x);y[a],y[b]=0,1
            j=positions[tuple(y)]
            result[np.ix_((i,j),(i,j))] = np.array([[1,sign],[-sign,1]])/math.sqrt(2)
    return result


def _literal(density, law, s=F(1), feedback=True):
    ref=len(density)//66
    clock=np.diag(np.exp(-1j*float(s)*np.log(list(map(float,c.local.reference())))))
    def stage(index,sign):
        r=_native_rotation(index,sign)
        return np.kron(r@clock@_native_rotation(index,-1),np.eye(ref))
    ready=stage(0,1)@density@stage(0,1).conj().T
    result={('abort',0):np.zeros((ref,ref),complex)}
    for record,matrix in _native_branches(law,c.SEAMS[0]).items():
        k=np.kron(matrix,np.eye(ref));rho=k@ready@k.conj().T
        u=stage(1,-1 if feedback and record[2]<0 else 1)
        rho=u@rho@u.conj().T
        result['abort',1,*record]=np.zeros((ref,ref),complex)
        for later,matrix in _native_branches(law,c.SEAMS[1]).items():
            k=np.kron(matrix,np.eye(ref))
            result['read',*record,*later]=k@rho@k.conj().T
    return result


def test_exact_bare_cap_membership_and_selected_exterior_countercontrol():
    certificate=c.exact_algebra()
    assert certificate['completed_net_operations_in_bare_cap_algebras']
    assert all(row['selected_exterior_generator_fails_bare_cap_membership'] for row in certificate['rows'])
    for index in (0,1):
        assert len(c.pairs(index))==10
        generator={address:sign for i,j,_,_ in c.pairs(index)
                   for address,sign in (((i,j),1),((j,i),-1))}
        assert c._in_algebra(generator,c.SEAMS[index])
        # Omitting only one spectator pair destroys the local algebra.
        i,j,_,_=c.pairs(index)[0]
        del generator[i,j];del generator[j,i]
        assert not c._in_algebra(generator,c.common.CAPS[index])


@pytest.mark.parametrize('law',('C','E'))
@pytest.mark.parametrize('index',(0,1))
def test_unread_local_operation_preserves_complete_outside_reference_state(law,index):
    rho=_random_density()
    before=_outside(rho,c.common.CAPS[index])
    branches,abort=c._read(rho,index,1,law,False,F(1))
    assert not np.any(abort)
    after=_outside(sum(branches.values()),c.common.CAPS[index])
    for charge in before:
        np.testing.assert_allclose(after[charge],before[charge],rtol=0,atol=2e-15)
    # A retained outcome can steer the correlated exterior; never claim
    # outcome-by-outcome decoupling from that exterior.
    assert any(sum(np.linalg.norm(a-np.trace(branch).real*before[q],'nuc')
                   for q,a in _outside(branch,c.common.CAPS[index]).items())>1e-3
               for branch in branches.values())


@pytest.mark.parametrize('law',('C','E'))
def test_remote_unitary_changes_neither_local_record_probabilities_nor_local_channel(law):
    rho=_random_density(ref=1)
    cap=c.common.CAPS[0]
    sectors=c.common.cut.cut_sectors(cap)
    u=np.eye(66,dtype=complex)
    # The SAME outside rotation on every cap configuration in q_cap=1.
    sector=next(x for x in sectors if x.charge==1)
    block=np.eye(len(sector.complement_states),dtype=complex)
    block[:2,:2]=np.array([[1,1j],[1j,1]])/math.sqrt(2)
    for row in sector.indices:u[np.ix_(row,row)]=block
    baseline,_=c._read(rho,0,1,law,False,F(1))
    changed,_=c._read(u@rho@u.conj().T,0,1,law,False,F(1))
    for record in baseline:
        np.testing.assert_allclose(changed[record],u@baseline[record]@u.conj().T,rtol=0,atol=3e-15)
        assert abs(np.trace(changed[record])-np.trace(baseline[record]))<3e-15


@pytest.mark.parametrize('s',(F(0),F(1,2),F(1),F(4),F(-4)))
def test_both_isolated_caps_have_the_counted_two_frequency_current_law(s):
    sigma=np.diag(list(map(float,c.local.reference()))).astype(complex)
    expected=c.current_law(s)
    for index in (0,1):
        for law in ('C','E'):
            branches,_=c._read(sigma,index,1,law,False,s)
            actual={j:sum(np.trace(r).real for key,r in branches.items() if key[2]==j) for j in (-1,0,1)}
            for j in expected:assert abs(actual[j]-expected[j])<2e-14
    assert abs(expected[-1]+expected[1]-181/900)<1e-15


def test_mean_clock_inverse_is_scoped_to_a_fixed_monotone_window():
    certificate=c.calibration_certificate()
    assert F(certificate['absolute_mean_derivative_lower'])>F(1,200)
    values=[]
    for n in range(51,101):
        s=F(n,100)
        values.append(F(5,33)*(c.current_law(s)[1]-c.current_law(s)[-1]))
    assert all(a>b for a,b in zip(values,values[1:]))
    assert all(abs(a-b)>.01/200 for a,b in zip(values,values[1:]))
    # Local frequency readouts alone do not choose a global time orientation.
    assert c.current_law(F(1))==c.current_law(F(-1))


@pytest.fixture(scope='module')
def source_budget():
    return c.source_budget()


def test_original_whole_source_rotations_include_every_spectator_and_padding(source_budget):
    for index in (0,1):
        for sign in (-1,1):
            array=c.rotation(index,sign)
            assert not array.flags.writeable
            assert c.certify_rotation(index,sign)<c.ROTATION_GATE
            np.testing.assert_allclose(array,c.ideal_rotation(index,sign),rtol=0,atol=2e-11)
        wrong=c.rotation(index,1).copy();wrong[511,511]+=.125
        assert c.certify_rotation(index,1,wrong)>F(1,9)
    assert source_budget['maximum_source_pulses']==45656
    assert not source_budget['exterior_configuration_selection']
    # A caller cannot poison later certificate values or costs.
    original=c.source_budget();original['channel_norm_upper']=F(0)
    assert c.source_budget()==source_budget


@pytest.mark.parametrize('law',('C','E'))
def test_continuing_reference_history_matches_independent_native_oracle(source_budget,law):
    rho=_random_density()
    original=rho.copy()
    expected=_literal(rho,law)
    for raw,tolerance in ((False,3e-15),(True,2e-11)):
        actual=c.execute(rho,law,raw=raw)
        assert actual.keys()==expected.keys()
        for key in expected:np.testing.assert_allclose(actual[key],expected[key],rtol=0,atol=tolerance)
        assert c.history.distance(actual,expected)<float(source_budget['complete_history_half_diamond_upper'])
    np.testing.assert_array_equal(rho,original)
    assert c.history.distance(expected,_literal(rho,law,feedback=False))>1e-3


@pytest.mark.parametrize('law',('C','E'))
def test_complete_preparation_and_native_feedback_keep_all_failures(source_budget,law):
    actual=c.prepared_history(law)
    expected=c.target_history(law)
    b=c.local.budget()
    bound=source_budget['channel_norm_upper']*b['complete_flagged_reference_halftrace_upper']+F(5,33)*source_budget['complete_history_half_diamond_upper']
    assert c.history.distance(actual,expected)<float(bound)
    assert abs(sum(c.history.probabilities(actual).values())-1)<2e-10
    assert all(('preparation',flag) in actual for flag in c.local.FLAGS if flag!='success')
    assert any(key[:2]==('abort',1) for key in actual)
    assert ('abort',0) in actual
    # The second marginal belongs to a continuing disturbed carrier.
    second_mean=sum(key[-1]*np.trace(rho).real for key,rho in expected.items() if key[0]=='read')
    first_mean=F(5,33)*(c.current_law(F(1))[1]-c.current_law(F(1))[-1])
    assert abs(second_mean-first_mean)>1e-4


@pytest.mark.parametrize('law',('C','E'))
def test_both_export_aborts_retain_the_untouched_reference(monkeypatch,law):
    # A deliberately escaping unitary tests the complete export instrument,
    # not the accuracy of the admitted source word. Every native amplitude
    # splits equally into native and distinct padding coordinates.
    codes=c.full.support()[0]
    padding=[i for i in range(512) if i not in set(codes)][:66]
    escape=np.eye(512,dtype=complex)
    for a,b in zip(codes,padding,strict=True):
        escape[np.ix_((a,b),(a,b))]=np.array([[1,1],[-1,1]])/math.sqrt(2)
    monkeypatch.setattr(c,'_stage',lambda *args:(escape.dtype.str,escape.tobytes()))
    rho=_random_density()
    outputs=c.execute(rho,law,raw=False)
    assert abs(np.trace(outputs['abort',0]).real-.5)<2e-15
    second=[value for key,value in outputs.items() if key[:2]==('abort',1)]
    assert abs(sum(np.trace(value).real for value in second)-.25)<2e-15
    total_reference=np.zeros((2,2),complex)
    for key,value in outputs.items():
        total_reference+=value if key[0]=='abort' else np.einsum('aras->rs',value.reshape(66,2,66,2))
    original_reference=np.einsum('aras->rs',rho.reshape(66,2,66,2))
    np.testing.assert_allclose(total_reference,original_reference,rtol=0,atol=3e-15)


@pytest.mark.parametrize('bad',((True,2,1),(0,512,1),(0,0,1),(0,1,True),(0.,1,1)))
def test_pair_compiler_validates_original_addresses_and_sign(bad):
    with pytest.raises(ValueError):c.pair_word(*bad)


def test_raw_history_parameter_and_density_contracts_fail_closed():
    rho=np.eye(66)/66
    for kw in ({'s':F(2)},{'raw':1},{'feedback':0}):
        with pytest.raises(ValueError):c.execute(rho,'C',**kw)
    with pytest.raises(ValueError):c.execute(np.ma.array(rho),'C')
    with pytest.raises(ValueError):c.current_law(float('nan'))


def test_receipt_rejects_negative_probabilities_and_bad_complete_mass(source_budget):
    receipt=c.strict_json.loads(c.full.preparation.canonical(c.run()))
    c.compare_receipt(receipt,receipt)
    for quantity in ('complete_prepared_history_halftrace_diagnostic','complete_prepared_history_mass_diagnostic'):
        wrong=deepcopy(receipt);wrong['histories'][0][quantity]=-1e-12
        with pytest.raises(ValueError):c.compare_receipt(wrong,wrong)
    wrong=deepcopy(receipt);wrong['histories'][0]['isolated_cap_current_laws'][0]['0']=-1e-12
    with pytest.raises(ValueError):c.compare_receipt(wrong,wrong)
