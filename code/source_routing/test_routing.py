"""Adversarial scientific controls for local-only destructive transport."""
import copy
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

import pytest

HERE=Path(__file__).resolve().parent


def load_module(name):
    spec=importlib.util.spec_from_file_location(name,HERE/f'{name}.py')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


verify=load_module('verify_routing')
producer=load_module('build_routing')


@pytest.fixture
def packet():
    return json.loads((HERE/'runtime/path_tomography_receipt.json').read_text())


def rehash(episode):
    """Rechain mutations so a hash mismatch cannot substitute for science checks."""
    previous='0'*64
    for eid,event in enumerate(episode['events']):
        event['id']=eid;event['previous_hash']=previous
        event.pop('event_hash',None)
        event['event_hash']=hashlib.sha256(verify.raw(event)).hexdigest()
        previous=event['event_hash']
    episode['final_event_hash']=previous


def test_complete_independent_replay(packet):
    assert len(verify.verify(packet))==12


def test_producer_reproducible(packet):
    assert producer.make()==packet


def test_deleted_calibration_rejected_even_if_rehashed(packet):
    ep=packet['episodes'][3]
    first=next(i for i,e in enumerate(ep['events']) if e['op']=='pair_mean')
    del ep['events'][first]
    rehash(ep)
    with pytest.raises(ValueError,match='event loss'):
        verify.verify(packet)


def test_wrong_law_rejected_even_if_rehashed(packet):
    ep=packet['episodes'][0]
    mean=next(e for e in ep['events'] if e['op']=='pair_mean')
    for w in mean['writes']:w['value']=mean['reads'][0]['value']
    rehash(ep)
    with pytest.raises(ValueError,match='wrong law'):
        verify.verify(packet)


def test_wrong_wiring_rejected_even_if_rehashed(packet):
    packet['episodes'][0]['path'][1]=15359
    with pytest.raises(ValueError,match='absent seam'):
        verify.verify(packet)


def test_hidden_remote_read_rejected_even_if_rehashed(packet):
    ep=packet['episodes'][0]
    observed=next(e for e in ep['events'] if e['op']=='local_readback')
    # A valid existing remote register still violates the destination firewall.
    observed['reads']=[copy.deepcopy(ep['events'][0]['writes'][0])]
    rehash(ep)
    with pytest.raises(ValueError,match='side-channel'):
        verify.verify(packet)


def test_stale_version_rejected_even_if_rehashed(packet):
    ep=packet['episodes'][3]
    means=[e for e in ep['events'] if e['op']=='pair_mean']
    means[-1]['reads'][0]['version']-=1
    rehash(ep)
    with pytest.raises(ValueError,match='version, writer or value'):
        verify.verify(packet)


def test_omitted_variant_rejected(packet):
    packet['episodes'].pop()
    with pytest.raises(ValueError,match='scheduled outcome'):
        verify.verify(packet)


def test_scope_promotion_rejected(packet):
    packet['scope']['q13_q21_complete_routing']=True
    with pytest.raises(ValueError,match='scope promotion'):
        verify.verify(packet)


def test_optimistic_noise_claim_rejected(packet):
    packet['episodes'][-1]['conditioning']['source_error_per_uniform_sample_error']='1'
    with pytest.raises(ValueError,match='error amplification'):
        verify.verify(packet)


def test_capacity_undercount_rejected(packet):
    packet['episodes'][-1]['cost']['protected_sample_words']=1
    with pytest.raises(ValueError,match='capacity accounting'):
        verify.verify(packet)


def test_forward_sweep_without_calibration_has_real_ambiguity():
    # These two genuine three-port histories have identical receiver-local
    # initial and final records, while their original source values differ.
    def one_sweep(initial):
        x,y,z=map(F,initial)
        m=(x+y)/2
        return z,(m+z)/2
    assert one_sweep([1,3,5])==one_sweep([2,2,5])
    assert producer.decode_local([F(5),F(4),F(13,4)])[0][0] == 1


def test_rounding_can_hide_source_intervention(packet):
    base=packet['episodes'][9]; intervention=packet['episodes'][10]
    a=F(base['local_samples'][-1]);b=F(intervention['local_samples'][-1])
    assert b-a==F(1,2**24)
    assert round(a,4)==round(b,4)


def test_exact_observation_error_norm_is_attained():
    observation,_=verify.observation_matrix(12)
    row=verify.invert(observation)[0]
    epsilon=F(1,10**8)
    errors=[epsilon if x>0 else -epsilon if x<0 else F(0) for x in row]
    assert sum(x*y for x,y in zip(row,errors))==47321*epsilon


def test_duplicate_or_inexact_json_rejected(tmp_path):
    path=tmp_path/'bad.json'
    for contents in ('{"a":1,"a":2}', '{"a":NaN}', '{"a":1.0}'):
        path.write_text(contents)
        with pytest.raises(ValueError):
            verify.load(path)


def test_boolean_cost_is_not_an_integer_count(packet):
    packet['episodes'][0]['cost']['seam_mean_operations']=True
    with pytest.raises(ValueError,match='capacity accounting'):
        verify.verify(packet)


def test_exact_conditioning_law_from_actual_mean_matrices():
    """Compare the analytic all-depth law with independently inverted histories."""
    gains = [1, 3]
    for d in range(2, 25):
        gains.append(2*gains[-1] + gains[-2])
    for d in range(25):
        observation, _ = verify.observation_matrix(d)
        assert all(x >= 0 for row in observation for x in row)
        assert all(sum(row) == 1 for row in observation)
        inverse = verify.invert(observation)
        norms = [sum(abs(x) for x in row) for row in inverse]
        assert norms == list(reversed(gains[:d+1]))
        # The source row alternates by sample delay, including exact zeros.
        assert all(x*(-1)**(d-j) >= 0 for j, x in enumerate(inverse[0]))
        errors = [F((-1)**(d-j), 10**8) for j in range(d+1)]
        assert sum(x*y for x, y in zip(inverse[0], errors)) == F(gains[d], 10**8)
