"""Plausible false-green rewrites of the exact population experiment."""
import copy,hashlib,importlib.util,json
from pathlib import Path
import pytest
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('population_independent',HERE/'verify_population.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
@pytest.fixture(scope='module')
def packet():return v.load(HERE/'runtime/population_receipt.json')
def reseal(p):
 previous='0'*64
 for e in p['events']:
  e['previous_hash']=previous
  e['event_hash']=hashlib.sha256(v.raw({k:x for k,x in e.items() if k!='event_hash'})).hexdigest();previous=e['event_hash']
 p['final_event_hash']=previous

def test_complete_exact_replay(packet):assert v.check(packet)['accepted']

def test_rehashed_wrong_mean(packet):
 p=copy.deepcopy(packet);e=next(e for e in p['events'] if e['op']=='pair_mean');key=next(iter(e['writes']));e['writes'][key]='0';reseal(p)
 with pytest.raises(ValueError,match='pair mean'):v.check(p)

def test_rehashed_missing_parent(packet):
 p=copy.deepcopy(packet);next(e for e in p['events'] if e['op']=='pair_mean')['parents']=[];reseal(p)
 with pytest.raises(ValueError,match='read-from'):v.check(p)

def test_rehashed_false_decrement(packet):
 p=copy.deepcopy(packet);next(e for e in p['events'] if e['op']=='pair_mean')['quadratic_decrement']='0';reseal(p)
 with pytest.raises(ValueError,match='potential'):v.check(p)

def test_forged_live_address(packet):
 p=copy.deepcopy(packet);p['sites'][1]['live_address_after']=p['sites'][1]['protected_address_after']
 with pytest.raises(ValueError,match='live readout'):v.check(p)

def test_erased_retained_address(packet):
 p=copy.deepcopy(packet);p['sites'][1]['protected_address_after']=[['0','0']]*3
 with pytest.raises(ValueError,match='protected write'):v.check(p)

def test_wrong_metric_census(packet):
 p=copy.deepcopy(packet);p['summary']['changed_unordered_metric_relations']=0
 with pytest.raises(ValueError,match='exact census'):v.check(p)

def test_physical_promotion(packet):
 p=copy.deepcopy(packet);p['scope']['physical_premise_discharged']=True
 with pytest.raises(ValueError,match='scope boundary'):v.check(p)

def test_missing_site(packet):
 p=copy.deepcopy(packet);p['sites'].pop()
 with pytest.raises(ValueError,match='census'):v.check(p)

def test_changed_source_pin(packet):
 p=copy.deepcopy(packet);p['source_reference']['sha256']='0'*64
 with pytest.raises(ValueError,match='custody'):v.check(p)

def test_duplicate_json_rejected(tmp_path):
 p=tmp_path/'bad.json';p.write_text('{"q":13,"q":21}')
 with pytest.raises(ValueError,match='duplicate'):v.load(p)

@pytest.mark.parametrize('mutate,match',[
 (lambda p:p.update(schema='invented.schema'),'schema version'),
 (lambda p:p['scope'].update(physical_premise_discharged=0),'boolean type'),
 (lambda p:p['events'][0].update(undeclared_write='0'),'event schema'),
 (lambda p:p['sites'][0].update(site=False),'unexpected boolean'),
])
def test_strict_schema(packet,mutate,match):
 p=copy.deepcopy(packet);mutate(p)
 with pytest.raises(ValueError,match=match):v.check(p)

def test_overflow_float_rejected(tmp_path):
 p=tmp_path/'bad.json';p.write_text('{"q":1e400}')
 with pytest.raises(ValueError,match='floating'):v.load(p)
