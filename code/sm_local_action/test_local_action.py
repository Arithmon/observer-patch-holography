"""Independent action replay and non-vacuous local-Ward/Grassmann controls."""
import copy
import json
from fractions import Fraction
from pathlib import Path
import sys
import pytest

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import jet_action as producer
import verify_local_action as verifier

@pytest.fixture(scope='module')
def receipt():
 return verifier.load()

def test_fresh_independent_replay(receipt):
 report=verifier.verify(receipt)
 assert report['action_sectors']==65
 assert report['action_monomials']==2295
 assert report['Grassmann_generators']==900
 assert report['local_Ward_directions']==60
 assert report['symmetric_second_parameter_jets']==120
 assert report['variational_matter_currents']==48
 assert not report['source_action_emission']

def test_cheap_custody_does_not_claim_mathematical_replay(receipt):
 result=verifier.validate_custody(receipt)
 assert result['accepted_custody'] and result['mathematical_replay'] is False

@pytest.mark.parametrize('field',[
 'source_action_emitted','physical_Spin_attached','coupling_or_mass_values_selected',
 'renormalized_quantum_theory','empirical_comparison','current_joined_to_source_population',
 'full_action_formalized_in_Lean'],ids=lambda x:x)
def test_no_scientific_promotion(receipt,field):
 bad=copy.deepcopy(receipt);bad['scope'][field]=True
 with pytest.raises(ValueError,match='scientific scope'):verifier.validate_custody(bad)

@pytest.mark.parametrize('case',['source_hash','source_size','parent_hash','parent_size','extra_source','missing_parent'])
def test_stale_or_resealed_parent_custody_rejected(receipt,case):
 bad=copy.deepcopy(receipt)
 if case=='extra_source':bad['source_pins']['invented.py']={'sha256':'0'*64,'bytes':0}
 elif case=='missing_parent':del bad['parent_pins'][next(iter(bad['parent_pins']))]
 else:
  kind='source_pins' if case.startswith('source') else 'parent_pins';p=next(iter(bad[kind]))
  key='sha256' if case.endswith('hash') else 'bytes';bad[kind][p][key]='0'*64 if key=='sha256' else True
 with pytest.raises(ValueError,match='pins'):verifier.validate_custody(bad)

@pytest.mark.parametrize('case',['premise','units','parameter','representation','fixture','family_bool','generator_label','missing_dtheta','dtheta_bool','missing_hessian','asymmetric_hessian','mass_sign','numeric_yukawa','missing_offdiagonal_yukawa','missing_current','claimed_zero_negative','wrong_negative_direction','extra_root'])
def test_scope_and_complete_jet_schema(receipt,case):
 bad=copy.deepcopy(receipt)
 if case=='premise':bad['premises'].pop()
 elif case=='units':bad['conventions']['units']='physical seconds and measured GeV'
 elif case=='parameter':bad['conventions']['parameter']='constant gauge parameter only'
 elif case=='representation':bad['representation_fields']['u_c']['Y']='2/3'
 elif case=='fixture':bad['fixture']['jets']['H'][0][0]='1'
 elif case=='family_bool':bad['fixture']['families']=True
 elif case=='generator_label':bad['generator_labels'][0][-1]=0
 elif case=='missing_dtheta':bad['ward_directions'].pop()
 elif case=='dtheta_bool':bad['ward_directions'][0]['mu']=False
 elif case=='missing_hessian':bad['symmetric_second_parameter_jets'].pop()
 elif case=='asymmetric_hessian':bad['symmetric_second_parameter_jets'][0]['nu']=1
 elif case=='mass_sign':bad['comparison']['mass_sign_dictionary']='mu2=m_squared'
 elif case=='numeric_yukawa':bad['sector_coefficients']['Y_up_0_1_re']={'kind':'measured','value':'1'}
 elif case=='missing_offdiagonal_yukawa':del bad['action_coefficients']['Y_up_0_1_re']
 elif case=='missing_current':del bad['matter_current_coefficients']['0/11']
 elif case=='claimed_zero_negative':bad['negative_controls'][0]['nonzero_sector_terms']={}
 elif case=='wrong_negative_direction':bad['negative_controls'][0]['kind']='theta'
 elif case=='extra_root':bad['new_claim']=True
 with pytest.raises(ValueError):verifier.validate_custody(bad)

@pytest.mark.parametrize('case',['bool_index','negative_index','large_index','duplicate_index','duplicate_word','unsorted_word','noncanonical_fraction','float_token','complex_token','zero_term'])
def test_strict_exterior_coefficient_schema(receipt,case):
 bad=copy.deepcopy(receipt);terms=bad['action_coefficients']['kinetic_Q'];term=terms[0]
 if case=='bool_index':term[0][0]=False
 elif case=='negative_index':term[0][0]=-1
 elif case=='large_index':term[0][-1]=900
 elif case=='duplicate_index':term[0]=[0,0]
 elif case=='duplicate_word':terms.append(copy.deepcopy(term))
 elif case=='unsorted_word':term[0]=[2,0]
 elif case=='noncanonical_fraction':term[1][0]='2/2'
 elif case=='float_token':term[1][0]=0.5
 elif case=='complex_token':term[1][0]='1+2j'
 elif case=='zero_term':term[1]=['0','0']
 with pytest.raises((ValueError,TypeError)):verifier.validate_custody(bad)

@pytest.mark.parametrize('case',['kinetic_sign','offdiagonal_yukawa','missing_derivative_term','current_sign','commuting_spin_pair','residual_digest'])
def test_mathematical_forgeries_rejected_by_full_verifier(receipt,case):
 bad=copy.deepcopy(receipt)
 if case in ('kinetic_sign','offdiagonal_yukawa'):
  k='kinetic_Q' if case=='kinetic_sign' else 'Y_up_0_1_im'
  coeff=bad['action_coefficients'][k][0][1];index=0 if coeff[0]!='0' else 1;coeff[index]=str(-Fraction(coeff[index]))
 elif case=='missing_derivative_term':
  terms=bad['action_coefficients']['kinetic_Q'];labels=bad['generator_labels']
  ix=next(i for i,(word,c) in enumerate(terms) if any(labels[n][4]!='value' for n in word));terms.pop(ix)
 elif case=='current_sign':
  term=next(t for t in bad['matter_current_coefficients']['0/11'] if t[0]);coeff=term[1];ix=0 if coeff[0]!='0' else 1;coeff[ix]=str(-Fraction(coeff[ix]))
 elif case=='commuting_spin_pair':bad['grassmann_controls']['spin_pair']=[]
 elif case=='residual_digest':bad['negative_controls'][0]['residual_digest']='0'*64
 # Self-consistent schemas are not mathematical certificates.
 assert verifier.validate_custody(bad)['mathematical_replay'] is False
 with pytest.raises(ValueError,match='coefficients|current formula|Grassmann sign|negative Ward'):verifier.verify(bad)

@pytest.mark.parametrize('raw',[b'{"a":1,"a":2}',b'{"x":NaN}',b'{"x":Infinity}',b'{"x":0.5}',b'{' ],ids=['duplicate','nan','infinity','float','truncated'])
def test_strict_json_loader(tmp_path,raw):
 path=tmp_path/'bad.json';path.write_bytes(raw)
 with pytest.raises((ValueError,json.JSONDecodeError)):verifier.load(path)

def test_maximum_packet_size(tmp_path):
 path=tmp_path/'large.json';path.write_bytes(b' '*10000000)
 with pytest.raises(ValueError,match='packet size'):verifier.load(path)

@pytest.mark.parametrize('x',[True,False,0.5,complex(1,0),'1',float('nan')])
def test_exact_producer_scalars_reject_lossy_casts(x):
 with pytest.raises(TypeError):producer.C(x)

@pytest.mark.parametrize('x',[True,False,-1,4,0.5])
def test_family_and_fixture_integer_inputs(x):
 with pytest.raises(ValueError):producer.Algebra(x)
 if x not in (4,):
  if type(x) is not int or x<0:
   with pytest.raises(ValueError):producer.fixture(x)

@pytest.mark.parametrize('args',[('theta',True,0),('theta',-1,0),('theta',12,0),('dtheta',0,False),('dtheta',0,-1),('dtheta',0,4),('invented',0,0)])
def test_direction_inputs_reject_aliases(args):
 with pytest.raises(ValueError):producer.direction(*args)

def test_exterior_adjoint_and_noncommutativity_independent_replay():
 a=producer.Algebra(1);x=producer.P.generator(0);y=producer.P.generator(2)
 assert not x*x and not (x*y+y*x)
 assert (x*y-y*x).encode()==[[[0,2],['2','0']]]
 assert (x*y).adjoint(a.partners).encode()==[[[1,3],['-1','0']]]
 z=verifier.W({1:verifier.ONE});w=verifier.W({4:verifier.ONE})
 assert (z*w-w*z).encode()==(x*y-y*x).encode()
 assert (z*w).star().encode()==(x*y).adjoint(a.partners).encode()
 assert (z*w).star().star()==z*w

def test_all_electroweak_current_components_equal_direct_formula():
 # One smaller compiler family is unsuitable here: the independent formula
 # intentionally sums all three declared families.
 a=producer.Algebra(3);j=producer.fixture();other=verifier.fixture()
 for mu,g in [(0,11),(1,8),(2,9),(3,10)]:
  computed=producer.compile_action(a,j,producer.direction('connection',g,mu))
  result=sum((v.delta for k,v in computed.items() if k.startswith('kinetic_') or k=='scalar_kinetic'),producer.P())
  assert result.encode()==verifier.direct_current(other,mu,g).encode()

def test_spacetime_parameter_derivative_is_essential():
 a=producer.Algebra(1);j=producer.fixture(1);d=producer.direction('dtheta',11,2)
 good=producer.compile_action(a,j,d)
 bad=producer.compile_action(a,j,d,'drop_dtheta_connection')
 assert not any(v.delta for v in good.values())
 assert bad['kinetic_Q'].delta and bad['scalar_kinetic'].delta
