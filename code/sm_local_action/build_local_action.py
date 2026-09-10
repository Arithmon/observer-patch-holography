"""Emit exact finite local-action, Ward and variational-current evidence."""
from pathlib import Path
import argparse,hashlib,json
import jet_action as a
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];OUTPUT=HERE/'local_action_receipt.json'
OWN=('jet_action.py','build_local_action.py','verify_local_action.py','test_local_action.py')
PARENTS=('code/a5_closure/exterior_sm_completion.json','code/e9_kinetic/gauge_kinetic_invariant_forms.certificate.json',
 'Lean/Screen/WeylYukawaConventions.lean','Lean/Screen/LocalGaugeJetAction.lean','paper/tex_fragments/LOCAL_SM_JET_ACTION.tex',
 'code/particles/calibration/wz_upstream_completion/outputs/sm_eft_action_1.json')
SCOPE={'conditional_local_classical_action':True,'three_families_declared':True,'fermions_anticommuting':True,'all_complex_Yukawa_entries_symbolic':True,
 'local_parameter_derivatives_retained':True,'same_action_variational_matter_current':True,'source_connection_selected':False,
 'source_action_emitted':False,'physical_Spin_attached':False,'coupling_or_mass_values_selected':False,'broken_vacuum_selected':False,
 'curved_spacetime_action':False,'renormalized_quantum_theory':False,'quantum_anomaly_from_measure_proved':False,
 'full_action_formalized_in_Lean':False,'empirical_comparison':False,'current_joined_to_source_population':False}
PREMISES=['flat oriented Minkowski chart and supplied two-component Spin convention',
 'declared local anti-Hermitian SU3/SU2/U1 connections on the registered fundamental/dual/singlet representations',
 'representation-to-classified-bracket identification and invariant trace normalization are declared; no response-history emission',
 'minimal first-derivative Hermitian Weyl and scalar kinetic law',
 'Higgs doublet H=(1,2)_(1/2), scalar degree<=4 polynomial, three families and all three Yukawa intertwiners declared',
 'positive free gauge trace coefficients; real mass-squared and quartic parameters; unrestricted complex Yukawa matrices',
 'smooth gauge parameter with symmetric second derivatives; compact support for the integrated Noether argument']
def pin(path):
 b=path.read_bytes();return {'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)}
def enc(x):return [enc(v) for v in x] if isinstance(x,list) else x.encode()
def scalar_jet_packet(j):return {k:enc(v) for k,v in j.items()}
def parameter(k):
 if k.startswith('kinetic_') or k=='scalar_kinetic':return {'kind':'declared_minimal_normalization','value':'1'}
 return {'kind':'free_symbol','name':k,'domain':'positive real' if k.startswith('kappa_') else 'real'}
def nonzero(action):return {k:len(v.delta.t) for k,v in sorted(action.items()) if v.delta}
def produce():
 algebra=a.Algebra(3);jets=a.fixture();base=a.compile_action(algebra,jets)
 if any(v.value!=v.value.adjoint(algebra.partners) for v in base.values()):raise ValueError('Hermitian density')
 ward=[]
 for kind in ('theta','dtheta'):
  for g in range(12):
   for mu in (range(4) if kind=='dtheta' else [0]):
    computed=a.compile_action(algebra,jets,a.direction(kind,g,mu));res=nonzero(computed)
    if res:raise ValueError('nonzero local Ward variation '+str((kind,g,mu,res)))
    ward.append({'kind':kind,'generator':g,'mu':mu,'nonzero_sectors':res})
 hessian=[]
 for g in range(12):
  for mu in range(4):
   for nu in range(mu,4):
    # The supplied second parameter jet is symmetric; its two signed terms cancel.
    _,_,d2,_=a.direction('d2theta',g,mu,nu)
    if any(any(x for row in a.madd(a.mscale(-1,d2[u][v][s]),d2[v][u][s]) for x in row) for u in range(4) for v in range(4) for s in range(3)):raise ValueError('symmetric Hessian cancellation')
    hessian.append({'generator':g,'mu':mu,'nu':nu,'curvature_variation_zero':True})
 currents={}
 for mu in range(4):
  for g in range(12):
   computed=a.compile_action(algebra,jets,a.direction('connection',g,mu))
   current=sum((v.delta for k,v in computed.items() if k.startswith('kinetic_') or k=='scalar_kinetic'),a.P())
   if not current:raise ValueError('vacuous matter current')
   currents[f'{mu}/{g}']=current.encode()
 controls=[]
 for variant,kind,g,mu in [('drop_dtheta_connection','dtheta',11,0),('drop_A_commutator','dtheta',0,1),('wrong_dual_color','theta',1,0),('mixed_up_hypercharge','theta',11,0)]:
  bad=a.compile_action(algebra,jets,a.direction(kind,g,mu),variant);r=nonzero(bad)
  if not r:raise ValueError('mutation did not fail '+variant)
  controls.append({'variant':variant,'kind':kind,'generator':g,'mu':mu,'nonzero_sector_terms':r,'residual_digest':a.digest(a.encode_action(bad,'delta'))})
 x=a.P.generator(0);y=a.P.generator(2)
 exterior=json.loads((ROOT/PARENTS[0]).read_text());external=json.loads((ROOT/PARENTS[-1]).read_text())
 return {'schema':'oph.local_sm_jet_action.v1','scope':SCOPE,'premises':PREMISES,
  'conventions':{'metric':[1,-1,-1,-1],'connection':'anti-Hermitian; D=partial+rho(A), F=dA+[A,A]','parameter':'delta psi=theta psi; delta A=-dtheta+[theta,A]',
   'fermions':'all left Weyl; independent Grassmann conjugate and first-derivative generators','adjoint':'antilinear, reverses order, exchanges adjacent value/conjugate generator pairs',
   'Pauli':'bar_sigma=(I,-sigma1,-sigma2,-sigma3); epsilon01=1','gauge_normalization':'kappa_s sum_(mu<nu) eta_mu eta_nu Tr(F_mu_nu squared); trace bases are not orthonormal',
   'scalar_potential':'V=m_squared Hdagger H+lambda(Hdagger H)^2; no broken branch selected',
   'Yukawa':'-Y T-conj(Y)Tdagger, represented by independent real and imaginary coefficient sectors',
   'units':'formal local action/jet conventions; no laboratory calibration'},
  'representation_fields':exterior['matter_package']['fields'],'representation_boundary':exterior['conditional_closure']['physical_status'],
  'fixture':{'kind':'exact algebraic test jets; not measured or source-emitted fields','seed':0,'families':3,'jets':scalar_jet_packet(jets)},
  'generator_labels':[list(x) for x in algebra.labels],'sector_coefficients':{k:parameter(k) for k in sorted(base)},
  'action_coefficients':a.encode_action(base),'ward_directions':ward,'symmetric_second_parameter_jets':hessian,'matter_current_coefficients':currents,
  'negative_controls':controls,'grassmann_controls':{'square':(x*x).encode(),'anticommutator':(x*y+y*x).encode(),'spin_pair':(x*y-y*x).encode(),'adjoint_product':(x*y).adjoint(algebra.partners).encode()},
  'comparison':{'external_action_role':external['claim_boundary'],'physical_operator_shapes_assembled':['gauge_kinetic_SU3','gauge_kinetic_SU2','gauge_kinetic_U1','fermion_kinetic','higgs_kinetic','higgs_mass','higgs_quartic','yukawa_up','yukawa_down','yukawa_lepton'],
   'not_assembled':['gauge_fixing','ghost_sector'],'no_external_parameter_values_consumed':True,
   'mass_sign_dictionary':'external +mu2 |H|^2 corresponds to local -m_squared |H|^2 with mu2=-m_squared',
   'Yukawa_comparison':'same neutral invariant channels under charge conjugation and Hermitian conjugation; external coefficient symbols are not numerically identified'},
  'source_pins':{x:pin(HERE/x) for x in OWN},'parent_pins':{x:pin(ROOT/x) for x in PARENTS}}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();packet=produce();raw=a.canonical(packet)
 if args.write:OUTPUT.write_bytes(raw)
 else:
  if OUTPUT.read_bytes()!=raw:raise SystemExit('local action receipt differs')
 print(json.dumps({'accepted':True,'sectors':len(packet['action_coefficients']),'Ward_directions':len(packet['ward_directions']),'second_parameter_jets':len(packet['symmetric_second_parameter_jets']),'currents':len(packet['matter_current_coefficients']),'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}))
