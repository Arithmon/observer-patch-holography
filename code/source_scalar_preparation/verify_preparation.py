"""Independent exact checks in the 1,sqrt(5) basis; no producer import.

The pinned parent certificates are identified, not rerun here. Their own
independent verifiers remain the authorities for the inherited instrument.
"""
from __future__ import annotations
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PARENTS={
 'code/source_scalar_execution/source_scalar_execution_receipt.json':'6981fcd4e13fbc33f371224228341a995f34b4733aa70c2ec0b3f6ee5975c7af',
 'code/source_scalar_finite_instrument/finite_instrument_receipt.json':'a9e922b62896419bb96ac24d9023c24f1c7f4f0cc40ca84c601ece12577a09ad',
}
ZERO=(F(0),F(0))

def add(x,y): return x[0]+y[0],x[1]+y[1]
def neg(x): return -x[0],-x[1]
def mul(x,y): return x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def scalar(x): return F(x),F(0)
def phi(x):
    a,b=map(F,x)
    return a+b/2,b/2
def total(xs):
    out=ZERO
    for x in xs: out=add(out,x)
    return out
def sign(x):
    a,b=x
    if a==0:return (b>0)-(b<0)
    if b==0:return (a>0)-(a<0)
    if a*b>0:return (a>0)-(a<0)
    square=a*a-5*b*b
    return ((square>0)-(square<0))*((a>0)-(a<0))
def require(test,label):
    if not test: raise ValueError(label)

def check(value=None,root=ROOT):
    if value is None:value=json.loads((HERE/'preparation_receipt.json').read_text())
    require(value['schema']=='oph.scalar-finite-local-preparation.v1','schema')
    require(value['parents']==PARENTS,'parent identities')
    parent=[]
    for path,digest in PARENTS.items():
        data=(root/path).read_bytes();require(sha256(data).hexdigest()==digest,'parent bytes')
        parent.append(json.loads(data))
    p,ins=parent
    mass=list(map(phi,p['mass_Qphi']));v=list(map(phi,p['initial_velocity_Qphi']))
    rows=[{j:phi(x) for j,x in row} for row in p['action_rows_Qphi']]
    require(len(rows)==len(mass)==len(v)==64,'full carrier')
    for i,row in enumerate(rows):
        require(sign(mass[i])>0,'positive mass')
        require(sign(add(total(row.values()),scalar(-1)))>=0,'A lower bound')
        for j,x in row.items():
            require(i==j or sign(x)<=0,'weighted graph sign')
            require(mul(mass[i],x)==mul(mass[j],rows[j].get(i,ZERO)),'mass self-adjointness')
    av=[total(mul(a,v[j]) for j,a in row.items()) for row in rows]
    norm0=total(mul(m,mul(x,x)) for m,x in zip(mass,v))
    norm2=total(mul(m,mul(x,x)) for m,x in zip(mass,av))
    mo=value['moments']; require(phi(mo['v_mass_squared_Qphi'])==norm0,'velocity norm')
    require(phi(mo['Av_mass_squared_Qphi'])==norm2,'action velocity norm')
    upper=F(mo['Av_norm_upper']);quantum=F(1,10**15)
    require(upper>0 and sign(add(scalar(upper**2),neg(norm2)))>=0,'norm enclosure upper')
    require(sign(add(scalar((upper-quantum)**2),neg(norm2)))<0,'norm enclosure tightness')
    energy=F(mo['ideal_added_energy_upper'])
    require(sign(add(scalar(energy),neg(mul(scalar(F(1,2)),norm0))))>=0,'ideal energy upper')
    require(sign(add(scalar(energy-quantum),neg(mul(scalar(F(1,2)),norm0))))<0,'ideal energy tightness')
    require(value['law']=={
      'equation':'q_ddot + A q = v/delta + e(t) during [-delta/2,delta/2]',
      'Hamiltonian':'H0 - sum_i m_i*v_i*q_i/delta',
      'reference':'free evolution from D(momentum=M*v)|vacuum> at virtual t=0',
      'effective_velocity':'sinc(delta*sqrt(A)/2)*v',
      'ideal_trace_bound':'delta^2*||A v||_M/(24*sqrt(2)) <= delta^2*||A v||_M/24',
      'force_error_contract':'integral ||e(t)||_M dt <= force_L1_error',
      'force_trace_bound':'force_L1_error/sqrt(2) <= force_L1_error',
      'initial_error_contract':'trace distance from original vacuum at pulse start <= vacuum_error',
      'baseline_force':'0, with same declared residual force error contract'},'declared force law')
    par=value['parameters']
    require(par=={'delta_min':'1/10000','delta_max':'1/500','vacuum_error':'1/1000000','force_L1_error':'1/1000000'},'frozen parameter region')
    delta=F(par['delta_max']);pulse=delta**2*upper/24
    baseline=F(par['vacuum_error'])+F(par['force_L1_error'])
    intervention=baseline+pulse
    bound=value['bounds'];budget=F(ins['parameters']['preparation_trace_distance_max'])
    first=ins['schedule'][0];offset=first['operations'][0]['start_offset']
    first_start=F(first['center_time_interval'][0])+F(offset[0])*F(ins['parameters']['duration_max'])+F(offset[1])
    require(bound=={'pulse_trace_error_upper':str(pulse),'baseline_trace_error_upper':str(baseline),
       'intervention_trace_error_upper':str(intervention),'inherited_per_branch_preparation_budget':str(budget),
       'both_branches_fit_inherited_budget':True,'first_pointer_preparation_start_lower':str(first_start),
       'preparation_end_to_first_pointer_gap_lower':str(first_start-delta/2)},'preparation bounds')
    require(intervention<=budget and first_start>delta/2,'composition acceptance')
    support=[i for i,x in enumerate(v) if sign(x)]
    require(support==list(range(16)),'local force support')
    ops=value['control_schedule'];require(len(ops)==16,'complete force schedule')
    for i,op in zip(support,ops):
        require(op['site']==i and phi(op['velocity_Qphi'])==v[i] and phi(op['force_coefficient_m_v_Qphi'])==mul(mass[i],v[i]),'force coefficient or site')
        require(op['start_over_delta']=='-1/2' and op['stop_over_delta']=='1/2','centered finite pulse')
        require(op['control_reads']==['mass','velocity','duration'],'consumed control inputs')
        require(op['on_write']=='force=-m*v/delta multiplying q in H' and op['off_write']=='force=0','force on/off')
    res=value['resource_model']
    require(res['parallel_force_ports']==16 and res['coefficient_reads']==48 and res['on_off_writes']==32,'resource accounting')
    require(res['elapsed_preparation_time']=='delta, not 16*delta' and res['phase_and_outcome_records']=='none synthesized','execution scope')
    require(res['coefficient_bit_precision']=='exact Q(phi) law plus declared integrated force error; no hardware precision claim','precision scope')
    comp=value['instrument_composition']
    require(comp=={'all21_unconditional_marginals_retained':True,'independent_shots_assumed':False,
      'resolved_steps':ins['summary']['resolved_steps'],'parent_error_budgets_remain_valid':True,
      'step16_response_lower':ins['rows'][15]['noisy_paired_response_abs_lower'],
      'step16_total_error_upper':ins['rows'][15]['total_paired_response_error_upper']},'instrument composition')
    scope=value['scope']
    true={'same_original64_mode_action','free_field_active_throughout_preparation','local_force_support_equals_original16_site_intervention',
          'original_vacuum_at_pulse_start_supplied','force_control_law_and_error_budget_supplied','trace_comparison_valid_after_pulse_end_only','all_bounded_later_instruments_share_trace_error_bound'}
    false={'virtual_center_state_is_actual_midpulse_state','source_selected_quantization_or_vacuum','physical_clock_or_quantum_outcome_custody','actual_quantum_hardware_execution','noisy_energy_bound_from_trace_distance'}
    require(set(scope)==true|false and all(scope[k] is True for k in true) and all(scope[k] is False for k in false),'scientific scope')
    expected={f'code/source_scalar_preparation/{f}' for f in ('preparation.py','verify_preparation.py','test_preparation.py','README.md')}
    require(set(value['source_pins'])==expected,'source inventory')
    for f,digest in value['source_pins'].items():require(sha256((root/f).read_bytes()).hexdigest()==digest,'source bytes')
    return {'verified':True,'preparation_error_upper':str(intervention),'pulse_ports':16,
            'resolved_steps':comp['resolved_steps'],'parent_receipts_replayed':False,
            'operator_lower_and_symmetry_reconstructed':True}

if __name__=='__main__':print(json.dumps(check(),indent=2))
