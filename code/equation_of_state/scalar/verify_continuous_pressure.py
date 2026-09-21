"""Independent exact sqrt(5)-arithmetic audit of continuous scalar pressure.

Reuses the parent's independent arithmetic and event verifier, but imports
neither continuous-pressure producer nor scalar-pressure producer. Rebuilds
both Taylor polynomials by Horner evaluation and all quadratic observables.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
RER = HERE.parents[2]
PARENT = 'code/source_scalar_execution/source_scalar_execution_receipt.json'
sys.path.insert(0, str(RER/'code/source_scalar_execution'))
import verify_source_scalar_execution as parent

R, parse = parent.R, parent.parse
require = parent.require
EXPECTED_SCOPE = {'same_supplied_finite_action': True, 'all_22_model_times_retained': True,
                  'continuous_pressure_rigorously_enclosed': True, 'native_repair_eos': False,
                  'equilibrium_thermodynamics': False, 'new_observer_history': False,
                  'recorded_split_values_identified_with_exact_continuous_values': False,
                  'physical_clock_calibration': False, 'continuum_limit': False}


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def remainder():
    q = F(4)*F(3,2)**163*614**81/factorial(163)
    v = F(4)*F(3,2)**162*614**81/factorial(162)
    p = F(1,2)*(10+v)*v + F(308,3)*(14+q)*q
    require(0 < q < v < F(1,10**28), 'remainder_budget')
    return q,v,p


def preflight(report):
    require(type(report) is dict and report.get('schema') == 'oph.continuous_scalar_pressure.v1', 'schema')
    require(report['scope'] == EXPECTED_SCOPE, 'scope_promotion')
    require(all(report['scope'][k] is v for k,v in EXPECTED_SCOPE.items()), 'scope_type')
    require(report['polynomial'] == {'q_degree':161,'v_degree':160,'terms':81,
                                     'all_operator_modes_retained':64}, 'all_mode_polynomial')
    dq,dv,dp = remainder()
    expected = {'operator_lower':'1','operator_upper':'614','gradient_operator_upper':'613',
                'initial_velocity_norm_upper':'4','time_upper':'3/2',
                'q_polynomial_norm_upper':'7','v_polynomial_norm_upper':'5',
                'q_tail':str(dq),'v_tail':str(dv),'pressure_tail':str(dp),
                'pressure_tail_formula':'5*dv+dv^2/2+(616/3)*(7*dq+dq^2/2)'}
    require(report['bounds'] == expected, 'spectral_or_quadratic_error_bound')
    require(len(report['samples']) == 22, 'missing_time')
    require(report['assumptions'] == 'supplied q5 Dirichlet scalar action, preparation, model time, V=1 and fixed-original-q,p metric work', 'assumptions')
    require(report['baseline_control'] == {'continuous_q_v_and_pressure':'identically zero',
                                          'w':None,'all_22_times_retained_in_parent':True}, 'zero_control')
    require(report['schedule_control'] == 'both independently replayed intervention schedules give identical layer states', 'schedule_control')
    for j,row in enumerate(report['samples']):
        require(row['step'] == j, 'time_order')
        require(row['pressure_tail_bound'] == str(F() if j==0 else dp), 'sample_tail_bound')
        for key in ('continuous_pressure_interval','continuous_w_interval',
                    'split_minus_continuous_pressure_interval'):
            lo,hi = map(parent.rational,row[key])
            require(lo <= hi, 'interval_order')


def bounds(value,radius=F()):
    return parent.widen(parent.value_bounds(value),radius)


def full_check(report):
    preflight(report)
    require(set(report) == {'schema','scope','assumptions','polynomial','bounds',
        'continuous_energy_Qphi','continuous_energy_interval','baseline_control',
        'schedule_control','samples','parent_replay','source_sha256','producer_sha256'},'top_fields')
    paths=(PARENT,'code/source_scalar_execution/source_scalar_execution.py',
        'code/source_scalar_execution/scalar_execution_algebra.py',
        'code/source_scalar_execution/verify_source_scalar_execution.py')
    require(report['source_sha256'] == {p:sha(RER/p) for p in paths},'parent_source_bytes')
    require(report['producer_sha256']==sha(HERE/'continuous_pressure.py'),'producer_bytes')
    packet=parent.load(RER/PARENT)
    replay=parent.verify(packet,root=RER)
    require(report['parent_replay']==replay,'parent_replay')
    mass=[parse(x) for x in packet['mass_Qphi']]
    rows=[[(j,parse(x)) for j,x in row] for row in packet['action_rows_Qphi']]
    velocity=[parse(x) for x in packet['initial_velocity_Qphi']]
    layers=[[parse(x) for x in layer] for layer in packet['traces']['ascending_intervention']['layers']]
    # Independently establish the self-adjoint interval [1,614]. The
    # lower bound follows from the explicit nonnegative edge-square form.
    matrix=[dict(row) for row in rows]
    for i in range(64):
        require(mass[i].sign()>0,'positive_mass')
        sumabs=R()
        rowsum=R(-1)
        for j,c in rows[i]:
            require(mass[i]*c==mass[j]*matrix[j].get(i,R()),'weighted_symmetry')
            if i != j:
                require(c.sign()<=0,'negative_offdiagonal')
            sumabs+=c if c.sign()>=0 else -c
            rowsum+=c
        require(rowsum.sign()>=0,'gradient_edge_square_boundary')
        require((614-sumabs).sign()>=0,'gershgorin_upper')
    norm0=parent.inner(mass,velocity,velocity)
    require((16-norm0).sign()>=0 and (R(F(3,2))-21*parent.TAU).sign()>0,
            'initial_norm_and_time')
    energy=norm0/2
    require(report['continuous_energy_Qphi']==energy.encode(),'exact_energy')
    require(report['continuous_energy_interval']==[str(x) for x in parent.value_bounds(energy)],'energy_interval')
    energy_lower=parent.value_bounds(energy)[0]
    require(energy_lower>0,'positive_energy')
    powers=[velocity]
    for k in range(80):
        powers.append(parent.apply(rows,powers[-1]))
    _,_,dp=remainder()
    for step,row in enumerate(report['samples']):
        q,v=[],[]
        square=F(step*step,245)
        for i in range(64):
            qq,vv=R(),R()
            for k in range(80,-1,-1):
                qq=qq*square+F((-1)**k,factorial(2*k+1))*powers[k][i]
                vv=vv*square+F((-1)**k,factorial(2*k))*powers[k][i]
            q.append(qq*(step*parent.TAU))
            v.append(vv)
        require((49-parent.inner(mass,q,q)).sign()>=0,'q_polynomial_norm')
        require((25-parent.inner(mass,v,v)).sign()>=0,'v_polynomial_norm')
        aq=parent.apply(rows,q)
        K=parent.inner(mass,v,v)/2
        U=parent.inner(mass,q,q)/2
        G=parent.inner(mass,q,aq)/2-U
        pressure=K-G/3-U
        previous,actual=layers[step:step+2]
        force=parent.apply(rows,actual)
        actual_v=[(a-b)/parent.TAU-parent.TAU*f/2 for a,b,f in zip(actual,previous,force)]
        actual_k=parent.inner(mass,actual_v,actual_v)/2
        actual_u=parent.inner(mass,actual,actual)/2
        actual_g=parent.inner(mass,actual,force)/2-actual_u
        actual_p=actual_k-actual_g/3-actual_u
        radius=F() if step==0 else dp
        expected={
            'step':step,'model_time_Qphi':(step*parent.TAU).encode(),
            'polynomial_q_Qphi':[x.encode() for x in q],
            'polynomial_v_Qphi':[x.encode() for x in v],
            'polynomial_K_Qphi':K.encode(),'polynomial_G_Qphi':G.encode(),
            'polynomial_U_Qphi':U.encode(),'polynomial_pressure_Qphi':pressure.encode(),
            'continuous_pressure_interval':bounds(pressure,radius),
            'continuous_w_interval':bounds(pressure/energy,radius/energy_lower),
            'split_pressure_Qphi':actual_p.encode(),
            'split_minus_continuous_pressure_interval':bounds(actual_p-pressure,radius),
            'pressure_tail_bound':str(radius)}
        require(row==expected,f'continuous_pressure_replay:{step}')


def verify(report):
    try:
        full_check(report)
        passed,errors=True,[]
    except (ValueError,KeyError,TypeError,ZeroDivisionError,OSError,OverflowError,RecursionError) as e:
        passed,errors=False,[str(e)]
    return {'schema':'oph.continuous_scalar_pressure.verification.v1',
            'verdict':'PASS' if passed else 'FAIL','receipt':passed,
            'independent_implementation':True,'producer_imported':False,
            'all_22_times_replayed':passed,'exact_q_phi_and_sqrt5_implementations_agree':passed,
            'rigorous_pressure_enclosures':passed,'native_eos':False,
            'parent_events_replayed':5888 if passed else 0,'errors':errors}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path',nargs='?',type=Path,default=HERE/'continuous_pressure_receipt.json')
    parser.add_argument('--out',type=Path,default=HERE/'continuous_pressure_verification.json')
    args=parser.parse_args()
    result=verify(json.loads(args.path.read_text()))
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['receipt'] else 1)
