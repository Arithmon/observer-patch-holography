"""Exact biased-core renewal and complete flagged quantum reset service.

The optimized quantity is the number of classical core-label resets in the
already declared sharp-code family.  Quantum recovery remains active after
every attempt.  The theorem does not identify apparatus work with energy.
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import math
import numpy as np
import sympy as sp
from ..moving_defect import transport as t


@lru_cache(None)
def exact_chain(fair=False):
    if type(fair) is not bool:
        raise ValueError('explicit fair-control flag required')
    p=t.local.Parameters() if fair else t.parameters()
    full=t.kernel(p,'C')
    if full!=t.kernel(p,'E'):
        raise ArithmeticError('tagged C/E core kernels differ')
    q=tuple(tuple(row[t.CORES.index(c)] for c in t.LIVE) for row in full)
    hazard=tuple(1-sum(row) for row in q)
    matrix=sp.eye(60)-sp.Matrix(q)
    h=tuple(F(v) for v in matrix.inv()*sp.ones(60,1))
    if any(h[i]!=1+sum(q[i][j]*h[j] for j in range(60)) for i in range(60)):
        raise ArithmeticError('first-hit potential fails its original equations')
    maximum=max(h)
    best=tuple(t.LIVE[i] for i,value in enumerate(h) if value==maximum)
    return p,full,q,hazard,h,maximum,best


def renewal_certificate(fair=False):
    p,full,q,hazard,h,maximum,best=exact_chain(fair)
    target=best[0];j=t.LIVE.index(target)
    # The occupation measure before the next boundary hit; no simulated fit.
    seed=sp.zeros(60,1);seed[j]=1
    visit=(sp.eye(60)-sp.Matrix(q).T).inv()*seed
    stationary=tuple(F(v)/maximum for v in visit)
    controlled=tuple(tuple(q[i][k]+hazard[i]*F(k==j) for k in range(60)) for i in range(60))
    if sum(stationary)!=1 or any(sum(stationary[i]*controlled[i][k] for i in range(60))!=stationary[k] for k in range(60)):
        raise ArithmeticError('continued full-core stationary law failed')
    rate=sum(stationary[i]*hazard[i] for i in range(60))
    if rate!=1/maximum:
        raise ArithmeticError('renewal/reset-rate equality failed')
    eps=min(hazard[i]+sum(q[i][k]*hazard[k] for k in range(60)) for i in range(60))
    if eps<=0:
        raise ArithmeticError('finite renewal tail not established')
    return {'parameters':p.description(),'target_core':list(target),
            'all_maximizing_cores':[list(c) for c in best],
            'core_order':[list(c) for c in t.LIVE],
            'expected_first_boundary_attempts':[str(v) for v in h],
            'maximum_mean_excursion_attempts':str(maximum),
            'minimum_two_attempt_boundary_probability':str(eps),
            'continued_kernel':[[str(v) for v in row] for row in controlled],
            'stationary_core_probabilities':[str(v) for v in stationary],
            'stationary_resets_per_native_attempt':str(rate),
            'asymptotic_reset_rate_globally_minimal_in_declared_code_family':True,
            'finite_horizon_universal_lower_bound':'E[N_n] >= max(0,(n-h(initial_core))/max(h))',
            'proof':'Extend h by zero on forbidden cores. Every native attempt lowers its conditional expectation by one. A counted reset raises h by at most max(h); telescope. Reset only at the boundary to an h-maximizer attains renewal mean max(h).',
            'domain':'The fixed60 supplied codes; uniform complete seam attempts; fixed endpoint weights; arbitrary history-dependent C/E choices and adaptive/randomized legal-core resets; every attempt restored to a legal code, no halting/postselection to evade attempts.',
            'optimizes_source_pulses_or_physical_work':False}


def bellman_certificate(horizon=64):
    if type(horizon) is not int or not 0<=horizon<=128:
        raise ValueError('integer horizon0..128 required')
    _,_,q,hazard,_,_,_=exact_chain()
    sparse=[[(j,v) for j,v in enumerate(row) if v] for row in q]
    values=[F(0)]*60
    rows=[]
    for n in range(1,horizon+1):
        boundary_cost=1+min(values)
        continuation=[sum(v*values[j] for j,v in row)+hazard[i]*boundary_cost
                      for i,row in enumerate(sparse)]
        immediate_reset=1+min(continuation)
        new=[min(v,immediate_reset) for v in continuation]
        rows.append({'attempts':n,'optimal_value_from_initial':str(new[t.LIVE.index(t.INITIAL)]),
                     'minimizing_initial_cores':[list(t.LIVE[i]) for i,v in enumerate(new) if v==min(new)],
                     'voluntary_reset_better_cores':[list(t.LIVE[i]) for i,v in enumerate(continuation) if v>immediate_reset],
                     'maximum_continuation_value_spread':str(max(continuation)-min(continuation))})
        values=new
    return {'horizon':horizon,'rows':rows,'terminal_values':[str(v) for v in values],
            'recurrence':'B_n(c)=sum_live Q(c,d)V_(n-1)(d)+r(c)[1+min V_(n-1)]; V_n(c)=min(B_n(c),1+min B_n)',
            'exact_dynamic_program_for_each_declared_finite_horizon':True,
            'same_boundary_target_proved_optimal_for_all_horizons':False}


def recover_boundary(center,law,seam,coin,raw,state,*,target=t.INITIAL):
    """Decode the retained forbidden-core image and re-encode its qubit.

    A full512-address emitted unitary acts on the unnormalized branch. No
    weight is divided out. Other branches already have their ordinary
    transported-core recovery and remain unchanged by this extra service.
    """
    t._core(target)
    value=t.routing._state(state)
    images,moved,_=t.branch(center,law,seam,coin,raw)
    if not images[0] or moved not in t.STOP:
        return value
    value=t.record_code._apply_permutation(t.branch_permutation(center,law,seam,coin,raw),value,True)
    value=t.compiler.execute(t.fixed._mixing_word(len(images[0]),True),9,value)
    return t._encode_register(target,value)


def boundary_resources(center,law,seam,coin,raw,*,target=t.INITIAL):
    t._core(target)
    images,moved,_=t.branch(center,law,seam,coin,raw)
    if not images[0] or moved not in t.STOP:
        return dict.fromkeys(t.routing.COUNT_KEYS,0)
    return t.fixed._counts(
        t.record_code._permutation_resources(t.branch_permutation(center,law,seam,coin,raw)),
        t.fixed._mixing_resources(len(images[0]),True),
        t.fixed._mixing_resources(len(t.supports(target)[0])),
        t.record_code._permutation_resources(t.encoder_permutation(target)))


def source_branches(center,law,seam,state,*,target=t.INITIAL):
    rows=t.source_branches(center,law,seam,state)
    for row in rows:
        row['native_core_after']=row['core_after']
        row['reset_performed']=row['status']=='opposite_core_stop'
        row['amplitude']=recover_boundary(center,law,seam,row['coin'],row['raw'],row['amplitude'],target=target)
        if row['reset_performed']:
            row['core_after']=target
            row['status']='reset_continue'
    return rows


def exact_service_certificate():
    census=Counter();schedules=set();costs=[]
    for center in t.LIVE:
        for law in ('C','E'):
            for seam in t.SEAMS:
                total=F()
                for coin in (0,1):
                    for raw in range(1<<t.circuit.RECORD_BITS[law]):
                        images,moved,abort=t.branch(center,law,seam,coin,raw)
                        total+=t.branch_probability(center,law,seam,coin,raw)
                        census['all_raw_slots']+=1
                        if not images[0]:
                            census['raw_abort' if abort else 'valid_zero']+=1
                        elif moved in t.STOP:
                            census['boundary_reset_positive']+=1
                            perm=t.branch_permutation(center,law,seam,coin,raw)
                            if perm not in schedules:
                                executed=t.record_code._integer_schedule_action(t.routing.permutation_edges(perm,inverse=True))
                                if executed!=tuple(perm.index(i) for i in range(512)):
                                    raise ArithmeticError('reset decode Gray schedule differs')
                                schedules.add(perm)
                            costs.append(boundary_resources(center,law,seam,coin,raw)['native_pulses'])
                        else:
                            census['ordinary_continue_positive']+=1
                if total!=1:
                    raise ArithmeticError('continued instrument loses complete raw weight')
    return {'census':dict(census),'distinct_new_inverse_schedules_checked':len(schedules),
            'addresses_checked_per_schedule':512,'maximum_extra_boundary_reset_pulses':max(costs),
            'ideal_flag_identity':'R_reset K_flag V_core = sqrt(p_flag) V_(0,1)',
            'arbitrary_logical_reference_supported':True,
            'reset_count_is_not_recovery_count':True,
            'fixed_format_record_bits_per_attempt':{'C':25,'E':26},
            'record_format':'law1,seam5,coin1,raw2or3,oldcore7,nativecore7,reset1,halt1; target(0,1) is fixed',
            'archive_growth_in_bits_is_linear':True}


def resource_summary():
    """Charge every old recovery and every added boundary service anew."""
    old=t.resource_certificate()
    stationary=list(map(F,renewal_certificate()['stationary_core_probabilities']))
    baseline=F();extra=F();peak=0
    for i,center in enumerate(t.LIVE):
        baseline+=stationary[i]*sum(F(old['rows'][i]['laws'][law]['uniform_seam_expected_source_pulses_per_attempt']) for law in ('C','E'))/2
        for law in ('C','E'):
            for seam in t.SEAMS:
                for coin in(0,1):
                    for raw in range(1<<t.circuit.RECORD_BITS[law]):
                        weight=t.branch_probability(center,law,seam,coin,raw)
                        if not weight:
                            continue
                        added=boundary_resources(center,law,seam,coin,raw)['native_pulses']
                        extra+=stationary[i]*weight*added/60
                        peak=max(peak,t.fixed.repair_resources(law,seam)['native_pulses']+
                                 t.recovery_resources(center,law,seam,coin,raw)['native_pulses']+added)
    return {'stationary_equal_C_E_source_pulses_per_attempt':str(baseline+extra),
            'stationary_existing_native_and_ordinary_recovery_pulses_per_attempt':str(baseline),
            'stationary_added_boundary_service_pulses_per_attempt':str(extra),
            'peak_source_pulses_per_positive_complete_attempt':peak,
            'initial_encoder':old['initial_encoder'],
            'coherent_recovery_after_each_attempt_supplied':True,
            'proper_code_hardware_failure_certificate_constructed':False,
            'storage_and_classical_control_physical_cost_established':False,
            'units':'source pulses per native attempt, not physical energy or time'}


def clock_boundary_certificate():
    continued=renewal_certificate()
    pi=list(map(F,continued['stationary_core_probabilities']))
    clocks=t.clock_compatibility()['rows']
    coefficients=[F(row['modular_code_leakage_sine_squared_coefficient']) for row in clocks]
    incompatible=[i for i,row in enumerate(clocks) if not row['code_invariant_under_all_modular_parameters']]
    best=t.LIVE.index(t.INITIAL)
    kernel=[[F(v) for v in row] for row in continued['continued_kernel']]
    return {'best_reset_core_is_fixed_clock_invariant':clocks[best]['code_invariant_under_all_modular_parameters'],
            'stationary_mass_on_clock_incompatible_codes':str(sum(pi[i] for i in incompatible)),
            'stationary_mean_clock_leakage_sine_squared_coefficient':str(sum(p*c for p,c in zip(pi,coefficients,strict=True))),
            'one_attempt_probability_to_leave_clock_invariant_codes_from_reset_core':str(sum(kernel[best][i] for i in incompatible)),
            'modular_leakage_parameter':'sin^2(s log(25/8)/2)',
            'joint_preserved_moving_clock_subspace_obtained':False}


def source_diagnostics():
    logical=np.array([[1,1j],[1j,1]],complex)/math.sqrt(2)
    state=t.encode(t.INITIAL,logical)
    weight=F(1);events=list(t.witness_path())+[t.witness_path()[0]]
    result=[];worst=0.;mass_error=0.
    for center,law,seam,selected_coin,selected_raw in events:
        rows=source_branches(center,law,seam,state)
        mass_error=max(mass_error,abs(sum(float(np.linalg.norm(r['amplitude'])**2) for r in rows)-float(np.linalg.norm(state)**2)))
        for row in rows:
            expected=(np.zeros_like(state) if not row['probability'] else
                      math.sqrt(float(weight*row['probability']))*(t.ideal_code(row['core_after'])@logical))
            worst=max(worst,float(np.linalg.norm(row['amplitude']-expected,2)))
        chosen=next(r for r in rows if(r['coin'],r['raw'])==(selected_coin,selected_raw))
        state=chosen['amplitude'];weight*=chosen['probability']
        result.append({'core_before':list(center),'law':law,'seam':list(seam),
                       'coin':selected_coin,'raw':selected_raw,
                       'native_core_after':list(chosen['native_core_after']),
                       'core_after':list(chosen['core_after']),'status':chosen['status'],
                       'complete_raw_outputs_retained':len(rows),
                       'unnormalized_prefix_weight_given_selected_seams':str(weight),
                       'extra_reset_resources':boundary_resources(center,law,seam,selected_coin,selected_raw)})
    for number in (worst,mass_error):
        if not math.isfinite(number) or not 0<=number<1e-8:
            raise ArithmeticError('reset continuation failed fixed numerical source gate')
    return {'path':result,'maximum_branch_operator_error':worst,'maximum_total_mass_error':mass_error}


def admitted_source_report():
    result=source_diagnostics()
    return {'path':result['path'],
            'maximum_branch_operator_error_checked_interval':['0','1/100000000'],
            'maximum_total_mass_error_checked_interval':['0','1/100000000'],
            'floating_diagnostics_are_outward_certified':False,
            'complete_branches_at_each_declared_prefix_executed':True,
            'all43200_floating_branch_words_executed':False}
