"""Exact action-measure reconstruction and deterministic density controls.

The authenticated scalar family supplies its response law, preparation,
protected records and model clock. The reconstruction does not consume the
parent mass vector. Metric identification requires the separately specified
conformal action model used only by the geometric control below.
"""
from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EVIDENCE = ROOT/'evidence/source_density_geometry'
PARENT = 'code/source_scalar_execution/source_scalar_execution_receipt.json'
PARENT_SHA = '6981fcd4e13fbc33f371224228341a995f34b4733aa70c2ec0b3f6ee5975c7af'
ALGEBRA = ROOT/'code/source_scalar_execution/scalar_execution_algebra.py'
spec = importlib.util.spec_from_file_location('_density_qphi', ALGEBRA)
algebra = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = algebra
spec.loader.exec_module(algebra)
Q, parse, bounds = algebra.Q, algebra.parse, algebra.rational_bounds
TAU = Q(F(-1,35),F(2,35))


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)+'\n').encode('ascii')


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def load_parent():
    data = (ROOT/PARENT).read_bytes()
    if hashlib.sha256(data).hexdigest() != PARENT_SHA:
        raise ValueError('Changed authenticated scalar parent')
    return json.loads(data)


def decode_rows(rows):
    return [{int(j): parse(value) for j,value in row} for row in rows]


def encode_rows(rows):
    return [[[j,value.encode()] for j,value in sorted(row.items())] for row in rows]


def reconstruct(rows, root=0):
    """Infer the normalized positive diagonal symmetrizer from A alone."""
    n=len(rows)
    if not n or not 0 <= root < n:
        raise ValueError('Nonempty finite operator and root required')
    weights=[None]*n;weights[root]=Q(1)
    parents=[None]*n;depths=[0]*n;queue=deque([root])
    for i,row in enumerate(rows):
        if i not in row or any(type(j) is not int or not 0 <= j < n for j in row):
            raise ValueError('Invalid operator support')
        for j,a in row.items():
            if j != i and (a.sign() >= 0 or i not in rows[j] or rows[j][i].sign() >= 0):
                raise ValueError('Off-diagonal support must be paired and strictly negative')
    while queue:
        i=queue.popleft()
        for j,a in sorted(rows[i].items()):
            if i==j:continue
            proposal=weights[i]*a/rows[j][i]
            if weights[j] is None:
                weights[j]=proposal;parents[j]=i;depths[j]=depths[i]+1;queue.append(j)
            elif weights[j] != proposal:
                raise ValueError('Inconsistent edge-ratio cycle')
    if any(w is None for w in weights):
        raise ValueError('Disconnected response support: relative component scales are free')
    total=sum(weights,Q())
    if any(w.sign() <= 0 for w in weights):raise ValueError('Nonpositive inferred action weight')
    edges=sum(len(row)-1 for row in rows)//2
    return [w/total for w in weights], {'root':root,'parents':parents,'depths':depths,
        'undirected_edges':edges,'independent_cycle_edges':edges-n+1,'max_tree_depth':max(depths)}


def floor_q(value):
    lo,hi=map(F,bounds(value))
    k=lo.numerator//lo.denominator
    while (value-(k+1)).sign() >= 0:k+=1
    while (value-k).sign() < 0:k-=1
    return k


def interval(trace, lower=(1,21), upper=(4,21)):
    ancestors={};events={}
    for event in trace['events']:
        key=tuple(event['id']);past=set()
        for read in event['reads']:
            parent=tuple(read[2]);past.add(parent);past.update(ancestors[parent])
        ancestors[key]=past;events[key]=event
    if lower not in ancestors[upper]:raise ValueError('Unordered interval anchors')
    return [events[key] for key in sorted(ancestors[upper]|{upper})
            if key==lower or lower in ancestors[key]]


def emission_panel(events, weights, source_strength, production_rate, budget=4096):
    """Local readback leaves; field-state writes and their order remain fixed.

    The deterministic production law is declared, not inferred: floor(K*r*v)
    copies of each selected cell's immutable local field record, with v=tau*m.
    Every leaf records its same-site writer/version/value and a local write.
    These extra operations are counted; no physical duration is assigned to them.
    """
    cells=[];leaves=[];exact=Q();previous='0'*64
    for cell,event in enumerate(events):
        site=event['id'][1];volume=TAU*weights[site];target=budget*production_rate*volume
        count=floor_q(target);exact+=volume
        if count < 0 or (target-count).sign()<0 or (target-count-1).sign()>=0:
            raise ValueError('Deterministic production discrepancy')
        cells.append({'field_event':event['id'],'field_hash':event['hash'],'site':site,
                      'action_cell_measure_Qphi':volume.encode(),'record_count':count})
        for version in range(1,count+1):
            identity=[source_strength,production_rate,cell,version]
            material={'id':identity,'site':site,'read':event['write'],
                      'write':[128+cell,version,identity,event['write'][3]],'ledger_parent':previous}
            previous=digest(material);leaves.append({**material,'hash':previous})
    estimate=F(len(leaves),budget*production_rate);error=F(len(events),budget*production_rate)
    if (exact-estimate).sign()<0 or (exact-estimate-error).sign()>=0:
        raise ValueError('Finite calibrated-volume enclosure')
    return {'source_strength':source_strength,'production_rate':production_rate,'budget':budget,
        'cells':cells,'observer_records':leaves,'final_observer_hash':previous,
        'raw_record_count':len(leaves),'calibrated_action_measure':str(estimate),
        'exact_action_measure_Qphi':exact.encode(),'one_sided_error_bound':str(error),
        'extra_local_read_write_operations':len(leaves)}


def conformal_operator(rows, base_mass, profile):
    """Declared conformal scalar action: kinetic u, potential u^2, u=Omega^2.

    Undirected gradient terms use the symmetric endpoint average of u;
    Dirichlet boundary terms use the interior endpoint. At u=1 this is the
    original action exactly. This is a control family, not source selection.
    """
    result=[]
    for i,row in enumerate(rows):
        mass=profile[i]*base_mass[i];new={};gradient=Q()
        boundary=base_mass[i]*row[i]-base_mass[i]
        for j,a in row.items():
            if i==j:continue
            conductance=-base_mass[i]*a
            boundary-=conductance
            changed=(profile[i]+profile[j])*conductance/2
            new[j]=-changed/mass;gradient+=changed
        if boundary.sign()<0:raise ValueError('Negative Dirichlet boundary energy')
        new[i]=(gradient+profile[i]*boundary+profile[i]**2*base_mass[i])/mass
        result.append(new)
    return result


def tensor_control(parent, rows, mass):
    """Directional variations of the very same finite scalar action.

    An x/y permutation commutes with the cubic source operator. It preserves
    action energy, action volume and scalar dilation response but interchanges
    directional responses. These are exact discrete-action derivatives, not a physical
    stress-source identification or an additional native-history claim.
    """
    site=[address['site'] for address in parent['source_addresses']]
    index={tuple(p):i for i,p in enumerate(site)}
    permutation=[index[(p[1],p[0],p[2])] for p in site]
    if any(mass[i]!=mass[permutation[i]] or
           {permutation[j]:a for j,a in rows[i].items()}!=rows[permutation[i]] for i in range(64)):
        raise ValueError('Coordinate permutation is not an action symmetry')
    q=list(map(parse,parent['traces']['ascending_intervention']['layers'][8]))
    previous=list(map(parse,parent['traces']['ascending_intervention']['layers'][7]))
    points=[Q()]+[parse(parent['source_addresses'][16*i]['coordinate_Qphi'][0]) for i in range(4)]+[Q(1)]
    gaps=[b-a for a,b in zip(points,points[1:])];widths=[(a+b)/2 for a,b in zip(gaps,gaps[1:])]
    controls=[]
    for label,field,old in [('x_preparation',q,previous),
            ('y_preparation',[q[p] for p in permutation],[previous[p] for p in permutation])]:
        gradient=[Q(),Q(),Q()]
        for i,row in enumerate(rows):
            for j,a in row.items():
                if i<j:
                    axes=[k for k in range(3) if site[i][k]!=site[j][k]]
                    if len(axes)!=1:raise ValueError('Nonaxial gradient edge')
                    gradient[axes[0]]+=-mass[i]*a*(field[i]-field[j])**2/2
            for axis,k in enumerate(site[i]):
                if k in (0,3):
                    gap=gaps[0] if k==0 else gaps[-1]
                    gradient[axis]+=mass[i]*field[i]**2/(2*widths[k]*gap)
        potential=sum((m*x*x for m,x in zip(mass,field)),Q())/2
        total=sum(gradient,Q())+potential
        direct=sum((mass[i]*field[i]*sum((a*field[j] for j,a in row.items()),Q()) for i,row in enumerate(rows)),Q())/2
        if total!=direct:raise ValueError('Directional energy omits an action term')
        derivative=[total-2*g for g in gradient]
        kinetic=sum((m*(a-b)**2 for m,a,b in zip(mass,field,old)),Q())/(2*TAU)
        lagrangian=kinetic-TAU*total
        action_derivative=[kinetic-TAU*value for value in derivative]
        controls.append({'label':label,'field_Qphi':[x.encode() for x in field],
            'previous_field_Qphi':[x.encode() for x in old],
            'directional_gradient_energy_Qphi':[x.encode() for x in gradient],
            'mass_potential_energy_Qphi':potential.encode(),'total_potential_energy_Qphi':total.encode(),
            'log_axis_stretch_potential_derivative_Qphi':[x.encode() for x in derivative],
            'discrete_kinetic_action_Qphi':kinetic.encode(),'discrete_lagrangian_Qphi':lagrangian.encode(),
            'log_axis_stretch_action_derivative_Qphi':[x.encode() for x in action_derivative],
            'isotropic_action_dilation_derivative_Qphi':sum(action_derivative,Q()).encode()})
    if controls[0]['log_axis_stretch_action_derivative_Qphi']==controls[1]['log_axis_stretch_action_derivative_Qphi']:
        raise ValueError('Tensor distinction is vacuous')
    return {'source_layer':8,'previous_source_layer':7,'potential_endpoint':'right',
            'x_y_site_permutation':permutation,'controls':controls,
            'interpretation':'Fixed-configuration discrete-action derivatives including kinetic and potential terms; no physical stress coupling selected.'}


def build():
    parent=load_parent();rows=decode_rows(parent['action_rows_Qphi'])
    weights,tree=reconstruct(rows)
    mass=list(map(parse,parent['mass_Qphi']));mass_total=sum(mass,Q())
    if weights != [m/mass_total for m in mass]:raise ValueError('Independent parent action weights disagree')
    parent_module_path=ROOT/'code/source_scalar_execution/source_scalar_execution.py'
    sys.path.insert(0,str(parent_module_path.parent))
    try:
        spec=importlib.util.spec_from_file_location('_density_scalar_driver',parent_module_path)
        driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
        _,_,_,driver_rows,velocity,_,_=driver.model()
        traces={str(s):driver.execute(driver_rows,[F(s,2)*v for v in velocity],list(range(64))) for s in (1,2)}
    finally:sys.path.remove(str(parent_module_path.parent))
    if traces['2'] != parent['traces']['ascending_intervention']:
        raise ValueError('Full-strength replay differs from authenticated parent')
    intervals={s:interval(trace) for s,trace in traces.items()}
    if [[e['id'] for e in interval_] for interval_ in intervals.values()][0] != [[e['id'] for e in interval_] for interval_ in intervals.values()][1]:
        raise ValueError('Source amplitude changed field ancestry')
    panels=[emission_panel(intervals[str(s)],weights,s,r) for s in (1,2) for r in (1,2)]
    x=[parse(address['coordinate_Qphi'][0]) for address in parent['source_addresses']]
    controls=[]
    for slope in (F(0),F(1,2)):
        profile=[1+slope*z for z in x]
        changed=conformal_operator(rows,mass,profile)
        recovered,_=reconstruct(changed)
        # An interior row has no Dirichlet boundary contribution. With the
        # declared unit mass potential its row sum is u, fixing the common
        # profile factor from the full operator rather than an imported u0.
        reference=21;calibration=sum(changed[reference].values(),Q())
        u=[(recovered[i]/weights[i])/(recovered[reference]/weights[reference])*calibration for i in range(64)]
        if u != profile:raise ValueError('Conformal profile reconstruction')
        volume=[weights[i]*value**2 for i,value in enumerate(u)]
        kinetic=[weights[i]*value for i,value in enumerate(u)]
        volume_sum=sum(volume,Q());kinetic_sum=sum(kinetic,Q())
        fitted_slope=(u[48]-u[0])/(x[48]-x[0]);fitted_intercept=u[0]-fitted_slope*x[0]
        if any(value!=fitted_intercept+fitted_slope*coordinate for value,coordinate in zip(u,x)):
            raise ValueError('Nonaffine conformal control')
        controls.append({'supplied_profile_slope':str(slope),'supplied_mass_potential_squared':'1',
            'calibration_site':reference,'interior_rowsum_Qphi':calibration.encode(),
            'identified_affine_intercept_Qphi':fitted_intercept.encode(),'identified_affine_slope_Qphi':fitted_slope.encode(),
            'operator_rows_Qphi':encode_rows(changed),'normalized_action_weights_Qphi':[m.encode() for m in recovered],
            'reconstructed_u_Qphi':[v.encode() for v in u],
            'normalized_metric_four_volume_Qphi':[(v/volume_sum).encode() for v in volume],
            'normalized_kinetic_measure_Qphi':[(v/kinetic_sum).encode() for v in kinetic],
            'scalar_curvature_Qphi':[(F(3,2)*fitted_slope**2/(v**3)).encode() for v in u],
            'raw_support_order_unchanged':all(set(a)==set(b) for a,b in zip(rows,changed))})
    epsilon=F(1,1000);ratio=(1+epsilon)/(1-epsilon);depth=tree['max_tree_depth']
    initial_energy=sum((m*parse(v)**2 for m,v in zip(mass,parent['initial_velocity_Qphi'])),Q())/2
    return {'schema':'oph.source-action-density-identification.v1',
        'parent':{'path':PARENT,'sha256':PARENT_SHA},
        'operator_rows_sha256':digest(parent['action_rows_Qphi']),
        'reconstruction':{'normalized_mass_Qphi':[m.encode() for m in weights],
            'parent_mass_total_Qphi':mass_total.encode(),'spanning_tree':tree,
            'parent_mass_used_by_reconstruction':False,
            'finite_edge_error_control':{'assumed_relative_error':str(epsilon),
                'edge_ratio_distortion_upper':str(ratio),'root_relative_distortion_upper':str(ratio**depth),
                'normalized_weight_distortion_upper':str(ratio**(2*depth))}},
        'source_controls':[{'strength':s,'amplitude_multiplier':str(F(s,2)),
            'trace_file':f'source_strength_{s}.json','trace_sha256':digest(traces[str(s)]),
            'field_writes':len(traces[str(s)]['events']),'initial_energy_Qphi':(F(s*s,4)*initial_energy).encode(),
            'action_and_recovered_measure_unchanged':True} for s in (1,2)],
        'interval':{'lower':[1,21],'upper':[4,21],'field_event_ids':[e['id'] for e in intervals['1']],
            'clock_Qphi':TAU.encode(),'measure':'model tick times normalized action mass; not an identified physical metric volume'},
        'emission_panels':panels,'conformal_action_controls':controls,
        'tensor_action_control':tensor_control(parent,rows,mass),
        'scope':{'source_population_action_read_grammar_and_clock_supplied':True,
            'same_family_local_observers_check_actual_immutable_values':True,
            'observer_operations_charged_but_physical_durations_not_identified':True,
            'aggregate_certificate_is_not_single_observer_remote_record_transport':True,
            'deterministic_production_law_selected_by_axioms':False,
            'global_action_scale_identified_by_reciprocity':False,
            'kinetic_density_identified_as_metric_volume_without_constitutive_inputs':False,
            'conformal_controls_are_explicitly_supplied_not_sourced_geometries':True,
            'curvature_from_unweighted_order_counts':False,
            'physical_stress_coupling_or_Einstein_response_derived':False}}, traces


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--check',action='store_true');args=ap.parse_args()
    packet,traces=build();outputs={'receipt.json':packet,**{f'source_strength_{s}.json':trace for s,trace in traces.items()}}
    EVIDENCE.mkdir(parents=True,exist_ok=True)
    for name,value in outputs.items():
        data=canonical(value);path=EVIDENCE/name
        if args.check:
            if path.read_bytes()!=data:raise ValueError('Changed generated evidence: '+name)
        else:path.write_bytes(data)
    print('SOURCE_ACTION_DENSITY_REPRODUCED')


if __name__=='__main__':main()
