"""Independent certificate checks in the parent's sqrt(5) arithmetic.

No new producer or shared new numerical helper is imported. Each source trace
is independently replayed by the original source verifier. Positivity, every
reciprocity equation, connectivity and normalization certify the recovered
measure; the spanning-tree algorithm is not rerun here.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
EVIDENCE=ROOT/'evidence/source_density_geometry'
PARENT='code/source_scalar_execution/source_scalar_execution_receipt.json'
PARENT_SHA='6981fcd4e13fbc33f371224228341a995f34b4733aa70c2ec0b3f6ee5975c7af'
VERIFIER='code/source_scalar_execution/verify_source_scalar_execution.py'
VERIFIER_SHA='c03892572732bbf7d019e721d19b6a249b445e2dc30f54d3c237be938cf88d42'


def need(condition,message):
    if not condition:raise ValueError(message)


def canonical(value):
    return (json.dumps(value,sort_keys=True,indent=2,ensure_ascii=True,allow_nan=False)+'\n').encode('ascii')


def digest(value):return hashlib.sha256(canonical(value)).hexdigest()


def context():
    raw=(ROOT/PARENT).read_bytes();need(hashlib.sha256(raw).hexdigest()==PARENT_SHA,'parent pin')
    code=(ROOT/VERIFIER).read_bytes();need(hashlib.sha256(code).hexdigest()==VERIFIER_SHA,'independent parent verifier pin')
    spec=importlib.util.spec_from_file_location('_density_independent_parent',ROOT/VERIFIER)
    v=importlib.util.module_from_spec(spec);sys.modules[spec.name]=v
    exec(compile(code,str(ROOT/VERIFIER),'exec'),v.__dict__)
    return v,json.loads(raw)


def certify_measure(v,rows,weights):
    n=len(rows);need(len(weights)==n and n>0,'measure dimension')
    need(all(w.sign()>0 for w in weights),'positive action measure')
    need(sum(weights,v.R())==v.R(1),'normalized action measure')
    components=[{i} for i in range(n)];edges=0
    for i,row in enumerate(rows):
        need(i in row,'operator diagonal')
        for j,a in row.items():
            if i==j:continue
            need(0<=j<n and i in rows[j] and a.sign()<0 and rows[j][i].sign()<0,'paired negative support')
            need(weights[i]*a==weights[j]*rows[j][i],'action reciprocity')
            if i<j:
                edges+=1
                left=next(group for group in components if i in group)
                right=next(group for group in components if j in group)
                if left is not right:left.update(right);components.remove(right)
    need(len(components)==1,'connected action support')
    return edges


def verify(packet=None, traces=None):
    v,parent=context();R,parse=v.R,v.parse
    packet=json.loads((EVIDENCE/'receipt.json').read_text()) if packet is None else packet
    need(packet['schema']=='oph.source-action-density-identification.v1','schema')
    need(packet['parent']=={'path':PARENT,'sha256':PARENT_SHA},'parent identity')
    mass,action,velocity,detector,addresses=v.model(ROOT)
    need(parent['action_rows_Qphi']==[[[j,a.encode()] for j,a in row] for row in action],'actual source operator')
    need(parent['mass_Qphi']==[m.encode() for m in mass],'actual source quadrature')
    rows=[dict(row) for row in action]
    need(packet['operator_rows_sha256']==digest(parent['action_rows_Qphi']),'operator identity')
    record=packet['reconstruction'];weights=list(map(parse,record['normalized_mass_Qphi']))
    edges=certify_measure(v,rows,weights)
    total=sum(mass,R());need(record['parent_mass_total_Qphi']==total.encode(),'independent total quadrature')
    need(weights==[m/total for m in mass],'independent normalized source quadrature')
    need(record['parent_mass_used_by_reconstruction'] is False,'reconstruction boundary')
    tree=record['spanning_tree'];need(tree['root']==0,'root')
    need(len(tree['parents'])==len(tree['depths'])==64,'tree dimensions')
    for i,(parent_index,depth) in enumerate(zip(tree['parents'],tree['depths'])):
        if i==0:need(parent_index is None and depth==0,'tree root');continue
        need(type(parent_index) is int and 0<=parent_index<64 and parent_index in rows[i],'tree edge')
        need(type(depth) is int and depth==tree['depths'][parent_index]+1,'acyclic rooted depths')
    maximum=max(tree['depths'])
    need(tree['max_tree_depth']==maximum and tree['undirected_edges']==edges and tree['independent_cycle_edges']==edges-63,'tree census')
    error=record['finite_edge_error_control'];epsilon=F(error['assumed_relative_error'])
    need(epsilon==F(1,1000),'edge error budget')
    ratio=(1+epsilon)/(1-epsilon)
    need(error=={'assumed_relative_error':str(epsilon),'edge_ratio_distortion_upper':str(ratio),
        'root_relative_distortion_upper':str(ratio**maximum),
        'normalized_weight_distortion_upper':str(ratio**(2*maximum))},'finite ratio propagation bound')
    need(len(packet['source_controls'])==2,'source strength census')
    controls=packet['source_controls'];interval_events={};replayed=0
    selected_ids=None;all_traces={}
    for row,strength in zip(controls,(1,2),strict=True):
        need(row['strength']==strength and row['amplitude_multiplier']==str(F(strength,2)),'source strengths')
        name=f'source_strength_{strength}.json';need(row['trace_file']==name,'bounded trace filename')
        trace=json.loads((EVIDENCE/name).read_text()) if traces is None else traces[str(strength)]
        all_traces[str(strength)]=trace
        need(row['trace_sha256']==digest(trace),'source trace bytes')
        _,ancestors=v.verify_trace(trace,action,[F(strength,2)*x for x in velocity],list(range(64)))
        replayed+=len(trace['events'])
        need(row['field_writes']==1472 and row['action_and_recovered_measure_unchanged'] is True,'field census and fixed action')
        energy=F(strength*strength,8)*sum((m*x*x for m,x in zip(mass,velocity)),R())
        need(row['initial_energy_Qphi']==energy.encode(),'source energy')
        lower,upper=(1,21),(4,21)
        selected=sorted(key for key in ancestors[upper]|{upper} if key==lower or lower in ancestors[key])
        need(selected_ids is None or selected==selected_ids,'same source field order');selected_ids=selected
        byid={tuple(event['id']):event for event in trace['events']}
        interval_events[strength]=[byid[key] for key in selected]
    need(all_traces['2']==parent['traces']['ascending_intervention'],'full-strength parent history')
    need(packet['interval']=={'lower':[1,21],'upper':[4,21],'field_event_ids':[list(key) for key in selected_ids],
        'clock_Qphi':v.TAU.encode(),'measure':'model tick times normalized action mass; not an identified physical metric volume'},'interval and clock interpretation')
    need(len(packet['emission_panels'])==4,'factorial source/production controls')
    observer_events=0
    for panel,(strength,rate) in zip(packet['emission_panels'],((1,1),(1,2),(2,1),(2,2)),strict=True):
        need(panel['source_strength']==strength and panel['production_rate']==rate and panel['budget']==4096,'emission parameters')
        selected=interval_events[strength];need(len(panel['cells'])==len(selected),'cell census')
        chain='0'*64;offset=0;exact=R()
        for index,(cell,source) in enumerate(zip(panel['cells'],selected,strict=True)):
            site=source['id'][1];volume=v.TAU*weights[site];exact+=volume
            count=cell['record_count'];need(type(count) is int and count>=0,'integer local record count')
            delta=4096*rate*volume-count
            need(delta.sign()>=0 and (delta-1).sign()<0,'exact production discrepancy')
            need(cell=={'field_event':source['id'],'field_hash':source['hash'],'site':site,
                'action_cell_measure_Qphi':volume.encode(),'record_count':count},'cell identity')
            for version in range(1,count+1):
                ident=[strength,rate,index,version]
                expected={'id':ident,'site':site,'read':source['write'],
                          'write':[128+index,version,ident,source['write'][3]],'ledger_parent':chain}
                chain=digest(expected);expected['hash']=chain
                need(offset<len(panel['observer_records']) and panel['observer_records'][offset]==expected,
                     'actual same-site observer read/write/version/value')
                offset+=1
        need(offset==len(panel['observer_records'])==panel['raw_record_count']==panel['extra_local_read_write_operations'],'complete observer bill')
        need(chain==panel['final_observer_hash'],'observer commitment')
        estimate=F(offset,4096*rate);bound=F(len(selected),4096*rate)
        need(panel['calibrated_action_measure']==str(estimate) and panel['one_sided_error_bound']==str(bound)
             and panel['exact_action_measure_Qphi']==exact.encode(),'calibrated measure/error')
        need((exact-estimate).sign()>=0 and (exact-estimate-bound).sign()<0,'finite action-measure enclosure')
        observer_events+=offset
    for first,second in ((0,2),(1,3)):
        a,b=packet['emission_panels'][first],packet['emission_panels'][second]
        need(a['raw_record_count']==b['raw_record_count'] and a['exact_action_measure_Qphi']==b['exact_action_measure_Qphi'],'source amplitude is not geometry')
    need(packet['emission_panels'][1]['raw_record_count']>packet['emission_panels'][0]['raw_record_count'],'sampling-only density control is nonvacuous')
    need(len(packet['conformal_action_controls'])==2,'geometric control census')
    for control,slope in zip(packet['conformal_action_controls'],(F(0),F(1,2)),strict=True):
        need(control['supplied_profile_slope']==str(slope),'declared geometric control')
        u=[1+slope*parse(address['coordinate_Qphi'][0]) for address in addresses]
        changed=[{j:parse(value) for j,value in row} for row in control['operator_rows_Qphi']]
        need(len(changed)==64,'conformal dimensions')
        # Assemble the complete symmetric conformal edge energy independently.
        stiffness=[{} for _ in mass]
        for i,row in enumerate(rows):
            boundary=mass[i]*(sum(row.values(),R())-1)
            stiffness[i][i]=u[i]**2*mass[i]+u[i]*boundary
        for i,row in enumerate(rows):
            for j,a in row.items():
                if i<j:
                    conductance=-(u[i]+u[j])*mass[i]*a/2
                    stiffness[i][i]+=conductance;stiffness[j][j]+=conductance
                    stiffness[i][j]=stiffness[j][i]=-conductance
        expected=[{j:k/(mass[i]*u[i]) for j,k in row.items()} for i,row in enumerate(stiffness)]
        need(changed==expected,'same conformal action stress/kinetic/gradient coefficients')
        recovered=list(map(parse,control['normalized_action_weights_Qphi']))
        certify_measure(v,changed,recovered)
        need(control['supplied_mass_potential_squared']=='1' and control['calibration_site']==21,'declared mass scale and interior calibration')
        calibration=sum(changed[21].values(),R())
        need(control['interior_rowsum_Qphi']==calibration.encode(),'operator-derived calibration')
        for i,address in enumerate(addresses):
            if all(k in (1,2) for k in address['site']):
                need(sum(changed[i].values(),R())==u[i],'all interior row-sum calibrations')
        reconstructed=[(recovered[i]*weights[21])/(recovered[21]*weights[i])*calibration for i in range(64)]
        need(reconstructed==u and control['reconstructed_u_Qphi']==[x.encode() for x in u],'conformal identification under its action model')
        x=[parse(address['coordinate_Qphi'][0]) for address in addresses]
        beta=(reconstructed[48]-reconstructed[0])/(x[48]-x[0]);alpha=reconstructed[0]-beta*x[0]
        need(control['identified_affine_intercept_Qphi']==alpha.encode() and control['identified_affine_slope_Qphi']==beta.encode()
             and all(value==alpha+beta*coordinate for value,coordinate in zip(reconstructed,x)),'identified conformal profile')
        kinetic=[m*z for m,z in zip(weights,u)];volume=[m*z*z for m,z in zip(weights,u)]
        need(control['normalized_kinetic_measure_Qphi']==[(x/sum(kinetic,R())).encode() for x in kinetic],'kinetic density exponent')
        need(control['normalized_metric_four_volume_Qphi']==[(x/sum(volume,R())).encode() for x in volume],'metric volume exponent')
        need(control['scalar_curvature_Qphi']==[(R(F(3,2)*slope*slope)/(z**3)).encode() for z in u],'specified conformal scalar curvature')
        need(control['raw_support_order_unchanged'] is True and all(set(a)==set(b) for a,b in zip(changed,rows)),'unchanged causal support control')
        if slope:need(control['normalized_kinetic_measure_Qphi']!=control['normalized_metric_four_volume_Qphi'],'action density is not metric four-volume')
        else:need(changed==rows,'flat control reproduces original action')
    tensor=packet['tensor_action_control']
    need(tensor['source_layer']==8 and tensor['previous_source_layer']==7 and tensor['potential_endpoint']=='right'
         and tensor['interpretation']==
         'Fixed-configuration discrete-action derivatives including kinetic and potential terms; no physical stress coupling selected.',
         'tensor interpretation')
    sites=[tuple(address['site']) for address in addresses];index={site:i for i,site in enumerate(sites)}
    permutation=[index[(p[1],p[0],p[2])] for p in sites]
    need(tensor['x_y_site_permutation']==permutation and len(tensor['controls'])==2,'tensor family symmetry')
    need(all(mass[i]==mass[permutation[i]] and {permutation[j]:a for j,a in rows[i].items()}==rows[permutation[i]]
             for i in range(64)),'permutation commutes with complete source action')
    original=list(map(parse,parent['traces']['ascending_intervention']['layers'][8]))
    previous=list(map(parse,parent['traces']['ascending_intervention']['layers'][7]))
    points=[R()]+[parse(addresses[16*i]['coordinate_Qphi'][0]) for i in range(4)]+[R(1)]
    def stretched_energy(field,axis,stretch):
        # Reassemble the differential operator on the stretched source grid,
        # instead of consuming the producer's directional edge energies.
        total_energy=R()
        for i,site in enumerate(sites):
            force=field[i]
            for direction,k in enumerate(site):
                scale=stretch if direction==axis else F(1)
                left=(points[k+1]-points[k])*scale;right=(points[k+2]-points[k+1])*scale
                width=(left+right)/2
                for step,gap in ((-1,left),(1,right)):
                    neighbor=list(site);neighbor[direction]+=step;neighbor=tuple(neighbor)
                    value=field[index[neighbor]] if neighbor in index else R()
                    force+=(field[i]-value)/(width*gap)
            total_energy+=mass[i]*stretch*field[i]*force/2
        return total_energy
    for control,label,field,old in zip(tensor['controls'],('x_preparation','y_preparation'),
            (original,[original[p] for p in permutation]),
            (previous,[previous[p] for p in permutation]),strict=True):
        need(control['label']==label and control['field_Qphi']==[x.encode() for x in field] and
             control['previous_field_Qphi']==[x.encode() for x in old],'actual tensor comparison field')
        energy=v.inner(mass,field,v.apply(action,field))/2
        potential=v.inner(mass,field,field)/2
        kinetic=v.inner(mass,[a-b for a,b in zip(field,old)],[a-b for a,b in zip(field,old)])/(2*v.TAU)
        gradient=[];derivative=[];action_derivative=[]
        for axis in range(3):
            enlarged=stretched_energy(field,axis,F(2));contracted=stretched_energy(field,axis,F(1,2))
            gradient.append(F(2,3)*(2*energy-enlarged))
            derivative.append(F(2,3)*(enlarged-contracted))
            action_derivative.append(F(2,3)*((2*kinetic-v.TAU*enlarged)-(kinetic/2-v.TAU*contracted)))
        need(control['directional_gradient_energy_Qphi']==[x.encode() for x in gradient],
             'independent directional gradient energies')
        need(control['mass_potential_energy_Qphi']==potential.encode() and
             control['total_potential_energy_Qphi']==energy.encode() and sum(gradient,R())+potential==energy,
             'complete tensor-action energy')
        need(control['log_axis_stretch_potential_derivative_Qphi']==[x.encode() for x in derivative],
             'independent potential-energy derivatives')
        need(control['discrete_kinetic_action_Qphi']==kinetic.encode() and
             control['discrete_lagrangian_Qphi']==(kinetic-v.TAU*energy).encode() and
             control['log_axis_stretch_action_derivative_Qphi']==[x.encode() for x in action_derivative] and
             control['isotropic_action_dilation_derivative_Qphi']==sum(action_derivative,R()).encode(),
             'independent axis-stretch derivatives')
    first,second=tensor['controls']
    need(first['total_potential_energy_Qphi']==second['total_potential_energy_Qphi'] and
         first['discrete_kinetic_action_Qphi']==second['discrete_kinetic_action_Qphi'] and
         first['isotropic_action_dilation_derivative_Qphi']==second['isotropic_action_dilation_derivative_Qphi'] and
         first['log_axis_stretch_action_derivative_Qphi']!=second['log_axis_stretch_action_derivative_Qphi'],
         'scalar energy and dilation response do not determine tensor response')
    need(packet['scope']=={'source_population_action_read_grammar_and_clock_supplied':True,
        'same_family_local_observers_check_actual_immutable_values':True,
        'observer_operations_charged_but_physical_durations_not_identified':True,
        'aggregate_certificate_is_not_single_observer_remote_record_transport':True,
        'deterministic_production_law_selected_by_axioms':False,
        'global_action_scale_identified_by_reciprocity':False,
        'kinetic_density_identified_as_metric_volume_without_constitutive_inputs':False,
        'conformal_controls_are_explicitly_supplied_not_sourced_geometries':True,
        'curvature_from_unweighted_order_counts':False,
        'physical_stress_coupling_or_Einstein_response_derived':False},'scientific interpretation')
    return {'verified':True,'sites':64,'undirected_reciprocity_checks':edges,
            'field_events_independently_replayed':replayed,'local_observer_records_replayed':observer_events,
            'factorial_density_controls':4,'conformal_action_controls':2,
            'matched_scalar_distinct_tensor_controls':2,
            'relative_action_measure_identified':True,'physical_source_law_selected':False}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.parse_args()
    print(json.dumps(verify(),indent=2,sort_keys=True))
