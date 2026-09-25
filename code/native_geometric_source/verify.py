"""Independent exact bit-mask enumeration; never imports the producer.

Uses a lazy edge transposition transition rather than nearest-integer averaging,
and checks all retained covariance entries, including separated coarse blocks.
The finite numerical chart-volume control is checked separately from exact
rational load moments. None of these controls measures physical curvature.
"""
from __future__ import annotations
import argparse
import ast
from fractions import Fraction as F
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import numpy as np

RER=Path(__file__).resolve().parents[2]
EVIDENCE=RER/'evidence/native_geometric_source_20260925'
VENDOR=RER/'evidence/observer_dynamics_20260925/oph-physics-sim'

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def require(condition,message):
    if not condition: raise ValueError(message)
def arr(value): return np.array([[F(x) for x in row] for row in value],dtype=object)
def eq(actual,expected,label):
    require(arr(actual).shape==expected.shape and np.array_equal(arr(actual),expected),label)
def gram(a,b,den=1):
    return np.array([[F(sum(int(a[i,j])*int(b[i,k]) for i in range(len(a))),den)
                      for k in range(b.shape[1])] for j in range(a.shape[1])],dtype=object)

def source_edges():
    tree=ast.parse((VENDOR/'oph_fpe/dynamics/self_readback_repair_closure.py').read_text())
    faces=None
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='ORIENTED_BASE_FACES' for t in node.targets):
            faces=ast.literal_eval(node.value)
    require(faces is not None,'source incidence missing')
    return sorted({tuple(sorted((face[i],face[(i+1)%3]))) for face in faces for i in range(3)})

def check_finite(row,edges,antipodes):
    r=row['raised']; masks=[s for s in range(1<<12) if s.bit_count()==r]
    require(row['configurations']==len(masks),'configuration census')
    require(row['preparation']=='uniform law on all fixed-occupancy configurations','source preparation')
    require(row['clock']=='one independent uniform seam attempt and fair endpoint coin','clock definition')
    where={s:i for i,s in enumerate(masks)}
    neighbor=[[] for _ in range(12)]
    for a,b in edges: neighbor[a].append(b);neighbor[b].append(a)
    occupation=np.array([[(s>>i)&1 for i in range(12)] for s in masks],dtype=np.int64)
    load=12*occupation-r
    drive=np.array([[sum(int(z[i]-z[j]) for j in neighbor[i]) for i in range(12)] for z in load])
    mismatch=np.array([[sum(int(n[i]!=n[j]) for j in neighbor[i]) for i in range(12)] for n in occupation])
    p=F(r,12); kappa=F(r*(12-r),132)
    require(F(row['kappa'])==kappa,'absolute kappa')
    cz=np.array([[kappa*(F(int(i==j))-F(1,12)) for j in range(12)] for i in range(12)],dtype=object)
    a,b=edges[0];c,d=antipodes[a],antipodes[b]
    require(row['coarse_blocks']==[[a,b],[c,d]] and len({a,b,c,d})==4,'separated coarse blocks')
    coarse=np.zeros((2,12),dtype=object);coarse[0,a]=coarse[0,b]=F(1,2);coarse[1,c]=coarse[1,d]=F(1,2)
    lap=np.array([[len(neighbor[i]) if i==j else -int(j in neighbor[i]) for j in range(12)] for i in range(12)],dtype=np.int64)
    normalized={}
    for name,(value,scale) in {'load':(load,12),'local_drive':(drive,12),'local_mismatch':(mismatch,1)}.items():
        record=row['readouts'][name]; n=len(value)
        mean=np.array([F(sum(int(x) for x in value[:,i]),n*scale) for i in range(12)],dtype=object)
        require([str(x) for x in mean]==record['mean'],'mean '+name)
        bmat=gram(value,load,n*scale*12)/kappa
        # E[z]=0, so no centering term in the density projection.
        eq(record['density_projection_B'],bmat,'projection '+name)
        moments=[];f=value.copy()
        for t in range(7):
            moments.append(gram(value,f,n*scale*scale*60**t)-np.outer(mean,mean))
            if t<6:
                next_f=30*f.copy()
                for a,b in edges:
                    indexes=[where[s ^ ((1<<a)|(1<<b))] if ((s>>a)^(s>>b))&1 else where[s] for s in masks]
                    next_f+=f[indexes]
                f=next_f
        require(len(record['lag_covariance'])==7,'lag census')
        for t in range(7): eq(record['lag_covariance'][t],moments[t],f'lag {name} {t}')
        eq(record['instant_covariance'],moments[0],'instant '+name)
        residual=moments[0]-bmat@cz@bmat.T
        eq(record['residual_covariance'],residual,'residual '+name)
        for stride,key in [(1,'record_average_covariance'),(2,'two_attempts_per_record_covariance')]:
            require(set(record[key])=={'1','2','4'},'window census')
            for count in (1,2,4):
                total=count*moments[0]
                for t in range(1,count): total+=(count-t)*(moments[stride*t]+moments[stride*t].T)
                eq(record[key][str(count)],total/(count*count),'clock/window '+name)
        eq(record['cross_scale_covariance'],moments[0]@coarse.T,'cross scale '+name)
        eq(record['coarse_covariance'],coarse@moments[0]@coarse.T,'separated coarse '+name)
        if name in ('load','local_drive'):
            normalized[name]=moments[0]/np.trace(moments[0])
            eq(row['unit_total_variance_covariances'][name],normalized[name],'normalized '+name)
        if r==6 and name=='local_mismatch':
            require(all(v==0 for v in bmat.flat),'half-filled mismatch must be wholly residual')
            require(np.trace(residual)>0,'nonlinear residual is nonzero')
    require(all(normalized['load'][i,i]==F(1,12)==normalized['local_drive'][i,i] for i in range(12)),'matched diagonal')
    require(normalized['load'][0,antipodes[0]]==F(-1,132),'load antipodal value')
    require(normalized['local_drive'][0,antipodes[0]]==0,'drive antipodal value')
    require(normalized['local_drive'][edges[0][0],edges[0][1]]==F(-1,45),'drive adjacent value')
    require(not np.array_equal(normalized['load'],normalized['local_drive']),'source shape remains unselected')
    for key in ('geometry_q_covariance','geometry_density_projection_B','geometry_residual_covariance'):
        eq(row[key],np.zeros((12,12),dtype=object),'fixed geometry '+key)
    return {'raised':r,'configurations':len(masks),'exact_matrices_checked':3*(1+7+1+1+6+2)+2+3}

def determinant(v):
    a,b,c=v
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))

def check_trace(row):
    n=row['ports'];require(n==240*4**row['level'],'native rung size')
    initial=row['initial_loads'];final=row['final_loads'];require(len(initial)==n==len(final),'trace dimension')
    require(set(initial)<={0,1} and sum(initial)==n//2,'fixed-total preparation')
    if row['preparation']=='clustered_fixed_total': require(initial==[1]*(n//2)+[0]*(n//2),'clustered law')
    else: require(row['preparation']=='uniform_fixed_total' and initial!=[1]*(n//2)+[0]*(n//2),'distinct uniform draw')
    x=list(initial);swaps=waits=0
    require(len(row['seams'])==len(row['coins'])==64,'attempt census')
    for (a,b),coin in zip(row['seams'],row['coins']):
        require(0<=a<n and 0<=b<n and a!=b and coin in (0,1),'valid attempt')
        s=x[a]+x[b];low=s//2;high=s-low
        target=(high,low) if coin else (low,high)
        if target==(x[a],x[b]): waits+=1
        else: swaps+=1
        x[a],x[b]=target
    require(x==final and row['swaps']==swaps and row['waits']==waits,'independent native replay')
    vertices=row['vertices'];faces=row['faces'];volumes=row['reference_cone_volumes']
    require(len(faces)==20*4**row['level']==len(volumes),'geometry cell census')
    for face,volume in zip(faces,volumes):
        require(volume>0 and abs(abs(determinant([vertices[i] for i in face]))/6-volume)<1e-14,'chart reference volume')
    require(row['geometry_max_log_ratio']==0 and row['geometry_is_supplied_reference'] is True,'geometry interpretation')
    require(row['physical_collar_volume_measured'] is False,'unmeasured physical volume')
    require(0<=row['rotation_volume_max_error']<1e-13 and 0<=row['dilation_volume_max_error']<1e-12,'chart controls')
    return {'level':row['level'],'preparation':row['preparation'],'attempts':64,'swaps':swaps}

def verify(path=EVIDENCE/'receipt.json'):
    row=json.loads(Path(path).read_text());spec=json.loads((EVIDENCE/'spec.json').read_text())
    require(row['schema']=='oph.native-geometric-source.receipt.v1','schema')
    require(row['spec_sha256']==sha(EVIDENCE/'spec.json'),'spec hash')
    require(row['producer_sha256']==sha(Path(__file__).with_name('build.py')),'producer hash')
    archive=json.loads((VENDOR.parent/'manifest.json').read_text())['files']
    require(row['source_sha256']==spec['source_sha256'],'source pin census')
    for relative,digest in row['source_sha256'].items():
        require(sha(VENDOR/relative)==digest==archive['oph-physics-sim/'+relative]['sha256'],'native source custody')
    require(spec['carrier_occupancies']==[1,6] and spec['record_attempt_counts']==[1,2,4] and spec['federation_levels']==[0,1],'frozen controls')
    require(row['observational_data_read'] is False,'target independence declaration')
    require(row['nonclaims']==spec['nonclaims'] and len(row['nonclaims'])==5,'nonclaims')
    require(row['decision']=='FIXED_GEOMETRY_SOURCE_EXCLUDED_AND_VOLUME_COUPLING_UNIDENTIFIED','decision boundary')
    require(row['remaining_interface']=='Native local volume/metric update plus independently specified reference slice and observable calibration; an attached load function alone is not that interface.','missing interface')
    edges=source_edges();require(row['seams']==[list(e) for e in edges],'native seam incidence')
    distance=[[0 if i==j else 99 for j in range(12)] for i in range(12)]
    for a,b in edges: distance[a][b]=distance[b][a]=1
    for k in range(12):
        for i in range(12):
            for j in range(12): distance[i][j]=min(distance[i][j],distance[i][k]+distance[k][j])
    anti=[next(j for j in range(12) if distance[i][j]==3) for i in range(12)]
    require(row['antipodes']==anti,'native antipodes')
    require([r['raised'] for r in row['finite_controls']]==[1,6],'two distinct finite laws')
    finite=[check_finite(r,edges,anti) for r in row['finite_controls']]
    require([(r['level'],r['preparation']) for r in row['geometry_traces']]==[(l,p) for l in (0,1) for p in ('uniform_fixed_total','clustered_fixed_total')],'full trace census')
    traces=[check_trace(r) for r in row['geometry_traces']]
    check=subprocess.run([sys.executable,str(Path(__file__).with_name('check_native_geometry.py')),str(Path(path).resolve()),str(EVIDENCE/'spec.json')],capture_output=True,text=True)
    require(check.returncode==0,'native geometry/preparation custody: '+check.stderr)
    return {'status':'PASS','receipt_sha256':sha(path),'verifier_sha256':sha(__file__),
            'native_custody_checker_sha256':sha(Path(__file__).with_name('check_native_geometry.py')),
            'finite_controls':finite,'trace_controls':traces,'scope':'Exact finite load moments and replayed fixed supplied geometry; no geometric source identification.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--receipt',type=Path,default=EVIDENCE/'receipt.json');p.add_argument('--write',action='store_true');args=p.parse_args()
    result=verify(args.receipt)
    if args.write: (EVIDENCE/'verification_receipt.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))
