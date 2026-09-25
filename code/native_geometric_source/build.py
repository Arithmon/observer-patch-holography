"""Small native fixed-geometry and source-identification controls.

All finite moments are exact rational numbers. The production trace supplies
only an implementation witness; the frozen-geometry statement is universal
under the explicitly unchanged geometry, as proved in Lean and the paper.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np

RER = Path(__file__).resolve().parents[2]
EVIDENCE = RER / 'evidence/native_geometric_source_20260925'
ARCHIVE = RER / 'evidence/observer_dynamics_20260925'
VENDOR = ARCHIVE / 'oph-physics-sim'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def encode(value):
    if isinstance(value, F): return str(value)
    if isinstance(value, np.ndarray): return encode(value.tolist())
    if isinstance(value, (tuple, list)): return [encode(v) for v in value]
    if isinstance(value, dict): return {k: encode(v) for k,v in value.items()}
    if isinstance(value, np.integer): return int(value)
    return value

def matrix(a, denominator=1):
    return np.array([[F(int(v), denominator) for v in row] for row in a], dtype=object)

def centered_cov(x, y):
    n = len(x)
    return matrix(x.T @ y, n) - np.outer([F(int(v),n) for v in x.sum(0)],
                                       [F(int(v),n) for v in y.sum(0)])

def native_modules(spec):
    for path, digest in spec['source_sha256'].items():
        if sha(VENDOR/path) != digest: raise ValueError('native source digest: '+path)
    sys.path.insert(0, str(VENDOR))
    from oph_exact import carrier, federation
    from oph_fpe.core.icosahedral import build_geodesic_icosahedral_tower
    return carrier, federation, build_geodesic_icosahedral_tower

def finite_control(r, carrier, windows):
    p = 12
    states = np.array([[int(i in s) for i in range(p)]
                       for s in itertools.combinations(range(p), r)], dtype=np.int64)
    where = {tuple(s): i for i,s in enumerate(states)}
    edges = list(carrier.seams())
    lap = carrier.laplacian()
    # Integer numerators: z=x/12. Keeping numerators integral avoids rounding.
    x = 12*states-r
    drive = x@lap
    mismatch = np.zeros_like(states)
    for a,b in edges:
        sq = (states[:,a]-states[:,b])**2
        mismatch[:,a] += sq; mismatch[:,b] += sq
    values = {'load': (x,12), 'local_drive':(drive,12), 'local_mismatch':(mismatch,1)}
    successors = []
    for state in states:
        row=[]
        for a,b in edges:
            for coin in (False,True):
                target=state.copy()
                target[a],target[b]=carrier.integer_nearest_agreement(int(state[a]),int(state[b]),ceiling_to_first=coin)
                row.append(where[tuple(target)])
        successors.append(row)
    successors = np.asarray(successors)
    kappa=F(r*(p-r),p*(p-1))
    cz = centered_cov(x,x)/144
    # Two disjoint two-port blocks. The second is antipodal to the first.
    a,b = edges[0]; c,d = carrier.antipode()[a],carrier.antipode()[b]
    coarse=np.zeros((2,p),dtype=object)
    coarse[0,a]=coarse[0,b]=F(1,2); coarse[1,c]=coarse[1,d]=F(1,2)
    out={}
    for name,(f,den) in values.items():
        bmat=centered_cov(f,x)/(den*12*kappa)
        lag=[]; propagated=f.copy()
        for t in range(2*(max(windows)-1)+1):
            mu=[F(int(v),len(f)*den) for v in f.sum(0)]
            cmat=matrix(f.T@propagated,len(f)*den*den*60**t)-np.outer(mu,mu)
            lag.append(cmat)
            if t<2*(max(windows)-1): propagated=propagated[successors].sum(axis=1)
        covs={}
        for count in windows:
            total=count*lag[0]
            for t in range(1,count): total=total+(count-t)*(lag[t]+lag[t].T)
            covs[str(count)]=total/(count*count)
        altered_clock={}
        for count in windows:
            total=count*lag[0]
            for t in range(1,count): total=total+(count-t)*(lag[2*t]+lag[2*t].T)
            altered_clock[str(count)]=total/(count*count)
        c0=lag[0]; residual=c0-bmat@cz@bmat.T
        out[name]={'mean':[F(int(v),len(f)*den) for v in f.sum(0)],
                   'density_projection_B':bmat,'residual_covariance':residual,
                   'instant_covariance':c0,'lag_covariance':lag,
                   'record_average_covariance':covs,
                   'two_attempts_per_record_covariance':altered_clock,
                   'cross_scale_covariance':c0@coarse.T,
                   'coarse_covariance':coarse@c0@coarse.T}
    normalized={name:out[name]['instant_covariance']/np.trace(out[name]['instant_covariance'])
                for name in ('load','local_drive')}
    return encode({'raised':r,'configurations':len(states),'kappa':kappa,
                   'preparation':'uniform law on all fixed-occupancy configurations',
                   'clock':'one independent uniform seam attempt and fair endpoint coin',
                   'coarse_blocks':[[a,b],[c,d]], 'readouts':out,
                   'unit_total_variance_covariances':normalized,
                   'geometry_q_covariance':np.zeros((12,12),dtype=np.int64),
                   'geometry_density_projection_B':np.zeros((12,12),dtype=np.int64),
                   'geometry_residual_covariance':np.zeros((12,12),dtype=np.int64)})

def geometry_traces(spec, carrier, federation, tower_builder):
    rows=[]
    tower=tower_builder(max(spec['federation_levels']))
    for level in spec['federation_levels']:
        mesh=tower.levels[level]; fed=federation.build_federation(level,'port_pair')
        vertices=mesh.vertices.copy(); faces=mesh.faces.copy()
        volumes=np.abs(np.linalg.det(vertices[faces]))/6
        if np.any(volumes<=0): raise AssertionError('positive reference cell volumes')
        for index,prep in enumerate(spec['preparations']):
            rng=np.random.default_rng(spec['seed']+level*10+index)
            loads=np.zeros(fed.ports,dtype=np.int64)
            if prep=='uniform_fixed_total': loads[rng.choice(fed.ports,fed.ports//2,replace=False)]=1
            else: loads[:fed.ports//2]=1
            initial=loads.copy()
            sequence=rng.integers(0,fed.seams,size=spec['trace_attempts'],dtype=np.int64)
            coins=rng.integers(0,2,size=spec['trace_attempts'],dtype=np.int64)
            tally=federation.IntegerTally()
            federation.integer_moves_python(loads,fed.graph,sequence,coins,
                                             int(initial@initial),fed.ports//2,tally)
            # Pure rotations and common coordinate dilation are chart controls,
            # not changes to the repair generator or physical identifications.
            rotate=np.array([[0,-1,0],[1,0,0],[0,0,1]],dtype=float)
            rotated=np.abs(np.linalg.det((vertices@rotate.T)[faces]))/6
            dilated=np.abs(np.linalg.det((2*vertices)[faces]))/6
            rows.append({'level':level,'preparation':prep,'ports':fed.ports,
                         'initial_loads':initial.tolist(),'final_loads':loads.tolist(),
                         'seams':[[int(fed.seam_a[s]),int(fed.seam_b[s])] for s in sequence],
                         'coins':coins.tolist(),'waits':tally.waits,'swaps':tally.swaps,
                         'geometry_hash':mesh.geometry_hash,
                         'vertices':vertices.tolist(),'faces':faces.tolist(),
                         'reference_cone_volumes':volumes.tolist(),
                         'geometry_max_log_ratio':0.0,
                         'rotation_volume_max_error':float(np.max(np.abs(rotated-volumes))),
                         'dilation_volume_max_error':float(np.max(np.abs(dilated-8*volumes))),
                         'geometry_is_supplied_reference':True,
                         'physical_collar_volume_measured':False})
    return rows

def build():
    spec=json.loads((EVIDENCE/'spec.json').read_text())
    carrier,federation,tower=native_modules(spec)
    result={'schema':'oph.native-geometric-source.receipt.v1',
            'spec_sha256':sha(EVIDENCE/'spec.json'),
            'producer_sha256':sha(__file__), 'observational_data_read':False,
            'source_sha256':spec['source_sha256'],
            'seams':[list(e) for e in carrier.seams()],
            'antipodes':list(carrier.antipode()),
            'finite_controls':[finite_control(r,carrier,spec['record_attempt_counts']) for r in spec['carrier_occupancies']],
            'geometry_traces':geometry_traces(spec,carrier,federation,tower),
            'decision':'FIXED_GEOMETRY_SOURCE_EXCLUDED_AND_VOLUME_COUPLING_UNIDENTIFIED',
            'remaining_interface':'Native local volume/metric update plus independently specified reference slice and observable calibration; an attached load function alone is not that interface.',
            'nonclaims':spec['nonclaims']}
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    output=json.dumps(build(),indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write: (EVIDENCE/'receipt.json').write_text(output)
    else: print(output,end='')
