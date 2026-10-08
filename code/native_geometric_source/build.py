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
from math import comb
from numbers import Integral
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

def reject_masked(value):
    """Inspect original containers before NumPy can discard an inner mask."""
    if np.ma.isMaskedArray(value):
        raise ValueError('integer matrix cannot contain masked data')
    entries = value.flat if isinstance(value, np.ndarray) else value
    if isinstance(value, (np.ndarray, list, tuple)):
        for entry in entries:
            reject_masked(entry)


def integers(value):
    """Preserve original integer scalars before any product or reduction."""
    reject_masked(value)
    a = np.asarray(value, dtype=object)
    if a.ndim != 2:
        raise ValueError('integer matrix must have two dimensions')
    if any(isinstance(v, (bool, np.bool_)) or not isinstance(v, Integral)
           for v in a.flat):
        raise ValueError('integer matrix requires exact integer entries')
    return np.array([int(v) for v in a.flat], dtype=object).reshape(a.shape)


def matrix(a, denominator=1):
    if (isinstance(denominator, (bool, np.bool_)) or
            not isinstance(denominator, Integral) or denominator <= 0):
        raise ValueError('positive integer denominator required')
    a = integers(a)
    return np.array([[F(v, int(denominator)) for v in row] for row in a], dtype=object)

def centered_cov(x, y):
    x, y = integers(x), integers(y)
    n = len(x)
    if n == 0 or len(y) != n:
        raise ValueError('covariance needs matching nonempty populations')
    return matrix(x.T @ y, n) - np.outer([F(int(v),n) for v in x.sum(0)],
                                       [F(int(v),n) for v in y.sum(0)])

def native_modules(spec):
    for path, digest in spec['source_sha256'].items():
        if sha(VENDOR/path) != digest: raise ValueError('native source digest: '+path)
    sys.path.insert(0, str(VENDOR))
    from oph_exact import carrier, federation
    from oph_fpe.core.icosahedral import build_geodesic_icosahedral_tower
    return carrier, federation, build_geodesic_icosahedral_tower

def moment_coordinates(r, carrier):
    """Exact degree-two closure of the original binary repair operation.

    These 79 monomials span the three readouts. They need not be independent
    on a fixed-occupancy slice; the (possibly singular) moment matrix is
    only multiplied, never inverted. Products use set UNION since n_i^2=n_i.
    """
    p = 12
    edges = list(carrier.seams())
    # Bind the simplification to the supplied primitive, not just its name:
    # for binary endpoints a fair coin is half identity, half transposition.
    for a, b in itertools.product((0, 1), repeat=2):
        targets = [tuple(carrier.integer_nearest_agreement(a, b, ceiling_to_first=coin))
                   for coin in (False, True)]
        if sorted(targets) != sorted(((a, b), (b, a))):
            raise ValueError('binary repair is not a fair lazy transposition')
    if (len(edges) != 30 or len(set(edges)) != 30 or
            any(not (0 <= a < b < p) for a, b in edges)):
        raise ValueError('expected thirty oriented native seams')
    lap = np.zeros((p, p), dtype=object)
    for a, b in edges:
        lap[a, a] += 1; lap[b, b] += 1
        lap[a, b] -= 1; lap[b, a] -= 1
    if not np.array_equal(integers(carrier.laplacian()), lap):
        raise ValueError('native Laplacian and seam operations disagree')
    features = [()] + [(i,) for i in range(p)] + list(itertools.combinations(range(p), 2))
    index = {s: i for i, s in enumerate(features)}
    # Counts, not rounded probabilities: exactly this many r-subsets contain S.
    counts = [comb(p-k, r-k) if k <= r else 0 for k in range(5)]
    moments = np.array([[counts[len(set(s) | set(t))] for t in features]
                        for s in features], dtype=object)
    x = np.zeros((len(features), p), dtype=object)
    x[0] = -r
    for i in range(p):
        x[index[(i,)], i] = 12
    drive = x @ lap
    mismatch = np.zeros_like(x)
    for a,b in edges:
        for port in (a, b):
            mismatch[index[(a,)], port] += 1
            mismatch[index[(b,)], port] += 1
            mismatch[index[(a, b)], port] -= 2
    values = {'load': (x,12), 'local_drive':(drive,12), 'local_mismatch':(mismatch,1)}
    permutations = np.array([
        [index[tuple(sorted(b if i == a else a if i == b else i for i in s))]
         for s in features] for a, b in edges])
    return edges, moments, values, permutations


def finite_control(r, carrier, windows):
    if type(r) is not int or not 1 <= r <= 11:
        raise ValueError('exact occupancy integer in 1..11 required')
    if (type(windows) not in (list, tuple) or not windows or
            any(type(count) is not int or count <= 0 for count in windows) or
            len(set(windows)) != len(windows)):
        raise ValueError('nonempty unique positive integer record windows required')
    p = 12
    edges, moments, values, permutations = moment_coordinates(r, carrier)
    configurations = comb(p, r)
    x = values['load'][0]
    kappa=F(r*(p-r),p*(p-1))
    cz = matrix(x.T @ moments @ x, configurations*144)
    # Two disjoint two-port blocks. The second is antipodal to the first.
    a,b = edges[0]; c,d = carrier.antipode()[a],carrier.antipode()[b]
    coarse=np.zeros((2,p),dtype=object)
    coarse[0,a]=coarse[0,b]=F(1,2); coarse[1,c]=coarse[1,d]=F(1,2)
    out={}
    for name,(f,den) in values.items():
        # E[z]=0 exactly. Integer counts are retained through the Gram product.
        weighted = f.T @ moments
        bmat=matrix(weighted @ x, configurations*den*12)/kappa
        mu=[F(v,configurations*den) for v in moments[0] @ f]
        lag=[]; propagated=f.copy()
        for t in range(2*(max(windows)-1)+1):
            cmat=matrix(weighted@propagated,configurations*den*den*60**t)-np.outer(mu,mu)
            lag.append(cmat)
            if t<2*(max(windows)-1):
                propagated=30*propagated+propagated[permutations].sum(axis=0)
        averages={}
        for stride,key in ((1,'record_average_covariance'),(2,'two_attempts_per_record_covariance')):
            averages[key]={}
            for count in windows:
                total=count*lag[0]
                for t in range(1,count):
                    total=total+(count-t)*(lag[stride*t]+lag[stride*t].T)
                averages[key][str(count)]=total/(count*count)
        c0=lag[0]; residual=c0-bmat@cz@bmat.T
        out[name]={'mean':mu,
                   'density_projection_B':bmat,'residual_covariance':residual,
                   'instant_covariance':c0,'lag_covariance':lag,
                   **averages,
                   'cross_scale_covariance':c0@coarse.T,
                   'coarse_covariance':coarse@c0@coarse.T}
    normalized={name:out[name]['instant_covariance']/np.trace(out[name]['instant_covariance'])
                for name in ('load','local_drive')}
    return encode({'raised':r,'configurations':configurations,'kappa':kappa,
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
    if args.write: (EVIDENCE/'receipt.json').write_bytes(output.encode('utf-8'))
    else: print(output,end='')
