"""Check finite trace inputs against the pinned native geometry in a fresh process."""
import json
from pathlib import Path
import sys
import numpy as np
RER=Path(__file__).resolve().parents[2]
VENDOR=RER/'evidence/observer_dynamics_20260925/oph-physics-sim'
sys.path.insert(0,str(VENDOR))
from oph_exact import federation
from oph_fpe.core.icosahedral import build_geodesic_icosahedral_tower

def require(condition, message):
    if not condition: raise ValueError(message)

def check(receipt,spec):
    tower=build_geodesic_icosahedral_tower(1)
    for row in receipt['geometry_traces']:
        level=row['level'];mesh=tower.levels[level];fed=federation.build_federation(level,'port_pair')
        require(row['geometry_hash']==mesh.geometry_hash,'native geometry hash')
        require(row['vertices']==mesh.vertices.tolist() and row['faces']==mesh.faces.tolist(),'native geometry arrays')
        index=spec['preparations'].index(row['preparation'])
        rng=np.random.default_rng(spec['seed']+10*level+index)
        initial=np.zeros(fed.ports,dtype=np.int64)
        if index==0: initial[rng.choice(fed.ports,fed.ports//2,replace=False)]=1
        else: initial[:fed.ports//2]=1
        sequence=rng.integers(0,fed.seams,size=spec['trace_attempts'],dtype=np.int64)
        coins=rng.integers(0,2,size=spec['trace_attempts'],dtype=np.int64)
        require(row['initial_loads']==initial.tolist(),'preparation seed')
        require(row['coins']==coins.tolist(),'native coin sequence')
        require(row['seams']==[[int(fed.seam_a[s]),int(fed.seam_b[s])] for s in sequence],'native routed seams')
    print('PASS: four finite traces use pinned native geometry, preparations and attempts')
if __name__=='__main__':
    check(json.loads(Path(sys.argv[1]).read_text()),json.loads(Path(sys.argv[2]).read_text()))
