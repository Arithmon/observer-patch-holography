"""Generate a compact candidate; the verifier independently reconstructs it."""

from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from . import model
from .check import digest
from .verify import pins, parent_claims


def graph(config,family,topology):
    masks,cuts,reads=model.graph_cuts(config,family,topology)
    q=config[0]
    vectors=model.stencil(config,family)
    rows={}
    for name,mask in masks.items():
        rows[name]=dict(volume_sites=int(mask.sum()),cut_pairs=cuts[name],
                        ordered_crossing_reads=2*cuts[name],
                        pi_times_normalized_entropy=str(F(4*cuts[name],q**4)),
                        mask_sha256=hashlib.sha256(mask.astype('uint8').tobytes()).hexdigest())
    changed=int((masks['localized']!=masks['connected']).sum())
    return dict(config=list(config),family=family,topology=topology,
                stencil_sha256=digest(vectors),degree=len(vectors),ordered_reads=reads,masks=rows,
                bridge=dict(changed_sites=changed,cut_change=cuts['connected']-cuts['localized'],
                            degree_bound=len(vectors)*changed,nearest_neighbor_connected=True))


def build():
    return dict(schema='oph-m1-interfaces-v1',sources=pins(),parent_claims=parent_claims(),
                scales={str(t):model.scale(t) for t in range(1,13)},
                moments={kind:model.moment_case(kind) for kind in ('raw','critical','repaired')},
                graphs={f'{c[0]}:{f}:{b}':graph(c,f,b) for c in model.CONFIGS
                        for f in ('T','U') for b in ('periodic','clipped')})


if __name__=='__main__':
    packet=build()
    Path(__file__).with_name('receipt.json').write_text(
        json.dumps(packet,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('Built 12 complete graph comparisons, 96 cut cases, 12 scale levels and 3 moment regimes')
