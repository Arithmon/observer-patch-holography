"""Capture the declared production W12 support, without running its campaigns.

This optional fixture producer needs the sibling oph-physics-sim checkout.
The scientific verifier consumes the captured finite combinatorics only.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
SIM = HERE.parents[2] / 'oph-physics-sim'
sys.path.insert(0, str(SIM))
from oph_exact.support_wiring import build_wiring, cell_geometry
from oph_exact.carrier import seams


def capture():
    wiring = build_wiring(3)
    geometry = cell_geometry(3)
    parents = ['oph_exact/support_wiring.py', 'oph_exact/carrier.py',
               'oph_exact/federation.py', 'oph_fpe/core/icosahedral.py',
               'oph_fpe/core/screen_ports.py']
    fixture = {
        'schema': 'oph.source_routing.w12_support.v1', 'level': 3,
        'carriers': wiring.carriers,
        'faces': geometry.faces.tolist(),
        'intra_carrier_seams': [list(pair) for pair in seams()],
        'glued_pairs': wiring.pairs.tolist(),
        'source': {'repository': 'oph-physics-sim',
                   'base_checkout_commit': subprocess.check_output(['git', '-C', str(SIM), 'rev-parse', 'HEAD'], text=True).strip(),
                   'source_identity_scope': 'The base checkout commit is contextual; file hashes identify the captured local sources, which need not be committed at that revision.',
                   'files': {p: hashlib.sha256((SIM / p).read_bytes()).hexdigest() for p in parents}},
        'scope': 'Declared W12 geometric port assignment captured from the production builder; no wiring selection is claimed.',
    }
    data = (json.dumps(fixture, sort_keys=True, separators=(',', ':')) + '\n').encode()
    (HERE / 'support_w12_l3.json').write_bytes(data)
    print(hashlib.sha256(data).hexdigest())


if __name__ == '__main__':
    capture()
