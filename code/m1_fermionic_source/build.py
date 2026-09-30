"""Emit compact evidence after independent replay, never by blessing hashes alone."""

import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = 'm1_fermionic_source'

from . import circuits, experiments, model, preparation, spatial, verify, walk


def candidate():
    return dict(small_graphs=model.small_graphs(), native_rotations=circuits.gadgets(),
                native_measurements=circuits.measurements(), weighted_measurements=circuits.weighted_measurements(),
                geometry=spatial.graphs(),
                preparation=preparation.evidence(), budgets=spatial.budgets(), walks=walk.evidence(),
                conversions=experiments.conversions(), resolved_reads=experiments.resolved_reads(),
                detector=experiments.detector(), interaction=experiments.interaction(),
                noise=experiments.noise(), clock_witness=spatial.clock_witness())


def build():
    evidence = candidate()
    verify.parent.verify()
    verify.verify_evidence(evidence)
    packet = dict(schema='oph-fermionic-source-v1', sources=verify.pins(),
                  claims=verify.claim_pins(), evidence=evidence)
    (verify.HERE/'receipt.json').write_text(json.dumps(packet, sort_keys=True, separators=(',', ':'),
                                                      allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print('Built independently replayed fermionic-source receipt')


if __name__ == '__main__':
    build()
