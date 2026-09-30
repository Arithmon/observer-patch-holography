"""Create a compact receipt only after independent replay."""

import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = 'm1_source_realization'

from . import circuits, model, topology, verify
from .check import verify_evidence


def build():
    evidence = model.candidate() | circuits.candidate() | {'topology': topology.candidate()}
    verify_evidence(evidence)
    packet = dict(schema='oph-source-realization-v1', sources=verify.pins(), claims=verify.claim_pins(), evidence=evidence)
    (verify.HERE/'receipt.json').write_text(json.dumps(packet, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n',
                                          encoding='utf-8', newline='\n')
    print('Built independently replayed source-realization receipt')


if __name__ == '__main__':
    build()
