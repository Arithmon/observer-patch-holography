"""Build only after full independent replay of the candidate evidence."""

import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = 'm1_operational_clocks'

from . import model, verify
from .check import verify_evidence


def build():
    evidence = model.candidate()
    verify_evidence(evidence)
    packet = dict(schema='oph-m1-operational-clocks-v1', sources=verify.pins(), claims=verify.claim_pins(), evidence=evidence)
    (verify.HERE/'receipt.json').write_text(json.dumps(packet, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n',
                                          encoding='utf-8', newline='\n')
    print('Built independently replayed operational clock receipt')


if __name__ == '__main__':
    build()
