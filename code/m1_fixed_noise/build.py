"""Create compact evidence only after independent physics replay."""

import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = 'm1_fixed_noise'

from . import algebra, noise, recovery, scaling, verify


def candidate():
    return dict(decoder=algebra.decoder(), interfaces=algebra.interfaces(), recovery=recovery.audit(),
                memory=algebra.memory(), continuous_noise=noise.examples(),
                accounting=algebra.accounting(), scaling=scaling.evidence())


def build():
    evidence = candidate()
    verify.parent.verify()
    verify.verify_evidence(evidence)
    packet = dict(schema='oph-fixed-rate-source-v1', sources=verify.pins(), claims=verify.claims(), evidence=evidence)
    (verify.HERE/'receipt.json').write_text(json.dumps(packet, sort_keys=True, separators=(',', ':'),
                                                      allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print('Built independently replayed fixed-rate source receipt')


if __name__ == '__main__':
    build()
