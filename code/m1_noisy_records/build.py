"""Write a compact receipt only after independent replay."""

import json
from .circuits import candidate
from . import verify


def build():
    evidence = candidate()
    verify.parent.verify()
    verify.verify_evidence(evidence)
    packet = dict(schema='oph-noisy-source-records-v1', sources=verify.pins(),
                  claims=verify.claims(), evidence=evidence)
    (verify.HERE/'receipt.json').write_text(
        json.dumps(packet, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n',
        encoding='utf-8', newline='\n')
    print('Built independently replayed noisy-record receipt')


if __name__ == '__main__':
    build()
