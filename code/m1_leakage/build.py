"""Produce compact evidence, then independently verify before serialization."""

import json
from . import interfaces, noise, correction, scaling, composition, verify


def candidate():
    return dict(interfaces=interfaces.candidate(), noise=noise.candidate(),
                correction=correction.candidate(), scaling=scaling.candidate(), composition=composition.candidate())


def build():
    evidence = candidate()
    verify.parent.verify()
    verify.verify_evidence(evidence)
    packet = dict(schema='oph-native-leakage-v1', sources=verify.pins(), claims=verify.claims(), evidence=evidence)
    (verify.HERE/'receipt.json').write_text(
        json.dumps(packet, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n',
        encoding='utf-8', newline='\n')
    print('Built independently replayed native-leakage receipt')


if __name__ == '__main__':
    build()
