"""Produce a compact candidate only after independent semantic replay."""

import json
from . import graphs, circuit, decoder, bounds, export, verify


def candidate():
    return dict(graphs=graphs.candidate(), circuit=circuit.candidate(), decoder=decoder.candidate(),
                bounds=bounds.candidate(), export=export.candidate())


def build():
    evidence = candidate()
    verify.parent.verify()
    verify.verify_evidence(evidence)
    packet = dict(schema='oph-fixed-strength-public-records-v1', sources=verify.pins(),
                  claims=verify.claims(), evidence=evidence)
    (verify.HERE/'receipt.json').write_text(
        json.dumps(packet, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n',
        encoding='utf-8', newline='\n')
    print('Built independently replayed fixed-strength public-record receipt')


if __name__ == '__main__':
    build()
