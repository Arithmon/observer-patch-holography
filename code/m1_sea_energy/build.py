"""Write compact evidence only after independent replay, never by a hash refresh."""

import json
from . import spectrum, preparation, clock, bounds, verify


def candidate():
    return dict(spectrum=spectrum.candidate(), preparation=preparation.candidate(),
                clock=clock.candidate(), bounds=bounds.candidate())


def build():
    evidence = candidate()
    verify.verify_evidence(evidence)
    verify.protected.verify()
    verify.clocks.verify()
    packet = dict(schema='oph-source-sea-energy-v1', sources=verify.pins(), claims=verify.claims(), evidence=evidence)
    (verify.HERE/'receipt.json').write_text(json.dumps(packet, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print('Built independently replayed filled-source receipt')


if __name__ == '__main__':
    build()
