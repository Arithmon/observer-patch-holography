"""Replay a complete candidate before replacing its public receipt."""
import json
from . import model, verify


def build():
    packet = model.build()
    verify.check.verify(packet)
    return dict(schema='oph-dark-source-lensing-v1', sources=verify.pins(),
                claim=verify.claim_digest(), evidence=packet)


if __name__ == '__main__':
    (verify.HERE/'receipt.json').write_text(json.dumps(build(), indent=2, sort_keys=True)+'\n',
                                         encoding='ascii', newline='\n')
    print('Built independently replayed dark-source dynamics/lensing certificate')
