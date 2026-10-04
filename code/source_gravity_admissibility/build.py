"""Recompute the complete compact primary experiment packet."""
import json
from .model import candidate
from .verify import HERE, pins, claim_pin


def build():
    return dict(schema='oph-source-gravity-admissibility-v1', sources=pins(),
                claim=claim_pin(), evidence=candidate())


def render(packet):
    # One line per experiment preserves readable coverage without huge matrix tapes.
    evidence = packet['evidence']
    head = {k: v for k, v in packet.items() if k != 'evidence'}
    text = json.dumps(head, indent=2, sort_keys=True)[:-2]+',\n  "evidence": {\n'
    parts = []
    for key, value in evidence.items():
        if key in ('line', 'cube', 'clocks', 'symbols', 'robustness', 'timelines', 'bands'):
            encoded = '[\n'+',\n'.join('      '+json.dumps(row, sort_keys=True, allow_nan=False) for row in value)+'\n    ]'
        else:
            encoded = json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)
        parts.append('    '+json.dumps(key)+': '+encoded)
    return text+',\n'.join(parts)+'\n  }\n}\n'


if __name__ == '__main__':
    packet = build()
    (HERE/'receipt.json').write_text(render(packet), encoding='ascii', newline='\n')
    print('Built complete admitted source process evidence')
