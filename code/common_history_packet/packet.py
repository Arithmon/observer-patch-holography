"""A finite scalar history whose local readout feedback is used by the next step.

Classical exact software operations and supplied model time; no quantum or
physical-clock interpretation is assigned to the execution.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'code/source_scalar_execution'))
from source_scalar_execution import model, TAU, Q, inner
from scalar_execution_algebra import rational_bounds

OUT = ROOT / 'evidence/common_history_packet_20260925/packet.json'
PARENT = 'code/source_scalar_execution/source_scalar_execution_receipt.json'
PARENT_SHA = '6981fcd4e13fbc33f371224228341a995f34b4733aa70c2ec0b3f6ee5975c7af'
STEPS = 17
EPS = F(1, 10**6)
PREP = F(1, 10**6)
DETECT = F(1, 10**7)
CASES = ('baseline', 'intervention', 'double_input', 'noisy_baseline',
         'noisy_intervention', 'disabled_feedback')


def canonical(x):
    return (json.dumps(x, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


def digest(x):
    return sha256(canonical(x)).hexdigest()


def model_identity(mass, rows, velocity, detector, addresses):
    return {'mass': [x.encode() for x in mass],
            'action': [[[j, x.encode()] for j, x in row] for row in rows],
            'velocity': [x.encode() for x in velocity],
            'detector': [x.encode() for x in detector], 'addresses': addresses,
            'tau': TAU.encode(), 'boundary': 'Dirichlet', 'dimensionless_mass': '1'}


class Recorder:
    def __init__(self, identity):
        self.state, self.writers, self.events = {}, {}, []
        self.chain = identity

    def emit(self, op, args, keys, function):
        if len(keys) != len(set(keys)):
            raise ValueError('duplicate read key')
        reads = {k: {'writer': self.writers[k], 'value': self.state[k].encode()} for k in keys}
        output = function({k: self.state[k] for k in keys})
        event = {'id': len(self.events), 'op': op, 'args': args, 'reads': reads,
                 'writes': {k: Q.of(v).encode() for k, v in output.items()}, 'parent': self.chain}
        event['hash'] = digest(event)
        self.events.append(event)
        self.chain = event['hash']
        for k, v in output.items():
            self.state[k] = Q.of(v)
            self.writers[k] = event['id']


def execute(case, identity, mass, rows, velocity, detector):
    selected = [i for i, g in enumerate(detector) if g.sign()]
    rec = Recorder(digest({'model': identity, 'case': case}))
    amplitude = 0 if 'baseline' in case else (2 if case == 'double_input' else 1)
    noisy = case.startswith('noisy')
    sign = -1 if case == 'noisy_baseline' else 1
    constants = {'tau': TAU, 'zero': Q(), 'clock': Q(), 'gain': Q(amplitude)}
    constants.update({f'm/{i}': x for i, x in enumerate(mass)})
    constants.update({f'g/{i}': x for i, x in enumerate(detector)})
    constants.update({f'v/{i}': x for i, x in enumerate(velocity)})
    constants.update({f'K/{i}/{j}': x for i, row in enumerate(rows) for j, x in row})
    constants.update({f'prep/{i}': Q(sign*PREP if noisy and i == 0 else 0) for i in range(64)})
    constants.update({f'eta/{j}/{i}': Q(sign*EPS*(-1)**(i+j) if noisy else 0)
                      for j in range(1, STEPS+1) for i in selected})
    constants.update({f'nu/{j}': Q(sign*DETECT*(-1)**j if noisy else 0) for j in range(1, STEPS+1)})
    rec.emit('inputs', [], [], lambda _: constants)
    for i in range(64):
        rec.emit('prepare', [i], ['zero', 'tau', 'gain', f'v/{i}', f'prep/{i}'],
                 lambda v, i=i: {f'u/0/{i}': v['zero'],
                                 f'u/1/{i}': -v['tau']*(v['gain']*v[f'v/{i}']+v[f'prep/{i}'])})
    for i in selected:
        rec.emit('pointer_reset', [i], ['zero'], lambda v, i=i: {f'p/{i}': v['zero']})
    for j in range(1, STEPS+1):
        previous, current = j % 2, (j-1) % 2
        for i in range(64):
            keys = [f'u/{previous}/{i}', 'tau']
            keys += [f'u/{current}/{k}' for k, _ in rows[i]]
            keys += [f'K/{i}/{k}' for k, _ in rows[i]]
            rec.emit('evolve', [j, i], keys, lambda v, i=i, previous=previous, current=current: {
                f'u/{previous}/{i}': 2*v[f'u/{current}/{i}']-v[f'u/{previous}/{i}']-
                v['tau']**2*sum((v[f'K/{i}/{k}']*v[f'u/{current}/{k}'] for k, _ in rows[i]), Q())})
        rec.emit('tick', [j], ['clock', 'tau'], lambda v: {'clock': v['clock']+v['tau']})
        rec.emit('accumulator_reset', [j], ['zero'], lambda v: {'acc': v['zero']})
        for i in selected:
            field, pointer = f'u/{previous}/{i}', f'p/{i}'
            rec.emit('average_probe', [j, i], [field, pointer],
                     lambda v, field=field, pointer=pointer: {
                         field: (v[field]+v[pointer])/2, pointer: (v[field]+v[pointer])/2})
            rec.emit('retain', [j, i], [pointer, f'eta/{j}/{i}'], lambda v, j=j, i=i, pointer=pointer: {
                f'r/{j}/{i}': v[pointer]+v[f'eta/{j}/{i}']})
            op = 'feedback' if case != 'disabled_feedback' else 'feedback_disabled'
            rec.emit(op, [j, i], [f'r/{j}/{i}', 'zero'],
                     lambda v, j=j, i=i, field=field, pointer=pointer, op=op: {
                         **({field: 2*v[f'r/{j}/{i}']} if op == 'feedback' else {}), pointer: v['zero']})
            rec.emit('route', [j, i], [f'r/{j}/{i}'],
                     lambda v, j=j, i=i: {f'inbox/{j}/{i}': v[f'r/{j}/{i}']})
            rec.emit('accumulate', [j, i], ['acc', f'inbox/{j}/{i}', f'm/{i}', f'g/{i}'],
                     lambda v, j=j, i=i: {'acc': v['acc']+2*v[f'm/{i}']*v[f'g/{i}']*v[f'inbox/{j}/{i}']})
        rec.emit('public', [j], ['acc', 'clock', f'nu/{j}'], lambda v, j=j: {
            f'D/{j}': v['acc']+v[f'nu/{j}'], f't/{j}': v['clock']})
    return {'case': case, 'model_sha256': identity, 'events': rec.events,
            'final_hash': rec.chain,
            'detector_intervals': {str(j): rational_bounds(rec.state[f'D/{j}']) for j in range(1, STEPS+1)}}


def build():
    parent = (ROOT/PARENT).read_bytes()
    if sha256(parent).hexdigest() != PARENT_SHA:
        raise ValueError('parent custody')
    _, _, mass, rows, velocity, detector, addresses = model()
    identity = model_identity(mass, rows, velocity, detector, addresses)
    return {'schema': 'oph.common-history-feedback.v1',
        'parent': {'path': PARENT, 'sha256': PARENT_SHA}, 'model': identity,
        'steps': STEPS, 'cases': [execute(c, digest(identity), mass, rows, velocity, detector) for c in CASES],
        'contract': {'record_error': str(EPS), 'preparation_velocity_mass_norm_error': str(PREP),
          'final_detector_error': str(DETECT), 'model_time_readback_error': '1/1000000',
          'window_step_endpoints': [15, 17], 'reconstruction': 'causal previous-sample hold',
          'observer_partition': '64 field patches; pointers at detector sites32..63; public collector owns inbox/acc/D; controller owns clock',
          'units': 'dimensionless parent action time c*t/L and parent field normalization',
          'clock': 'one supplied tau per completed action step; auxiliary events counted separately',
          'routing': 'supplied authenticated classical record channels to one detector accumulator',
          'auxiliary_dynamics': 'action state advances only at evolve events; no physical duration assigned',
          'arithmetic': 'exact Q(phi); coefficient and finite localization errors zero on the declared finite model',
          'quantum_claim': False, 'physical_clock_claim': False, 'source_selected_action_claim': False}}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args()
    raw = canonical(build())
    if args.check:
        if OUT.read_bytes() != raw:
            raise SystemExit('packet drift')
    else:
        OUT.write_bytes(raw)
    print(json.dumps({'receipt': str(OUT.relative_to(ROOT)), 'bytes': len(raw)}))
