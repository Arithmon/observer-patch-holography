"""Independent replay and composed error certificate, with no producer import.

The inherited verifier assembles undirected edge energies in Q(sqrt(5)).
Here a separate operation interpreter checks actual reads, restoration and
the next action step. Generic error bounds cover adversarial bounded errors,
not just the deterministic signed fixture errors.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'code/source_scalar_execution'))
import verify_source_scalar_execution as parent

R, parse, need, equal = parent.R, parent.parse, parent.require, parent.equal
TAU = parent.TAU
PATH = ROOT/'evidence/common_history_packet_20260925/packet.json'
EVIDENCE = ROOT/'evidence/common_history_packet_20260925'
PIN_PATHS = (
    'code/common_history_packet/packet.py', 'code/common_history_packet/verify.py',
    'code/common_history_packet/test_common_history_packet.py', 'code/common_history_packet/README.md',
    'Lean/Screen/CommonHistoryFeedback.lean',
    'paper/tex_fragments/COMMON_HISTORY_FEEDBACK_PACKET.tex',
    'evidence/common_history_packet_20260925/packet.json',
    'evidence/common_history_packet_20260925/verification.json',
)
PARENT = 'code/source_scalar_execution/source_scalar_execution_receipt.json'
PARENT_SHA = '6981fcd4e13fbc33f371224228341a995f34b4733aa70c2ec0b3f6ee5975c7af'
CASES = ['baseline', 'intervention', 'double_input', 'noisy_baseline', 'noisy_intervention', 'disabled_feedback']
STEPS = 17
EPS, PREP, DETECT, CLOCK = F(1, 10**6), F(1, 10**6), F(1, 10**7), F(1, 10**6)
CONTRACT = {'record_error': str(EPS), 'preparation_velocity_mass_norm_error': str(PREP),
    'final_detector_error': str(DETECT), 'model_time_readback_error': str(CLOCK),
    'window_step_endpoints': [15, 17], 'reconstruction': 'causal previous-sample hold',
    'observer_partition': '64 field patches; pointers at detector sites32..63; public collector owns inbox/acc/D; controller owns clock',
    'units': 'dimensionless parent action time c*t/L and parent field normalization',
    'clock': 'one supplied tau per completed action step; auxiliary events counted separately',
    'routing': 'supplied authenticated classical record channels to one detector accumulator',
    'auxiliary_dynamics': 'action state advances only at evolve events; no physical duration assigned',
    'arithmetic': 'exact Q(phi); coefficient and finite localization errors zero on the declared finite model',
    'quantum_claim': False, 'physical_clock_claim': False, 'source_selected_action_claim': False}


def digest(x):
    return sha256(parent.canonical(x)).hexdigest()


def manifest(root=ROOT):
    return {'schema': 'oph.common-history-feedback.custody.v1',
        'files': {p: {'sha256': sha256((root/p).read_bytes()).hexdigest(), 'bytes': (root/p).stat().st_size}
                  for p in PIN_PATHS}}


def check_manifest(root=ROOT):
    folder = root/'evidence/common_history_packet_20260925'
    need({p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file()} == {'packet.json', 'verification.json', 'manifest.json'},
         'closed evidence file inventory')
    equal(json.loads((folder/'manifest.json').read_text()), manifest(root), 'source/evidence byte custody')


def load(path=PATH):
    raw = Path(path).read_bytes()
    need(len(raw) < 80_000_000, 'bounded packet bytes')
    return json.loads(raw, object_pairs_hook=parent.strict_pairs,
                      parse_float=parent.forbidden_number, parse_constant=parent.forbidden_number)


def enclosure(x):
    scale = 10**12
    lo, hi = -scale, scale
    while (x-R(F(lo, scale))).sign() < 0:
        lo *= 2
    while (x-R(F(hi, scale))).sign() > 0:
        hi *= 2
    while hi-lo > 1:
        middle = (hi+lo)//2
        if (x-R(F(middle, scale))).sign() < 0:
            hi = middle
        else:
            lo = middle
    return [str(F(lo, scale)), str(F(hi, scale))]


def context(full_parent=True):
    raw = (ROOT/PARENT).read_bytes()
    need(sha256(raw).hexdigest() == PARENT_SHA, 'immutable parent bytes')
    packet = parent.load(ROOT/PARENT)
    report = parent.verify(packet) if full_parent else None
    m, a, v, g, addresses = parent.model()
    identity = {'mass': [x.encode() for x in m],
        'action': [[[j, x.encode()] for j, x in row] for row in a],
        'velocity': [x.encode() for x in v], 'detector': [x.encode() for x in g],
        'addresses': addresses, 'tau': TAU.encode(), 'boundary': 'Dirichlet', 'dimensionless_mass': '1'}
    return m, a, v, g, identity, packet, report


def expected_program(detector):
    yield 'inputs', []
    for i in range(64):
        yield 'prepare', [i]
    for i in detector:
        yield 'pointer_reset', [i]
    for j in range(1, STEPS+1):
        for i in range(64):
            yield 'evolve', [j, i]
        yield 'tick', [j]
        yield 'accumulator_reset', [j]
        for i in detector:
            for op in ('average_probe', 'retain', 'feedback', 'route', 'accumulate'):
                yield op, [j, i]
        yield 'public', [j]


def replay(trace, ctx):
    m, a, v0, g, identity, base, _ = ctx
    case = trace['case']
    need(case in CASES, 'known case')
    need(set(trace) == {'case', 'model_sha256', 'events', 'final_hash', 'detector_intervals'}, 'trace schema')
    equal(trace['model_sha256'], digest(identity), 'history action/state/readout identity')
    selected = [i for i, x in enumerate(g) if x.sign() != 0]
    program = list(expected_program(selected))
    need(len(trace['events']) == len(program), 'complete auxiliary/event inventory')
    noisy = case in ('noisy_baseline', 'noisy_intervention')
    amplitude = 0 if case in ('baseline', 'noisy_baseline') else (2 if case == 'double_input' else 1)
    direction = -1 if case == 'noisy_baseline' else 1
    values, writers, outputs, snapshots = {}, {}, {}, {}
    chain = digest({'model': digest(identity), 'case': case})
    counts, read_count, write_count, max_bits = Counter(), 0, 0, 0
    restored_reads, peak_slots, mismatches = 0, 0, 0
    for event, (op, args) in zip(trace['events'], program, strict=True):
        if op == 'feedback' and case == 'disabled_feedback':
            op = 'feedback_disabled'
        event_id = sum(counts.values())
        need(set(event) == {'id', 'op', 'args', 'reads', 'writes', 'parent', 'hash'}, 'event schema')
        equal([event['id'], event['op'], event['args']], [event_id, op, args], 'event program order')
        equal(event['parent'], chain, 'audit chain predecessor')
        material = {k: x for k, x in event.items() if k != 'hash'}
        equal(event['hash'], digest(material), 'event commitment')
        chain = event['hash']
        used = set()
        def read(k):
            nonlocal restored_reads
            need(k in values and k in event['reads'], 'missing consumed record '+k)
            equal(event['reads'][k], {'writer': writers[k], 'value': values[k].encode()}, 'consumed writer/value '+k)
            used.add(k)
            if op == 'evolve' and trace['events'][writers[k]]['op'] == 'feedback':
                restored_reads += 1
            return values[k]
        out = {}
        if op == 'inputs':
            out = {'tau': TAU, 'zero': R(), 'clock': R(), 'gain': R(amplitude)}
            out.update({f'm/{i}': x for i, x in enumerate(m)})
            out.update({f'g/{i}': x for i, x in enumerate(g)})
            out.update({f'v/{i}': x for i, x in enumerate(v0)})
            out.update({f'K/{i}/{k}': x for i, row in enumerate(a) for k, x in row})
            out.update({f'prep/{i}': R(direction*PREP if noisy and i == 0 else 0) for i in range(64)})
            out.update({f'eta/{j}/{i}': R(direction*EPS*(-1)**(i+j) if noisy else 0)
                        for j in range(1, STEPS+1) for i in selected})
            out.update({f'nu/{j}': R(direction*DETECT*(-1)**j if noisy else 0) for j in range(1, STEPS+1)})
        elif op == 'prepare':
            i, = args
            out = {f'u/0/{i}': read('zero'), f'u/1/{i}': -read('tau')*(read('gain')*read(f'v/{i}')+read(f'prep/{i}'))}
        elif op == 'pointer_reset':
            out = {f'p/{args[0]}': read('zero')}
        elif op == 'evolve':
            j, i = args
            old, current = f'u/{j%2}/{i}', (j-1)%2
            # First-order kick in the difference variable, independently of
            # the producer's direct second-order assignment.
            q = read(f'u/{current}/{i}')
            difference = q-read(old)
            force = sum((read(f'K/{i}/{k}')*read(f'u/{current}/{k}') for k, _ in a[i]), R())
            out = {old: q+difference-read('tau')**2*force}
        elif op == 'tick':
            out = {'clock': read('clock')+read('tau')}
        elif op == 'accumulator_reset':
            out = {'acc': read('zero')}
        elif op == 'average_probe':
            j, i = args
            field, pointer = f'u/{j%2}/{i}', f'p/{i}'
            need(read(pointer) == R(), 'reset pointer before probe')
            average = (read(field)+read(pointer))/2
            out = {field: average, pointer: average}
        elif op == 'retain':
            j, i = args
            out = {f'r/{j}/{i}': read(f'p/{i}')+read(f'eta/{j}/{i}')}
        elif op in ('feedback', 'feedback_disabled'):
            j, i = args
            reconstructed, zero = 2*read(f'r/{j}/{i}'), read('zero')
            out = {f'p/{i}': zero}
            if op == 'feedback':
                out[f'u/{j%2}/{i}'] = reconstructed
        elif op == 'route':
            j, i = args
            out = {f'inbox/{j}/{i}': read(f'r/{j}/{i}')}
        elif op == 'accumulate':
            j, i = args
            q = 2*read(f'inbox/{j}/{i}')
            out = {'acc': read('acc')+read(f'm/{i}')*read(f'g/{i}')*q}
            if q != values[f'u/{j%2}/{i}']:
                mismatches += 1
        elif op == 'public':
            j, = args
            out = {f'D/{j}': read('acc')+read(f'nu/{j}'), f't/{j}': read('clock')}
            need(out[f't/{j}'] == j*TAU, 'model clock/action correspondence')
            outputs[j] = out[f'D/{j}']
            snapshots[j] = [values[f'u/{j%2}/{i}'] for i in range(64)]
        else:
            raise ValueError('unknown operation')
        need(set(event['reads']) == used, 'complete exact consumed read set')
        equal(event['writes'], {k: x.encode() for k, x in out.items()}, 'operation output')
        for k, x in out.items():
            values[k], writers[k] = x, event_id
            for encoded in x.encode():
                rational = F(encoded)
                max_bits = max(max_bits, abs(rational.numerator).bit_length(), rational.denominator.bit_length())
        peak_slots = max(peak_slots, len(values))
        counts[op] += 1
        read_count += len(used)
        write_count += len(out)
    equal(trace['final_hash'], chain, 'complete chain endpoint')
    equal(trace['detector_intervals'], {str(j): enclosure(x) for j, x in outputs.items()}, 'detector enclosures')
    need(max_bits <= 4096, 'declared exact register coefficient capacity')
    need((mismatches > 0) == (case == 'disabled_feedback'), 'restoration semantics')
    if case != 'disabled_feedback':
        need(restored_reads > 0, 'future action consumes feedback writers')
    # Compare every full restored state, not merely its public contraction.
    for j, actual in snapshots.items():
        reference = [amplitude*parse(x) for x in base['traces']['ascending_intervention']['layers'][j+1]]
        difference = [x-y for x, y in zip(actual, reference, strict=True)]
        norm2 = sum((w*x*x for w, x in zip(m, difference, strict=True)), R())
        if noisy:
            bound = j*TAU*PREP+EPS*j*(j+1)
            need((norm2-bound*bound).sign() <= 0, 'generic propagated perturbation bound')
        elif case != 'disabled_feedback':
            need(norm2 == R(), 'exact same-history action reconstruction')
    return outputs, {'events': len(program), 'reads': read_count, 'writes': write_count,
        'operation_counts': dict(counts), 'restored_writer_reads_in_future_action': restored_reads,
        'retained_scalar_slots': peak_slots, 'maximum_rational_coefficient_bits': max_bits,
        'declared_rational_coefficient_capacity_bits': 4096,
        'reconstruction_mismatches': mismatches}


def verify(packet, ctx=None):
    ctx = context() if ctx is None else ctx
    m, a, v, g, identity, base, parent_report = ctx
    need(set(packet) == {'schema', 'parent', 'model', 'steps', 'cases', 'contract'}, 'packet schema')
    equal(packet['schema'], 'oph.common-history-feedback.v1', 'schema')
    equal(packet['parent'], {'path': PARENT, 'sha256': PARENT_SHA}, 'parent identity')
    equal(packet['model'], identity, 'same carrier/action/preparation/readout')
    equal(packet['steps'], STEPS, 'finite duration')
    equal(packet['contract'], CONTRACT, 'scientific contract and error budget')
    equal([x['case'] for x in packet['cases']], CASES, 'complete control cases')
    norms = [sum((mi*x*x for mi, x in zip(m, w, strict=True)), R()) for w in (v, g)]
    need((norms[0]-F(9,16)).sign() <= 0 and (norms[1]-F(1,16)).sign() <= 0, 'sharpened velocity/detector norms')
    need((sum(m, R())-1).sign() <= 0 and (TAU-F(8,125)).sign() < 0, 'mass and clock bounds')
    source = [i for i, x in enumerate(v) if x.sign()]
    detector = [i for i, x in enumerate(g) if x.sign()]
    need(not set(source)&set(detector), 'disjoint source and detector')
    for s in source:
        for d in detector:
            distance = parse(identity['addresses'][d]['coordinate_Qphi'][0])-parse(identity['addresses'][s]['coordinate_Qphi'][0])
            need((distance-F(3,8)).sign() > 0, 'finite separated support')
    outputs, resources = {}, {}
    for trace in packet['cases']:
        outputs[trace['case']], resources[trace['case']] = replay(trace, ctx)
    for j in range(1, STEPS+1):
        need(outputs['baseline'][j] == R(), 'zero baseline')
        need(outputs['double_input'][j] == 2*outputs['intervention'][j], 'altered-input linear response')
    need(outputs['disabled_feedback'][17] != outputs['intervention'][17], 'disabled feedback changes the history')
    # Generic two-record comparison: both preparations and both feedback
    # histories may carry any errors inside the declared pointwise budgets.
    time_upper = F(8,125)
    feedback = F(1,4)*2*EPS*17*18
    preparation = F(1,4)*2*17*time_upper*PREP
    detector_error = 2*DETECT
    measurement = feedback+preparation+detector_error
    discretization = F(1,4)*max(F(x['full_mass_norm_discretization_error'][1])
        for x in base['diagnostics']['comparisons'] if 15 <= x['step'] <= 17)
    # The underlying same finite continuous solution satisfies
    # |D'(t)| <= ||g|| ||v|| <= 3/16. Causal holding uses no future record.
    hold_error, clock_error = F(3,16)*time_upper, F(3,16)*CLOCK
    total = measurement+discretization+hold_error+clock_error
    ideal_lower = min(-F(enclosure(outputs['intervention'][j])[1]) for j in (15,16,17))
    response_lower = ideal_lower-measurement
    need(0 < total < response_lower, 'total error smaller than admitted observer response')
    for j in (15,16,17):
        delta = outputs['noisy_intervention'][j]-outputs['noisy_baseline'][j]-outputs['intervention'][j]
        need((delta-R(measurement)).sign() <= 0 and (delta+R(measurement)).sign() >= 0, 'paired noisy readback bound')
    return {'schema': 'oph.common-history-feedback.verification.v1', 'verdict': 'PASS',
        'parent_full_independent_replay': parent_report is not None,
        'model_sha256': digest(identity), 'sites': 64, 'source_sites': len(source), 'detector_sites': len(detector),
        'protected_model_serialized_bytes': len(parent.canonical(identity)),
        'separation_lower': '3/8', 'window_Qphi': [(j*TAU).encode() for j in (15,17)],
        'window_intervals': [enclosure(j*TAU) for j in (15,17)],
        'paired_response_lower': str(response_lower), 'comparison_total_error_upper': str(total),
        'strict_margin': str(response_lower-total), 'continuous_signal_lower': str(response_lower-total),
        'error_budget': {'feedback_accumulation': str(feedback), 'preparation': str(preparation),
          'final_detector': str(detector_error), 'finite_action_discretization': str(discretization),
          'causal_hold_reconstruction': str(hold_error), 'clock_readback': str(clock_error),
          'exact_arithmetic_on_declared_coefficients': '0', 'finite_support_localization': '0'},
        'resources': resources,
        'disabled_feedback_final_detector': enclosure(outputs['disabled_feedback'][17]),
        'scope': {'classical_software_observer': True, 'same_history_feedback_consumed': True,
          'uniform_adversarial_bounded_record_errors': True, 'continuous_same_finite_action_reference': True,
          'finite_spatial_continuum_error_claim': False, 'native_law_or_physical_clock_selected': False,
          'quantum_outcomes_or_nature_identification': False}}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('path', nargs='?', type=Path, default=PATH)
    ap.add_argument('--write', type=Path)
    args = ap.parse_args()
    if args.write:
        need(args.write.resolve() == EVIDENCE/'verification.json', 'canonical verification output')
    else:
        check_manifest()
    result = verify(load(args.path))
    result['packet_sha256'] = sha256(args.path.read_bytes()).hexdigest()
    if args.write:
        args.write.write_bytes(parent.canonical(result))
        (EVIDENCE/'manifest.json').write_bytes(parent.canonical(manifest()))
    else:
        equal(json.loads((EVIDENCE/'verification.json').read_text()), result,
              'independently recomputed retained verification summary')
    print(json.dumps({k: result[k] for k in ('verdict', 'parent_full_independent_replay',
        'paired_response_lower', 'comparison_total_error_upper', 'strict_margin')}, indent=2))
