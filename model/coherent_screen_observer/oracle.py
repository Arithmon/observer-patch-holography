"""Small full-configuration native oracle and independent runtime control.

Ports 0 and 1 are the observed pair. The minor tensor factor is an untouched
reference qubit. The stopping cap counts external native attempts; only the
incident-edge mask, or timeout, is exposed. No dual evolution is used here.
"""
from collections import defaultdict
from fractions import Fraction as F
import hashlib
from itertools import product
from math import gcd, sqrt

import numpy as np

from ..primitive_observer import runtime


RECORDS = (('timeout', ()), *tuple(('timeout', (m,)) for m in range(4)),
           *tuple(('complete', word) for word in product(range(4), repeat=2)))


def _parameters(n, edges, root, incident, law, cap, feedback, rho):
    if type(n) is not int or not 2 <= n <= 5:
        raise ValueError('two through five original integer ports required')
    if (type(edges) is not tuple or not edges or any(type(e) is not tuple or len(e) != 2
            or any(type(v) is not int for v in e) or not 0 <= e[0] < e[1] < n for e in edges)
            or tuple(sorted(set(edges))) != edges):
        raise ValueError('a sorted unique tuple of oriented graph edges is required')
    if (type(root) is not tuple or any(type(v) is not int or not 0 <= v < n for v in root)
            or tuple(sorted(set(root))) != root or not {0, 1} <= set(root)):
        raise ValueError('the clock region must contain observed ports zero and one')
    if incident is None:
        incident = tuple(e for e in edges if set(e) & set(root))
    if (type(incident) is not tuple or not incident
            or any(type(e) is not tuple or len(e) != 2 or any(type(v) is not int for v in e)
                   for e in incident)
            or tuple(sorted(set(incident))) != incident or not set(incident) <= set(edges)
            or any((0 in e or 1 in e) and e not in incident for e in edges)):
        raise ValueError('incident edges must be a complete clock menu covering the observed pair')
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('an original C or E law is required')
    if type(cap) is not int or not 1 <= cap <= 4 or type(feedback) is not bool:
        raise ValueError('an original cap from one through four and Boolean feedback required')
    if (type(rho) is not F or type(rho.numerator) is not int or type(rho.denominator) is not int
            or rho.denominator <= 0 or gcd(rho.numerator, rho.denominator) != 1):
        raise ValueError('a canonical rational Bernoulli population is required')
    rho = F(rho.numerator, rho.denominator)
    if not 0 <= rho <= 1:
        raise ValueError('Bernoulli population must lie in [0,1]')
    return incident, rho


def _image(mask, edge, coin):
    a, b = edge
    left, right = (mask >> a) & 1, (mask >> b) & 1
    total = left+right
    out_a = (total+coin)//2
    out_b = total-out_a
    output = (mask & ~((1 << a) | (1 << b))) | (out_a << a) | (out_b << b)
    return output, left-out_a, left-right


def _initial(n, rho):
    result = {}
    for remote in range(1 << (n-2)):
        count = remote.bit_count()
        mass = rho**count*(1-rho)**(n-2-count)/2
        a, b = 2*((remote << 2) | 1), 2*((remote << 2) | 2)+1
        for left, right in product((a, b), repeat=2):
            result[left, right] = mass
    return {key: value for key, value in result.items() if value}


def _stopped(rho, edges, incident, law, cap):
    output = {m: defaultdict(F) for m in range(4)}
    pending = rho
    weight = F(1, 2*len(edges))
    for _ in range(cap):
        remote = defaultdict(F)
        for edge in edges:
            destination = output[int(0 in edge)+2*int(1 in edge)] if edge in incident else remote
            for coin in (0, 1):
                for (left, right), value in pending.items():
                    x, y = _image(left//2, edge, coin), _image(right//2, edge, coin)
                    if x[1 if law == 'C' else 2] == y[1 if law == 'C' else 2]:
                        destination[2*x[0]+left % 2, 2*y[0]+right % 2] += weight*value
        pending = {key: value for key, value in remote.items() if value}
    return {m: dict(block) for m, block in output.items()}, pending


def _full_history(n, edges, incident, law, cap, feedback, rho):
    first, timeout = _stopped(_initial(n, rho), edges, incident, law, cap)
    result = {('timeout', ()): timeout}
    for mask, state in first.items():
        if feedback and mask == 1:
            state = {(a, b): value*(-1)**(((a//2) & 1)+((b//2) & 1))
                     for (a, b), value in state.items()}
        second, timeout = _stopped(state, edges, incident, law, cap)
        result['timeout', (mask,)] = timeout
        result.update({('complete', (mask, following)): block for following, block in second.items()})
    return {record: result[record] for record in RECORDS}


def _reduced(block):
    result = np.full((8, 8), F(0), dtype=object)
    for (a, b), value in block.items():
        if a//8 == b//8:
            result[a % 8, b % 8] += value
    return result


def exact_history(n, edges, root=(0, 1), *, incident=None, law='C', cap=3,
                  feedback=True, rho=F(1, 3)):
    """All 21 exact pair/reference outputs; explicit incident overrides root."""
    incident, rho = _parameters(n, edges, root, incident, law, cap, feedback, rho)
    return {record: _reduced(block) for record, block in
            _full_history(n, edges, incident, law, cap, feedback, rho).items()}


def _native_patterns(n, edges, incident, law):
    """Independent original Kraus grouping on every binary basis column."""
    rows = []
    for edge in edges:
        for coin in (0, 1):
            groups = {}
            for mask in range(1 << n):
                target, current, contrast = _image(mask, edge, coin)
                key = current if law == 'C' else contrast
                groups.setdefault(key, [-1]*(1 << n))[mask] = target
            label = int(0 in edge)+2*int(1 in edge) if edge in incident else None
            rows.extend((label, tuple(pattern)) for pattern in groups.values())
    return rows


def _stopped_patterns(n, edges, incident, law, cap):
    # Compositions remain partial permutations. Merge identical Kraus maps
    # by adding their exact weights, never by coherently adding alternatives.
    pending = {tuple(range(1 << n)): F(1)}
    output = {m: defaultdict(F) for m in range(5)}
    native = _native_patterns(n, edges, incident, law)
    for _ in range(cap):
        remote = defaultdict(F)
        for before, weight in pending.items():
            for label, after in native:
                combined = tuple(after[i] if i >= 0 else -1 for i in before)
                if any(i >= 0 for i in combined):
                    (remote if label is None else output[label])[combined] += weight/(2*len(edges))
        pending = remote
    output[4] = pending
    return output


def _matrix(block, size):
    result = np.zeros((size, size), complex)
    for (a, b), value in block.items():
        result[a, b] = float(value)
    return result


def _distance(first, second):
    return sum(float(np.linalg.svd(first[key]-second[key], compute_uv=False).sum())/2 for key in first)


def execute_history(n, edges, root=(0, 1), *, incident=None, law='C', cap=3,
                    feedback=True, rho=F(1, 3)):
    """Execute the native stopped instrument and replay its retained checkpoint.

    Numerical Kraus amplitudes are square roots of exact native path weights.
    This is CP-map semantics, not an emitted physical pulse implementation.
    Reported residuals are diagnostics rather than roundoff certificates.
    """
    incident, rho = _parameters(n, edges, root, incident, law, cap, feedback, rho)
    size = 1 << n
    patterns = _stopped_patterns(n, edges, incident, law, cap)
    groups = []
    for rows in patterns.values():
        group = []
        for pattern, weight in rows.items():
            matrix = np.zeros((size, size), complex)
            for column, row in enumerate(pattern):
                if row >= 0:
                    matrix[row, column] = sqrt(float(weight))
            group.append(matrix)
        groups.append(tuple(group) or (np.zeros((size, size), complex),))
    if any(len(group) > 256 for group in groups):
        raise ValueError('this small runtime control exceeds its explicit Kraus allocation')
    cost = (('native_attempts_upper', cap), ('stopped_intervals', 1))
    binding = {'n': n, 'edges': edges, 'incident': incident, 'law': law,
               'cap': cap, 'feedback': feedback, 'rho': str(rho)}
    ordinary = runtime.Instrument(tuple(groups), True, resources=cost, binding=binding)
    phase = np.diag([(-1)**(mask & 1) for mask in range(size)])
    controlled = runtime.Instrument(tuple(tuple(k@phase for k in group) for group in groups),
                                    True, resources=cost+(('logical_phase_controls', 1),),
                                    binding={**binding, 'control': 'Z0'})
    idle = runtime.Instrument(tuple((np.eye(size) if m == 4 else np.zeros((size, size)),)
                                    for m in range(5)), True)
    initial = _matrix(_initial(n, rho), 2*size)
    work = np.trace(initial.reshape(size, 2, size, 2), axis1=1, axis2=3)
    program = runtime.Program('coherent stopped native control', size, work,
        (runtime.Event('first stopped interval', (ordinary,)),
         runtime.Event('second stopped interval', tuple(idle if m == 4 else
             controlled if feedback and m == 1 else ordinary for m in range(5)), control=0)),
        binding=binding)
    whole = runtime.execute(program, initial=initial)
    prefix = runtime.execute(program, initial=initial, stop_after=1)
    payload = runtime.checkpoint_to_json(prefix)
    resumed = runtime.execute(program, checkpoint=runtime.checkpoint_from_json(program, payload))
    replay = runtime.complete_distance(whole, resumed)
    if replay != 0 or whole.resources != resumed.resources:
        raise ArithmeticError('serialized continuation changed the complete runtime output')
    full, reduced, padding = {}, {}, 0.
    for branch in whole.branches:
        a, b = branch.record
        if a == 4 and b != 4:
            padding += float(np.linalg.norm(branch.density, 'nuc'))
            continue
        record = ('timeout', ()) if a == 4 else ('timeout', (a,)) if b == 4 else ('complete', (a, b))
        full[record] = branch.density
        reduced[record] = np.trace(branch.density.reshape(1 << (n-2), 8, 1 << (n-2), 8),
                                   axis1=0, axis2=2)
    exact = _full_history(n, edges, incident, law, cap, feedback, rho)
    error = _distance(full, {key: _matrix(value, 2*size) for key, value in exact.items()})
    reduced_error = _distance(reduced, {key: np.asarray(_reduced(value), dtype=complex)
                                       for key, value in exact.items()})
    if padding != 0 or error > 2e-12 or abs(whole.probability_mass-1) > 2e-12:
        raise ArithmeticError('native runtime differs from the full exact stopped history')
    q = F(len(edges)-len(incident), len(edges))
    resource_names = {name for branch in whole.branches for name, _ in branch.resources}
    maximum_resources = {name: max(dict(branch.resources).get(name, 0) for branch in whole.branches)
                         for name in sorted(resource_names)}
    return {record: reduced[record] for record in RECORDS}, {
        'complete_runtime_records': 25, 'comparison_records': 21,
        'zero_padding_nuclear_norm': padding, 'checkpoint_distance': replay,
        'checkpoint_bytes': len(payload.encode()), 'checkpoint_json': payload,
        'checkpoint_sha256': hashlib.sha256(payload.encode()).hexdigest(),
        'full_output_distance_numeric': error,
        'reduced_output_distance_numeric': reduced_error,
        'probability_mass': whole.probability_mass,
        'maximum_native_attempts': 2*cap, 'declared_joint_preparations': 1,
        'expected_native_attempts': sum((q**i for i in range(cap)), F(0))*(2-q**cap),
        'timeout_probability': 1-(1-q**cap)**2,
        'maximum_branch_resources': maximum_resources,
        'expanded_tree_resources': dict(whole.resources),
        'kraus_per_outcome': tuple(map(len, groups)), 'physical_pulse_execution': False}


def run():
    """Bounded source-checked runtime summaries suitable for a JSON receipt."""
    from ..observer import native
    checked = 0
    for bits, coin in product(range(4), (0, 1)):
        loads = tuple((bits >> i) & 1 for i in range(12))
        image, record = native.native(loads, (0, 1), coin)
        expected = (sum(value << i for i, value in enumerate(image)),
                    record[1], loads[0]-loads[1])
        if _image(bits, (0, 1), coin) != expected:
            raise ArithmeticError('binary oracle differs from the original native endpoint rule')
        checked += 1
    edges = ((0, 1), (0, 2), (1, 2), (2, 3), (3, 4))
    rows = []
    for law in ('C', 'E'):
        _, diagnostics = execute_history(5, edges, root=(0, 1, 2), law=law, cap=2)
        for key in ('checkpoint_json', 'checkpoint_bytes', 'checkpoint_sha256'):
            diagnostics.pop(key)
        diagnostics['checkpoint_replayed'] = True
        diagnostics['kraus_per_outcome'] = list(diagnostics['kraus_per_outcome'])
        rows.append({'law': law, 'ports': 5, 'cap_per_interval': 2,
                     **{key: str(value) if type(value) is F else value
                        for key, value in diagnostics.items()}})
    return {'original_native_endpoint_checks': checked, 'executions': rows}
