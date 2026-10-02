"""Independent direct schedule construction and vectorized fault experiments."""

from functools import lru_cache
import numpy as np
from .format import keys, need, integer, digest
from .graph_check import expected_maps


def reference_schedule(m):
    n, d = m*m, 8
    width = n*(d+2)
    layers = [[['idle', i] for i in range(n)]+[['reset', i] for i in range(n, width)]]
    for j, permutation in enumerate(expected_maps(m)):
        layers.append([['copy', source, n+d*voter+j] for source, voter in enumerate(permutation)])
    for end in range(d-1, 0, -1):
        for j in range(end):
            layers.append([['sort', n+d*v+j, n+d*v+j+1] for v in range(n)])
    layers.append([['copy', n+d*v+3, n*(d+1)+v] for v in range(n)])
    for t, layer in enumerate(layers):
        used = {wire for op in layer for wire in op[1:]}
        layer.extend(['idle', i] for i in range(width) if i not in used)
        layer.sort(key=lambda op: min(op[1:]))
        need(len([w for op in layer for w in op[1:]]) == width, 'no parallel wire conflict')
    return dict(m=m, power=1, width=width, layers=layers)


def validate(program):
    keys(program, 'm power width layers')
    integer(program['m'], 2, 3)
    integer(program['power'], 1, 1)
    expected = reference_schedule(program['m'])
    need(type(program['layers']) is list and len(program['layers']) == len(expected['layers']),
         'complete layer catalog')
    integer(program['width'], expected['width'], expected['width'])
    for layer, ref in zip(program['layers'], expected['layers']):
        need(type(layer) is list and len(layer) == len(ref), 'all live locations')
        for op, target in zip(layer, ref):
            need(type(op) is list and op == target and all(type(w) is int for w in op[1:]),
                 'complete native service and ownership')
    return [op for layer in program['layers'] for op in layer]


@lru_cache(maxsize=2)
def census(m):
    p = reference_schedule(m)
    flat = validate(p)
    n, width = m*m, p['width']
    fault_rows = [(j, mask) for j, op in enumerate(flat) for mask in range(1, 2**(len(op)-1))]
    faults = np.array(fault_rows, dtype=int)
    output, live = np.zeros(n+1, dtype=int), np.zeros(n+1, dtype=int)
    for label in (0, 1):
        for old_error in [-1]+list(range(n)):
            bits = np.zeros((width, len(fault_rows)), dtype=bool)
            bits[:n] = bool(label)
            if old_error >= 0:
                bits[old_error] ^= True
            counts = np.full(len(fault_rows), int(old_error >= 0))
            maximum = counts.copy()
            for j, op in enumerate(flat):
                name, a, *tail = op
                if name == 'reset':
                    bits[a] = False
                elif name == 'copy':
                    bits[tail[0]] ^= bits[a]
                elif name == 'sort':
                    b = tail[0]
                    both = bits[a] & bits[b]
                    bits[b] |= bits[a]
                    bits[a] = both
                active = np.flatnonzero(faults[:, 0] == j)
                for operand, wire in enumerate(op[1:]):
                    selected = active[(faults[active, 1] >> operand & 1) != 0]
                    if wire < n:
                        counts[selected] += np.where(bits[wire, selected] == bool(label), 1, -1)
                        maximum[selected] = np.maximum(maximum[selected], counts[selected])
                    bits[wire, selected] ^= True
            errors = np.count_nonzero(bits[-n:] != bool(label), axis=0)
            output += np.bincount(errors, minlength=n+1)
            live += np.bincount(maximum, minlength=n+1)
    return dict(m=m, power=1, schedule_sha256=digest(p), width=width, depth=len(p['layers']),
                locations=len(flat), cases=2*(n+1)*len(fault_rows),
                output_histogram=output.tolist(), live_histogram=live.tolist())


def verify_indexed(rows):
    need(type(rows) is list and len(rows) == 10, 'complete full-power index catalog')
    m, n, d = 17, 289, 8**12
    depth, width = 2+d+d*(d-1)//2, n*(d+2)
    points = [(0, 0), (0, n), (1, 13), (d//2, 211), (d, n-1),
              (d+1, n), (d+2, n+1), (depth-2, n),
              (depth-1, n+d//2-1), (depth-1, width-1)]
    maps = expected_maps(m)
    for row, (t, w) in zip(rows, points):
        keys(row, 'm power layer wire operation')
        for key, x in [('m', m), ('power', 12), ('layer', t), ('wire', w)]:
            integer(row[key], x, x)
        if t == 0:
            op = ['idle' if w < n else 'reset', w]
        elif t <= d:
            v, code = w, t-1
            for _ in range(12):
                v, code = maps[code % 8][v], code//8
            op = ['copy', w, n+d*v+t-1]
        elif t == depth-1:
            voter = (w-n)//d if w < n*(d+1) else w-n*(d+1)
            op = ['copy', n+d*voter+d//2-1, n*(d+1)+voter]
        else:
            # Frozen first, second and last comparator positions.
            pos = {d+1: 0, d+2: 1, depth-2: 0}[t]
            op = ['sort', n+pos, n+pos+1]
        need(type(row['operation']) is list and row['operation'] == op
             and all(type(x) is int for x in row['operation'][1:]), 'indexed full-power service')


def verify(row):
    keys(row, 'censuses indexed')
    need(type(row['censuses']) is list and len(row['censuses']) == 2, 'complete native census')
    for supplied, m in zip(row['censuses'], (2, 3)):
        expected = census(m)
        keys(supplied, ' '.join(expected))
        for field, value in expected.items():
            if type(value) is int:
                integer(supplied[field], value, value)
            elif type(value) is list:
                need(type(supplied[field]) is list and supplied[field] == value
                     and all(type(x) is int for x in supplied[field]), 'exhaustive fault '+field)
            else:
                need(supplied[field] == value, 'schedule custody')
    verify_indexed(row['indexed'])
