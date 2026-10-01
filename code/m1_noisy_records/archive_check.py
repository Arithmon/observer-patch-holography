"""Independent Boolean truth-table replay and arbitrary-r structural checks."""

import math
from .independent import keys, need, integer


def validate(program):
    keys(program, 'length width outputs layers')
    n = program['length']
    integer(n, 5, 65)
    need((n-1) % 4 == 0 and type(program['width']) is int
         and program['width'] == n*(n+2), 'refresh dimensions')
    need(type(program['outputs']) is list
         and all(type(q) is int for q in program['outputs'])
         and program['outputs'] == list(range(n*(n+1), n*(n+2))), 'output wires')
    need(type(program['layers']) is list and len(program['layers']) == 2+n+n*(n-1)//2,
         'complete refresh depth')
    flat, copied, sorts = [], set(), [[] for _ in range(n)]
    for t, layer in enumerate(program['layers']):
        need(type(layer) is list and layer, 'nonempty layer')
        touched = set()
        for op in layer:
            need(type(op) is list and len(op) in (2, 3), 'instruction')
            name, args = op[0], op[1:]
            need(name in ('reset', 'copy', 'sort', 'idle') and
                 len(args) == (2 if name in ('copy', 'sort') else 1), 'operation arity')
            for q in args:
                integer(q, 0, program['width']-1)
                need(q not in touched, 'one service per wire per layer')
                touched.add(q)
            if t == 0:
                need(name == ('idle' if args[0] < n else 'reset'), 'blank supply and initial idles')
            elif t <= n and name != 'idle':
                need(name == 'copy' and args[0] < n and n <= args[1] < n*(n+1), 'private input copy')
                voter, position = divmod(args[1]-n, n)
                need(position == args[0] and (voter, position) not in copied, 'one copy per voter/input')
                copied.add((voter, position))
            elif t == len(program['layers'])-1 and name != 'idle':
                need(name == 'copy' and args[1] in program['outputs'], 'output copy')
                voter = args[1]-n*(n+1)
                need(args[0] == n+n*voter+n//2, 'majority position')
            elif name != 'idle':
                need(name == 'sort' and n <= args[0] < args[1] < n*(n+1), 'private comparator')
                v, i = divmod(args[0]-n, n)
                need(args[1] == args[0]+1 and i+1 < n, 'same private lane')
                sorts[v].append(i)
            flat.append(op)
        need(touched == set(range(program['width'])), 'every live wire, including idles')
    need(len(copied) == n*n, 'all voter inputs')
    ideal_sorts = [i for end in range(n-1, 0, -1) for i in range(end)]
    need(all(s == ideal_sorts for s in sorts), 'complete sorting networks')
    # Check every output has exactly one non-idle writer in the final layer.
    need(sum(op[0] == 'copy' for op in program['layers'][-1]) == n, 'all outputs written')
    return flat


def run(program, inputs, faults=()):
    flat = validate(program)
    need(type(inputs) is list and len(inputs) == program['length'] and
         all(type(x) is int and x in (0, 1) for x in inputs), 'binary input')
    need(type(faults) in (tuple, list), 'fault list')
    by_location = {}
    for pair in faults:
        need(type(pair) in (list, tuple) and len(pair) == 2, 'fault pair')
        i, mask = pair
        integer(i, 0, len(flat)-1)
        integer(mask, 1, (1 << (len(flat[i])-1))-1)
        need(i not in by_location, 'duplicate fault location')
        by_location[i] = mask
    state = inputs+[0]*(program['width']-len(inputs))
    for j, op in enumerate(flat):
        name, a, *rest = op
        if name == 'reset':
            state[a] = 0
        elif name == 'copy':
            b = rest[0]
            state[b] = (state[a]+state[b]) % 2
        elif name == 'sort':
            b = rest[0]
            state[a], state[b] = min(state[a], state[b]), max(state[a], state[b])
        mask = by_location.get(j, 0)
        for bit, q in enumerate(op[1:]):
            state[q] = (state[q]+(mask >> bit & 1)) % 2
    return [state[q] for q in program['outputs']]


def verify(row):
    keys(row, 'program census')
    p = row['program']
    ops = validate(p)
    need(p['length'] == 5, 'frozen exhaustive census size')
    keys(row['census'], 'locations single_fault_cases output_error_histogram')
    histogram = [0]*6
    for logical in (0, 1):
        for error in range(-1, 5):
            bits = [logical ^ (i == error) for i in range(5)]
            need(run(p, bits) == [logical]*5, 'no-fault recovery')
            for j, op in enumerate(ops):
                for mask in range(1, 1 << (len(op)-1)):
                    result = run(p, bits, [(j, mask)])
                    histogram[sum(x != logical for x in result)] += 1
    need(not any(histogram[2:]), 'one old error plus one fresh fault must remain correctable')
    expected = dict(locations=len(ops), single_fault_cases=sum(histogram),
                    output_error_histogram=histogram)
    need(row['census'] == expected and all(type(row['census'][k]) is int
         for k in ('locations', 'single_fault_cases')) and
         all(type(x) is int for x in row['census']['output_error_histogram']), 'exact fault census')


def log_failure_majorant(q, rate_constant=1.):
    """Log upper bound only: no invented full-library threshold constant."""
    need(type(q) is int and q >= 2, 'integer refinement q >= 2')
    need(type(rate_constant) in (int, float) and math.isfinite(rate_constant)
         and rate_constant > 0, 'positive finite illustrative constant')
    r = (2*q-1).bit_length()  # ceil(log2(2q)), without floating rounding.
    n = 4*r+1
    width = n*(n+2)
    depth = n*(n-1)//2+n+2
    volume = width*depth  # bounds primitive service locations, idles included
    # Union over q^5 log^8(2q) records/cycles; constants fixed in the proof.
    return 5*math.log(q)+8*math.log(math.log(2*q))+(r+1)*math.log(volume*rate_constant/q)
