"""Explicit noisy classical archive refresh, including all live-wire idles."""


def schedule(n):
    if type(n) is not int or n < 5 or (n-1) % 4:
        raise ValueError('archive length must be 4r+1, r >= 1')
    width = n*(n+2)
    layers = [[['reset', j] for j in range(n, width)]]
    # Each source bit is copied once per round. Every voter receives a private
    # copy of every source; gates never touch two shared source inputs.
    for step in range(n):
        layers.append([['copy', i, n+n*((i+step) % n)+i] for i in range(n)])
    # Bubble-sort each private binary word; the central bit is its majority.
    for end in range(n-1, 0, -1):
        for i in range(end):
            layers.append([['sort', n+n*v+i, n+n*v+i+1] for v in range(n)])
    outputs = list(range(n+n*n, width))
    layers.append([['copy', n+n*v+n//2, outputs[v]] for v in range(n)])
    for layer in layers:
        used = {q for op in layer for q in op[1:]}
        layer.extend(['idle', q] for q in range(width) if q not in used)
    return dict(length=n, width=width, outputs=outputs, layers=layers)


def locations(program):
    return [op for layer in program['layers'] for op in layer]


def run(program, input_mask, fault=None, include_live=False):
    n = program['length']
    if program != schedule(n):
        raise ValueError('run expects the complete generated refresh schedule')
    if (type(program['width']) is not int or any(type(q) is not int for q in program['outputs'])
            or any(type(q) is not int for layer in program['layers'] for op in layer for q in op[1:])):
        raise ValueError('integer wire indices required')
    if type(input_mask) is not int or not 0 <= input_mask < 1 << n:
        raise ValueError('input word outside archive length')
    if type(include_live) is not bool:
        raise ValueError('Boolean live-trace option required')
    ops = locations(program)
    if fault is not None:
        if type(fault) not in (tuple, list) or len(fault) != 2:
            raise ValueError('fault must be a location and nonzero operand mask')
        i, mask = fault
        if (type(i) is not int or not 0 <= i < len(ops) or type(mask) is not int
                or not 1 <= mask < 1 << (len(ops[i])-1)):
            raise ValueError('fault location or mask outside circuit')
    state = input_mask
    word_mask = (1 << n)-1
    label_word = word_mask if input_mask.bit_count() > n//2 else 0
    max_live_wrong = (input_mask ^ label_word).bit_count()
    for index, op in enumerate(ops):
        name, a, *rest = op
        if name == 'reset':
            state &= ~(1 << a)
        elif name == 'copy':
            state ^= ((state >> a & 1) << rest[0])
        elif name == 'sort':
            b = rest[0]
            if state >> a & 1 and not state >> b & 1:
                state ^= (1 << a) | (1 << b)
        if fault is not None and index == fault[0]:
            for k, bit in enumerate(op[1:]):
                state ^= ((fault[1] >> k & 1) << bit)
        if include_live:
            max_live_wrong = max(max_live_wrong, ((state & word_mask) ^ label_word).bit_count())
    output = sum((state >> q & 1) << j for j, q in enumerate(program['outputs']))
    return (output, max_live_wrong) if include_live else output


def census(program):
    n = program['length']
    histogram = [0]*(n+1)
    live_histogram = [0]*(n+1)
    inputs = [0]+[1 << i for i in range(n)]
    for logical in (0, 1):
        target = ((1 << n)-1)*logical
        for error in inputs:
            for loc, op in enumerate(locations(program)):
                for fault in range(1, 1 << (len(op)-1)):
                    output, live_wrong = run(program, target ^ error, (loc, fault), include_live=True)
                    wrong = (output ^ target).bit_count()
                    histogram[wrong] += 1
                    live_histogram[live_wrong] += 1
    return dict(locations=len(locations(program)), single_fault_cases=sum(histogram),
                output_error_histogram=histogram, live_input_error_histogram=live_histogram)


def evidence():
    program = schedule(5)
    return dict(program=program, census=census(program))
