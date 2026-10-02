"""Indexed native schedule; finite controls execute the unpowered graph."""

from .format import digest


def dimensions(m, power):
    if type(m) is not int or m < 2 or type(power) is not int or not 1 <= power <= 12:
        raise ValueError('integer graph size and bounded power')
    n, d = m*m, 8**power
    return n, d, n*(d+2), 2+d+d*(d-1)//2


def destination(m, source, word, power, inverse=False):
    labels = [(word // 8**j) % 8 for j in range(power)]
    if inverse:
        labels = [v ^ 1 for v in reversed(labels)]
    x, y = divmod(source, m)
    for label in labels:
        offset, sign = (label % 4)//2, 2*(label % 2)-1
        if label < 4:
            x = (x+sign*(2*y+offset)) % m
        else:
            y = (y+sign*(2*x+offset)) % m
    return m*x+y


def comparator_position(d, t):
    # Invert the triangular row count without enumerating D^2 comparators.
    lo, hi = 0, d-1
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid*(2*d-mid-1)//2 <= t:
            lo = mid
        else:
            hi = mid
    return t-lo*(2*d-lo-1)//2


def operation(m, power, layer, wire):
    """The service touching one wire, with its complete operand list.

    Random access supports the full degree 8^12 without materializing its
    enormous circuit. This is a schedule compiler, not an execution claim.
    """
    n, d, width, depth = dimensions(m, power)
    if (type(layer) is not int or not 0 <= layer < depth
            or type(wire) is not int or not 0 <= wire < width):
        raise ValueError('layer or wire outside schedule')
    if layer == 0:
        return ['idle' if wire < n else 'reset', wire]
    if 1 <= layer <= d:
        word = layer-1
        if wire < n:
            voter = destination(m, wire, word, power)
            return ['copy', wire, n+d*voter+word]
        voter, position = divmod(wire-n, d)
        if voter < n and position == word:
            source = destination(m, voter, word, power, inverse=True)
            return ['copy', source, wire]
        return ['idle', wire]
    if layer == depth-1:
        if wire >= n*(d+1):
            voter = wire-n*(d+1)
            return ['copy', n+d*voter+d//2-1, wire]
        if n <= wire < n*(d+1):
            voter, position = divmod(wire-n, d)
            if position == d//2-1:
                return ['copy', wire, n*(d+1)+voter]
        return ['idle', wire]
    position = comparator_position(d, layer-d-1)
    if n <= wire < n*(d+1):
        voter, j = divmod(wire-n, d)
        if j in (position, position+1):
            return ['sort', n+d*voter+position, n+d*voter+position+1]
    return ['idle', wire]


def schedule(m):
    if type(m) is not int or m not in (2, 3):
        raise ValueError('full finite execution is limited to m=2,3 and power=1')
    _, _, width, depth = dimensions(m, 1)
    layers = []
    for t in range(depth):
        covered, layer = set(), []
        for wire in range(width):
            if wire not in covered:
                op = operation(m, 1, t, wire)
                covered.update(op[1:])
                layer.append(op)
        layers.append(layer)
    return dict(m=m, power=1, width=width, layers=layers)


def histogram(words, case_mask):
    bins = [case_mask]+[0]*len(words)
    for word in words:
        new = [0]*len(bins)
        for k, mask in enumerate(bins):
            new[k] |= mask & ~word
            if k+1 < len(bins):
                new[k+1] |= mask & word
        bins = new
    return bins


def census(m):
    program = schedule(m)
    n, d, width, _ = dimensions(m, 1)
    flat = [op for layer in program['layers'] for op in layer]
    cases = sum((1 << (len(op)-1))-1 for op in flat)
    all_cases = (1 << cases)-1
    out_hist, live_hist = [0]*(n+1), [0]*(n+1)
    # Bits of each Python integer are independent fault experiments.
    for logical in (0, 1):
        for initial in [-1]+list(range(n)):
            wires = [0]*width
            for i in range(n):
                wires[i] = all_cases if logical ^ (i == initial) else 0
            maximum = [0]*(n+1)
            maximum[int(initial >= 0)] = all_cases
            cursor = 0
            for op in flat:
                name, a, *rest = op
                if name == 'reset':
                    wires[a] = 0
                elif name == 'copy':
                    wires[rest[0]] ^= wires[a]
                elif name == 'sort':
                    b = rest[0]
                    wires[a], wires[b] = wires[a] & wires[b], wires[a] | wires[b]
                for operand, wire in enumerate(op[1:]):
                    flip = sum(1 << (cursor+mask-1) for mask in range(1, 1 << (len(op)-1))
                               if mask >> operand & 1)
                    wires[wire] ^= flip
                cursor += (1 << (len(op)-1))-1
                if a < n:
                    current = histogram([w ^ (all_cases*logical) for w in wires[:n]], all_cases)
                    new = [0]*(n+1)
                    for i, mask in enumerate(maximum):
                        for j, other in enumerate(current):
                            new[max(i, j)] |= mask & other
                    maximum = new
            outputs = wires[n*(d+1):]
            for i, bits in enumerate(histogram([x ^ (all_cases*logical) for x in outputs], all_cases)):
                out_hist[i] += bits.bit_count()
            for i, bits in enumerate(maximum):
                live_hist[i] += bits.bit_count()
    return dict(m=m, power=1, schedule_sha256=digest(program), width=width,
                depth=len(program['layers']), locations=len(flat),
                cases=2*(n+1)*cases, output_histogram=out_hist, live_histogram=live_hist)


def indexed_cases():
    rows = []
    m, power = 17, 12
    n, d, width, depth = dimensions(m, power)
    points = [(0, 0), (0, n), (1, 13), (d//2, 211), (d, n-1),
              (d+1, n), (d+2, n+1), (depth-2, n),
              (depth-1, n+d//2-1), (depth-1, width-1)]
    for t, w in points:
        rows.append(dict(m=m, power=power, layer=t, wire=w,
                         operation=operation(m, power, t, w)))
    return rows


def candidate():
    return dict(censuses=[census(m) for m in (2, 3)], indexed=indexed_cases())
