"""Independent complete-event and clock-budget replay of the RG reference.

Imports no producer. It validates the generated primitive process, not the
additional assertion that a physical source implements this process theory.
"""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json


def require(test, message):
    if not test:
        raise ValueError(message)


def rational(value):
    require(type(value) is str and len(value) < 128, "rational type/size")
    result = Q(value)
    require(str(result) == value, "canonical rational")
    return result


def check(q, layers, intervention, events):
    require(type(q) is int and 2 <= q <= 9, "cutoff")
    require(type(layers) is int and 1 <= layers <= 4, "layer count")
    count = q**3
    require(intervention is None or type(intervention) is int and 0 <= intervention < count,
            "intervention")
    # Independent coordinate decoder, and an integer budget inequality.
    xyz = lambda i: (i//(q*q), (i//q) % q, i % q)
    squared_integer = lambda i, j: sum((a-b)**2 for a, b in zip(xyz(i), xyz(j)))
    allowed = {(i, j) for i in range(count) for j in range(count)
               if (2*q+1)*squared_integer(i, j) <= 2*q*q}
    offdiag = {pair for pair in allowed if pair[0] != pair[1]}
    incoming_count = Counter(j for i, j in allowed)
    # A finite positive timing margin exists: equality would require the
    # odd integer 2q+1, coprime to q, to divide 2.
    require(all((2*q+1)*squared_integer(i, j) != 2*q*q
                for i in range(count) for j in range(count)), "clock-boundary clearance")
    store, wires, delivered, waited, consumed = {}, {}, set(), set(), set()
    sums = [0]*count
    commits = set()
    read_count = Counter()
    layer = 0
    phase = -1
    last_arrival = -1
    counters = Counter()
    digest = hashlib.sha256()
    inputs = []
    responses = []
    parents = {}
    for event in events:
        require(type(event) is list and event and type(event[0]) is str, "event schema")
        name = event[0]
        lengths = {"prepare": 4, "fork": 5, "flight": 6, "wait": 6,
                   "accumulate": 7, "commit": 4, "checkpoint": 3}
        require(name in lengths and len(event) == lengths[name], "event type/length")
        # Numbers in the event language are actual integers, never booleans.
        integer_fields = {"prepare": [1, 3], "fork": [1, 2, 3, 4],
                          "flight": [1, 2, 3, 5], "wait": [1, 2, 3],
                          "accumulate": [1, 2, 3, 4, 5, 6], "commit": [1, 2, 3],
                          "checkpoint": [1]}
        require(all(type(event[i]) is int for i in integer_fields[name]), "integer event field")
        require(layer < layers or name == "prepare", "events after final checkpoint")
        if name == "prepare":
            _, i, address, value = event
            require(phase == -1 and layer == 0 and i == len(inputs) and i < count, "preparation coverage")
            require(type(address) is list and all(type(a) is int for a in address) and
                    address == list(xyz(i)), "source address")
            require(value == i+1+int(intervention == i), "source intervention")
            inputs.append(value)
            store[-1, i] = value
        else:
            require(len(inputs) == count and event[1] == layer, "stage lineage")
            next_phase = {"fork": 0, "flight": 1, "wait": 2, "accumulate": 3,
                          "commit": 3, "checkpoint": 4}[name]
            require(phase <= next_phase, "noncausal operation ordering")
            phase = next_phase
            if name in {"fork", "flight", "wait"}:
                i, j = event[2:4]
                pair = i, j
                require(pair in offdiag, "unsupported clock-budget flight")
                if name == "fork":
                    require(pair not in wires and event[4] == store[layer-1, i], "immutable writer fork")
                    wires[pair] = event[4]
                elif name == "flight":
                    require(len(wires) == len(offdiag) and pair not in delivered, "complete distinct forks/flights")
                    require(event[5] == wires[pair], "flight payload")
                    squared = rational(event[4])
                    require(squared == Q(squared_integer(i, j), q*q), "unit-flight duration")
                    require(squared >= last_arrival, "flight arrival chronology")
                    last_arrival = squared
                    delivered.add(pair)
                else:
                    require(len(delivered) == len(offdiag) and pair not in waited, "arrival before wait")
                    budget, used = map(rational, event[4:6])
                    require(budget == Q(2, 2*q+1) and
                            used == Q(squared_integer(i, j), q*q) and used <= budget,
                            "wait clock/nonnegativity")
                    waited.add(pair)
            elif name == "accumulate":
                _, _, i, j, sample, before, after = event
                require(len(waited) == len(offdiag) and (i, j) in allowed and (i, j) not in consumed,
                        "complete deadline and read coverage")
                require(j not in commits and sample == store[layer-1, i], "read version/value")
                if i != j:
                    require(sample == wires[i, j], "receiver-local delivered copy")
                require(before == sums[j] and after == before+sample, "local accumulation")
                sums[j] = after
                consumed.add((i, j))
                read_count[j] += 1
                parents.setdefault((layer, j), []).append((layer-1, i))
            elif name == "commit":
                _, _, j, value = event
                require(0 <= j < count and j not in commits and value == sums[j], "commit value")
                require(read_count[j] == incoming_count[j],
                        "complete inferred parent menu")
                store[layer, j] = value
                commits.add(j)
            else:
                require(type(event[2]) is list and len(event[2]) == count and
                        all(type(v) is int for v in event[2]), "checkpoint schema")
                require(commits == set(range(count)) and consumed == allowed, "whole-layer completion")
                require(event[2] == [store[layer, i] for i in range(count)], "checkpoint values")
                responses.append(event[2])
                layer += 1
                phase, last_arrival = -1, -1
                wires, delivered, waited, consumed, commits = {}, set(), set(), set(), set()
                sums = [0]*count
                read_count = Counter()
        digest.update((json.dumps(event, separators=(",", ":"))+"\n").encode())
        counters[name] += 1
    require(layer == layers, "truncated process")
    require(counters == {"prepare": count, "fork": layers*len(offdiag),
                         "flight": layers*len(offdiag), "wait": layers*len(offdiag),
                         "accumulate": layers*len(allowed), "commit": layers*count,
                         "checkpoint": layers}, "whole-operation census")
    return {"sha256": digest.hexdigest(), "events": sum(counters.values()),
            "operation_counts": dict(sorted(counters.items())),
            "reads_per_layer": len(allowed), "source_population": count,
            "minimum_squared_clock_margin": str(min(Q(2, 2*q+1)-Q(squared_integer(i, j), q*q)
                                                       for i, j in allowed)),
            "layers": responses, "parents": parents}
