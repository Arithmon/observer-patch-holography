"""Producer: expand the native write/read controller into complete histories."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
if __package__:
    from . import codec
else:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_native_updates import codec

PORTS = [0,9,7,4,8,11,1,5,14,40,23,38,45,39]
REQUESTS = [2,0,2,1]
LEVELS = [F(-3,2), F(-1,2), F(1,2), F(3,2)]


def program():
    ops = []
    def copy(s, t):
        ops.extend(((2*s,2*t), (2*s+1,2*t+1)))
    copy(0, 3)
    copy(1, 2)
    copy(3, 2)
    ops.append((6,7))
    for record in REQUESTS:
        copy(record, 3)
        for cell in range(3, 6):
            copy(cell, cell+1)
        for _ in range(80):
            ops.append((6,7))
            for cell in range(3, 6):
                copy(cell, cell+1)
    return ops


def history(a, b):
    q = 2**20
    values = [2*q]*14
    values[:4] = [int((2+a)*q), int((2-a)*q), int((2+b)*q), int((2-b)*q)]
    initial = values[:]
    writers = [-i-1 for i in range(14)]
    tape, reads = [], []
    uses = [0,0,0]
    for k, (u, v) in enumerate(program()):
        output = round(F(values[u]+values[v], 2))
        tape.append([u,v,values[u],values[v],writers[u],writers[v],output])
        values[u] = values[v] = output
        writers[u] = writers[v] = k
        if k >= 7 and (k-7) % 568 == 7:
            cycle = (k-7)//568
            record = REQUESTS[cycle]
            uses[record] += 1
            exponent = (1,1,2)[record]+uses[record]+3
            alphabet = LEVELS if record < 2 else list(map(F, range(-3,4)))
            observation = F(values[12]-values[13], 2*q)*2**exponent
            ranked = sorted((abs(observation-c), c) for c in alphabet)
            if ranked[0][0] == ranked[1][0]:
                raise ValueError("ambiguous receiver")
            reads.append({"after_mean": k, "record": record, "use": uses[record],
                          "exponent": exponent, "local_units": values[12:14],
                          "writers": writers[12:14], "decoded": str(ranked[0][1])})
        if k == 6:
            committed = values[:]
    return {"payloads": [str(a),str(b)], "initial_units": initial,
            "commit_units": committed, "reads": reads, "final_units": values,
            "event_count": len(tape), "tape_sha256": codec.digest(tape)}, tape


def build():
    cases, retained = [], None
    for a, b in product(LEVELS, repeat=2):
        case, tape = history(a,b)
        cases.append(case)
        if (a,b) == (F(1,2),F(-3,2)):
            retained = tape
    return {"schema": "oph-native-updates-v1", "pins": codec.pins(),
            "ports": PORTS, "cases": cases, "representative_payloads": ["1/2","-3/2"],
            "representative_tape": retained}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=codec.HERE/"controls.json")
    args = parser.parse_args()
    args.output.write_bytes(codec.canonical(build()))


if __name__ == "__main__":
    main()
