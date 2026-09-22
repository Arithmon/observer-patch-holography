"""Execute every two-carrier central-record intervention and deduplicate tapes."""
from fractions import Fraction as Q
from pathlib import Path
import json


def central_execution(labels, steps):
    deviation = [Q(1, 4) for _ in labels]
    logs = [[] for _ in labels]
    for tick in range(steps):
        for carrier, label in enumerate(labels):
            query = tick % 12
            before = deviation[carrier]
            if label == query:
                deviation[carrier] = Q(0)
            logs[carrier].append((tick, label, query, str(before),
                                  str(deviation[carrier]), label, 1))
    return logs


def capture():
    tapes, indices, cases = [], {}, []
    for source in range(12):
        for receiver in range(12):
            executed = central_execution([source, receiver], 24)
            refs = []
            for tape in executed:
                key = tuple(tape)
                if key not in indices:
                    indices[key] = len(tapes)
                    tapes.append(tape)
                refs.append(indices[key])
            cases.append({"labels": [source, receiver], "tapes": refs})
    return {"schema": "oph.scalar-seam-records.v1", "tapes": tapes, "cases": cases}


if __name__ == "__main__":
    target = Path(__file__).with_name("records.json")
    target.write_text(json.dumps(capture(), separators=(",", ":"))+"\n",
                      encoding="utf-8", newline="\n")
    print(target.name, target.stat().st_size)
