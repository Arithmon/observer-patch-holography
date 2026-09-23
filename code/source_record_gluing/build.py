"""Capture compact commitments to every unfiltered reference-process event."""
from collections import Counter
import hashlib
import json
from pathlib import Path
from .process import execute


def capture():
    cases = []
    for q in (2, 3, 5):
        probes = [(-1, None)]+[(stage, i) for stage in (-1, 0) for i in range(q**3)]
        for stage, intervention in probes:
            digest = hashlib.sha256()
            counts = Counter()
            checkpoints = []
            for row in execute(q, 2, intervention, stage):
                digest.update((json.dumps(row, separators=(",", ":"))+"\n").encode())
                counts[row[0]] += 1
                if row[0] == "checkpoint":
                    checkpoints.append(row[2])
            cases.append({"q": q, "layers": 2, "intervention": intervention,
                          "intervention_stage": stage,
                          "stream_sha256": digest.hexdigest(),
                          "operation_counts": dict(sorted(counts.items())),
                          "checkpoint_sha256": hashlib.sha256(json.dumps(
                              checkpoints, separators=(",", ":")).encode()).hexdigest()})
    return {"schema": "oph.record-gluing-reference.v2",
            "scope": "generated process theory of proposed RG; no native implementation asserted",
            "cases": cases}


if __name__ == "__main__":
    path = Path(__file__).with_name("capture.json")
    path.write_text(json.dumps(capture(), separators=(",", ":"))+"\n", encoding="utf-8", newline="\n")
    print(path.name, path.stat().st_size)
