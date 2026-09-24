"""Construct a candidate receipt and require independent replay before writing it."""

import json
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "m1_quantum_transport"

from . import model, verify
from .check import verify_evidence


def build():
    evidence = model.candidate()
    verify_evidence(evidence)
    packet = dict(schema="oph-m1-quantum-transport-v1", sources=verify.pins(),
                  parent_claims=verify.parent_claims(), evidence=evidence)
    path = verify.HERE/"receipt.json"
    path.write_text(json.dumps(packet, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n",
                    encoding="utf-8", newline="\n")
    print("Built independently replayed native quantum transport evidence")


if __name__ == "__main__":
    build()
