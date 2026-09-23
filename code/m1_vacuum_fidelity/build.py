"""Build a compact candidate; verify.py reconstructs its evidence independently."""

import json
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "m1_vacuum_fidelity"

from . import model
from .verify import HERE, parent_claims, pins


def build():
    return dict(schema="oph-m1-vacuum-fidelity-v1", sources=pins(),
                parent_claims=parent_claims(), evidence=model.candidate())


if __name__ == "__main__":
    (HERE/"receipt.json").write_text(json.dumps(build(), sort_keys=True, separators=(",", ":"))+"\n",
                                    encoding="utf-8", newline="\n")
    print("Built full-population vacuum-fidelity evidence")
