"""Cross-package checks against the audited issue 776 diagnostic.

Shared wiring does not identify the two operation grammars. The q13 comparison
is applied to logical commits only after the native verifier has authenticated
every routed primitive and its consumed writers.
"""
from __future__ import annotations

from contextlib import redirect_stdout
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
from support_wiring import verify as diagnostic


def require(condition, message):
    if not condition:
        raise ValueError(message)


def compare_support(packet, geometry, level):
    require(packet["level"] == level and packet["carriers"] == len(geometry["faces"]),
            "776/777 support census")
    for name, field in (("faces", "faces"), ("intra_carrier_seams", "seams"), ("glued_pairs", "w12")):
        require(np.array_equal(packet[name], geometry[field]), "776/777 support differs: "+name)


def check_shared_support():
    diagnostic.verify_inventory(diagnostic.DEFAULT)
    result = {}
    for level in (3, 4, 5):
        path = (ROOT/"code/source_routing/support_w12_l3.json" if level == 3
                else HERE/f"support_w12_l{level}.json")
        geometry_path = diagnostic.HERE/f"geometry/geometry_l{level}.npz"
        with np.load(geometry_path) as geometry:
            packet = diagnostic.read_json(path)
            compare_support(packet, geometry, level)
        result[str(level)] = {"carriers": packet["carriers"], "glued_seams": len(packet["glued_pairs"]),
                              "routing_capture_sha256": diagnostic.digest(path),
                              "diagnostic_capture_sha256": diagnostic.digest(geometry_path)}
    return result


def compare_projection_arrays(reference, offsets, indices, logical):
    require(np.array_equal(offsets, reference["indptr"])
            and np.array_equal(indices, reference["indices"]), "776/777 q13 read menus differ")
    values = reference["values"]
    require(logical.dtype == np.dtype("<i8") and logical.shape == values.shape
            and values.dtype == np.dtype("<i8"), "776/777 exact logical value format")
    require(np.array_equal(logical, values), "776/777 q13 logical values differ")


def check_q13_projection(info, offsets, indices, logical):
    """Called only after full baseline tape replay and logical refinement."""
    require(info["q"] == 13, "776 comparison is the q13 baseline")
    supports = check_shared_support()
    # Keep the routing CLI's single JSON result on stdout machine-readable.
    with redirect_stdout(sys.stderr), np.load(diagnostic.HERE/"geometry/geometry_l3.npz") as geometry:
        row = diagnostic.verify_q13(diagnostic.DEFAULT, geometry)
        # Scientific controls remain independently executed diagnostic families;
        # this does not assert routed 1D/2D executions.
        diagnostic.verify_q13_controls(diagnostic.DEFAULT, geometry)
    require((info["sites"], info["rounds"], info["centre_site"], info["logical_reads"])
            == (row["sites"], row["rounds"], row["centre"], row["reads"]), "776/777 q13 census")
    with np.load(diagnostic.DEFAULT/"q13_reads.npz") as reference:
        compare_projection_arrays(reference, offsets, indices, logical)
    return {"scope": "q13 baseline logical projection; not the full primitive-event order",
            "support_captures": supports,
            "diagnostic_receipt_sha256": diagnostic.digest(diagnostic.DEFAULT/"support_wiring_receipt.json"),
            "logical_events_compared": int(logical.size), "logical_reads": row["reads"],
            "logical_values_sha256": hashlib.sha256(logical.tobytes()).hexdigest(),
            "q13_interval_interior_flags": [r["interior"] for r in row["intervals"]],
            "same_support_and_logical_read_menu": True,
            "canonical_only_routing_derived": False}


if __name__ == "__main__":
    print(json.dumps({"shared_support": check_shared_support()}, sort_keys=True))
