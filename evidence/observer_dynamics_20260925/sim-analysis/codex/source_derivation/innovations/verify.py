"""Fast read-only custody checks, without fresh RNG reconstruction."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(directory=HERE, sim=None):
    directory = Path(directory)
    sim = directory.parents[2] if sim is None else Path(sim)
    receipt = json.loads((directory / "receipt.json").read_text())
    nested = json.loads((directory / "nested_receipt.json").read_text())
    spec = json.loads((directory / "spec.json").read_text())
    pins = [
        (receipt["producer_sha256"], directory / "extract.py"),
        (receipt["spec_sha256"], directory / "spec.json"),
        (nested["producer_sha256"], directory / "nested.py"),
        (nested["spec_sha256"], directory / "nested_spec.json"),
        (nested["input_receipt_sha256"], directory / "receipt.json"),
    ]
    pins.extend((entry["sha256"], sim / entry["path"]) for entry in spec["inputs"])
    for expected, path in pins:
        if sha(path) != expected:
            raise ValueError(f"digest mismatch: {path}")
    expected_inputs = sorted(entry["path"] for entry in spec["inputs"]
                             if entry["role"] == "refinement_chain")
    actual_inputs = sorted(path.relative_to(sim).as_posix()
                           for path in (sim / "data/refine_ensemble").glob("chain_*.json"))
    if actual_inputs != expected_inputs:
        raise ValueError("archived chain inventory mismatch")
    for data in [receipt, nested]:
        if sorted(c["source_path"] for c in data["chains"]) != expected_inputs:
            raise ValueError("artifact chain inventory mismatch")
    if len(receipt["chains"]) != spec["chain_count"]:
        raise ValueError("chain count mismatch")
    for chain in receipt["chains"]:
        if [r["cells_per_group"] for r in chain["by_scale"]] != spec["group_cells"]:
            raise ValueError("artifact scale inventory mismatch")
    for chain in nested["chains"]:
        if chain["group_cells"] != spec["group_cells"]:
            raise ValueError("nested scale inventory mismatch")
    return {"hashes_checked": len(pins), "chains": len(receipt["chains"]),
            "scale_rows": sum(len(c["by_scale"]) for c in receipt["chains"])}


if __name__ == "__main__":
    print(json.dumps(check(), sort_keys=True))
