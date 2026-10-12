"""Verify artifact bindings and replay the runnable model receipts."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PACKETS = ("model/observer", "model/observer/p_response", "model/local_tick",
           "model/observer/scalar_response", "model/observer/source_scalar",
           "model/predictions", "model/common_source")
REPLAYS = (
    ("model.observer", "native", "model/observer/results/native.json"),
    ("model.observer", "capacity", "model/observer/results/capacity.json"),
    ("model.observer", "run", "model/observer/results/run.json"),
    ("model.observer.p_response", None, "model/observer/p_response/results.json"),
    ("model.local_tick", None, "model/local_tick/results.json"),
    ("model.observer.scalar_response", None, "model/observer/scalar_response/results.json"),
    ("model.observer.scalar_response.local_evolution", None,
     "model/observer/scalar_response/local_results.json"),
    ("model.observer.source_scalar", None, "model/observer/source_scalar/results.json"),
    ("model.predictions", None, "model/predictions/results.json"),
    ("model.common_source", None, "model/common_source/results.json"),
    ("model.common_source.bridge", None, "model/common_source/bridge_results.json"),
    ("model.common_source.smooth_bridge", None, "model/common_source/smooth_results.json"),
    ("model.common_source.thermal_bridge", None, "model/common_source/thermal_results.json"),
    ("model.common_source.causality", None, "model/common_source/causality_results.json"),
)


def verify_artifacts():
    covered = set()
    for name in PACKETS:
        folder = ROOT / name
        manifest = json.loads((folder / "manifest.json").read_text())
        files = manifest["files"]
        if not files:
            raise ValueError(f"Empty artifact manifest: {name}")
        for relative, expected in files.items():
            target = (folder / relative).resolve()
            if not target.is_relative_to(folder.resolve()):
                raise ValueError(f"Artifact outside its package: {name}/{relative}")
            actual = hashlib.sha256(target.read_bytes()).hexdigest()
            if actual != expected:
                raise ValueError(f"Artifact digest mismatch: {name}/{relative}")
            covered.add(target)
    for name in PACKETS:
        for target in (ROOT / name).rglob("*"):
            if (target.is_file() and target.suffix in {".py", ".md", ".json"}
                    and target.name != "manifest.json"
                    and target.resolve() not in covered):
                raise ValueError(f"Unbound artifact: {target.relative_to(ROOT)}")
    return len(covered)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hashes-only", action="store_true",
                        help="check artifact integrity without recalculation")
    args = parser.parse_args()
    count = verify_artifacts()
    print(f"Verified {count} model artifacts", flush=True)
    if args.hashes_only:
        return
    environment = dict(os.environ)
    environment.setdefault("OPENBLAS_NUM_THREADS", "1")
    environment.setdefault("OMP_NUM_THREADS", "1")
    for module, command, receipt in REPLAYS:
        argv = [sys.executable, "-m", module]
        if command is not None:
            argv.append(command)
        argv.extend(["--check", receipt])
        print(f"Replaying {receipt}", flush=True)
        subprocess.run(argv, cwd=ROOT, env=environment, check=True)
    print(f"All {len(REPLAYS)} model receipts reproduced", flush=True)


if __name__ == "__main__":
    main()
