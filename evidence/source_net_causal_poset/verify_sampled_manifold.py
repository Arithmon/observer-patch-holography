"""Verify pinned sampled-manifold sources, retained arithmetic and small controls.

Large q55/q89 sample draws and their covariance were not retained here. This
checks their derived arithmetic and the original producer against brute-force
small orders; it does not reconstruct missing large-sample standard errors.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE / "manifold_sampled"

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", action="store_true")
    parser.add_argument("--replay-references", action="store_true", help="also regenerate the five million seeded continuum reference pairs")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "manifest.json").read_text())
    for path, row in manifest["files"].items():
        blob = (ROOT / path).read_bytes()
        assert len(blob) == row["bytes"] and hashlib.sha256(blob).hexdigest() == row["sha256"], path
    receipt = HERE / "source_net_manifold_sampled_q55_q89_2026-09-25.json"
    assert hashlib.sha256(receipt.read_bytes()).hexdigest() == manifest["receipt_sha256"]
    env = dict(os.environ, PYTHONPATH=str(ROOT), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    commands = [[sys.executable, "-m", "oph_exact.verify_source_net_manifold_sampled_independent", str(receipt)] + (["--reference-pairs", "-1"] if args.replay_references else []),
                [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", str(ROOT / "tests")]]
    rows=[]
    for command in commands:
        completed = subprocess.run(command, env=env, text=True, capture_output=True)
        print(completed.stdout, end="")
        if completed.returncode:
            print(completed.stderr, file=sys.stderr)
            raise SystemExit(completed.returncode)
        rows.append({"command": [str(Path(x).relative_to(HERE)) if x.startswith(str(HERE)) else x for x in command[1:]], "passed": True, "output": completed.stdout.strip()})
    value={"schema":"oph.source-net-manifold-sampled.verification.v1", "verified":True, "manifest_sha256":hashlib.sha256((ROOT/'manifest.json').read_bytes()).hexdigest(), "checks":rows, "continuum_reference_draws_replayed":args.replay_references, "verifier_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "scope":"Retained derived arithmetic, seeded continuum reference draws, original source pins, and small exact-order controls. The q55/q89 event samples, joint chain-estimator covariance and weighted-CDF standard errors are not reconstructed."}
    if args.write_receipt:
        (ROOT / "verification.json").write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    print("SAMPLED_MANIFOLD_PACKAGE_PASS")

if __name__ == '__main__': main()
