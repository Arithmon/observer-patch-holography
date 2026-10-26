"""Replay the complete retained mass and charged-field experiments for CI."""
from __future__ import annotations

import argparse
import importlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from .hierarchy_diagnosis.audit import (
    REPLAY_ABSOLUTE_TOLERANCE,
    REPLAY_RELATIVE_TOLERANCE,
    compare_replay,
)


PACKAGES = ("mass_response", "particle_identification", "chiral_candidate")


def compare_report(retained, current, *, host_timing=False):
    """Compare every scientific field; validate host duration/version metadata."""
    current = json.loads(json.dumps(current, allow_nan=False))
    retained = json.loads(json.dumps(retained, allow_nan=False))
    if host_timing:
        for document in (retained, current):
            environment = document.pop("environment")
            if (not isinstance(environment, dict)
                    or set(environment) != {"python", "numpy", "scipy"}
                    or any(not isinstance(value, str) or not value.strip()
                           for value in environment.values())):
                raise ValueError("invalid execution environment metadata")
        if retained["cases"].keys() != current["cases"].keys():
            raise ValueError("complete experiment case set changed")
        for name in retained["cases"]:
            for document in (retained, current):
                duration = document["cases"][name].pop("elapsed_host_seconds")
                if (type(duration) not in (int, float)
                        or not math.isfinite(duration) or duration < 0):
                    raise ValueError("invalid host duration")
    compare_replay(retained, current)


def compare_arrays(retained, current):
    """Check quantum outputs, probabilities and complete record labels."""
    if set(retained) != set(current):
        raise ValueError("retained quantum array keys changed")
    for key, value in current.items():
        before, after = np.asarray(retained[key]), np.asarray(value)
        if before.shape != after.shape or before.dtype.kind != after.dtype.kind:
            raise ValueError(f"retained quantum array shape/type changed: {key}")
        if after.dtype.kind in "fc":
            equal = (np.isfinite(before).all() and np.isfinite(after).all()
                     and np.allclose(before, after, rtol=REPLAY_RELATIVE_TOLERANCE,
                                     atol=REPLAY_ABSOLUTE_TOLERANCE))
        else:
            equal = np.array_equal(before, after)
        if not equal:
            raise ValueError(f"retained quantum array differs from execution: {key}")


def replay(package, path):
    if package not in PACKAGES:
        raise ValueError("unknown experiment")
    module = importlib.import_module(f"model.common_source.{package}.run")
    path = Path(path)
    module.verify(path)
    retained = strict_json.load(path)
    recomputed = module.run()
    if package == "mass_response":
        compare_report(retained, recomputed, host_timing=True)
    else:
        current, arrays = recomputed
        # verify() checks this adjacent artifact's checksum and custody before
        # the independent execution; compare its contents, not ZIP timestamps.
        current["artifact"] = retained["artifact"]
        compare_report(retained, current)
        with np.load(path.with_name(retained["artifact"]["filename"]),
                     allow_pickle=False) as stored:
            compare_arrays(dict(stored), arrays)
    print(f"Recomputed complete {package} history and checked every retained output", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", choices=PACKAGES)
    parser.add_argument("--check", type=Path, required=True)
    args = parser.parse_args()
    replay(args.package, args.check)


if __name__ == "__main__":
    main()
