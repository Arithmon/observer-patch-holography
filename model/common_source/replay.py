"""Replay the complete retained mass and charged-field experiments for CI."""
from __future__ import annotations

import argparse
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from tools import strict_json
from .hierarchy_diagnosis.audit import (
    REPLAY_ABSOLUTE_TOLERANCE,
    REPLAY_RELATIVE_TOLERANCE,
    compare_replay,
)


PACKAGES = ("mass_response", "particle_identification", "chiral_candidate")


def _validate_preparation_zero_diagnostic(case):
    """Validate the exact-zero angle tally without treating it as a cost.

    LAPACK roundoff can change whether a computed preparation angle is
    exactly zero. The compiler executes and charges every request, including
    zeros; native pulse, service and integrated-norm costs remain compared.
    This exception applies only to the mass preparation's named diagnostic.
    """
    if "resources_per_realized_nonabort_history" not in case:
        return
    try:
        word = case["resources_per_realized_nonabort_history"]["preparation_source_word"]
        zero_count = word.pop("zero_angle_requests")
        pulse_count = word["native_pulses"]
    except (KeyError, TypeError, AttributeError) as error:
        raise ValueError("missing preparation zero-angle diagnostic") from error
    if (type(zero_count) is not int or type(pulse_count) is not int
            or not 0 <= zero_count <= pulse_count):
        raise ValueError("invalid preparation zero-angle diagnostic")


def compare_report(retained, current, *, host_timing=False):
    """Compare science and costs; validate host metadata and exact-zero tally."""
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
                _validate_preparation_zero_diagnostic(document["cases"][name])
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
    if package == "mass_response":
        # This producer records every loaded repository module. Execute its
        # original CLI in a fresh interpreter so verifier imports cannot alter
        # its computational source closure. Scratch outputs live outside the
        # repository and are removed even if a comparison fails.
        with tempfile.TemporaryDirectory(prefix="oph-mass-replay-") as scratch:
            output = Path(scratch) / "results.json"
            subprocess.run([sys.executable, "-m", "model.common_source.mass_response.run",
                            "--output", str(output)], cwd=module.ROOT, check=True)
            compare_report(retained, strict_json.load(output), host_timing=True)
    else:
        current, arrays = module.run()
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
