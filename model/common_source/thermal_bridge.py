"""Replay the fixed-reference modular interpretation of the common dynamics."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import thermal_clock


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path, help="create a receipt, refusing overwrite")
    group.add_argument("--check", type=Path, help="reexecute and compare a retained receipt")
    args = parser.parse_args()
    result = thermal_clock.report()
    if args.write:
        serialized = canonical(result)
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(serialized)
    if args.check:
        compare_receipt(json.loads(args.check.read_text(encoding="utf-8")), result)
    print(json.dumps({"all_pass": result["all_pass"],
                      "receipt": "fixed Gibbs reference separated from actual recorded state",
                      "reference_selected_by_repair": False,
                      "physical_clock_calibration_established": False}, indent=2))


if __name__ == "__main__":
    main()
