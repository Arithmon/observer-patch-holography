"""Replay the separate finite-resolution instrument and its error certificate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import calibration, sanity, smooth_protocol, smooth_refinement, units


def run(steps=512):
    execution = smooth_protocol.run(steps)
    readouts = calibration.report(steps=steps)
    unit_check = units.report()
    sanity_checks = sanity.report(steps=steps)
    return {
        "schema": "oph.common_source.smooth_bridge.v1",
        "execution": execution,
        "refinement": smooth_refinement.report(steps=steps),
        "internal_readout_checks": readouts,
        "conditional_unit_calibration": unit_check,
        "physical_relation_sanity_checks": sanity_checks,
        "all_pass": all(row["all_pass"] for row in
                        (execution, readouts, unit_check, sanity_checks)),
        "scope": {
            "detector_meaning": "Ramsey effects (I +/- sin(Q0))/2",
            "distinct_from_original_sharp_sign_readout": True,
            "old_modular_clock_transferred_to_this_protocol": False,
            "source_executed_occupation_dimensions": [2, 4],
            "physical_units_calibrated": False,
            "spatial_continuum_or_empirical_agreement_established": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path, help="create a receipt, refusing overwrite")
    group.add_argument("--check", type=Path, help="reexecute and compare a retained receipt")
    args = parser.parse_args()
    result = run()
    if args.write:
        serialized = canonical(result)
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(serialized)
    if args.check:
        compare_receipt(json.loads(args.check.read_text(encoding="utf-8")), result)
    print(json.dumps({"all_pass": result["all_pass"],
                      "refinement": result["refinement"]["status"],
                      "scope": result["scope"]}, indent=2))


if __name__ == "__main__":
    main()
