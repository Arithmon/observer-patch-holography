"""Run and replay the conditional scalar observer experiment."""
import argparse
import json
from pathlib import Path

from . import experiment


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=experiment.TRIALS)
    parser.add_argument("--seed", type=int, default=experiment.SEED)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path, help="create a receipt, refusing overwrite")
    group.add_argument("--check", type=Path, help="replay exact fields and floating diagnostics within 2e-10")
    args = parser.parse_args()
    result = experiment.run(args.trials, args.seed)
    if args.write:
        with args.write.open("x") as stream:
            stream.write(experiment.canonical(result))
    if args.check:
        experiment.compare_receipt(json.loads(args.check.read_text()), result)
    print(json.dumps({"all_pass": result["all_pass"], "grids": list(experiment.GRIDS),
                      "trials_per_grid": args.trials,
                      "finest_current": result["executions"][-1]["sample"]["statistics"]["current"],
                      "finest_mathematical_bound": result["refinement_bounds"][-1]}, indent=2))


if __name__ == "__main__":
    main()
