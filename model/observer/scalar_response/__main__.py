"""Run and replay the conditional scalar observer experiment."""
import argparse
import json
from pathlib import Path

from . import experiment
from tools.strict_json import load


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--continue-field", action="store_true",
                        help="two reads of one retained field and three native repairs")
    parser.add_argument("--trials", type=int, default=None)
    parser.add_argument("--seed", type=int, default=experiment.SEED)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path, help="create a receipt, refusing overwrite")
    group.add_argument("--check", type=Path, help="replay exact fields and floating diagnostics within 2e-10")
    args = parser.parse_args()
    runner = experiment
    if args.continue_field:
        from . import continuation
        runner = continuation
    trials = runner.TRIALS if args.trials is None else args.trials
    result = runner.run(trials, args.seed)
    if args.write:
        with args.write.open("x") as stream:
            stream.write(experiment.canonical(result))
    if args.check:
        experiment.compare_receipt(load(args.check), result)
    if args.continue_field:
        print(json.dumps({"all_pass": result["all_pass"], "grids": list(runner.GRIDS),
                          "trials_per_grid": trials,
                          "protocol": "two retained field reads, three native repairs"}, indent=2))
        return
    print(json.dumps({"all_pass": result["all_pass"], "grids": list(experiment.GRIDS),
                      "trials_per_grid": trials,
                      "finest_current": result["executions"][-1]["sample"]["statistics"]["current"],
                      "finest_mathematical_bound": result["refinement_bounds"][-1]}, indent=2))


if __name__ == "__main__":
    main()
