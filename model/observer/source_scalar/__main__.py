"""Run the emitted source scalar apparatus and replay its retained receipt."""
import argparse
import json
from pathlib import Path

from . import experiment
from tools.strict_json import load


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=experiment.TRIALS)
    parser.add_argument("--seed", type=int, default=experiment.SEED)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path, help="create a receipt, refusing overwrite")
    group.add_argument("--check", type=Path, help="replay exact fields and floats within 2e-10")
    parser.add_argument("--emit", type=Path, help="also emit complete d=8,n=128 routines and repeat recipe, refusing overwrite")
    args = parser.parse_args()
    result = experiment.run(args.trials, args.seed)
    if args.write:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(experiment.canonical(result))
    if args.check:
        experiment.compare_receipt(load(args.check), result)
    if args.emit:
        with args.emit.open("x", encoding="utf-8") as stream:
            stream.write(experiment.canonical(experiment.program(8, 128)))
    finest = result["levels"][-1]["step_refinement"][-1]
    print(json.dumps({"all_pass": result["all_pass"], "dimensions": list(experiment.DIMENSIONS),
                      "steps": list(experiment.STEPS), "trials_per_dimension": args.trials,
                      "finest_plus_probability": .5+finest["response"],
                      "finite_split_certificate": finest["finite_split_certificate"],
                      "resources": finest["histories"][-1]["resource_counts"],
                      "scope": result["scope"]}, indent=2))


if __name__ == "__main__":
    main()
