"""Execute and verify the source-compiled C/E discrimination experiment."""
import argparse
import json
from pathlib import Path

from tools.strict_json import load
from . import experiment


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=experiment.SEED)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path, help="create a receipt, refusing overwrite")
    group.add_argument("--check", type=Path, help="recompute all fields; 2e-10 floating tolerance")
    args = parser.parse_args()
    expected = load(args.check) if args.check else None
    result = experiment.run(args.seed)
    if args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(experiment.canonical(result))
    if args.check:
        experiment.compare_receipt(expected, result)
    print(json.dumps({"all_pass": result["all_pass"],
                      "protocol_sha256": result["protocol_sha256"],
                      "plus_probabilities": {
                          name: {row["setting"]: row["plus_probability"] for row in candidate["executions"]}
                          for name, candidate in result["candidates"].items()},
                      "empirical_confirmation": False,
                      "scope": result["scope"]}, indent=2))


if __name__ == "__main__":
    main()
