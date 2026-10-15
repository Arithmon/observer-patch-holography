"""Replay both empirical packets while sharing parameter-keyed source caches.

Each existing CLI performs its complete prediction, scoring and receipt check.
Only the same process's canonical word and executed-operation caches are reused;
no prediction, trajectory, fitted preparation or previous acceptance is reused.
The fixed module/receipt allowlist does not accept caller-selected code or files.
"""
from __future__ import annotations

import argparse
import importlib
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
EMPIRICAL_REPLAYS = (
    ("model.anharmonic_transfer", None, "model/anharmonic_transfer/results.json"),
    ("model.anharmonic_observations", None, "model/anharmonic_observations/results.json"),
)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args(argv)
    for module, _, receipt in EMPIRICAL_REPLAYS:
        print(f"Replaying {receipt}", flush=True)
        canonical_main = importlib.import_module(f"{module}.__main__").main
        previous_argv = sys.argv
        try:
            sys.argv = [module, "--check", str(ROOT / receipt)]
            canonical_main()
        finally:
            sys.argv = previous_argv


if __name__ == "__main__":
    main()
