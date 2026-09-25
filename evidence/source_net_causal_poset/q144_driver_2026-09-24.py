"""Record-metric causal-poset lane at q = 144 (Fibonacci index 12) with sampled pair counting.

Writes a standalone receipt; the committed family receipt (q <= 89, exact pairs) is untouched.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

from oph_exact import source_net as S


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=62)
    parser.add_argument("--exact-support-limit", type=int, default=500_000)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    t0 = time.time()
    receipt = S.build(levels=(12,), processes=args.processes, exact_support_limit=args.exact_support_limit)
    receipt["extension_note"] = ("sampled extension of the committed exact family receipt: level 12 (q = 144) with stratified "
                                 "sampled pair counting on both intervals; the committed receipt's claim boundary is unchanged")
    receipt["wall_seconds"] = round(time.time() - t0, 1)
    data = S.canonical(receipt)
    args.out.write_bytes(data)
    print("Q144_RECEIPT", len(data), hashlib.sha256(data).hexdigest(), f"{time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
