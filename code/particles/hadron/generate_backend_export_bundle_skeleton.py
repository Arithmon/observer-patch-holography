#!/usr/bin/env python3
"""Emit a fillable backend-export skeleton for hadron production.

This uses the documented backend-export contract. Instead
of asking an external RHMC/HMC code to hand-build a giant JSON file, it emits a
manifest that records the exact dataset paths and per-source coordinates. When
h5py is available it can also emit a placeholder HDF5 file with the canonical
dataset tree.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.hadron.backend_export_bundle import build_backend_export_skeleton


def _load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _array_file_for_dataset(dataset_path: str) -> str:
    safe = dataset_path.strip("/").replace("/", "__")
    return f"arrays/{safe}.txt"


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate a raw backend export bundle skeleton.")
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--payload", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--profile-id", default="oph_reference_backend_v1")
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="Write only backend_run_manifest.json and dataset_index.json; do not require h5py or create correlators.h5.",
    )
    parser.add_argument(
        "--dataset-index-output",
        default=None,
        help="Optional explicit path for the dataset index JSON. Defaults to <out-dir>/dataset_index.json.",
    )
    args = parser.parse_args()

    receipt = _load_json(args.receipt)
    payload = _load_json(args.payload)
    manifest, datasets = build_backend_export_skeleton(receipt, payload, profile_id=args.profile_id)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "backend_run_manifest.json"
    dataset_index_path = Path(args.dataset_index_output) if args.dataset_index_output else out_dir / "dataset_index.json"
    h5_path = out_dir / "correlators.h5"

    if not args.manifest_only:
        try:
            import h5py  # type: ignore
            import numpy as np  # type: ignore
        except Exception as exc:  # pragma: no cover - dependency error path
            raise RuntimeError(
                "h5py and numpy are required to create correlators.h5; rerun with --manifest-only "
                "to emit the production handoff without an HDF5 placeholder"
            ) from exc

        with h5py.File(h5_path, "w") as h5:
            for dset_path, length in datasets:
                dset = h5.create_dataset(dset_path, shape=(length,), dtype="<f8")
                dset[...] = np.full((length,), math.nan, dtype=np.float64)
                dset.attrs["status"] = "fill_with_real_backend_output"

    dataset_index = {
        "artifact": "oph_hadron_backend_dataset_index",
        "format_version": 1,
        "manifest": str(manifest_path.name),
        "correlators_hdf5": manifest["files"]["correlators_hdf5"],
        "placeholder_hdf5_written": not args.manifest_only,
        "datasets": [
            {
                "path": dset_path,
                "length": length,
                "dtype": "float64_le",
                "array_file": _array_file_for_dataset(dset_path),
                "fill_status": "requires_real_backend_output",
            }
            for dset_path, length in datasets
        ],
        "notes": [
            "Every listed dataset must be filled by production backend execution before ingestion.",
            "NaN placeholders or synthetic arrays must not be normalized as production hadron data.",
        ],
    }

    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    dataset_index_path.parent.mkdir(parents=True, exist_ok=True)
    dataset_index_path.write_text(json.dumps(dataset_index, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {manifest_path}")
    print(f"wrote {dataset_index_path}")
    if not args.manifest_only:
        print(f"wrote {h5_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
