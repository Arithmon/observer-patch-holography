"""Nested same-location covariance from complete ancestral regression rows.

This is a mathematical projection identity, not a spatial two-point function.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def covariance_from_nested_residuals(variances):
    variances = np.asarray(variances, dtype=float)
    if np.any(~np.isfinite(variances)) or np.any(variances < 0):
        raise ValueError("invalid variance")
    if np.any(np.diff(variances) > 2e-9 * variances.max()):
        raise ValueError("nested residual variances must decrease with coarsening")
    index = np.arange(len(variances))
    covariance = variances[np.maximum(index[:, None], index[None, :])]
    # Do not clamp numerical or physical discrepancies away.
    bands = np.concatenate([-np.diff(variances), variances[-1:]])
    return covariance, bands


def build(directory=HERE):
    directory = Path(directory)
    path = directory / "receipt.json"
    source = json.loads(path.read_text())
    spec = directory / "nested_spec.json"
    chains = []
    for chain in source["chains"]:
        result = {key: chain[key] for key in ["source_path", "parent_level", "mode", "parent_seed", "fresh_seed"]}
        result["group_cells"] = [r["cells_per_group"] for r in chain["by_scale"]]
        for label, field in [("inherited_linear", "inherited_linear_residual_variance"),
                             ("joint_linear", "joint_linear_residual_variance")]:
            covariance, bands = covariance_from_nested_residuals([r[field] for r in chain["by_scale"]])
            result[label] = {"same_location_nested_covariance": covariance.tolist(),
                             "orthogonal_scale_increment_variances_plus_coarsest": bands.tolist()}
        chains.append(result)
    result = {
        "schema": "oph.native-nested-linear-residual-covariance.v1",
        "input_receipt_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "spec_sha256": hashlib.sha256(spec.read_bytes()).hexdigest(),
        "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "interpretation": "Covariance of nested group means lifted to the same uniformly selected finest group; not covariance at two separated sky positions.",
        "identity": "Cov(e_M,e_N)=Var(e_max(M,N)); nested residual projection spaces are subspaces of finer ones.",
        "chains": chains,
    }
    return result


def main(argv=None, directory=HERE):
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true",
                        help="Compare the deterministic transform without writing any files")
    args = parser.parse_args(argv)
    result = build(directory)
    destination = Path(directory) / "nested_receipt.json"
    if args.check:
        if json.loads(destination.read_text()) != result:
            raise SystemExit("nested receipt differs from deterministic recomputation")
        print("nested receipt verified")
    else:
        destination.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
