#!/usr/bin/env python3
"""Retrospective fixed-tilt comparison on official ACT DR6 posterior chains.

The retained input is a lossless two-column projection of the four ordered
chains, including burn-in. No likelihood is fitted or new chain generated.
Posterior positions are not p-values, model probabilities or FZ-15 verdicts.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import mpmath

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
from cosmology.fixed_capacity_wlaw.official_desi_dr2_chain_audit import weighted_quantiles

PROJECTION = HERE / "act_dr6_ns_weight.json.gz"
RECEIPT = HERE / "results.json"
SOURCE_ARCHIVE_SHA256 = "8f2801a69663ef1142bb05660e0118ef7b1e0bdf3d05afba434a0569023b98ab"
PROJECTION_SHA256 = "13dab540b652300bd000cf1ac3b197af55c7c411de215682482ea55993f40485"
DECIMAL_TOKEN = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
P = F(16309682094039593, 10**16)
EDGE_TILT = 1 - P / 48
QUANTILES = (F(1, 1000), F(25, 1000), F(16, 100), F(1, 2),
             F(84, 100), F(975, 1000), F(999, 1000))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def predictions():
    # The historical alternative is explicitly a diagnostic, not a replacement
    # of the frozen edge-center branch. Its numerical precision is disclosed.
    mp = mpmath.mp.clone()
    mp.dps = 60
    clock = 1 - mp.e * (mp.mpf(P.numerator) / P.denominator - (1 + mp.sqrt(5)) / 2)
    return {"edge_center": EDGE_TILT,
            "historical_clock_diagnostic": F(str(clock))}


def parse_rows(rows):
    """Preserve decimal boundaries before comparisons; no float coercion."""
    parsed = []
    for row in rows:
        if (not isinstance(row, list) or len(row) != 2
                or any(not isinstance(token, str) or not DECIMAL_TOKEN.fullmatch(token)
                       for token in row)):
            raise ValueError("each projected row needs two original decimal tokens")
        try:
            weight, ns = map(F, row)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError("finite weight and ns tokens required") from exc
        if weight <= 0 or weight.denominator != 1:
            raise ValueError("this ordered Cobaya projection requires positive integer multiplicities")
        if not 0 < ns < 2:
            raise ValueError("ns is outside the declared scalar-index domain")
        parsed.append((weight, ns))
    if not parsed:
        raise ValueError("empty chain")
    return parsed


def discard_burn_in(rows, fraction=F(1, 2)):
    """Drop the official fraction of stored rows, separately in each chain."""
    if not isinstance(fraction, F) or not 0 <= fraction < 1:
        raise ValueError("burn-in must be an exact fraction in [0,1)")
    drop = len(rows) * fraction.numerator // fraction.denominator
    kept = rows[drop:]
    if not kept:
        raise ValueError("no samples remain after burn-in")
    return kept, drop


def ordering_diagnostics(rows, targets):
    """Count actual changes and tail entries; these are not ESS estimates."""
    changes = sum(left[1] != right[1] for left, right in zip(rows, rows[1:]))
    entries = {}
    for name, target in targets.items():
        previous, count = False, 0
        for _, value in rows:
            current = value <= target
            count += current and not previous
            previous = current
        entries[name] = count
    return {"ns_changes_between_stored_rows": changes,
            "contiguous_tail_excursions": entries,
            "effective_sample_size_estimator": False}


def summary(rows, targets=None):
    if not rows or any(w <= 0 for w, _ in rows):
        raise ValueError("nonempty positive-weight posterior required")
    targets = predictions() if targets is None else targets
    # ACT oversamples nuisance updates: many consecutive rows have identical
    # ns. Aggregate only for marginal arithmetic, retaining order in the input.
    atoms = defaultdict(F)
    for w, x in rows:
        atoms[x] += w
    weight = sum(atoms.values())
    center = sum(w * x for x, w in atoms.items()) / weight
    variance = sum(w * (x - center)**2 for x, w in atoms.items()) / weight
    positions = {}
    for name, target in targets.items():
        below = sum(w for x, w in atoms.items() if x < target)
        equal = atoms.get(target, F(0))
        above = sum(w for x, w in atoms.items() if x > target)
        assert below + equal + above == weight
        positions[name] = {
            "target": str(target), "target_decimal": float(target),
            "weight_strictly_below": str(below), "weight_equal": str(equal),
            "weight_strictly_above": str(above),
            "posterior_mass_le_target_exact": str((below + equal) / weight),
            "posterior_mass_le_target": float((below + equal) / weight),
            "posterior_mass_gt_target": float(above / weight),
        }
    return {
        "stored_rows": len(rows), "distinct_ns_values": len(atoms),
        "total_multiplicity": str(weight),
        "ns_mean": float(center), "ns_posterior_sd": math.sqrt(float(variance)),
        "step_cdf_quantiles": weighted_quantiles(list(atoms.items()), QUANTILES),
        "positions": positions,
    }


def analyze(payload):
    if payload.get("schema") != "oph.act_dr6.ns_projection.v1":
        raise ValueError("unexpected projection schema")
    if payload.get("archive_sha256") != SOURCE_ARCHIVE_SHA256:
        raise ValueError("different source archive; review before comparison")
    if str(payload.get("burn_in_fraction")) != "0.5":
        raise ValueError("the official baseline burn-in is 0.5")
    chains = payload.get("chains", [])
    if len(chains) != 4 or len({c["member"] for c in chains}) != 4:
        raise ValueError("four distinct original ordered chains required")
    targets = predictions()
    selected, per_chain = [], []
    for chain in chains:
        if len(chain.get("source_sha256", "")) != 64:
            raise ValueError("original member digest missing")
        rows, drop = discard_burn_in(parse_rows(chain["rows"]))
        selected.append(rows)
        per_chain.append({"member": chain["member"],
                          "source_sha256": chain["source_sha256"],
                          "discarded_stored_rows": drop,
                          "ordering_diagnostics": ordering_diagnostics(rows, targets),
                          "summary": summary(rows, targets)})
    pooled = summary([r for rows in selected for r in rows], targets)
    # No normal error model is inferred from this finite-chain sensitivity.
    leave_one_out = []
    for omitted, chain in enumerate(chains):
        parts = [c["summary"] for index, c in enumerate(per_chain) if index != omitted]
        total = sum(F(part["total_multiplicity"]) for part in parts)
        leave_one_out.append({
            "omitted_member": chain["member"],
            "posterior_mass_le_target": {
                key: float(sum(F(part["positions"][key]["weight_strictly_below"])
                               + F(part["positions"][key]["weight_equal"])
                               for part in parts) / total)
                for key in targets},
        })
    return {
        "schema": "oph.act_dr6.tilt_retrospective.v1",
        "classification": "retrospective conditional-branch posterior-position audit",
        "source": {key: value for key, value in payload.items() if key != "chains"},
        "prediction": {
            "P_exact": str(P), "P_role": "inherited CODATA-alpha-calibrated coordinate",
            "edge_center_ns_exact": str(EDGE_TILT),
            "physical_primordial_embedding": "additional conditional hypothesis",
            "native_observer_model_derivation_established": False,
            "historical_clock_alternative": "exposed diagnostic; 60-digit evaluation; no branch substitution",
        },
        "analysis": {
            "burn_in": "discard floor(0.5 * stored_rows) separately in each ordered chain",
            "weights": "original positive integer Cobaya multiplicities",
            "quantiles": "inverse empirical weighted CDF; no interpolation",
            "tail_arithmetic": "exact rational comparison of original decimal tokens",
            "mc_uncertainty": "per-chain and leave-one-out sensitivity only; no tail-ESS or confidence certificate",
            "posterior_is_likelihood": False, "tail_mass_is_p_value": False,
            "posterior_probability_OPH_true": None,
            "new_independent_data": False, "FZ15_prospective_verdict": None,
            "joint_calibration_and_observation_uncertainty_propagated": False,
        },
        "pooled": pooled, "per_chain": per_chain, "leave_one_chain_out": leave_one_out,
        "total_tail_excursions_across_chains": {
            key: sum(c["ordering_diagnostics"]["contiguous_tail_excursions"][key]
                     for c in per_chain) for key in targets},
    }


def build(projection=PROJECTION):
    if digest(projection) != PROJECTION_SHA256:
        raise ValueError("retained projection hash differs; reconstruct from the pinned archive")
    with gzip.open(projection, "rt", encoding="utf-8") as stream:
        result = analyze(json.load(stream))
    result["custody"] = {
        "projection_sha256": digest(projection),
        "analysis_source_sha256": digest(__file__),
        "quantile_source_sha256": digest(ROOT / "code/cosmology/fixed_capacity_wlaw/official_desi_dr2_chain_audit.py"),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projection", type=Path, default=PROJECTION)
    parser.add_argument("--output", type=Path, default=RECEIPT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build(args.projection)
    if args.check:
        if result != json.loads(args.output.read_text()):
            raise SystemExit("retained tilt posterior receipt differs from replay")
        print("ACT DR6 retrospective tilt receipt matches exact weighted replay")
    else:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps(result["pooled"], indent=2))


if __name__ == "__main__":
    main()
