#!/usr/bin/env python3
"""Compare-only distances between the criticality boundary family and data.

This module consumes the Higgs reference and the cross-section top pole-mass
reference from the frozen PDG table, together with the boundary-scale scan.
It reports arithmetic distances to those measurements.  The Higgs output is a
tree coordinate and the top conversion lacks a complete electroweak matching
budget, so distances in experimental-sigma units are not physical pulls.
Nothing here is a prediction ancestor; the scan itself reads no reference.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SCAN = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_criticality_boundary_scan.json"
)
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_criticality_comparison.json"
)


def _distance_row(
    row: dict[str, Any],
    mt_ref: float,
    mt_sigma: float,
    mh_ref: float,
    mh_sigma: float,
    mh_band: float,
    mt_band: float,
) -> dict[str, Any]:
    mt = float(row["mt_pole_gev"])
    mh = float(row["mh_tree_gev"])
    return {
        "mt_pole_gev": mt,
        "mh_tree_gev": mh,
        "mt_relative": mt / mt_ref - 1.0,
        "mh_relative": mh / mh_ref - 1.0,
        "mt_offset_in_experimental_sigma_units": (mt - mt_ref) / mt_sigma,
        "mh_offset_in_experimental_sigma_units": (mh - mh_ref) / mh_sigma,
        "mt_within_matching_band": abs(mt / mt_ref - 1.0) <= mt_band,
        "mh_within_matching_band": abs(mh / mh_ref - 1.0) <= mh_band,
    }


def _relation_test(
    scan: dict[str, Any], mt_ref: float, mh_ref: float, mh_sigma: float
) -> dict[str, Any]:
    """Evaluate m_H on the criticality curve at the measured top mass.

    The curve is fit-free; the measured top selects the point on it, so this
    block tests the m_t to m_H relation independently of the boundary-scale
    selection theorem.  Compare-only by construction.
    """

    result: dict[str, Any] = {}
    for label, curve in (
        ("one_loop", scan.get("one_loop_curve", [])),
        ("two_loop", scan.get("two_loop_curve", [])),
    ):
        rows = [r for r in curve if isinstance(r, dict) and "mt_pole_gev" in r]
        rows.sort(key=lambda r: r["mt_pole_gev"])
        mts = [r["mt_pole_gev"] for r in rows]
        mhs = [r["mh_tree_gev"] for r in rows]
        scales = [r["boundary_scale_gev"] for r in rows]
        if len(rows) < 2 or not (mts[0] <= mt_ref <= mts[-1]):
            result[label] = {
                "status": "measured_top_outside_curve_range",
                "curve_mt_range_gev": [mts[0], mts[-1]] if mts else None,
            }
            continue
        for i in range(len(mts) - 1):
            if mts[i] <= mt_ref <= mts[i + 1]:
                frac = (mt_ref - mts[i]) / (mts[i + 1] - mts[i])
                mh_at = mhs[i] + frac * (mhs[i + 1] - mhs[i])
                scale_at = math.exp(
                    math.log(scales[i])
                    + frac * (math.log(scales[i + 1]) - math.log(scales[i]))
                )
                result[label] = {
                    "status": "interpolated",
                    "mh_on_curve_at_measured_mt_gev": mh_at,
                    "mh_relative": mh_at / mh_ref - 1.0,
                    "mh_offset_in_experimental_sigma_units": (
                        mh_at - mh_ref
                    ) / mh_sigma,
                    "implied_boundary_scale_gev": scale_at,
                }
                break
    return result


def build_artifact(scan: dict[str, Any], references: dict[str, Any]) -> dict[str, Any]:
    entries = references["entries"]
    higgs = entries["higgs"]
    top = entries["top_quark"]
    if top["source"]["summary_id"] != "Q007TP4":
        raise ValueError("top comparison requires the cross-section pole-mass entry")
    mh_ref = float(higgs["value_gev"])
    mh_sigma = float(higgs["error_plus_gev"])
    mt_ref = float(top["value_gev"])
    mt_sigma = float(top["error_plus_gev"])
    mh_band = float(scan["matching_bands"]["mh_tree_to_pole_relative_band"])
    mt_band = float(scan["matching_bands"]["mt_electroweak_matching_relative_band"])

    named = {}
    for name, row in scan["one_loop_named_boundaries"].items():
        named[name] = _distance_row(
            row, mt_ref, mt_sigma, mh_ref, mh_sigma, mh_band, mt_band
        )

    two_loop = {}
    if scan["two_loop_named_boundaries"].get("available"):
        for name, row in scan["two_loop_named_boundaries"].items():
            if isinstance(row, dict) and "mt_pole_gev" in row:
                two_loop[name] = _distance_row(
                    row, mt_ref, mt_sigma, mh_ref, mh_sigma, mh_band, mt_band
                )

    curve = scan["one_loop_curve"]
    best = min(
        curve,
        key=lambda r: (r["mt_pole_gev"] / mt_ref - 1.0) ** 2
        + (r["mh_tree_gev"] / mh_ref - 1.0) ** 2,
    )
    best_row = _distance_row(best, mt_ref, mt_sigma, mh_ref, mh_sigma, mh_band, mt_band)
    best_row["boundary_scale_gev"] = best["boundary_scale_gev"]

    relation_test = _relation_test(scan, mt_ref, mh_ref, mh_sigma)

    archived = named["mu_U_gauge_unification"]
    single_anchor = named["E_cell_pixel_energy"]

    return {
        "artifact": "oph_d11_criticality_comparison",
        "schema_version": 1,
        "status": "compare_only_distance_surface",
        "row_class": "compare_only_never_a_prediction_ancestor",
        "promotion_allowed": False,
        "physical_fit_evaluable": False,
        "comparison_convention": (
            "PDG_2026_top_cross_section_pole_mass_vs_model_QCD_converted_top; "
            "PDG_2026_Higgs_mass_vs_model_tree_Higgs_coordinate"
        ),
        "distance_status": (
            "experimental_sigma_offsets_are_arithmetic_only; model_matching_and_"
            "truncation_bands_are_declared_tolerances_not_statistical_uncertainties"
        ),
        "references": {
            "mh_gev": mh_ref,
            "mh_sigma": mh_sigma,
            "mh_source_id": higgs["source"]["summary_id"],
            "mt_gev": mt_ref,
            "mt_sigma": mt_sigma,
            "mt_source_id": top["source"]["summary_id"],
            "note": "consumed here only; the scan artifact reads no reference",
        },
        "one_loop_named_distances": named,
        "two_loop_named_distances": two_loop,
        "one_loop_curve_nearest_point": best_row,
        "mt_mh_relation_test": relation_test,
        "reading": {
            "archived_mu_U_boundary": {
                "mt_relative": archived["mt_relative"],
                "mh_relative": archived["mh_relative"],
            },
            "single_anchor_E_cell_boundary": {
                "mt_relative": single_anchor["mt_relative"],
                "mh_relative": single_anchor["mh_relative"],
            },
            "statement": (
                "The one-loop curve does not reach the measured top reference. "
                "The two-loop relation uses that top datum as an input and "
                "returns a tree Higgs coordinate under supplied SM running "
                "and matching assumptions. Declared truncation/matching bands "
                "are not statistical uncertainties. The adopted boundary-scale "
                "formula was selected after exposure to the implied scale; "
                "source selection, physical pole matching and theory uncertainty "
                "remain open."
            ),
        },
        "checks": {
            "scan_checks_pass": bool(scan["checks_pass"]),
            "distances_finite": all(
                math.isfinite(v)
                for row in named.values()
                for v in row.values()
                if isinstance(v, float)
            ),
        },
        "checks_pass": bool(scan["checks_pass"]),
    }


def build() -> dict[str, Any]:
    scan = json.loads(SCAN.read_text(encoding="utf-8"))
    references = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))
    return build_artifact(scan, references)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "one_loop_named_distances": {
                    name: {
                        "mt_rel": round(row["mt_relative"], 5),
                        "mh_rel": round(row["mh_relative"], 5),
                    }
                    for name, row in artifact["one_loop_named_distances"].items()
                },
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
