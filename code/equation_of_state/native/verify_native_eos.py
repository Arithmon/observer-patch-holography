"""Independent exact verifier; imports neither the producer nor the simulator.

Rebuilds the carrier over Q(phi), where phi²=phi+1, and replays every native
pair mean using rational arithmetic. Work-pressure entries are checked as
analytic jets of a declared extension, never as a thermodynamic measurement.
"""

from __future__ import annotations

import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re


def require(condition: bool, message: str):
    if not condition:
        raise ValueError(message)


def fraction(text) -> F:
    require(isinstance(text, str), "rational_not_string")
    value = F(text)
    require(str(value) == text, "noncanonical_rational")
    return value


def digest(data: dict) -> str:
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def golden_edges() -> list[list[int]]:
    """Exact icosahedral squared distances; coordinates are pairs a+b phi."""
    zero, one, minus_one, phi, minus_phi = (0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)
    vertices = [
        (minus_one, phi, zero), (one, phi, zero),
        (minus_one, minus_phi, zero), (one, minus_phi, zero),
        (zero, minus_one, phi), (zero, one, phi),
        (zero, minus_one, minus_phi), (zero, one, minus_phi),
        (phi, zero, minus_one), (phi, zero, one),
        (minus_phi, zero, minus_one), (minus_phi, zero, one),
    ]
    edges = []
    for i in range(12):
        for j in range(i + 1, 12):
            constant = golden = 0
            for vi, vj in zip(vertices[i], vertices[j]):
                a, b = vi[0] - vj[0], vi[1] - vj[1]
                constant += a * a + b * b
                golden += 2 * a * b + b * b
            if (constant, golden) == (4, 0):
                edges.append([i, j])
    require(len(edges) == 30, "internal_exact_carrier_count")
    require([sum(i in edge for edge in edges) for i in range(12)] == [5] * 12,
            "internal_exact_carrier_degree")
    return edges


def check_jets(rows: list, q: F):
    alphas = [F(-1), F(0), F(1, 3), F(1, 2), F(1)]
    require(len(rows) == len(alphas), "work_jet_count")
    for row, alpha in zip(rows, alphas):
        require(set(row) == {"alpha", "candidate_energy", "candidate_density",
                            "volume_derivative", "candidate_pressure", "candidate_w"},
                "work_jet_fields")
        values = {key: fraction(value) for key, value in row.items()}
        require(values["alpha"] == alpha, "work_jet_alpha")
        require(values["candidate_energy"] == q > 0, "work_jet_energy")
        require(3 * values["candidate_density"] == q, "work_jet_density")
        # This is the analytic logarithmic derivative of the declared power law.
        require(3 * values["volume_derivative"] == -alpha * q, "work_jet_derivative")
        require(values["candidate_pressure"] == -values["volume_derivative"], "work_sign")
        require(values["candidate_w"] * values["candidate_density"] ==
                values["candidate_pressure"], "work_jet_ratio")
        require(values["candidate_w"] == alpha, "work_jet_w")


def check_body(report: dict):
    require(set(report) == {"schema", "claim_boundary", "contract", "seams", "work_extension",
                           "energy_offset_control", "linear_direction_readout_control", "cases",
                           "provenance", "payload_sha256"}, "top_level_fields")
    require(report["schema"] == "oph.native_work_nonidentifiability.v1", "schema")
    payload = copy.deepcopy(report)
    received_digest = payload.pop("payload_sha256")
    require(received_digest == digest(payload), "payload_hash")
    expected_boundary = {
        "primitive": "one twelve-port scalar mean-repair carrier",
        "native_thermodynamic_pressure": None, "native_thermodynamic_energy_density": None,
        "native_thermodynamic_w": None, "native_physical_energy_selected": False,
        "quadratic_drop_identified_as_heat": False, "thermodynamic_equilibrium_established": False,
        "physical_momentum_no_go": False, "all_OPH_work_laws_excluded": False,
        "independent_laboratory_measurement": False,
    }
    require(report["claim_boundary"] == expected_boundary, "scope_promotion_or_missing_boundary")
    # Python 0==False must not allow a numeric substitution for a scope boolean.
    for key, expected in expected_boundary.items():
        require(type(report["claim_boundary"][key]) is type(expected), "scope_type")
    require(report["contract"] == {
        "ports": 12, "seam_count": 30, "attempts_per_case": 60,
        "stopping_rule": "two complete lexicographic or reverse sweeps, with all four cases retained",
        "quadratic_definition": "Q(x)=sum_i x_i^2/2",
        "loss_identity": "Q(x)-Q(R_ij x)=(x_i-x_j)^2/4",
        "bookkeeping_record": "B_next=B+Q_before-Q_after; B_initial=0",
        "bookkeeping_interpretation": "an added diagnostic record; no native bath or heat unit",
    }, "contract")
    edges = golden_edges()
    require(report["seams"] == edges, "exact_carrier_seams")
    work = report["work_extension"]
    require(work == {
        "domain": "V>0, Q(x)>0", "reference_volume": "3",
        "family": "H_alpha(x,V)=Q(x)*(V/V0)^(-alpha)",
        "convention": "p=-partial_V H_alpha at fixed x; rho=H_alpha/V",
        "analytic_consequence": "p/rho=alpha for every V>0",
        "identifiability_scope": "all extensions agree in H, all x derivatives, and repair Q losses at V0",
        "alphas": ["-1", "0", "1/3", "1/2", "1"],
        "physical_volume_supplied": True, "work_energy_identification_supplied": True,
        "work_law_supplied": True, "equilibrium_not_claimed": True,
    }, "work_contract")
    for key in ("physical_volume_supplied", "work_energy_identification_supplied",
                "work_law_supplied", "equilibrium_not_claimed"):
        require(work[key] is True, "work_scope_type")
    definitions = [
        ("pulse_ascending", [12] + [0] * 11, True, edges * 2),
        ("pulse_descending", [12] + [0] * 11, True, list(reversed(edges)) * 2),
        ("constant_control", [1] * 12, True, edges * 2),
        ("no_repair_control", [12] + [0] * 11, False, edges * 2),
    ]
    require(len(report["cases"]) == len(definitions), "case_count")
    for case, (name, initial, enabled, schedule) in zip(report["cases"], definitions):
        require(set(case) == {"name", "repair_enabled", "initial_state", "events", "initial_quadratic",
                              "final_quadratic", "cumulative_quadratic_drop", "quadratic_floor_at_fixed_load",
                              "reference_work_jets_initial", "reference_work_jets_final"}, "case_fields")
        require(case["name"] == name and case["repair_enabled"] is enabled, "case_contract")
        require(case["initial_state"] == list(map(str, initial)), "prepared_state")
        state = list(map(F, initial))
        initial_q = sum(v * v for v in state) / 2
        load = sum(state)
        cumulative = F()
        require(len(case["events"]) == 60, "missing_or_extra_events")
        for number, (event, seam) in enumerate(zip(case["events"], schedule), 1):
            i, j = seam
            before = state.copy()
            before_q = sum(v * v for v in before) / 2
            if enabled:
                state[i] = state[j] = (before[i] + before[j]) / 2
            q = sum(v * v for v in state) / 2
            drop = (before[i] - before[j]) ** 2 / 4 if enabled else F()
            require(before_q - q == drop >= 0, "internal_loss_identity")
            cumulative += drop
            expected = {
                "attempt": number, "seam": seam,
                "readback_before": [str(before[i]), str(before[j])],
                "state_after": list(map(str, state)), "total_scalar_load": str(load),
                "quadratic_before": str(before_q), "quadratic_after": str(q),
                "quadratic_drop": str(drop), "cumulative_quadratic_drop": str(cumulative),
                "bookkeeping_total": str(initial_q),
            }
            require(event == expected, f"event_replay:{name}:{number}")
            require(sum(state) == load and q + cumulative == initial_q, "internal_conservation")
        require(fraction(case["initial_quadratic"]) == initial_q, "initial_quadratic")
        require(fraction(case["final_quadratic"]) == q, "final_quadratic")
        require(fraction(case["cumulative_quadratic_drop"]) == cumulative, "cumulative_drop")
        require(fraction(case["quadratic_floor_at_fixed_load"]) == load * load / 24,
                "fixed_load_floor")
        require(q >= load * load / 24, "quadratic_below_cauchy_schwarz_floor")
        if name.startswith("pulse_"):
            require(cumulative > 0, "nontrivial_repair_control")
        else:
            require(cumulative == 0, "null_control")
        check_jets(case["reference_work_jets_initial"], initial_q)
        check_jets(case["reference_work_jets_final"], q)

    offset = report["energy_offset_control"]
    require(offset == {
        "description": "Adding a V-independent constant changes rho and w but no x force or pressure",
        "reference_Q": "72", "reference_volume": "3", "alpha": "1/3",
        "constant_offset": "72", "pressure_before_and_after": "8",
        "rho_before": "24", "rho_after": "48", "w_before": "1/3", "w_after": "1/6",
    }, "offset_control")
    v, q0, constant, alpha = (fraction(offset[k]) for k in
                               ("reference_volume", "reference_Q", "constant_offset", "alpha"))
    p = alpha * q0 / v
    require(fraction(offset["pressure_before_and_after"]) == p, "offset_pressure")
    for suffix, energy in (("before", q0), ("after", q0 + constant)):
        require(fraction(offset[f"rho_{suffix}"]) == energy / v, "offset_density")
        require(fraction(offset[f"w_{suffix}"]) == p * v / energy, "offset_w")

    moment = report["linear_direction_readout_control"]
    require(moment == {
        "meaning": "candidate x-direction moment of unnormalized icosahedral port coordinates",
        "initial_state": ["12"] + ["0"] * 11, "seam": [0, 1],
        "endpoint_x_coordinates": ["-1", "1"],
        "candidate_direction_moment_before": "-12", "candidate_direction_moment_after": "0",
        "conclusion": "this fixed direction-weighted scalar moment is not conserved by all seam means",
        "scope": "does not exclude other fields, stored momenta, protected memory or enriched dynamics",
    }, "linear_direction_control")
    require([0, 1] in edges, "direction_control_seam")
    require(-12 != (-1 * 6 + 1 * 6), "direction_control_not_conserved")
    provenance = report["provenance"]
    require(set(provenance) == {"repository", "head", "source_sha256", "producer_sha256",
                                "native_primitive_imported"}, "provenance_fields")
    require(provenance["repository"] == "https://github.com/muellerberndt/oph-physics-sim",
            "source_repository")
    require(isinstance(provenance["head"], str) and
            re.fullmatch(r"[0-9a-f]{40}", provenance["head"]) is not None, "source_commit_format")
    require(provenance["native_primitive_imported"] is True, "source_import_boundary")
    require(set(provenance["source_sha256"]) == {"oph_fpe/dynamics/canonical_seam_repair.py",
                                               "oph_fpe/core/icosahedral.py"}, "source_paths")
    for value in [*provenance["source_sha256"].values(), provenance["producer_sha256"]]:
        require(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None,
                "source_hash_format")


def verify(report, sim_root: Path | None = None) -> dict:
    try:
        require(isinstance(report, dict), "report_not_object")
        check_body(report)
        if sim_root is not None:
            for relative, expected in report["provenance"]["source_sha256"].items():
                require(hashlib.sha256((sim_root / relative).read_bytes()).hexdigest() == expected,
                        f"source_bytes_changed:{relative}")
        passed, errors = True, []
    except (AssertionError, KeyError, TypeError, ValueError, ZeroDivisionError, OSError,
            OverflowError, RecursionError) as error:
        passed, errors = False, [str(error)]
    return {"schema": "oph.native_work_nonidentifiability.verification.v1",
            "status": "PASS" if passed else "FAIL", "receipt": passed,
            "independent_implementation": True, "producer_imported": False,
            "simulator_imported": False, "exact_topology_reconstructed": passed,
            "exact_events_replayed": 240 if passed else 0,
            "candidate_work_jets_checked": 40 if passed else 0,
            "source_bytes_checked": sim_root is not None and passed,
            "native_thermodynamic_w": None,
            "scope": "finite native scalar replay and supplied work-law ambiguity; no thermodynamic realization",
            "errors": errors}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", nargs="?", type=Path,
                        default=Path(__file__).with_name("native_eos_receipt.json"))
    parser.add_argument("--sim-root", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = verify(json.loads(args.receipt.read_text()), args.sim_root)
    serialized = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.out:
        args.out.write_text(serialized)
    print(serialized, end="")
    raise SystemExit(0 if result["receipt"] else 1)


if __name__ == "__main__":
    main()
