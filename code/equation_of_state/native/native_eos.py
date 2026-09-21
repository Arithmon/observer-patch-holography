r"""Exact native repair traces and supplied-work nonidentifiability witnesses.

This diagnostic observes a bounded twelve-port software patch with local scalar
state, seam readback and mean repair, retaining every read/write in a public
evidence bundle. The recorder is external; it is not a native thermal bath.
No thermodynamic pressure, energy, temperature or equation of state is inferred.


# Scoped work-law ambiguity for a fixed scalar repair primitive

Let `x ∈ R^n` with `n ≥ 2`. The native seam update `R_ij` replaces `x_i,x_j`
with their arithmetic mean and leaves all other coordinates fixed. Let
`S(x)=Σ_i x_i` and `Q(x)=½Σ_i x_i²`.

Direct expansion gives

\[
S(R_{ij}x)=S(x),\qquad
Q(x)-Q(R_{ij}x)=\frac{(x_i-x_j)^2}{4}\geq0.
\]

Thus an added record satisfying `B_0=0` and
`B_(t+1)=B_t+Q(x_t)-Q(x_(t+1))` obeys
`Q(x_t)+B_t=Q(x_0)` exactly for every finite schedule. The record is a
bookkeeping extension. Calling it heat would require additional physical
energy, bath and unit identifications. At fixed total load,
`Q(x) ≥ S(x)²/(2n)` by Cauchy–Schwarz. None of these identities mentions volume,
stress, temperature or thermodynamic equilibrium.

For any reference volume `V0>0`, any real `α`, any `V>0` and any `x` with
`Q(x)>0`, define the additional work-energy function

\[
H_\alpha(x,V)=Q(x)(V/V_0)^{-\alpha}.
\]

This function is smooth and positive on its domain. At `V=V0` it agrees with
`Q(x)` for **every** state, not only the recorded trajectory. Every derivative
with respect to scalar state coordinates at that volume is therefore the same
for every `α`. All native repair maps, schedules, load values and quadratic
decreases at `V0` are unchanged. Nevertheless, under the explicitly supplied
work convention and energy-density identification

\[
p_\alpha=-\left.\frac{\partial H_\alpha}{\partial V}\right|_x,
\qquad \rho_\alpha=H_\alpha/V,
\]

ordinary differentiation gives

\[
p_\alpha=\alpha H_\alpha/V,
\qquad p_\alpha/\rho_\alpha=\alpha.
\]

In particular `α=1/3` and `α=1/2` are positive-pressure work extensions that
agree on all native state-only observations at the reference volume but have
different pressure/density ratios. Fixed-volume native traces, even augmented
by exact knowledge of `Q` on all states, do not select such a ratio within this
extension class. This is an identifiability statement for this declared
primitive and these absent work data. It is not a no-go theorem for OPH,
nor a construction of an equilibrium equation of state.

The energy reference also matters. Adding a volume-independent constant `C`
to a supplied `H` changes no scalar gradient or pressure derivative, but changes
`ρ` to `(H+C)/V`. Wherever both energies are positive, the ratio consequently
changes from `pV/H` to `pV/(H+C)`. The executable control uses `Q=72`, `V0=3`,
`α=1/3` and `C=72`: pressure stays `8`, while the ratio changes from `1/3` to
`1/6`. No assumption of thermodynamic extensivity is used or established.

Finally, a fixed linear readout `m(x)=Σ_i c_i x_i` changes under one seam by

\[
m(R_{ij}x)-m(x)=\tfrac12(c_i-c_j)(x_j-x_i).
\]

It is invariant for every state on that seam exactly when `c_i=c_j`.
On a connected allowed-seam graph all such coefficients are constant, so only
total-load multiples survive in this fixed linear scalar class. The twelve-port
certificate uses the actual edge `(0,1)` with unnormalized icosahedral
x-coordinates `(-1,1)` and initial port loads `(12,0)`. Its candidate direction
moment changes from `-12` to `0` after averaging. This distinguishes that scalar
candidate from a momentum-conserving collision gas. It excludes neither stored
momentum coordinates nor tensor fields, non-scalar readouts, protected records,
driving or other dynamics. The existing general native linear-invariant proof
is `Lean/Geometry/SourceRecordProtection.lean` in the research repository.

To select a physical EoS, the missing input must constrain the work/strain
response and the physical energy reference, or supply independently justified
energy and stress observables with a thermodynamic preparation and measurement
protocol. More fixed-volume scalar replay cannot distinguish the extensions
constructed here.

"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import sys


SCHEMA = "oph.native_work_nonidentifiability.v1"
REFERENCE_VOLUME = F(3)
ALPHAS = (F(-1), F(0), F(1, 3), F(1, 2), F(1))


def payload_hash(data: dict) -> str:
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def quadratic(state) -> F:
    return sum((x * x for x in state), F()) / 2


def jet(alpha: F, q: F) -> dict:
    """Analytic derivative of supplied H=Q (V/V0)^(-alpha), at V=V0.

    The returned numbers do not estimate a derivative from native execution.
    They instantiate the declared work convention p=-partial_V H at fixed x.
    """
    density = q / REFERENCE_VOLUME
    derivative = -alpha * density
    return {"alpha": str(alpha), "candidate_energy": str(q),
            "candidate_density": str(density), "volume_derivative": str(derivative),
            "candidate_pressure": str(-derivative), "candidate_w": str(alpha)}


def produce(sim_root: Path) -> dict:
    sys.path.insert(0, str(sim_root))
    from oph_fpe.dynamics.canonical_seam_repair import (
        reference_edges, edge_conditional_expectation, apply_fraction_matrix,
    )

    edges = reference_edges()
    matrices = {edge: edge_conditional_expectation(edge) for edge in edges}
    cases = []
    definitions = (
        ("pulse_ascending", [12] + [0] * 11, True, False),
        ("pulse_descending", [12] + [0] * 11, True, True),
        ("constant_control", [1] * 12, True, False),
        ("no_repair_control", [12] + [0] * 11, False, False),
    )
    for name, initial, enabled, descending in definitions:
        state = tuple(map(F, initial))
        initial_q = quadratic(state)
        initial_load = sum(state)
        cumulative = F()
        events = []
        schedule = tuple(reversed(edges)) if descending else edges
        for attempt, edge in enumerate(schedule * 2, 1):
            before = state
            state = apply_fraction_matrix(matrices[edge], state) if enabled else state
            old_q, q = quadratic(before), quadratic(state)
            loss = old_q - q
            expected = (before[edge[0]] - before[edge[1]]) ** 2 / 4 if enabled else F()
            assert loss == expected >= 0
            assert sum(state) == initial_load
            cumulative += loss
            assert q + cumulative == initial_q
            events.append({
                "attempt": attempt, "seam": list(edge),
                "readback_before": [str(before[i]) for i in edge],
                "state_after": list(map(str, state)),
                "total_scalar_load": str(sum(state)),
                "quadratic_before": str(old_q), "quadratic_after": str(q),
                "quadratic_drop": str(loss), "cumulative_quadratic_drop": str(cumulative),
                "bookkeeping_total": str(q + cumulative),
            })
        cases.append({
            "name": name, "repair_enabled": enabled,
            "initial_state": list(map(str, initial)), "events": events,
            "initial_quadratic": str(initial_q), "final_quadratic": str(quadratic(state)),
            "cumulative_quadratic_drop": str(cumulative),
            "quadratic_floor_at_fixed_load": str(initial_load ** 2 / 24),
            "reference_work_jets_initial": [jet(alpha, initial_q) for alpha in ALPHAS],
            "reference_work_jets_final": [jet(alpha, quadratic(state)) for alpha in ALPHAS],
        })

    source_paths = ("oph_fpe/dynamics/canonical_seam_repair.py", "oph_fpe/core/icosahedral.py")
    result = {
        "schema": SCHEMA,
        "claim_boundary": {
            "primitive": "one twelve-port scalar mean-repair carrier",
            "native_thermodynamic_pressure": None,
            "native_thermodynamic_energy_density": None,
            "native_thermodynamic_w": None,
            "native_physical_energy_selected": False,
            "quadratic_drop_identified_as_heat": False,
            "thermodynamic_equilibrium_established": False,
            "physical_momentum_no_go": False,
            "all_OPH_work_laws_excluded": False,
            "independent_laboratory_measurement": False,
        },
        "contract": {
            "ports": 12, "seam_count": 30, "attempts_per_case": 60,
            "stopping_rule": "two complete lexicographic or reverse sweeps, with all four cases retained",
            "quadratic_definition": "Q(x)=sum_i x_i^2/2",
            "loss_identity": "Q(x)-Q(R_ij x)=(x_i-x_j)^2/4",
            "bookkeeping_record": "B_next=B+Q_before-Q_after; B_initial=0",
            "bookkeeping_interpretation": "an added diagnostic record; no native bath or heat unit",
        },
        "seams": [list(edge) for edge in edges],
        "work_extension": {
            "domain": "V>0, Q(x)>0",
            "reference_volume": str(REFERENCE_VOLUME),
            "family": "H_alpha(x,V)=Q(x)*(V/V0)^(-alpha)",
            "convention": "p=-partial_V H_alpha at fixed x; rho=H_alpha/V",
            "analytic_consequence": "p/rho=alpha for every V>0",
            "identifiability_scope": "all extensions agree in H, all x derivatives, and repair Q losses at V0",
            "alphas": list(map(str, ALPHAS)),
            "physical_volume_supplied": True,
            "work_energy_identification_supplied": True,
            "work_law_supplied": True,
            "equilibrium_not_claimed": True,
        },
        "energy_offset_control": {
            "description": "Adding a V-independent constant changes rho and w but no x force or pressure",
            "reference_Q": "72", "reference_volume": "3", "alpha": "1/3",
            "constant_offset": "72", "pressure_before_and_after": "8",
            "rho_before": "24", "rho_after": "48", "w_before": "1/3", "w_after": "1/6",
        },
        "linear_direction_readout_control": {
            "meaning": "candidate x-direction moment of unnormalized icosahedral port coordinates",
            "initial_state": ["12"] + ["0"] * 11, "seam": [0, 1],
            "endpoint_x_coordinates": ["-1", "1"],
            "candidate_direction_moment_before": "-12", "candidate_direction_moment_after": "0",
            "conclusion": "this fixed direction-weighted scalar moment is not conserved by all seam means",
            "scope": "does not exclude other fields, stored momenta, protected memory or enriched dynamics",
        },
        "cases": cases,
        "provenance": {
            "repository": "https://github.com/muellerberndt/oph-physics-sim",
            "head": subprocess.check_output(["git", "-C", str(sim_root), "rev-parse", "HEAD"], text=True).strip(),
            "source_sha256": {p: hashlib.sha256((sim_root / p).read_bytes()).hexdigest() for p in source_paths},
            "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "native_primitive_imported": True,
        },
    }
    result["payload_sha256"] = payload_hash(result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sim-root", type=Path,
                        default=Path(__file__).resolve().parents[3].parent / "oph-physics-sim")
    parser.add_argument("--out", type=Path, default=Path(__file__).with_name("native_eos_receipt.json"))
    args = parser.parse_args()
    report = produce(args.sim_root.resolve())
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.out), "native_w": None,
                      "supplied_work_w_values": report["work_extension"]["alphas"],
                      "case_count": len(report["cases"]), "attempt_count": 240,
                      "payload_sha256": report["payload_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
