"""Collar-matrix realization probe for the THERMO-REALIZATION receipt.

The receipt requires a source-derived collar transition matrix that
either equals the conditional-resampling kernel or passes
stochasticity, stationarity, and protected-charge preservation, with
detailed balance where microscopic reversibility is claimed. This
probe consumes the empirical repair transition matrix that the
simulator's own transition-clock builder counts from stored
observer histories of the earned 64k universe run, on the declared
visible-packet quotient. The matrix is counted from the run's
repair dynamics; no entry is constructed from the target kernel
formula. The probe measures the receipt properties on that matrix and
records the result with content pins to the source run; it does not
close the receipt.

Measured here, with the protected datum declared as the record
family: row stochasticity, off-fibre leakage of the protected datum,
the equal-fibre-row comparison in both the pairwise-row and the
stationary-profile form, the detailed-balance error of the
reversibilized chain, and relative-entropy descent to the stationary
law of the reversibilized chain. The raw chain's reducibility, which
the simulator's own eligibility gate names as a blocker, is recomputed
from exact count support. The probe also exhausts all fifteen coordinate projections
obtained by retaining a nonempty subset of the four committed packet fields.
It isolates the eight-state repair-load count aggregation as ergodic but
nonreversible and checks every closed communicating class of the fine chain.
These projected count kernels are not automatically Markov quotients of the
fine chain; the probe checks exact count-kernel strong lumpability for each declared
coordinate map. It does not enumerate arbitrary partitions or nonlinear,
statistical, stochastic, weakly lumpable, or history-dependent quotient maps.
This is a bounded audit of the pinned 20-state table, without an additional
observer-patch simulation.

Run with --write to refresh the committed probe receipt.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

import exact_count_chain as exact

from conditional_repair_certificate import (
    ThermoError,
    canonical_json_bytes,
    require,
    tagged_sha256,
)

HERE = Path(__file__).resolve().parent
ARTIFACT_DIR = HERE / "runtime" / "realization_probe_64k"
MATRIX_PATH = ARTIFACT_DIR / "finite_repair_transition_matrix.npz"
REPORT_PATH = ARTIFACT_DIR / "finite_repair_transition_matrix_report.json"
PROBE_PATH = HERE / "runtime" / "collar_matrix_realization_probe.json"

# Content pins established at production time. The matrix and report
# were produced by the simulator CLI from the earned 64k run snapshot;
# the source pins identify the exact inputs consumed.
PINS: dict[str, str] = {
    "matrix_npz_sha256": (
        "sha256:314114f430a0d9a93b6a573ec61bebec213193d5"
        "ee7ff9a8bbbd7d57390b16bf"
    ),
    "report_json_sha256": (
        "sha256:01f880b8b1d8bc77ac497efb8c6bdb1c88c489df"
        "15b7528179f9bb04ac988c54"
    ),
    "source_observer_views_sha256": (
        "sha256:3568f030412ac3db259f438ff6c36e5cee51214d"
        "5004984569814905a6abec44"
    ),
    "source_manifest_sha256": (
        "sha256:11ba89c21f612da1579f4aa24242a4e0933b85b8"
        "c4c91c5a421065e82f312d24"
    ),
    "producer_checkout_commit": "4aa2ce703b5cd172f687c7c9bc3d2d4aa04ed11e",
    "source_run": (
        "oph-workspace://oph-physics-sim/data/earned_runs/"
        "oph_universe_64k_3p1d_reearned"
    ),
    "producer_command": (
        "oph-fpe finite-repair-transition-clock --run-dir "
        "data/earned_runs/oph_universe_64k_3p1d_reearned --out <out> "
        "--packet-fields "
        "record_family,checkpoint_class,s3_sector_class,"
        "repair_load_bucket"
    ),
}

FIBRE_FIELD = "record_family"
DB_TOL = 1e-12
KL_STEPS = 16
LUMPABILITY_TOL = 1e-12


def file_sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def load_inputs() -> tuple[dict[str, Any], dict[str, Any]]:
    """Load the pinned artifacts.

    Preserve stored scalar values. NumPy decodes the container and replays
    its original binary64 export; exact rational algebra then owns the
    normalized count kernels. Logarithmic readouts use a private context.
    """
    require(MATRIX_PATH.exists(), "matrix artifact missing")
    require(REPORT_PATH.exists(), "report artifact missing")
    require(
        file_sha256(MATRIX_PATH) == PINS["matrix_npz_sha256"],
        "matrix artifact drifted from its pin",
    )
    require(
        file_sha256(REPORT_PATH) == PINS["report_json_sha256"],
        "report artifact drifted from its pin",
    )
    with np.load(MATRIX_PATH, allow_pickle=True) as z:
        require(all(z[name].dtype == np.dtype('float64')
                    for name in ('counts', 'raw_empirical', 'reversible_empirical')),
                "pinned export must retain its binary64 source representation")
        data = {
            "counts": z["counts"].tolist(),
            "raw": z["raw_empirical"].tolist(),
            "rev": z["reversible_empirical"].tolist(),
            "labels": [json.loads(str(s)) for s in z["state_labels"]],
        }
    report = json.loads(REPORT_PATH.read_text())
    fields = report['packet_fields']
    require(fields and len(set(fields)) == len(fields), "duplicate source packet fields")
    require(canonical_json_bytes(data['labels']) == canonical_json_bytes(
        [json.loads(label) for label in report['state_labels']]),
        "matrix and report label order differ")
    require(type(report['state_count']) is int
            and len(data['labels']) == report['state_count'], "source state count mismatch")
    require(all([name for name, _ in label] == fields for label in data['labels']),
            "source labels must contain each declared field exactly once in order")
    require(len({canonical_json_bytes(label) for label in data['labels']}) == len(data['labels']),
            "source labels must identify distinct states")
    return data, report


def fibre_of(label: list[list[Any]]) -> Any:
    return label_value(label, FIBRE_FIELD)


def push_vec(mu: list[float], matrix: list[list[float]]) -> list[float]:
    n = len(mu)
    return [
        sum(mu[x] * matrix[x][y] for x in range(n)) for y in range(n)
    ]


def label_value(label: list[list[Any]], field: str) -> Any:
    """Read one named value from the simulator's canonical packet label."""
    values = [value for name, value in label if name == field]
    require(len(values) == 1, f"label must contain exactly one quotient field {field}")
    return values[0]


def coarsen_counts(
    counts: list[list[float]],
    labels: list[list[list[Any]]],
    fields: tuple[str, ...],
) -> tuple[list[list[list[Any]]], list[list[float]], list[list[float]]]:
    """Aggregate the pinned count table through a coordinate map.

    This is a deterministic aggregation of the source-counted table, not a
    fitted transition law. Row normalization produces an empirical stochastic
    kernel but does not by itself prove that the coordinate map is lumpable.
    Keys are sorted by canonical JSON so the resulting matrix has a stable row
    order.
    """
    counts = exact.weight_matrix(counts)
    require(len(labels) == len(counts), "count/label state count mismatch")
    blocks = coordinate_partition_blocks(labels, fields)
    source_index = [0]*len(labels)
    for i, block in enumerate(blocks):
        for j in block:
            source_index[j] = i
    n = len(blocks)
    coarse_counts = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i, row in enumerate(counts):
        for j, value in enumerate(row):
            coarse_counts[source_index[i]][source_index[j]] += value
    matrix = exact.kernel_from_weights(coarse_counts)
    coarse_labels = [
        [[field, label_value(labels[block[0]], field)] for field in fields]
        for block in blocks
    ]
    return coarse_labels, coarse_counts, matrix


# Keep the public probe entry points; one exact implementation owns graph logic.
support_reachability = exact.support_reachability
closed_communicating_classes = exact.closed_classes
stationary_of = exact.stationary_distribution


def pairwise_row_tv_max(matrix: list[list[float]]) -> float:
    """Largest total-variation distance between two rows."""
    out = Fraction(0)
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            out = max(
                out,
                Fraction(1, 2)
                * sum(abs(x - y) for x, y in zip(matrix[i], matrix[j])),
            )
    return out


def coordinate_partition_blocks(
    fine_labels: list[list[list[Any]]], fields: tuple[str, ...]
) -> list[list[int]]:
    """Fine-state blocks induced by a coordinate map, canonically ordered."""
    require(bool(fields) and len(set(fields)) == len(fields),
            "a coarsening needs distinct nonempty fields")
    encoded_blocks: dict[str, list[int]] = {}
    for i, label in enumerate(fine_labels):
        key = tuple(label_value(label, field) for field in fields)
        encoded = json.dumps(key, separators=(",", ":"), sort_keys=True, allow_nan=False)
        encoded_blocks.setdefault(encoded, []).append(i)
    return [encoded_blocks[key] for key in sorted(encoded_blocks)]


def strong_lumpability_max_err(
    fine_matrix: list[list[float]],
    fine_labels: list[list[list[Any]]],
    fields: tuple[str, ...],
) -> float:
    """Maximum strong-lumpability defect for one coordinate map.

    Strong lumpability requires two fine states in the same source block to
    have identical total transition probability into every target block. A
    row-normalized aggregation of counted transitions is a valid stochastic
    kernel even when this condition fails, but in that case it is not a
    certified Markov quotient of the fine chain.
    """
    return exact.lumpability_defect(
        fine_matrix, coordinate_partition_blocks(fine_labels, fields))


def audit_irreducible_chain(matrix, mp) -> dict[str, Any]:
    """Exact finite-chain classifications, with separate numerical KL samples."""
    matrix = exact.stochastic_matrix(matrix)
    n = len(matrix)
    irreducible = all(all(row) for row in support_reachability(matrix))
    period = exact.irreducible_period(matrix) if irreducible else None
    out = {
        "state_count": n,
        "row_sum_max_err": 0.0,
        "irreducible": irreducible,
        "period": period,
        "aperiodic": period == 1 if irreducible else None,
        "aperiodic_self_loop_witness": next(
            (i for i in range(n) if matrix[i][i] > 0), None),
    }
    if not irreducible:
        out.update(dict.fromkeys((
            "stationary_residual_max_err", "stationary_distribution",
            "stationary_distribution_exact", "stationary_min",
            "detailed_balance_max_err", "detailed_balance_defect_exact",
            "reversible", "kl_to_stationary_initial", "kl_to_stationary_final",
            "kl_min_stepwise_descent")))
        return out
    pi = stationary_of(matrix)
    db = exact.detailed_balance_defect(matrix, pi)
    out.update({
        "stationary_residual_max_err": 0.0,
        "stationary_distribution": [float(value) for value in pi],
        "stationary_distribution_exact": [str(value) for value in pi],
        "stationary_min": float(min(pi)),
        "detailed_balance_max_err": float(db),
        "detailed_balance_defect_exact": str(db),
        "reversible": db == 0,
        **kl_diagnostics(matrix, pi, mp),
    })
    return out


def kl(p, q, mp):
    """KL on exact normalized laws; positive Bregman terms avoid cancellation.

    Sum q*((1+u)*log(1+u)-u), u=(p-q)/q. For small u use
    its alternating Taylor series, retaining exact supplied differences before
    evaluating logarithms. This is numerical evaluation, not an interval proof.
    """
    p, q = [exact.rational(x) for x in p], [exact.rational(x) for x in q]
    require(len(p) == len(q) and sum(p) == sum(q) == 1
            and min(p) >= 0 and min(q) >= 0, "KL requires normalized probability laws")
    terms = []
    for a, b in zip(p, q):
        require(b > 0 or a == 0, "stationary law misses support of the iterate")
        if not b:
            continue
        weight = mp.mpf(b.numerator)/b.denominator
        delta = (a-b)/b
        u = mp.mpf(delta.numerator)/delta.denominator
        if not a:
            phi = mp.mpf(1)
        elif abs(delta) <= Fraction(1, 2):
            term = phi = u*u/2
            degree = 2
            while term:
                term *= -u*(degree-1)/(degree+1)
                phi += term
                degree += 1
                if abs(term) <= mp.eps*abs(phi):
                    break
        else:
            ratio = a/b
            r = mp.mpf(ratio.numerator)/ratio.denominator
            phi = r*mp.log(r)-r+1
        terms.append(weight*phi)
    return mp.fsum(terms)


def kl_diagnostics(matrix, pi, mp):
    # Neither equilibrium nor qualitative classifications depend on these
    # sampled, finite-precision logarithms. Preserve the caller's context.
    matrix = exact.stochastic_matrix(matrix)
    pi = [exact.rational(value) for value in pi]
    require(len(pi) == len(matrix) and min(pi) > 0 and sum(pi) == 1
            and push_vec(pi, matrix) == pi,
            "entropy diagnostics require a faithful stationary reference")
    ctx = mp.clone()
    ctx.dps = max(90, mp.dps)
    mu = [Fraction(1, len(matrix))]*len(matrix)
    values = [kl(mu, pi, ctx)]
    descents = []
    for _ in range(KL_STEPS):
        following = push_vec(mu, matrix)
        # Chain rule, without subtracting two almost equal entropies:
        # A_ij=mu_i K_ij, B_ij=(mu K)_j pi_i K_ij/pi_j are probability
        # laws. Stationarity implies D(A||B)=D(mu||pi)-D(mu K||pi).
        # This holds without detailed balance, including genuine zero loss.
        forward, reconstructed = [], []
        for i, row in enumerate(matrix):
            for j, transition in enumerate(row):
                forward.append(mu[i]*transition)
                reconstructed.append(following[j]*pi[i]*transition/pi[j])
        descents.append(kl(forward, reconstructed, ctx))
        mu = following
        values.append(kl(mu, pi, ctx))
    return {
        "kl_to_stationary_initial": ctx.nstr(values[0], 12),
        "kl_to_stationary_final": ctx.nstr(values[-1], 12),
        "kl_min_stepwise_descent": ctx.nstr(min(descents), 12),
    }


def source_kernels(data):
    """Bind both stored arrays to the archived weighted-count constructions.

    Replay the binary64 export exactly before deriving rational kernels from
    the supplied binary64 masses. The latter are the explicitly normalized
    count models, not a claim that rounded exported rows sum exactly to one.
    """
    counts = exact.weight_matrix(data["counts"])
    array = np.asarray(data["counts"], dtype=float)
    require(all(exact.rational(float(array[i, j])) == counts[i][j]
                for i in range(len(counts)) for j in range(len(counts))),
            "source weight conversion lost supplied precision")
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        sym = (array+array.T)/2
        for name, weights in (("raw", array), ("rev", sym)):
            totals = weights.sum(axis=1)
            require(np.isfinite(weights).all() and np.isfinite(totals).all(),
                    "source export normalization is not finite")
            expected = np.zeros_like(weights)
            occupied = totals > 0
            expected[occupied] = weights[occupied]/totals[occupied, None]
            for i in np.flatnonzero(~occupied):
                expected[i, i] = 1
            stored = exact.weight_matrix(data[name])
            require(len(stored) == len(counts) and all(
                exact.rational(float(expected[i, j])) == stored[i][j]
                for i in range(len(counts)) for j in range(len(counts))),
                name+" matrix is not the recorded weighted-count export")
    symmetric = [[counts[i][j]+counts[j][i] for j in range(len(counts))]
                 for i in range(len(counts))]
    raw = exact.kernel_from_weights(counts, empty_rows='absorbing')
    rev = exact.kernel_from_weights(symmetric, empty_rows='absorbing')
    for name, kernel in (('raw', raw), ('rev', rev)):
        require(all((kernel[i][j] > 0) == (data[name][i][j] > 0)
                    for i in range(len(counts)) for j in range(len(counts))),
                name+" source export lost count-kernel support")
    return counts, raw, rev


def build_probe() -> dict[str, Any]:
    import mpmath

    data, report = load_inputs()
    counts, raw, rev = source_kernels(data)
    labels = data["labels"]
    n = len(labels)
    require(
        len(raw) == n and all(len(row) == n for row in raw),
        "matrix shape drift",
    )
    fibres = [canonical_json_bytes(fibre_of(lab)) for lab in labels]
    row_mass = [sum(counts[i]) for i in range(n)]
    total_mass = sum(row_mass)
    require(total_mass > 0, "empty transition count table")
    visit = [m / total_mass for m in row_mass]

    # stochasticity of counted rows
    counted = [i for i in range(n) if row_mass[i] > 0]
    row_sum_err = max(
        abs(sum(raw[i]) - 1) for i in counted
    )
    require(row_sum_err == 0, "normalized count rows are not stochastic")

    # protected-datum leakage: off-fibre mass of counted rows
    off_masses = []
    for i in counted:
        off = sum(
            raw[i][j] for j in range(n) if fibres[j] != fibres[i]
        )
        off_masses.append((off, visit[i]))
    off_max = max(m for m, _ in off_masses)
    off_visit_weighted = sum(m * w for m, w in off_masses) / sum(
        w for _, w in off_masses
    )

    # equal-fibre-row comparison, pairwise form: within a fibre, the
    # conditional-resampling kernel has identical rows after
    # restriction to the fibre
    pair_tv_max = Fraction(0)
    pair_tv_count = 0
    profile_tv_max = Fraction(0)
    for fib in sorted(set(fibres)):
        members = [
            i for i in counted if fibres[i] == fib
        ]
        if len(members) < 2:
            continue
        restricted = []
        for i in members:
            block = [
                raw[i][j] for j in range(n) if fibres[j] == fib
            ]
            s = sum(block)
            if s > 0:
                restricted.append([v / s for v in block])
        for a in range(len(restricted)):
            for b in range(a + 1, len(restricted)):
                tv = Fraction(1, 2) * sum(
                    abs(x - y)
                    for x, y in zip(restricted[a], restricted[b])
                )
                pair_tv_max = max(pair_tv_max, tv)
                pair_tv_count += 1
        # profile form: restricted rows against the visit-weight
        # profile of the fibre
        fibre_visit = [
            visit[i] for i in range(n) if fibres[i] == fib
        ]
        fv_total = sum(fibre_visit)
        if fv_total > 0 and restricted:
            profile = [v / fv_total for v in fibre_visit]
            for row in restricted:
                tv = Fraction(1, 2) * sum(
                    abs(x - y) for x, y in zip(row, profile)
                )
                profile_tv_max = max(profile_tv_max, tv)

    # reversibilized chain: detailed balance and relative-entropy
    # descent to its stationary law
    rev_entry = report["matrices"]["reversible_empirical"]
    rev_summary = rev_entry.get("summary", rev_entry)
    require(
        bool(rev_summary["irreducible"]) and bool(rev_summary["aperiodic"]),
        "reversibilized chain lost ergodicity",
    )
    mp = mpmath.mp.clone()
    mp.dps = 90
    rev_audit = audit_irreducible_chain(rev, mp)
    require(rev_audit["irreducible"] and rev_audit["aperiodic"],
            "count-derived reversibilized chain lost ergodicity")
    db_err = rev_audit["detailed_balance_max_err"]
    require(rev_audit["reversible"], "count-derived reversibilized detailed balance fails")
    kl_initial = rev_audit["kl_to_stationary_initial"]
    kl_final = rev_audit["kl_to_stationary_final"]
    kl_min_descent = rev_audit["kl_min_stepwise_descent"]

    inherited_blockers = list(report.get("blockers", []))

    # Exhaust every nonempty coordinate projection of the four committed
    # packet fields. This is only 15 maps of a 20-state table; it is a bounded
    # audit of the pinned artifact, not an enumeration of arbitrary quotient
    # maps, a new simulation, or a parameter search.
    available_fields = tuple(str(field) for field in report["packet_fields"])
    coarsening_rows: list[dict[str, Any]] = []
    # Fifteen syntactic maps induce only four distinct partitions. Calculate
    # each kernel once, while preserving every declared map in the receipt.
    partition_audits = {}
    for size in range(1, len(available_fields) + 1):
        for fields in itertools.combinations(available_fields, size):
            blocks = coordinate_partition_blocks(labels, fields)
            signature = tuple(tuple(block) for block in blocks)
            if signature not in partition_audits:
                _, _, coarse_matrix = coarsen_counts(counts, labels, fields)
                chain = audit_irreducible_chain(coarse_matrix, mp)
                defect = exact.lumpability_defect(raw, blocks)
                chain.update({
                    "equal_row_pairwise_tv_max": float(pairwise_row_tv_max(coarse_matrix)),
                    "equal_row_pairwise_tv_max_exact": str(pairwise_row_tv_max(coarse_matrix)),
                    "fine_chain_strong_lumpability_max_err": float(defect),
                    "fine_chain_strong_lumpability_defect_exact": str(defect),
                    "fine_chain_strongly_lumpable": defect == 0,
                    # Retained display diagnostic; decisions use the exact field.
                    "fine_chain_strongly_lumpable_at_tolerance": defect <= Fraction(str(LUMPABILITY_TOL)),
                })
                partition_audits[signature] = chain
            coarse_labels, _, _ = coarsen_counts(counts, labels, fields)
            coarsening_rows.append({**partition_audits[signature],
                "packet_fields": list(fields), "labels": coarse_labels,
                "induced_fine_partition_blocks": blocks})

    selected_fields = ["repair_load_bucket"]
    selected = next(
        row
        for row in coarsening_rows
        if row["packet_fields"] == selected_fields
    )
    require(selected["state_count"] == 8, "repair quotient state count drift")
    require(
        bool(selected["irreducible"]) and bool(selected["aperiodic"]),
        "repair quotient lost its raw-chain ergodicity",
    )
    require(
        selected["reversible"] is False,
        "repair quotient unexpectedly became reversible",
    )
    nontrivial_irreducible = [
        row
        for row in coarsening_rows
        if row["state_count"] > 1 and row["irreducible"]
    ]
    nontrivial_reversible = [
        row
        for row in nontrivial_irreducible
        if row["reversible"]
    ]
    require(
        not nontrivial_reversible,
        "a nontrivial reversible raw coarsening appeared",
    )
    distinct_partition_signatures = {
        tuple(tuple(block) for block in row["induced_fine_partition_blocks"])
        for row in coarsening_rows
    }
    require(
        len(distinct_partition_signatures) == 4,
        "constant-field projection deduplication drifted",
    )

    # The alternative recurrent-class route named by issue #688 is also
    # decided on the pinned fine quotient.  Its only closed class is the
    # absorbing freezeout state, so restriction would be mathematically
    # valid but physically degenerate.
    recurrent_classes = closed_communicating_classes(raw)
    require(
        len(recurrent_classes) == 1 and len(recurrent_classes[0]) == 1,
        "fine-quotient recurrent-class structure drifted",
    )
    recurrent_labels = [
        [labels[index] for index in cls] for cls in recurrent_classes
    ]
    protected_cardinality = len(set(fibres))

    body: dict[str, Any] = {
        "schema": "oph.collar_matrix_realization_probe.v4",
        "arithmetic": {
            "source_weights": "exact represented values of the pinned binary64 accumulated weighted counts",
            "source_export_binding": "both stored arrays exactly replay row-normalized C and (C+C.T)/2 in binary64",
            "kernel": "exact rational row normalization of those supplied masses; zero rows explicitly absorbing",
            "stationarity": "unique exact balance solve; no mixing-time or residual stopping rule",
            "classifications": "exact positive support, cycle period, detailed balance and coordinate lumpability",
            "display": "binary64 summaries may round; exact rational fields determine classifications",
            "entropy": "90-digit numerical logarithms on exact finite iterates; contraction loss from forward/reconstructed joint KL, without entropy subtraction; not an interval or continuum certificate",
            "pre_accumulation_precision_recovered": False,
        },
        "status": (
            "MEASURED_PROBE__SOURCE_PRODUCED_MATRIX_ATTAINED__"
            "DECLARED_15_SYNTACTIC_COORDINATE_PROJECTIONS_AUDITED__"
            "FOUR_DISTINCT_INDUCED_PARTITIONS__"
            "OTHER_QUOTIENT_MAPS_OPEN__"
            "EIGHT_STATE_COUNT_AGGREGATED_NONREVERSIBLE_H_THEOREM_PROBE__"
            "FINE_20_STATE_CHAIN_ONLY_CLOSED_CLASS_IS_SINGLETON_FREEZEOUT__"
            "REVERSIBILIZED_KL_DESCENT_VERIFIED__"
            "RAW_CHAIN_REDUCIBLE__RECEIPT_OPEN"
        ),
        "receipt_target": "THERMO-REALIZATION",
        "noncircularity": (
            "the transition matrix is counted by the simulator's "
            "transition-clock builder from stored observer repair "
            "histories of the earned run; no entry is constructed "
            "from the conditional-resampling formula, and the "
            "comparison profiles are built from the chain's own "
            "visit weights"
        ),
        "pins": PINS,
        "quotient": {
            "packet_fields": report["packet_fields"],
            "fibre_field": FIBRE_FIELD,
            "state_count": n,
            "counted_state_count": len(counted),
            "observer_count": report["observer_count"],
            "transition_count": report["transition_count"],
            "weight_field": report["weight_field"],
        },
        "measurements": {
            "row_sum_max_err": float(row_sum_err),
            "off_fibre_mass_max": float(off_max),
            "off_fibre_mass_max_exact": str(off_max),
            "off_fibre_mass_visit_weighted": float(off_visit_weighted),
            "equal_fibre_row_pairwise_tv_max": float(pair_tv_max),
            "equal_fibre_row_pairwise_tv_max_exact": str(pair_tv_max),
            "equal_fibre_row_pairwise_pairs": pair_tv_count,
            "fibre_profile_tv_max": float(profile_tv_max),
            "fibre_profile_tv_max_exact": str(profile_tv_max),
            "protected_datum_cardinality": protected_cardinality,
            "reversibilized_detailed_balance_max_err": db_err,
            "reversibilized_spectral_gap": (
                1.0 - float(rev_summary["lambda_2"])
            ),
            "reversibilized_spectral_gap_source": "archived report; not recomputed or used by the exact equilibrium audit",
            "kl_to_stationary_initial": kl_initial,
            "kl_to_stationary_final": kl_final,
            "kl_steps": KL_STEPS,
            "kl_min_stepwise_descent": kl_min_descent,
            "raw_chain_irreducible": all(all(row) for row in support_reachability(raw)),
            "raw_chain_aperiodic": None,
            "raw_chain_aperiodic_scope": "global period omitted for reducible chain; closed-class periods below",
            "reversibilized_stationary_distribution_exact": rev_audit["stationary_distribution_exact"],
            "reversibilized_reversible": rev_audit["reversible"],
        },
        "raw_coarsening_audit": {
            "audit_bound": (
                "all 15 syntactic coordinate maps that retain a nonempty "
                "subset of the four fields in the pinned 20-state table, "
                "inducing four distinct partitions; 1,028 "
                "observers and 31,744 counted transitions; no additional "
                "observer-patch simulation"
            ),
            "enumeration_scope": {
                "map_class": (
                    "coordinate maps x -> tuple(field(x) for field in S), "
                    "for every nonempty subset S of the four committed fields"
                ),
                "aggregation_rule": (
                    "sum the pinned transition counts over each source and "
                    "target coordinate block, then row-normalize"
                ),
                "exhaustive_within_declared_coordinate_grammar": True,
                "syntactic_maps_may_induce_duplicate_partitions": True,
                "duplicate_reason": (
                    "record_family and s3_sector_class are constant on the "
                    "pinned table, so the 15 syntactic maps induce only four "
                    "distinct fine-state partitions"
                ),
                "arbitrary_set_partitions_enumerated": False,
                "nonlinear_maps_enumerated": False,
                "statistical_or_learned_maps_enumerated": False,
                "stochastic_maps_enumerated": False,
                "history_dependent_maps_enumerated": False,
                "weak_lumpability_tested": False,
                "arbitrary_strongly_lumpable_partition_search_performed": False,
                "coordinate_strong_lumpability_diagnostic_recorded": True,
                "exact_lumpability_proof_emitted": True,
                "interpretation": (
                    "the row-normalized count aggregations are stochastic "
                    "kernels of the supplied weighted masses. Exact rational "
                    "block sums decide strong lumpability of their normalized "
                    "fine kernel; this does not recover an unobserved process "
                    "or precision lost before the counts were stored"
                ),
            },
            "strong_lumpability_tolerance": LUMPABILITY_TOL,
            "syntactic_coordinate_map_count": len(coarsening_rows),
            "distinct_induced_partition_count": len(
                distinct_partition_signatures
            ),
            "nontrivial_irreducible_syntactic_map_count": len(
                nontrivial_irreducible
            ),
            "nontrivial_irreducible_reversible_syntactic_map_count": len(
                nontrivial_reversible
            ),
            "constant_fields": ["record_family", "s3_sector_class"],
            "selected_raw_equilibrium_probe": {
                **selected,
                "microscopic_reversibility_claimed": False,
                "protected_record_family_cardinality": protected_cardinality,
                "protected_charge_test_nontrivial": (
                    protected_cardinality > 1
                ),
                "common_reference_with_state_optimizer_identified": False,
                "interpretation": (
                    "the source-counted repair-load aggregation is an "
                    "eight-state irreducible aperiodic stochastic kernel, "
                    "but its projected process is not a certified Markov "
                    "quotient because the coordinate map fails the recorded "
                    "exact count-kernel strong-lumpability test. "
                    "Its computed full-support stationary law and sampled "
                    "relative-entropy descent support a nonreversible "
                    "finite H-theorem branch. It does not identify that law "
                    "with the state optimizer's source reference, and the "
                    "record-family conservation check is vacuous because "
                    "only one record-family value occurs in this run"
                ),
            },
            "rows": coarsening_rows,
        },
        "recurrent_class_audit": {
            "fine_quotient_closed_class_count": len(recurrent_classes),
            "closed_class_periods": [exact.irreducible_period(
                [[raw[i][j] for j in cls] for i in cls]) for cls in recurrent_classes],
            "closed_class_sizes": [
                len(cls) for cls in recurrent_classes
            ],
            "closed_class_state_indices": recurrent_classes,
            "closed_class_labels": recurrent_labels,
            "nontrivial_recurrent_restriction_available": any(
                len(cls) > 1 for cls in recurrent_classes
            ),
            "verdict": (
                "the fine 20-state chain's only closed communicating class "
                "is one absorbing freezeout state. Its recurrent-class "
                "restriction is a trivial one-state equilibrium, not a "
                "nontrivial thermodynamic relaxation realization; this does "
                "not classify recurrent classes of other quotient maps"
            ),
        },
        "inherited_blockers": inherited_blockers,
        "verdict": {
            "stochasticity": "pass",
            "reversibilized_detailed_balance": "pass",
            "reversibilized_kl_descent": "pass",
            "protected_datum_preservation": (
                "zero off-fibre leakage, but only one record-family value "
                "occurs, so this conservation check is exact and vacuous"
            ),
            "equal_fibre_row": (
                "measured deviation; the kernel is Metropolis repair, "
                "and coincidence with conditional resampling is a "
                "hypothesis to test, never an assumption"
            ),
            "receipt_state": (
                "open; the exhaustive declared coordinate-projection audit "
                "finds a nontrivial ergodic repair-load count kernel, but it "
                "is nonreversible, fails the fine-chain strong-lumpability "
                "test for the exact normalized counts, lacks a "
                "nontrivial protected charge, and has no source "
                "identification with the state optimizer's reference. The "
                "fine chain's only closed recurrent class is a singleton "
                "freezeout state. This rules out closure through the "
                "declared coordinate grammar and the fine-chain recurrent "
                "restriction only; nonlinear, statistical, stochastic, "
                "history-dependent, weakly lumpable, and other "
                "partition-based quotient maps are not excluded"
            ),
        },
    }
    body["receipt_sha256"] = tagged_sha256(canonical_json_bytes(body))
    return body


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    probe = build_probe()
    if args.write:
        PROBE_PATH.parent.mkdir(parents=True, exist_ok=True)
        PROBE_PATH.write_bytes(canonical_json_bytes(probe) + b"\n")
        print(f"wrote {PROBE_PATH}")
    print(probe["status"])
    print(
        "off-fibre visit-weighted mass:",
        probe["measurements"]["off_fibre_mass_visit_weighted"],
    )
    print(
        "pairwise fibre-row TV max:",
        probe["measurements"]["equal_fibre_row_pairwise_tv_max"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
