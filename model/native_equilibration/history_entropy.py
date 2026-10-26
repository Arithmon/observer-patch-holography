"""Counted entropy and charge readout on one continuing native history.

The source has binary load two on twelve ports. A bounded, record-adaptive
prefix may cross the protected core. A two-edge census then resolves the core
from E contrasts with certainty, or C currents on a state-independent success
flag. Every native
coin, current/contrast, disturbed endpoint and timeout remains in the instrument.

For a complete terminal tape b, native branches are injective translations.
If b resolves the initial core charge q, its d_b compatible endpoints have
reference path weight w_q 2**(-T_b). Counting endpoints independently therefore
gives, for a state-independent retained flag of probability beta,

 F(p) = H(B|flag) + E_p log(d_B) - E_p[T] log(2) - log(beta),
 F(p)-S(sigma) = h1 delta<m> + h2 delta<c> - D(B_p || B_sigma).

The fixed-moment first variation vanishes; its Hessian is the negative
classical Fisher form of the actual record response. The cost subtraction is
essential even for bounded stopping times. It accounts for retained fair-coin
entropy, not a physical tape erasure. The successful census is an operational
readout, not supplied access to source occupations. Timeout fibres need not
resolve q; their full-instrument counting excess is bounded explicitly.

Completion kernels for different caps are marginals of the SAME conditional
endpoint distribution. Their overlap agreement is a finite consensus result.
Unmeasured cap occupations in that kernel are predictions, not extra records.
No area, proper volume, stress tensor, geometric boost or Einstein stationarity
is identified by this finite population theorem.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import math

from . import common_variation as v
from ..native_modular_reference import local_preparation as local
from ..observer import native, quantum_mixing

TREE = ((0, 1), (0, 5))
COMPONENT = (0, 1, 5)


@dataclass(frozen=True)
class Path:
    initial: int
    endpoint: tuple[int, ...]
    tape: tuple[tuple[int, int, int, int], ...]
    chance: F
    terminal: str


def current(record, law):
    """Decode the signed transfer from the declared primitive record."""
    _, _, coin, label = record
    if law == "C":
        return label
    if law == "E":
        return (label-(2*coin-1)*(abs(label) % 2))//2
    raise ValueError("original C or E record law required")


def _step(state, tape, seam, coin, law):
    a, b = seam
    out, (_, transfer) = native.native(state, seam, coin)
    label = transfer if law == "C" else state[a]-state[b]
    return out, (*tape, (a, b, coin, label))


def adaptive_seam(tape, law):
    """A declared causal controller, using records only, not hidden loads."""
    if not tape:
        return (0, 5)
    if len(tape) % 2:
        return (0, 1) if current(tape[-1], law) else (1, 8)
    return (0, 7) if current(tape[-1], law) else (2, 10)


def edge_words(deadline):
    if type(deadline) is not int or not 2 <= deadline <= 16:
        raise ValueError("integer census deadline from 2 through 16 required")
    return tuple(((*((coin,)*(length-1)), 1-coin), "success")
                 for length in range(2, deadline+1) for coin in (0, 1)) + tuple(
                     ((coin,)*deadline, "timeout") for coin in (0, 1))


def prefix_length(tape, law, limit, stop_on_current=False):
    """Replay the declared bounded stopping decision from committed records."""
    if stop_on_current:
        for index, record in enumerate(tape[:limit]):
            if current(record, law):
                return index+1
    return limit


@lru_cache(maxsize=12)
def histories(law="C", prefix=2, deadline=3, stop_on_current=False):
    """Enumerate the complete finite instrument; no terminal branch is dropped."""
    if law not in ("C", "E") or type(prefix) is not int or not 0 <= prefix <= 5:
        raise ValueError("C/E law and a bounded zero-to-five attempt prefix required")
    words = edge_words(deadline) if law == "C" else (((0,), "success"), ((1,), "success"))
    allowed = set(quantum_mixing.seams())
    terminal = []
    for initial, state in enumerate(native.states(binary=True)):
        live, stopped = [(state, ())], []
        for _ in range(prefix):
            continuing = []
            for end, tape in live:
                if stop_on_current and tape and current(tape[-1], law):
                    stopped.append((end, tape))
                    continue
                for coin in (0, 1):
                    seam = adaptive_seam(tape, law)
                    if seam not in allowed:
                        raise ArithmeticError("controller selected a non-source seam")
                    continuing.append(_step(end, tape, seam, coin, law))
            live = continuing
        for end, tape in (*stopped, *live):
            live = [(end, tape)]
            for edge in TREE:
                continuing = []
                for end, tape in live:
                    for word, flag in words:
                        out, record = end, tape
                        for coin in word:
                            out, record = _step(out, record, edge, coin, law)
                        if flag == "timeout":
                            terminal.append(Path(initial, out, record,
                                                 F(1, 2**len(record)), "timeout"))
                        else:
                            continuing.append((out, record))
                live = continuing
            terminal.extend(Path(initial, end, tape, F(1, 2**len(tape)), "success")
                            for end, tape in live)
    return tuple(terminal)


def decode_core(tape, prefix, law="C"):
    """Recover initial core bits from a successful native tape, without a state."""
    position = prefix
    displacement = {i: 0 for i in COMPONENT}
    differences = []
    for edge in TREE:
        if position >= len(tape) or tape[position][:2] != edge:
            raise ValueError("incomplete native census")
        start = position
        if law == "E":
            contrast = tape[start][3]
        else:
            first_coin = tape[start][2]
            while position < len(tape) and tape[position][2] == first_coin:
                if tape[position][:2] != edge:
                    raise ValueError("missing native coin flip")
                position += 1
            if position >= len(tape) or tape[position][:2] != edge:
                raise ValueError("census timeout does not resolve the core")
            first = current(tape[start], law)
            flip = current(tape[position], law)
            direction = 2*first_coin-1
            contrast = first if first else (direction if flip else 0)
        a, b = edge
        differences.append((a, b, contrast-displacement[a]+displacement[b]))
        transfer = sum(current(record, law) for record in tape[start:position+1])
        displacement[a] -= transfer
        displacement[b] += transfer
        position += 1
    if position != len(tape):
        raise ValueError("unexpected records after the completed census")
    potential = {0: 0}
    for a, b, contrast in differences:
        potential[b] = potential[a]-contrast
    # A three-port component cannot be all occupied on the binary total-two
    # source. Its lowest occupation is therefore zero, fixing the offset.
    offset = -min(potential.values())
    initial = {i: potential[i]+offset for i in COMPONENT}
    if any(x not in (0, 1) for x in initial.values()) or sum(initial.values()) > 2:
        raise ValueError("records do not describe the binary total-two source")
    core = [initial[0], initial[1]]
    for record in tape[:prefix]:
        a, b = record[:2]
        transfer = current(record, law)
        for i in (0, 1):
            core[i] -= transfer*((b == i)-(a == i))
    if any(x not in (0, 1) for x in core):
        raise ValueError("prefix records do not recover a binary initial core")
    return tuple(core)


def fibres(paths, cap=()):
    """Group actual histories, optionally refine by a predicted endpoint cap."""
    if tuple(sorted(set(cap))) != cap or any(type(i) is not int or not 0 <= i < 12 for i in cap):
        raise ValueError("sorted distinct cap indices required")
    result = defaultdict(list)
    for path in paths:
        result[(path.tape, tuple(path.endpoint[i] for i in cap))].append(path)
    for row in result.values():
        if len({p.endpoint for p in row}) != len(row) or len({p.chance for p in row}) != 1:
            raise ArithmeticError("a native tape lost its injective, equal-coin fibre")
    return dict(result)


def flag_probability(paths):
    by_input = [F() for _ in range(66)]
    for path in paths:
        by_input[path.initial] += path.chance
    if min(by_input) <= 0 or len(set(by_input)) != 1:
        raise ValueError("retained flag must have positive state-independent probability")
    return by_input[0]


def _sum_forms(*forms):
    result = defaultdict(F)
    for form in forms:
        for prime, coefficient in form.items():
            result[int(prime)] += coefficient
    return v._clean(result)


def source_entropy(weights):
    v.moments(weights)
    result = defaultdict(F)
    for weight in weights:
        if weight:
            v._log_add(result, -weight, weight)
    return v._clean(result)


def ledger(paths, weights, cap=()):
    """Exact prime-log entropy ledger; edge dimensions are counted first."""
    v.moments(weights)
    beta = flag_probability(paths)
    rows = fibres(paths, cap)
    bulk, edge, cost, normalization = (defaultdict(F) for _ in range(4))
    probabilities = {}
    mean_attempts = F()
    for key, row in rows.items():
        probability = sum((weights[p.initial]*p.chance for p in row), F())/beta
        probabilities[key] = probability
        if probability:
            v._log_add(bulk, -probability, probability)
            v._log_add(edge, probability, len(row))
        mean_attempts += probability*len(row[0].tape)
    v._log_add(cost, -mean_attempts, 2)
    v._log_add(normalization, -1, beta)
    return {"probabilities": probabilities, "bulk": v._clean(bulk),
            "edge": v._clean(edge), "coin_cost": v._clean(cost),
            "flag_normalization": v._clean(normalization),
            "functional": _sum_forms(bulk, edge, cost, normalization),
            "flag_probability": beta, "mean_native_attempts": mean_attempts}


def relative_entropy(probabilities, reference):
    result = defaultdict(F)
    if probabilities.keys() != reference.keys():
        raise ValueError("the same record outcomes are required")
    for key, probability in probabilities.items():
        if probability:
            v._log_add(result, probability, probability/reference[key])
    return v._clean(result)


def conditional_cap(row, cap):
    """The common reference's uniform conditional endpoint marginal."""
    if len({v.charges()[p.initial] for p in row}) != 1:
        raise ValueError("counted reference completion needs a charge-resolved fibre")
    result = defaultdict(F)
    for path in row:
        result[tuple(path.endpoint[i] for i in cap)] += F(1, len(row))
    return dict(result)


def value(form):
    return sum(float(coefficient)*math.log(int(prime)) for prime, coefficient in form.items())


def run():
    weights = local.reference()
    target = source_entropy(weights)
    results = {}
    for law in ("C", "E"):
        paths = histories(law)
        kept = tuple(p for p in paths if p.terminal == "success")
        beta = flag_probability(kept)
        if flag_probability(paths) != 1 or beta != (F(9, 16) if law == "C" else 1):
            raise ArithmeticError("native census terminal probability changed")
        for path in kept:
            if decode_core(path.tape, 2, law) != native.states(binary=True)[path.initial][:2]:
                raise ArithmeticError("native record core decoder failed")
        rows = fibres(kept)
        observed = ledger(kept, weights)
        if observed["functional"] != target:
            raise ArithmeticError("record-only counted entropy envelope failed")
        # All cap comparisons refine this same history measure. They neither
        # prepare a fresh source nor independently resample an overlap state.
        caps = ((), (0, 1), (0, 1, 5), (0, 1, 2, 5), (0, 1, 5, 8), tuple(range(12)))
        for cap in caps:
            if ledger(kept, weights, cap)["functional"] != target:
                raise ArithmeticError("same-history cap envelope changed")
        complete = ledger(paths, weights)
        excess = _sum_forms(complete["functional"], {k: -x for k, x in target.items()})
        bound = (1-beta)*math.log(float(max(weights)/min(weights)))
        if (law == "C" and not 0 < value(excess) < bound) or (law == "E" and excess):
            raise ArithmeticError("timeout ambiguity control lost its strict excess or bound")
        results[law] = {"terminal_paths": len(paths), "successful_paths": len(kept),
            "record_only_success_fibres": len(rows),
            "success_probability": beta, "timeout_probability": 1-beta,
            "mean_success_attempts": observed["mean_native_attempts"],
            "mean_complete_attempts": complete["mean_native_attempts"],
            "exact_source_entropy": target, "exact_success_functional": observed["functional"],
            "complete_counting_excess": excess,
            "complete_excess_upper_coefficient": 1-beta,
            "complete_excess_upper_log_argument": max(weights)/min(weights),
            "successful_completion_dimensions": sorted({len(row) for row in rows.values()}),
            "same_history_cap_refinements_checked": len(caps)}
    stopped = histories("E", prefix=3, stop_on_current=True)
    for path in stopped:
        count = prefix_length(path.tape, "E", 3, True)
        if decode_core(path.tape, count, "E") != native.states(binary=True)[path.initial][:2]:
            raise ArithmeticError("record-stopped prefix decoder failed")
    stopped_ledger = ledger(stopped, weights)
    if stopped_ledger["functional"] != target:
        raise ArithmeticError("record-dependent native stopping cost was lost")
    per_input_cost = [sum((p.chance*len(p.tape) for p in stopped if p.initial == i), F())
                      for i in range(66)]
    h = 24
    alpha = 1-F(1, 2**(h-1))
    finite = v.stationarity_budget("local")
    return {"laws": results, "prefix_attempts": 2, "census_deadline_tested": 3,
        "bounded_C_readout_service": {"deadline": h, "native_attempts_upper": 2+2*h,
            "timeout_probability": 1-alpha**2,
            "mean_attempts_upper": 8,
            "timeout_counting_excess_upper_log_argument": max(weights)/min(weights),
            "timeout_counting_excess_upper_coefficient": 1-alpha**2},
        "unconditional_E_readout_service": {"native_attempts": 4, "timeout_probability": 0,
            "core_census_attempts": 2, "prefix_attempts": 2},
        "record_stopped_E_service": {"prefix_attempts_upper": 3, "total_attempts_upper": 5,
            "conditional_input_mean_attempts_range": (min(per_input_cost), max(per_input_cost)),
            "mean_native_attempts": stopped_ledger["mean_native_attempts"],
            "exact_cost_corrected_functional": stopped_ledger["functional"]},
        "finite_preparation_stationarity": {
            "operator_sandwich_radius": finite["operator_sandwich_radius"],
            "fixed_moment_tangent_error_per_trace_norm_upper":
                finite["entropy_tangent_error_per_trace_norm_upper"],
            "local_reference_native_attempts_upper": finite["native_total_attempts_upper"],
            "complete_flagged_reference_halftrace_upper": finite["complete_flagged_reference_halftrace_upper"],
            "adaptive_record_channel_preserves_the_reference_sandwich": True,
            "input_dependent_stopping_cost_must_be_varied": True},
        "scope": {"readout_uses_only_declared_native_records": True,
            "prefix_has_record_feedback_and_core_crossings": True,
            "all_timeout_disturbed_outputs_retained": True,
            "cap_completions_share_one_conditional_endpoint_measure": True,
            "unmeasured_cap_predictions_are_additional_records": False,
            "coherent_reference_entropy_theorem": False,
            "proper_volume_area_stress_dictionary": False}}
