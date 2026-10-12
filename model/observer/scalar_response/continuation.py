"""Two field reads and three native repairs on one retained joint state.

The main experiment prepares the field and source once. Classical controller
records persist while the conditional field, including its first measurement
backaction, supplies the second observation. Explicit replacement-field and
disconnected-actuator controls are separate experiments.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
from itertools import product
import json
import math
import random

from model.observer import native
from .experiment import canonical, decoder_seam, mean_summary, sources


GRIDS = (13, 34, 89, 233)
TIMES = (Q(19, 20), Q(1))
GAMMA = 4
TRIALS = 50_000
SEED = 20261012


@dataclass(frozen=True)
class Record:
    """A bounded append-only tape: native, field, native, field, native."""

    prefix: tuple = ()

    def __post_init__(self):
        if type(self.prefix) is not tuple or len(self.prefix) > 5:
            raise ValueError("a tuple of at most five committed records is required")
        for index, value in enumerate(self.prefix):
            if index % 2:
                valid = type(value) is int and value in (-1, 1)
            else:
                valid = (type(value) is tuple and len(value) == 2
                         and all(type(x) is int for x in value)
                         and value in native.LETTERS)
            if not valid:
                raise ValueError("invalid alternating native/field record")

    @property
    def phase(self):
        return len(self.prefix)

    @property
    def halted(self):
        return self.phase == 5

    def append(self, value):
        if self.halted:
            raise ValueError("the full record tape cannot be appended")
        return Record(self.prefix + (value,))

    def checkpoint(self):
        """The controller checkpoint contains no source or field state."""
        return {"phase": self.phase,
                "prefix": [list(x) if type(x) is tuple else x for x in self.prefix]}

    @classmethod
    def restore(cls, payload):
        if (type(payload) is not dict or set(payload) != {"phase", "prefix"}
                or type(payload["phase"]) is not int
                or type(payload["prefix"]) is not list):
            raise ValueError("invalid controller checkpoint")
        record = cls(tuple(tuple(x) if type(x) is list else x for x in payload["prefix"]))
        if record.phase != payload["phase"]:
            raise ValueError("checkpoint phase disagrees with the tape")
        return record


def selected_seam(record):
    """Actuator access: phase and the latest committed detector bit only."""
    if record.phase == 0:
        return native.E
    if record.phase in (2, 4):
        return decoder_seam(record.prefix[-1])
    return None


def native_step(loads, record, coin):
    """Advance the same source; the controller receives no loads as input."""
    if (type(loads) is not tuple or len(loads) != 12
            or any(type(x) is not int or x < 0 for x in loads) or sum(loads) != 2):
        raise ValueError("a total-two native source is required")
    if record.halted:
        return loads, record
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError("a binary native coin is required")
    seam = selected_seam(record)
    if seam is None:
        raise ValueError("the next event is a field observation")
    after, letter = native.native(loads, seam, coin)
    return after, record.append(letter)


@dataclass(frozen=True)
class Path:
    initial: tuple
    after_first: tuple
    after_second: tuple
    final: tuple
    record: Record


@lru_cache(maxsize=1)
def native_paths():
    """All labelled paths; field probabilities are supplied by actual branches."""
    rows = []
    for initial in native.states(binary=True):
        for c1, a, c2, b, c3 in product((0, 1), (-1, 1), (0, 1), (-1, 1), (0, 1)):
            first, r1 = native_step(initial, Record(), c1)
            second, r2 = native_step(first, r1.append(a), c2)
            final, r3 = native_step(second, r2.append(b), c3)
            rows.append(Path(initial, first, second, final, r3))
    return tuple(rows)


@lru_cache(maxsize=1)
def record_codes():
    """Lossless eight-bit representation of the reached controller prefixes."""
    prefixes = {path.record.prefix[:n] for path in native_paths() for n in range(6)}
    return tuple(sorted(prefixes, key=lambda prefix: (len(prefix), prefix)))


def encode_record(record):
    try:
        return record_codes().index(record.prefix)
    except ValueError as exc:
        raise ValueError("record is outside the declared preparation support") from exc


def decode_record(code):
    if type(code) is not int or not 0 <= code < len(record_codes()):
        raise ValueError("invalid reached controller code")
    return Record(record_codes()[code])


def exact_statistics():
    """Exact native coefficients for an arbitrary parity-symmetric field law.

    p_j(a,b) = [1+j*mu1*a+j*mu2*b+kappa_abs(j)*a*b]/4.
    Each complete path has probability p_j(a,b)/528.
    """
    functions = {
        "second_current": lambda row: row.record.prefix[2][1],
        "third_current": lambda row: row.record.prefix[4][1],
        "final_port0": lambda row: row.final[0],
    }
    coefficients = {}
    for name, function in functions.items():
        values = dict.fromkeys(("constant", "mu1", "mu2", "kappa0", "kappa1"), Q())
        for row in native_paths():
            r1, a, _, b, _ = row.record.prefix
            j, value = r1[1], function(row)
            values["constant"] += Q(value, 2112)
            values["mu1"] += Q(value*a*j, 2112)
            values["mu2"] += Q(value*b*j, 2112)
            values["kappa0" if j == 0 else "kappa1"] += Q(value*a*b, 2112)
        coefficients[name] = {key: str(value) for key, value in values.items()}
    expected = {
        "second_current": ("0", "5/264", "0", "0", "0"),
        "third_current": ("0", "-5/352", "5/352", "1/264", "-1/264"),
        "final_port0": ("1/6", "-5/1056", "-5/352", "-1/264", "1/264"),
    }
    if any(tuple(coefficients[key].values()) != value for key, value in expected.items()):
        raise ArithmeticError("native continuation coefficients disagree")
    counts = [sum(len(prefix) == phase for prefix in record_codes()) for phase in range(6)]
    if counts != [1, 4, 8, 24, 48, 124] or len(native_paths()) != 2112:
        raise ArithmeticError("native continuation reachability changed")
    allocated = [1, 6, 12, 72, 144, 864]
    words = {row.record.prefix for row in native_paths()}
    return {
        "labelled_histories": 2112,
        "path_probability": "p_j(y1,y2)/528",
        "moment_coefficients": coefficients,
        "reachable_prefix_counts": counts,
        "reachable_prefix_total": len(record_codes()),
        "controller_bits": 8,
        "terminal_record_words": 124,
        "terminal_fixed_width_bits": 7,
        "terminal_record_log_count_nats": math.log(124),
        "allocated_record_sectors": allocated,
        "source_times_record_dimension_excluding_field": 78*sum(allocated),
        "coarse_terminal_record_counts": {
            "native_letters_only": len({(w[0], w[2], w[4]) for w in words}),
            "latest_detector_and_native_letter": len({w[3:] for w in words}),
            "detector_pair": len({(w[1], w[3]) for w in words}),
            "delete_first_detector_only": len({(w[0], w[2], w[3], w[4]) for w in words}),
            "complete_erasure": 1,
        },
        "prefix_encoding_sha256": hashlib.sha256(canonical(record_codes()).encode()).hexdigest(),
    }


def field_laws(field):
    """Execute both instruments, including checkpointed conditional continuation.

    The independent closed two-read formula is deliberately not used here.
    The reset control returns to the unmeasured preparation in the same
    interaction-picture kernel; only this explicitly named control replaces it.
    """
    from .kernel import WeylKernel, WeylState

    joint, reset, donor, independent, commuting = {}, {}, {}, {}, {}
    descriptors, norm_errors, checkpoint_errors = [], [], []
    base = WeylKernel.from_field(field, times=tuple(map(float, TIMES)), current=1,
                                gamma=GAMMA, couplings=(1, 1))
    for current in (-1, 0, 1):
        kernel = base.for_current(current)
        classical = WeylKernel(base.covariance.real, base.unit_means, base.times,
                               current, base.mode_count)
        original = kernel.vacuum()
        first_probabilities, second_marginal = {}, defaultdict(float)
        free_second = {b: original.instrument(1, b)[0] for b in (-1, 1)}
        for a in (-1, 1):
            p1, branch = original.instrument(0, a)
            classical_p1, classical_branch = classical.vacuum().instrument(0, a)
            first_probabilities[a] = p1
            checkpoint = json.loads(json.dumps(branch.checkpoint(), allow_nan=False))
            restored = WeylState.restore(kernel, checkpoint)
            norm_errors.extend((abs(branch.norm_squared()-1), abs(restored.norm_squared()-1)))
            for b in (-1, 1):
                p2, final = branch.instrument(1, b)
                restored_p2, restored_final = restored.instrument(1, b)
                checkpoint_errors.append(abs(p2-restored_p2))
                norm_errors.extend((abs(final.norm_squared()-1), abs(restored_final.norm_squared()-1)))
                joint[current, a, b] = p1*p2
                second_marginal[b] += p1*p2
                reset[current, a, b] = p1*free_second[b]
                independent[current, a, b] = p1/2
                commuting[current, a, b] = classical_p1*classical_branch.instrument(1, b)[0]
                terminal = final.checkpoint()
                # This execution binding is platform-dependent floating data.
                # Its round trip was tested above; do not freeze its digest.
                terminal.pop("kernel")
                descriptors.append({"first_current": current, "first_detector": a,
                                    "second_detector": b, "first_probability": p1,
                                    "conditional_second_probability": p2,
                                    "joint_probability": p1*p2,
                                    "conditional_field_descriptor": terminal})
        for a, b in product((-1, 1), repeat=2):
            donor[current, a, b] = first_probabilities[a]*second_marginal[b]
    laws = {"retained": joint, "conditional_field_donor": donor,
            "field_reset": reset, "independent_second_bit": independent,
            "commuting_gaussian": commuting}
    for name, law in laws.items():
        for current in (-1, 0, 1):
            probabilities = [law[current, a, b] for a, b in product((-1, 1), repeat=2)]
            if min(probabilities) <= 0 or abs(math.fsum(probabilities)-1) > 2e-12:
                raise ArithmeticError(f"invalid complete positive outcome law: {name}")
    if max(norm_errors+checkpoint_errors) > 2e-11:
        raise ArithmeticError("conditional field checkpoint or normalization failed")
    return laws, {
        "maximum_conditional_norm_error": max(norm_errors),
        "maximum_checkpoint_continuation_probability_error": max(checkpoint_errors),
        "kernel_statistics": base.statistics(),
        "terminal_field_branches": descriptors,
    }


def evaluate(law, disconnected=False):
    """Enumerate every path and forecasts issued before the final native event."""
    totals = defaultdict(float)
    full, coarse = defaultdict(lambda: [0., 0.]), defaultdict(lambda: [0., 0.])
    for row in native_paths():
        r1, a, r2, b, r3 = row.record.prefix
        weight = law[r1[1], a, b]/528
        final, last = (native.native(row.after_second, native.E, r3[0])
                       if disconnected else (row.final, r3))
        event = int(last[1] > 0)
        for key, value in (("mass", 1), ("second_current", r2[1]),
                           ("third_current", last[1]), ("final_port0", final[0]),
                           ("second_detector_plus", b == 1),
                           ("detector_product", a*b),
                           ("zero_current_detector_product", (r1[1] == 0)*a*b)):
            totals[key] += weight*value
        for target, key in ((full, row.record.prefix[:4]), (coarse, (r2, b))):
            target[key][0] += weight
            target[key][1] += weight*event
    full_forecast = {key: event/mass for key, (mass, event) in full.items()}
    coarse_forecast = {key: event/mass for key, (mass, event) in coarse.items()}
    full_risk = math.fsum(event*(1-event/mass) for mass, event in full.values())
    coarse_risk = math.fsum(event*(1-event/mass) for mass, event in coarse.values())
    # Conditional on a fixed first current, both predictors know the same
    # preparation. Deleting y1 changes information, not the physical field.
    field_full_risk = field_coarse_risk = 0.
    for j, prior in ((-1, 5/66), (0, 28/33), (1, 5/66)):
        plus = math.fsum(law[j, a, 1] for a in (-1, 1))
        field_coarse_risk += prior*plus*(1-plus)
        for a in (-1, 1):
            mass = law[j, a, -1]+law[j, a, 1]
            event = law[j, a, 1]
            field_full_risk += prior*event*(1-event/mass)
    totals["zero_current_conditional_detector_product"] = totals["zero_current_detector_product"]/(28/33)
    totals.update({"native_full_history_brier": full_risk,
                   "native_latest_pair_brier": coarse_risk,
                   "native_record_deletion_brier_gain": coarse_risk-full_risk,
                   "native_within_latest_pair_shuffle_brier": 2*coarse_risk-full_risk,
                   "field_full_record_brier": field_full_risk,
                   "field_first_bit_deleted_brier": field_coarse_risk,
                   "field_record_deletion_brier_gain": field_coarse_risk-field_full_risk})
    return dict(totals), full_forecast, coarse_forecast


def _second_plus(law, current, first_bit):
    return law[current, first_bit, 1]/sum(law[current, first_bit, b] for b in (-1, 1))


def sample(laws, predictions, trials=TRIALS, seed=SEED):
    if type(trials) is not int or trials < 2 or type(seed) is not int:
        raise ValueError("at least two trials and an integer seed are required")
    rng, totals, squares, words = random.Random(seed), Counter(), Counter(), Counter()
    configurations, traces = native.states(binary=True), []
    actual = laws["retained"]
    _, full_forecasts, coarse_forecasts = evaluate(actual)
    for trial in range(trials):
        loads, record = native_step(rng.choice(configurations), Record(), rng.randrange(2))
        j = record.prefix[0][1]
        plus1 = sum(actual[j, 1, b] for b in (-1, 1))
        a = 1 if rng.random() < plus1 else -1
        loads, record = native_step(loads, record.append(a), rng.randrange(2))
        u, final_coin = rng.random(), rng.randrange(2)
        bits = {name: 1 if u < _second_plus(law, j, a) else -1
                for name, law in laws.items()}
        b = bits["retained"]
        before_final = record.append(b)
        # Both forecasts are read before the third source update.
        full_forecast = full_forecasts[before_final.prefix]
        coarse_forecast = coarse_forecasts[(record.prefix[2], b)]
        final, complete = native_step(loads, before_final, final_coin)
        values = {"third_current": complete.prefix[4][1], "final_port0": final[0],
                  "second_detector_plus": int(b == 1), "detector_product": a*b,
                  "zero_current_detector_product": (j == 0)*a*b}
        for name in ("conditional_field_donor", "field_reset", "independent_second_bit", "commuting_gaussian"):
            _, control = native_step(loads, record.append(bits[name]), final_coin)
            values[name] = control.prefix[4][1]
        values["disconnected"] = native.native(loads, native.E, final_coin)[1][1]
        values["paired_current_minus_field_donor"] = values["third_current"]-values["conditional_field_donor"]
        event = int(complete.prefix[4][1] > 0)
        values["native_full_history_brier"] = (event-full_forecast)**2
        values["native_latest_pair_brier"] = (event-coarse_forecast)**2
        for name, value in values.items():
            totals[name] += value
            squares[name] += value*value
        if native_step(final, complete, 1-final_coin) != (final, complete):
            raise ArithmeticError("halt changed a completed history")
        if decode_record(encode_record(complete)) != complete:
            raise ArithmeticError("compact controller lost a record")
        words[complete.prefix] += 1
        if trial < 8:
            traces.append({"record": complete.checkpoint(),
                           "controller_code": encode_record(complete),
                           "forecasts_issued_before_event_3": {"full_history": full_forecast,
                                                                "latest_pair": coarse_forecast},
                           "terminal_field_branch": {"first_current": j, "outcomes": [a, b]},
                           "halted": True})
    expected = {key: predictions["retained"][key] for key in
                ("third_current", "final_port0", "second_detector_plus", "detector_product",
                 "zero_current_detector_product", "native_full_history_brier", "native_latest_pair_brier")}
    expected.update({key: predictions[key]["third_current"] for key in
                     ("conditional_field_donor", "field_reset", "independent_second_bit", "commuting_gaussian", "disconnected")})
    expected["paired_current_minus_field_donor"] = expected["third_current"]-expected["conditional_field_donor"]
    return {"trials": trials, "seed": seed,
            "statistics": {key: mean_summary(totals[key], squares[key], trials, value)
                           for key, value in expected.items()},
            "record_count_encoding": "[controller_code, count]; exact prefix generator and its digest define the code",
            "record_counts": sorted([encode_record(Record(word)), count] for word, count in words.items()),
            "first_eight_histories": traces}


def certificate_checks(q, kernel_statistics, moments, certificate):
    """Compare independently executed diagnostics with ideal-model enclosures.

    The small comparison tolerance only accommodates binary64 evaluation;
    it does not promote that evaluation to a directed roundoff certificate.
    """
    level = next((row for row in certificate["levels"] if row["q"] == q), None)
    if level is None:
        raise ValueError("the numerical grid has no matching mathematical certificate")
    witnesses = level["finite_field_witnesses"]
    real, imaginary = kernel_statistics["covariance_real"], kernel_statistics["covariance_imag"]
    def contains(interval, value):
        low, high = map(lambda x: float(Q(x)), interval)
        return math.isfinite(value) and low-2e-12 <= value <= high+2e-12
    def above(bound, value):
        return math.isfinite(value) and value+2e-12 >= float(Q(bound))
    checks = {
        "kernel_matches_certified_grid_and_times": (
            kernel_statistics["mode_count"] == (q-1)**3
            and kernel_statistics["times"] == list(map(float, TIMES))),
        "both_variances_inside_rational_enclosure": all(
            contains(witnesses["variance_interval"], real[index][index]) for index in (0, 1)),
        "cross_covariance_inside_rational_enclosure": contains(
            witnesses["real_cross_covariance_interval"], real[0][1]),
        "commutator_inside_rational_enclosure": contains(
            witnesses["commutator_over_i_interval"], 2*imaginary[0][1]),
        "vacuum_memory_above_certified_lower_bound": above(
            witnesses["vacuum_bit_correlation_lower"], moments["kappa0"]),
        "second_plus_backaction_above_certified_lower_bound": above(
            witnesses["second_plus_disturbance_given_positive_current_lower"],
            moments["commutator_backaction_signed_mean_difference"]/2),
    }
    for name in ("vacuum_bit_correlation_lower", "second_plus_disturbance_given_positive_current_lower"):
        if q in (89, 233):
            checks[f"{name}_strictly_positive"] = Q(witnesses[name]) > 0
    return checks


def run(trials=TRIALS, seed=SEED, grids=GRIDS):
    from .field import ScalarField
    from . import continuation_bounds

    if type(trials) is not int or trials < 2 or type(seed) is not int:
        raise ValueError("at least two trials and an integer seed are required")
    if (type(grids) not in (tuple, list) or not grids
            or any(type(q) is not int or q not in GRIDS for q in grids)
            or len(set(grids)) != len(grids)):
        raise ValueError("choose distinct grids from the fixed certified family")
    sources()
    bounds = continuation_bounds.run()
    exact, executions, checks = exact_statistics(), [], {}
    for q in grids:
        field = ScalarField(q)
        laws, continuation_checks = field_laws(field)
        predictions = {name: evaluate(law)[0] for name, law in laws.items()}
        predictions["disconnected"] = evaluate(laws["retained"], disconnected=True)[0]
        law = laws["retained"]
        moments = {"mu1": sum(a*law[1, a, b] for a, b in product((-1, 1), repeat=2)),
                   "mu2": sum(b*law[1, a, b] for a, b in product((-1, 1), repeat=2)),
                   "kappa0": sum(a*b*law[0, a, b] for a, b in product((-1, 1), repeat=2)),
                   "kappa1": sum(a*b*law[1, a, b] for a, b in product((-1, 1), repeat=2))}
        moments["second_signed_mean_commuting_gaussian"] = sum(
            b*laws["commuting_gaussian"][1, a, b] for a, b in product((-1, 1), repeat=2))
        moments["commutator_backaction_signed_mean_difference"] = (
            moments["mu2"]-moments["second_signed_mean_commuting_gaussian"])
        checks.update({f"q{q}_{name}": value for name, value in certificate_checks(
            q, continuation_checks["kernel_statistics"], moments, bounds).items()})
        for name, coefficients in exact["moment_coefficients"].items():
            expected = float(Q(coefficients["constant"])) + math.fsum(
                float(Q(coefficients[key]))*moments[key] for key in ("mu1", "mu2", "kappa0", "kappa1"))
            checks[f"q{q}_{name}_matches_exact_native_coefficients"] = abs(expected-predictions["retained"][name]) < 2e-12
        for name, prediction in predictions.items():
            checks[f"q{q}_{name}_all_history_mass"] = abs(prediction["mass"]-1) < 2e-12
        checks[f"q{q}_earlier_record_improves_field_forecast"] = predictions["retained"]["field_record_deletion_brier_gain"] > 0
        checks[f"q{q}_earlier_records_improve_native_forecast"] = predictions["retained"]["native_record_deletion_brier_gain"] > 0
        checks[f"q{q}_field_donor_removes_zero_current_memory"] = abs(predictions["conditional_field_donor"]["zero_current_conditional_detector_product"]) < 2e-12
        checks[f"q{q}_same_second_action_marginal"] = all(abs(predictions[name]["second_detector_plus"]-.5) < 2e-12 for name in laws)
        sampled = sample(laws, predictions, trials, seed)
        for name, item in sampled["statistics"].items():
            checks[f"q{q}_{name}_sample_agrees"] = abs(item["mean"]-item["forecast"]) <= 6*item["standard_error"]+1/trials
        executions.append({"q": q, "mode_count": field.mode_count, "field_moments": moments,
                           "predictions": predictions, "conditional_continuation": continuation_checks,
                           "sample": sampled})
    result = {"schema": "oph.observer.scalar_continuation.v1", "grids": list(grids),
              "model_times": list(map(str, TIMES)), "gamma": GAMMA,
              "exact_native_statistics": exact, "refinement_bounds": bounds, "executions": executions,
              "resources_per_history": {"source_preparations": 1, "field_preparations": 1,
                                        "native_events": 3, "field_observations": 2,
                                        "classical_record_writes": 5, "source_resets": 0,
                                        "field_resets": 0, "maximum_coherent_terms": 4,
                                        "controller_bits": 8, "halt_after_native_event": 3},
              "scope": {
                  "state": "same native loads and conditional field continue; one initial preparation per complete history",
                  "controller": "phase and latest detector bit select the native seam; forecasts read records only",
                  "memory": "8 bits encode 209 reached classical prefixes; program, source and field are separate resources",
                  "field": "infinite-dimensional oscillator state represented by a bounded four-term Weyl word; no Fock cutoff",
                  "time": "supplied scalar model coordinate; native operation count is not physical duration",
                  "controls": "conditional-field donor and unmeasured-field reset are explicit replacement services; neither occurs in the main process",
                  "commuting_control": "same real Gaussian covariance and means with zero commutator; temporal correlation alone does not distinguish this classical alternative",
                  "deletion": "forecast-only record deletion leaves the actual source, field and action unchanged",
                  "numerics": "floating full-mode execution checked against a separate directed rational two-read ideal-model certificate; numerical roundoff is not certified",
                  "sampling": "50,000-trial default checks broad distribution agreement and does not resolve the small commutator backaction contrast",
                  "refinement": "the attached certificate covers the fixed four-grid family; the grids list identifies the levels actually executed",
                  "halt": "bounded continuing history, followed by halt and reread; no unbounded memory or repeated reset cycle",
              }, "checks": checks, "all_pass": all(checks.values())}
    if not result["all_pass"]:
        raise ArithmeticError("continuation check failed: "+", ".join(key for key, value in checks.items() if not value))
    return result
