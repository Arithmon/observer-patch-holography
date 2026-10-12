#!/usr/bin/env python3
"""Execute reproducible fresh two-event trials of the compact native observer."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent


def dependencies():
    from . import capacity
    source, _ = capacity.dependencies()
    _, _, prefixes, transitions, histories = capacity.enumerate_source(source, True)
    encoding = capacity.encode_prefixes(prefixes, transitions)
    compiled = capacity.compile_controller(source, histories, encoding)
    return capacity, source, compiled



def measured(total, trials, reference):
    mean = Q(total, trials)
    exact = Q(reference)
    return {"sample_mean": str(mean), "sample_decimal": float(mean),
            "exact_ensemble_reference": str(exact), "reference_decimal": float(exact),
            "sample_minus_reference": str(mean - exact)}


def run(trials=20_000, seed=20261012):
    if type(trials) is not int or trials < 1 or type(seed) is not int:
        raise ValueError("Trials must be positive and seed must be an integer")
    compact, source, compiled = dependencies()
    exact = source.run()
    states = source.states()
    binary = [x for x in states if max(x) == 1]
    # These 132 first-event paths implement the exact prior for independent donors.
    donors = [first.prefix[0] for x in binary for _, first in source.tick(source.State(x))]
    donors_by_action = {"E": [r for r in donors if r[1] <= 0],
                        "F": [r for r in donors if r[1] > 0]}
    oracle = {(states[row["initial_state"]], row["first_record"][0], row["second_record"][0]): row
              for row in exact["main_histories"]}
    if len(oracle) != 264 or len(compiled["decode"]) != 16:
        raise ValueError("Exact binary experiment or compact code inventory changed")
    rng = random.Random(seed)
    counts = Counter()
    loss = {key: Q(0) for key in ("native_record", "action_known_prior", "within_action_shuffle")}
    traces = []
    feedback_example = None
    for trial in range(trials):
        initial = rng.choice(binary)  # fresh uniform diagonal preparation
        initial_code = compiled["encode"][()]
        coin1 = rng.randrange(2)
        weight1, code1, loads1 = compact.compact_tick(compiled, initial_code, initial)[coin1]
        view1 = compact.compact_view(compiled, code1)
        record1 = view1["prefix"][0]
        action = "F" if view1["seam"] == source.F else "E"
        # All three forecasts are issued now, before coin2 or event2 exists.
        forecast = view1["forecast"]
        action_forecast = Q(0) if action == "F" else Q(5, 61)
        shuffle_record = rng.choice(donors_by_action[action])
        shuffled_forecast = source.FORECAST[shuffle_record]
        selector_record = rng.choice(donors)
        selector_seam = source.F if selector_record[1] > 0 else source.E
        coin2 = rng.randrange(2)
        weight2, code2, loads2 = compact.compact_tick(compiled, code1, loads1)[coin2]
        view2 = compact.compact_view(compiled, code2)
        record2 = view2["prefix"][1]
        halt = compact.compact_tick(compiled, code2, loads2)
        if weight1 != Q(1, 2) or weight2 != Q(1, 2) or halt != ((Q(1), code2, loads2),):
            raise ValueError("Native branch weights or terminal halt changed")
        expected = oracle[initial, coin1, coin2]
        if (record1 != tuple(expected["first_record"])
                or record2 != tuple(expected["second_record"])
                or loads2 != states[expected["final_state"]]
                or view1["seam"] != tuple(expected["second_seam"])
                or forecast != Q(expected["issued_before_event_2"])):
            raise ValueError(f"Trial {trial} disagrees with the exact path reference")
        # Paired controls share the same first source state and second fair coin.
        disconnected, _ = source.native(loads1, source.E, coin2)
        donor_control, _ = source.native(loads1, selector_seam, coin2)
        outcome = int(record2[1] > 0)
        counts["event"] += outcome
        counts["action_F"] += action == "F"
        counts["donor_action_F"] += selector_seam == source.F
        counts["native_port0"] += loads2[0]
        counts["disconnected_port0"] += disconnected[0]
        counts["donor_port0"] += donor_control[0]
        for name, predicted in (("native_record", forecast), ("action_known_prior", action_forecast),
                                ("within_action_shuffle", shuffled_forecast)):
            loss[name] += (predicted - outcome) ** 2
        trace = {
            "trial": trial + 1,
            "controller_codes": [initial_code, code1, code2],
            "first_record": list(record1), "second_action": action,
            "forecasts_issued_before_event_2": {
                "native_record": str(forecast), "action_known_prior": str(action_forecast),
                "within_action_shuffle": str(shuffled_forecast)},
            "second_record": list(record2), "outward_current_event": bool(outcome),
            "halted": True,
            "external_source_audit": {"initial_occupied_ports": [i for i, x in enumerate(initial) if x],
                                      "final_port0": loads2[0], "disconnected_final_port0": disconnected[0],
                                      "donor_control_final_port0": donor_control[0]},
        }
        if trial < 8:
            traces.append(trace)
        if feedback_example is None and action == "F":
            feedback_example = trace
    risks = exact["forecast"]
    return {
        "schema": "oph.compact_observer.seeded_run.v1", "seed": seed, "trials": trials,
                "execution": {"fresh_preparations": trials, "native_events": 2 * trials,
                      "explicit_halt_checks": trials, "exact_path_reference_matches": trials,
                      "controller_bits": 4, "trace_capacity_per_trial": 2,
                      "parameter": "two sequential native event invocations per freshly prepared trial"},
        "next_outward_current": measured(counts["event"], trials, risks["event_rate"]),
        "second_seam_F": measured(counts["action_F"], trials, "5/66"),
        "independent_donor_second_seam_F": measured(counts["donor_action_F"], trials, "5/66"),
        "brier": {
            "native_record": measured(loss["native_record"], trials, risks["brier_risk"]),
            "action_known_prior": measured(loss["action_known_prior"], trials, risks["action_known_prior_brier_risk"]),
            "within_action_shuffle": measured(loss["within_action_shuffle"], trials, risks["within_action_shuffle_brier_risk"]),
        },
        "final_port0_audit": {
            "native_record_feedback": measured(counts["native_port0"], trials, "35/264"),
            "record_disconnected": measured(counts["disconnected_port0"], trials, "1/6"),
            "independent_donor_record": measured(counts["donor_port0"], trials, "1/6"),
            "paired_difference_vs_disconnected": measured(counts["native_port0"] - counts["disconnected_port0"], trials, "-3/88"),
            "paired_difference_vs_donor": measured(counts["native_port0"] - counts["donor_port0"], trials, "-3/88"),
        },
        "first_eight_trials": traces, "first_record_selected_F_example": feedback_example,
        "scope": {
            "protocol": "Independent fresh two-event trials, each followed by halt; not a continuing life.",
            "control": "Binary seam decision reads first current sign; retained record contains both coin/current letters.",
            "forecast": "Supplied exact table under known preparation; no online learning.",
            "time": "Native operation count, not calibrated physical duration.",
            "readout": "Native current letters are observer inputs; full loads appear only in external verification.",
            "controls": "Common second coin; independent marginal donor for selector and action-conditioned independent donor for prediction.",
            "sampling": "Sample means are compared with exact finite-ensemble references; no confidence interval or new physics inference.",
        },
    }


def canonical(result):
    return json.dumps(result, indent=2, sort_keys=True) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20261012)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run(args.trials, args.seed)
    output = canonical(result)
    if args.write:
        with args.write.open("x") as stream:
            stream.write(output)
    if args.check and args.check.read_text() != output:
        raise SystemExit("Seeded run receipt mismatch")
    print(json.dumps({"trials": result["trials"], "oracle_matches": result["execution"]["exact_path_reference_matches"],
                      "brier": {key: value["sample_decimal"] for key, value in result["brier"].items()},
                      "receipt_sha256": hashlib.sha256(output.encode()).hexdigest()}))
