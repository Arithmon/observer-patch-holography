"""Independent basis replay and exact verification of both certificate alternatives."""
import argparse
from fractions import Fraction as F
import hashlib
import re
from itertools import product
from pathlib import Path

from . import codec, check_tomography


def require(condition, label):
    if not condition:
        raise ValueError(label)


def rational(value):
    require(type(value) is str and len(value) <= 300, "rational encoding")
    require(re.fullmatch(r"-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?", value) is not None,
            "rational syntax")
    result = F(value)
    require(str(result) == value, "noncanonical rational")
    return result


def check_certificate(rows, target, cert):
    """A primal decoder or a normalized null witness exhausts both alternatives."""
    require(type(target) in (list, tuple) and all(type(v) in (int, F) for v in target),
            "exact target vector")
    require(type(rows) in (list, tuple) and all(
        type(row) in (list, tuple) and len(row) == len(target)
        and all(type(v) in (int, F) for v in row) for row in rows),
        "exact rectangular observation matrix")
    require(type(cert) is dict, "certificate object")
    if cert.get("kind") == "decoder":
        require(set(cert) == {"kind", "coefficients", "sample_error_gain"}, "decoder fields")
        require(type(cert["coefficients"]) is list and len(cert["coefficients"]) == len(rows),
                "decoder sample count")
        coefficients = [rational(v) for v in cert["coefficients"]]
        for j, wanted in enumerate(target):
            require(sum(c*row[j] for c, row in zip(coefficients, rows)) == wanted,
                    "decoder does not reconstruct target")
        require(rational(cert["sample_error_gain"]) == sum(abs(c) for c in coefficients),
                "undercharged sample error")
        return True
    require(cert.get("kind") == "ambiguity", "certificate alternative")
    require(set(cert) == {"kind", "perturbation"}, "ambiguity fields")
    require(type(cert["perturbation"]) is list and len(cert["perturbation"]) == len(target),
            "perturbation dimension")
    vector = [rational(v) for v in cert["perturbation"]]
    require(all(sum(a*b for a, b in zip(row, vector)) == 0 for row in rows),
            "visible ambiguity witness")
    require(sum(F(a)*b for a, b in zip(target, vector)) == 1, "trivial ambiguity witness")
    return False


def identifies_two(rows, target):
    # Independent two-column rank test, with no producer elimination code.
    pivot = next((row for row in rows if any(row)), None)
    if pivot is None:
        return not any(target)
    if any(pivot[0]*row[1] != pivot[1]*row[0] for row in rows):
        return True
    return pivot[0]*target[1] == pivot[1]*target[0]


def toy_rows():
    outputs = []
    for horizon in range(9):
        totals = [0]*4
        continuable = [0]*4
        work = [0]*4
        hit_counts = [[0]*(horizon+1) for _ in range(4)]
        stream = hashlib.sha256()
        # Scalar experiments on separate basis preparations construct each column.
        for number in range(3**horizon):
            digits = [0]*horizon
            remainder = number
            for k in reversed(range(horizon)):
                remainder, digits[k] = divmod(remainder, 3)
            columns = []
            final_columns = []
            for basis in range(2):
                state = [F(i == basis) for i in range(4)]
                samples = [state[3]]
                for edge in digits:
                    state[edge] = state[edge+1] = (state[edge]+state[edge+1])/2
                    samples.append(state[3])
                columns.append(samples)
                final_columns.append(state)
            rows = list(zip(*columns))
            remaining_state = list(zip(*final_columns))
            first = []
            for t, target in enumerate(codec.TARGETS):
                found = None
                for prefix in range(horizon+1):
                    if identifies_two(rows[:prefix+1], target):
                        found = prefix
                        break
                first.append(found)
                continuable[t] += identifies_two(rows+remaining_state, target)
                totals[t] += found is not None
                work[t] += horizon if found is None else found
                if found is not None:
                    hit_counts[t][found] += 1
            stream.update(codec.canonical([digits, [[str(x) for x in row] for row in rows], first]))
        outputs.append({"horizon": horizon, "attempt_words": 3**horizon,
                        "accepted_counts": totals, "first_counts": hit_counts,
                        "continuable_counts": continuable,
                        "irreversible_counts": [3**horizon-a for a in continuable],
                        "pending_completable_counts": [a-b for a, b in zip(continuable, totals)],
                        "abort_counts": [3**horizon-a for a in totals],
                        "stopped_mean_totals": work,
                        "stopped_mean_read_totals": [2*w for w in work],
                        "stopped_mean_write_totals": [2*w for w in work],
                        "preparation_writes": 4*3**horizon,
                        "expected_means": [str(F(w, 3**horizon)) for w in work],
                        "expected_samples": [str(F(w+3**horizon, 3**horizon)) for w in work],
                        "response_sha256": stream.hexdigest()})
    return outputs


def adjacency():
    support = codec.load(codec.ROOT/codec.SUPPORT)
    require(support["carriers"] == 1280 and support["level"] == 3, "support level")
    edges = set()
    for c in range(1280):
        for a, b in support["intra_carrier_seams"]:
            edges.add(frozenset((12*c+a, 12*c+b)))
    for c, a, d, b in support["glued_pairs"]:
        edges.add(frozenset((12*c+a, 12*d+b)))
    require(len(edges) == 46050 and all(len(e) == 2 for e in edges), "support edge census")
    return edges


def guarded_trace(word):
    # Two separate scalar basis experiments, no producer row elimination.
    columns = [[F(v == j) for v in range(4)] for j in range(2)]
    observations = [[F(0)], [F(0)]]
    decisions = []
    hit = None

    def rank_two(rows):
        return identifies_two(rows, (1, 0)) and identifies_two(rows, (0, 1))

    for step, move in enumerate(word):
        trial = [state[:] for state in columns]
        for state in trial:
            state[move] = state[move+1] = F(1, 2)*(state[move]+state[move+1])
        joint = list(zip(*observations))+list(zip(*trial))
        # The proposed receiver sample is already a row of the trial state.
        accept = rank_two(joint)
        decisions.append(accept)
        if accept:
            columns = trial
            for j in range(2):
                observations[j].append(trial[j][3])
        require(rank_two(list(zip(*observations))+list(zip(*columns))), "guard lost record")
        if rank_two(list(zip(*observations))):
            hit = step+1
            break
    return {"proposals": list(word), "admitted": decisions, "first_complete": hit,
            "means": sum(decisions), "stopped_proposals": len(decisions),
            "receiver_rows": [list(map(str, r)) for r in zip(*observations)]}


def check_guarded(packet):
    codec.equal(packet["controls"], [guarded_trace(w) for w in ((0, 0, 0, 0, 0, 0),
                (0, 1, 2, 0, 1, 2), (2, 1, 2, 0, 1, 2))], "guard control trace")
    rows = []
    for h in range(9):
        complete = means = proposals = 0
        stream = hashlib.sha256()
        for number in range(3**h):
            digits, rest = [0]*h, number
            for i in reversed(range(h)):
                rest, digits[i] = divmod(rest, 3)
            trace = guarded_trace(digits)
            complete += trace["first_complete"] is not None
            means += trace["means"]
            proposals += trace["stopped_proposals"]
            stream.update(codec.canonical(trace))
        count = 3**h
        rows.append({"horizon": h, "attempt_words": count, "complete_count": complete,
                     "pending_count": count-complete, "irreversible_count": 0,
                     "stopped_proposal_total": proposals, "executed_mean_total": means,
                     "rejected_proposal_total": proposals-means, "mean_read_total": 2*means,
                     "mean_write_total": 2*means, "preparation_write_total": 4*count,
                     "sample_total": count+means, "expected_proposals": str(F(proposals, count)),
                     "expected_means": str(F(means, count)), "trace_sha256": stream.hexdigest()})
    expected = {"scope": "maximal_joint_kernel_guard_uniform_full_proposals",
                "protected_records": [0, 1], "horizons": rows,
                "controls": [guarded_trace(w) for w in ((0, 0, 0, 0, 0, 0),
                                                       (0, 1, 2, 0, 1, 2), (2, 1, 2, 0, 1, 2))],
                "universal_completion_word": [2, 1, 2, 0, 1, 2],
                "block_length": 6, "block_words": 3**6,
                "failure_bound": "(728/729)^k after 6*k proposals",
                "expected_proposal_upper_bound": 6*3**6,
                "physical_guard_derived": False, "full_axiom_schedule_selected": False}
    codec.equal(packet, expected, "guard decisions, full attempt law, work or scope")
    return rows[-1]


def check_captured(controls):
    require(type(controls) is list and len(controls) == 3, "captured inventory")
    names = ("aggregate_only", "separated_records", "wrong_temporal_order")
    expected = (
        [[0, 1], [1, 14], [14, 23], [23, 45]],
        [[1, 14], [14, 23], [23, 45], [0, 1], [1, 14], [14, 23], [23, 45]],
        [[23, 45], [14, 23], [1, 14], [0, 1]],
    )
    graph = adjacency()
    decisions = []
    for control, name, word in zip(controls, names, expected):
        codec.equal(control.get("name"), name, "control name")
        codec.equal(control.get("word"), word, "complete temporal word")
        require(all(frozenset(edge) in graph for edge in word), "unsupported seam")
        observations = []
        for j in range(2):
            state = {p: F(p == j) for p in (0, 1, 14, 23, 45)}
            samples = [state[45]]
            for a, b in word:
                state[a] = state[b] = F(1, 2)*(state[a]+state[b])
                samples.append(state[45])
            observations.append(samples)
        rows = list(zip(*observations))
        codec.equal(control["response_rows"], [list(map(str, r)) for r in rows], "source responses")
        certs = control["prefix_certificates"]
        require(type(certs) is list and len(certs) == len(rows), "complete prefix inventory")
        accepted = []
        for k, group in enumerate(certs):
            require(type(group) is list and len(group) == 4, "complete target inventory")
            accepted.append([check_certificate(rows[:k+1], t, c) for t, c in zip(codec.TARGETS, group)])
        decisions.append(accepted)
        cases = []
        for x, y in product((-1, 0, 1), repeat=2):
            state = {p: F(3)+(x if p == 0 else y if p == 1 else 0)
                     for p in (0, 1, 14, 23, 45)}
            writer = {p: -p-1 for p in state}
            tape = []
            samples = [state[45]]
            for k, (a, b) in enumerate(word):
                mean = (state[a]+state[b])/2
                tape.append([a, b, writer[a], writer[b], str(state[a]), str(state[b]), str(mean)])
                state[a] = state[b] = mean
                writer[a] = writer[b] = k
                samples.append(state[45])
            require(all(sample-3 == row[0]*x+row[1]*y for sample, row in zip(samples, rows)),
                    "basis response does not match native execution")
            cases.append({"payload": [x, y], "samples": list(map(str, samples)),
                          "tape_sha256": codec.digest(tape)})
        expected_control = {"name": name, "word": word, "receiver": 45,
                            "preparation_ports": [0, 1], "baseline": "3",
                            "response_rows": [list(map(str, r)) for r in rows],
                            "prefix_certificates": certs, "cases": cases,
                            "means_per_case": len(word), "mean_reads_per_case": 2*len(word),
                            "mean_writes_per_case": 2*len(word), "samples_per_case": len(rows),
                            "full_support_preparation_writes": 15360}
        codec.equal(control, expected_control, "native state, sample, work or schema")
    require(decisions[0][-1] == [False, False, True, False], "aggregate boundary")
    require(decisions[1][-1] == [True]*4, "complete record recovery")
    require(decisions[2][-1] == [False]*4, "temporal-order negative control")
    return decisions


def verify(packet):
    require(type(packet) is dict, "packet object")
    codec.equal([packet.get(key) for key in ("m1_derived", "full_axiom_instantiation",
                                            "physical_precision_derived")],
                [False, False, False], "scope promotion")
    codec.equal(packet.get("pins"), codec.pins(), "source/proof custody")
    decisions = check_captured(packet["captured"])
    completion = check_tomography.check(packet["tomography"], adjacency())
    guarded = check_guarded(packet["guarded"])
    # Check loss on the entire captured state, not merely on the current receiver.
    difference = [F(0)]*15360
    difference[0], difference[1] = F(1), F(-1)
    initial_receiver = difference[45]
    difference[0] = difference[1] = (difference[0]+difference[1])/2
    require(all(v == 0 for v in difference), "global erasure identity")
    erasure = {"first_seam": [0, 1], "initial_difference": [[0, "1"], [1, "-1"]],
               "receiver_initial_difference": str(initial_receiver),
               "post_mean_difference": [[i, str(v)] for i, v in enumerate(difference) if v],
               "target_differences": [str(F(a-b)) for a, b in codec.TARGETS],
               "uniform_first_seam_probability": str(F(1, len(adjacency()))),
               "sound_individual_abort_lower_bound": str(F(1, len(adjacency()))),
               "scope": "same_continuation_no_retained_discriminating_sample"}
    expected = {"schema": "oph-temporal-acceptance-v1", "pins": codec.pins(),
                "scope": "finite_linear_observation_grammar_and_conditional_classical_selection",
                "m1_derived": False, "full_axiom_instantiation": False,
                "physical_precision_derived": False, "targets": [list(t) for t in codec.TARGETS],
                "toy": toy_rows(), "captured": packet["captured"], "erasure": erasure,
                "tomography": packet["tomography"]}
    expected["guarded"] = packet["guarded"]
    codec.equal(packet, expected, "full packet, stopped law or scope")
    return {"schema": "oph-temporal-acceptance-receipt-v1", "pins": codec.pins(),
            "controls_sha256": codec.digest(packet), "verified": True,
            "toy_word_count": sum(r["attempt_words"] for r in packet["toy"]),
            "captured_histories": 27, "captured_means": 135,
            "native_completion": completion,
            "guarded_horizon_eight": guarded,
            "prefix_target_certificates": sum(len(c["prefix_certificates"])*4 for c in packet["captured"]),
            "captured_first_reads": [[next((k for k, row in enumerate(d) if row[t]), None)
                                      for t in range(4)] for d in decisions],
            "m1_derived": False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--controls", type=Path, default=codec.HERE/"controls.json")
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    require(args.controls.stat().st_size < 2_000_000, "oversized evidence packet")
    packet = codec.load(args.controls)
    require(args.controls.read_bytes() == codec.canonical(packet), "noncanonical packet")
    receipt = verify(packet)
    if args.write_receipt:
        (codec.HERE/"receipt.json").write_bytes(codec.canonical(receipt))
    else:
        codec.equal(receipt, codec.load(codec.HERE/"receipt.json"), "committed receipt")
    print("VERIFIED", codec.digest(packet))


if __name__ == "__main__":
    main()
