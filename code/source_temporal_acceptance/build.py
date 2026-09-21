"""Producer: source response matrices, primal/dual certificates and stopped laws."""
from fractions import Fraction as F
import hashlib
from itertools import product

from . import codec, tomography


def certificate(rows, target):
    """Eliminate the observation rows while retaining their sample coefficients."""
    n, m = len(target), len(rows)
    basis = {}
    for i, row in enumerate(rows):
        v = list(map(F, row))
        weights = [F(i == j) for j in range(m)]
        for p, (b, w) in sorted(basis.items()):
            scale = v[p]
            v = [x-scale*y for x, y in zip(v, b)]
            weights = [x-scale*y for x, y in zip(weights, w)]
        pivot = next((j for j in range(n) if v[j]), None)
        if pivot is not None:
            scale = v[pivot]
            basis[pivot] = ([x/scale for x in v], [x/scale for x in weights])
    residual = list(map(F, target))
    decoder = [F(0)]*m
    for p, (b, w) in sorted(basis.items()):
        scale = residual[p]
        residual = [x-scale*y for x, y in zip(residual, b)]
        decoder = [x+scale*y for x, y in zip(decoder, w)]
    if not any(residual):
        return {"kind": "decoder", "coefficients": list(map(str, decoder)),
                "sample_error_gain": str(sum(map(abs, decoder)))}
    free = next(j for j, value in enumerate(residual) if value)
    witness = [F(0)]*n
    witness[free] = 1
    for p, (b, _) in sorted(basis.items(), reverse=True):
        witness[p] = -sum(b[j]*witness[j] for j in range(p+1, n))
    scale = sum(F(a)*b for a, b in zip(target, witness))
    witness = [x/scale for x in witness]
    return {"kind": "ambiguity", "perturbation": list(map(str, witness))}


def evolution(word, size, receiver):
    state = [[F(i == j) for j in range(2)] for i in range(size)]
    rows = [state[receiver][:]]
    for a, b in word:
        value = [(x+y)/2 for x, y in zip(state[a], state[b])]
        state[a] = state[b] = value
        rows.append(value[:] if receiver in (a, b) else state[receiver][:])
    return rows, state


def responses(word, size, receiver):
    return evolution(word, size, receiver)[0]


def toy():
    alphabet = ((0, 1), (1, 2), (2, 3))
    result = []
    for horizon in range(9):
        accepted = [0]*4
        continuable = [0]*4
        stopped_means = [0]*4
        first_counts = [[0]*(horizon+1) for _ in range(4)]
        stream = hashlib.sha256()
        for word in product(range(3), repeat=horizon):
            rows, state = evolution([alphabet[e] for e in word], 4, 3)
            first = []
            for t, target in enumerate(codec.TARGETS):
                hit = next((k for k in range(horizon+1)
                            if certificate(rows[:k+1], target)["kind"] == "decoder"), None)
                first.append(hit)
                continuable[t] += certificate(rows+state, target)["kind"] == "decoder"
                accepted[t] += hit is not None
                stopped_means[t] += horizon if hit is None else hit
                if hit is not None:
                    first_counts[t][hit] += 1
            stream.update(codec.canonical([word, [[str(x) for x in row] for row in rows], first]))
        result.append({"horizon": horizon, "attempt_words": 3**horizon,
                       "accepted_counts": accepted, "first_counts": first_counts,
                       "continuable_counts": continuable,
                       "irreversible_counts": [3**horizon-a for a in continuable],
                       "pending_completable_counts": [a-b for a, b in zip(continuable, accepted)],
                       "abort_counts": [3**horizon-a for a in accepted],
                       "stopped_mean_totals": stopped_means,
                       "stopped_mean_read_totals": [2*w for w in stopped_means],
                       "stopped_mean_write_totals": [2*w for w in stopped_means],
                       "preparation_writes": 4*3**horizon,
                       "expected_means": [str(F(w, 3**horizon)) for w in stopped_means],
                       "expected_samples": [str(1+F(w, 3**horizon)) for w in stopped_means],
                       "response_sha256": stream.hexdigest()})
    return result


def captured():
    # Both preparations vary ports 0 and 1; every other port has the same baseline.
    words = {
        "aggregate_only": [(0, 1), (1, 14), (14, 23), (23, 45)],
        "separated_records": [(1, 14), (14, 23), (23, 45), (0, 1),
                              (1, 14), (14, 23), (23, 45)],
        "wrong_temporal_order": [(23, 45), (14, 23), (1, 14), (0, 1)],
    }
    result = []
    for name, word in words.items():
        rows = responses(word, 46, 45)
        cases = []
        for payload in product((-1, 0, 1), repeat=2):
            values = [F(3)]*46
            values[0] += payload[0]
            values[1] += payload[1]
            samples = [values[45]]
            tape = []
            writers = [-i-1 for i in range(46)]
            for event, (a, b) in enumerate(word):
                value = (values[a]+values[b])/2
                tape.append([a, b, writers[a], writers[b], str(values[a]), str(values[b]), str(value)])
                values[a] = values[b] = value
                writers[a] = writers[b] = event
                samples.append(values[45])
            cases.append({"payload": list(payload), "samples": list(map(str, samples)),
                          "tape_sha256": codec.digest(tape)})
        result.append({"name": name, "word": [list(e) for e in word], "receiver": 45,
                       "preparation_ports": [0, 1], "baseline": "3",
                       "response_rows": [[str(x) for x in row] for row in rows],
                       "prefix_certificates": [[certificate(rows[:k+1], t) for t in codec.TARGETS]
                                               for k in range(len(rows))],
                       "cases": cases, "means_per_case": len(word),
                       "mean_reads_per_case": 2*len(word), "mean_writes_per_case": 2*len(word),
                       "samples_per_case": len(word)+1,
                       "full_support_preparation_writes": 15360})
    return result


def guarded_word(word):
    state = [[F(i == j) for j in range(2)] for i in range(4)]
    rows = [state[3][:]]
    admitted = []
    first = None
    for k, edge in enumerate(word, 1):
        proposed = [r[:] for r in state]
        proposed[edge] = proposed[edge+1] = [
            (x+y)/2 for x, y in zip(state[edge], state[edge+1])]
        trial_rows = rows+[proposed[3][:]]
        safe = all(certificate(trial_rows+proposed, t)["kind"] == "decoder"
                   for t in codec.TARGETS[:2])
        admitted.append(safe)
        if safe:
            state, rows = proposed, trial_rows
        if all(certificate(rows, t)["kind"] == "decoder" for t in codec.TARGETS[:2]):
            first = k
            break
    return {"proposals": list(word), "admitted": admitted, "first_complete": first,
            "means": sum(admitted), "stopped_proposals": len(admitted),
            "receiver_rows": [list(map(str, r)) for r in rows]}


def guarded_census():
    census = []
    for horizon in range(9):
        complete = work = proposals = 0
        stream = hashlib.sha256()
        for word in product(range(3), repeat=horizon):
            row = guarded_word(word)
            complete += row["first_complete"] is not None
            work += row["means"]
            proposals += row["stopped_proposals"]
            stream.update(codec.canonical(row))
        count = 3**horizon
        census.append({"horizon": horizon, "attempt_words": count,
                       "complete_count": complete, "pending_count": count-complete,
                       "irreversible_count": 0, "stopped_proposal_total": proposals,
                       "executed_mean_total": work, "rejected_proposal_total": proposals-work,
                       "mean_read_total": 2*work, "mean_write_total": 2*work,
                       "preparation_write_total": 4*count, "sample_total": count+work,
                       "expected_proposals": str(F(proposals, count)),
                       "expected_means": str(F(work, count)), "trace_sha256": stream.hexdigest()})
    return {"scope": "maximal_joint_kernel_guard_uniform_full_proposals",
            "protected_records": [0, 1], "horizons": census,
            "controls": [guarded_word(w) for w in ((0, 0, 0, 0, 0, 0),
                                                     (0, 1, 2, 0, 1, 2), (2, 1, 2, 0, 1, 2))],
            "universal_completion_word": [2, 1, 2, 0, 1, 2],
            "block_length": 6, "block_words": 729,
            "failure_bound": "(728/729)^k after 6*k proposals",
            "expected_proposal_upper_bound": 4374,
            "physical_guard_derived": False, "full_axiom_schedule_selected": False}


def build():
    return {"schema": "oph-temporal-acceptance-v1", "pins": codec.pins(),
            "scope": "finite_linear_observation_grammar_and_conditional_classical_selection",
            "m1_derived": False, "full_axiom_instantiation": False,
            "physical_precision_derived": False, "targets": [list(t) for t in codec.TARGETS],
            "toy": toy(), "captured": captured(), "tomography": tomography.build(),
            "guarded": guarded_census(),
            "erasure": {"first_seam": [0, 1], "initial_difference": [[0, "1"], [1, "-1"]],
                        "receiver_initial_difference": "0", "post_mean_difference": [],
                        "target_differences": ["1", "-1", "0", "2"],
                        "uniform_first_seam_probability": "1/46050",
                        "sound_individual_abort_lower_bound": "1/46050",
                        "scope": "same_continuation_no_retained_discriminating_sample"}}


if __name__ == "__main__":
    (codec.HERE/"controls.json").write_bytes(codec.canonical(build()))
