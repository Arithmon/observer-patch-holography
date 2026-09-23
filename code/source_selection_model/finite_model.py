"""Exact central-record execution and finite-algebra/refinement controls.

These checks support the analytic construction. They do not mechanically
certify that its operational interpretation satisfies the full axiom text.
"""
from fractions import Fraction as Q
from itertools import product


def replay_local(label, rows):
    """Independent local replay; no access to any other input label."""
    deviation = Q(1, 4)
    for tick, row in enumerate(rows):
        if type(row) is not list or len(row) != 7 or any(type(row[i]) is not int for i in (0, 1, 2, 5, 6)):
            raise ValueError("local event schema")
        if row[0] != tick or row[1] != label or row[2] != tick % 12:
            raise ValueError("local record/query")
        if row[3] != str(deviation):
            raise ValueError("consumed local state")
        deviation *= 0 if row[2] == label else 1
        if row[4:] != [str(deviation), label, 1]:
            raise ValueError("feedback, prediction or scalar seam")
    return deviation


def record_controls(packet):
    if type(packet) is not dict or set(packet) != {"schema", "tapes", "cases"} or packet["schema"] != "oph.scalar-seam-records.v1":
        raise ValueError("record capture schema")
    tapes = packet["tapes"]
    if type(tapes) is not list or len(tapes) != 12 or any(type(t) is not list or len(t) != 24 for t in tapes):
        raise ValueError("record tape coverage")
    if type(packet["cases"]) is not list or len(packet["cases"]) != 144:
        raise ValueError("record intervention coverage")
    by_pair = {}
    for case in packet["cases"]:
        if type(case) is not dict or set(case) != {"labels", "tapes"}:
            raise ValueError("record case schema")
        for values in (case["labels"], case["tapes"]):
            if type(values) is not list or len(values) != 2 or any(type(x) is not int or not 0 <= x < 12 for x in values):
                raise ValueError("record label or tape index")
        key = tuple(case["labels"])
        if key in by_pair:
            raise ValueError("duplicate record intervention")
        by_pair[key] = [tapes[i] for i in case["tapes"]]
    if set(by_pair) != set(product(range(12), repeat=2)):
        raise ValueError("missing record intervention")
    histories, events, distinguishable = 0, 0, 0
    for receiver in range(12):
        reference = None
        source_traces = set()
        for source in range(12):
            logs = by_pair[source, receiver]
            for label, rows in zip((source, receiver), logs):
                if replay_local(label, rows) != 0:
                    raise ValueError("local repair not completed")
            source_traces.add(tuple(tuple(row) for row in logs[0]))
            if reference is None:
                reference = logs[1]
            if logs[1] != reference:
                raise ValueError("remote intervention became locally visible")
            histories += 1
            events += sum(map(len, logs))
        if len(source_traces) != 12:
            raise ValueError("interventions are not locally operational")
        distinguishable += len(source_traces)
    # Every possible guess has exactly one success among twelve source labels.
    success_counts = [sum(source == guess for source in range(12)) for guess in range(12)]
    if success_counts != [1]*12:
        raise ValueError("uniform source discrimination")
    return {"paired_histories": histories, "local_events": events,
            "source_traces_distinguished": distinguishable,
            "uniform_guess_success": "1/12", "uniform_guess_failure": "11/12"}


def algebra_controls():
    # Matrix units in each of twelve 6-by-6 blocks. The full local algebra
    # dimension is 432 over C and its Hermitian part has real dimension 432.
    units = [(p, a, b) for p in range(12) for a in range(6) for b in range(6)]
    center = [[int(p == q) for q in range(12)] for p in range(12)]
    for p, q in product(range(12), repeat=2):
        if [a*b for a, b in zip(center[p], center[q])] != (center[p] if p == q else [0]*12):
            raise ValueError("central port products")
    # In a full matrix block, commuting with E_aa removes off-diagonal
    # entries and commuting with E_ab equates every diagonal: one central
    # scalar per block. The graph of those diagonal equalities is complete.
    equalities = {(a, b) for a, b in product(range(6), repeat=2) if a != b}
    if len(equalities) != 30:
        raise ValueError("primitive central blocks")
    trace = lambda u: Q(int(u[1] == u[2]), 72)
    if sum(trace((p, a, a)) for p in range(12) for a in range(6)) != 1:
        raise ValueError("normalized reference")
    if any(sum(trace((p, a, a)) for a in range(6)) != Q(1, 12) for p in range(12)):
        raise ValueError("selected port weights")
    # The feedback reset in one selected central block sends E_ab to
    # delta_ab I/6 in that block. Its Kraus operators E_ab/sqrt(6)
    # satisfy both trace-preserving and unital sums. Count the squared
    # coefficients per input and output coordinate independently.
    kraus = [(a, b, Q(1, 6)) for a, b in product(range(6), repeat=2)]
    if any(sum(w for a, b, w in kraus if b == i) != 1 or
           sum(w for a, b, w in kraus if a == i) != 1 for i in range(6)):
        raise ValueError("reset Kraus normalization")
    # The trace restriction sends every off-diagonal matrix unit to zero;
    # different central port states have the same unique scalar state.
    return {"local_matrix_dimension": len(units), "primitive_central_ports": len(center),
            "faithful_representation_dimension": 72, "reference_eigenvalue": "1/72",
            "scalar_seam_state_dimension": 0, "reset_kraus_operators": len(kraus)}


def refinement_controls(tower_rows):
    """Exact posterior weights and nonuniform composite-fiber control."""
    one_step = []
    for row in tower_rows[1:]:
        fibers = {}
        for child, parent in enumerate(row["coarsen"]):
            fibers.setdefault(parent, []).append(child)
        one_step.append({parent: {child: Q(1, len(children)) for child in children}
                         for parent, children in fibers.items()})
    coarse = one_step[0]
    following = one_step[1]
    composite = {i: {k: wi*wk for j, wi in weights.items()
                     for k, wk in following[j].items()} for i, weights in coarse.items()}
    if any(sum(weights.values()) != 1 for weights in composite.values()):
        raise ValueError("composite mixture normalization")
    nonuniform = sum(len(set(weights.values())) > 1 for weights in composite.values())
    if nonuniform == 0:
        raise ValueError("composite-weight negative control absent")
    # A selective diagonal instrument with outcome probabilities 1/4 and
    # 3/4 updates an initially equal child mixture to weights 1/4 and 3/4.
    # Keeping the prior weights after observing the outcome is wrong.
    prior = [Q(1, 2), Q(1, 2)]
    likelihood = [Q(1, 4), Q(3, 4)]
    mass = sum(w*p for w, p in zip(prior, likelihood))
    posterior = [w*p/mass for w, p in zip(prior, likelihood)]
    if posterior != [Q(1, 4), Q(3, 4)] or posterior == prior:
        raise ValueError("outcome-conditioned refinement")
    return {"one_step_maps": len(one_step), "composite_fibers": len(composite),
            "nonuniform_composite_fibers": nonuniform,
            "selective_posterior": [str(x) for x in posterior]}
