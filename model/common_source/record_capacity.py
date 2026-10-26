"""Recovery bounds for the complete, uniformly addressed native instrument.

The bound concerns an arbitrary k-dimensional subspace of the binary total-two
carrier and arbitrary record-conditioned CPTP recovery. It uses entanglement
fidelity on the maximally mixed logical input, not a selected successful branch.
The exact source branch ranks and their complete probabilities supply the bound.
"""
from __future__ import annotations

from fractions import Fraction as F

from ..observer import native, record_recovery


def _law(value):
    if type(value) is not str or value not in ("C", "E"):
        raise ValueError("native record law C or E required")
    return value


def branch_partition(law):
    """Return effective binary branches as (edge, coin, label, input indices).

    Probability 1/(2*number of edges) multiplies each input projector. Empty
    binary labels remain included; raw padding/abort labels have zero binary
    projector and do not change this bound. This independent endpoint enumeration
    does not use the protected code or a target recovery map.
    """
    law = _law(law)
    states = native.states(binary=True)
    rows = []
    for edge in record_recovery.seams():
        a, b = edge
        for coin in (0, 1):
            groups = {label: [] for label in (-1, 0, 1)}
            for index, state in enumerate(states):
                end, (_, current) = native.native(state, edge, coin)
                if (sum(end) != 2 or any(value not in (0, 1) for value in end)
                        or end[a] != (state[a]+state[b]+coin)//2):
                    raise ArithmeticError("native binary endpoint rule changed")
                label = current if law == "C" else state[a]-state[b]
                groups[label].append(index)
            rows.extend((edge, coin, label, tuple(indices))
                        for label, indices in groups.items())
    return tuple(rows)


def recovery_bound(law, dimension):
    """Exact all-code entanglement-fidelity ceiling on this fixed instrument.

    For a branch K=sqrt(w)UP and code isometry V, set M=V*PV. Even the
    optimal decoder contributes at most w*(Tr sqrt(M))^2/k^2. Cauchy--Schwarz
    bounds this by w*min(k,rank(P))*Tr(M)/k^2. Summing the deficit gives a
    diagonal source operator A. The variational trace bound gives
    Tr(V*AV)/k >= the average of its k smallest diagonal entries. Therefore
    the resulting ceiling holds uniformly over V and every CPTP decoder. The same
    deficit lower-bounds half trace distance on a maximally entangled input.
    It is a lower bound on error, not an attainable-error claim.
    """
    law = _law(law)
    states = native.states(binary=True)
    if type(dimension) is not int or not 1 <= dimension <= len(states):
        raise ValueError("integer logical dimension from 1 through 66 required")
    rows = branch_partition(law)
    edge_count = len(record_recovery.seams())
    weight = F(1, 2*edge_count)
    total = [F(0) for _ in states]
    deficit = [F(0) for _ in states]
    ranks = {}
    for _, _, _, indices in rows:
        rank = len(indices)
        ranks[rank] = ranks.get(rank, 0)+1
        loss = max(F(0), F(dimension-rank, dimension))
        for index in indices:
            total[index] += weight
            deficit[index] += weight*loss
    if any(value != 1 for value in total):
        raise ArithmeticError("incomplete native record partition")
    lower = sum(sorted(deficit)[:dimension], F(0))/dimension
    return {
        "law": law,
        "logical_dimension": dimension,
        "carrier_dimension": len(states),
        "uniform_seams": edge_count,
        "branch_rank_counts": {str(rank): ranks[rank] for rank in sorted(ranks)},
        "entanglement_fidelity_upper": str(1-lower),
        "reference_trace_distance_lower": str(lower),
        "maximum_basis_deficit": str(max(deficit)),
        "scope": "one uniformly selected seam; complete records; arbitrary subspace and CPTP recovery",
    }


def certificate():
    return {law: recovery_bound(law, 15) for law in ("C", "E")}
