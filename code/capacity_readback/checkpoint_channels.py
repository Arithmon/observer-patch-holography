"""Exact finite checkpoint kernels and deterministic worst-input decoding.

Near-unit input row masses retain the historical 1e-12 acceptance boundary.
Their relative weights are normalized exactly, explicitly disclosed in each
decoder certificate. No support entry or decision threshold is rounded.
"""
from decimal import Decimal
from fractions import Fraction
import math
from typing import Mapping, Sequence


ROW_MASS_TOLERANCE = Fraction(1, 10**12)


def probability(value, name="probability"):
    if isinstance(value, bool) or not isinstance(value, (int, float, Fraction, Decimal)):
        raise ValueError(f"{name} must be a finite real number, not a Boolean or string")
    try:
        result = Fraction(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not 0 <= result <= 1:
        raise ValueError(f"{name} must lie in [0,1]")
    return result


def alphabet(values, name="record alphabet"):
    if (not isinstance(values, Sequence) or isinstance(values, (str, bytes))
            or not values or any(type(x) is not str or not x for x in values)
            or len(set(values)) != len(values)):
        raise ValueError(f"{name} must contain distinct nonempty string labels")
    return tuple(sorted(values))


def channel_rows(channel, reachable):
    """Return exact stochastic rows and the supplied pre-normalization masses."""
    records = alphabet(reachable)
    if not isinstance(channel, Mapping):
        raise ValueError("checkpoint channel must be a mapping")
    rows = channel.get("rows")
    if not isinstance(rows, Mapping) or set(rows) != set(records):
        raise ValueError("every checkpoint kernel needs one row per reachable record")
    normalized, masses = {}, {}
    for source in records:
        row = rows[source]
        if (not isinstance(row, Mapping) or not row
                or any(type(out) is not str or not out for out in row)):
            raise ValueError("checkpoint rows require nonempty string output labels")
        exact = {out: probability(p) for out, p in row.items()}
        mass = sum(exact.values(), Fraction(0))
        if mass == 0 or abs(mass-1) > ROW_MASS_TOLERANCE:
            raise ValueError("checkpoint rows must have unit mass within 1e-12")
        normalized[source] = {out: p/mass for out, p in exact.items()}
        masses[source] = mass
    return normalized, masses


def channel_family(channels, reachable):
    if (not isinstance(channels, Sequence) or isinstance(channels, (str, bytes))
            or not channels):
        raise ValueError("a nonempty complete checkpoint channel family is required")
    return [channel_rows(channel, reachable) for channel in channels]


def ratio(value):
    value = Fraction(value)
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def directed_float(value, *, upward):
    """Outward display bound; exact ratios remain authoritative."""
    value = Fraction(value)
    result = float(value)
    if (upward and Fraction(result) < value) or (not upward and Fraction(result) > value):
        result = math.nextafter(result, math.inf if upward else -math.inf)
    return result


class SearchBudget:
    def __init__(self, limit):
        if type(limit) is not int or limit < 1:
            raise ValueError("max_decoder_nodes must be a positive integer")
        self.remaining = limit

    def charge(self):
        self.remaining -= 1
        if self.remaining < 0:
            raise ValueError("exact decoder search budget exhausted; no capacity certified")


def optimal_decoder(rows, code, budget):
    """Exact minimax decoder after eliminating dominated output assignments.

    An output possible for only one codeword is assigned to that word. For a
    shared output, assignments to words having zero probability are dominated.
    The current partial errors are lower bounds for every completion. An
    explicit stack avoids dependence on Python's recursion depth.
    """
    outputs = sorted({out for source in code for out in rows[source]})
    fixed, ambiguous = {}, []
    for out in outputs:
        owners = tuple(i for i, source in enumerate(code) if rows[source].get(out, 0) > 0)
        if len(owners) <= 1:
            fixed[out] = code[owners[0]] if owners else code[0]
        else:
            ambiguous.append((out, owners, tuple(rows[source].get(out, Fraction(0)) for source in code)))
    best_error, best_decoder = None, None
    # Keep one path of owner indices, not a copied decoder dictionary at
    # every pending branch. Budget allocations as well as visited leaves.
    choices = [None]*len(ambiguous)
    budget.charge()
    stack = [(0, (Fraction(0),)*len(code), None)]
    while stack:
        index, errors, owner = stack.pop()
        if index:
            choices[index-1] = owner
        lower = max(errors)
        if best_error is not None and lower >= best_error:
            continue
        if index == len(ambiguous):
            best_error = lower
            best_decoder = {**fixed, **{item[0]: code[chosen]
                                        for item, chosen in zip(ambiguous, choices)}}
            continue
        out, owners, masses = ambiguous[index]
        for owner in reversed(owners):
            updated = tuple(error+(p if i != owner else 0)
                            for i, (error, p) in enumerate(zip(errors, masses)))
            budget.charge()
            stack.append((index+1, updated, owner))
    return best_error, best_decoder


def decoder_certificate(rows, masses, code, decoder):
    errors = {source: sum((p for out, p in rows[source].items()
                           if decoder[out] != source), Fraction(0)) for source in code}
    return {"decoder": decoder, "error_by_input": {x: ratio(e) for x, e in errors.items()},
            "worst_input_error": ratio(max(errors.values())),
            "input_row_masses": {x: ratio(masses[x]) for x in code}}
