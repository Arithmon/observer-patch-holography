"""Record-only reconstruction and fair stopping on the existing native carrier.

The fixed forests minimize E's predetermined event count at total load two.
C's doubled-forest and stopping protocols are constructions, not optimality
claims. All coin histories retain their raw probability, including prefixes
that have not stopped. No protected quantum register or coherent feedback is
included in these classical seam-selection protocols.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import comb

from . import record_recovery as recovery


_FULL_FOREST = ((0, 1), (0, 5), (2, 3), (2, 4), (6, 8), (8, 9), (7, 10), (10, 11))
_BINARY_FOREST = ((1, 5), (2, 3), (2, 4), (6, 8), (8, 9), (7, 10), (10, 11))


def forest(*, binary=False):
    recovery.carrier(binary)  # Validate before bool/int dictionary aliases.
    return _BINARY_FOREST if binary else _FULL_FOREST


def doubled_forest(*, binary=False):
    return tuple(edge for edge in forest(binary=binary) for _ in range(2))


def components(edges):
    """Connected components of the visited graph, including isolated ports."""
    if not isinstance(edges, (tuple, list)):
        raise ValueError("a finite seam sequence required")
    groups = [{i} for i in range(12)]
    for edge in edges:
        a, b = recovery.seam(edge)
        left = next(group for group in groups if a in group)
        right = next(group for group in groups if b in group)
        if left is not right:
            left.update(right)
            groups.remove(right)
    return tuple(sorted(tuple(sorted(group)) for group in groups))


def separates(edges, *, binary=False):
    """Whether E's records determine every initial load for this fixed word.

    Equal histories differ by a constant on each visited component. With
    total load two, full-carrier ambiguity needs two components of size at
    most two; binary ambiguity needs two singletons or two size-two components.
    """
    recovery.carrier(binary)
    sizes = [len(group) for group in components(edges)]
    return ((sizes.count(1) <= 1 and sizes.count(2) <= 1) if binary
            else sum(size <= 2 for size in sizes) <= 1)


def initial_constraints(law, history):
    """Derive original-load contrast constraints without replaying any input.

    A C record gives d_before in {2*j, 2*j+(2*c-1)}; E gives d_before
    directly. Subtracting the accumulated, known current translations gives
    constraints on the original load. This also treats reversed seam labels.
    """
    law = recovery.candidate(law)
    if not isinstance(history, (tuple, list)):
        raise ValueError("a complete history sequence required")
    displacement = [0] * 12
    constraints = {}
    for record in history:
        current = recovery.record_current(law, record)
        a, b, coin, value = record
        shift = displacement[a]-displacement[b]
        values = ({value-shift} if law == "E"
                  else {2*current-shift, 2*current+(2*coin-1)-shift})
        edge = (a, b) if a < b else (b, a)
        if a > b:
            values = {-value for value in values}
        constraints[edge] = constraints.get(edge, values) & values
        if not constraints[edge]:
            raise ValueError("inconsistent initial-load constraints")
        displacement[a] -= current
        displacement[b] += current
    return tuple((a, b, tuple(sorted(values)))
                 for (a, b), values in sorted(constraints.items()))


def decode(law, history, *, binary=False):
    """All compatible initial loads, using only committed classical records."""
    constraints = initial_constraints(law, history)
    result = tuple(state for state in recovery.carrier(binary)
                   if all(state[a]-state[b] in values for a, b, values in constraints))
    if not result:
        raise ValueError("record history has no compatible initial load")
    return result


def reconstruct(law, history, *, binary=False):
    result = decode(law, history, binary=binary)
    if len(result) != 1:
        raise ValueError("record history does not uniquely identify the initial load")
    return result[0]


def next_seam(history, *, binary=False):
    """Stopped C controller: repeat an edge until its first opposite coin.

    Only the supplied record prefix is read. Completed edge blocks reveal
    their initial contrast as 2*j_first+j_last. No load register is read.
    """
    decode("C", history, binary=binary)
    edges = forest(binary=binary)
    position, first_coin = 0, None
    for a, b, coin, _ in history:
        if position == len(edges) or (a, b) != edges[position]:
            raise ValueError("history does not follow the stopped forest controller")
        if first_coin is None:
            first_coin = coin
        elif coin != first_coin:
            position += 1
            first_coin = None
    return None if position == len(edges) else edges[position]


def execute(law, initial, coins, *, binary=False):
    """Execute a supplied finite coin tape, retaining unsuccessful stop prefixes.

    E visits its forest once. C follows ``next_seam``. Extra coins after a
    completed run are returned as unused; insufficient tapes remain explicit
    incomplete outcomes. Quantum branches are never renormalized.
    """
    law, initial = recovery.candidate(law), recovery.load(initial)
    if initial not in recovery.carrier(binary):
        raise ValueError("initial load is outside the requested carrier")
    if (not isinstance(coins, (tuple, list))
            or any(type(coin) is not int or coin not in (0, 1) for coin in coins)):
        raise ValueError("a finite sequence of integer fair-coin values required")
    state, history = initial, []
    edges = forest(binary=binary)
    for coin in coins:
        edge = (next_seam(history, binary=binary) if law == "C"
                else edges[len(history)] if len(history) < len(edges) else None)
        if edge is None:
            break
        state, record = recovery.recorded_step(law, state, edge, coin)
        history.append(record)
    complete = (next_seam(history, binary=binary) is None if law == "C"
                else len(history) == len(edges))
    decoded = decode(law, history, binary=binary)
    if initial not in decoded or (complete and decoded != (initial,)):
        raise ArithmeticError("record-only reconstruction differs from native execution")
    return {"history": tuple(history), "terminal_load": state,
            "compatible_initial_loads": decoded, "complete": complete,
            "attempts": len(history), "coin_prefix_probability": F(1, 2**len(history)),
            "unused_coins": tuple(coins[len(history):])}


def stopping_mass(attempts, *, binary=False):
    """Exact P(T=t), where T=m+NegBin(m,1/2) and m is the forest size."""
    attempts = recovery.attempts(attempts)
    count = len(forest(binary=binary))
    return (F(0) if attempts < 2*count
            else F(comb(attempts-count-1, count-1), 2**(attempts-count)))


def stopping_tail(attempts, *, binary=False):
    """Exact unconditional P(T>t); unfinished runs remain in this tail."""
    attempts = recovery.attempts(attempts)
    count = len(forest(binary=binary))
    if attempts < 2*count:
        return F(1)
    trials = attempts-count
    return F(sum(comb(trials, successes) for successes in range(count)), 2**trials)


def common_c_history(policy, attempts):
    """Construct a rank-at-least-two history of a bounded record-only policy.

    The fixed binary inputs differ at two nonadjacent ports. Every chosen
    seam admits a matching coin; their signed difference therefore remains
    unchanged. The policy receives only the common classical history and may
    return None to halt. This does not cover coherent operations between events.
    """
    attempts = recovery.attempts(attempts)
    if not callable(policy):
        raise ValueError("a record-only seam policy required")
    edges = {frozenset(edge) for edge in recovery.seams()}
    u, v = next((u, v) for u in range(12) for v in range(u+1, 12)
                if frozenset((u, v)) not in edges)
    spectator = next(i for i in range(12) if i not in (u, v))
    left, right = [0]*12, [0]*12
    left[u] = left[spectator] = right[v] = right[spectator] = 1
    initial = (tuple(left), tuple(right))
    left, right = initial
    history = []
    for _ in range(attempts):
        edge = policy(tuple(history))
        if edge is None:
            break
        edge = recovery.seam(edge)
        for coin in (0, 1):
            out_left, record_left = recovery.recorded_step("C", left, edge, coin)
            out_right, record_right = recovery.recorded_step("C", right, edge, coin)
            if record_left == record_right:
                left, right = out_left, out_right
                history.append(record_left)
                break
        else:
            raise ArithmeticError("a nonadjacent difference has no matching coin")
    if left == right or any(state not in decode("C", history, binary=True) for state in initial):
        raise ArithmeticError("the bounded-policy entangled branch lost its two inputs")
    return {"initial_pair": initial, "terminal_pair": (left, right),
            "history": tuple(history), "coin_prefix_probability": F(1, 2**len(history))}


def audit():
    # Exhaust the local identity used by every stopped or doubled edge block.
    for contrast in range(-2, 3):
        for coin in (0, 1):
            current = (contrast-(2*coin-1)*(abs(contrast) % 2))//2
            after = contrast-2*current
            opposite = (after+(2*coin-1)*(abs(after) % 2))//2
            if 2*current+opposite != contrast:
                raise ArithmeticError("opposite-coin contrast reconstruction failed")
    reports = {}
    for binary in (False, True):
        edges = forest(binary=binary)
        groups = components(edges)
        count = len(edges)
        if not separates(edges, binary=binary) or count != (7 if binary else 8):
            raise ArithmeticError("the source forest misses its separation certificate")
        reports["binary" if binary else "full"] = {
            "forest": [list(edge) for edge in edges],
            "component_sizes": sorted(map(len, groups)),
            "minimum_predetermined_E_attempts": count,
            "doubled_C_attempts": 2*count,
            "rank_one_C_effect_sum": str(F(1, 2**count)),
            "C_stopping_mean_attempts": 3*count,
            "C_stopping_tail_at_mean": str(stopping_tail(3*count, binary=binary)),
            "C_stopping_tail_at_twice_mean": str(stopping_tail(6*count, binary=binary)),
            "C_stopping_time_law": "T=m+NegBin(m,1/2); P(T=t)=choose(t-m-1,m-1)/2^(t-m)",
        }
    return {"schema": "oph.native_record_schedule.v1", "carriers": reports,
            "C_bounded_record_only_policy_entanglement_breaking": False,
            "C_stopped_controller_almost_surely_finite": True,
            "coherent_controls_or_protected_registers_covered": False}
