"""Exact deterministic-scan preparation on the existing twelve-port carrier.

One sweep visits every native seam once, in an immutable, record-independent
order.  There is no random seam scheduler.  The placement coin and the C/E
Kraus records are those of the existing native primitive.  Warm-up records
and any law-selecting controller are discarded at the comparison boundary.

The certificate compares the joint output, for every untouched reference,
with I_binary/66 tensor the ORIGINAL reference marginal.  It does not assert
that a retained warm-up tape decouples, select a conserved particle sector,
or generate a nontrivial thermal/modular reference.  ``full=True`` admits the
whole 78-dimensional total-two carrier; one covering sweep removes collision
support exactly, after which the binary mixing certificate applies.

No eigenvalue theorem for random scans is used.  The diagonal kernel is
multiplied in the actual scan order using integers.  Its centered row and
column sums bound its singular norm, including for nonreversible scans.
The nonzero signed-difference blocks are checked seam by seam: C contracts
by at most 2**-10 per sweep; E kills them.  For every admitted current-record
law, positivity of the current Gram blocks gives the C0 envelope, whose
bound is 2**-5.  This last certificate also permits classical law-only
feedback from earlier coarse currents; the scan itself stays fixed and no
other operation acts on the source during preparation.  Such a controller
starts in a fixed state independent of the input source and reference.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

import numpy as np

from model.observer import native, quantum_mixing


FULL_STATES = native.states()
BINARY_STATES = native.states(binary=True)
SEAMS = quantum_mixing.seams()
DEFAULT_ORDER = tuple(range(len(SEAMS)))
MAX_SWEEPS = 64
D = len(BINARY_STATES)


@dataclass(frozen=True)
class ClassicalPower:
    """Output-row/input-column integer matrix divided by denominator."""

    counts: tuple[tuple[int, ...], ...]
    denominator: int


def _count(value):
    if type(value) is not int or not 0 <= value <= MAX_SWEEPS:
        raise ValueError(f"an integer sweep count in [0,{MAX_SWEEPS}] required")
    return value


def _full(value):
    if type(value) is not bool:
        raise ValueError("full carrier switch must be Boolean")
    return value


def _law(value, *, execution=False):
    if type(value) is not str or value not in (("C", "E") if execution else ("C", "E", "all")):
        raise ValueError("native law C/E required" if execution else "law C/E/all required")
    return value


def _order(value):
    if value is None:
        return DEFAULT_ORDER
    if (type(value) not in (tuple, list) or len(value) != len(SEAMS)
            or any(type(i) is not int for i in value)
            or set(value) != set(DEFAULT_ORDER)):
        raise ValueError("a permutation of all thirty seam indices required")
    return tuple(value)


@lru_cache(maxsize=2)
def _images(full):
    basis = FULL_STATES if full else BINARY_STATES
    positions = {state: i for i, state in enumerate(basis)}
    events = []
    for a, b in SEAMS:
        coins = []
        for coin in (0, 1):
            image = []
            for state in basis:
                output, record = native.native(state, (a, b), coin)
                if record != (coin, state[a]-output[a]) or output not in positions:
                    raise ArithmeticError("native repair left its carrier/record contract")
                image.append((positions[output], record[1], state[a]-state[b]))
            coins.append(tuple(image))
        events.append(tuple(coins))
    return tuple(events)


@lru_cache(maxsize=128)
def _classical_power(sweeps, order):
    if not sweeps:
        return ClassicalPower(tuple(tuple(int(i == j) for j in range(D)) for i in range(D)), 1)
    previous = _classical_power(sweeps-1, order)
    counts = np.array(previous.counts, dtype=object)
    denominator = previous.denominator
    for edge in order:
        updated = np.zeros((D, D), dtype=object)
        for coin in _images(False)[edge]:
            for source, (destination, _, _) in enumerate(coin):
                updated[destination, :] += counts[source, :]
        counts = updated
        denominator *= 2
    result = tuple(tuple(int(x) for x in row) for row in counts)
    if (any(sum(row) != denominator for row in result)
            or any(sum(result[i][j] for i in range(D)) != denominator for j in range(D))):
        raise ArithmeticError("the native binary sweep is not doubly stochastic")
    return ClassicalPower(result, denominator)


def classical_power(sweeps, order=None):
    """Exact, immutable diagonal kernel after the specified complete sweeps."""
    return _classical_power(_count(sweeps), _order(order))


def _sqrt_upper(value, bits=80):
    """A rational upward square root, verified by its defining inequality."""
    scale = 1 << bits
    root = isqrt(value.numerator * scale * scale // value.denominator)
    if root*root*value.denominator < value.numerator*scale*scale:
        root += 1
    return F(root, scale)


@lru_cache(maxsize=128)
def _classical_bounds(sweeps, order):
    power = _classical_power(sweeps, order)
    centered = tuple(tuple(D*x-power.denominator for x in row) for row in power.counts)
    denominator = D*power.denominator
    row = F(max(sum(abs(x) for x in line) for line in centered), denominator)
    column = F(max(sum(abs(centered[i][j]) for i in range(D)) for j in range(D)), denominator)
    # Every local heat bath is a contraction, independently of this norm bound.
    return row, column, min(F(1), _sqrt_upper(row*column))


def _swap(state, edge):
    a, b = edge
    result = list(state)
    result[a], result[b] = result[b], result[a]
    return tuple(result)


@lru_cache(maxsize=1)
def _structural_audit():
    if len(FULL_STATES) != 78 or D != 66 or len(SEAMS) != 30:
        raise ArithmeticError("unexpected native total-two carrier")
    if len(set(SEAMS)) != 30 or any(a == b or not 0 <= a < b < 12 for a, b in SEAMS):
        raise ArithmeticError("invalid native seam set")
    degree = Counter(v for edge in SEAMS for v in edge)
    if any(degree[v] != 5 for v in range(12)):
        raise ArithmeticError("the declared carrier is not five-regular")
    reached = {0}
    for _ in range(12):
        reached |= {b for a, b in SEAMS if a in reached} | {a for a, b in SEAMS if b in reached}
    if len(reached) != 12:
        raise ArithmeticError("the scan carrier is disconnected")
    index = {state: i for i, state in enumerate(BINARY_STATES)}
    differences = set()
    local_units = 0
    for edge_index, edge in enumerate(SEAMS):
        a, b = edge
        events = _images(False)[edge_index]
        # The same current has one repaired image.  Forgetting the coin
        # gives the exact coarse-current probabilities used by C0.
        coarse = []
        for i in range(D):
            outcomes = Counter((event[i][1], event[i][0]) for event in events)
            table = {}
            for (current, output), count in outcomes.items():
                if current in table:
                    raise ArithmeticError("one current has multiple repaired images")
                table[current] = (output, F(count, 2))
            coarse.append(table)
        for i, left in enumerate(BINARY_STATES):
            for j, right in enumerate(BINARY_STATES):
                z = tuple(x-y for x, y in zip(left, right, strict=True))
                differences.add(z)
                delta = abs(z[a]-z[b])
                if delta == 0:
                    common = Counter({(i, j): 1})
                    swapped = (index[_swap(left, edge)], index[_swap(right, edge)]) if z[a] == z[b] == 0 else (i, j)
                    common[swapped] += 1
                else:
                    common = Counter()
                for law, record_index in (("C", 1), ("E", 2)):
                    actual = Counter()
                    for event in events:
                        if event[i][record_index] == event[j][record_index]:
                            actual[event[i][0], event[j][0]] += 1
                    expected = common if delta == 0 else Counter({(i, j): 1}) if law == "C" and delta == 1 else Counter()
                    if actual != expected:
                        raise ArithmeticError("native signed-difference reduction failed")
                # All coefficients are nonnegative.  At delta=0 the exact
                # coefficient is rational; at delta=1/2 its square is 1/2 or
                # 1/4.  Compare squares so this audit uses no square roots.
                squared = {}
                for current in coarse[i].keys() & coarse[j].keys():
                    out_i, p_i = coarse[i][current]
                    out_j, p_j = coarse[j][current]
                    pair = out_i, out_j
                    if pair in squared:
                        raise ArithmeticError("unexpected duplicate C0 image")
                    squared[pair] = p_i*p_j
                expected_squared = ({pair: F(count, 2)**2 for pair, count in common.items()}
                                    if delta == 0 else {(i, j): F(1, 2**delta)})
                if squared != expected_squared:
                    raise ArithmeticError("native C0 squared-coefficient reduction failed")
                local_units += 1
    nonzero = tuple(z for z in differences if any(z))
    cuts = tuple((sum(abs(z[a]-z[b]) == 1 for a, b in SEAMS),
                  sum(abs(z[a]-z[b]) == 2 for a, b in SEAMS)) for z in nonzero)
    c_surviving_minimum = min(c1 for c1, c2 in cuts if not c2)
    envelope_minimum = min(c1+2*c2 for c1, c2 in cuts)
    if c_surviving_minimum != 10 or envelope_minimum != 10 or any(c1+c2 == 0 for c1, c2 in cuts):
        raise ArithmeticError("the full-sweep off-diagonal cut bound failed")
    # The first incident repair takes a collision (2,0) to (1,1), and all
    # subsequent binary states stay binary.  This covers every seam order.
    for edge_index, (a, b) in enumerate(SEAMS):
        for event in _images(True)[edge_index]:
            for i, state in enumerate(FULL_STATES):
                output = FULL_STATES[event[i][0]]
                if max(state) <= 1 and max(output) > 1:
                    raise ArithmeticError("native repair recreated a collision")
                if state[a] == 2 or state[b] == 2:
                    if max(output) != 1:
                        raise ArithmeticError("an incident repair did not split a collision")
    return (local_units, len(nonzero), c_surviving_minimum, envelope_minimum)


def structural_audit():
    """Exact local identities and graph inequalities supporting every bound."""
    units, sectors, c_cut, all_cut = _structural_audit()
    return {"binary_dimension": D, "full_dimension": len(FULL_STATES),
            "seams_per_sweep": len(SEAMS), "checked_local_matrix_units": units,
            "nonzero_difference_sectors": sectors, "C_surviving_cut_minimum": c_cut,
            "C0_weighted_cut_minimum": all_cut, "C_offdiagonal_per_sweep": F(1, 1024),
            "E_offdiagonal_per_sweep": F(0), "all_offdiagonal_per_sweep": F(1, 32),
            "collision_absorption_sweeps": 1}


def certificate(sweeps, law="C", *, full=False, order=None):
    """Reference-uniform half trace distance, entirely rational.

    Write Q=P**m-J.  The diagonal contribution is at most ||Q||_1/2:
    its input reference blocks are positive and their traces sum to one.
    On off-diagonals the signed-difference block contraction is c**m.
    For pure inputs the reference Schmidt rank is at most 66, hence trace
    norm <=66 times Hilbert--Schmidt norm; convexity handles all mixed inputs
    and arbitrary reference dimensions.  This gives ||Q||_1/2+33*c**m.
    The independent full traceless-operator bound is 33*max(||Q||_2,c**m).
    """
    sweeps, law, full, order = _count(sweeps), _law(law), _full(full), _order(order)
    structural_audit()
    mixing_sweeps = max(0, sweeps-int(full))
    row, column, diagonal_hs = _classical_bounds(mixing_sweeps, order)
    per_sweep = {"C": F(1, 1024), "E": F(0), "all": F(1, 32)}[law]
    offdiagonal = per_sweep**mixing_sweeps
    hs = max(diagonal_hs, offdiagonal)
    trace_bound = min(F(1), F(D, 2)*hs, column/2+F(D, 2)*offdiagonal)
    return {"law": law, "full_carrier": full, "sweeps": sweeps,
            "binary_mixing_sweeps": mixing_sweeps,
            "native_attempts": len(SEAMS)*sweeps, "seam_order": order,
            "collision_absorption_sweeps": int(full and sweeps > 0),
            "classical_row_l1_upper": row, "classical_column_l1_upper": column,
            "diagonal_hilbert_schmidt_upper": diagonal_hs,
            "offdiagonal_hilbert_schmidt_upper": offdiagonal,
            "hilbert_schmidt_upper_on_binary_traceless": hs,
            "trace_distance_upper": trace_bound}


def budget(tolerance, law="C", *, full=False, order=None):
    """First certified sweep count, not an optimal physical mixing time.

    Only exact integer/Fraction tolerances are accepted.  At most 64 sweeps
    are searched; failure to reach a requested accuracy is explicit.
    """
    if type(tolerance) not in (int, F) or not 0 < tolerance <= 1:
        raise ValueError("an exact rational tolerance in (0,1] required")
    law, full, order = _law(law), _full(full), _order(order)
    for sweeps in range(MAX_SWEEPS+1):
        report = certificate(sweeps, law, full=full, order=order)
        if report["trace_distance_upper"] <= tolerance:
            return dict(report, requested_tolerance=F(tolerance))
    raise ValueError("requested tolerance not certified within the 64-sweep search")


@lru_cache(maxsize=120)
def _kraus_groups(law, full, edge):
    record_index = 1 if law == "C" else 2
    groups = []
    for event in _images(full)[edge]:
        grouped = defaultdict(list)
        for source, image in enumerate(event):
            grouped[image[record_index]].append((source, image[0]))
        for pairs in grouped.values():
            if len({destination for _, destination in pairs}) != len(pairs):
                raise ArithmeticError("a native hidden-record Kraus map is not injective")
            groups.append(tuple(pairs))
    return tuple(groups)


def execute(rho, law="C", sweeps=1, *, reference_dimension=1, full=False, order=None):
    """Numerically execute actual per-seam C/E Kraus maps on source tensor R.

    The input may be subnormalized (including zero); its mass is preserved,
    never renormalized.  For such an input the certificate's state bound is
    multiplied by its trace.  The source index is major.  The full carrier
    remains 78-dimensional,
    including the zero collision coordinates after a covering sweep.  The
    numerical execution is a diagnostic; exact certificate values do not
    depend on floating execution, eigenvalues or rounding tolerances.
    """
    sweeps, law, full, order = _count(sweeps), _law(law, execution=True), _full(full), _order(order)
    if type(reference_dimension) is not int or reference_dimension < 1:
        raise ValueError("a positive integer reference dimension required")
    dimension = (len(FULL_STATES) if full else D)*reference_dimension
    if type(rho) is not np.ndarray or rho.dtype.kind not in "fc":
        raise ValueError("a plain floating or complex ndarray required")
    value = np.asarray(rho, dtype=complex)
    if value.shape != (dimension, dimension) or not np.isfinite(value).all():
        raise ValueError("a finite density matrix on the declared source/reference required")
    if (np.max(abs(value-value.conj().T)) > 1e-10
            or abs(np.trace(value).imag) > 1e-10
            or not -1e-10 <= np.trace(value).real <= 1+1e-10
            or np.linalg.eigvalsh(value).min() < -1e-10):
        raise ValueError("a Hermitian, positive input with trace at most one required")
    value = value.copy()
    reference = np.arange(reference_dimension)
    for _ in range(sweeps):
        for edge in order:
            updated = np.zeros_like(value)
            for pairs in _kraus_groups(law, full, edge):
                source = (np.array([i for i, _ in pairs])[:, None]*reference_dimension+reference).ravel()
                destination = (np.array([j for _, j in pairs])[:, None]*reference_dimension+reference).ravel()
                updated[np.ix_(destination, destination)] += value[np.ix_(source, source)]/2
            value = updated
    if not np.isfinite(value).all():
        raise ArithmeticError("nonfinite numerical native execution")
    return value
