"""Directed exact enclosures for a prepared logical-qubit observer channel.

Only degree-one/two exclusion observables and one fixed signed-difference
operator are needed. No full-screen density matrix or independent-particle
approximation is formed. See README.md for the channel and error proof.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import gcd
from typing import NamedTuple

import numpy as np
from scipy import sparse

from ..observer_time_fields import dual


class Output(NamedTuple):
    ports: tuple
    edges: tuple
    incident: tuple
    law: str
    cap: int
    bits: int
    rho: F
    feedback: bool
    records: tuple
    columns: tuple
    coherence: tuple
    arithmetic_error: F
    timeout: F

    @property
    def diagonal_bounds(self):
        return dict(zip(self.records, self.columns, strict=True))

    @property
    def coherences(self):
        return dict(zip(self.records, self.coherence, strict=True))


def records():
    return (tuple(('complete', w) for w in product(range(4), repeat=2))
            + tuple(('timeout', (m,)) for m in range(4))
            + (('timeout', ()),))


def category(edge):
    return int(0 in edge)+2*int(1 in edge)


def _fraction(value):
    try:
        return (type(value) is F and type(value.numerator) is int
                and type(value.denominator) is int and value.denominator > 0
                and gcd(value.numerator, value.denominator) == 1)
    except AttributeError:
        return False


def _parameters(ports, edges, incident, law, cap, bits, rho, feedback):
    dual.build(ports, edges, incident)
    if not {0, 1} <= set(ports) or any(
            set(e).intersection((0, 1)) and e not in incident for e in edges):
        raise ValueError('both logical ports and all their clock seams are required')
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('native law C or E required')
    if (type(cap) is not int or not 1 <= cap <= 1024
            or type(bits) is not int or not 8 <= bits <= 52
            or 2*len(edges)*(1 << bits) >= 1 << 63):
        raise ValueError('bounded original cap/grid with safe integer products required')
    if not _fraction(rho) or not 0 < rho < 1 or type(feedback) is not bool:
        raise ValueError('exact interior density and Boolean feedback required')


@lru_cache(maxsize=12)
def _maps(ports, edges, incident, degree):
    """Reuse the existing exact exclusion lift; no new transition formula."""
    op = dual.build(ports, edges, incident, degree)
    # Repository-internal sparse kernels are read-only, never exposed by this API.
    remote = dual._matrices(op)[0]
    marked = []
    for mask in range(4):
        menu = tuple(e for e in incident if category(e) == mask)
        marked.append(dual._matrices(dual.build(ports, edges, menu, degree))[1]
                      if menu else sparse.csr_matrix(remote.shape, dtype=np.int64))
    return remote, tuple(marked)


def _step(vector, remote, marked, denominator, cap):
    """A marked backward step sum_{j<h} A^j B_m, all divisions downward."""
    waiting = (marked@vector)//denominator
    completed = waiting.copy()
    for _ in range(1, cap):
        waiting = (remote@waiting)//denominator
        completed += waiting
    if np.any(completed < 0) or int(completed.sum()) > int(vector.sum()):
        raise ArithmeticError('the marked pullback exceeded its integer mass')
    return completed


def _interval(vector, labels, initial, rho, grid, mass):
    """Bernoulli expectation; the missing positive coefficient mass is charged."""
    sums = {}
    for label, weight in zip(labels, vector, strict=True):
        mask = int(0 in label)+2*int(1 in label)
        degree = len(label)-mask.bit_count()
        sums[mask, degree] = sums.get((mask, degree), 0)+int(weight)
    lo = sum((F(weight, grid)*rho**degree
              for (mask, degree), weight in sums.items()
              if mask & initial == mask), F(0))
    deficit = mass-F(sum(map(int, vector)), grid)
    if deficit < 0 or lo < 0 or lo > mass:
        raise ArithmeticError('directed coefficients violated their exact mass bound')
    return lo, min(mass, lo+deficit)


def _population(mass, first, second, joint):
    a, b, c = first, second, joint
    raw = ((mass-a[1]-b[1]+c[0], mass-a[0]-b[0]+c[1]),
           (a[0]-c[1], a[1]-c[0]), (b[0]-c[1], b[1]-c[0]), c)
    answer = tuple((max(F(0), lo), min(mass, hi)) for lo, hi in raw)
    if any(lo > hi for lo, hi in answer):
        raise ArithmeticError('inconsistent complete population enclosure')
    return answer


@lru_cache(maxsize=8)
def _population_data(ports, edges, incident, cap, bits, density_ratio):
    # Builtin integer keys cannot be mutated through a returned Fraction.
    rho = F(*density_ratio)
    grid, denominator = 1 << bits, 2*len(edges)
    q = 1-(1-F(len(incident), len(edges)))**cap
    weights = tuple(q*F(sum(category(e) == m for e in incident), len(incident))
                    for m in range(4))
    words = ((), *tuple((m,) for m in range(4)), *tuple(product(range(4), repeat=2)))
    mass = {word: product_fraction(weights[m] for m in word) for word in words}
    moments = {word: [[], []] for word in words}
    for degree, label in ((1, (0,)), (1, (1,)), (2, (0, 1))):
        op = dual.build(ports, edges, incident, degree)
        labels = dual.states(op)
        remote, marked = _maps(ports, edges, incident, degree)
        initial = np.zeros(len(labels), dtype=np.int64)
        initial[labels.index(label)] = grid
        coefficients = {(): initial}
        for word in words[1:]:
            # The last read acts first backwards; prepend one earlier read.
            coefficients[word] = _step(coefficients[word[1:]], remote,
                                       marked[word[0]], denominator, cap)
        for word in words:
            for r, initial_pair in enumerate((1, 2)):
                moments[word][r].append(_interval(coefficients[word], labels,
                    initial_pair, rho, grid, mass[word]))
    columns = []
    for status, word in records():
        factor = 1-q if status == 'timeout' else F(1)
        columns.append(tuple(tuple((factor*lo, factor*hi) for lo, hi in
            _population(mass[word], *moments[word][r])) for r in range(2)))
    # Midpoint linear map: half trace-norm error <= 1/4 max column width sum,
    # also after tensoring any reference. No claim that midpoints form a CP map.
    error = max(sum((hi-lo for column in columns for lo, hi in column[r]), F(0))
                for r in range(2))/4
    return tuple(columns), error, q, weights


def product_fraction(values):
    result = F(1)
    for value in values:
        result *= value
    return result


def _clone(value):
    # Fraction's protected slots are writable; keep cached certificates private.
    if type(value) is F:
        return F(value.numerator, value.denominator)
    if type(value) is tuple:
        return tuple(_clone(item) for item in value)
    return value


def output(ports, edges, incident, *, law='C', cap=384, bits=52,
           rho=F(1, 3), feedback=True):
    """Complete prepared-qubit channel, with exact rational entry enclosures.

    columns[record][input logical bit][output pair mask] is a population
    interval. coherence[record] multiplies |10><01| (and its adjoint).
    The initial environment is the declared independent Bernoulli product;
    the logical qubit can have any untouched reference. No input is disclosed.
    """
    _parameters(ports, edges, incident, law, cap, bits, rho, feedback)
    rho = F(rho.numerator, rho.denominator)
    columns, error, q, weights = _population_data(ports, edges, incident, cap, bits,
                                                (rho.numerator, rho.denominator))
    chi = (F(1), F(1, 2), F(1, 2), F(0)) if law == 'C' else (F(1), F(0), F(0), F(0))
    coherence = []
    for status, word in records():
        value = product_fraction(weights[m]*chi[m] for m in word)
        if status == 'timeout':
            value *= 1-q
        if feedback and word and word[0] == 1:
            value = -value
        coherence.append(value)
    return Output(ports, edges, incident, law, cap, bits, rho, feedback, records(),
                  _clone(columns), _clone(tuple(coherence)), _clone(error), 1-q*q)


def validate(result):
    """Recompute from original source parameters; a digest is not validation."""
    if type(result) is not Output:
        raise ValueError('a typed complete observer output required')
    expected = output(result.ports, result.edges, result.incident, law=result.law,
        cap=result.cap, bits=result.bits, rho=result.rho, feedback=result.feedback)
    if (type(result.records) is not tuple or type(result.columns) is not tuple
            or type(result.coherence) is not tuple
            or not _fraction(result.arithmetic_error) or not _fraction(result.timeout)
            or any(not _fraction(x) for x in result.coherence)
            or any(type(record) is not tuple or len(record) != 2
                or type(record[0]) is not str or type(record[1]) is not tuple
                or any(type(m) is not int for m in record[1]) for record in result.records)
            or any(type(c) is not tuple or len(c) != 2 or any(
                type(row) is not tuple or len(row) != 4 or any(
                    type(pair) is not tuple or len(pair) != 2
                    or any(not _fraction(x) for x in pair) for pair in row) for row in c)
                for c in result.columns) or result != expected):
        raise ValueError('complete output differs from the source calculation')
    return expected


def bell_blocks(result):
    """Numerical midpoints of the normalized Choi probe, for diagnostics only."""
    result = validate(result)
    blocks = {}
    for record, columns, coherence in zip(result.records, result.columns, result.coherence, strict=True):
        matrix = np.zeros((8, 8), dtype=np.complex128)
        for r in range(2):
            for state, (lo, hi) in enumerate(columns[r]):
                matrix[2*state+r, 2*state+r] = float((lo+hi)/4)
        matrix[2, 5] = matrix[5, 2] = float(coherence/2)
        blocks[record] = matrix
    return blocks
