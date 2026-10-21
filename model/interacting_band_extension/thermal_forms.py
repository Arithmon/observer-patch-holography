"""Finite lower forms for the fixed two-mode interacting Gibbs operator.

The existing source code uses frequency-two Hermite functions in both modes.
This independent partition bracket uses the natural harmonic frequencies
(1, 2), only to make the complement trace explicit.  It changes neither that
source code nor its occupation projection.  Rational Bargmann coefficients
and their factorial metric are exact.  Conversion to an orthonormal scipy
matrix is a floating diagnostic, not an eigenvalue or interval certificate.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial, sqrt

from scipy import sparse


def _integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('an original integer in the admitted range is required')
    return value


def _rational(value, *, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and not value):
        raise ValueError('an original exact nonnegative rational is required')
    return F(value)


def _orthonormal_matrix(entries, metric, dimension):
    rows, columns, values = [], [], []
    for row, column, value in entries:
        rows.append(row)
        columns.append(column)
        values.append(float(value)*sqrt(float(F(metric[row], metric[column]))))
    return sparse.csr_matrix((values, (rows, columns)), shape=(dimension, dimension))


@dataclass(frozen=True)
class FiniteLowerForm:
    """Exact sparse coefficients in the factorial-metric Bargmann basis.

    Index ``n1 + occupation*n2`` denotes the natural-frequency occupation
    pair.  If ``M`` is the rational matrix and ``G=diag(metric)``, then
    ``G M = M.T G`` and the orthonormal matrix is ``sqrt(G) M / sqrt(G)``.
    """
    shell: int
    coupling: F
    occupation: int
    entries: tuple[tuple[int, int, F], ...]
    metric: tuple[int, ...]
    basis_frequencies: tuple[int, int] = (1, 2)

    @property
    def dimension(self):
        return self.occupation**2

    def orthonormal_matrix(self):
        """Fresh floating matrix; its construction carries no roundoff bound."""
        return _orthonormal_matrix(self.entries, self.metric, self.dimension)


@dataclass(frozen=True)
class FixedCodeForm:
    """Full-polynomial compression in the existing frequency-two source code."""
    coupling: F
    occupation: int
    entries: tuple[tuple[int, int, F], ...]
    metric: tuple[int, ...]
    basis_frequencies: tuple[int, int] = (2, 2)

    @property
    def dimension(self):
        return self.occupation**2

    def orthonormal_matrix(self):
        """Fresh floating diagnostic, with the same exact metric convention."""
        return _orthonormal_matrix(self.entries, self.metric, self.dimension)


def _word_column(occupation, power, *, momentum=False):
    result = {occupation: F(1)}
    for _ in range(power):
        following = defaultdict(F)
        for index, value in result.items():
            following[index+1] += value*(1 if momentum else F(1, 2))
            if index:
                following[index-1] += value*index*(-1 if momentum else F(1, 2))
        result = following
    return {n: -v if momentum else v for n, v in result.items() if v}


def fixed_code_form(occupation, coupling=F(1)):
    """Independent exact H_d in the original two frequency-two Hermite modes.

    Form Q**2, P**2 and Q**4 on the infinite algebra first, then retain
    their output rows.  This is not the power of a compressed coordinate.
    For boundary residuals, build occupation d+4 and select the columns
    n1+(d+4)*n2 with n1,n2<d; all nonzero H J_d rows are then included.
    """
    occupation, coupling = _integer(occupation, 1), _rational(coupling)
    facts = tuple(factorial(n) for n in range(occupation))
    metric = tuple(facts[n1]*facts[n2]
                   for n2 in range(occupation) for n1 in range(occupation))
    entries = defaultdict(F)
    for n2 in range(occupation):
        for n1 in range(occupation):
            column = n1+occupation*n2
            q2 = (_word_column(n1, 2), _word_column(n2, 2))
            for mode, (index, frequency) in enumerate(((n1, 1), (n2, 2))):
                terms = ((_word_column(index, 2, momentum=True), F(1, 2)),
                         (q2[mode], F(frequency**2, 2)),
                         (_word_column(index, 4), coupling/16))
                for polynomial, coefficient in terms:
                    for row_index, value in polynomial.items():
                        if row_index < occupation:
                            row = row_index+occupation*n2 if mode == 0 else n1+occupation*row_index
                            entries[row, column] += coefficient*value
            for i, first in q2[0].items():
                for j, second in q2[1].items():
                    if i < occupation and j < occupation:
                        entries[i+occupation*j, column] += coupling*first*second/4
    nonzero = tuple((i, j, value) for (i, j), value in sorted(entries.items()) if value)
    return FixedCodeForm(coupling, occupation, nonzero, metric)


def _ladder_row(row, power):
    """Nonzero columns of (a+a*)**power in unnormalized Hermite coordinates."""
    if power == 1:
        result = {row+1: row+1}
        if row:
            result[row-1] = 1
        return result
    result = {row: 2*row+1, row+2: (row+1)*(row+2)}
    if row >= 2:
        result[row-2] = 1
    return result


def natural_lower_form(shell, coupling=F(1)):
    """Return H0 + sum c_j K_j* P_shell K_j on its finite support.

    H0 has energies n1+2*n2+3/2, K=(Q1**2,Q2**2,Q1*Q2), and
    c=(g/16,g/16,g/4).  P_shell retains 0 <= n1,n2 < shell.
    The perturbation is supported inside occupation ``shell+2``; on its
    orthogonal complement the infinite lower operator is exactly H0.
    Shell zero is the harmonic lower bound.  Increasing shell increases
    the infinite quadratic form, not necessarily every matrix entry.
    """
    shell, coupling = _integer(shell), _rational(coupling)
    occupation = shell+2
    facts = tuple(factorial(n) for n in range(occupation))
    metric = tuple(facts[n1]*facts[n2]
                   for n2 in range(occupation) for n1 in range(occupation))
    entries = defaultdict(F)
    for n2 in range(occupation):
        for n1 in range(occupation):
            address = n1+occupation*n2
            entries[address, address] = F(2*n1+4*n2+3, 2)

    def gram_row(row_metric, coefficients, weight):
        for left, first in coefficients.items():
            for right, second in coefficients.items():
                entries[left, right] += weight*first*second*F(row_metric, metric[left])

    # Squaring removes all irrational natural-coordinate normalization:
    # Q1**2=(a1+a1*)**2/2, Q2**2=(a2+a2*)**2/4,
    # Q1*Q2=(a1+a1*)(a2+a2*)/sqrt(8).
    for n2 in range(shell):
        for n1 in range(shell):
            row_metric = facts[n1]*facts[n2]
            first = {k+occupation*n2: v for k, v in _ladder_row(n1, 2).items()}
            second = {n1+occupation*k: v for k, v in _ladder_row(n2, 2).items()}
            mixed = {k+occupation*l: v*w
                     for k, v in _ladder_row(n1, 1).items()
                     for l, w in _ladder_row(n2, 1).items()}
            gram_row(row_metric, first, coupling/64)
            gram_row(row_metric, second, coupling/256)
            gram_row(row_metric, mixed, coupling/32)
    nonzero = tuple((i, j, value) for (i, j), value in sorted(entries.items()) if value)
    return FiniteLowerForm(shell, coupling, occupation, nonzero, metric)


def harmonic_outside_data(occupation, beta=F(1, 4)):
    """Exact formula data for Tr and Tr(H0 .) outside a natural-mode box.

    Let Z0=e**(-3*beta/2)/[(1-e**(-beta))*(1-e**(-2*beta))]
    and E0=3/2+1/(e**beta-1)+2/(e**(2*beta)-1).
    For each (sign, shift) in ``terms``, the outside trace contribution is
    sign*Z0*e**(-beta*shift), and the outside energy contribution is
    sign*Z0*e**(-beta*shift)*(E0+shift).  These inclusion-exclusion formulas
    charge all states outside the box, not an ordered spectral-rank tail.
    Evaluating the exponentials and diagonalizing the finite lower block
    require separate arithmetic certificates.
    """
    occupation = _integer(occupation, 1)
    beta = _rational(beta, positive=True)
    return {'occupation': occupation, 'inverse_temperature': beta,
            'frequencies': (1, 2), 'ground_energy': F(3, 2),
            'terms': ((1, F(occupation)), (1, F(2*occupation)),
                      (-1, F(3*occupation))),
            'finite_matrix_or_exponential_roundoff_certified': False}
