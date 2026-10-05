"""Information without subtracting binary64 entropies.

Reductions retain exact sums of the supplied binary64 entries. Classical
information is a sum of scalar relative-entropy remainders. Quantum entropy
combinations use private multiprecision arithmetic, with exact product/rank
checks for zero cases and explicit refusal of unresolved precision.
"""

from fractions import Fraction
from itertools import product
from math import prod

import mpmath
import numpy as np

from .gibbs import _numeric
from .states import _parts, density_matrix


_ZERO = (Fraction(0), Fraction(0))


def _add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _multiply(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def _mp_real(ctx, value):
    return ctx.mpf(value.numerator)/value.denominator


def _as_float(ctx, value):
    result = float(value)
    if not np.isfinite(result) or (value != 0 and (
            result == 0 or abs(ctx.mpf(result)-value) > 32*np.finfo(float).eps*abs(value))):
        raise ValueError("information underflow or insufficient output precision")
    return result


class _Reductions:
    """Exact linear reductions; no intermediate density gets renormalized."""

    def __init__(self, rho, dims):
        raw = _numeric(rho, "information state")
        a = density_matrix(raw)
        if len(a) != prod(dims):
            raise ValueError("subsystem dimensions do not match density matrix")
        self.dims = dims
        self.indices = list(product(*(range(d) for d in dims)))
        # Apply the validator's Hermitian-roundoff convention exactly, so
        # averaging cannot erase a nonzero subnormal component first.
        self.full = [[((Fraction(float(raw[i, j].real))+Fraction(float(raw[j, i].real)))/2,
                       (Fraction(float(raw[i, j].imag))-Fraction(float(raw[j, i].imag)))/2)
                      for j in range(len(a))] for i in range(len(a))]
        self.cache = {tuple(range(len(dims))): self.full}
        self.trace = sum(self.full[i][i][0] for i in range(len(a)))
        # Information inequalities do not apply to indefinite matrices
        # admitted by a roundoff tolerance. Validate the actual rational
        # Hermitian representative, even if a partial trace hides a defect.
        _positive_rank(self.full)

    def index(self, coordinate, keep):
        result = 0
        for k in keep:
            result = result*self.dims[k]+coordinate[k]
        return result

    def get(self, keep):
        keep = tuple(sorted(keep))
        if keep not in self.cache:
            size = prod(self.dims[k] for k in keep)
            result = [[_ZERO for _ in range(size)] for _ in range(size)]
            omitted = tuple(k for k in range(len(self.dims)) if k not in keep)
            mapped = [self.index(index, keep) for index in self.indices]
            for i, left in enumerate(self.indices):
                for j, right in enumerate(self.indices):
                    if all(left[k] == right[k] for k in omitted):
                        r, c = mapped[i], mapped[j]
                        result[r][c] = _add(result[r][c], self.full[i][j])
            self.cache[keep] = result
        return self.cache[keep]

    def factorizes(self, x, y):
        """Exact rho_XY Tr(rho) = rho_X tensor rho_Y, in original order."""
        x, y = tuple(sorted(x)), tuple(sorted(y))
        selected = tuple(sorted(x+y))
        joint, left, right = self.get(selected), self.get(x), self.get(y)
        coordinates = list(product(*(range(self.dims[k]) for k in selected)))
        maps = []
        for coord in coordinates:
            full = [0]*len(self.dims)
            for k, value in zip(selected, coord):
                full[k] = value
            maps.append((self.index(full, x), self.index(full, y)))
        for i, (a, b) in enumerate(maps):
            for j, (c, d) in enumerate(maps):
                if tuple(z*self.trace for z in joint[i][j]) != _multiply(left[a][c], right[b][d]):
                    return False
        return True

    def tensor_reference(self, x, y):
        """rho_X tensor rho_Y, reordered to the joint subsystem order."""
        x, y = tuple(sorted(x)), tuple(sorted(y))
        selected = tuple(sorted(x+y))
        left, right = self.get(x), self.get(y)
        maps = []
        for coord in product(*(range(self.dims[k]) for k in selected)):
            full = [0]*len(self.dims)
            for k, value in zip(selected, coord):
                full[k] = value
            maps.append((self.index(full, x), self.index(full, y)))
        return [[_multiply(left[a][c], right[b][d]) for c, d in maps] for a, b in maps]


def _diagonal(matrix):
    return all(value == _ZERO for i, row in enumerate(matrix)
               for j, value in enumerate(row) if i != j)


def _conditional_cells(data, a, b, c):
    selected = tuple(sorted(a+b+c))
    joint, ab, bc, middle = (data.get(part) for part in (selected, a+b, b+c, b))
    coordinates = product(*(range(data.dims[k]) for k in selected))
    for i, coord in enumerate(coordinates):
        full = [0]*len(data.dims)
        for k, value in zip(selected, coord):
            full[k] = value
        ix, iy, iz = (data.index(full, sorted(part)) for part in (a+b, b+c, b))
        p = joint[i][i][0]
        mass = middle[iz][iz][0]
        if p < 0 or mass < 0:
            raise ValueError("information requires a positive semidefinite state")
        if mass == 0:
            continue
        q = ab[ix][ix][0]*bc[iy][iy][0]/mass
        yield p, q


def _classical(data, a, b, c):
    """D(p_ABC || p_AB p_BC / p_B), using exact conditional products."""
    ctx = mpmath.mp.clone()
    ctx.dps = 80
    terms = []
    for p, q in _conditional_cells(data, a, b, c):
        if p == q:
            continue
        if p == 0:
            terms.append(_mp_real(ctx, q))
        elif q <= 0:
            raise ValueError("inconsistent classical information support")
        else:
            ratio = (p-q)/q
            r, mq = _mp_real(ctx, ratio), _mp_real(ctx, q)
            if abs(r) < ctx.mpf('0.125'):
                # ((1+r)log(1+r)-r)/r^2; all differences precede rounding.
                g = ctx.mpf(0)
                for k in reversed(range(100)):
                    g = -r*g + ctx.mpf(1)/((k+1)*(k+2))
                terms.append(mq*r*r*g)
            else:
                mp = _mp_real(ctx, p)
                terms.append(mp*(ctx.log(mp)-ctx.log(mq))-mp+mq)
    return _as_float(ctx, ctx.fsum(terms))


def _characteristic(matrix):
    import sympy as sp
    exact = sp.Matrix([
        [sp.Rational(r.numerator, r.denominator)
         + sp.I*sp.Rational(i.numerator, i.denominator) for r, i in row]
        for row in matrix])
    return [Fraction(int(x.p), int(x.q)) for x in exact.charpoly().all_coeffs()]


def _positive_rank(matrix):
    """Exact PSD/rank for a Hermitian rational matrix, without a cutoff.

    The coefficients of det(t I + A) are nonnegative iff A is PSD: a
    negative eigenvalue would give a positive root of that polynomial.
    Trailing zero coefficients give the exact zero-eigenvalue multiplicity.
    """
    coefficients = _characteristic(matrix)
    if any((-1)**k*value < 0 for k, value in enumerate(coefficients)):
        raise ValueError("information requires a positive semidefinite state")
    return max(k for k, value in enumerate(coefficients) if value)


def _power_traces(matrix):
    """Exact Newton sums, without building large tensor-product matrices."""
    coefficients = _characteristic(matrix)[1:]
    n, powers = len(matrix), [Fraction(len(matrix))]
    k = 1
    while True:
        value = -sum(coefficients[j-1]*powers[k-j] for j in range(1, min(k, n+1)))
        if k <= n:
            value -= k*coefficients[k-1]
        powers.append(value)
        yield value
        k += 1


def _isospectral_entropy_identity(data, a, b, c):
    """Sufficient exact CMI-zero certificate, including recharted products.

    Compare rho_AB tensor rho_BC and rho_B tensor rho_ABC by the first N
    power traces. Newton identities then give identical spectra; equal
    traces of all reductions imply the desired entropy identity. This is
    not a necessary test for every possible multi-sector Markov state.
    """
    matrices = [data.get(part) for part in (a+b, b+c, b, a+b+c)]
    count = len(matrices[0])*len(matrices[1])
    sequences = [_power_traces(matrix) for matrix in matrices]
    for _ in range(count):
        p, q, r, s = (next(sequence) for sequence in sequences)
        if p*q != r*s:
            return False
    return True


def _relative_modular_moments(source, reference):
    """Exact Tr(rho^k sigma^(1-k)); pseudoinverse only on exact support."""
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix

    def matrix(entries):
        return sp.Matrix([[sp.Rational(r.numerator, r.denominator)
                           + sp.I*sp.Rational(i.numerator, i.denominator)
                           for r, i in row] for row in entries])

    rho, sigma = matrix(source), matrix(reference)
    try:
        inverse = sigma.inv()
    except sp.matrices.exceptions.NonInvertibleMatrixError:
        inverse = sigma.pinv()
    if sigma*inverse*rho != rho:
        raise ValueError("relative modular reference does not contain the state support")
    r = DomainMatrix.from_Matrix(rho).convert_to(sp.QQ_I)
    inverse = DomainMatrix.from_Matrix(inverse).convert_to(sp.QQ_I)
    left = r
    right = DomainMatrix.eye(len(source), sp.QQ_I)
    while True:
        value = left.matmul(right)
        yield sum(value[i, i].element for i in range(len(source)))
        left = left.matmul(r)
        right = right.matmul(inverse)


def _exact_markov(data, a, b, c):
    """Complete finite CMI-zero test on exact PSD input; see proof note.

    Relative modular operators have at most n^2 spectral points. Equality
    of N positive moments, for N=n_ABC^2+n_AB^2, identifies their weighted
    positive spectra by a Vandermonde system, hence their t log(t) values.
    HJPW factorization gives the converse, including singular supports.
    """
    full, reduced = data.get(a+b+c), data.get(a+b)
    full_moments = _relative_modular_moments(full, data.tensor_reference(a, b+c))
    reduced_moments = _relative_modular_moments(reduced, data.tensor_reference(a, b))
    for _ in range(len(full)**2+len(reduced)**2):
        if next(full_moments) != next(reduced_moments):
            return False
    return True


def _entropy(ctx, matrix):
    if _diagonal(matrix):
        values = [_mp_real(ctx, matrix[i][i][0]) for i in range(len(matrix))]
    else:
        a = ctx.matrix([[ctx.mpc(_mp_real(ctx, r), _mp_real(ctx, i))
                         for r, i in row] for row in matrix])
        values = list(ctx.eighe(a, eigvals_only=True))
        threshold = ctx.eps*len(matrix)**3*100
        uncertain = [i for i, value in enumerate(values) if abs(value) <= threshold]
        if uncertain:
            if _positive_rank(matrix) != len(matrix)-len(uncertain):
                raise ArithmeticError("information spectrum is unresolved at this precision")
            for i in uncertain:
                values[i] = ctx.mpf(0)
    if any(value < 0 for value in values):
        raise ValueError("information requires a positive semidefinite state")
    return -ctx.fsum(value*ctx.log(value) for value in values if value)


def _quantum(data, a, b, c):
    # A product over any declared split of B gives zero CMI exactly. The
    # test is sufficient, not a complete recognition of the HJPW structure.
    exact_product = False
    for mask in range(1 << len(b)):
        left = a+tuple(k for j, k in enumerate(b) if mask & (1 << j))
        right = c+tuple(k for j, k in enumerate(b) if not mask & (1 << j))
        if data.factorizes(left, right):
            exact_product = True
            break
    parts = ((1, a+b), (1, b+c), (-1, b), (-1, a+b+c))
    # 400 digits leave a wide margin below every representable binary64
    # result. A cancellation below that resolution is refused, not floored.
    for precision in (80, 400):
        ctx = mpmath.mp.clone()
        ctx.dps = precision
        try:
            # Positivity is checked even for an exact product: an indefinite
            # product is not evidence of independent physical subsystems.
            entropies = {tuple(sorted(a+b+c)): _entropy(ctx, data.get(a+b+c))}
            if exact_product:
                return 0.0
            for _, part in parts:
                key = tuple(sorted(part))
                if key not in entropies:
                    entropies[key] = _entropy(ctx, data.get(part))
            terms = [sign*entropies[tuple(sorted(part))] for sign, part in parts]
        except ArithmeticError:
            if precision == 80:
                continue
            raise ValueError("information spectrum remains unresolved at available precision")
        value = ctx.fsum(terms)
        resolution_scale = ctx.eps*max(1, sum(abs(term) for term in terms))*len(data.full)**4*1000
        if value > 10**20*resolution_scale:
            # A second precision must agree; the scale guard alone is not
            # an a posteriori interval certificate for an eigensolver.
            refined = mpmath.mp.clone()
            refined.dps = precision+40
            try:
                check = refined.fsum(sign*_entropy(refined, data.get(part)) for sign, part in parts)
            except ArithmeticError:
                if precision == 80:
                    continue
                raise ValueError("information refinement has unresolved spectral precision")
            if check > 0 and abs(refined.mpf(value)-check) <= abs(check)*refined.mpf('1e-25'):
                return _as_float(refined, check)
            if precision == 80:
                continue
            raise ValueError("information evaluation fails the precision agreement check")
        if value < -10**20*resolution_scale:
            raise ValueError("negative information: state positivity or precision is unresolved")
    if _isospectral_entropy_identity(data, a, b, c) or _exact_markov(data, a, b, c):
        return 0.0
    raise ValueError("information cancellation is unresolved at available precision")


def weighted_information(weights, values):
    """Sum positive weighted diagnostics without silently underflowing terms."""
    weights = _numeric(weights, "information weights", real=True)
    values = _numeric(values, "information values", real=True)
    if (weights.ndim != 1 or values.shape != weights.shape or not len(weights)
            or np.any(weights < 0) or np.any(values < 0)):
        raise ValueError("paired nonnegative information weights and values required")
    exact = sum(Fraction(float(p))*Fraction(float(value)) for p, value in zip(weights, values))
    ctx = mpmath.mp.clone()
    ctx.dps = 80
    return _as_float(ctx, _mp_real(ctx, exact))


def conditional_information(rho, dims, part_a, part_b, part_c):
    """Tr(rho) I(A:C|B) of rho/Tr(rho), with exact linear reductions.

    This equals the usual four-entropy combination, also for accepted trace
    roundoff. With empty B it gives homogeneous mutual information, including
    the +Tr(rho) log Tr(rho) term. No state entries are normalized or clipped.
    """
    dims, (a, b, c) = _parts(dims, (part_a, part_b, part_c))
    data = _Reductions(rho, dims)
    selected = data.get(a+b+c)
    if _diagonal(selected):
        return _classical(data, a, b, c)
    return _quantum(data, a, b, c)


def is_markov_exact(rho, dims, part_a, part_b, part_c):
    """Decide zero CMI algebraically for the supplied finite rational state.

    Exact means the validated binary64 entries and their exactly symmetrized
    Hermitian representative. This is neither a noise tolerance nor a claim
    about an unknown physical state. Nontrivial relative-modular moment
    checks can be expensive; the interface targets small finite witnesses.
    """
    dims, (a, b, c) = _parts(dims, (part_a, part_b, part_c))
    data = _Reductions(rho, dims)
    if _diagonal(data.get(a+b+c)):
        return all(p == q for p, q in _conditional_cells(data, a, b, c))
    for mask in range(1 << len(b)):
        left = a+tuple(k for j, k in enumerate(b) if mask & (1 << j))
        right = c+tuple(k for j, k in enumerate(b) if not mask & (1 << j))
        if data.factorizes(left, right):
            return True
    return _isospectral_entropy_identity(data, a, b, c) or _exact_markov(data, a, b, c)
