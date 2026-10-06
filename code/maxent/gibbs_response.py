"""Gibbs Hessians with population and observable precision kept until output.

The pairwise Gram formula avoids subtracting two large mean products.
Independent precision agreement is a numerical check, not an interval
certificate. Inputs are validated by information_projection's public APIs.
"""

import math
from fractions import Fraction
from collections import defaultdict

import mpmath
import numpy as np

from .hamiltonian_assembly import _sum_components


def _diagonal_hamiltonian(operators, coefficients):
    real = _sum_components(operators, coefficients, 'real')
    imag = _sum_components(operators, coefficients, 'imag')
    real, imag = (real+real.T)/2, (imag-imag.T)/2
    if any(real[i, j] or imag[i, j] for i in range(len(real)) for j in range(i)):
        return None
    return list(real.diagonal())


def _rational_observables(operators):
    matrices = []
    for a in operators:
        r = np.array([Fraction(float(x)) for x in a.real.flat], dtype=object).reshape(a.shape)
        s = np.array([Fraction(float(x)) for x in a.imag.flat], dtype=object).reshape(a.shape)
        matrices.append(((r+r.T)/2, (s-s.T)/2))
    return matrices


def _uniform_covariance(operators):
    """Exact trace Gram form at I/d, including its exact zero entries.

    Arithmetic on the supplied binary64 entries is rational. Evaluate
    Tr(A B)/d-Tr(A)Tr(B)/d^2 before rounding, rather than comparing tiny
    numerical cancellation residues at two precisions.
    """
    matrices = _rational_observables(operators)
    size = len(operators[0])
    result = np.zeros((len(operators), len(operators)))
    for i, (r, s) in enumerate(matrices):
        for j, (u, v) in enumerate(matrices[:i+1]):
            pairing = sum((r*u+s*v).flat, Fraction(0))
            value = (size*pairing-sum(r.diagonal())*sum(u.diagonal()))/size**2
            try:
                rounded = float(value)
            except OverflowError as exc:
                raise ValueError("Duhamel covariance exceeds finite numerical range") from exc
            if value and (rounded == 0 or abs(Fraction(rounded)-value)
                          > Fraction(32*np.finfo(float).eps)*abs(value)):
                raise ValueError("Duhamel covariance underflow or insufficient output precision")
            result[i, j] = result[j, i] = rounded
    return result


def _spectral_zero_entries(operators, coefficients):
    """Sufficient zero certificate from rational eigenvalues/projectors.

    Z^2 C_ab is an exponential polynomial with rational coefficients when
    H has a rational spectrum. Group its coefficients exactly. Projectors
    are polynomials in the supplied H, so the certificate survives complex
    basis changes and degeneracy. An unsupported spectrum yields no zeros.
    """
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix

    def domain(r, s):
        return DomainMatrix.from_Matrix(sp.Matrix([
            [sp.Rational(x.numerator, x.denominator)
             + sp.I*sp.Rational(y.numerator, y.denominator) for x, y in zip(rr, ss)]
            for rr, ss in zip(r, s)])).convert_to(sp.QQ_I).to_dense()

    real = _sum_components(operators, coefficients, 'real')
    imag = _sum_components(operators, coefficients, 'imag')
    real, imag = (real+real.T)/2, (imag-imag.T)/2
    size = len(real)
    real -= real[0, 0]*np.eye(size, dtype=int)
    h = domain(real, imag)
    roots = h.to_Matrix().charpoly().as_poly().ground_roots()
    if sum(roots.values()) != size:
        return []
    energies = sorted(roots)
    identity = DomainMatrix.eye(size, sp.QQ_I).to_dense()
    projectors = []
    for energy in energies:
        projector = identity
        for other in energies:
            if other != energy:
                projector = projector.matmul(h-identity.scalarmul(sp.QQ_I.convert(other)))
                projector = projector.scalarmul(sp.QQ_I.convert(1/(energy-other)))
        projectors.append(projector)
    matrices = [domain(r, s) for r, s in _rational_observables(operators)]
    blocks = [[p.matmul(a) for a in matrices] for p in projectors]

    def trace_real(matrix):
        value = sp.QQ_I.to_sympy(sum(matrix[k, k].element for k in range(size)))
        real = sp.re(value)
        return Fraction(int(real.p), int(real.q))

    means = [[trace_real(a) for a in block] for block in blocks]
    zeros = []
    for a in range(len(operators)):
        for b in range(a+1):
            terms = defaultdict(Fraction)
            for r, hr in enumerate(energies):
                for s, hs in enumerate(energies):
                    terms[hr+hs] -= means[r][a]*means[s][b]
                diagonal = trace_real(blocks[r][a].matmul(blocks[r][b]))
                for hk, count in roots.items():
                    terms[hr+hk] += int(count)*diagonal
                for s in range(r):
                    hs = energies[s]
                    coefficient = 2*trace_real(blocks[r][a].matmul(blocks[s][b]))
                    gap = Fraction(int((hs-hr).p), int((hs-hr).q))
                    for hk, count in roots.items():
                        terms[hr+hk] += int(count)*coefficient/gap
                        terms[hs+hk] -= int(count)*coefficient/gap
            if not any(terms.values()):
                zeros.append((a, b))
    return zeros


class _PrecisionDisagreement(ValueError):
    pass


def _precision(operators, coefficients=()):
    spread = 0.
    for values in ([x for a in operators for part in (a.real, a.imag)
                    for x in part.flat if x], list(coefficients)):
        exponents = [math.log10(abs(float(x))) for x in values if x]
        if exponents:
            spread += max(exponents)-min(exponents)
    # A basis-change error is squared in a variance. A 550-digit floor
    # provides headroom for binary64 responses with large observables and
    # populations far below binary64 range; it is not an error certificate.
    # Input dynamic range may demand still more; agreement is checked below.
    return max(550, 2*math.ceil(spread)+80)


def _matrix(ctx, a):
    return ctx.matrix([[ctx.mpc(float(z.real), float(z.imag)) for z in row] for row in a])


def _observables(ctx, operators):
    result = []
    for a in operators:
        matrix = _matrix(ctx, a)
        matrix = (matrix+matrix.H)/2
        # Remove the scalar in this precision, before changing basis. A
        # huge energy origin must not create response via V*V roundoff.
        result.append(matrix-ctx.re(matrix[0, 0])*ctx.eye(len(a)))
    return result


def _gram(ctx, operators, probabilities, vectors, gaps):
    total = ctx.fsum(probabilities)
    p = [x/total for x in probabilities]
    rotated = [vectors.H*a*vectors for a in operators]
    pairs = [(i, j) for i in range(len(p)) for j in range(i)]
    kernels = []
    for i, j in pairs:
        gap = abs(gaps[i]-gaps[j])
        kernels.append(max(p[i], p[j])*(-ctx.expm1(-gap)/gap if gap else 1))
    result = ctx.zeros(len(operators))
    for a, left in enumerate(rotated):
        for b in range(a+1):
            right = rotated[b]
            terms = []
            for (i, j), kernel in zip(pairs, kernels):
                terms.append(p[i]*p[j]*ctx.re(left[i, i]-left[j, j])
                             *ctx.re(right[i, i]-right[j, j]))
                terms.append(2*kernel*ctx.re(left[i, j]*ctx.conj(right[i, j])))
            result[a, b] = result[b, a] = ctx.fsum(terms)
    return result


def _checked(evaluate, precision):
    contexts, values = [], []
    for digits in (precision, precision+40):
        ctx = mpmath.mp.clone()
        ctx.dps = digits
        contexts.append(ctx)
        values.append(evaluate(ctx))
    ctx = contexts[-1]
    coarse, refined = values
    result = np.zeros((refined.rows, refined.cols))
    for i in range(refined.rows):
        for j in range(i+1):
            value = refined[i, j]
            if abs(ctx.mpf(coarse[i, j])-value) > abs(value)*ctx.mpf('1e-25'):
                raise _PrecisionDisagreement("Duhamel covariance precision refinement disagrees")
            rounded = float(value)
            if not math.isfinite(rounded):
                raise ValueError("Duhamel covariance exceeds finite numerical range")
            if value and (rounded == 0 or abs(ctx.mpf(rounded)-value)
                          > 32*np.finfo(float).eps*abs(value)):
                raise ValueError("Duhamel covariance underflow or insufficient output precision")
            result[i, j] = result[j, i] = rounded
    return result


def covariance_from_spectrum(operators, probabilities, vectors):
    """Same Gram calculation on the optimizer's supplied binary64 spectrum.

    This does not recover the eigensolver's lost accuracy. The public
    response API instead recomputes the thermal family before rounding.
    """
    if np.all(probabilities == probabilities[0]):
        return _uniform_covariance(operators)

    def evaluate(ctx):
        p = [ctx.mpf(float(x)) for x in probabilities]
        return _gram(ctx, _observables(ctx, operators), p, _matrix(ctx, vectors),
                     [ctx.log(x) for x in p])
    return _checked(evaluate, _precision(operators))


def gibbs_covariance(operators, coefficients):
    """Hessian of log Tr exp(-sum lambda_a T_a) in supplied observable units.

    Reconstruct from original entries, not already-centered binary64
    matrices or rounded Gibbs populations. Only the final Hessian entries
    are converted to binary64. A state can be unrepresentable in binary64
    while its weighted response remains accurately representable.
    """
    diagonal = _diagonal_hamiltonian(operators, coefficients)
    if diagonal is not None and len(set(diagonal)) == 1:
        return _uniform_covariance(operators)
    zero_entries = []

    def evaluate(ctx):
        matrices = _observables(ctx, operators)
        hamiltonian = ctx.zeros(len(operators[0]))
        for coefficient, matrix in zip(coefficients, matrices):
            hamiltonian += ctx.mpf(float(coefficient))*matrix
        if all(hamiltonian[i, j] == 0 for i in range(hamiltonian.rows)
               for j in range(i)):
            energies = [ctx.re(hamiltonian[i, i]) for i in range(hamiltonian.rows)]
            vectors = ctx.eye(hamiltonian.rows)
        else:
            energies, vectors = ctx.eighe(hamiltonian)
        origin = min(energies)
        weights = [ctx.exp(-(e-origin)) for e in energies]
        result = _gram(ctx, matrices, weights, vectors, energies)
        for i, j in zero_entries:
            result[i, j] = result[j, i] = 0
        return result
    precision = _precision(operators, coefficients)
    try:
        return _checked(evaluate, precision)
    except _PrecisionDisagreement:
        zero_entries = _spectral_zero_entries(operators, coefficients)
        if not zero_entries:
            raise
        return _checked(evaluate, precision)
