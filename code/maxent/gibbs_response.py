"""Gibbs Hessians with population and observable precision kept until output.

The pairwise Gram formula avoids subtracting two large mean products.
Independent precision agreement is a numerical check, not an interval
certificate. Inputs are validated by information_projection's public APIs.
"""

import math

import mpmath
import numpy as np


def _precision(operators, coefficients=()):
    spread = 0.
    for values in ([x for a in operators for part in (a.real, a.imag)
                    for x in part.flat if x], list(coefficients)):
        exponents = [math.log10(abs(float(x))) for x in values if x]
        if exponents:
            spread += max(exponents)-min(exponents)
    # A basis-change error is squared in a variance. A 550-digit floor
    # separates that error from binary64 responses even when observables
    # approach 1e308 and rare populations are far below binary64 range.
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
                raise ValueError("Duhamel covariance precision refinement disagrees")
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
        return _gram(ctx, matrices, weights, vectors, energies)
    return _checked(evaluate, _precision(operators, coefficients))
