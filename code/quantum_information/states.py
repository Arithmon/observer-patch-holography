"""Validated finite-dimensional density matrices, reductions and entropy.

Natural logarithms and the ordinary matrix trace are used. ATOL bounds
floating-point Hermiticity, normalization and negative-eigenvalue roundoff;
it is not a cutoff on positive eigenvalues or a physical uncertainty.
No positive eigenvalue is discarded or replaced by a logarithm floor.
These numerical diagnostics do not certify exact rank from approximate data.
"""

from math import frexp, fsum, nextafter, prod

import numpy as np

ATOL = 1e-12


def finite_real_scalar(value, name="value"):
    if np.ma.is_masked(value):
        raise ValueError(f"{name} must not be masked or missing")
    raw = np.asarray(value)
    if raw.ndim != 0 or raw.dtype.kind not in "iuf" or not np.isfinite(raw):
        raise ValueError(f"{name} must be a finite real scalar")
    result = float(raw)
    if not np.isfinite(result):
        raise ValueError(f"{name} exceeds finite numerical range")
    return result


def dimensions(dims):
    if not isinstance(dims, (list, tuple)) or not dims:
        raise ValueError("nonempty subsystem dimensions required")
    if any(isinstance(d, (bool, np.bool_)) or not isinstance(d, (int, np.integer))
           or d <= 0 for d in dims):
        raise ValueError("subsystem dimensions must be positive integers")
    return tuple(int(d) for d in dims)


def probabilities(values):
    raw = np.asarray(values)
    if raw.ndim != 1 or raw.size == 0 or raw.dtype.kind not in "iuf":
        raise ValueError("nonempty real probability vector required")
    p = raw.astype(float)
    if not np.all(np.isfinite(p)) or np.any(p < 0):
        raise ValueError("probabilities must be finite and nonnegative")
    if abs(float(p.sum()) - 1) > ATOL:
        raise ValueError("probabilities must sum to one; no implicit normalization")
    return p


def _spectrum(value):
    raw = np.asarray(value)
    if (raw.ndim != 2 or raw.shape[0] == 0 or raw.shape[0] != raw.shape[1]
            or raw.dtype.kind not in "iufc"):
        raise ValueError("nonempty square numeric density matrix required")
    a = raw.astype(complex)
    if not np.all(np.isfinite(a)):
        raise ValueError("density matrix must be finite")
    if np.linalg.norm(a - a.conj().T, ord="fro") > ATOL:
        raise ValueError("density matrix must be Hermitian")
    if abs(np.trace(a) - 1) > ATOL:
        raise ValueError("density matrix must have trace one")
    # Symmetrization is restricted to the declared roundoff tolerance.
    a = (a + a.conj().T) / 2
    eigenvalues, vectors = np.linalg.eigh(a)
    if eigenvalues[0] < -ATOL:
        raise ValueError("density matrix must be positive semidefinite")
    return a, np.maximum(eigenvalues, 0), vectors


def density_matrix(value):
    """Validate a normalized state; never normalize an invalid input."""
    return _spectrum(value)[0]


def _unresolved_positive_spectrum(a, eigenvalues):
    # A coordinate-diagonal input directly supplies its spectrum. For a
    # dense eigensolve, positive values on the roundoff scale do not resolve
    # rank. Refusing that case does not remove or floor an eigenvalue.
    roundoff = 8 * np.finfo(float).eps * len(a) * np.linalg.norm(a, ord="fro")
    return (np.any(a != np.diag(np.diag(a)))
            and np.any((eigenvalues > 0) & (eigenvalues <= roundoff)))


def _require_faithful(a, eigenvalues):
    if eigenvalues[0] <= 0:
        raise ValueError("state is not faithful; full matrix logarithm undefined")
    if _unresolved_positive_spectrum(a, eigenvalues):
        raise ValueError("faithful-state support is numerically unresolved")


def faithful_density_matrix(value):
    """Validate a state whose full support is resolved at this precision."""
    a, eigenvalues, _ = _spectrum(value)
    _require_faithful(a, eigenvalues)
    return a


def _indices(indices, count):
    if not isinstance(indices, (list, tuple)):
        raise ValueError("subsystem indices must be a list or tuple")
    if any(isinstance(i, (bool, np.bool_)) or not isinstance(i, (int, np.integer))
           or not 0 <= i < count for i in indices):
        raise ValueError("invalid subsystem index")
    if len(set(indices)) != len(indices):
        raise ValueError("duplicate subsystem index")
    return tuple(sorted(indices))


def partial_trace(rho, dims, keep):
    """Retain subsystems in their original order, including the empty set."""
    dims = dimensions(dims)
    keep = _indices(keep, len(dims))
    a = density_matrix(rho)
    if a.shape[0] != prod(dims):
        raise ValueError("subsystem dimensions do not match density matrix")
    tensor = a.reshape(dims + dims)
    current = len(dims)
    for i in reversed(range(len(dims))):
        if i not in keep:
            tensor = np.trace(tensor, axis1=i, axis2=i + current)
            current -= 1
    size = prod(dims[i] for i in keep)
    return tensor.reshape(size, size)


def von_neumann_entropy(rho):
    _, eigenvalues, _ = _spectrum(rho)
    positive = eigenvalues[eigenvalues > 0]
    return float(-np.sum(positive * np.log(positive)))


def shannon_entropy(ps):
    p = probabilities(ps)
    positive = p[p > 0]
    return float(-np.sum(positive * np.log(positive)))


def faithful_log(rho):
    """Logarithm on the full algebra; a singular state has no such log."""
    a, eigenvalues, vectors = _spectrum(rho)
    _require_faithful(a, eigenvalues)
    return (vectors * np.log(eigenvalues)) @ vectors.conj().T


def _exact_hermitian(value):
    """Rational Hermitian representative of validated binary64 entries."""
    from fractions import Fraction

    raw = np.asarray(value, dtype=complex)
    return [[((Fraction(float(raw[i, j].real))
               + Fraction(float(raw[j, i].real)))/2,
              (Fraction(float(raw[i, j].imag))
               - Fraction(float(raw[j, i].imag)))/2)
             for j in range(len(raw))] for i in range(len(raw))]


def _resolved_subnormal(ctx, value, resolution):
    """Round an MP component only when its subnormal midpoint is resolved.

    Resolution is the caller's numerical estimate, not an interval bound.
    Rational conversion avoids mpf-to-float double rounding at this boundary.
    """
    from fractions import Fraction
    from mpmath.libmp import to_rational

    if value == 0:
        return 0.0
    rounded = float(Fraction(*to_rational(value._mpf_)))
    if rounded == 0:
        raise ArithmeticError("nonzero component is below output resolution")
    lower = (ctx.mpf(rounded)+ctx.mpf(nextafter(rounded, -float("inf"))))/2
    upper = (ctx.mpf(rounded)+ctx.mpf(nextafter(rounded, float("inf"))))/2
    if min(abs(value-lower), abs(value-upper)) <= resolution:
        raise ArithmeticError("subnormal rounding midpoint is unresolved")
    return rounded


def _relative_entropy_refined(rho, sigma, *, subtract_trace=False):
    """Reevaluate range-sensitive contractions from the supplied entries.

    Exact Hermitian representatives determine positivity, rank and support.
    Two private numerical precisions must resolve the final cancellation and
    agree on its binary64 result; this is not an interval certificate.
    Consumers can request D-Tr(rho-sigma) before either term is rounded.
    """
    import mpmath
    from .information import _positive_rank

    a, b = _exact_hermitian(rho), _exact_hermitian(sigma)
    if a == b:
        return 0.0
    ranks = [_positive_rank(matrix) for matrix in (a, b)]
    if ranks[1] < len(b):
        import sympy as sp
        from sympy.polys.matrices import DomainMatrix
        matrices = [DomainMatrix.from_Matrix(sp.Matrix([
            [sp.Rational(r.numerator, r.denominator)
             + sp.I*sp.Rational(i.numerator, i.denominator) for r, i in row]
            for row in matrix])).convert_to(sp.QQ_I) for matrix in (a, b)]
        if not matrices[0].matmul(matrices[1].nullspace().transpose()).is_zero_matrix:
            return float("inf")

    def evaluate(precision):
        ctx = mpmath.mp.clone()
        ctx.dps = precision
        matrices = [ctx.matrix([[ctx.mpc(ctx.mpf(r.numerator)/r.denominator,
                                         ctx.mpf(i.numerator)/i.denominator)
                                 for r, i in row] for row in matrix])
                    for matrix in (a, b)]
        spectra = []
        for matrix, rank in zip(matrices, ranks):
            values, vectors = ctx.eighe(matrix)
            floor = 64*len(a)*ctx.eps*max(1, ctx.norm(matrix))
            positive = list(range(len(a)-rank, len(a)))
            if any(values[i] <= floor for i in positive):
                raise ArithmeticError("unresolved positive spectrum")
            spectra.append((values, vectors, positive))
        values, _, positive = spectra[0]
        terms = [values[i]*ctx.log(values[i]) for i in positive]
        values, vectors, positive = spectra[1]
        for i in positive:
            v = vectors[:, i]
            mass = ctx.re((v.H*matrices[0]*v)[0])
            terms.append(-mass*ctx.log(values[i]))
        if subtract_trace:
            terms.extend(ctx.re(matrices[1][i, i]-matrices[0][i, i])
                         for i in range(len(a)))
        result = ctx.fsum(terms)
        scale = max(1, ctx.fsum(abs(term) for term in terms))
        resolution = 256*len(a)*ctx.eps*scale
        if abs(result) <= resolution:
            raise ArithmeticError("unresolved entropy cancellation")
        if abs(result) < np.finfo(float).tiny:
            return _resolved_subnormal(ctx, result, resolution)
        return float(result)

    for precision in (80, 400, 800):
        try:
            result, check = evaluate(precision), evaluate(precision+40)
        except (ArithmeticError, ValueError):
            continue
        if result == check and np.isfinite(result) and result != 0:
            return result
    raise ValueError("relative entropy range or cancellation is numerically unresolved")


def _relative_entropy_remainder(a, b, eig_a, vec_a, eig_b, vec_b):
    """Stable trace Bregman remainder plus the supplied trace difference.

    Sum |<u_i,v_j>|^2 [p_i log(p_i/q_j) - p_i + q_j], then add
    Tr(a-b). Each bracket is nonnegative, so nearby normalized states do
    not require subtracting two entropy terms of order one.
    """
    difference = a - b
    # BLAS need not report an underflow inside the two contractions. Bound
    # the exponent of every nonzero real product before evaluating them.
    exponents = []
    for matrix in (vec_a, difference, vec_b):
        parts = np.concatenate((np.abs(matrix.real).ravel(),
                                np.abs(matrix.imag).ravel()))
        positive = parts[parts != 0]
        if positive.size:
            exponents.append(frexp(float(np.min(positive)))[1])
    if len(exponents) == 3 and sum(exponents) < -1018:
        raise FloatingPointError("relative entropy contraction range")
    action = vec_a.conj().T @ difference @ vec_b
    # This sum is exact up to its final rounding, unlike the spectral
    # terms below. Do not count the cancelling unit traces as numerical
    # uncertainty when both supplied traces coincide.
    trace_difference = fsum([float(x.real) for x in np.diag(a)]
                            + [-float(x.real) for x in np.diag(b)])
    terms = [trace_difference]
    for i, p in enumerate(eig_a):
        for j, q in enumerate(eig_b):
            if q == 0:
                # The caller has already checked containment of support.
                continue
            if p > 0 and abs(p - q) <= .5 * q:
                r = (p - q) / q
                # g(r) = ((1+r) log(1+r) - r)/r^2, with g(0)=1/2.
                # The series avoids cancellation in the scalar bracket.
                if abs(r) < .125:
                    g = 0.0
                    for k in reversed(range(18)):
                        g = -r * g + 1 / ((k + 1) * (k + 2))
                else:
                    g = ((1 + r) * np.log1p(r) - r) / (r * r)
                # (p-q)<u,v> = <u|(a-b)|v>. Evaluating the small
                # difference first avoids eigensolver trace and overlap
                # roundoff dominating a divergence below machine epsilon.
                # Divide before squaring to retain subnormal spectral mass.
                term = abs(action[i, j] / np.sqrt(q)) ** 2 * g
            else:
                # Complex division may form an overflowing reciprocal for
                # a subnormal gap even when both real quotients are bounded.
                weight = np.hypot(action[i, j].real / (p - q),
                                  action[i, j].imag / (p - q)) ** 2
                # Keep the distant bracket's terms separate for fsum. If
                # q >> p, forming the bracket first could lose p*log(p/q),
                # then cancellation with the trace term would expose the loss.
                if p > 0:
                    terms.extend((float(weight * (p * (np.log(p) - np.log(q)))),
                                  float(-weight * p)))
                term = weight * q
            terms.append(float(term))
    # Include the actual trace difference even inside the validation
    # tolerance: this API does not normalize its inputs or clip D to zero.
    result = fsum(terms)
    magnitude = fsum(abs(term) for term in terms)
    # The individual spectral terms are approximate even when fsum is
    # accurate. Cancellation can expose their roundoff, including against
    # an accepted trace offset. Reevaluate before that amplification
    # exceeds 2^10; the guard is not an interval error certificate.
    if magnitude and abs(result)*1024 <= magnitude:
        raise FloatingPointError("relative entropy cancellation")
    return result


def relative_entropy(rho, sigma):
    """Umegaki D(rho || sigma), with +infinity on detected support escape.

Ambiguous reference eigenvalues trigger exact positivity, rank and support
checks on the supplied Hermitian representative. Certified support escape
gives infinity, even if small; unresolved finite evaluation raises instead
of inventing either a finite answer or infinity. No leakage is rounded away.
These input-specific checks do not certify arbitrary numerical eigensolves.
Nearby PSD states use a compensated trace Bregman remainder. If the source
eigensolve clips spectral roundoff without an exactly annihilated kernel,
the legacy entropy/reference pairing is retained; the stable PSD identity
does not apply to an indefinite matrix admitted within the input tolerance.
"""
    from .gibbs import _numeric
    rho = _numeric(rho, "relative entropy source state")
    sigma = _numeric(sigma, "relative entropy reference state")
    a, eig_a, vec_a = _spectrum(rho)
    b, eig_b, vec_b = _spectrum(sigma)
    if a.shape != b.shape:
        raise ValueError("relative entropy requires the same algebra")
    if np.array_equal(a, b):
        if np.array_equal(np.asarray(rho), np.asarray(sigma)):
            return 0.0
        # Accepted Hermitian averaging can erase a half-subnormal entry.
        return _relative_entropy_refined(rho, sigma)
    # Either sign of a roundoff-sized eigenvalue can conceal full support.
    # Even a rounded zero matrix/kernel product is not an exact certificate:
    # certify the original representative before returning infinity.
    if np.any(eig_b == 0) or _unresolved_positive_spectrum(b, eig_b):
        return _relative_entropy_refined(rho, sigma)
    source_kernel = vec_a[:, eig_a == 0]
    if source_kernel.size and np.any(a @ source_kernel != 0):
        # The input tolerance also admits clipped negative spectral
        # roundoff. The PSD Bregman identity does not apply to that matrix;
        # retain the existing spectral-roundoff convention in this case.
        positive_a = eig_a[eig_a > 0]
        positive_b = eig_b > 0
        v = vec_b[:, positive_b]
        diagonal = np.real(np.diag(v.conj().T @ a @ v))
        return float(np.sum(positive_a * np.log(positive_a))
                     - np.dot(diagonal, np.log(eig_b[positive_b])))
    try:
        with np.errstate(over="raise", under="raise", invalid="raise", divide="raise"):
            return _relative_entropy_remainder(a, b, eig_a, vec_a, eig_b, vec_b)
    except FloatingPointError:
        return _relative_entropy_refined(rho, sigma)


def _parts(dims, parts):
    dims = dimensions(dims)
    groups = [_indices(part, len(dims)) for part in parts]
    flat = [i for group in groups for i in group]
    if len(flat) != len(set(flat)):
        raise ValueError("information-theoretic subsystems must be disjoint")
    return dims, groups


def mutual_information(rho, dims, part_x, part_y):
    """Homogeneous mutual information; preserve small correlation defects."""
    from .information import conditional_information
    return conditional_information(rho, dims, part_x, [], part_y)


def conditional_mutual_information(rho, dims, part_a, part_b, part_c):
    """Conditional information with exact reductions and guarded precision."""
    from .information import conditional_information
    return conditional_information(rho, dims, part_a, part_b, part_c)


def direct_sum_state(weights, states):
    """Retain a classical sector label; this is not an unlabelled mixture."""
    p = probabilities(weights)
    if len(p) != len(states):
        raise ValueError("one normalized state is required per sector weight")
    blocks = [density_matrix(state) for state in states]
    out = np.zeros((sum(len(block) for block in blocks),) * 2, dtype=complex)
    offset = 0
    for weight, block in zip(p, blocks):
        end = offset + len(block)
        scaled = weight * block
        if weight > 0 and (
                np.any((block.real != 0) & (scaled.real == 0))
                or np.any((block.imag != 0) & (scaled.imag == 0))):
            raise ValueError("direct-sum block underflow; precision is insufficient")
        out[offset:end, offset:end] = scaled
        offset = end
    return out


def one_sided_projection(operator, d_left, d_right):
    """HS projection onto M_left tensor 1 + 1 tensor M_right.

    This is a linear operator calculation, not a density-matrix operation.
    A multi-sector collar applies it separately inside each central block.
    """
    left, right = dimensions((d_left,d_right))
    a = np.asarray(operator)
    if (a.shape != (left*right,)*2 or a.dtype.kind not in "iufc"
            or not np.all(np.isfinite(a))):
        raise ValueError("finite operator with matching tensor dimensions required")
    a = a.astype(complex)
    tensor = a.reshape(left,right,left,right)
    l = np.trace(tensor,axis1=1,axis2=3)/right
    r = np.trace(tensor,axis1=0,axis2=2)/left
    return np.kron(l,np.eye(right))+np.kron(np.eye(left),r)-np.trace(a)/(left*right)*np.eye(left*right)
