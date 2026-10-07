"""Separate matrix-exponential, analytic and covariance-law controls."""

import mpmath
import numpy as np
import pytest
from scipy.linalg import expm, expm_frechet

from maxent.information_projection import duhamel_covariance, gibbs_state


X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])


def frechet_control(operators, multipliers):
    minus_h = -sum(c*a for c, a in zip(multipliers, operators))
    exponential = expm(minus_h)
    partition = np.trace(exponential).real
    result = np.empty((len(operators), len(operators)))
    for b, direction in enumerate(operators):
        derivative = expm_frechet(minus_h, -direction, compute_expm=False)
        response = derivative/partition-exponential*np.trace(derivative)/partition**2
        for a, observable in enumerate(operators):
            result[a, b] = -np.trace(response@observable).real
    return result


@pytest.mark.parametrize("size", [2, 3, 4])
@pytest.mark.parametrize("seed", range(4))
def test_noncommuting_response_matches_frechet_derivatives(size, seed):
    rng = np.random.default_rng(920+seed)
    operators = []
    for _ in range(3):
        a = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
        operators.append((a+a.conj().T)/4)
    multipliers = [.31, -.27, .46]
    result = duhamel_covariance(operators, multipliers)
    np.testing.assert_allclose(result, frechet_control(operators, multipliers),
                               rtol=3e-13, atol=3e-14)
    assert np.linalg.eigvalsh(result).min() > 0


@pytest.mark.parametrize("axis", [X, X+2*Z, X+Y+2*Z])
@pytest.mark.parametrize("scale", [1e100, 1e150, 1e300])
def test_rare_dense_thermal_response_keeps_the_commuting_variance(axis, scale):
    a = scale*axis
    beta = 1000/(2*scale*np.linalg.norm(axis, 2))
    ctx = mpmath.mp.clone()
    ctx.dps = 700
    # Independent two-level partition function: exact supplied components.
    e = ctx.sqrt(ctx.mpf(float(a[0, 0].real))**2
                 +ctx.mpf(float(a[0, 1].real))**2+ctx.mpf(float(a[0, 1].imag))**2)
    p = 1/(1+ctx.exp(2*ctx.mpf(beta)*e))
    expected = 4*e*e*p*(1-p)
    result = duhamel_covariance([a], [beta])
    assert result[0, 0] == pytest.approx(float(expected), rel=3e-13, abs=0)
    with pytest.raises(ValueError, match="faithful|population|underflow"):
        gibbs_state([a], [beta])


@pytest.mark.parametrize("seed", range(5))
def test_observable_recharts_and_complex_basis_changes(seed):
    rng = np.random.default_rng(180+seed)
    unitary, _ = np.linalg.qr(rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2)))
    transform = np.array([[2., -.5], [.25, 1.5]])
    shifts = [3., -5.]
    changed = [unitary@(sum(c*a for c, a in zip(row, [X, Z]))+s*np.eye(2))@unitary.conj().T
               for row, s in zip(transform, shifts)]
    lam = np.array([.2, -.4])
    expected = transform@frechet_control([X, Z], transform.T@lam)@transform.T
    np.testing.assert_allclose(duhamel_covariance(changed, lam), expected, rtol=3e-13, atol=3e-14)


def test_thermal_hessian_distinguishes_transverse_and_commuting_response():
    beta = 200.
    result = duhamel_covariance([Z, X, Y], [beta, 0., 0.])
    assert result[0, 0] == pytest.approx(1/np.cosh(beta)**2, rel=2e-14, abs=0)
    assert result[1, 1] == result[2, 2] == pytest.approx(np.tanh(beta)/beta, rel=2e-14)
    assert np.count_nonzero(result-np.diag(np.diag(result))) == 0


def test_identity_and_duplicate_observables_have_the_predicted_nullspace():
    result = duhamel_covariance([X, Z, 2*X-3*Z, 1e200*np.eye(2)], [.2, -.1, .3, 0.])
    np.testing.assert_allclose(result@[-2., 3., 1., 0.], 0, atol=2e-15)
    assert np.count_nonzero(result[3, :]) == np.count_nonzero(result[:, 3]) == 0
    assert np.linalg.matrix_rank(result, tol=1e-12) == 2


def test_mixed_observable_units_do_not_overflow_intermediate_products():
    ops = [1e150*X, 1e-150*Z]
    result = duhamel_covariance(ops, [0., 0.])
    np.testing.assert_allclose(np.diag(result), [1e300, 1e-300], rtol=2e-15, atol=0)
    assert result[0, 1] == 0


def test_large_scalar_origin_does_not_enter_the_response():
    actual = duhamel_covariance([1e200*np.eye(2)+X, Y], [.4, -.3])
    np.testing.assert_allclose(actual, frechet_control([X, Y], [.4, -.3]), rtol=2e-14, atol=0)


def test_private_multiprecision_does_not_change_global_settings():
    before = mpmath.mp.dps
    duhamel_covariance([X, Z], [.4, -.3])
    assert mpmath.mp.dps == before


@pytest.mark.parametrize("size", [4, 6, 8, 10])
def test_long_coherent_paths_do_not_become_false_zero_responses(size):
    scale, coefficient = 1e150, 1e-230
    chain = scale*(np.eye(size, k=1)+np.eye(size, k=-1))
    edge = np.zeros((size, size))
    edge[0, 1] = edge[1, 0] = scale
    end = np.zeros((size, size))
    end[0, -1] = end[-1, 0] = scale
    if size == 10:
        with pytest.raises(ValueError, match="precision|underflow"):
            duhamel_covariance([chain, edge, end], [coefficient, 0., 0.])
        return
    ctx = mpmath.mp.clone()
    ctx.dps = 700
    g = ctx.mpf(scale)*ctx.mpf(coefficient)
    # Independent walk expansion of exp(-g*path adjacency). The shortest
    # end-to-end path has size-1 edges, each used once. The next terms have
    # two extra steps, giving relative O(g^2) < 1e-150 here.
    leading = 2*ctx.mpf(scale)**2*g**(size-2)/(size*ctx.factorial(size-2))
    result = duhamel_covariance([chain, edge, end], [coefficient, 0., 0.])
    assert result[0, 2] == pytest.approx(float(leading), rel=2e-13, abs=0)
    assert result[1, 2] == pytest.approx(float(leading/(size-1)), rel=2e-13, abs=0)
