"""Maintainer audit: certify uniform zeros without a cancellation tolerance."""

from fractions import Fraction

import numpy as np
import pytest

from maxent.information_projection import duhamel_covariance, _covariance
from maxent.gibbs_response import _spectral_zero_entries


def orthogonal_integer_pair(size, seed):
    rng = np.random.default_rng(seed)
    pair = []
    for _ in range(2):
        raw = rng.integers(-3, 4, (size, size))+1j*rng.integers(-3, 4, (size, size))
        a = raw+raw.conj().T
        pair.append(size*a-np.trace(a)*np.eye(size))
    a, c = pair
    # Exact integer Gram-Schmidt. Entries and every contraction are below
    # 2**53; this is an independent zero certificate, not a small tolerance.
    b = np.trace(a@a).real*c-np.trace(a@c).real*a
    assert np.trace(a) == np.trace(b) == np.trace(a@b) == 0
    return a, b


@pytest.mark.parametrize("size", [3, 5, 6])
@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("kind", ["zero_parameters", "scalar_sum", "newton"])
def test_exact_uniform_cross_zero_is_accepted(size, seed, kind):
    a, b = orthogonal_integer_pair(size, seed)
    if kind == "newton":
        result = _covariance([a, b], np.full(size, 1/size), np.eye(size))
    elif kind == "scalar_sum":
        result = duhamel_covariance([a, b, a+b, 1e30*np.eye(size)], [1., 1., -1., .25])
        assert np.count_nonzero(result[3]) == 0
    else:
        result = duhamel_covariance([a, b], [0., 0.])
    assert result[0, 1] == result[1, 0] == 0.
    assert result[0, 0] == pytest.approx(np.trace(a@a).real/size, rel=1e-15)
    assert result[1, 1] == pytest.approx(np.trace(b@b).real/size, rel=1e-15)


@pytest.mark.parametrize("small", [1e-20, 1e-150, 1e-300])
def test_uniform_exact_arithmetic_does_not_floor_a_small_coherence(small):
    a = np.diag([1., -1., 0.]).astype(complex)
    a[0, 2], a[2, 0] = small*1j, -small*1j
    b = np.diag([1., 1., -2.]).astype(complex)
    b[0, 2], b[2, 0] = 1j, -1j
    expected = float(2*Fraction(small)/3)
    result = duhamel_covariance([a, b], [0., 0.])
    assert result[0, 1] == pytest.approx(expected, rel=2e-15, abs=0)


def test_nonuniform_family_is_not_replaced_by_its_uniform_covariance():
    a = np.diag([1., -1.])
    beta = .75
    result = duhamel_covariance([a], [beta])
    assert result[0, 0] == pytest.approx(1/np.cosh(beta)**2, rel=2e-15)


@pytest.mark.parametrize("scale", [1e-200, 1e200])
def test_exact_uniform_result_still_obeys_output_range(scale):
    with pytest.raises(ValueError, match="precision|underflow|numerical range"):
        duhamel_covariance([scale*np.diag([1., -1.])], [0.])


@pytest.mark.parametrize("size", [3, 5, 6])
@pytest.mark.parametrize("seed", range(3))
def test_nonuniform_degenerate_thermal_sectors_keep_exact_zeros(size, seed):
    a, b = orthogonal_integer_pair(size, seed)
    aa = np.zeros((size+1, size+1), complex)
    bb = aa.copy()
    aa[:size, :size], bb[:size, :size] = a, b
    h = np.diag([0.]*size+[1.])
    beta = .7
    result = duhamel_covariance([h, aa, bb], [beta, 0., 0.])
    # Each observable has zero trace within each energy eigenspace and
    # Tr(a b)=0. Their weighted trace products vanish for every beta.
    assert result[0, 1] == result[0, 2] == result[1, 2] == 0
    partition = size+np.exp(-beta)
    assert result[1, 1] == pytest.approx(np.trace(a@a).real/partition, rel=2e-15)
    assert result[2, 2] == pytest.approx(np.trace(b@b).real/partition, rel=2e-15)


@pytest.mark.parametrize("small", [1e-20, 1e-150, 1e-300])
def test_exact_diagonal_zero_certificate_preserves_a_nearby_response(small):
    h = np.diag([0., 0., 0., 1.])
    a = np.diag([1., -1., small, 0.])
    b = np.diag([1., 1., -2., 0.])
    result = duhamel_covariance([h, a, b], [.7, 0., 0.])
    expected = -2*small/(3+np.exp(-.7))
    assert result[1, 2] == pytest.approx(expected, rel=3e-15, abs=0)


@pytest.mark.parametrize("size", [4, 8])
@pytest.mark.parametrize("imaginary", [False, True])
def test_exact_zero_certificate_survives_a_dense_dyadic_basis_change(size, imaginary):
    from scipy.linalg import hadamard
    a, b = orthogonal_integer_pair(size-1, 0)
    aa = np.zeros((size, size), complex)
    bb = aa.copy()
    aa[:-1, :-1], bb[:-1, :-1] = a, b
    h = np.diag([0.]*(size-1)+[1.])
    u = np.kron(np.eye(size//4), hadamard(4)/2).astype(complex)
    if imaginary:
        u[:, 0] *= 1j
    operators = [u@matrix@u.conj().T for matrix in [h, aa, bb]]
    result = duhamel_covariance(operators, [.7, 0., 0.])
    assert result[0, 1] == result[0, 2] == result[1, 2] == 0
    assert result[1, 1] == pytest.approx(np.trace(a@a).real/(size-1+np.exp(-.7)), rel=2e-15)
    assert result[2, 2] == pytest.approx(np.trace(b@b).real/(size-1+np.exp(-.7)), rel=2e-15)


@pytest.mark.parametrize("small", [1e-20, 1e-150, 1e-300])
def test_dense_certificate_does_not_promote_nearby_nonzero_to_zero(small):
    a, b = orthogonal_integer_pair(3, 0)
    aa, bb, h = (np.zeros((8, 8), complex) for _ in range(3))
    aa[:3, :3], bb[:3, :3] = a, b
    h[6:, 6:] = .5  # exact rank-one spectral projector, eigenvalues 0 and 1
    aa[0, 6:] = aa[6:, 0] = small
    bb[0, 6:] = bb[6:, 0] = 1.
    certificates = _spectral_zero_entries([h, aa, bb], [.7, 0., 0.])
    assert (1, 0) in certificates and (2, 0) in certificates
    assert (2, 1) not in certificates
    result = duhamel_covariance([h, aa, bb], [.7, 0., 0.])
    expected = 4*small*(-np.expm1(-.7))/(.7*(7+np.exp(-.7)))
    assert result[1, 2] == pytest.approx(expected, rel=3e-15, abs=0)


def test_unsupported_exact_spectrum_does_not_certify_a_numerical_zero():
    h = np.array([[1., 1.], [1., -1.]], complex)  # irrational eigenvalues
    assert _spectral_zero_entries([h, np.diag([1., -1.])], [1., 0.]) == []


def test_zero_certificate_retains_interenergy_terms_and_multiplicities():
    # Explicit two-level coherences in energy spaces of dimensions 2,1,1.
    # All three gaps contribute; complex conjugation cannot be omitted.
    from scipy.linalg import expm, expm_frechet
    h = np.diag([0., 0., 1., 3.]).astype(complex)
    a = np.array([[1, 1j, 2, -1j], [-1j, -1, 1j, 3],
                  [2, -1j, 0, 1], [1j, 3, 1, 0]], complex)
    b = np.array([[0, 1, 1j, 2], [1, 0, -2j, -1],
                  [-1j, 2j, 1, 1j], [2, -1, -1j, -1]], complex)
    certificates = _spectral_zero_entries([h, a, b], [1., 0., 0.])
    assert (2, 1) not in certificates
    exponential, derivative = expm(-h), expm_frechet(-h, -b, compute_expm=False)
    partition = np.trace(exponential).real
    expected = -np.trace((derivative/partition-exponential*np.trace(derivative)/partition**2)@a).real
    assert abs(expected) > .01
    assert duhamel_covariance([h, a, b], [1., 0., 0.])[1, 2] == pytest.approx(expected, rel=3e-14)


def test_product_thermal_factors_need_the_mean_product_subtraction():
    left = np.diag([0., 0., 1., 1.]).astype(complex)
    right = np.diag([0., 1., 0., 1.]).astype(complex)
    certificates = _spectral_zero_entries([left, right], [1., 2.])
    assert (1, 0) in certificates
    result = duhamel_covariance([left, right], [1., 2.])
    assert result[0, 1] == 0
    for i, energy in enumerate([1., 2.]):
        p = 1/(1+np.exp(energy))
        assert result[i, i] == pytest.approx(p*(1-p), rel=2e-15)


def test_interenergy_variance_cannot_be_certified_as_zero():
    z = np.diag([1., -1.]).astype(complex)
    x = np.array([[0., 1.], [1., 0.]], complex)
    assert (1, 1) not in _spectral_zero_entries([z, x], [.7, 0.])
