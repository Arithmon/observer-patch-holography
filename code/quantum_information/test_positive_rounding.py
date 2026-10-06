"""Independent characteristic-polynomial checks of the returned Gram entries."""
import numpy as np
import pytest
import sympy as sp
from sympy.polys.matrices import DomainMatrix

from quantum_information.positive_rounding import _is_exact_psd, round_positive_gram


def exact(a):
    return sp.Matrix([[sp.Rational(float(z.real))+sp.I*sp.Rational(float(z.imag))
                       for z in row] for row in np.asarray(a, complex)])


def positive(a):
    a = DomainMatrix.from_Matrix(a.applyfunc(sp.expand)).convert_to(sp.QQ_I).to_Matrix()
    assert a == a.H
    return all((-1)**k*x >= 0 for k, x in enumerate(a.charpoly().all_coeffs()))


@pytest.mark.parametrize("shape", [(2, 1), (3, 2), (5, 4)])
@pytest.mark.parametrize("scale", [1e-200, 1e-100, 1., 1e150])
@pytest.mark.parametrize("imaginary", [False, True])
def test_returned_matrix_and_rounding_error_are_positive_with_certified_trace(shape, scale, imaginary):
    rng = np.random.default_rng(31)
    factor = rng.normal(size=shape).astype(complex)
    if imaginary:
        factor += 1j*rng.normal(size=shape)
    factor *= scale
    matrix, bound = round_positive_gram(factor)
    b = DomainMatrix.from_Matrix(exact(factor)).convert_to(sp.QQ_I).to_dense()
    adjoint = DomainMatrix.from_Matrix(exact(factor).H).convert_to(sp.QQ_I).to_dense()
    gram = b.matmul(adjoint).to_Matrix()
    difference = exact(matrix)-gram
    assert positive(exact(matrix))
    assert positive(difference)
    assert 0 <= sp.trace(difference) <= sp.Rational(bound)
    assert sp.Rational(bound)-sp.trace(difference) <= sp.Rational(float(np.spacing(bound)))
    # Independent row inequalities imply PSD and catch missing compensation
    # even where an eigensolver would round a negative direction to zero.
    for i in range(shape[0]):
        assert difference[i, i] >= sum(abs(sp.re(difference[i, j]))
                                      + abs(sp.im(difference[i, j]))
                                      for j in range(shape[0]) if j != i)


@pytest.mark.parametrize("size", [2, 3, 5])
def test_exact_schur_check_matches_an_independent_polynomial_control(size):
    rng = np.random.default_rng(52)
    for rank in range(1, size+1):
        factor = (rng.integers(-4, 5, (size, rank))
                  + 1j*rng.integers(-4, 5, (size, rank)))/8
        gram = factor@factor.conj().T  # Exact dyadic sums in this input range.
        assert _is_exact_psd(gram) and positive(exact(gram))
        changed = gram.copy()
        changed[0, 0] = -1/8
        assert not _is_exact_psd(changed) and not positive(exact(changed))
        # Permuting zero pivots exercises singular Schur complements.
        zero = np.zeros((size+1, size+1), complex)
        zero[1:, 1:] = gram
        assert _is_exact_psd(zero) and positive(exact(zero))


@pytest.mark.parametrize("small", [1e-20, 1e-200, np.nextafter(0., 1.)])
def test_zero_pivot_cannot_hide_an_indefinite_complex_direction(small):
    a = np.array([[0., 1j*small], [-1j*small, 1.]])
    assert not _is_exact_psd(a)
    assert not positive(exact(a))


def test_exactly_representable_gram_has_zero_rounding_bound():
    factor = np.array([[1., 1j], [0., .5], [-1j, 0.]])/4
    actual, bound = round_positive_gram(factor)
    np.testing.assert_array_equal(actual, factor@factor.conj().T)
    assert bound == 0


@pytest.mark.parametrize("bad", [[], [1.], [[True]], [[float("nan")]], [[float("inf")]],
                               np.ma.array([[1.]], mask=[[True]])])
def test_gram_rounding_rejects_invalid_input(bad):
    with pytest.raises(ValueError):
        round_positive_gram(bad)


def test_gram_rounding_refuses_output_overflow():
    with pytest.raises(ValueError, match="finite numerical range"):
        round_positive_gram([[1e308]])
