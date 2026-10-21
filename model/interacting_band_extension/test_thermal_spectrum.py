"""Independent arithmetic and corrupted-candidate controls for thermal traces."""
from dataclasses import replace
from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest

from . import thermal_forms as forms
from . import thermal_spectrum as spectrum


@pytest.mark.parametrize('value', (F(-500), F(-17, 3), F(-1, 100), F(0), F(1, 3), F(20)))
def test_rational_exponential_encloses_independent_high_precision(value):
    lower, upper = spectrum.exponential_interval(value)
    with mp.workdps(100):
        expected = mp.exp(mp.mpf(value.numerator)/value.denominator)
        assert mp.mpf(lower.numerator)/lower.denominator <= expected
        assert expected <= mp.mpf(upper.numerator)/upper.denominator


def test_limb_gram_matches_python_integer_oracle_with_wide_signed_entries():
    rng = np.random.default_rng(2188)
    matrix = rng.integers(-(1 << 50), 1 << 50, size=(13, 9), dtype=np.int64)
    expected = np.array([[sum(int(matrix[k, i])*int(matrix[k, j]) for k in range(13))
                          for j in range(9)] for i in range(9)], dtype=object)
    np.testing.assert_array_equal(spectrum.exact_gram(matrix), expected)
    with pytest.raises(ValueError, match='overflow'):
        spectrum.exact_gram(np.ones((16, 2), dtype=np.int64), limb_bits=30)
    with pytest.raises(ValueError, match='absolute'):
        spectrum.exact_gram(np.array([[np.iinfo(np.int64).min]], dtype=np.int64))


@pytest.mark.parametrize('coupling', (F(0), F(1)))
def test_certified_full_partition_contains_independent_finite_trace(coupling):
    form = forms.natural_lower_form(4, coupling)
    certificate = spectrum.certify_spectrum(form)
    matrix = form.orthonormal_matrix().toarray()
    energies = np.linalg.eigvalsh(matrix)
    for value, (lower, upper) in zip(energies, certificate['eigenvalue_intervals'], strict=True):
        assert float(lower)-1e-13 <= value <= float(upper)+1e-13
    lower, upper = certificate['partition_function_interval']
    assert float(lower) <= np.exp(-energies/4).sum() <= float(upper)
    assert float(upper-lower) < 1e-7
    if not coupling:
        harmonic = sorted(F(2*n+4*m+3, 2) for m in range(6) for n in range(6))
        assert all(lo <= energy <= hi for energy, (lo, hi) in
                   zip(harmonic, certificate['eigenvalue_intervals'], strict=True))
    assert not certificate['eigensolver_accuracy_assumed']


def _candidates(form):
    d = form.occupation
    matrix = form.orthonormal_matrix().toarray()
    return {parity: np.linalg.eigh(matrix[np.ix_(addresses, addresses)])
            for parity in ((0, 0), (0, 1), (1, 0), (1, 1))
            for addresses in [[i for i in range(d*d) if (i % d % 2, i//d % 2) == parity]]}


@pytest.mark.parametrize('mutation', ('zero_basis', 'nonfinite', 'boolean', 'shape', 'huge'))
def test_untrusted_eigensolver_candidates_fail_closed(mutation):
    form = forms.natural_lower_form(2)
    candidate = _candidates(form)
    values, vectors = candidate[0, 0]
    if mutation == 'zero_basis':
        vectors = np.zeros_like(vectors)
    elif mutation == 'nonfinite':
        vectors[0, 0] = np.nan
    elif mutation == 'boolean':
        vectors = vectors.astype(bool)
    elif mutation == 'shape':
        values = values[:-1]
    else:
        values[0] = 1e308
    candidate[0, 0] = values, vectors
    with pytest.raises((ValueError, ArithmeticError)):
        spectrum.certify_spectrum(form, candidate=candidate)


def test_wrong_eigenvalues_are_charged_instead_of_trusted():
    form = forms.natural_lower_form(2, F(0))
    candidate = _candidates(form)
    values, vectors = candidate[0, 0]
    candidate[0, 0] = values+F(1, 8), vectors
    candidate[0, 0] = np.asarray(candidate[0, 0][0], dtype=float), vectors
    certificate = spectrum.certify_spectrum(form, candidate=candidate)
    block = next(row for row in certificate['blocks'] if row['parity'] == (0, 0))
    assert block['eigenvalue_error_upper'] >= F(1, 8)
    exact_trace = sum((spectrum.exponential_interval(-F(2*n+4*m+3, 8))[0]
                       for m in range(4) for n in range(4)), F())
    lower, upper = certificate['partition_function_interval']
    assert lower <= exact_trace <= upper


def test_exact_metric_symmetry_and_sparse_addresses_are_checked():
    form = forms.natural_lower_form(2)
    with pytest.raises(ValueError, match='duplicate'):
        spectrum.certify_spectrum(replace(form, entries=form.entries+(form.entries[0],)))
    entries = list(form.entries)
    index = next(i for i, (row, column, _) in enumerate(entries) if row != column)
    row, column, value = entries[index]
    entries[index] = row, column, value+1
    with pytest.raises(ValueError, match='self-adjoint'):
        spectrum.certify_spectrum(replace(form, entries=tuple(entries)))


@pytest.mark.parametrize('coupling', (F(0), F(1)))
def test_original_frequency_two_code_has_a_certified_partition(coupling):
    form = forms.fixed_code_form(6, coupling)
    certificate = spectrum.certify_spectrum(form, bits=48)
    matrix = form.orthonormal_matrix().toarray()
    values = np.linalg.eigvalsh(matrix)
    lower, upper = certificate['partition_function_interval']
    assert float(lower) <= np.exp(-values/4).sum() <= float(upper)
    assert float(upper-lower) < 1e-8


def test_candidate_export_preserves_unsorted_column_pairing_and_is_opt_in():
    form = forms.fixed_code_form(4)
    candidate = _candidates(form)
    values, vectors = candidate[0, 0]
    candidate[0, 0] = values[::-1].copy(), vectors[:, ::-1].copy()
    result = spectrum.certify_spectrum(form, candidate=candidate, return_candidates=True)
    row = next(block for block in result['blocks'] if block['parity'] == (0, 0))
    assert row['dyadic_eigenvalues'][0] > row['dyadic_eigenvalues'][-1]
    assert row['addresses'] == (0, 2, 8, 10)
    assert type(row['dyadic_basis']) is tuple
    basis = np.array(row['dyadic_basis'], dtype=float)/row['dyadic_scale']
    diagonal = np.array(row['dyadic_eigenvalues'], dtype=float)/row['dyadic_scale']
    matrix = form.orthonormal_matrix().toarray()[np.ix_(row['addresses'], row['addresses'])]
    assert np.linalg.norm(matrix@basis-basis*diagonal, 'fro') < 1e-9
    plain = spectrum.certify_spectrum(form)
    assert all('dyadic_basis' not in block and 'addresses' not in block for block in plain['blocks'])
    with pytest.raises(ValueError, match='Boolean'):
        spectrum.certify_spectrum(form, return_candidates=1)
