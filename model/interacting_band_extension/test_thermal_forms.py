"""Independent padded-oscillator controls for the Gibbs lower forms."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import thermal_forms as forms


def _padded_operators(size):
    a = np.diag(np.sqrt(np.arange(1, size)), 1)
    x = a+a.T
    q1, q2 = x/np.sqrt(2), x/2
    return q1, q2


@pytest.mark.parametrize('shell', (0, 1, 2, 5))
def test_exact_lower_form_matches_independent_padded_gram(shell):
    result = forms.natural_lower_form(shell, F(3, 2))
    n = result.occupation
    q1, q2 = _padded_operators(n+3)
    q1square, q2square = (q1@q1)[:n, :n], (q2@q2)[:n, :n]
    operators = (np.kron(np.eye(n), q1square),
                 np.kron(q2square, np.eye(n)),
                 np.kron(q2[:n, :n], q1[:n, :n]))
    rows = [i+n*j for j in range(shell) for i in range(shell)]
    expected = np.diag([i+2*j+1.5 for j in range(n) for i in range(n)])
    for coefficient, operator in zip((3/32, 3/32, 3/8), operators, strict=True):
        projected = operator[rows, :]
        expected += coefficient*projected.T@projected
    np.testing.assert_allclose(result.orthonormal_matrix().toarray(), expected,
                               atol=2e-13, rtol=2e-15)
    exact = {(i, j): value for i, j, value in result.entries}
    assert all(result.metric[i]*value == result.metric[j]*exact[j, i]
               for (i, j), value in exact.items())
    assert result.basis_frequencies == (1, 2)


def test_increasing_lower_forms_and_compression_have_distinct_boundary_terms():
    small, large = forms.natural_lower_form(2), forms.natural_lower_form(4)
    n = large.occupation
    embedded = np.diag([i+2*j+1.5 for j in range(n) for i in range(n)])
    addresses = [i+n*j for j in range(small.occupation) for i in range(small.occupation)]
    embedded[np.ix_(addresses, addresses)] = small.orthonormal_matrix().toarray()
    difference = large.orthonormal_matrix().toarray()-embedded
    assert np.linalg.eigvalsh(difference)[0] > -1e-13
    q1, q2 = _padded_operators(n+4)
    full = (np.kron(np.eye(n), (q1@q1@q1@q1)[:n, :n])/16
            + np.kron((q2@q2@q2@q2)[:n, :n], np.eye(n))/16
            + np.kron((q2@q2)[:n, :n], (q1@q1)[:n, :n])/4)
    full += np.diag([i+2*j+1.5 for j in range(n) for i in range(n)])
    assert np.linalg.eigvalsh(full-large.orthonormal_matrix().toarray())[0] > -1e-12
    assert np.linalg.norm(full-large.orthonormal_matrix().toarray()) > 1


def test_outside_partition_and_energy_formula_matches_explicit_sum():
    data = forms.harmonic_outside_data(5, F(1, 4))
    beta, n = float(data['inverse_temperature']), data['occupation']
    z0 = np.exp(-1.5*beta)/(-np.expm1(-beta)*-np.expm1(-2*beta))
    e0 = 1.5+1/np.expm1(beta)+2/np.expm1(2*beta)
    trace = z0*sum(sign*np.exp(-beta*float(shift)) for sign, shift in data['terms'])
    energy = z0*sum(sign*np.exp(-beta*float(shift))*(e0+float(shift))
                    for sign, shift in data['terms'])
    levels = [(i+2*j+1.5) for j in range(180) for i in range(180) if i >= n or j >= n]
    assert trace == pytest.approx(sum(np.exp(-beta*e) for e in levels), abs=5e-13)
    assert energy == pytest.approx(sum(e*np.exp(-beta*e) for e in levels), abs=1e-11)
    assert not data['finite_matrix_or_exponential_roundoff_certified']


def test_lower_partition_can_bound_original_frequency_two_compression():
    from .construction import ProductBand
    low = forms.natural_lower_form(8)
    n, beta = low.occupation, .25
    data = forms.harmonic_outside_data(n)
    z0 = np.exp(-1.5*beta)/(-np.expm1(-beta)*-np.expm1(-2*beta))
    upper = np.exp(-beta*np.linalg.eigvalsh(low.orthonormal_matrix().toarray())).sum()
    upper += z0*sum(sign*np.exp(-beta*float(shift)) for sign, shift in data['terms'])
    original = ProductBand(occupation=4).hamiltonian.toarray()
    lower = np.exp(-beta*np.linalg.eigvalsh(original)).sum()
    assert 3.64 < lower < 3.65 < upper
    # The new bracket does not relabel the original finite source matrix.
    assert np.linalg.norm(original-forms.natural_lower_form(2).orthonormal_matrix()) > 1


@pytest.mark.parametrize('occupation,coupling', ((2, F(0)), (4, F(1)), (7, F(3, 2))))
def test_original_code_form_matches_existing_production_without_changing_its_basis(occupation, coupling):
    from .construction import ProductBand
    form = forms.fixed_code_form(occupation, coupling)
    expected = ProductBand(coupling=coupling, occupation=occupation).hamiltonian.toarray()
    np.testing.assert_allclose(form.orthonormal_matrix().toarray(), expected, atol=3e-13, rtol=2e-15)
    assert form.basis_frequencies == (2, 2)
    entries = {(i, j): value for i, j, value in form.entries}
    assert all(form.metric[i]*value == form.metric[j]*entries[j, i]
               for (i, j), value in entries.items())


def test_four_extra_occupations_contain_the_entire_H_boundary_residual():
    d = 4
    little, big = forms.fixed_code_form(d+4), forms.fixed_code_form(d+7)
    columns = [i+little.occupation*j for j in range(d) for i in range(d)]
    bigger_columns = [i+big.occupation*j for j in range(d) for i in range(d)]
    rows = [i+big.occupation*j for j in range(little.occupation) for i in range(little.occupation)]
    full = big.orthonormal_matrix()[:, bigger_columns].toarray()
    expected = np.zeros_like(full)
    expected[rows] = little.orthonormal_matrix()[:, columns].toarray()
    np.testing.assert_array_equal(full, expected)
    assert np.linalg.norm(full[[i for i in range(big.dimension) if i not in bigger_columns]]) > 0


@pytest.mark.parametrize('value', (True, -1, 2., np.int64(2), '2'))
def test_shell_rejects_nonoriginal_or_negative_integers(value):
    with pytest.raises(ValueError):
        forms.natural_lower_form(value)
    with pytest.raises(ValueError):
        forms.fixed_code_form(value)


@pytest.mark.parametrize('value', (True, -1, .25, np.float64(.25), '1/4'))
def test_parameters_reject_inexact_or_negative_values(value):
    with pytest.raises(ValueError):
        forms.natural_lower_form(2, value)
    with pytest.raises(ValueError):
        forms.fixed_code_form(2, value)
    with pytest.raises(ValueError):
        forms.harmonic_outside_data(4, value)
    with pytest.raises(ValueError):
        forms.harmonic_outside_data(4, 0)
