"""Exact dyadic conditional products and genuine support escapes (#1049)."""

from fractions import Fraction

import numpy as np
import pytest

from quantum_information import conditional_mutual_information, is_markov_exact
from quantum_information.information import _relative_modular_moments
from collar_alignment.msa_characterizations import (
    alignment_information_budget, collar_cmi,
)


def conditional_product(singular=False):
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    left = [np.eye(2)/2 + (x+y)/8, np.diag([.75, .25])]
    right = [np.diag([.75, .25]), np.diag([.625, .375])]
    if singular:
        left[0] = (np.eye(2)+y)/2
    rho = np.zeros((8, 8), complex)
    for b in range(2):
        ix = [4*a+2*b+c for a in range(2) for c in range(2)]
        rho[np.ix_(ix, ix)] = np.kron(left[b], right[b])/2
    return rho


@pytest.mark.parametrize("singular", [False, True])
@pytest.mark.parametrize("api", ["exact", "cmi", "collar", "budget"])
def test_complex_conditional_products_have_zero_cmi(singular, api):
    rho = conditional_product(singular)
    # B is a classical flag; each conditional A,C state is a product.
    # This certificate uses the construction, not an entropy subtraction.
    if api == "exact":
        assert is_markov_exact(rho, [2, 2, 2], [0], [1], [2])
    elif api == "cmi":
        assert conditional_mutual_information(rho, [2, 2, 2], [0], [1], [2]) == 0
    elif api == "collar":
        assert collar_cmi([(1., rho, (2, 2, 1, 2))]) == 0
    else:
        result = alignment_information_budget([(1., rho, (2, 2, 1, 2))])
        assert result["weighted"]["collar"] == 0
        assert result["weighted"]["alignment"] > 0
        assert abs(result["sectors"][0]["chain_rule_residual"]) < 1e-15


def test_complex_nonmarkov_state_is_false_without_support_error():
    rho = conditional_product()
    # A within-sector diagonal correlation keeps both marginals fixed.
    for i, sign in zip([0, 1, 4, 5], [1, -1, -1, 1]):
        rho[i, i] += sign/128
    assert np.linalg.eigvalsh(rho).min() > .03
    assert not is_markov_exact(rho, [2, 2, 2], [0], [1], [2])
    assert conditional_mutual_information(rho, [2, 2, 2], [0], [1], [2]) > 0


def entries(matrix):
    return [[(Fraction(float(z.real)), Fraction(float(z.imag))) for z in row]
            for row in np.asarray(matrix, complex)]


@pytest.mark.parametrize("imaginary", [False, True])
def test_exact_support_check_still_rejects_escape(imaginary):
    vector = np.array([1, 1j if imaginary else 1])/2
    reference = 2*np.outer(vector, vector.conj())
    with pytest.raises(ValueError, match="does not contain"):
        next(_relative_modular_moments(entries(np.eye(2)/2), entries(reference)))
