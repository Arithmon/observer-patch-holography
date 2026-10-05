"""Maintainer-audit boundaries: exact zeros and irrelevant tensor factors."""

import numpy as np
import pytest

from quantum_information import conditional_mutual_information as cmi
from quantum_information import is_markov_exact, mutual_information as mi
from quantum_information.information import _Reductions


def _tiny_positive_gram():
    # Integer Gram entries fit exactly in binary64. The smallest eigenvalue
    # is positive but below 1e-400 after scaling; no entry is unrepresentable.
    upper = np.eye(6)+2.**26*np.eye(6, k=1)
    return np.nextafter(0., 1.)*(upper.T@upper)


def test_exact_product_zero_does_not_require_resolved_eigenvalues():
    left = np.zeros((7, 7))
    left[0, 0] = 1.
    left[1:, 1:] = _tiny_positive_gram()
    rho = np.kron(left, np.diag([1., 0.]))
    assert is_markov_exact(rho, [7, 2], [0], [], [1])
    assert mi(rho, [7, 2], [0], [1]) == 0.


def test_exact_conditional_zero_survives_unresolved_positive_spectrum():
    # B=C is a classical record. Conditional on it, A is independent of C.
    # Neither A:BC nor AB:C factorizes because both A and C identify B.
    rho = np.zeros((28, 28))
    rho[0, 0] = 1.
    indices = [4*a+3 for a in range(1, 7)]
    rho[np.ix_(indices, indices)] = _tiny_positive_gram()
    data = _Reductions(rho, (7, 2, 2))
    assert not data.factorizes((0,), (1, 2))
    assert not data.factorizes((0, 1), (2,))
    assert is_markov_exact(rho, [7, 2, 2], [0], [1], [2])
    assert cmi(rho, [7, 2, 2], [0], [1], [2]) == 0.


@pytest.mark.parametrize("evaluate", (cmi, is_markov_exact))
def test_trivial_conditioning_factors_do_not_repeat_physical_splits(monkeypatch, evaluate):
    x = np.array([[0., 1.], [1., 0.]])
    rho = np.eye(4)/4+.125*np.kron(x, x)
    calls = 0
    original = _Reductions.factorizes

    def counted(data, left, right):
        nonlocal calls
        calls += 1
        assert calls <= 4, "dimension-one factors repeated an identical physical split"
        return original(data, left, right)

    monkeypatch.setattr(_Reductions, "factorizes", counted)
    result = evaluate(rho, [2]+[1]*24+[2], [0], list(range(1, 25)), [25])
    if evaluate is is_markov_exact:
        assert result is False
    else:
        # The two joint eigenvalues each occur twice; marginals are tracial.
        expected = .75*np.log(1.5)+.25*np.log(.5)
        assert result == pytest.approx(expected, abs=1e-15)
