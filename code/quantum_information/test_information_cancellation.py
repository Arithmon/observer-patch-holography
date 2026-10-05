"""Independent counterexamples to entropy-subtraction information scores."""

import mpmath as mp
import numpy as np
import pytest

from quantum_information import conditional_mutual_information, mutual_information
from collar_alignment.msa_characterizations import is_ec_aligned


def binary_information(delta, size):
    """Exact scalar spectrum: size/2 copies of 1/size +/- delta."""
    ctx = mp.mp.clone()
    ctx.dps = 750
    x = size * ctx.mpf(delta)
    return float(((1+x)*ctx.log1p(x)+(1-x)*ctx.log1p(-x))/2)


@pytest.mark.parametrize("power", (20, 30, 40, 50))
def test_small_classical_mutual_information_is_not_erased(power):
    delta = 2.**-power
    rho = np.diag([.25+delta, .25-delta, .25-delta, .25+delta])
    expected = binary_information(delta, 4)
    assert mutual_information(rho, [2, 2], [0], [1]) == pytest.approx(
        expected, rel=2e-12, abs=0)


@pytest.mark.parametrize("power", (20, 30, 40))
def test_classical_conditional_correlation_survives_entropy_cancellation(power):
    delta = 2.**-power
    parity = np.array([1, -1, -1, 1, -1, 1, 1, -1])
    rho = np.diag(.125 + delta*parity)
    expected = binary_information(delta, 8)
    assert conditional_mutual_information(rho, [2, 2, 2], [0], [1], [2]) == pytest.approx(
        expected, rel=2e-12, abs=0)


@pytest.mark.parametrize("delta", (1e-10, 1e-50, 1e-150))
def test_quantum_coherences_cannot_become_a_false_markov_state(delta):
    x = np.array([[0., 1.], [1., 0.]])
    y = np.array([[0., -1j], [1j, 0.]])
    rho = np.eye(8)/8 + delta*np.kron(np.kron(x, y), x)
    expected = binary_information(delta, 8)
    assert conditional_mutual_information(rho, [2, 2, 2], [0], [1], [2]) == pytest.approx(
        expected, rel=2e-12, abs=0)


@pytest.mark.parametrize("power", (30, 40))
def test_small_correlation_cannot_falsely_pass_collar_alignment(power):
    delta = 2.**-power
    rho = np.diag([.25+delta, .25-delta, .25-delta, .25+delta])
    assert not is_ec_aligned([(1., rho, (2, 1, 1, 2))], tol=delta*delta)


def test_unrepresentable_positive_information_is_not_reported_as_zero():
    x = np.array([[0., 1.], [1., 0.]])
    rho = np.eye(4)/4 + 1e-200*np.kron(x, x)
    with pytest.raises(ValueError, match="precision|underflow|represent"):
        mutual_information(rho, [2, 2], [0], [1])


def test_missing_density_entries_cannot_supply_an_information_score():
    rho = np.ma.array(np.eye(4)/4, mask=False)
    rho.mask[0, 1] = True
    with pytest.raises(ValueError, match="masked|missing"):
        mutual_information(rho, [2, 2], [0], [1])
