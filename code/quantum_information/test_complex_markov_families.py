"""Complex conditional products and correlated neighbours for review #1049."""
import numpy as np
import pytest

from quantum_information import conditional_mutual_information, is_markov_exact
from collar_alignment.msa_characterizations import collar_cmi


@pytest.mark.parametrize("seed", range(24))
def test_complex_conditional_product_and_its_nonmarkov_neighbour(seed):
    rng = np.random.default_rng(seed)
    pauli = [np.array([[0, 1], [1, 0]], complex),
             np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.])]
    left = [np.eye(2)/2 + sum(int(k)*p for k, p in zip(
        rng.integers(-2, 3, 3), pauli))/16 for _ in range(2)]
    right = [np.eye(2)/2 + sum(int(k)*p for k, p in zip(
        rng.integers(-2, 3, 3), pauli))/16 for _ in range(2)]
    rho = np.zeros((8, 8), complex)
    for b, mass in enumerate([.25, .75]):
        ix = [4*a+2*b+c for a in range(2) for c in range(2)]
        rho[np.ix_(ix, ix)] = mass*np.kron(left[b], right[b])
    # Every entry and product is exactly dyadic. B retains the sector label;
    # conditional independence is certified by construction, independently
    # of either the entropy calculation or the relative-modular moments.
    assert np.linalg.eigvalsh(rho).min() > .015
    assert is_markov_exact(rho, [2, 2, 2], [0], [1], [2])
    assert conditional_mutual_information(rho, [2, 2, 2], [0], [1], [2]) == 0
    assert collar_cmi([(1., rho, (2, 2, 1, 2))]) == 0

    # A Z_A Z_C correlation in the B=0 block leaves both marginals fixed.
    # Their product is unique, so the changed state is strictly non-Markov.
    # The perturbation is smaller than the analytic positivity margin.
    for i, sign in zip([0, 1, 4, 5], [1, -1, -1, 1]):
        rho[i, i] += sign/1024
    assert np.linalg.eigvalsh(rho).min() > .014
    assert not is_markov_exact(rho, [2, 2, 2], [0], [1], [2])
    value = conditional_mutual_information(rho, [2, 2, 2], [0], [1], [2])
    assert value > 0
    assert collar_cmi([(1., rho, (2, 2, 1, 2))]) == value
