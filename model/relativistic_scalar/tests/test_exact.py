"""Independent checks of oscillator compression and exact certificate bounds."""
from fractions import Fraction as F
from math import factorial

import numpy as np
import pytest
from scipy import sparse

from model.relativistic_scalar import exact


def test_cube_hamiltonian_and_witness_derivative():
    assert len(exact.EDGES) == 12
    vector = exact.apply_hamiltonian(exact.vacuum())
    assert vector.get(exact.VACUUM, 0) == 0
    four = (4,) + (0,) * 7
    assert vector[four] == F(1, 384)
    # <e|H|vac> = i*4!*coefficient/5, so p'(0)=2*(1/5)*(4!*coefficient/5).
    derivative = F(2, 25) * factorial(4) * vector[four]
    assert derivative == F(1, 200)
    assert exact.norm_squared(exact.vacuum()) == 1
    assert exact.apply_hamiltonian(exact.vacuum(), center=False)[exact.VACUUM] == F(129, 16)


def test_quartic_compression_retains_paths_outside_local_cutoff():
    # Site0 in |1>, all other sites in |0>. Compare the quartic difference
    # against an independent oscillator matrix computed above the cutoff.
    n = 10
    annihilation = np.diag(np.sqrt(np.arange(1, n)), 1)
    b = annihilation + annihilation.T
    b4 = np.linalg.matrix_power(b, 4)
    key = (1,) + (0,) * 7
    on = exact.apply_hamiltonian({key: F(1)}, coupling=F(1), local_cutoff=2)
    off = exact.apply_hamiltonian({key: F(1)}, coupling=F(0), local_cutoff=2)
    # Ground expectations at every other site cancel under centering.
    assert on[key]-off[key] == F(1, 32)
    assert float(on[key]-off[key]) == pytest.approx((b4[1, 1]-3)/384)
    b_truncated = b[:2, :2]
    assert not np.isclose(b4[1, 1], np.linalg.matrix_power(b_truncated, 4)[1, 1])


def test_untruncated_moments_and_exact_coefficient_agreement():
    free = exact.vacuum_moments(F(0), 2)
    interacting = exact.vacuum_moments(F(1), 2)
    assert free[1]["norm_squared"] == F(3, 4)
    assert interacting[1]["norm_squared"] == F(145, 192)
    assert interacting[2]["maximum_local_occupation"] == 8
    assert all(row["coefficients_match"] for row in exact.coefficient_comparison(9, 2))
    comparison = exact.coefficient_comparison((8,) + (2,)*7, 2)
    assert comparison[0]["coefficients_match"]
    assert not comparison[1]["coefficients_match"]
    assert comparison[1]["difference_norm_squared"] > 0


def test_sparse_operator_is_hermitian_in_factorial_metric():
    first = (2, 1) + (0,)*6
    second = (3, 0) + (0,)*6
    hfirst = exact.apply_hamiltonian({first: F(1)})
    hsecond = exact.apply_hamiltonian({second: F(1)})
    lhs = hfirst.get(second, 0) * factorial(3)
    rhs = hsecond.get(first, 0) * factorial(2)*factorial(1)
    assert lhs == rhs != 0


def test_compressed_polynomial_matches_canonical_q_p_hamiltonian():
    dims = (5,) + (2,)*7
    full_dimension = exact.tensor_dimension(dims)
    q, p2, q2, q4 = [], [], [], []
    for d in dims:
        a = np.diag(np.sqrt(np.arange(1, d+4)), 1)
        position = (a+a.T)/2
        momentum = 1j*(a.T-a)
        q.append(sparse.csr_matrix(position[:d, :d]))
        p2.append(sparse.csr_matrix((momentum@momentum).real[:d, :d]))
        q2.append(sparse.csr_matrix((position@position)[:d, :d]))
        q4.append(sparse.csr_matrix(np.linalg.matrix_power(position, 4)[:d, :d]))

    def embed(operators):
        matrix = sparse.csr_matrix([[1.0]])
        for site, d in enumerate(dims):
            matrix = sparse.kron(matrix, operators.get(site, sparse.eye(d)), format="csr")
        return matrix

    # Canonical p²/2 + m²q²/2 + sum_edges(qx-qy)²/2 + sum q⁴/24.
    matrix = sparse.csr_matrix((full_dimension, full_dimension))
    for site in range(8):
        matrix += embed({site: p2[site]/2 + q2[site]/2 + q4[site]/24})
    for x, y in exact.EDGES:
        matrix += embed({x: q2[x]/2}) + embed({y: q2[y]/2}) - embed({x: q[x], y: q[y]})
    matrix -= float(F(129, 16)) * sparse.eye(full_dimension)
    state = {exact.VACUUM: F(1, 3), (4,)+(0,)*7: F(-1, 7), (2, 1)+(0,)*6: F(2, 5)}

    def normalized_amplitudes(vector):
        out = np.zeros(full_dimension)
        for key, coefficient in vector.items():
            out[np.ravel_multi_index(key, dims)] = float(coefficient)*np.sqrt(
                np.prod([factorial(n) for n in key]))
        return out

    expected = matrix @ normalized_amplitudes(state)
    actual = normalized_amplitudes(exact.apply_hamiltonian(state, local_cutoff=dims))
    np.testing.assert_allclose(actual, expected, atol=3e-14, rtol=2e-14)


def test_directed_roots_and_nonvacuous_short_time_certificates():
    for value in (F(0), F(1, 4), F(2), F(1234567, 999)):
        upper = exact.sqrt_upper(value, 12)
        assert upper*upper >= value
        if upper:
            assert (upper-F(1, 10**12))**2 < value or value == 0
    time = F(1, 100000)
    infinite = exact.interaction_witness_certificate(time)
    finite = exact.interaction_witness_certificate(time, (8,) + (2,)*7)
    assert infinite["strictly_positive"]
    assert finite["strictly_positive"]
    truncation = exact.truncation_certificate(time, (8,) + (4,)*7, order=2)
    assert 0 < truncation["readout_tv_upper"] < F(1, 100000)
    assert exact.tensor_dimension((8,) + (2,)*7) == 1024


def test_invalid_inputs_fail_closed():
    for bad in (True, 0, (2, 2), (2,)*7 + (False,), 2.0):
        with pytest.raises(ValueError):
            exact.cutoffs(bad)
    with pytest.raises(ValueError):
        exact.apply_hamiltonian({(5,)+(0,)*7: F(1)}, local_cutoff=4)
    with pytest.raises(ValueError):
        exact.interaction_witness_certificate(F(1, 1000), 4)
    with pytest.raises(ValueError):
        exact.sqrt_upper(F(-1))
    with pytest.raises(ValueError):
        exact.vacuum_moments(F(1, 2))
