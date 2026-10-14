from fractions import Fraction as F
from math import factorial

import numpy as np
import pytest
from scipy.linalg import expm

from model.relativistic_scalar import bounds
from model.relativistic_scalar.lattice import ScalarLattice, ScalarParameters, free_dispersion


@pytest.mark.parametrize('d', [2, 4, 8])
@pytest.mark.parametrize('coupling', [0, 1])
def test_exact_onsite_oracle_and_directed_norm(d, coupling):
    lattice = ScalarLattice((d,)+(1,)*7, ScalarParameters(coupling=coupling))
    exact = bounds.local_centered_matrix(d, coupling)
    normalized = np.array([[float(exact[i][j])*np.sqrt(factorial(i)/factorial(j))
                            for j in range(d)] for i in range(d)])
    target = lattice.onsite(0)-(1+coupling/128)*np.eye(d)
    assert np.max(abs(normalized-target)) < 2e-13
    assert np.linalg.norm(target, 2) <= float(bounds.onsite_norm_upper(d, coupling))+2e-13


def test_strang_bound_against_independent_small_cell_exponential():
    lattice = ScalarLattice((2,)*8)
    onsite = sum((lattice.embed_local(lattice.onsite(i), [i])
                  for i in range(8)))
    onsite = onsite.toarray()
    h = lattice.hamiltonian.toarray()
    t, n = F(1, 16), 2
    half = expm(-0.5j*float(t)/n*onsite)
    product = np.linalg.matrix_power(half@expm(-1j*float(t)/n*(h-onsite))@half, n)
    assert np.linalg.norm(product-expm(-1j*float(t)*h), 2) <= float(bounds.product_error((2,)*8, t, n))


@pytest.mark.parametrize('a', [F(1), F(1, 2), F(1, 4)])
def test_free_window_controls_all_directions_and_phase_error(a):
    rng = np.random.default_rng(417)
    vectors = rng.normal(size=(100, 3))
    vectors *= rng.uniform(0, 1, 100)[:, None]/np.linalg.norm(vectors, axis=1)[:, None]
    continuum = np.sqrt(1+np.sum(vectors*vectors, axis=1))
    spatial = free_dispersion(vectors, spacing=float(a))
    bound = bounds.free_window(a, F(1), F(3))
    assert np.min(continuum-spatial) > -1e-14
    assert np.max(continuum-spatial) <= float(bound['frequency_error_upper'])
    assert np.max(abs(np.exp(-3j*continuum)-np.exp(-3j*spatial))) <= float(bound['free_evolution_vector_error_upper'])
    assert not bound['interaction_included']


def test_spatial_error_decreases_quadratically_without_claiming_interacting_limit():
    coarse = bounds.free_window(F(1), F(1), F(1))
    fine = bounds.free_window(F(1, 2), F(1), F(1))
    assert coarse['frequency_error_upper'] == 4*fine['frequency_error_upper']
    assert not fine['eight_site_boundary_included']
    gap = bounds.gap_statement()
    assert not gap['finite_occupation_gap_claimed']
    assert not gap['gap_measured_by_protocol']


@pytest.mark.parametrize('args', [(0, 1, 1), (1, 2, 1), (1., 1, 1), (1, 1, -1)])
def test_undeclared_spatial_domains_rejected(args):
    with pytest.raises(ValueError):
        bounds.free_window(*args)
