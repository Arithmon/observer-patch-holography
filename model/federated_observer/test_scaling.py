"""Independent source, sector, spectral and propagation controls."""
from fractions import Fraction as F
from itertools import combinations, product
import math

import numpy as np
import pytest

from . import scaling as s


def test_local_identity_extends_to_every_number_sector_on_four_cycle():
    assert len(s.local_density_identity()) == 8
    edges = ((0, 1), (1, 2), (2, 3), (0, 3))
    weights = (F(2, 3), F(-7, 2), F(5), F(1, 4))
    lf = [F() for _ in weights]
    for a, b in edges:
        lf[a] += weights[a]-weights[b]
        lf[b] += weights[b]-weights[a]
    for number in range(5):
        for occupied in combinations(range(4), number):
            state = tuple(int(i in occupied) for i in range(4))
            value = sum(f*x for f, x in zip(weights, state, strict=True))
            actual = F()
            for edge in edges:
                for coin in (0, 1):
                    output, _ = s.native.native(state, edge, coin)
                    actual += value-sum(f*x for f, x in zip(weights, output, strict=True))
            assert actual == sum(f*x for f, x in zip(lf, state, strict=True))
    # This is the same normalization as NativeRepairFourCycle.coordinate_drift:
    # per-edge rateone has generator Q=−G/2, not−G.


def test_triangle_gram_is_full_binary_space_certificate():
    b = np.asarray(s.triangle_gram(), dtype=np.int64)
    np.testing.assert_array_equal(b.T@b, 3*b)
    a = (3*np.eye(8)+b)/2
    np.testing.assert_array_equal(a-1.5*np.eye(8), b.T@b/6)
    assert np.linalg.eigvalsh(a).min() >= 1.5-1e-14
    # Endpoint swap has a negative eigenvalue; positivity is in the sum,
    # not a false assertion that each swap is positive.
    swap = np.array([[0, 1], [1, 0]])
    assert np.linalg.eigvalsh(swap)[0] == -1


@pytest.mark.parametrize('level', [0, 1, 2])
def test_identical_triangle_cells_give_uniform_all_sector_log_bound(level):
    graph = s.source.federation(level)
    assert s.transfer_floor(level) == F(10, 21)
    assert len(graph.seams) == 63*graph.cells//2
    assert max(map(len, graph.neighbours)) == 6
    assert s.log_remainder(level, 0) == F(11, 10)
    for order in (4, 16):
        bound = s.log_remainder(level, order)
        d = F(11, 21)
        partial_tail = sum(d**k/k for k in range(order+1, order+70))
        assert partial_tail <= bound
        assert s.log_remainder(level, order+1) < bound


def test_sparse_source_spectrum_matches_dense_and_the_density_eigenrelation():
    graph = s.source.federation(0)
    lap = s.laplacian(0).toarray()
    np.testing.assert_array_equal(lap, lap.T)
    np.testing.assert_array_equal(lap@np.ones(graph.ports), np.zeros(graph.ports))
    values, vectors = np.linalg.eigh(lap)
    receipt = s.spectral_diagnostic(0)
    np.testing.assert_allclose(receipt['graph_laplacian_positive_eigenvalues_diagnostic'], values[1:9], atol=2e-11, rtol=0)
    # Any field eigenvector lifts to a fixed-number configuration eigenmode;
    # a nontrivial two-occupied preparation checks the exact source update.
    f = vectors[:, 1]
    state = np.zeros(graph.ports, dtype=int); state[[0, 17]] = 1
    mean = 0.
    for edge in graph.seams:
        for coin in (0, 1):
            output, _ = s.native.native(tuple(state.tolist()), edge, coin)
            mean += np.dot(f, output)/(2*len(graph.seams))
    assert abs(mean-(1-values[1]/(2*len(graph.seams)))*np.dot(f, state)) < 2e-14


@pytest.mark.parametrize('size', [8, 16, 32, 64, 128, 1024])
def test_cycle_exact_ratio_enclosures_exclude_gapless_linear_dispersion(size):
    row = s.cycle_diagnostic(size)
    assert F(row['ratio_lower']) > 2
    assert F(row['ratio_lower']) <= row['energy_ratio_diagnostic'] <= F(row['ratio_upper'])
    # Independent cycle matrix at the small declared controls.
    if size <= 32:
        lap = 2*np.eye(size)
        for i in range(size):
            lap[i, (i+1) % size] -= 1; lap[i, (i-1) % size] -= 1
        eig = np.linalg.eigvalsh(lap)
        e1, e2 = [-math.log1p(-v/(2*size)) for v in (eig[1], eig[3])]
        assert abs(e2/e1-row['energy_ratio_diagnostic']) < 5e-13
    assert s.cycle_bounds(2*size)[1]-s.cycle_bounds(2*size)[0] < F(row['ratio_upper'])-F(row['ratio_lower'])


def test_poisson_tail_enclosure_and_read_union_are_exact():
    time = F(1, 4)
    for radius in (1, 2, 4, 8):
        bound = s.locality_bound(6, time, radius)
        assert bound == min(F(1), F(3, 4)**radius/math.factorial(radius))
        # Independent direct Poisson summation is a diagnostic of the bound.
        tail = sum(math.exp(-.75)*.75**k/math.factorial(k) for k in range(radius, 100))
        assert tail <= float(bound)
        assert s.locality_bound(6, time, radius, reads=3) == min(F(1), 3*bound)
    assert s.locality_bound(6, F(0), 3) == 0
    assert s.locality_bound(0, F(9), 3) == 0
    assert s.locality_bound(6, F(100), 3) == 1


def test_poissonized_source_propagation_keeps_distant_nonzero_tails():
    row = s.propagation_diagnostic(0)
    assert row['fixed_word_dependency_ball_sizes']['0'] == 1
    assert row['fixed_word_dependency_ball_sizes']['1'] >= 6
    assert all(0 <= check['probability_at_distance_at_least_radius_diagnostic']
               <= float(F(check['uniform_exact_probability_bound'])) for check in row['tail_checks'])
    assert row['tail_checks'][0]['probability_at_distance_at_least_radius_diagnostic'] > 0
    assert not row['physical_light_cone_claimed']
    assert not row['native_CP_coherent_propagation_claimed']


def test_fixed_word_dependency_cone_is_not_a_zero_tail_continuous_clock():
    # Source rule on a path: an initially distant occupied endpoint cannot
    # influence vertexzero in fewer than its graphdistance edge operations.
    edges = ((0, 1), (1, 2), (2, 3))
    for word in product(edges, repeat=2):
        for coins in product((0, 1), repeat=2):
            state = (0, 0, 0, 1)
            for edge, coin in zip(word, coins, strict=True):
                state, _ = s.native.native(state, edge, coin)
            assert state[0] == 0
    state = (0, 0, 0, 1)
    for edge in reversed(edges):
        state, _ = s.native.native(state, edge, 1)
    assert state[0] == 1


@pytest.mark.parametrize('bad', [True, 0, 7, 65537, 8., [], None])
def test_cycle_domain_rejects_undeclared_settings(bad):
    with pytest.raises(ValueError):
        s.cycle_bounds(bad)


@pytest.mark.parametrize('args', [(True, F(1), 2), (6, 1., 2), (6, F(-1), 2),
                                 (6, F(1), True), (6, F(1), 0), (-1, F(1), 2)])
def test_tail_certificate_rejects_ambiguous_or_negative_inputs(args):
    with pytest.raises(ValueError):
        s.locality_bound(*args)


def test_log_and_origin_boundaries():
    for order in (True, -1, 129, 1.):
        with pytest.raises(ValueError):
            s.log_remainder(0, order)
    for origin in (True, -1, s.source.federation(0).ports):
        with pytest.raises(ValueError):
            s.distances(0, origin)
