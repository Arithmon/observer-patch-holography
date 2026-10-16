"""Compare the closed detector law with independent native path enumeration."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations

import pytest

from ..observer import native
from . import correlations as c


def graph(vertices):
    edges = tuple((i, i+1) for i in range(vertices-1))
    return edges, tuple(tuple(j for j in range(vertices)
                             if tuple(sorted((i,j))) in edges) for i in range(vertices))


def complete_paths(vertices, particles, steps):
    edges, _ = graph(vertices)
    prior = tuple(tuple(int(i in chosen) for i in range(vertices))
                  for chosen in combinations(range(vertices), particles))
    paths = {(x, x, ()): F(1, len(prior)) for x in prior}
    for _ in range(steps):
        following = defaultdict(F)
        for (initial, state, increments), weight in paths.items():
            for edge in edges:
                for coin in (0, 1):
                    out, record = native.native(state, edge, coin)
                    signed = record[1]*(int(edge[1] == 0)-int(edge[0] == 0))
                    following[initial, out, increments+(signed,)] += weight/F(2*len(edges))
        paths = following
    assert sum(paths.values()) == 1
    return paths


@pytest.mark.parametrize('vertices,particles', [(3,1),(3,2),(4,1),(4,2),(4,3),(5,2)])
def test_closed_field_covariance_matches_native_many_body_histories(vertices, particles):
    _, neighbours = graph(vertices)
    density = F(particles, vertices)
    for steps in (0,1,2,3):
        paths = complete_paths(vertices, particles, steps)
        for other in range(vertices):
            expected = sum(weight*(initial[0]-density)*(state[other]-density)
                           for (initial,state,_), weight in paths.items())
            assert c.covariance(neighbours, particles, 0, other, steps) == expected
        current_variance = sum(weight*sum(word)**2 for (_,_,word), weight in paths.items())
        assert c.increment_covariance(neighbours, particles, 0, 0,
                                      (0,steps), (0,steps)) == current_variance


def test_successive_current_windows_are_direct_record_reads_with_negative_covariance():
    _, neighbours = graph(4)
    paths = complete_paths(4, 2, 3)
    direct = sum(weight*word[0]*sum(word[1:]) for (_,_,word), weight in paths.items())
    assert direct < 0
    assert c.increment_covariance(neighbours, 2, 0, 0, (0,1), (1,3)) == direct
    assert all(sum(word) == state[0]-initial[0] for initial,state,word in paths)


def test_uniform_particle_conservation_and_empty_full_sectors():
    _, neighbours = graph(5)
    for particles in range(6):
        for gap in (0,1,8):
            assert sum(c.covariance(neighbours, particles, 0, j, gap)
                       for j in range(5)) == 0
    for particles in (0,5):
        assert c.increment_covariance(neighbours, particles, 0, 1, (0,8),(2,4)) == 0


def test_finite_native_attempt_propagation_has_exact_graph_support():
    _, neighbours = graph(5)
    for steps in range(4):
        row, denominator = c.one_point_counts(neighbours, 0, steps)
        assert sum(row) == denominator == 8**steps
        assert all(row[i] == 0 for i in range(steps+1,5))
        assert row[steps] == 1


@pytest.mark.parametrize('bad', [((1,),()), ((0,),(0,)), ((1,1),(0,)), ((True,),(0,)), ((),())])
def test_bad_graphs_fail_closed(bad):
    with pytest.raises(ValueError):
        c.one_point_counts(bad,0,1)


@pytest.mark.parametrize('start,steps', [(True,1),(0,True),(-1,1),(3,1),(0,-1),(0,65)])
def test_bad_detector_or_attempt_parameters_fail_closed(start,steps):
    with pytest.raises(ValueError):
        c.one_point_counts(graph(3)[1],start,steps)


def test_three_federation_levels_produce_exact_current_predictions():
    result = c.run()
    rows = result['integer_arithmetic_predictions']
    assert len(rows) == 24
    assert {row['elementary_observers'] for row in rows} == {20,80,320}
    assert all(F(row['net_current_variance']) > 0 for row in rows)
    assert all(F(row['opposite_endpoint_net_current_covariance']) < 0 for row in rows)
    assert all(F(row['cross_endpoint_covariance_over_local_variance']) == F(-1,6)
               for row in rows if row['uniform_seam_attempts'] == 1)
    assert not result['initial_occupation_access_needed_for_increment_reads']
