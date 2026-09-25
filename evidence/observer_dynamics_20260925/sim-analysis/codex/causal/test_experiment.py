import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location("causal_experiment", Path(__file__).with_name("experiment.py"))
E = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(E)


def test_spherical_caps_against_independent_cross_sections():
    for r, s in [(1, 1), (.3, 1), (.8, .4), (.2, .7)]:
        for d in [0, .1, .4, .9, 1.1, 2.0]:
            assert E.lens_volume(r, s, d) == pytest.approx(E.lens_volume_cross_section(r, s, d), rel=2e-12, abs=2e-13)


def test_continuum_lorentz_volume_and_temporal_parity_control():
    for k in [2, 3, 4, 5, 6, 8, 12]:
        delta = 1 / 7
        for beta in [0, .2, .5, .75, .9]:
            row = E.geometry_control(k, delta, beta * k * delta)
            assert row["calibration_relative_error"] < 1e-10
            assert row["proper_time_ratio"] == pytest.approx(math.sqrt(1 - beta ** 2))
        rest = E.geometry_control(k, delta, 0)
        b = 1 + 4 / k ** 2 if k % 2 == 0 else (1 - 1 / k ** 2) ** 2
        assert rest["finite_layer_over_continuum_volume"] == pytest.approx(b, abs=1e-13)


def test_actual_read_history_retains_wait_events_and_distinct_layers():
    # Three spatial sites with reads 0<->1<->2 and waiting at every site.
    ptr = np.array([0, 2, 5, 7])
    indices = np.array([0, 1, 0, 1, 2, 1, 2])
    to_last = np.array([2, 1, 0])
    assert E.longest_read_history(ptr, indices, 0, 4).tolist() == [4, 4, 4]
    history = E.witnessed_history(ptr, indices, 0, 2, to_last, 4)
    assert [e["site"] for e in history] == [0, 0, 0, 1, 2]
    assert len({(e["layer"], e["site"]) for e in history}) == 5
    assert E.witnessed_history(ptr, indices, 0, 2, to_last, 1) == []


def test_mutation_without_waiting_breaks_rank_prediction():
    # Removing waits genuinely changes the law. A bipartite graph cannot
    # return to its start after an odd number of one-layer reads.
    ptr = np.array([0, 1, 2])
    indices = np.array([1, 0])
    scores = E.longest_read_history(ptr, indices, 0, 3)
    assert scores[0] < 0
    assert scores[1] == 3


def test_axis_control_changes_graph_but_keeps_population_and_waits():
    sites = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
    ptr = np.arange(0, 17, 4)
    indices = np.tile(np.arange(4), 4)
    aptr, aidx = E.axis_graph(ptr, indices, sites)
    assert aidx[aptr[0]:aptr[1]].tolist() == [0, 1, 2]
    assert all(u in aidx[aptr[u]:aptr[u + 1]] for u in range(4))


def test_counts_use_both_endpoint_cones_and_keep_endpoints():
    dx, dy = np.array([0, 1, 2, -1]), np.array([2, 1, 0, -1])
    assert E.interval_counts(dx, dy, 2) == [1, 1, 1]
    assert E.interval_counts(dx, dy, 1) == [0, 0]
