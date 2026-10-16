"""Audit controls for what the existing finite records can identify.

Run from the repository root with ``python -m
model.identifiability_checks.architecture_counterexamples``.  These are independent
mathematical/reference calculations, not additional native source executions.
The occupation ambiguity is an exact identity; the continued-history contrast
is an ordinary floating sparse-Hamiltonian diagnostic, not an interval proof.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product
import json

import numpy as np

from model.interacting_reconstruction.prediction import RecoveredField
from model.observer_reconstruction.recovery import cube_geometry
from model.relativistic_scalar.lattice import CUBE_EDGES
from model.scalar_observations.lattice import GalerkinParameters


def occupation_ambiguity():
    """A local positive term invisible to the entire four-level calibration."""
    def extra(n):
        return F(n*(n-1)*(n-2)*(n-3), 24)

    diagonal = [extra(n) for n in range(8)]
    assert diagonal[:4] == [0]*4
    assert diagonal[4] == 1
    assert all(value >= 0 for value in diagonal)
    bonds = np.zeros((8, 8))
    for left, right in CUBE_EDGES:
        bonds[left, right] = bonds[right, left] = 1.
    dimensions = (8,)+(2,)*7
    base = RecoveredField(dimensions, [4.]*8, [1/24]*8, bonds)
    changed = RecoveredField(dimensions, [4.]*8, [1/24]*8, bonds)
    # The same nonnegative local addition is specified at every site.  At
    # this asymmetric withheld cutoff only the origin can reach n >= 4.
    for site, dimension in enumerate(dimensions):
        changed.hamiltonian += changed.embed({site: np.diag(
            [float(extra(n)) for n in range(dimension)])})
    contrasts = []
    for sign in (-1, 1):
        reference = base.history(1/128, sign=sign)
        alternative = changed.history(1/128, sign=sign)
        p = np.array([row['probability'] for row in reference])
        q = np.array([row['probability'] for row in alternative])
        assert abs(p.sum()-1) < 2e-12 and abs(q.sum()-1) < 2e-12
        tv = float(np.abs(p-q).sum()/2)
        assert tv > 1e-6
        contrasts.append({'pointer_sign': sign, 'base_probabilities': p.tolist(),
                          'alternative_probabilities': q.tolist(), 'tv': tv})
    return {
        'added_term_per_site': 'N(N-1)(N-2)(N-3)/24',
        'eight_level_diagonal_exact': [str(x) for x in diagonal],
        'four_level_calibration_operator_difference_exact': '0',
        'calibration_identity_scope': 'every input and duration in the complete four-level-per-site code',
        'withheld_dimensions': list(dimensions),
        'model_time_per_segment': 1/128,
        'continued_history_diagnostics': contrasts,
        'source_execution': False,
        'numerical_contrast_is_rigorous_bound': False,
    }


def embedding_ambiguity():
    """Same cube graph and unit edges, but eight distinct planar vertices."""
    directions = ((F(1), F(0)), (F(3, 5), F(4, 5)),
                  (F(-5, 13), F(12, 13)))
    labels = tuple(product((0, 1), repeat=3))
    planar = tuple(tuple(sum((bits[j]*directions[j][axis] for j in range(3)), F(0))
                          for axis in range(2)) for bits in labels)
    assert len(set(planar)) == 8
    squared_distance = lambda i, j: sum((planar[i][axis]-planar[j][axis])**2
                                      for axis in range(2))
    assert all(squared_distance(i, j) == 1 for i, j in CUBE_EDGES)
    # Two independent centered edge differences prove the affine rank is two.
    assert directions[0][0]*directions[1][1]-directions[0][1]*directions[1][0] != 0
    adjacency = np.zeros((8, 8))
    for left, right in CUBE_EDGES:
        adjacency[left, right] = adjacency[right, left] = 1
    recovered = cube_geometry(4*np.eye(8)-adjacency)
    assert recovered['spatial_dimension'] == 3
    mismatches = [abs(squared_distance(i, j)-sum(abs(x-y) for x, y in zip(labels[i], labels[j])))
                  for i in range(8) for j in range(8)]
    assert max(mismatches) > 0
    return {
        'planar_vertices_exact': [[str(x) for x in point] for point in planar],
        'all_twelve_graph_edge_squared_lengths_exact': '1',
        'planar_affine_rank_exact': 2,
        'canonical_cube_inverse_reported_rank': recovered['spatial_dimension'],
        'max_difference_from_graph_distance_exact': str(max(mismatches)),
        'hamiltonian_and_records_change': False,
        'scope': 'weighted graph dynamics do not alone select a Euclidean embedding',
    }


def fixed_band_grid_check():
    """Only one gradient coefficient changes in the retained seven-mode grids."""
    exact = GalerkinParameters(grid=None)
    rows = []
    for grid in (5, 9, 17, 33, 65):
        finite = GalerkinParameters(grid=grid)
        for mode in range(7):
            assert finite.onsite_quartic(mode) == exact.onsite_quartic(mode)
            for other in range(mode+1, 7):
                assert finite.pair_coefficient(mode, other) == exact.pair_coefficient(mode, other)
        rows.append({'grid': grid, 'field_mode_count': 7, 'kappa': finite.kappa,
                     'quadratic_frequency_squared_error': 1-finite.kappa})
    return {'exact_gradient_already_supported': True, 'grids': rows}


def run():
    return {'scope': 'architecture audit; no empirical outcomes or frozen receipts changed',
            'occupation_extension': occupation_ambiguity(),
            'geometry': embedding_ambiguity(),
            'fixed_band': fixed_band_grid_check()}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
