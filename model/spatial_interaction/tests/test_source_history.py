"""Source-word, independent evolution, and complete-record regression checks."""
from fractions import Fraction as Q
import math

import numpy as np
from scipy.linalg import expm

from model.spatial_interaction import bounds, experiment as e, source
from model.spatial_interaction.lattice import CubicBosonLattice, local_sector_leakage


def test_source_bond_on_all_local_inputs_and_reference_columns():
    gate, info = source.bond(.071)
    expected = expm(-.071j*source.bond_generator())
    rng = np.random.default_rng(419)
    reference = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    assert np.linalg.norm((gate-expected)@reference) < 3e-10
    assert local_sector_leakage(gate) < 2e-11
    assert info['resources']['proper_code_service_events'] == 2*info['resources']['native_pulses']


def test_compiled_pointer_has_both_outcomes_and_retains_backreaction():
    kraus, info = source.instrument('ramsey')
    detector = np.zeros(16, complex)
    detector[2], detector[5] = 1/math.sqrt(2), 1j/math.sqrt(2)
    projector = np.outer(detector, detector.conj())
    assert np.linalg.norm(kraus[1]-projector) < 2e-11
    assert np.linalg.norm(kraus[0]-(np.eye(16)-projector)) < 2e-11
    assert info['completeness_residual'] < 2e-12
    # The no outcome retains remote-particle configurations and the vacuum.
    assert abs(kraus[0][0, 0]-1) < 2e-11
    assert abs(kraus[1][0, 0]) < 2e-11
    assert np.linalg.norm(kraus[0]@kraus[1]) < 2e-11


def test_source_product_converges_to_independent_interacting_hamiltonian():
    lattice = CubicBosonLattice(1, .2, 1.3)
    initial = lattice.initial_state()
    expected = lattice.reference_evolve(initial, .3)
    coarse = e.evolve(lattice, initial, Q(3, 10), 4)
    fine = e.evolve(lattice, initial, Q(3, 10), 8)
    coarse_error, fine_error = (np.linalg.norm(v-expected) for v in (coarse, fine))
    assert 3.9 < coarse_error/fine_error < 4.1
    bound = bounds.strang_bound(8, Q(3, 10), Q(1, 5), Q(13, 10))
    assert fine_error < bound['operator_error_upper_float']


def test_full_feedback_history_compares_the_same_instruments_and_preparation():
    lattice = CubicBosonLattice(1, float(e.HOPPING), 1.)
    for phase in (0, 1):
        for connected in (False, True):
            actual = e.history(lattice, 8, phase, connected)
            reference = e.history(lattice, 8, phase, connected, source_execution=False)
            assert [r['bits'] for r in actual] == [[0, 0], [0, 1], [1, 0], [1, 1]]
            assert abs(sum(e.probabilities(actual))-1) < 2e-9
            assert e.total_variation(actual, reference) < 2e-7
    laws = [e.history(lattice, 8, phase) for phase in (0, 1)]
    assert e.total_variation(*laws) > .02  # actual current-controlled preparation
    records = e.native_records(laws)
    detached = e.native_records(laws, connected=False)
    assert len(records) == 66*2*4*2
    assert abs(math.fsum(r['probability'] for r in records)-1) < 2e-9
    assert any(a['records'] != b['records'] for a, b in zip(records, detached))
    assert all(r['probability'] >= 0 for r in records)


def test_matched_interaction_deletion_has_certified_positive_margin():
    row = e.witness(1, 8)
    assert row['strictly_positive_certified_margin']
    assert row['interaction_contrast'] > row['ideal_source_product_contrast_lower_float'] > .004


def test_boundary_comparison_and_spectral_threshold():
    small = CubicBosonLattice(1, float(e.HOPPING), 1.)
    large = CubicBosonLattice(2, float(e.HOPPING), 1.)
    law_small = e.history(small, 1, source_execution=False)
    law_large = e.history(large, 1, source_execution=False)
    budget = bounds.boundary_bound(1, 2*e.INTERVAL, e.HOPPING, e.INTERACTION)
    assert 0 < e.total_variation(law_small, law_large) < budget['complete_history_tv_upper_float']
    # Constant one-particle function has zero open-Laplacian kinetic energy;
    # this is the finite-volume upper half of the separately proved gap m.
    laplacian = np.diag(small.degree.astype(float))
    for x, y in small.edges:
        laplacian[x, y] = laplacian[y, x] = -1
    assert np.array_equal(laplacian@np.ones(small.volume), np.zeros(small.volume))
    assert np.linalg.eigvalsh(laplacian).min() > -1e-12
