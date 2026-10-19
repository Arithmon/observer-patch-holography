"""Native oracle, dynamic elimination, clock units and coherent-wave controls."""
from fractions import Fraction as F
import copy
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import eigh
import sympy as sp

from model.observer.scalar_response.experiment import compare_receipt
from . import spatial_dynamics as d


@pytest.mark.parametrize('attempts', (0, 1, 2, 4, 8))
def test_exact_native_configuration_oracle_preserves_every_field_pair(attempts):
    covariance = d.covariance(attempts)
    assert covariance.shape == (12, 12)
    assert covariance == covariance.T
    assert covariance*sp.ones(12, 1) == sp.zeros(12, 1)
    assert covariance[0, 0] > 0
    if attempts:
        assert covariance[0, 0] < d.covariance(0)[0, 0]
    assert d.native_density()['branches'] == 3960


def test_readout_is_canonical_exclusion_not_independent_site_preparation():
    covariance = d.covariance(0)
    assert covariance[0, 0] == sp.Rational(5, 36)
    assert covariance[0, 1] == -sp.Rational(5, 396)
    assert covariance != sp.Rational(5, 36)*sp.eye(12)


@pytest.mark.parametrize('z', (F(1, 2), 1, 4))
def test_native_discounted_correlations_are_exact_resolvents(z):
    result = d.native_discount(z)
    source = d.native_density()
    assert (source['laplacian']+z*sp.eye(12))*result == source['chi']*source['center']
    # A wrong discount placement retains a measurable normalization error.
    assert result*(1+F(z, 60)) != result


def test_finite_discrete_discount_tail_is_charged_including_the_constant_mode():
    source = d.native_density()
    matrix = np.eye(12)-np.asarray(source['laplacian'], dtype=float)/60
    z, terms = F(3, 2), 7
    discount = F(60, 60+z)
    partial = sum(float(discount**(n+1)/60)*np.linalg.matrix_power(matrix, n) for n in range(terms))
    full = np.linalg.inv(np.asarray(source['laplacian'], dtype=float)+float(z)*np.eye(12))
    error = np.linalg.norm(full-partial, 2)
    assert abs(error-float(discount**terms/z)) < 2e-15


@pytest.mark.parametrize('z', (F(1, 2), 1, 4))
def test_actual240port_dynamic_Schur_identity_and_its_exact_residual(z):
    value = d.base_resolvent(z)
    assert value['row_bound'] < F(1, 8)
    difference = value['difference']
    # Exact positivity on the constant-free principal coordinate chart.
    _, pivots = difference[:19, :19].LDLdecomposition(hermitian=False)
    assert all(pivots[i, i] > 0 for i in range(19))
    assert difference*sp.ones(20, 1) == sp.zeros(20, 1)
    assert value['compressed']*sp.ones(20, 1) == sp.ones(20, 1)/sp.Rational(z)


def test_independent_full_port_solve_matches_exact_compression():
    graph, lap, _, _, _ = d.source_green.matrices(0)
    v = np.repeat(np.eye(20)/np.sqrt(12), 12, axis=0)
    expected = v.T@np.linalg.solve(lap.toarray()+np.eye(graph.ports), v)
    np.testing.assert_allclose(expected, np.asarray(d.base_resolvent()['compressed'], dtype=float),
                               rtol=2e-14, atol=2e-15)


@pytest.mark.parametrize('permutation', ((0, 1, 2), (11, 7, 3)))
def test_bound_does_not_depend_on_geometric_choice_of_cross_ports(permutation):
    graph = d.source_green.source.federation(0)
    local = np.asarray(d.native_density()['laplacian'], dtype=float)
    full = np.kron(np.eye(20), local)
    used = [0]*20
    for left, right in graph.inter_seams:
        a, b = left//12, right//12
        u, v = 12*a+permutation[used[a]], 12*b+permutation[used[b]]
        used[a] += 1; used[b] += 1
        full[u, u] += 1; full[v, v] += 1; full[u, v] -= 1; full[v, u] -= 1
    u = np.repeat(np.eye(20)/np.sqrt(12), 12, axis=0)
    # Independently eliminate using a complete orthogonal within-cell basis.
    _, basis = np.linalg.eigh(np.eye(12)-np.ones((12, 12))/12)
    q = np.kron(np.eye(20), basis[:, 1:])
    a, c, block = u.T@full@u, u.T@full@q, q.T@full@q
    schur = a-c@np.linalg.solve(block, c.T)
    assert np.linalg.eigvalsh(block).min() > 2
    for z in (.25, 1, 4):
        effective = np.linalg.inv(schur+z*np.eye(20))
        actual = u.T@np.linalg.solve(full+z*np.eye(240), u)
        difference = effective-actual
        assert np.linalg.eigvalsh(difference).min() > -2e-13
        assert np.linalg.norm(difference, 2) < 1/8


def test_cell_mean_reset_is_not_the_actual_native_two_time_law():
    value = d.memory_control()
    assert value['two_attempt_compressed_transfer_defect'] == F(13, 19051200)
    assert value['defect_trace_before_attempt_normalization'] == F(25, 3)
    assert not value['cell_average_native_process_is_closed_Markov_chain']


@pytest.mark.parametrize('level', (0, 1, 3, 10, 100))
def test_all_level_dynamic_bounds_are_exact_and_use_distinct_clock_counts(level):
    value = d.bounds(level)
    assert value['dynamic_Schur_resolvent_error'] == F(1, 8*4**level)
    assert value['native_attempt_normalization'] == 1260*16**level
    assert value['root_incident_normalization'] == 66*4**level
    assert value['scaled_attempt_log_error_on_band'] < value['scaled_incident_ordered_log_eigenvalue_error_on_band']
    assert d.bounds(level+1 if level < 100 else level)['dynamic_Schur_resolvent_error'] <= value['dynamic_Schur_resolvent_error']


def test_log_generator_band_error_does_not_claim_global_linearization():
    for level in range(4):
        value = d.bounds(level, band=16)
        scale, denominator = 4**level, value['native_attempt_normalization']
        for a in (1, 4, 16):
            error = -denominator*np.log1p(-a/denominator)-a
            assert 0 <= error <= float(value['scaled_attempt_log_error_on_band'])+2e-14
    # The finite-density Laplacian's high band is outside a fixed scaled band.
    with pytest.raises(ValueError):
        d.bounds(0, band=66)


def test_ordered_incident_eigenvalue_bound_on_actual_native_density_module():
    value = d.incident_spectrum()
    laplacian = np.asarray(d.native_density()['laplacian'], dtype=float)
    density = np.linalg.eigvalsh(laplacian)
    decay = eigh(laplacian, np.asarray(value['D'], dtype=float), eigvals_only=True)
    a = 2*value['incident_seams']
    assert np.all(decay >= density/(a+density)-2e-15)
    assert np.all(decay <= density/a+2e-15)
    # The actual local-clock kernel is nonsymmetric: no scalar holding model.
    assert value['kernel'] != value['kernel'].T


def test_three_native_bands_cannot_share_one_relativistic_mass_and_speed():
    value = d.dispersion()
    y = value['squared_frequency_over_first_band']
    root = np.sqrt(5.)
    eigen = np.array([5-root, 6, 5+root])
    frequencies = -np.log1p(-eigen/60)
    numeric = (frequencies/frequencies[0])**2
    for actual, (lower, upper) in zip(numeric, y, strict=True):
        assert float(lower) <= actual <= float(upper)
    k = (eigen[2]-eigen[0])/(eigen[1]-eigen[0])
    defect = numeric[2]-1-k*(numeric[1]-1)
    assert float(value['affine_fit_defect_interval'][0]) <= defect <= float(value['affine_fit_defect_interval'][1])
    assert value['minimum_maximum_normalized_squared_frequency_error'] > F(3, 10)
    # A genuine massive relativistic comparison has exactly affine squares.
    for mass_squared, speed_squared in ((0, 1), (3, 2)):
        target = mass_squared+speed_squared*eigen
        target /= target[0]
        assert abs(target[2]-1-k*(target[1]-1)) < 1e-14


def test_native_coherent_wave_executes_amplitudes_and_both_Ramsey_outcomes():
    value = d.coherent_field()
    assert value['native_arc_dimension'] == 3960
    for row in value['all_return_and_complementary_Ramsey_outcomes']:
        assert row['Ramsey_plus']+row['Ramsey_minus'] == 1
        assert 0 <= row['Ramsey_plus'] <= 1
    assert value['all_return_and_complementary_Ramsey_outcomes'][2]['Ramsey_plus'] < 1
    assert value['same_first_event_opposite_reflection_two_step_Ramsey_plus'] == 1
    assert value['same_committed_instrument_coherent_trace_cube_difference'] == -F(1, 900)


@pytest.mark.parametrize('level', (0, 1, 2, 4))
def test_same_spatial_band_has_a_controlled_positive_coherent_wave_limit(level):
    value = d.wave_bounds(level, band=16, parameter=4)
    normalization = value['squared_coherent_clock_normalization']
    for a in (.25, 1, 4, 16):
        # arcsin is numerically stable as level grows; arccos(1-epsilon) is not.
        theta = 2*np.arcsin(np.sqrt(a/(4*normalization)))
        omega = np.sqrt(normalization)*theta
        assert 0 <= omega/np.sqrt(a)-1 <= float(value['relative_wave_frequency_error'])+2e-15
        steps = int(4*np.sqrt(normalization))
        t = steps/np.sqrt(normalization)
        error = abs(np.cos(steps*theta)-np.cos(t*np.sqrt(a)))
        assert error <= float(value['cosine_propagator_error_on_band'])+2e-14
    assert d.wave_bounds(level+1)['relative_wave_frequency_error'] < d.wave_bounds(level)['relative_wave_frequency_error']


@pytest.mark.parametrize('value', (True, -1, 1., '1'))
def test_ambiguous_source_level_or_rational_parameter_rejected(value):
    with pytest.raises(ValueError):
        d.bounds(value)
    with pytest.raises(ValueError):
        d.wave_bounds(value)
    with pytest.raises(ValueError):
        d.native_discount(value)


def test_native_matrices_are_not_mutable_shared_certificates():
    with pytest.raises(TypeError):
        d.native_density()['chi'] = 1
    with pytest.raises(TypeError):
        d.base_resolvent()['compressed'][0, 0] = 1


def test_retained_receipt_replays_and_rejects_promoted_scope():
    path = Path(__file__).with_name('spatial_dynamics_results.json')
    result = d.run()
    compare_receipt(json.loads(path.read_text()), result)
    changed = copy.deepcopy(result)
    changed['scope']['physical_relativistic_QFT_or_CFT_identified'] = True
    with pytest.raises(ValueError):
        compare_receipt(changed, result)


def test_fresh_runtime_import_closure_is_pinned():
    command = '''
import sys
from pathlib import Path
from model.observer_time_fields import spatial_dynamics as d
d.run()
missing=[]
for module in tuple(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if path is None:continue
    try:name=Path(path).resolve().relative_to(d.ROOT).as_posix()
    except ValueError:continue
    if name.endswith('.py') and name not in d.SOURCE_PATHS:missing.append(name)
assert not missing,missing
'''
    subprocess.run([sys.executable, '-c', command], cwd=d.ROOT, check=True)
