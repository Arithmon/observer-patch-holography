"""Independent and adversarial finite-source checks on bounded chains."""
import copy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg as la
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
MODULE_SPEC = importlib.util.spec_from_file_location("native_finite_source", HERE/"run.py")
n = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(n)


@pytest.fixture
def chain_data():
    states, q, lap = n.chain(5, 2, [(0,1,.5),(1,2,2.),(2,3,1.),(3,4,3.),(4,0,1.)])
    return states, q, lap, states-.4, .3


def test_response_matches_direct_double_time_integral():
    for u in [0, 1e-9, 1e-3, .5, 30, 300]:
        actual = 2*quad(lambda x:(1-x)*np.exp(-u*x), 0, 1, epsabs=1e-13)[0]
        assert float(n.average_response(u)) == pytest.approx(actual, abs=1e-12)


def test_absolute_variance_is_never_normalized_away(chain_data):
    states, q, lap, z, kappa = chain_data
    cov = n.graph_covariance(lap, kappa, 1e-9)
    np.testing.assert_allclose(cov, z.T@z/len(z), atol=1e-9)
    assert cov[0,0] == pytest.approx(.4*.6, abs=1e-9)
    assert n.graph_covariance(lap, kappa, 2000)[0,0] < .001


def test_independent_moment_equation_matches_general_readout(chain_data):
    states, q, lap, z, kappa = chain_data
    f = states[:,0]*states[:,1]+.3*z[:,3]
    f -= f.mean()
    eigen, vec = la.eigh(-q)
    for window in [.05, .5, 12.]:
        actual = n.moment_variance(q, f, window)
        predicted = n.spectral_covariance(eigen, vec, f, window)[0,0]
        assert actual == pytest.approx(predicted, abs=1e-11)


def test_finite_nonlinear_residual_is_psd_and_no_cross_term(chain_data):
    states, q, lap, z, kappa = chain_data
    f = np.column_stack([states[:,0]*states[:,1], states[:,2]*states[:,3]]).astype(float)
    f -= f.mean(axis=0)
    b, residual = n.density_split(f, z, kappa)
    eigen, vec = la.eigh(-q)
    for window in [.001, .2, 3, 100]:
        total = n.spectral_covariance(eigen, vec, f, window)
        rest = n.spectral_covariance(eigen, vec, residual, window)
        density = b@n.graph_covariance(lap, kappa, window)@b.T
        np.testing.assert_allclose(total, rest+density, atol=1e-11)
        assert la.eigvalsh(rest)[0] > 0
        assert np.linalg.norm(total-density) > 1e-5  # Reject silently discarding g.


def test_unit_rescaling_cancels_but_memory_change_changes_power(chain_data):
    states, q, lap, z, kappa = chain_data
    original = n.graph_covariance(lap, kappa, 4., rate=2.)
    changed_units = n.graph_covariance(lap, kappa, 4.*37, rate=2./37)
    np.testing.assert_allclose(original, changed_units, atol=1e-14)
    longer_memory = n.graph_covariance(lap, kappa, 8., rate=2.)
    assert np.trace(longer_memory) < .6*np.trace(original)


def test_clock_and_gain_cannot_be_renamed_curvature_selection(chain_data):
    states, q, lap, z, kappa = chain_data
    eigen, vec = la.eigh(-q)
    base = n.spectral_covariance(eigen, vec, z, 2)
    changed_gain = n.spectral_covariance(eigen, vec, 3*z, 2)
    np.testing.assert_allclose(changed_gain, 9*base, atol=1e-12)
    # Replacing T^-1 with T^-1/2 changes variance by a factor T.
    assert not np.allclose(base, 2*base)


def test_equal_area_coefficients_retain_quadrature_factor(chain_data):
    states, q, lap, z, kappa = chain_data
    eigen, vec = la.eigh(lap)
    mass = 4*np.pi/len(lap)
    coefficients = vec[:,1:].T@z[0]
    continuum_coefficients = np.sqrt(mass)*coefficients
    assert continuum_coefficients@continuum_coefficients == pytest.approx(mass*(z[0]@z[0]))
    assert (eigen[1:]/mass)@(continuum_coefficients**2) == pytest.approx(z[0]@lap@z[0])
    assert not np.isclose(coefficients@coefficients, continuum_coefficients@continuum_coefficients)


def test_ordinary_average_cannot_preserve_continuum_power_by_clock_choice():
    # Every duration/rate enters only through u>=0, and h(u)<=1.
    for u in np.concatenate([[0.], np.geomspace(1e-12, 1e12, 101)]):
        assert 0 < float(n.average_response(u)) <= 1
    rows = json.loads((HERE/"receipt.json").read_text())["fixed_area_refinement_bound"]
    bounds = [x["all_windows_all_rates_mode_variance_upper_bound"] for x in rows]
    assert all(a > b for a,b in zip(bounds,bounds[1:]))
    for row in rows:
        assert row["area_weight"]*row["kappa"] == pytest.approx(np.pi/(row["ports"]-1))


def test_receipt_mutations_fail():
    stored = json.loads((HERE/"receipt.json").read_text())
    n.validate(stored)
    bad = copy.deepcopy(stored)
    bad["native_carrier"]["windows"][0]["readback_covariance"][0][0] *= 2
    with pytest.raises(AssertionError):
        n.validate(bad)
    bad = copy.deepcopy(stored)
    bad["native_carrier"]["windows"][0]["mode_variances_native"][0] *= 2
    with pytest.raises(AssertionError):
        n.compare_reproduction(bad, stored)
    bad = copy.deepcopy(stored); bad["observational_data_read"] = True
    with pytest.raises(AssertionError):
        n.validate(bad)


def test_harmonic_bound_needs_no_exact_quadrature_or_distinct_directions():
    directions = np.array([[1,2,3],[-4,2,.1],[1,0,0],[1,0,0],[1,0,0]], dtype=float)
    for row in n.harmonic_coefficient_bounds(directions, [0,1,2,5,9], .3, gain_l1=2.):
        assert row["addition_theorem_absolute_error"] < 1e-13
        assert row["constant_projection_norm_increase"] < 1e-13
        assert row["universal_pseudo_Cl_upper_bound"] == pytest.approx(.3*4*4*np.pi/5)
        assert row["projected_pseudo_Cl_upper_bound"] <= row["universal_pseudo_Cl_upper_bound"]+1e-13


def test_invalid_inputs_fail(chain_data):
    states, q, lap, z, kappa = chain_data
    for value in [-1., np.nan, np.inf]:
        with pytest.raises(ValueError):
            n.average_response(value)
    with pytest.raises(ValueError):
        n.chain(3, 1, [(0,1,1.)])
    with pytest.raises(ValueError):
        n.graph_covariance(lap, kappa, 0)
    with pytest.raises(ValueError):
        n.spectral_covariance(*la.eigh(-q), states, 1.)
