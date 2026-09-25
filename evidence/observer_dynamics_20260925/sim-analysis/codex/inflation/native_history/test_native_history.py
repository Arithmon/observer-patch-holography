import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg as la
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
module_spec = importlib.util.spec_from_file_location("codex_native_history_tested", HERE/"run.py")
native = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(native)


@pytest.mark.parametrize("u", [0, 1e-12, 1e-6, .001, .1, 1, 100])
def test_finite_window_formula_against_dimensionless_quadrature(u):
    window, variance = 3.7, .23
    expected = 2*variance*window*quad(lambda s:(1-s)*np.exp(-u*s),0,1,
                                   epsabs=1e-13,epsrel=1e-13)[0]
    assert native.occupation_variance(u/window,window,variance) == pytest.approx(expected,rel=1e-10)


def test_two_port_chain_independent_moment_evolution_and_discrete_clock():
    states,q = native.configuration_generator(2,1,[(0,1)])
    kappa = native.exact_checks(states,q,[(0,1)])
    q = np.array(q,dtype=float)
    f = np.array(states) @ np.array([1,-1])/np.sqrt(2)
    for window in [.1,1,10]:
        moments = native.occupation_moments(q,f,window)
        assert moments[1] == pytest.approx(0,abs=1e-12)
        assert moments[2] == pytest.approx(native.occupation_variance(1,window,float(kappa)),rel=1e-11)
    # At every discrete attempt this two-port state is resampled independently.
    # Poisson holding intervals change the limiting occupation variance by2.
    assert native.discrete_occupation_variance(0,100,1,float(kappa)) == .5
    assert native.occupation_variance(1,1e8,float(kappa)) == pytest.approx(1,rel=2e-8)


def test_disconnected_graph_exposes_missing_ergodicity_hypothesis():
    states,q = native.configuration_generator(4,2,[(0,1),(2,3)])
    q = np.array(q,dtype=float)
    f = np.array(states) @ np.array([.5,.5,-.5,-.5])
    assert np.allclose(q@f,0)
    variance = float(np.mean(f*f))
    for window in [1,10]:
        moments = native.occupation_moments(q,f,window,order=2)
        assert moments[2] == pytest.approx(variance*window)
    # Component-total randomness cannot be replaced by the connected-graph L+ result.
    assert variance > 0


def test_leaky_formula_against_independent_configuration_resolvent():
    edges=[(0,1),(1,2),(2,3)]
    states,q=native.configuration_generator(4,2,edges)
    kappa=float(native.exact_checks(states,q,edges))
    q=np.array(q,dtype=float)
    eigen,vectors=la.eigh(native.laplacian(4,edges).astype(float))
    for j in range(1,4):
        f=np.array(states)@vectors[:,j]
        for leak in [.03,.3,3.]:
            expected=2*f@la.solve(leak*np.eye(len(q))-q,f)/len(q)
            assert native.leaky_variance(eigen[j]/2,leak,kappa)==pytest.approx(expected,rel=1e-12)


def test_receipt_checks_independent_methods_and_non_gaussian_finite_windows():
    receipt=json.loads((HERE/"receipt.json").read_text())
    for name,expected in receipt["inputs_sha256"].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected
    for chain in receipt["chains"]:
        assert chain["exact_stationarity_drift_and_covariance_checks"]
        assert chain["max_semigroup_covariance_absolute_error"]<1e-11
        assert chain["clt_covariance_max_absolute_error"]<1e-11
        for row in chain["windows"]:
            assert row["variance_moment_ode"]==pytest.approx(row["variance_analytic"],rel=1e-10)
            assert row["variance_semigroup_quadrature"]==pytest.approx(row["variance_analytic"],rel=1e-10)
        assert chain["windows"][0]["kurtosis"]<2.6
        assert 2.95<chain["windows"][-1]["kurtosis"]<3
        for row in chain["discrete_controls"]:
            assert row["variance_chain"]==pytest.approx(row["variance_formula"],rel=1e-11)
        for row in chain["leak_controls"]:
            assert row["variance_configuration_resolvent"]==pytest.approx(row["variance_analytic"],rel=1e-11)


def test_actual_graph_predictions_preserve_geometry_and_finite_window_cutoff():
    receipt=json.loads((HERE/"receipt.json").read_text())
    previous=json.loads((HERE.parents[1]/"dynamics/receipt.json").read_text())
    # Retained graph source hashes make this more than a synthetic lattice control.
    assert receipt["production_source_commit"]=="14d1699"
    assert len(receipt["production_source_sha256"])==7
    for graph in receipt["actual_graphs"]:
        assert graph["max_eigenpair_relative_residual"]<1e-6
        assert graph["ports"]==12*graph["carriers"]
        assert sum(graph["port_attempt_rates_per_sweep"].values())==graph["ports"]
        for window in graph["windows"]:
            ratio=np.array(window["mode_variance_over_gff_limit"])
            assert ((ratio>0)&(ratio<1)).all()
            assert np.diff(ratio).min()>-1e-10
        assert graph["windows"][-1]["mode_variance_over_gff_limit"][0]==pytest.approx(.99)
    prior={row["level"]:row for row in previous["levels"]}
    for graph in receipt["actual_graphs"]:
        assert graph["seam_endpoints_sha256"]==prior[graph["level"]]["seam_endpoints_sha256"]
