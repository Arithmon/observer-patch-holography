import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import pytest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("codex_weighted_universality_tested",HERE/"certify.py")
model=importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)


def fractions(rows):
    return [[F(x) for x in row] for row in rows]


@pytest.fixture(scope="module")
def receipt():
    return json.loads((HERE/"receipt.json").read_text())


def test_retained_source_hashes(receipt):
    for name,expected in receipt["inputs_sha256"].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected
    _,pin=model.production_rule("14d1699")
    assert pin==receipt["production_rule_source"]


@pytest.mark.parametrize("name",["path4","weighted_path4","cycle5"])
def test_full_configuration_certificate_recomputed_exactly(name,receipt):
    design=json.loads((HERE/"spec.json").read_text())
    graph=next(g for g in design["graphs"] if g["name"]==name)
    rule,_=model.production_rule(design["source_commit"])
    states,q,lap,incidence,z,kappa,pi=model.configuration_data(
        graph["ports"],graph["raised"],[(a,b,F(w)) for a,b,w in graph["edges"]],rule)
    gplus=model.pseudoinverse_connected_laplacian(model.scale(F(-1),q))
    record=next(g for g in receipt["graphs"] if g["name"]==name)
    for data in record["records"].values():
        b=fractions(data["B"])
        observed=model.multiply(z,model.transpose(b))
        direct=model.scale(F(2,len(states)),model.congruence(model.transpose(observed),gplus))
        assert direct==fractions(data["continuous_occupation_limit"])


def test_hand_computed_path_load_and_native_seam_countercontrol(receipt):
    graph=next(g for g in receipt["graphs"] if g["name"]=="path4")
    expected=fractions([["7/6","1/6","-1/2","-5/6"],
                        ["1/6","1/2","-1/6","-1/2"],
                        ["-1/2","-1/6","1/2","1/6"],
                        ["-5/6","-1/2","1/6","7/6"]])
    assert fractions(graph["records"]["load"]["continuous_occupation_limit"])==expected
    assert fractions(graph["records"]["seam"]["continuous_occupation_limit"])==model.scale(F(4,3),model.eye(3))
    assert fractions(graph["records"]["drive"]["continuous_occupation_limit"])==model.scale(F(4,3),fractions(graph["laplacian"]))


def test_production_rule_check_rejects_biased_remainder_and_excludes_unbalanced_state():
    rule,_=model.production_rule("14d1699")
    biased=lambda a,b,ceiling_to_first:rule(a,b,ceiling_to_first=True)
    with pytest.raises(AssertionError):
        model.configuration_data(3,1,[(0,1,F(1)),(1,2,F(1))],biased)
    assert rule(2,4,ceiling_to_first=True)==(3,3)
    assert rule(2,4,ceiling_to_first=True)!=(4,2)


def test_successful_event_clock_has_different_exact_stationary_covariance(receipt):
    graph=next(g for g in receipt["graphs"] if g["name"]=="path4")
    jump=graph["successful_event_clock"]
    assert jump["stationary_configuration_weights"]==["1/12","1/4","1/6","1/6","1/4","1/12"]
    assert F(jump["stationary_load_covariance"][0][1])==F(-1,6)
    assert F(graph["records"]["load"]["stationary_covariance"][0][1])==F(-1,12)
    assert all(g["successful_event_clock"]["differs_from_attempt_clock_covariance"] for g in receipt["graphs"])


def test_discrete_clock_correction_is_preserved_for_all_records(receipt):
    for graph in receipt["graphs"]:
        for data in graph["records"].values():
            continuous=fractions(data["continuous_occupation_limit"])
            static=fractions(data["stationary_covariance"])
            expected=model.add(continuous,model.scale(-1/F(graph["total_attempt_rate"]),static))
            assert expected==fractions(data["discrete_attempt_occupation_limit"])


def test_finite_window_bounds_independent_quadrature_and_initial_state_control(receipt):
    for graph in receipt["graphs"]:
        for row in graph["finite_window_checks"]:
            assert row["quadrature_relative_matrix_error"]<1e-10
            assert row["lower_bound_min_eigenvalue"]>-1e-11
            assert row["upper_bound_min_eigenvalue"]>-1e-11
        for row in graph["nonstationary_initial_checks"]:
            assert row["mean_moment_ode"]==pytest.approx(row["mean_graph_formula"],rel=1e-10,abs=1e-12)
    path=next(g for g in receipt["graphs"] if g["name"]=="path4")
    window=next(r for r in path["nonstationary_initial_checks"] if r["window"]==1)
    assert window["variance_deterministic_initial"]<.2*window["variance_stationary_formula"]


def test_uniformly_comparable_rates_give_inverse_covariance_bounds(receipt):
    graphs={g["name"]:g for g in receipt["graphs"]}
    base=np.array(fractions(graphs["path4"]["laplacian_pseudoinverse"]),dtype=float)
    weighted=np.array(fractions(graphs["weighted_path4"]["laplacian_pseudoinverse"]),dtype=float)
    assert np.linalg.eigvalsh(weighted-base/4).min()>-1e-12
    assert np.linalg.eigvalsh(4*base-weighted).min()>-1e-12


def test_covariance_symmetry_requires_equivariant_record(receipt):
    cycle=next(g for g in receipt["graphs"] if g["name"]=="cycle5")
    sigma=fractions(cycle["records"]["load"]["continuous_occupation_limit"])
    permutation=model.zeros(5,5)
    for i in range(5):permutation[(i+1)%5][i]=F(1)
    assert model.congruence(permutation,sigma)==sigma
    gain=model.eye(5);gain[0][0]=F(2)
    changed=model.congruence(gain,sigma)
    assert model.congruence(permutation,changed)!=changed


def test_fixed_window_refinement_loses_gff_limit(receipt):
    rows=receipt["refinement_countercontrol"]
    assert all(a["fixed_window_gff_fraction"]>b["fixed_window_gff_fraction"] for a,b in zip(rows,rows[1:]))
    assert rows[0]["fixed_window_gff_fraction"]>.67
    assert rows[-1]["fixed_window_gff_fraction"]<.0016
    assert all(row["scaled_window_gff_fraction"]==pytest.approx(.9000045399929762) for row in rows)


def test_nonlinear_projection_and_poisson_certificate_independently():
    result=json.loads((HERE/"nonlinear_receipt.json").read_text())
    for name,expected in result["inputs_sha256"].items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected
    rule,_=model.production_rule("14d1699")
    states,q,lap,incidence,z,kappa,pi=model.configuration_data(
        4,2,[(0,1,F(1)),(1,2,F(1)),(2,3,F(1))],rule)
    gplus=model.pseudoinverse_connected_laplacian(model.scale(F(-1),q))
    for row in result["observables"].values():
        observable=fractions(row["observable"])
        b=model.scale(1/(len(states)*kappa),model.multiply(model.transpose(observable),z))
        assert b==fractions(row["B"])
        residual=model.add(observable,model.scale(F(-1),model.multiply(z,model.transpose(b))))
        assert model.multiply(model.transpose(residual),z)==model.zeros(1,4)
        direct=model.scale(F(2,len(states)),model.congruence(model.transpose(observable),gplus))[0][0]
        assert direct==F(row["total_asymptotic_variance"])
        density=model.scale(4*kappa,model.congruence(b,model.pseudoinverse_connected_laplacian(lap)))[0][0]
        rest=model.scale(F(2,len(states)),model.congruence(model.transpose(residual),gplus))[0][0]
        assert direct==density+rest
        assert rest>=0
    even=result["observables"]["even_pair"]
    assert fractions(even["B"])==model.zeros(1,4)
    assert F(even["total_asymptotic_variance"])==F(1,6)
    assert F(result["observables"]["mixed"]["total_asymptotic_variance"])==F(4,3)
