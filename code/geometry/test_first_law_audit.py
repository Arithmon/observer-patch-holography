"""Regressions for step-size false positives in the entropy receipts (#1033)."""

import numpy as np
import pytest

from geometry.einstein_closure_receipts import (
    central_z, first_law_receipt, maxent_multiplier_receipt,
)


def test_reversing_a_variation_cannot_make_normalization_errors_negative():
    report = first_law_receipt(z_weights=[0., 0.], eps=-1e-6)
    for name in ["first_law_defect", "edge_identification_defect",
                 "bulk_identity_defect", "predicted_edge_defect"]:
        assert report[name] >= 0
    assert report["edge_identification_defect"] == pytest.approx(np.log(3/2), abs=1e-12)


def test_a_rounded_away_variation_cannot_certify_wrong_edge_normalization():
    with pytest.raises(ValueError, match="resolved"):
        first_law_receipt(z_weights=[0., 0.], eps=1e-20)


def test_wrong_normalization_has_a_finite_split_error_not_nan():
    report = first_law_receipt(z_weights=[0., 0.])
    assert np.isfinite(report["split_identity_defect"])
    assert report["split_identity_defect"] == pytest.approx(np.log(3/2), abs=1e-12)


def test_small_multiplier_is_not_replaced_by_entropy_subtraction_noise():
    report = maxent_multiplier_receipt(lam=1e-20)
    assert report["ds_dt"] == pytest.approx(1e-20, rel=1e-13, abs=0)


@pytest.mark.parametrize("step", [0., 1e-20])
def test_unresolved_gibbs_secant_is_an_explicit_input_error(step):
    with pytest.raises(ValueError, match="step|resolved"):
        maxent_multiplier_receipt(lam=1.3, dlam=step)


@pytest.mark.parametrize("dims,weights", [([True, 1], [0., 0.]),
                                        ([1, 1], [0., np.nan])])
def test_central_generator_requires_dimensions_and_finite_real_weights(dims, weights):
    with pytest.raises(ValueError):
        central_z(dims, [1, 1], weights)


@pytest.mark.parametrize("step", [1e-2, 1e-6, -1e-6])
def test_fixed_weights_do_not_certify_the_unsampled_center(step):
    report = first_law_receipt(z_weights=[0., 0.], move_weights=False, eps=step)
    assert report["edge_identification_defect"] == 0
    assert report["all_sector_transfer_defect"] == pytest.approx(np.log(3/2), abs=1e-15)
    assert report["sector_transfer_witness"] == [1., -1.]


@pytest.mark.parametrize("step", [1e-2, 1e-6, -1e-6])
def test_tangent_law_is_separate_from_a_nonzero_finite_remainder(step):
    report = first_law_receipt(eps=step)
    for field in ["first_law_defect", "bulk_identity_defect", "split_identity_defect"]:
        assert report[field] < 1e-13
    assert report["finite_relative_entropy"] > 0
    assert report["finite_remainder_decomposition_defect"] < 3e-15


def test_declared_diagonal_fixture_has_the_correct_finite_entropy_change(monkeypatch):
    import mpmath as mp
    import geometry.einstein_closure_receipts as receipts

    initial = [np.diag([.75, .25]), np.diag([.5, .25, .25])]
    target = [np.diag([.25, .75]), np.diag([.25, .5, .25])]
    inputs = iter(initial+target)
    monkeypatch.setattr(receipts, "random_faithful", lambda *_: next(inputs))
    step = .125
    report = receipts.first_law_receipt(eps=step)
    # Ordinary scalar Shannon entropies of all repeated full-state eigenvalues.
    def scalar_entropy(ps, states):
        result = mp.mpf(0)
        for p, state, edge in zip(ps, states, [2, 3]):
            for v in np.diag(state):
                eigenvalue = mp.mpf(float(p*v/edge))
                result -= edge*eigenvalue*mp.log(eigenvalue)
        return result
    with mp.workdps(80):
        varied = [s+step*(t-s) for s, t in zip(initial, target)]
        expected = float(scalar_entropy([.6+step, .4-step], varied)
                         - scalar_entropy([.6, .4], initial))
    assert report["finite_entropy_change"] == pytest.approx(expected, abs=2e-15)
    assert abs(report["finite_modular_change"]-expected) > .01


def test_gibbs_secant_and_analytic_slope_are_not_conflated():
    report = maxent_multiplier_receipt(lam=1.3, dlam=.1)
    assert report["ds_dt"] == pytest.approx(1.3, abs=1e-15)
    assert report["secant_multiplier_defect"] > 1e-3
