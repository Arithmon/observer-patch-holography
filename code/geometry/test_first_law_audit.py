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
