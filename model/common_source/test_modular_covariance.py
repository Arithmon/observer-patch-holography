"""Full-algebra bounds, compressed boundary controls and native witness."""
import math

import numpy as np
import pytest

from . import modular_covariance as covariance


def test_ward_identity_noncommuting_pair_and_sharp_commuting_control():
    g = np.array([[1, 2j], [-2j, -3]])
    p = np.array([[2, 1], [1, 1]])
    for rate in (-2., .25, 3.):
        report = covariance.ward_certificate(g, p, rate)
        assert report["squared_identity_residual"] < 3e-13
        assert report["global_relative_defect"] >= 1
        sharp = covariance.ward_certificate(g, np.eye(2), rate)
        assert sharp["global_relative_defect"] == pytest.approx(1.)


def test_finite_parameter_norm_bound_both_orientations():
    g = np.diag([0., 1., 4.])
    p = np.array([[1, .3j, .2], [-.3j, 2, .5], [.2, .5, 3]])
    for s in (-1., .25, 2.):
        result = covariance.finite_dilation_certificate(g, p, s, .5)
        assert result["operator_defect"]+1e-13 >= result["operator_defect_lower_bound"]
    identity = covariance.finite_dilation_certificate(g, np.eye(3), .75, .5)
    assert identity["operator_defect"] == pytest.approx(identity["operator_defect_lower_bound"])


def test_spectral_compression_preserves_obstruction():
    g = np.diag([0., 1., 4.])
    p = np.array([[1, .3j, .2], [-.3j, 2, .5], [.2, .5, 3]])
    result = covariance.compressed_certificate(g, p, np.diag([1., 1., 0.]), .5)
    assert result["boundary_term_hs_norm"] == 0
    assert result["decomposition_hs_error"] < 1e-13
    assert result["compressed_defect_hs_norm"] >= result["necessary_boundary_plus_defect_lower_bound"]


def test_native_derivative_match_requires_boundary_and_fails_finite_flow():
    report = covariance.native_certificate()
    assert report["dimension"] == 66
    compressed = report["compressed_ward"]
    assert compressed["compressed_defect_hs_norm"] < 2e-14
    assert compressed["internal_derivative_hs_norm"] < 2e-14
    assert compressed["boundary_term_hs_norm"] == pytest.approx(math.log(5)/2)
    assert compressed["generator_leakage_hs_norm"] == pytest.approx(math.log(5)/2)
    assert compressed["translation_leakage_hs_norm"] == pytest.approx(.5)
    assert report["global_ward"]["global_relative_defect"] > 1
    assert report["compressed_finite_parameter"]["difference"] < -.22
    assert report["physical_boost_or_null_momentum_identified"] is False
    assert report["continuum_no_go_claimed"] is False


@pytest.mark.parametrize("bad", [np.array([[1, 1], [0, 1]]),
                                 np.array([[np.nan]]), np.zeros((2, 3)),
                                 np.ma.array([[1]], mask=[[True]])])
def test_bad_generators_rejected(bad):
    with pytest.raises(ValueError):
        covariance.ward_certificate(bad, np.eye(2), .5)


def test_nonprojector_and_inconsistent_algebras_rejected():
    with pytest.raises(ValueError, match="orthogonal projector"):
        covariance.compressed_certificate(np.eye(2), np.eye(2), np.eye(2)/2, 1.)
    with pytest.raises(ValueError, match="one algebra"):
        covariance.ward_certificate(np.eye(2), np.eye(3), 1.)
    with pytest.raises(ValueError, match="finite real"):
        covariance.ward_certificate(np.eye(2), np.eye(2), float("nan"))
