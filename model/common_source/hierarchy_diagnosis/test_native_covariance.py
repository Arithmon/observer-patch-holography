"""Independent moment, normalization and mixing controls."""
import json

import sympy as s

from . import native_covariance as c


def test_complete_connected_reference_covariance_matches_direct_pair_weights():
    covariance = c.port_covariance()
    means, second, *_ = c._data()
    assert means[0] == means[1] == s.Rational(29, 90)
    assert all(means[i] == s.Rational(61, 450) for i in range(2, 12))
    assert second[0, 1] == s.Rational(2, 45)
    assert second[0, 2] == s.Rational(1, 36)
    assert second[2, 3] == s.Rational(2, 225)
    assert covariance == second-means*means.T
    assert covariance == covariance.T
    assert covariance*s.ones(12, 1) == s.zeros(12, 1)
    assert covariance.rank() == 11


def test_full_spectrum_and_family_covariance_positive_ranks():
    covariance = c.port_covariance()
    assert covariance.eigenvals() == {
        s.Integer(0): 1, s.Rational(19, 150): 9,
        s.Rational(644, 3375): 1, s.Rational(5, 18): 1}
    compressed = c.family_covariance()
    assert compressed.rank() == 3
    assert all(value.is_positive is True for value in compressed.diagonal())
    assert compressed == s.diag(s.Rational(37, 225)-17*s.sqrt(5)/2250,
                                s.Rational(3283, 22500)+433*s.sqrt(5)/112500,
                                s.Rational(19, 150))


def test_canonical_family_normalization_and_F_commutation():
    _, _, covariance, basis, gram, feature, compressed = c._data()
    assert gram == (10+2*s.sqrt(5))*s.eye(3)
    normalized = basis/s.sqrt(10+2*s.sqrt(5))
    assert (normalized.T*normalized-s.eye(3)).applyfunc(s.simplify) == s.zeros(3)
    assert (normalized.T*covariance*normalized-compressed).applyfunc(s.simplify) == s.zeros(3)
    assert feature*compressed == compressed*feature
    # Even distinct nonlinear spectral maps of these two matrices do not
    # change their common left mixing basis.
    up = feature**2 + s.I*compressed
    down = compressed**3 - 2*feature
    assert (up*up.H*down*down.H-down*down.H*up*up.H).applyfunc(s.expand) == s.zeros(3)


def test_report_preserves_statistical_scope_and_strict_spread_bound():
    packet = c.report()
    json.dumps(packet, allow_nan=False)
    compression = packet["family_compression"]
    spread = s.sympify(compression["largest_over_smallest"])
    assert (s.Rational(5, 4)-spread).is_positive is True
    assert packet["comparison_to_existing_feature"]["commutes_with_F_exact"]
    assert not packet["scope"]["new_coherent_interaction_or_mass_selected"]
    assert not packet["scope"]["covariance_and_F_alone_supply_nontrivial_quark_mixing"]
    assert not packet["scope"]["arbitrary_interacting_native_histories_excluded"]
