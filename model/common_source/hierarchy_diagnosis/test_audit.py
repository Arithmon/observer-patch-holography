"""Independent matrix and retained-record checks for hierarchy diagnosis."""
import json

import numpy as np
import pytest
import sympy as s

from . import audit


def test_native_affine_relation_is_exact_for_complex_coefficients():
    feature, values, projectors = audit.spectral_data()
    assert s.simplify((values[2]-values[0])/(values[1]-values[0])-(3+s.sqrt(5))/2) == 0
    assert sum(projectors, s.zeros(3)) == s.eye(3)
    assert sum((x*p for x, p in zip(values, projectors, strict=True)), s.zeros(3)) == feature
    u, v = 3+5*s.I, -2+7*s.I
    phi = (1+s.sqrt(5))/2
    z = [u+v*x for x in values]
    assert s.simplify(z[2]-phi**2*z[1]+phi*z[0]) == 0


def test_affine_singular_bound_survives_phases_and_independent_field_bases():
    # Independent numerical SVD of complex affine matrices, with arbitrary
    # unitary left/right field bases; includes cancellation at each eigenaxis.
    feature, values, _ = audit.spectral_data()
    f = np.array(feature, dtype=complex)
    rng = np.random.default_rng(20261026)
    phi = (1+np.sqrt(5))/2
    choices = [(-float(x), 1) for x in values]
    choices += [(complex(*rng.normal(size=2)), complex(*rng.normal(size=2))) for _ in range(32)]
    for u, v in choices:
        left, _ = np.linalg.qr(rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)))
        right, _ = np.linalg.qr(rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)))
        masses = np.sort(np.linalg.svd(left @ (u*np.eye(3)+v*f) @ right, compute_uv=False))
        assert masses[2] <= phi**2*masses[1]+phi*masses[0]+1e-12


def test_affine_exclusion_does_not_exclude_arbitrary_single_feature_functions():
    # A synthetic counterexample confines the no-go to affine laws. It is
    # deliberately not a fit to an empirical mass spectrum.
    feature, values, projectors = audit.spectral_data()
    targets = (1, 2, 100)
    polynomial = sum((x*p for x, p in zip(targets, projectors, strict=True)), s.zeros(3))
    assert sorted(polynomial.eigenvals()) == list(targets)
    phi = (1+s.sqrt(5))/2
    assert (100-2*phi**2-phi).is_positive
    assert polynomial*feature == feature*polynomial
    # A distinct complex function changes spectra but shares all axes.
    other = sum((x*p for x, p in zip((s.I, 4-2*s.I, 3), projectors, strict=True)), s.zeros(3))
    a, b = polynomial*polynomial.H, other*other.H
    assert a*b == b*a


def test_current_records_recover_the_same_bad_hierarchy():
    retained = json.loads(audit.BASELINE.read_text())
    result = audit.record_audit(retained)
    assert result["record_count"] == 64
    assert result["maximum_absolute_record_action_residual"] < 1e-11
    assert 1 < result["decoded_ascending_ratios"][1] < result["decoded_ascending_ratios"][2] < 1.5


def test_benchmark_is_separate_and_complex_affine_gap_is_large():
    source = audit.source_algebra()
    retained = json.loads(audit.BASELINE.read_text())
    benchmark = audit.benchmark_exclusion(source, retained)
    assert source["empirical_inputs"] == []
    assert 543 < benchmark["maximum_affine_heavy_ratio_at_admitted_middle_ratio"]["float_diagnostic"] < 544
    assert benchmark["affine_class_excluded_exact"]
    controls = audit.nonlinear_shapes(retained)
    lo, hi = controls["benchmark_Q_window_diagnostic"]
    assert .52 < lo < hi < .54
    assert all(row["adjacent_log_spacing_diagnostic"] > 1.08 for row in controls["fixed_source_controls"])
    assert not any(row["both_ratios_pass"] for row in controls["original_exact_interval_comparisons"])


def test_kinetic_bound_is_saturated_by_an_adversarial_diagonal_example():
    # Do not merely retest the formula: diagonal kinetic matrices produce
    # the extremal spectral ratio directly after canonical normalization.
    m = np.diag([1., 1.2, 1.5])
    eps = .9
    z = np.diag([1+eps, 1., 1-eps])
    invroot = np.diag(1/np.sqrt(np.diag(z)))
    masses = np.linalg.svd(invroot @ m @ invroot, compute_uv=False)
    assert masses[0]/masses[-1] == pytest.approx(1.5*(1+eps)/(1-eps))
    bounds = audit.correction_bounds(s.Rational(3, 2), 3470)
    assert s.sympify(bounds["required_epsilon_if_both_Z_within_epsilon_of_identity_at_least"]["exact"]) == s.Rational(6937, 6943)
    assert s.sympify(bounds["required_additive_norm_over_original_minimum_at_least"]["exact"]) == s.Rational(6937, 6942)
    with pytest.raises(ValueError):
        audit.correction_bounds(2, 1)


def test_additive_bound_is_attained_without_dropping_or_relabeling_modes():
    original = np.diag([1., 1.2, 1.5])
    d = .9
    changed = original + np.diag([-d, 0, d])
    masses = np.linalg.svd(changed, compute_uv=False)
    assert masses[0]/masses[-1] == pytest.approx((1.5+d)/(1-d))
