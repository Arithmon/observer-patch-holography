"""Independent numerical and boundary checks of the rational gap certificate."""

from fractions import Fraction as F
import json

import mpmath as mp
import pytest

from .gap_certificate import (
    RationalInterval as Interval,
    certificate,
    exp_negative_bounds,
    hubbard_charge_gap_bounds,
    log_bounds,
    native_reference_data,
    pi_bounds,
    sqrt_bounds,
)


def _mp(value: F):
    return mp.mpf(value.numerator) / value.denominator


@pytest.fixture(scope="module")
def source():
    return native_reference_data()


@pytest.fixture(scope="module")
def gap_reference():
    # This is an independent diagnostic, not a proof step or certificate input.
    # Lieb--Wu's original Bessel formula differs from the positive K1 bounds.
    with mp.workdps(65):
        omega = mp.log(5)
        t, u = omega / 2, mp.log(mp.mpf(125) / 64)
        b = u / (2 * t)

        def integrand(x):
            return mp.mpf(1) / 4 if x == 0 else mp.besselj(1, x) / (x * (1 + mp.exp(b * x)))

        integral = mp.quad(integrand, list(range(0, 201, 5)))
        value = u - 4 * t + 8 * t * integral
        # |J1(x)| <= 1 gives this truncation bound; it does not certify quadrature.
        tail = 8 * t * mp.exp(-200 * b) / (200 * b)
        return +value, +tail


@pytest.mark.parametrize("value", [F(1, 5), F(1), F(5), F(125, 64)])
def test_exact_log_bounds_contain_independent_high_precision_value(value):
    result = log_bounds(value)
    with mp.workdps(80):
        reference = mp.log(_mp(value))
        assert _mp(result.lower) <= reference <= _mp(result.upper)
        assert _mp(result.upper - result.lower) < mp.mpf("1e-45")


def test_pi_and_square_root_bounds_contain_independent_values():
    with mp.workdps(80):
        pi = pi_bounds()
        assert _mp(pi.lower) < mp.pi < _mp(pi.upper)
        for value in (F(0), F(2), F(25, 4)):
            result = sqrt_bounds(Interval.enclose(value))
            reference = mp.sqrt(_mp(value))
            assert _mp(result.lower) <= reference <= _mp(result.upper)


@pytest.mark.parametrize("value", [F(0), F(1, 3), F(8), F(40)])
def test_range_reduced_reciprocal_taylor_encloses_exponential(value):
    result = exp_negative_bounds(Interval.enclose(value))
    with mp.workdps(80):
        reference = mp.exp(-_mp(value))
        assert _mp(result.lower) <= reference <= _mp(result.upper)
        assert result.lower >= 0


def test_outward_interval_operations_retain_exact_endpoint_examples():
    a = Interval(F(-2, 3), F(5, 7))
    b = Interval(F(1, 9), F(7, 5))
    for x in (F(-2, 3), F(5, 7)):
        for y in (F(1, 9), F(7, 5)):
            assert (a + b).lower <= x + y <= (a + b).upper
            assert (a * b).lower <= x * y <= (a * b).upper
            assert (a / b).lower <= x / y <= (a / b).upper
    with pytest.raises(ZeroDivisionError):
        _ = a / Interval(F(-1), F(1))


def test_native_coefficients_come_from_all_charge_sector_weights(source):
    assert source["weights"] == (F(2, 225), F(1, 36), F(2, 45))
    assert source["multiplicities"] == (45, 20, 1)
    assert source["clock_ratio"] == 5
    assert source["connected_ratio"] == F(125, 64)
    assert sum(n * w for n, w in zip(source["multiplicities"], source["weights"], strict=True)) == 1


def test_certificate_encloses_independent_original_lieb_wu_formula(source, gap_reference):
    omega = log_bounds(source["clock_ratio"])
    u = log_bounds(source["connected_ratio"])
    gap = hubbard_charge_gap_bounds(omega / 2, u)
    reference, tail = gap_reference
    with mp.workdps(65):
        assert _mp(gap.lower) < reference - tail
        assert reference + tail < _mp(gap.upper)
    assert gap.lower > F(1, 1000)
    assert gap.upper < F(11, 10000)
    assert gap.upper - gap.lower < F(1, 10**9)


def test_higher_global_taylor_bounds_tighten_the_source_enclosure(source, gap_reference):
    omega = log_bounds(source["clock_ratio"])
    u = log_bounds(source["connected_ratio"])
    coarse = hubbard_charge_gap_bounds(omega / 2, u, bessel_even=8)
    fine = hubbard_charge_gap_bounds(omega / 2, u, bessel_even=14)
    assert coarse.lower < fine.lower < fine.upper < coarse.upper
    with mp.workdps(65):
        assert _mp(fine.lower) < gap_reference[0] < _mp(fine.upper)


def test_zero_interaction_control_and_small_repulsion_envelopes():
    t = log_bounds(5) / 2
    zero = hubbard_charge_gap_bounds(t, Interval.enclose(0))
    assert zero.lower == zero.upper == 0
    half = hubbard_charge_gap_bounds(t, t / 2)
    quarter = hubbard_charge_gap_bounds(t, t / 4)
    eighth = hubbard_charge_gap_bounds(t, t / 8)
    assert 0 < eighth.lower < eighth.upper < quarter.lower < quarter.upper < half.lower
    assert eighth.upper < F(1, 10**19)


def test_json_safe_certificate_keeps_threshold_and_interpretation_distinct(monkeypatch):
    def forbidden_quadrature(*_args, **_kwargs):
        raise AssertionError("quadrature must not enter the certificate")

    monkeypatch.setattr(mp, "quad", forbidden_quadrature)
    result = certificate()
    json.dumps(result, allow_nan=False)
    gap = result["charge_gap"]
    exact_lo, exact_hi = F(gap["lower_exact"]), F(gap["upper_exact"])
    assert F.from_float(gap["lower"]) <= exact_lo <= exact_hi <= F.from_float(gap["upper"])
    assert F(result["holon_energy"]["lower_exact"]) <= exact_lo / 2
    assert F(result["holon_energy"]["upper_exact"]) >= exact_hi / 2
    assert result["proof"]["strictly_positive_certified_lower_endpoint"] is True
    assert result["proof"]["quadrature_used"] is False
    assert result["source"]["empirical_mass_inputs"] == []
    assert result["architecture"]["architecture_assignment_forced_by_native_repair"] is False
    assert result["interpretation"]["positive_global_many_body_gap_claimed"] is False
    assert result["interpretation"]["isolated_electron_pole_claimed"] is False
    assert result["interpretation"]["three_plus_one_identification"] is False
    assert result["interpretation"]["observed_Standard_Model_mass"] is False


def test_invalid_domains_are_rejected():
    with pytest.raises(ValueError):
        log_bounds(0)
    with pytest.raises(ValueError):
        log_bounds(-1)
    with pytest.raises(TypeError):
        log_bounds(1.5)
    with pytest.raises(ValueError):
        sqrt_bounds(Interval(F(-1), F(1)))
    with pytest.raises(ValueError):
        exp_negative_bounds(Interval(F(-1), F(0)))
    with pytest.raises(ValueError):
        hubbard_charge_gap_bounds(Interval.enclose(0), Interval.enclose(1))
    with pytest.raises(ValueError):
        hubbard_charge_gap_bounds(Interval.enclose(1), Interval.enclose(-1))
    with pytest.raises(ValueError):
        hubbard_charge_gap_bounds(Interval.enclose(1), Interval(F(0), F(1)))
    with pytest.raises(ValueError):
        hubbard_charge_gap_bounds(Interval.enclose(1), Interval.enclose(1), bessel_even=3)
    with pytest.raises(ValueError):
        Interval(F(2), F(1))
    with pytest.raises(ValueError):
        log_bounds(5, precision_bits=32)
