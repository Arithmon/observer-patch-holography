"""Original-input controls for exponent cancellation in Kerr observables.

Expected values use only ordinary-sized constants at 200 digits and exact
Decimal tuple exponent shifts. No oracle forms 10**MAX_EMAX as an integer,
or reads a producer-generated mass, root, frequency, or scale.
"""
from decimal import Context, Decimal, MAX_EMAX, MIN_EMIN, localcontext

import pytest

import integer_k_comb_template as producer


D = Decimal
N = MAX_EMAX
HALF_N = N // 2
PI = D("3.1415926535897932384626433832795028841971693993751")


def _shift(value, exponent):
    """Shift a supplied finite Decimal exactly, without context arithmetic."""
    sign, digits, original_exponent = value.as_tuple()
    return D((sign, digits, original_exponent + exponent))


def _unit_constants():
    with localcontext(Context(prec=200)):
        c, gm = D(299792458), D("1.3271244e20")
        # chi=0, m=0, k=2: R=2GM/c² and f=c³ ln(2)/(16 pi² GM).
        radius = 2 * gm / c ** 2
        frequency_times_pi_squared = c ** 3 * D(2).ln() / (16 * gm)
        frequency = frequency_times_pi_squared / PI ** 2
    return radius, frequency, frequency_times_pi_squared


def _check(actual, expected):
    assert isinstance(actual, D) and actual.is_finite() and actual > 0
    # The reference carries about 150 guard digits beyond the tested output.
    with localcontext(Context(prec=250, Emax=MAX_EMAX, Emin=MIN_EMIN)):
        error = abs(actual - expected)
        ulp = D((0, (1,), max(actual.adjusted(), MIN_EMIN) - 49))
        assert error < ulp
        assert error <= abs(expected) * D("1e-49")


@pytest.mark.parametrize("case", ["radius", "detector_mass", "pi_and_mass"])
def test_representable_final_result_survives_overflowing_intermediate_scales(case):
    radius, frequency, frequency_pi2 = _unit_constants()
    with localcontext(Context(prec=50, Emax=MAX_EMAX, Emin=MIN_EMIN)):
        if case == "radius":
            expected = _shift(radius, N - 5)
            actual = producer.r_plus_si(D(f"1e{N - 5}"), D(0))
        elif case == "detector_mass":
            # M_det=10*10**N cannot fit, but the requested inverse observable can.
            expected = _shift(frequency, -N - 1)
            actual = producer.detector_frame_tooth_frequency_hz(
                D(f"1e{N}"), D(9), D(0), 0, 2, PI)
        else:
            # pi²*M has exponent N-6-2*floor(N/2), an ordinary small integer.
            expected = _shift(frequency_pi2, 2 * HALF_N - N + 6)
            actual = producer.tooth_frequency_hz(
                D(f"1e{N - 6}"), D(0), 0, 2, D(f"1e{-HALF_N}"))
    _check(actual, expected)


@pytest.mark.parametrize("case", ["radius", "detector_mass", "pi_and_mass"])
def test_nearby_resolved_scales_remain_supported(case):
    radius, frequency, frequency_pi2 = _unit_constants()
    with localcontext(Context(prec=50, Emax=MAX_EMAX, Emin=MIN_EMIN)):
        if case == "radius":
            expected = _shift(radius, N - 30)
            actual = producer.r_plus_si(D(f"1e{N - 30}"), D(0))
        elif case == "detector_mass":
            expected = _shift(frequency, -N)
            actual = producer.detector_frame_tooth_frequency_hz(
                D(f"1e{N}"), D(0), D(0), 0, 2, PI)
        else:
            expected = frequency_pi2
            actual = producer.tooth_frequency_hz(D("1e600"), D(0), 0, 2, D("1e-300"))
    _check(actual, expected)


@pytest.mark.parametrize("case", ["radius", "detector_mass", "pi_and_mass"])
def test_reversed_scale_cancellation_survives_changes_to_evaluation_order(case):
    radius, frequency, frequency_pi2 = _unit_constants()
    with localcontext(Context(prec=50, Emax=MAX_EMAX, Emin=MIN_EMIN)):
        if case == "radius":
            expected = _shift(radius, -N)
            actual = producer.r_plus_si(D(f"1e{-N}"), D(0))
        elif case == "detector_mass":
            # Original M_det = 10**(-N)*(1+10**N) = 1+10**(-N).
            # Thus f differs from the unit-mass reference by less than 10**(-N)
            # relatively, far below the 200-digit reference's truncation error.
            expected = frequency
            actual = producer.detector_frame_tooth_frequency_hz(
                D(f"1e{-N}"), D(f"1e{N}"), D(0), 0, 2, PI)
        else:
            expected = _shift(frequency_pi2, N - 2 * HALF_N)
            actual = producer.tooth_frequency_hz(
                D(f"1e{-N}"), D(0), 0, 2, D(f"1e{HALF_N}"))
    _check(actual, expected)
