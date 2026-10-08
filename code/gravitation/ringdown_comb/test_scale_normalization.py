"""Original-input controls for exponent cancellation in Kerr observables.

Expected values use only ordinary-sized constants at 200 digits and exact
Decimal tuple exponent shifts. No oracle forms 10**MAX_EMAX as an integer,
or reads a producer-generated mass, root, frequency, or scale.
"""
from decimal import Context, Decimal, MAX_EMAX, MIN_EMIN, Rounded, localcontext

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


@pytest.mark.parametrize("pi,chi,m", [
    (D(".125"), D(".67"), -2),
    (D(1), D(".67"), -1),
    (D("0." + "9" * 80), D("-.67"), -1),
    (D("1." + "0" * 79 + "1"), D(".67"), 2),
    (D(10), D("-.67"), 2),
])
def test_both_pi_scalings_preserve_signed_original_input_frequencies(pi, chi, m):
    import mpmath as mp

    with mp.workdps(250):
        x, supplied_pi = mp.mpf(str(chi)), mp.mpf(str(pi))
        root = mp.sqrt(1 - x * x)
        kappa = mp.mpf(299792458) ** 3 * root / (2 * 62 * mp.mpf("1.3271244e20") * (1 + root))
        omega = mp.mpf(299792458) ** 3 * x / (2 * 62 * mp.mpf("1.3271244e20") * (1 + root))
        expected = D(mp.nstr(m * omega / (2 * supplied_pi) + kappa * mp.log(2) / (4 * supplied_pi ** 2), 200))
    with localcontext(Context(prec=50)):
        actual = producer.tooth_frequency_hz(D(62), chi, m, 2, pi)
    assert (actual < 0) == (expected < 0)
    _check(actual.copy_abs(), expected.copy_abs())


def test_true_output_overflow_refuses_without_mutating_caller_flags():
    with localcontext(Context(prec=50, Emax=MAX_EMAX, Emin=MIN_EMIN)) as context:
        context.flags[Rounded] = True
        flags, traps = dict(context.flags), dict(context.traps)
        with pytest.raises(producer.NumericalResolutionError):
            producer.omega_h_si(D(f"1e{-N}"), D(".67"))
        assert dict(context.flags) == flags
        assert dict(context.traps) == traps


@pytest.mark.parametrize("upper", [False, True])
def test_signed_cancellation_survives_combined_implementation_scale_exponents(upper):
    import mpmath as mp

    with mp.workdps(250):
        log_three = mp.log(3)
        coefficient = int(mp.floor(log_three / 2 * 10 ** 60)) + int(upper)
        a = D((0, tuple(map(int, str(coefficient))), -60))
        chi = _shift(a, -HALF_N)
        numerator, denominator = a.as_integer_ratio()
        a_exact = mp.mpf(numerator) / denominator
        # Original chi=a*10**(-H), pi=10**H and M=10**(-N).
        # Replacing sqrt(1-chi²) by1 changes f by <10**(-N+10),
        # negligible compared with the retained 200-digit oracle and its
        # nonzero ~10**(-55) answer. All finite supplied coefficients remain.
        kappa_at_unit_mass = mp.mpf(299792458) ** 3 / (4 * mp.mpf("1.3271244e20"))
        expected = D(mp.nstr(kappa_at_unit_mass * (log_three - 2 * a_exact)
                            / 4 * 10 ** (N - 2 * HALF_N), 200))
    with localcontext(Context(prec=50, Emax=MAX_EMAX, Emin=MIN_EMIN)):
        actual = producer.tooth_frequency_hz(D(f"1e{-N}"), chi, -1, 3, D(f"1e{HALF_N}"))
    assert (actual < 0) == upper
    _check(actual.copy_abs(), expected.copy_abs())
