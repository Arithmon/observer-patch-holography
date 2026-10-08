"""Original-input arithmetic controls; no event data or physical derivation.

The independent oracle evaluates the supplied finite Decimal values directly
at 200 digits. It never reads the producer's Kerr roots, spacings, rounded
rotation line, or other intermediate results.
"""
from __future__ import annotations

from decimal import (
    Context, Decimal, Inexact, MAX_EMAX, MIN_EMIN, ROUND_DOWN,
    ROUND_HALF_EVEN, ROUND_UP, Rounded, localcontext,
)
from fractions import Fraction
import json
from pathlib import Path

import pytest

import integer_k_comb_template as producer
import verify_integer_k_comb_independent as verifier


D = Decimal
PI = D("3.1415926535897932384626433832795028841971693993751")
CANCEL_CHI = D("0.10965258099938507963811633773723788073636340857259")
NEAR_CHIS = (D("0." + "9" * 50 + "7"), D("0." + "9" * 60))
ABOVE_ONE = D("1." + "0" * 59 + "1")
REFUSAL = (TypeError, ValueError, ArithmeticError)


def _oracle(mass, chi, m=2, k=2, pi=PI, a=D(1)):
    with localcontext(Context(prec=200, rounding=ROUND_HALF_EVEN,
                              Emax=MAX_EMAX, Emin=MIN_EMIN)):
        c, gm = D(299792458), D(mass) * D("1.3271244e20")
        root = (D(1) - chi * chi).sqrt()
        horizon = D(1) + root
        omega = c ** 3 * chi / (2 * gm * horizon)
        kappa = c ** 3 * root / (2 * gm * horizon)
        g = 4 * gm * kappa / c ** 3
        base = kappa / (4 * pi ** 2)
        offset = base * D(k).ln()
        rotation = D(m) * omega / (2 * pi)
        result = dict(root=root, radius=gm * horizon / c ** 2,
                      omega=omega, kappa=kappa, g=g, base=base,
                      offset=offset, rotation=rotation,
                      frequency=rotation + offset)
        if root:
            result["linewidth"] = 64 * pi ** 2 * D("2e-4") / (a * g ** 2 * D(k).ln())
        return result


def _assert_close(actual, expected, precision):
    assert isinstance(actual, Decimal)
    assert actual.is_finite(), actual
    if expected == 0:
        assert actual == 0
        return
    with localcontext(Context(prec=220)):
        tolerance = D(10) ** (2 - precision)
        assert abs(actual - expected) <= abs(expected) * tolerance, (actual, expected)


PUBLIC_COMPONENTS = {
    "root": lambda mass, chi, m, k, pi: producer.sqrt_one_minus_chi_squared(chi),
    "radius": lambda mass, chi, m, k, pi: producer.r_plus_si(mass, chi),
    "omega": lambda mass, chi, m, k, pi: producer.omega_h_si(mass, chi),
    "kappa": lambda mass, chi, m, k, pi: producer.kappa_si(mass, chi),
    "g": lambda mass, chi, m, k, pi: producer.g_of_chi(chi),
    "base": lambda mass, chi, m, k, pi: producer.base_spacing_hz_per_nat(mass, chi, pi),
    "offset": lambda mass, chi, m, k, pi: producer.tooth_offset_hz(mass, chi, k, pi),
    "rotation": lambda mass, chi, m, k, pi: producer.rotation_line_hz(mass, chi, m, pi),
    "frequency": lambda mass, chi, m, k, pi: producer.tooth_frequency_hz(mass, chi, m, k, pi),
    "linewidth": lambda mass, chi, m, k, pi: producer.linewidth_fraction(D(1), chi, k, pi),
}


@pytest.mark.parametrize("chi", NEAR_CHIS)
@pytest.mark.parametrize("negative", [False, True])
@pytest.mark.parametrize("component", PUBLIC_COMPONENTS)
def test_original_near_extremal_inputs_preserve_each_component(chi, negative, component):
    chi = chi.copy_negate() if negative else chi
    expected = _oracle(D(62), chi)
    with localcontext(Context(prec=50)):
        actual = PUBLIC_COMPONENTS[component](D(62), chi, 2, 2, PI)
    _assert_close(actual, expected[component], 50)


@pytest.mark.parametrize("mass", [D(62), D("6.2e-49")])
@pytest.mark.parametrize("negative", [False, True])
def test_signed_rotation_cancellation_is_evaluated_from_original_inputs(mass, negative):
    chi = CANCEL_CHI.copy_negate() if negative else CANCEL_CHI
    m = 1 if negative else -1
    expected = _oracle(mass, chi, m=m)["frequency"]
    assert expected > 0
    with localcontext(Context(prec=50)):
        actual = producer.tooth_frequency_hz(mass, chi, m, 2, PI)
    _assert_close(actual, expected, 50)


@pytest.mark.parametrize("chi", [D(0), D(".67"), D("-.67"), D(1), D(-1)])
@pytest.mark.parametrize("m", [-3, 0, 2])
def test_resolved_signed_inputs_and_exact_kerr_endpoints_remain_supported(chi, m):
    expected = _oracle(D("10.5"), chi, m=m, k=3)
    with localcontext(Context(prec=50)):
        for component, evaluate in PUBLIC_COMPONENTS.items():
            if component == "linewidth" and chi.copy_abs() == 1:
                continue
            actual = evaluate(D("10.5"), chi, m, 3, PI)
            _assert_close(actual, expected[component], 50)


@pytest.mark.parametrize("negative", [False, True])
@pytest.mark.parametrize("component", ["root", "frequency"])
def test_original_superextremal_input_is_refused_before_context_rounding(negative, component):
    chi = ABOVE_ONE.copy_negate() if negative else ABOVE_ONE
    with localcontext(Context(prec=50)), pytest.raises(REFUSAL):
        PUBLIC_COMPONENTS[component](D(62), chi, 2, 2, PI)


@pytest.mark.parametrize("precision", [4, 6, 50])
def test_entropy_loss_uses_the_exact_integer_transition(precision):
    before, k = 2 * 10 ** 1000, 2
    assert producer.integer_division_after(before, k) == 10 ** 1000
    with localcontext(Context(prec=200)):
        expected = D(k).ln()
    with localcontext(Context(prec=precision)):
        signed, loss = producer.transition_entropy_nats(before, k)
    assert signed < 0 < loss
    _assert_close(signed, expected.copy_negate(), precision)
    _assert_close(loss, expected, precision)


@pytest.mark.parametrize("before,k", [(60, 3), (2 ** 4096, 2), (3 ** 200, 3)])
def test_exact_integer_division_preserves_arbitrary_size(before, k):
    result = producer.integer_division_after(before, k)
    assert isinstance(result, int) and not isinstance(result, bool)
    assert k * result == before


def test_float_divisor_cannot_round_an_odd_original_dimension_to_even():
    before = 2 ** 53 + 1
    assert before % 2 == 1
    with pytest.raises(REFUSAL):
        producer.integer_division_after(before, 2.0)


@pytest.mark.parametrize("bad", [2.0, D(2), Fraction(2), True, "2", None])
@pytest.mark.parametrize("argument", ["dimension", "divisor", "tooth", "azimuthal"])
def test_integer_coordinates_reject_noninteger_types(bad, argument):
    with localcontext(Context(prec=50)), pytest.raises(REFUSAL):
        if argument == "dimension":
            producer.integer_division_after(bad, 2)
        elif argument == "divisor":
            producer.integer_division_after(6, bad)
        elif argument == "tooth":
            producer.tooth_frequency_hz(D(62), D(".67"), 2, bad, PI)
        else:
            producer.rotation_line_hz(D(62), D(".67"), bad, PI)


@pytest.mark.parametrize("bad", [D(0), D(-1), D("Infinity"), D("-Infinity"), D("NaN")])
@pytest.mark.parametrize("argument", ["mass", "pi", "linewidth_a"])
def test_nonpositive_or_nonfinite_required_positive_scalars_are_refused(bad, argument):
    with localcontext(Context(prec=50)), pytest.raises(REFUSAL):
        if argument == "mass":
            producer.tooth_frequency_hz(bad, D(".67"), 2, 2, PI)
        elif argument == "pi":
            producer.tooth_frequency_hz(D(62), D(".67"), 2, 2, bad)
        else:
            producer.linewidth_fraction(bad, D(".67"), 2, PI)


@pytest.mark.parametrize("bad", [D(-1), D("Infinity"), D("-Infinity"), D("NaN")])
def test_detector_frame_refuses_invalid_redshift(bad):
    with localcontext(Context(prec=50)), pytest.raises(REFUSAL):
        producer.detector_frame_tooth_frequency_hz(D(62), bad, D(".67"), 2, 2, PI)


@pytest.mark.parametrize("bad", [D("Infinity"), D("-Infinity"), D("NaN")])
def test_kerr_spin_must_be_finite(bad):
    with localcontext(Context(prec=50)), pytest.raises(REFUSAL):
        producer.sqrt_one_minus_chi_squared(bad)


@pytest.mark.parametrize("a", [D(".125"), D(1), D(10), D(32)])
def test_positive_linewidth_parameter_is_not_limited_to_receipt_display_range(a):
    expected = _oracle(D(62), D(".67"), k=3, a=a)["linewidth"]
    with localcontext(Context(prec=50)):
        actual = producer.linewidth_fraction(a, D(".67"), 3, PI)
    _assert_close(actual, expected, 50)


@pytest.mark.parametrize("chi", [D(1), D(-1)])
def test_linewidth_refuses_singular_extremal_limit(chi):
    with localcontext(Context(prec=50)), pytest.raises(REFUSAL):
        producer.linewidth_fraction(D(1), chi, 2, PI)


@pytest.mark.parametrize("value,root", [("1e-600", "1e-300"), ("1e600", "1e300"), ("4", "2"), (".04", ".2"), ("0", "0")])
def test_independent_verifier_sqrt_preserves_original_decimal_scale(value, root):
    with localcontext(Context(prec=60)):
        actual = verifier.v_sqrt(D(value))
    _assert_close(actual, D(root), 60)


def _context_state(context):
    return (context.prec, context.rounding, context.Emin, context.Emax,
            context.capitals, context.clamp, dict(context.flags), dict(context.traps))


def _ambient(case):
    context = Context(prec=13, rounding=ROUND_HALF_EVEN)
    if case == "round_down":
        context.rounding = ROUND_DOWN
    elif case == "round_up":
        context.rounding = ROUND_UP
    elif case == "exponents":
        context.Emin, context.Emax, context.clamp = -9, 9, 1
    elif case == "traps":
        context.traps[Inexact] = context.traps[Rounded] = True
    return context


BUILDERS = (producer.build_receipt, verifier.v_build_receipt)


@pytest.mark.parametrize("builder", BUILDERS, ids=["producer", "independent_verifier"])
@pytest.mark.parametrize("case", ["precision", "round_down", "round_up", "exponents", "traps"])
def test_canonical_receipts_are_independent_of_ambient_decimal_policy(builder, case):
    expected = (Path(__file__).parent / "runtime" / "integer_k_comb_template_receipt.json").read_bytes()
    with localcontext(_ambient(case)):
        receipt = builder()
    actual = (json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n").encode("ascii")
    assert actual == expected


@pytest.mark.parametrize("builder", BUILDERS, ids=["producer", "independent_verifier"])
def test_receipt_builder_preserves_the_callers_complete_decimal_context(builder):
    with localcontext(_ambient("precision")) as context:
        context.flags[Rounded] = True
        before = _context_state(context)
        builder()
        assert _context_state(context) == before


def test_large_mass_intermediate_does_not_hide_representable_surface_gravity():
    mass = D("1e999980")
    expected = _oracle(mass, D(0))["kappa"]
    assert expected > 0
    with localcontext(Context(prec=50)):
        actual = producer.kappa_si(mass, D(0))
    _assert_close(actual, expected, 50)


def test_unrepresentable_nonzero_result_is_refused_instead_of_rounded_to_zero():
    # At precision 9 and Emin=-3, the smallest Decimal subnormal is 1e-11.
    # The original operation is exactly (1+0)*1e-12 = 1e-12, not zero.
    with localcontext(Context(prec=9, Emin=-3, Emax=99)), pytest.raises(REFUSAL):
        producer.detector_frame_mass_solar(D("1e-12"), D(0))


def test_exact_representable_subnormal_remains_supported():
    with localcontext(Context(prec=9, Emin=-3, Emax=99)):
        assert producer.detector_frame_mass_solar(D("1e-10"), D(0)) == D("1e-10")


def test_public_arithmetic_uses_half_even_and_preserves_flags_and_traps():
    expected = _oracle(D(62), D(".67"), m=-1, k=3)["frequency"]
    values = []
    for rounding in (ROUND_DOWN, ROUND_UP):
        with localcontext(Context(prec=35, rounding=rounding)) as context:
            context.traps[Inexact] = context.traps[Rounded] = True
            context.flags[Rounded] = True
            before = _context_state(context)
            values.append(producer.tooth_frequency_hz(D(62), D(".67"), -1, 3, PI))
            assert _context_state(context) == before
    assert values[0] == values[1]
    _assert_close(values[0], expected, 35)


def test_detector_frame_cancellation_uses_original_source_parameters():
    mass, redshift = D("6.2e-49"), D(".123456789123456789123456789123456789123456789123456789")
    with localcontext(Context(prec=200)):
        original_detector_mass = mass * (1 + redshift)
    expected = _oracle(original_detector_mass, CANCEL_CHI, m=-1)["frequency"]
    with localcontext(Context(prec=50)):
        actual = producer.detector_frame_tooth_frequency_hz(mass, redshift, CANCEL_CHI, -1, 2, PI)
    _assert_close(actual, expected, 50)


def test_compute_pi_supports_the_requested_high_precision():
    import mpmath

    with mpmath.workdps(160):
        expected = D(mpmath.nstr(mpmath.pi, 150))
    with localcontext(Context(prec=100)):
        actual = producer.compute_pi()
    _assert_close(actual, expected, 100)


def test_near_cancellation_retains_a_genuinely_negative_frequency():
    chi = D("0.10965258099938507963811633773723788073636340857260")
    expected = _oracle(D("6.2e-49"), chi, m=-1)["frequency"]
    assert expected < 0
    with localcontext(Context(prec=50)):
        actual = producer.tooth_frequency_hz(D("6.2e-49"), chi, -1, 2, PI)
    _assert_close(actual, expected, 50)


def test_independent_verifier_sqrt_supports_the_requested_high_precision():
    with localcontext(Context(prec=250)):
        expected = D(2).sqrt()
    with localcontext(Context(prec=200)):
        actual = verifier.v_sqrt(D(2))
    _assert_close(actual, expected, 200)


@pytest.mark.parametrize("value", [D(-1), D("NaN"), D("Infinity"), True, 2.0])
def test_independent_verifier_sqrt_refuses_invalid_original_scalars(value):
    with pytest.raises(REFUSAL):
        verifier.v_sqrt(value)


@pytest.mark.parametrize("value", [0, -1, 13, True, 2.5, D(2), Fraction(2)])
def test_independent_verifier_log_refuses_invalid_original_integers(value):
    with pytest.raises(REFUSAL):
        verifier.v_ln_nat(value)


@pytest.mark.parametrize("k", [1, 2.0])
@pytest.mark.parametrize("consumer", ["ladder", "position", "kms", "offset", "linewidth"])
def test_every_public_tooth_consumer_enforces_the_integer_domain(k, consumer):
    calls = {
        "ladder": lambda: producer.ladder_ratio(k),
        "position": lambda: producer.universal_position(k, PI),
        "kms": lambda: producer.kms_weight(k),
        "offset": lambda: producer.tooth_offset_hz(D(62), D(".67"), k, PI),
        "linewidth": lambda: producer.linewidth_fraction(D(1), D(".67"), k, PI),
    }
    with localcontext(Context(prec=50)), pytest.raises(REFUSAL):
        calls[consumer]()


@pytest.mark.parametrize("value", [0.5, True])
@pytest.mark.parametrize("consumer", ["dec", "spin", "mass"])
def test_scalar_conversion_does_not_silently_accept_float_or_boolean(value, consumer):
    calls = {
        "dec": lambda: producer.dec(value),
        "spin": lambda: producer.g_of_chi(value),
        "mass": lambda: producer.gm_si(value),
    }
    with localcontext(Context(prec=50)), pytest.raises(REFUSAL):
        calls[consumer]()
