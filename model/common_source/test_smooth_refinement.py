"""Independent smooth instrument, energy and full cutoff arithmetic checks."""
from fractions import Fraction as F
import json
import math

import mpmath as mp
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import expm
import sympy as sp

from . import smooth_refinement as s
from . import refinement as r
from ..observer.source_scalar import field


def _mp_fraction(value):
    value = F(value)
    return mp.mpf(value.numerator) / value.denominator


@pytest.mark.parametrize("d", (2, 4))
def test_ancilla_unitary_gives_declared_smooth_kraus_maps_and_effect_sign(d):
    q = np.diag(np.sqrt(np.arange(1, d) / 2), 1)
    q += q.T
    identity = np.eye(d)
    hadamard = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
    phase_dagger = np.diag([1., -1j])
    unitary = (np.kron(hadamard @ phase_dagger, identity)
               @ expm(-.5j * np.kron(np.diag([1., -1.]), q))
               @ np.kron(hadamard, identity))
    expected = ((expm(-.5j * q) - 1j * expm(.5j * q)) / 2,
                (expm(-.5j * q) + 1j * expm(.5j * q)) / 2)
    sine = (expm(1j * q) - expm(-1j * q)) / (2j)
    for bit in (0, 1):
        actual = unitary[bit * d:(bit + 1) * d, :d]
        np.testing.assert_allclose(actual, expected[bit], atol=3e-15)
        np.testing.assert_allclose(actual.conj().T @ actual,
                                   (identity + (-1)**bit * sine) / 2, atol=3e-15)
    np.testing.assert_allclose(sum(m.conj().T @ m for m in expected), identity, atol=3e-15)


def test_exact_kraus_energy_identities_and_coefficient_bound():
    q = sp.Symbol("q", real=True)
    maps = [sp.cos(q / 2 - sp.pi / 4), sp.cos(q / 2 + sp.pi / 4)]
    assert sp.trigsimp(sum(m**2 for m in maps)) == 1
    assert sp.trigsimp(sum(m * sp.diff(m, q) for m in maps)) == 0
    assert sp.trigsimp(sum(sp.diff(m, q)**2 for m in maps)) == sp.Rational(1, 4)
    assert sp.Rational(1, 2) * sp.Rational(1, 4) == sp.Rational(1, 8)
    assert sp.Rational(1, 2) / sp.sqrt(2) < sp.Rational(1, 2)
    assert s.POINTER_MOMENT_RATE_UPPER == 4 * s.POINTER_LINEAR_COEFFICIENT_UPPER


@pytest.mark.parametrize("center,momentum", ((0., 0.), (.7, 1.1), (-.4, -.6)))
def test_infinite_continuum_kinetic_energy_injection_on_nontrivial_wavefunctions(center, momentum):
    def wave(x):
        value = np.pi**(-.25) * np.exp(-.5 * (x - center)**2 + 1j * momentum * x)
        return value, (-(x - center) + 1j * momentum) * value

    before = quad(lambda x: abs(wave(x)[1])**2 / 2, -np.inf, np.inf,
                  epsabs=1e-13, epsrel=1e-13)[0]
    after, mass, plus = 0., 0., 0.
    for sign in (-1, 1):
        def derivative_integrand(x):
            value, derivative = wave(x)
            angle = x / 2 + sign * np.pi / 4
            return abs(np.cos(angle) * derivative - .5 * np.sin(angle) * value)**2 / 2
        branch_mass = quad(lambda x: np.cos(x / 2 + sign * np.pi / 4)**2 * abs(wave(x)[0])**2,
                           -np.inf, np.inf, epsabs=1e-13, epsrel=1e-13)[0]
        after += quad(derivative_integrand, -np.inf, np.inf, epsabs=1e-13, epsrel=1e-13)[0]
        mass += branch_mass
        if sign == -1:
            plus = branch_mass
    assert mass == pytest.approx(1, abs=3e-12)
    assert plus == pytest.approx((1 + np.exp(-.25) * np.sin(center)) / 2, abs=3e-12)
    assert after - before == pytest.approx(1 / 8, abs=3e-12)


def test_infinite_energy_identity_is_not_falsely_transferred_to_finite_ccr_cutoff():
    q = field.quadratures(2)[0]
    kinetic = field.projected_squares(2)[1] / 2
    maps = [(expm(-.5j * q) + sign * 1j * expm(.5j * q)) / 2 for sign in (-1, 1)]
    increment = sum(m.conj().T @ kinetic @ m for m in maps) - kinetic
    assert np.trace(increment) == pytest.approx(0, abs=3e-15)
    assert np.linalg.norm(increment - np.eye(2) / 8, 2) > .1
    assert not s.status()["finite_cutoff_canonical_energy_identity_claimed"]


@pytest.mark.parametrize("d", (2, 4, 16, 64, 256))
def test_full_pointer_duhamel_bound_is_directed_above_independent_high_precision_formula(d):
    certificate = s.smooth_occupation_bound(d)
    with mp.workdps(90):
        initial = mp.sqrt(mp.mpf(143748821) / 8192000)
        growth = mp.exp(mp.mpf(39) / 16)
        before_pointer = initial * growth
        interaction_integral = initial * (growth - 1) / (mp.mpf(13) / 8)
        pair = d // 2
        tail = mp.mpf(81) / 80 * math.comb(2 * pair, pair) / (mp.mpf(4)**pair * 81**pair)
        projected = mp.sqrt(2 * tail) + interaction_integral * (
            mp.mpf(9) / (32 * (d + 2)) + 1 / (4 * mp.mpf(d + 2)**mp.mpf("1.5")))
        pointer_residual = before_pointer * (mp.e - 1) / (2 * mp.mpf(d + 2)**mp.mpf("1.5"))
        outside = before_pointer * mp.e / (d + 2)**2
        expected = mp.sqrt((projected + pointer_residual)**2 + outside**2)
        actual = _mp_fraction(certificate["complete_stinespring_vector_error_upper"])
        assert expected <= actual
        assert actual - expected < mp.mpf("1e-27")


def test_full_bound_decreases_in_cutoff_and_does_not_erase_pointer_cost():
    values = []
    for d in (4, 16, 64, 256, 1024):
        smooth = s.smooth_occupation_bound(d)
        value = F(smooth["complete_stinespring_vector_error_upper"])
        old = r.occupation_evolution_bound(d)
        assert value > F(old["pre_readout_stinespring_vector_error_upper"])
        assert F(smooth["pointer_projected_duhamel_residual_integral_upper"]) > 0
        values.append(value)
    assert all(a > b for a, b in zip(values, values[1:]))
    assert values[-1] < F(1, 100)


@pytest.mark.parametrize("d", (2, 4))
def test_complete_budget_adds_both_temporal_segments_with_explicit_vacuity(d):
    coarse = s.complete_error_budget(d, 128)
    fine = s.complete_error_budget(d, 256)
    for a, b in zip(coarse["ideal_temporal_segment_certificates"], fine["ideal_temporal_segment_certificates"]):
        assert F(a["unclipped_operator_error_upper"]) == 4 * F(b["unclipped_operator_error_upper"])
    cutoff = F(fine["occupation_certificate"]["complete_instrument_trace_distance_upper"])
    temporal = sum(F(item["ideal_channel_half_diamond_error_upper"])
                   for item in fine["ideal_temporal_segment_certificates"])
    assert F(fine["complete_ideal_product_to_infinite_instrument_error_upper"]) == min(F(1), cutoff + temporal)
    # The existing tiny source registers do not yet earn a precise infinite
    # detector prediction merely because the asymptotic rate is controlled.
    assert not fine["nonvacuous_combined_bound"]
    assert not fine["floating_execution_or_physical_service_error_included"]


def test_new_smooth_scope_does_not_overwrite_old_sharp_limit_or_claim_diamond_bound():
    report = s.report((4, 64))
    json.dumps(report)
    assert not report["status"]["old_sharp_readout_gap_closed_by_this_result"]
    assert not report["status"]["spatial_continuum_or_empirical_completion"]
    assert r.sharp_readout_status()["sharp_readout_rate"] == "not supplied"
    for item in report["occupation_bounds"]:
        assert item["includes_final_smooth_pointer_and_all_unnormalized_branches"]
        assert not item["bound_is_an_arbitrary_input_diamond_norm_claim"]
        assert not item["applies_to_old_sharp_sign_instrument_or_modular_clock"]


@pytest.mark.parametrize("d,time", ((True, F(1)), (3, F(1)), (4, True), (4, .5), (4, F(2))))
def test_rejects_undeclared_cutoff_or_model_time(d, time):
    with pytest.raises(ValueError):
        s.smooth_occupation_bound(d, time)


@pytest.mark.parametrize("d,steps", ((16, 512), (2, True), (2, 0), (2, 1.5)))
def test_budget_requires_actual_source_dimension_and_integer_steps(d, steps):
    with pytest.raises(ValueError):
        s.complete_error_budget(d, steps)
