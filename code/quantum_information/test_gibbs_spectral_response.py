"""Original-spectrum controls for Gibbs response and its returned state."""

import mpmath
import numpy as np
import pytest

from quantum_information.entropy_response import gibbs_entropy_response


def partition(energies, beta):
    """Independent scalar Boltzmann moments, without a matrix eigensolve."""
    ctx = mpmath.mp.clone()
    ctx.dps = 400
    e, b = [ctx.mpf(float(x)) for x in energies], ctx.mpf(float(beta))
    weights = [ctx.exp(-b*x) for x in e]
    total = ctx.fsum(weights)
    probabilities = [w/total for w in weights]
    mean = ctx.fsum(p*x for p, x in zip(probabilities, e))
    variance = ctx.fsum(p*(x-mean)**2 for p, x in zip(probabilities, e))
    return ctx, probabilities, variance, -b*variance


@pytest.mark.parametrize("power, units, sign, origin", [
    (40, 0, 1, 0.), (60, 0, 1, 0.), (60, 400, 1, 0.), (60, -400, 1, 0.),
    (60, 0, -1, 0.), (60, 400, -1, 0.), (60, 0, 1, 2.**-10),
])
def test_small_original_splitting_controls_each_response_component(power, units, sign, origin):
    scale = 2.**units
    energies = [origin, origin+sign*scale*2.**-power, origin+sign*scale]
    beta = sign*700./scale
    _, probabilities, variance, entropy_derivative = partition(energies, beta)
    report = gibbs_entropy_response(np.diag(energies), beta)
    for name, expected in (("energy_variance", variance), ("dt_dlambda", -variance),
                           ("ds_dlambda", entropy_derivative)):
        assert report[name] == pytest.approx(float(expected), rel=3e-14, abs=0)
    np.testing.assert_allclose(np.diag(report["state"]).real,
                               [float(p) for p in probabilities], rtol=3e-14, atol=0)
    assert report["ds_dt"] == pytest.approx(beta, rel=3e-15, abs=0)


@pytest.mark.parametrize("phase", [1., 1j])
def test_exact_complex_doublet_keeps_small_state_coherence(phase):
    gap = 2.**-60
    observable = np.array([[gap/2, phase*gap/2, 0],
                           [np.conj(phase)*gap/2, gap/2, 0], [0, 0, 1.]])
    _, probabilities, variance, ds = partition([0., gap, 1.], 20.)
    report = gibbs_entropy_response(observable, 20.)
    # The exact two spectral projectors are (I +/- [[0,phase],[phase*,0]])/2.
    expected_coherence = complex(phase)*float((probabilities[1]-probabilities[0])/2)
    assert report["state"][0, 1] == pytest.approx(expected_coherence, rel=3e-14, abs=0)
    assert report["state"][1, 0] == report["state"][0, 1].conjugate()
    assert report["state"][2, 2].real == pytest.approx(float(probabilities[2]), rel=3e-14, abs=0)
    assert report["energy_variance"] == pytest.approx(float(variance), rel=3e-14, abs=0)
    assert report["dt_dlambda"] == pytest.approx(float(-variance), rel=3e-14, abs=0)
    assert report["ds_dlambda"] == pytest.approx(float(ds), rel=3e-14, abs=0)


def test_zero_multiplier_keeps_original_uniform_variance():
    observable = np.array([[0, .25j, 0], [-.25j, 0, 0], [0, 0, 1.]])
    _, _, variance, _ = partition([-.25, .25, 1.], 0.)
    report = gibbs_entropy_response(observable, 0.)
    assert report["energy_variance"] == pytest.approx(float(variance), rel=3e-15, abs=0)
    assert report["dt_dlambda"] == pytest.approx(float(-variance), rel=3e-15, abs=0)
    assert report["ds_dlambda"] == report["ds_dt"] == 0
    np.testing.assert_array_equal(report["state"], np.eye(3)/3)


@pytest.mark.parametrize("power", [500, 600])
@pytest.mark.parametrize("phase", [1., 1j])
def test_normal_response_does_not_silently_erase_a_genuine_state_component(power, phase):
    tiny = 2.**-power
    observable = np.array([[0., phase*tiny, 0],
                           [np.conj(phase)*tiny, 1., 0], [0, 0, 2.]])
    ctx = mpmath.mp.clone()
    ctx.dps = 850
    t = ctx.mpf(tiny)
    # Exact 2x2 characteristic roots, with the small root rationalized.
    root = ctx.sqrt(1+4*t*t)
    energies = [-2*t*t/(1+root), (1+root)/2, ctx.mpf(2)]
    weights = [ctx.exp(-t*e) for e in energies]
    total = ctx.fsum(weights)
    coherence = t*(weights[1]-weights[0])/(root*total)
    if float(coherence) == 0:
        assert coherence != 0
        with pytest.raises(ValueError, match="component|precision|resolved"):
            gibbs_entropy_response(observable, tiny)
    else:
        mean = ctx.fsum(w*e for w, e in zip(weights, energies))/total
        variance = ctx.fsum(w*(e-mean)**2 for w, e in zip(weights, energies))/total
        report = gibbs_entropy_response(observable, tiny)
        assert report["state"][0, 1] == pytest.approx(phase*float(coherence), rel=3e-14, abs=0)
        assert report["energy_variance"] == pytest.approx(float(variance), rel=3e-14, abs=0)
        assert report["ds_dlambda"] == pytest.approx(float(-t*variance), rel=3e-14, abs=0)


def test_representable_subnormal_coherence_is_not_discarded():
    tiny = np.nextafter(0., 1.)
    observable = np.array([[0., 1j*tiny, 0], [-1j*tiny, 0., 0], [0, 0, 1.]])
    _, probabilities, variance, ds = partition([-tiny, tiny, 1.], 20.)
    expected = 1j*float((probabilities[1]-probabilities[0])/2)
    assert expected != 0
    report = gibbs_entropy_response(observable, 20.)
    assert report["state"][0, 1] == expected
    assert report["energy_variance"] == pytest.approx(float(variance), rel=3e-14, abs=0)
    assert report["ds_dlambda"] == pytest.approx(float(ds), rel=3e-14, abs=0)


def test_structural_complex_zero_survives_spectral_roundoff():
    observable = np.array([[0, 1+1j, 1-1j], [1-1j, 0, 0], [1+1j, 0, 0]])
    beta = .3
    ctx, _, variance, ds = partition([-2., 0., 2.], beta)
    b = ctx.mpf(beta)
    # Direct source identity T^3=4T gives this quadratic exponential;
    # no production spectral data or source-power zero decision is reused.
    z = 1+2*ctx.cosh(2*b)
    expected = (np.eye(3)/float(z)-float(ctx.sinh(2*b)/(2*z))*observable
                +float((ctx.cosh(2*b)-1)/(4*z))*(observable @ observable))
    report = gibbs_entropy_response(observable, beta)
    np.testing.assert_allclose(report["state"], expected, rtol=3e-14, atol=0)
    assert report["state"][2, 1].real == 0
    assert report["state"][2, 1].imag != 0
    assert report["energy_variance"] == pytest.approx(float(variance), rel=3e-14, abs=0)
    assert report["ds_dlambda"] == pytest.approx(float(ds), rel=3e-14, abs=0)


@pytest.mark.parametrize("observable, beta", [
    (np.eye(3)*1e300, .7),
    (np.diag([0., 2.**-60, 1.]), 710.),
    (np.array([[2.**-61, 1j*2.**-61, 0], [-1j*2.**-61, 2.**-61, 0], [0, 0, 1.]]), 700.),
])
def test_scalar_or_unresolved_population_and_dense_support_still_refused(observable, beta):
    with pytest.raises(ValueError):
        gibbs_entropy_response(observable, beta)


@pytest.mark.parametrize("phase", [1., 1j])
@pytest.mark.parametrize("family", ["qubit", "doublet_with_spectator"])
@pytest.mark.parametrize("sign", [-1., 1.])
@pytest.mark.parametrize("middle, offset", [(n, step) for n in (3, 5, 7) for step in (-1, 0, 1)])
def test_subnormal_component_rounding_resolves_both_sides_of_midpoints(
        phase, family, sign, middle, offset):
    tiny = np.nextafter(0., 1.)
    beta = sign*(middle+offset)*tiny
    if family == "qubit":
        observable = np.array([[0, phase], [np.conj(phase), 0]])
        variance = 1.
    else:
        observable = np.array([[0, 1.5*phase, 0],
                               [1.5*np.conj(phase), 0, 0], [0, 0, 1.5]])
        variance = 2.
    ctx = mpmath.mp.clone()
    ctx.dps = 1100
    b = ctx.mpf(beta)
    if family == "qubit":
        coherence = -ctx.tanh(b)/2
        # For b>0, tanh(b)<b: an odd multiple lies just BELOW its
        # half-subnormal midpoint. The cubic correction matters for rounding.
    else:
        w = ctx.expm1(3*b)
        coherence = -w/(2*(3+w))
        # Here |coherence|=b/2+b²/4+O(b³) for b>0: the opposite
        # midpoint side. The source spectrum is exactly {-1.5,1.5,1.5}.
    # Round in integer units of the smallest binary64 subnormal. Direct
    # float(mpf) can itself double-round onto a midpoint, so it cannot be
    # the oracle for the production conversion being challenged here.
    units = abs(coherence)/ctx.mpf(tiny)
    lower = int(ctx.floor(units))
    fraction = units-lower
    rounded_units = lower+(fraction > ctx.mpf('.5')
                           or (fraction == ctx.mpf('.5') and lower % 2))
    expected = phase*int(ctx.sign(coherence))*rounded_units*tiny
    if offset == 0:
        # Agreement of two values already rounded onto a midpoint cannot
        # decide its side. Accurate output or explicit resolution refusal
        # are both honest; silently selecting the midpoint's even neighbor is not.
        try:
            report = gibbs_entropy_response(observable, beta)
        except ValueError as error:
            assert "component" in str(error) and ("rounding" in str(error)
                                                  or "resolved" in str(error)
                                                  or "precision" in str(error))
            return
    else:
        # Adjacent even multipliers put the coherence near a representable
        # subnormal, away from the midpoint. A blanket range refusal fails.
        report = gibbs_entropy_response(observable, beta)
    assert report["state"][0, 1] == expected
    assert report["state"][1, 0] == expected.conjugate()
    assert report["energy_variance"] == variance
    assert report["ds_dlambda"] == -beta*variance
    assert report["ds_dt"] == beta


@pytest.mark.parametrize("phase", [1., 1j])
def test_positive_component_below_zero_midpoint_is_not_accepted_as_zero(phase):
    observable = np.array([[0, phase], [np.conj(phase), 0]])
    with pytest.raises(ValueError, match="component.*(resolved|rounding|precision)"):
        gibbs_entropy_response(observable, np.nextafter(0., 1.))
    report = gibbs_entropy_response(observable, 0.)
    np.testing.assert_array_equal(report["state"], np.eye(2)/2)
    assert report["energy_variance"] == 1.
    assert report["ds_dlambda"] == report["ds_dt"] == 0.


@pytest.mark.parametrize("phase", [1., 1j])
@pytest.mark.parametrize("beta_units", [5, -3])
def test_resolved_subnormal_midpoint_side_is_rounded_once(phase, beta_units):
    tiny, spectator = np.nextafter(0., 1.), 2.**400
    beta = beta_units*tiny
    observable = np.array([[0, 1.5*phase, 0],
                           [1.5*np.conj(phase), 0, 0], [0, 0, spectator]])
    ctx = mpmath.mp.clone()
    ctx.dps = 900
    b, e = ctx.mpf(beta), ctx.mpf(spectator)
    low, high, third = ctx.expm1(ctx.mpf('1.5')*b), ctx.expm1(-ctx.mpf('1.5')*b), ctx.expm1(-e*b)
    total = 3+low+high+third
    coherence = (high-low)/(2*total)
    units = abs(coherence)/ctx.mpf(tiny)
    if beta_units > 0:
        assert ctx.mpf('2.5') < units < 3
        rounded_units = -3
    else:
        assert 1 < units < ctx.mpf('1.5')
        rounded_units = 1
    # The spectator moves the coherence off the midpoint by O(beta^2 E),
    # resolvable at the producer's precision. The exact integer cell decides
    # its binary64 rounding without converting the MP component to float.
    energies, weights = [-ctx.mpf('1.5'), ctx.mpf('1.5'), e], [1+low, 1+high, 1+third]
    mean = ctx.fsum(w*x for w, x in zip(weights, energies))/total
    variance = ctx.fsum(w*(x-mean)**2 for w, x in zip(weights, energies))/total
    report = gibbs_entropy_response(observable, beta)
    assert report["state"][0, 1] == phase*rounded_units*tiny
    assert report["state"][1, 0] == report["state"][0, 1].conjugate()
    assert report["energy_variance"] == pytest.approx(float(variance), rel=3e-14, abs=0)
    assert report["dt_dlambda"] == pytest.approx(float(-variance), rel=3e-14, abs=0)
    assert report["ds_dlambda"] == pytest.approx(float(-b*variance), rel=3e-14, abs=0)
