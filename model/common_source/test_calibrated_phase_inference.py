"""Exact rational source controls, independent of every matter packet."""
from fractions import Fraction as F
import math
import random

import pytest

from . import calibrated_phase_inference as p


def _sincos(theta):
    """Independent alternating-series enclosure on |theta|<=1."""
    assert abs(theta) <= 1
    x = abs(theta)
    cosine = sum(((-1)**k*x**(2*k)/math.factorial(2*k) for k in range(24)), F())
    sine = sum(((-1)**k*x**(2*k+1)/math.factorial(2*k+1) for k in range(24)), F())
    cb = cosine, cosine+x**48/math.factorial(48)
    sb = sine, sine+x**49/math.factorial(49)
    return cb, sb if theta >= 0 else (-sb[1], -sb[0])


def _observations(phase, prior, level, amplitude, noise=F(1, 10**12)):
    offset = sum(prior, F())*2**level/2
    c, s = _sincos(2**level*phase-offset)
    def expand(box):
        lo, hi = amplitude*box[0]-noise, amplitude*box[1]+noise
        scale = 2**80
        return F(math.floor(lo*scale), scale), F(math.ceil(hi*scale), scale)
    return expand(c), expand(s)


@pytest.mark.parametrize("phase", [F(-3, 7), F(0), F(1, 10**9), F(2, 9), F(10000)+F(5, 17)])
def test_generic_adaptive_phase_coverage_tiny_and_thousand_turns(phase):
    # The initial window supplies only a broad branch, not a target fitted
    # at each round. Every later offset depends solely on the last box.
    center = F(10000) if phase > 1000 else F()
    prior = center-F(1), center+F(1)
    rng = random.Random(691)
    for level in range(25):
        visibility = F(rng.randrange(1, 1001), 1000)
        cosine, sine = _observations(phase, prior, level, visibility,
                                     F(1, 10**12*65536**level))
        update = p.phase_update(prior, level, cosine, sine)
        assert update["offset_exact"] == sum(prior, F())*2**level/2
        assert update["certified"], update["reason"]
        new = update["phase_base_interval"]
        assert prior[0] <= new[0] <= phase <= new[1] <= prior[1]
        assert new[1]-new[0] < prior[1]-prior[0]
        assert update["source_phase_interval"][0] <= 2**level*phase-update["offset_exact"] <= update["source_phase_interval"][1]
        assert not update["unknown_amplitude_fitted"]
        prior = new
    assert prior[1]-prior[0] < F(1, 10**15)


def test_ratio_cancels_random_unknown_amplitude_without_a_packet():
    prior, phase = (F(-1), F(1)), F(1, 3)
    for visibility in (F(1, 10000), F(1, 2), F(3), F(1000)):
        c, s = _observations(phase, prior, 0, visibility, F(1, 10**20))
        out = p.phase_update(prior, 0, c, s)
        assert out["certified"]
        assert out["phase_base_interval"][0] <= phase <= out["phase_base_interval"][1]
    assert "calibrated_matter" not in p.__dict__


def test_large_exact_record_denominators_are_outward_bounded_before_atan():
    prior, phase = (F(-1), F(1)), F(1, 3)
    cosine, sine = _observations(phase, prior, 0, F(3, 7))
    tiny = F(1, 2**10000+1)
    cosine = cosine[0]-tiny, cosine[1]+tiny
    sine = sine[0]-tiny, sine[1]+tiny
    out = p.phase_update(prior, 0, cosine, sine)
    assert out["certified"]
    assert out["phase_base_interval"][0] <= phase <= out["phase_base_interval"][1]
    low, high = out["sine_cosine_ratio_interval"]
    corners = [y/x for y in sine for x in cosine]
    assert low <= min(corners) <= max(corners) <= high
    assert low.denominator <= 2**80 and high.denominator <= 2**80


@pytest.mark.parametrize("cosine,sine", [((F(-1), F(1)), (F(-1), F(1))),
    ((F(0), F(0)), (F(0), F(0))), ((F(-2), F(-1)), (F(0), F(1)))])
def test_noise_zero_visibility_and_wrong_orientation_unresolved(cosine, sine):
    prior = F(-1), F(1)
    out = p.phase_update(prior, 0, cosine, sine)
    assert not out["certified"]
    assert out["phase_base_interval"] == prior
    assert "positive cosine" in out["reason"]


def test_positive_cosine_but_uninformative_noise_is_not_certification():
    prior = F(-1, 10), F(1, 10)
    out = p.phase_update(prior, 0, (F(1, 10), F(1)), (F(-1), F(1)))
    assert not out["certified"]
    assert "does not contract" in out["reason"]
    assert out["phase_base_interval"] == prior


def test_inconsistent_quadratures_and_alias_window_refuse():
    prior = F(1), F(11, 10)
    # Enacted offset1.05, but atan(1)=pi/4 is outside residual[-.05,.05].
    out = p.phase_update(prior, 0, (F(1), F(1)), (F(1), F(1)))
    assert not out["certified"]
    assert "inconsistent" in out["reason"]
    out = p.phase_update((F(-2), F(2)), 0, (F(1), F(1)), (F(0), F(0)))
    assert not out["certified"]
    assert "right half-plane" in out["reason"]
    # An independently narrow phase window also becomes aliased if a
    # dyadic gain is chosen without respecting its inherited width.
    out = p.phase_update((F(-1), F(1)), 1, (F(1), F(1)), (F(0), F(0)))
    assert not out["certified"]


@pytest.mark.parametrize("bad", [(0, F(1)), (F(0), 1.), [F(0), F(1)],
                                  (F(1), F(0)), (F(0), F(0)), (False, F(1))])
def test_original_prior_fraction_domain(bad):
    with pytest.raises(ValueError):
        p.phase_update(bad, 0, (F(1), F(1)), (F(0), F(0)))


@pytest.mark.parametrize("bad", [False, -1, 1., 257])
def test_original_gain_level_domain(bad):
    with pytest.raises(ValueError):
        p.phase_update((F(-1), F(1)), bad, (F(1), F(1)), (F(0), F(0)))


@pytest.mark.parametrize("bad", [[F(0), F(1)], (0., F(1)), (F(2), F(1)), (True, F(1))])
def test_quadrature_original_interval_domain(bad):
    with pytest.raises(ValueError):
        p.phase_update((F(-1), F(1)), 0, bad, (F(0), F(0)))
