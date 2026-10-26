"""Packet-independent exact interval inversion for adaptive dyadic phases.

The declared source means are C=A*cos(2**level*phi-offset) and
S=A*sin(2**level*phi-offset), with one UNKNOWN common A>0. The enacted
offset is the midpoint of the previous base-phase interval, multiplied by
2**level. The prior must already place every residual phase strictly in
the right half-plane. A positive-cosine confidence rectangle then supplies
an oriented atan interval; the common amplitude cancels without a fit.

This is a deterministic coverage implication conditional on both supplied
mean intervals and that source model. It supplies no sampling/confidence
claim, target phase, calibrated mass, source execution, or authenticity of
records. The caller must charge all trials and establish the source and
joint interval error budget, including adaptive conditioning.
"""
from __future__ import annotations

from fractions import Fraction as F

from .halfline import pi_interval
from ..native_modular_reference.phase_inference import atan_interval

RATIO_ROUND_BITS = 80


def _interval(value, *, prior=False):
    if (type(value) is not tuple or len(value) != 2
            or any(type(x) is not F or type(x.numerator) is not int
                   or type(x.denominator) is not int or x.denominator <= 0 for x in value)):
        raise ValueError("original two-endpoint Fraction intervals required")
    low, high = value
    if low > high or (prior and low == high):
        raise ValueError("ordered intervals and a positive-width prior required")
    return low, high


def phase_update(prior, level, cosine, sine):
    """Return an exact nested phase box, or a reason no update is certified.

    phase_base_interval is the prior on refusal; a caller must check
    certified before treating it as an informative new observation. No
    hidden calibration packet or point estimator enters this function.
    """
    prior = _interval(prior, prior=True)
    cosine, sine = _interval(cosine), _interval(sine)
    if type(level) is not int or not 0 <= level <= 256:
        raise ValueError("an original dyadic level integer in0..256 required")
    gain = F(2**level)
    offset = sum(prior, F())*gain/2
    radius = (prior[1]-prior[0])*gain/2
    result = {"schema": "oph.calibrated_phase_inference.v1", "certified": False,
        "status": "unresolved", "reason": "", "dyadic_level": level,
        "phase_gain_exact": gain, "offset_exact": offset,
        "prior_phase_base_interval": prior, "phase_base_interval": prior,
        "residual_source_phase_window": (-radius, radius),
        "cosine_mean_interval": cosine, "sine_mean_interval": sine,
        "common_positive_amplitude_required": True, "unknown_amplitude_fitted": False,
        "sampling_coverage_or_source_execution_proved_here": False}
    if radius >= pi_interval()[0]/2:
        result["reason"] = "prior residual window is not certified in one right half-plane"
        return result
    if cosine[0] <= 0:
        result["reason"] = "quadrature rectangle does not certify positive cosine/visibility"
        return result
    corners = tuple(y/x for y in sine for x in cosine)
    scale = 2**RATIO_ROUND_BITS
    low, high = min(corners)*scale, max(corners)*scale
    # Exact outward rounding bounds the size of atan's rational powers,
    # even when complete tensor-output norms have enormous denominators.
    ratio = F(low.numerator//low.denominator, scale), F(-(-high.numerator//high.denominator), scale)
    residual = atan_interval(ratio[0])[0], atan_interval(ratio[1])[1]
    candidate = (offset+residual[0])/gain, (offset+residual[1])/gain
    clipped = max(candidate[0], prior[0]), min(candidate[1], prior[1])
    result.update(sine_cosine_ratio_interval=ratio, source_phase_interval=residual,
                  candidate_phase_base_interval=candidate, ratio_outward_round_bits=RATIO_ROUND_BITS)
    if clipped[0] > clipped[1]:
        result["reason"] = "observed quadratures are inconsistent with the declared prior/source phase"
        return result
    if clipped == prior:
        result["reason"] = "quadrature uncertainty does not contract the declared prior"
        return result
    result.update(certified=True, status="certified_interval", reason="oriented exact atan intersection",
                  phase_base_interval=clipped)
    return result
