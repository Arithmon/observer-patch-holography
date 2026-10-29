"""A sharper full-field spectral/clock bound from the contact-zero compression.

Unlike the seven-state control, this applies to the COMPLETE bounded Hilbert
space whenever H=T+U P, T=T*, P=P*=P^2, PTP=0 and ||T||^2<=b_squared.
These operator hypotheses are premises, not verified from scalar inputs.
"""
from fractions import Fraction as F


def enclosure(b_squared, interaction, negative_magnitude, positive_edge):
    """Certify rational outward edges and the centered generator norm bound.

The lower edge is -b. The upper edge is no larger than
(U+sqrt(U^2+4b^2))/2. The supplied rational bounds are checked without
floating square roots; the return values do not round the certificate.
"""
    if any(type(v) is not F for v in
           (b_squared, interaction, negative_magnitude, positive_edge)):
        raise ValueError('exact Fraction bounds required')
    if b_squared < 0 or interaction <= 0 or negative_magnitude < 0:
        raise ValueError('nonnegative squared norm/negative magnitude and positive U required')
    if negative_magnitude**2 < b_squared:
        raise ValueError('negative edge fails to enclose the free spectrum')
    if (positive_edge < interaction or positive_edge**2 < b_squared
            or positive_edge*(positive_edge-interaction) < b_squared):
        raise ValueError('positive edge fails the contact resolvent bound')
    return {'lower': -negative_magnitude, 'upper': positive_edge,
            'center': (positive_edge-negative_magnitude)/2,
            'centered_norm_upper': (positive_edge+negative_magnitude)/2}


def canonical_clock_bound(parameter=F(1), beta_error=F(3,80)):
    """Optional exact rate-error bound for the complete canonical N=2 field.

For unitary channels the half-diamond distance is bounded by the centered
operator difference, uniformly with untouched references. This covers ONLY
the rate-error contribution. Statistical failure, splitting and source errors
remain separate, as in band_clock; no existing retained receipt is changed.
"""
    if type(parameter) is not F or abs(parameter)>1:
        raise ValueError('exact Fraction modular parameter with |s|<=1 required')
    if type(beta_error) is not F or beta_error < 0:
        raise ValueError('nonnegative exact Fraction beta error required')
    result=enclosure(F(3,16),F(1),F(7,16),F(7,6))
    return {**result, 'rate_error_halfdiamond_upper': min(
        F(1),result['centered_norm_upper']*abs(parameter)*beta_error)}
