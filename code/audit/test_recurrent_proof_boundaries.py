"""Counterexamples to missing assumptions in the neural interpretation paper."""

from fractions import Fraction as Q
import itertools


def test_approximate_descent_allows_nonstationary_periodic_limit_points():
    # On r=1, Phi(r,theta)=sin(theta)/2+(r-1)*sqrt(2-x/2-x*x/4),
    # x=cos(theta), and the smooth flow is theta'=1, r'=0.
    # Its tangential and radial gradient components are orthogonal.
    for x in (Q(-1), Q(-3, 5), Q(0), Q(3, 5), Q(1)):
        tangent_squared = x*x/4
        radial_squared = 2-x/2-x*x/4
        gradient_squared = tangent_squared + radial_squared
        derivative = x/2
        assert radial_squared > 0
        assert derivative == -gradient_squared + 2
    # At theta=pi each revolution, the gradient exceeds epsilon/c=2.
    assert Q(5, 2) > 2


def test_energy_nonincrease_does_not_supply_scheduler_fairness():
    # Independent spins with positive fields: (-1) at coordinate 1 is
    # unstable forever if the scheduler only revisits coordinate 0.
    state = (1, -1)
    def update(s, i):
        return tuple(1 if j == i else value for j, value in enumerate(s))
    assert update(state, 0) == state
    assert update(state, 1) != state


def test_tie_preserving_hopfield_flips_strictly_lower_energy():
    weights = ((0, 2, -1), (2, 0, 1), (-1, 1, 0))
    thresholds = (1, 0, -1)
    def energy(s):
        return -sum(weights[i][j]*s[i]*s[j] for i in range(3)
                    for j in range(i)) + sum(t*x for t, x in zip(thresholds, s))
    saw_tie = False
    for s in itertools.product((-1, 1), repeat=3):
        for i in range(3):
            field = sum(w*x for w, x in zip(weights[i], s)) - thresholds[i]
            spin = s[i] if field == 0 else (1 if field > 0 else -1)
            after = s[:i] + (spin,) + s[i+1:]
            if field == 0:
                saw_tie = True
                assert after == s
            if after != s:
                assert energy(after) - energy(s) == -2*abs(field) < 0
    assert saw_tie


def test_equilibrium_propagation_sign_on_exact_quadratic_branch():
    # E=(s-theta)^2/2, C=(s-y)^2/2, s_beta=(theta+beta*y)/(1+beta).
    theta, y = Q(3), Q(1)
    for beta in (Q(1, 2), Q(1, 10), Q(1, 100)):
        s_beta = (theta+beta*y)/(1+beta)
        contrast = (theta-s_beta)/beta
        assert contrast == (theta-y)/(1+beta)
        assert 0 < contrast < theta-y
        # Finite-nudge contrast is biased; the limit has the positive sign.
        assert (theta-y)-contrast == (theta-y)*beta/(1+beta)
