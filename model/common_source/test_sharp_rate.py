"""Exact independent checks of the vacuum finite-code resource theorem.

These tests use the vacuum only as an analytic control. They do not run an
interacting history or infer its convergence rate from this control.
"""
from fractions import Fraction as F
import math

import pytest
import sympy as sp

from . import halfline, sharp_output


def _central(k):
    return F(math.comb(2*k, k), 4**k)


def test_recurrence_and_symbolic_induction_certify_every_integer_index():
    k = sp.symbols("k", integer=True, positive=True)
    # These nonnegative polynomial remainders prove the induction steps,
    # rather than testing the desired inequalities only at sample indices.
    assert sp.expand((2*k+1)**2-4*k*(k+1)) == 1
    assert sp.expand(4*(k+1)**2*(3*k+1)-(2*k+1)**2*(3*k+4)) == k
    assert _central(1)**2 == F(1, 4)  # Lower-bound induction base.
    assert _central(0)**2 == 1       # Upper-bound induction base.
    for n in range(65):
        assert _central(n+1) == _central(n)*F(2*n+1, 2*n+2)
        assert _central(n)-_central(n+1) == _central(n)/F(2*n+2)
        assert _central(n)**2 <= F(1, 3*n+1)
        if n:
            assert _central(n)**2 >= F(1, 4*n)


@pytest.mark.parametrize("m", (1, 2, 3, 8, 32, 128))
def test_finite_tail_telescopes_with_the_required_uniform_ratio(m):
    end = m+37
    differences = sum((_central(k)-_central(k+1) for k in range(m, end)), F(0))
    tail = sum((_central(k)/F(2*k+1) for k in range(m, end)), F(0))
    ratio = F(2*m+2, 2*m+1)
    assert differences == _central(m)-_central(end)
    assert differences <= tail <= ratio*differences
    for k in range(m, end):
        assert 1 <= F(2*k+2, 2*k+1) <= ratio
    # c_end -> 0 follows from the proved c_end^2 <= 1/(3 end+1).


@pytest.mark.parametrize("m", (1, 2, 4, 8, 32, 128))
def test_vacuum_tail_witness_and_constructive_surrogate_have_certified_bounds(m):
    vacuum = {(0, 0, 0): (F(1), F(0))}
    alpha_low, alpha_high = halfline.inverse_sqrt_two_pi_interval()
    pi_low, pi_high = halfline.pi_interval()
    series = sum((_central(k)/F(2*k+1) for k in range(m)), F(0))
    tail = (F(1, 4)-series*alpha_high**2,
            F(1, 4)-series*alpha_low**2)
    ratio = F(2*m+2, 2*m+1)
    assert _central(m)*alpha_low**2 <= tail[0] <= tail[1]
    assert tail[1] <= ratio*_central(m)*alpha_high**2
    assert 0 < tail[0] <= tail[1] < F(1, 8)

    root_m = sharp_output.root_interval((F(m), F(m)))
    root_3m = sharp_output.root_interval((F(3*m), F(3*m)))
    assert F(1)/(4*pi_high*root_m[1]) <= tail[0]
    assert tail[1] <= F(2)/(3*pi_low*root_3m[0])

    for outcome in (0, 1):
        masses = sharp_output.projected_masses(vacuum, 2*m, outcome)
        assert masses["total"] == (F(1, 2), F(0), F(0))
        assert masses["inside"] == (F(1, 4), F(0), series)
        assert masses["outside"] == (F(1, 4), F(0), -series)

    witness = sharp_output.support_witness(vacuum, 2*m, 1)
    witness_interval = tuple(2*x for x in witness["half_trace_distance_witness_interval"])
    # L^2=4*(1/4-series*alpha^2)*(1/4+series*alpha^2).
    analytic = sharp_output.root_interval((F(1, 4)-4*series**2*alpha_high**4,
                                            F(1, 4)-4*series**2*alpha_low**4))
    assert max(analytic[0], witness_interval[0]) <= min(analytic[1], witness_interval[1])
    assert witness_interval[0]**2 >= F(3)/(8*pi_high*root_m[1])

    # The surrogate restores each projected branch's original mass 1/2.
    # Its squared overlap with the true branch is a/2, hence the complete
    # two-block distance is sqrt(1-2*a)=sqrt(2*T), not the finite B_d law.
    surrogate = sharp_output.root_interval((2*tail[0], 2*tail[1]))
    assert witness_interval[1] < surrogate[0]
    assert surrogate[1]**2 <= F(4)/(3*pi_low*root_3m[0])


def test_surrogate_rank_one_identity_and_epsilon_fourth_power_obligation():
    t = sp.symbols("t", real=True)
    inside, branch_mass = sp.Rational(1, 2)-t, sp.Rational(1, 2)
    overlap_squared = branch_mass*inside
    assert sp.expand((2*branch_mass)**2-4*overlap_squared) == 2*t
    # L >= sqrt(3/(8*pi))*m^(-1/4), d=2m, implies
    # d*L^4 >= 9/(32*pi^2). Check the algebra before applying a tolerance.
    m, pi = sp.symbols("m pi", positive=True)
    assert sp.simplify(2*m*(3/(8*pi*sp.sqrt(m)))**2) == 9/(32*pi**2)
    _, pi_upper = halfline.pi_interval()
    for tolerance in (F(1, 10), F(1, 100), F(1, 1000)):
        required = F(9)/(32*pi_upper**2*tolerance**4)
        assert required > 0
        assert required*tolerance**4 == F(9)/(32*pi_upper**2)
