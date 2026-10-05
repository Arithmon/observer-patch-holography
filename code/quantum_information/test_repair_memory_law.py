"""Exact maps and independent evolution controls for the false-memory witness.

The analytic spectrum and finite-time response are derived in
REPAIR_FIXED_POINTS.md. These controls construct the declared rational
dephasing maps independently of the GNS-projection implementation.
"""

from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm

from quantum_information.algebras import FiniteAlgebra, algebra_intersection
from quantum_information.expectations import apply_map, choi_matrix, repair_generator


I = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.diag(1, -1)


def exact_maps(epsilon):
    n = Z+epsilon*X
    return ((sp.eye(4)+sp.kronecker_product(Z, Z))/2,
            (sp.eye(4)+sp.kronecker_product(n, n)/(1+epsilon**2))/2)


@pytest.mark.parametrize("epsilon", (sp.Rational(0), sp.Rational(3, 4),
                                     sp.Rational(1, 10**14), sp.Rational(1, 10**100)))
@pytest.mark.parametrize("rates", ((sp.Rational(1, 2), sp.Rational(1, 2)),
                                   (sp.Rational(2, 7), sp.Rational(3, 5))))
def test_full_exact_channel_has_the_claimed_kernel_and_complete_spectrum(epsilon, rates):
    p, q = exact_maps(epsilon)
    for channel in (p, q):
        assert channel*channel == channel
        assert channel == channel.conjugate().T
        assert channel*sp.Matrix([1, 0, 0, 1]) == sp.Matrix([1, 0, 0, 1])
        # Independently assemble Choi blocks from the action on all units.
        choi = sp.zeros(4)
        for i in range(2):
            for j in range(2):
                unit = sp.zeros(4, 1); unit[2*i+j] = 1
                choi[2*i:2*i+2, 2*j:2*j+2] = (channel*unit).reshape(2, 2)
        assert choi == choi.conjugate().T
        assert choi.eigenvals() == {sp.Integer(0): 2, sp.Integer(1): 2}
    a, b = rates
    decay = a*(sp.eye(4)-p)+b*(sp.eye(4)-q)
    t = sp.Symbol('t')
    expected = t*(t-a-b)*(t**2-(a+b)*t+a*b*epsilon**2/(1+epsilon**2))
    assert sp.expand(decay.charpoly(t).as_expr()-expected) == 0
    assert len(decay.nullspace()) == (2 if epsilon == 0 else 1)


@pytest.mark.parametrize("epsilon", (0., .125, .75, 2., 16.))
@pytest.mark.parametrize("rates", ((.5, .5), (.25, 2.)))
def test_record_evolution_matches_bloch_law_and_full_quantum_semigroup(epsilon, rates):
    e = sp.Rational(epsilon)
    p, q = (np.array(s).astype(complex) for s in exact_maps(e))
    identity, x, z = (np.array(s).astype(complex) for s in (I, X, Z))
    algebras = [FiniteAlgebra([identity, z]), FiniteAlgebra([identity, z+epsilon*x])]
    a, b = rates
    generator, report = repair_generator([p, q], algebras, identity/2, [a, b])
    gamma = a+b
    determinant = a*b*epsilon**2/(1+epsilon**2)
    root = np.sqrt((a-b)**2+4*a*b/(1+epsilon**2))
    slow, fast = 2*determinant/(gamma+root), (gamma+root)/2
    # If e=0 the zero slow mode is protected; the gap on its complement
    # is gamma. If e!=0 that same mode belongs in the decay spectrum.
    assert report['gap'] == pytest.approx(gamma if epsilon == 0 else slow, rel=1e-11)
    assert report['intersection_dimension'] == (2 if epsilon == 0 else 1)
    bzz = b*epsilon**2/(1+epsilon**2)
    slow_weight, fast_weight = (fast-bzz)/root, (bzz-slow)/root
    initial = (identity+.5*z)/2
    for time in (0., .25, 3.):
        channel = expm(time*generator)
        assert np.linalg.eigvalsh(choi_matrix(channel, 2))[0] >= -2e-14
        evolved = apply_map(channel, initial, dual=True)
        expected = .5*(slow_weight*np.exp(-slow*time)+fast_weight*np.exp(-fast*time))
        assert np.trace(evolved@z).real == pytest.approx(expected, abs=2e-14)
        assert np.trace(evolved).real == pytest.approx(1., abs=2e-14)
        assert np.linalg.eigvalsh(evolved)[0] > 0


@pytest.mark.parametrize("epsilon", ('1e-14', '1e-100', '1e-300'))
def test_nearby_repairs_keep_a_record_temporarily_but_erase_it_on_the_true_timescale(epsilon):
    ctx = mp.mp.clone(); ctx.dps = 80
    e = ctx.mpf(epsilon)
    length = ctx.sqrt(1+e*e)
    # Rationalized eigenvalue: no subtraction of two numbers near one.
    slow = e*e/(2*length*(length+1))
    fast = 1-slow
    slow_weight = (1+1/length)/2
    fast_weight = e*e/(2*length*(length+1))
    # At ordinary time, evaluate lost information with expm1 instead of
    # subtracting nearly equal record values. At the true relaxation time
    # the surviving Z signal is order one below its initial value.
    lost_at_one = (-slow_weight*ctx.expm1(-slow)
                   - fast_weight*ctx.expm1(-fast))
    at_one_lifetime = slow_weight*ctx.exp(-1)+fast_weight*ctx.exp(-fast/slow)
    assert lost_at_one > 0
    assert abs(lost_at_one/(e*e)-(2-ctx.exp(-1))/4) < ctx.mpf('1e-25')
    assert abs(at_one_lifetime-ctx.exp(-1)) < ctx.mpf('1e-25')
    assert abs(slow/(e*e)-ctx.mpf('.25')) < ctx.mpf('1e-25')
    # The source entries are still distinct finite binary64 data at 1e-300.
    identity = np.eye(2); x = np.array([[0., 1.], [1., 0.]]); z = np.diag([1., -1.])
    assert Fraction(float(e)) > 0
    assert algebra_intersection([FiniteAlgebra([identity, z]),
                                 FiniteAlgebra([identity, z+float(e)*x])]).dimension == 1


def test_dirichlet_record_bound_holds_for_noncommuting_repairs():
    # A complete independent 3x3 Bloch evolution checks the bound on every
    # coordinate combination tested, including a nearly preserved record.
    for epsilon in (.01, .5, 2.):
        n = np.array([epsilon, 0., 1.])/np.hypot(1., epsilon)
        z = np.array([0., 0., 1.])
        decay = .25*(np.eye(3)-np.outer(z, z))+2*(np.eye(3)-np.outer(n, n))
        for value in (z, np.array([.3, .2, -.7]), n):
            energy = value@decay@value
            assert energy > 0
            for time in (.01, .5, 10.):
                evolved = expm(-time*decay)@value
                assert np.linalg.norm(evolved-value)**2 <= time*energy+2e-14
                assert value@value-value@evolved <= time*energy+2e-14
