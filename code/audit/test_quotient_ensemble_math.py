"""Independent finite witnesses for the quotient-ensemble paper corrections.

These test mathematical distinctions, not a prose substring or a sampled
physical realization. Exact rational controls expose selection and gap errors;
the quantum support example has a closed-form matrix exponential.
"""

from fractions import Fraction as Q
import math

import pytest


def test_settlement_pushforward_is_not_uniform_normal_form_counting():
    normalizer = (0, 1, 1)
    prior = (Q(1, 3),) * 3
    pushforward = tuple(sum(prior[x] for x in range(3) if normalizer[x] == q)
                        for q in range(2))
    quotient_counting = (Q(1, 2), Q(1, 2))
    assert pushforward == (Q(1, 3), Q(2, 3))
    # The indicator of the second normal form is an observable prediction.
    assert pushforward[1] - quotient_counting[1] == Q(1, 6)
    # A lift can reproduce the target, but its fiber masses are extra data.
    matching_prior = (Q(1, 2), Q(1, 4), Q(1, 4))
    assert tuple(sum(matching_prior[x] for x in range(3) if normalizer[x] == q)
                 for q in range(2)) == quotient_counting


def test_conserved_record_annuls_globally_centered_gap():
    # Conditional averaging inside two protected record fibers.
    def expectation(f):
        return ((f[0] + f[1]) / 2,) * 2 + ((f[2] + f[3]) / 2,) * 2

    record = tuple(map(Q, (1, 1, -1, -1)))
    assert sum(record) == 0
    assert expectation(record) == record
    # Thus the infimum over f perpendicular to constants is zero.
    # On the full fixed-space complement, I-E is exactly the identity.
    for f in ((Q(1), Q(-1), Q(0), Q(0)), (Q(0), Q(0), Q(1), Q(-1))):
        ef = expectation(f)
        assert ef == (0, 0, 0, 0)
        assert sum((x - y) ** 2 for x, y in zip(f, ef)) / sum(x*x for x in f) == 1


def test_boundary_quantum_maxent_compresses_log_reference():
    # H couples basis vectors 0 and 2, leaving 1 fixed. sigma=exp(H)/Z.
    # Restrict feasible states to P=diag(1,1,0). P H P=0, so the minimizer
    # of D(rho||sigma) is I_P/2. Compressing sigma first gives a wrong bias.
    z = 2 * math.cosh(1) + 1
    wrong_p = math.cosh(1) / (math.cosh(1) + 1)

    def relative_entropy(p):
        return p * math.log(p) + (1-p) * math.log(1-p) + math.log(z)

    assert wrong_p > 0.6
    assert relative_entropy(wrong_p) > relative_entropy(0.5) + 0.02
    # Strict convexity gives the unique optimum on this diagonal family;
    # off-diagonal coherence only decreases entropy at fixed diagonal.
    for p in (0.1, 0.25, 0.75, 0.9):
        assert relative_entropy(p) > relative_entropy(0.5)


@pytest.mark.parametrize("theta", [-10, -1, 0, 1, 10])
def test_finite_quantum_multiplier_cannot_reach_pure_boundary(theta):
    # Constraint Tr(rho |1><1|)=0 forces a pure state. No finite Gibbs
    # multiplier with a faithful reference can realize that constraint.
    excited_probability = 1 / (1 + math.exp(theta))
    assert 0 < excited_probability < 1
