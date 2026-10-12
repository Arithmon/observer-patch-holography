"""Independent controls for the interacting occupation certificate."""
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import occupation as oc


def test_concrete_refinement_and_exact_characteristic():
    rows = [oc.occupation_cutoff_bound(d) for d in (32, 64, 128)]
    errors = [F(row["pre_final_readout_trace_distance_upper"]) for row in rows]
    assert errors[0] < F(11408, 10**6)
    assert errors[1] < F(2498, 10**9)
    assert errors[2] < F(6175, 10**17)
    assert errors[0] > errors[1] > errors[2] > 0
    for row in rows:
        moment = row["number_moment_certificate"]
        q = F(moment["exponential_weight"])
        growth = F(moment["characteristic_exponential_upper"])
        assert q == F(113, 64)
        assert q >= F(7, 4)
        assert F(moment["number_moment_upper"]) < 13
        assert (q - 1) / (q + 1) * growth <= F(7, 9)
        tail = F(moment["number_moment_upper"]) / q**row["dimension"]
        assert F(row["unclipped_rectangular_outside_probability_upper"]) == tail
        assert not row["final_scalar_sign_readout_comparison_included"]
        assert not row["arbitrary_scalar_input_diamond_bound"]


@pytest.mark.parametrize("time", [F(0), F(1, 10), F(1, 2), F(1), F(3, 2)])
def test_directed_rational_moment_factors(time):
    certificate = oc.number_moment_bound(time)
    q = F(certificate["exponential_weight"])
    initial = F(certificate["initial_number_moment_upper"])
    power = F(certificate["characteristic_power_upper"])
    # Integer powers verify the fractional-power majorant independently
    # of whichever concavity tangent the implementation selects.
    assert initial**2 >= F(80, 17)
    assert power**22 >= (8/q)**25
    assert F(certificate["number_moment_upper"]) == initial*power
    # Check the exponential enclosure against a separate fixed-order
    # Taylor enclosure, tighter than the exported rounding precision.
    x = oc.CHARACTERISTIC_RATE*time
    terms = [x**k/F(math.factorial(k)) for k in range(101)]
    independent_upper = sum(terms, F()) + terms[-1]*x/101/(1-x/102)
    assert F(certificate["characteristic_exponential_upper"]) >= independent_upper


def test_cross_boundary_constant_against_independent_dense_matrices():
    # Build H on a larger Fock square. Every matrix element crossing into
    # P_d is present already at d+2; no infinite operator is substituted.
    for d in (2, 4, 8):
        size = d + 2
        a = np.diag(np.sqrt(np.arange(1, size)), 1)
        identity = np.eye(size)
        a0, a1 = np.kron(identity, a), np.kron(a, identity)
        n0, n1 = a0.T @ a0, a1.T @ a1
        kappa, g = float(oc.KAPPA), float(oc.MAXIMUM_COUPLING)
        h = ((1 + kappa / 2) * (n0 + n1 + np.eye(size * size))
             + kappa / 4 * (a0 @ a0 + a0.T @ a0.T + a1 @ a1 + a1.T @ a1.T)
             - kappa / 2 * (a0 + a0.T) @ (a1 + a1.T)
             + g * (a0 + a0.T))
        inside = [n0 + size * n1 for n1 in range(d) for n0 in range(d)]
        outside = [i for i in range(size * size) if i not in inside]
        actual = np.linalg.norm(h[np.ix_(inside, outside)], 2)
        certified = F(oc.occupation_cutoff_bound(d)["projected_hamiltonian_cross_norm_upper"])
        assert actual <= float(certified)


def test_weighted_commutator_inequality_for_every_compressed_state():
    # Finite number compression retains the exact weighted commutator.
    # A noncommuting matter Hamiltonian is included and drops out because
    # it commutes with N; the force remains fully matter-entangled.
    d = 7
    a = np.diag(np.sqrt(np.arange(1, d)), 1)
    a0, a1 = np.kron(np.eye(d), a), np.kron(a, np.eye(d))
    n = a0.T @ a0 + a1.T @ a1
    kappa, g = float(oc.KAPPA), float(oc.MAXIMUM_COUPLING)
    scalar = ((1 + kappa / 2) * (n + np.eye(d*d))
              + kappa / 4 * (a0 @ a0 + a0.T @ a0.T + a1 @ a1 + a1.T @ a1.T)
              - kappa / 2 * (a0 + a0.T) @ (a1 + a1.T))
    z = np.diag([1., -1.])
    x = np.array([[0., 1.], [1., 0.]])
    h = np.kron(scalar, np.eye(2)) + np.kron(np.eye(d*d), x) + g*np.kron(a0+a0.T, z)
    numbers = np.repeat(np.diag(n), 2)
    for q in (F(5, 4), F(2), F(8)):
        sqrt_weights = float(q)**(numbers/2)
        sinh = (float(q) - 1/float(q))/2
        # Congruence by exp(-sN/2) avoids exponentially large weights.
        # Positive semidefiniteness checks all compressed states, and
        # remains true after tensoring with any untouched reference.
        conjugated = h*sqrt_weights[None, :]/sqrt_weights[:, None]
        derivative = 1j*(conjugated-conjugated.T)
        upper = np.diag(sinh*(float(oc.CHARACTERISTIC_RATE)*numbers
                             + float(oc.MOMENT_RATE)))
        assert np.linalg.eigvalsh(upper-derivative).min() >= -1e-11
        assert np.linalg.eigvalsh(upper+derivative).min() >= -1e-11


def test_initial_tail_and_no_rounding_floor():
    for d in (2, 8, 32):
        tail = oc.preparation_tail(d)
        # Independently sum the omitted positive squeezed-vacuum series.
        first_pair = d//2
        omitted = sum(F(math.comb(2*k, k), 4**k * 81**k)
                      for k in range(first_pair, first_pair + 30))
        # 4*sqrt(5)/9 < 1, so no irrational computation enters this lower bound.
        assert F(99, 100)*omitted <= tail
    assert 0 < oc.preparation_tail(128) < F(1, 10**120)


def test_initial_vacuum_and_actual_compressed_hamiltonian():
    from ..observer.source_scalar import field
    from . import scalar_coupling as model

    assert oc.KAPPA == field.KAPPA
    assert oc.MAXIMUM_COUPLING == F(model.COUPLING)
    for d in (2, 4):
        a = np.diag(np.sqrt(np.arange(1, d)), 1)
        a0, a1 = np.kron(np.eye(d), a), np.kron(a, np.eye(d))
        n = a0.T@a0 + a1.T@a1
        scalar = ((73/64)*(n+np.eye(d*d))
                  + (9/128)*(a0@a0+a0.T@a0.T+a1@a1+a1.T@a1.T)
                  - (9/64)*(a0+a0.T)@(a1+a1.T))
        terms = model.hamiltonian_terms(d)
        expected = (np.kron(scalar, np.eye(4)) + terms["hopping"]
                    + terms["electric"]
                    + np.kron(a0+a0.T, np.diag([1., -1., 1., -1.]))/4)
        assert np.linalg.norm(expected-model.ideal_hamiltonian(d), 2) < 1e-13
        _, preparation = field.vacuum_vector(d)
        actual_tail = 1-preparation["retained_probability"]
        assert abs(actual_tail-float(oc.preparation_tail(d))) < 1e-14


@pytest.mark.parametrize("dimension", [True, 0, 1, 3, 4.0, "4"])
def test_bad_dimensions(dimension):
    with pytest.raises(ValueError):
        oc.occupation_cutoff_bound(dimension)


@pytest.mark.parametrize("time", [True, 0, 1.0, "1", F(-1), F(2)])
def test_time_requires_exact_fraction(time):
    with pytest.raises(ValueError):
        oc.occupation_cutoff_bound(32, time)


@pytest.mark.parametrize("dimensions", [(), [], {}, {2, 4}, "32", [True], [3]])
def test_bad_table_parameters(dimensions):
    with pytest.raises(ValueError):
        oc.certificate_table(dimensions)
