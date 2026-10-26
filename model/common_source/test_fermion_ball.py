"""Independent CAR, Fourier, source and entropy controls for the same field."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest
from scipy.linalg import eigh

from . import fermion_ball as f


def test_actual_full_action_and_exact_local_particle_hole_attachment():
    assert f.exact_attachment()["intertwining_and_CAR_identity_exact"]
    result = f.full_action_attachment()
    assert result["charged_color_family_blocks_checked"] == 21
    assert result["maximum_intertwining_error"] < 2e-12


@pytest.mark.parametrize("mass", [0., 1e-7, 0.8, 9.])
def test_analytic_spinors_solve_the_independent_dirac_eigenproblem(mass):
    p = (0.7, -0.4, 0.2)
    energy, negative, positive = f.spinors(p, mass)
    alpha, beta = f.matrices()
    h = np.einsum("i,ijk->jk", p, alpha)+mass*beta
    frame = np.column_stack((negative, positive))
    assert np.max(abs(frame.conj().T@frame-np.eye(4))) < 2e-14
    assert np.max(abs(h@negative+energy*negative)) < 2e-14
    assert np.max(abs(h@positive-energy*positive)) < 2e-14


def test_fft_restriction_matches_direct_twisted_plane_waves():
    n, length, mass = 4, F(6), 0.3
    labels = np.array([[0, 1, 2], [3, 0, 1], [1, 3, 3]])
    got = f.restricted_covariance(n, length, mass, labels)
    expected = np.zeros_like(got)
    wave = 2*np.pi*(np.fft.fftfreq(n)*n+.5)/float(length)
    for index in np.ndindex(n, n, n):
        p = wave[list(index)]
        _, negative, _ = f.spinors(p, mass)
        phases = np.exp(1j*(float(length)/n*labels)@p)/n**1.5
        columns = (phases[:, None, None]*negative).reshape(4*len(labels), 2)
        expected += columns@columns.conj().T
    assert np.max(abs(got-expected)) < 2e-15
    assert np.max(abs(got-got.conj().T)) < 2e-15
    assert np.array_equal(got, f.restricted_covariance(n, length, mass, labels.astype(np.uint64)))


def test_massless_dirac_region_is_two_equal_entropy_complementary_weyl_blocks():
    labels, _ = f.region(4, F(6), F(3, 2))
    covariance = f.restricted_covariance(4, F(6), 0., labels)
    slots = np.arange(4*len(labels)).reshape(-1, 4)
    left, right = slots[:, :2].ravel(), slots[:, 2:].ravel()
    c_left, c_right = covariance[np.ix_(left, left)], covariance[np.ix_(right, right)]
    tolerance = 64*len(covariance)*np.finfo(float).eps*max(1., np.linalg.norm(covariance, 2))
    assert np.linalg.norm(c_left+c_right-np.eye(len(c_left)), 2) <= tolerance
    assert np.linalg.norm(covariance[np.ix_(left, right)], 2) <= tolerance
    dirac_entropy = f.modular(covariance)[0]
    left_entropy = f.modular(c_left)[0]
    right_entropy = f.modular(c_right)[0]
    entropy_tolerance = tolerance*max(1., abs(dirac_entropy))
    assert abs(dirac_entropy-left_entropy-right_entropy) <= entropy_tolerance
    assert abs(left_entropy-right_entropy) <= entropy_tolerance


def _slater(columns):
    modes, occupied = columns.shape
    state = np.zeros(1 << modes, complex)
    for rows in combinations(range(modes), occupied):
        state[sum(1 << j for j in rows)] = np.linalg.det(columns[list(rows)])
    return state


def _reduce_first_two(state):
    # Little-endian mode labels: retained low two bits, traced high two.
    tensor = state.reshape(4, 4).T
    return tensor@tensor.conj().T


def _entropy(rho):
    values = np.linalg.eigvalsh(rho)
    assert np.min(values) > 0
    return float(-np.dot(values, np.log(values)))


def test_correlation_entropy_and_positive_particle_variation_against_full_fock_trace():
    generator = np.random.default_rng(431)
    raw = generator.normal(size=(4, 4))+1j*generator.normal(size=(4, 4))
    frame, _ = np.linalg.qr(raw)
    base = _slater(frame[:, :2])
    excited = _slater(frame[:, :3])
    reduced = _reduce_first_two(base)
    covariance = (frame[:, :2]@frame[:, :2].conj().T)[:2, :2]
    entropy, modular, _, _ = f.modular(covariance)
    assert abs(entropy-_entropy(reduced)) < 2e-14
    psi = frame[:2, 2]
    derivative = float(np.vdot(psi, modular@psi).real)
    epsilon = 1e-5
    after = (1-epsilon)*reduced+epsilon*_reduce_first_two(excited)
    assert abs((_entropy(after)-entropy)/epsilon-derivative) < 2e-4
    gaussian = covariance+epsilon*np.outer(psi, psi.conj())
    assert abs(f.modular(gaussian)[0]-_entropy(after)) < 2e-14


def test_region_is_actual_partial_trace_not_a_new_isolated_ground_state():
    labels, _ = f.region(12, F(6), F(1))
    covariance = f.restricted_covariance(12, F(6), 0.4, labels)
    entropy, _, gap, _ = f.modular(covariance)
    assert entropy > 40 and gap > 0
    # A newly selected pure regional Slater vacuum fails the modular domain.
    wrong = np.diag([1.]*(len(covariance)//2)+[0.]*(len(covariance)//2))
    with pytest.raises(ValueError, match="purity"):
        f.modular(wrong)


def test_actual_source_prepares_original_all_left_vacuum():
    result = f.source_preparation(mass=0.7)
    assert result["local_particle_hole_phases_and_FSWAP_executed"]
    assert result["all_four_local_Fock_input_columns_executed"]
    assert max(result["numerical_errors"].values()) < 1e-11
    assert result["actual_resource_counts"]["native_pulses"] == 13616
    assert result["actual_resource_counts"]["proper_code_service_events"] == 27232
    assert not result["physical_service_error_certified"]


def test_pair_source_corruption_is_not_replaced_by_target_covariance(monkeypatch):
    actual = f.source.local_gate

    def broken(kind, coefficient=0., duration=0.):
        gate = actual(kind, coefficient, duration)
        if kind == "pairing":
            return replace(gate, orthogonal=np.eye(4))
        return gate

    monkeypatch.setattr(f.source, "local_gate", broken)
    with pytest.raises(ArithmeticError, match="preparation failed"):
        f.source_preparation(mass=0.7)


@pytest.mark.parametrize("mass", [0., 0.8, 1.7])
def test_regional_entropy_first_law_uses_independent_entropy_secant(mass):
    result = f.ball(n=12, mass=mass, secant=True)
    control = result["entropy_first_law_control"]
    value = control["independent_entropy_difference_over_step"]
    assert control["analytic_secant_lower"] <= value+control["floating_secant_admission"]
    assert value <= control["analytic_secant_upper"]+control["floating_secant_admission"]
    assert not result["massive_conformal_identity_asserted"]


@pytest.mark.parametrize("n", [True, 0, 3, 41, 12.0])
def test_grid_rejects_ambiguous_or_unbounded_size(n):
    with pytest.raises(ValueError):
        f.vacuum_grid(n, F(6), 0.3)
    with pytest.raises(ValueError):
        f.region(n, F(6), F(1))


@pytest.mark.parametrize("mass", [True, -1., float("nan"), float("inf"), "1"])
def test_mass_domain_refuses_nonphysical_values(mass):
    with pytest.raises(ValueError):
        f.ball(n=12, mass=mass)


def test_cached_kernel_cannot_be_mutated_or_renormalized():
    kernel, _ = f.vacuum_grid(4, F(6), 0.3)
    with pytest.raises(ValueError):
        kernel.setflags(write=True)
    with pytest.raises(ValueError):
        f.restricted_covariance(4, F(6), 0.3, np.array([[0, 0, 0], [0, 0, 0]]))
    with pytest.raises(ValueError):
        f.modular(np.ma.array(np.eye(2)/2, mask=False))
    with pytest.raises(ValueError):
        f.restricted_covariance(4, F(6), 0.3, np.ma.array([[0, 0, 0]], mask=False))
    f.vacuum_grid(4, F(6), 1.)
    with pytest.raises(ValueError):
        f.vacuum_grid(4, F(6), True)


def test_relative_entropy_monotonicity_uses_the_same_global_state():
    result = f.nested_relative_entropy_control()
    a, b = result["regions"]
    assert 0 < a["relative_entropy"] < b["relative_entropy"]
    assert a["sites"] < b["sites"]


def test_numerical_metadata_and_exact_contracts_cannot_drift():
    with pytest.raises(ValueError):
        f.compare({"floating_secant_admission": 2e-7}, {"floating_secant_admission": 2.01e-7})
    with pytest.raises(ValueError):
        f.compare({"derived": False}, {"derived": 0})
    with pytest.raises(ValueError):
        f.compare({"purity_gap": 1e-4}, {"purity_gap": 1.1e-4})
