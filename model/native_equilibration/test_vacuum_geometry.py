"""Native spectrum, actual partial trace, shape sign and local-stress controls."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
import scipy.sparse as sparse
from scipy.sparse.linalg import eigsh
import sympy as s

from . import vacuum_geometry as v
from ..observer import native, quantum_mixing


def _entropy(position, momentum):
    values, vectors = np.linalg.eigh(position)
    root = (vectors*np.sqrt(values))@vectors.T
    symplectic = np.sqrt(np.linalg.eigvalsh(root@momentum@root))
    assert min(symplectic) > .5
    occupation = symplectic-.5
    return float(np.sum(np.log1p(occupation)+occupation*np.log1p(1/occupation)))


@pytest.fixture(scope="module")
def report():
    return v.run()


def test_actual_native_endpoint_mean_supplies_the_entire_angular_operator():
    lap = np.asarray(v.native_angular(), dtype=int)
    states, seams = native.states(binary=True), quantum_mixing.seams()
    assert len(states) == 66 and len(seams) == 30
    reconstructed = np.zeros((12, 12), dtype=int)
    for a, b in seams:
        reconstructed[a, a] += 1
        reconstructed[b, b] += 1
        reconstructed[a, b] -= 1
        reconstructed[b, a] -= 1
    np.testing.assert_array_equal(lap, reconstructed)
    for state in states:
        total = np.zeros(12, dtype=int)
        for edge in seams:
            for coin in (0, 1):
                total += np.asarray(native.native(state, edge, coin)[0])
        np.testing.assert_array_equal(total, (60*np.eye(12, dtype=int)-lap)@state)
    expected = [0.]+[5-math.sqrt(5)]*3+[6.]*5+[5+math.sqrt(5)]*3
    np.testing.assert_allclose(np.linalg.eigvalsh(lap), expected, atol=5e-15, rtol=0)


def test_symbolic_monotonicity_and_independent_geometry_determinants(report):
    identities = v.exact_identities()
    assert identities["identities_use_independent_metric_determinants"]
    assert identities["two_face_area_log_rate"] == "-2*u/3 - 2*v/3"
    assert identities["proper_radial_lattice_spacing_log_rate"] == "2*u/3 + 2*v/3"
    assert identities["proper_tangential_lattice_spacing_log_rate"] == "-u/3 - v/3"
    assert identities["proper_regulator_aspect_ratio_log_rate"] == "u + v"
    assert report["fixed_slab_proper_volume_derivative_exact"] == "0"
    assert report["regional_vacuum_entropy_derivative"] < 0
    assert report["two_face_area_derivative"] < 0
    assert sum(row["multiplicity"] for row in report["angular_modes"]) == 12
    assert report["angular_modes"][0]["entropy_derivative"] == 0
    assert all(row["entropy_derivative"] < 0 for row in report["angular_modes"][1:])
    for newton in (F(1, 100), F(1), F(100)):
        derivative = report["regional_vacuum_entropy_derivative"]+report["two_face_area_derivative"]/(4*float(newton))
        assert derivative < 0


def test_gaussian_reduction_matches_independent_integration_of_exterior_positions(report):
    cov = v.covariance()
    potential = cov["potential"]
    values, vectors = np.linalg.eigh(potential)
    ground_kernel = (vectors*np.sqrt(values))@vectors.T
    middle = np.arange(12, 24)
    outside = np.r_[0:12, 24:36]
    aa = ground_kernel[np.ix_(middle, middle)]
    ab = ground_kernel[np.ix_(middle, outside)]
    bb = ground_kernel[np.ix_(outside, outside)]
    # Direct Gaussian integration of psi(x,y) psi(x',y), not an assumed
    # reduced thermal covariance or a hidden measurement of exterior y.
    traced_position = np.linalg.inv(aa-ab@np.linalg.solve(bb, ab.T))/2
    traced_momentum = aa/2
    np.testing.assert_allclose(traced_position, cov["regional_position"], atol=2e-15, rtol=0)
    np.testing.assert_allclose(traced_momentum, cov["regional_momentum"], atol=2e-14, rtol=0)
    assert abs(_entropy(traced_position, traced_momentum)-report["regional_vacuum_entropy_nats"]) < 3e-13
    # Setting the unobserved exterior coordinate to zero would instead
    # prepare a pure conditional wavefunction and erase the real entropy.
    conditioned_position = np.linalg.inv(aa)/2
    np.testing.assert_allclose(conditioned_position@traced_momentum, np.eye(12)/4, atol=2e-15, rtol=0)
    assert np.linalg.norm(conditioned_position-traced_position) > .01


def _fock_partial_trace(native_eigenvalue, cutoff):
    kappa = 4.
    local_frequency = math.sqrt(2*kappa+native_eigenvalue)
    annihilate = sparse.diags(np.sqrt(np.arange(1, cutoff)), 1, shape=(cutoff, cutoff), format="csr")
    local_coordinate = annihilate+annihilate.T
    identity = sparse.eye(cutoff, format="csr")
    number = sparse.diags(np.arange(cutoff)+.5, format="csr")
    kron = lambda a, b, c: sparse.kron(sparse.kron(a, b, format="csr"), c, format="csr")
    hamiltonian = local_frequency*(kron(number, identity, identity)
        +kron(identity, number, identity)+kron(identity, identity, number))
    hamiltonian -= kappa/(2*local_frequency)*(kron(local_coordinate, local_coordinate, identity)
        +kron(identity, local_coordinate, local_coordinate))
    _, vector = eigsh(hamiltonian, k=1, which="SA", tol=1e-13,
                     v0=np.ones(cutoff**3)/math.sqrt(cutoff**3))
    wave = vector[:, 0].reshape(cutoff, cutoff, cutoff)
    reduced = np.einsum("abc,adc->bd", wave, wave.conj())
    eigenvalues = np.linalg.eigvalsh(reduced)
    assert min(eigenvalues) > -1e-14 and abs(np.trace(reduced)-1) < 1e-13
    return -math.fsum(float(x)*math.log(float(x)) for x in eigenvalues if x > 0)


@pytest.mark.parametrize("value", [0., 5-math.sqrt(5), 6., 5+math.sqrt(5)])
def test_literal_fock_ground_state_partial_trace_agrees_for_every_actual_native_band(value, report):
    target = next(row["entropy_nats"] for row in report["angular_modes"]
                  if row["native_angular_eigenvalue"] == value)
    coarse = _fock_partial_trace(value, 4)
    fine = _fock_partial_trace(value, 10)
    assert abs(fine-target) < abs(coarse-target)
    assert abs(fine-target) < 2e-6


@pytest.mark.parametrize("rates", [(F(1), F(0)), (F(2), F(-1)), (F(0), F(1)), (F(-1), F(0)), (F(1), F(-1))])
def test_independent_radial_response_has_no_nonisometric_fixed_volume_escape(rates):
    u, radial_rate = rates
    current = v.run(angular_log_rate=u, radial_length_log_rate=radial_rate)
    epsilon = 1e-4
    values = []
    for sign in (-1, 1):
        z = sign*epsilon
        cov = v.covariance(length=F(2*math.exp(float(radial_rate)*z)),
            angular_scale=F(math.exp(2*float(u)*z)), density=F(math.exp(-2*float(u+radial_rate)*z/3)))
        values.append(_entropy(cov["regional_position"], cov["regional_momentum"]))
    secant = (values[1]-values[0])/(2*epsilon)
    assert abs(secant-current["regional_vacuum_entropy_derivative"]) < 5e-7
    assert current["bare_entropy_plus_positive_area_coefficient_stationary"] == (u+radial_rate == 0)
    if u+radial_rate == 0:
        assert abs(values[1]-values[0]) < 5e-13
        assert current["two_face_area_derivative"] == 0


def test_action_shape_response_is_independent_Hamiltonian_variation(report):
    cov = v.covariance()
    lap = np.asarray(v.native_angular(), float)
    shape_derivative = 2*np.kron(np.eye(3), lap)
    # At the ground state kinetic/potential equipartition cancels the mu'
    # terms. The remaining force is <q.B'.q>/2, independently of geometry.
    force = float(np.trace(cov["position"]@shape_derivative)/2)
    assert abs(force-report["independent_action_shape_response"]) < 5e-13
    epsilon = 1e-5
    energy = []
    for z in (-epsilon, epsilon):
        potential = cov["potential"]+(math.expm1(2*z))*np.kron(np.eye(3), lap)
        energy.append(float(np.sqrt(np.linalg.eigvalsh(potential)).sum()/2))
    assert abs((energy[1]-energy[0])/(2*epsilon)-force) < 1e-8


def test_actual_regional_modular_first_law_distinguishes_equal_port_local_stresses(report):
    cov = v.covariance()
    _, modular_p = v.regional_modular()
    off = modular_p-np.diag(np.diag(modular_p))
    a, b = np.unravel_index(np.argmax(abs(off)), off.shape)
    assert a != b and abs(modular_p[a, b]) > .001
    slopes = []
    epsilon = 1e-6
    for sign in (-1, 1):
        direction = np.zeros(12)
        direction[a], direction[b] = 1/math.sqrt(2), sign/math.sqrt(2)
        covariance_change = np.outer(direction, direction)
        # A classical Gaussian momentum displacement gives this positive
        # covariance increment, so it is an admissible full quantum state.
        changed = cov["regional_momentum"]+epsilon*covariance_change
        actual = (_entropy(cov["regional_position"], changed)-report["regional_vacuum_entropy_nats"])/epsilon
        predicted = float(direction@modular_p@direction)/2
        assert abs(actual-predicted) < 2e-5
        slopes.append(predicted)
        # Both signs have exactly the same local p_i^2 and unchanged q,
        # hence the same pairing with ANY port-local canonical stress kernel.
        expected_diagonal = (np.eye(12)[a]+np.eye(12)[b])*(1/math.sqrt(2))**2
        assert np.array_equal(np.diag(covariance_change), expected_diagonal)
    assert abs(slopes[1]-slopes[0]-modular_p[a, b]) < 2e-15


@pytest.mark.parametrize("kwargs", [{"length": True}, {"length": 2.}, {"density": 0},
    {"angular_scale": -1}, {"angular_log_rate": .5}, {"radial_length_log_rate": False}])
def test_comparison_domain_refuses_nonoriginal_or_wrong_sign_parameters(kwargs):
    with pytest.raises(ValueError):
        v.run(**kwargs)


def test_nonzero_exact_deformations_cannot_underflow_into_false_stationarity():
    tiny = F(1, 10**1000)
    with pytest.raises(ValueError, match="execution precision"):
        v.run(angular_log_rate=tiny)
    with pytest.raises(ValueError, match="shape rate"):
        v.run(angular_log_rate=F(1), radial_length_log_rate=-1+tiny)
