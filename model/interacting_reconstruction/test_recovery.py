"""Independent finite systems test the record-only interacting inverse."""
from copy import deepcopy
from itertools import combinations, product
from pathlib import Path
import ast

import numpy as np
import pytest
from scipy.linalg import expm

from model.interacting_reconstruction import recovery as r


def operators(dimension, omega):
    # Assemble the infinite-polynomial matrix entries in ample independent
    # workspace, rather than borrowing any source executor or action helper.
    a = np.diag(np.sqrt(np.arange(1, dimension+8)), 1)
    q = (a+a.T)/np.sqrt(2*omega)
    p = 1j*np.sqrt(omega/2)*(a.T-a)
    return [x[:dimension, :dimension] for x in (q, q@q, (p@p).real, q@q@q@q)]


def law(generator, a, b, time):
    state = expm(-1j*time*generator)[:, a]
    records = []
    for sign in (-1, 1):
        ket = np.zeros(len(generator), complex)
        ket[a], ket[b] = 1/np.sqrt(2), sign*1j/np.sqrt(2)
        yes = abs(np.vdot(ket, state))**2
        records.append([1-yes, yes])
    return np.array(records)


def clock(time):
    return [(1+np.cos(time))/2, (1-np.sin(time))/2]


def cube():
    positions = np.array(list(product((0, 1), repeat=3)))
    adjacency = (np.sum(abs(positions[:, None]-positions[None, :]), axis=-1) == 1).astype(int)
    permutation = np.array([7, 2, 0, 5, 3, 1, 6, 4])
    return adjacency[np.ix_(permutation, permutation)], positions[permutation]-.5


def records(time, *, mass_squared=1., quartic=1/24, coupling=1., omega=2., adjacency=None):
    if adjacency is None:
        adjacency, _ = cube()
    q, q2, p2, q4 = operators(4, omega)
    mass_squared = np.broadcast_to(mass_squared, (8,))
    quartic = np.broadcast_to(quartic, (8,))
    onsite = []
    for site in range(8):
        a = mass_squared[site]+coupling*adjacency[site].sum()
        h = .5*p2+.5*a*q2+quartic[site]*q4
        onsite.append([law(h, 0, 2, time), law(h, 1, 3, time)])
    q = operators(2, omega)[0]
    pairs = list(combinations(range(8), 2))
    bonds = [law(-coupling*adjacency[i, j]*np.kron(q, q), 0, 3, time) for i, j in pairs]
    return {"onsite_probabilities": np.array(onsite), "bond_probabilities": np.array(bonds),
            "pairs": pairs}


def recover(first=None, second=None, *, parameters=None, **kwargs):
    parameters = parameters or {}
    u, v = 1/512, 1/256
    first = records(u, **parameters) if first is None else first
    second = records(v, **parameters) if second is None else second
    return r.recover(first, second, clock(u), clock(v), omega=2,
                     norm_bound_onsite=16, norm_bound_bond=2, **kwargs)


def test_hidden_interacting_cube_recovers_geometry_and_action_within_intervals():
    result = recover()
    adjacency, positions = cube()
    assert np.array_equal(result["adjacency"], adjacency)
    assert result["spatial_dimension"] == 3
    assert np.allclose(result["gram"], positions@positions.T)
    assert np.allclose(result["coordinates"]@result["coordinates"].T, result["gram"])
    for name, value in (("quartic", 1/24), ("edge_weight", 1.), ("bare_mass_squared", 1.)):
        interval = result["homogeneous_intervals"][name]
        assert interval["lower"] < value < interval["upper"]
        assert abs(interval["value"]-value) < 1e-8
        assert interval["lower"] > 0
    assert np.all(abs(result["onsite_quartic"]-1/24) <= result["onsite_quartic_radius"])
    assert np.all(abs(result["site_bare_mass_squared"]-1) <= result["site_bare_mass_squared_radius"])
    assert result["edge_gap_is_supplied"]
    assert result["bare_mass_is_not_measured_interacting_gap"]
    assert result["bare_mass"] == pytest.approx(1., abs=1e-8)
    assert result["bare_speed_per_unit_edge"] == pytest.approx(1., abs=1e-8)
    assert result["bond_kappa"] == pytest.approx(1., abs=1e-8)
    assert np.allclose(result["bare_quadratic_stiffness"], 4*np.eye(8)-adjacency)
    assert not {"mass", "speed_per_unit_edge", "coupling", "stiffness", "frequencies"}.intersection(result)


@pytest.mark.parametrize("omega,mass_squared,quartic", [(1.3, .7, .08), (2., 1., .04), (3.1, 4., .11)])
def test_two_transitions_resolve_infinite_polynomials_compressed_afterwards(omega, mass_squared, quartic):
    _, q2, p2, q4 = operators(4, omega)
    h = .5*p2+.5*mass_squared*q2+quartic*q4
    result = r.onsite_parameters(h[2, 0], h[3, 1], 0., 0., omega=omega)
    assert result["quartic"] == pytest.approx(quartic, abs=1e-14)
    assert result["onsite_quadratic"] == pytest.approx(mass_squared, abs=1e-14)
    # Squaring compressed q would silently change H31 at the boundary.
    q = operators(4, omega)[0]
    assert not np.allclose(q4, q@q@q@q)


def test_two_level_onsite_space_cannot_resolve_mass_from_quartic():
    omega, change = 2., .17
    _, q2, _, q4 = operators(2, omega)
    difference = change*q4-3*change/omega*q2
    assert np.allclose(difference, -3*change/(4*omega**2)*np.eye(2))
    # A scalar energy offset drops out of every state probability.  The
    # four-level protocol supplies two independent transitions instead.


def test_richardson_recovers_random_real_generator_at_unequal_times():
    random = np.random.default_rng(157)
    for dimension in (3, 4, 7):
        matrix = random.normal(size=(dimension, dimension))
        matrix = (matrix+matrix.T)/2
        matrix /= np.linalg.norm(matrix, 2)
        first, second = .013, .041
        p = np.array([law(matrix, 0, 2, t) for t in (first, second)])
        estimate = r.coefficient_from_records(p, first, second, norm_bound=1)
        assert abs(estimate["matrix_element"]-matrix[2, 0]) < estimate["radius"]
        raw_derivative = (p[1, 0, 1]-p[1, 1, 1])/(2*second)
        assert abs(estimate["matrix_element"]-matrix[2, 0]) < abs(raw_derivative-matrix[2, 0])


def test_sign_order_is_fixed_by_the_projector_protocol():
    matrix = np.array([[.4, .8], [.8, -.3]])
    p = np.array([law(matrix, 0, 1, t) for t in (.01, .02)])
    result = r.coefficient_from_records(p, .01, .02, norm_bound=2)
    reversed_result = r.coefficient_from_records(p[:, ::-1], .01, .02, norm_bound=2)
    assert result["matrix_element"] == pytest.approx(.8, abs=result["radius"])
    assert reversed_result["matrix_element"] == pytest.approx(-.8, abs=result["radius"])


def test_scalar_energy_offsets_do_not_change_the_inverse():
    matrix = np.array([[.4, .8], [.8, -.3]])
    for offset in (0, 500):
        p = np.array([law(matrix+offset*np.eye(2), 0, 1, t) for t in (.01, .02)])
        estimate = r.coefficient_from_records(p, .01, .02, norm_bound=2)
        assert abs(estimate["matrix_element"]-.8) < estimate["radius"]


def test_probability_error_budget_covers_adversarial_law_perturbations():
    matrix = np.array([[.4, .8], [.8, -.3]])
    u, v, error = .01, .02, 1e-5
    p = np.array([law(matrix, 0, 1, t) for t in (u, v)])
    for i, j in product(range(2), repeat=2):
        perturbation = error*(-1)**(i+j)
        p[i, j] += [-perturbation, perturbation]
    result = r.coefficient_from_records(p, u, v, norm_bound=2, probability_error=error)
    assert abs(result["matrix_element"]-.8) < result["radius"]
    assert result["probability_radius"] == pytest.approx(3*error/v)


def test_dyadic_finite_time_radius_matches_closed_form():
    p = np.full((2, 2, 2), .5)
    result = r.coefficient_from_records(p, 1/512, 1/256, norm_bound=16)
    assert result["finite_time_radius"] == pytest.approx(16**5*(1/256)**4/18)
    assert result["finite_time_radius_exact"] == "1/73728"


def test_data_without_quartic_cannot_be_reported_as_positive_interaction():
    with pytest.raises(ValueError, match="quartic interaction.*strictly positive"):
        recover(parameters={"quartic": 0.})


def test_matched_free_control_can_report_an_interval_containing_zero():
    result = recover(parameters={"quartic": 0.}, require_positive_interaction=False)
    assert not result["positive_quartic_interaction_resolved"]
    assert result["homogeneous_intervals"]["quartic"]["lower"] == 0
    assert result["homogeneous_intervals"]["quartic"]["upper"] > 0
    assert result["quartic"] == pytest.approx(0, abs=1e-12)


def test_optional_free_control_does_not_allow_an_attractive_action():
    with pytest.raises(ValueError, match="repulsive action"):
        recover(parameters={"quartic": -.03}, require_positive_interaction=False)


def test_negative_quartic_cannot_be_reported_as_stable_interaction():
    with pytest.raises(ValueError, match="quartic interaction.*strictly positive"):
        recover(parameters={"quartic": -.03})


def test_inhomogeneous_mass_does_not_get_averaged_into_a_homogeneous_fit():
    with pytest.raises(ValueError, match="bare mass squared.*homogeneous"):
        recover(parameters={"mass_squared": np.linspace(.8, 1.2, 8)})


def test_inhomogeneous_quartic_does_not_get_averaged_into_a_homogeneous_fit():
    with pytest.raises(ValueError, match="quartic interaction.*homogeneous"):
        recover(parameters={"quartic": np.linspace(.02, .07, 8)})


def test_weak_nonzero_edges_violate_the_supplied_separation_family():
    with pytest.raises(ValueError, match="edge gap"):
        recover(parameters={"coupling": .05})


def test_large_readout_error_leaves_topology_unresolved():
    with pytest.raises(ValueError, match="edge gap"):
        recover(probability_error=.0001)


def test_other_graph_is_not_silently_assigned_cube_geometry():
    ring = np.zeros((8, 8), dtype=int)
    for i in range(8):
        ring[i, (i+1) % 8] = ring[(i+1) % 8, i] = 1
    with pytest.raises(ValueError, match="cube|couplings"):
        recover(parameters={"adjacency": ring})


@pytest.mark.parametrize("change", ["pair_missing", "pair_repeated", "bad_shape", "nan", "complex", "bad_normalization", "missing_outcomes"])
def test_malformed_calibration_is_rejected(change):
    first = records(1/512)
    if change == "pair_missing":
        first["pairs"] = first["pairs"][:-1]
    elif change == "pair_repeated":
        first["pairs"][0] = first["pairs"][1]
    elif change == "bad_shape":
        first["onsite_probabilities"] = first["onsite_probabilities"][:7]
    elif change == "nan":
        first["bond_probabilities"][0, 0, 0] = np.nan
    elif change == "complex":
        first["onsite_probabilities"] = first["onsite_probabilities"].astype(complex)
    elif change == "bad_normalization":
        first["bond_probabilities"][0, 0] = [.4, .4]
    elif change == "missing_outcomes":
        first["bond_probabilities"] = first["bond_probabilities"][..., 1]
    with pytest.raises(ValueError, match="records|pairs"):
        recover(first)


@pytest.mark.parametrize("name,value", [("norm_bound", 0), ("norm_bound", True), ("norm_bound", np.nan),
                                        ("probability_error", -.1), ("probability_error", .5)])
def test_invalid_coefficient_hypotheses_are_rejected(name, value):
    options = {"norm_bound": 2}
    options[name] = value
    with pytest.raises(ValueError):
        r.coefficient_from_records(np.full((2, 2, 2), .5), .01, .02, **options)


@pytest.mark.parametrize("times", [(0, .01), (.02, .01), (.01, .01), (True, 2), (.01, np.inf)])
def test_invalid_clock_order_is_rejected(times):
    with pytest.raises(ValueError):
        r.coefficient_from_records(np.full((2, 2, 2), .5), *times, norm_bound=2)


def test_input_arrays_are_not_mutated():
    first, second = records(1/512), records(1/256)
    backup = deepcopy((first, second))
    recover(first, second)
    for actual, original in zip((first, second), backup):
        for key in actual:
            assert np.array_equal(actual[key], original[key])


def test_inverse_has_no_source_executor_or_hidden_action_import():
    tree = ast.parse(Path(r.__file__).read_text())
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
    assert not any("source" in name or "lattice" in name or "experiment" in name for name in imports)
