"""Native collective transport, local collapse records and statistics controls."""
from fractions import Fraction as F

import pytest
import sympy as s
import numpy as np
from scipy.linalg import expm

from . import collective_response as collective


def test_collective_charge_and_matrix_unit_statistics_are_exact():
    modes = collective.fields()
    core = (0, 1)
    charge = s.diag(*[sum(state[p] for p in core) for state in collective.STATES])
    operators = [s.SparseMatrix(66, 66, {(x, y): 1}) for _, x, y in modes]
    for i, operator in enumerate(operators):
        assert charge*operator-operator*charge == operator
        assert operator**2 == s.zeros(66)
        projector = operator*operator.H+operator.H*operator
        assert projector**2 == projector
        assert s.trace(projector) == 2
        assert projector != s.eye(66)
        for other in operators[i+1:]:
            assert operator*other.H == other.H*operator == s.zeros(66)
    assert collective.algebra()["conserved_core_number_charge"] == 1
    assert not collective.algebra()["ten_independent_CAR_modes_realized"]


def test_source_charge_density_and_current_are_actual_native_correlations():
    certificate = collective.charge_continuity()
    assert certificate["source_branch_input_addresses_checked"] == 66*20*2
    modes = collective.fields()
    for j, x, y in modes:
        field = s.SparseMatrix(66, 66, {(x, y): 1})
        for location in range(1, 12):
            density = s.diag(*[state[0]*state[location] for state in collective.STATES])
            assert density*field-field*density == (field if location == j else s.zeros(66))


@pytest.mark.parametrize("law,loss", [("C", s.Rational(3, 64)), ("E", s.Rational(3, 32))])
def test_actual_channel_is_the_induced_common_occupation_graph(law, loss):
    data = collective.transfer(law)
    assert data["survival"] == s.Rational(27, 32)
    assert data["killing"] == loss
    assert data["accepted"] == (data["survival"]-loss)*s.eye(10)-data["laplacian"]/2
    assert data["accepted"]*s.ones(10, 1) == (data["survival"]-loss)*s.ones(10, 1)
    common = [j for j, _, _ in data["fields"]]
    # The common protected port is an isolated response; the nine outer
    # coordinates contain actual transport, including 2 -> 6.
    assert data["accepted"][common.index(6), common.index(2)] == s.Rational(1, 64)
    assert all(data["accepted"][common.index(1), k] == 0 for k in range(10) if k != common.index(1))
    ratios = dict(zip(common, data["modular_weight_ratios"], strict=True))
    assert ratios[1] == F(8, 5)
    assert {ratios[j] for j in common if j != 1} == {F(25, 8)}
    assert all(value == 1 for value in collective.transfer(law, (0, 1))["modular_weight_ratios"])
    assert collective.transfer(law, (0, 1))["killing"] == 0


def test_native_record_choice_changes_damping_without_changing_transport_or_frequency():
    c, e = collective.transfer("C"), collective.transfer("E")
    assert c["laplacian"] == e["laplacian"]
    assert c["modular_weight_ratios"] == e["modular_weight_ratios"]
    assert c["accepted"]-e["accepted"] == s.eye(10)*s.Rational(3, 64)


@pytest.mark.parametrize("law", ["C", "E"])
def test_clocked_propagation_requires_the_two_local_repairs(law):
    common = [j for j, _, _ in collective.fields()]
    before, after = common.index(2), common.index(7)
    one = collective.incident_kernel(law, deadline=1)
    two = collective.incident_kernel(law, deadline=2)
    assert one["success"][after, before] == 0
    assert two["success"][after, before] == s.Rational(1, 2048)
    # First incident and remote transfers do not define one scalar time
    # rescaling; their noncommuting order belongs to the measured clock law.
    data = collective.transfer(law)
    assert data["incident"]*data["remote"] != data["remote"]*data["incident"]


@pytest.mark.parametrize("law,paths", [("C", 4881), ("E", 10921)])
def test_complete_same_history_records_match_the_clock_kernel_and_phase_control(law, paths):
    plus = collective.clock_experiment(law)
    minus = collective.clock_experiment(law, phase_sign=-1)
    assert plus["enumerated_terminal_paths"] == minus["enumerated_terminal_paths"] == paths
    assert plus["zero_probability_paths_retained"] > 0
    assert plus["terminal_probabilities_exact"] == minus["terminal_probabilities_exact"] == {
        "deadline_timeout": "121/256", "selector_timeout": "135/512", "success": "135/512"}
    assert plus["transported_site_7_signed_success_response_exact"] == "1/2048"
    assert minus["transported_site_7_signed_success_response_exact"] == "-1/2048"
    for flag, row in plus["shared_local_detector_joint_law_exact"].items():
        opposite = minus["shared_local_detector_joint_law_exact"][flag]
        assert row["+"] == opposite["-"] and row["-"] == opposite["+"]
        assert row["rest"] == opposite["rest"]
        assert sum(map(F, row.values())) == F(plus["terminal_probabilities_exact"][flag])
        # The equal mixture of opposite phases is the diagonal input control.
        assert (F(row["+"])+F(opposite["+"]))/2 == (F(row["-"])+F(opposite["-"]))/2


def test_overlapping_local_detectors_have_the_same_complete_postcollapse_law():
    certificate = collective.detector_consensus()
    assert certificate["shared_region"] == [0, 5, 7]
    assert certificate["local_embeddings_identical_on_actual_source_exact"]
    first, second = map(tuple, certificate["observer_regions"])
    a = [collective.embed_local(effect, first) for effect in collective.local_effects(first)]
    b = [collective.embed_local(effect, second) for effect in collective.local_effects(second)]
    # An independent nonuniform state has support both in the Ramsey pair
    # and outside it. This checks normalization and all complementary bins.
    vector = s.zeros(66, 1)
    vector[collective.address(0, 7)] = s.Rational(2, 3)
    vector[collective.address(5, 7)] = s.Rational(1, 3)
    vector[collective.address(2, 3)] = s.Rational(2, 3)
    rho = vector*vector.H
    assert s.trace(rho) == 1
    probabilities = []
    for i, effect in enumerate(a):
        probability = s.trace(effect*rho)
        probabilities.append(probability)
        assert probability > 0
        selected = effect*rho*effect/probability
        assert s.trace(selected) == 1
        for j, other in enumerate(b):
            assert s.trace(other*effect*rho*effect*other) == (probability if i == j else 0)
            assert s.trace(other*selected) == (1 if i == j else 0)
    assert sum(probabilities) == 1


def test_detector_rejects_an_observer_region_without_shared_access():
    with pytest.raises(ValueError):
        collective.local_effects((0, 1, 2, 8))


@pytest.mark.parametrize("law", ["C", "E"])
def test_flat_modular_frequency_cannot_give_the_declared_free_field_dispersion(law):
    result = collective.conditional_dispersion_test(law)
    assert result["graph_zero_modes_exact"] == 1
    assert result["graph_positive_modes_exact"] == 8
    assert s.Rational(result["graph_laplacian_trace_exact"]) > 0
    assert not result["tested_sharp_free_dispersion_attachment_passes"]
    assert not result["interacting_thermal_or_alternative_continuum_fields_excluded"]


@pytest.mark.parametrize("law", ["C", "E"])
def test_exact_log_spread_bounds_cover_both_actual_clock_maps(law):
    data = collective.transfer(law)
    bounds = collective.hierarchy_bounds(law)["conventions"]
    success = s.Rational(135, 512)
    maps = {
        "accepted_proposal": data["accepted"],
        "surviving_proposal": data["accepted"]/data["survival"],
        "finite_first_incident_success": collective.incident_kernel(law, deadline=2)["success"],
        "conditional_finite_first_incident_success": collective.incident_kernel(law, deadline=2)["success"]/success,
    }
    for label, matrix in maps.items():
        row = bounds[label]
        lower, upper = map(s.Rational, (row["transfer_lower_exact"], row["transfer_upper_exact"]))
        bound = s.Rational(row["positive_log_rate_spread_upper_exact"])
        assert bound == (1-lower)/(lower*(1-upper))
        assert 0 < lower <= upper < 1 and bound < 100
        values = np.linalg.eigvals(np.array(matrix, float))
        assert max(abs(values.imag)) < 1e-12
        assert min(values.real) >= float(lower)-1e-12
        assert max(values.real) <= float(upper)+1e-12
        rates = -np.log(values.real)
        assert max(rates)/min(rates) <= float(bound)


@pytest.mark.parametrize("copies", [1, 2, 3])
def test_fluctuation_moments_against_full_tensor_reference(copies):
    # Keep an explicit inactive third state: renormalizing the two active
    # states would give the wrong polarization, finite-copy error and KMS state.
    wa, wb = s.Rational(1, 36), s.Rational(2, 225)
    delta = wa-wb
    rho_one = s.diag(wa, wb, 1-wa-wb)
    x = s.SparseMatrix(3, 3, {(0, 1): 1})
    rho = s.kronecker_product(*([rho_one]*copies))
    total = s.zeros(3**copies)
    for index in range(copies):
        factors = [s.eye(3)]*copies
        factors[index] = x
        total += s.kronecker_product(*factors)
    commutator = (total*total.H-total.H*total)/(copies*delta)
    certificate = collective.fluctuation_certificate(copies)
    assert s.trace(rho*commutator) == 1
    assert s.trace(rho*total.H*total)/(copies*delta) == s.Rational(certificate["mean_occupation_exact"])
    assert s.trace(rho*(commutator-s.eye(3**copies))**2) == s.Rational(certificate["commutator_mean_square_error_exact"])
    assert s.trace(rho*total.H**2*total**2)/(copies*delta)**2 == s.Rational(certificate["normal_ordered_fourth_moment_exact"])


def test_ordered_weyl_formula_against_tensor_exponentials_and_order_control():
    wa, wb = 1/36, 2/225
    delta = wa-wb
    x = np.array([[0., 1., 0.], [0., 0., 0.], [0., 0., 0.]], complex)
    total = np.kron(x, np.eye(3))+np.kron(np.eye(3), x)
    q, p = (total+total.conj().T)/np.sqrt(4*delta), (total-total.conj().T)/(1j*np.sqrt(4*delta))
    rho = np.kron(np.diag([wa, wb, 1-wa-wb]), np.diag([wa, wb, 1-wa-wb]))
    vectors = [(0.2, 0.), (0., 0.3)]
    direct = np.trace(rho@expm(0.2j*q)@expm(0.3j*p))
    assert abs(direct-collective.ordered_weyl(vectors, 2)[0]) < 1e-13
    errors = []
    for copies in (4, 64, 1024, 16384):
        observed, limit = collective.ordered_weyl(vectors, copies)
        reverse, reverse_limit = collective.ordered_weyl(vectors[::-1], copies)
        assert abs(reverse-observed.conjugate()) < 1e-13
        assert abs(reverse_limit-limit.conjugate()) < 1e-13
        assert abs(observed.imag) > 0.02
        errors.append(abs(observed-limit))
    assert all(a > b for a, b in zip(errors, errors[1:]))
    assert errors[-1] < 3e-6


def test_local_return_response_exposes_all_transfer_poles():
    for law in ("C", "E"):
        transfer = collective.transfer(law)["accepted"]
        eigenvalues = np.linalg.eigvalsh(np.array(transfer, float))
        for proposals in range(1, 5):
            # Averaging the ten actual local Ramsey return contrasts gives
            # positive modal weight, so a dark coordinate cannot hide a pole.
            actual = s.trace(transfer**proposals)/10
            assert abs(float(actual)-sum(eigenvalues**proposals)/10) < 1e-13
