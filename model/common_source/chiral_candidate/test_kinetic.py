"""Independent algebraic/negative controls for the finite chiral symbol."""
import json
import math

import numpy as np
import pytest
import sympy as s

from . import action, kinetic


def numeric(value):
    return np.asarray(value.evalf(17) if hasattr(value, "evalf") else value, complex)


@pytest.mark.parametrize("p", [(0, 0, 0), (.3, 0, 0), (0, -.4, 0),
                               (0, 0, .5), (.23, -.31, .17)])
def test_default_symbol_has_exact_mass_singular_value_dispersion(p):
    k = numeric(action.mass_kernel())
    h = kinetic.bdg_symbol(p)
    expected = np.repeat(np.sqrt(np.linalg.svd(k, compute_uv=False)**2
                                 + kinetic.coordinate_speed()**2*np.dot(p, p)), 2)
    expected = np.sort(np.r_[-expected, expected])
    np.testing.assert_allclose(h, h.conj().T, atol=1e-14)
    np.testing.assert_allclose(np.linalg.eigvalsh(h), expected, atol=2e-14)


def test_pairing_matches_action_and_handles_complex_symmetric_kernel():
    np.testing.assert_allclose(kinetic.hamiltonian_pairing(), numeric(action.pairing_kernel()), atol=0)
    k = np.array([[.3+.2j, .11-.07j], [.11-.07j, -.5j]])
    delta = kinetic.hamiltonian_pairing(k)
    np.testing.assert_allclose(delta, -delta.T, atol=0)
    # H_mass=-L_mass and CAR reordering give Delta=-(epsilon K)^dagger.
    a = np.kron(kinetic.spin_epsilon(), k)
    np.testing.assert_allclose(delta, -a.conj().T, atol=0)
    p = (.2, -.3, .7)
    h = kinetic.bdg_symbol(p, k)
    n = 4
    z = np.zeros((n, n))
    cp2 = kinetic.coordinate_speed()**2*np.dot(p, p)
    expected = np.block([[cp2*np.eye(n)+np.kron(np.eye(2), k.conj()@k), z],
                         [z, cp2*np.eye(n)+np.kron(np.eye(2), k@k.conj())]])
    np.testing.assert_allclose(h@h, expected, atol=1e-14)


def test_hole_spin_basis_exposes_left_right_mass_mixing():
    p = (.2, -.3, .7)
    k = np.array([[.3+.2j, .11-.07j], [.11-.07j, -.5j]])
    left = kinetic.weyl_kinetic(p, 2)
    expected = np.block([[left, np.kron(np.eye(2), k.conj())],
                         [np.kron(np.eye(2), k), -left]])
    np.testing.assert_allclose(kinetic.spin_conjugated_symbol(p, k), expected, atol=1e-14)


def test_nambu_conjugation_has_opposite_momentum_and_is_redundancy():
    k = np.array([[.1j, .4], [.4, -.2j]])
    p = (.3, -.2, .7)
    h = kinetic.bdg_symbol(p, k)
    d = len(h)//2
    exchange = np.block([[np.zeros((d, d)), np.eye(d)],
                         [np.eye(d), np.zeros((d, d))]])
    np.testing.assert_allclose(exchange@h.conj()@exchange,
                               -kinetic.bdg_symbol(tuple(-x for x in p), k), atol=1e-14)
    # Using h(p) instead of h(-p) in the lower Nambu block destroys dispersion.
    hp, delta = kinetic.weyl_kinetic(p, 2), kinetic.hamiltonian_pairing(k)
    wrong = np.block([[hp, delta], [delta.conj().T, -hp.T]])
    assert np.linalg.norm(wrong@wrong-h@h) > .1


def test_opposite_momentum_actual_car_source_splits_into_two_nambu_blocks():
    k, p = np.array([[0, .4], [.4, 0]]), (.2, -.3, .7)
    hp, delta = kinetic.momentum_pair_source(p, k)
    np.testing.assert_allclose(hp, hp.conj().T, atol=0)
    np.testing.assert_allclose(delta, -delta.T, atol=0)
    full = np.block([[hp, delta], [delta.conj().T, -hp.T]])
    d = len(hp)//2
    # Original order (p,-p,p†,-p†); compare (p,-p†) and (-p,p†).
    order = np.r_[np.arange(d), np.arange(3*d, 4*d),
                  np.arange(d, 2*d), np.arange(2*d, 3*d)]
    permuted = full[np.ix_(order, order)]
    np.testing.assert_allclose(permuted[:2*d, :2*d], kinetic.bdg_symbol(p, k), atol=0)
    np.testing.assert_allclose(permuted[2*d:, 2*d:],
                               kinetic.bdg_symbol(tuple(-x for x in p), k), atol=0)
    np.testing.assert_allclose(permuted[:2*d, 2*d:], 0, atol=0)
    with pytest.raises(ValueError, match="self-conjugate"):
        kinetic.momentum_pair_source((0, 0, 0), k)


def test_spin_half_and_helicity_are_correct_in_raw_hole_convention():
    spin = kinetic.nambu_spin(1)
    np.testing.assert_allclose(sum(s@s for s in spin), .75*np.eye(4), atol=0)
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        np.testing.assert_allclose(spin[i]@spin[j]-spin[j]@spin[i], 1j*spin[k], atol=0)
    p = (.2, -.3, .7)
    h = kinetic.bdg_symbol(p, np.zeros((1, 1)))
    helicity = kinetic.helicity(p, 1)
    positive = (np.eye(4)+h/(kinetic.coordinate_speed()*math.hypot(*p)))/2
    particle = positive @ np.diag([1, 1, 0, 0])
    antiparticle = positive @ np.diag([0, 0, 1, 1])
    np.testing.assert_allclose(np.trace(particle), 1, atol=1e-14)
    np.testing.assert_allclose(np.trace(helicity@particle), -.5, atol=1e-14)
    np.testing.assert_allclose(np.trace(helicity@antiparticle), .5, atol=1e-14)
    massive = kinetic.bdg_symbol(p, np.array([[.4]]))
    np.testing.assert_allclose(helicity@massive, massive@helicity, atol=1e-14)
    chirality = np.diag([-1, -1, 1, 1])
    assert np.linalg.norm(chirality@massive-massive@chirality) > .5
    with pytest.raises(ValueError, match="undefined"):
        kinetic.helicity((0, 0, 0))


def test_chiral_internal_currents_kinetic_and_unbroken_mass_charge():
    p = (.2, -.3, .7)
    h = kinetic.weyl_kinetic(p)
    for t in action.generators():
        generator = np.kron(np.eye(2), numeric(t))
        np.testing.assert_allclose(generator@h, h@generator, atol=1e-14)
    full = kinetic.bdg_symbol(p)
    charge = kinetic.nambu_charge()
    np.testing.assert_allclose(charge@full, full@charge, atol=1e-14)
    # A fixed neutral Higgs preserves Q but breaks hypercharge separately.
    hypercharge = kinetic.nambu_charge(action.identification_operators()["Y"])
    assert np.linalg.norm(hypercharge@full-full@hypercharge) > .1
    k = numeric(action.mass_kernel())
    derivatives = [numeric(t).T@k+k@numeric(t) for t in action.generators()]
    constraints = np.column_stack([np.r_[d.real.ravel(), d.imag.ravel()] for d in derivatives])
    assert np.linalg.matrix_rank(constraints, tol=1e-10) == 3


def test_no_right_neutrino_added_and_nambu_zero_multiplicity_not_species():
    k = numeric(action.mass_kernel())
    zero_rows = np.flatnonzero(np.max(abs(k), axis=1) == 0)
    np.testing.assert_array_equal(zero_rows, [30, 31, 32])
    assert np.linalg.matrix_rank(k) == 42
    zero_values = np.linalg.eigvalsh(kinetic.bdg_symbol((0, 0, 0)))
    assert np.count_nonzero(abs(zero_values)<1e-12) == 12
    assert np.count_nonzero(zero_values>1e-12) == 84
    moving_values = np.linalg.eigvalsh(kinetic.bdg_symbol((.2, -.3, .7)))
    assert np.count_nonzero(moving_values>1e-12) == 90


def test_naive_lattice_mirrors_are_explicit_negative_control():
    nodes = kinetic.naive_weyl_nodes()
    assert len(nodes) == 8
    assert sum(node["jacobian_chirality"] for node in nodes) == 0
    assert sum(node["jacobian_chirality"]<0 for node in nodes) == 4
    for node in nodes:
        k = tuple(math.pi*x for x in node["pi_coordinates"])
        np.testing.assert_allclose(kinetic.naive_lattice_symbol(k, 1), 0, atol=2e-16)
        # chirality=sign det[d(-c sin k)/dk], not a relabeling of the node.
        jacobian = -kinetic.coordinate_speed()*np.diag(np.cos(k))
        assert np.sign(np.linalg.det(jacobian)) == node["jacobian_chirality"]
    # Deleting exact pi from an odd torus does not suppress near-mirror modes.
    norms = [np.linalg.norm(kinetic.naive_lattice_symbol((math.pi*(1-1/n), 0, 0), 1))
             for n in (9, 19, 39)]
    assert norms[2] < norms[1] < norms[0]


def test_actual_source_weak_raising_current_has_doublet_and_singlet_action():
    coefficients, current = kinetic._weak_raising_exact()
    matrices = action.matrices()
    fundamental = sum((c*g for c, g in zip(coefficients, matrices["generators5"], strict=True)),
                      s.zeros(5)).applyfunc(s.expand)
    expected = s.zeros(5)
    expected[3, 4] = 1
    assert fundamental == expected
    j = numeric(current)
    operators = action.identification_operators()
    t3 = numeric(operators["T3"])
    np.testing.assert_allclose(t3@j-j@t3, j, atol=0)
    # Each native family charged L is raised to its neutral L partner, whereas
    # the all-left charge-conjugated singlet carries no weak raising vertex.
    for f in range(3):
        assert j[30+f, 33+f] == 1
        np.testing.assert_allclose(j[:, 27+f], 0, atol=0)


def test_weak_helicity_discriminator_sums_complete_output_space():
    packet = kinetic.weak_chirality_report((.13, -.21, .34))
    assert packet["complete_final_space_used_in_JdaggerJ"]
    assert packet["unprojected_spin_family_strength_trace"] == {
        "charged_left_doublet": 6., "charge_conjugated_singlet_holes": 0.}
    for row in packet["zero_Higgs_control"]:
        np.testing.assert_allclose(row["total_weak_current_weight"]["negative_helicity"], 1, atol=1e-14)
        np.testing.assert_allclose(row["total_weak_current_weight"]["positive_helicity"], 0, atol=1e-14)
    faster = kinetic.weak_chirality_report((1.3, -2.1, 3.4))
    for row, fast in zip(packet["massive_native_trial"], faster["massive_native_trial"], strict=True):
        actual, expected = row["total_weak_current_weight"], row["predicted_chiral_mixing_weight"]
        for label in ("negative_helicity", "positive_helicity"):
            np.testing.assert_allclose(actual[label], expected[label], atol=1e-14)
        np.testing.assert_allclose(sum(actual.values()), 1, atol=1e-14)
        assert 0 < actual["positive_helicity"] < actual["negative_helicity"] < 1
        assert fast["total_weak_current_weight"]["positive_helicity"] < actual["positive_helicity"]
    with pytest.raises(ValueError, match="nonzero momentum"):
        kinetic.weak_chirality_report((0, 0, 0))


@pytest.mark.parametrize("momentum", [(0, 0), (0, 0, 0, 0), (True, 0, 0),
                                       (float("nan"), 0, 0), (1j, 0, 0)])
def test_bad_momentum_rejected(momentum):
    with pytest.raises(ValueError):
        kinetic.weyl_kinetic(momentum)


def test_report_scopes_physical_and_regulator_identification():
    packet = kinetic.report()
    json.dumps(packet, allow_nan=False)
    assert packet["dimensions"]["charged_Dirac_fields_counting_color"] == 21
    assert packet["dimensions"]["independent_right_neutrinos"] == 0
    assert packet["fixed_neutral_Higgs_current_span"]["unbroken_dimension_numeric"] == 9
    assert sum(packet["positive_energy_electric_charge_multiplicities"].values()) == 90
    assert not packet["local_gauge_covariant_chiral_regulator_established"]
    assert not packet["measured_Standard_Model_particles_identified"]
    assert not packet["native_trial_masses_match_observed_hierarchy"]
    assert max(packet["identities"].values()) < 1e-12
