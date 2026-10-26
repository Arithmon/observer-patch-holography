import json
import math

import numpy as np
import pytest

from . import clifford as c
from ..mass_response.dispersion import bloch_hamiltonian, refinement_coefficients


def test_continuum_lorentz_and_spin_algebras_are_constructed():
    matrices = c.continuum_matrices()
    gamma = matrices["lorentz_gamma"]
    spin = matrices["spin"]
    for i, a in enumerate(gamma):
        for j, b in enumerate(gamma):
            expected = 2 * (1 if i == 0 else -1) * np.eye(8) if i == j else np.zeros((8, 8))
            np.testing.assert_allclose(a @ b + b @ a, expected, atol=2e-15)
    np.testing.assert_allclose(sum(s @ s for s in spin), .75 * np.eye(8), atol=2e-15)
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        np.testing.assert_allclose(spin[i] @ spin[j] - spin[j] @ spin[i], 1j * spin[k], atol=2e-15)


def test_full_operator_basis_proves_two_dirac_copies_and_complete_commutant():
    matrices = c.continuum_matrices()
    generators = (*matrices["alpha"], matrices["beta"])
    flavor = c.commutant_basis()
    words = c.clifford_words(generators)
    basis = np.array([(a @ t).ravel() for a in words for t in flavor])
    np.testing.assert_allclose(basis.conj() @ basis.T / 8, np.eye(64), atol=3e-15)
    assert c.commutant_nullity(generators) == 4
    f0, f1 = c.flavor_projectors()
    assert np.trace(f0) == pytest.approx(4.)
    assert np.trace(f1) == pytest.approx(4.)
    np.testing.assert_allclose(f0 @ f1, np.zeros((8, 8)), atol=2e-15)
    t1 = matrices["flavor"][0]
    np.testing.assert_allclose(t1 @ f0 @ t1, f1, atol=2e-15)


@pytest.mark.parametrize("p", [(0., 0., 0.), (.2, -.3, .4), (1., .5, -.7)])
def test_each_energy_sign_has_two_spins_for_each_of_two_flavors(p):
    h = c.continuum_hamiltonian(p)
    positive, negative = c.energy_projectors(p)
    np.testing.assert_allclose(positive @ negative, np.zeros((8, 8)), atol=2e-15)
    np.testing.assert_allclose(positive + negative, np.eye(8), atol=2e-15)
    scales = c.source_scales()
    energy = math.sqrt(scales["candidate_rest_energy"] ** 2
                       + scales["candidate_speed"] ** 2 * sum(x * x for x in p))
    np.testing.assert_allclose(h @ positive, energy * positive, atol=2e-15)
    np.testing.assert_allclose(h @ negative, -energy * negative, atol=2e-15)
    for flavor in c.flavor_projectors():
        for projector in (positive, negative):
            joint = flavor @ projector
            assert np.trace(joint) == pytest.approx(2.)
            np.testing.assert_allclose(joint @ joint, joint, atol=2e-15)


def test_charge_conjugation_reverses_energy_and_momentum_but_is_not_em_identification():
    data = c.continuum_matrices()
    b = data["charge_conjugation_matrix"]
    p = (.31, -.27, .43)
    np.testing.assert_allclose(b @ b.conj(), np.eye(8), atol=2e-15)
    np.testing.assert_allclose(b @ c.continuum_hamiltonian(p).conj() @ b.conj().T,
                               -c.continuum_hamiltonian(tuple(-x for x in p)), atol=2e-15)


def test_continuum_rotation_covariance_and_massive_chirality_mixing():
    data = c.continuum_matrices()
    angle = .37
    sz = data["spin"][2]
    u = math.cos(angle / 2) * np.eye(8) - 2j * math.sin(angle / 2) * sz
    p = (.3, -.4, .7)
    rotated = (math.cos(angle) * p[0] - math.sin(angle) * p[1],
               math.sin(angle) * p[0] + math.cos(angle) * p[1], p[2])
    np.testing.assert_allclose(u @ c.continuum_hamiltonian(p) @ u.conj().T,
                               c.continuum_hamiltonian(rotated), atol=2e-15)
    chirality = data["chirality"]
    commutator = c.continuum_hamiltonian(p) @ chirality - chirality @ c.continuum_hamiltonian(p)
    assert np.linalg.norm(commutator, 2) == pytest.approx(2 * c.source_scales()["candidate_rest_energy"])


def test_finite_lattice_does_not_inherit_continuum_flavor_commutant():
    assert c.commutant_nullity(c.lattice_probe_matrices()) == 1
    h = bloch_hamiltonian((.2, -.5, .9))
    assert max(np.linalg.norm(h @ t - t @ h) for t in c.commutant_basis()[1:]) > .01


def test_regulator_hamiltonian_approaches_the_classified_continuum_operator():
    p = (.3, -.4, .7)
    errors = []
    for level in (4, 8, 16, 32):
        a, b = refinement_coefficients(level)
        h = bloch_hamiltonian(tuple(x / level for x in p), intra=a, inter=b)
        errors.append(np.linalg.norm(h - c.continuum_hamiltonian(p), 2))
    assert errors[-1] < errors[0] / 7


def test_report_is_json_ready_and_keeps_physical_species_open():
    record = c.report()
    json.dumps(record, allow_nan=False)
    assert record["kinematic_classification"]["continuum_dirac_multiplicity"] == 2
    assert record["kinematic_classification"]["continuum_commutant_nullity_numeric"] == 4
    assert record["regulator_boundary"]["lattice_commutant_nullity_numeric"] == 1
    assert record["standard_model_boundary"]["observed_species_identified"] == []
    assert not record["standard_model_boundary"]["photon_or_higgs_identified"]
    assert not record["regulator_boundary"]["massless_single_particle_pole_in_this_free_gapped_continuation"]
    assert record["particle_antiparticle_contract"]["zero_temperature_neutral_bilinear_threshold_lower_bound"] > 0
    assert not record["physical_particle_identification_complete"]
    assert max(record["numerical_identity_residuals"].values()) < 5e-14


@pytest.mark.parametrize("p", [(1., 2.), (1., 2., float("nan")), (1., 2., True), (1., 2., 1j)])
def test_invalid_momenta_are_rejected(p):
    with pytest.raises(ValueError):
        c.continuum_hamiltonian(p)
