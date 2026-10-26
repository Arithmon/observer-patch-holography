"""Independent small-Fock-space checks of the Hubbard diagnostic."""
import json
import math

import numpy as np
import pytest

from . import hubbard, model


def _full_fock_hamiltonian(length, t, u):
    """Independent CAR construction, using tensor-product Pauli strings."""
    modes = 2 * length
    identity = np.eye(2)
    parity = np.diag([1., -1.])
    lower = np.array([[0., 1.], [0., 0.]])
    annihilators = []
    for mode in range(modes):
        operator = np.ones((1, 1))
        for position in reversed(range(modes)):
            local = lower if position == mode else parity if position < mode else identity
            operator = np.kron(operator, local)
        annihilators.append(operator)
    number = [a.T @ a for a in annihilators]
    h = np.zeros((1 << modes, 1 << modes))
    for spin in (0, 1):
        for site in range(length):
            a = annihilators[spin * length + site]
            b = annihilators[spin * length + (site + 1) % length]
            h -= t * (a.T @ b + b.T @ a)
    for site in range(length):
        h += u * number[site] @ number[length + site]
    return h


@pytest.fixture(scope="module")
def controls():
    return hubbard.finite_volume_controls()


def test_sparse_fixed_spin_matches_independent_full_car_and_dense_ground_energy():
    length, t, u = 4, .73, .91
    full = _full_fock_hamiltonian(length, t, u)
    for n_up, n_down in ((2, 2), (3, 2), (1, 2)):
        indices = [up | (down << length)
                   for up in hubbard.spin_basis(length, n_up)
                   for down in hubbard.spin_basis(length, n_down)]
        expected = full[np.ix_(indices, indices)]
        actual = hubbard.sector_hamiltonian(length, n_up, n_down,
                                             hopping=t, repulsion=u)
        np.testing.assert_allclose(actual.toarray(), expected, rtol=0, atol=1e-14)
        assert np.linalg.norm(actual.toarray() - actual.toarray().T) == 0
        diagnostic = hubbard.sector_ground_state(length, n_up, n_down,
                                                   hopping=t, repulsion=u)
        assert diagnostic["energy"] == pytest.approx(np.linalg.eigvalsh(expected)[0], abs=2e-12)
        assert diagnostic["eigenpair_residual_norm"] < 2e-11


def test_report_uses_actual_preparation_and_preserves_particle_hole_identity(controls):
    source = model.parameters()
    assert controls["source"] == source.record()
    assert controls["architecture"]["hopping"] == source.clock / 2
    assert controls["architecture"]["repulsion"] == source.interaction
    assert controls["architecture"]["hopping"] == pytest.approx(math.log(5) / 2)
    assert controls["architecture"]["repulsion"] == pytest.approx(math.log(125 / 64))
    assert [row["length"] for row in controls["rows"]] == [4, 6, 8]
    for row in controls["rows"]:
        assert row["charge_gap"] > 0
        assert abs(row["particle_hole_energy_residual"]) < 2e-10
        assert row["charge_gap"] == pytest.approx(
            row["Eplus"] + row["Eminus"] - 2 * row["E0"], abs=1e-14)
        for sector in row["sectors"].values():
            assert sector["dimension"] <= 4900
            assert sector["eigenpair_residual_norm"] < 2e-10
            assert sector["eigenvector_norm_error"] < 2e-14
            assert not sector["eigenvalue_ordering_certified"]
    assert len(json.dumps(controls, allow_nan=False)) > 100


def test_free_chain_shell_formula_matches_single_particle_filling_and_ed(controls):
    t = controls["architecture"]["hopping"]
    for row in controls["rows"]:
        length = row["length"]
        levels = np.sort([-2 * t * math.cos(2 * math.pi * j / length)
                          for j in range(length) for _spin in range(2)])
        filling_gap = levels[length] - levels[length - 1]
        assert row["free_charge_gap_analytic"] == pytest.approx(filling_gap, abs=2e-15)
        check = row["free_numerical_check"]
        if check is not None:
            assert check["charge_gap"] == pytest.approx(filling_gap, abs=2e-11)
            assert abs(check["particle_hole_energy_residual"]) < 2e-11
    assert controls["rows"][0]["free_charge_gap_analytic"] == 0
    assert controls["rows"][1]["free_charge_gap_analytic"] == pytest.approx(2 * t)
    assert controls["rows"][2]["free_charge_gap_analytic"] == 0


def test_small_volume_shell_effect_prevents_monotone_mass_interpretation(controls):
    gaps = [row["charge_gap"] for row in controls["rows"]]
    assert gaps[1] > gaps[0] > gaps[2]
    interpretation = controls["interpretation"]
    assert interpretation["numerical_diagnostic"]
    assert not interpretation["certified_spectrum"]
    assert not interpretation["thermodynamic_extrapolation_performed"]
    assert not interpretation["finite_addition_gap_is_physical_mass"]
    assert not interpretation["three_spatial_dimensions_established"]


def test_empty_and_full_spin_sectors_have_exact_energies():
    empty = hubbard.sector_ground_state(4, 0, 0, hopping=1., repulsion=.7)
    full = hubbard.sector_ground_state(4, 4, 4, hopping=1., repulsion=.7)
    assert empty["energy"] == 0
    assert full["energy"] == pytest.approx(2.8)
    assert empty["eigenpair_residual_norm"] == full["eigenpair_residual_norm"] == 0


@pytest.mark.parametrize("length", [True, 3, 5, 10, 4.0])
def test_invalid_volume_rejected(length):
    with pytest.raises(ValueError):
        hubbard.sector_hamiltonian(length, 2, 2)


@pytest.mark.parametrize("population", [True, -1, 5, 2.0])
def test_invalid_spin_population_rejected(population):
    with pytest.raises(ValueError):
        hubbard.sector_hamiltonian(4, population, 2)


@pytest.mark.parametrize("t,u", [(0, 1), (-1, 1), (1, -.1), (True, 1), (1, float("nan")),
                               (float("inf"), 1)])
def test_invalid_couplings_rejected(t, u):
    with pytest.raises(ValueError):
        hubbard.sector_hamiltonian(4, 2, 2, hopping=t, repulsion=u)


def test_duplicate_and_empty_volumes_rejected():
    for lengths in ((), (4, 4)):
        with pytest.raises(ValueError):
            hubbard.finite_volume_controls(lengths)
