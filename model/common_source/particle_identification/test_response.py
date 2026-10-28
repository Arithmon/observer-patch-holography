"""Independent correlator, CAR, grading and complete-source readout checks."""
import math

import numpy as np
import pytest
from scipy.linalg import expm

from ..mass_response import dispersion, model
from . import response as r


@pytest.mark.parametrize("interaction", [False, True])
def test_charged_residues_are_positive_complete_and_hit_the_measured_threshold(interaction):
    result = r.finite_response(interaction=interaction)
    integrated = sum(p.residue for p in (*result.addition, *result.removal))
    np.testing.assert_allclose(integrated, np.eye(4), atol=3e-15)
    for channel, charge in ((result.addition, 1), (result.removal, -1), (result.density, 0)):
        for pole in channel:
            assert pole.charge_change == charge
            assert pole.target_number == 2+charge
            np.testing.assert_allclose(pole.residue, pole.residue.conj().T, atol=1e-15)
            assert np.linalg.eigvalsh(pole.residue).min() > -2e-15
    gap = dispersion.ring_gap_predictions(coupling=model.parameters().interaction if interaction else 0.)
    assert result.addition[0].energy_cost == pytest.approx(gap["charged_addition_gap"], abs=6e-15)
    assert result.removal[0].energy_cost == pytest.approx(gap["charged_removal_gap"], abs=6e-15)
    assert np.trace(result.addition[0].residue).real > .1
    assert np.trace(result.removal[0].residue).real > .1


@pytest.mark.parametrize("time", [0., .23, -.41])
def test_all_Lehmann_time_functions_equal_independent_operator_evolution(time):
    result = r.finite_response()
    h, g, energy = result.hamiltonian, result.ground, result.ground_energy
    ladders = model.annihilators()
    forward = expm(-1j*(h-energy*np.eye(16))*time)
    plus = np.array([[np.vdot(g, a@forward@b.conj().T@g) for b in ladders] for a in ladders])
    minus = np.array([[np.vdot(g, b.conj().T@forward.conj().T@a@g) for b in ladders] for a in ladders])
    centered = [a.conj().T@a-mean*np.eye(16) for a, mean in zip(ladders, result.density_means, strict=True)]
    density = np.array([[np.vdot(g, a@forward@b@g) for b in centered] for a in centered])
    np.testing.assert_allclose(r.charged_correlation(result, time), plus, atol=5e-15)
    np.testing.assert_allclose(r.charged_correlation(result, time, kind="removal"), minus, atol=5e-15)
    np.testing.assert_allclose(r.density_correlation(result, time), density, atol=5e-15)


def test_neutral_density_has_distinct_selection_rules_and_no_charge_pole():
    result = r.finite_response()
    assert sum(p.multiplicity for p in result.density) == 6
    assert any(np.trace(p.residue).real < r.VISIBILITY_TOLERANCE for p in result.density)
    assert all(p.charge_change == 0 for p in result.density)
    for pole in result.density:
        assert abs(pole.residue.sum()) < 2e-15  # total N is conserved
    report = r.report()
    assert report["lowest_visible_neutral_density_cost"] != pytest.approx(
        report["lowest_CAR_visible_addition_cost"], abs=1e-8)
    assert max(abs(value) for value in report["sum_rules"].values()) < 1e-13
    assert not report["Standard_Model_species_identified"]
    assert not report["Lorentzian_interacting_particle_pole_established"]


def test_common_chemical_potential_moves_odd_costs_and_cancels_from_transfer():
    base, shifted = r.finite_response(), r.finite_response(chemical_potential=.37)
    for left, right in zip(base.addition, shifted.addition, strict=True):
        assert right.energy_cost-left.energy_cost == pytest.approx(.37, abs=1e-14)
        np.testing.assert_allclose(right.residue, left.residue, atol=3e-15)
    for left, right in zip(base.removal, shifted.removal, strict=True):
        assert right.energy_cost-left.energy_cost == pytest.approx(-.37, abs=1e-14)
        np.testing.assert_allclose(right.residue, left.residue, atol=3e-15)
    assert r.transfer_visibility(base)["pair_frequency"] == pytest.approx(
        r.transfer_visibility(shifted)["pair_frequency"], abs=2e-14)
    assert r.transfer_visibility(base)["particle_hole_gap_identification_admitted"]
    assert not r.transfer_visibility(shifted)["particle_hole_gap_identification_admitted"]
    assert r.transfer_visibility(shifted)["single_gap_over_clock_under_particle_hole_symmetry"] is None


def test_retarded_resolvent_has_positive_spectral_measure_and_CAR_normalization():
    result = r.finite_response()
    z = .24+.17j
    green = r.green_function(result, z)
    spectral = -(green-green.conj().T)/(2j)
    assert np.linalg.eigvalsh(spectral).min() > 0
    expected = np.zeros((4, 4), complex)
    g, h, energy = result.ground, result.hamiltonian, result.ground_energy
    plus_inverse = np.linalg.inv(z*np.eye(16)-(h-energy*np.eye(16)))
    minus_inverse = np.linalg.inv(z*np.eye(16)+(h-energy*np.eye(16)))
    for i, a in enumerate(model.annihilators()):
        for j, b in enumerate(model.annihilators()):
            expected[i, j] = np.vdot(g, a@plus_inverse@b.conj().T@g+b.conj().T@minus_inverse@a@g)
    np.testing.assert_allclose(green, expected, atol=3e-15)


@pytest.mark.parametrize("quadrature", ["X", "Y", "Z"])
@pytest.mark.parametrize("sites", [(0, 0), (3, 1)])
def test_actual_source_transfer_instrument_matches_positive_Kraus_on_entangled_references(quadrature, sites):
    rng = np.random.default_rng(4106)
    # An existing clock spectator, with three untouched reference columns.
    state = rng.normal(size=(512, 3))+1j*rng.normal(size=(512, 3))
    state /= 3*np.linalg.norm(state)  # An unnormalized pre-existing branch.
    observable = r.transfer_operator(*sites, quadrature=quadrature)
    result = r.transfer_read(state, *sites, quadrature=quadrature)
    targets = r.transfer_kraus(*sites, quadrature=quadrature)
    for bit, sign in enumerate((1, -1)):
        kraus = (np.eye(256)/math.sqrt(2)+sign*observable/2
                  +(.5-1/math.sqrt(2))*(observable@observable))
        expected = model.apply_subsystem(state, kraus, range(8), 9)
        np.testing.assert_allclose(targets[bit], kraus, atol=2e-16)
        np.testing.assert_allclose(result.outputs[bit], expected, atol=3e-14, rtol=0)
    assert all(output.shape == state.shape for output in result.outputs)
    assert sum(np.vdot(output, output).real for output in result.outputs) == pytest.approx(1/9, abs=3e-14)
    assert result.resources["work_preparations"] == 0
    assert result.resources["fresh_pointer_blanks"] == 1
    assert result.resources["native_pulses"] == len(result.word)
    assert not result.resources["physical_service_error_certified"]
    charge = np.diag([mask.bit_count() for mask in range(256)])
    assert np.max(abs(charge@observable-observable@charge)) == 0


def test_two_arm_visibility_has_correct_global_fermion_sign_and_phase():
    result = r.finite_response()
    e1, g1 = model.sector_ground(result.hamiltonian, 1)
    e3, g3 = model.sector_ground(result.hamiltonian, 3)
    ground = np.kron(result.ground, result.ground)
    charged = np.kron(g1, g3)
    x, y = (r.transfer_operator(quadrature=axis) for axis in ("X", "Y"))
    visibility = r.transfer_visibility(result)
    assert visibility["transfer_weight"] > .01
    assert visibility["factorization_residual"] < 1e-16
    assert visibility["global_occupation_U1_invariant"]
    assert not visibility["local_electromagnetic_gauge_invariance_established"]
    m = complex(*visibility["transfer_matrix_element"])
    gap = e1+e3-2*result.ground_energy
    psi, specification = model.charge_transfer_preparation()
    np.testing.assert_allclose(psi, (ground+charged)/math.sqrt(2), atol=1e-15)
    assert gap == pytest.approx(specification["pair_gap"], abs=1e-15)
    for time in (0., .25, .61):
        evolved = (ground+np.exp(-1j*gap*time)*charged)/math.sqrt(2)
        assert np.vdot(evolved, x@evolved).real == pytest.approx((m*np.exp(1j*gap*time)).real, abs=2e-15)
        assert np.vdot(evolved, y@evolved).real == pytest.approx((m*np.exp(1j*gap*time)).imag, abs=2e-15)
    # The Jordan-Wigner string affects an occupied intervening mode.
    probe = np.zeros(256, complex)
    probe[(1 << 4)+(1 << 1)] = 1
    transferred = (x+1j*y)@probe/2
    assert transferred[(1 << 0)+(1 << 1)] == -1
    assert np.count_nonzero(transferred) == 1


def test_actual_simple_charged_sectors_have_only_an_abelian_internal_commutant():
    result = r.finite_response()
    certificate = r.charged_sector_symmetry(result)
    assert certificate["strict_order_proof"] == {"exp(b)": "8/5", "exp(a-b)": "125/64",
                                                "both_exact_ratios_exceed_one": True}
    for number in (1, 3):
        indices = [mask for mask in range(16) if mask.bit_count() == number]
        h = result.hamiltonian[np.ix_(indices, indices)]
        energies, basis = np.linalg.eigh(h)
        np.testing.assert_allclose(energies, certificate["unshifted_energies"], atol=2e-15)
        commutator_map = np.kron(np.eye(4), h)-np.kron(h.T, np.eye(4))
        assert 16-np.linalg.matrix_rank(commutator_map, tol=1e-12) == 4
        projectors = [np.outer(column, column.conj()) for column in basis.T]
        assert max(np.max(abs(a@b-b@a)) for a in projectors for b in projectors) < 1e-15
    assert all(rank == 1 for ranks in certificate["charged_residue_ranks_at_visibility_tolerance"].values()
               for rank in ranks)
    assert not certificate["nontrivial_energy_and_charge_preserving_color_triplet_by_relabeling"]
    assert not certificate["composite_or_refined_constructions_excluded"]


def test_ten_work_bits_four_references_and_fixed_total_charge_continue_without_reset():
    rng = np.random.default_rng(1064)
    state = rng.normal(size=(1024, 4))+1j*rng.normal(size=(1024, 4))
    allowed = np.array([(mask & 255).bit_count() == 4 for mask in range(1024)])
    state[~allowed] = 0
    state /= np.linalg.norm(state)
    result = r.transfer_read(state, quadrature="Y")
    assert all(output.shape == (1024, 4) for output in result.outputs)
    assert all(np.linalg.norm(output[~allowed]) < 5e-14 for output in result.outputs)
    np.testing.assert_allclose(sum(output.conj().T@output for output in result.outputs),
                               state.conj().T@state, atol=5e-14)
    # No old outcome or untouched reference is conditioned away.
    assert sum(np.vdot(output, output).real for output in result.outputs) == pytest.approx(1, abs=5e-13)
    assert result.resources["work_qubits"] == 10
    assert result.resources["reference_columns"] == 4
    assert result.resources["work_preparations"] == 0


def test_zero_branch_retained_invalid_inputs_rejected_and_data_immutable():
    outputs = r.transfer_read(np.zeros((1024, 4), complex)).outputs
    assert len(outputs) == 2 and all(np.count_nonzero(value) == 0 for value in outputs)
    result = r.finite_response()
    with pytest.raises(ValueError):
        result.ground.setflags(write=True)
    for invalid in (True, math.nan, math.inf, "0.2", [0.2]):
        with pytest.raises(ValueError):
            r.finite_response(chemical_potential=invalid)
    for invalid in (0., 1.-.1j, math.nan+1j):
        with pytest.raises(ValueError):
            r.green_function(result, invalid)
    for state, site_a, site_b, axis in ((np.zeros(16), 0, 0, "X"),
                                       (np.zeros(256), True, 0, "X"),
                                       (np.zeros(256), 0, 4, "X"),
                                       (np.zeros(256), 0, 0, "Q")):
        with pytest.raises(ValueError):
            r.transfer_read(state, site_a, site_b, quadrature=axis)
