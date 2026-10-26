"""Charged spectral weight, clocked dispersion and Euclidean calibration."""
import json
import math

import numpy as np
import pytest

from ..chiral_candidate import action
from . import mass_observable as m


def test_existing_source_clock_supplies_the_inverse_frequency_denominator():
    clock = m._clock_record()
    omega = float(action.source_coefficients()["Omega"])
    assert clock["decoded_phase"] == pytest.approx(omega*float(m.INTERVAL), abs=1e-14)
    assert clock["purification_reference_columns_retained"] == 4
    assert clock["clock_is_actual_source_execution"]
    assert clock["preparation_resources"]["native_pulses"] > 0


def test_actual_charge_projected_CAR_residues_give_clocked_mass_shell_on_three_axes():
    rows = m.charged_pole_diagnostics()
    assert len(rows) == 9
    squared = {family: [] for family in range(3)}
    for row in rows:
        assert row["electric_charge_of_Nambu_branch"] == -1
        assert 0 < row["positive_pole_residue"] < 1
        assert 0 < row["negative_pole_residue"] < 1
        assert row["positive_pole_residue"]+row["negative_pole_residue"] == pytest.approx(1, abs=1e-14)
        assert row["two_signed_pole_response_residual"] < 1e-14
        assert row["invariant_residual"] < 2e-14
        assert row["clock_decoded_energy_over_Omega"] == pytest.approx(row["pole_energy_over_Omega_diagnostic"], abs=1e-13)
        squared[row["native_family_axis"]].append(row["invariant_E_over_Omega_squared_minus_p_squared_over4"])
    for values in squared.values():
        np.testing.assert_allclose(values, values[0], atol=2e-14, rtol=0)
    # Spin projection produces unequal signed residues on z; testing only the
    # transverse axes would miss a mistaken equal-residue implementation.
    assert any(abs(row["positive_pole_residue"]-.5) > .1 for row in rows)
    # The correct observable reads the existing modest hierarchy; it does not
    # change the chosen action's mass ratios by a new interpretation.
    rest = [math.sqrt(values[0]) for values in squared.values()]
    assert max(rest)/min(rest) < 1.5


def test_literal_ground_charged_Euclidean_decay_has_mass_not_blank_quench_factor():
    report = m.rest_pair_comparison()
    mass = report["source_mass_parameter"]
    np.testing.assert_allclose(report["Fock_energies"], [-mass, 0, 0, mass], atol=1e-14)
    assert report["ground_electric_charge"] == pytest.approx(0, abs=1e-14)
    assert report["addition_electric_charge"] == pytest.approx(-1, abs=1e-14)
    assert report["removal_electric_charge"] == pytest.approx(1, abs=1e-14)
    assert report["addition_residue"] == pytest.approx(.5, abs=1e-14)
    assert report["removal_residue"] == pytest.approx(.5, abs=1e-14)
    for row in report["ground_Euclidean_correlations"]:
        assert row["addition_correlator"] == pytest.approx(row["same_H_prediction"], abs=1e-14)
        assert row["removal_correlator"] == pytest.approx(row["same_H_prediction"], abs=1e-14)
    assert report["Euclidean_rate_equals_charged_gap_residual"] < 1e-14
    assert report["source_Lorentzian_anticommutator_cos_mt_residual"] < 1e-13
    assert report["blank_quench_number_contrast"] == pytest.approx(report["blank_quench_prediction_cos_2mt"], abs=1e-13)
    assert abs(report["blank_quench_number_contrast"]-math.cos(mass*float(m.INTERVAL))) > 1e-3


def test_same_positive_H_neutral_Euclidean_rate_is_different_from_charged_mass():
    report = m.rest_pair_comparison()
    for row in report["ground_Euclidean_correlations"]:
        assert row["neutral_connected_number_correlator"] == pytest.approx(row["neutral_same_H_prediction"], abs=1e-14)
    assert report["neutral_number_commutes_electric_charge_residual"] == 0
    assert report["neutral_decay_equals_twice_charged_gap_residual"] < 1e-14
    assert report["same_H_neutral_number_decay_rate"] == pytest.approx(2*report["Euclidean_rate_of_same_H_charged_correlator"], abs=1e-14)
    assert not report["Euclidean_contraction_executed_as_native_repair"]


def test_propagating_CAR_response_is_recovered_from_completed_native_words():
    packet = m.source_propagating_check()
    coarse, fine = packet["native_source_convergence"]
    assert packet["work_CAR_modes_at_p_and_minus_p"] == 8
    assert fine["full_CAR_linear_map_max_residual"] < coarse["full_CAR_linear_map_max_residual"]*.26
    assert fine["rest_squared_residual"] < coarse["rest_squared_residual"]*.26
    assert fine["rest_squared_residual"] < 1e-6
    orthogonal = np.array(packet["final_source_Majorana_map"])
    np.testing.assert_allclose(orthogonal@orthogonal.T, np.eye(16), atol=3e-12, rtol=0)
    source = packet["final_source_execution"]
    assert source["resources"]["routing_fswap_calls"] > 0
    assert source["resources"]["retained_work_fermion_modes"] == 8
    assert source["diagnostics"]["all_local_Fock_columns_contracted"]
    assert not source["target_BdG_exponential_used_by_source"]
    assert not packet["charged_probe_instrument_for_this_new_momentum_menu_executed"]


def test_report_preserves_physical_and_native_operator_boundaries():
    packet = m.report()
    json.dumps(packet, allow_nan=False)
    assert packet["largest_invariant_residual"] < 2e-14
    assert packet["smallest_charged_pole_residue"] > 0
    assert not packet["native_repair_operator_boundary"]["same_operator_or_same_history_intertwiner_established"]
    assert not packet["native_repair_operator_boundary"]["native_CP_relaxation_eigenvalue_equals_mass_proved"]
    assert not packet["scope"]["new_native_interacting_mass_generation_derived"]
    assert not packet["scope"]["observed_Standard_Model_mass_spectrum_identified"]
