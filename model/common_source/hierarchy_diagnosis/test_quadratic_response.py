"""Independent finite pole, residue, blank-quench and reference checks."""
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import quadratic_response as q
from ..chiral_candidate import gaussian_source


def test_exact_finite_CAR_identities_pin_the_pole_and_frequency_distinction():
    proof = q.exact_identities()
    assert proof["equal_time_CAR_exact"]
    assert proof["Fock_energies"] == ["-|d|", "0", "0", "+|d|"]
    assert proof["ground_addition_and_removal_weights_exact"] == ["1/2", "1/2"]
    assert proof["normal_local_residues_exact"] == ["1/2", "1/2"]
    assert proof["blank_quench_frequency_is_twice_CAR_gap"]


@pytest.mark.parametrize("coefficient", [.173, -.21, .17+.23j])
@pytest.mark.parametrize("time", [.37, -.62])
def test_retarded_anticommutator_is_same_operator_for_distinct_reference_states(coefficient, time):
    result = q.state_independence_check(coefficient, time)
    assert result["operator_anticommutator_residual"] < 2e-15
    assert set(result["reference_expectation_residuals"]) == {
        "encoded_blank", "illustrative_pair_ground", "maximally_mixed", "nonGaussian_even_pair_mixture"}
    assert max(result["reference_expectation_residuals"].values()) < 2e-15
    assert result["reference_stationarity_residuals"]["encoded_blank"] > .1
    assert result["reference_stationarity_residuals"]["illustrative_pair_ground"] < 1e-15
    # Reference independence is specific to the retarded CAR response.
    greater = result["normal_greater_correlations_do_depend_on_reference"]
    assert np.linalg.norm(np.array(greater["encoded_blank"])-greater["illustrative_pair_ground"]) > .1


def test_actual_native_four_column_pair_word_has_the_same_retarded_operator():
    coefficient, time = .17+.23j, .37
    h, nambu, xi = q.pair_matrices(coefficient)
    source = gaussian_source.local_gate("pairing", coefficient, time)
    prediction = expm(-1j*time*nambu)
    for i, left in enumerate(xi):
        evolved = source.fock_unitary.conj().T@left@source.fock_unitary
        for j, right in enumerate(xi):
            anticommutator = evolved@right.conj().T+right.conj().T@evolved
            np.testing.assert_allclose(anticommutator, prediction[i, j]*np.eye(4), atol=3e-14, rtol=0)
    np.testing.assert_allclose(source.fock_unitary, expm(-1j*time*h), atol=3e-14, rtol=0)


def test_pair_report_pins_ground_charged_residues_and_quench_frequency():
    row = q.pair_report(.17+.23j)
    magnitude = abs(.17+.23j)
    np.testing.assert_allclose(row["Fock_energy_eigenvalues"], [-magnitude, 0, 0, magnitude], atol=1e-15)
    assert row["ground_to_odd_gap"] == pytest.approx(magnitude, abs=1e-15)
    assert row["ground_to_even_pair_gap"] == pytest.approx(2*magnitude, abs=1e-15)
    np.testing.assert_allclose(row["normal_pole_residues"], [.5, .5], atol=2e-15)
    assert row["blank_pair_occupation_cosine_residual"] < 2e-15


def test_report_uses_actual_all_three_lepton_pairs_without_an_empirical_target():
    report = q.report()
    rows = report["actual_lepton_pair_rows"]
    assert {row["native_family_axis"] for row in rows} == {0, 1, 2}
    assert len({row["charged_CAR_gap_exact"] for row in rows}) == 3
    gaps = [row["charged_CAR_gap"] for row in rows]
    assert 1 < max(gaps)/min(gaps) < 1.5
    assert max(row["maximum_normal_residue_error_from_one_half"] for row in rows) < 2e-15
    assert not report["interpretation"]["same_quadratic_generator_reference_change_moves_linear_retarded_poles"]
    assert not report["interpretation"]["physical_interacting_vacuum_selected"]
    assert report["empirical_mass_inputs"] == []
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("invalid", [True, 0., math.nan, math.inf, "0.2", [0.2]])
def test_invalid_pair_coefficient_is_rejected(invalid):
    with pytest.raises(ValueError):
        q.pair_matrices(invalid)
