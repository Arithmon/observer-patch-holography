"""Actual stochastic current law and its restricted mean-field interpretation."""
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest

from ..observer import native
from . import native_fields as f


@pytest.fixture(scope="module")
def result():
    return f.run()


def test_every_actual_branch_current_and_conditional_mean_checked(result):
    assert result["primitive"]["primitive_branches_checked"] == 78 * 30 * 2
    assert result["primitive"]["conditional_current_means_checked"] == 78 * 30
    assert result["primitive"]["markov_rows_checked"] == 78
    assert result["all_pass"]


def test_full_markov_transition_is_independently_built_from_all_primitive_events():
    states = native.states()
    for state, row in zip(states, f.transition()):
        dense = [F(0)] * 78
        for seam in f.geometry.source().seams():
            for coin in (0, 1):
                final, _ = native.native(state, seam, coin)
                dense[states.index(final)] += F(1, 60)
        assert tuple((i, p) for i, p in enumerate(dense) if p) == row
        assert sum(dense) == 1


def test_markov_multi_step_means_match_separate_matrix_powers():
    kernel = np.zeros((78, 78))
    for i, row in enumerate(f.transition()):
        for j, probability in row:
            kernel[i, j] = float(probability)
    states = np.array(native.states())
    operator = np.array(f.mean_operator(), dtype=float)
    for n in (1, 2, 5, 11):
        np.testing.assert_allclose(np.linalg.matrix_power(kernel, n) @ states,
                                   states @ np.linalg.matrix_power(operator, n).T,
                                   rtol=0, atol=2e-15)


def test_incidence_factorization_dissipation_and_constant_mean():
    edges = f.geometry.source().seams()
    incidence = np.zeros((12, 30), dtype=int)
    for index, (a, b) in enumerate(edges):
        incidence[a, index], incidence[b, index] = 1, -1
    np.testing.assert_array_equal(incidence @ incidence.T, f.laplacian())
    operator, lap = f.mean_operator(), f.laplacian()
    for m in native.states() + ((F(1, 6),) * 12,):
        evolved = f._matvec(operator, m)
        lm = f._matvec(lap, m)
        drop = sum(x*x for x in m) - sum(x*x for x in evolved)
        expected = sum(m[i]*lm[i] for i in range(12))/30 - sum(x*x for x in lm)/3600
        assert drop == expected >= 0
        assert sum(evolved) == sum(m)


def test_conditional_path_surprisal_is_exact_additive_and_separates_start_prior(result):
    audit = result["stochastic_action"]
    assert audit["conditional_path_probability"] == "product_k K[x_k,x_(k+1)]"
    assert audit["local_surprisal"] == "-log K[x,y]=log(60/m[x,y]) nats"
    assert "separate boundary term" in audit["start_prior"]
    assert audit["exact_transition_multiplicity_histogram"] == {"1": 600, "2": 60, "50": 48, "52": 30}
    assert audit["local_surprisal_classes_nats"] == {
        "1": "log(60)", "2": "log(30)", "50": "log(6/5)", "52": "log(15/13)"}
    assert audit["fully_labeled_seam_coin_branch_probability"] == "1/60"
    assert audit["fully_labeled_branch_surprisal_nats"] == "log(60), constant per event"
    assert audit["not_proportional_to_quadratic_load_drop"]
    assert audit["not_the_logarithm_of_the_kernel_as_a_linear_operator"]
    assert not audit["physical_action_units_derived"]
    most_likely = audit["one_step_most_likely_vs_mean"]
    assert most_likely["most_likely_next_state"] == most_likely["initial_state"]
    assert most_likely["most_likely_transition_probability"] == "13/15"
    assert most_likely["conditional_mean_after_one_event"] == [
        "14/15", "14/15", "0", "0", "0", "1/30", "0", "1/30", "1/60", "1/60", "1/60", "1/60"]
    assert most_likely["greedy_surprisal_path_is_not_the_conditional_mean_update"]
    examples = audit["state_only_path_examples"]
    repair = examples["concentration_repair"]["local_terms"][0]
    assert repair["conditional_probability"] == "1/30"
    assert repair["surprisal_nats"] == "log(30)"
    assert repair["quadratic_load_potential_drop"] == "1"
    stationary = examples["split_state_self_transition"]["local_terms"][0]
    rearrangement = examples["split_state_zero_drop_move"]["local_terms"][0]
    assert stationary["quadratic_load_potential_drop"] == rearrangement["quadratic_load_potential_drop"] == "0"
    assert stationary["surprisal_nats"] == "log(15/13)"
    assert rearrangement["surprisal_nats"] == "log(60)"
    assert not examples["split_state_self_transition"]["initial_state_prior_included"]


def test_conditional_path_surprisal_rejects_impossible_or_malformed_paths():
    x = (2,) + (0,) * 11
    impossible = (0, 2) + (0,) * 10
    with pytest.raises(ValueError, match="impossible native transition"):
        f.conditional_path_surprisal((x, impossible))
    with pytest.raises(ValueError, match="at least two"):
        f.conditional_path_surprisal((x,))


def test_mean_equilibrium_is_not_a_static_native_state(result):
    iteration = result["iteration"]
    assert iteration["uniform_78_mean_is_preserved"]
    assert not iteration["uniform_78_distribution_is_stationary"]
    assert iteration["uniform_66_binary_distribution_is_stationary"]
    assert iteration["stationary_binary_fluctuation_norm_squared"] == "5/3"
    first_chain = iteration["chains"][0]
    norms = [F(row["centered_mean_norm_squared"]) for row in first_chain["samples"]]
    assert all(a > b for a, b in zip(norms, norms[1:]))
    assert all(row["squared_mean_norm"] == "1/3" for row in iteration["chains"][1]["samples"])


def test_spectrum_and_oscillation_boundary_are_exact(result):
    spectrum = result["spectrum"]
    q = (55 + math.sqrt(5))/60
    assert spectrum["slow_multiplier_numeric"] == pytest.approx(q)
    assert spectrum["slow_symmetric_second_difference_coefficient_numeric"] == pytest.approx(q + 1/q - 2)
    assert spectrum["slow_symmetric_second_difference_coefficient_numeric"] > 0
    expected = np.linalg.eigvalsh(np.asarray(f.mean_operator(), dtype=float))
    actual = sorted(value["mean_multiplier_numeric"] for value in spectrum["bands"] for _ in range(value["multiplicity"]))
    np.testing.assert_allclose(expected, actual, atol=3e-15, rtol=0)
    assert all(0 < value < 1 for value in actual[:-1]) and actual[-1] == 1
    for angle in (.1, 1., math.pi, 2*math.pi):
        assert -4*math.sin(angle/2)**2 <= 0


def test_changed_schedule_or_coin_breaks_the_declared_uniform_fair_formula(result):
    control = result["premise_controls"]
    assert control["biased_conditional_current"] == "1/4"
    assert control["fair_conditional_current"] == "1/2"
    assert control["fixed_seam_fair_mean"] != control["uniform_seams_fair_mean"]
    assert not result["scope"]["canonical_uniform_path_law_derived_from_OPH_axioms"]
    assert not result["scope"]["other_native_history_readouts_excluded"]
    assert not result["scope"]["physical_clock_or_spatial_continuum_derived"]
    json.dumps(result, allow_nan=False)


def test_existing_actual_feedback_runs_are_not_uniform_all_seam_calibration(result):
    control = result["premise_controls"]
    assert control["actual_E_F_feedback_two_event_mean"] == ["1/4", "3/4", "1"] + ["0"] * 9
    assert control["actual_E_F_feedback_two_event_mean"] != control["uniform_calibration_two_event_mean"]
    paths = control["actual_E_F_feedback_paths"]
    assert len(paths) == 4 and all(row["probability"] == "1/4" for row in paths)
    assert {tuple(row["second_selected_seam"]) for row in paths} == {native.E, native.F}
    assert all(row["second_selected_seam"] == list(native.F if row["first_record"][1] > 0 else native.E) for row in paths)
    assert not result["scope"]["uniform_calibration_is_actual_E_F_feedback_protocol"]


@pytest.mark.parametrize("distribution", ((F(1, 77),) * 77, (0,) * 78, (1.0,) + (0,) * 77, (-1, 2) + (0,) * 76))
def test_nonexact_or_invalid_distribution_rejected(distribution):
    with pytest.raises(ValueError):
        f.advance(distribution)


def test_cli_round_trip_rejects_changed_current_result(tmp_path, monkeypatch, capsys):
    path = tmp_path / "native_fields.json"
    monkeypatch.setattr("sys.argv", ["native_fields", "--output", str(path)])
    f.main()
    monkeypatch.setattr("sys.argv", ["native_fields", "--check", str(path)])
    f.main()
    result = json.loads(path.read_text())
    result["premise_controls"]["biased_conditional_current"] = "1/2"
    path.write_text(json.dumps(result))
    with pytest.raises(ValueError, match="Receipt mismatch"):
        f.main()
    capsys.readouterr()
