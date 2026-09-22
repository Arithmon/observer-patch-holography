"""Tests reject changed statistics, hidden top cutoffs and promoted premises."""
import copy
import importlib.util
from pathlib import Path
import sys
import pytest

HERE = Path(__file__).resolve().parent


def imported(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


producer = imported("pauli_stability_producer", HERE/"build.py")
verifier = imported("pauli_stability_verifier", HERE/"verify.py")


def test_exact_independent_replay_and_committed_receipt():
    receipt = producer.build()
    assert receipt == verifier.load(HERE/"receipt.json")
    result = verifier.verify(receipt)
    assert result["stable_signed_full_ladder_q"] == -1
    assert result["finite_dimensional_q"] == -1
    assert not result["source_premises_derived"]


def test_word_normal_order_is_not_assuming_CAR():
    assert producer.normal("ac") == {"": (1,), "ca": (0, 1)}
    assert producer.expectation("aacc") == (1, 1)
    assert producer.expectation("aacacc") == (1, 2, 1)
    assert producer.value(producer.expectation("aacc"), producer.F(1)) == 2


def test_intermediate_q_cannot_keep_exact_number_generator():
    receipt = producer.build()
    receipt["word_certificate"]["admissible_real_q_roots"] = [-1, 0, 1]
    receipt["word_certificate"]["q_controls"][2]["generator_residual_expectation"] = "0"
    with pytest.raises(ValueError, match="scalar-q"):
        verifier.verify(receipt)


def test_wrong_residual_cannot_select_q_by_renaming_roots():
    receipt = producer.build()
    receipt["word_certificate"]["exact_generator_residual_expectation"] = [0, -1, 1]
    receipt["word_certificate"]["admissible_real_q_roots"] = [0, 1]
    with pytest.raises(ValueError, match="scalar-q"):
        verifier.verify(receipt)


def test_positive_shifted_generator_is_not_same_source_hopping():
    receipt = producer.build()
    edge = receipt["source_edge"]
    for i in range(4):
        edge["hamiltonian"][i][i] = ["1", "0"]
    for mode in edge["orthonormal_eigenmodes"]:
        mode["energy"] = str(producer.F(mode["energy"])+1)
    with pytest.raises(ValueError, match="source edge"):
        verifier.verify(receipt)


def test_wrong_eigenmode_phase_rejected():
    receipt = producer.build()
    receipt["source_edge"]["orthonormal_eigenmodes"][0]["vector"][2][1] = "-1/2"
    with pytest.raises(ValueError, match="source edge"):
        verifier.verify(receipt)


def test_neutral_species_census_cannot_be_truncated():
    receipt = producer.build()
    receipt["neutral_bosonic_descent"]["multiplets"].pop()
    with pytest.raises(ValueError, match="neutral ladder"):
        verifier.verify(receipt)


def test_charge_chemical_potential_cannot_stabilize_neutral_sequence():
    receipt = producer.build()
    receipt["neutral_bosonic_descent"]["charge_chemical_potential_changes_energy"] = True
    with pytest.raises(ValueError, match="neutral ladder"):
        verifier.verify(receipt)


def test_finite_occupation_samples_are_not_all_n_proof():
    receipt = producer.build()
    receipt["scope"]["finite_samples_prove_unboundedness"] = True
    with pytest.raises(ValueError, match="scientific scope"):
        verifier.verify(receipt)


def test_filled_sea_shift_cannot_change_relative_excitation():
    receipt = producer.build()
    receipt["CAR_two_mode_sea"]["states"][0]["sea_shifted_energy"] = "0"
    with pytest.raises(ValueError, match="CAR sea"):
        verifier.verify(receipt)


def test_truncated_boson_does_not_satisfy_one_scalar_q():
    receipt = producer.build()
    control = receipt["truncated_boson_control"]
    assert control["double_creation_vacuum_norm_squared"] == "2"
    control["no_single_real_q_satisfies_scalar_relation"] = False
    with pytest.raises(ValueError, match="truncated boson"):
        verifier.verify(receipt)


def test_truncated_wrong_gram_cannot_hide_boundary():
    receipt = producer.build()
    receipt["truncated_boson_control"]["gram_diagonal"] = ["1", "1", "1"]
    with pytest.raises(ValueError, match="truncated boson"):
        verifier.verify(receipt)


def test_finite_dimensional_CCR_trace_must_fail():
    receipt = producer.build()
    receipt["finite_dimension_alternative"]["finite_trace_controls"][1]["commutator_trace"] = 2
    with pytest.raises(ValueError, match="finite matter dimension"):
        verifier.verify(receipt)


def test_A1_response_dimension_cannot_be_promoted_to_matter_dimension():
    receipt = producer.build()
    receipt["finite_dimension_alternative"]["A1_finite_response_implies_this_matter_realization"] = True
    with pytest.raises(ValueError, match="finite matter dimension"):
        verifier.verify(receipt)


@pytest.mark.parametrize("key", ["physical_spin_statistics_theorem", "A1_A3_derive_CAR",
    "all_statistics_or_all_observer_models_excluded", "native_preparation_or_physical_energy_identified",
    "full_gauge_constrained_source_Hamiltonian"])
def test_scope_promotions_rejected(key):
    receipt = producer.build()
    receipt["scope"][key] = True
    with pytest.raises(ValueError, match="scientific scope"):
        verifier.verify(receipt)


@pytest.mark.parametrize("key", ["one_common_real_scalar_q_for_all_modes",
    "number_is_a_dagger_a_with_exact_unit_phase_generator",
    "unrestricted_ladder_domain_and_fixed_signed_generator"])
def test_indispensable_premise_cannot_be_erased(key):
    receipt = producer.build()
    receipt["premises"][key] = False
    with pytest.raises(ValueError, match="selector premises"):
        verifier.verify(receipt)


def test_number_shift_cannot_claim_same_Cayley_factor():
    receipt = producer.build()
    receipt["escapes"]["uniform_number_shift"]["same_step_Cayley_map_inherited"] = True
    with pytest.raises(ValueError, match="escape boundaries"):
        verifier.verify(receipt)


def test_duplicate_receipt_keys_rejected(tmp_path):
    path = tmp_path/"duplicate.json"
    path.write_text('{"schema":"bad","schema":"oph.pauli_stability.v1"}')
    with pytest.raises(ValueError, match="duplicate"):
        verifier.load(path)
