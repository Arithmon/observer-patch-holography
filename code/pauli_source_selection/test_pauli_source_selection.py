"""Adversarial checks of missing source-selection conditions."""
import copy
import importlib.util
from pathlib import Path
import sys
import pytest

HERE = Path(__file__).resolve().parent


def imported(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE/filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


producer = imported("pauli_source_producer", "build.py")
verifier = imported("pauli_source_verifier", "verify.py")


@pytest.fixture(scope="module")
def original():
    return producer.build()


def test_independent_labeled_tensor_verification(original):
    assert original == verifier.load(HERE/"receipt.json")
    result = verifier.verify(original)
    assert result["verified"] and result["dimension"] == 15
    assert not result["source_premises_derived"]


@pytest.mark.parametrize("mode", [0, 1, 2, 3])
def test_erasing_double_creation_fails_tensor_map(original, mode):
    receipt = copy.deepcopy(original)
    for row in range(5, 15):
        for col in range(1, 5):
            receipt["creation"][mode][row][col] = ["0", "0"]
    with pytest.raises(ValueError, match="creation intertwiner"):
        verifier.verify(receipt)


def test_replacing_factorial_gram_with_identity_fails(original):
    receipt = copy.deepcopy(original)
    receipt["gram_diagonal"] = ["1"]*15
    with pytest.raises(ValueError, match="tensor Gram"):
        verifier.verify(receipt)


def test_empty_hamiltonian_cannot_claim_same_source_lift(original):
    receipt = copy.deepcopy(original)
    receipt["quadratic_hamiltonian"] = [[["0", "0"] for _ in range(15)] for _ in range(15)]
    with pytest.raises(ValueError, match="quadratic generator"):
        verifier.verify(receipt)


def test_no_bosonic_scalar_q_promotion(original):
    receipt = copy.deepcopy(original)
    receipt["scalar_q_failure"]["one_common_q_exists"] = True
    with pytest.raises(ValueError, match="scalar-q"):
        verifier.verify(receipt)


def test_idle_outcome_cannot_be_erased(original):
    receipt = copy.deepcopy(original)
    receipt["three_outcome_instrument"]["idle_effect_diagonal"] = ["0"]*15
    receipt["three_outcome_instrument"]["exhaustive_two_adjoint_outcomes"] = True
    with pytest.raises(ValueError, match="three-outcome"):
        verifier.verify(receipt)


def test_scaled_kraus_cannot_inherit_unit_phase(original):
    receipt = copy.deepcopy(original)
    receipt["contractive_phase_route"]["capped_half_scaled_number_phase_coefficient"] = "1"
    with pytest.raises(ValueError, match="contractive phase"):
        verifier.verify(receipt)


def test_source_history_rank_cannot_be_reversed(original):
    receipt = copy.deepcopy(original)
    receipt["source_history_number_obstruction"]["available_target_dimension"] = 4
    with pytest.raises(ValueError, match="source-history"):
        verifier.verify(receipt)


def test_basis_only_records_cannot_select_relative_exchange_sign(original):
    receipt = copy.deepcopy(original)
    receipt["basis_only_control"]["occupation_record_probabilities_select_relative_sign"] = True
    with pytest.raises(ValueError, match="basis-only"):
        verifier.verify(receipt)


def test_individual_raisers_do_not_certify_superposition_phase(original):
    receipt = copy.deepcopy(original)
    row = receipt["basis_only_control"]["models"]["commuting_bits"]
    row["mixed_double_creation_vacuum_amplitude"] = ["0", "0"]
    row["mixed_phase_defect_on_mask1_amplitude"] = ["0", "0"]
    with pytest.raises(ValueError, match="basis-only"):
        verifier.verify(receipt)


@pytest.mark.parametrize("key", ["all_modes_and_all_unitaries_exhaustively_numerically_tested",
    "physical_spin_statistics_theorem",
    "scalar_q_law_derived_from_source", "full_A1_A3_countermodel", "native_physical_matter_or_instrument",
    "local_quantum_Gauss_or_relativistic_field_locality", "same_parent_fifteen_particle_state"])
def test_scope_promotions_rejected(original, key):
    receipt = copy.deepcopy(original)
    receipt["scope"][key] = True
    with pytest.raises(ValueError, match="scope"):
        verifier.verify(receipt)


def test_duplicate_json_keys_rejected(tmp_path):
    path = tmp_path/"duplicate.json"
    path.write_text('{"schema":0,"schema":1}')
    with pytest.raises(ValueError, match="duplicate"):
        verifier.load(path)


@pytest.mark.parametrize("malformed", [["0", "0", "7"], [False, "0"], ["0/2", "0"]])
def test_malformed_complex_entries_fail_closed(original, malformed):
    receipt = copy.deepcopy(original)
    receipt["creation"][0][0][0] = malformed
    with pytest.raises(ValueError, match="canonical"):
        verifier.verify(receipt)


@pytest.mark.parametrize("location", ["creation", "quadratic_hamiltonian", "lifted_cayley", "source_cayley"])
def test_zero_padded_matrix_shape_fails_closed(original, location):
    receipt = copy.deepcopy(original)
    matrix = (receipt["creation"][0] if location == "creation" else
              receipt["source_edge"]["cayley"] if location == "source_cayley" else receipt[location])
    dimension = len(matrix)
    for row in matrix:
        row.append(["0", "0"])
    matrix.append([["0", "0"] for _ in range(dimension+1)])
    with pytest.raises(ValueError, match="matrix dimension"):
        verifier.verify(receipt)
