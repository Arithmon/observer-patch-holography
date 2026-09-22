"""Adversarial controls for statistics selection and retained premises."""
import copy
import importlib.util
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parent


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


producer = module("spin_exchange_producer", HERE/"build.py")
verifier = module("spin_exchange_verifier", HERE/"verify.py")


def test_independent_tensor_replay_and_committed_receipt():
    receipt = producer.build()
    assert receipt == verifier.load(HERE/"receipt.json")
    assert verifier.verify(receipt)["coincidence_gap"] == "576/625"


@pytest.mark.parametrize("wrong_model", ["boson", "distinguishable"])
def test_coherent_normalized_wrong_statistics_rejected(wrong_model):
    receipt = producer.build()
    receipt["results"]["fermion"] = copy.deepcopy(receipt["results"][wrong_model])
    # All wrong-model probabilities remain normalized; the 1-RDM is unchanged.
    assert sum(map(producer.F, receipt["results"]["fermion"]["probabilities"].values())) == 1
    assert receipt["results"]["fermion"]["one_body_density"] == receipt["results"]["boson"]["one_body_density"]
    with pytest.raises(ValueError, match="tensor probabilities"):
        verifier.verify(receipt)


def test_boson_factorial_normalization_cannot_be_dropped():
    receipt = producer.build()
    p = receipt["results"]["boson"]["probabilities"]
    p["20"] = p["02"] = "144/625"
    p["11"] = "337/625"  # Renormalize to the distinguishable answer.
    with pytest.raises(ValueError, match="tensor probabilities"):
        verifier.verify(receipt)


def test_two_particle_preparation_cannot_be_replaced_by_parent_single_channel_state():
    receipt = producer.build()
    receipt["preparation"]["two_particles_in_same_channel"] = False
    with pytest.raises(ValueError, match="preparation"):
        verifier.verify(receipt)


def test_unbound_but_unitary_identity_splitter_rejected():
    receipt = producer.build()
    receipt["splitter"] = [[["1", "0"], ["0", "0"]], [["0", "0"], ["1", "0"]]]
    with pytest.raises(ValueError, match="parent hopping"):
        verifier.verify(receipt)


def test_missing_parent_hash_rejected():
    receipt = producer.build()
    receipt["source_pins"].pop("code/sm_fermion_current/quantum_link.py")
    with pytest.raises(ValueError, match="source pins"):
        verifier.verify(receipt)


def test_spin_central_sign_does_not_select_exchange():
    receipt = producer.build()
    receipt["spin_rotation"]["two_particle_sign_in_both_statistics"] = -1
    with pytest.raises(ValueError, match="spin rotation"):
        verifier.verify(receipt)


@pytest.mark.parametrize("field", ["physical_spin_statistics_theorem", "A1_A3_no_go",
    "source_selected_preparation_or_statistics", "native_or_laboratory_detector_outcomes",
    "parent_quantum_link_operator_Gauss_inherited", "continuum_or_interacting_QFT_result"])
def test_scope_promotion_rejected(field):
    receipt = producer.build()
    receipt["scope"][field] = True
    with pytest.raises(ValueError, match="scientific scope"):
        verifier.verify(receipt)


def test_widened_error_budgets_cannot_keep_separation_claim():
    receipt = producer.build()
    receipt["error_contract"]["preparation_trace_distance_at_most"] = "1/2"
    receipt["error_contract"]["total_probability_error_at_most"] = "51/100"
    with pytest.raises(ValueError, match="error contract"):
        verifier.verify(receipt)


def test_duplicate_json_keys_rejected(tmp_path):
    path = tmp_path/"duplicate.json"
    path.write_text('{"schema": "bad", "schema": "oph.spin_exchange.v1"}')
    with pytest.raises(ValueError, match="duplicate"):
        verifier.load(path)


def test_boolean_cannot_be_replaced_with_numeric_one():
    receipt = producer.build()
    receipt["scope"]["finite_declared_statistics_models"] = 1
    with pytest.raises(ValueError, match="scientific scope"):
        verifier.verify(receipt)


def test_partial_distinguishability_controls():
    receipt = producer.build()
    controls = receipt["partial_distinguishability"]
    assert [x["coincidence_gap"] for x in controls] == ["0", "144/625", "576/625"]
    assert controls[0]["results"]["fermion"] == controls[0]["results"]["boson"]
    assert verifier.verify(receipt)["verified"]


def test_orthogonal_hidden_labels_cannot_keep_full_exchange_gap():
    receipt = producer.build()
    controls = receipt["partial_distinguishability"]
    controls[0]["results"] = copy.deepcopy(controls[2]["results"])
    controls[0]["coincidence_gap"] = "576/625"
    with pytest.raises(ValueError, match="partial distinguishability"):
        verifier.verify(receipt)


def test_hidden_overlap_is_squared_not_amplitude():
    receipt = producer.build()
    receipt["partial_distinguishability"][1]["squared_hidden_state_overlap"] = "1/2"
    with pytest.raises(ValueError, match="partial distinguishability"):
        verifier.verify(receipt)
