"""Custody checks reject damaged receipts; they do not replace a full replay."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from . import run


@pytest.fixture(scope="module")
def receipt():
    return json.loads((Path(__file__).with_name("results.json")).read_text())


def test_retained_receipt_passes_custody_and_record_checks(capsys):
    run.verify(Path(__file__).with_name("results.json"))
    assert "does not rerun the experiment" in capsys.readouterr().out


def damage(document, kind):
    case = document["cases"]["interacting_n64"]
    if kind == "missing_control":
        del document["cases"]["pinned_work"]
    elif kind == "missing_source":
        del document["source_bindings"]["model/common_source/mass_response/model.py"]
    elif kind == "changed_source":
        document["source_bindings"]["model/common_source/mass_response/model.py"] = "0"*64
    elif kind == "escaped_source":
        document["source_bindings"]["../outside.py"] = "0"*64
    elif kind == "missing_normal_outcome":
        del case["record_probabilities"]["00000"]
    elif kind == "missing_abort":
        del case["source_abort_probabilities"][next(iter(case["source_abort_probabilities"]))]
    elif kind == "negative_probability":
        case["record_probabilities"]["00000"] = -1e-4
    elif kind == "lost_mass":
        case["record_probabilities"]["00000"] += .01
    elif kind == "nonfinite":
        case["complete_quantum_output_half_trace_distance"] = float("nan")
    elif kind == "wrong_readout":
        case["receiver_occupation_probability"] += .01
    elif kind == "wrong_decoded_gap":
        case["observer_gap_inverse"]["finite_single_charged_gap_over_clock"] += .01
    elif kind == "wrong_signal":
        document["matched_intervention_signals"]["no_kick"]["separation_after_ideal_split_budgets"] += .01
    elif kind == "over_budget":
        case["complete_quantum_output_half_trace_distance"] = 1.
    elif kind == "particle_promotion":
        document["claim_scope"]["standard_model_particle_identification"] = True
    elif kind == "missing_joint_path":
        document["same_history_native_feedback"]["actual"]["normal_paths"] -= 1
    elif kind == "wrong_native_selector":
        document["same_history_native_feedback"]["kick_probability_exact"] = "1/2"
    elif kind == "no_feedback":
        document["same_history_native_feedback"]["final_native_classical_total_variation"] = 0.
    else:
        raise AssertionError(kind)


@pytest.mark.parametrize("kind", ["missing_control", "missing_source", "changed_source", "escaped_source",
                                 "missing_normal_outcome", "missing_abort", "negative_probability", "lost_mass",
                                 "nonfinite", "wrong_readout", "wrong_decoded_gap", "wrong_signal", "over_budget",
                                 "particle_promotion", "missing_joint_path", "wrong_native_selector", "no_feedback"])
def test_damaged_receipt_is_rejected(receipt, tmp_path, kind):
    altered = deepcopy(receipt)
    damage(altered, kind)
    path = tmp_path/"damaged.json"
    path.write_text(json.dumps(altered))
    with pytest.raises(ValueError):
        run.verify(path)


def test_duplicate_json_keys_are_rejected(tmp_path):
    path = tmp_path/"duplicate.json"
    path.write_text('{"schema":"x", "schema":"y"}')
    with pytest.raises(ValueError, match="duplicate"):
        run.verify(path)
