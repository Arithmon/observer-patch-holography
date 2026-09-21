"""False-green controls for source binding, algebra scope and chronology."""
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import pytest

HERE = Path(__file__).resolve().parent


def module(name):
    spec = importlib.util.spec_from_file_location(f"source_operator_join_{name}", HERE / f"{name}.py")
    obj = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = obj
    spec.loader.exec_module(obj)
    return obj


producer = module("join")
consumer = module("verify")


@pytest.fixture
def packet():
    return producer.produce()


def test_committed_packet_reproduces_and_independent_verifier_checks(packet):
    assert producer.canonical(packet) == (HERE / "operator_join_packet.json").read_bytes()
    result = consumer.verify(packet)
    assert result["transition_carrier_points"] == 31 * 2366
    assert result["chronological_intervals"] == 496
    assert result["physical_adequacy"] is False
    assert result["observed_quantum_outcomes"] is False
    assert result["lean_kernel_replay_performed"] is False
    assert result["universal_operator_theorems_verified_by_python"] is False


def test_utf8_source_and_packet_replay_with_cp1252_default(monkeypatch, tmp_path):
    original_read_text = Path.read_text

    def windows_read_text(path, encoding=None, errors=None, **kwargs):
        return original_read_text(path, encoding=encoding or "cp1252", errors=errors, **kwargs)

    monkeypatch.setattr(Path, "read_text", windows_read_text)
    # Confirm that the simulation reaches the actual Windows failure mode.
    with pytest.raises(UnicodeDecodeError):
        (producer.ROOT / producer.PARENTS[0]).read_text()
    packet = producer.produce()
    assert producer.canonical(packet) == (HERE / "operator_join_packet.json").read_bytes()
    assert consumer.verify(consumer.load(HERE / "operator_join_packet.json"))["source_rows"] == 32
    path = tmp_path / "unicode.json"
    value = {"label": "∀ ψ ∈ ℂ"}
    path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")
    assert consumer.load(path) == value


def test_missing_adjacent_transition(packet):
    del packet["transitions"][8]
    with pytest.raises(consumer.Rejected, match="31 actual"):
        consumer.verify(packet)


def test_invented_terminal_edge(packet):
    packet["transitions"].append({"row": 31, "next_row": 0, "permutations": []})
    with pytest.raises(consumer.Rejected, match="31 actual"):
        consumer.verify(packet)


def test_wrong_observer_factor_transport(packet):
    packet["transitions"][0]["permutations"][0] = list(range(13))
    with pytest.raises(consumer.Rejected, match="transition0"):
        consumer.verify(packet)


def test_reordered_transition_history(packet):
    packet["transitions"][0], packet["transitions"][1] = packet["transitions"][1], packet["transitions"][0]
    with pytest.raises(consumer.Rejected, match="transition0"):
        consumer.verify(packet)


def test_identity_is_not_a_selected_spectator_state(packet):
    packet["embeddings"]["unused_factor"] = "projector_onto_zero"
    with pytest.raises(consumer.Rejected, match="unused_factor"):
        consumer.verify(packet)


def test_wrong_hinge_axis(packet):
    packet["embeddings"]["common_axes"] = [1]
    with pytest.raises(consumer.Rejected, match="common_axes"):
        consumer.verify(packet)


def test_record_diagonal_cannot_replace_full_hinge(packet):
    packet["carrier"]["common_complex_vector_dimension"] = 13
    with pytest.raises(consumer.Rejected, match="common_complex_vector_dimension"):
        consumer.verify(packet)


def test_pair_marginal_cannot_be_replaced_by_uniform_state(packet):
    packet["state"]["pair88_counts"][0][-1] += 1
    with pytest.raises(consumer.Rejected, match="pair88_counts"):
        consumer.verify(packet)


@pytest.mark.parametrize("field", [
    "overlapping_pair_algebras_commute", "physical_regions", "physical_clock",
    "observed_quantum_outcomes", "continuum_refinement", "scalar64_instrument_same_carrier",
    "universal_algebraic_pushout",
])
def test_forbidden_scope_promotion(packet, field):
    packet["interpretation"][field] = True
    with pytest.raises(consumer.Rejected, match=field):
        consumer.verify(packet)


def test_single_path_projector_cannot_replace_source_generators(packet):
    packet["interpretation"]["admitted_generation"] = "singleton_path_projector"
    with pytest.raises(consumer.Rejected, match="admitted_generation"):
        consumer.verify(packet)


def test_tower_state_cannot_be_identified_with_counted_state(packet):
    packet["interpretation"]["tower_state"] = "same_as_empirical"
    with pytest.raises(consumer.Rejected, match="tower_state"):
        consumer.verify(packet)


def test_boolean_cannot_spoof_zero_row(packet):
    packet["transitions"][0]["row"] = False
    with pytest.raises(consumer.Rejected, match="row: type"):
        consumer.verify(packet)


def test_integer_cannot_spoof_false_scope_flag(packet):
    packet["interpretation"]["physical_clock"] = 0
    with pytest.raises(consumer.Rejected, match="physical_clock: type"):
        consumer.verify(packet)


def test_rehashed_changed_parent_is_not_accepted(packet, tmp_path):
    for path in consumer.PINS:
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(consumer.ROOT / path, target)
    path = "Lean/QFT/SourceOperatorGeneration.lean"
    target = tmp_path / path
    target.write_text(target.read_text(encoding="utf-8").replace(
        "path := ![1,2,2,7", "path := ![0,2,2,7", 1), encoding="utf-8")
    with pytest.raises(consumer.Rejected, match="pinned source changed"):
        consumer.verify(packet, tmp_path)
    packet["parents"][path] = hashlib.sha256(target.read_bytes()).hexdigest()
    with pytest.raises(consumer.Rejected, match="parents"):
        consumer.verify(packet, tmp_path)


def test_duplicate_json_key_is_rejected(tmp_path):
    path = tmp_path / "duplicate.json"
    path.write_text('{"schema": "a", "schema": "b"}', encoding="utf-8")
    with pytest.raises(consumer.Rejected, match="duplicate key"):
        consumer.load(path)


def test_sparse_matrix_unit_controls_are_noncommutative():
    assert consumer.sparse_controls() == 18


def test_spectator_trace_is_not_multiplicative():
    # Normalized trace over observer 247: E01 and E10 separately vanish,
    # while their product has normalized trace 1/13. This precludes treating
    # either observable retraction or state marginal as a star homomorphism.
    x = consumer.sparse_unit((2,), (0,), (1,))
    y = consumer.sparse_unit((2,), (1,), (0,))
    diagonal = lambda op: sum(value for (a, b), value in op.items() if a == b)
    assert diagonal(x) == diagonal(y) == 0
    assert diagonal(consumer.multiply(x, y)) == 182
