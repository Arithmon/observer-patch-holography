"""Check retained scientific evidence and fail-closed provenance boundaries."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from model.observer_particles import run


def test_retained_computations_and_premise_views_agree():
    result = run.verify()
    assert result == json.loads((run.HERE / "results/verification.json").read_text())
    witnesses = run.premise_witnesses()
    assert witnesses == json.loads((run.HERE / "results/premise_witnesses.json").read_text())
    assert len(result["producers"]) == 16
    assert result["physical_mass_predictions_promoted"] is False
    for row in witnesses["rows"]:
        assert row["producer"] in result["producers"]
        assert (run.ROOT / row["receipt"]).is_file()
        assert row["physical_premise_removed"] is False


@pytest.mark.parametrize("filename", [
    "final_end_to_end_predictions.json", "blind_prediction_provenance.json",
    "derivation_chain_closure_matrix.json", "theorem_gap_register.json",
    "postdiction_ledger.json",
])
def test_particle_ledgers_share_the_verified_finite_witnesses(filename):
    record = json.loads((run.ROOT / "code/particles/runs/status" / filename).read_text())
    assert record["finite_model_witnesses"] == run.premise_witnesses()


@pytest.fixture
def isolated_receipt(tmp_path, monkeypatch):
    """Use real files to exercise changed inputs without editing evidence."""
    here = tmp_path / "model/observer_particles"
    (here / "results").mkdir(parents=True)
    producer = here / "entangled_alpha.py"
    source = tmp_path / "model/source.py"
    artifact = here / "results/operator.bin"
    producer.write_text("# isolated producer\n")
    source.write_text("# isolated computational input\n")
    artifact.write_bytes(b"finite operator")
    receipt = {
        "status": "PASS_FINITE_CONTRACT",
        "producer_sha256": run.sha(producer),
        "source_sha256": {"model/source.py": run.sha(source)},
        "artifacts": {"operator.bin": run.sha(artifact)},
        "evidence_contract": {"physical_particle_masses_source_selected": False},
    }
    record = here / "results/entangled_alpha.json"
    record.write_text(json.dumps(receipt))
    monkeypatch.setattr(run, "HERE", here)
    monkeypatch.setattr(run, "ROOT", tmp_path)
    monkeypatch.setattr(run, "PRODUCERS", ("entangled_alpha",))
    assert run.verify()["artifact_bindings"] == 1
    return producer, source, artifact, record


@pytest.mark.parametrize("changed", [0, 1, 2])
def test_changed_producer_input_or_operator_is_rejected(isolated_receipt, changed):
    isolated_receipt[changed].write_bytes(b"changed")
    with pytest.raises(AssertionError):
        run.verify()


def test_receipt_cannot_promote_physical_masses(isolated_receipt):
    record = isolated_receipt[3]
    receipt = json.loads(record.read_text())
    receipt["evidence_contract"]["physical_particle_masses_source_selected"] = True
    record.write_text(json.dumps(receipt))
    with pytest.raises(AssertionError):
        run.verify()


def test_receipt_cannot_bind_a_file_outside_repository(isolated_receipt):
    record = isolated_receipt[3]
    receipt = json.loads(record.read_text())
    receipt["source_sha256"] = {"../outside.py": "0" * 64}
    record.write_text(json.dumps(receipt))
    with pytest.raises(AssertionError):
        run.verify()
