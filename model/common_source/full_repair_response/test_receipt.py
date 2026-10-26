"""Receipt tampering must not turn a negative physical result into a claim."""
import copy
import json

import numpy as np
import pytest

from . import run


def test_portable_roundoff_does_not_relax_physical_or_exact_claims():
    original = {"physical_mass_identified": False, "rank": 14,
                "exact_bound": "1/100", "source_sha256": "a" * 64,
                "residual": 1e-14, "spread": 25.0}
    replay = dict(original, residual=2e-14, spread=25.0 + 1e-10)
    run.compare_replay(original, replay)
    for key, value in (("physical_mass_identified", True), ("rank", 15),
                       ("exact_bound", "1/101"), ("source_sha256", "b" * 64),
                       ("spread", 206.0), ("rank", 14.0),
                       ("residual", float("nan"))):
        changed = dict(replay, **{key: value})
        with pytest.raises(ValueError):
            run.compare_replay(original, changed)
    with pytest.raises(ValueError):
        run.compare_replay(original, {key: value for key, value in replay.items()
                                      if key != "physical_mass_identified"})


@pytest.fixture
def receipt(tmp_path, monkeypatch):
    arrays = {"work": np.array([[1.0, 0.0], [0.0, 1.0j]]),
              "terminal_minus": np.zeros((2, 2), dtype=complex)}
    artifact = tmp_path / "response.npz"
    np.savez_compressed(artifact, **arrays)
    report = {"schema": "oph.full_repair_response.v1",
              "physical_mass_identified": False}
    retained = {**report, "artifact": {"filename": artifact.name,
                                       "sha256": run._hash(artifact)}}
    path = tmp_path / "response.json"
    path.write_text(json.dumps(retained))
    monkeypatch.setattr(run, "run", lambda: (copy.deepcopy(report), arrays))
    return path, artifact, arrays, retained


def test_receipt_accepts_the_complete_recomputed_work(receipt):
    run.verify(receipt[0])


def test_changed_work_bytes_fail_before_recomputation(receipt, monkeypatch):
    path, artifact, _, _ = receipt
    artifact.write_bytes(artifact.read_bytes() + b"altered")
    def forbidden_replay():
        pytest.fail("altered artifact should be rejected before recomputation")
    monkeypatch.setattr(run, "run", forbidden_replay)
    with pytest.raises(ValueError, match="artifact changed"):
        run.verify(path)


def test_rehashed_missing_terminal_outcome_is_still_rejected(receipt):
    path, artifact, arrays, retained = receipt
    np.savez_compressed(artifact, work=arrays["work"])
    retained["artifact"]["sha256"] = run._hash(artifact)
    path.write_text(json.dumps(retained))
    with pytest.raises(ValueError, match="missing retained work"):
        run.verify(path)


def test_rehashed_work_cannot_change_the_recomputed_history(receipt):
    path, artifact, arrays, retained = receipt
    changed = dict(arrays, work=arrays["work"] * 0.9)
    np.savez_compressed(artifact, **changed)
    retained["artifact"]["sha256"] = run._hash(artifact)
    path.write_text(json.dumps(retained))
    with pytest.raises(ValueError, match="differs from recomputation"):
        run.verify(path)


def test_receipt_cannot_promote_the_physical_verdict(receipt):
    path, _, _, retained = receipt
    retained["physical_mass_identified"] = True
    path.write_text(json.dumps(retained))
    with pytest.raises(ValueError, match="claim or exact value changed"):
        run.verify(path)
