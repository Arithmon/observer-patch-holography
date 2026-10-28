"""Check retained scientific evidence and fail-closed provenance boundaries."""
from __future__ import annotations

import json
import ast
import copy
from pathlib import Path, PureWindowsPath

import pytest

from model.observer_particles import run


@pytest.mark.parametrize('producer', (*run.PRODUCERS, 'run'))
def test_receipt_and_operator_writes_are_canonical_lf(producer, tmp_path, monkeypatch):
    """Exercise the actual write calls with a Windows-like default translation."""
    original = Path.open

    def windows_open(self, *args, **kwargs):
        if kwargs.get('mode', args[0] if args else 'r') == 'w':
            if kwargs.get('newline') is None:
                kwargs['newline'] = '\r\n'
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Path, 'open', windows_open)
    control = tmp_path / 'default.txt'
    control.write_text('first\nsecond\n')
    assert control.read_bytes() == b'first\r\nsecond\r\n'
    tree = ast.parse((run.HERE / (producer + '.py')).read_text(encoding='utf-8'))
    writes = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
              and isinstance(node.func, ast.Attribute) and node.func.attr == 'write_text']
    assert writes
    for index, node in enumerate(writes):
        call = copy.deepcopy(node)
        call.func.value = ast.Name(id='target', ctx=ast.Load())
        call.args = [ast.Constant('first\nsecond\n')]
        expression = ast.fix_missing_locations(ast.Expression(call))
        target = tmp_path / f'{index}.json'
        eval(compile(expression, str(run.HERE / (producer + '.py')), 'eval'),
             {'target': target})
        assert target.read_bytes() == b'first\nsecond\n'


def test_windows_paths_preserve_complete_public_source_closure(monkeypatch):
    from model.observer_particles import evidence_contract as contract

    sources = ('model/native_modular_reference/local_preparation.py',
               'model/observer_particles/geometric_stationarity.py')
    bindings = {str(PureWindowsPath(path)): run.sha(run.ROOT / path)
                for path in sources}
    original = Path.relative_to
    with monkeypatch.context() as patch:
        patch.setattr(Path, 'relative_to',
                      lambda self, *args, **kwargs: PureWindowsPath(
                          original(self, *args, **kwargs)))
        result = contract.public_receipt({'source_sha256': bindings})
    expected = {path: run.sha(run.ROOT / path) for path in sources}
    helper = 'model/observer_particles/evidence_contract.py'
    expected[helper] = run.sha(run.ROOT / helper)
    assert result['source_sha256'] == expected


def test_windows_source_path_still_checks_computational_hash():
    from model.observer_particles import evidence_contract as contract

    path = str(PureWindowsPath('model/native_modular_reference/local_preparation.py'))
    with pytest.raises(ValueError, match='Changed computational input'):
        contract.public_receipt({'source_sha256': {path: '0' * 64}})


def test_dot_source_alias_is_canonicalized_and_hash_checked():
    from model.observer_particles import evidence_contract as contract

    path = 'model/native_modular_reference/local_preparation.py'
    expected = run.sha(run.ROOT / path)
    result = contract.public_receipt({'source_sha256': {'./' + path: expected}})
    assert result['source_sha256'][path] == expected
    assert './' + path not in result['source_sha256']
    with pytest.raises(ValueError, match='Changed computational input'):
        contract.public_receipt({'source_sha256': {'./' + path: '0' * 64}})


@pytest.mark.parametrize('bad_alias_first', [False, True])
def test_every_duplicate_source_alias_hash_is_checked(bad_alias_first):
    from model.observer_particles import evidence_contract as contract

    path = 'model/native_modular_reference/local_preparation.py'
    entries = [(path, run.sha(run.ROOT / path)), ('./' + path, '0' * 64)]
    if bad_alias_first:
        entries.reverse()
    with pytest.raises(ValueError, match='Changed computational input'):
        contract.public_receipt({'source_sha256': dict(entries)})


def test_geometric_source_collection_uses_portable_windows_paths(monkeypatch):
    monkeypatch.syspath_prepend(str(run.HERE))
    from model.observer_particles import geometric_stationarity as producer

    original = Path.relative_to
    with monkeypatch.context() as patch:
        patch.setattr(Path, 'relative_to',
                      lambda self, *args, **kwargs: PureWindowsPath(
                          original(self, *args, **kwargs)))
        result = producer.source_hashes()
    assert all('\\' not in key for key in result)
    for path in ('model/observer_particles/geometric_stationarity.py',
                 'model/native_modular_reference/local_preparation.py'):
        assert result[path] == run.sha(run.ROOT / path)


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
