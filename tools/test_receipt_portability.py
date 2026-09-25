from __future__ import annotations

import json
from pathlib import Path

import pytest

import check_receipt_portability as portability

from check_receipt_portability import check, find_violations


@pytest.mark.parametrize(
    "leaked",
    [
        "/Users/alice/work/oph/receipt.json",
        "/home/alice/work/oph/receipt.json",
        r"C:\Users\Alice\work\oph\receipt.json",
    ],
)
def test_developer_home_paths_fail_closed(tmp_path: Path, leaked: str) -> None:
    receipt = tmp_path / "receipt.json"
    receipt.write_text(
        json.dumps({"artifact": "mutation", "source": {"path": leaked}}),
        encoding="utf-8",
    )
    violations = find_violations([receipt])
    assert len(violations) == 1
    assert violations[0].json_pointer == "/source/path"
    with pytest.raises(ValueError, match="developer-home path leaked"):
        check([receipt])


def test_clone_stable_references_pass(tmp_path: Path) -> None:
    receipt = tmp_path / "receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "repo_relative": "code/particles/runs/flavor/receipt.json",
                "sibling": "oph-workspace://oph-physics-sim/runs/e1/receipt.json",
                "external": "external-file://measurement-pack.json",
            }
        ),
        encoding="utf-8",
    )
    check([receipt])


def archived_fixture(tmp_path, monkeypatch):
    archive = portability.ROOT / portability.OBSERVER_ARCHIVE
    source = next((archive / 'sim-analysis/data/refine_ensemble').glob('chain_*.json'))
    target = tmp_path / source.relative_to(portability.ROOT)
    target.parent.mkdir(parents=True)
    target.write_bytes(source.read_bytes())
    manifest = tmp_path / portability.OBSERVER_ARCHIVE / 'manifest.json'
    manifest.write_bytes((archive / 'manifest.json').read_bytes())
    monkeypatch.setattr(portability, 'ROOT', tmp_path)
    return target, manifest


def test_exact_archived_provenance_is_portable_metadata(tmp_path, monkeypatch):
    target, _ = archived_fixture(tmp_path, monkeypatch)
    check([target])


@pytest.mark.parametrize('mutation', ['value', 'field', 'name', 'manifest', 'manifest_and_file'])
def test_original_provenance_exception_rejects_changed_identity(tmp_path, monkeypatch, mutation):
    target, manifest = archived_fixture(tmp_path, monkeypatch)
    if mutation in ('value', 'field', 'manifest_and_file'):
        data = json.loads(target.read_bytes())
        data['parent_run' if mutation != 'field' else 'runtime_path'] = '/home/changed/not-an-original'
        target.write_text(json.dumps(data))
    if mutation == 'name':
        replacement = target.with_name('new_receipt.json')
        target.rename(replacement)
        target = replacement
    if mutation in ('manifest', 'manifest_and_file'):
        data = json.loads(manifest.read_bytes())
        data['description'] = 'changed custody'
        if mutation == 'manifest_and_file':
            import hashlib
            rel = target.relative_to(tmp_path / portability.OBSERVER_ARCHIVE).as_posix()
            data['files'][rel]['sha256'] = hashlib.sha256(target.read_bytes()).hexdigest()
            data['files'][rel]['bytes'] = target.stat().st_size
        manifest.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='developer-home path leaked'):
        check([target])
