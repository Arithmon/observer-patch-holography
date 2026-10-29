"""Independent receipt custody and hostile-output checks for exchange statistics."""
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from tools import strict_json
from . import statistics_packet as packet


@pytest.fixture(scope="module")
def replay():
    return packet.report()


def test_canonical_receipt_replays_every_field(replay):
    packet.compare(strict_json.load(packet.OUTPUT), replay)
    assert replay["physical_particle_identification"] is False
    assert replay["physical_comparison_count_added"] == 0
    assert replay["integrated_gate_accepted"] is False


def _leaves(value, path=()):
    if type(value) is dict:
        for key, item in value.items():
            yield from _leaves(item, path+(key,))
    elif type(value) is list:
        for i, item in enumerate(value):
            yield from _leaves(item, path+(i,))
    else:
        yield path, value


def test_every_receipt_leaf_and_missing_subtree_are_rejected(replay):
    checked = 0
    for path, value in _leaves(replay):
        changed = deepcopy(replay)
        parent = changed
        for key in path[:-1]:
            parent = parent[key]
        parent[path[-1]] = (not value if type(value) is bool else
                            value+1 if type(value) is int else "altered")
        with pytest.raises(ValueError):
            packet.compare(changed, replay)
        checked += 1
    assert checked > 100
    for key in replay:
        changed = deepcopy(replay)
        del changed[key]
        with pytest.raises(ValueError):
            packet.compare(changed, replay)


@pytest.mark.parametrize("text", ('{"x":0,"x":0}', '{"x":NaN}', '{"x":Infinity}',
                                 '{"x":-Infinity}', '{"x":false}', '{"x":0.0}', '[]'))
def test_malformed_numbers_duplicate_keys_and_types_cannot_pass(text):
    with pytest.raises((TypeError, ValueError)):
        packet.compare(strict_json.loads(text), {"x": 0})


def test_required_geometry_and_package_initializers_are_bound(replay):
    bound = replay["source_sha256"]
    loaded = set()
    for name, module in tuple(sys.modules.items()):
        if not name.startswith(("model", "tools", "oph_exact", "oph_fpe")):
            continue
        filename = getattr(module, "__file__", None)
        if filename and (path := Path(filename).resolve()).is_relative_to(packet.ROOT):
            # Other pytest test modules may already be imported. Bind runtime
            # implementation modules, not an unrelated test collection order.
            if not path.name.startswith("test_"):
                loaded.add(path.relative_to(packet.ROOT).as_posix())
    # The focused CLI below independently tests a clean interpreter; a broad
    # suite may legitimately have loaded unrelated implementation packages.
    assert packet.OWN_SOURCES[0] in bound
    for name in ("model/__init__.py", "model/exchange_sector/__init__.py",
                 "evidence/native_geometric_source_20260925/spec.json",
                 "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
                 "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/__init__.py"):
        assert name in bound
    assert any(name.endswith("oph_exact/carrier.py") for name in loaded)


def test_clean_process_replay_binds_its_entire_loaded_import_closure(tmp_path):
    script = """
import json, pathlib, sys
from model.exchange_sector import statistics_packet as p
value = p.report()
p.compare(json.loads(p.OUTPUT.read_text()), value)
loaded = []
for name, module in tuple(sys.modules.items()):
    filename = getattr(module, '__file__', None)
    if filename and name.startswith(('model', 'tools', 'oph_exact', 'oph_fpe')):
        path = pathlib.Path(filename).resolve()
        if path.is_relative_to(p.ROOT):
            loaded.append(path.relative_to(p.ROOT).as_posix())
missing = set(loaded) - value['source_sha256'].keys()
assert not missing, sorted(missing)
print('complete clean-process exchange replay passed')
"""
    env = dict(os.environ, PYTHONUTF8="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    result = subprocess.run([sys.executable, "-X", "utf8", "-c", script],
                            cwd=packet.ROOT, env=env, capture_output=True, text=True,
                            encoding="utf-8", timeout=300, check=False)
    assert result.returncode == 0, result.stdout+result.stderr
    assert "complete clean-process exchange replay passed" in result.stdout


def test_contract_scope_cannot_be_promoted(tmp_path, monkeypatch):
    value = strict_json.load(packet.HERE / "statistics_contract.json")
    value["physical_particle_identification"] = True
    (tmp_path / "statistics_contract.json").write_text(json.dumps(value), encoding="utf-8")
    monkeypatch.setattr(packet, "HERE", tmp_path)
    with pytest.raises(ValueError, match="contract changed"):
        packet.contract()


def test_source_changes_during_execution_are_rejected(replay, monkeypatch):
    pins = replay["source_sha256"]
    calls = []

    def changing():
        calls.append(1)
        return pins if len(calls) == 1 else {**pins, "changed.py": "0"*64}

    monkeypatch.setattr(packet, "source_bindings", changing)
    monkeypatch.setattr(packet.statistics_topology, "certificate", lambda: replay["topology"])
    monkeypatch.setattr(packet.statistics, "certificate", lambda: replay["instruments"])
    monkeypatch.setattr(packet.statistics_service, "certificate", lambda: replay["execution"])
    with pytest.raises(ArithmeticError, match="changed during replay"):
        packet.report()


def test_warm_geometry_cache_cannot_accept_an_invalid_source_manifest(monkeypatch):
    geometry = packet.slow_geometry
    geometry.source()
    assert geometry.source.cache_info().currsize == 1
    spec = geometry.SPEC
    altered = strict_json.load(spec)
    altered['source_sha256']['oph_exact/carrier.py'] = '0'*64
    altered_bytes = json.dumps(altered).encode('utf-8')
    original_bytes, original_text = Path.read_bytes, Path.read_text

    def changed_bytes(path):
        return altered_bytes if path == spec else original_bytes(path)

    def changed_text(path, *args, **kwargs):
        return altered_bytes.decode('utf-8') if path == spec else original_text(path, *args, **kwargs)

    def premature_calculation():
        raise AssertionError('geometry must be validated before scientific calculation')

    monkeypatch.setattr(Path, 'read_bytes', changed_bytes)
    monkeypatch.setattr(Path, 'read_text', changed_text)
    for owner in (packet.statistics_topology, packet.statistics, packet.statistics_service):
        monkeypatch.setattr(owner, 'certificate', premature_calculation)
    with pytest.raises(ValueError, match='pinned native geometry source has changed'):
        packet.report()


@pytest.mark.parametrize("name", (
    "evidence/native_geometric_source_20260925/spec.json",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/__init__.py",
    "code/source_selection_model/response.json",
))
def test_actual_source_bytes_control_standalone_custody(name, replay, monkeypatch):
    original = Path.read_bytes
    target = packet.ROOT / name

    def changed(path):
        value = original(path)
        return value+b"\n" if path == target else value

    monkeypatch.setattr(Path, "read_bytes", changed)
    bindings = packet.source_bindings()
    assert bindings[name] != replay["source_sha256"][name]
    with pytest.raises(ValueError):
        packet.compare({**replay, "source_sha256": bindings}, replay)
