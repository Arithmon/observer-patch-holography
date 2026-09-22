"""Frozen source identity, independent semantics and hostile artifact checks."""
from copy import deepcopy
import hashlib
from pathlib import Path
import subprocess
import sys

import pytest

from . import simulator_verifier as check
from .verify import HERE, verify
from .check_live import canonicalize_source_aliases


@pytest.fixture(scope="module")
def packet():
    return check.strict_load(HERE/"capture.json")


def reseal(value):
    value["sha256"] = hashlib.sha256(check.canonical({k: v for k, v in value.items()
                                                    if k != "sha256"})).hexdigest()
    return value


def test_retained_receipt_and_observation_limits(packet):
    result = verify(packet)
    check.equal(result, check.strict_load(HERE/"receipt.json"), "retained receipt")
    cases = result["native_derivation"]["cases"]
    assert [r["commits"] for r in cases] == [8, 16, 24, 24, 48, 96, 24]
    assert [cases[i]["ideal_real_read_algebra"]["hidden_dimension"] for i in (2, 4, 5)] == [21, 65, 153]
    assert all(r["repair_read_after_write_edges"] == 0 for r in cases)
    assert result["M1_derived"] is result["complete_A1_A3_model"] is False


@pytest.mark.parametrize("path,value", [
    (("source", "revision"), "0"*40),
    (("source", "files", "oph_fpe/core/echosahedral_dynamics.py"), "0"*64),
    (("cases", 0, "commits"), 0),
    (("cases", 2, "observer_records", 0, 4, 0), "0x0.0p+0"),
    (("cases", 2, "observer_records", 0, 3), 11),
    (("cases", 0, "order"), []),
    (("quantum", "unitary_hex", 0, 0), ["0x0.0p+0", "0x0.0p+0"]),
    (("quantum", "uniform_snapshot_feedback", "native_next_port"), 0),
    (("scope", "M1_derived"), True),
    (("cases",), []),
])
def test_resealed_semantic_and_custody_forgery(packet, path, value):
    bad = deepcopy(packet)
    target = bad
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        verify(reseal(bad))


def test_verifier_is_bound_to_pinned_simulator_bytes(packet, monkeypatch, tmp_path):
    import source_operation_reads.verify as wrapper
    (tmp_path/"source_snapshot.json").write_bytes((HERE/"source_snapshot.json").read_bytes())
    (tmp_path/"simulator_verifier.py").write_bytes(b"def verify(x): return True\n")
    monkeypatch.setattr(wrapper, "HERE", tmp_path)
    with pytest.raises(ValueError, match="verifier source pin"):
        wrapper.verify(packet)


@pytest.mark.parametrize("raw", ['{}', '[]', '{"x":0,"x":1}', '{"x":NaN}', '{"x":1.0}'])
def test_real_cli_rejects_trash(raw, tmp_path):
    path = tmp_path/"bad.json"
    path.write_text(raw, encoding="ascii")
    run = subprocess.run([sys.executable, "-m", "source_operation_reads.verify", "--packet", str(path)],
                         capture_output=True, text=True, timeout=30)
    assert run.returncode != 0
    assert '"verified":true' not in run.stdout


def test_independent_core_has_no_producer_or_numerical_imports(packet, tmp_path):
    path = tmp_path/"capture.json"
    path.write_bytes(check.canonical(packet))
    script = '''
import importlib.abc, runpy, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'oph_fpe','numpy','scipy'}:
            raise AssertionError('forbidden import: '+fullname)
sys.meta_path.insert(0, Block())
module=runpy.run_path(sys.argv[1],run_name='independent')
assert module['verify'](module['strict_load'](sys.argv[2]))['M1_derived'] is False
'''
    run = subprocess.run([sys.executable, "-c", script, str(Path(check.__file__)), str(path)],
                         capture_output=True, text=True, timeout=30)
    assert run.returncode == 0, run.stderr


@pytest.fixture
def source_alias(tmp_path):
    canonical = "oph_fpe/core/patchnet.py"
    alias = "oph_fpe/core/PatchNet.py"
    path = tmp_path/canonical
    path.parent.mkdir(parents=True)
    path.write_bytes(b"class PatchNet: pass\n")
    alias_path = tmp_path/alias
    if not alias_path.exists():
        alias_path.hardlink_to(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    source = {"repository": "pinned", "revision": "a"*40, "files": {canonical: digest}}
    packet = reseal({"source": {**source, "files": {canonical: digest, alias: digest}}})
    return packet, source, tmp_path, canonical, alias


def test_live_source_alias_requires_same_file_and_preserves_packet(source_alias):
    packet, source, root, _, _ = source_alias
    original = deepcopy(packet)
    result = canonicalize_source_aliases(packet, source, root)
    assert result["source"] == source
    assert result == reseal(deepcopy(result))
    assert packet == original


@pytest.mark.parametrize("mutation", ["missing_pin", "unknown_import", "wrong_hash", "bad_seal"])
def test_live_alias_does_not_relax_source_custody(source_alias, mutation):
    packet, source, root, canonical, alias = source_alias
    if mutation == "missing_pin":
        del packet["source"]["files"][canonical]
    elif mutation == "unknown_import":
        packet["source"]["files"]["oph_fpe/core/other.py"] = source["files"][canonical]
    elif mutation == "wrong_hash":
        packet["source"]["files"][alias] = "0"*64
    if mutation != "bad_seal":
        reseal(packet)
    else:
        packet["sha256"] = "0"*64
    with pytest.raises(ValueError):
        canonicalize_source_aliases(packet, source, root)


def test_live_case_variant_separate_file_is_rejected(source_alias, monkeypatch):
    packet, source, root, canonical, alias = source_alias
    # Exercise the case-sensitive filesystem branch on every CI platform.
    # Identical bytes alone cannot authorize deletion from the source closure.
    original = Path.samefile
    def separate(path, other):
        return False if path == root/alias else original(path, other)
    monkeypatch.setattr(Path, "samefile", separate)
    with pytest.raises(ValueError, match="same pinned file"):
        canonicalize_source_aliases(packet, source, root)


def test_live_source_alias_rejects_changed_bytes(source_alias):
    packet, source, root, canonical, _ = source_alias
    (root/canonical).write_bytes(b"class PatchNet: changed = True\n")
    with pytest.raises(ValueError, match="same pinned file"):
        canonicalize_source_aliases(packet, source, root)
