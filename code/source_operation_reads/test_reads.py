"""Frozen source identity, independent semantics and hostile artifact checks."""
from copy import deepcopy
import hashlib
from pathlib import Path
import subprocess
import sys

import pytest

from . import simulator_verifier as check
from .verify import HERE, verify


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
