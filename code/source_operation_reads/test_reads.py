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
    assert [r['terminal_real_projection_rank'] for r in cases] == [40,32,24,24,48,96,24]


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
    (("instruments", "extensions_selected_by_source"), True),
    (("instruments", "threshold_output_hex", 1, 0), "0x1.0000000000000p-1"),
    (("instruments", "phase_lift_coherence_hex", 0, 0), "0x0.0p+0"),
    (("instruments", "twirl_output_sparse"), []),
    (("instruments", "next_read_hex", 0), "0x0.0p+0"),
    (("instruments", "measured_recurrence", 2, "steps", 1, "positive_slot_entries"), 0),
    (("instruments", "measured_recurrence", 0, "erasure", "inverse_flow_hex"), ["0x0.0p+0"]*48),
    (("instruments", "retained_flags"), []),
    (("instruments", "retained_flags", 1, "recovered_digest"), "0"*64),
    (("instruments", "deterministic_archive", 2, "chord_differences"), []),
    (("instruments", "deterministic_archive", 2, "tree_seams"), [0]*15),
    (("instruments", "deterministic_archive", 2, "recovered_sha256"), "0"*64),
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


def test_constructed_channel_family_keeps_the_erasure_and_selection_boundary(packet):
    result=verify(packet)['native_derivation']
    assert result['instruments']['extensions_selected_by_source'] is False
    assert result['quantum']['coordinate_ports_central_in_reversible_closure'] is False
    for row in result['instruments']['measured_recurrence']:
        assert row['normalized_first_cycle_kernel_dimension']==5*row['carriers']+1
        assert row['normalized_first_cycle_image_dimension']==6*row['carriers']-1
        assert row['carrier_totals_are_Markov_state'] is False


def test_constructed_recovery_keeps_record_cost_and_native_boundary(packet):
    result=verify(packet)['native_derivation']['instruments']
    assert result['all_flag_words_recover_input_analytically'] is True
    assert result['retained_flag_bits_per_full_sweep']==24
    for row in result['deterministic_archive']:
        assert row['retained_real_coordinates']==row['linear_archive_minimum']==5*row['carriers']+1
        assert row['exact_input_recovered'] is True
        assert row['native_record_interface'] is False


@pytest.mark.parametrize('field', ['initial_hex','final_hex','observer_records'])
def test_live_bridge_rejects_different_values_even_with_identical_summaries(packet, field):
    from .check_live import compare_replays
    changed = deepcopy(packet)
    values = (changed['cases'][0]['observer_records'][0][4]
              if field == 'observer_records' else changed['cases'][0][field])
    values[0] = (float.fromhex(values[0])+2**-30).hex()
    # Pass the same summary deliberately: numeric reproduction is an
    # additional gate, not a consequence of rank/provenance summary equality.
    summary = verify(packet)['native_derivation']
    with pytest.raises(ValueError,match='live numeric reproduction'):
        compare_replays(changed,packet,summary,summary)


def test_live_bridge_accepts_last_bit_variation_and_checks_source_basis(packet):
    import math
    from .check_live import compare_replays
    summary = verify(packet)['native_derivation']
    changed = deepcopy(packet)
    values = changed['cases'][0]['initial_hex']
    values[0] = math.nextafter(float.fromhex(values[0]),math.inf).hex()
    compare_replays(changed,packet,summary,summary)
    changed['quantum']['phase_ports'] = [0,11]
    with pytest.raises(ValueError,match='live quantum source'):
        compare_replays(changed,packet,summary,summary)
