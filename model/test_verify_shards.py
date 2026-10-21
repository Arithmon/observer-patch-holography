"""CI partitions execution, never the canonical receipt or artifact coverage."""
from collections import Counter
import hashlib
import json
import sys

import pytest
import yaml

from . import verify


def flatten(groups):
    return tuple(replay for group in groups for replay in group)


@pytest.mark.parametrize('count', [1, 2, 3, 7])
def test_shards_are_exhaustive_disjoint_and_keep_empirical_worker_atomic(count):
    shards = [verify.selected_replay_groups(index, count) for index in range(count)]
    assert Counter(replay for groups in shards for replay in flatten(groups)) == Counter(verify.REPLAYS)
    assert sum(group == verify.EMPIRICAL_REPLAYS for groups in shards for group in groups) == 1
    assert all(group == verify.EMPIRICAL_REPLAYS or len(group) == 1 for groups in shards for group in groups)
    if count > 1:
        assert shards[0] == (verify.EMPIRICAL_REPLAYS,)
    assert all(shards)


def test_default_and_single_shard_keep_every_original_replay_in_order():
    assert flatten(verify.selected_replay_groups()) == verify.REPLAYS
    assert flatten(verify.selected_replay_groups(0, 1)) == verify.REPLAYS


def test_new_replays_are_included_without_a_parallel_allowlist():
    extra = ('model.example_new_certificate', None, 'model/example_new_certificate/results.json')
    declaration = (extra,)+verify.REPLAYS
    shards = [verify.selected_replay_groups(index, 3, replays=declaration) for index in range(3)]
    assert Counter(replay for groups in shards for replay in flatten(groups)) == Counter(declaration)
    assert extra in flatten(shards[1])


@pytest.mark.parametrize('index,count', [(None,3),(0,None),(-1,3),(3,3),(0,0),
                                       (0,10000),(True,3),(0,True),(0.,3),(0,3.)])
def test_invalid_shard_arguments_fail_closed(index, count):
    with pytest.raises(ValueError):
        verify.selected_replay_groups(index, count)


def test_malformed_empirical_group_cannot_silently_run_split():
    first = verify.REPLAYS.index(verify.EMPIRICAL_REPLAYS[0])
    changed = list(verify.REPLAYS)
    changed[first], changed[first+1] = changed[first+1], changed[first]
    with pytest.raises(ValueError, match='contiguous'):
        verify.replay_groups(changed)
    with pytest.raises(ValueError):
        verify.replay_groups(tuple(x for x in verify.REPLAYS if x not in verify.EMPIRICAL_REPLAYS))
    with pytest.raises(ValueError, match='unique'):
        verify.replay_groups(verify.REPLAYS+(verify.REPLAYS[0],))


@pytest.mark.parametrize('index', [0,1,2])
def test_every_shard_hashes_all_artifacts_and_executes_its_complete_groups(monkeypatch, index):
    monkeypatch.setattr(sys, 'argv', ['verify', '--shard-index', str(index), '--shard-count', '3'])
    hashed, expanded = [], []
    def hashes():
        hashed.append(True)
        return 534
    def run(argv, **kwargs):
        assert kwargs['check'] is True
        assert kwargs['cwd'] == verify.ROOT
        if argv == [sys.executable, '-m', 'model.replay_empirical']:
            expanded.extend(verify.EMPIRICAL_REPLAYS)
        else:
            position = argv.index('--check')
            expanded.append((argv[2], None if position == 3 else argv[3], argv[position+1]))
    monkeypatch.setattr(verify, 'verify_artifacts', hashes)
    monkeypatch.setattr(verify.subprocess, 'run', run)
    verify.main()
    assert hashed == [True]
    assert tuple(expanded) == flatten(verify.selected_replay_groups(index, 3))


def test_invalid_cli_shard_stops_before_hashes_or_execution(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['verify', '--shard-index', '3', '--shard-count', '3'])
    monkeypatch.setattr(verify, 'verify_artifacts', lambda: pytest.fail('invalid arguments reached hashes'))
    with pytest.raises(SystemExit) as error:
        verify.main()
    assert error.value.code == 2


def test_hash_only_shard_checks_all_artifacts_without_replays(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['verify', '--hashes-only', '--shard-index', '2', '--shard-count', '3'])
    calls = []
    monkeypatch.setattr(verify, 'verify_artifacts', lambda: calls.append('all') or 534)
    monkeypatch.setattr(verify.subprocess, 'run', lambda *a, **k: pytest.fail('hash-only executed a replay'))
    verify.main()
    assert calls == ['all']


@pytest.mark.parametrize('field', ['sources', 'source_sha256', 'source_bindings',
    'input_sha256', 'source_binding_sha256', 'source_control_bindings', 'source_pins',
    'implementation_sha256', 'source_interface_hashes'])
def test_hash_only_rejects_stale_inner_source_even_with_current_outer_manifest(
        monkeypatch, tmp_path, field):
    packet = tmp_path/'model/example'
    packet.mkdir(parents=True)
    source = packet/'source.py'
    source.write_text('value = 1\n', encoding='utf-8')
    retained = hashlib.sha256(source.read_bytes()).hexdigest()
    source_key = {'implementation_sha256': 'source.py',
                  'source_interface_hashes': 'example/source.py'}.get(field, 'model/example/source.py')
    receipt = packet/'results.json'
    receipt.write_text(json.dumps({'nested': [{field: {source_key: retained}}]}),
                       encoding='utf-8')
    source.write_text('value = 2\n', encoding='utf-8')
    # Both files have current outer hashes; the receipt's dependency binding
    # is stale and must fail before any scientific replay.
    manifest = {'files': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (source, receipt)}}
    (packet/'manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'PACKETS', ('model/example',))
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, 'model/example/results.json'),))
    monkeypatch.setattr(verify, 'selected_replay_groups', lambda *args: ())
    monkeypatch.setattr(sys, 'argv', ['verify', '--hashes-only'])
    monkeypatch.setattr(verify.subprocess, 'run', lambda *a, **k: pytest.fail('a replay executed'))
    with pytest.raises(ValueError, match=f'root.nested\\[0\\].{field}'):
        verify.main()
    receipt.write_text(json.dumps({'nested': [{field: {
        source_key: hashlib.sha256(source.read_bytes()).hexdigest()}}]}),
        encoding='utf-8')
    manifest['files']['results.json'] = hashlib.sha256(receipt.read_bytes()).hexdigest()
    (packet/'manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    verify.main()


def test_implementation_bindings_resolve_package_parent_without_root_fallback(monkeypatch, tmp_path):
    folder = tmp_path/'model/example'
    folder.mkdir(parents=True)
    source = folder.parent/'source.py'
    source.write_text('actual implementation\n', encoding='utf-8')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    receipt = folder/'results.json'
    receipt.write_text(json.dumps({'implementation_sha256': {'../source.py': digest}}),
                       encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, 'model/example/results.json'),))
    verify.verify_receipt_source_bindings()
    source.write_text('changed implementation\n', encoding='utf-8')
    # A matching file at another plausible root cannot hide the local drift.
    (tmp_path/'source.py').write_text('actual implementation\n', encoding='utf-8')
    with pytest.raises(ValueError, match='digest mismatch'):
        verify.verify_receipt_source_bindings()


@pytest.mark.parametrize('binding', [{}, [], False, {'missing.py': '0'*64},
    {'../../../outside.py': '0'*64}, {'https://example.org/source.py': '0'*64}])
@pytest.mark.parametrize('field', ['implementation_sha256', 'source_interface_hashes',
    'source_sha256', 'source_bindings', 'input_sha256', 'source_binding_sha256',
    'source_control_bindings', 'source_pins'])
def test_local_bindings_cannot_skip_malformed_or_unresolved_inputs(
        monkeypatch, tmp_path, binding, field):
    folder = tmp_path/'model/example'
    folder.mkdir(parents=True)
    receipt = folder/'results.json'
    receipt.write_text(json.dumps({field: binding}), encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, 'model/example/results.json'),))
    with pytest.raises(ValueError, match=field):
        verify.verify_receipt_source_bindings()


def test_endpoint_pins_are_checked_at_the_declared_nested_location(monkeypatch, tmp_path):
    source = tmp_path/'source.py'
    source.write_text('original\n', encoding='utf-8')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    receipt = tmp_path/'results.json'
    receipt.write_text(json.dumps({'retrospective_endpoint_audit': {
        'source': {'pins': {'source.py': digest}}}, 'unrelated': {'pins': {'alias': 'owner'}}}),
        encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, 'results.json'),))
    verify.verify_receipt_source_bindings()
    source.write_text('changed\n', encoding='utf-8')
    with pytest.raises(ValueError, match='root.retrospective_endpoint_audit.source.pins'):
        verify.verify_receipt_source_bindings()


def test_live_source_gate_leaves_external_aliases_and_frozen_captures_to_owners(
        monkeypatch, tmp_path):
    source = tmp_path/'source.py'
    source.write_text('value = 1\n', encoding='utf-8')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    receipt = tmp_path/'results.json'
    receipt.write_text(json.dumps({'sources': {'source.py': digest,
        'external/raw.csv': 'external custody', 'https://example.org/data': 'external',
        '../outside.py': 'not repository-relative'},
        'external_interface_hashes': {'unknown_alias.py': 'owning verifier'}}), encoding='utf-8')
    (tmp_path/'frozen_capture.json').write_text(
        json.dumps({'source_sha256': {'source.py': 'capture-time hash'}}), encoding='utf-8')
    monkeypatch.setattr(verify, 'ROOT', tmp_path)
    monkeypatch.setattr(verify, 'REPLAYS', (('model.example', None, 'results.json'),))
    verify.verify_receipt_source_bindings()


def test_ci_covers_all_three_shards_both_platforms_and_complete_regression():
    workflow = yaml.safe_load((verify.ROOT/'.github/workflows/observer-model.yml').read_text())
    observer, regression = (workflow['jobs'][key] for key in ('observer', 'regression'))
    assert observer['strategy']['matrix']['os'] == ['ubuntu-latest', 'windows-latest']
    assert observer['strategy']['matrix']['shard'] == [0,1,2]
    assert any(step.get('run') == 'python -m model.verify --shard-index ${{ matrix.shard }} --shard-count 3'
               for step in observer['steps'])
    scalar = [step for step in observer['steps'] if 'verify_source_scalar_time_refinement.py' in step.get('run','')]
    assert len(scalar) == 1 and scalar[0]['if'] == 'matrix.shard == 1'
    assert regression['strategy']['matrix']['os'] == ['ubuntu-latest', 'windows-latest']
    assert any(step.get('run') == 'python -m pytest -q -W error model' for step in regression['steps'])
