"""Mutation tests for archive identity, simulator pinning and receipt comparisons."""
import hashlib
import importlib.machinery
import json
import os
from pathlib import Path, PureWindowsPath
import subprocess
import sys
import pytest
from .archive_adapter import COMMIT, VENDOR, SIM, check_manifest, frozen_sources, load
from .verify import compare, subprocess_pass_count


def test_changed_evidence_fails_before_scientific_replay(tmp_path):
    (tmp_path / "input").write_bytes(b"valid")
    data = {"files": {"input": {"bytes": 5, "sha256": hashlib.sha256(b"valid").hexdigest()}}}
    (tmp_path / "manifest.json").write_text(json.dumps(data))
    assert check_manifest(tmp_path)["files"] == 1
    (tmp_path / "input").write_bytes(b"wrong")
    with pytest.raises(ValueError, match="digest mismatch"):
        check_manifest(tmp_path)


def test_archive_cannot_pin_outside_its_directory(tmp_path):
    (tmp_path / "manifest.json").write_text(json.dumps({"files": {"../outside": {"bytes": 0, "sha256": ""}}}))
    with pytest.raises(ValueError, match="escapes"):
        check_manifest(tmp_path)


def test_pinned_source_reads_exact_blob_and_rejects_wrong_revision():
    with frozen_sources():
        actual = subprocess.check_output(["git", "-C", str(VENDOR), "show", f"{COMMIT}:oph_exact/carrier.py"])
        assert actual == (VENDOR / "oph_exact/carrier.py").read_bytes()
        with pytest.raises(ValueError, match="revision"):
            subprocess.check_output(["git", "-C", str(VENDOR), "show", "0000000:oph_exact/carrier.py"])
        with pytest.raises(ValueError, match="Unarchived"):
            subprocess.check_output(["git", "-C", str(VENDOR), "show", f"{COMMIT}:missing.py"])


def test_receipt_comparison_preserves_exact_claims_and_rejects_large_errors():
    compare({"matrix": ["1/3", 2], "numeric": 1.0}, {"matrix": ["1/3", 2], "numeric": 1.0+1e-9})
    with pytest.raises(AssertionError): compare({"matrix": ["1/3"]}, {"matrix": ["1/2"]})
    with pytest.raises(AssertionError): compare({"numeric": 1.0}, {"numeric": 1.1})
    with pytest.raises(AssertionError): compare({"a": 1}, {"a": 1, "b": 2})


def test_crosspins_reject_wrong_terminal_state():
    from .archive_adapter import SIM
    from .retained_measurements import scale_ensemble
    value=json.loads((SIM/'data/codex_audit_20260925/scale_ensemble.json').read_text())
    assert scale_ensemble(value)['terminal_digest_crosspins']==48
    value['levels'][0]['inputs'][0]['array_sha256']='0'*64
    with pytest.raises(ValueError,match='crosspin'):
        scale_ensemble(value)


def test_schedule_power_identity_rejects_changed_summary():
    from .archive_adapter import SIM
    from .retained_measurements import scale_ensemble
    value=json.loads((SIM/'data/codex_audit_20260925/scale_ensemble.json').read_text())
    value['levels'][0]['scales'][0]['cross_schedule_common_power_unbiased']+=.1
    with pytest.raises(ValueError,match='arithmetic'):
        scale_ensemble(value)



def test_coherently_rescaled_power_summary_cannot_change_retained_absolute_amplitude():
    from .archive_adapter import SIM
    from .retained_measurements import scale_ensemble
    value=json.loads((SIM/'data/codex_audit_20260925/scale_ensemble.json').read_text())
    scale=value['levels'][0]['scales'][0]
    for name in ('total_terminal_variance', 'conditional_schedule_variance_unbiased',
                 'cross_schedule_common_power_unbiased'):
        scale[name] *= 2
    # Both aggregate identities and the schedule-noise fraction still hold;
    # the independent schedule variances must expose the amplitude change.
    with pytest.raises(ValueError,match='arithmetic'):
        scale_ensemble(value)


def test_changed_exchangeable_normalization_is_rejected():
    from .archive_adapter import SIM
    from .retained_measurements import scale_ensemble
    value=json.loads((SIM/'data/codex_audit_20260925/scale_ensemble.json').read_text())
    value['levels'][0]['scales'][0]['per_schedule'][0]['variance_over_exchangeable_equilibrium'] *= 2
    with pytest.raises(ValueError,match='arithmetic'):
        scale_ensemble(value)


class WindowsRelativePaths(type(Path())):
    """Keep native file access, but serialize relative paths as Windows does."""
    def relative_to(self, *args, **kwargs):
        return PureWindowsPath(super().relative_to(*args, **kwargs).as_posix())


def test_windows_spectrum_paths_preserve_strict_receipt_comparison(monkeypatch):
    with frozen_sources():
        module = load('codex/observer_cmb/observer_spectrum.py')
        for key in ('SPECTRUM', 'RECEIPT', 'INI'):
            monkeypatch.setattr(module, key, WindowsRelativePaths(getattr(module, key)))
        raw = module.make_output.__wrapped__()
        expected = json.loads((SIM/'codex/observer_cmb/observer_spectrum.json').read_text())
        with pytest.raises(ValueError, match='value mismatch'):
            module._compare(raw['inputs_relative_to_codex'], expected['inputs_relative_to_codex'])
        adapted = module.verify_output()
        adapted['inputs_relative_to_codex']['spectrum'] = 'wrong/spectrum.txt'
        with pytest.raises(ValueError, match='value mismatch'):
            module._compare(adapted['inputs_relative_to_codex'], expected['inputs_relative_to_codex'])


def test_windows_measurement_paths_match_immutable_input_pins(monkeypatch):
    with frozen_sources():
        module = load('codex/observer_cmb/measured_comparison.py')
        monkeypatch.setattr(module, 'CODEX', WindowsRelativePaths(module.CODEX))
        expected = json.loads((SIM/'codex/observer_cmb/measured_comparison.json').read_text())['measurement_inputs']
        assert module.load_measurements.__wrapped__()[1] != expected
        assert module.load_measurements()[1] == expected


def test_windows_manual_module_loader_preserves_planck_audit(monkeypatch):
    with frozen_sources():
        module = load('codex/boltzmann_data/audit_inputs.py')
        for key in ('HERE', 'OBSERVATION'):
            monkeypatch.setattr(module, key, WindowsRelativePaths(getattr(module, key)))
        expected = json.loads((SIM/'codex/boltzmann_data/input_audit.json').read_text())
        assert module.verify_downloads.__wrapped__() != expected['inputs']
        assert module.make_receipt() == expected


def test_archive_context_restores_loader_and_subprocess_hooks():
    loader = importlib.machinery.SourceFileLoader.exec_module
    command = subprocess.check_output
    with frozen_sources():
        assert importlib.machinery.SourceFileLoader.exec_module is not loader
        assert subprocess.check_output is not command
    assert importlib.machinery.SourceFileLoader.exec_module is loader
    assert subprocess.check_output is command


def test_windows_forward_summary_retains_original_failure_log(monkeypatch):
    with frozen_sources():
        module = load('codex/boltzmann/summarize.py')
        monkeypatch.setattr(module, 'HERE', WindowsRelativePaths(module.HERE))
        expected = json.loads((SIM/'codex/boltzmann/summary.json').read_text())
        assert module.build.__wrapped__()['observational_input_hashes'] != expected['observational_input_hashes']
        compare(expected, module.build())


@pytest.mark.parametrize('skip_body', [
    'import pytest\npytest.skip("missing control", allow_module_level=True)\n',
    'import pytest\ndef test_control():\n    pytest.skip("missing control")\n',
])
def test_real_pytest_collection_and_runtime_skips_fail_closed(tmp_path, skip_body):
    (tmp_path/'test_present.py').write_text('def test_present():\n    assert True\n')
    (tmp_path/'test_missing.py').write_text(skip_body)
    script = '''import json, runpy, sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(sys.argv[1]).parent))
verifier = runpy.run_path(sys.argv[1])
counts = verifier["Counts"]()
code = pytest.main(["-q", "-p", "no:cacheprovider", sys.argv[2]], plugins=[counts])
print("CONTROL_COUNTS=" + json.dumps(vars(counts)))
assert code == 0, "The false-green fixture itself must return pytest success"
verifier["require_complete_controls"](counts.skipped)
'''
    result = subprocess.run([sys.executable, '-c', script,
                             str(Path(__file__).with_name('verify.py')), str(tmp_path)],
                            capture_output=True, text=True,
                            env={**os.environ, 'PYTEST_DISABLE_PLUGIN_AUTOLOAD': '1'})
    assert result.returncode != 0
    counts = json.loads(next(line.split('=', 1)[1] for line in result.stdout.splitlines()
                             if line.startswith('CONTROL_COUNTS=')))
    assert counts == {'passed': 1, 'failed': 0, 'skipped': 1}
    assert 'Mandatory archived controls were skipped' in result.stderr
    with pytest.raises(RuntimeError, match='controls were skipped'):
        subprocess_pass_count(result.stdout)
