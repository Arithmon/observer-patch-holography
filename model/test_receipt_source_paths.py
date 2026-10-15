"""Exercise actual receipt producers with Windows path semantics on every OS."""
import hashlib
from pathlib import Path, PureWindowsPath

from model.empirical_diagnostics import experiment as diagnostic
from model.primitive_observer import experiment as primitive
from model.primitive_observer import action, completion_report, observer_report, selection

REAL_ROOT=Path(__file__).resolve().parents[1]
WINDOWS_ROOT=PureWindowsPath('C:/portable-repository')


class WindowsSourcePath(PureWindowsPath):
    def real(self):
        return REAL_ROOT/self.relative_to(WINDOWS_ROOT).as_posix()

    def read_bytes(self):
        return self.real().read_bytes()

    def glob(self,pattern):
        return (WindowsSourcePath(WINDOWS_ROOT/p.relative_to(REAL_ROOT).as_posix())
                for p in self.real().glob(pattern))


def assert_portable_bindings(bindings):
    assert bindings and all('\\' not in name for name in bindings)
    for name,digest in bindings.items():
        assert hashlib.sha256((REAL_ROOT/name).read_bytes()).hexdigest()==digest


def test_diagnostic_producer_uses_portable_names_with_windows_paths(monkeypatch):
    monkeypatch.setattr(diagnostic,'ROOT',WindowsSourcePath(WINDOWS_ROOT))
    monkeypatch.setattr(diagnostic,'FOLDER',WindowsSourcePath(WINDOWS_ROOT/'model/empirical_diagnostics'))
    assert_portable_bindings(diagnostic._source_bindings())


def test_primitive_producer_uses_portable_names_with_windows_glob(monkeypatch):
    # Only the expensive scientific calculations are stubbed. The actual
    # source discovery, naming and byte hashing execute unchanged.
    monkeypatch.setattr(primitive,'ROOT',WindowsSourcePath(WINDOWS_ROOT))
    monkeypatch.setattr(primitive,'FOLDER',WindowsSourcePath(WINDOWS_ROOT/'model/primitive_observer'))
    for module in (action,completion_report,observer_report,selection):
        monkeypatch.setattr(module,'run',lambda:{'all_pass':True})
    monkeypatch.setattr(primitive,'predict_laboratory',lambda **kw:{})
    monkeypatch.setattr(primitive,'evaluate_laboratory',lambda prediction:{
        'scores':{'metrics':{'combined':{'observations':125}}},
        'independent_existing_executor_comparison':{'maximum_joint_probability_residual':0.},
        'prediction':{'positive_time_outcomes_used':False}})
    bindings=primitive.run()['source_bindings']
    assert_portable_bindings(bindings)
    assert 'model/primitive_observer/transfer.py' in bindings
