"""Run actual receipt producers with Windows path semantics on every host."""
from importlib import import_module
from pathlib import PureWindowsPath

import pytest


@pytest.mark.parametrize('package', ['bound_pair', 'field_tastes'])
def test_source_bindings_have_the_same_keys_under_windows_paths(package, monkeypatch):
    producer = import_module(f'model.{package}.experiment')
    root = producer.ROOT
    expected = producer.source_bindings()

    class WindowsSourcePath(PureWindowsPath):
        """Windows path operations, with file bytes supplied by this checkout."""
        def read_bytes(self):
            return (root/self.relative_to(windows_root).as_posix()).read_bytes()

        def glob(self, pattern):
            actual = root/self.relative_to(windows_root).as_posix()
            # Deliberately reverse discovery order as well as changing platform.
            return iter(reversed([self/path.name for path in actual.glob(pattern)]))

    windows_root = WindowsSourcePath('C:/receipt-source-checkout')
    monkeypatch.setattr(producer, 'ROOT', windows_root)
    monkeypatch.setattr(producer, 'HERE', windows_root/'model'/package)
    actual = producer.source_bindings()
    assert actual == expected
    assert all('\\' not in name and not PureWindowsPath(name).is_absolute()
               for name in actual)
