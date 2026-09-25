"""Compiler-identity mutations guard against platform false failures and wrong compilers."""
from __future__ import annotations

from copy import deepcopy
import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location('canonical_lean_verifier', Path(__file__).with_name('verify.py'))
assert SPEC is not None and SPEC.loader is not None
VERIFIER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFIER)

MAC = 'Lean (version 4.29.1, arm64-apple-darwin24.6.0, commit f72c35b3f637c8c6571d353742168ab66cc22c00, Release)'
LINUX = MAC.replace('arm64-apple-darwin24.6.0', 'x86_64-unknown-linux-gnu')


def fixture(report: str) -> dict:
    return {'lean_version': report, 'compiler_identity': VERIFIER.compiler_identity(report),
            'sources': [{'source': 'Lean/Example.lean', 'sha256': '1234'}]}


def test_same_compiler_on_another_architecture_preserves_original_provenance() -> None:
    expected, actual = fixture(MAC), fixture(LINUX)
    before = deepcopy((expected, actual))
    VERIFIER.check_receipt(expected, actual)
    assert (expected, actual) == before


@pytest.mark.parametrize('report', [LINUX.replace('4.29.1', '4.29.2'),
                                   LINUX.replace('f72c35b3', 'e72c35b3')])
def test_other_version_or_source_commit_is_rejected(report: str) -> None:
    with pytest.raises(ValueError, match='version or source commit'):
        VERIFIER.check_receipt(fixture(MAC), fixture(report))


def test_forged_identity_cannot_mask_a_wrong_commit() -> None:
    actual = fixture(LINUX.replace('f72c35b3', 'e72c35b3'))
    actual['compiler_identity'] = VERIFIER.compiler_identity(MAC)
    with pytest.raises(ValueError, match='disagrees with its provenance'):
        VERIFIER.check_receipt(fixture(MAC), actual)


def test_architecture_exception_cannot_hide_a_changed_proof_pin() -> None:
    actual = fixture(LINUX)
    actual['sources'][0]['sha256'] = 'changed'
    with pytest.raises(ValueError, match='differs from fresh verification'):
        VERIFIER.check_receipt(fixture(MAC), actual)


def test_malformed_compiler_report_is_rejected() -> None:
    with pytest.raises(ValueError, match='Unrecognized'):
        VERIFIER.compiler_identity('Lean 4.29.1 with no source commit')
