#!/usr/bin/env python3
"""Fail closed when committed scientific receipts leak developer-home paths."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECEIPT_ROOTS = (
    ROOT / "code" / "particles" / "runs",
    ROOT / "claims",
    ROOT / "evidence",
)
# Attested bytes that carry a developer-home path and cannot be rewritten: the
# anchored FZ-01 registry records where its pinned tracker section lived on the
# author's machine at freeze time, and its sha256 is fixed by the 2026-07-17
# OpenTimestamps attestation. The pinned bytes are vendored beside it as
# tracker_section8_pinned_text.md, and PROVENANCE.md in that directory records
# the repository reference. Each entry is (repo-relative path, JSON pointer).
ATTESTED_PATH_ALLOWLIST: frozenset[tuple[str, str]] = frozenset(
    {
        (
            "evidence/custody/falsification/frozen_targets/fz01_2026-07-17/"
            "fz01_freeze_registry_2026-07-17.json",
            "/pinned_source/file",
        ),
    }
)
# Original remote-machine provenance is immutable, while the canonical replay
# uses the vendored relative paths. Permit only these two historical fields,
# in original files whose bytes and containing manifest match the reviewed
# archive. A changed file/manifest, new field, or new file still fails closed.
OBSERVER_ARCHIVE = "evidence/observer_dynamics_20260925/"
OBSERVER_MANIFEST_SHA256 = "c01616121d3cba371633f396b5b78c67c48a6dfa1e487004ea7f3ef165375915"


def _original_observer_provenance(rel: str, pointer: str, data: bytes) -> bool:
    if not rel.startswith(OBSERVER_ARCHIVE):
        return False
    local = rel[len(OBSERVER_ARCHIVE):]
    chain = pointer == "/parent_run" and re.fullmatch(
        r"sim-analysis/data/refine_ensemble/chain_\d{3}_L[78]_s\d+_f\d+_(?:chain|shuffled)\.json", local
    )
    geometry = pointer == "/geometry/cache_dir" and local in {
        f"sim-analysis/data/codex_audit_20260925/inputs/data/receipts/L{level}/receipt.json"
        for level in (8, 9, 10)
    }
    if not (chain or geometry):
        return False
    manifest = ROOT / OBSERVER_ARCHIVE / "manifest.json"
    if not manifest.is_file():
        return False
    original = manifest.read_bytes()
    if hashlib.sha256(original).hexdigest() != OBSERVER_MANIFEST_SHA256:
        return False
    entry = json.loads(original)["files"].get(local, {})
    return entry.get("bytes") == len(data) and entry.get("sha256") == hashlib.sha256(data).hexdigest()
DEVELOPER_HOME_PATTERNS = (
    re.compile(r"/Users/[^/]+/"),
    re.compile(r"/home/[^/]+/"),
    re.compile(r"[A-Za-z]:[\\/]Users[\\/][^\\/]+[\\/]", re.IGNORECASE),
)


@dataclass(frozen=True)
class PortabilityViolation:
    path: Path
    json_pointer: str
    value: str


def _escape_pointer(token: object) -> str:
    return str(token).replace("~", "~0").replace("/", "~1")


def _is_developer_home_path(value: str) -> bool:
    return any(pattern.search(value) for pattern in DEVELOPER_HOME_PATTERNS)


def _walk(value: Any, pointer: str = "") -> Iterable[tuple[str, str]]:
    if isinstance(value, dict):
        for key, item in value.items():
            yield from _walk(item, f"{pointer}/{_escape_pointer(key)}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from _walk(item, f"{pointer}/{index}")
    elif isinstance(value, str):
        yield pointer or "/", value


def receipt_paths(roots: Iterable[Path] = DEFAULT_RECEIPT_ROOTS) -> list[Path]:
    paths: set[Path] = set()
    for root in roots:
        if root.is_file() and root.suffix == ".json":
            paths.add(root)
        elif root.exists():
            paths.update(root.rglob("*.json"))
    return sorted(paths)


def find_violations(paths: Iterable[Path]) -> list[PortabilityViolation]:
    violations: list[PortabilityViolation] = []
    for path in paths:
        data = path.read_bytes()
        try:
            payload = json.loads(data.decode("utf-8"))
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid JSON receipt {path}: {exc}") from exc
        try:
            rel = path.resolve().relative_to(ROOT).as_posix()
        except ValueError:
            rel = path.as_posix()
        for pointer, value in _walk(payload):
            if (rel, pointer) in ATTESTED_PATH_ALLOWLIST:
                continue
            if _is_developer_home_path(value):
                if _original_observer_provenance(rel, pointer, data):
                    continue
                violations.append(PortabilityViolation(path, pointer, value))
    return violations


def check(paths: Iterable[Path]) -> None:
    violations = find_violations(paths)
    if not violations:
        return

    def display_path(path: Path) -> str:
        try:
            return path.relative_to(ROOT).as_posix()
        except ValueError:
            return path.as_posix()

    rendered = "\n".join(
        f"{display_path(item.path)}{item.json_pointer}: {item.value}"
        for item in violations
    )
    raise ValueError(
        "developer-home path leaked into committed scientific receipts:\n"
        f"{rendered}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        help="optional receipt files/directories (defaults to particle runs, claims and evidence)",
    )
    args = parser.parse_args()
    roots = tuple(path.resolve() for path in args.paths) or DEFAULT_RECEIPT_ROOTS
    paths = receipt_paths(roots)
    check(paths)
    print(f"receipt portability OK ({len(paths)} JSON artifacts)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
