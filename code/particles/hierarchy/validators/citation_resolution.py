#!/usr/bin/env python3
"""Resolve the path-like citations carried by a hierarchy certificate payload.

Each certificate in this bundle names its premises by repository path. A
citation that names a file the repository does not contain is a dangling
pointer: the published certificate then quotes a document no reader can open.

The scan below walks a whole certificate payload, collects every path-like
token in every string, and resolves each token against three roots: the
repository root, the bundle root, and the bundle ``certificates`` directory.
A token that resolves against none of them is reported, except in a string
that declares an external result (a theorem imported from outside this
repository, for example a Mathlib lemma), because such a citation has no
repository path to resolve. Tokens that do resolve are checked in every
string, external declaration or not.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Iterator

CITED_EXTENSIONS = (
    "tex",
    "json",
    "py",
    "lean",
    "md",
    "yaml",
    "yml",
    "csv",
    "txt",
    "sh",
    "pdf",
)

CITATION_PATTERN = re.compile(
    r"(?<![\w./-])"
    r"[A-Za-z0-9_][A-Za-z0-9_./-]*"
    r"\.(?:" + "|".join(CITED_EXTENSIONS) + r")"
    r"(?![\w/])"
)

EXTERNAL_RESULT_MARKERS = (
    "external result",
    "external theorem",
    "imported result",
    "external mathematical result",
)

ROOT_MARKERS = ("paper", "code", "Lean")


def declares_external_result(text: str) -> bool:
    """Report whether a citation string marks its source as external."""

    lowered = text.lower()
    return any(marker in lowered for marker in EXTERNAL_RESULT_MARKERS)


def iter_strings(node: Any, trail: str = "") -> Iterator[tuple[str, str]]:
    """Yield ``(location, text)`` for every string in a JSON-like payload."""

    if isinstance(node, dict):
        for key, value in node.items():
            child = f"{trail}.{key}" if trail else str(key)
            yield from iter_strings(value, child)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from iter_strings(value, f"{trail}[{index}]")
    elif isinstance(node, str):
        yield trail, node


def repository_root(start: Path) -> Path:
    """Walk up from ``start`` to the directory carrying the repository roots."""

    start = start.resolve()
    for candidate in (start, *start.parents):
        if all((candidate / marker).is_dir() for marker in ROOT_MARKERS):
            return candidate
    return start


def citation_tokens(text: str) -> list[str]:
    """Return the path-like tokens of one string, anchors stripped."""

    return [match.group(0) for match in CITATION_PATTERN.finditer(text)]


def unresolved_citations(
    payload: Any,
    *,
    repo_root: Path,
    bundle_root: Path,
) -> list[dict[str, str]]:
    """Return every citation token that resolves against none of the roots."""

    roots = (repo_root, bundle_root, bundle_root / "certificates")
    findings: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for location, text in iter_strings(payload):
        external = declares_external_result(text)
        for token in citation_tokens(text):
            relative = token[2:] if token.startswith("./") else token
            if any((root / relative).exists() for root in roots):
                continue
            if external:
                continue
            key = (location, token)
            if key in seen:
                continue
            seen.add(key)
            findings.append(
                {
                    "citation": token,
                    "location": location,
                    "quoted_string": text if len(text) <= 200 else text[:197] + "...",
                }
            )
    return findings


def scan_certificate(payload: Any, validator_file: str) -> list[dict[str, str]]:
    """Scan one certificate payload from a validator inside ``validators/``."""

    bundle_root = Path(validator_file).resolve().parents[1]
    return unresolved_citations(
        payload,
        repo_root=repository_root(bundle_root),
        bundle_root=bundle_root,
    )
