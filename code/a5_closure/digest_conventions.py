#!/usr/bin/env python3
"""Tag-agnostic resolver for the three sha256 conventions used by this package.

Three distinct serializations produce the digests pinned across
``code/a5_closure/manifests/`` and ``code/a5_closure/receipts/``:

``raw_bytes``
    sha256 over the file bytes exactly as committed, including the trailing
    newline.  This is ``baryon_dimension_six_census.sha256_file``.

``canonical_json``
    sha256 over ``json.dumps(value, sort_keys=True, separators=(",", ":"),
    ensure_ascii=False)`` encoded as UTF-8, with no trailing newline.  This is
    ``echosahedral_selector_certificate.sha256_json``, reused by
    ``axis_center_descent_certificate`` and the other producers that import it.

``canonical_json_self_digest_stripped``
    the same canonical serialization applied to an artifact body whose own
    self-digest field has been removed, which is how a semantic artifact signs
    itself and how a downstream manifest pins that signature.

The ``sha256:`` string prefix names none of the three.  It marks raw-byte
digests in ``receipts/baryon_dimension_six_census.receipt.json`` and
self-digest-stripped digests in ``manifests/axis_center_descent_reference.json``
alike, and the canonical-JSON pins in
``manifests/common_ew_order_unit_carrier_reference.json`` carry no prefix at
all.  The resolver therefore ignores the tag and the field name, tries all
three serializations, and reports which ones reproduce a given pin.

Run as a script to classify every path-paired pin in the package:

    python3 digest_conventions.py
    python3 digest_conventions.py --json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

MODULE_DIR = Path(__file__).resolve().parent
REPO_ROOT = MODULE_DIR.parents[1]

PIN_TREES = ("manifests", "receipts")

RAW_BYTES = "raw_bytes"
CANONICAL_JSON = "canonical_json"
CANONICAL_JSON_SELF_DIGEST_STRIPPED = "canonical_json_self_digest_stripped"

CONVENTIONS = (RAW_BYTES, CANONICAL_JSON, CANONICAL_JSON_SELF_DIGEST_STRIPPED)

DIGEST_TAG = "sha256:"
DIGEST_RE = re.compile(r"^(?:sha256:)?[0-9a-f]{64}$")

# Suffixes that mark a field as carrying a digest, so that ``<stem>_sha256``
# can be paired with a sibling ``<stem>_path``.
DIGEST_FIELD_SUFFIXES = ("_sha256", "_hash", "_digest")

# Suffixes that mark a string as naming a file, so that a mapping key such as
# "local_domain/stage1_event_complex.py" is read as a path.
PATH_SUFFIXES = frozenset(
    {".json", ".py", ".lean", ".tex", ".md", ".txt", ".csv", ".npz", ".gz", ".yaml", ".yml"}
)


class DigestConventionError(RuntimeError):
    """Raised when a digest convention cannot be applied to a file."""


def normalize_digest(digest: str) -> str:
    """Strip the uninformative ``sha256:`` tag and lowercase the hex."""
    text = digest.strip()
    if text.startswith(DIGEST_TAG):
        text = text[len(DIGEST_TAG) :]
    return text.lower()


def looks_like_digest(value: Any) -> bool:
    return isinstance(value, str) and bool(DIGEST_RE.match(value.strip()))


def looks_like_path(value: Any) -> bool:
    if not isinstance(value, str) or not value or " " in value:
        return False
    if looks_like_digest(value):
        return False
    return "/" in value or Path(value).suffix in PATH_SUFFIXES


def canonical_json_bytes(value: Any) -> bytes:
    """Compact, key-sorted UTF-8 JSON with no trailing newline."""
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")


def sha256_hex(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def digest_raw_bytes(path: Path) -> str:
    """Convention 1: sha256 over the committed file bytes."""
    try:
        return sha256_hex(path.read_bytes())
    except OSError as exc:  # pragma: no cover - surfaced through the CLI
        raise DigestConventionError(f"cannot read {path}: {exc}") from exc


def load_json_document(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise DigestConventionError(f"cannot parse {path} as JSON: {exc}") from exc


def digest_canonical_json(path: Path) -> str:
    """Convention 2: sha256 over the canonical JSON serialization."""
    return sha256_hex(canonical_json_bytes(load_json_document(path)))


def self_digest_fields(document: Any) -> tuple[str, ...]:
    """Top-level fields whose value signs the rest of the document.

    A field qualifies when removing it and canonically serializing what is left
    reproduces that field's own value.  The test is on the bytes, so it names
    the self-digest field without trusting its name.
    """
    if not isinstance(document, Mapping):
        return ()
    found: list[str] = []
    for field, value in document.items():
        if not looks_like_digest(value):
            continue
        body = {key: item for key, item in document.items() if key != field}
        if sha256_hex(canonical_json_bytes(body)) == normalize_digest(value):
            found.append(field)
    return tuple(found)


def digest_canonical_json_self_digest_stripped(path: Path) -> dict[str, str]:
    """Convention 3: canonical JSON of the body without its self-digest field.

    Returns a mapping from each self-digest field to the digest of the body that
    field signs.  An artifact that signs itself has exactly one such field; a
    document that signs nothing returns an empty mapping.
    """
    document = load_json_document(path)
    digests: dict[str, str] = {}
    for field in self_digest_fields(document):
        body = {key: item for key, item in document.items() if key != field}
        digests[field] = sha256_hex(canonical_json_bytes(body))
    return digests


def convention_digests(path: Path) -> dict[str, tuple[str, ...]]:
    """Every digest each convention produces for one file.

    ``raw_bytes`` always yields one value.  The two JSON conventions yield
    nothing for a file that is not JSON, and the stripped convention yields
    nothing for a JSON document that carries no self-digest field.
    """
    digests: dict[str, tuple[str, ...]] = {RAW_BYTES: (digest_raw_bytes(path),)}
    try:
        document = load_json_document(path)
    except DigestConventionError:
        digests[CANONICAL_JSON] = ()
        digests[CANONICAL_JSON_SELF_DIGEST_STRIPPED] = ()
        return digests
    digests[CANONICAL_JSON] = (sha256_hex(canonical_json_bytes(document)),)
    stripped = digest_canonical_json_self_digest_stripped(path)
    digests[CANONICAL_JSON_SELF_DIGEST_STRIPPED] = tuple(sorted(set(stripped.values())))
    return digests


def resolve_convention(digest: str, path: Path) -> tuple[str, ...]:
    """Which of the three conventions reproduce ``digest`` for ``path``.

    The field name and the ``sha256:`` tag are ignored.  An empty result means
    no convention in this package produces that digest; more than one result
    means the pin cannot be read back to a single serialization.
    """
    target = normalize_digest(digest)
    produced = convention_digests(path)
    return tuple(name for name in CONVENTIONS if target in produced[name])


@dataclass(frozen=True)
class PathPairedPin:
    """A digest field paired, inside one JSON object, with a file path."""

    source: Path
    pointer: str
    field: str
    digest: str
    declared_path: str

    @property
    def location(self) -> str:
        return f"{self.source.name}{self.pointer}"


def _pair_path(container: Mapping[str, Any], field: str) -> str | None:
    """The path field that a digest field in the same object refers to."""
    if looks_like_path(field):
        # Mapping form: the key is the path and the value is its digest.
        return field
    for suffix in DIGEST_FIELD_SUFFIXES:
        if not field.endswith(suffix):
            continue
        stem = field[: -len(suffix)]
        for candidate in (f"{stem}_path", f"{stem}path", stem):
            if looks_like_path(container.get(candidate)):
                return str(container[candidate])
    if looks_like_path(container.get("path")):
        return str(container["path"])
    return None


def iter_path_paired_pins(source: Path) -> Iterator[PathPairedPin]:
    """Every digest in ``source`` that names the file it pins."""

    def walk(node: Any, pointer: str) -> Iterator[PathPairedPin]:
        if isinstance(node, Mapping):
            for field, value in node.items():
                if not looks_like_digest(value):
                    continue
                declared = _pair_path(node, field)
                if declared is None:
                    continue
                yield PathPairedPin(
                    source=source,
                    pointer=f"{pointer}/{field}",
                    field=field,
                    digest=str(value),
                    declared_path=declared,
                )
            for field, value in node.items():
                yield from walk(value, f"{pointer}/{field}")
        elif isinstance(node, Sequence) and not isinstance(node, (str, bytes)):
            for index, value in enumerate(node):
                yield from walk(value, f"{pointer}/{index}")

    yield from walk(load_json_document(source), "")


def locate(declared_path: str, roots: Sequence[Path]) -> Path | None:
    """Resolve a declared pin path against the package's two path bases.

    Pins are written either repo-root-relative (``code/a5_closure/...``) or
    package-relative (``receipts/...``).  A declared path that resolves under
    neither names a file outside this tree, such as a module in an upstream
    producer checkout; guessing a third base would let a basename collision
    classify the wrong file.
    """
    for root in roots:
        candidate = root / declared_path
        if not candidate.is_file():
            continue
        resolved = candidate.resolve()
        if resolved.is_relative_to(REPO_ROOT):
            return resolved
    return None


@dataclass(frozen=True)
class Classification:
    pin: PathPairedPin
    resolved: Path | None
    conventions: tuple[str, ...]

    @property
    def located(self) -> bool:
        return self.resolved is not None


def classify_pins(
    trees: Sequence[Path] | None = None,
    roots: Sequence[Path] | None = None,
) -> list[Classification]:
    """Classify every path-paired pin under the given pin trees."""
    if trees is None:
        trees = [MODULE_DIR / name for name in PIN_TREES]
    if roots is None:
        roots = (REPO_ROOT, MODULE_DIR)
    results: list[Classification] = []
    for tree in trees:
        for source in sorted(Path(tree).rglob("*.json")):
            for pin in iter_path_paired_pins(source):
                resolved = locate(pin.declared_path, roots)
                conventions = (
                    resolve_convention(pin.digest, resolved) if resolved else ()
                )
                results.append(Classification(pin, resolved, conventions))
    return results


def summarize(results: Sequence[Classification]) -> dict[str, Any]:
    located = [row for row in results if row.located]
    split = {name: 0 for name in CONVENTIONS}
    unreproduced: list[Classification] = []
    ambiguous: list[Classification] = []
    for row in located:
        if not row.conventions:
            unreproduced.append(row)
        elif len(row.conventions) > 1:
            ambiguous.append(row)
        else:
            split[row.conventions[0]] += 1
    outside = [row for row in results if not row.located]
    return {
        "pins": len(located),
        "unreproduced": unreproduced,
        "ambiguous": ambiguous,
        "split": split,
        "outside_tree": outside,
    }


def summary_line(summary: Mapping[str, Any]) -> str:
    return (
        f"{summary['pins']} path-paired pins, "
        f"{len(summary['unreproduced'])} unreproduced, "
        f"{len(summary['ambiguous'])} ambiguous"
    )


def report(results: Sequence[Classification]) -> tuple[str, int]:
    summary = summarize(results)
    lines = [summary_line(summary)]
    width = max(len(name) for name in CONVENTIONS)
    for name in CONVENTIONS:
        lines.append(f"  {name:<{width}}  {summary['split'][name]}")
    outside = summary["outside_tree"]
    if outside:
        paths = sorted({row.pin.declared_path for row in outside})
        lines.append(
            f"  {len(outside)} pin sites name {len(paths)} paths outside this tree; "
            "not classified"
        )
    for row in summary["unreproduced"]:
        lines.append(
            f"  UNREPRODUCED {row.pin.location} -> {row.pin.declared_path} "
            f"({row.pin.digest})"
        )
    for row in summary["ambiguous"]:
        lines.append(
            f"  AMBIGUOUS {row.pin.location} -> {row.pin.declared_path} "
            f"({', '.join(row.conventions)})"
        )
    status = 1 if summary["unreproduced"] or summary["ambiguous"] else 0
    return "\n".join(lines), status


def as_json(results: Sequence[Classification]) -> str:
    summary = summarize(results)
    payload = {
        "pins": summary["pins"],
        "split": summary["split"],
        "unreproduced": [row.pin.location for row in summary["unreproduced"]],
        "ambiguous": [row.pin.location for row in summary["ambiguous"]],
        "outside_tree": sorted({row.pin.declared_path for row in summary["outside_tree"]}),
        "classified": [
            {
                "source": row.pin.source.relative_to(REPO_ROOT).as_posix(),
                "pointer": row.pin.pointer,
                "field": row.pin.field,
                "digest": row.pin.digest,
                "path": row.pin.declared_path,
                "conventions": list(row.conventions),
            }
            for row in results
            if row.located
        ],
    }
    return json.dumps(payload, indent=2, sort_keys=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Classify every path-paired sha256 pin in this package by the "
            "serialization that reproduces it, ignoring the field name and the "
            "sha256: tag."
        )
    )
    parser.add_argument(
        "--json", action="store_true", help="emit the full classification as JSON"
    )
    parser.add_argument(
        "path",
        nargs="*",
        type=Path,
        help="files or digests to classify instead of the whole package",
    )
    parser.add_argument(
        "--digest",
        help="one pin value to resolve against each path given",
    )
    args = parser.parse_args(argv)

    if args.digest:
        status = 0
        for path in args.path:
            conventions = resolve_convention(args.digest, path)
            print(f"{path}: {', '.join(conventions) if conventions else 'none'}")
            if len(conventions) != 1:
                status = 1
        return status

    results = classify_pins()
    if args.json:
        print(as_json(results))
        summary = summarize(results)
        return 1 if summary["unreproduced"] or summary["ambiguous"] else 0
    text, status = report(results)
    print(text)
    return status


if __name__ == "__main__":
    sys.exit(main())
