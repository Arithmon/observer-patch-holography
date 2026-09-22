#!/usr/bin/env python3
"""Lock the compiler-trusted ``native_decide`` inventory in the Lean tree.

The proofs in this inventory are admission-free, but Lean reports generated
native-code evaluation axioms for them. Any change must be reviewed together
with the public trust-boundary wording.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LEAN_ROOT = REPO_ROOT / "Lean"

EXPECTED = Counter(
    {
        "Screen/PortGramA5Isometry.lean": 6,
        "Screen/PortGramRepairCovariance.lean": 7,
    }
)

NATIVE_TOKEN_RE = re.compile(r"\bnative_decide\b")


def code_view(text: str) -> str:
    """Lean source with block comments, line comments and string literals blanked.

    Block comments nest. Blanking keeps every other byte in place, so a match
    inside a tactic block, a term-mode ``(by native_decide)`` or a one-line
    ``exact`` counts, while prose that merely names the tactic does not.
    """

    out: list[str] = []
    i = 0
    depth = 0
    n = len(text)
    while i < n:
        two = text[i : i + 2]
        if depth:
            if two == "/-":
                depth += 1
                i += 2
            elif two == "-/":
                depth -= 1
                i += 2
            else:
                i += 1
            continue
        if two == "/-":
            depth = 1
            i += 2
            continue
        if two == "--":
            end = text.find("\n", i)
            i = n if end < 0 else end
            continue
        if text[i] == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            i = min(j + 1, n)
            continue
        out.append(text[i])
        i += 1
    return "".join(out)


def inventory() -> Counter[str]:
    found: Counter[str] = Counter()
    for path in sorted(LEAN_ROOT.rglob("*.lean")):
        if any(part.startswith(".") for part in path.relative_to(LEAN_ROOT).parts):
            continue
        count = len(NATIVE_TOKEN_RE.findall(code_view(path.read_text(encoding="utf-8"))))
        if count:
            found[str(path.relative_to(LEAN_ROOT))] = count
    return found


def main() -> int:
    found = inventory()
    total = sum(found.values())
    print(f"Lean native_decide proof inventory: {total}")
    for path, count in sorted(found.items()):
        print(f"  {path}: {count}")
    if found != EXPECTED:
        print("FAIL: native_decide trust inventory changed", file=sys.stderr)
        print(f"expected: {dict(sorted(EXPECTED.items()))}", file=sys.stderr)
        print(f"found:    {dict(sorted(found.items()))}", file=sys.stderr)
        return 1
    print("OK: compiler-trusted proof inventory matches the reviewed set")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
