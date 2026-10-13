#!/usr/bin/env python3
"""Render the claim matrices and the dependency graph from the claim registry.

``claims/claim_registry.yaml`` is the single source of the claims ledger. Its
``claims`` rows hold each claim's statement, scope and evidence, and its
``claim_context`` section holds, per claim, the prior-work record with its
citations, the owner section, the falsifiers and the claims it depends on. The
rows stay byte-stable because several lane receipts hash complete rows. This
tool renders the three tabular views that tests, lane verifiers and readers
consume:

* ``claims/novelty_matrix.csv``: prior work, overlap, imported step, OPH delta,
  novelty type, owner section, falsifier and the prior-work citation keys
  (entries of ``claims/prior_work.json``) per claim;
* ``claims/falsification_matrix.csv``: mathematical, physical-identification and
  phenomenological falsifiers and the scope if false;
* ``claims/dependency_graph.json``: claim-to-claim edges with their roles.

Rows follow registry order. ``--check`` fails when a committed view differs from
its rendering, so the views cannot drift from the registry again.

Usage::

    python3 tools/build_claim_views.py           # write the views
    python3 tools/build_claim_views.py --check   # compare, write nothing
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "claims" / "claim_registry.yaml"
NOVELTY = ROOT / "claims" / "novelty_matrix.csv"
FALSIFICATION = ROOT / "claims" / "falsification_matrix.csv"
GRAPH = ROOT / "claims" / "dependency_graph.json"

NOVELTY_COLUMNS = (
    "claim_id",
    "closest_prior_work",
    "exact_overlap",
    "imported_step",
    "oph_specific_delta",
    "novelty_type",
    "owner_theorem_or_section",
    "falsifier",
    "prior_work_citations",
)
FALSIFICATION_COLUMNS = (
    "claim_id",
    "mathematical_falsifier",
    "physical_identification_falsifier",
    "phenomenological_falsifier",
    "scope_if_false",
)


def load_claims(registry: Path = REGISTRY) -> list[dict]:
    """Claim rows merged with their context entries (a view, not the stored form)."""
    data = json.loads(registry.read_text(encoding="utf-8"))
    context = data["claim_context"]
    return [{**claim, "_context": context[claim["claim_id"]]} for claim in data["claims"]]


def _csv(columns: tuple[str, ...], rows: list[dict]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(columns), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def render_novelty(claims: list[dict]) -> str:
    rows = []
    for claim in claims:
        context = claim["_context"]
        prior = context["prior_work"]
        summary = context.get("matrix_summary", {})
        rows.append(
            {
                "claim_id": claim["claim_id"],
                "closest_prior_work": prior["summary"],
                "exact_overlap": prior["overlap"],
                "imported_step": prior["imported_step"],
                "oph_specific_delta": summary.get("oph_specific_delta", claim["oph_specific_delta"]),
                "novelty_type": claim["novelty_type"],
                "owner_theorem_or_section": context["owner_section"],
                "falsifier": summary.get("falsifier", claim["falsifier"]),
                "prior_work_citations": "; ".join(
                    dict.fromkeys(key for record in prior["citations"] for key in record["keys"])
                ) or "internal",
            }
        )
    return _csv(NOVELTY_COLUMNS, rows)


def render_falsification(claims: list[dict]) -> str:
    rows = []
    for claim in claims:
        falsifiers = claim["_context"]["falsifiers"]
        rows.append(
            {
                "claim_id": claim["claim_id"],
                "mathematical_falsifier": falsifiers.get("mathematical", claim["falsifier"]),
                "physical_identification_falsifier": falsifiers["physical_identification"],
                "phenomenological_falsifier": falsifiers["phenomenological"],
                "scope_if_false": claim["scope_if_false"],
            }
        )
    return _csv(FALSIFICATION_COLUMNS, rows)


def render_graph(claims: list[dict]) -> str:
    edges = [
        {"from": dependency["claim_id"], "to": claim["claim_id"], "role": dependency["role"]}
        for claim in claims
        for dependency in claim["_context"]["depends_on"]
    ]
    graph = {
        "schema_version": 1,
        "nodes": [claim["claim_id"] for claim in claims],
        "edges": edges,
    }
    return json.dumps(graph, indent=2, ensure_ascii=False) + "\n"


def renderings(root: Path = ROOT) -> dict[Path, str]:
    claims = load_claims(root / "claims" / "claim_registry.yaml")
    return {
        root / "claims" / "novelty_matrix.csv": render_novelty(claims),
        root / "claims" / "falsification_matrix.csv": render_falsification(claims),
        root / "claims" / "dependency_graph.json": render_graph(claims),
    }


def stale_views(root: Path = ROOT) -> list[str]:
    stale = []
    for path, text in renderings(root).items():
        if not path.is_file() or path.read_text(encoding="utf-8") != text:
            stale.append(str(path.relative_to(root)))
    return stale


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="compare the committed views; write nothing")
    args = parser.parse_args(argv)
    if args.check:
        stale = stale_views()
        if stale:
            print("claim views differ from the registry: " + ", ".join(stale), file=sys.stderr)
            print("run `python3 tools/build_claim_views.py` and commit the result", file=sys.stderr)
            return 1
        print("claim views OK: novelty matrix, falsification matrix and dependency graph match the registry")
        return 0
    for path, text in renderings().items():
        path.write_text(text, encoding="utf-8", newline="")
        print(f"wrote {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
