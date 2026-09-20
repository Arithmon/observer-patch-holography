# OPH Claim Registry

The papers are the standalone source for every theorem, assumption, falsifier, and claim boundary.
This directory contains the machine-readable scientific registry that keeps those standalone
statements synchronized across the public stack.
It is not public reading-path material. Do not link the registry
from the top-level public README files; point readers to the papers, falsifiability map, and public
explainers instead.
README numeric summaries should distinguish source-only rows, empirical closures, compare-only
rows, and SI convention/display rows.

`claims/claim_registry.yaml` is written as JSON, which YAML 1.2 also accepts.
The validator parses it as strict JSON and rejects duplicate keys; a YAML loader
reads the same data but lets a duplicate key through silently, so consumers
should parse it as JSON.

The registry is part of the working process:

- `claims/axiom_registry.yaml` records the normative three-axiom identities,
  interfaces, exclusions, and scientific realization boundaries.
- `claims/claim_registry.yaml` records top-level claim IDs, owner papers, claim tiers, imported
  mathematics, OPH-specific deltas, assumptions, evidence, falsifiers, and survival rules.
- `claims/novelty_matrix.csv` maps each claim against prior work.
- `claims/falsification_matrix.csv` records mathematical, physical-identification, and
  phenomenological failure modes.
- `claims/dependency_graph.json` records cross-claim dependencies.
- `claims/assumption_dictionary.md` gives stable names to recurring assumptions.
- `claims/frozen_prediction_register.json` records frozen and pending
  prediction contracts; `tools/build_fz_registry.py` validates it and renders
  `docs/FROZEN_PREDICTION_LADDER.md`.
- `claims/emergent_instrument_register.json` records scientific simulation and
  measurement instruments; `tools/build_instrument_register.py` validates it
  and renders `docs/INSTRUMENT_REGISTER_V3.md`.
- `claims/selection_ledger.json` and
  `claims/physical_identification_registry.json` record scientific selection
  classes and physical-identification boundaries; `tools/build_selection_ledger.py`
  validates them and renders `docs/SELECTION_LEDGER.md`.
- `claims/gravity_premise_ladder.json` records the gravity premise-elimination
  rungs; `tools/build_gravity_ladder.py` validates it and renders
  `docs/GRAVITY_PREMISE_LADDER.md`.
- `claims/public_surface_quantitative_claims.json` controls quantitative
  statements on public summary surfaces and is checked by
  `tools/check_public_surface_claims.py`.
- `claims/active_surface_inventory.json` is a generated reachability
  projection written by `tools/check_axiom_consistency.py --inventory`; it is
  not a hand-edited registry or project-status surface.

The validator is:

```bash
python3 tools/check_claim_registry.py
```

It checks that the registry release ID matches `paper/release_info.tex`, that every claim has an
owner file and falsifier, that the novelty/falsification matrices and dependency graph contain
every canonical claim ID with no unknown IDs, that the one-row-per-claim novelty and DAG node
projections have no duplicates, and that paper sources do not depend on direct paths to this
registry. The falsification matrix may keep several independently scoped rows for one claim.

Three row-level contracts carry their own machine-checked declaration:

- Evidence medium. The artifact medium of every evidence path is read from its suffix, and an
  unlisted suffix fails closed. A theorem-asserting row (`conditional_implication`,
  `branch_entry`, `empirical_implementation`, `physical_establishment`) whose entire evidence
  list is prose declares `proof_medium: paper_prose`; the validator rejects that field on every
  other row, so a Lean-backed or run-backed row cannot acquire the prose label. An empty
  evidence list fails.
- Topical gate owners. Each row that the V3 topical-owner policy in
  `tools/check_claim_registry.py` names declares `required_topical_gate_owners`, a sorted subset
  of that row's own `gates`. The validator compares the declaration against the policy in both
  directions and rejects a policy key that names no registered claim, so a one-sided edit to
  either surface fails. This field is the one gate-named key admitted outside `gates`; every
  other gate-named key is rejected as a side channel.
- Owner medium. An owner file under `paper/`, `extra/`, `cosmology/`, or `flagship/` is a paper
  source and declares no medium. An owner outside those roots declares
  `owner_medium: protocol_record`, which is admissible only for
  `claim_class: emitted_artifact`.

The GitHub workflow runs the validator on registry changes and on public claim-surface changes.
When a pull request changes paper TeX or the README claim narrative, it must also touch this
registry/check surface. That rule keeps the registry from becoming a stale snapshot.

## Native-record theorem-count review (#918)

The README theorem-count floor is mechanically checked by
`tools/check_lean_theorem_count.py`; declaration counts are explicitly excluded
from physical quantitative claims by `public_surface_quantitative_claims.json`.
The 39 declarations in `Lean/Geometry/SourceNativeRecords.lean` and
`Lean/Geometry/SourceBusScaling.lean` bring the reviewed library to 11,022
declarations, requiring a public floor of 11,000 in both README languages.

The scientific scope review finds a finite sum-write construction under the
supplied scalar mean law and a polynomial cleanup bound on a supplied paired
chain. Preparation, rail embedding, controller, codebook and physical error
bounds remain inputs. The range obstruction applies to unscaled growing raw
values; the two-schedule witness is not a full A1--A3 countermodel. These
results leave the complete native q=13/q=21 compiler and M1 read-law selection
open. No claim-registry status or physical-identification class is promoted
by the corrected declaration count. The detailed evidence and limits are in
`code/source_native_updates/CONTRACT.md` and `code/source_native_updates/AUDIT.md`.
