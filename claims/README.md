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
  and renders `docs/registers/INSTRUMENT_REGISTER_V3.md`.
- `claims/selection_ledger.json` and
  `claims/physical_identification_registry.json` record scientific selection
  classes and physical-identification boundaries; `tools/build_selection_ledger.py`
  validates them and renders `docs/registers/SELECTION_LEDGER.md`.
- `claims/gravity_premise_ladder.json` records the gravity premise-elimination
  rungs; `tools/build_gravity_ladder.py` validates it and renders
  `docs/registers/GRAVITY_PREMISE_LADDER.md`.
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

## Numerical audit consolidation, 2026-10-08

The review of [#1052](https://github.com/FloatingPragma/observer-patch-holography/pull/1052),
[#1057](https://github.com/FloatingPragma/observer-patch-holography/pull/1057),
[#1058](https://github.com/FloatingPragma/observer-patch-holography/pull/1058),
[#1059](https://github.com/FloatingPragma/observer-patch-holography/pull/1059) and
[#1060](https://github.com/FloatingPragma/observer-patch-holography/pull/1060)
checks their combined effect on the existing claims. All 326 claim payloads,
including their assumptions, dependencies, gates, statuses and falsifiers,
remain unchanged. Frozen prediction targets, decision rules and custody bytes
are preserved. These repairs add no physical-model prerequisite.

| Evidence family | Sharper supported result | Public-surface disposition |
| --- | --- | --- |
| Petz recovery and tomography | Complete counts survive a later analysis failure; valid CMI remains available when recovery support is unresolved. Strict support validation and the squared-fidelity convention remain in force. | The Standard Model paper and IBM overview label old Stage 1 metrics as archived and unrevalidated because complete counts/states are missing. Historical values and the separate issue-509 claim are preserved. |
| Scalar Bose readout | The equilibrium identity `p V = 2 sum(G) / 3` retains small pressure; nonzero free energies and positive mode occupations cannot disappear behind an absolute verifier tolerance. | Existing equations and supplied-action/ensemble scope remain correct. The live receipt is regenerated and independently checked. |
| DESI posterior accounting | Exact decimal weights, endpoint tests, centered moments and separately accumulated subset masses preserve narrow spread and positive excluded tails. | The retrospective receipt and its downstream ledger are regenerated. Existing comparison classifications and the prospective FZ-13/FZ-15 commitments remain unchanged. |
| Neutral quantum packet | Stable scalar invariants preserve norm suppression and normalized radius, including when only the radius is reportable. | Existing packet formulas, supplied preparation and unproved quantum propagation status remain correct. The live preparation receipt is independently replayed. |
| Finite collar certificates | All declared influences and multiplicities contribute; the existing `1/2` calibration and `3/8` witness floors remain exact. | The Yang--Mills paper states the shared checker boundary: supplied finite arithmetic cannot authenticate physical-source or continuum evidence. |

The paper, flagship and book passages state the applicable thermal,
packet and cosmology premises. They need no new qualifications from these
numerical fixes. Local replay of the repaired recovery producer is software
evidence and supplies no new hardware measurement. The older Stage 1 summaries
remain archived, without a claim that their measured numbers survived a
recalculation. Validation and exact reviewed heads are recorded in the
[standing audit](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).
