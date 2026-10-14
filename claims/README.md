# OPH Claim Registry

The papers are the standalone source for their theorems, assumptions,
falsifiers and claim boundaries. Executed constructions may instead name a
self-contained protocol record as their owner. This directory contains the
machine-readable scientific registry that keeps those statements synchronized
across the public stack.
It is not public reading-path material. Do not link the registry
from the top-level public README files; point readers to the papers, falsifiability map, and public
explainers instead.
README numeric summaries should distinguish source-only rows, empirical closures, compare-only
rows, and SI convention/display rows.

`claims/claim_registry.yaml` is written as JSON, which YAML 1.2 also accepts.
The validator parses it as strict JSON and rejects duplicate keys; a YAML loader
reads the same data but lets a duplicate key through silently, so consumers
should parse it as JSON.

## The five ledgers

The repository keeps its scientific bookkeeping in five ledgers. Each ledger
has one machine source; the other files are views rendered from it or
sub-registers validated with it, and the checks fail when a view drifts from
its source.

| Ledger | Machine source | Views and sub-registers |
| --- | --- | --- |
| Claims | `claims/claim_registry.yaml`: the claim rows and a `claim_context` section with each claim's prior work and citations, owner section, falsifiers and dependencies; `claims/prior_work.json`, the verified bibliography of that prior work; `claims/assumption_dictionary.md`, the assumption names | `claims/novelty_matrix.csv`, `claims/falsification_matrix.csv` and `claims/dependency_graph.json` (rendered by `tools/build_claim_views.py`); the credited-prior-work section of each paper (`tools/build_prior_work.py`); the selection and physical-identification registers `claims/selection_ledger.json` and `claims/physical_identification_registry.json` (validated and rendered to `docs/registers/SELECTION_LEDGER.md` by `tools/build_selection_ledger.py`); the public quantitative-display policy `claims/public_surface_quantitative_claims.json` |
| Premises | `tracking/premise_register.json`: the premise rows, with the model-supply crosswalk, the Standard Model coverage and the gravity rungs as its `model_supply`, `sm_coverage` and `gravity_rungs` sections; the theorem-library premises A1-A3 in `claims/axiom_registry.yaml`, stated in `docs/AXIOM_REFERENCE.md` | `docs/registers/PREMISE_REGISTER_V3.md`, whose closing sections show the Standard Model coverage and the gravity rungs; `model/PREMISES.md`, the model view of the crosswalk |
| Predictions and tests | `claims/frozen_prediction_register.json`, the custody record of frozen predictions, and `tracking/observation_ledger.json`, whose `instruments` and `constants` sections hold the emergent-instrument register and the constants ancestry | `docs/FROZEN_PREDICTION_LADDER.md`, `docs/registers/OBSERVATION_LEDGER_V3.md`, `docs/POSTDICTION_LEDGER.md`, `tracking/null_model_scorecard.md`; `docs/OPH_FALSIFICATION_PROGRAM.md` summarizes the falsification rules for readers |
| Evidence | the per-package manifests and verifiers under `evidence/`, `code/` and `model/`, the receipt-promotion and external-data registries under `code/audit/`, and `Lean/docs/PROOF_INDEX.md` | `evidence/README.md` |
| Papers | `paper/content_inventory.yml` and the paper list in `tools/build_tex_papers.py` | `paper/paper_release_manifest.json`, `docs/PAPER_SET.md` |

Model construction uses the existing claim identities and premises through
[Model construction](../model/CONSTRUCTION.md) and the generated
[model–theorem premise map](../model/PREMISES.md). Candidate hypotheses and
issue completion do not change a registered claim's scientific status. Record
exact same-model attachments and retained inputs there; change this registry
when the established claim or its evidence changes.

`claims/active_surface_inventory.json` is a generated reachability projection
written by `tools/check_axiom_consistency.py --inventory`; it is not a ledger.

The validator is:

```bash
python3 tools/check_claim_registry.py
```

It checks that the registry release ID matches `paper/release_info.tex`, that every claim has an
owner file and falsifier, that the novelty/falsification matrices and dependency graph contain
every canonical claim ID with no unknown IDs, that the one-row-per-claim novelty and DAG node
projections have no duplicates, that those three views match their rendering from the registry,
that every established result a claim uses is cited from a verified entry of
`claims/prior_work.json` or marked as a result of the project itself, that the papers'
credited-prior-work sections are current, and that paper sources do not depend on direct paths
to this registry. The falsification matrix may keep several independently scoped rows for one claim.

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

## Finite model constructions, 2026-10-12

Six `emitted_artifact` rows register the runnable model and its retained
evidence. Their owners are the self-contained implementation documents under
`model/`, with `owner_medium: protocol_record`; no paper theorem is promoted
or changed by these additions.

| Claim | Executed scope |
| --- | --- |
| `OPH-MODEL-NATIVE-OBSERVER-EXECUTION` | [Native current recorder](../model/observer/README.md), exact classical preservation, supplied forecasts and feedback |
| `OPH-MODEL-OBSERVER-RECORD-CAPACITY` | [Finite correctable capacity](../model/observer/CAPACITY.md) and lossless reached-prefix controller encoding |
| `OPH-MODEL-P-RESPONSE-EXECUTION` | [P-dependent binary detector](../model/observer/p_response/README.md), with accepted closure input and declared running/matching map |
| `OPH-MODEL-N4-LOCAL-TICK-CERTIFICATE` | [Incident-attempt memory certificate](../model/local_tick/README.md), with supplied occupation preparation and readout |
| `OPH-MODEL-SCALAR-RESPONSE-EXECUTION` | [Scalar detector and native continuation](../model/observer/scalar_response/README.md), with supplied field action/vacuum and inherited history refinement control |
| `OPH-MODEL-COMMON-SOURCE-EXECUTION` | [Common interacting observer](../model/common_source/README.md), complete scalar/matter/flux records and feedback, observer readouts and conditional error bounds |

The detector consumes PR-17 and PR-37. Its separate physical same-quantity
and Thomson readback attachments are undischarged PR-38 and PR-39. The other
five rows certify their supplied finite contracts without consuming the
canonical physical premises. Local tape capacity is not cosmic N closure;
the native current recorder does not implement the N4 occupation apparatus.
Uniqueness of a designed rule is not required for these implementation claims.
Existing claim payloads, physical statuses and frozen predictions are preserved.
The scalar adapter's conditional detector protocol is independently testable
in software. Its scalar action, vacuum, apparatus and model clock remain
declared; it establishes no new physical prediction or laboratory verdict.

## Continued-observer maintainer audit, 2026-10-14

PRs [#1099](https://github.com/FloatingPragma/observer-patch-holography/pull/1099),
[#1100](https://github.com/FloatingPragma/observer-patch-holography/pull/1100) and
[#1101](https://github.com/FloatingPragma/observer-patch-holography/pull/1101)
extend `OPH-MODEL-COMMON-SOURCE-EXECUTION` in dependency order:

- [Composition](../model/common_source/COMPOSITION.md) executes four reads,
  later evolution and a gauge clock probe on one prepared work state; its
  fixed-program smooth occupation bound controls the complete output.
- [Quantum-native execution](../model/common_source/QUANTUM_NATIVE.md)
  compiles the surrounding complete C/E instruments and feedback-selected
  routing, including all abort and padding outputs. It retains quantum
  states within record branches, without preserving cross-record coherence.
- [Active encoded memory](../model/common_source/ENCODED_MEMORY.md) compiles
  encoding and record-selected recovery and preserves a logical qubit and
  its reference correlations before the optional destructive read. This is
  a supplied active protocol, not protection by passive repair dynamics.

The existing claim identity remains an `emitted_artifact`. Its protocol
records own the detailed statements; the construction and flagship papers
summarize them in the refactored paper set. Canonical physical-premise
statuses, frozen prediction bytes and observation dispositions are unchanged.
The two executed occupation cutoffs have complete-output occupation
budgets of one. Physical services, observer-clock calibration, interacting
spatial refinement and independent observations remain open. The old sharp
history, its infinite signal certificate and its later clock probes keep
their distinct scopes; the new smooth bound does not certify those probes.

## Executable-model premise audit, 2026-10-12

The compact public README summarizes the registered results for agreement,
quantum observations, carrier geometry, gauge symmetry, thermodynamics,
record capacity and fields, alongside the common-source execution. Its
computational, observation-first interpretation states the programme's
proposal and motivation. It promotes no claim or physical identification.
Physical identification here means agreement with observers' measurement
statistics and interventions. It does not require an independently existing
classical world or fundamental fields and spacetime.
A continuum certificate concerns convergence and controlled errors of those
operational readouts, without asserting a fundamental continuous substrate.
The existing input assumptions of `OPH-MODEL-COMMON-SOURCE-EXECUTION` retain
their names and dependencies. The finite continuation, coherent native
attachment and active encoded-memory variants add their own explicit
contracts. The distinction is between consequences
proved from a construction and inputs selected to define that construction:

| Derived or executed consequence | Inputs that the result retains |
| --- | --- |
| Exact continuity and conditional current mean and variance | Integer seam rule and fair placement coin |
| Full-kernel mean diffusion and rank-three slow response | Fixed twelve-port incidence, uniform seam calibration and geometric readout interpretation |
| Complete record preservation and certified coupled feedback | Quantum algebra, specified action, instruments and coherent source services |
| Schedule-independent completed history | Fixed event graph, boundary, source maps and fair constraint-repair schedule |
| Observer Hilbert representations and modular flow | Specified observer algebra and its state restriction |
| Alignment of Hamiltonian and modular flow | Chosen full-work Gibbs reference |

The source-discrimination construction compares two complete finite
measurement channels and implements their difference as a recorded parity
bit. Its executed trial counts are synthetic, and coherent preparation,
readout and load-to-code identification remain declared inputs. This result
adds no independent physical data or canonical premise disposition.

The consequences in this table are registered conclusions, not extra input
premises to delete. Stationarity of the uniform binary distribution does not
derive its preparation. Setting the thermal scale parameter to one is a
unit convention; it does not select the Gibbs reference or derive the
interacting Hamiltonian. Compiling a chosen apparatus into admitted source
controls establishes implementation while retaining the chosen action and
the physical service contract.

The [native field derivation](../model/common_source/NATIVE_FIELDS.md) excludes
undamped wave modes for the stated uniform mean readout and forgotten-record
quantum channel. It does not identify repair-attempt count with subjective
modular time or exclude oscillatory modular correlations on a restricted
history state. A wave claim requires an actual observable/clock pairing and
its propagation law. No physical wave, continuum or empirical premise is
discharged by that restriction.

## Audit #1078 attachment corrections

The [construction guide](../model/CONSTRUCTION.md#retained-audit-controls)
records the dispositions of issues #1087–#1090. Their exact controls reject
the whole-simple-algebra central-port identification, the original
same-marked-action triplet identification, unrestricted extension of the
single-product feedback exclusion, and selection of inertia from the
configuration Gibbs marginal alone.

The [addressed native/N4 construction](../code/native_response/README.md)
supplies a separate twelve-block carrier and a fixed-full-reference abstract
kinetic counterpair under a restricted collective response grammar. Its
occupied-port instrument makes the information loss and retained-record
requirements explicit. These results are linked in the model–theorem map;
they change no existing claim ID, canonical inertia hypothesis, physical
premise disposition, observation rung or frozen prediction.

## Numerical audit consolidation, 2026-10-08

The review of [#1052](https://github.com/FloatingPragma/observer-patch-holography/pull/1052),
[#1057](https://github.com/FloatingPragma/observer-patch-holography/pull/1057),
[#1058](https://github.com/FloatingPragma/observer-patch-holography/pull/1058),
[#1059](https://github.com/FloatingPragma/observer-patch-holography/pull/1059) and
[#1060](https://github.com/FloatingPragma/observer-patch-holography/pull/1060)
checks their combined effect on the existing claims. All 326 claim payloads of that date, including their assumptions, dependencies, gates, statuses and falsifiers,
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
