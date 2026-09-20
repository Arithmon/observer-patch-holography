# INS-03: source-bound phase-sensitive readout instrument (design)

This is an unfrozen instrument
design, not a reverse-engineering-reality result document and not a simulator
receipt.

Status: design only. This document specifies an instrument; it arms none. No
seed is drawn, no freeze event exists, no run is authorized, no register row
is written, and no ledger or premise row changes. Execution requires an
immutable preregistration, a target-blind freeze, and explicit authorization
under the owning simulation-instrument lane
[#737](https://github.com/FloatingPragma/observer-patch-holography/issues/737).

Sources of record for every statement below: the deep audit
the 2026-08-20 post-r2020 deep audit (findings F1, F2, F4,
the fastest-discriminating-route list, and the residual-limitations list), the
V3.17 post-r2020 deep-audit correction, the RER Lean modules and
premise-register rows
cited by name in each section, the simulator surfaces
`oph_fpe/quantum/phase_operation.py`, `oph_fpe/qm_observer/DESIGN.md`, and
`oph_fpe/dynamics/source_operator_inventory.py`, and the campaign record
the 2026-08-20 simulator alignment record. RER content is read at committed head
`609f88ad`.

## A. Identity and scope

- Instrument id: INS-03. Owning lane: issue #737, which owns the instrument
  surface and the register `claims/emergent_instrument_register.json` with its
  generated view `docs/INSTRUMENT_REGISTER_V3.md`.
- Bound ledger row: OL-C5, "Phase completion and full operator tomography",
  lane issue 730, status `partial`, premises PR-02 and PR-04, open premises
  PR-03, PR-64, PR-65 (`tracking/observation_ledger.json`, generated view
  `docs/OBSERVATION_LEDGER_V3.md`). The register binds one instrument to
  exactly one ledger row; INS-03 binds OL-C5 and does not enter the
  OL-A1 lineage.
- Register status a later step would record: `SPECIFIED`, which the register
  policy defines as mutable design with no freeze artifacts, no freeze time,
  and no verdict receipts. This document creates no such row.
- Lineage predecessors on the register are INS-01 (`FAILED`, controlling for
  OL-A1) and INS-02 (`SPECIFIED`). INS-03 shares their format and their
  discipline and supersedes neither.

### Division of labour with the PR-64 Lean instrument

A separate session builds the PR-64 completely-positive instrument in Lean.
INS-03 does not duplicate that construction and states no channel theorem.
INS-03 specifies the simulator-side export that the PR-64 instrument's checks
consume. The export contract below is the whole interface obligation of this
instrument toward that work.

Candidate export schema id: `oph.sim.ins03_phase_instrument_export.v1`. One
JSON object, canonicalized as `json.dumps(..., sort_keys=True, indent=2)` over
content free of paths, timestamps, and environment data, hashed by the
committed `canonical_sha256` of `oph_fpe/core/charged_response.py` over the
JSON text without a trailing line break, matching the convention recorded in
the 2026-08-20 simulator alignment record.

Number encoding: every exact rational is a two-integer array
`[numerator, denominator]`, following `oph.sim.qm_observer_viz.v1` of
`oph_fpe/qm_observer/DESIGN.md` section 7; every element of `Q(sqrt(3), i)` is
one exact `C3` scalar in the encoding of `oph_fpe/quantum/phase_operation.py`,
the four-string list `[re_rational, re_sqrt3, im_rational, im_sqrt3]` written
by `C3.encode`, whose components are the two `Q3` parts `a + b*sqrt(3)` of the
real and imaginary halves and whose four strings are rational literals. A
matrix is written by `encode_matrix` as nested lists of those entries.
`encode_matrix` ranges over the two-by-two shape only, so a carrier of
dimension other than two is exported under the explicit extension named here,
`encode_matrix_n`: the same per-entry four-string encoding over `n` rows and
`n` columns in row-major order, with the integer `n` carried beside the
matrix. That extension is a producer-side addition specified by this design
and written by no committed module. No floating-point value carries receipted
content; float renderings, if present, sit inside a block labeled
`derived_for_display`.

Required fields, per context `c` in the frozen context list and per outcome
`i` of that context:

1. `outcome_maps[c][i]`: the Kraus family of the outcome map, as an ordered
   list of exact matrices. The Lean side obtains complete positivity from
   Kraus form through `OPH.Dynamics.isCompletelyPositive_of_kraus`
   (`Lean/Dynamics/ChoiCPTP.lean`), so the export supplies Kraus data rather
   than a positivity assertion.
2. `outcome_maps[c][i].effect_from_kraus`: the exact matrix
   `sum_k K_k^dagger K_k` for that outcome, and
   `outcome_maps[c][i].declared_effect`: the entry of the frozen effect list
   for that context and outcome, which for the committed eight-context table
   is the transcription of `phase_operation.named_effects()` taken at freeze
   time and for effects outside that table is the frozen effect literal. The
   count-producing path reads the frozen effect list and not that module,
   which section C.2 (c) excludes. Both matrices appear, plus the exact
   residual matrix of their difference, required to be exactly zero under the
   residual gate of section C.4. The PR-64 predicate that consumes this is the
   effect-agreement clause of the register statement of PR-64.
3. `outcome_maps[c][i].trace_nonincreasing`: for each element of a frozen
   spanning set of the matrix space (the four matrix units suffice by
   linearity), the exact values `Tr(Phi_{c,i}(X))` and `Tr(X)` with their
   exact difference. The Lean predicate is
   `OPH.Dynamics.IsTracePreserving` applied to the summed channel and its
   trace-nonincreasing weakening applied outcome-wise.
4. `summed_channel[c]`: the exact Kraus normalization identity
   `sum_i sum_k K_{i,k}^dagger K_{i,k} = 1` with its residual matrix, and the
   exact trace values `Tr(sum_i Phi_{c,i}(X))` against `Tr(X)` on the same
   spanning set. This is the trace-preserving summed-channel field named in
   PR-64 and in the V3.17 plan section.
5. `readback[c][i]`: the public outcome symbol the observer records, the
   integer count, the context mass, and the exact compatibility residual
   `Tr(rho E_{c,i}) - count_i / mass`, where `rho` is the produced
   preparation of field 6 and `E_{c,i}` is `effect_from_kraus`. This is the
   effect/readback compatibility field of PR-64. The residual is exported as
   an exact rational and is not an algebraic identity: it compares an exact
   Born weight against an integer count over a finite mass, so a produced run
   realizes it as a sampling residual. Its gate is the preregistered exact
   rational tolerance `TOL_READBACK` of section C.4, not exact zero.
6. `preparation`: the produced state in the general coordinatization of a
   Hermitian trace-one two-by-two matrix, that is, two real diagonal entries
   `rho_00` and `rho_11` summing to one, exported as exact rationals, and one
   off-diagonal coordinate `rho_01`, exported as one exact `C3` element, with
   the lower entry its conjugate; the exact positivity certificate
   `|rho_01|^2 <= rho_00 * rho_11`; a `record_diagonal` boolean certificate
   with the exact off-diagonal coordinate that justifies it; the carrier id;
   the source record ids consumed; `operations`, the ordered list of source
   operations applied to reach the state, each entry naming the operation, its
   carrier before and after, and an `in_class` boolean flagged against the
   constructor list of `EventAlgebra.SourceReachability.Reachable` (`base247`,
   `base88`, `baseTriple`, `step247`, `step88`, `stepTriple`, `swap`,
   `anchor247`, `anchor88`, `margSnd247`, `margFst247`, `margSnd88`,
   `margFst88`, `margSndTriple`, `margMiddle`, `expectL247`, `expectR247`,
   `expectL88`, `expectR88`, `scal`, `prodMarg`, `mix`), an operation matching
   no constructor carrying `in_class = false` with its own description; and
   `preparation_content_sha256` over the canonical preparation body with no
   context field inside it. The two coordinate facts this field exports are
   general: every Hermitian trace-one two-by-two matrix carries that
   coordinatization, and positive semidefiniteness of such a matrix is the
   stated modulus bound. `EventAlgebra.prep_coordinates` and
   `EventAlgebra.prep_offdiag_normSq_le` of
   `Lean/EventAlgebra/PhaseInstrumentDetermination.lean` are the
   committed-instance specializations of those two facts: both run through the
   `diagonal_counts_run` field of the ideal model, which pins the diagonal to
   `111/179` and `68/179` and the modulus bound to `7548/32041`. The export
   states the general form and imports no pinned diagonal.
7. `provenance`: producer module paths with sha256, simulator commit,
   canonical repository URL, RER commit pin, run id, the declared input
   inventory, the runtime read log, and the import-graph independence report
   in the style of `oph_fpe/qm_observer/DESIGN.md` section 6.
8. `labels`: `exploratory`/`evidential` booleans and the verbatim claim
   boundary, following the labeling of `oph_fpe/qm_observer/DESIGN.md`
   section 7 and the boundary paragraph of `oph_fpe/quantum/phase_operation.py`.

An export that omits any of fields 1 through 5 is outside the PR-64 interface
and returns the INCONCLUSIVE verdict of section C rather than a scientific
outcome. An export that carries fields 1 through 5 and no field 6
preparation and no field 7 provenance carries channel structure without
provenance, which is exactly the PR-64 without PR-65 case named in audit
finding F1.

## B. The preparation gap, stated exactly

### B.1 The committed source operation class produces record-diagonal states

The committed source operation class is the inductive predicate
`EventAlgebra.SourceReachability.Reachable` of
`Lean/EventAlgebra/SourceReachabilityDelimitation.lean`. Its base elements are
the three committed counted states `correlationState`, `correlationState88`,
and `tripleCorrelationState`. Its closure operations are the walk-step
transports `stepEvolve`, `pair88StepEvolve`, `tripleStepEvolve`, the slot
exchange `slotSwap`, the ambient anchorings `ambientEquiv` and
`ambientEquiv247`, the marginalizations `ptraceSnd`, `ptraceFst`, and
`ptraceMiddle`, the regional conditional expectations `leftSlotExpectation`,
`rightSlotExpectation`, and `scalarExpectation`, marginal-product formation on
the hinge alphabets, and convex mixtures. The module states which operations
are excluded: general Kraus or CPTP application, conjugation by arbitrary
unitaries, the declared Bell-state preparation and its transport, and any
operation writing off-diagonal entries.

The closure theorem is `reachable_isRecordLaw`, with direction
`Reachable c M -> IsRecordLaw M`, where `IsRecordLaw` is diagonality in the
record basis with nonnegative real weights summing to one. Its corollaries are
`reachable_recordDiagonal` and `reachable_state_certificate`. There is no
converse theorem and none is claimed. Audit finding F2 records the exact
error that would follow from reading the converse: the enumerated carriers of
`Carrier.Index` are `Fin 13`, `Fin 14`, `Fin 13 x Fin 13`, `Fin 13 x Fin 14`,
the triple carrier, `Fin 169`, and `Fin 182`, and they do not include `Fin 2`,
so no typed reachability witness exists for the two-by-two committed run
matrix. `EventAlgebra.committedRunState` in
`Lean/EventAlgebra/OperationalPhaseAttainment.lean` carries that statement in
its own docstring: the matrix has the diagonal shape of the delimited states,
and the implication runs from reachability to diagonality, not back.

On the simulator side the same shape holds by construction.
`oph_fpe/quantum/phase_operation.record_diagonal_state` builds exactly
`diag(first/mass, second/mass)`. `committed_core_state(x, y)` accepts one
off-diagonal coordinate under the positivity bound `x^2 + y^2 <= 111*68/179^2`,
and it is a declared constructor with no source producer behind it. Gap 3 of
the scout report in the 2026-08-20 simulator alignment record records that the
simulator state space is entirely real and record-diagonal.
`oph_fpe/dynamics/source_operator_inventory.py` carries the constant
`STATUS = "NO_REGISTERED_ACCEPTED_VERTEX12_BRIDGE_PACKET_ON_TRACKED_SERIALIZED_DATA_SURFACE"`,
so the tracked serialized-data surface holds no registered accepted
source-operator bridge packet.

### B.2 The committed phase context returns 1/2 on record-diagonal preparations

The exact Lean facts, with theorem names and module paths. The three
`EventAlgebra` modules cited here open the namespace `EventAlgebra` and no
module-named namespace inside it, so the identifiers carry no module segment
and the file path is the separate pointer.

- `EventAlgebra.bornWeight_rhoY_sub_of_isHermitian`
  (`Lean/EventAlgebra/PhaseInstrumentDetermination.lean`): for Hermitian `E`,
  `bornWeight rhoYPlus E - bornWeight rhoYMinus E = -(2 * (E 0 1).im)`.
- `EventAlgebra.phase_sensitivity_iff_offdiag_im` (same module): a Hermitian
  effect satisfies the phase clause exactly when `(E 0 1).im != 0`. Phase
  sensitivity is a property of the off-diagonal imaginary coordinate of the
  effect and of nothing else.
- `EventAlgebra.bornWeight_fin_two` (same module): the Born pairing in entry
  coordinates, `bornWeight rho P = rho 0 0 * P 0 0 + rho 0 1 * P 1 0 +
  rho 1 0 * P 0 1 + rho 1 1 * P 1 1`. On a trace-one record-diagonal `rho`
  the two off-diagonal terms vanish, so `Tr(rho E)` pairs the diagonal of `E`
  with the diagonal of `rho` and reads nothing else. The committed phase
  effect `OPH.QFT.sourcePhaseLift`, equal to the Pauli +Y projector by
  `OPH.QFT.sourcePhaseLift_eq_rhoYPlus`, carries the constant diagonal
  `(1/2, 1/2)`, so `Tr(rho E) = 1/2` exactly for every trace-one
  record-diagonal `rho`. This is the general Born-weight fact of this section
  and it assumes no count model.
- `EventAlgebra.frequency_affine` (same module): for every inhabitant of
  `IdealPhasePOVMCountModel` and every instrument context there are real
  `u, v, w` with
  `binaryFrequency (I.counts c) = u + v * (I.prep 0 1).re + w * (I.prep 0 1).im`.
  Every count frequency of that model is affine in the single free state
  coordinate.
- `EventAlgebra.phase_frequency_eq` (same module): over the committed core the
  phase-count frequency of an inhabitant equals `1/2 - (I.prep 0 1).im`. This
  theorem and `frequency_affine` quantify over inhabitants of the ideal model
  and run through its assumed `born_matches` fit field and its
  `diagonal_counts_run` field, which pins the diagonal to the committed run
  literals. They state what the ideal model forces, not what a produced run
  returns.
- `EventAlgebra.committedRunState_offdiag_zero`
  (`Lean/EventAlgebra/OperationalPhaseAttainment.lean`): the committed
  record-diagonal preparation has free coordinate zero.
- `EventAlgebra.rhoY_agree_on_diagonal` and
  `EventAlgebra.diagonal_receipts_do_not_identify`
  (`Lean/EventAlgebra/OperationalPhaseInstrument.lean`): the dual statement on
  the effect side, that diagonal effects read only the shared diagonal.

Composing them at the Born-weight level: on a trace-one record-diagonal
preparation the off-diagonal coordinate is zero, so by `bornWeight_fin_two`
every context Born weight is the diagonal pairing, and for the committed phase
effect it is exactly `1/2`. The value is independent of `(E 0 1).im`, which is
the entire content of the effect's phase sensitivity by
`phase_sensitivity_iff_offdiag_im`. On one record-diagonal preparation a
phase-sensitive effect and a phase-blind effect with the same diagonal carry
identical Born weights, so they are indistinguishable in the produced counts.
That same-diagonal pair is the scope of the statement: an effect whose
diagonal is not constant does separate two record-diagonal preparations with
different diagonals, and what this section states is blindness to the
off-diagonal coordinate, not blindness to every state distinction. Inside the
ideal model the same collapse appears as `frequency_affine` at free coordinate
zero and as `phase_frequency_eq` returning `1/2`. Since `reachable_isRecordLaw`
sends the committed source class into `IsRecordLaw`, and
`IsRecordLaw.recordDiagonal` sends it into record diagonality, any preparation
produced by that class carries free coordinate zero on any
record-diagonality-preserving two-dimensional compression.

The effect side of the same gap is registered separately.
`EventAlgebra.QuantumSurface.real_closure_blind` and
`EventAlgebra.QuantumSurface.phase_lift_outside_real_closure`
(`Lean/EventAlgebra/QuantumAdequacySurface.lean`) prove that the generous
phase-free real closure of the committed source effects is tomographically
incomplete and does not contain the committed lift. That is the reason PR-04
stands as a register row with disposition `axiomatize` rather than as a
consequence of the committed payload.

Therefore a source-bound phase instrument requires a preparation producer
outside the enumerated class. Deterministic recalculation of declared Born
weights supplies no such producer, which is the statement of PR-65 and of
audit finding F1.

### B.3 Candidate routes the audit leaves open

Route R1: contextual or non-jointly-diagonal readouts. Audit finding F4
records that the CHSH result covers four jointly record-diagonal readouts on
reachable record-diagonal states and one declared Bell-state target, and that
the contextual and non-jointly-diagonal readout routes stay outside the
theorem. The nonclaims paragraph of
`Lean/EventAlgebra/SourceReachabilityDelimitation.lean` states the same scope.
What must be produced: a source-produced state on an enumerated carrier, a
readout family that is not jointly diagonal in the record basis, context
selection driven by a source record rather than by the analyst, and public
outcome counts per context. The discriminating quantity is the failure of a
single noncontextual valuation on the produced effect set, measured through
the frozen estimator. What falsifies it: an exhibited noncontextual assignment
reproducing every produced frequency; or a produced readout family that
commutes with the record projector family, in which case the run has produced
jointly diagonal readouts and tests nothing outside the committed bound; or a
produced frequency set consistent with the committed countermodel
`producedCubicValuation` of `Lean/EventAlgebra/OperationalAdditivityBoundary.lean`.

Route R2: the PR-44 product split under slot-local settings that are not
jointly record-diagonal. The jointly record-diagonal case is decided and sits
outside this route: `reachable_diagonal_chsh_le_two` of
`Lean/EventAlgebra/SourceReachabilityDelimitation.lean` bounds every reachable
state on every enumerated carrier by CHSH value 2 against every four jointly
record-diagonal readouts valued in `[-1, 1]`, so a run confined to that case
tests nothing this route is for. R2 is therefore scoped to slot-local settings
that fail joint diagonality in the record basis. Its explicit prerequisite is
the theorem that item 3 of the audit's fastest scientifically discriminating
completion route names and that the nonclaims paragraph of
`SourceReachabilityDelimitation` records as unproved: separability of source
record-diagonal product-basis states, with the CHSH consequence for arbitrary
slot-local settings. Until that theorem exists R2 carries no bound of its own
to test against, and its non-jointly-diagonal readout question is the question
of route R1. PR-44 is the finite subsystem split and normalized
local-operation interface, consumed by lanes 728, 730, and 740. What must be
produced: the prerequisite theorem; a source-produced bipartite state on a
committed pair carrier (`Fin 13 x Fin 13` or `Fin 13 x Fin 14`); slot-local
settings selected by source records on each slot, with the joint-diagonality
check on the produced setting family recorded; joint outcome counts; and the
separability certificate of the produced state. What falsifies it: produced
joint statistics exceeding CHSH value 2 against jointly record-diagonal
readouts contradict `reachable_diagonal_chsh_le_two` and indicate a producer
defect or an operation outside the enumerated class, which the run must then
name; produced statistics bounded by the prerequisite theorem's value for
every slot-local setting falsify the hope that the product split alone
supplies phase-sensitive discrimination, and the phase question returns
SOURCE_PRODUCER_MISSING under section C.

Route R3: the declared lift used as readback rather than preparation. PR-04
stays a declared effect; the declared Pauli +Y projector enters the instrument
on the readback side as a Lueders outcome map (`EventAlgebra.luedersUpdate`,
`Lean/EventAlgebra/Lueders.lean`), with complete positivity through
`OPH.Dynamics.isCompletelyPositive_of_kraus` and summed-channel trace
preservation in the shape of `OPH.Dynamics.partitionPinching_isCPTP`
(`Lean/Dynamics/ChoiCPTP.lean`). The PR-64 register note names a declared
Lueders instrument for the Pauli-Y projector as a useful conditional
construction that leaves source selection and provenance under PR-65. What
must be produced: CP trace-nonincreasing outcome maps, the trace-preserving
summed channel, effect/readback compatibility, a source-produced preparation,
and public outcomes. What falsifies it as a PR-65 discharge: the blindness of
section B.2. With a record-diagonal preparation the phase context returns
exactly `1/2`, so the run cannot separate the declared phase readback from any
phase-blind readback with the same diagonal. This route can satisfy the PR-64
interface of section A while returning SOURCE_PRODUCER_MISSING on the PR-65
question, and the decision rule states that outcome explicitly.

Route R4: a typed source-to-`Fin 2` compression with an off-diagonal-writing
source operation. Audit finding F2 states the required next theorem: construct
a typed source-to-`Fin 2` map and prove that the target matrix is its image.
What must be produced: that typed map, a source operation admitted into the
class that writes a nonzero off-diagonal record coordinate (the exclusion list
of `SourceReachabilityDelimitation` names such operations as outside the
committed class), and the produced state's off-diagonal coordinate as one
exact `C3` element under the general positivity bound
`|rho_01|^2 <= rho_00 * rho_11` of the produced diagonal, of which
`EventAlgebra.prep_offdiag_normSq_le` is the committed-diagonal specialization
at bound `7548/32041`. What falsifies it: a proof that the candidate operation
preserves `IsRecordLaw`, which extends `reachable_isRecordLaw` and pins the
compression's free coordinate at zero; the route then returns
SOURCE_PRODUCER_MISSING for that operation set, and the falsifier is the
record-law preservation proof itself.

## C. The preregistered decision rule

The rule is target-blind and freezes before any run and before any seed draw.
The rule text below is the design of the rule; the numeric placeholders bind
at freeze time. This document binds none of them.

### C.1 Verdict grammar

Four verdict words are recorded verbatim in the run receipt: PASS, FAIL,
SOURCE_PRODUCER_MISSING, INCONCLUSIVE. The register requires a decision rule
naming REPLICATED, FAILED, and INCONCLUSIVE, so the frozen text also carries
the mapping: PASS maps to register status `REPLICATED`, FAIL maps to `FAILED`,
and both SOURCE_PRODUCER_MISSING and INCONCLUSIVE map to `INCONCLUSIVE`, with
the distinguishing word retained in the receipt and in every public sentence
about the run.

Evaluation order, fixed before the run and enforced by the analysis code:

1. Preparation and provenance gates (C.2 conjuncts a, b, c) and the
   discrimination gate (C.2 conjunct d) execute and are recorded before any
   outcome count is exposed to the analysis. A failure here returns
   SOURCE_PRODUCER_MISSING and the run reports no endpoint value.
2. PR-64 interface checks on the export fields of section A. A failure returns
   INCONCLUSIVE.
3. Controls of section D, excluding the engineering-conformance comparison of
   control D.3, which is not gating at this step. A control failure returns
   INCONCLUSIVE.
4. Endpoints: the phase endpoint of C.2 (e) and the PR-03 additivity endpoint
   of C.3, each reported separately, then the overall verdict.
5. The engineering-conformance comparison of control D.3, run after the
   endpoint values are recorded so that it cannot steer the analysis.
   Agreement supports no verdict, per section C.5; disagreement returns
   INCONCLUSIVE with a defect report and leaves the separately reported
   endpoint values standing in the receipt.

No component pass overwrites a failed control or a failed endpoint, and no
component failure erases a separately reported component pass. This is the
outcome-separation discipline of the INS-02 design and of the register policy.

### C.2 PR-65 criteria, quantitative

PR-65 states that one source construction supplies the common preparation used
across the phase and real contexts, selects the effect or instrument, produces
the public outcomes through that instrument, and binds preparation, context,
operation, readback, and receipt custody to one run.

(a) Common preparation identity across contexts by content hash. The producer
emits one preparation record per context in the frozen context list. Each
record carries `preparation_content_sha256` over a canonical body containing
the source record ids, the carrier id, the produced state coordinates in the
exact encoding of section A, and the producer commit, with no context field
inside the hashed body. Criterion: the count of distinct preparation hashes
over the frozen context list equals exactly 1. Any other count returns
SOURCE_PRODUCER_MISSING, because the contexts were then read on different
preparations and the cross-context comparison has no common preparation.

(b) Public outcome custody. Every outcome count in the export is a tally over
outcome record files emitted by the run before any analysis code executes.
Each record file carries its own sha256; the set of record hashes is written
into a custody manifest; the manifest hash is recorded in the run log before
analysis starts. Criteria, all exact integer or byte equalities: the export
count for each context and outcome equals the integer tally over the record
files; the manifest hash read by the analysis equals the hash recorded before
analysis; no outcome record file is added, removed, or modified between the
two readings. Any mismatch fails closed and returns SOURCE_PRODUCER_MISSING.

(c) Provenance chain from source records to outcome counts with no access to
the committed weights. The producer declares its complete input inventory. The
frozen excluded-input set contains, at minimum,
`code/phase_operation_producer/PHASE_OPERATION_RECEIPT.v1.json` of the
research repository, the whole module `oph_fpe/quantum/phase_operation.py`
with every committed run literal and declared state constructor it carries
(`RUN_COUNTS`, `RUN_MASS`, the reference-receipt constants,
`committed_core_state`, and `record_diagonal_state`), and the declared branch
tables of `oph_fpe/qm_observer/tables.py`. Criteria: the intersection of the
excluded-input set with the runtime read log is empty; the import-graph check,
built by AST parsing in the style of the independence receipt of
`oph_fpe/qm_observer/DESIGN.md` section 6, reports zero edges from the
count-producing path into the excluded module set and zero dynamic-import
calls on that path; no committed run literal appears as a constant of the
count-producing path; the preparation coordinates of field 6 are
reconstructible from the produced source records named in that field alone,
and the receipt carries that reconstruction as an exact recomputation over
those records; the source hashes of the inspected modules match the audited
bytes. A nonempty intersection, a nonzero edge count, a committed literal on
the path, or a preparation coordinate that no produced source record
reconstructs returns SOURCE_PRODUCER_MISSING, because the produced numbers
could then be a replay of declared weights rather than a production from
source records. The exclusion covers the whole module, so a run needing that
exact arithmetic carries an independent implementation of it on the
count-producing path; `phase_operation.py` stays a control surface under
section C.5 and enters the analysis only at the post-endpoint comparison of
control D.3.

(d) Discrimination gate. If (a) through (c) hold and the produced preparation
carries `record_diagonal = true` with off-diagonal coordinate exactly zero,
the rule returns SOURCE_PRODUCER_MISSING without inspecting the phase counts,
because the committed phase effect carries the constant diagonal `(1/2, 1/2)`,
so `Tr(rho E) = 1/2` exactly on every trace-one record-diagonal `rho` by the
entry expansion of `bornWeight_fin_two`, and no count pattern separates the
declared phase readback from a phase-blind readback with the same diagonal.
The gate further requires the produced off-diagonal coordinate to satisfy
`|Im(rho_01)| >= DELTA_PHASE`, the frozen minimum phase content that keeps the
Born prediction of (e) separated from `1/2` at the declared mass; a produced
preparation below that threshold returns SOURCE_PRODUCER_MISSING with the
produced coordinate recorded.

(e) Phase endpoint. The endpoint is not a restatement of gate (d). It is the
block-level residual between the produced phase frequency and the Born
prediction computed from the preparation coordinates of field 6, which gate
(c) requires to be produced from source records independently of the phase
counts: for the committed phase effect that prediction is `1/2 - Im(rho_01)`,
by the same entry expansion. The residual
`f_phase - (1/2 - Im(rho_01))` carries a multiplicity-adjusted interval at
familywise level `ALPHA_FAMILY` and a frozen equivalence half width
`EPS_PHASE`. PASS on this endpoint requires the adjusted interval to lie
wholly inside `[-EPS_PHASE, EPS_PHASE]` while gate (d) holds, which is a
produced phase-sensitive readout on a source-produced preparation carrying
phase content. FAIL on this endpoint requires the adjusted interval to exclude
zero by more than `EPS_PHASE` with every gate of C.2 and section D passing;
that is the kill band, a produced disagreement between the phase readback and
the independently produced preparation, and it is a scientific result about
the produced architecture, reported with the prominence of a positive
outcome.

### C.3 PR-03 criterion, quantitative

PR-03 states that a public frame valuation is additive on every required
coexistent effect sum, across measurement contexts.
`Lean/EventAlgebra/OperationalAdditivityBoundary.lean` locates the row
exactly: by `committed_coexistent_sum_iff` a coexistent sum formed inside the
committed effect set is an effect precisely when it equals the sure effect, and
by `committed_sums_carry_only_normalisation` additivity on those sums is
equivalent to per-context normalization and forces nothing.

Consequences for the frozen design:

- The frozen effect list must contain at least one required coexistent pair
  `(A, B)` with `A + B` an effect and `A + B` not equal to the sure effect.
  The committed eight-context table of `phase_operation.named_effects()` does
  not supply such a pair, so the instrument must produce outcome maps for
  effects outside that table. The committed witness pair is
  `halfWitness = (1/4) rhoYPlus + (1/8) 1` with
  `witnessSum = halfWitness + halfWitness`, both proved effects by
  `halfWitness_isEffect` and `witnessSum_isEffect`. By `witnessSum_coords`
  the lower off-diagonal coordinate is `(witnessSum 1 0).im = 1/4`, so the
  upper coordinate is `(witnessSum 0 1).im = -1/4`, which is nonzero, and by
  `phase_sensitivity_iff_offdiag_im` the pair is inside the phase-sensitive
  family this instrument concerns.
- The three quantities `v(A)`, `v(B)`, and `v(A + B)` must be produced on
  separately produced record blocks with disjoint source record ids and block
  seeds derived from the master seed by the frozen collision-free derivation,
  with execution order randomized within block. The test is the paired
  block-level contrast `v(A + B) - v(A) - v(B)` against zero, with a
  multiplicity-adjusted interval at `ALPHA_FAMILY` and an equivalence half
  width `EPS_ADD`.
- Rejection requirement against the committed countermodel. The committed
  non-Born valuation `producedCubicValuation` reproduces every fixture count
  frequency, maps every effect into `[0, 1]`, sends the sure effect to one,
  obeys the complement rule, and is additive on every committed coexistent
  sum, and it fails additivity exactly at the witness pair:
  `producedCubicValuation_halfWitness` gives `143/512`,
  `producedCubicValuation_witnessSum` gives `35/64`, so the additivity residual
  it produces is `35/64 - 143/256 = -3/256`. The same
  `producedCubicValuation_witnessSum` value stands against the run Born value
  `1/2` at that effect, a difference of `3/64`;
  `producedCubicValuation_deviates` states the inequality of those two values
  and not its size. The frozen rule must be able to reject that
  countermodel, so `EPS_ADD` binds strictly below `3/256` and the interval
  construction must resolve a residual of that size at the declared block
  count. A rule with `EPS_ADD >= 3/256` cannot separate the countermodel from
  Born additivity at the committed witness pair and is not a valid freeze.
- PASS on this endpoint requires the adjusted interval for the additivity
  residual to lie wholly inside `[-EPS_ADD, EPS_ADD]` while every gate of C.2
  passes. FAIL requires the adjusted interval to exclude zero by more than
  `EPS_ADD` with every gate passing, which reports a produced cross-context
  non-additivity.

### C.4 Overall verdict

- PASS: the PR-65 conjuncts (a), (b), (c) hold; the produced preparation is
  not record-diagonal and clears the discrimination gate (d); the phase
  endpoint residual of (e) lies inside `EPS_PHASE`; the PR-03 endpoint lies
  inside `EPS_ADD` with the countermodel rejectable; every PR-64 interface
  field of section A is present, every algebraic residual of fields 2, 3, and
  4 is exactly zero, and every field 5 readback residual is within
  `TOL_READBACK`; every control of section D passes.
- FAIL: every gate and control passes and at least one endpoint lands in its
  frozen kill band.
- SOURCE_PRODUCER_MISSING: any PR-65 conjunct fails, or the produced
  preparation is record-diagonal or fails the discrimination gate (d), or the
  `operations` list of field 6 carries no entry flagged `in_class = false`
  against the constructor list of `EventAlgebra.SourceReachability.Reachable`,
  so that every exercised source operation sits inside the enumerated class.
  The run states that it did not test the question. It discharges nothing and
  demotes nothing.
- INCONCLUSIVE: a control fails, a PR-64 interface field is absent, an
  algebraic residual of field 2, 3, or 4 differs from zero, a field 5 readback
  residual exceeds `TOL_READBACK`, a fail-closed error fires on the producer,
  estimator, or control path (a typed `PhaseOperationError` raised through
  `phase_operation.require`, including the `OUTCOME_POSITIVITY` code of the
  outcome table and the `WINDOW_SHAPE`, `WINDOW_POSITIVITY`, and
  `WINDOW_VIOLATION` codes of `phase_operation.check_phase_window`, which
  enforces the committed integer window `32041 (a - b)^2 <= 30192 (a + b)^2`,
  or the same-coded error of the independent count-producing implementation
  that reimplements those checks; the Lean statement of that window over
  committed-core inhabitants is `phase_counts_receipt_window`), or an endpoint
  interval neither lies inside its frozen equivalence band nor lands in its
  kill band.

Residual gate, stated separately for the two residual kinds. Fields 2, 3, and
4 carry algebraic identities of the exported matrices: `effect_from_kraus`
against `declared_effect`, the outcome-wise trace-nonincrease differences on
the frozen spanning set, and the summed-channel Kraus normalization
`sum_i sum_k K_{i,k}^dagger K_{i,k} = 1`. Each holds entry by entry in exact
arithmetic on any correct export of a produced run, so each residual matrix is
required to be exactly zero and one nonzero entry returns INCONCLUSIVE.

Field 5 is a different kind. Its residual compares an exact Born weight
against an integer count over a finite mass, so on a produced run it cannot be
exactly zero except by deterministic recomputation of the declared weight,
which section C.5 rules a control rather than evidence; a rule demanding exact
zero there would select for replay and reject production. Its gate is the
preregistered exact rational `TOL_READBACK`, applied as the exact comparison
`|Tr(rho E_{c,i}) - count_i / mass| <= TOL_READBACK` per context and outcome.
`TOL_READBACK` binds target-blind: it is the half width of the frozen interval
construction at level `ALPHA_INTERFACE` on the declared context mass, rounded
up to an exact rational whose denominator is that declared mass. Its inputs
are the declared context mass, `ALPHA_INTERFACE`, and the frozen interval
construction, each fixed before the seed draw. No produced count, no
preparation coordinate, no endpoint value, no control outcome, and no
committed count literal enters that computation, and the value binds under
section C.6 with the other placeholders.

### C.5 Controls, not evidence

Per audit finding F1 and item 1 of the audit's fastest-route list, a
deterministic Born-table calculator is only a control. The following surfaces
are CONTROLS only under this rule and carry no promotion authority:

- `oph_fpe/quantum/phase_operation.py`, which replays the committed count
  package in exact `Q(sqrt(3), i)` arithmetic and states on its face that the
  counts are produced by the declared semantics and never measured, together
  with the research-repository producer and verifier under
  `code/phase_operation_producer/`;
- the `oph_fpe/qm_observer/` record-ensemble consistency package, whose
  `DESIGN.md` records that the branch tables of section 3.3 are declared
  transcriptions of the conditional weights `Tr(p F)` and that the three-way
  identity is a consistency receipt over the transcription, the refinement
  rule, and the enumeration, not an independent confirmation.

Agreement between either surface and a run supports no verdict. Disagreement
is an engineering defect and returns INCONCLUSIVE with a defect report.

### C.6 Placeholder binding rule

`DELTA_PHASE`, `EPS_PHASE`, `EPS_ADD`, `ALPHA_FAMILY`, `ALPHA_INTERFACE`,
`ALPHA_CONTROL`, `TOL_READBACK`, the block count, the context list,
the frozen effect list, the excluded-input set, the interval construction, the
multiplicity family, and the missing-data rule each bind to a value before any
seed draw. A seed drawn while any placeholder is unbound voids the
registration under the owning lane's freeze discipline, and a frozen band
admits no post-draw revision. The block count is set by a target-blind
prospective precision calculation from declared minimum effects, interval
widths, multiplicity, and attrition, published with its inputs; no convenient
count is presumed sufficient. Runtime-only pilots may measure cost and may not
expose preparation coordinates, outcome counts, endpoint values, or control
outcomes.

## D. Controls and negative controls

Every control runs through the same observable code path and the same fixed
estimator as its matched main cell. A structural verifier may check the
invariants of a transform and may not replace the estimator in an inferential
gate.

1. Record-diagonal source run, checked at two levels. Born-weight level: for
   a run whose produced preparation carries the record-diagonal certificate,
   the analysis recomputes `Tr(rho E)` in exact arithmetic from the exported
   preparation coordinates of field 6 and the exported `effect_from_kraus` of
   field 2, in every context. For the committed phase effect, whose diagonal
   is the constant `(1/2, 1/2)`, that value is exactly `1/2` on every
   trace-one record-diagonal `rho` by the entry expansion of
   `bornWeight_fin_two`; in every other context it is the diagonal pairing of
   the exported effect with the exported preparation. Any deviation of the
   recomputed Born weight from those exact values is an implementation defect
   and returns INCONCLUSIVE for the campaign until repaired. Count level: the
   produced frequency of such a run is a sampling quantity, so the check is
   interval coverage rather than exact equality. The frozen interval
   construction at preregistered level `ALPHA_CONTROL`, applied to the
   produced counts of a context, must cover the Born weight recomputed for
   that context, which is `1/2` in the phase context. A coverage failure
   returns INCONCLUSIVE for the campaign until repaired. `phase_frequency_eq`
   and `frequency_affine` are cited here for the ideal-model identity only:
   they quantify over inhabitants of `IdealPhasePOVMCountModel`, whose assumed
   `born_matches` field equates frequency with Born weight by construction and
   whose `diagonal_counts_run` field pins the diagonal to the committed run
   literals, so neither states what a produced run returns. This control tests
   the section B.2 blindness fact against the running code.
2. Shuffled-context control. The context labels attached to outcome records
   are permuted by a frozen permutation after production and before analysis.
   The analysis must return no phase endpoint residual inside `EPS_PHASE` and
   no additivity residual outside `EPS_ADD`. Either outcome under shuffling
   shows the estimator reads a label artifact and returns INCONCLUSIVE.
3. Declared-table control. The deterministic calculator of
   `oph_fpe/quantum/phase_operation.py` and the committed reference receipt
   `code/phase_operation_producer/PHASE_OPERATION_RECEIPT.v1.json` are
   compared against the run counts as an engineering conformance check.
   Agreement supports no verdict, per section C.5; disagreement returns
   INCONCLUSIVE with a defect report. This control is not gating at step 3 of
   the evaluation order: the comparison runs at step 5, after the endpoint
   values are recorded, so it cannot steer the analysis.
4. Classical reveal-only mock. The reveal-only mock of the
   `oph_fpe/qm_observer` tests carries mediated counts equal to its direct
   counts, so its own interference gap is exactly `0`, while the receipted gap
   of the record-ensemble package is `129/1432` at common mass 2864
   (`DESIGN.md` section 5); the interference check rejects the mock because
   `0` differs from that receipted value. The mock passes through the same
   analysis path under a frozen mapping: its mediating context maps to the
   frozen context list entry for `web_conjugated_3` and its base context to
   the entry for `web_diagonal`, and its two outcome labels map to the two
   effects of each of those entries in the frozen effect list, so the
   estimator reads mock counts in the same fields it reads run counts. The
   analysis must reject it. Failure to reject invalidates the estimator and
   returns INCONCLUSIVE.

## E. Custody

Hashed at freeze, as `{path, sha256}` pairs under the register's artifact
rule, with simulator paths carrying the canonical repository URL and full
commit pin:

- the verbatim decision-rule text with every placeholder bound;
- the frozen context list, the frozen effect list including the coexistent
  pair, the excluded-input set, the interval construction, the multiplicity
  family, and the missing-data rule;
- the producer modules, the estimator module, the control modules, and the
  independent verifier, each by source sha256;
- this design document at its frozen revision, as the specification pointer;
- the simulator commit and the RER commit pin, with the audited RER head
  `609f88ad` recorded as the reading baseline of this design;
- the environment record.

Hashed at run time and retained: the preparation bodies with their content
hashes, every outcome record file, the custody manifest, the runtime read log,
the import-graph independence report, the export body, the receipt body with
its own `receipt_sha256` over its canonical body, and the raw captures
sufficient for an independent implementation to recompute every reported
observable.

Seeds: none is drawn by this document. If a run is ever authorized, one fresh
master seed is drawn at freeze time from a recorded source, with a frozen
collision-free derivation to block and cell seeds, paired blocks spanning every
context, every effect of the coexistent pair, and every control, execution
order randomized within block, no redraw, no cell substitution, no optional
stopping, and one recorded runtime environment.

Evidential versus exploratory: a run is evidential only when the register row
stands at `FROZEN` with a non-future UTC freeze time, hashed freeze artifacts,
and a decision rule naming REPLICATED, FAILED, and INCONCLUSIVE; when the run
configuration is reproducible from committed artifacts; and when raw captures
are retained for independent recomputation. The reproducibility boundary
recorded for INS-01 is the failure mode this design avoids: that campaign
retained the manifest, summary, and fifteen per-arm receipt files and no raw
feature matrices or fit captures, so no independent implementation could
reconstruct every observable, and the simulator's validate-only mode rechecks
receipt fields through the producer implementation rather than recomputing
observables independently. Every run outside those conditions is exploratory
and non-evidential, in the sense that `oph_fpe/qm_observer/DESIGN.md` and the
boundary paragraph of `oph_fpe/quantum/phase_operation.py` state on their face.

Standing verdicts: INS-01 remains the controlling completed verdict for OL-A1
and its pins are immutable, as recorded in the ledger-control lineage of
`docs/INSTRUMENT_REGISTER_V3.md` and in the 2026-08-20 simulator alignment record.
INS-03 binds OL-C5 and does not enter that lineage. This design changes no
ledger row, no premise row, and no register row: OL-C5 stays `partial` with
open premises PR-03, PR-64, and PR-65, and PR-04 stays consumed under its
`axiomatize` disposition.

## F. Relation to the ledger

- Ledger row that INS-03 could touch: OL-C5 only. A PASS could support OL-C5
  for the exact frozen run configuration, after an independent verifier
  recomputes the declared observables from retained captures, and after the
  register's ledger-control lineage explicitly selects INS-03 as controlling
  for that row. A controlling FAIL demotes the affected clause of OL-C5 with
  equal prominence. SOURCE_PRODUCER_MISSING and INCONCLUSIVE move the row in
  no direction.
- PR-04, declared phase-sensitive effect, disposition `axiomatize`: consumed,
  never discharged. INS-03 uses the declared effect and derives it from
  nothing. No verdict of this instrument removes the row.
- PR-64, operational phase-instrument channel realization, disposition
  `remove`: owned by the companion Lean instrument. INS-03's export supplies
  the sim-side outcome maps, summed-channel trace, and effect/readback
  compatibility fields those checks consume. Satisfying the section A
  interface discharges nothing by itself on PR-65.
- PR-65, source-produced phase preparation and readback provenance,
  disposition `remove`: only a source-bound run under this rule, after a
  target-blind freeze and explicit authorization, and only through a
  preparation carrying a nonzero off-diagonal coordinate, that is through
  route R4, can discharge it. No design, no export, no interface check, and no
  deterministic replay discharges it.
- PR-03, operational effect additivity, disposition `remove`: the section C.3
  endpoint on independently produced records can falsify additivity at the
  witness pair. Removal requires the operational derivation named in the
  PR-03 register note, which records that removal "needs a source-produced
  complex-tomographically-complete effect/instrument system with an
  operational theorem deriving noncontextual additivity on all required
  coexistent sums" and that "removal needs an operational derivation
  beyond context-wise normalization". No endpoint of this instrument supplies
  that derivation. Until it exists, `producedCubicValuation` keeps the row
  load-bearing, per the PR-03 register notes and
  `OperationalAdditivityBoundary.lean`.
- PR-44 and PR-52 are outside this instrument. PR-44 is named only as route R2
  of section B.3. PR-52 physical attachment stays open under every verdict;
  INS-03 makes no physical-region, spacelike-separation, or spacetime claim.

## G. Non-claims

- No instrument is armed. No seed is drawn. No freeze event exists. No run is
  authorized by this document.
- No ledger row, premise row, register row, or claim row changes. This
  document writes to one path and to nothing else.
- No premise is discharged. PR-03, PR-64, and PR-65 stay open, and PR-04 stays
  a declared row under its recorded decision of 2026-08-18, lane issue 730.
- Nothing here is a measurement of a physical system. The instrument concerns
  what simulated observers' records exhibit inside the architecture, which is
  the register's stated separation from the frozen-prediction ladder.
- The deterministic count package stays a static semantic-conformance fixture
  and the record-ensemble package stays a consistency receipt over declared
  transcriptions. Neither is evidence under this rule.
- No source reachability is claimed for any two-by-two preparation. The Lean
  direction is `Reachable -> IsRecordLaw`, and audit finding F2 records that
  the enumerated carriers exclude `Fin 2`.
- No Bell, locality, contextuality, or no-signalling result is claimed by this
  design. Routes R1 and R2 are stated as candidate routes with their own
  falsifiers.
- INS-01 stays the controlling verdict for OL-A1, and OL-A1 stays owed. INS-03
  supersedes no instrument by existing as a design.
- Nothing here establishes OPH as the correct fundamental theory, which the
  audit's executive verdict states directly.
