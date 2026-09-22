# Maintainer-style scope and correctness review

This review covers the local-clock reduction and its finite controls against
`main` commit `992486deb21dd64fd58a4afb2056a431327e2e59`, including the
maintainer's registry-scope and CI-cache corrections. The work is independent
of the unmerged implementation in PR #968.

## Claims examined

* **Cone selection:** finite symmetry fixes an invariant quadratic form;
  it does not create local directional coverage. Uniform quadratic regularity,
  operational attachment of the response ideal and complete time-faithful
  refinement are explicit premises. The proof includes all control/timing
  dependencies. The nonlinear-clock example removes exact scaling without
  claiming finite null-boundary attainment.
* **Coherence:** the twelve-central-port tensor model and copy channel are
  finite witnesses. An explicit joint observer resolves state compatibility
  for that finite family. Neither a dormant support complex nor a prescribed
  response subalgebra proves complete A1--A3 membership. These conclusions
  are excluded in the paper, contract and claim registry.
* **A3 join:** actual gate equations define a classical history state problem.
  Conditional relative entropy reproduces the input ensemble. A deliberately
  defective program remains feasible and fails the separate truth-table
  check. No successful-read condition defines the feasible family.
* **Programmability:** NOT/Toffoli suffice. The constant-control bit used to
  lower CNOT is prepared, cleaned and charged. Inputs and all work registers
  are checked, including nonblank output registers. This is a constructive
  all-size analytic argument with a bounded executable reference suite.
* **Writer identity:** the cached-version control agrees on every initial
  input and fails isolated intermediate interventions. Initial equivalence
  is not treated as evidence of correct writer ancestry.
* **Capacity:** the bound uses a complete finite transcript and independent
  interventions while caches are fixed. Timing side channels are counted.
  Incoming bandwidth excludes the receiver's own locally available record;
  the audited bound is c*q^2/(2L), not c*q^2/L. Aggregate read incidences are
  not unicast transmissions. Native event lower bounds require an additional
  per-event capacity bound. A logical receiver needs a growing composite
  register bank; one C^12 central record has only twelve alternatives.

## Verification coverage

The reference receipt recomputes all 263 circuits: 442,496 forward, inverse,
intervention and boundary-comparison gate executions, 7,880 isolated
intermediate interventions and 40 paired lowering probes,
all 1,728 basis units of the copy channel, exact clock/entropy calculations,
the complete 64-history Toffoli state problem and integer resource controls.
It stores hashes and counts, not the raw tapes. The checker does not import
the compiler. Its feasibility predicate does not test the desired truth table.

Hostile tests reject malformed types, out-of-range/duplicate operands,
incomplete tables, dropped/retargeted gates, dirty work, cached substitution,
stale source pins, missing cases, forged receipt quantities and a forged
entropy prior. JSON duplicate keys and nonfinite constants are rejected;
optimized Python cannot bypass the certificate checks.

The twelve Lean statements have a transitive axiom whitelist audit. Negative
controls reject both `sorryAx` and compiler-trust proof shortcuts. The Lean
coverage is algebraic/path/limit/capacity mathematics; the differentiability,
general compiler induction, matrix-state and entropy arguments are analytic
proofs in the derivation, not claimed as fully formalized Lean results.

## Exit assessment

The contract's finite coherence, simplification, countercontrols, consequences
and reproducibility deliverables are supplied. The remaining physical work
is the native source completion and testing of the proposed clock/refinement
law, not a withheld step of the stated conditional proofs. No empirical RG
validation, canonical axiom amendment or complete source-selection theorem
is asserted by this PR.

## Follow-up audit corrections

The reference receipt previously checked each producer-supplied truth table
but did not independently enforce the meaning and completeness of its named
catalog. The checker now specifies all 263 target functions independently,
rejects duplicates/missing names and rejects coordinated circuit/table
substitution during receipt recomputation. The retained-copy entry also
checks its named intermediate writer against a separate boundary contract;
initial truth-table agreement cannot bypass this check. The additional
comparison gates and probes are counted.

The copy-channel validator accepts actual submitted sparse Kraus operators
and rejects dropped, duplicated, source-erasing and wrong-buffer operators.
An independent producer remains responsible for constructing the candidate.

The tracial-minimum proof requires local references to be restrictions of
the same joint reference. Faithfulness and traciality without compatibility
do not suffice: an exact two-reference counterexample has nonzero entropy
derivative at the proposed joint minimum. The paper and assumption dictionary
state the compatible-reference premise. The latency argument states its
nonnegative grading, zero-displacement waits and inheritance of quadratic
symmetry explicitly. Inner-menu resource bounds are distinguished from the
all-n bound for the full-radius experiment.

CI exposed two inherited integration failures. The source-operation receipt
still pinned the old `ChoiCPTP.lean` comment bytes after main's documentation
correction. A full recomputation changes only that source hash; the native
derivation and every other receipt field remain identical. The Windows CSV
fixture translated embedded LF text to CRLF, so the strict registry-scope
comparison correctly rejected its changed content. The fixture now preserves
bytes and tests both LF and CRLF cells; the production validator is unchanged.

Local follow-up validation passes 232 affected tests on Windows and 167 on
Linux, including 67 RG tests. Both the record-gluing and source-operation
Lean audit targets build. Registry, axiom consistency, reader/Lean style,
CI-budget and receipt-portability gates pass. The changed paper rebuilds
with no overfull boxes or reference/citation/glyph/font problems; its preview
manifest validates and appendix pages 183--186 were visually checked.
