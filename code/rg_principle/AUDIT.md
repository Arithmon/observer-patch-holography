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

The reference receipt recomputes all 263 circuits: 442,320 forward, inverse
and intervention gate executions, 7,880 isolated intermediate interventions,
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
