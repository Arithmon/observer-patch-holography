# Self-audit of the passive-memory derivation

Base: main `c0857e3bcd2b529d266c77d263377dd6ffb8aa7d`, fetched before this work.
This branch is independent of the unmerged PR #913. It imports existing main
proofs and uses the already retained W12 capture. No simulator modification,
new source-population capture or q=13/q=21 execution is claimed.

## Contract review

| Deliverable | Evidence | Exit and boundary |
| --- | --- | --- |
| Arbitrary state-only record classification | `SourceNonlinearRecord.lean`, 11 theorems; six retained equal-image witnesses | Factors through total load on connected finite support under invariance on **all** nonnegative states; no linearity/continuity assumption. No claim about restricted codes or observation histories. |
| Reset obstruction and matching finite-error construction | `SourcePassiveReset.lean`, 14 theorems; five star controls, reused-zero control, five captured neighbours | Positive raw loads never become exactly zero in finite words. Fresh zeros attain the lower bound with every used ancilla retained. Preparation and growing star support remain supplied. |
| Closed reuse and balanced-copy budget | `SourcePassiveMemoryBudget.lean`, 15 theorems; complete prepared-copy/cleanup tape | Exact quadratic loss, only trivial closed cycles, and no nonzero catalytic balanced copy on the whole finite register set. The preparation is payload-dependent; no unknown-input copier or physical energy is derived. |
| Kernel audit and independently checked evidence | `SourcePassiveMemoryAxiomAudit.lean`, `verify.py`, retained JSON and mutation suite | All 40 theorems audited transitively, all eight specified executions replayed from disk; finite controls illustrate the theorems rather than proving their universal scope. |

The bounded contract is met. Full M1 and the source-selection/common-world
obligations in #779/#740 remain open. Neither a broad issue closure nor a
physical selection claim is justified by this work.

## Mathematical review

- The record proof moves *readout equalities* across equal-image fibres.
  Each background and moved addend is nonnegative; no inverse mean is inserted
  as an operation. Path direction is legitimate because equality is symmetric.
  Its converse is actual total-load conservation. A total can itself carry
  arbitrarily much ideal-real information; the theorem is not a bit-capacity
  bound or a claim that no memory whatsoever survives.
- The native-zero theorem counts only means incident on the tested register
  and includes changing ancillas. All loads are initially nonnegative and
  remain so. For any particular adaptive execution its realized finite word
  satisfies the theorem. This does not assert a uniform adaptive schedule or
  a bound for another primitive operation. Natural-number indexing in the
  constructive proof has a separate finite-support theorem.
- Sharpness uses supplied fresh zero neighbours; the retained five-neighbour
  W12 subprogram does not establish arbitrarily many neighbours on that fixed
  port. Used positive zeros cannot be refreshed exactly by finite means.
  Allocations and preparation assignments are counted, without claiming to
  derive their physical cost or a universal minimum-space reset algorithm.
- The convexity defect is exact for each mean, including the identical-endpoint
  case in Lean. The global proof works over arbitrary finite real states.
  Equality forces every operation to be trivial, so restoring a subsystem
  while changing its environment is not excluded as a closed-system cycle.
- Balanced copy explicitly includes all ancillary coordinates and has the
  same total before and after. The missing resource is quantified by the
  quadratic budget. Common-baseline covariance permits nonnegative raw-load
  examples. Logical reset to equal nonzero rails is allowed, as the retained
  cleanup demonstrates. The unexecuted catalytic target is labelled as such.
- The existing feedback error identity gives residual-reset error directly:
  exact export/mean/readback with residual `r` returns `payload+r`. The new
  reset construction bounds this one contribution. It does not derive those
  other operations, their preparation, noise, repetition or clock.

## Evidence and prior-audit lessons

The verifier reads the retained JSON, with independent operation semantics
and independent experiment recipes. It checks actual operands, both consumed
writers, port identities, allowed menus, full final states and counts. A
fully valid alternative star order with the same source residual and Q is
rejected as a different experiment. Ancilla deletion, invented zero resets,
missing events, stale writers, free preparation and undercounted resources
are tested. The graph census is independently reconstructed by breadth-first
search rather than reusing the builder's union-find.

All 14 source/proof/support/toolchain pins are recomputed. JSON duplicate keys,
floating/nonfinite numbers, Boolean event indices and noncanonical rational
strings are rejected. The verifier does not import the builder. Unique package
imports avoid the previous `verify`/`build` collision problem in full collection;
the CLI also avoids Python's standard-library `code` module. Windows testing
caught and corrected a missing explicit UTF-8 read in the audit-coverage test.
The initial Linux test attempt needed its temporary parent directory created;
no scientific receipt was changed to fix that environment setup.

The mandatory runner is checked against its existing frozen source-projection
hash. The new suite has dedicated Windows/Linux CI. Existing scientific
receipts, claim registries, paper text and publication artifacts are unchanged.

## Follow-up audit of published head `1d93ef54`

The requested second review reread all three proof modules, the contract,
independent replay, retained records, CI selection and the current #779/#740
boundaries. No invalid mathematical implication was found in the 40 theorems.
Their statements and proofs, all 14 source pins, and both scientific JSON
artifacts are unchanged by the corrections below.

One enforcement gap was found and corrected: Lean CI previously built changed
modules and two fixed targets. If only a passive-memory proof or upstream
dependency changed, an unchanged downstream axiom audit was not necessarily
run. Importing that audit in an unchanged `Geometry.lean` was insufficient.
`Geometry.SourcePassiveMemoryAxiomAudit` is now an explicit default CI target;
the control workflow also runs when this CI wiring changes. A regression guard
checks this obligation. The existing build timeouts and frozen runner remain
unchanged.

The actual old/new Bash selector was executed with controlled changed-file
lists for `SourcePassiveReset.lean`, `ScalarSeamRepair.lean`, and the toolchain.
The old selector omitted the audit in all three cases; the corrected selector
included it in all three. Separately, an isolated temporary Lean module with
an extra axiom compiled normally, while the real audit rejected its dependent
theorem with `Passive memory axiom audit rejected [UntrustedProbe.extra]`.
That deliberately untrusted fixture is outside the repository and is not part
of the mathematical evidence. The clean 40-theorem audit passes.

Four rational-matrix scope tests were added without importing either tape
engine. They demonstrate why signed raw loads, encoded logical zero,
disconnected components and a connected restricted code must not be covered
by stronger readings of the results. The restricted-code test exercises 121
words, alongside the original 121 replay controls. These are test-time
counterexample controls, not new retained production executions.

The README now states a narrow conclusion about the actual record/reset/copy
interfaces rather than suggesting a no-go for every implementation of
feedback. It lists the remaining premises of the current scalar-feedback
route and makes the scalar cost model's exclusions explicit: writer metadata,
controller/schedule state, audit history, graph work and exact-arithmetic bits
are not physical memory resources accounted for by the register count.

## Validation performed

- `lake build Geometry`: passed, 8,427 build jobs; all 40 new theorem dependencies
  are within `propext`, `Classical.choice`, `Quot.sound`. The audit's negative
  controls reject `sorryAx` and compiler-trusted Boolean reduction.
- Windows Python 3.13: 110 tests passed in the final combined run: 58
  passive-memory tests, 34 invariant-mining freeze tests, and 18 affected
  CI-budget/classification regressions.
- Linux Python 3.12: 76 tests passed in the final combined run: the same 58
  passive-memory tests and 18 affected CI regressions. Exact retained byte
  regeneration and the 14 input pins pass on both platforms.
- The original validation also ran the 31 existing feedback-transport tests
  successfully; no feedback implementation changed during this follow-up.
- Repository-wide Python collection: 5,222 tests collected successfully.
  This is an import/collection check, not execution of all 5,222 tests.
- `git diff --check`: passed.
- The existing Lean `native_decide` inventory check passes on Linux with its
  reviewed 13 declarations unchanged. On Windows the same unmodified script
  compares backslash paths against slash-form expected paths and reports a
  mismatch; the counts and modules agree. Lean CI runs this check on Linux.
  This pre-existing platform limitation is not a new proof-trust failure.

These are local validation results. The reviewed mathematical assumptions and
bounded exit are stated above; no broader CI or scientific certification is
inferred from these checks.
