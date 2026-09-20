# Constrained read selection: maintainer-style audit

Audit date: 2026-09-20. Base: #918 at `bafa6378`. The original #777 acceptance
contract and the canonical A1--A3 reference were re-read. The conditional
repair-word exit and maintainer corrections of #663 were reviewed. The
current main at the start was `2b03a95c`, already contained in this branch's
base. No full M1 derivation or administrative closure is claimed.

## Objective and scientific deliverables

The final [contract](CONTRACT.md) covers the general constrained selection
criterion, its observer-cover boundary, its native-read bridge and exact
executable evidence. It preserves the larger source-derivation obligation.

| Deliverable | Evidence | Audit conclusion |
| --- | --- | --- |
| Arbitrary convex constrained selection | `SourceConstrainedSelection.lean`, weighted KL expansion, mixture gap and support theorem | No full-simplex, IID, uniform-reference or reference-feasibility assumption. An attained minimum, convexity, nonnegative coordinates, faithful reference and positive weights are explicit. |
| Guaranteed-read criterion | Zero-event and positive-readout equivalences, `SourceConstrainedRead.lean` | A complete-history KL guarantee forces the entire feasible family to exclude cut-avoiding histories. A coarse cover requires a nonnegative representation of the failure readout. |
| Actual scalar transition bridge | Arbitrary real-valued cut-locality and complete receiver-transcript theorems | Both interventions share schedule, exterior preparation and decoder metadata. Every prefix is covered, including the initial state; no intermediate observation escapes the argument. |
| Exact correlated constrained optimizer | Four-word affine family, KL/Pythagorean identity, group-ratio certificate | The universal minimum is kernel checked. Rational producer/verifier certificates reconstruct it independently; it is not a fitted or sampled optimizer. |
| Actual cover boundary | Same convex three-history family and same faithful compatible reference, two state-determining covers | Both unique minima and their support difference are kernel checked. Cover injectivity does not imply history-atom support preservation. The negative affine offset is explicit. |
| Captured probability specialization | Entire W12 level-three seam census and exact rational enclosures | 46,050 distinct seams, 11 cut edges. No assumed regular degree, stochastic sampling, or binary64 tolerance. The 29% bound belongs only to the stated uniform specialization. |
| Executable controls | 2,186 complete toy histories; four retained captured histories; constrained and cover certificates | All toy words and operations are reconstructed on every verification. Positive transport succeeds. Nonconvex and assumed-success-face controls retain their distinct scopes. |

## Proof review

The support proof is an independently formalized classical I-projection
result; it is not presented as a novel entropy principle. It uses the exact
identity at a newly populated coordinate and chooses a small positive mixing
parameter whose t log t term defeats any finite objective difference. The
mixture is feasible by convexity. There is no optimality oracle or additional
project axiom. The KL expression is Mathlib's `klFun` f-divergence; its
entropy-plus-linear form is proved. Normalized classical observer states
with one positive weight per observer instantiate the weighted local KL use.

Forty-eight theorem declarations are checked transitively using the existing
audit command. Only `propext`, `Classical.choice` and `Quot.sound` are allowed.
The audit's guarded sorry/compiler-trust controls run, and a separate injected
axiom hidden behind an intermediate theorem was rejected. All fourteen
checked Lean/toolchain inputs were byte-compared with the WSL build tree.
The new audit is an unconditional per-change CI root and is also imported by
`Geometry.lean`. Tests reject deleting, commenting out or substituting its
CI root and require every local proof dependency, including the imported
audit-command implementation, to be pinned.

## Inferences rejected or strengthened during review

* A full-support uniform word argument alone would leave correlated and
  constrained selection untouched. The general convex theorem and direct
  native-read bridge remove those restrictions in the finite classical
  complete-history branch.
* Extending that theorem from directly scored atoms to every history
  reconstructed by an injective observer cover would be false. The proved
  same-reference cover example rejects this extension. A nonnegative
  failure-readout representation is now an explicit sufficient condition.
* Comparing only the final receiver value leaves an avoidable observation
  loophole. The final theorem covers the complete prefix transcript and even
  grants the decoder knowledge of the whole word.
* The cut calculation must use the captured graph. Carrier 0 has eleven
  crossing seams, not an assumed twelve. The verifier independently builds
  adjacency and rejects duplicate primitive seams.
* The universal error bound applies to an independent balanced input pair,
  hence to at least one input. A decoder always guessing one payload refutes
  an alleged identical lower bound for both inputs.
* A successful route proves possibility, not schedule selection. Conditioning
  on success changes the feasible family and requires an accepted-record
  policy and cost accounting. The successful-face control explicitly labels
  this as an added constraint.
* The uniform all-horizon obstruction concerns every prescribed finite
  deadline. It does not exclude eventual almost-sure delivery, an unbounded
  stopping policy, high-probability service or another source-derived model.
* A nonconvex feasible menu can select a vertex despite a full-support
  competitor. The exact two-point log-score comparison demonstrates why
  convexity cannot be omitted.

## Independent executable verification

The producer updates exact payload states; the verifier propagates exact
source coefficients and reconstructs consumed values, writers and local
samples. The producer builds an edge set; the verifier builds symmetric
adjacency, checks ranges/duplicates and enumerates each edge once. Toy words
are generated by Cartesian products in one implementation and base-three
index expansion in the other. Probability intervals use Fraction powers in
the producer and integer powers/Euclidean division in the verifier.

The constrained producer constructs group projections. Its verifier instead
checks faithful normalized probabilities, the complete partition, every
group constraint and the exact p/r ratio within each group. Those conditions
give the moment annihilator used by the proved Pythagorean identity. A
feasible but nonoptimal candidate fails the ratio certificate. The coarse
cover's constrained injectivity is checked by an independently expanded
nonzero determinant, while its universal mathematical claim is proved in
Lean. The nonconvex comparison uses exact rational exponentiated KL scores,
not floating-point logarithms.

Tests disable all producer entry points while verification continues, inspect
imports, mutate all seven event-field types at six positions across both programs,
alter source files in real copied checkouts, corrupt proof pins, falsify
probability enclosures and resource counts, omit histories, forge constrained
optimizers and cover certificates, and invoke the actual CLI on resealed
artifacts and malformed JSON. Canonical null/empty objects, duplicate keys,
floats, booleans, nonfinite and noncanonical rationals are rejected. Forged
receipt digests do not replace semantic replay.

The two artifacts total 16,025 bytes and 535 lines. Compact toy commitments
do not bypass execution: every committed history is expanded and checked.
No frozen runner, old proof, evidence receipt, tolerance, paper mathematics
or scientific claim status was modified. The published theorem-count floor
remains correct at 11,000 for 11,070 declarations.

## Remaining M1 obligations

The actual A1-generated history/observable grammar, all A2-visible constraints,
same-source response and endogenous transport, reference/cover selection and
the accepted-record optimizer-to-output map remain unconstructed for M1.
The finite classical controls do not instantiate all those clauses, and the
noncommutative A3 case is not proved here. No theorem in this packet supplies
the complete native q=13/q=21 compiler or its operational refinement.

What is now determined is the selection test to meet: in a positively
failure-visible classical cover, guaranteed reads require the entire feasible
family to exclude the corresponding failure histories. A merely injective
coarse cover instead needs its own source-selected optimizer/readout theorem.
The source must justify that choice; declaring a desired successful-read menu
does not discharge it. M1 and #777/#779/#740 remain open.

## Final local validation

* Windows and Linux: **544 focused tests passed on each**, including all 185
  selection tests and the existing native-update, reusable-bus, encoded-memory,
  Lean-CI-budget and mandatory-suite-wiring regressions.
* Full `Geometry` build and the 48-declaration transitive audit passed.
  The indirect injected-axiom probe was rejected, and the exact built inputs
  matched the worktree.
* Standalone evidence verification, claim registry, public-surface registry,
  theorem-count floor, Lean CI budget and whitespace checks passed.
* Frozen mandatory runner SHA256 remains
  `04bab737fa2d1b7b51c5545242e40c20375eea98d1a41eb0645ba32238360291`.

An intermediate Windows run overlapped the final pinned-file edits and was
discarded: its stale-fixture failures correctly reported source-pin mismatch.
The final full run used a stable tree and passed without weakening those
checks. Future final validation must likewise freeze all pinned inputs first.
