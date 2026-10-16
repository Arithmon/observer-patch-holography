# Auditing finite numerical evidence

This procedure applies to numerical changes under the [standing evidence
audit, #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).
A passing test suite records execution on one tree. It does not establish
the accuracy, domain or scientific interpretation of the calculation.

## Lessons from the October 2026 reviews

The maintainer found defects after the author's audits and passing CI:

| Review | Missed failure | Required independent check |
| --- | --- | --- |
| [#1048](https://github.com/FloatingPragma/observer-patch-holography/pull/1048#pullrequestreview-5423083959) | Subnormal Gibbs populations were rounded before conversion to exact fractions. Variance and entropy derivative were both wrong by 75% while their ratio was correct. | Compute each response component from the original supplied parameters using a separate scalar partition function. Cross population precision boundaries at ordinary and rescaled observable units. |
| [#1049](https://github.com/FloatingPragma/observer-patch-holography/pull/1049#pullrequestreview-5423084144) | Structural symbolic equality rejected an ordinary full-rank complex Markov state whose exact support residual was zero. Real examples and many existing tests missed it. | Compare operators in a canonical exact domain. Exercise complex conditional products, singular supports, nearby non-Markov states and genuine support escape through the public consumers. |
| [#1052](https://github.com/FloatingPragma/observer-patch-holography/pull/1052#pullrequestreview-5439504720) | Correct Petz support refusal aborted a valid count-to-CMI run before its acquired evidence was saved. Default circuit selection hid the separate seed-1 failure. | Save acquired evidence before analysis; test recognized partial results and unexpected errors through serialization, including a nondefault valid circuit. |

The retained reproductions are
[`test_gibbs_population_response.py`](../../code/quantum_information/test_gibbs_population_response.py)
in #1048 and
[`test_complex_markov_support.py`](https://github.com/MarioPoneder/observer-patch-holography/blob/ca79333e2420ee0d4c04dfc366ba87257112a0e8/code/quantum_information/test_complex_markov_support.py)
in #1049. The first repair explicitly refuses insufficient population precision; it
does not claim to evaluate that regime. The second preserves valid complex
states and rejects actual support violations.

## Before changing the implementation

1. Record the reviewed head, exact inputs, observed result and expected
   result or justified refusal. Retain a regression that fails against the
   defective implementation for the intended reason, rather than an import,
   timeout or unavailable dependency.
2. State the numerical contract: units, normalization, tensor order,
   allowed rank/support, output range and any accepted roundoff. Identify
   every conversion, centering, eigensolve, reduction and normalization
   between supplied data and the reported observable.
3. Derive an independent control before reading its answer from the
   producer. Preserve the original input values: high precision applied to
   a damaged intermediate cannot validate its accuracy. An oracle
   must not reuse the producer's thermal states, support decisions or
   cancellation-prone formula when those are the operations under test.

## Challenge the domain and the formula

Select cases that apply to the changed calculation and explain omissions.
Keep the cases and assertions reviewable; a large random sample alone is
not a domain argument.

- **Components and identities:** independently check each component of a
  ratio, difference or normalized diagnostic. Correlated errors can leave
  identities, positivity and normalization intact.
  Derive retained-statistic null denominators from the original population
  and estimator weights. Check that responses, powers and schedule noise can
  describe one common population. Propagate subtraction uncertainty before a
  square root; a small erased residual can have a much larger effect on a
  joint bound. Exercise valid zero-noise and negative unbiased-common-power
  controls through the public consumers as well as their numerical helpers.
  For nearly identical responses, independently centered variance and variance
  about a rounded mean differ by a quadratic centering term. Bound that term
  from the mean's arithmetic error; a tolerance on the full response scale
  can hide a fabricated spread. Retain controls for both accurate and rounded
  centering, erased spread, inflated spread and inconsistent standard errors.
  For a transfer logarithm, test the full Hamiltonian as well as its gap,
  and every relevant stationary-population ratio as well as normalization.
  Resolving small transfer eigenvalues and resolving the leading eigenspace
  are different obligations. Where the source supplies a positive factor,
  preserve that structure before forming a rounded Gram matrix; a successful
  eigensolver on the damaged matrix cannot recover discarded information.
  For thermal responses, test small splittings within the populated energy
  levels as well as the full spectral width. Subtracting a common energy
  origin in finite precision can merge a doublet and destroy its response
  while every population and the derivative ratio pass their checks.
  Preserve the original entries through centering and spectral evaluation;
  challenge each returned state component and derivative separately.
  For operator bases, `Q*Q=I` checks orthonormality, not preservation of
  the supplied span. Check its projector and observable action against the
  original entries, including exactly representable ill-conditioned shears.
  A failed chart is not an intrinsically unresolved span: compare alternate
  exact presentations before attributing an intermediate underflow to the
  source. Check small retained projector entries on their own scale.
  For repeated near-identity operations, check the zero-duration limit and
  accumulated observables: a small unitary phase error can repeat at every
  step. Cancel inverse words only where their source identity is exact;
  retain other phase remainders, raw contractions and correction diagnostics.
  An affine inverse-word evaluation such as `I+B(D-I)C` can suppress phase
  cancellation error while introducing norm drift when the rounded `C` is
  not exactly `B^-1`. Passing a single-operation target-residual check does
  not guarantee probability conservation after repetition. Test complete
  process mass and retained output laws at both coarse and fine time steps,
  including every retained refinement cutoff. When exact source unitarity
  justifies an operation-level polar representation, retain the raw and
  corrected phase-sensitive target checks, fixed polar guards and a bound
  on the total representation change. Do not normalize evolving states or
  selected outcomes to conceal the defect. Include nearby unitary inverse
  perturbations that expose accumulated norm drift, and verify that wrong
  unitary phases fail after restoration.
  A trace-preserving channel can have the wrong action. Check it on
  independent input operators and its reference on each nonzero component's
  own scale. An exactly positive source can have both a computed zero
  eigenvalue and a rounded zero kernel residual; neither determines rank.
  Correctly rounded matrix-function components need not remain accurate when
  multiplied into a channel. Retain a valid neighboring tiny-output control
  alongside any precision refusal. Check the final output against the source
  formula before rounding intermediate factors, and form observable differences
  before rounding the separate expectations. A correct rounded intermediate
  does not guarantee a correct rounded final result.
- **Combined conditioning:** challenge metric whitening and constraint
  inversion together. Separate precision gates can each pass while their
  compounded error spoils the minimum. Compare with the original-input
  solution, and check reported residuals on their own scale so a false zero
  cannot hide inside the much larger target's roundoff budget.
  Test the action of a validated map as well as its matrix identities.
  Measure residuals before rounding a projection: a correctly rounded
  projected operator can give an incorrect distance after subtraction.
  Pair a large conserved component with a tiny transverse perturbation.
- **Stationarity and graph structure:** a small stationary residual need not
  bound population error in a slowly mixing chain. Check an independently
  solved law, periodic chains and multiple closed classes. Decide support,
  reversibility and lumpability from the supplied masses without an absolute
  edge floor; check that a rounded export has not erased a positive edge.
- **Intermediate and output scales:** cross normal/subnormal/zero
  boundaries before and after weighting. Use ordinary units, very small
  and large units, both signs where allowed, and cancelling energy origins.
  A normal output need not have resolved intermediate data.
  A final float conversion can double-round an exact value at a subnormal
  midpoint. Check original-value inequalities against adjacent floats and
  exercise ties in both directions, including the zero boundary.
  Equal results at two precisions can both miss the higher-order term that
  selects a midpoint's side. Compare midpoint distance with numerical
  resolution before accepting the conversion; retain nearby successful
  subnormals. Test cancellation against accepted trace offsets separately
  from normalized-state controls, including unequal inputs with exact zero.
- **Representations:** include real and complex inputs, exact dyadic basis
  changes, repeated eigenvalues, full-rank and singular states, and
  one-dimensional tensor factors. Check algebraic equality in a canonical
  domain, not by expression-tree equality. Check original scalars in mixed
  numeric containers before a common array dtype erases integer precision,
  Boolean types or missingness; validating the coerced array is too late.
  Array subclasses can override `copy`, `nonzero` and indexing while retaining
  ordinary numeric entries. A matrix certificate must reject such subclasses
  or detach them to a trusted base representation before inspecting support
  or entries. Test an off-support perturbation hidden by method dispatch.
  Compare against the exact domain's own zero; equality with a Python
  integer is not guaranteed by a canonical representation. Exercise a
  successful case of each exact fallback as well as its refusals.
  A rounded zero eigenvalue and a rounded zero matrix/kernel product can
  coexist for an exactly full-rank matrix. Check source support before
  reporting infinite divergence, with a finite determinant-based control
  and genuine complex singular containment and escape controls.
- **Arithmetic context:** exercise exact Decimal inputs longer than the
  working precision, inherited rounding modes, exponent limits and traps.
  Deterministic receipt builders must isolate their arithmetic policy and
  preserve the caller's context, including flags. More precision applied after
  a rounded subtraction or coercion cannot recover the original input.
- **Both sides of a decision:** pair an exact zero with a nearby nonzero
  input; valid support with actual support escape; complete evidence with
  incomplete, empty and malformed evidence. A blanket rejection is not a
  successful repair of a valid-input regression. A blanket acceptance or
  zero floor must fail the corresponding negative controls.
- **Shared uncertainty:** when several observations use one unknown calibration,
  state or parameter, require a single admissible witness for the complete
  population. Per-row or pairwise compatibility does not establish joint
  feasibility. Exercise contradictions across controls, executions and held-out
  blocks, alongside feasible boundary cases. Keep a consistency check separate
  from fitting: observations must not silently narrow their own frozen bounds.
- **Callers:** follow the public entry point through its wrappers and
  downstream diagnostics. Test the internal helper and the returned
  observable or classification. Exercise serialization and replay where
  the result is retained as evidence. Inject a failure after acquisition
  and verify that raw evidence survives. A recognized unavailable diagnostic
  must not discard independent valid results or become a successful check;
  unrelated errors must propagate. Exercise reused output locations so
  replacement evidence cannot be paired with a stale success report.
  Passing original inputs to a repaired scalar routine does not repair a
  wrapper's separate pairings: exercise a small linear component alongside
  an ordinary quadratic response, as well as a completely small receipt.
  Reject nonfinite diagnostics at the scientific acceptance gate itself;
  a later JSON serialization failure does not make an earlier `all_pass`
  result valid. For continuing quantum instruments, include wrong-but-CPTP
  controls (wrong tensor order, feedback or hidden-draw correlation) alongside
  trace-preservation checks. Check complete reference/record allocation before
  executing source services, and test mutable checkpoint envelopes as well
  as their array buffers.
  Round-trip actual producer outputs through the strict JSON comparator;
  a CLI test returning a serialized fixture misses runtime type
  mismatches such as tuples becoming lists.
  Compare freshly generated source-binding maps with retained receipt maps
  on each supported platform. Repository-relative keys use `/`; checking
  existing file hashes alone misses platform-dependent path separators.
  Exercise validation with cold and warm caches. Python cache keys can equate
  Boolean, integer, floating and exact rational values. Validation inside a
  cached function may therefore be skipped for an invalid argument after a
  valid equivalent key was stored. Validate original arguments before entering
  the cache; a passing cold-input refusal does not establish this boundary.
- **Numerical versus exact claims:** repeated-precision agreement is not
  an interval bound or proof of a zero. An exact certificate must establish
  the property of the supplied matrix, without silently projecting it onto
  a preferred state family.

## Verify the repair and its integration

Run the retained failure, independent controls, relevant callers and
affected suites with warnings treated as errors. Restore the defective
operation in an isolated copy: the appropriate scientific assertion must
fail. Choose mutations from plausible failure mechanisms, including
constant success, blanket refusal and missing validation when applicable;
report surviving mutations and controls added because of them. A mutation
count is not a completeness claim.

Before reporting readiness, fetch current main and inspect every published
review and conversation thread on the affected PRs. Recheck the exact final
heads, including any fixes made after approval. Combine overlapping open
PRs in a temporary tree and test shared consumers; do not infer integration
from separate green branches. Record the source commits and test commands.

For each changed consumer, state whether numerical results, scientific
classification, live receipts, frozen evidence, registry payloads or public
claims change. Regenerate live artifacts with their canonical producers.
Never refresh frozen evidence simply to make a hash check pass. Resolve a
generated-PDF merge conflict by rebuilding from the merged source and
current release metadata, then regenerate and validate the manifest and
visually inspect the changed pages; choosing either old binary is not enough.

Claim-registry edits also affect receipts that bind complete downstream claim
rows. In particular, `code/m1_necessity` follows the registered dependency
closure even when none of its numerical code changes. Regenerate its live
receipt with `PYTHONPATH=code python -m m1_necessity.build`, inspect the exact
field changes, and run its producer-free verifier and hostile controls.
Confirm that any custody-only update leaves the numerical evidence unchanged;
a passing model manifest or registry validator does not check that consumer.

## Reporting and exit

Report each maintainer finding with its corrective commit and retained
control. Distinguish a code correction from formal maintainer approval;
unresolved review status remains visible until the maintainer revisits it.
Separate tests actually run on the final tree from historical results.
Check all triggered CI on the exact pushed head when a complete audit is
requested, and state expected skips separately. Recheck mergeability after
main moves. Avoid saying an audit is complete merely because CI is green,
and never describe numerical repair as evidence for a new physical law.
