# Captured reusable-bus audit

## Objective and exit review

The objective is a finite-precision reusable record protocol on captured W12
seams. The [contract](CONTRACT.md) fixes the four-level alphabet, two prepared
records, request sequence, grid, cleanup count and assumed error allowances.
It does not permit closure of M1 or transfer a supplied scalar archive into
the native operation grammar by renaming it.

| Deliverable | Evidence | Limit |
| --- | --- | --- |
| Native finite operations and cleanup | `program_native`, `scrub_native`, `scrub_contraction`, `clean_bound`, count and archive-preservation theorems | Code, schedule, isolation and preparation are supplied |
| Distinct-carrier transport | Injective 12-port map; all 1704 means per history checked against the captured seams | One route on the declared capture; no full-family compiler |
| Reusable finite-precision records | 16 complete initial-payload histories, each reading records `[0,1,0]`; all 48 reads decode correctly | No record creation, arbitrary overwrite, immutable scalar value or unlimited reuse |
| Error and precision | Exact ideal coefficient bounds, kernel-checked retention/contrast/rounding lemmas, positive decoding margins | Added disturbance, initialization and comparator bounds are hypotheses |
| Evidence and enforcement | Independent rational/integers replay, complete consumed writers, mutation controls, two-platform workflow and explicit Lean audit target | External audit storage and supplied controller are not source-generated memory |

The bounded exit is met by the retained protocol and its stated error model.
The receiver uses two local samples and declared codebook/use information.
Initial payloads enter the verifier to check results; they are not arguments
to either receiver decoder. The 16 histories exhaust the specified two-record
alphabet. The coefficient certificate independently bounds every real input
amplitude in `[-3/2,3/2]`; successful decoding is asserted for the four levels.

## Mathematical and topology checks

The source and receiver are on different actual carriers, 0 and 3; these
carriers are glued neighbours through seam `5–40`. The positive route crosses
carrier 1. The two intermediate paired cells span
carriers, so their rails cannot be treated as local reset pairs. The program
uses only the verified source-side reset edge, followed by paired forward
sweeps. Removing the physical seam `23–45` makes verification fail.

The seven-mean cleanup matrix has nonnegative bus coefficients and maximum
row sum `7/8`. It has identity rows on the two source records and no source
coefficient in a bus row. The independent full-basis matrix replay matches
the Lean formula. Repeated cleanup bounds the ideal residual by
`(3/2)*(7/8)^80`, without a supplied ancillary refresh or assigned zero.

The retained-source error is relative to the attenuated value `a/2^k`.
The receiver reference adds another factor `1/8` for the three-link bus.
No theorem restores a source amplitude or turns the decoder's scaling into
a native amplifier. Uniform transfer errors are bounded by the exact sum of
absolute coefficient errors, multiplied by the initial amplitude bound.

Integer nearest/even rounding is distinct from exact native pair averaging.
The formal integer bound gives an error of at most `1/(2Q)` per mean and
proves that rounded values lie between their integer inputs. The theorem
`program_readout_error` composes native compilation, initial and per-step
errors, both local readout errors and the ideal transfer estimate. Its
per-step hypothesis can include rounding plus the additional disturbance.
The analytical initial, additional-disturbance and
readout allowances are supplied; the retained execution implements rounding
only. No experimental hardware validation is inferred.

The smallest certified decoding margin exceeds `0.0150217`; this includes
the cleanup residual, all means preceding that read, initial error and both
sampled-rail readout errors through the half-contrast bound. A two-fractional-bit
control produces equal receiver rails on the first read and is not decodable.
Exact decision-midpoint controls likewise provide no unique payload.

## Custody, interventions and resource boundaries

Every tape row contains the two input values and the actual prior writers;
the verifier derives subsequent writers instead of trusting declared parents.
Cleanup, equal-input means and the final bus state are retained. The request
after a different record inherits that record's preparation and mean ancestry.
The ideal reread's coefficient from the other record is nonzero. Thus the
decoded payload intervention result must not be restated as exact raw-value
independence or absence of physical causal interaction. Some rounded states
can become exactly blank without changing this ideal-model distinction.

There are 12 prepared scalar registers and 1704 means per history, including
729 cross-carrier means and all 240 cleanup sweeps. Input/output traffic and
six receiver samples are counted. Scalar storage, intermediate arithmetic,
expanded mean schedule, port mapping and read metadata have explicit encoding
bounds. These are not a complete physical-machine resource model: controller
realization, timing, calibration/codebook representation, framing and external
logging are excluded explicitly. No omitted ancillary register carries payload.
The widths describe the retained integer model, including a signed integer
threshold decoder, rather than Python object storage or arbitrary real noise.
The 12 active registers do not count the spectator ports of the captured
15,360-port support or supply its isolation. Individual rounded means can
change total load by a grid unit; exact native conservation is not claimed
for those approximations.

All addresses and versions refer to two records supplied at initialization.
The protocol derives their finite read/reuse behavior under an admitted word;
it does not produce a population, choose the word, write a new arbitrary
record, or establish that every allowed repair preserves stored information.
This is consistent with the scoped passive-memory obstructions in #914.

## Verification and trust

The builder and verifier share serialization and hashing only. The verifier
constructs its word by event-index arithmetic, rounds via rational distance
and parity, and separately builds ideal linear maps. Altered operations,
inputs, writers, codewords, declared errors, precision, resource counts and
histories are rejected. An identity-valued event cannot be deleted to reduce
cost. Both generated JSON files are compared byte for byte on each platform.

All 25 public Lean theorems have transitive dependency audits. The audit
accepts only `propext`, `Classical.choice` and `Quot.sound`; built-in controls
reject `sorryAx` and compiler-trusted Boolean evaluation. An isolated theorem
depending on an injected axiom was also rejected. The audit is an explicit
Lean CI target, so dependency-only changes cannot omit it. Test-time parsing
checks executable target words rather than accepting commented-out entries.

## Maintainer-style follow-up findings

The follow-up review checked proof hypotheses, the actual captured edges,
operation order, every arithmetic and storage width, receiver inputs and
error propagation. It found the following evidence gaps and corrected them:

- Resource totals in the verifier were matching literals. The totals were
  correct, but this was weaker than an independent ledger. Counts now come
  from the validated tape and topology, widths from the declared encodings,
  and observed arithmetic maxima are retained. A control deletes a crossing
  event and checks the resulting traffic counts; an overflow control rejects
  an oversized receiver intermediate.
- The aggregate receiver bound was explained by combining separate lemmas.
  The new `program_readout_error` checks that composition in Lean, including
  arbitrary signed errors on idle registers and two local readout bounds.
  The finite route coefficients are independently checked rational evidence;
  the Python replay itself is not a Lean-verified program.
- The standalone CLI accepted noncanonical bytes although the CI regeneration
  test rejected them. The CLI now enforces canonical controls and receipts.
  Isolated invocations reject altered proof bytes, tapes, promoted receipts
  and noncanonical bytes after first accepting each untouched copy.
- CI tests now parse the actual workflow selector and both trigger lists,
  cover every pinned input, and reject deleted or commented-out audit targets.

Six exact adversarial trajectories choose disturbances by backward receiver
coefficients, maximizing each read's signed error in both directions. They
check every intermediate coordinate bound, idle-archive drift and final
decoding. These are tests of the abstract allowance, not added production
histories or a physical noise source. A separate control demonstrates that
cleanup preserves a common-mode offset; that error is not contracted away.

The follow-up does not change the initial payloads, 27,264 retained means,
read outcomes, declared resource totals, error budgets or decoding margins.
The receipt gains independently computed arithmetic maxima; source digests
change with the added theorem and verification code.

The prior encoded-memory, feedback-transport and captured-routing JSON files
are byte-identical to base `6eca2c8c`. No existing scientific classification,
paper source, PDF or retained M1 parent result is changed. The frozen mandatory
runner retains SHA-256
`04bab737fa2d1b7b51c5545242e40c20375eea98d1a41eb0645ba32238360291`.

Local validation:

- `lake build Geometry` succeeds (8427 jobs), including the registered bus
  audit. The direct bus audit build succeeds with all 25 declarations and
  no warnings in the added modules. The isolated extra-axiom probe fails
  with the expected rejection.
- All 24 stages of mandatory shard 0 pass on Windows and Linux, including
  the complete Lean comment/axiom-coverage gate and the subsequent register
  checks. The frozen runner is invoked directly. A Linux-only wrapper maps
  Windows worktree Git metadata to Linux paths for historical-source checks;
  it does not change any repository file or acceptance rule.
  Windows reports the existing `custody_parser_unsupported_on_platform`
  result; Linux verifies external custody. Source-history checks pass on both.
- **215 tests pass on each platform**: 97 reusable-bus controls, 64 encoded
  memory controls, and 54 mandatory-workflow/Lean-budget regressions. The
  bus suite includes strict byte regeneration and independent receipt replay.
- The scientific Python suite collects **5348 tests** successfully. This is
  a collection result, not a claim that all 5348 tests were executed.
- `git diff --check` passes. A virtual merge with the green #914 head
  `f317cfb4` succeeds without conflicts.

The remaining M1 obligations include general payload/write/commit services,
source-selected control and preparation, complete q=13/q=21 routing with the
required intervention-preserving refinement, and the corresponding resource
accounting. Physical clock/count and common-world joins are separate. The
bounded result closes none of #777, #779 or #740.
