# Native reusable accumulator contract

## Objective

Replace the repeatedly supplied `START = 1`, zero reset and host `ADD`
primitives of the conditional M1 recurrence with an actual scalar-mean word.
For every finite supplied list of input requests, the result is
`1 + sum(requested input values)`. Repeated reads preserve the inputs'
decoded values. The previous mutable accumulator is retired at the next
start. Its old value is not an immutable archive.

This is an implementation of a conditional arithmetic interface. It does
not select metric neighbours, generate requests from A1--A3, establish the
M1 scale laws, or promote M1 from its retained status.

## Deliverables

1. `SourceNativeAccumulator.lean` proves scale alignment, addition, native
   rectangle reset, seed reuse and arbitrary finite request composition.
2. `SourceAccumulatorProgram.lean` proves whole-program execution, retained
   input semantics, polynomial mean-work and scale bounds, and composition
   with the final native read and local integer publication under an explicit
   whole-word error margin. `SourceAccumulatorAxiomAudit.lean` checks every
   listed theorem transitively against the standard Lean axioms.
3. Independent producer and verifier replay the captured W12 L3 seams,
   including every current writer and rounded mean. Complete deterministic
   regeneration supplies the large tapes; small examples and commitments
   are checked in. All five integer payloads from -2 through 2 and all request
   words of lengths 0 through 7 are covered, plus five 15-episode reuse runs.
4. Adversarial controls reject malformed and semantically forged evidence,
   unsupported means, incorrect scales, stale writers, uncharged work,
   premature publication and promotions of the retained assumptions.

## Exit

The kernel accepts all declarations without admissions or nonstandard
axioms. The independent replay reproduces every committed history, checks
the exact ideal semantics after every operation, and checks every reported
local publication, input retention and cost. Producer regeneration is byte
identical. Tests include actual rounding and signed disturbances. Passing
this exit closes the stated conditional arithmetic interface, not M1 or a
general native compiler for the q=13 and q=21 populations.

## State and operation interface

Logical value `v_i` is represented by the pair
`(b + g*v_i/2^e_i, b - g*v_i/2^e_i)`. Indices 0 and 1 are the mutable
accumulator and workspace; index 2 is the unit seed and index 3 the payload.
The generic Lean theorem permits any finite list of inputs with indices at
least 2. Its abstract star and reset rectangle must have supported seams.
The captured instance supplies the following seven distinct pairs:

| Index | Role | Physical ports |
| --- | --- | --- |
| 0 | accumulator | 0, 9 |
| 1 | workspace | 1, 5 |
| 2 | unit input | 7, 4 |
| 3 | payload input | 8, 11 |
| 4 | bus | 14, 40 |
| 5 | bus | 23, 38 |
| 6 | receiver | 45, 39 |

One preparation writes 14 scalar values. Afterwards only the scalar means
write these records. Scales and request positions are public controller
metadata. They depend on requests and prior use counts, not on payloads.
Reset is four means across the accumulator/workspace rectangle; seed export
adds two. Alignment halves a record through three means, and an aligned
addition uses five. The final four paired copies cost eight more means.
The receiver is on a neighbouring carrier, not a demonstrated long route.

Mean reads and writes cost two each. Every declared checkpoint samples the
accumulator and both retained inputs (six scalar samples); final inspection
samples these three records and the receiver (eight samples). Instrumented
diagnostics are included in those sample counts. Controller arithmetic,
metadata, physical isolation, inactive support, evidence storage and a
physical publication instrument are not included in the native mean bound.

## Precision and local publication

The captured replay uses `Q = 2^256`, `g = (Q-3)/Q`, `b = 2`, tie-to-even
rounding on the grid, and one initial unit seed. Negative, zero and positive
payloads are included; prepared scalar values are nonnegative. The public
error model is preparation at most `1/(4Q)`, error at every mean at most
`5/(8Q)` per coordinate (rounding plus `1/(8Q)` disturbance), and readout
error at most `1/(4Q)`. These are conditional error allowances, not measured
hardware performance. The ordinary evidence executes rounding; separate
controls inject nonzero preparation, update and readout disturbances.

After `W` total means the normalized radius for a record of scale `e` is
`R = (4+5W)*2^e / (8*(Q-3))`. The counter never resets when the accumulator
or workspace is cleared. The local decoder takes only its two samples,
`e`, `W` and a public integer magnitude bound. It publishes exactly when
the resulting closed interval has one feasible integer. It abstains on
empty and ambiguous intervals. Availability requires `2R < 1`; finite
precision is not claimed to support indefinitely many episodes.

For `S` episodes, `N` requests and initial scale bound `M`, the proved
program bounds are `e <= M+S+2N` and
`W <= 6S + N*(8+3*(M+S+2N))`. A final remote read adds four to its decoding
scale and eight means. Choosing finite precision with
`(4+5W)*2^e < 4*(Q-3)` suffices for integer identification. Thus sufficient
precision in bits grows at most linearly in `M+S+N`, with logarithmic work
overhead. These are bounds for this supplied controller, not a derivation
of physical capacity or M1 coefficients.

## Retained source obligations

The input preparation, one-time unit contrast, valid requests, supported
placement, time ordering, error limits, scale metadata, checkpoints and
publication instrument are supplied. Logical input values survive by
changing representation scales; their raw amplitudes are not immutable.
The result removes repeated fresh unit writes and host data addition from
this interface. It does not derive initial preparation, immutable committed
layer storage, general routing, a source-selected temporal grammar, its A3
reference/cover, or the complete M1 family.

The source-selection criterion and cover obstruction established
in PR #919 remain applicable. This construction supplies arithmetic words
for a supplied source-derived controller; existence of the words
does not imply that A1--A3 select them.
