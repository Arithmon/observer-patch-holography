# Native reset, addition and finite reuse

Write the amplitude of record `i` as `a_i = g*v_i/2^e_i`, encoded on the two
rails `b+a_i` and `b-a_i`. A paired copy replaces both source and destination
amplitudes by their mean. The construction from the native-record package
halves a record into a zero workspace and clears that workspace using three
scalar means. Repeating it `k` times divides the amplitude by `2^k`.

## Addition without writing the sum

Let `s` be an input distinct from the accumulator and workspace, and put
`t = max(e_0, e_s+1)`. Halve the accumulator `t-e_0` times and the input
`t-e_s-1` times. Copy the input into the workspace, copy the accumulator
into the workspace, then clear the workspace. The surviving amplitudes are

```
accumulator: g*(v_0+v_s)/2^(t+1)
input s:     g*v_s/2^t
workspace:   0
```

The five means in the final three instructions perform the addition.
Host arithmetic updates public exponents and request positions; it never
writes the data sum. `add_formula` proves the full state equality, including
every untouched input. `accumulation_formula` and `accumulated_sum` compose
the actual words over arbitrary finite request lists, including repeats.

## A native start from a retained unit

The accumulator's two ports need not share a seam. The captured pairs
`(0,9)` and `(1,5)` have the four cross-pair edges `(0,1)`, `(9,5)`, `(0,5)`
and `(9,1)`. The first two means give both pairs the same amplitude; the
last two mean opposite rails. Both amplitudes become zero. This works for
arbitrary initial amplitudes of both records with a common baseline.
`reset_rectangle` proves that four-mean identity.

Two further means export the unit seed to the cleared accumulator. The seed
and accumulator amplitudes both become half the seed's initial amplitude;
the public scales increase accordingly. `start_native` proves the complete
state transformation. No blank or unit record is allocated between episodes.
No noisy-state error budget is refreshed by the ideal clearing identity.

## Composition and local publication

`program_native` identifies execution of the concatenated scalar word with
the encoded program semantics. `program_input_records` preserves each input
value across all episodes. `program_last_sum` identifies the final value as
the original seed plus the requested original inputs, independently of the
number of preceding episodes.

The final bus uses four paired copies, costing eight means. The receiver's
amplitude is the final accumulator amplitude divided by sixteen. Both the
bus word and that attenuation appear in `program_remote_native`.

For initial error `E`, per-mean error `delta`, sample error `rho` and `W`
executed means, the contrast error is at most `E+W*delta+rho`. It includes
every start, reset, alignment, addition and transport mean. Division by the
record gain multiplies the error allowance by the same factor as the
observation. `integer_identification` proves uniqueness when the normalized
radius is strictly below one half; `integer_publication` connects this to the
canonical singleton policy.

`native_recurrence_publication` composes the entire construction. From a unit
seed and arbitrary integer inputs, valid finite requests, the actual noisy
mean word, and the stated strict error margin, the final local policy returns
`1 + sum(requested original inputs)`. The arithmetic answer is a conclusion.
The theorem's hypotheses include the supplied requests and supported-use
contract; they do not include a source-selection derivation.

The executable decoder intersects its closed error interval with the public
integer magnitude range. `candidate_iff` and `decode_is_canonical` identify
its clipped ceiling/floor endpoints with the canonical policy for that exact
candidate set, including empty and ambiguous sets. `recurrence_bound` obtains
the public cap from the sum of input magnitude bounds. The composed
`native_bounded_publication` theorem applies this exact decoder to the noisy
native recurrence; its magnitude cap does not depend on the desired answer.

## Resources and scientific boundary

With `S` episodes, `N` additions and initial scale bound `M`, the actual native
word has at most `6S + N*(8+3*(M+S+2N))` means, before the final eight transport
means. Scales are at most `M+S+2N`, before the receiver's additional four.
The sufficient grid inequality in the contract therefore requires a finite
number of precision bits growing linearly in the finite program size, with
logarithmic work overhead. Neither indefinite fixed-precision reuse nor
physical capacity calibration follows.

The construction closes the specified repeated-start/add/reset arithmetic
interface. Input records preserve logical values by changing amplitudes and
scales. Immutable committed layer archives and a native compiler for the
entire q=13/q=21 family are different interfaces. Selecting the temporal
grammar, feasible words, reference and observer cover from A1--A3 is also
outside the result. The selection package's constrained-KL support criterion
and cover obstruction do not supply these data.
