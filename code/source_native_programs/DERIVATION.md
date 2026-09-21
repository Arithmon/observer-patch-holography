# Finite native stored programs

## Representation and supported layout

A record of integer `v` and public exponent `e` has reference loads
`b ± g v / 2^e`. Its logical value survives source attenuation; raw amplitudes
are mutable. The preparation has `|v| <= P`, a retained unit of value one,
common blank loads elsewhere, `A=max(1,P)`, and baseline `(A+1)Q` in grid units.
The gain is `g=(Q-3)/Q`. All `12C` scalar loads are charged as preparation writes.

The hub carrier uses accumulator `(2,1)`, workspace `(4,9)`, unit `(6,7)`, and
operand `(3,5)` or `(5,3)`. All four operand/workspace rectangle edges exist,
so either incoming polarity supports addition and operand retirement. The
accumulator has no reset rung. Repeating workspace-clear then paired
accumulator/workspace means `k` times, followed by workspace-clear, leaves
accumulator amplitude divided by `2^k`; this residual is charged. Copying the
retained unit starts the accumulator. No unit-valued dynamic write occurs.

Each bank slot occupies ports `(0,1)` of its own carrier. Its export scratch
is `(5,9)`. The local seams support export and both reset rungs. Two banks
require `2N+1 <= C`, including the core. Routes avoid every other bank record
and every live core port. Their scalar paths are simple, disjoint and equal
in length. An incoming polarity is public metadata, not a sign chosen by
examining the payload. The path finder is an untrusted finite search with no
all-level completeness theorem. Its complete witnesses at the three captured
levels are independently validated.

## Shuttle and composition

Write a route as `m=n+2` paired cells, numbered `0,...,n+1`, and keep the source
outside it. Clear a retired destination if needed, copy the source into blank
cell zero, and execute paired means along the route. The source amplitude is
halved and the destination receives `source/2^m`. Cells through `n` form the
bus. Only its first cell needs a rung. Repeat root-clear plus a forward sweep
`(n+1)^2 k` times, excluding the destination. The bound is

```text
bus residual <= A (n+1)^2 / (2n+1) 2^(-k).
read word length  = 4 + 2n + (1+2n)(n+1)^2 k
store word length = 5 + 2n + (1+2n)(n+1)^2 k.
```

The reference state projects only the cleaned bus to the baseline. This is
an error comparison, never a physical reset. A read into a reference-blank
operand omits its nonexistent reset rung. Physical residuals remain in the
global initial discrepancy for the next step. Different routes can pair
otherwise blank ports differently because their reference load is the same.

`SourceNativeShuttle` proves the exact source/destination formulas, outside
record preservation, scalar word lengths and bus error. `SourceNativeCore`
proves restart and its residual. `SourceNativeStoredProgram.Construction`
admits start, aligned addition, rectangle retirement, store, read, scalar-port
permutations and concatenation. `correct`, `noisy_correct` and `publication`
derive the composed certificate by induction over those constructors. There
is no constructor accepting an arbitrary unproved local error certificate.
The lower-level `SourceNativeProgramError` also exposes the stage composition
lemma used to express the global error ledger.

For a captured finite route, its distinct ports and all other live record
pairs form a finite injective assignment. Extend it to a permutation of the
countable scalar index set used by the theorem; unused indices have baseline
loads. Only its finitely supported native word acts. The route checker proves
the assignment and support conditions as finite computations. This scalar
extension/re-pairing argument remains analytic. The general paired-cell bank
compiler below has its own kernel-checked induction; it does not identify
every captured scalar embedding or verify the Python implementation generally.

## Layer and lifetime induction

At layer `t`, the input bank contains the values of layer `t-1`. Every other
bank slot is retired or contains a completed output of layer `t`. For each
site, start at one, shuttle the requested input versions in list order, perform
scale-aligned addition, retire the operand through its supported rectangle,
then store into the opposite bank. A source export changes only its public
exponent, not its decoded value. A store excludes the entire input bank and
all other output records. Thus the input-bank invariant holds throughout the
layer, even when a later site reads an input after earlier sites have stored
their results. Swap banks only after all sites finish. Induction gives

```text
v[t,i] = 1 + sum(v[t-1,j] for j in menu[t-1,i]).
cap[t,i] = 1 + sum(cap[t-1,j] for j in menu[t-1,i]).
```

The cap is public and derived from `P`; it does not reveal the answer.
`certify_plan` reconstructs every edge block, repetition, bank version,
orientation, exponent and cap independently. Native replay uses only scalar
means after preparation. A separate logical oracle checks publications; it
never supplies native values or enters the local decoder.

## Kernel-checked bank compiler

`SourceBankMachine` emits `start`, `transfer`, `add` and `retire` instructions
from finite typed menus. Logical registers 0, 1, 2 and 3 hold the accumulator,
workspace, retained unit and operand. Two disjoint banks begin at register 4.
These are logical pair coordinates, not the captured graph's carrier numbers.
`program_recurrence` derives the iterated recurrence from the instruction
semantics for arbitrary width, signed initial values, list lengths and layer
count. `recurrence_bound` propagates public caps even with signed cancellation.

`SourceBankInvariant` proves clean workspace/operand conditions, protected
input-bank lifetimes and finite support preservation from the compiler.
`SourceBankLowering` proves each instruction's native implementation and
reference-state identity. A route annotation is a cell permutation with
source/destination images and a depth of at least three. Its scratch cells
must be workspace or outside the supported records. It has no correctness,
desired-output, arbitrary-word or error-certificate field.

`SourceBankCompiler.native_construction` derives the whole native word and
error ledger from those geometric annotations. `compiled_publication`
combines native noisy execution with the actual compiled recurrence and its
derived cap. The caller does not supply the answer or prove a completed
`Construction`. `native_word_length` identifies the scalar count exactly with
the computable `work` function, including alignment and all cleanup repeats.

`SourceBankBusWitness` constructs every required route at depth three using
explicit finite permutations. Its read bus uses two scratch cells; its write
bus uses the workspace and one scratch cell. `every_program_native` therefore
has no route-existence premise. This is a declared star bus, not an embedding
into W12 and not a replacement for the captured-route certificates.

`SourceBankExecution.finite_execution` combines this witness with the precision
theorem: for every finite program and public amplitude bound there exist a
cleanup count and grid exponent such that every execution within the stated
preparation, per-mean and sampling allowances publishes the requested layered
recurrence. The choices precede and are independent of the payload. No final
correctness or strict-margin hypothesis remains. Selection and realization of
the controller, prepared records and physical instrument are supplied.

The compact generated `SourceBankControl` checks the captured signed control
against the same general instruction stream. Kernel reduction checks all 42
stages, every captured observed scale, each native stage length, their total
of 1,798,154 means, and the nine signed final output pairs. Six negative
certificates reject an empty trace, missing/extra stage, wrong instruction,
undercharged work and wrong scale. `lean_control.py --check` ties the source
to the complete independently validated plan; changed data cannot silently
leave a stale theorem. This finite check does not generalize itself to q=3,
q=13 or q=21, or prove captured scalar re-pairing correct in Lean.

## Global precision and resources

Reference amplitudes remain bounded by `A`: means preserve the amplitude
interval and comparison projections set amplitudes to zero. Let `J` count
starts and shuttles, `D` bound bus cells, and `E` bound every observed exponent.
Choose

```text
k = E + ceil(log2(64 A J D^2))
R = A J D^2 / 2^k
B = E + ceil(log2(64 (W+1)))
Q = 2^B,
```

where `W` is the complete scalar word length with that `k`. Each start adds at
most `A/2^k` and each shuttle at most `A D^2/2^k`, so `R` bounds their sum.
Preparation and sampling allowances are `1/(4Q)` each, and signed per-mean
disturbances including grid rounding are bounded by `5/(8Q)`, even on idle
ports. The normalized local observation and radius are

```text
z = (plus-minus) 2^e / (2(Q-3))
r = ((4+5 Wprefix)/8 + Q R) 2^e / (Q-3).
```

The decoder returns the unique integer in `[-cap,cap]` intersected with
`[z-r,z+r]`, or returns no value. `SourceNativeProgramBudget.precision_margin`
proves a strict singleton margin for the sufficient bounds above. Neither
the error clock nor `R` resets at a version boundary.

The general Lean existence proof uses deliberately coarse, finite choices.
Let `C=sum(cleanupCharge)` with charge one for a start and `(d-1)^2` for a
transfer, and let `E` be the initial scale bound plus the sum of public
instruction scale charges. `SourceBankCompiler.error_bound` proves
`error <= A C / 2^k`, while `scales_bound` proves every final scale is at most
`E`. Choosing `k=ceil(64 * 2^E * A * C)` and then
`B=ceil(64 * (work(k)+1) * 2^E)` gives the strict whole-word margin through
`SourceBankPrecision.exists_global_margin`. The power-of-two bounds, exact
work dependence on `k`, and nonnegative error sum are proved in Lean. These
existence choices are much larger than the logarithmic Python choices above;
they establish finite conditional execution, not hardware attainability.

For a stationary menu with `T` layers and `N` sites, let `S=TN`, `H` be the
total requested reads, `d_in` the largest number of reads of one input in a
layer, `d_out` the largest request-list length, and `L` the maximum route cells.
The unit exponent is at most `S`. In a layer, input exponents grow by at most
`d_in`, a transfer adds at most `L`, and successive additions increase the
maximum aligned exponent by at most `d_out+1`. Storing adds at most `L`.
Layer induction therefore permits

```text
E <= S + T(d_in+d_out+2L+1)
J = 2S+H, D=L-1
W <= S(3k+3) + (S+H)(2L+1+(2L-3)(L-1)^2 k) + H(3E+12).
```

The last term includes alignment, addition and the four-mean operand
retirement. These are explicit sufficient finite bounds, not measured large
executions or physical efficiency claims. They count means, `2W` mean reads,
`2W` mean writes and `12C` prepared scalar registers. Public exponent storage,
route search, the controller, geometric address data, comparison circuitry and
physical isolation are separate supplied resources.

Precision also costs storage. A sufficient unsigned scalar word has
`B+ceil(log2(2A+2))` bits: the prepared load is at most `(2A+1)Q`, means do not
increase the scalar range, and the global signed grid error is less than `Q`.
Thus scalar storage alone has the sufficient bound
`12C (B+ceil(log2(2A+2)))` bits. Counting registers does not make those bits free.

## Evidence and trust boundary

The finite controls include signed cancellation, duplicate and empty requests,
three-layer bank reuse, complete golden metric menus at q=3, and a centre-site
intervention. Every checkpoint compares local samples with its immutable
logical version and public scale. Stationary-block acceleration checks the
whole touched state after an executed block; an unchanged deterministic state
has identical repeated images. `fixed_point` proves that implication and
`repeated_word_length` counts all repetitions. The receipts separate evaluated
means from stationary means whose effect is certified this way.

At q=13 and q=21 the evidence checks complete two-bank route inventories and
every metric decision, with sufficient finite word/precision bounds. It does
not execute those large native words. No source-selection, observer-cover,
physical clock or continuum conclusion follows from the existence of a
compiled word. The finite controller is supplied, as specified in the contract.

The finite log-recovery result in `evidence/federation_recovery` also applies:
replaying an observed program does not identify unvisited support topology.
Its two indistinguishable complete gluings are an identification boundary,
not a way to infer a W12 embedding from the native histories. The declared
bus witness and the captured-route evidence therefore remain separate.

| Population | Complete route witnesses | Maximum paired cells | Native execution evidence |
| --- | ---: | ---: | --- |
| 27 | 108 | 24 | Two full two-layer histories, baseline and centre intervention |
| 2,197 | 8,788 | 102 | Analytic sufficient bound |
| 9,261 | 37,044 | 197 | Analytic sufficient bound |

Each 27-site native word contains 704,826,613 scalar means and 197,636 scalar
samples. The baseline explicitly evaluates 146,277,951 means and certifies
558,548,662 stationary repetitions. The intervention evaluates 146,279,344 and
certifies 558,547,269. Its `+1` change affects seven first-layer outputs and 26
second-layer outputs; the centre second-layer value changes from 1,416 to 1,423.
The nine signed three-layer controls account for 16,183,386 scalar means.

The large-family sufficient bounds are expensive: 27,466,293,578,568,192 means
and 11,309 grid bits at q=13; 11,563,547,398,990,878,570 means and 52,720 grid
bits at q=21. These upper bounds establish finite compilability under the
supplied premises. They do not establish physical efficiency or feasibility.
