# A native reference state supplies a modular frequency

This packet connects the existing native event statistics to a programmed,
observer-relative modular readout, including a source-compiled flow on the
full native `M66` algebra. The faithful Palm reference is
`sigma = diag(degree/600)`, derived in
[`event_circulation`](../event_circulation/README.md). Its two eigenvalues have
ratio `5/4`. The reference therefore supplies a nonzero modular frequency
`log(5/4)` without choosing a Gibbs state from an externally supplied
Hamiltonian. This is a finite reference-state attachment, not a recovery of
interacting spacetime or a new comparison with experimental data.

For the active-corner Ramsey experiment below, the reference is **not the
probe preparation**. That probe begins at
`I66/66`, uses the existing fixed-seam signed-current herald, and prepares a
coherent mixture with the compiled native balanced mixer. Each herald sign
has probability `5/66`; failed preparation has probability `28/33` and remains
in the raw record. The ten orientation pairs have degree census
`(8,10):2`, `(10,8):2`, `(8,8):2`, `(10,10):4`. Thus four pairs see a
nontrivial relative modular phase and six do not.

## Exact relation and executed computation

With Heisenberg convention `sigma^(is) A sigma^(-is)`, the dual state action
uses `sigma^(-is)`. On the full native algebra,

```
-log(sigma) = log(60) I + log(5/4) P_degree8.
```

The twenty-dimensional active corner commutes with the reference. There the
implemented diagonal target is `D(theta)=exp(i theta P_degree8)`, up to the
common scalar phase, with `theta=s log(5/4)`. Normalizing the restricted
reference changes only that scalar. Its full-reference weight is `23/75`.
That experiment's extension assigns zero phase to the remaining native addresses;
it omits 22 other degree-eight states. **It does not implement the modular
flow on the full 66-dimensional native algebra.**

Each observer program executes:

1. Native signed-current herald, retaining failure and both signs.
2. On success, the 34-pulse compiled balanced mixer `B`.
3. A fresh compiled diagonal source word `D` on seven work bits, lifted by
   two spectator bits to all 512 encoded and padding coordinates.
4. The same 34-pulse mixer `B`, then the native pointer with sign-aware
   decoding. Failure receives its explicit read record and no phase cost.

The diagonal compiler acts on eight declared degree-eight addresses among
the twenty active addresses. All source matrices are raw binary pulse
products, with no polar projection. Padding remains in the complete runtime;
native herald padding is failure and pointer padding has an unbiased bit.
There are six joint herald/read records, including the zero-probability
failed-preparation/read-one record. The C and E candidate laws are both run.
The reference preparation itself reuses the separately executed Palm-herald
construction; this packet does not rerun that experiment as a probe.

All settings consume the actual `runtime.Program` events, controls and
Kraus arrays with a factor executor: the 66 initial basis columns retain
their original binary weight `1/66`, and every hidden Kraus alternative
adds a column. Final densities are weighted sums of outer products on all
512 coordinates. This avoids subtractive cancellation at exactly zero
ideal probabilities without discarding small amplitudes, projecting states
or normalizing outcomes. At phase `pi/2`, the generic dense observer runtime
independently reproduces the full output and serialized checkpoint
continuation. Its relative positivity guard can reject cancellation noise
at nominally zero outcomes; it is not relaxed for this packet.

The exact ideal conditional sign-match law is

```
P(match | successful herald, s) = 4/5 + cos(s log(5/4))/5.
```

| Programmed phase | Modular parameter | Conditional match | Joint match, each sign |
|---|---|---|---|
| 0 | 0 | 1 | 5/66 |
| pi/2 | pi/(2 log(5/4)) | 4/5 | 2/33 |
| pi | pi/log(5/4) | 3/5 | 1/22 |
| 2 pi | 2 pi/log(5/4) | 1 | 5/66 |

These settings use exact rational-complex target entries `1`, `i` and `-1`.
No binary logarithm or exponential is treated as an exact source target.
The displayed decimal frequency is diagnostic only. Uniform reference
`I66/66` gives identity modular action and conditional match one at every
setting. Conversely, the native unread dissipative attempt has a different,
law-dependent contraction. Programmed modular control is not asserted to be
that natural waiting evolution.

## Error boundary and operational meaning

The certificate converts every original binary component of the local
128-by-128 phase map to an exact rational, compares every coordinate against
the exact diagonal target, and encloses its Frobenius error `epsilon`.
Consequently its operator error is at most `epsilon`, its half-diamond error
at most `delta=epsilon(2+epsilon)/2`, and its channel norm at most
`g=(1+epsilon)^2`. The spectator lift and arbitrary untouched quantum
references preserve these operator bounds. The existing eight-dimensional
mixer certificate is reused. For the whole raw `B D B` sequence the bound is

```
delta_B (1 + g_D g_B) + delta_D g_B.
```

Both mixer uses and numerical amplification are charged. Exact ideal native
instruments and adaptive failed branches contract this comparison. This
bound covers the compiled gates against their stated complete extensions;
it does **not** certify binary native Kraus production, subsequent floating
state propagation, or hardware service errors. Full terminal density and
probability comparisons and checkpoint replay are separately labeled
numerical diagnostics. Physical service discrepancy remains unbounded.
For portable replay the fresh exact phase comparison must pass the fixed
operator-error envelope `1e-10`; the mixer must pass `1e-11`. The exported
whole-sequence half-diamond allowance is outward rounded from those fixed
envelopes, rather than depending on platform-specific last bits.

The law is an ensemble clock candidate on the monotone window
`0 <= s log(5/4) <= pi`; one bit is not an exact clock reading, and sign and
period aliases remain. A noncommuting prepared probe is necessary: the
reference itself is stationary under its modular flow. Coherent services,
the probe design and selected phase settings remain supplied. Calling `s`
an observer's candidate time needs no universal clock or conversion to SI
seconds. What remains is to establish consistency with continuing clock
records, other accessible correlations and interventions, and effective
field dynamics. This packet supplies a concrete reference-derived frequency
and testable finite readout law, rather than that larger identification.

The finite modular-inner attachment is the existing
[`PrivateInner.lean`](../../Lean/Dynamics/PrivateInner.lean); the current
construction uses explicit finite spectral calculus. The prior native
interface supplies the exact herald and pointer laws, while
[`compiled_native_interface`](../compiled_native_interface/README.md) supplies
the complete source mixer and its numerical comparison theorem. No additional
axiom or newly claimed abstract theorem is needed.

## Finite preparation of the faithful reference

[`preparation.py`](preparation.py) removes exact equilibrium preparation as
a premise for this reference approximation. For any initial binary-sector
state, including one entangled with an untouched reference, the existing
[native warm-up theorem](../record_laws/README.md) gives
`D(rho_nR, I66/66 tensor rho_R) <= delta = min(1, 33(19079/20000)^n)`.
The preparation tape and controller are unavailable to subsequent work;
the bound retains the theorem's exogenous seam-coverage hypotheses.

One subsequent **global random-seam movement herald**, the same CP success
map for C and E, sends `I66/66` to `q sigma`, with `q=5/33`. The complete
success/failure instrument gives `|p-q|<=delta` and conditional work-reference
distance at most `delta/q`. Its source effect also bounds `2/15<=p<=1/6`.
Failure remains part of the experiment and its resource cost. This full
Palm herald is distinct from the fixed-seam probe herald, whose successful
support has dimension 20 and cannot prepare a faithful full-algebra state.

Put `t=66 delta`. Positivity gives the operator sandwich
`(1-t)sigma <= M(rho_n)/q <= (1+t)sigma`. For `t<1`, the actual successful
work state `sigma_n=M(rho_n)/p` is faithful, with minimum eigenvalue at least
`(1-t)q/[75 min(1/6,q+delta)]`. Operator monotonicity of the logarithm and
Duhamel's identity then give

```
||log(sigma_n)-log(sigma)+log(p/q) I|| <= t/(1-t),
half-diamond(modular_n(s), modular(s)) <= min(1, |s| t/(1-t)).
```

The scalar phase cancels; no commutation assumption is needed. At **369
native warm-up attempts plus one herald**, the eigenvalue floor exceeds
`0.01333244` and the channel error is below `0.000971` uniformly for
`|s|<=16`, including arbitrary probe references. Repeating the complete
preparation on failure costs at most `2442.015` expected native attempts.
Faithfulness concerns the successful work algebra, not a possibly singular
external reference or every classical record block. The explicit replay
also retains exact finite diagonal histories and failure masses.

Each such state admits the existing finite inner-flow theorem and
[`gibbsState_kms`](../../Lean/QFT/StructuralInheritance.lean) at
`H=-log(sigma_n), beta=1`; its Gibbs Heisenberg parameter is `-s` in the
modular convention used here. These are finite algebraic conclusions.
The full-support construction below approximates the actual `sigma_n` flow
without computing its unknown logarithm. Passive repair is not identified
with modular evolution, and neither attempt counts nor `s` are calibrated
to physical seconds.

## Full native modular flow from the same source services

[`full_flow.py`](full_flow.py) compiles
`U(s)=exp(i s log(5/4) P_degree8)` on **all 30 degree-eight native addresses**,
including the 22 omitted by the active-corner experiment. This is the dual
Palm modular action on every native state and coherence, up to the common
scalar `exp(i s log(60))`. The declared extension acts as identity on all
446 nonnative/padding addresses; it supplies no faithful modular state there.

The complete nine-bit Walsh expansion is verified at every one of the 512
addresses. Closed Gray parity walks execute it using 512 existing `I`/`Z`
phase requests and 510 existing phase-exact CNOTs: **2042 source pulses**.
Each CNOT uses the existing `H J H` word from the same proper-code source
controls. No new global phase primitive is required. Intermediate words
retain every encoded address: projecting onto the native 66-dimensional
shell during execution loses amplitude and is explicitly rejected by a
control. The certificate covers the complete terminal operation, not an
experiment that measures a partial word.

For any requested exact rational `|s|<=16`, rational series enclose
`log(5/4)` and its complex phase. Every entry of the actual raw 512-by-512
source map is compared to this target; the Frobenius enclosure must pass a
portable operator-error gate `epsilon=1e-9`. Thus the half-diamond bound is
`delta_word=epsilon(2+epsilon)/2`, with arbitrary untouched references.
The retained executions use `s=0,1,16`; other requested words must pass their
own numerical source certificate. The actual prepared-reference comparison is

```
half-diamond(source(s), modular_sigma_n(s))
    <= delta_word + min(1, |s| 66 delta_n / (1-66 delta_n)).
```

At 369 native warm-up attempts and one full movement herald this is below
`0.000971` throughout the declared modular window, for each certified word.
The complete preparation-and-source experiment, retaining failure, instead
has half-trace error at most `delta_n + p_max delta_word` against its ideal
Palm experiment. These bounds concern different comparisons; successful
preparation is never treated as free or certain.

The Ramsey witness uses native addresses 16 and 144, outside the active
corner, and predicts `(1+cos(s log(5/4)))/2`. Its initialization and readout
are complete native protocols:

1. A [current-only controller](../native_clock/configuration.py) resets the
   work to code 16. A separate native auxiliary is reset to the
   [pointer's ready configuration](../native_interface/pointer.py).
   Each reset permits arbitrary total-two inputs, including collisions;
   both public timeout flags remain. Private reset tapes are discarded.
2. Two existing source `Y` pulses surround the 2042-pulse modular word.
3. Six source pulses couple work bit 7 to the auxiliary orientation. One
   native coin/current repair reads it. Tracing the auxiliary gives the
   continuing work-bit Lüders instrument on all 512 work coordinates.

The maximum cost is **5357 native attempts, 2050 source pulses and 4100
proper-code service events**. The native work and pointer preparation needs
no external blank or random-seam service. With four seed scans and 32 retries per
mixed edge, the two complete reset errors sum to `140/2^32`. The ideal
coupling/readout theorem covers arbitrary work/reference inputs. Its raw
six-pulse map has operator error at most `1e-11`, half-diamond error
`delta_pointer` and norm bound `g_pointer=(1+1e-11)^2`. Thus

```
complete history error <= 140/2^32 + g_pointer*delta_Ramsey + delta_pointer
                       < 3.4e-8.
```

The norm factor accounts for possible amplification by a raw numerical
map. This compares the complete history with the ideal Palm Ramsey
experiment, including empty ideal timeout blocks. It is distinct from the
`0.000971` comparison to an unknown finitely prepared reference above.
The retained execution averages the exact native reset laws over independent
`I78/78` inputs; the uniform bound does not require that input or independence.
On a correlated input, successful reset can steer the reference. The proof
uses the complete failure branch and never assumes otherwise.

All four public outcomes are retained for both C and E: work-reset timeout,
auxiliary-reset timeout, bit zero and bit one. At `s=16`, the successful
match law is about `0.04524740944`. Padding is retained in both bit outcomes;
nonnative weight is a diagnostic, not a measured membership flag. Raw
outputs are not renormalized. The independent source-array, rank-one reset
and complete pointer-instrument checks establish composition on this common
finite apparatus.

Coherent proper-code load/export, their service ancillas and routing remain
architecture services; the native resets do not construct those service blanks.
Their physical errors and durations are not authenticated by source-array
arithmetic. Committed native repair from a joint diagonal accessible input
with fresh independent diagonal ancillas cannot generate the needed coherent
controls. No passive-repair/modular-time or
laboratory-clock identification follows from implementing this experiment.

## Replay

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference --check model/native_modular_reference/results.json
python -m model.native_modular_reference.preparation --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.full_flow --check model/native_modular_reference/full_flow_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/native_modular_reference
```

Controls cover the faithful reference and 4+6 pair census, omitted native
coordinates, all padding columns, inverse/common phase, scalar amplification
that defeats an unweighted telescope, retained preparation failures, uniform
reference, native waiting, checkpoint continuation and hostile ndarray
subclasses. The retained receipt is software evidence, not measured data.
