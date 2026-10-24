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
   native coin/current repair reads it. The ideal reduced target is the
   continuing work-bit Lüders instrument on all 512 work coordinates;
   the literal source word retains its separate numerical error bound.

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

The retained `full_flow_results.json` witness contains the four probability
marginals and trace/leakage diagnostics for the declared `I78/78` inputs.
The complete quantum-output and arbitrary-reference bound follows by
composition of the independently checked reset, source-array and pointer
instruments. The pointer also has separate reference-entangled executions.
Failure quantum channels are specified by the finite native controller;
the full-flow receipt does not enumerate their dense output matrices.
Replay checks the probability and diagnostic domains and the complete
record mass as well as the comparison bound. Portable numerical tolerances
do not permit negative probabilities or relax execution error gates.

Coherent proper-code load/export, their service ancillas and routing remain
architecture services; the native resets do not construct those service blanks.
Their physical errors and durations are not authenticated by source-array
arithmetic. Committed native repair from a joint diagonal accessible input
with fresh independent diagonal ancillas cannot generate the needed coherent
controls. No passive-repair/modular-time or
laboratory-clock identification follows from implementing this experiment.

## Regional algebras under one native modular action

[`regional.py`](regional.py) constructs regional algebras directly on the
native total-two binary shell. For a set of ports `R`, the bare algebra is
`A_R = direct_sum_q M(local q-configurations) tensor I(outside)`. The local
occupation records lie in its diagonal part. Every such classical record
commutes with the Palm reference, so its modular action is trivial. A moving
modular clock needs a coherent probe and the admitted control services.

The full Palm action does not generally preserve `A_R`: the phase on a local
matrix unit depends on the outside configuration. For the adjacent ports
`R=(0,1)`, one coherence has degree ratio `4/5` in one outside context and
`1` in another. Its reduced faithful state is tracial in each charge block,
so its own modular flow is identity. Faithful regional restriction alone
does not establish restriction of the ambient clock.

There is an exact, canonical completion. Group outside configurations when
their positive degree vectors over local configurations are proportional.
Each resulting block has `sigma=p mu tensor nu`, with both factors faithful.
The algebra `C_R = direct_sum_blocks M(local) tensor I(outside context)` is
modular invariant, and its state-preserving expectation is

```
E_R(X) = direct_sum_blocks [(id tensor nu)(P_block X P_block) tensor I].
```

Its Kraus formula gives complete positivity; the displayed factorization
gives unitality, the bimodule identity, idempotence and preservation of
`sigma`. It need not preserve the ordinary matrix trace. The Schrödinger
restriction to all local factors is CPTP and retains every context.

This is the **smallest modular-invariant algebra containing `A_R`**. Spectral
components of a local matrix unit separate its distinct degree ratios;
products with the reverse unit give the corresponding outside-context
projectors. Their joint partition is exactly proportionality of the degree
columns. Minimality proves `C_R subset C_S` whenever `R subset S`; the one
ambient modular group agrees on every common observable. For the cap
`R=(0,1,5,7)`, the algebra dimension grows from 53 to 149. This supplies
finite regional modular compatibility, with explicit additional context.
Its implementation by local records and coherent controls remains an access
obligation. Disjoint completed regions can fail to commute; the receipt
exhibits a nonzero commutator for `(0,1)` and `(2,8)`.

Nested expectations compose: `E_R E_S=E_S E_R=E_R` for `R subset S`.
The bimodule identity and state preservation make each expectation the
orthogonal projection for the faithful inner product `Tr(sigma X*Y)`;
ordinary nested-projection composition then applies.

These explicit finite maps instantiate the relevant algebraic formulas,
not the entire [`CPRegionalNet`](../../Lean/QFT/CPRestrictionNet.lean)
interface: that interface additionally requires trace-preserving expectations,
locality and its other net laws. There is no claimed spacetime causal net.
The existing full source word implements the ideal action simultaneously on
all `C_R`. Regional partial trace preserves the finite-preparation operator
sandwich, giving the same modular-error bound for the actual regional states;
exact invariance of the actual finite-prepared ambient state is not assumed.

## One native entropy completion and its regional clock

[`completed_clock.py`](completed_clock.py) joins the
[native complement-repair entropy completion](../native_equilibration/cut_completion.py)
to a clock on the **same prepared state and bare cap algebra**. For the cap
above, run the admitted native C/E repairs on complementary seams. Their
derived tracial completion `E_R^trace` preserves every cap observable and
selects `bar_sigma=E_R^trace(sigma)`. The source degree counts give its five
eigenvalues `5/336, 19/1200, 37/2400, 1/75, 1/60`. It is the unique entropy
maximum at the retained Palm cap marginal. Its distance from Palm is `8/175`;
the tracial entropy completion and Palm-preserving expectation are distinct.

Now `bar_sigma` factors with a uniform complement in each charge block, so
its modular flow **does preserve the original 53-dimensional cap algebra**.
The same counted state gives its cap modular gaps:

| Cap charge | Nonzero modular frequency magnitude |
| --- | --- |
| 0 | None |
| 1 | `log(38/37)` |
| 2 | `log(5/4)` |

The dimensionless frequency ratio `log(38/37)/log(5/4)` has no fitted rate.
It is an internal prediction of this preparation and access protocol.
No physical temperature, atomic transition or SI time conversion is attached.

Two complete 2042-pulse words use the existing source controls: the direct
state action `U_bar(s)=bar_sigma^(-is)` up to the common native scalar, and
the relative action `V_s=bar_sigma^(-is) sigma^(is)`. The exact diagonal
identity is `U_bar(s)=V_s U_Palm(s)` in the same gauge. The Heisenberg relative
cocycle is `V_s*`; this relates two reference-dependent clocks without
identifying their frequencies. Rational logarithm/phase enclosures certify
every original binary source-array entry against its target. Target centers
are rounded to dyadic rationals with both component errors charged; the
executed arrays are unchanged. Each requested rational `|s|<=16` word must
pass operator error `1e-9`; the retained words use `s=16`.

Finite preparation also joins. Put `t=66 delta`, where `delta` is the Palm
warm-up error, and `u=75 eta`, where `eta` is the complement-repair error.
Complete positivity and `bar_sigma>=I/75` give the actual unnormalized
successful work state a relative sandwich of radius `r=t+u+tu`, up to its
herald normalization scalar. Logarithm stability gives

```
half-diamond(source_bar(s), actual_completed_modular(s))
    <= delta_word + min(1, |s| r/(1-r)).
```

At **369 warm-up attempts, one global movement herald, then 7440 complement
attempts on success**, this is below `0.001` for each certified word in the
window. The cap-only comparison uses `t` instead of `r`: every complement
attempt already preserves the complete cap marginal. The complete flagged
single-trial preparation/source error is `delta + (5/33) eta + p_max delta_word`.
All failed herald outcomes remain records. Repeating failed Palm trials and
then running one successful complement stage costs fewer than 9882.015
expected native attempts. Preparation records/controllers are unavailable to
later work; attempts and `s` retain their different meanings.

The same clock also respects the continuing cut protocol. Every actual C/E
complement Kraus branch preserves the cap configuration, and therefore
commutes with `bar_sigma` and its modular unitary. This proves covariance
of the **complete labelled instrument**, including arbitrary untouched
references and any finite concatenation of its records. Ideal modular words
can be moved between complete complement attempts without changing those
current-record probabilities. Numerical source words carry their certified
control error; measurements during a partial word are outside this claim.

An executed cap Ramsey probe uses native codes 328 and 73, with the same
outside configuration and cap charge one. Seven existing source pulses
prepare their coherence and seven analyze it. The complete work-register
bit-0 PVM reads `(1+cos(s log(38/37)))/2`; all 512 coordinates contribute to
its two outcomes. At `s=16` the match probability lies between `0.95` and
`0.96`. The 2056-pulse experiment retains padding throughout; neither a
membership projection nor outcome normalization is inserted. Configuration
preparation and the projective bit read remain admitted apparatus, with
their physical costs outside this pulse count.

This supplies the finite faithful-state/inner-flow hypotheses used in
[`PrivateInner`](../../Lean/Dynamics/PrivateInner.lean) and the finite Gibbs/KMS
identity in [`StructuralInheritance`](../../Lean/QFT/StructuralInheritance.lean)
for the selected native entropy completion. The cut schedule, discarded
records, coherent services and cap choice remain explicit. Passive cut repair
is not the modular unitary. Exact regional invariance belongs to the ideal
completed reference; actual finite-prepared flows have the stated error
bound. No boost, area, smooth spacetime or laboratory-clock identification
is inferred.

## One reference for overlapping caps

[`common_caps.py`](common_caps.py) supplies one native-selected reference
for a specified family of overlapping bare cap algebras. Fairly select
`R=(0,1,5,7)` or `S=(0,1,2,8)`, then a uniform complementary native seam
and a fair placement coin. The 20-seam union is connected on the ten ports
outside `I=R intersection S=(0,1)`. For both C and E, the exact native
channel has fixed algebra `A_I` and quantum gap at least `1/600` per
attempt. Its unread marginal converges to the tracial completion `E_I`.
Records may remain stored if subsequent work does not access them or use
them for feedback.

Applied to the same globally heralded Palm reference, this protocol selects
`sigma_common=E_I(Palm)`. Its per-configuration weights are

| Occupation of the protected core I | Weight | Configurations |
| --- | --- | --- |
| 0 | `68/4500` | 45 |
| 1 | `69/4500` | 20 |
| 2 | `60/4500` | 1 |

This single density lies in `A_I`, hence its modular group preserves every
bare `A_R` containing `I`; all 1024 such cap algebras are checked. It is
central in `A_I`, so **the shared intersection's own modular action is
the identity**. Nontrivial clocks require larger-cap coherences transferring
occupation between the protected core and another cap port. The mixed
preparation preserves the core marginal, not the individual larger-cap
marginals. The state and clock are common across these caps without
independently selecting a reference for each observer.

One 2042-pulse source word implements the common modular action. Two
cap-supported selected-pair Ramsey probes use it unchanged: codes 328/73
for `R`, and codes 0/73 for `S`. With their preparation/analyzer mixers,
they cost 2056 and 2068 pulses and both read
`(1+cos(s log(69/68)))/2` through the complete 512-coordinate bit-0 PVM.
At the retained `s=16`, the ideal probability is approximately
`0.986421907373`. The complete binary source arrays and both mixer errors
are certified; padding is retained without projection or renormalization.
Configuration preparation, coherent source services and the projective bit
read remain admitted apparatus. The selected pairs have the same cap
exterior; this does not assert that each source mixer is a local operation
on arbitrary exterior inputs.

With 369 Palm warm-up attempts, one global movement herald, and 18,600
mixed-cap attempts after success, the same relative-sandwich proof bounds
the source/common-actual-state modular half-diamond discrepancy by
`0.000989302` at `s=16`. The preparation costs fewer than 21,042.015 expected
native attempts. The complete flagged single-trial bound retains every
herald failure. Each other requested rational `|s|<=16` requires its own
word certificate; a continuous family of floating implementations is not
certified by the retained word. Regional invariance is exact for the ideal
common reference; finite-prepared modular actions have the stated error.

The all-region boundary is also exact. On the binary two-occupation shell over
`N>=3` ports, a faithful diagonal reference with weights `w_ij` preserves
**every bare two-port algebra** iff `w_ij=lambda*r_i*r_j` with positive
factors. Equivalently, every four distinct ports satisfy
`w_ij*w_kl=w_ik*w_jl=w_il*w_jk`. Local matrix units require their modular
ratio to be independent of the spectator port; factorization gives the
converse. The generator is a scalar minus `sum_i log(r_i) N_i`.
Native Palm, the one-cap completion and this common-core completion fail
that all-pair condition. A connected symmetric native repair preserves a
factorized reference only when all factors agree. A concrete unprotected
two-cap schedule has a connected 28-seam union and uniquely selects
`I66/66`, whose stationary modular action is trivial.

These are a constructive protected-core family and an exact boundary for
the indicated regional algebras. The protected core and preparation protocol
remain specified choices. No arbitrary spacetime region net, physical boost,
laboratory time or continuum geometry is inferred.

## A common observer from two local stopping clocks

[`local_preparation.py`](local_preparation.py) connects the regional stopping
theorem to the common-cap completion. On the same binary total-two carrier,
count incidents at the shared core `I=(0,1)`, including zero currents. After
the declared warm-up count, retain both outcomes of the **next specified**
incident movement flag. This is not the first successful movement. On
success, run the mixed complementary repair of the two caps above, stopping
after the declared number of incidents at the first cap. The core has no
incidents during this complementary repair, so it cannot serve as that
second counter.

The first observer has incident fraction `3/10`; the second has fraction
`5/24` during mixed completion. The complete ideal herald succeeds with
probability `5/33`. Its completion selects one common density:

| Core occupation | Per-configuration weight | Multiplicity |
| --- | --- | --- |
| 0 | `2/225` | 45 |
| 1 | `1/36` | 20 |
| 2 | `2/45` | 1 |

The nonzero elementary modular gaps are `log(25/8)` and `log(8/5)`. These
follow from this local preparation; the global-Palm gaps are not substituted.
The density lies in the protected-core algebra and is faithful. Its one
inner modular action therefore preserves all 1024 core-containing bare cap
algebras. The common core's own action remains trivial. This instantiates
the finite faithful-state, inner-action and Gibbs/KMS contracts with
`H=-log(sigma)` and inverse temperature one in modular units; it does not
turn the unread repair semigroup into that unitary action.

Both stopping clocks have finite controllers. Every tick has a bounded
remote-wait deadline. The existing reusable native auxiliary supplies each
scheduler choice through its coin-only instrument, whose effect is `I/2`
even when the auxiliary is correlated with the work. No external random
service, auxiliary reset or independent auxiliary state is assumed. A
bounded rejection selector retains its failure. The two complementary
menus have different sizes, so their finite acceptance factors differ.
They are included in the actual channel rather than canceled by conditioning.

For either phase let `A` be the accepted incident subchannel, `R` the
accepted remote subchannel, and `t` the state-independent selector failure
probability. For a per-tick deadline `h`, the executed success and timeout
maps are

```
T_h = A (Id + R + ... + R^(h-1)),
F_h = t (Id + R + ... + R^(h-1)) + R^h.
```

They form a complete instrument. Repeating `n` ticks gives success `T_h^n`
and terminal failure `F_h sum_(j=0)^(n-1) T_h^j`, preserving the actual
disturbed failure state. The code evaluates these maps by binary powering,
without enumerating exponentially many fine histories. In the uncapped
comparison, the positive-contraction similarity proof applies on the
orthogonal complement of the common fixed algebra, replacing the tracial
projection by `E_I`. Finite selector and waiting losses enter the complete
error **before** any conditioning or logarithm stability bound.

The retained flags are warm-up timeout, herald failure, completion timeout,
and success. Preparation seam/coin/current tapes may be stored but are not
accessible to later feedback in this comparison. The controller needs
private bounded loop/deadline counters; the public preparation clocks count
local incidents. This removes the global warm-up timer and separately
selected cap references, not the need for controller hardware.

The retained precision uses **496 warm-up ticks, one specified herald tick
and 19,200 completion ticks**, with a 256-slot deadline per tick and 32
selector rounds. Its complete preparation distance is below `5.682e-8`;
conditional success distance is below `7.500e-7`. Expected work plus auxiliary
cost is below **105,445 native calls**. The conservative finite maximum is
659,460,352 calls. The reusable controller-register subtotal is 51 bits;
program tables, interpreter storage, private tapes and quantum registers
are excluded from that subtotal. These are bounds for the stated controller,
not optimal costs or physical durations.

## Continuing source-clock reads of that prepared state

[`local_history.py`](local_history.py) executes two coherent interventions
and native current reads on the **same prepared mixed carrier**. It does
not reset to code 328 or code 0 between probes. Each exact two-address
rotation fixes every other native and padding address; the cheaper
unconditional XOR mixers valid for selected pure inputs would not implement
this mixed-input instrument. Both source words use the same local reference
and modular parameter. The second rotation's sign depends only on the first
retained signed current. Both actual C/E coin/current instruments retain
their zero outcomes and continuing quantum states.

The complete 512-address source arrays are certified entry by entry against
rational enclosures of their exact targets. Every native export also has
an explicit complement abort, on which the work is discarded and an
untouched reference is retained. No escaped probability is renormalized
away. The receipt joins the finite preparation maps to this history and
compares the complete quantum outputs, including all three preparation
failures, both export aborts and every native read label. It reports the
full source-pulse/service budget separately from preparation's work and
auxiliary attempts. The all-input/reference guarantee is compositional;
the retained numerical preparation witness starts from one explicit basis
state and does not replace that proof.

At retained modular parameter one, the full preparation/history bound is
below `6.32e-8`. The two-probe source word uses **31,768 pulses and 63,536
proper-code service events**, plus two native current reads and two complete
native exports. Its single-clock comparison with the actual finitely prepared
state's modular action has a separate bound below `8.437e-5`. The source-array
enclosures and analytic preparation bounds are certificates; intermediate
floating state propagation and its reported residuals are diagnostics.

Coherent proper-code load/export and source controls remain specified
architecture. Selecting these two address pairs requires exterior context;
this is not a derivation of cap-local coherent control on arbitrary
spectators. No projective work-bit read or independently selected pure probe
is needed. The comparison supplies neither physical service fidelity nor a
laboratory time, temperature, atom or field dictionary.

## The boundary includes noncommuting references

[`compatibility.py`](compatibility.py) removes a possible loophole in the
all-region test. Let a faithful state on the binary two-occupation shell over
`N>=3` ports preserve every bare two-port algebra under its modular group,
and put `H=-log(sigma)`. For each port `i`, its occupation `N_i` belongs to
every `A_{ij}`. Thus `[H,N_i]` lies in the intersection over `j != i` of
those algebras. Their off-diagonal supports are disjoint, so that intersection
is diagonal; `[H,N_i]` also has zero diagonal and hence vanishes. All joint
occupation eigenspaces on this shell are one-dimensional, forcing `H` and
`sigma` to be diagonal. This includes `N=3`, where the intersection need
not be only `span(I,N_i)`.

Spectator-independent modular ratios now give exactly the positive product
weights `w_ij=lambda*r_i*r_j` already characterized above. Consequently
**every faithful reference**, not only a preselected diagonal one, preserves
all bare two-port algebras precisely in that product family. The executable
certificate supplies forbidden-entry witnesses for every native
off-diagonal Hamiltonian entry and exact independent constraints for the
remaining diagonal energies. Connected symmetric repair still selects only
the tracial member at equilibrium. Coherences cannot evade this specified
all-region boundary; protected regional families and changed observable
algebras are different contracts.

## Replay

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference --check model/native_modular_reference/results.json
python -m model.native_modular_reference.preparation --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.full_flow --check model/native_modular_reference/full_flow_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.regional --check model/native_modular_reference/regional_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.completed_clock --check model/native_modular_reference/completed_clock_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.common_caps --check model/native_modular_reference/common_caps_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.local_preparation --check model/native_modular_reference/local_preparation_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_modular_reference.local_history --check model/native_modular_reference/local_history_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/native_modular_reference
```

Controls cover finite local-selector/wait timeouts, initially entangled
references, noncommuting compatibility controls, the faithful reference and
4+6 pair census, omitted native
coordinates, all padding columns, inverse/common phase, scalar amplification
that defeats an unweighted telescope, retained preparation failures, uniform
reference, native waiting, checkpoint continuation and hostile ndarray
subclasses. The retained receipt is software evidence, not measured data.
