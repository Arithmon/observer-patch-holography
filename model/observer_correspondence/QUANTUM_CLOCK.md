# A quantum instrument in the observer's incident clock

The [diagonal correspondence](README.md) transports complete finite records
and terminal populations through a source-derived change of representation.
The extension here constructs an ordinary completely positive effective
instrument. It preserves complete local quantum outputs, their correlations
with an untouched reference, and adaptive local interventions. Its finite
implementation retains timeouts and has a uniform error and work budget.

For the specified C and E laws the unread effective transfer is positive on
the Hilbert space of operators. A finite logarithm requires a separate
strict-positivity certificate. That certificate depends on both the chosen
clock region and the operator domain. A permitted signed record law provides
an explicit negative-mode obstruction despite identical classical repair
data. None of these statements selects a physical record law.

## The same source and observation contract

The source is the original twelve-port carrier with all thirty native seams.
Scheduled seams are independent uniform choices. Every repair uses the
original fair native coin and its declared C or E record instrument. Hidden
coin/record alternatives are summed as CP maps, not as coherent amplitudes.
The [record-instrument distinction](../record_laws/INSTRUMENT_CONTRACT.md)
remains part of the premise.

| Interface | Retained realization |
| --- | --- |
| Quantum domain | The twelve-dimensional total-one sector for executed histories, with an optional untouched reference. The theorem below applies to the full finite binary tensor algebra and to invariant sectors. |
| Observers | Region A is `(0,1)` and region B is `(0,5)`. Their union O is `(0,1,5)`. |
| Clock | A scheduled repair ticks each region whose ports its seam meets, even when its current is zero. Stop after two incidents on O. The union meets twelve of the thirty seams, so `p=2/5`. |
| Committed marks | `(tick_A,tick_B,delta_n0)`. There are five success marks: A-only with zero increment, B-only with zero increment, and a shared tick with increment `-1`, `0` or `+1`. |
| Preparation | A coherent root superposition on ports zero and one, or its root/reference Bell version, with the remote ports in vacuum. No initial occupation disclosure dephases this input. |
| Feedback | When enabled, the first B-only mark requests the local phase `Z_0`. It acts on the continuing conditional state. |
| Clock and control cost | `Z_0` is a supplied coherent control, not a native repair or an additional incident-clock tick. Its gate cost is recorded separately. |
| Analyzer | After two successful hits, the local phase-plus, phase-minus and complementary outcomes form a complete quantum measurement, tensored with identity on the reference. The complement is retained. Timeout branches keep their conditional quantum states without running the analyzer. |
| Finite execution | Each requested interval has a cap of `h` native operations. A timeout aborts the remaining requests and remains an output branch with its preceding records and unnormalized state. |
| Continuation | The common runtime retains the conditional quantum state, classical record and control memory through serialization and subsequent work. |

Coherent preparation, the phase gate, the analyzer and faithful record access
are declared local quantum services already used by the native quantum
controls. They are not derived from the diagonal repair grammar. Their
service errors are idealized here; no hardware accuracy or duration is
inferred. Attempt caps are external execution services and do not supply
an additional clock observable to the feedback policy.
When working in a fixed-number sector, the declared controls must preserve
that sector. The executed phase gate and analyzer satisfy this condition;
the broader arbitrary-local-instrument theorem uses the full tensor
domain when an intervention changes total number.

## A CP square-root waiting channel

Write the scheduled unread superoperator on source density matrices as

```text
P = A+B,                 B = sum_m B_m,
A = remote-seam contribution,
B_m = union-incident contribution with committed mark m,
p = incident seam fraction,     r = 1-p,
D = I-A.
```

The letter `I` here is the identity superoperator. Remote attempts act
entirely on the complement of O. They have constant probability `r`, so
`A^*(1)=r 1`; similarly `B^*(1)=p 1`. For `0<p<1`, let `Q=A/r`, a
channel using only remote native repairs. The case `p=1` has no waiting
and follows by setting `A=0`.

Set

```text
c_n = binom(2n,n)/4^n,
R = (I-A)^(-1/2) = sum_(n>=0) c_n A^n,
J = sqrt(p) R.
```

The series converges in completely bounded trace norm because
`||A^n||_diamond=r^n`. Every coefficient is nonnegative, so `R` is CP.
Since `sum_n c_n r^n=(1-r)^(-1/2)=p^(-1/2)`, J is trace preserving.
It is a genuine remote channel, with the operational mixture

```text
J = sum_(n>=0) w_n Q^n,
w_n = sqrt(p) c_n r^n,             sum_n w_n = 1.
```

This is a half-geometric, or negative-binomial shape-one-half, waiting
distribution. Its mean number of remote native repairs is `r/(2p)`.
Its definition uses the source seam fraction, not a fitted time scale.

The native completed incident instrument and the effective instrument are

```text
K_m = B_m (I-A)^(-1),
T_m = R B_m R = J (B_m/p) J.
```

Every `T_m` is CP and their sum T is trace preserving. Individual marked
maps need not be self-adjoint or positive as Hilbert-space operators.
Complete positivity of the instrument is a separate assertion from the
spectral positivity of its unread sum.

## Complete local histories and arbitrary references

The exact intertwining identity is

```text
T_m J = J K_m.
```

It follows directly from `R^2=(I-A)^(-1)`. Every local intervention on O,
including its selective branches and any local controller memory, commutes
with J because J acts only on the remote tensor factors. The same holds
after adjoining an arbitrary untouched reference.

Begin the source experiment in any source/reference density matrix `rho`.
Begin the effective experiment in `(J tensor id)(rho)`. For every retained
adaptive record branch, induction gives

```text
effective terminal state = (J tensor id)(source terminal state).
```

The induction applies the same record-dependent local intervention and
then the corresponding marked clock map on each branch. It retains all
outcomes and all their unnormalized states. Since J is a trace-preserving
remote channel, tracing out the remote factors removes it. Thus the two
experiments have exactly the same joint local records, terminal state on O,
controller memory and reference. Their complete local output channels
agree, not just a chosen expectation value.

The initial encoding matters. An arbitrary other input need not satisfy
`J rho=rho`; replacing its effective input by the unencoded state can
change predictions. The executed root superposition and Bell inputs have
remote vacuum, which each remote native map fixes. They therefore satisfy
`J rho=rho` exactly and use the same initial state in both representations.

Full terminal source states are not asserted equal. Their relation uses
J, whose inverse need not be positive. The local-output result needs no
physical inverse decoder. Likewise, a control that touches the remote
region, reads an inaccessible record or consults global waiting counts
does not inherit this theorem merely by being a finite operation.

## A finite ordinary instrument with rational weights

Two independent half-waits have a geometric total because

```text
sum_(a=0)^n c_a c_(n-a) = 1.
```

Thus T has a direct sampling construction. Draw a total remote count N
with probability `p r^N`; conditional on N, draw a split a with probability
`c_a c_(N-a)`. Execute the right remote half, a uniformly selected incident
native operation with its committed mark, and the left remote half.
All conditional split probabilities are exact rationals. The total
expected number of native operations per completed interval is
`E[N+1]=1/p`, equal to the native scheduler's expectation.

For a cap of h operations, retain totals `N<h`. The successful maps are

```text
T_(m,h) = sum_(a+b<h) c_a c_b A^a B_m A^b.
```

If the draw reaches h remote operations without success, execute the
remote prefix and retain the timeout map

```text
Q_h = A^h.
```

The sum of success traces is `s=1-r^h`, and the timeout trace is `r^h`.
This is an ordinary complete CP instrument on every input and reference.
It uses only rational mixtures of source-native maps. It is not obtained
by conjugating a finite native cap with a generally nonpositive inverse.

The successful difference `T_m-T_(m,h)` is CP term by term. Its aggregate
trace is exactly `r^h`. Embed the completed output in the capped alphabet
with a zero timeout branch. For every input/reference density matrix the
complete CQ trace distance is then at most `r^h`: the lost success mass
and retained timeout mass each contribute `r^h` to the trace norm.

For k requested intervals with local adaptive trace-preserving controls
and abort after the first timeout, the successful capped histories are
termwise subchannels of their completed counterparts and have total mass
`s^k`. Hence the uniform complete-output allowance is

```text
epsilon_cap(k,h) = 1-(1-r^h)^k.
```

The native before-hit cap has maps
`K_(m,h)=B_m sum_(n=0)^(h-1) A^n` and the same timeout `A^h`.
It has the same allowance to its own completed target. Combining both
allowances with the exact local correspondence gives the conservative
native-capped versus effective-capped local-output bound

```text
D_trace(output_native,h, output_effective,h)
    <= min(1, 2 epsilon_cap(k,h)).
```

This comparison assumes the declared J-related preparations, or the
remote-fixed preparation executed here. A separately truncated preparation
of J adds its own complete-output error; it cannot be silently omitted
for an arbitrary input. For a preparation cutoff h that error is the
half-wait tail `1-sqrt(p) sum_(n<h) c_n r^n`.

Neither complete-output ensemble renormalizes its successful branches.
The sampled trajectories below use ordinary normalized conditional states;
their sampling probabilities remain part of the complete instrument.
A phase analyzer or any later local output processing is applied after
these budgets and cannot increase the complete-output trace distance.
A bounded observable of norm one has
expectation error at most twice the trace-distance allowance; a POVM
outcome probability has error at most the allowance itself.

## Spectral positivity and its limits

The existing [native Liouville result](../native_action/README.md) proves
that every binary C or E unread seam superoperator is self-adjoint and
positive in the Hilbert--Schmidt inner product, with spectrum in `[0,1]`.
Thus `0<=A<I`, `B>=0` and `D-B=I-P>=0`. It follows that

```text
0 <= T = D^(-1/2) B D^(-1/2) <= I.
```

This is a source-derived positive Euclidean transfer on the Hilbert space
of operators. It is also the CPTP channel constructed above. If a separate
certificate gives `B>=alpha I>0`, then `D<=I` implies `T>=alpha I`, and

```text
H = -log T,
exp(-n H)=T^n                         for integer n>=0,
||H-sum_(j=1)^L (I-T)^j/j||
    <= (1-alpha)^(L+1)/((L+1)alpha).
```

H is a Liouville operator, not a Hamiltonian on source state vectors.
Its integer parameter counts union incidents. The two observers' own
clock increments remain correlated marks; substituting their mean rates
for the joint clock record would change the experiment.
Neither arbitrary fractional powers `T^t` nor the unitary continuation
`exp(-itH)` are asserted to be native quantum channels or physical waiting.
The finite capped instrument is a numerical implementation with its own
timeout branch; no logarithm of that augmented instrument is used.

Strictness has several distinct domains:

* On the full binary operator algebra, an incident native triangle gives
  the C bound `B>=I/M`. The exact 64-dimensional triangle count matrix
  `C_triangle=2 sum_e Phi_(C,e)` is symmetric and obeys the annihilator
  with roots `{2,3,4,6}`. It follows that the three-seam sum is at least I.
  The actual carrier triangle `(0,1,5)` is incident to the stated union,
  yielding the full-binary floor `1/30` for C.
* On the executed twelve-dimensional total-one sector, direct source
  count matrices give incident floors `7/30` for C and `1/10` for E.
  The population block is symmetric and its lower bound follows from
  its diagonal minus its off-diagonal row sum. Every off-diagonal
  matrix unit is a scalar block checked from the original endpoint and
  record rule. These checks cover all 144 operator matrix units.
* The finite certificate on the actual 66-dimensional binary total-two
  sector supplies incident floors `2/15` for C and `1/60`
  for E under the union clock. These also give conservative transfer
  floors on that domain. All 4,356 operator matrix units are partitioned
  into the existing signed-difference blocks. Where the row bound alone
  is insufficient, exact positive-pivot elimination certifies the
  shifted E block. This is not an inference from the total-one result.
* For the smaller clock region `(0,1)`, E has a zero mode already on
  total-one operators: `X=|particle at 0><particle at 1|`. Remote maps fix
  X, while every incident E map kills it. Therefore `T_E X=0` and a finite
  full-domain `-log T_E` does not exist for that clock. For C the same
  mode has eigenvalue `4/9`.

The classified [signed C record law](../record_laws/README.md) makes the
boundary sharper. It keeps the same classical transitions and definite
current labels but changes the sign of the waiting coherence. Under the
two-port clock its same X mode has eigenvalue `-4/9`. Its effective clock
maps remain a CP instrument by the half-wait construction, but the unread
transfer cannot equal `exp(-H)` for a self-adjoint H. Classical current
probabilities alone therefore do not establish positive Euclidean
correspondence. The zero- and negative-mode controls are not failures of
the stated union-clock C/E construction.

## A phase-sensitive prediction in the same local clock

Let X be the root off-diagonal between a particle at port zero and a
particle at port one, with remote vacuum. Of the twelve union-incident
seams, nine touch zero or one. For C, their edge `(0,1)` kills X and the
other eight retain one half of it. The remaining three seams are B-only
and preserve X. Thus a completed union tick multiplies X by `7/12` for C
and `1/4` for E.

The first B-only mark selects `Z_0`, changing the sign of precisely the
three-seam contribution. With feedback, the first factors become `1/12`
for C and `-1/4` for E. The second tick is unchanged. For the unreferenced
root-plus preparation, the completed two-tick local phase contrast,
`Pr(plus)-Pr(minus)`, is

| Record law | Baseline | B-only-controlled phase |
| --- | ---: | ---: |
| C | `49/144` | `7/144` |
| E | `1/16` | `-1/16` |

All diagonal populations and pre-analyzer classical mark probabilities
still agree between C and E. For this plus input, the phase contrast
accesses coherence that those classical data do not contain.
The capped contrast is also known exactly.
The remote factor on X is r, so every capped marked coherence is s times
its completed value, where `s=1-r^h`. The unconditioned success-only
analyzer contrast is therefore exactly `s^2` times the corresponding
table entry, for both native and effective chronologies. This gives an
independent finite-cap observable check. The complete local quantum output
still uses the uniform cap allowance, and actual arithmetic rounding
requires its separately measured budget.

The Bell preparation has a different local readout. Its reduced source
input is the incoherent mixture of the two root occupations. The executed
analyzer acts only on the source, so its plus/minus contrast is zero for
both candidates, with or without feedback. The table does not describe
this Bell probe. Its complete root/reference output retains different
conditional reference coherences under C and E, even when all local
analyzer probabilities agree. That comparison tests quantum disturbance
and reference correlations; no joint Bell analyzer is added to the
declared measurement service.

The retained Bell state also supplies a direct entanglement witness before
the destructive analyzer. Retain the public prefix consisting of two shared
zero-increment marks, `((1,1,0),(1,1,0))`. For C, each such capped mark
multiplies the root off-diagonal by `5s/24`. The unnormalized root/reference
state therefore has cross entry `25s^2/1152`. The corresponding two crossed
occupation diagonal entries vanish. Its partial transpose has the exact
two-by-two principal minor

```text
det_minor = -625 s^4 / 1327104 < 0.
```

This proves entanglement of that retained conditional branch for both
chronologies, without discarding its probability from the complete output.
The first mark does not request feedback, so the same witness applies under
both policies. At cap sixteen the minor is approximately
`-0.0004704190635525181`. E gives zero for this minor; that zero alone is
not a separability certificate. This is an intermediate-state witness,
not an extra measurement service or a claim that the subsequent destructive
local analyzer preserves entanglement.

## Resources, verification and remaining obligations

For the union clock, `p=2/5` and `r=3/5`. Each requested capped interval
uses at most h native operations and has expected work `s/p`, with
`s=1-r^h`. Two requests with abort on timeout therefore use at most `2h`
native operations and expected work `s(1+s)/p`. In the completed limit
the expectation is five. The controlled phase is used with probability
`s/4`, because three of the twelve incident seams give the B-only mark.
It is a separate coherent gate, not a native clock increment.
The phase analyzer is used with probability `s^2`, and preparation is one
declared plus or Bell service. These control and analyzer counts are not
folded into the native-attempt expectation.
The default plus preparation is charged by the common runtime. The Bell
state is supplied as a joint source/reference input to that runtime, so its
one joint-preparation service is recorded separately; it is not a derived
native pulse cost or a second preparation inside the program.

There are 25 successful two-mark words, each with three analyzer outcomes,
and six timeout prefixes: one before the first hit and five after it.
The complete declared local output therefore has 81 classical labels,
each carrying its unnormalized state on O and the untouched reference.
The runtime additionally stores the full conditional source state for
continuation; that larger state is not asserted equal between native and
effective representations. Internal fixed-program
not-run placeholders are collapsed to these same timeout prefixes; they
do not create a second timeout event or an uncharged analyzer call.
The retained runtime stores 108 padded terminal branches, including zero
branches, at work/reference dimension at most 24. Its terminal state arrays
therefore contain at most 62,208 complex entries, with at most three public
symbols per padded branch. These are terminal representation counts, not
a bound on all temporary storage during factorization or verification.

The source-native and effective executions must retain every mark prefix,
timeout, analyzer outcome and continuing quantum state. Exact source-map
checks on matrix units establish reference extensions. Independent
endpoint/coin recursion checks the native instrument; split-wait histories
check the effective mixture. Complete CQ comparison must include the
trace norm of local output states, not just a comparison of traces. Each
runtime may separately compare its full state with its own exact native
or effective oracle. Actual runtime rounding is a further budget alongside
the analytic cap tail.
Serialized continuation must resume the saved quantum state and memory,
with fresh record environments for subsequent native events.

The retained [quantum-clock producer](quantum_experiment.py) uses `h=16`.
Writing `s=1-(3/5)^16`, the exact completed-target allowance is

```text
delta = 1-s^2
      = 0.00056414239488010053599115608064.
```

The common runtime receives a separate complete-output allowance
`eta=10^-10` for each executed preparation, candidate, feedback choice
and representation. It converts the actual binary64 output entries to
their exact rational values and sums their absolute entry errors against
the exact full-state oracle, thereby bounding the complete trace distance.
This is an a-posteriori bound for those executed probes, not an all-input
diamond-norm certificate for the floating Kraus implementation. The
analytic cap bound itself is uniform over arbitrary references on its
stated preparation domain.

| Cap-16 quantity | Certified expression or outward upper bound |
| --- | ---: |
| Per-interval timeout | `(3/5)^16 = 0.0002821109907456` |
| Exact capped to completed complete output | `delta < 0.000564142395` |
| One runtime output to its completed local target | `delta+eta < 0.000564142495` |
| Native runtime versus effective runtime, complete local output | `2delta+2eta < 0.001128284990` |
| Maximum native calls | `32` |
| Expected native calls | `s(1+s)/(2/5)`, approximately `4.997884366536` |
| Expected enabled phase controls | `s/4 = 0.2499294722523136` |
| Expected terminal analyzers | `s^2 = 0.99943585760511989946400884391936` |

The retained [receipt](quantum_results.json) evaluates both C/E laws, both
feedback choices and both plus/Bell preparations, with native and effective
representations for each case: sixteen common-runtime executions. Every
execution is compared with its own exact full-state oracle within eta.
Its saved prefix is serialized and restored, and continuation reproduces
the fresh complete execution's records, resources and output arrays.
For each of the eight native/effective pairs, direct exact comparison of
all local CQ entries gives the sharper trace-distance upper bound
of approximately `1.7343294049995123e-7` before runtime arithmetic. The exact rational value
in the receipt is used for verification; an outward decimal upper bound is
`1.73432941e-7`. Adding both numerical allowances gives the executed
runtime-pair bound `1.73632941e-7`. This probe-specific improvement does
not replace the uniform cap theorem.

These work counts describe the source-native sampling construction. The
matrix-factorized common-runtime execution separately charges each
requested aggregate interval its declared cap budget. A sampled native
history retains its actual seam, coin, mark and operation count; it is
an execution witness, not an estimate of the complete-law error.
The receipt includes 32 such exact-sampling histories: C/E, native/effective,
plus/Bell, caps one/sixteen and fixed seeds zero/one, with feedback enabled.
They retain timeouts, terminal-state digests and the private native
coin/record audit. The controller receives only the declared public marks;
private record labels do not become feedback inputs. These fixed samples
demonstrate execution coverage and supply no confidence estimate for the
complete quantum law.

Finite precision, interval cap, tested operator domains and source bindings
are part of the certificate. Increasing cap and precision is a mathematical
approximation scheme; a bounded implementation does not thereby certify
arbitrarily small error.

This result gives a finite quantum clock and field correspondence with
an explicit access contract. It does not derive the coherent service
apparatus, select C or E, identify physical seconds, prove a relativistic
continuum limit, or establish a physical quantum field theory. Extending
the observation region, changing the record basis, admitting new controls,
or claiming global terminal-state equivalence requires the corresponding
premises and error bounds.
