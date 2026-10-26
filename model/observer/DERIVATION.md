# Finite current-record controller: mathematical account

## Native source and distinct quantum chart

Let `X` be all nonnegative twelve-tuples of total two. For seam `e=(a,b)`
and coin `c`, the native endpoint repair is `F_e,c`, with current
`j=x_a-(F_e,c(x))_a`. Define

`L^e_(c,j) = 2^(-1/2) sum_{x: current(x,c)=j} |F_e,c(x)><x|`.

Each fixed-current map is an injective translation: its inverse adds `j` at
`a` and subtracts `j` at `b`. Each source basis vector has one branch for
each coin. Consequently `sum_r L_r† L_r = I`, and
`W_e=sum_r |r> tensor L^e_r` is an isometry. This is the current law C
(legacy label N1). The **contrast** instrument E (legacy N1-E) groups branches by pre-event contrast
instead. Its classical endpoint probabilities agree with C, but its coherent
sums differ. Classical grouping of E records is not an implementation of C.

All C labels have `c in {0,1}` and `j in {-1,0,1}` on X. All six occur on
the full 78-state source; four occur on the binary shell. The code generates X directly, checks all 312 seam/coin/source outputs
against independent floor/ceiling endpoint averaging, and constructs C's
coherent branches from the recorded currents.

## Fixed space, adaptive write and complete positivity

Let `R=C^6`, `H=C^78`, and
`H_cap=H ⊕ (R tensor H) ⊕ (R tensor R tensor H)`. Its sectors have
dimensions 78, 468 and 2808. They all exist before operation. Put `e=(0,1)`
and `f=(0,5)`, with `s(r)=f` if its current is positive and `s(r)=e` otherwise.
The two nonfull-sector maps are

```
V0 |psi> = W_e |psi>
V1 (|r> tensor |psi>) = |r> tensor W_s(r) |psi>.
```

The orthogonal retained labels and `W_s†W_s=I` give `V0†V0=I` and
`V1†V1=I`, even on arbitrary unreachable first-letter/source combinations.
Extend V0 and V1 by zero to other input sectors, calling them A0 and A1.
Let P2 project onto the full sector. Then

`T(rho)=A0 rho A0† + A1 rho A1† + P2 rho P2`

is CPTP because `A0†A0+A1†A1+P2=I`. It is one fixed map, with all
record-dependent selection internal to its definition. Cross-phase coherence
is dephased. From phase zero, two applications compose the isometry V1 V0;
at full capacity every operator is fixed. Neither statement implies global
unitarity, unlimited operation or uniqueness of this controller.

These isometry statements concern the source/tape tick without classical
readout. The executed forecast interface reads the spectral projections of
`sum_r q(r) Q_r` (with abstention on the two off-prior letters) before the
second tick. Its forgotten-output action is pinching in the forecast classes.
It is CPTP and preserves every Q_r, but can destroy quantum coherences.
The first coherent write can create record coherences even when the initial
working density is diagonal. V1 retains each old-record block without mixing
distinct labels, and all scored final effects are diagonal in those labels.
They therefore receive no contribution from off-diagonal record blocks;
forecast pinching leaves their complete history probabilities unchanged.
The full observed execution is therefore not claimed isometric. The quantum
coherence witness below concerns the coherent feedback map before that readout.

The checker verifies every nonfull isometry column exactly: each has two
entries of amplitude `1/sqrt(2)`, and different columns have disjoint output
supports. Thus the Gram identity holds for all input coherences, not merely
on basis probabilities. The 78 two-event columns likewise have four
disjoint-support entries of amplitude `1/2`. Earlier input loads are also
reconstructed by reversing both recorded current translations.

## Which memory algebra is protected

On the inhabited subspace `K=H1 ⊕ H2`, let Q_r be the projector that the
first stored letter equals r, independent of the current phase and loads.
Every subsequent map intertwines Q_r, so `T*Q_r=Q_r` on K. All linear
combinations of these projectors form the protected classical first-record
algebra. Repeated projective reads have the same value with probability one.
This is a statement **after that record is written**. Extending Q_r by zero
onto phase zero does not make it a global fixed observable: a first write
creates the record. New second-slot cylinders become protected only once
committed, when this machine halts.

For off-diagonal first-record operators, partial tracing the next write gives
cross factors `W_s(t)† W_s(r)`, which need not be identity. An explicit
reachable witness starts with the normalized equal superposition of binary
occupations `{0,2}` and `{1,2}`. After V0, the first-record reduced-density
entry `(r,t)=((0,1),(0,0))` is `1/4`. Its associated working state has
occupations `{1,2}`. The f and e isometries on that state overlap by `1/2`,
so V1 changes that entry to `1/8`. Classical diagonals are unchanged.

This illustrates a general exact constraint. Suppose a CPTP channel from
`R tensor W` to `R tensor W'` preserves every observable in the full matrix
algebra of R: `Lambda*(A tensor I)=A tensor I` for all A. The unital CP map
Lambda* then has these operators in its multiplicative domain: equality
holds in its Schwarz inequality for A and A†. Hence

`Lambda*((A tensor I) X (B tensor I))
 = (A tensor I) Lambda*(X) (B tensor I)`.

Since `I tensor O` commutes with all A, its image also commutes with every
`A tensor I`; it must be `I tensor Phi*(O)`. Therefore changing only the
R preparation in a product input cannot change any W' output statistic.
Genuine record-to-working-state control is incompatible with universal
preservation of the **entire** unknown quantum record state. Classical
diagonal preservation allows the controlled isometry above. This is an
application of the standard CP multiplicative-domain argument, not a new
physical no-go and not an obstruction to observers with classical records.

The [postcommit observer realization](../common_source/DERIVATION.md#a-native-postcommit-operational-observer)
checks every `ConsensusTower` and `OperationalObserver` field on this K,
with its actual first-event state, `e`/`f` seam controls and declared blank-buffer
refinements. Its accessible algebra is the classical first-letter algebra;
its typed prediction advances issuance to retention of the same evidence.
This is an exact mathematical realization, without a compiled Lean
inhabitant or inter-observer overlap claim. The next-current forecast below
is a separate function of the record: equal record probabilities need not
give equal forecasts.

## Prediction and independent controls

The supplied prior is uniform over the 66 binary configurations. The first
repair and fair coin give the four record masses in README.md. The controller
forecasts the indicator `Y=1{j2>0}` from its stored r before V1 executes.
Direct conditional enumeration gives `q(r) = 0, 5/56, 0, 1/2` in the stated
row order. This is a fixed supplied model/table, with no learning claim.

For squared error, `E[(Y-q(R))^2]=E[Var(Y|R)]`. Let S be the chosen seam,
a deterministic function of R. The best S-only forecast has excess risk
`E[Var(q(R)|S)]`. If R' is independently sampled from the conditional law
of R given S, using q(R') against the actual Y adds twice that variance
above the true record forecast. Thus the within-action shuffle preserves
action information while testing the additional native-record information.
All fractions are also verified directly by independent branch enumeration;
no conditional branch or timeout is discarded. Because the forecast is the
model's own conditional law and the controller is a design choice, these
are consistency checks of the declared construction, not evidence of
learning. There is no delivery service
in this two-event local protocol.

For the feedback audit, a positive first current occurs on ten binary inputs
with coin zero, hence mass `5/66`. On that branch after the first event,
`x0=0`, `x1=1`, and the conditional mean of x5 is `1/10`. A second e event
has port-0 mean `1/2`; a second f event has mean `1/20`. Therefore replacing
e by f on precisely this branch changes the unconditional final mean by
`(5/66)(1/20-1/2)=-3/88`, giving `1/6-3/88=35/264`.
An independently chosen selector of the same marginal leaves the first
event's unconditional uniform shell state uniform under either next seam,
so its mean remains `1/6`. The script retains every branch as a separate
exact contribution.

## Interpretation and limits

The source, tape and controller form one declared finite state/channel law.
The tape, prior, forecast decoder, conditional seam selector, isolation and
halt policy are declared apparatus choices. The construction implements
those choices; it does not require self-consistency to select a unique
apparatus. It does not claim to manufacture the tape using unrestricted
bare C updates on the original working loads. The global all-seam stationary state
is not the halted recorder state, and the local event count is not calibrated
physical duration. No N4 geometry, clock, gauge or field result is silently
transferred onto this C observer construction.
