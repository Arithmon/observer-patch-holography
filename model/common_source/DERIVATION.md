# One finite interacting source: derivation and theorem attachments

This construction implements a specified joint quantum instrument. Its
Hamiltonian, state preparation, local code services and model-coordinate
schedule are declared inputs. It connects existing finite source, CAR,
operator-Gauss and protected-record results on the same evolving work
register. It proves a nonzero finite feedback response with a rational
interval. It does not establish the physical completion target in
[the construction guide](../CONSTRUCTION.md#completion-target-and-acceptance-criteria).

## Existing source and what this construction adds

The [coherent-source construction](../../docs/research/COHERENT_SOURCE_CLOCKS.md)
supplies a common abstract regional algebra: a twelve-address
center, a six-dimensional processor, privately owned quantum buffers and
coherent proper-code services. Independent buffer factors are introduced
once and retained under exact hardware-inventory refinement. The
[source compiler](../observer/source_scalar/compiler.py) executes the
processor's actual generators through these services. We reuse that source;
we do not infer a new existence theorem from a numerical simulation.

The addition here is an explicit interacting gauge/matter/scalar word,
its two complete quantum instruments and its native record-controlled
continuation. Separate scalar and charged-matter calculations become a
single nonproduct process. Every quantum outcome probability is the norm
of a retained branch of that process. Native integer repair remains its
existing declared channel; no compilation of the entire 78-dimensional
native channel into the six-dimensional processor is claimed.


## A common finite quantum gauge and matter code

This result concerns the declared finite observable sector, with a supplied
CAR representation, integer electric link, charge pair, preparation and
source-control services. It gives a common interacting process and readable
records. It does not supply a spatial refinement, physical clock, empirical
identification or the entire Standard Model.

### Physical carrier and the Gauss constraints

Order the four fermion modes as `L0,R0,L1,R1`. The occupation basis is
`(c0*)^n0 (c1*)^n1 (c2*)^n2 (c3*)^n3 |vacuum>` with increasing mode order.
The channels carry charges `q0=+1`, `q1=-1`; each contains exactly one
fermion. Let `b_j=0` mean left occupation and `b_j=1` right occupation.
The logical integer index is `b0+2*b1`, matching the existing compiler's
little-endian convention.

The electric link has basis `|ell>`, `ell in Z`, with
`E|ell>=ell|ell>` and `U|ell>=|ell+1>`. Define the isometry

```text
W |b0,b1> = |occupied modes b0, 2+b1> tensor |-b0+b1>.
```

Its four image columns are distinct orthonormal occupation/flux basis
vectors. On every image column,

```text
Q_L=(1-b0)-(1-b1)=-b0+b1=E,
Q_R=b0-b1=-E,
G_L=E-Q_L=0,   G_R=-E-Q_R=0.
```

Thus `G_L W=G_R W=0`, an operator constraint on the whole physical code,
including arbitrary superpositions and entanglement with spectators.
The supported flux set is `{-1,0,1}`. This is a dynamically invariant
sector of the infinite link, not replacement of the link by a cyclic
three-state rotor. In particular, no finite-dimensional global equation
`[E,U]=U` is asserted.

### The same action couples matter and electric flux

The two dressed hopping terms are

```text
h_j = cLj* U^(qj) cRj + cRj* U^(-qj) cLj.
```

The CAR sign of a ladder at mode `k` is
`(-1)^(sum_(i<k) n_i)`. The producer checks the full four-mode CAR algebra
on all 16 occupation states, not only the code. The two ladders in an
active adjacent-mode hop have equal signs, giving positive code amplitude.
Independent Jordan-Wigner matrices verify these signs; nonadjacent moves
retain the negative fermionic controls. Right-to-left motion increases
electric flux by `q_j`, and left-to-right motion decreases it by `q_j`.
Both moves preserve both ambient Gauss values, even away from their zero
sector. Their images of `W` remain in the same four-dimensional code.

Consequently

```text
W* (h_0+h_1) W = X0+X1 = H_hop,
W* E W = diag(0,-1,+1,0),
W* E^2 W = diag(0,1,1,0) = (I-Z0 Z1)/2 = H_E.
```

The chosen dimensionless action is `H=H_hop+H_E`. Its terms do not
commute: `||[H_hop,H_E]||=2`. Electric evolution therefore changes later
matter observations. Moreover `H_E` contains a joint `Z0 Z1` term, so
it is not merely a sum of two independent channel energies. The current
`J=-i[H_hop,E]` satisfies `i[H,E]=-J`, since `E^2` commutes with `E`.
This is the same electric and matter action; a separately read classical
electric parameter is not substituted for the quantum link.

The scalar interaction in the combined construction may use the first
channel's imbalance `D0=nL0-nR0=Z0`. It commutes with both Gauss generators.
Its later coupling to a neutral scalar preserves the physical code while
acting on the same matter factors as `H_hop` and `H_E`.

### Exact electric-to-matter response

Let `theta=atan2(4,3)`, so `cos(theta)=3/5` and `sin(theta)=4/5`.
Use the finite echo

```text
u=(3I-4iX)/5,
U=u tensor u=exp(-i theta H_hop),
z=(3-4i)/5,
V=diag(1,z,z,1)=exp(-i theta H_E),
T=U* V U.
```

Here `U*` denotes the adjoint, not complex conjugation alone. Every entry
is a Gaussian rational and the three code unitaries are checked exactly
on all columns. On the initial gauge/matter code state `|00>`, the echo
gives, in integer-index order,

```text
T|00> = (1/3125) *
  (2549-1152i, 336-168i, 336-168i, -576-1152i).
```

For the right-occupation detector of channel 0, the outcome-one branch
contains indices 1 and 3, giving exactly

```text
p(nR0=1) = 576/3125 = 0.18432.
```

Omitting `V` leaves `U*U=I` and the same detector returns one with
probability zero. This is a controlled change of the electric operation
with the same state and readout. It is an internal model discrimination,
not comparison with experimental counts. On every other basis input,
the probability that either selected occupation bit flips is also
`576/3125`; both outcomes are retained.

### Native execution, records and composition

The inherited compiler implements actual exponentials of the twelve M6
source generators through its two- and four-dimensional proper codes.
`source_hopping_word()` compiles the two Pauli rotations. The diagonal
electric word includes its identity/global-phase pulse. The echo requires
22 native pulses and 44 coherent load/export service events. No arbitrary
four-dimensional unitary is inserted into the execution path. Analogue
control precision and code preparation, transfer, storage and measurement
remain the inherited apparatus services, with complementary failure
outcomes retained by their full-channel contract.

The exact ideal matrix is reconstructed from CAR/flux operations without
calling the native compiler. Execution checks all four columns and all
16 channel matrix units; tests include a spectator quantum reference and
remapped interleaved registers. Numerical residuals diagnose that replay;
they are not outward floating-point or physical service error bounds.

Each occupation read has complementary projectors `P0,P1` satisfying
`P0+P1=I`, `P_r^2=P_r`, and `P0 P1=0`. The recorded map keeps each
unnormalized block `P_r rho P_r`. It preserves the physical Gauss sector.
Record-controlled subsequent operations also preserve it whenever their
individual code operations do. The retained record label is a separate
classical slot; neither its value nor a zero-probability branch is erased
to restore unwanted interference. With reliable committed records, block
control fixes all previous record projectors even when it changes quantum
coherence between different record values.

The tracial A3/conditional-expectation theorem applies to nonselective
binary readout with retained block algebra
`N={nR0}'=M2 direct-sum M2` and faithful reference `I4/4`:
`E_N(A)=P0 A P0+P1 A P1` is its trace-preserving conditional expectation.
Retaining only the commutative algebra generated by the bit gives a
different entropy completion which additionally averages within each
block. These two operations must not be identified.

`source_step_word(dt)` separately emits the symmetric factor
`exp(-i dt H_hop/2) exp(-i dt H_E) exp(-i dt H_hop/2)` in chronological
order. It uses dimensionless model duration and is not asserted equal to
`exp(-i dt H)` at a finite step. This interface is suitable for composition
with scalar and detector factors and their separately established error
bounds. The exact echo above is a chosen finite protocol, not an error
estimate for that summed-Hamiltonian evolution.

### Existing result and remaining obligations

The nearest existing quantum-link construction is
`code/sm_fermion_current/quantum_link.py`: its flux-dressed CAR sectors
satisfy operator Gauss, but the retained scope explicitly omits electric
energy. The present smaller, differently prepared opposite-charge pair
adds noncommuting electric energy and an actual native observer-compatible
word. It does not inherit that construction's full one-generation census.
The coherent native source primitives are inherited from
`code/m1_source_realization` through the live scalar compiler.

A spatial continuum of the combined interacting process, empirical charge
and coupling identification, physical clocks and responsive geometry are
separate open obligations. Neither this finite code nor its record
conditional expectations discharge them.

## A scalar field and charged matter in the same finite source program

The working register has gauge-sector bits 0 and 1, the first scalar site's
occupation bits beginning at bit 2, and the second site's bits after those.
An index is therefore `b0 + 2*b1 + 4*n0 + 4*d*n1`. We execute `d=2,4`.
The gauge isometry attaches one fermion in each charged channel, with
`b_c=0` meaning left occupation and `b_c=1` right occupation. Charges are
`(+1,-1)` and the electric flux is `E=-b0+b1`. Thus
`D0=n_L0-n_R0=Z0`, not `-Z0`. Every encoded basis vector obeys both endpoint
operator Gauss constraints and has even fermion parity.

With `kappa=9/32` and declared coupling `g=1/4`, define

```
H = H_hop + E^2 + H_scalar + H_int,
H_hop = X0 + X1,
H_scalar = sum_j [P_j^2 + (1+kappa) Q_j^2]/2 - kappa Q0 Q1,
H_int = g sqrt(2) Q0 (n_L0-n_R0).
```

In `H_scalar`, each onsite quadratic is compressed from the infinite
oscillator: `Pi Q^2 Pi` and `Pi P^2 Pi`. Squaring the finite quadratures
instead would remove the boundary contribution `d/2` and change the model.
The bond and interaction use compressed linear `Q`. The interaction is
neutral under the gauge group and commutes with both Gauss generators and
fermion parity. It does not commute with the charged hopping:
`[H_int,H_hop]=2 i g sqrt(2) Q0 Y0`. This is a coupling between the sectors,
not a product of separately supplied response probabilities.

For `d=2`, the entire common Hamiltonian has rational entries:

```
H = X0+X1+(I-Z0 Z1)/2
    +(73/64)(n_scalar0+n_scalar1+1)
    -(9/64) X_scalar0 X_scalar1
    +(1/4) X_scalar0 Z0.
```

The gauge initial state is `|00>`. The scalar preparation uses the retained
normalized projection of the two-oscillator squeezed vacuum: total occupation
through `d-2`. At `d=4`, relative amplitudes are 1 on `|00>`,
`-1/(18 sqrt(2))` on each of `|20>,|02>`, and `1/18` on `|11>`, followed by
normalization. The existing real-state preparation compiler emits this state
from blank registers; its executed unitary is applied to the blank vector.
The target state is not inserted into the evolution.

The common time step is the following chronological palindrome:

```
gauge_step(h/2), onsite(h/2), bond(h/2), interaction(h),
bond(h/2), onsite(h/2), gauge_step(h/2).
```

Here `gauge_step(s)=hop(s/2),electric(s),hop(s/2)` is itself symmetric,
and the two onsite evolutions commute. Scalar quadrature and onsite bases
are compiled with the existing real-orthogonal source compiler; diagonal
phases use its complete commuting Walsh expansion. For the interaction,
the quadrature eigenbasis is transformed on site 0, and the diagonal phase
is `-h*g*sqrt(2)*q_eigenvalue*(1-2*b0)` on that same register. No controlled
matrix gate is granted separately. All operations resolve to the existing
single/pair proper codes and twelve native M6 generators. Zero coupling
retains the identical sequence of service requests and zero-angle pulses.

The implementation contracts each actual native pulse on every input column
of the full register. It then uses the existing guarded polar correction of
that executed matrix, independently of the target Hamiltonian, and repeats
the executed step by matrix powering. Counts include every repeated pulse
and both proper-code service events per pulse. The numerical guard and
independent full-matrix agreement are rounding diagnostics, not certified
physical service or floating-point error budgets.

### Controlled finite-time limit

This bound concerns the ideal finite-dimensional words with exact target
coefficients. Let `N` denote the onsite occupation operator. Positivity of
the compressed infinite quadratic forms gives

```
Q_d^2 <= Pi Q^2 Pi <= 2N+1 <= (2d-1) I,
Pi(P^2+Q^2)Pi/2 = N+1/2.
```

Consequently the total norm weight of the step's component generators is

```
S <= 3 + (1+2*kappa)(2d-1) + |g| sqrt(4d-2).
```

The 3 is the gauge bound `||H_hop||+||E^2||=2+1`. Both onsite terms together
contribute at most `(1+kappa)(2d-1)`, the bond contributes at most
`kappa(2d-1)`, and the scalar/charge coupling contributes at most
`|g|sqrt(4d-2)`. All coefficients of step duration in the displayed
palindrome are nonnegative. Their weighted norms sum to this same `S`,
including both inner gauge decompositions.

For any product `F(h)=prod_j exp(-i h a_j A_j)` with Hermitian generators,
the multinomial Leibniz expansion bounds its third derivative in operator
norm by `(sum_j |a_j| ||A_j||)^3`: all undifferentiated exponentials are
unitary. Palindromic reversal gives `F(-h)=F(h)^*`; its first derivative is
`-iH`, and differentiating `F(h)F(-h)=I` twice at zero gives
`F''(0)=-H^2`. Thus the product and `exp(-ihH)` agree through order two.
Their integral Taylor remainder norms are individually at most
`S^3 |h|^3/6`, hence

```
||F(h)-exp(-ihH)|| <= S^3 |h|^3/3,
||F(t/n)^n-exp(-itH)|| <= S^3 |t|^3/(3 n^2).
```

The second line follows by unitary telescoping. The operator difference can
be capped at 2, and the associated channel half-diamond error at 1. Tensoring
an arbitrary reference system changes neither norm bound. The implementation
computes an exact rational upper bound for the square root and subsequent
formula; the certificate accepts only integer/Fraction time and coupling.
It rejects inexact floating inputs. At `n=512,t=1,g=1/4`, the half-diamond
upper bounds are below `0.000728` for `d=2` and `0.004184` for `d=4`.

Inserted complete measurements, conditional feedback and record retention
are channels and contract trace distance. Replacing two evolution periods
by their continuous finite-Hamiltonian limits therefore changes the full
joint quantum/classical history by at most `epsilon(t1)+epsilon(t2)`. This
bounds the entire unnormalized outcome instrument and every final outcome
probability. It does not grant the same error on a normalized rare branch
after dividing by its probability. With periods 1 and 1/2 at equal n, the
sum is `9/8` of the first bound. Comparing two differently controlled
experiments charges both experiments' errors.

### Scope

This is one declared finite interacting Hamiltonian and a controlled temporal
limit. It does not establish an occupation or spatial continuum limit,
interacting quantum-field identification, empirical masses or couplings,
responsive geometry, or a calibrated physical clock. In particular the
earlier free Gaussian occupation bound cannot be transferred through this
new interaction. Comparing d=2 with d=4 is a numerical cutoff diagnostic,
not a cutoff certificate. The scalar action, interaction strength, initial
state and proper-code capabilities remain stated design choices. These
choices require no proof of unique selection before implementing them, but
they do not become measured properties of nature by being consistently
implemented.

## Checkpoints, recorded instruments and the same quantum working state

This adapter executes the existing record proof in
`model/observer/DERIVATION.md` on the actual native source/tape carrier and
on the common gauge/matter carrier. It reuses the dual fixed-observable
criterion in `Lean/Dynamics/ProtectedCharge.lean`; it does not claim a new
abstract conservation theorem or a complete `OperationalObserver` inhabitant.

Let the committed record space have basis `|r>` and let W be the
entire current quantum work space. For every r supply a complete instrument
`{K_(r,y)}`, including outcomes of zero probability, with

```
sum_y K_(r,y)^dagger K_(r,y) = I_W.
B_y = sum_r |r><r| tensor K_(r,y).
Lambda_y(X) = B_y X B_y^dagger.
```

This acts on a single work state. It does not draw independent marginals of
separately simulated components. For any off-diagonal record block,
`Lambda_y(X)_(r,s)=K_(r,y) X_(r,s) K_(s,y)^dagger`. Outcome y is retained as
an orthogonal classical output block; no branch normalization or selection
occurs. Completeness gives `sum_y B_y^dagger B_y=I`, and hence a CPTP
instrument, including on a reference-entangled input.

Put `Q_r=|r><r| tensor I_W`. Every B_y commutes with Q_r, so

```
Lambda*(Q_r) = sum_y B_y^dagger Q_r B_y
            = Q_r sum_y B_y^dagger B_y = Q_r.
```

Consequently every old-record observable `sum_r a_r Q_r` has exactly
preserved expectation for every input operator. This is the existing
dual-channel criterion, bound to the executed common-source instrument.
It also proves repeatability of a read of a committed old record. The proof
is unaffected by tensoring an arbitrary untouched reference. It does not
assert that every work observable or every off-diagonal quantum record
coherence is fixed. Individual unnormalized outcome branches need not
preserve the old probabilities separately; their complete sum does.

There are two useful and different readback maps:

```
E_block(X) = sum_r Q_r X Q_r,
E_public(X) = sum_r Tr(Q_r X Q_r)/dim(W) * Q_r.
```

`E_block` is the pinching conditional expectation onto
`direct_sum_r B(W)`. It preserves all within-record quantum operators,
including off-diagonal work matrix units, fixes that whole algebra, is
idempotent and obeys its bimodule law. Its fixed algebra is generally
noncommutative. `E_public` is the conditional expectation onto the
commutative algebra spanned by the old-record projectors. It is unital,
completely positive, trace preserving and idempotent, fixes each Q_r, and
satisfies `E_public E_block=E_public`. The latter is the relevant classical
self-readback algebra. The instrument's preservation of Q_r does not imply
that it fixes the larger block algebra pointwise.

### Native checkpoint attachment and the blank-sector control

Use the actual native finite space

```
H_cap = H_0 direct_sum H_1 direct_sum H_2,
dim(H_0,H_1,H_2) = (78,468,2808),
T(X)=A_0 X A_0^dagger+A_1 X A_1^dagger+P_2 X P_2.
```

The code reconstructs every native write column from `native.tick`.
Nonfull columns each have two amplitudes `1/sqrt(2)` with disjoint output
supports; full columns are the original identity halt. It checks all
78+468+2808 columns. The separate phase Kraus branches give the full Gram
identity without incorrectly retaining cross-phase coherences.

On `K=H_1 direct_sum H_2`, each old first-letter projector is intertwined
by A_1 and P_2. Therefore the restricted channel T_K is CPTP and its dual
fixes every first-letter projector on the 3276-dimensional checkpoint
space, including unreachable first-letter/load combinations. Sparse exact
Gaussian-rational execution verifies arbitrary operator matrix units and
the original reachable coherence witness: the first-letter off-diagonal
entry `((0,1),(0,0))` changes `1/4 -> 1/8`, while all first-letter diagonal
probabilities are unchanged.

Writing is a distinct transition into a checkpoint. For the blank input
with occupations `(1,1,0,...,0)`, extend the first-letter `(0,0)` projector
by zero on H_0. Its expectation is zero before a tick and exactly 1/2
after it. It is therefore false that `T*-I` annihilates that projector on
the entire blank-and-written space. `OperationalObserver.record_durable`
requires a fixed-observable law for its declared persistence generator;
that generator must describe the postcommit persistence/control regime if
these raw record projectors are used. Writing transitions must be attached
separately. This package does not silently identify the writing generator
with the persistence generator, and it does not instantiate all other
fields of `OperationalObserver` or its cross-observer extension.

## Exact hardware refinement and approximate decoded physics

The declared coherent-source construction has an exact inventory tower:
retain every old privately owned processor and buffer once, add independent
factors, embed old observables by `X -> X tensor I`, and restrict states by
partial trace. Lift an old operation by acting identically on added factors.
These exact maps support the corresponding state/operation compatibility
identities. Adding duplicate copies of an unknown old quantum buffer would
violate this construction.

An occupation cutoff is a different map. The isometry V_d from a finite
oscillator code into Fock space gives the embedding
`X -> V_d X V_d^dagger`, whose image of identity is the cutoff projector,
not the ambient identity. The reverse compression is unital as an operator
system map but need not preserve multiplication. Normalized cutoff vacua
and compressed Hamiltonians need not obey exact state or generator
naturality across cutoffs. Thus the quantitative occupation bounds cannot
be relabeled as the exact unital-star-algebra tower fields in
`Lean/Tower/ConsensusTower.lean`.

The operational comparison can instead use specified finite experiments.
For each decoded preparation/intervention/readout protocol a and regulator
d, compare complete record laws p_d(.|a) and p(.|a), with an established
bound

```
TV(p_d(.|a), p(.|a)) <= epsilon_d(a),   epsilon_d(a) -> 0.
```

When a state-specific complete-output trace-distance bound is available,
every common subsequent CPTP continuation and readout contracts it. The
two-site occupation certificate has that scope on its declared prepared
histories; it is not an all-input diamond-norm estimate. If the continuation
itself changes with d, its own comparison error must be added. A common
finite experiment menu permits a uniform maximum of its finitely many
bounds. It does not produce a bound over all physical interventions or a
spatial continuum from occupation refinement alone. The hardware tower,
decoded approximation, and empirical identification remain distinct typed
objects, with explicit maps required between them.

## Exact detector response for the common interacting finite history

The following certificate concerns the same finite source construction with
onsite scalar dimension `d=2`. It proves a detector-response difference for
the declared model and gives a bound for ideal native source words. It does
not use floating-point predictions as premises, identify a physical
apparatus, or establish an occupation or spatial continuum limit.

### Exact common generator and preparation

Index the 16 basis vectors by little-endian bits `(b0,b1,n0,n1)`.
The first two bits are right occupations of the charge `(+1,-1)` CAR
channels, with the same physical isometry and electric flux `E=-b0+b1`.
The final two bits are the scalar occupations at the two declared sites.

At `d=2`, the projected scalar coordinate is `Q=X/sqrt(2)`. The projected
squares are the matrix elements of the original oscillator squares, rather
than the squares of the truncated coordinate/momentum matrices. With
`kappa=9/32`, the two onsite terms are
`(73/64)(n0+n1+1)`. Consequently the complete declared Hamiltonian is

```text
H = X0+X1+(-b0+b1)^2 +(73/64)(n0+n1+1)
    -(9/64)X2 X3 +(1/4)X2 Z0.
```

This is a real symmetric rational matrix. It is independently compared
with the scalar quadrature/physical-square and gauge-matter assemblers.
The normalized `d=2` preparation is exactly `|0000>`: only the scalar
total-occupation-zero term survives the declared preparation projection.
It is not identified with the full infinite interacting vacuum.

The five term-norm upper bounds are, respectively,

```text
2, 1, 219/64, 9/64, 1/4,
```

whose sum is `S=109/16`. Thus `||H||<=S`. The norm bound includes the
scalar boundary correction and the noncommuting scalar/matter interaction.

### Complete retained instrument

Let `U1=exp(-i H)` and `U2=exp(-i H/2)`. The first binary read has
projectors `P_a` onto `b0=a`, for `a=0,1`. The committed first record
controls `F_0=I` and, when feedback is enabled,
`F_1=(3I-4iX0)/5`. With feedback disabled both `F_a` are identity.
The second read concerns the sign of `Q0`: `B_1=(I+X2)/2` and
`B_0=(I-X2)/2`.

For each setting, all four unnormalized branch vectors are

```text
v_(a,b) = B_b U2 F_a P_a U1 |0000>.
```

They retain the first record even when feedback changes the working
occupation bit. Neither conditioning on a successful first result nor
coherent recombination of different record values is used. The
unconditional positive-coordinate probability is
`p_plus=sum_a ||v_(a,1)||^2`. The observable contrast is its value with
feedback minus its value without feedback.

### Rational Taylor enclosure

Define `T_m(t)=sum_(k=0)^m (-itH)^k/k!`. The integral Taylor remainder
for the selfadjoint finite generator gives

```text
||T_m(t)-exp(-itH)|| <= epsilon_t
epsilon_t = (S t)^(m+1)/(m+1)!.
```

There is no `exp(S t)` factor: the integral remainder contains the
unitary `exp(-isH)`, of norm one. The production calculation uses only
exact Fraction pairs for Gaussian-rational vector arithmetic. Its norm
bound is an operator bound, so when the input is unnormalized the vector
error includes the input norm.

Writing `epsilon1=epsilon_1` and `epsilon2=epsilon_(1/2)`, replace both
exponentials by their degree-`m` polynomials. Projectors are contractions
and feedback is unitary. The error of each branch is at most

```text
eta = epsilon2 + (1+epsilon2)*epsilon1.
```

Indeed, telescope at the first and second evolution, use
`||T_m(1/2)||<=1+epsilon2`, and use the unit norm of the exact initial
state. Each exact branch has norm at most one, hence its squared norm
differs from the polynomial squared norm by at most `2 eta+eta²`.
Summing the two first-record branches gives probability error at most
`4 eta+2 eta²`. Subtracting the two configurations gives contrast error
at most `8 eta+4 eta²`. These conservative estimates use no normalization
of approximate branches.

At `m=64`, the unconditional probability radius is less than
`7.091e-37`. Exact outward rounding to denominator `10^15` gives

```text
feedback enabled:  [0.370697726608684, 0.370697726608685],
feedback disabled: [0.402359206973491, 0.402359206973492],
contrast:          [-0.031661480364808, -0.031661480364807].
```

The implementation retains exact rational endpoint strings. The decimal
intervals here are those same outward rational bounds, not rounded
floating-point central values. A separate numerical exponential check is
corroborating evidence, not part of the proof. The contrast excludes zero
for this complete continuously evolving finite-model instrument.

### Ideal source-word response

The common scalar/gauge/matter word is a symmetric product whose generator
sum is exactly `H`. Its finite splitting certificate bounds the
half-diamond error in each evolution segment. Inserting the same
projective instruments and feedback channels cannot enlarge trace
distance. The unconditional record-probability error for one complete
configuration is therefore at most the sum `delta1+delta2` of those two
segment channel bounds. The feedback contrast has additional error at
most `2(delta1+delta2)`.

At 512 source steps per segment, combining those rational bounds with the
exact Taylor enclosure gives

```text
ideal source-word contrast:
[-0.033297300667839, -0.030025660061775].
```

This interval also excludes zero. It certifies the ideal source word
implementing the declared finite action, not the numerical contraction or
the physical availability/accuracy of the source services. The certificate
names those exclusions. Coarse splitting or a deliberately inadequate
Taylor degree instead yields an interval containing zero; the evaluator
does not manufacture a resolved signal in those cases.

This result proves a common interacting observer process can produce a
nonzero record-conditioned field response with an exact quantitative
bound. It supplies no comparison with independent measurements in nature
and no transfer to another occupation cutoff or a spatial continuum.

## Native observation, interacting quantum response and native feedback

The common experiment uses one declared state space for the native loads,
the gauge/scalar quantum work registers and the committed history. The native
integer-load repair is its existing instrument primitive. Gauge, scalar,
pointer and conditional preparation/control operations are emitted and
executed through the M6 proper-code source compiler. This joining result
does not assert that the entire native M78 instrument has itself been
compiled into those M6 instructions, or that its prepared state is the
unconstrained tracial A3 optimizer.

Let `L_r^E` be the first native current-record Kraus operator, where
`r=(c1,j)` and `E=(0,1)`. Let the native current control the actual emitted
gauge pulse

```
V_j = exp(-i*j*theta*X_0),    theta=atan2(4,3),    j in {-1,0,+1}.
```

The common quantum source prepares its declared finite scalar vacuum and
gauge state. It then applies V_j, its first interacting evolution U_1,
the species-0 right-occupation instrument P_a, the recorded-a-controlled
hopping F_a, its second interacting evolution U_2 and the scalar-coordinate
sign instrument R_b. Here a,b are ordinary committed bits, with b=1
meaning positive scalar coordinate. Both quantum outcomes are retained.
The quantum history operator, apart from its initial V_j, is

```
K_(a,b) = R_b U_2 F_a P_a U_1.
sum_(a,b) K_(a,b)^dagger K_(a,b) = I.
```

All matrices act on the same gauge/scalar register. In particular, P_a
changes the actual state passed to U_2; it is not an independent draw from
a separately evolved marginal. The first classical record a survives F_a
and the second read, while the quantum state can change.

Only the committed final b enters the native actuator:
`e(b)=(0,5)` if b=1, and E otherwise. For the second native record s,
the entire experiment has Kraus operators

```
M_(r,a,b,s) = (L_s^e(b) L_r^E) tensor (K_(a,b) V_j).
```

Write the outcome `(r,a,b,s)` into orthogonal classical record blocks.
Summing M^dagger M first over s, then (a,b), then r gives identity.
This is one CPTP instrument on arbitrary jointly prepared native/quantum
inputs, including entangled inputs. The first native record is classical
before it controls V_j. This permits blockwise conditional composition;
it is not a promise to preserve arbitrary quantum coherence between
distinct committed record values.

### The executed prepared law and its retained quantum outputs

The reported benchmark prepares the native system uniformly on the 66
binary total-two configurations and prepares the independently declared
quantum vector psi. Each fixed placement coin contributes a native Kraus
amplitude 1/sqrt(2). Therefore for each initial native x, first coin c1,
quantum records (a,b) and second coin c2, the complete history probability is

```
p(a,b|j)/264,    p(a,b|j) = ||K_(a,b) V_j psi||^2.
```

There are exactly `66*2*4*2=1056` rows. All rows are retained, including
zero-probability outcomes. No parity identity or fair-current-zero detector
law is assumed: the current-zero quantum response need not be unbiased.
The joining code checks every reported p against the norm of the actual
retained unnormalized quantum vector before enumerating histories.

Each row identifies its final native basis state, complete classical
records and the referenced quantum vector `v_(j,a,b)=K_(a,b)V_j psi`.
Its quantum terminal density contribution is exactly
`|v_(j,a,b)><v_(j,a,b)|/264`, tensored with the stated native state and
record projectors. The twelve vectors are stored once, indexed by j and
(a,b); retaining their references preserves these terminal density
operators rather than reducing the experiment to a probability lookup.
The initial-source and final-source integer indexes are audit labels,
not extra inputs available to the observer's bit-only controller. For
general coherent native inputs the displayed Kraus formula supplies the
channel; the 1056-row enumeration specifically describes this prepared
classical native mixture.

Completeness preserves the original first-native-record probabilities:
`14/33,14/33,5/66,5/66` for `(0,0),(1,0),(0,1),(1,-1)` respectively.
Later states and seam decisions depend on the retained quantum record;
no later event rewrites those old native values. The checkpoint adapter's
fixed-projector proof applies to this record-controlled continuation.

### Controls isolate the effect of the recorded scalar bit

The fixed-seam control always applies E to the actual native state after
the first repair. It retains all prior native and quantum records and
terminal quantum states. The stronger control draws an independent
actuator bit with the actual experiment's overall b=1 frequency. Its
choice is independent of the native state, first native record and both
quantum records. Thus it has the same seam frequencies as the actual
record-controlled protocol. Both placement-coin branches and both donor
bit choices remain in the calculation. A difference in second native
current or final port-0 load measures the effect of their correlation
with the actual record, rather than only a change in action frequency.
The code computes these differences and does not require them to be
nonzero for arbitrary supplied model parameters.

### Transfer of the same finite temporal comparison

For each admitted j, the source experiment supplies a half-diamond bound
epsilon_j comparing its ideal split quantum instrument with ideal
continuous evolution of the same finite Hamiltonian, preparation pulse
and intermediate measurements. Taking `epsilon=max_j epsilon_j` bounds
the corresponding block-controlled channel on the committed first-record
algebra. The unchanged exact native instrument before it and the complete
native continuation after it are CPTP and contract the bound. Thus the
complete joint instrument has half-diamond error at most epsilon, and
the full recorded benchmark law has total-variation error at most epsilon.
Since the binary-shell second native current lies in [-1,1], its mean
error is at most 2 epsilon. The code transfers only a supplied complete
set of per-current certificates; it does not invent one for a test
fixture lacking those certificates.

This is temporal approximation of a specified finite interacting model.
Floating compilation and physical service errors require additional
terms. The older free two-oscillator occupation estimate does not cover
the present gauge/scalar interaction, and no spatial continuum theorem
or empirical identification follows from this composition.

## Compatible theorem use and open physical connections

| Existing result | Attachment here | Boundary |
| --- | --- | --- |
| Proper-code source control | Actual M6 pulse contraction on all work columns, including electric, matter and scalar factors | Coherent load/export, blanks, storage and pointer readout remain apparatus services |
| CAR and operator Gauss | Explicit flux-dressed isometry; every action and read preserves its physical sector | Opposite-charge four-state code, not the full observed particle census |
| A3 and conditional expectation | Faithful tracial reference with block and public-record algebras distinguished | Does not select the prepared pure state or prove one global repair channel is its stationary optimizer |
| Protected records | Fixed old-record projectors on committed checkpoints and complete controlled instruments | Blank-to-written transitions are separate; no complete Lean `OperationalObserver` inhabitant is asserted |
| Hardware refinement | Old privately owned factors embed by tensoring identities | Oscillator cutoff embeddings are not exact unital algebra/state/generator tower maps |
| Finite temporal approximation | Rational all-input bounds for this same interacting Hamiltonian, with arbitrary spectator references | Numerical/service errors remain separate; the new occupation analysis has a state-specific scope |
| Interacting occupation refinement | [Moment/Duhamel estimate and strong sharp-readout limit](OCCUPATION_REFINEMENT.md) for this same two-site history | No finite-cutoff sharp-readout rate or spatial field limit |
| Observer representation and modular flow | [GNS and modular readouts](OBSERVER_READOUTS.md) of the actual observation state; emitted clock probe retains all records | Quantum algebra and access supplied; full conditioned scalar states may be nonfaithful; physical clock calibration open |
| Quantum regional restriction | CP conditional expectations and Hamiltonian-derived interaction support on the same carrier | Supplied factorization; no multiplicative restriction, metric or growing spatial family |

The [model–theorem premise map](../PREMISES.md) retains the broader
conditional results. They must be applied with compatible hypotheses;
classical and noncommuting repair, exact and approximate Markov structure,
terminating and continuing schedules, and different native record choices
are not one simultaneous set of assumptions.

The next physical obligation is a common refining family of these
**observer statistics**, with controlled errors and identified preparations,
interventions, states and clocks. It need not make a fundamental spacetime
or a fundamental QFT. It must make their observed joint statistics and
responses effective descriptions of this same process. Two occupation
cutoffs on a fixed two-site carrier do not supply spatial refinement.
Nonabelian/chiral matter, the physical coupling map, dynamical stress and
geometry, and independent data comparison remain open. The free scalar
Gaussian occupation bound cannot be reused after adding this interaction.
Neither a tracial reference's trivial modular flow nor a solver counter
is a calibrated observer clock. P/N closure assumptions do not determine
the declared coupling g or the action/preparation in this packet.
