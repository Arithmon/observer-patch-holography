# Coherent continuation: native attachment and nonselection

The executable [exact experiment](coherent_continuation.py) attaches a
reflection walk to the actual N1 repair branches on the 66-state binary
total-two shell. It derives a wave recurrence for compressed **amplitudes**
and gives two exact countercontrols against selecting that wave law from
the repair instrument alone. This is a law-selection experiment, separate
from the observer that commits a classical tape after each event. The finite
[compiled density walk](#source-compilation-of-the-finite-density-walk) realizes
its density-sector operation with the existing source controls.

```sh
python3 -m model.observer.coherent_continuation --check model/observer/results/coherent_continuation.json
python3 -m pytest -q -W error model/observer/test_coherent_continuation.py
```

## Actual native branches

Let `x` run over the binary shell, `e` over the thirty pinned source seams,
and `c` over the two placement coins. There are 3,960 arc basis vectors
`|x,e,c>`. Apply the existing native rule to obtain `y=F_c^e(x)` and current
`j`. Define the reversal

```text
S|x,e,c> = |y,e,c'>,   c'=c if j=0, and c'=1-c otherwise.
```

The executable checks every arc: `S²=I`. There are 300 exchanged pairs and
3,360 fixed arcs. The output record `(y,e,c,j)` and the reversed arc are
bijective: applying the reverse event recovers `x` and `-j`. This is a
coherent re-encoding, not erasure of an observed classical record.

The uniform coherent event embedding is
`A|x> = sum_(e,c)|x,e,c>/sqrt(60)`. Its compression obeys `A†SA=K`, where
`K` is exactly the independently enumerated symmetric native transition
kernel under uniform seam attempts and fair coins. All entries are rational.
The implementation stores arc amplitudes multiplied by `sqrt(60)`; it
handles real rational columns, with complex and reference columns extending
by linearity. The norm includes the compensating factor `1/60`.

The binary restriction matters. On the full 78-state carrier the event
`(2,0,...)->(1,1,...)` on seam `(0,1)` has no native reverse event. Its
transition kernel is asymmetric, whereas `A†SA` for an involutive unitary
`S` is Hermitian. This excludes this same-domain reversal construction on
that carrier; it does not exclude a larger history-space dilation.

## Imported theorem and source-specific consequence

Set `Pi=AA†` and `W=S(2Pi-I)`. The actual objects satisfy the coisometry and
involution hypotheses of Higuchi, Segawa and Suzuki's
[abstract quantum-walk spectral mapping theorem](https://arxiv.org/abs/1506.06457):
their `d_A` is `A†`, their discriminant is `K`, and their unitary is `W`.
The equivalent product-of-reflections construction in
[Szegedy](https://arxiv.org/abs/quant-ph/0401053) corresponds to `W²`;
its eigenphase has an extra factor of two. No new general walk theorem is
claimed here.

On the embedded sector, the standard reflection identity gives
`A†W^n A=T_n(K)`, with `T_n` the Chebyshev polynomial. For the centered native
load matrix `X_(x,p)=x_p-1/6`, direct enumeration gives

```text
K X = X (I-L/60),
U_(n+1) - 2 U_n + U_(n-1) = - U_n L/30,
U_n = A† W^n A X,
U_0 = X,  U_1 = X (I-L/60).
```

Here `L` is the twelve-port source graph Laplacian. A Laplacian mode has
`cos(omega)=1-lambda/60`. The parameter `n` counts applications of the chosen
coherent walk. It is neither the committed observer's repair counter nor a
calibrated physical clock. Squared amplitudes and measured mean loads do
not obey this linear recurrence. The existing forgotten-record mean-load
law remains diffusive.

The new result is this exact attachment to the native branch and record
dictionary, together with its selection obstructions. Existing
[quantum transport](../../code/m1_quantum_transport/) and
[operational clocks](../../code/m1_operational_clocks/) already supply
walks and interferometric readings for their declared dynamics.

## Same initial event, different continuations

For any phase `z` of modulus one, define
`W_z=S[Pi+z(I-Pi)]`. Every choice has the **same complete first coherent
event**, `W_z A=SA`. Nevertheless,

```text
A† W_z² A = K² + z(I-K²).
z=+1: W_z²=I;             z=-1: A† W_z² A=2K²-I.
```

The implementation executes the two real endpoints with exact fractions.
An added Ramsey control compares the walk arm with an identity arm and
retains both outputs `(A f +/- W_z^n A f)/2`. For its normalized basis input
at index zero, the two-step plus probabilities are `1` and `113/150`,
respectively; the latter return **amplitude** is `38/75`. The minus outcomes
complete each probability to one. These are specified interferometric
services, not an emitted source apparatus or a physical clock.

A second control changes a scalar phase on one committed record label:
`eta=-1` on seam `(0,1)` with current `+1`, and `+1` elsewhere. This leaves
each committed C CP branch unchanged, including on arbitrary references;
the E label also determines that current, so the same statement holds for E.
For the phase-twisted embedding `A_eta`, however,

```text
D = A_eta† S A_eta,
Tr(D³)-Tr(K³) = -1/900.
```

Ten undirected kernel transitions change sign. The trace difference rules
out merely a basis change of the compressed kernel. Thus even a complete
committed quantum instrument does not determine its coherent continuation.

## Premises and research consequence

The binary preparation, uniform coherent seam/coin control, uncommitted
record access, reversible workspace re-encoding, reflection, continuation
phase and complete interferometric readout are supplied. There is no
derivation of a spatial refinement family, relativistic dispersion,
calibrated clock, physical service accuracy or empirical fit.

The next source-selection task must constrain the cross-record phases and
continuation action, or exhibit an observable insensitive to that freedom.
Recovering a wave recurrence after choosing the reflection does not discharge
that task. The exact countercontrols are acceptance tests for any proposed
selection argument.

## Selection under fixed structural constraints

The extended exact receipt checks the two real continuation endpoints under
all sixty proper rotations of the pinned carrier, including orientation
reversal of an edge's coin convention. It verifies 237,600 arc covariance
identities per reversal. The two endpoints share conservative one-seam
support, real orthogonality, a fixed uniform coherent ket and the same
algebraic antiunitary reversal `Theta=S conjugation`. Even an additional
full event-permutation coin symmetry leaves identity and the reflection;
this stronger symmetry is a hypothetical restriction, not a source premise.
Both preserve amplitude disagreement with the event-uniform subspace.

An all-seam moving-arc sign control removes the distinguished seam from the
earlier phase counterexample. Its 600 negative moving arcs preserve every
committed C/E CP branch and proper covariance, while changing the compressed
trace cube by `-1/90`. It has uniform-return outcomes `529/1089,560/1089`
and fixed-Ramsey outcomes `28/33,5/33`. Unlike the endpoint pair, this signed
control does not fix the uniform coherent ket or share the first coherent
event. The [source-selection guide](../SOURCE_SELECTION.md) states the
three primitive obligations and the appropriate interface equivalence.

The [source-executed record test](source_selection.py) separately realizes
an observable phase distinction with existing C/E source words. It compares
coherent record access with commit-first access; it is not a source-pulse
realization of the entire 3,960-arc walk. Replay it with
`python3 -m model.observer.source_selection --check
model/observer/results/source_selection.json`.


## Source compilation of the finite density walk

The same source Laplacian admits a smaller representation of this precise
density sector. On twelve one-particle port addresses, take sixty events
`(e,sigma)`: `sigma=0` fixes the port and `sigma=1` exchanges the endpoints
of seam `e`. This is the native one-particle event under the reversible
coin relabelling `sigma=1-c` at the first endpoint and `sigma=c` elsewhere.
All 720 native branches and inverse coins are checked. The event embedding
`A_12` has uniform amplitudes and `A_12* S_12 A_12=K_12=I-L/60`.

For the original 66-state total-two source, set
`F[x,p]=(x_p-1/6)/sqrt(10)`. Exact enumeration gives

```text
F*F = I-J/12,                    K_66 F = F K_12.
```

Thus `F` is isometric on the eleven-dimensional mean-zero port space.
The two pairs `(A_12 f,S_12 A_12 f)` and
`(A_66 Ff,S_66 A_66 Ff)` have the same block Gram matrix
`[[I,K_12],[K_12,I]]`. Sending the first pair to the second therefore
extends to an isometry of their cyclic subspaces. The identity
`W A=SA`, `W SA=2SAK-A` shows that it intertwines their walks and every
density compression. This proves a finite density-sector representation;
it does not implement the full 3,960-arc unitary or its higher-configuration
observables.

[`density_walk_program`](coherent_continuation.py) emits the following word
using the existing [source compiler](source_scalar/README.md). Bits `0..3`
encode the port, `4..9` encode the event, and five initially blank bits
supply reversible scratch. A source word `U_60` prepares the uniform first
sixty event addresses. The event reflection is
`U_60 diag(1,-1,...,-1) U_60*`. A controlled permutation exchanges each
seam's endpoint addresses only for its odd event. Reversible AND ladders
and Gray paths compile this permutation using the compiler's two-dimensional
and four-dimensional proper-code services; no new logical matrix is a
primitive.

One walk has **114,910 native pulses**: 73,414 in the reflection, 1,010
forty-pulse Toffoli routines and 548 two-pulse bit flips. Every phase and
zero-angle request is charged, with 114,910 loads and 114,910 exports.
Preparing the event register costs a separate 36,288 pulses. The complete
chronological word is generated on demand and its SHA-256 is retained;
the uniform-state angles use exact integer multiplicities and the pinned
arbitrary-precision backend at eighty decimal digits. Emitted angles are
rounded to multiples of `2^-44`. No expanded pulse
array is stored in the repository.

All 32,768 reversible basis inputs are executed exactly. Every one of the
1,024 clean-scratch address inputs returns clean scratch. Of these, 720
addresses are valid and 304 are retained padding. The shift fixes padding;
the reflection acts on the event register at every site and equals `-I`
on event padding. This specifies a unitary extension, including the padding
block, without postselection. The declared proper-code export channel
retains each processor-complement failure branch; source generators preserve
the valid code exactly, so these failure probabilities are zero for the
stated valid inputs. Invalid processor inputs have explicit failure outputs.
Application-address padding and processor failure are different branches.

Each routine is executed pulse by pulse from the full six-level source
exponentials, on all its input columns. The 15-qubit replay reuses those raw
executed matrices and retains every output component, including scratch
and padding. It checks a complex reference input and an event-prepared
two-step density input. This is **factorized verified execution**, not
literal replay of every expanded pulse on the entire 15-qubit array.
No polar projection, target-matrix substitution or output normalization is
used. An independent test compares literal pulses and factor reuse on a
complete smaller register, including overlapping controls.

The numerical guards bound the reflection's Frobenius error by `10^-8`
and each Toffoli/flip error by `10^-10`; these also bound operator errors.
For one walk their sum is `s=829/5,000,000,000`. Multiplicative error
accounting gives `delta<=s/(1-s)` on the complete register, including any
spectator reference. The induced unnormalized output trace-norm error is
at most `2 delta+delta^2`. Event preparation has a `10^-8` first-column
guard; its two-step bound uses `s=10^-8+2(829/5,000,000,000)` and requires
an initially blank event register and scratch. These are checked floating
bounds, not interval or physical precision certificates. If physical pulse
implementation errors are at most `epsilon` in operator norm, `P` pulses
add `P epsilon/(1-P epsilon)` when `P epsilon<1`; load/export errors need
their own allowance. The existing source-pulse and proper-code apparatus,
blank scratch, coherent storage and analogue control remain explicit.

This constructs the finite wave service for the density representation.
Its [incident-clock attachment](../relative_observer_geometry/README.md#active-coherent-waves-in-a-native-incident-clock)
separates a timing carrier from a protected wave register. A source-pulse
Ramsey controller and a compiler with vanishing error along a growing
refinement are not executed here. The identity/reflection and record-phase
controls continue to show that committed repair alone does not select this
coherent law. No physical vacuum, observed Lorentzian field or laboratory
clock identification follows from compiling the finite service.

## A finite edge-generator representation from native repair

The native twelve-port density Laplacian has an exact finite-group form.
Let `H` be the sixty proper rotations of the recovered carrier and `U_g`
their vertex permutations. Exactly one of its two twelve-element order-five
conjugacy classes satisfies

```text
L = 6 I - (1/2) sum_(g in C) U_g,
K_12 = I-L/60 = (9/10) I + (1/120) sum_(g in C) U_g.
```

`C` is inverse symmetric, so this is a positive central Cayley Laplacian
represented on the actual density module. The complete native kernel
retains `K_66 F=F K_12`. The other order-five class has the same spectrum
but predicts zero one-attempt transfer across a native seam, where the
source predicts `1/60`. Thus the local record dictionary distinguishes
these isospectral choices.

There are four orbits of ordered vertex pairs, so the representation's
commutant has dimension four. The four distinct eigenvalues of `L` have
multiplicities `1,3,5,3`, exhausting that commutant: each eigenspace is an
irreducible representation occurring once. The group image is therefore
`M_1 ⊕ M_3 ⊕ M_5 ⊕ M_3`, of dimension 44; independent exact matrix rank
checks this count. Its sector projectors have ranks `d_R`. This realizes
the generator and rank clauses of the
[finite edge-law theorem](../../paper/tex_fragments/PAPER.tex) on this
specific representation, without importing a compact gauge group.

Assigning `rho_t=exp(-t L)/Tr(exp(-t L))` then gives
`p_R=d_R exp(-t lambda_R)/Z`. Poissonizing the native attempt count gives
`E[K_12^N]=exp(-t L)` for `N~Poisson(60t)`, but an averaged density
transition operator is not a physical preparation of `rho_t`. Normalized
trace on this twelve-dimensional density representation gives sector
weights `1/12,1/4,5/12,1/4`; this is not the stationary occupation
distribution or density covariance of the 66-state carrier.
The missing four-dimensional irreducible block also prevents identifying
this 44-dimensional image with the full 60-dimensional group algebra.
No native gauge collar, Gibbs preparation, complete group sewing, compact
group refinement or worldsheet limit is established by this representation.
