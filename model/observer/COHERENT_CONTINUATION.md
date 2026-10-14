# Coherent continuation: native attachment and nonselection

The executable [exact experiment](coherent_continuation.py) attaches a
reflection walk to the actual N1 repair branches on the 66-state binary
total-two shell. It derives a wave recurrence for compressed **amplitudes**
and gives two exact countercontrols against selecting that wave law from
the repair instrument alone. This is a law-selection experiment, separate
from the observer that commits a classical tape after each event.

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
