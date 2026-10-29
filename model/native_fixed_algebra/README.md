# Native conservation and a common-reference contraction

The native binary two-load carrier has a stronger conservation boundary than
the [branch-commutant certificate](../completion/defect_algebra.py) alone.
For positive all-seam exploration and fixed positive endpoint weights, every
source-only Heisenberg conserved observable is scalar on H66. This conclusion
covers arbitrary classical record feedback when conservation is required on
every admitted source/controller input. It is not merely constancy of an
expectation in one stationary state.

The same source maps have an explicit positive gap in the weighted
Hilbert--Schmidt norm. This gives a quantitative test for nearly conserved
observables and convergence under arbitrary **exogenous** varying selectors
with the same endpoint weights. It does not establish mixing of a growing
record-feedback process or a physical particle, mass or time calibration.
The finite matrix arguments and executable certificates here are not new
Lean theorems.

## Carrier, channel and reference

Fix the existing H66 occupation basis and positive endpoint weights `theta`.
The declared reference is `pi_ij=theta_i theta_j/Z`, selected independently
of the update. On seam `(a,b)`, the endpoint heat-bath coin puts the load at
`a` with probability `theta_a/(theta_a+theta_b)`. C and E retain their complete
raw partitions. Each native branch is `K=sqrt(w) M`, with `M` a partial
injection. Every seam's C and E channel preserves the same faithful `pi`.

An attempt chooses any positive seam distribution, any C/E mixture separately
on each seam, and a positive repair probability; its complement is an identity
hold. These choices are independent of the source conditional on the declared
classical controller state. The quantitative uniform constants below use the
bounded [local family](../local_repair_family/README.md). The conservation
argument only needs positivity and all-seam support, with fixed `theta`.

Write the Heisenberg pullback as `Phi(X)=sum_b K_b* X K_b`, and use
`<X,Y>_pi=Tr(pi X*Y)`. Forming this pullback ignores the new output flag only
for this observable; the complete physical instrument and its archive are
not erased or assumed independent. `faithful_channel_certificate` checks
completeness and faithful invariance exactly from the actual raw maps.

## Every source-only fixed observable is scalar

For any matrix `X`, completeness and `Phi_*(pi)=pi` give the identity

```text
Re <X,(I-Phi)X>_pi = (1/2) sum_b ||X K_b-K_b X||_pi^2.
```

For `Phi(X)=X`, the right side vanishes. Faithfulness makes every
commutator zero. Since `Phi` preserves adjoints, applying the same argument
to `X*` also gives commutation with every `K_b*`. The existing exact effect
separation and transition-tree certificate shows that the supported C
branches generate the full H66 matrix algebra. E refines C and contains its
integer branch sums. Consequently any positive per-seam choice of C, E or
their mixture has the same scalar commutant. No numerical rank test is used.

For a classical controller state `h`, denote its conditional source pullback
by `Phi_h`. A source-only observable on the full instrument is `X tensor I`
on source and controller/archive. Its Heisenberg conservation on **all**
admitted source/controller inputs requires `Phi_h(X)=X` at every active `h`.
The preceding argument applies at each such state, without any stationary
density on the complete growing tape. Thus classical record feedback alone
does not provide a nonconstant universally conserved source observable.

An expectation in the single invariant state `pi` is constant for *every*
observable, including those whose Heisenberg pullback changes. The tests
separate these two assertions. Source/controller-correlated charges such as
the [co-moving relational charge](../conserved_sector_action/README.md),
history observables, restricted input classes, active quantum recovery,
coherent feedback, different carriers and refining limits are outside this
source-only theorem. H78 has a nonfaithful stationary support and is not
included in this faithful-reference argument. Removing all seams incident
to one port restores its nonconstant occupation charge, providing a positive
control for the all-seam hypothesis.

## An exact weighted spectral bound

The single-seam channels admit an additional simplification. In the matrix
unit basis, each map splits into blocks of size one or two. The nontrivial
two-unit blocks are the endpoint heat-bath block

```text
[[p_a,p_b], [p_a,p_b]],
```

or the corresponding relabeling, with reference weights in ratio
`pi_a/pi_b=p_a/p_b`. Multiplying rows by their reference weights gives a
symmetric positive matrix of determinant zero. Its weighted identity-minus-
map matrix is also positive. The one-unit blocks multiply by a number in
`[0,1]`. This description includes spectator coherences: when both sides of
a matrix unit move together they form the same two-unit block; otherwise
the common-record rule retains a stationary unit or kills it.

Therefore each seam's C/E map is a positive self-adjoint contraction in
`<.,.>_pi`. Here positive quadratic form is a Hilbert-space property, distinct
from the map's complete positivity. `seam_GNS_certificate` verifies the
weighted symmetry and every one- and two-unit principal minor directly.
For the canonical biased law this checks 261,360 units across all 60 maps.

For the fair uniform C comparison channel, reuse the exact population
minorization `P^20 >= beta 1 nu^T` from the
[source preparation proof](../local_repair_family/mixing.py), where

```text
beta = 142390556543804527497099183523 / 743847338880000000000000000000.
g0   = min(beta/20, 1/6) = beta/20 > 0.00957122.
```

The minorization contracts zero-mass population differences by `1-beta`.
Symmetry then bounds every nonconstant population eigenvalue away from one
by at least `beta/20`. The actual off-diagonal C blocks are symmetric,
entrywise nonnegative and have row sums at most `5/6`, giving the other
term. These are rational checks, not fitted eigenvalues or an imported
exclusion-process gap theorem.

Replacing C by E on any individual seam only removes unchanged matrix-unit
terms with nonnegative coefficients. Thus `C-E` is positive as a
Hilbert--Schmidt *quadratic form*, although it is not asserted to be a
positive map. The exact per-seam/per-coin check covers all 261,360 units.
Arbitrary mixtures selected separately on different seams retain the fair
comparison gap `g0`.

For repair probability `r`, seam probabilities `s_e` and endpoint coin
probabilities `p_(e,c)`, comparison in the commutator identity gives

```text
Gamma = g0 (pi_min/pi_max) 60 r min_(e,c)(s_e p_(e,c)),
Re <X,(I-Phi)X>_pi >= Gamma ||X-Tr(pi X) I||_pi^2.
```

Indeed, compare each branch coefficient with its fair value `1/60`, replace
`pi` below by `pi_min I`, and bound the centered weighted variance above
by `pi_max` times the unweighted centered variance. In the bounded family,
`pi_min/pi_max>=1/16`, `r>=1/16`, `s_e>=1/480` and `p_(e,c)>=1/5`.
Hence the common bound is `Gamma>=g0/10240>9.3468e-7` per native attempt.
The canonical biased selector has the stronger bound greater than
`0.00047519`. These deliberately conservative bounds require no fitted
mass, physical clock or stationary archive.

## Consequences and limits

Let `X0=X-Tr(pi X)I`. Cauchy--Schwarz gives the robust conservation test

```text
||X0||_pi <= ||X-Phi(X)||_pi / Gamma.
```

Self-adjoint positivity also gives `||Phi(X0)||_pi <= (1-Gamma)||X0||_pi`.
Products of different exogenously chosen attempt maps obey the product of
these factors, provided they use the same endpoint weights and hence the
same reference. They need not commute or be identically distributed. For a
fixed selector, all nontrivial eigenmodes have eigenvalues in `[0,1-Gamma]`:
the source channel is primitive and its unique stationary state is `pi`.
The supported Kraus span independently contains a strictly positive diagonal
matrix; its intertwining relation rules out a nontrivial peripheral phase.
That span witness is mathematical, not a new coherent source service.

An untouched reference factor extends the weighted inequalities by tensoring
the source map with identity and using the product metric `pi tensor sigma_R`
for a faithful reference density `sigma_R`. Centering then projects onto `I tensor B(R)`,
not onto scalar matrices on the whole joint system. These norm statements do
not by themselves specify physical preparation or a diamond-distance budget.

With record feedback, selector choices can correlate with the source. The
exogenous product contraction must not be applied to the joint feedback
process. Only the conditional-map source conservation conclusion above is
claimed there. No finite spectral gap here is a particle mass, a Hamiltonian
energy gap or an exclusion of actively protected logical information.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.native_fixed_algebra.fixed
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/native_fixed_algebra
```
