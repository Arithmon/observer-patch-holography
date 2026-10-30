# Conservation and contraction for either exchange character

Both complete [native exchange lifts](../exchange_sector/STATISTICS.md),
`chi=+1` and `chi=-1`, have only scalar universally conserved source
observables on H66 under positive all-seam exploration. They also have the
same rigorous lower contraction bound in the declared weighted norm.
Changing the exchange character therefore supplies no source-only
protection mechanism in this finite carrier.

One step means a **complete selected-character channel application**.
For the negative character this includes the declared record-selected
phase service after the native attempt. Equal pulse cost, physical duration,
or autonomous selection of the character is not asserted. The character
is a fixed classical choice for each application; a coherent character
register or an additional internal multiplicity changes the carrier.

## Faithful reference and exact source algebra

Use the same admitted positive endpoint weights `theta` and reference
`pi_ij=theta_i theta_j/Z` as the [native contraction theorem](README.md).
The complete signed branch has

```
K_h = sqrt(w_h) sum_(i in D_h) s_h(i) |f_h(i)><i|,
s_h(i) in {+1,-1},  f_h injective.
```

Its effect is still `w_h P_Dh`, and its action on the diagonal reference
is unchanged. Every seam and both C/E laws preserve the same faithful
`pi`, including arbitrary per-seam classical law mixtures, identity holds,
zero slots and reserved abort addresses. No prepared reference is supplied
by this invariance argument.

The effects separate every occupation address. Products of those effects
and their complements therefore construct each rank-one diagonal `P_i`.
For each supported transition,

```
P_j K_h P_i / sqrt(w_h) = s_h(i) |j><i|.
```

The nonzero sign can be removed by scalar multiplication *inside the
generated algebra*. A checked spanning tree of transitions generates all
matrix units by multiplication and adjoints. This proves the full source
matrix algebra for either character and every supported per-seam C/E
choice. It does not equate the complete positive and negative instruments
by scalar phases on their whole raw branches.

Completeness and faithful invariance give the same exact Dirichlet identity
for arbitrary complex `X`:

```
Re Tr[pi X* (X-Phi(X))] = (1/2) sum_h ||X K_h-K_h X||_pi^2,
Phi(X) = sum_h K_h* X K_h,   ||X||_pi^2 = Tr(pi X*X).
```

If `Phi(X)=X`, every commutator vanishes. Apply the argument to `X*` as
well to obtain commutation with the adjoints. The full algebra then makes
`X` scalar. With classical record feedback the same conclusion holds when
the source-only observable must be conserved on every admitted
source/controller input: it must be fixed by each conditional map.

## An exact transfer of spectra and bounds

There is a further simplification on this carrier. For occupation pairs
`i,j`, let `C=i intersect j` and `D=i symmetric_difference j`, and define

```
S(|i><j|) = (-1)^(sum_(q in C,d in D) 1[d<q]) |i><j|.
```

On a surviving moving matrix unit, both configurations move the same common
particle; their two exclusive spectators remain fixed. The product of their
fermionic signs counts exactly the exclusive spectators crossed by that
common particle. Hence it equals the initial S sign times the final S sign.
On holds and diagonal units the identity is immediate; zero branches remain
zero. This proves, separately for every complete labeled branch,
`Phi_minus=S Phi_plus S`, checked on all surviving matrix units.

S is an involution, fixes every diagonal observable, preserves adjoints and
is an isometry for every diagonal reference metric. Thus every fixed
selector's two channels have exactly the same spectrum, and their centered
GNS bounds transfer without an estimate. The same isometry acts on arbitrary
untouched reference blocks. It is a mathematical change on operator space,
not a physical identification of the instruments: S is not even positive.
Applied to the equal coherent superposition of `{0,1}`, `{0,2}`, `{1,2}`,
its output has expectation `-1/3` in the normalized vector `(1,-1,1)`.
This control prevents using S as a preparation, channel or physical gauge.

## The quantitative bound retains the signs

An independent entrywise check also verifies the transferred bounds. Each
actual signed seam channel splits the matrix-unit space into blocks
of size at most two. One-unit coefficients lie in `[0,1]`. The nontrivial
two-unit block has the endpoint heat-bath form with a common off-diagonal
sign `epsilon=+1` or `-1`. Weighted detailed balance still holds. Conjugating
this two-unit block by `diag(1,epsilon)` gives the positive heat-bath block
and preserves its diagonal reference metric. Both the weighted channel
matrix and its identity complement have nonnegative principal minors.
The implementation checks the actual signed entries and these minors
exactly on all `60*66^2` units. These individual observable-space block
changes are not one occupation gauge for every source seam.

For the fair uniform C channel, the population block is exactly unchanged.
For each off-diagonal matrix-unit entry, the triangle inequality bounds the
absolute signed sum by the corresponding positive-lift entry. The latter
has row sum at most `5/6`. The signed matrix is itself symmetric, so its
operator norm on that subspace is at most `5/6`. This is a verified norm
bound, not replacement of a signed quantum channel by its entrywise absolute
values. The existing exact population minorization therefore supplies

```
g0 = min(beta/20, 1/6) = beta/20 > 0.00957122,
beta = 142390556543804527497099183523 / 743847338880000000000000000000.
```

E refines C with the same signs on their shared columns. The only extra C
matrix-unit terms are unchanged positive diagonal terms in this observable
basis. Thus `C-E` is a nonnegative quadratic form for either character,
checked separately for every seam and coin. No positive-map claim about
`C-E` is made. The same coefficient comparison in the Dirichlet identity
now gives

```
Gamma = g0 (pi_min/pi_max) 60 r min_(e,c)(s_e p_(e,c)),
Re <X,(I-Phi)X>_pi >= Gamma ||X-Tr(pi X)I||_pi^2.
```

Here `r` is repair probability, `s_e` seam probability and `p_(e,c)` the
endpoint coin probability. In the admitted bounded family,
`Gamma >= g0/10240 > 9.3468e-7` per complete lift application. The canonical
biased parameters give `Gamma > 0.00047519`. Hence

```
||X-Tr(pi X)I||_pi <= ||X-Phi(X)||_pi / Gamma,
||Phi(X-Tr(pi X)I)||_pi <= (1-Gamma)||X-Tr(pi X)I||_pi.
```

Exogenously varying selectors and even varying classical character choices
obey the product of these contraction factors when the endpoint weights,
and therefore `pi`, remain fixed. Arbitrary untouched references extend the
proof with the product metric `pi tensor sigma_R`. The fixed algebra then
is `I_source tensor B(R)`, and centering projects onto it. A reference is
not falsely classified as lost information.

The norm bound does not establish a physical clock, Hamiltonian gap, mass,
or diamond-distance preparation budget. It does not assert contraction of
source-correlated adaptive feedback or of a growing archive. Active recovery,
additional coherent controls beyond the specified character lift, relational
charges, restricted input sets, other carriers and refining limits remain
outside the conservation exclusion. Removing seams incident to one port
restores its conserved occupation and provides an explicit control.

```
OPENBLAS_NUM_THREADS=1 python -m model.native_fixed_algebra.signed
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/native_fixed_algebra/test_signed.py
```
