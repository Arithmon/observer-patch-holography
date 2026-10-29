# Protection and transport of the same finite bound pair

This component advances [CP-03 / #1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139)
and [CP-04 / #1140](https://github.com/FloatingPragma/observer-patch-holography/issues/1140)
within the declared [Weyl-contact action](README.md). It uses the same
`c=1/8`, `U=1` as the energy, current and clock calculations. The advance
is an explicit localized composite wavepacket that moves while remaining
inside one isolated band of that action. The spatial lattice, CAR fields,
action and packet preparation remain supplied.

The nearest existing result is the complete two-fermion action and local
source-word implementation in [`interacting_field`](../interacting_field/DERIVATION.md).
Its open cubes are replaced here by a periodic cube for exact momentum
resolution. This is a declared boundary change, not an assertion of equal
finite observer histories. Native active core transport, including its
forbidden destinations and reset services, is a different construction.
The following closed-Hamiltonian result does not supply the all-branch
native protection required to close CP-03.

## The spectral gap does not depend on a fitted mass

On the complete antisymmetric two-fermion space, write `H=T+UP`, where
`T=dGamma(h)` and `P` projects onto the on-site spin-singlet pairs. Fourier
diagonalization of the one-body symbol gives

```text
||T|| <= b = 2 sqrt(3) c.
```

A single hopping term moves one constituent off the shared site; therefore
`PTP=0`. Let `V` be the number of sites, so `rank P=V`. The min–max principle
gives at least V eigenvalues at least U, since the Rayleigh quotient on P
is identically U. On `P-perp`, the interaction vanishes and the Rayleigh
quotient is at most b, so there are at most V eigenvalues above b.
Consequently, if `U>b`, exactly V eigenvalues form a separated upper band,
and its separation from the entire remaining spectrum is at least `U-b`.
There is no extra factor of two in this gap: the identity `PTP=0` supplies
the stronger lower band edge U.

For the canonical action the gap is at least `1-sqrt(3)/4`, which is
strictly greater than the exact rational lower bound `9/16`. That lower
bound uses `sqrt(3)/4 < 7/16`. The band is invariant under the closed
Hamiltonian, including arbitrary superpositions of its momenta. These
statements do not constrain a driven detector, an arbitrary error operator,
unconditioned native repair, or sectors with additional particles. No
finite lifetime under those enlarged dynamics is inferred.

## Relative localization has a uniform tail bound

In each momentum fiber, let d be the normalized contact vector. Its isolated
eigenvector is

```text
psi = U sqrt(Z) (E-T)^-1 d,
E >= U > b,      0 < Z <= 1.
```

Let `Q_R` project onto relative lattice Manhattan distance **at least R**,
using shortest periodic distance on a finite torus. For every integer
`R>=1`, a sequence of fewer than R one-particle hops cannot take a contact
pair into that region. Hence `Q_R T^n d=0` for `n<R`. The convergent Neumann
series gives

```text
||Q_R psi|| <= [U sqrt(Z)/(E-b)] (b/E)^R
            <= [U/(U-b)] (b/U)^R.
```

Squaring and capping at one bounds the complete tail probability. For the
canonical action, [`tail_bound`](protection.py) uses the outward rational
bound

```text
Pr(relative distance >= R) <= min(1, (16/9)^2 (7/16)^(2R)).
```

The index is inclusive: the R=1 tail contains every separated pair. There
is no unmentioned discarded boundary shell. The relative-distance projector
commutes with simultaneous translations, so momentum fibers are orthogonal
also after applying it. The same bound therefore holds for every normalized
wavepacket in this band, without a factor for the number of momenta. The
Neumann argument also applies to the infinite lattice at fixed spacing;
it is not a spatial continuum limit.

## Complete spatial transport replay

[`PeriodicPair`](protection.py) independently assembles every oriented
periodic bond and its Hermitian conjugate, then lifts the one-body operator
through the CAR exterior square. The retained side-four cube has 64 sites,
128 fermion modes and the full 8,128-dimensional N=2 space. All 64 momentum
fibers are solved. Their four spin components at every relative momentum
are Fourier-transformed into that full position space. No on-site-pair-only
effective hopping Hamiltonian or relative-coordinate truncation is used.
The Fourier map rejects vectors that violate fermion antisymmetry.

For all 64 eigenvectors, the independently assembled position Hamiltonian
acts with their computed energies. Their Gram matrix is the identity. This
checks the spin ordering, momentum convention, periodic seam signs, CAR
lift and contact normalization together. The finite band runs from 1 to
approximately 1.08859449; every contact residue is positive.

The initial packet coefficients are declared before propagation:

```text
a_K proportional to exp{(1/2)[cos(K_x-pi/2)
                           +cos(K_y-pi)+cos(K_z-pi)]}.
```

They are normalized over all 64 reciprocal-grid momenta. The positive
contact gauge fixes the relative phases of the eigenvectors. This is a
supplied state preparation, not a source-compiled native field encoder.
Evolution multiplies every coefficient by `exp(-i E(K) t)`, where t is
Hamiltonian **model time**. A separate sparse exponential of the complete
position Hamiltonian agrees with the fiber evolution.

The recorded observables include the particle number at every site, every
site contact probability `Pr(both constituents at x)`, and the single
complementary separated-pair outcome. The complete distribution keeps unit
probability and total particle number two. At model duration 16 its total
variation from the initial contact/separated distribution is approximately
0.0898600. The signed difference between the +x and -x neighboring contact
probabilities changes from zero to approximately -0.0709936. Thus this is
directed spatial motion, not merely changing internal phase or a periodic
center-of-mass labeling convention.

The packet remains exactly in the same band in the ideal mathematics. Its
relative separation distribution is constant under this one-band evolution;
only the center-of-mass distribution changes. The complete replay retains
the non-contact weight, approximately 0.0335570, throughout. Replacing the
band energies by a constant removes the motion. Rephasing each eigenvector
and inversely rephasing its coefficient leaves the entire state unchanged.
Removing U produces a valid free action but loses this isolated repulsive
band premise. These controls are tested independently of the reported
successful trajectory.

Numerical matrix residuals are replay diagnostics, not outward-certified
floating-point or hardware error bounds. The gap and tail bounds above are
analytic statements; their rational relaxations do not rely on numerical
eigenvalue fitting.

## The band is topologically trivial

The positive contact overlap never vanishes. Choosing
`psi_K=U sqrt(Z(K))(E(K)-T_K)^-1 d_K` therefore gives one globally defined,
periodic, nonzero section of the isolated complex line bundle. Smoothness
follows from the uniform gap and the implicit resolvent equation. The
infinite-lattice statement follows from this actual operator family, not
from sampling a finite grid of eigenvectors.

More directly, contract the hopping by `H_lambda=lambda T+UP`, with
`0<=lambda<=1`. The separating gap is bounded below by `U-b` throughout,
and at lambda zero the band consists of strictly on-site pairs. This
gap-preserving contraction makes the band trivial. The tests verify the
same positive contact section at several contraction parameters; the
all-parameter conclusion is the stated uniform-bound argument.

This pair consequently has **no assigned nontrivial real-space defect
charge**. A Weyl node's momentum-space chirality is a different invariant
and is not renamed as the pair's charge. The proof establishes a useful
protected moving composite in the supplied action while preserving the
missing native invariant, reference-like exterior and full-source
implementation obligations of CP-03.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/bound_pair/test_protection.py
```
