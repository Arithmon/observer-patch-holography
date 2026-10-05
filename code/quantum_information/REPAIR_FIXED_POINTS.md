# Exact conservation and slow repair are different conclusions

This repair under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
starts from main `3d848b7f`. It fixes a false conserved-observable verdict
in the existing finite repair implementation. The associated convergence
theorem remains valid. The implementation had silently removed its slowest
decaying mode by classifying that mode as an exact fixed point.

## Reproduction and changed verdict

Let `X,Z` be Pauli matrices, `rho=I/2`, and take the two observable algebras
`span(I,Z)` and `span(I,Z+e X)`. Give each reference-preserving expectation
rate `1/2`. Comparing off-diagonal entries proves that their intersection
is the scalars for every nonzero `e`. The two-dimensional diagonal algebra
survives only at `e=0`.

| Input | Previous result | Correct conclusion |
| --- | --- | --- |
| `e=1e-12` | Intersection construction rejected the identity | The intersection is exactly `span(I)` |
| `e=1e-14` | Intersection dimension 2; reported gap approximately 1 | Dimension 1; true gap approximately `2.5e-29` |
| `e=1e-100` | Intersection dimension 2; reported gap 1 | Dimension 1; true gap approximately `2.5e-201` |
| Masked algebra, channel or reference entries | Hidden payload could pass validation | Missing evidence must raise |

Commit `3abccf1e` retains nine failures before the repair. Three additional
masked-reference failures were reproduced during implementation. The new
intersection calculation gives the correct dimension down to the smallest
positive binary64 tilt. The existing gap routine now refuses an unresolved
positive mode instead of excluding it and reporting the next eigenvalue.
It does not replace a supplied rounded channel by an ideal channel.

## Exact intersection, without a numerical rank cutoff

If the independent vectorized basis columns are U and V, solve

```text
[U -V] (a,b)^T = 0.
```

Sending `(a,b)` to `U a` maps this kernel onto the intersection. It is
injective: `U a=0` forces `a=0`, and then `V b=0` forces `b=0`.
Consequently a full exact kernel basis supplies an independent basis of
the entire intersection. Iterating gives the common span of every member
of a finite family. No alignment angle is declared zero by a threshold.

Every real and imaginary component of a supplied binary64 entry is rational.
The implementation retains these entries before QR and performs this
calculation over `Q(i)`. It then converts an appropriately scaled basis to
a numerical chart for the existing projections and channel diagnostics.
The exact result is retained internally across subsequent intersections.
Without that retention, rounding a new rational basis can change the next
intersection; a three-dimensional, integer-input regression demonstrates it.

This gives a complete intersection calculation for the supplied spans.
It does not assert that a noisy presentation equals an unknown ideal
algebra. In particular, `basis` remains a floating QR chart: exporting it,
rounding it and constructing a new object supplies new data. No exact
round-trip guarantee applies to those changed entries. Numerical algebra
closure, modular invariance, channel checks and spectral-gap evaluation
retain their stated tolerances. An unresolved numerical output basis or
an intersection without a common identity is rejected. Exact elimination
is intended for the existing small finite witnesses and can cost more
than a floating singular-value calculation.

The shared numeric conversion now checks operators and references before
array coercion. Masks, mixed Boolean entries and lossy conversions cannot
supply evidence for the conserved algebra. Supplied basis arrays are copied
and retained read-only so caller mutation cannot change the exact problem
after the numerical object has been validated.

## Complete relaxation law for the counterexample

The finite source maps themselves can be supplied exactly. Put

```text
N = Z+e X,
P(A) = (A+Z A Z)/2,
Q(A) = (A+N A N/(1+e^2))/2,
L = a(P-I)+b(Q-I),                  a,b > 0.
```

Since `N^2=(1+e^2)I`, the normalized N is a Hermitian unitary. Both maps
are averages of unitary conjugations, hence completely positive, unital
and trace preserving. They are idempotent, self-adjoint for the tracial
GNS inner product, and preserve `rho=I/2`. They have precisely the two
displayed algebras as their ranges. For rational e and rates, all entries
of these superoperators and their generator are rational.

On the Bloch coordinates `(X,Y,Z)`, the positive decay operator `-L` is

```text
       [ a+b/(1+e^2)       0       -b e/(1+e^2) ]
B  =   [      0          a+b             0      ].
       [ -b e/(1+e^2)      0        b e^2/(1+e^2)]
```

The full operator space also contains the identity with eigenvalue zero.
Set `g=a+b`, `d=a b e^2/(1+e^2)` and
`r=sqrt((a-b)^2+4 a b/(1+e^2))`. Direct determinant calculation gives

```text
det(t I-(-L)) = t (t-g) (t^2-g t+d),
lambda_minus = 2d/(g+r),
lambda_plus  = (g+r)/2.
```

These are the complete four eigenvalues, including multiplicities. The
rationalized expression for `lambda_minus` avoids subtracting two nearly
equal positive numbers. For every nonzero finite e, its value is strictly
positive, the only fixed observables are scalars, and every initial state
converges to `I/2`. At e=0, Z is exactly conserved and the gap on the
remaining complement is g. This discontinuity in the dimension of the
fixed algebra is why omitting the small mode changes the conclusion.

For the initial state `(I+z0 Z)/2`, its later Z expectation is exactly

```text
z(t)/z0 = w_minus exp(-lambda_minus t) + w_plus exp(-lambda_plus t),
w_minus = (lambda_plus-B_ZZ)/r,
w_plus  = (B_ZZ-lambda_minus)/r.
```

The weights are nonnegative and sum to one; the formula is also valid at
e=0, where the slow rate is zero and its weight is one. For equal half-rates,

```text
lambda_minus = e^2 / [2 sqrt(1+e^2)(sqrt(1+e^2)+1)]
             = e^2/4 + O(e^4),
relaxation time = 1/lambda_minus ~ 4/e^2.
```

At time one the fractional Z loss is
`e^2 (2-exp(-1))/4 + O(e^4)`. At one slow relaxation time the retained
fraction tends to `exp(-1)` as e tends to zero. Thus an arbitrarily stable
finite-time record can still be erased at later model times. The limits
`e -> 0` and `t -> infinity` do not commute for this observable. No
continuum fit, selected successful instance or extra dynamical premise
enters this finite conclusion.

## Connection to general primitive repair geometry

For any finite family of the reference-preserving expectations already
specified in the [repair theorem](ALGEBRAS_AND_REPAIR.md), let `A=-L` and
use its GNS inner product. Their Dirichlet identity is

```text
E(x) = <x,A x>_rho = sum_m gamma_m ||x-E_m x||_rho^2.
```

It decides conservation: `E(x)=0` exactly when every primitive fixes x.
The spectral theorem also gives the finite-time bounds

```text
||exp(-t A)x-x||_rho^2 <= t E(x),
||x||_rho^2-<x,exp(-t A)x>_rho <= t E(x),            t >= 0.
```

Indeed, integrate `(1-exp(-t u))^2 <= t u` and
`1-exp(-t u) <= t u` against the positive spectral measure of x.
Small Dirichlet energy therefore bounds temporary drift. It never supplies
an exact zero merely because the drift is below a measurement threshold.

For two repairs the entire decay spectrum follows from their principal
angles in the same GNS space. On a nontrivial angle block, the projections
are `P=diag(1,0)` and `Q=[[c^2,cs],[cs,s^2]]`, with `c^2+s^2=1`.
The block of `a(I-P)+b(I-Q)` has trace `a+b` and determinant `a b s^2`,
hence eigenvalues `[(a+b) +/- sqrt((a-b)^2+4 a b c^2)]/2`.
The four common range/kernel subspaces give rates `0,b,a,a+b`.
This is the standard two-projection decomposition; see the deterministic
spectral reduction in [Kargin, Lemma 2.1](https://arxiv.org/pdf/1205.0993).
Its use here connects the existing finite repair convergence theorem to
the complete slow-memory witness above. It is not a new physical law.

## Impact on OPH evidence

- The common-algebra dimension used by `repair_generator` is now computed
  from the retained exact spans. A small positive mode stays in the
  complement and must pass the existing spectral-resolution check.
- The finite canonical-repair convergence theorem is unchanged. Its
  numerical recognizer can no longer infer extra protected observables
  from almost coincident declared algebras.
- The exact protected-record/Green-Kubo witness has explicitly preserved
  fibres and its own rational verifier. Its protected label and published
  conclusions are unchanged. The counterexample explains why a measured
  small drift alone would not establish that witness's exact protection.
- Existing nontracial, complex recharting, clock-unit, collar, regional
  separation and MaxEnt controls pass. Signatures and ordinary trace
  conventions are retained. No current committed receipt calls this
  intersection routine; no pinned data, claim payload, paper, book or
  frozen registration is regenerated. The mandatory runner is unchanged.
- Rates and model time remain supplied. This audit does not select a
  physical clock, realize a simulator law, or satisfy the model-selection
  tasks #1025/#1026. It corrects an existing finite evidence claim without
  adding a prerequisite queue.

## Reproduction and hostile controls

```text
PYTHONPATH=code python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent -W error
```

The PR adds **45 cases**. Independent controls compare all four exact
superoperator eigenvalues, construct Choi matrices directly from matrix-unit
actions, evolve positive states with a separate matrix exponential, and
check the near-zero lifetime formula at private 80-digit precision.
Other controls cover tensor intersections in every input order, genuinely
complex common observables, changes of basis units, rational intersections
whose numerical output cannot be reused as exact input, and missing data.

Full affected suites: **827 passed on Linux; 821 passed on Windows with six
existing extended-precision skips**, with warnings treated as errors.
Ten isolated mutations were rejected by the 45 new cases:

| Deliberate defect | Failing cases |
| --- | ---: |
| Restore the old rank-threshold intersection | 12 |
| Ignore the final algebra | 24 |
| Conjugate the supplied complex span | 1 |
| Drop entries below `1e-12` | 12 |
| Conjugate the nullspace coefficients | 3 |
| Discard exact output before a later intersection | 1 |
| Strip operator masks and conversion checks | 7 |
| Strip reference masks | 2 |
| Infer the kernel by deleting small eigenvalues | 3 |
| Disable the unresolved positive-gap check | 3 |

Each mutant ran against its own copied source outside the repository import
path. The retained failures demonstrate false conservation, wrong rates,
incorrect complex spans and accepted missing evidence; receipt checksums
alone are not the verification mechanism.
