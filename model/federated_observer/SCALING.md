# What propagates on a federation of identical elementary carriers

The pinned federation supports an exact density field on every binary
fixed-number sector. Its native transfer is uniformly positive, so its
logarithm exists at every size. The same local rule also gives a controlled
propagation estimate independent of the number of observers. These results
use the declared geodesic-screen wiring in [source.py](source.py); they do
not select that wiring or prove a Lorentzian continuum limit.

There are `C=20*4^level` elementary twelve-port carriers, `V=12C` ports and
`M=63C/2` seams. Each member has its thirty internal icosahedral seams and
three connections to neighbouring members. Each port has degree at most
six. A seam attempt applies the same native endpoint operation with a fair
coin. On binary loads this is identity or endpoint swap with equal
probabilities. Uniform selection among all `M` seams is the declared
scheduler. A federation's role as a larger observer changes its accessible
records; it does not require a different elementary update law.

## Exact density closure in every number sector

For a configuration `x`, let `n_i(x)=x_i` and
`F_f(x)=sum_i f_i n_i(x)`. On the configuration transfer Hilbert space set

```
G = sum_edges (I - Swap_edge),    K = I - G/(2M).
```

If `L` is the port-graph Laplacian, then

```
G F_f = F_(L f),
K F_f = F_((I-L/(2M))f).
```

Proof: for an edge `(a,b)`, swapping its endpoints changes `F_f` by
`(f_a-f_b)(x_a-x_b)`. Sum over edges and collect each `x_i`. This is an
identity for every configuration, not a closure based on independent
particles or an assumed Gaussian state. Fixed total number is preserved,
so the same identity holds on every sector, including interacting exclusion
sectors. The executable checks all four endpoint inputs and both coefficient
basis vectors against the actual native function. Tests additionally
enumerate every number sector of a four-cycle independently.

This generalizes the committed exact identity
`Q*Z=-(1/2) Z*L` in
[`NativeRepairFourCycle.lean`](../../Lean/Thermodynamics/NativeRepairFourCycle.lean).
That theorem uses rate-one attempts on each edge, whereas `K` here means
one globally selected seam attempt. The factor `1/(2M)` is therefore
essential. No new Lean theorem is claimed by this executable attachment.

## A uniform positive transfer and controlled logarithm

An edge-average operator is `P_e=(I+Swap_e)/2`, an orthogonal projection.
For the three edges of a triangle on binary ports, put
`B=Swap_01+Swap_02+Swap_12`. Direct exact integer arithmetic on the full
eight-dimensional three-port space gives

```
B* B = 3B,
sum_triangle P_e - (3/2)I = B* B/6 >= 0.
```

Every elementary icosahedron has twenty triangles, and every internal seam
belongs to two. Summing these inequalities therefore gives
`sum_internal P_e >= 15I` per carrier. Operators on other ports are identity
extensions, and every inter-carrier `P_e` is positive. Consequently,

```
(15C/M) I <= K <= I,       15C/M = 10/21.
```

This bound holds on the full binary configuration space and on every
invariant fixed-number sector. It does not enumerate a `2^V` matrix or
assume that the seam projections commute. In particular `H=-log K` is
well defined with a size-independent positive lower bound on `K`.
Writing `D=I-K`, exact spectral calculus gives

```
0 <= D <= (11/21) I,
H_m = sum_(k=1..m) D^k/k,
0 <= H-H_m <= [(11/21)^(m+1)/((m+1)(10/21))] I.
```

The last step is the positive logarithm series with its denominators bounded
by `m+1` in the tail. It attaches the spectral-calculus method used in
[`TransferMatrixLift.lean`](../../Lean/SourceLaw/TransferMatrixLift.lean)
to this actual finite exchange transfer. It does not import that module's
separately supplied Gaussian field action.

For an eigenvector `L f=lambda f`, the source density mode has the exact
transfer eigenvalue and reconstructed energy

```
r = 1-lambda/(2M),     E = -log(1-lambda/(2M)).
```

For a nonconstant `f` and a sector `0<q<V`, this is a nontrivial mode. The
faithful mathematical reference `sigma_K=K/Tr K` supplies the same modular
energy differences: `-log sigma_K=H+log(Tr K)I`. The scalar term cancels in
modular automorphisms. This is a reference on the configuration transfer
Hilbert space, not an asserted native density preparation. The native
recorded CP instrument, its classical history transfer and the reconstructed
unitary group are different objects.

## Density poles and a controlled dispersion counterexample

The receipt diagonalizes only the sparse `V`-dimensional port Laplacian at
levels 0, 1 and 2. Exact closure transports these density modes into every
nontrivial number sector, avoiding exponential configuration enumeration.
The eight lowest positive Laplacian values, their transfer energies and
eigenpair residuals are numerical diagnostics. A sweep means `M` native
attempts, so `M E` is also reported as a rescaled energy. The gap ratios
between successive twofold screen resolutions are useful finite checks;
three levels establish no scaling limit or spatial dimension.

A periodic cycle is included as an explicit comparison graph, **not** as
the selected federation. Its density modes obey

```
x_k = (1-cos(2 pi k/N))/N,
E_k = -log(1-x_k),
E_2/E_1 -> 4.
```

For every declared `N>=8`, exact rational bounds are

```
4(1-10/N^2)(1-20/N^3) <= E_2/E_1 <= 4/(1-80/N^3).
```

They follow from `x <= -log(1-x) <= x/(1-x)`,
`x_2/x_1=4 cos^2(pi/N)`, and `pi^2<10`. Both bounds tend to four and the
lower bound exceeds two. Thus this gapless density sector has quadratic
long-wavelength dispersion (`z=2`), rather than the ratio two of a
gapless linearly dispersing mode. Multiplying the time coordinate by one
positive common factor changes both energies equally and cannot change
the ratio. Passing from `K` to the modular reference `K/Tr K` also leaves
it unchanged.

This refutes an automatic passage from exchange dynamics plus transfer or
modular reconstruction to linear relativistic dispersion. It is not a
no-go for other OPH sectors, coherent protocols or federation laws.
Diffusive conserved modes can also appear in thermal states of relativistic
theories. A massive relativistic dispersion is quadratic above its rest
energy at small momentum, so the control is explicitly about these
**gapless** modes, not all quadratic spectral features.

The square-root frequency in
[`ContinuumConnection.lean`](../../Lean/SourceLaw/ContinuumConnection.lean)
belongs to a specified Gaussian slab action. Substituting `sqrt(L)` for
the actual `-log(I-L/(2M))` would change the model's source correlations;
it is not an available consequence of the transfer reconstruction.

## A uniform operational propagation estimate

For a fixed word of `n` endpoint operations, a port value can depend only
on initial values in its closed graph-distance-`n` neighbourhood. This is
the direct endpoint instance of
[`DependencyCone.lean`](../../Lean/ObserverPatchHolography/Locality/DependencyCone.lean).
It is a bound in committed operation count, with no duration attached.
The uniform state-independent scheduler used here introduces no additional
state consultation. An adaptive scheduler would require the consultation
region in
[`AdaptiveScheduler.lean`](../../Lean/ObserverPatchHolography/Locality/AdaptiveScheduler.lean).

A second, explicitly declared construction places independent rate-one
Poisson attempt clocks on each seam. It is equivalent to Poissonizing the
global uniform attempt process at rate `M`. Thinning by the fair swap coin
gives rate-`1/2` swaps on each seam. Its density expectation obeys

```
u(t) = exp(-t L/2) u(0).
```

Here `t` is the **per-edge attempt Poisson parameter**, not a recovered
physical time, native first-event counter, or modular parameter. The
coefficient agrees with the existing four-cycle drift theorem. The
semigroup is the density law of this stochastic construction; it is not
the reconstructed unitary group `exp(-itH)`.

Trace one occupation label backwards through the swaps. Its jump rate is
at most `Delta/2`, where `Delta<=6` is the maximum port degree. To reach an
initial perturbation at distance at least `R`, the label must jump at least
`R` times. Its jump count is dominated by `Poisson(Delta*t/2)`. The factorial
moment inequality gives the exact bound

```
P(influence across distance R by t)
    <= min(1, (Delta*t/2)^R/R!).
```

For `r` retained port-value reads at specified times at most `t`, a union
bound multiplies the allowance by `r`. Net changes formed from those
readings inherit the same comparison. This coupling sums out the visible
native coin: its conversion to a swap bit depends on the endpoint
orientation. It therefore does not certify the joint visible native
coin/current history. Nor is it a bound for arbitrary coherent native CP
observables.
The executable evaluates the sparse heat propagation at `t=1/4` and compares
the probability beyond radii 2, 4 and 8 against the rational allowances.
It retains the nonzero tails; Poissonized evolution has no strict finite
light cone.

The bound is independent of the number of carriers. If a family of
federations has matching initial data and increasingly large identical
rooted neighbourhoods, couple the shared seam clocks and hidden swap bits.
The same tail bound makes each finite occupation-read law Cauchy. This is a concrete
conditional local-limit criterion. The pinned level refinements are not
asserted to satisfy it: local port-frame routing and state/coarse-graining
compatibility need their own check. A spherical drawing alone supplies
neither that convergence nor a bulk spacetime dimension.

For a continuum density limit on an independently specified manifold,
[van Ginkel and Redig's hydrodynamic theorem](https://link.springer.com/article/10.1007/s10955-019-02420-2)
reduces the exclusion process to convergence of its underlying graph
walks, given the required approximating-grid and initial-profile hypotheses.
Its limit is a heat equation. It is a possible reuse target once those
hypotheses are proved for the selected federation; the theorem does not
select the screen, its wiring or a relativistic clock.

The useful next attachment is therefore precise: establish compatible
local access and refinement for this growing family, then determine which
source-derived sector and observer clock support the target joint
propagation statistics. Positivity, a noncommutative finite algebra and a
modular parameter do not on their own discharge that propagation test.
