# Preserving the crossing observable as well as spacetime limits

The S and A counterexamples establish that the causal/count limit alone does
not fix crossing entropy. There is also an explicit sparse family that
preserves the existing leading crossing count with exactly one nat per
unordered crossing pair. It requires neither a direction-dependent entropy
weight nor a fitted scalar recalibration.

## Construction and path proof

For every integer t>=1 define

    q=2^(8t), R0=2^(4t), m=2^(4t), K=2^(3t), R=m*K=2^(7t), h=L/q.

The comparison B0 is the original dense integer ball of radius R0, with tick
h*R0/c and a0=L/sqrt(q). Define T to contain every m*u with u in Z^3 and
|u|<=K, together with +/-2^b e_i for 0<=b<4t. Include zero once. T has
exact sign/permutation symmetry and radius R. Use the same q^3 sites,
previous-layer reads clipped to the box, tick delta=hR/c and event weight
delta*h^3. The T radius is a=L*q^(-1/8), not the B0 radius. These are
fully specified families; no successful propagation or entropy condition is
used to select individual edges.

To route an integer displacement d, round d/m coordinatewise to z in Z^3,
with nearest-integer ties away from zero, and write d=m*z+r. Each residual
coordinate has magnitude at most m/2. If z is nonzero, put
k=ceil(|z|/(K-2)); otherwise use no coarse steps. Round j*z/k coordinatewise
at j=0,...,k. Successive differences have norm at most
|z|/k+sqrt(3)<=K-2+sqrt(3)<K, so their m-scaled versions are legal coarse
reads. Correct r by its signed binary digits, using at most 12t extra steps.
Every step is charged one full T tick, and

    |d|/R <= N <= (|d|/m+2)/(K-2)+1+12t.               (B1)

The first inequality is the range bound. For the second use
|z|<=|d|/m+sqrt(3)/2 and ceil(x)<=x+1. Each rounded coarse position lies
within sqrt(3)*m of the straight segment from 0 to d: one half-cell error
comes from rounding the path and one from rounding its endpoint. During the
binary corrections the distance to d is at most |r|_1<=3m/2. The whole
route therefore lies in the Euclidean tube of radius 2m about that segment.

In physical units the leading term in (B1) is
K/(K-2)*distance/c; its additive term is at most
(2/(K-2)+1+12t)*a/c. Both the time excess and the tube width 2mh tend to
zero. Zero edges provide integer waiting. The same strict-margin and
equal-cell argument as in DERIVATION.md proves the round causal limit,
actual generated-diamond and strict-pair limits, 1/10 ordering fraction,
fourth-root count-clock ratios and the weighted mass statements for T.
Finite cell weights and ticks are not equal between T and B0; it is their
normalized volume/count and clock-ratio limits that coincide.

The degree is |Z^3 intersect B_K|+24t=Theta(K^3)=Theta(q^(9/8)), versus
Theta(q^(3/2)) for B0. Over horizon L/c there are q/R=q^(1/8) transitions,
so T uses Theta(q^(17/4)) logical read incidences versus Theta(q^5) for B0.
These comparisons fix q; the T radius and convergence errors shrink more
slowly. They are not equal-accuracy runtime or native resource savings.

## Exact coordinate-cut count and the unchanged entropy coefficient

Let

    C_K = sum_(u in Z^3, |u|<=K, u_x>0) u_x.

Every dyadic x-axis step contributes its length times q^2 crossing pairs,
and the other added axis steps do not cross this cut. Therefore

    X_T = sum_(|u|<=K, u_x>0)
          m*u_x*(q-m*|u_y|)*(q-m*|u_z|) + q^2*(m-1).  (B2)

In particular,

    m*(q-mK)^2*C_K <= X_T <= m*q^2*C_K + q^2*(m-1).

To evaluate the leading coefficient without a fit, unite the unit cubes
centered at the integer points in B_K. That union contains B_(K-1) and is
contained in B_(K+1). The integral of x_+ over a cube is its center's x_+
except at x=0, where the integral is 1/8. There are at most (2K+1)^2 such
central cubes. The spherical integral is pi*r^4/4. Hence

    pi*(K-1)^4/4-(2K+1)^2/8 <= C_K <= pi*(K+1)^4/4.

Since m*K^4=q^2, (B2) gives the explicit finite enclosure

    (1-mK/q)^2 * [pi*(1-1/K)^4/4-(2K+1)^2/(8K^4)]
        <= X_T/q^4
        <= pi*(1+1/K)^4/4+(m-1)/q^2.                (B3)

Both ends tend to pi/4. B0 has X_B0/q^4 -> pi/4, so X_T/X_B0 ->1.
This preserves the original golden-family leading value as well: its
complete-neighbor theorem has the same pi*q^4/4 asymptotic. The exact finite
golden values, Fibonacci sequence and finite box errors are not identified.

## Every planar orientation, with the same scalar normalization

For any fixed positive-area regular plane section Sigma of the cube, let
X_T(Sigma) count unordered crossing pairs. Here a regular section is a cut
through the interior, with positive-volume portions of the cube on both
sides and a piecewise smooth section boundary. A supporting boundary face
is excluded: it has positive area but no crossing pairs. The coarse displacement mesh is
epsilon=mh; its ratio to the read radius is epsilon/a=1/K. Each allowed
coarse displacement samples the same radius ball with density epsilon^-3.
The receiver lattice has density h^-3. Replacing each receiver by its h-cell
changes a fixed-displacement crossing slab count by O(h^-2 L^2): its bounded
faces have an h-wide collar. Summing over O(K^3) displacements gives a
relative error O(1/(mK)) against the leading m*K^4*q^2 count.

For the displacement sum, cell boundary collars and the Lipschitz plane-flux
function give relative error O(epsilon/a)=O(1/K). Except within an a-wide
collar of the section boundary, clipping to the cube does not change the
crossing slab; its relative contribution is O(a/L)=O(mK/q). These are the
same cap-count estimates as the existing dense theorem, now with distinct
receiver and displacement meshes. The added dyadic edges contribute at most
O(m*q^2), which is o(q^4), for each fixed orientation. Thus

    X_T(Sigma) = (pi/4) * [a^4/(m^3*h^6)] * area(Sigma) * (1+o(1))
               = (pi*q^4/(4L^2)) * area(Sigma) * (1+o(1)).       (B4)

The coefficient is rotation independent because
integral_(|v|<=a) (v dot normal)_+ dv=pi*a^4/4 for every unit normal.
The relative error is O(1/K+mK/q+(m-1)/q^2), with constants depending on the
fixed section. Nearby offset sections used in the integer fixtures have
uniform bounds. This proves the same leading one-nat crossing entropy and
Planck-area dictionary as B0, for every such orientation, without assigning
different entropy to different reads. Crossing entropy per layer is the
observable being preserved; crossing traffic per unit model time is not.

## Scalar action and an explicit free quantum detector comparison

T has the same exact cubic symmetries, positive second moment and range bound
used in DERIVATION.md. For the final replacement action take one half of
that exact-moment operator and one half of the ordinary nearest-neighbor
Laplacian, whose edges are already present. This preserves every positive
coupling and closes the uniform-weight action's high-frequency alias problem,
as proved in SPECTRAL_COMPLETION.md. Its spatial error is at most a^2*M4/4
and its finite-time
Klein-Gordon/detector bound is c^2*tau^2*a^2*L^(3/2)*M4/8. Here a=L*q^(-1/8),
so the error tends to zero. The positive periodic action also satisfies the
same centered-leapfrog obstruction at its own tick a/c.

A bounded free quantum comparison follows without identifying the
layered read program with that evolution. Use a periodic cube and positive
mass mu, and canonically quantize each defined positive action. For a fixed
finite real Fourier basis of smooth modes, sufficiently fine lattices sample
these modes with exactly their continuum inner products. Translation
invariance makes each sampled mode an exact eigenvector. For wavevector p,

    lambdahat_V(p) = (alpha * sum_v [1-cos(h*p dot v)] + lambda_nn(p))/2,
    |lambdahat_V(p)-|p|^2| <= a^2*|p|^4/4.

Put omega_V(p)=sqrt(c^2*lambdahat_V(p)+mu^2). Then
|omega_V(p)-omega(p)| <= c^2*a^2*|p|^4/(8*mu).
Each vacuum mode has position variance hbar/(2*omega_V), momentum variance
hbar*omega_V/2 and zero symmetrized cross covariance. With explicitly
sampled finite-Fourier coherent initial means, exact harmonic evolution
converges mode by mode, uniformly on every bounded time interval. Therefore
the joint Weyl characteristic functions of every finite collection of these
position/momentum smearings converge to the same continuum Gaussian values
for S, A, B0 and T. This is an explicit free quantum observable result, not
an inference from the classical expectation alone. The reciprocal-frequency
bound uses mu>0; no massless vacuum, interacting theory, local physical
instrument, source-selected quantization or full vacuum-state convergence
is asserted by this finite-mode argument alone. SPECTRAL_COMPLETION.md gives
the stronger completed result: uniform whole-spectrum control, the full
thermal partition/energy/entropy limit, smooth compact Weyl observables and
the causal commutator on the same clipped box. The clipped-box result also
allows zero mass because its boundary supplies a proved frequency gap.

## Why a radius change is needed for this sparsification

There is a sharp obstruction within translation-invariant sub-stencils of
the original radius ball. Let V_n be any symmetric subset of B_(R0), at the
same q sites and with one nat per unordered crossing pair. Suppose its
crossings agree with B0 to leading order at each of the three coordinate
midplanes. Let E_n=B_(R0) minus V_n be the removed vectors. Exact finite
crossing formulas give

    sum_(i=1..3)(X_B0,i-X_V,i)
        >= (q-R0)^2/2 * sum_(v in E_n) |v|_1.

The left side is o(R0^4*q^2), so sum_E |v|_1=o(R0^4). For every fixed
eta>0, separate removed vectors inside the cube |v|_infinity<=eta*R0
from the rest. Then

    |E_n|/R0^3 <= (2*eta+1/R0)^3
                   + sum_E |v|_1/(eta*R0^4).

Taking limsup and then eta down to zero gives |E_n|=o(R0^3). Since the
ball has Theta(R0^3) lattice points, |V_n|/|B_(R0)| tends to one.
Thus no asymptotically substantial stencil thinning at the original radius
can retain these three leading crossing counts under this fixed dictionary.
This does not constrain arbitrary position-dependent read graphs. T evades
the proved obstruction by its explicitly larger, still vanishing radius.

## Decision

The completed comparison is stronger than a conditional proposal to replace
M1. T is an actual defined family for which the named causal/count, scalar,
free quantum and thermal, and leading planar crossing observables are proved to
coincide. S and A show why kinematic agreement alone is insufficient; the
fixed-radius bound proves a real obstruction and T supplies its constructive
escape. Exact finite tuples, native operations and additional interacting or
gravitational attachments retain their original contracts. No physical
premise is discharged by calling a regulator an OPH source.
