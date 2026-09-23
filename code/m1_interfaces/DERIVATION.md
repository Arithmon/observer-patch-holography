# From fixed sections to variational area

## 1. The finite object and its continuum dictionary

First set L=1; physical lengths multiply by L and areas by L^2. Work either
on the flat unit torus, or on the unit cube with every occupied cell confined
to the fixed interior cube C=[eta,1-eta]^3, where 0<eta<1/2. Reads in the
cube are clipped to the actual site population. The interior restriction
makes their cut count equal to the full zero-extended count once the read
radius is smaller than eta. No outside records or edges are then added.

There are q^3 sites with spacing h=1/q. A site partition is a binary vector
u; identify it with its equal-cell, piecewise-constant indicator u_h. A
symmetric integer stencil V contains zero. Each site reads the available
site at every offset v in the preceding layer. For the crossing observable,
count each unordered pair once, not twice for its two layer directions:

    X_V(u) = (1/2) sum_x sum_(v in V) |u(x+v)-u(x)|,
    E_V(u) = 4 X_V(u)/(pi q^4).                         (1)

All the stencils below have radius strictly less than q/2, so periodic
offsets identify no distinct read vectors. Zero contributes no crossing.
The physical normalized entropy is L^2 E_V. Every crossing has the same
one-nat weight; E is a fixed conversion of that statistic to area units.

Put, for integer n>=2,

    q=2^(8n), m=2^(4n), K=2^(3n), R=mK,
    epsilon=m/q, a=R/q, and m K^4=q^2.                 (2)

The spaced-ball stencil T contains every m v with |v|<=K, and every signed
axis step 2^b e_i for 0<=b<4n. Include zero once. Its degree is Theta(K^3),
its physical radius a=q^(-1/8), and its layer tick is a/c.

The exact translation form of its coarse contribution is

    E_coarse(u_h) = (2/(pi q)) sum_(|v|<=K)
                         ||u_h(.+epsilon v)-u_h||_1.  (3)

Here a^4=q epsilon^3. Consequently the formal continuum integral has the
normalization 2/(pi a^4), since

    integral_(|z|<=a) |z dot nu| dz = pi a^4/2

for every unit normal nu. This fixes the perimeter coefficient to one.
The dyadic contribution is lower order for every fixed finite-perimeter
set. Neither the choice of a surface nor its observed entropy defines T.

## 2. The missing compactness and a false minimum

On the torus take u(x)=1 exactly when x_1 mod m < m/2. Its volume is exactly
one half. Every coarse step preserves this label. A positive dyadic x-step
s<=m/2 changes a fraction 2s/m of the labels; y and z steps change none.
The two signs and the factor one half in (1) give the exact count

    X_T(u) = 2 q^3 (m-1)/m,
    E_T(u) = 8(m-1)/(pi m q) -> 0.                   (4)

The indicators converge weakly to 1/2 and have no strongly L1 convergent
subsequence. A strong limit of binary indicators is binary almost everywhere,
whereas their unique weak limit is not. Thus bounded crossing entropy alone
does not produce macroscopic regions. In particular the minimum at half
volume tends to zero, although the torus perimeter minimum is positive.

This is not only a periodic-boundary example. Let F be a ball strictly
inside C, of volume 2v, and retain the same active residue half inside F.
Coarse displacements preserve the residue. Their crossing differences are
the crossing differences of F multiplied by the active-residue indicator.
In every epsilon-wide residue period this indicator has mean one half.
For each displacement, the crossing slab has uniformly bounded boundary
area. Replacing the residue census by its mean costs O(epsilon) times that
area; after summing displacements and normalizing, the error is O(epsilon/a).
The fine cell error is O(h/a), and the smooth surface error is O(a).
The dyadic contribution is at most O(log(q)/q) even for arbitrary labels.
Therefore these actual clipped-box partitions satisfy

    |u_h| -> v,     E_T(u_h) -> (1/2) Per(F)
                       = 2^(-1/3) (36 pi)^(1/3) v^(2/3).       (5)

The last value is strictly below the Euclidean isoperimetric value at
volume v. Volume rounding does not rescue the minimum: the cell count
error is O(m q^2). Correct it to round(v q^3) by changing that many cells
inside C. A changed cell affects at most D_T crossing pairs, so the
normalized error is O(m D_T/q^2)=O(1/K)+o(1), which vanishes.

Connectivity alone also does not rescue compactness. For the ball example,
join the occupied slabs by a one-cell-wide axial path through the ball's
center. The cross-sections of each nonempty slab contain that path, so the
union is connected through nearest-neighbor cells. This changes O(q) labels,
volume O(q^-2), and normalized entropy O(D_T/q^3), all vanishing. The
connected sequence still has weak limit (1/2)1_F and no strong subsequence.
This connectedness statement uses volume convergence; exact-volume minimum
statements above range over all site partitions, without a connectivity
constraint. No regularity or shape is inferred from mere connectivity.

The fixed smooth-section area law is not contradicted. The failure concerns
varying microscopic partitions, which that law did not quantify over.

## 3. An explicit repair and its sharp scale window

Add every integer vector in the short ball B_r to T, retaining set union:

    U = T union B_r,       r<=m/2.                    (6)

No edge receives a special entropy weight. Suppose r tends to infinity,

    r^4/(m q) -> infinity,      r^4/q^2 -> 0.         (7)

These are analyzed parameter schedules, not assumed outcomes. Our final
choice is completely explicit: take n=2t, t>=1, and

    q=2^(16t), m=2^(8t), K=2^(6t), R=2^(14t), r=2^(7t).
    m q/r^4=2^(-4t),       r^4/q^2=2^(-4t).         (8)

The fine-ball degree is Theta(r^3); its contribution on a fixed interface
is O(r^4/q^2) Per(F), and hence vanishes. It nevertheless penalizes the
residue partitions. For the periodic half-residue partition in section 2,

    X_(B_r)(u) = (q^3/m) sum_(|w|<=r) |w_1|
               = (2q^3/m) C_r,
    C_r = sum_(|w|<=r,w_1>0) w_1 ~ pi r^4/4.         (9)

Thus E_(B_r)(u) ~ 2 r^4/(m q). Adding the dyadic overlaps only once
changes this by at most O(log(q)/q). If r^4/(m q) has a bounded subsequence,
these noncompact partitions have bounded repaired entropy on that
subsequence. The same localized partitions in F have bounded entropy and
no strong subsequence in the clipped box. Compactness therefore requires
r^4/(m q)->infinity within this repair class. Section 4 proves sufficiency.
The other condition in (7) retains the original leading area coefficient.

For (8), the degree is

    D_U = |B_K intersect Z^3| + |B_r intersect Z^3|-1+6(t-1)
        = Theta(q^(21/16)).                          (10)

The final term counts precisely the dyadic axis lengths above r and below m.
All other dyadic vectors already occur in B_r. Over a horizon 1/c there
are q/R=q^(1/8) transitions, so logical read incidences are
Theta(q^(71/16)), compared with Theta(q^5) for the original radius-q^(1/2)
dense family. This fixes population q. The radius a=q^(-1/8) has slower
convergence; these counts are not equal-accuracy or native resource savings.

## 4. Compactness proved from the actual reads

Write d(w)=||u_h(.+h w)-u_h||_1 and let E_fine be the ball contribution.
Even though the balls overlap the dyadic part, E_coarse<=E_U and
E_fine<=E_U separately; no double-counted energy is used in the conclusion.

For every |s|<=r/2, average the triangle inequality
d(s)<=d(w)+d(s+w) over integer |w|<=r/2. Both arguments lie in B_r, giving

    d(s) <= C r^-3 sum_(|w|<=r) d(w)
         <= C q r^-3 E_fine.                        (11)

Integer rounding changes only the universal constant for r>=4. Split any
integer translation of length O(m) into O(m/r) such steps. Fractional cell
translations of a piecewise-constant function are convex combinations of
adjacent integer translation differences, coordinate by coordinate. If
P_epsilon is convolution with the normalized epsilon-cube, this proves

    ||u_h-P_epsilon u_h||_1 <= C (m q/r^4) E_fine.   (12)

To control the remaining coarse oscillations, choose a fixed nonnegative
smooth bump phi supported in |y|<1/2 and positive near zero. Normalize
w_v=phi(v/K)/sum_z phi(z/K). Its support and first differences satisfy
w_v<=C K^-3 and |w_v-w_(v-e_i)|<=C K^-4, with support inside B_K for large K.
Let Q u=sum_v w_v u(.+epsilon v), and f_h=P_epsilon Q u_h. Then

    ||u_h-Q u_h||_1 <= C (q/K^3) E_coarse = C a E_coarse.

The derivative of an epsilon-box average is an epsilon-step difference
divided by epsilon, with averaging in the other coordinates. Move that
difference onto w_v by shifting its index. Since the differences of the
weights sum to zero, subtract the unshifted u_h before estimating. Therefore

    ||grad f_h||_1 <= C/(epsilon K^4) sum_(|v|<=K) d(mv)
                    <= C q/(epsilon K^4) E_coarse = C E_coarse. (13)

Also ||f_h||_1<=||u_h||_1, and

    ||u_h-f_h||_1 <= C [a + m q/r^4] E_U.            (14)

The identities in (13)--(14) use m K^4=q^2, not an assumed coercivity
constant. Bounded entropy and (7) give L1 approximation by a uniformly
bounded BV family. BV compactness on the torus, or on a fixed bounded
support in the zero-extended cube, gives strongly convergent subsequences.
Their limits are binary and belong to BV. This proves the promised
compactness, including the stated sharp threshold in the short-ball class.

## 5. Sharp lower bound and binary recovery

Suppose u_h->1_E strongly in L1 and energies have finite lower limit.
Compactness identifies E as a set of finite perimeter. For y in the unit
ball, assign the displacement v=round(K y); retain its displacement cell
only when |v|<=K. Except on a null limiting boundary, v/K->y. Define

    g_h(y) = a^-1 ||u_h(.+epsilon v)-u_h||_1.

For any compact smooth scalar test psi of sup norm at most one, changing
variables in the difference quotient shows its distributional limit is
y dot D1_E. The test-function difference quotient converges uniformly,
so strong L1 convergence suffices; no convergence rate relative to a is
needed. Consequently liminf g_h(y)>=|D_y 1_E|. Fatou on displacement
cells and (3) yield

    liminf E_U(u_h) >= (2/pi) integral_(|y|<=1) |D_y 1_E| dy
                     = Per(E).                     (15)

The last equality is the polar decomposition of the BV derivative and
the same isotropic flux integral as in section 1. The auxiliary edges are
nonnegative and cannot spoil this lower bound.

For recovery, begin with any finite-perimeter E in the admissible domain.
For its continuum indicator, the translation bound
||1_E(.+z)-1_E||_1<=|z| Per(E) makes displacement-cell quadrature errors
O(K^-1) Per(E). The directional BV translation formula and dominated
convergence give E_coarse(1_E)->Per(E). The fine ball contributes at most
C(r^4/q^2) Per(E), and the dyadic edges at most C(m/q^2) Per(E).

Replace E by majority occupancy in every h-cell. The cube BV Poincare
inequality gives L1 error at most C h Per(E). In the clipped case restrict
to cells lying wholly in C; the additional error is O(h) times the fixed
area of its boundary. The binary rounding error changes energy by at most

    C (D_U/q) ||u_h-1_E||_1 = O(D_U/q^2) ->0.        (16)

This follows directly by bounding the change in each translation norm by
twice the L1 error. Thus the recovery is a sequence of actual binary site
partitions, not fractional occupations or hidden edge weights. Equations
(15)--(16) prove the full perimeter Gamma limit in the strong L1 topology,
with the compactness already established rather than assumed.

## 6. Exact volume, finite minima and limiting shape

Fix 0<v<|C| in the clipped case, or 0<v<1 on the torus, and impose exactly
M_q=round(v q^3) occupied sites. Binary recovery differs from this count by
O(q^2). Add or remove that many admissible cells. There are enough cells
for large q because v is strictly between zero and the available volume.
One changed cell changes X_U by at most D_U, so its total normalized cost
is O(D_U/q^2)->0 and its volume error vanishes. The volume-constrained
recovery therefore obeys the exact finite constraint.

Every finite problem has a minimum because its admissible set is finite.
Recovery bounds the minimum above by the continuum perimeter infimum.
Compactness and (15) bound it below by the same value. Every sequence of
minimizers, or of partitions with normalized energy within o(1) of the
finite minimum, has subsequences converging strongly to continuum
perimeter minimizers of volume v. This is a consequence of the constructed
functional, not a presumed success criterion for choosing U.

In the clipped cube, if a ball of volume v fits strictly inside C, the
Euclidean isoperimetric inequality and its equality case identify the
limiting minimum and shape:

    min E_U -> (36 pi)^(1/3) v^(2/3),               (17)

and every limiting minimizer is a ball, up to translation and null sets.
In physical units, for volume V=L^3 v, the value is
(36 pi)^(1/3) V^(2/3). In particular a sphere of radius s has limiting
crossing entropy pi q^4 (4 pi s^2)/(4L^2) before the fixed area conversion.
This is an interface-selection theorem for a declared volume-constrained
cut problem. It does not identify an event horizon or a gravitational law.

## 7. The spacetime and free-field properties remain available

All T edges remain in U, and its maximum radius remains R. Round d/m to z
coordinatewise; take ceil(|z|/(K-2)) coarse steps obtained by rounding the
straight interpolation to z, then correct the residual by signed binary
axis steps. Each coarse increment has norm below K, every residual step
is present, and every step costs one tick. The route has

    |d|/R <= N <= (|d|/m+2)/(K-2)+1+12n,

and lies within 2m of its straight segment. Its physical time excess and
tube width vanish. Zero offsets supply waiting. Strict timelike margins
give eventual reachability; the radius bound excludes spacelike pairs.
Equal cells of weight (a/c)h^3 and almost-everywhere convergence of generated
interval and strict-pair indicators then give the flat diamond volume,
ordering fraction 1/10 and fourth-root count-clock ratios. These measure
conclusions do not depend on the new short edges. They are not physical
clock or native execution statements.

The field consequence is stronger than the availability of a separately
stabilized action. Remove zero, put M2=sum_(v in U)|v|^2 and use exactly

    A_U f = alpha sum_v [f(x)-f(x+hv)],
    alpha=6/(h^2 M2).                                (18)

Every read has the same positive action coefficient. There is no separate
nearest-neighbor coupling. Cubic symmetry gives exact isotropic second
moments; the Taylor error is at most a^2 B4/4, where B4 bounds the fourth
directional derivatives of the test function. Positivity includes the
fixed-zero exterior boundary terms.
The entropy weights in (1) remain one, distinct from the action normalization.
Smooth consistency alone would not justify a thermal limit. The full
spectral analysis is given next.

## 8. One threshold for cut compactness and hidden field modes

Retain (2), let 8<=r<=m/2, r/m->0, and set s=r/K. Observe the exact identity

    r^4/(m q) = s^4.                                (19)

The compactness condition is therefore s->infinity. The same parameter
controls the uniform action's extra modes.

### An elementary ball-symbol bound

For integer rho>=8 and principal theta in [-pi,pi]^3, define
D_rho(theta)=sum_(|w|<=rho) [1-cos(theta dot w)]. Then

    D_rho(theta) >= (rho^3/49152) min(rho^2 |theta|^2,1). (20)

Here is a uniform proof, including frequencies away from zero. The ball
contains the cube [-b,b]^3, b=floor(rho/2), whose side N=2b+1 satisfies
N>=rho and N>=9. For the one-dimensional Dirichlet sum d_N(x),

    N^2-|d_N(x)|^2 = 2 sum_(d=1..N-1) (N-d)[1-cos(d x)].

If |x|<=4/N, restrict to d<=floor(N/4), use
1-cos(d x)>=2d^2 x^2/pi^2, and sum d^2>=N^3/1536. Dividing by
N(N+|d_N|)<=2N^2 gives 1-|d_N|/N>=N^2 x^2/(1024 pi^2).
If |x|>=4/N, the geometric sum and sin(|x|/2)>=|x|/pi give
|d_N|/N<=pi/(N|x|)<=pi/4. Thus in both ranges
1-|d_N|/N>=min(N^2 x^2,1)/16384. Sum over the other two cube coordinates,
and choose a coordinate with |theta_i|>=|theta|/sqrt(3). This proves (20).
It does not replace the entire spectrum by a few sampled modes.

The moment of U obeys M2<=56 m^2 K^5: bound each integer ball by its
containing cube, use r^5<=m^2 K^5, and bound all dyadic second moments by
2m^2. Hence alpha>=q^2/(10m^2K^5). Its moment also has the asymptotic

    M2 ~ (4pi/5)m^2 K^5,                            (21)

since the short-ball and dyadic moments are lower order for r/m->0.

### Resolving the alias cells

For a principal periodic Fourier index k, decompose each coordinate uniquely
as k=m ell+j with -m/2<=j<m/2. The allowed ell range can depend on the
endpoint j; sums may be enlarged to all integer ell and j in estimates.
The coarse phase depends only on j. For ell nonzero,
|2pi k/q|>=pi |ell|/m. Apply (20) to both balls and use
q=m^2, K^4=m^3. When s>=1 this yields, with a fixed positive constant C,

    lambda_U(k) >= C min(|j|^2+s^5 |ell|^2,a^-2).    (22)

Indeed the coarse term bounds min(|j|^2,a^-2). The fine term bounds
min(s^5 |ell|^2, q^2 r^3/(m^2K^5)); its saturation value is
a^-2 s^3>=a^-2. For ell=0 the coarse bound alone suffices.
Constants absorb the fixed factors of pi. The proof is unchanged if s is
bounded below by any fixed positive number, with a constant depending on
that lower bound. The special choice (8) has s^5=2^(5t)>=a^-2=2^(4t).
In its central alias cell the coarse estimate gives ordinary |p|^2 control;
outside that cell the fine estimate gives at least C a^-2. Thus

    C min(|p|^2,a^-2) <= lambda_U(p) <= |p|^2         (23)

for every principal physical wave vector p=2pi k on the unit torus. The
upper inequality is exact second-moment
normalization and 1-cos(x)<=x^2/2.

### The critical limit is an extra oscillator spectrum

Suppose s->ell_0 with 0<ell_0<infinity. For each fixed pair of integer
triples (j,ell), the Fourier index k=m ell+j exists for all large n.
Second-moment expansion using (21) gives

    lambda_U(m ell+j) -> 4pi^2 (|j|^2+ell_0^5 |ell|^2). (24)

For the coarse ball, the integer ell phases cancel and K/m->0 controls
the Taylor remainder. For the short ball, r/m->0 controls that remainder,
and the second-moment ratio gives s^5 |ell+j/m|^2. The dyadic symbol is
at most O(n alpha)->0. These are all the terms of the defined action.

At positive mass, (24) is precisely the oscillator spectrum of a scalar
field on a product of two three-tori of side lengths 1 and ell_0^(-5/2).
This describes spectral multiplicities of the regulator, not six physical
spacetime dimensions. Smooth fixed spatial probes see only ell=0 and would
miss the additional sector.

The complete thermal limit also follows. Set omega=sqrt(c^2 lambda+mu^2),
with c,mu,beta,hbar>0. The thermal log-partition summand is
-log(1-exp(-beta hbar omega)); the normal-ordered energy summand is
hbar omega/(exp(beta hbar omega)-1). A positive mass gap bounds both by
C exp(-gamma omega), after using omega exp(-beta hbar omega/2)<=
2/(beta hbar e) for the energy. Below the saturation level a^-2, (22)
gives summable bounds C exp[-gamma'(|j|+|ell|)] when s stays above a
positive constant. Above it, all q^3 modes together contribute at most
C q^3 exp(-gamma'/a)->0, because a=q^(-1/8). Dominated convergence proves
the full six-index thermal partition, energy and Gibbs entropy limit,
not merely (24) for selected modes.

If s->infinity, every ell nonzero instead escapes to infinite frequency.
The same double-index domination shows its full contribution vanishes;
the ell=0 sector gives exactly the ordinary three-dimensional free thermal
limit. If s->0, each fixed (0,ell) approaches the mass frequency. Keeping
arbitrarily many such distinct ell gives arbitrarily large limiting energy,
so normal-ordered thermal energy diverges. In particular the original
dyadic-only T action has this defect. The absence of the fine ball can be
checked directly: each fixed alias has eigenvalue at most O(n alpha)->0.

Thus the explicit sparse repair suppresses both the noncompact entropy
partitions and the hidden field modes. At critical scale a finite extra
spectral sector survives. Fixed smooth-section and fixed-mode agreements
do not establish either full equivalence.

### The explicit uniform action on the clipped cube

For (8), extend grid vectors by zero to a periodic cube of side two.
The padded periodic form equals the zero-exterior form on the original
cube. Estimate (23) holds on this doubled reciprocal lattice too, with
p=pi k and alias spacing 2m in k. Bounded
energy then bounds the squared L2 mass outside any fixed Fourier radius P by
C/P^2 for all sufficiently large n. The cellwise versions of any fixed
finite collection of discrete Fourier modes converge in L2 to the
corresponding continuum modes. Their orthogonal complement has the stated
tail bound, so finite low-frequency truncation gives L2 compactness. Fatou on
fixed Fourier modes identifies an H1 limit of the zero extensions, hence
an H0^1 limit on the cube.

Consistency on compact smooth tests, the positive-form dual inequality and
their H0^1 density give the full Dirichlet lower form limit; smooth recovery
gives the upper limit. Min-max gives each ordered eigenvalue limit with
multiplicity. Passing eigenvector equations against sampled smooth tests
identifies the limiting eigenspaces and their spectral projections at gaps. The periodic
symbol has only its constant zero mode. A vector supported in one eighth
of the doubled cube has at most one eighth of its squared norm in that
constant mode. Estimate (23) therefore gives a uniform Dirichlet gap.
Compression and the doubled-box eigenvalue count give, for a fixed C'>0,

    lambda_(n,j) >= C' min(j^(2/3),a^-2),           (25)

including j=1 through that gap. Below the saturation level, thermal terms
are bounded by a summable exp(-gamma j^(1/3)); above it, their total is
at most C q^3 exp(-gamma/a)->0. This proves the full Dirichlet thermal
partition, normal-ordered energy and Gibbs entropy limits for mu>=0.

Canonical Gaussian quantization gives the smooth compact vacuum/coherent
Weyl limits as well. The position covariance is bounded by the inverse
gap. For the unbounded momentum covariance, consistency bounds ||A_U f_h||
uniformly for compact smooth f; its tail above eigenvalue Lambda is
O(Lambda^-3/2). Bounded continuous spectral functions and projections at
gaps then give the joint finite smearing limits and the Klein-Gordon
commutator. Rounding each observation time to a layer tick changes its
coefficient by at most hbar (a/c)||f_h||||g_h||, which vanishes. Finite
propagation of the continuum equation gives the limiting response cone.

The exact finite ODE is not an implementation on the layer tick. In fact
the positive normalized action still has trace mean at least 6/a^2, so
centered leapfrog at dt=a/c remains unstable. No native execution or exact
finite causal support is claimed. Thermal Gibbs entropy and the crossing
statistic in (1) remain different observables.

## 9. Imported analysis and specific contribution

BV compactness, the directional translation formula, cube BV Poincare and
the Euclidean isoperimetric inequality are standard mathematical inputs.
The continuum nonlocal perimeter and Gamma-limit setting is developed,
for example, in A. C. Ponce, *A new approach to Sobolev spaces and connections
to Gamma-convergence*, Calc. Var. PDE 19 (2004), 229--255,
https://perso.uclouvain.be/augusto.ponce/05.pdf.

The discrete displacement kernels here are atomic and have a growing residue
structure, so continuum compactness cannot be imported without checking
that structure. Sections 2--4 supply that missing analysis: an exact alias
cut, a clipped-box isoperimetric violation, the sharp r^4/(m q) threshold,
and an explicit sparse repair. Sections 5--6 give the binary recovery and
exact cardinality argument for these actual reads. Section 8 derives the
critical product-torus spectrum and its disappearance under the same repair,
using the uniform action rather than added special coupling weights.
Universal claims do not
rest on finite simulations or the limited Lean reductions.
