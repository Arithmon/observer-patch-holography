# A spatial fluctuation field from one native history

Long occupation averages of the native repair process have a Gaussian field
limit with a calculable spatial action. On the actual twelve-port carrier,
this field has an exact projection onto its record-recovered rank-three
spatial Gram. A port's incident clock gives different axial and transverse
variances on that same band. On the existing joined spherical screens,
eliminating internal carrier variables gives a native Green operator with
a uniform, quantitative approximation bound. Neither construction inserts
a replacement Hamiltonian or identifies native evolution with a relativistic
wave. The supplied torus at the end is a comparison example, separate from
the actual-source results.

This extends the existing [passage geometry](README.md),
[native replica fluctuations](../native_fluctuations/README.md), and
[one-history circulation](../native_circulation/README.md). The new observable
is the time sum of occupations on **one continuing history**. It does not use
independent replicas, impose a Gaussian microscopic state, or replace repair
by a chosen Hamiltonian. The
[incident-field correspondence](../observer_time_fields/README.md) supplies the
stopping construction; the calculation below extracts its occupation
fluctuation action and its spatial convergence conditions.

```sh
python -m model.physical_space_time.occupation_field --check model/physical_space_time/occupation_field_results.json
python -m pytest -q -W error model/physical_space_time
```

## Source, observer access and compatibility

On the original twelve-port carrier, prepare the uniform diagonal law on the
66 binary configurations of total two. Use the existing native endpoint rule,
uniform independent selection of its thirty seams, and fair placement coins.
Both C and E have this diagonal law. Their coherent outputs are not equated.
Stationarity, the preparation, scheduler and readout access are premises.

The vector observer is given the initial occupations and every seam/current
record. It computes each subsequent occupation from the native continuity
equation, then samples the vector either after each attempt or at attempts
incident to port zero. The second clock ignores remote attempts in its count;
the vector decoder receives their currents. A port-only observer with
its initial occupation and incident-current tape can read the scalar at that
port, but cannot reconstruct the complete vector. The execution checks a
counterexample with the same initial state and identical port-zero incident
tape, but different remote occupations at the first tick. Full-field access
is a supplied observation interface, not a hidden-state oracle relabeled as
a local measurement. No intervention or feedback occurs during these averages.

The independent compatibility requirements are conservation of the disclosed
total, continuity-based readback, normalized complete stopped laws, and
agreement between the native configuration process and field predictions.
All 3,960 native branches are checked. An independent 66-state transition
calculation checks the field closure, stationary distribution and lag laws.
All counter moments include every terminal configuration; no event is
postselected. The large torus example reuses the proved local binary identity
and checks 16,384 actual single-particle native branches on its 64-port graph.

## Finite-graph theorem and the source-derived action

The algebra applies to any finite connected simple graph with `V` vertices,
`M` seams and a uniform fixed total `0<q<V`. Its binary native step is the
average of identity and endpoint exchange. Let

```text
eta = occupations - (q/V) 1,
P0 = I - 11^T/V,
chi = q(V-q)/[V(V-1)],
L = graph Laplacian,                  E[eta eta^T] = chi P0.
```

Choose a nonempty fixed set of `m` clock seams. An attempt in that set ticks,
including attempts with zero current. Let `L_i,L_r` be the incident and remote
Laplacians. On density vectors the remote and incident subtransfers are

```text
A = (1-m/M) I - L_r/(2M),
B = (m/M) I - L_i/(2M),
D = I-A,               K = B D^(-1).
```

Every remote word is summed before the next incident attempt. Its tail mass
after `h` attempts is `(1-m/M)^h`, independently of the configuration.
`K` is the exact completed density map, not a mean waiting-time substitution.
The full configuration stopped kernel is bistochastic and keeps the uniform
fixed-total preparation stationary. It is irreducible: each incident exchange
is a one-tick transition; a remote exchange followed by an incident no-change
branch is also a possible one-tick transition. Connected exchanges generate
the fixed-number sector. Positive self-loops give aperiodicity.

On the constant-free space, write `L^+` for the inverse of `L`, zero on the
constant mode. Since `I-K=(L/(2M))D^(-1)`,

```text
H = (I-K)^(-1) = (2m I+L_r)L^+,
Sigma = chi [H+H^T-P0]
      = chi [4m L^+ + L_r L^+ + L^+ L_r - P0].                 (1)
```

The inverses in these formulas are restricted inverses, not inverses on the
conserved constant direction. For the vector sum
`Y_n(u)=n^(-1/2) sum_(j<floor(nu)) eta_j`, the finite-chain functional CLT gives

```text
Y_n -> W,        Cov(W(u),W(v)) = min(u,v) Sigma.              (2)
```

Here `u` is rescaled count of the specified observer events. To check the
imported theorem's premises, use `g(eta)=H eta`. Its Poisson equation is
`g-E[g_next|eta]=eta`; the sum differs from the martingale with increments
`H(eta_next-K eta)` by a bounded endpoint term. The finite state space bounds
these increments. The ergodic theorem makes their averaged conditional
covariance converge. Its stationary value is exactly
`H chi(P0-KK^T) H^T=Sigma`, independently checked in the executable. This
satisfies the bounded-jump case of
[Whitt, Theorem 2.1](https://www.columbia.edu/~ww2040/WWproofs2007pub.pdf).
The probabilistic limit theorem is reused; finite tests do not prove a CLT.

At `u=1`, the field density is proportional to `exp(-y^T Sigma^+ y/2)` on
the constant-free subspace. The receipt contains its exact covariance,
precision and positive rational LDL pivots. On the twelve-port source,
`chi=5/33`. The two clock formulas become

```text
Sigma_attempt  = (5/33)(120 L^+ - P0),
Sigma_incident = (5/33)(20 L^+ + L_r L^+ + L^+ L_r - P0).
```

This is an occupation-history fluctuation action, not the action of the
instantaneous occupation state. At a finite averaging count the state remains
discrete and non-Gaussian. The exact first four moments, computed independently
from complete 66-state histories, retain nonzero fourth cumulants. No finite
discrete law is asserted close to a continuous Gaussian in total variation.

The mean wait is six attempts per port-zero incident tick. Nevertheless,
`Sigma_incident != Sigma_attempt/6`. The difference after division by clock
exposure, `Sigma_incident/5-Sigma_attempt/30`, has entries
`(0,0)=-5/216`, `(1,1)=23/1080`, `(0,1)=61/11880`. Its anisotropy rules out
repairing this discrepancy by a scalar amplitude or clock conversion.
The action precision at this finite source is not a scalar multiple of `L`.

## The actual rank-three spatial field and its clock

Let `P` be the projector recovered from native current calibration. Its unit
Gram is `4P`, with `P^2=P`, `rank(P)=3`, `P_ii=1/4`, and
`LP=lambda P`, `lambda=5-sqrt(5)`. The fluctuation readout is **`P eta` on
the same continuing history**, with the disclosure contract above. No new
preparation, observer replica or coordinate frame is introduced. Define

```text
R = 4 (P e_0)(P e_0)^T,       T = P-R,
P L_i P = [(3-sqrt(5))/2] R + (1/2) T.
```

`R` is the axis through port zero in the finite Gram embedding; `T` is its
two-dimensional transverse plane. This use of an axis does not construct a
bulk radial coordinate. Exact multiplication and the independent 66-state
lag calculation give the marginal Gaussian covariances:

| Native clock | Variance on `R` | Variance on each direction of `T` |
| --- | --- | --- |
| Every attempt | `(145+30sqrt(5))/33` | `(145+30sqrt(5))/33` |
| Port-zero incident attempts | `(5+sqrt(5))/6` | `(115+19sqrt(5))/132` |

If these variances are `v_r,v_t`, the marginal action is
`(1/2)y^T(R/v_r+T/v_t)y` on `range(P)`. Inverting the projected covariance
is essential: projecting the full precision instead would condition the
discarded modes and gives a different incident-clock law.

The attempt density transfer preserves this band. The incident transfer
does not have an autonomous linear three-mode description: its exact
two-tick defect is `P K^2 P-(P K P)^2`, whose `(0,0)` entry equals
`893/492840-31sqrt(5)/54760`, a nonzero value. The corresponding attempt
defect is zero. This proves failure of the proposed projected linear
transition, not a blanket claim about every nonlinear reduced description.

Rank three specifies the dimension of this finite Euclidean embedding and
of the selected fluctuation coefficients. A Gaussian distribution over
three coefficients is not a continuum of spatial sites. These results give
joint operational spatial, field and clock readouts on the actual carrier;
they do not alone determine spatial topology, a physical length, or a
Lorentzian metric.

There is also an exact finite-state boundary. The existing coordinate
readout has 73 distinct values on all 78 total-two states, and 61 on the
66 binary states; the six antipodal binary pairs share the zero value.
All 4,680 native total-two branches satisfy
`delta(P n)=current * P(e_b-e_a)`. Thus cumulative coordinate displacement
telescopes to the final minus initial readout for every native history.
Waiting longer in this fixed domain does not realize the dense additive
completion of unconstrained port directions. The history-average CLT
constructs a field distribution, not additional spatial sites.

The [causal-placement control](causal_placement.py) uses the same native
seams and Gram. Assign each attempt its seam midpoint, and take its parents
to be the last uses of its two endpoints. For an exogenous scheduler, the
complete seam tape defines a dependency height `h`: zero for roots and one
plus the largest parent height otherwise. All 270 ordered pairs of seams
sharing an endpoint obey the sharp squared-distance bound
`c^2=(1+1/sqrt(5))/2`. The triangle inequality therefore gives
`||x_b-x_a|| <= c(h_b-h_a)` whenever `a` is an ancestor of `b`.
This is an exact forward cone placement with the declared time coordinate
`t=c h`; it is not a calibrated incident clock.

The converse fails on actual native histories. For the tape
`A,B,A,A`, with `A=(0,1)` and `B=(2,3)`, event `B` is not an ancestor of
the final event, but their forward interval is null at height difference
two. The receipt executes all sixteen coin branches, retaining zero-current
events. Consequently the derived forward bound does not supply a faithful
identification of causal order with a Lorentz cone. The endpoint DAG also
requires extra edges for adaptive controllers or preparation dependencies;
those are outside this placement certificate.

The execution certificate emits native before/after endpoint values and
writer versions, checks floor/ceiling payloads independently, and
authenticates the endpoint parents. Local support is checked on all 4,680
native branches and is irredundant within the admitted two-endpoint read
interface. The [Lean witness](../../Lean/Geometry/NativeRepairSemanticPlacement.lean)
constructs the executed semantic log for the four-event tape and all sixteen
coin branches, proves canonical heights `(0, 0, 1, 2)`, and applies the
registered forward placement theorem; the converse countercontrol is also
proved. Python compares all current/terminal rows with the literal table
proved in Lean and checks its source-midpoint/scale chart. The general
all-tape cone bound has its separate analytic certificate. Complete
classical value/version logging and zero-current commit stamps are
instrumentation, not intrinsic global read access, arbitrary quantum
histories or a physical influence law. Adaptive scheduling requires its
additional dependency edges.

## Finite counts, caps and observation errors

Stationarity and density closure give an exact finite-count covariance:

```text
Sigma_n = chi [P0 + sum_(h=1..n-1) (1-h/n)(K^h+(K^T)^h)P0].   (3)
```

If `||K^h P0|| <= c r^h`, with `r<1`, summing the geometric series gives

```text
||Sigma_n-Sigma|| <= 2 chi c r/[n(1-r)^2].                    (4)
```

The exact source bounds `2P0<L<10P0` are checked by rational LDL certificates.
For the attempt clock use `c=1,r=29/30`. For the incident clock, `D` lies
between `I/6` and `I/3`. Its similarity to the symmetric positive matrix
`D^(-1/2) B D^(-1/2)` gives `c=2,r=9/10`. These are conservative operator
bounds. Equation (4) is a covariance error, **not a distributional CLT rate**.

A finite implementation that caps each requested interval at `h` attempts
must keep timeout. For `n` requested ticks its total timeout probability is
at most `n(1-m/M)^h`. Coupling to the completed experiment bounds the
complete-output total variation by that probability; expectations of bounded
characteristic functions differ by at most twice it. Fixed-cap success laws
must not be renormalized and treated as completed laws. Increasing caps with
`n(1-m/M)^h -> 0` preserves (2). No auxiliary Poisson clock is required.

## Controlled Green elimination on the actual source screens

The unchanged source builder joins twelve-port carriers on its declared
spherical cell graph. Write its actual port Laplacian as
`L=L_0+D D^T`, where `L_0` is the direct sum of internal carrier Laplacians
and columns of `D` are the oriented cross-seam differences. Every cross-seam
endpoint is distinct, so `D^T D=2I`. The executable checks this matching
condition and every full port-Laplacian row on source levels zero, one and
two: 240, 960 and 3,840 ports. The native two-endpoint exchange identity
binds these matrices to the occupation process.

Let `U` repeat each cell value divided by `sqrt(12)`, `Q=I-UU^T`, and
`R=D^T(sqrt(12)U)`. The exact internal gap is `lambda=5-sqrt(5)>2` on
`range(Q)`. Minimizing the native quadratic energy over those internal
variables gives

```text
A = U^T L U = L_cell/12,
F = I + D^T L_0^+ D,
S = (1/12) R^T F^(-1) R,
U^T L^+ U = S^+,
[lambda/(lambda+2)] A <= S <= A.
```

The identity for `S` follows by solving the internal normal equation
`L_0 q + QD D^T(Uu+q)=0`; substituting its solution leaves the matrix `F`.
The Green identity is block elimination on the global constant-free space.
Since `0<=L_0^+<=Q/lambda` and `D^T D=2I`,
`I<=F<2I`. These statements give the uniform Green comparison

```text
A^+ <= S^+ <= [(lambda+2)/lambda] A^+.
```

Thus replacing each carrier by its arithmetic mean is a bounded
approximation, but it is not exact. A fully rational matrix polynomial
provides arbitrarily small error without discarding the actual routing.
For `p>=1`, set

```text
F_p^(-1) = (1/2) sum_(j=0..p) (I-F/2)^j,
S_p = (1/12) R^T F_p^(-1) R,
0 <= S_p^+ - S^+ <= S^+/(2^(p+1)-1).
```

The tail is `F^(-1)(I-F/2)^(p+1)`, bounded by
`2^(-p-1) F^(-1)`. Sandwiching by `R` and inverting on the mean-free
cell space proves the bound. At order 16 its relative Green error is at most
`1/131071`, uniformly over graphs satisfying these source hypotheses.
The proof bounds exact arithmetic. Floating evaluations at orders 4, 8 and
16, checked against independent full-port solves, are diagnostics; their
roundoff is not included in this analytic error budget. An independent
24-port two-carrier test verifies the elimination exactly over the rationals.
For a declared cell observation matrix `W`, the covariance error is bounded
by `||W^T S^+ W||/(2^(p+1)-1)`. The corresponding limiting Gaussian
characteristic-function error is at most half that bound times `||theta||^2`.

The same history-averaged cell field, divided by `sqrt(4 chi M)`, has
attempt-clock covariance `S^+-P_cell/(4M)`. The executable also evaluates
the exact formula (1) for the clock counting the root carrier's 33 incident
seams, on the disclosed contrast `(N_cell0-N_cell1)/sqrt(24)`. This keeps
the clock correction and the native Green readout together. These cell
labels are not asserted to specify the same continuum smearing at different
levels; no convergence conclusion is drawn from that contrast table.

The local rank-three carrier band and the refining global low band are
different. For any normalized global eigenvector `L psi=omega psi`,
positivity gives

```text
omega = psi^T L psi >= lambda ||Q psi||^2.
```

Consequently, low global modes concentrate on cell averages as `omega->0`;
their nonconstant internal weight is at most `omega/lambda`. The internal
rank-three modes lie at a fixed nonzero carrier energy. The three lowest
positive eigenvalues on the executed screens are approximately
`0.04414..0.04455`, `0.01032..0.01046`, and `0.002555..0.002615`.
Their within-cell nonconstant weights satisfy the bound. These are finite
spectral diagnostics, not a proof of a target spectrum or physical dimension.

### One ancestry readout across the actual source levels

The [ancestry calculation](ancestry_green.py) fixes a common observation
dictionary before comparing levels. Each source cell inherits one of the
twenty base faces from the source's declared refinement tree. At level `l`,
write `s=4^l`, so each base face contains `12s` ports. Define the integer
matrix `B` to have entry `19` on a port's own ancestral face and `-1`
otherwise, and set `W=B/(20 sqrt(12s))`. Thus `W^T W=P20`, where
`P20=I-11^T/20`. The twenty readouts are centered ancestral occupation
counts divided by `sqrt(12s)`. They require the disclosed initial loads and
complete seam/current transcript used throughout this field construction;
they are not inferred from an individual port's tape. The ancestry tree
and the combinatorial operator normalization `s` are declared source data,
not an operational radial coordinate or a selected continuum metric.

Use the stationary binary preparation with total `V/3` on the graph's `V`
ports, hence `chi=2V/[9(V-1)]`. The observation-count limit of the field,
normalized by `1/sqrt(4 chi m s)`, has the following exact covariance on
these twenty observations:

```text
G_l = W^T (s L)^+ W = B^T L^+ B/(4800 s^2),
C_attempt = G_l - P20/(4 M s),
C_incident = G_l + P20/(132 s)
             - W^T (L_i L^+ + L^+ L_i) W/(132 s).
```

Here `M` is the total seam count; `L_i` contains the root carrier's 33
incident seams, including its 30 internal seams. Both clocks sample the
same native law and preparation. The observation-count limit is taken
separately at each finite graph before any spatial comparison.

The executable checks every seam, ancestry assignment and full Laplacian
row at levels zero through three: 240, 960, 3,840 and 15,360 ports. A sparse
solve proposes `X/2^40` as an approximation to `L^+B`; exact integer
arithmetic then checks its residual `E=LX-2^40 B`. If `e` is the cell-graph
eccentricity of its root, the port diameter is at most `D=8e+3` because
each carrier has diameter three. The elementary range/variance inequality
gives `||L^+||<=VD/4=:beta`. Consequently the symmetrized covariance
approximation has certified operator error at most

```text
epsilon = sqrt(4800s) beta ||E/2^40||_F /(4800s^2).
```

The implementation rounds each square root upwards over the rationals.
The incident covariance's error is at most `(13/11)epsilon`, using
`||L_i||<=12`; the attempt covariance has error `epsilon`. Conversion to
retained floating matrices and the receipt's existing entrywise replay
tolerance are included in a common `1e-6` operator enclosure. The numerical
solver is therefore a candidate generator, not the error certificate.

Consecutive Green differences on this fixed dictionary have operator norms
approximately `2.02562`, `0.616605` and `0.182154`. The receipt encloses each
with an exact rational Rayleigh lower bound and row-sum upper bound,
expanded for the retained covariance error. These finite decreases do not
bound the unexecuted refinement tail.

The clock correction has an all-level bound.
Only three cross-carrier seams can affect a cell-constant observation, so
`||W^T L_i W||<=1/(2s)`. Since `L_i<=L`, Cauchy–Schwarz in the `L_i`
quadratic form gives

```text
||C_incident-G_l|| <= [1+2 sqrt(||G_l||/2)]/(132s).
```

On the four executed levels,
`||G_l||<=36` is certified and the incident correction is at most
`19/(264s)`. The following proof supplies boundedness at every level and
compatible subsequential limits on an enlarging observation dictionary.
It does not supply a Cauchy tail or uniqueness of the spatial limit.

### A uniform bound on the actual spherical source

The source refines each triangular face into four children using shared
edge-midpoint identifiers. Normalizing the midpoint coordinates changes
their placement, but not this combinatorial subdivision or its dual cell
graph. Compare that graph with the affine midpoint subdivisions of the
flat regular icosahedron of circumradius one. This auxiliary surface is
used to prove an inequality; it replaces neither native repair nor the
source readout by a geometric Laplacian. Its side length `a`, face area
`a0` and inradius `r` satisfy

```text
a^2=16/(10+2sqrt(5)),  a0=sqrt(3)a^2/4,
r^2=(5+2sqrt(5))/15,  r>3/4.
```

At level `l` each flat triangle has area `a0/s` and side `a/sqrt(s)`.
For a set of `k<=10s` cells with `b` dual edges leaving it, its flat
boundary consists of exactly those `b` primal edges. Its complement has
at least the same flat area. Radial projection of the convex polyhedron
onto the unit sphere is one-to-one except for shared edge descriptions.
On each face its area Jacobian is `r/|x|^3>=r`, and its differential has
norm at most `1/|x|<=1/r`. Thus both spherical areas are at least
`r a0 k/s`, while the spherical boundary length is at most
`a b/(r sqrt(s))`. The finite-perimeter spherical isoperimetric inequality
`P^2>=A(4pi-A)>=2pi min(A,4pi-A)` applies also to disconnected unions
and sets with holes. We import this classical round-sphere case of
[Cavalletti–Mondino, Theorem 1.1 and Section 2.2](https://arxiv.org/html/1610.05044v2);
the cap profile follows from `A=2pi(1-cos(theta))`, `P=2pi sin(theta)`.
The cited theorem uses probability area: substitute `v=A/(4pi)` and
normalized perimeter `P/(4pi)` to obtain the displayed unnormalized formula.
Consequently the cell graph's edge Cheeger constant obeys

```text
b^2 >= (pi sqrt(3) r^3/2) k,
h^2 := min_(0<k<=10s) (b/k)^2 >= pi sqrt(3) r^3/(20s)
                                           >= 27/(256s).
```

For completeness, the graph inequality used here follows from layer cake.
If `g>=0` is supported on at most half the vertices, integration of its
level-set boundaries gives `h sum g_i^2<=sum_edges |g_i^2-g_j^2|`.
Cauchy–Schwarz bounds the right side by
`sqrt(2d sum g_i^2 sum_edges(g_i-g_j)^2)` for maximum degree `d`.
Apply this to the positive and negative parts of a mean-zero vector after
subtracting a median. Their energies sum to at most the original energy;
their squared norms sum to at least the original squared norm. Since the
actual dual graph has degree three,

```text
lambda_2(L_cell)>=h^2/6>=9/(512s).
```

The internal gap `5-sqrt(5)>2` and distinct-port routing give the native
Schur bound `S>=L_cell/24`. Hence for **any** constant-free cell observation
`Z` of operator norm at most one, with port lift `W=UZ`,

```text
G_l[Z]=W^T(sL)^+W=Z^T S^+ Z/s,    0<=G_l[Z]<= (4096/3) I.   (U)
```

This holds at every combinatorial level and for every routing using
distinct ports at a carrier. The source's three endpoint assignments
use distinct columns of its local assignment; no convergence of its
floating routing choices is assumed. The executable checks the base
geometry and subdivision identities exactly, binds the independent
combinatorial reconstruction to the source, and checks the finite native
matrices. The all-level conclusion follows by subdivision induction and
the inequalities above, rather than by enumerating finitely many levels.
It concerns the mathematical continuation of the declared source rule,
not the ability of a finite-precision constructor to allocate every level.

For the twenty-face dictionary, (U) gives the rational clock bounds

```text
||C_incident-G_l|| <= 5/(12s),
||C_attempt-G_l|| <= 1/(2520s^2).
```

Thus these normalized observer-clock fields share every spatial
subsequential limit. No replacement of an event clock by its mean wait
was made. Gaussian characteristic-function differences are at most half
these bounds times the squared argument norm.

### Compatible subsequential field laws and their limits

Fix any coarse level `k`. At each `l>=k`, center the occupation counts
over its `20*4^k` ancestral cells and divide by
`sqrt(12*4^(l-k))`; call this port observation matrix `W_lk`.
The isometry `J_k` repeats each coarse coefficient on its four children
and divides it by two. Exactly,

```text
W_lk=W_l,k+1 J_k,
G_l,k=J_k^T G_l,k+1 J_k,    ||G_l,k||<=4096/3.
```

For each fixed `k`, the covariance matrices lie in a compact finite
dimensional set. A diagonal subsequence of levels therefore converges
for every `k` simultaneously and preserves these identities. This gives
a positive bounded covariance form on the countable union of the
ancestral step-function spaces, hence compatible centered Gaussian
finite-dimensional laws. One can construct the countable Gaussian law
by a Gram representation in a separable Hilbert space and independent
normal variables; no claim of a square-integrable random field follows
from bounded covariance alone. The observation-count CLT is taken first
at each finite graph, followed by this common spatial subsequence.

The attempt correction is unchanged. For the root incident clock,
`||W_lk^T L_i W_lk||<=4^k/(2s)`, giving
`||C_incident,l,k-G_l,k||<=[1+54*2^k]/(132s)`. It vanishes for every
fixed coarse observation space, so both clocks have the same compatible
limiting Gaussian laws along the selected subsequence. This is not a
uniform finite-observation-count distribution bound.

The limit cannot be identically zero. Let `C` contain the twenty flat
base-face centroid vectors, so `C^T C=(20r^2/3)I` and `1^T C=0`.
At level `l`, use the fine flat face centroids as a three-component
cell-constant variational test field. Each adjacent centroid separation
is at most `a/sqrt(3s)`, including edges crossing a base-face hinge.
There are `30s` dual edges, so its energy matrix is at most `10a^2 I`.
Its average over descendants of each base face is exactly that face's
centroid. The variational formula for `L^+` then gives, for every unit
`q` in `range(C)`,

```text
q^T G_l,0 q >= 8r^2/a^2 = (7+3sqrt(5))/3 > 4.
```

These three nonvanishing comparison smearings are not additional spatial
sites or an identified physical coordinate chart. Compactness does not
select a unique covariance, identify a spherical Green operator, establish
a conformal theory, or reconstruct causal Lorentzian dynamics. A whole
sequence convergence estimate and that physical reconstruction remain
separate obligations.

### A common positive spatial generator from the native Green field

The [compactness calculation](local_continuum.py) strengthens the common
subsequence result from weak covariance convergence to **operator-norm
convergence**. Its inverse is one positive self-adjoint continuum spatial
generator shared by the field and the temporal constructions. The source
architecture, routing contract and ancestry observations are unchanged.
This does not select a unique whole-sequence limit or identify a physical
metric. The full-port lifting argument below proves that its limiting form
is a conservative strongly local Dirichlet form.

Let `P_k` be conditional averaging on the level-`k` ancestral cells in
`H=L2_0(mu)`, with the normalized area measure defined below. The level-`l`
cell space has dimension `20*4^l-1`. Write `B_l=4^l S_l` and extend
`G_l=B_l^+` by zero on its orthogonal complement. The quantitative new
ingredient is, for every cell vector `u` and every `l>=k`,

```text
||(I-P_k)u||^2 <= 624*4^(-k) <u,B_l u>.                       (K1)
```

Here is an elementary proof using the actual subdivision. A coarse triangle
has `m^2` descendant cells, `m=2^(l-k)`. Label the upright cells `U(i,j)`,
`i,j>=0, i+j<m`, and downward cells `D(i,j)`, `i+j<m-1`.
Each downward cell joins precisely `U(i,j),U(i+1,j),U(i,j+1)`.
Reflect upright values into an `m` by `m` square by mapping an index with
`i+j>=m` to `(m-1-j,m-1-i)`. Each upright value occurs once or twice;
every horizontal or vertical upright edge occurs exactly twice in the
square energy. The square variance is at most `m^2` times its edge energy:
apply the one-dimensional telescoping inequality to each row and then to
the row means, using Jensen's inequality on the latter differences.

Each upright edge is a two-edge path through a downward cell. The sum of
upright edge energies is at most four times the dual patch energy `E`.
Thus upright squared deviations about the square mean total at most
`8m^2 E`. Pair each downward cell with its unique `U(i,j)`; its squared
deviation is at most twice that upright deviation plus twice their edge
energy. The total variance, minimized over the constant, is consequently
at most `(24m^2+2)E<=26m^2 E`. Summing over coarse triangles discards only
nonnegative boundary-edge energies. Finally `S_l>=L_cell/24` gives (K1).
The exact executable controls bind every descendant adjacency to the
independent triangular stencil, check square edge multiplicities, and
independently certify the finite patch inequality by rational LDL factors.
These controls audit the stencil; the argument proves all refinement levels.

Using `||G_l||<=C=4096/3`, substitute `u=G_l f` into (K1):

```text
||(I-P_k)G_l|| <= sqrt(624 C)*2^(-k),
||G_l-P_k G_l P_k|| <= 2 sqrt(624 C)*2^(-k).                 (K2)
```

For fixed `k`, the middle operator acts on a fixed finite-dimensional space
and converges along the existing common covariance subsequence. Given an
error tolerance, choose `k` first in (K2), then take the subsequence indices
large. This proves that the full operators are norm Cauchy. Their limit `G`
is compact, positive and self-adjoint. The finite-dimensional min-max
principle also gives a quantitative spectral bound: positive covariance
eigenvalue number `20*4^k` and every later one is at most `624*4^(-k)`.
Indeed (K1) bounds the energy below on the subspace perpendicular to `P_k`,
whose codimension is `20*4^k-1`; norm convergence preserves that bound.
The matching lower spectral-count estimate below identifies the spectral
dimension of this particular source limit; it does not identify the observed
spatial dimension.

The limiting covariance has no further null directions. For a Lipschitz
mean-zero test `g` on the auxiliary flat polyhedron, adjacent cell averages
differ by at most `2a*2^(-l) Lip(g)`. Counting the `30*4^l` dual edges and
normalizing the `20*4^l` cell indicators yields
`<P_l g,B_l P_l g><=a^2 Lip(g)^2/2<Lip(g)^2`.
This upper bound also holds for the full cell-constant port lift, whose
energy is `<P_l g,4^l L_cell,l P_l g>/12`. For any unit-axis flat coordinate
`f=e.x`, `Lip(f)<=1` and
`||P_l f||^2>=||P_0 f||^2=r^2/3>3/16`, by the exact base-centroid Gram
above and the inradius bound `r^2>9/16`.
The Green variational identity gives
`<f,G_l f> >= |<f,P_l g>|^2/Lip(g)^2` for nonconstant `g`.
Lipschitz tests are dense in `H`, so for every nonzero `f` one can choose
such a test with nonzero pairing; `P_l g->g` then gives `<f,Gf>>0`.
Consequently `A=G^(-1)` on `range(G)` is densely defined, positive and
self-adjoint. This is an infinite-dimensional generator from the native
source, with compact inverse. Its closed quadratic form has domain
`range(G^(1/2))`. No round-sphere Laplacian was inserted into this argument.
The [joint dynamics](../observer_time_fields/README.md) uses this same
operator limit for temporal correlations and the admitted coherent wave.

There is also a quantitative spatial range statement for the finite action.
Set `T=I-F/2`, so `0<=T<=I/2`. In the exact Schur series
`S_p=(1/24) R^T sum_(j=0..p) T^j R`, the term at power `j` joins only cells
of graph distance at most `j+1`: `F` joins only seams sharing a carrier.
Its remainder obeys `0<=S-S_p<=2^(-p-1)S<=2^(-p-2)I`.
For vectors supported on cell sets at graph distance `d>=1`, choose
`p=d-2` when `d>=2`; use `||S||<=1/2` when `d=1`. Then

```text
|<u,B_l v>| <= 4^l 2^(-d) ||u|| ||v||.                     (K3)
```

Adjacent affine centroids have intrinsic distance at most
`a/(sqrt(3)*2^l)<2^(-l)`, and cell diameters are less than `2*2^(-l)`.
If two fixed test supports have comparison separation `delta>0`, the
supports of their cell averages therefore have graph distance at least
`ceil(delta*2^l-4)` once that expression is positive. Bound (K3) tends to
zero faster than every power of the cell size. The original action thus
has vanishing couplings between separated continuum regions. This proof
alone does not identify the relaxed limiting form: convergence of the energies
of raw cell averages is a separate claim from convergence of their inverses.
In particular, eliminating **cell means** is not vertex Kron reduction.
The exact native level-zero Schur entry `(5,9)` is positive, so its negative
is not a Markov jump generator. Neither an M-matrix premise nor finite
propagation for the limiting wave follows from (K3) alone.

### A local conductance graph with the same continuum Green limit

The actual three routed ports per carrier admit a different elimination.
It uses ordinary terminal Kron reduction and star–triangle equivalence
([Dörfler–Bullo, Theorem III.8 and Lemma III.10](https://arxiv.org/abs/1102.2950)),
with the native constants and observation error derived here. For carrier
Laplacian `L_ico`, its exact resistance
`R_ab=(e_a-e_b)^T L_ico^+ (e_a-e_b)` is `11/30`, `7/15` or `1/2` at native
hop distance one, two or three. For distinct routed ports `i,j,k`, set
`r_i=(R_ij+R_ik-R_jk)/2` and cyclically. Eliminating the other nine ports
gives precisely the terminal energy
`min_z sum_i (x_i-z)^2/r_i`. All 220 carrier triples are checked by an
independent exact Schur elimination. Their sorted star arms are

```text
(11/60,11/60,11/60), (2/15,7/30,7/30), (7/30,7/30,7/30),
(11/60,11/60,17/60), (1/6,1/5,3/10).
```

Thus `2/15<=r_i<=3/10`. Across each actual cross seam `(c,i)-(d,j)`,
eliminating its two terminal values at fixed star centers puts the two
arms and the unit seam resistance in series:

```text
w_cd = 1/(1+r_ci+r_dj),                 5/8 <= w_cd <= 15/19,
W_l = (4^l/12) L_w,l,                  G_l^w = W_l^+ .       (W1)
```

`L_w,l` is a positive local weighted Laplacian on the existing cell dual;
the coefficients use the actual routed port labels, with no fitted action.
Cell and port `L2` norms below use uniform weights `1/C` and `1/(12C)`,
where `C=20*4^l`. Write `s=4^l` and
`E_port(x)=s/(12C) sum_native_edges (x_a-x_b)^2`.

Two variational maps compare these fields with native cell means. Given
any native port vector, set each star center to
`z_c=sum_i x_ci/r_ci / sum_i 1/r_ci`. Terminal energy minimization and
then series minimization give `E_w(z)<=E_port(x)`. Since this center is
a convex combination of port values, the internal gap `5-sqrt(5)>2` gives
`||z-mean_cell(x)||^2<=6 E_port(x)/s`. Conversely, for any center vector
`z`, use seam current `I_cd=w_cd(z_c-z_d)`, terminal values
`x_ci=z_c-r_ci I_cd`, `x_dj=z_d+r_dj I_cd`, and the harmonic lift to
the nine other ports. Its native energy is at most `E_w(z)`. The harmonic
maximum principle and `r_i<=3/10` give
`||mean_cell(x)-z||^2<=(18/5) E_w(z)/s`. Subtracting global means keeps
both bounds and all energies. These are maps of quadratic energies and
readouts, not an equality of finite native transition or observation laws.

The existing dual Poincare bound and (W1) give `||G_l^w||<=16384/15`;
the native bound is `||G_l||<=4096/3`. In
`<f,Gf>=sup_u(2<f,u>-E(u))`, apply either map to the maximizing vector
and bound its readout error by the square root of its energy. For unit `f`,
the two directional errors are at most
`2 sqrt(6*4096/3)/sqrt(s)` and `2 sqrt((18/5)*16384/15)/sqrt(s)`.
Polarization for self-adjoint operators therefore gives

```text
||G_l-G_l^w|| <= 182*2^(-l).                                (W2)
```

The finite source controls check both energy maps in rational arithmetic
and audit the distinct Green operators numerically. Consequently the
explicit local conductance graph and native cell-mean field have exactly
the same operator-norm subsequential limits; whole-sequence uniqueness
for either would establish it for both. The conductances are not constant:
the actual source already uses inequivalent terminal stars and several
edge weights. Carrier symmetry alone does not remove this dependence.
For a concrete control, assigning a distance-one equilateral port triple
at every cell gives `w=30/41`; a distance-two triple gives `w=15/22` on
the same dual graph. Their weighted Green operators differ by the exact
factor `44/41`, which (W2) preserves in any common nonzero continuum limit.
Both retain the carrier and repair law but change the declared routing;
this distinguishes topology-only inference from the fixed source question.
Convergence of this routing-dependent conductance family to one effective
form is unproved. The reduction supplies no extra spatial direction,
physical metric, changed repair rule or finite Markov approximation claim.

### Locality and spectral dimension of the source continuum

The continuum form inherits contraction and locality from the **full native
port energy**, even though finite cell-mean elimination is not Markov.
Extend the mean-zero operator `A` to `L2(mu)` by adjoining the constant zero
mode. Write `E(f)=||A^(1/2)f||^2`; the reduced inverse is the preceding `G`.
The finite forms are `E_l(f)=<f,B_l f>` on the cell step functions, and
infinite outside them. Constants have zero energy at every level.

First these forms converge in the Mosco sense along the selected subsequence.
For weakly convergent `f_l->f`, their variational representation gives

```text
E_l(f_l) >= 2<f_l,g>-<g,G_l g>  for every mean-zero g.
```

Norm Green convergence and then the supremum over `g` give the energy
liminf inequality. For `f=Gg`, the recovery sequence `f_l=G_l g` converges
strongly and has energy `<g,G_l g>->E(f)`. Finite spectral truncations are
dense in the form norm; a diagonal choice supplies recovery for every
finite-energy `f`. Add its constant part unchanged. Thus convergence of
the inverses controls the relaxed energies, not just covariance pairings.

For a finite cell field, take its full-port energy minimizer
`v_l=U_l f_l+q_l`, where `U_l` is the normalized cell isometry and `q_l`
has zero mean within each carrier. The native internal gap gives

```text
4^l <v_l,L_l v_l> = E_l(f_l),
||q_l||^2 <= E_l(f_l)/(2*4^l).                              (L1)
```

Consequently every bounded-energy recovery has negligible internal
fluctuation in `L2`. For any real normal contraction `C` (`C(0)=0`,
`|C(x)-C(y)|<=|x-y|`), apply it to the scalar port values:
`v_p -> C(sqrt(V_l)*v_p)/sqrt(V_l)`. Each native edge energy decreases.
This contraction differs from applying `C` to the cell values by at most
`||q_l||` in norm. Average the contracted ports back to cells and use the
Mosco liminf inequality. It follows that `E(C(f))<=E(f)` on the full
`L2(mu)` space. The closed symmetric form is therefore a Dirichlet form;
`1` belongs to its domain with zero energy, so its heat semigroup is
conservative. Explicitly, the resolvent minimizes
`E(u)+t^(-1)||u-f||^2`; for `0<=f<=1`, clipping `u` to `[0,1]` lowers
both terms, so the minimizer lies in that interval. Euler products of
these resolvents give the same interval-preservation for the heat
semigroup, and the constant zero mode gives conservation.
The finite positive Schur entry is compatible with this
limit, because the internal correction in (L1) vanishes.

Locality needs more than this contraction argument. Let `chi` be a fixed
Lipschitz cutoff on the auxiliary flat polyhedron, with `0<=chi<=1`,
sampled at cell centroids and repeated at all twelve ports of each cell.
For each native seam, expand the difference of `chi*v` into the difference
of `v` plus the difference of `chi`. Internal seams have no cutoff
difference. Each cross-seam endpoint occurs once, and adjacent centroids
are less than `2^(-l)` apart. For every `epsilon>0`,

```text
4^l <chi*v,L_l chi*v>
 <= (1+epsilon) 4^l <v,L_l v>
    +(1+1/epsilon) Lip(chi)^2 ||v||^2_transition.            (L2)
```

The transition norm uses only endpoints of edges on which `chi` changes.
For a cutoff equal to one on a neighborhood of `supp(f)` and changing
only at positive distance from that support, this norm tends to zero on
any full-port recovery for `f`, by (L1) and strong convergence. Choose
`epsilon_l->0` slowly enough that its last term also tends to zero.
This gives recovery sequences supported in any fixed neighborhood of
`supp(f)`, with asymptotically unchanged energy.

For two form-domain functions with supports at positive distance, choose
such recoveries with disjoint support neighborhoods. No native seam joins
them at sufficiently fine levels. The **full-port trial energies** of
both their sum and difference are exactly additive. The projected cell
energies are bounded above by those trial energies; no exact additivity
of finite Schur energies is assumed. Mosco liminf therefore gives
`E(f+g)<=E(f)+E(g)` and `E(f-g)<=E(f)+E(g)`, hence `E(f,g)=0` by
polarization. Subtract a global constant to cover the case that `g` is
constant on a neighborhood of `supp(f)`. This proves **strong locality**
on the comparison topology.

The form also has a regular Sobolev core. On each fine primal triangle,
assign its cell value to the centroid, the mean of the two adjacent cell
values to each edge midpoint, and the mean of the incident cell values to
each primal vertex. Affine interpolation on the six resulting triangles
defines a continuous function `I_l f_l`, including across original face
hinges. Every primal vertex has valence five or six: midpoint refinement
preserves old valences and creates valence-six vertices. Around either
cycle, `sum_i (u_i-mean(u))^2<=sum_i (u_i-u_(i+1))^2`. The executable
proves these two finite inequalities by exact rational positive-semidefinite
factorizations, and binds the vertex stars to the source refinement.

For one small triangle with vertex, edge midpoint and cell centroid
values `v,m,c`, exact flat integration gives
`integral |grad I_l f_l|^2=(sqrt(3)/6)(v-m)^2+(sqrt(3)/2)(c-m)^2`.
Summing all small triangles gives `sum(c-m)^2=E_dual` and
`sum(v-m)^2<=4 E_dual`, where `E_dual` is the unnormalized dual-edge
energy of the scalar cell values. Dividing by total area `20a0`, with
`a0=sqrt(3)a^2/4`, and using the normalized cell dictionary gives

```text
integral |grad I_l f_l|^2 dmu
 <= 14/(3a^2) <f_l,4^l L_cell,l f_l> <= 112/a^2 E_l(f_l),
||I_l f_l-f_l||^2 <= 5/(18*4^l) <f_l,4^l L_cell,l f_l>
                  <= 20/(3*4^l) E_l(f_l).                  (L5)
```

For the second bound, convexity bounds each interpolated squared error by
the barycentric average of its three nodal squared errors. Summed nodal
errors contribute at most `E_dual+4 E_dual`; each small triangle has
normalized area `1/(120*4^l)` and each barycentric coordinate averages
`1/3`, giving `E_dual/(72*4^l)`. No boundary value is pinned.

Take `H1` on the glued flat surface to be the completion of continuous
piecewise-affine functions in the norm `||f||^2+integral |grad f|^2 dmu`.
The gradient is closable in `L2`, as integration against tests in the
interior of each face verifies, so this completion embeds in `L2(mu)`.
Applying (L5) to bounded-energy recovery sequences and taking a weak
`H1` subsequence gives `D(E)` contained in `H1`. Conversely, the
piecewise-affine energy upper bound proved below extends to its `H1`
completion by closedness of `E`. Therefore

```text
D(E)=H1,
(a^2/112) integral |grad f|^2 dmu <= E(f)
 <= (a^2/48) integral |grad f|^2 dmu.                       (L6)
```

The equivalent norms make continuous piecewise-affine functions a form
core. Their vertex interpolants approximate every continuous function
uniformly as mesh diameter vanishes, so the Dirichlet form is **regular**.
Its exact intrinsic metric, routing-independent uniqueness and physical
length/time identification have not been selected by these comparisons.

There is a quantitative two-dimensional spectral conclusion independent of
those further claims. On each of the `N_k=20*4^k` coarse triangles, let
`b=3 min(lambda_1,lambda_2,lambda_3)` in barycentric coordinates, and zero
outside that triangle. This continuous tent has normalized triangle
moments `mean(b)=1/3`, `mean(b^2)=1/6`, and gradient squared
`12*4^k/a^2` almost everywhere. For a fixed continuous piecewise-affine
field `f`, the dual triangular mesh satisfies

```text
lim_l 4^l <f_l,L_cell,l f_l> = (a^2/4) integral |grad f|^2 dmu,  (L3)
```

where `f_l` is its centroid-sampled cell step function with the normalized
indicator coefficients. Within each affine piece, the three neighbor
centroid vectors obey `sum_j d_j d_j^T=(a^2/(2*4^l)) I_2`; counting each
edge at both ends gives (L3). Only `O(2^l)` fine cells meet the finitely
many affine break lines or base-face hinges. Each such edge contributes
`O(4^(-l))` before the fixed factor `4^l/(20*4^l)`, so their total
contribution vanishes. No round-sphere Laplacian is substituted for `S_l`.
The upper native comparison `S_l<=L_cell,l/12` and the energy liminf
give `E(f)<=(a^2/48) integral |grad f|^2 dmu` for these tests.

For a linear combination of the `N_k` tents with coefficients summing to
zero, its squared norm is `sum_i c_i^2/(6N_k)` and its energy is at most
`4^k sum_i c_i^2/(4N_k)`. This is an `N_k-1` dimensional mean-zero
subspace, so min-max gives `lambda_(20*4^k-1)(A)<=(3/2)*4^k`.
Together with (K1), monotonicity and `lambda_1(A)>=3/4096`, every positive
eigenvalue, counted with multiplicity from `n=1`, satisfies

```text
n/49920 <= lambda_n(A) <= 3n/2,
1/(exp(3t/2)-1) <= Tr_H0 exp(-t A) <= 1/(exp(t/49920)-1).    (L4)
```

Here `t>0` is the common generator's scaled native-decay parameter, not
laboratory time.
Thus the heat-trace spectral dimension
`2 lim_(t downarrow 0) log Tr_H0(exp(-t A))/log(1/t)` is exactly **two**.
This concerns the native spherical observation domain and normalization.
It excludes unitary identification, up to a fixed positive scale and
bounded mass shift, with a second-order positive elliptic operator on a
smooth compact three-dimensional manifold without boundary and of positive
volume: the latter
has counting growth `N(Lambda)~c Lambda^(3/2)` by the classical Weyl law
([Rivera, §1, recalling Hörmander's theorem](https://arxiv.org/html/1801.07598#S1)).
The exclusion does not cover a different radial observation algebra,
noncompact volume limit, different operator order or different scaling.
Rank-three internal coefficients alone do not change this spectral count.

### A boundary trace test for a radial extension

A three-dimensional **classical static Dirichlet-energy bulk** has a
different boundary action from its bulk action. Minimizing
`integral_bulk |grad u|^2` at fixed boundary value `f` gives the harmonic
extension and boundary energy `<f,Lambda f>`, where `Lambda` is the
Dirichlet-to-Neumann operator. On a ball of radius `R`, a degree-`j`
spherical harmonic extends as `(r/R)^j Y_j`, with boundary precision `j/R`
and multiplicity `2j+1`. Thus its positive eigenvalues are exactly
`mu_n=floor(sqrt(n))/R`, indexed from one. On any smooth compact
three-dimensional Riemannian bulk with smooth boundary, the Steklov Weyl
law gives `mu_n=O(sqrt(n))`; see
[Girouard–Polterovich, Example 1.3.2 and §2.1](https://arxiv.org/abs/1411.6567).

This supplies a precise test of a proposed boundary interpretation. Suppose
the centered native boundary readout and this scalar boundary trace are
related by a boundedly invertible `L2` map `T`, with the **form congruence**

```text
E_A(f)=c E_Lambda(Tf)+<f,Bf>,    c>0, B bounded self-adjoint.
```

Take the inverse image under `T` of the first `n` Steklov eigenfunctions.
On that `n`-dimensional subspace the Rayleigh quotient is at most
`c ||T||^2 mu_n+||B||`; min-max gives the same upper bound on `lambda_n(A)`.
It contradicts `lambda_n(A)>=n/49920` for large `n`. Consequently this
ordinary scalar boundary-trace identification fails, including bounded
calibration corrections. This argument is about energy congruence, not an
arbitrary similarity of operators. It excludes neither quantum or interacting
bulk reconstructions nor unbounded, differently scaled observation maps.

One permissible mathematical change illustrates the distinction: on a finite
spectral cutoff, applying `A^(1/4)` to the native Gaussian changes covariance
`A^-1` to `A^-1/2`. The distributional spectral limit has the order appropriate
to the harmonic half-cylinder extension `U(r)=exp(-r sqrt(A))f`, whose
minimized action is `<f,sqrt(A)f>`. This is unbounded `L2` processing; the
additional radial domain and its action are construction inputs. It supplies
no native bulk fluctuations or physical radial readout by itself.

The internal carrier directions also cannot supply such a radius at the
present scale. Set `E_l^port(v)=4^l<v,L_l v>`. The internal gap gives
`||Q v_l||^2<=E_l^port(v_l)/[(5-sqrt(5))*4^l]`. Every bounded-energy family
therefore loses its within-carrier component, including the internal
rank-three band, as `O(2^(-l))` in norm. This statement concerns the current
full-port energy scaling, not all possible native observable limits.

### Tilted shell spectra under the same comparison

The min-max argument of the boundary trace test applies to every isotropic
shell precision. Let `K_theta` be a positive operator on the round sphere
whose degree-`l` harmonics have eigenvalue proportional to
`(l(l+1))^(1+theta/2)`, so its counting function satisfies
`mu_n ~ c_theta n^(1+theta/2)`. Suppose the centered native form and this
precision are related by a boundedly invertible `L2` map `T`:

```text
E_A(f)=c E_K(Tf)+<f,Bf>,    c>0, B bounded self-adjoint.
```

Min-max in both directions gives `lambda_n(A)<=c ||T||^2 mu_n+||B||` and
`mu_n<=||T^-1||^2 (lambda_n(A)+||B||)/c`. With `n/49920<=lambda_n(A)<=3n/2`
from (L4), `mu_n` grows linearly, so `theta=0`. A field whose angular
precision is boundedly congruent with the native form therefore has the
scale-invariant angular spectrum of a two-dimensional Gaussian free field:
`l(l+1)C_l` is bounded above and below. In particular, a shell field with
tilt `theta=1-n_s` different from zero cannot be a bounded congruent image of
the native long-average field; a tilt needs unbounded spectral processing or
a separately derived scale map. This is a statement about the comparison,
conditional on the identification it names, and no sky prediction.

### A continuum random distribution on the carrier sphere

The compatible law extends beyond the countable ancestry dictionary.
Let `mu` be normalized area on the auxiliary flat icosahedron, pushed onto
the unit sphere by radial projection, and let `omega=dOmega/(4pi)`.
Every level-`l` cell has `mu`-mass `1/(20*4^l)`. Its normalized indicators
are therefore an orthonormal dictionary, with the same child isometry
`J` used above. Their union is dense in `L2(mu)`: their diameters tend to
zero and continuous functions are dense. Centering removes the constant.
The radial Jacobian gives

```text
dmu/domega = pi |x|^3/(5 a0 r),       1/2 < dmu/domega < 2.
```

Indeed `4/9<a0<1/2`, `r>3/4`, `r<=|x|<=1` and `3<pi<22/7` give the
stronger lower and upper bounds `27/40` and `66/35`. These are comparison
measures for expressing native readouts, not a selected physical metric.
The finite Green covariance, extended by the ancestry projection to
`L2(mu)`, has norm at most `C=4096/3`. The bound for the **full** fine
dictionary is obtained by setting `k=l` above:

```text
||C_incident,l-G_l|| <= (1+54*2^l)/(132*4^l) -> 0.
```

Thus all three normalized covariance forms have the common bound
`M=C+5/12=5463/4`; their constant mode vanishes. Uniform boundedness
extends the common subsequential covariance uniquely from ancestral step
functions to all `L2(mu)` tests. For example, approximating `f,g` by their
ancestry projections changes the covariance by at most
`M (||f-P_k f|| ||g|| + ||P_k f|| ||g-P_k g||)`, uniformly in finer levels.
This proves convergence on every fixed pair of square-integrable tests.

Use real spherical harmonics normalized in `L2(omega)` to define the
distribution norm `||X||_-2^2=sum_(j,m) [1+j(j+1)]^-2 |X(Y_jm)|^2`.
Their orthogonality and eigenvalues are the standard
[spherical-harmonic identities](https://dlmf.nist.gov/14.30).
Each coefficient has variance at most `2M`. For `j>=1`,

```text
(2j+1)/[1+j(j+1)]^2 <= 1/j^2 - 1/(j+1)^2,
E ||(I-P_<=J) X_l||_-2^2 <= (5463/2)/(J+1)^2.             (D)
```

The omitted `j=0` coefficient is zero by conservation. Finite harmonic
blocks have convergent Gaussian covariances along the same subsequence.
Bound (D) controls the remaining norm uniformly. To see tightness directly,
choose cutoffs `J_k` so that the sum of the Markov bounds for tails larger
than `1/k` is arbitrarily small, and bound the first finite block. The
intersection of these constraints has compact closure in `H^-2(S2)`.
Consequently the finite Gaussian fields converge weakly, along that
subsequence, to a centered **Gaussian random distribution in `H^-2(S2)`**.
The Green, attempt-clock and root-incident-clock fields share its law.
Its covariance on `L2(mu)` is the preceding extension. In particular the
three nonzero coarse smearings persist as mean-square limits of smooth
test observations; they need not be pointwise pairings of a distribution
with discontinuous cell indicators.

There is also an exact spatial approximation budget. Transport the three
flat coordinate functions to the sphere. Uniform triangle integration
gives `sum_i ||f_i-P_l f_i||_L2(mu)^2=a^2/(12*4^l)`.
Their summed mean-square field error is at most `1821/(14*4^l)` under the
common bound `M`. The executable independently checks triangle moments and
full native port-Laplacian covariances for these continuous readouts.
Neither these auxiliary coordinates nor (D) establish three physical
spatial dimensions. This result supplies a continuum distribution on the
given carrier topology, without a pointwise field, an `L2` random function,
a unique continuum covariance, conformal invariance or Lorentzian dynamics.
The observation-count limit precedes spatial refinement throughout.

### One continuum observation path for both native clocks

The static limit extends to a joint observation process, and the two clocks
can be coupled on the **same native history**. This uses the stationary
fixed-total binary preparation, IID uniform seam selection and fair placement
coins. Each seam acts as `(identity+endpoint swap)/2`, so it separately
preserves the uniform sector law. The full field retains the initial-load
and complete current-tape access contract above; an incident tape alone
does not disclose remote occupations. No feedback is included in this result.

Write `X_j` for the post-repair configuration, `eta_j` for its centered
occupation, and `I_j` for membership of the selected seam in a fixed clock
set of size `m`. Let `p=m/M`. Conditional on `X_j`, the mark has probability
`p`: this follows from the separate stationary preservation of each seam.
Future repairs are independent of that mark given `X_j`. Consequently

```text
Z_j=(I_j/p-1) eta_j
```

is a reverse martingale difference. For a finite history ending at `N`,
use the decreasing information sets generated by `X_j` and all later
states and marks. `Z_j` is measurable one step earlier and has zero
conditional mean given the current set. Orthogonality gives the exact,
finite-count identity

```text
Cov(sum_(j=1..N) Z_j) = N chi (1-p)/p P0.                    (J1)
```

This is stronger than agreement of two marginal covariances. On the same
attempt prefix, the weighted incident sum is
`sum I_j eta_j/p = sum eta_j + sum Z_j`. Let `R=B/p-P`, where `P` is the
attempt density map and `B` its incident part. For `h>=1`,
`E[eta_h Z_0^T]=0` and `E[Z_h eta_0^T]=chi R P^(h-1)`.
Thus the asymptotic cross covariance is `2 chi M L^+ R`, and adding it,
its transpose and (J1) to `Sigma_attempt` gives exactly
`Sigma_incident/p` from (1). The executable checks this identity and all
two-step moments against all 237,600 actual endpoint/coin histories of the
twelve-port source. Reusing one random seam at the second attempt preserves
the uniform one-time state but gives `Tr E[Z_2 Z_1^T]=250/33`, instead of
zero. Stationarity without the stated scheduling hypothesis is insufficient.
A fixed nonstationary initial configuration also breaks conditional centering.

Let `s=4^l` and project onto the normalized cell dictionary. At level `l`,
define piecewise-linear cumulative paths, with sums evaluated at integer
attempt prefixes and linearly interpolated between them:

```text
A_(l,n)(t) = U^T sum_(j<=nt) eta_j / sqrt(4 chi M s n),
J_(l,n)(t) = U^T sum_(j<=nt) I_j eta_j / [p sqrt(4 chi M s n)].
```

Here `t` is rescaled observation count, not physical duration. Both paths
refer to the same history boundary. The canonical root clock has `m=33`,
`M=630s`. By (J1), the difference covariance per unit observation time on
the mean-free cell space is

```text
d_l P_cell,       d_l=(1-p)/(4 m s)<1/(132s).                 (J2)
```

There is a finite-count path bound. Reverse martingale **suffix** sums obey
Doob's inequality. A prefix is the total minus a suffix; the L2 triangle
inequality gives `E max_prefix ||sum Z||^2 <= 9N E||Z_1||^2`, rather than
incorrectly applying the forward constant four to prefixes. Apply this to
the finite-dimensional image in `H^-2`. The density bound `dmu/domega<2`
and the harmonic sum used in (D) give, for `T>0`,

```text
E sup_(t<=T) ||J_(l,n)(t)-A_(l,n)(t)||_-2^2
 <= 18 ceil(nT)/n d_l.                                      (J3)
```

Doob's inequality here is the ordinary nonnegative-submartingale result
applied to a Hilbert norm; see
[Acciaio et al., Theorem 1.1](https://arxiv.org/pdf/1202.0447).
The bound includes every native branch; there is no success conditioning.

At each fixed level, the finite-chain functional CLT applies jointly to
the bounded transition rewards `(eta_j,Z_j)`. For `H=2M L^+`, add endpoint
differences `(H-P0)(eta_j-eta_(j-1))` to the attempt reward and
`R H(eta_j-eta_(j-1))` to `Z_j`. Both resulting innovations have zero
conditional mean given the previous state. Their summed endpoint corrections
are bounded at each fixed level. Ergodicity supplies the joint predictable
covariance limit, so the bounded-jump martingale theorem used in (2) applies.
The executable checks every conditional mean and the complete joint
innovation covariance on all 3,960 native branches. The limiting pair
`(A_l,J_l)` is a jointly
Gaussian Wiener process. In particular, its difference is a martingale,
and the sharper limiting bound is

```text
E sup_(t<=T) ||J_l(t)-A_l(t)||_-2^2 <= 8T d_l -> 0.           (J4)
```

This also compares the actual incident-index clock. If `tau_k` is the
attempt containing its kth tick, the Bernoulli strong law gives
`tau_floor(pnt)/n -> t`, uniformly on bounded intervals, at each fixed
level. Joint functional convergence and continuity make this random time
change negligible for the centered sums. The incident-index normalization
is `1/sqrt(4 chi m s p n)`, equal to the prefactor in `J_(l,n)`.
Centering matters: a constant uncentered readout has nonzero count noise
at a fixed attempt prefix but none at a fixed tick count. The limit
`n -> infinity` therefore precedes refinement; no uniform Anscombe or
distributional CLT rate as `p_l -> 0` is asserted. Finite waiting caps retain
their timeout and use the complete-output allowance stated above.

Finally, every level's Gaussian path has covariance
`min(t,u) C_l`. Doob's inequality and (D) bound its angular path tail by

```text
E sup_(t<=T) ||(I-P_<=J) A_l(t)||_-2^2
 <= 4T harmonic_tail(J) = 10926 T/(J+1)^2.                  (J5)
```

For a centered Hilbert Gaussian increment with covariance `Q`,
`E||G||^4=(Tr Q)^2+2Tr(Q^2)<=3(Tr Q)^2`. Hence the fourth temporal
increment moment is at most `3(5463/2)^2 |t-u|^2`. Finite harmonic blocks
converge in continuous-path space: their covariance square roots converge
and can multiply the same finite standard Brownian path. Bound (J5) uniformly
approximates the full paths by those blocks. Along the same spatial
subsequence as the static field, they converge in
`C([0,T],H^-2(S2))` to a Gaussian Wiener field with covariance
`min(t,u) Q`, where `Q` is the native limiting spatial covariance. By (J4)
the two normalized native clocks converge **jointly to the same path**.
Any fixed finite set of admitted clock-seam sets has the same conclusion,
using its own `d_l<1/(4m_l s)` and the sum of the comparison bounds.

The finite-time Gaussian action is derived as well. On any nondegenerate
finite observation subspace with covariance `C`, observations at
`0<t1<...<tr` have action
`(1/2) sum_i (y_i-y_(i-1))^T C^-1 (y_i-y_(i-1))/(t_i-t_(i-1))`,
with `y0=0`. Its precision is the temporal difference precision tensor
the native spatial precision. The executable checks a 33-dimensional
three-time instance exactly. Native fields at finite counts have correlated
increments; this action belongs to the proved long-observation limit.

This gives a common space-and-observation-time continuum on the source's
carrier topology, including the nonzero three coarse smearings. It is a
Gaussian process of accumulated observations. It does not identify the
instantaneous density dynamics, a physical radial coordinate, a unique
whole-sequence spatial covariance, Lorentzian propagation, modular flow
or an interacting quantum vacuum. These require their respective dynamical
and physical correspondence conditions, rather than a reinterpretation of
the observation parameter.

## A quantitative spatial-limit criterion

For a graph family indexed by `a`, supply a positive spatial scale `s_a`,
fixed-total preparation `chi_a>0`, and a declared observation matrix `W_a`
whose columns lie in the constant-free space. These columns specify spatial
smearings, including their normalization; neither a plot nor a graph rank
chooses them. Normalize the history average by `1/sqrt(4 chi_a m_a s_a)`.
Its limiting covariance on those observations is exactly

```text
C_a = W_a^T (s_a L_a)^+ W_a + E_a,
E_a = W_a^T (L_r L_a^+ + L_a^+ L_r - P0) W_a/(4 m_a s_a).    (5)
```

If the first term differs from a target covariance `Gamma` by at most
`epsilon_G`, and `||E_a||<=epsilon_clock`, then
`||C_a-Gamma||<=epsilon_G+epsilon_clock`. For a Fourier test vector `theta`,
the corresponding Gaussian characteristic functions differ by at most
`||theta||^2 (epsilon_G+epsilon_clock)/2`. Equation (4), multiplied by
`||W_a||^2/(4 chi_a m_a s_a)`, also bounds the finite-count **covariance**
error. These are observable matrix quantities that can be measured or
certified; they do not assume that agreement of means establishes a field law.

Thus Green convergence and vanishing clock correction give convergence of
the limiting Gaussian finite-dimensional laws. The order is explicitly
`lim_(graph) lim_(observation count)`. No uniform simultaneous limit or
uniform finite-count distributional error has been proved. A fixed total
two on growing graphs has `chi_a -> 0`; its amplitude normalization cannot
be carried over from the twelve-port value unchanged.

## Comparison control: unchanged repair on a supplied two-dimensional torus

Choose the `N by N` periodic square graph, a chart of period `2pi` in each
direction, `s_N=N^2/(4pi^2)`, and canonical total `floor(N^2/3)`. This is a
declared graph/preparation family, **not the existing spherical federation**.
The unchanged native binary rule gives every transition. No oscillator,
Hamiltonian or different update law is inserted. The field observer reads
normalized spatial Fourier sums from its disclosed occupation/current data.

For any fixed finite set of distinct nonzero integer wavevectors, retain
normalized real sine and cosine modes with no duplicate signs or lattice
aliasing. Their Laplacian energies are

```text
kappa_N(k) = (N/pi)^2 [sin^2(pi k_1/N)+sin^2(pi k_2/N)],
0 <= |k|^2-kappa_N(k) <= pi^2 (k_1^4+k_2^4)/(3N^2).          (6)
```

The canonical susceptibility is
`q_N(N^2-q_N)/[N^2(N^2-1)]`, and is divided out in (5).
For the attempt clock (`m=M=2N^2`), the covariance on this band is diagonal:
`1/kappa_N(k)-pi^2/(2N^4)` for each real mode.

For the clock counting the four seams at the origin, let
`G_N=diag(1/kappa_N(k))` and `J_N=W_N^T L_i W_N`. Its covariance is

```text
C_N,incident = G_N + I/(16s_N) - (J_N G_N+G_N J_N)/16.       (7)
```

This retains finite-mesh cross-mode correlations. For the sum `S` of
`|k|^2` over real modes, their normalized endpoint differences give
`||J_N|| <= 16pi^2 S/N^4`. Equations (6) and (7) bound the clock correction by
`pi^2/(4N^2) + 2pi^2 S/(kappa_min N^4)`. Hence **both actual clocks converge
to covariance `diag(1/|k|^2)` on every fixed band**, with an `O(N^-2)` bound.

The execution retains wavevectors `(1,0),(0,1),(1,1),(2,0)` and their eight
real modes. It uses `pi^2<10` to store rational error bounds:

| N | Incident covariance error bound | Computed error diagnostic |
| --- | ---: | ---: |
| 8 | 0.187270 | 0.096215 |
| 16 | 0.028450 | 0.022855 |
| 32 | 0.006046 | 0.005647 |
| 64 | 0.001446 | 0.001408 |

The inequalities certify the exact mathematical matrices. Floating band
evaluations check them independently; their residuals are diagnostics.
The theorem applies to arbitrary fixed finite bands. Their limiting Gaussian
laws are compatible under deleting modes, and define the zero-mode-removed
massless Gaussian generalized field on the declared torus. This identifies
its finite-dimensional laws, not uniform convergence in a spatial path-space
topology. The target is a random distribution: for every `epsilon>0`, its
expected squared `H^(-epsilon)` norm is proportional to the convergent sum
`sum_(k!=0) (1+|k|^2)^(-epsilon)/|k|^2`. This target construction does not
prove tightness of the native field sequence. The formal action is the Dirichlet energy
`(1/2) integral |grad phi|^2` in the declared normalized spatial measure.
For the relation between inverse-Laplacian covariance, generalized fields and
two-dimensional conformal Dirichlet energy, see
[Sheffield, Sections 2.1–2.4](https://arxiv.org/pdf/math/0312099).
This is a Euclidean free-field connection under the supplied graph/chart.
No interacting CFT, stress tensor or Lorentzian continuation is certified here.

The same native mean density has the dissipative multiplier
`1-lambda/(2M)` per attempt. Its long-wavelength energy ratio at momenta
`2k` and `k` tends to four. A positive deterministic rescaling of this clock
does not turn decay into a propagating wave. The Gaussian spatial action
therefore does not identify the original repair dynamics with relativistic
field evolution. Other observables or selected coherent continuations need
their own dynamical correspondence.

## Requirements for the actual source screens

The [coherent screen result](../coherent_screen_observer/README.md) compares
complete local quantum outputs on two matched finite neighborhoods. The
[compiled clock history](../compiled_clock_history/README.md) implements the
finite local preparation/control/readout contract. Both are compatible with
the conditional interpretation here; neither supplies the spatial maps and
Green convergence in (5). A fixed-radius local history bound does not control
the inverse-Laplacian contribution from arbitrarily slow global modes.

The ancestry dictionary supplies certified finite covariances and an
enlarging observation algebra. The source subdivision and native Schur
comparison give an all-level uniform bound, vanishing normalized clock
corrections and a nonzero compatible Gaussian law along a common spatial
subsequence. A quantitative whole-sequence tail, a unique target Green
operator and control of preparation/readout errors under that limit require
further proof. Physical dimensionality, Lorentzian propagation, interacting
continuum fields, and identification with observed physics are separate
obligations. No federation axiom is needed for the finite single-carrier result.
