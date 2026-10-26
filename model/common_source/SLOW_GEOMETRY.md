# The existing three-dimensional response carrier on actual observer records

[slow_geometry.py](slow_geometry.py) attaches the existing native twelve-port
response geometry to the complete native/quantum observer history. It reuses
the pinned exact [carrier](../../evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py)
and the actual integer primitive in [native.py](../observer/native.py).
No spatial interpretation of the two-oscillator interaction graph is added.

## One carrier, two different protocols

On a seam `(a,b)`, the primitive replaces the endpoint loads by the floor and
ceiling of half their sum. Its fair placement coin averages to the endpoint
mean. Choosing uniformly among all 30 seams therefore gives the exact
one-attempt mean

\[
 T=I-L/60,\qquad L=5I-A.
\]

The executable checks this identity on all 78 total-two states, all 30 seams
and both coins: 4,680 actual primitive branches. This is a calibration
protocol. The original observer instead uses `E=(0,1)` first, followed by
`E` or `F=(0,5)` according to its committed quantum bit. That controller does
not have uniform mean `T`.

For the supplied icosahedral incidence, the nonconstant Laplacian bands are
`5-sqrt(5)` (multiplicity 3), `6` (5), and `5+sqrt(5)` (3). Put
`C_n=Q T^(2n) Q`, where `Q` removes the constant mode. The normalized response
`12 C_n / trace(C_n)` tends to `G=4 P_slow`, with `G` having unit diagonal;
the unnormalized `C_n` itself tends to zero. At every finite attempt count all eleven nonconstant modes
remain present: the rank-three statement concerns the limiting response band.
The limit follows because the other two positive response eigenvalues have
strictly smaller magnitude than `1-(5-sqrt(5))/60`.

The program checks the exact projector identities over `Q(sqrt(5))`, its trace
three, its slow eigenvalue and its coordinate Gram. It also checks all 144
entries against the registered Lean incidence after the exact relabel
`[8,10,1,3,7,11,0,4,6,9,2,5]`. The coordinate witness is the existing
[`rawAxis_dot_eq_sourceGram`](../../Lean/Screen/PrimitivePortFrameQuotient.lean)
construction; the registered identities are in
[PortFrameGram](../../Lean/Screen/PortFrameGram.lean),
[PortGramRepairBand](../../Lean/Screen/PortGramRepairBand.lean) and
[PortGramRepairCovariance](../../Lean/Screen/PortGramRepairCovariance.lean).
The covariance file uses its documented native-decision checks for large
finite tables; this adapter does not add a new Lean proof.

## Relative displacement from committed records

Let `u_p` be the normalized source frame directions. The supplied response
readout of a native load is `x(n)=sum_p n_p u_p`. Since the total load is two,
its centroid is `x(n)/2`. If a primitive records signed current
`j=n_a-n'_a`, conservation gives exactly

\[
 x(n')-x(n)=j(u_b-u_a).
\]

The observer can compute this **relative** displacement from the admitted
seam and committed current. It need not read hidden loads. Absolute load
coordinates appear only as an external implementation audit. They do not
determine the complete state: six distinct antipodal load pairs have the
same zero coordinate. No records or quantum outputs are discarded into
position.

Every one of the 1,056 original joint histories retains its actual native
states, both native records, quantum records, branch-vector reference and
probability `||v_(a,b|j)||^2/264`. Both recorded displacements are checked
against the executed load change exactly. The original `E/F` controller's
active displacement span is only two-dimensional, even though its response
carrier has dimension three.

## A binary controller that also exercises the third direction

The separate `three_direction_continuation` protocol keeps the same first
`E` event and the same quantum histories. Only its second native actuator
changes: final bit zero selects `F=(0,5)`, and bit one selects `G=(0,7)`.
All 1,056 second events are executed again on their actual after-first loads.
Their new currents, final loads and displacements are retained; their earlier
records, quantum outputs and probabilities are unchanged. Quantum references
address the complete vectors retained under `joint_readouts.quantum_histories`.

In the unnormalized source frame the three displacement directions are

\[
 (2,0,0),\quad(1,1-\varphi,\varphi),\quad
 (1,1-\varphi,-\varphi),\qquad \varphi=(1+\sqrt5)/2.
\]

Their determinant is exactly four. The execution checks that all three occur
with positive probability; tests also check the span of the actual total
history displacements. This closes the missing finite controlled direction
with a binary decision. The protocol remains different from the uniform
30-seam calibration and does not establish isotropic or autonomous motion.

## The native directions suffice for finite null-tensor tomography

The existing [Einstein tensor theorem](../../Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean)
uses nine supplied null directions to determine a symmetric four-tensor
up to a metric multiple. [null_tomography.py](null_tomography.py) gives an exact
executable replacement for that directional-completeness premise using the
same twelve source ports. This is an algebraic certificate over `Q(sqrt(5))`,
separate from the registered Lean proof.

Write the raw port axes as `v_p`, with `r²=phi+2`. On the **declared Lorentz
lift** `k_p=(1,v_p)`, the metric is `g=diag(-r²,1,1,1)`. The coordinate change
`t'=rt`, followed by rescaling each null ray by `1/r`, gives `(1,v_p/r)`
for `diag(-1,1,1,1)`; the corresponding quadratic charges scale by `1/r²`.
The twelve rows `q_p=B(k_p,k_p)` have rank nine and kernel exactly `span(g)`.
For the representative with `B00=0`, put `b_i=B0i`, `S_ij=Bij` and let
angle brackets denote the arithmetic mean over the twelve ports. The inverse is

\[
 b_i=\frac{3}{2r^2}\langle v_i q\rangle,\qquad
 S_{ij}=\frac{15}{2r^4}\langle v_i v_j q\rangle
        -\frac{3\delta_{ij}}{2r^2}\langle q\rangle.
\]

The exact matrix product of this decoder with the nine-column `B00=0`
design is the identity. Its product with the full ten-column design gives
`B+(B00/r²)g`, proving both rank nine and the full metric kernel without
floating rank thresholds. Equal charges for two symmetric tensors on these
twelve directions therefore imply metric ambiguity on every direction.

The compatible charge space is exactly the native constant band, the slow
rank-three band and the antipodal-even rank-five band. Three constraints
exclude the fast rank-three band with Laplacian eigenvalue `5+sqrt(5)`.
`decode_null_charges` checks all twelve equations and rejects incompatible
readings; a single-port reading is an explicit rejection control. Combining
opposite directions lowers the rank to six and erases all three mixed
time/space components. Direction labels must retain their orientation.

The exact decoder also bounds measurement error: a charge perturbation of
maximum magnitude `epsilon` changes every component of its canonical
tensor estimate by at most `(1-sqrt(5)/5)*epsilon`. Its compatibility
residual has maximum error at most `phi*epsilon`. These are exact induced
matrix norms; noisy incompatible readings remain rejected by the exact
decoder. No observational error distribution is assumed.

This closes the finite spanning requirement **if** a common symmetric
stress/geometry charge is measured on the declared lifted directions.
The source Gram does not construct that detector, select the Lorentz lift,
or establish small-ball geometry, entropy stationarity, Ward/Bianchi laws
or the physical Einstein equation. The clock and charge dictionary remain
the physical composition requirements.

## Tensor charges from actual local records

`native_record_tensor_certificate()` attaches the same decoder to measured
occupancy bits. Each port receives its initial bit and retains its own
incident signed-current records. It reconstructs every later bit by
`n_i=n_i(0)+sum j*(1_{i=b}-1_{i=a})`; no later hidden load is an input.
All 3,960 branches of the 66-state binary sector verify this readout against
the actual primitive. A complete two-event test retains all 3,600
seam/coin/current histories and reproduces their joint endpoint law.
Initial-bit disclosure, complete local tapes and collection at the common
attempt horizon are declared access conditions.

For the uniform all-seam, fair-coin protocol the twelve ensemble means obey
`m_n=(I-L/60)^n m_0`. This source-derived mean law gives the exact decomposition

\[
 m_n=D\widehat B_n+r_n,\qquad
 r_n=\left(\frac{55-\sqrt5}{60}\right)^n P_{\rm fast}m_0,
 \qquad
 \|r_n\|_\infty\le
 \left(\frac14+\frac{\sqrt5}{20}\right)
 \left(\frac{55-\sqrt5}{60}\right)^n.
\]

Here `D` is the same null-charge design, `Bhat` its canonical tensor
estimate and `P_fast` the native fast-band projector. The sharp initial
bound holds for every binary preparation and every mixture, and is below
`1e-4` after 64 attempts. Exact decoding continues to reject nonzero
residuals. The approximation has an explicit source error; it does not
replace incompatible data by an asserted exact tensor.

A disclosed preparation occupying one antipodal pair gives nonconstant
compatible means at every horizon. The nine ports `0,1,2,4,5,6,8,9,10`
then predict the three withheld means `3,7,11` exactly, with an exact inverse
and no fit to those readings. An adjacent-pair preparation fails this test
at every finite tested horizon, with the residual decay above. For rational
preparations, conjugating `sqrt(5)` exchanges the slow and fast projectors;
zero fast component therefore forces an antipodal-even mean. The six
antipodal binary preparations are exactly the compatible deterministic
ones. Selecting such a preparation is an experimental input, not a unique
state selected by repair. The mean statement is also distinct from a
history statement: after one antipodal trial, the mean passes but individual
readings fail exact compatibility with probability `1/6`.

`interval_tensor_estimate` also accepts twelve unconditional occupation
intervals from complete record laws with unresolved outcomes. Their
midpoints define an estimate; maximum halfwidth `rho` gives coefficient
error `(1-sqrt(5)/5)*rho` and residual error `phi*rho`. It reports the
midpoint's exact compatibility separately and can exclude a declared
native-attempt horizon when every allowed residual exceeds its source
bound. Passing that necessary test establishes neither interval feasibility
nor exact compatibility of the unknown mean. Conditioning away unresolved
outcomes is outside this input contract; sampling uncertainty is separate.

For `N` independent complete trials, binary-readout variance and a union
bound give probability at most `min(1,3/(N*eta²))` that any sampled mean
misses by `eta`; correlations between ports are allowed. Serial samples
from one trajectory do not meet that independence contract. On the
complementary event, canonical tensor error is at most
`(1-sqrt(5)/5)*eta`, and compatibility residual is at most the source bound
above plus `phi*eta`. The withheld-direction errors are bounded by
`(2+2*sqrt(5))*(source_bound+eta)`. These are finite occupancy statistics on the declared
Lorentz lift. No stress-energy normalization, Ward identity, physical
boost, elapsed clock duration or continuum field is inferred.

## What this establishes about space and momentum

This is a concrete common-carrier attachment of native records to an existing
rank-three response quotient. The incidence, calibration schedule and use
of its Gram as an operational readout remain declared construction inputs.
The Cartesian frame is a coordinate witness, determined only up to isometry.
No physical length, `P`-dependent scale, `N`-dependent capacity or calibrated
clock-to-displacement relation is inferred here.

The quotient/completion theorem provides an abstract three-dimensional
Euclidean carrier with its declared topology. Identifying that carrier with
measured position still requires a detector dictionary, translation dynamics
and a refinement/overlap attachment. Momentum additionally requires the
appropriate translation representation and its generators. An exact
canonical relation `[X,P]=iI` cannot hold on this finite load space, since a
commutator has trace zero. The finite observer, its readable modular clock and
this carrier can be joined without asserting that those remaining physical
identifications are already supplied.

Run `python -m pytest -q model/common_source/test_slow_geometry.py`. The common
source receipt retains `report()` together with the rest of the experiment.

<a id="modular-boosts-require-a-controlled-continuum-comparison"></a>
## Modular boosts require a controlled continuum comparison

[modular_covariance.py](modular_covariance.py) supplies a finite test of the
geometric modular-flow premise. Let $G=\log\sigma$ on a faithful finite
observer algebra, so the Heisenberg modular flow is
$U(s)=\sigma^{is}=e^{isG}$, with dimensionless modular parameter $s$.
A boost acting on a null-translation generator would require, after choosing
its orientation and normalization, $U(s)PU(s)^*=e^{\kappa s}P$ for some
nonzero real $\kappa$. No physical null-translation generator is assigned
to the native code by this test.

For every finite Hermitian $G,P$, cyclicity of the trace gives
$\mathrm{Tr}(P\,i[G,P])=0$. Consequently

```math
\|i[G,P]-\kappa P\|_{\rm HS}^2
=\|[G,P]\|_{\rm HS}^2+\kappa^2\|P\|_{\rm HS}^2.
```

The full-algebra relative Ward defect is at least one whenever
$\kappa P\ne0$. This is an exact identity; the executable checks its
attachment to the native reference without using a numerical rank decision.
There is also a finite-parameter operator bound,

```math
\|U(s)PU(s)^*-e^{\kappa s}P\|
\ge |e^{\kappa s}-1|\,\|P\|,
```

from unitary invariance of the operator norm and the reverse triangle
inequality. Positivity of $P$ is unnecessary. The operator bound applies
to bounded operators more generally; unbounded continuum translation
generators require their own domain and convergence arguments. Thus it
does not exclude a physical continuum limit of finite observers.
Finite Lorentz matrices acting on Weyl field indices are a different
representation; this bound concerns unitary modular conjugation on the
observer's full finite Hilbert algebra.

Compression changes the comparison, but its boundary term must be included.
For an orthogonal projection $Q$, put $G_Q=QGQ$, $P_Q=QPQ$ and

```math
B_Q=i\bigl(QG(1-Q)PQ-QP(1-Q)GQ\bigr).
```

Direct multiplication gives

```math
Q\,i[G,P]Q=i[G_Q,P_Q]+B_Q,
\qquad
\|Q(i[G,P]-\kappa P)Q\|_{\rm HS}+\|B_Q\|_{\rm HS}
\ge |\kappa|\,\|P_Q\|_{\rm HS}.
```

A subspace invariant under $G$ or $P$ has $B_Q=0$ and inherits the same
finite-dimensional obstruction. In particular, selecting a spectral band
of the modular generator alone cannot close the geometric Ward premise.
For a useful refinement, the permitted states, observable domain and
finite modular-parameter window must be specified, and their actual
compressed response must converge with the boundary contribution retained.
Requiring both a vanishing compressed defect and a vanishing boundary term
while keeping a nonzero normalized $P_Q$ contradicts the displayed bound.

The source example uses the same faithful 66-state reference as the particle
response. Choose native configurations $a,b$ with protected-core occupations
zero and two. Their weight ratio is five. On their two-state corner choose
$P=\begin{pmatrix}1&i/2\\-i/2&1\end{pmatrix}$, zero outside, and let
$Q$ project onto $(|a\rangle+|b\rangle)/\sqrt2$.
These are declared diagnostic probe choices. $P$ is positive, with nonzero
eigenvalues $1/2,3/2$. For $\kappa=\log(5)/2$, its compressed derivative
satisfies $Q\,i[G,P]Q=\kappa P_Q$ exactly. The internal compressed
commutator is zero; the whole derivative comes from $B_Q$.
Yet its finite-parameter probe response is

```math
\langle U(s)PU(s)^*\rangle_Q
=1+\tfrac12\sin(s\log 5),
```

which differs from the dilation target $e^{s\log(5)/2}$. At
$s=1/\log5$ the values are approximately $1.420735$ and $1.648721$.
A matched derivative at one parameter value therefore does not establish
geometric modular covariance. The computation is a matrix diagnostic on
the ideal reference, not a newly executed clock or momentum instrument.

The entropy premise has a different status. The existing
[native entropy refinement](../native_equilibration/entropy_refinement.py)
derives all-partition counting-gradient error at most $3\ell^5/32$,
quantum-gradient error at most $\ell^5/16$, and entropy deficit at most
$\ell^{10}/512$, after its dimension-dependent discarded warm-up schedule.
Here $\ell=2^{-l}$ is a declared refinement coordinate, and the permitted
tangent norms and source sector are specified. The
[common variation theorem](../native_equilibration/common_variation.py)
also supplies constrained MaxEnt stationarity for the nonuniform protected
core reference. These results remove exact entropy preparation as a finite
premise. They do not identify $\ell$ with a proper small-ball radius,
fixed native capacity with proper volume, entropy variation with
$\delta A/(4G)$, or the modular generator with the same physical stress
that occurs in the Einstein theorem. Those identifications and the
restricted continuum covariance comparison remain required on one family.

Run `python -m model.common_source.modular_covariance` for the native
diagnostic and `python -m pytest -q model/common_source/test_modular_covariance.py`
for the sharp, spectral-compression and finite-parameter controls.

## Entropy stationarity from continuing native records

[history_entropy.py](../native_equilibration/history_entropy.py) extends the
finite entropy envelope to a continuing native population history with
record-dependent seam choices, core-crossing currents and bounded stopping.
The experiment starts from the same faithful 66-state reference
$\sigma(x)=w_{q(x)}$, with $q=x_0+x_1$ and $c=q(q-1)/2$.
More generally, the proof holds for every faithful two-moment family in
[common_variation.py](../native_equilibration/common_variation.py).
The reference is prepared before the experiment; the continuing history
does not reset or re-equilibrate between its measurements.

The declared controller chooses seams from prior committed records. Every
fair coin and current (C law), or pre-update contrast (E law), remains in
the tape. If $j$ is the signed transfer on $(a,b)$, the actual native rule
gives, for every cap $R$,

```math
N_R(x')-N_R(x)=j\bigl(1_{b\in R}-1_{a\in R}\bigr).
```

For E records, the transfer itself is decoded as
$j=[d-(2u-1)(|d|\bmod2)]/2$, where $u$ is the recorded fair coin.
Thus endpoint plus complete tape reconstructs the initial configuration
by reversing these translations. This is an exact finite continuity law,
not an identification with a physical stress tensor.

The core readout uses native records only. After the prefix, a census on
the source seams $(0,1),(0,5)$ determines occupation differences on three
connected ports. At total binary load two, these three ports cannot all
be occupied, so their minimum occupation is zero; this fixes the common
offset. Reversing the prefix currents then determines the initial core
charge. Under E, the two contrasts are recorded in **two native attempts**,
with no failure or postselection. Under C, each edge waits for the first
coin change, with deadline $h$. Its first and final currents determine the
pre-edge contrast. The complete success probability is
$\beta=\alpha^2$, $\alpha=1-2^{1-h}$, independent of the input state.
All timeouts retain their disturbed endpoints and complete tapes.
For a two-attempt prefix and $h=24$, the C service uses at most 50 native
attempts and has timeout probability $1-(1-2^{-23})^2$; its mean is at
most eight attempts. These are native attempt counts, not observer proper
times. The two-attempt E prefix plus census uses four attempts exactly.
The already declared preparation and controller services remain additional
resources; no free occupation measurement is assumed.

For a complete successful tape $b$, let $d_b$ be the **integer number of
compatible native endpoints**, counted before evaluating any entropy.
All members of this fibre have the same reconstructed initial charge and
the same coin word of length $T_b$. Consequently

```math
\pi_\sigma(b)=\frac{d_b w_{q(b)}2^{-T_b}}{\beta}.
```

For E, $\beta=1$ and every terminal tape is included. For C, the displayed
law is the explicitly flagged successful subensemble; its timeout mass is
kept separately. For another initial population $p$, using the same
controller and readout, define

```math
\mathcal F(p)=H(\pi_p)+\sum_b\pi_p(b)\log d_b
-\mathbb E_p[T]\log2-\log\beta.
```

The counted completion and the coin costs are independently defined by
the executed instrument. They are not chosen to satisfy a gravitational
equation. Substituting the reference law gives the exact identity

```math
\mathcal F(p)-S(\sigma)
=h_1\Delta\langle q\rangle+h_2\Delta\langle c\rangle
-D(\pi_p\Vert\pi_\sigma),
\qquad
h_1=\log(w_0/w_1),\quad h_2=\log\frac{w_1^2}{w_0w_2}.
```

At fixed initial moments, this is a maximum and its first variation
vanishes. Its Hessian is $-\sum_b(\delta\pi_b)^2/\pi_\sigma(b)$.
It is semidefinite: an unobserved spectator perturbation is a tested flat
direction, so the record functional does not uniquely determine the full
source state. A record-dependent stopping rule changes the expected coin
cost under a perturbation. The executable checks this explicitly by ending
the prefix at its first nonzero current or its third attempt. Keeping the
old mean cost instead of varying $\mathbb E_p[T]$ breaks the identity.
Subtracting this known coin entropy does not entail physically erasing a
tape, and it does not identify the uncorrected entropy with geometric
generalized entropy.

The finite preparation error also transfers to this instrument. If the
prepared density lies between $r^{-1}\sigma$ and $r\sigma$, the actual
record probabilities obey the same sandwich, including normalized C
success because its probability is state-independent. On every normalized
fixed-moment tangent $X$, the entropy derivative error is at most
$\log(r)\|X\|_1$. With the existing preparation radius $v$,
$r=(1+v)/(1-v)$ gives the retained bound $2v\|X\|_1/(1-v)$.
Feedback and the properly varied stopping cost introduce no new factor.
This is a bound for record entropy; it does not assert exact finite
preparation or a von Neumann entropy theorem for the coherently continuing
work and its untouched quantum reference.

Timeouts show why the charge-resolution condition matters. If a complete C
timeout fibre mixes initial charges, replacing its conditional reference
law by uniform counting gives a positive excess. For the whole terminal
instrument, this excess is

```math
\sum_b\pi_\sigma(b)
D\bigl(\sigma(\cdot\mid b)\Vert\mathrm{Unif}(d_b)\bigr)
\le (1-\beta)\log(w_{\max}/w_{\min}).
```

The local reference has $w_{\max}/w_{\min}=5$. The exact C audit at $h=3$
finds an excess of approximately $0.0149548033$ nats. E's unconditioned
record readout has zero excess. Dropping the informative current/contrast
labels destroys charge resolution, even when the fair-coin word is kept.
These controls rule out treating a coarse or incomplete tape as the same
entropy observable.

Cap agreement uses a single history measure. In a charge-resolved fibre,
the common reference assigns the uniform conditional distribution to its
$d_b$ compatible endpoints. Every cap completion is a marginal of that
same distribution. Nested restrictions and the intersection of overlapping
caps therefore agree exactly, conditioned on the same committed record;
there is no independent state resampling for each observer. Refining a
record by a predicted cap occupation reallocates its bulk and counted
completion entropies while preserving $\mathcal F(\sigma)=S(\sigma)$.
An unmeasured cap occupation is a conditional prediction, not an additional
observed record. The tests compare overlapping caps directly on every
successful fibre.

This discharges a finite operational first-law and constrained-stationarity
interface **with adaptive recorded histories**. The source parameter is
the declared initial-state variation; the constraints are its two
occupation moments. It does not identify them with fixed proper volume,
$\log d_b$ with $A/(4G)$, the current with physical stress, or the stopping
cost subtraction with a gravitational counterterm. In particular, the
Einstein theorem still needs its independent small-ball bulk kernel,
fixed-volume area expansion and common continuum error control. The
destructive census is also a different protocol from the coherent Ramsey
particle experiment: a single quantum matter-plus-gravity history has not
been certified by combining their separate finite results.
The finite proof is analytic with exact rational source checks; it is not
a new Lean formalization of the Einstein stationarity premise.

The exact source audit is retained in
[variation_results.json](../native_equilibration/variation_results.json),
under `continuing_record_entropy`. Run
`python -m pytest -q model/native_equilibration/test_history_entropy.py`.
