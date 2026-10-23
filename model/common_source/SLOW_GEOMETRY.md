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
