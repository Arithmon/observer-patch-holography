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
