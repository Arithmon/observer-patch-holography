# Recovering a finite field description from observer records

The construction identifies a supplied finite canonical free field using
executed clock and field readouts. Its inverse receives probabilities and
intervention settings, and recovers a quadratic generator, a three-dimensional
cube geometry, and predictions for withheld readouts. A separate exact inverse
recovers the native twelve-port response geometry from committed currents.
The two geometries are separately identified; no identification of the native
carrier with the eight field sites is assumed or established.

## Relation to existing constructions

The nearest existing results are the [common-source calibration](../common_source/CALIBRATION.md),
its [native response geometry](../common_source/SLOW_GEOMETRY.md), the finite
[spatial attachment](../common_source/spatial.py), and the
[canonical scalar regulator](../relativistic_scalar/DERIVATION.md). Those
constructions already supply field execution, intrinsic response geometry,
or calibration comparisons. The addition here is an explicit boundary
between execution and inference: hidden coordinates and the backend matrix
are unavailable to the inverse, and predictions are made before their
withheld field circuits execute. The present free normal-mode cutoff is a
different regulator from the earlier local-site occupation cutoff. Its
vacuum, pointer and occupation proof are specified below; earlier
interacting-history certificates are not transferred to it.

The harmonic Weyl evolution used here is standard. Nachtergaele, Raz,
Schlein and Sims derive its normal-mode form in
[Section 3.1 of *Lieb-Robinson Bounds for Harmonic and Anharmonic Lattice Systems*](https://arxiv.org/abs/0712.3820).
Their Hamiltonian has kinetic coefficient one and mode evolution
`exp(-2 i gamma(k) t)`. Our kinetic coefficient is one half and the
corresponding evolution is `exp(-i omega_j t)`. Their periodic lattice and
distance hypotheses also differ from this finite cube. We use the harmonic
evolution identity, derived directly below; we do not import their
Lieb-Robinson velocity or claim a strict light cone for this execution.

## Finite execution contract

In units with `hbar=1`, the unrestricted comparison is

```text
H = (p^T p + q^T K q)/2,
[q_x,p_y] = i delta_xy,
K = m^2 I + kappa L_cube,
m = 5/4,   kappa = 3/4.
```

The cube has eight sites and twelve unit-weight edges. A fixed permutation
hides its coordinate labels from the inverse. The real Walsh basis `V`
diagonalizes `K`; every entry has magnitude `1/sqrt(8)` and

```text
omega_j^2 = (25 + 24 popcount(j))/16,    j=0,...,7.
```

The field Hilbert space is the tensor product of eight oscillator spaces.
The supplied normal-mode vacuum is the ground state of this free action.
Its gap is `min omega_j=m` in the specified finite volume. For executed
occupation dimension `d=8` or `d=16`, each mode retains levels `0,...,d-1`.
Writing `P_d` for that projection, its coordinate gate uses
`Q_j,d=P_d Q_j P_d`, where

```text
Q_j = (a_j+a_j^dagger)/sqrt(2 omega_j),
q_x = sum_j V_xj Q_j.
```

The finite free evolution is diagonal, `exp(-i t omega_j N_j)`; the omitted
vacuum energy contributes a common phase. The canonical commutation relation
is exact in the unrestricted comparison; finite occupation matrices have
boundary corrections. The certificate below controls the stated vacuum
histories rather than arbitrary states near that boundary. A field experiment
prepares all eight vacua, executes `exp(i gamma q_x)`, evolves for model interval `t`, and
applies `exp(i eta q_y)` on one arm of a common Ramsey pointer. The same
pointer visits every mode. Its control basis stays fixed until terminal
readout. Separate fresh preparations give the pointer laws at phases
`0` and `pi/2`; both binary outcomes are retained in each experiment.

[source.py](source.py) compiles and executes the local gates through the
existing source compiler. Independent matrix exponentials check the full
phase of each executed matrix; they do not supply its probabilities. The
numerical contraction uses every block of each executed controlled matrix,
including off-block roundoff, while tracing out visited modes. This is the
tensor contraction of the common-pointer circuit. Cached gates reduce host
work; the resource ledger still counts the requested source words for each
fresh trial. Ideal preparation, coherent services and readout remain supplied
apparatus. Numerical stabilization diagnostics are separate from the
mathematical occupation certificate.

## What crosses the observer boundary

[recovery.py](recovery.py) imports neither the source executor nor its target
matrix, coordinate permutation or eigenvectors. Its inputs are field and
clock probabilities, addressed preparation/readout labels, kick strengths,
and an explicit frequency alias bound. Both complete binary outcomes remain
in the receipt; the inverse needs only bit-zero marginals because the
complement is determined by normalization. These are floating evaluations of
ensemble laws, not counts from a finite experiment. No confidence interval
for sampled observations is claimed. Inverse arithmetic runs as a classical
Python calculation on accessible records, not as a source-compiled native
hardware program.

The separate two-level clock uses a supplied positive unit gap and orientation
`exp(+i t N)`. Its readouts are

```text
p0(0) = (1+cos(t))/2,
p0(pi/2) = (1-sin(t))/2.
```

Thus its complex characteristic is `exp(i t)`. The inverse reads its argument
on the declared window `0<t<pi`. The backend supplies model intervals
`tau=0.24` and `2 tau=0.48`, but their numerical labels are not passed to the
inverse: it receives the clock probabilities. The recovered parameter is in
this clock's gap units, rather than seconds, host execution time or a universal
physical clock. `t` and `t+2 pi` give identical clock records; the receipt
executes that alias control. The unwrapped window is a protocol premise.

The native inverse has a separate boundary. For each preparation `x=2 e_j`
it receives only the labelled seam `(a,b)`, committed coin, current `r` and
exact rational probability. The dictionary `delta x=r(e_b-e_a)` gives

```text
T_ij = delta_ij + (1/2) sum_records_given_j probability * r * (e_b-e_a)_i,
L = 60 (I-T).
```

It checks complete mass, fair retained coins and a common scheduler before
forming this operator. Under the declared uniform scheduler, exact arithmetic
in `Q(sqrt(5))` verifies the spectral polynomial with roots
`0, 5-sqrt(5), 6, 5+sqrt(5)`. The projector `P` onto `5-sqrt(5)` has rank three
and diagonal `1/4`; `G=4P` is the recovered unit-diagonal Gram matrix. Output
loads, the source incidence table and source coordinates do not enter this
inverse. Changed-scheduler and relabeling controls test the dependence on
those supplied interface choices.

## Recovering the quadratic generator and spatial relations

Put `Omega=sqrt(K)`. Hamilton's equations, or the canonical Heisenberg
equations, give

```text
q(t) = cos(t Omega) q(0) + R(t) p(0),
R(t) = Omega^(-1) sin(t Omega).
```

A kick `exp(i gamma q_x)` shifts the mean momentum by `gamma e_x`.
For the supplied vacuum, with `C=Omega^(-1)/2`, the ideal Ramsey
characteristic at read site `y` is therefore

```text
chi_yx(t;gamma,eta) = exp(-eta^2 C_yy/2 + i eta gamma R_yx(t)),
p0(theta) = [1 + Re(exp(i theta) chi_yx)]/2,
chi = (2 p0(0)-1) + i (1-2 p0(pi/2)).
```

Opposite kicks recover `R_yx(t)` by subtracting their characteristic phases
and dividing by `2 gamma eta`. Visibility and phase-window checks reject
ill-resolved or observed branch-cut cases. As with the clock, those checks
cannot detect an extra unseen winding without the declared phase premise.

For two calibrated intervals with ratio two,

```text
R(2 tau) = 2 cos(tau Omega) R(tau),
B = (1/2) R(2 tau) R(tau)^(-1),
K = [arccos(B)/tau]^2.
```

For this principal-arccos algorithm, the positive-frequency hypothesis and
supplied bound `tau*omega_max < pi` choose its branch and make `R(tau)`
invertible. This band is an algorithm contract, not a necessary condition
for every possible inverse: exact canonical amplitudes and invertible
`R(tau)` also give `K=(I-B^2) R(tau)^(-2)`. The implementation uses the
arccos route with `omega_max=3`, checks the
response condition number and reciprocity, checks that the ratio has a
symmetric spectrum strictly inside `(-1,1)`, and tests the reconstructed
response amplitude. It does not clip an inadmissible spectrum onto a
physical model. These are numerical consistency checks, not propagated
outward error intervals for `K`.

The geometry fit tests whether this recovered `K` has a common positive
massive cube-Laplacian form. Negative off-diagonal entries supply candidate
edges; reciprocal degree three, common edge weights, row sums and distance
relations must pass. For a unit-edge cube the graph distance `D_xy` equals
the squared Euclidean distance. Consequently

```text
J = I - 11^T/8,
G = -(1/2) J D J
```

must be positive semidefinite of rank three. Its factorization gives
coordinates up to a Euclidean isometry. The edge length is a declared unit.
The recovered `sqrt(kappa)` is a coefficient in those edge and clock units;
it is not an independently measured causal speed. Geometry is identified
within this tested finite cube hypothesis, with the supplied channel and
intervention dictionary.

## Quantum characterization and withheld observations

The free vacuum model gives the conditional prediction

```text
W_xy(t) = <q_x(t) q_y(0)> = [exp(-i t Omega)/(2 Omega)]_xy,
Im W(t) = -R(t)/2.
```

The receipt calculates this two-point function after recovery. It is a
prediction under the supplied canonical commutation relations and vacuum
state, rather than a directly measured continuum Wightman distribution.
The field characteristic above predicts both Ramsey outcomes at the
predeclared withheld intervals `0.31` and `0.61`, four addressed site pairs,
and kick/probe pairs `(0.3,0.4)` and `(0.2,0.3)`. Each prediction is computed
before its field circuit executes. Only after fitting and these predictions
does the audit open the hidden backend stiffness matrix.

Linear response does not establish quantum necessity: a classical harmonic
system with the same `K` obeys the same response equation. A classical Gaussian
phase-space distribution with covariances `C_qq=Omega^(-1)/2`,
`C_pp=Omega/2`, `C_qp=0` also reproduces the symmetrized vacuum covariance
and single Weyl characteristic. The added ordered-loop control checks an
operation beyond these shared predictions. With `U=exp(i a Q)` and
`V=exp(i b P)`, chronological execution of `U,V,U^-1,V^-1` yields

```text
V^-1 U^-1 V U = exp(i a b) I,
```

while reversal yields `exp(-i a b) I`. The source implements the momentum
kicks using quarter-period number rotations of coordinate kicks. The loop
acts on one pointer arm, so its relative phase is observable; a global field
phase without that control would not be. At `a=3/10`, `b=2/5`, the target
phase is `3/25`. Commuting `Q,P` instead predict the identity loop. This
rejects that specified commuting control, without certifying an untrusted
device as quantum or excluding all classical simulations with internal state.

A further control adds `0.25 Q_j^4` to the separate normal-mode generators.
It changes withheld records relative to the recovered free prediction. This
is a deliberate change of action; it is not a local site `phi^4` interaction,
and the free coherent-state cutoff certificate does not apply to it.

## A complete-output occupation bound

The following argument concerns ideal compressed gates against unrestricted
oscillators. Let one mode start in vacuum and undergo displacements `D(delta)`
and number rotations. In every fixed pointer arm suppose the sum of all
absolute displacement amplitudes is at most `A`, with `A^2<d+1`. Number
rotations preserve occupation and the modulus of the coherent amplitude.
Every intermediate unrestricted state is a coherent state up to phase, with
amplitude of modulus at most `A`.

During a unit-parameter displacement segment, the anti-Hermitian generator is
`G=delta a^dagger-conj(delta) a`. For the embedded finite state `psi_d` and
unrestricted state `psi`, the projected difference
`e=psi_d-P_d psi` has a norm-preserving homogeneous evolution plus residual
`-P_d G (I-P_d) psi`. Only lowering occupation `d` crosses this boundary, so

```text
||P_d G (I-P_d) psi||
    = |delta| sqrt(d) |<d|psi>|
    <= |delta| sqrt(d) A^d/sqrt(d!).
```

Integrating each displacement segment and summing its absolute amplitude
gives `||e||<=L`, with

```text
L^2 = d A^(2d+2)/d!.
```

There is no residual on number-rotation segments. For the terminal coherent
state the high-occupation tail obeys

```text
||(I-P_d) psi||^2
  <= sum_(n>=d) A^(2n)/n!
  <= T = min(1, A^(2d)/[d! (1-A^2/(d+1))]).
```

The second inequality uses the ratio bound `A^2/(d+1)` for successive
terms; dropping the coherent Gaussian factor only enlarges the result.
Projected and tail errors are orthogonal, hence the full single-mode vector
error is at most `sqrt(L^2+T)`. Because the controlled isometries are diagonal
in the fixed pointer basis, this estimate also holds when that pointer is
entangled with previously visited modes. Telescoping over `M` modes gives

```text
joint vector error <= M sqrt(L^2+T),
complete output trace distance <= min(1, M sqrt(L^2+T)).
```

The trace distance uses the convention `||rho-sigma||_1/2`. Terminal pointer
readout and retention of its classical record are completely positive
trace-preserving maps, so the same bound covers both outcome probabilities
and the complete output with its remaining fields. A pointer spin mean has
error at most twice the probability bound. No branch is postselected or
renormalized.

For every main calibration and withheld read, `|V_xj|=1/sqrt(8)` and
`omega_j>=5/4`. Therefore

```text
A_j <= (|gamma|+|eta|)/sqrt(16 omega_j),
A_j^2 <= (7/10)^2/20 = 49/2000.
```

[bounds.py](bounds.py) evaluates `L^2+T` with exact fractions and outward
square roots, then sums the eight mode budgets. The resulting probability
bounds are below `1.571645839e-8` at `d=8` and `2.680115492e-19` at `d=16`.
They compare the ideal declared rational settings; they do not certify
binary floating coefficients, numerical matrix contractions or hardware.

For the one-mode Weyl loop, the two coordinate kicks and two momentum kicks
instead have total displacement budget

```text
A = 2 (|a|+omega |b|)/sqrt(2 omega),
omega=5/4,  a=3/10,  b=2/5,
A^2 = 128/125.
```

The same proof applies through the quarter-period rotations and gives a
`d=16` probability bound below `1.103986796e-6`. The occupation estimates are
forward output bounds. Phase extraction, response inversion and spectral
geometry fitting can amplify perturbations; their certified error propagation
and finite-shot confidence regions are separate obligations.

## Criteria for the wider recovery

Full time, space and quantum-field recovery requires connecting the native
observable dynamics to a selected field action and its apparatus, then
controlling observation laws over common enlarging spatial, energy and clock
domains. An interacting relativistic claim additionally needs matched
renormalized parameters, a nontrivial interaction limit, locality and
Lorentz-covariant correlation or scattering predictions, with any mass-gap
statement proved in the same limiting representation. An observer account
must propagate finite-record, regulator and execution errors through its
inference and compare the resulting predictions with independent physical
data. The present contribution supplies the finite record-based inverse and
withheld tests on which those extensions can build.
