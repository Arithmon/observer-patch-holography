# Identifying a local interacting field from executed records

This construction recovers coefficients and spatial relations within a supplied
finite canonical scalar family, then predicts complete coupled measurement
histories. The recovery inputs are addressed preparation and instrument labels,
binary detector laws, and a separately calibrated clock. The hidden action,
coordinate labels and backend matrices do not enter the inverse. The action
family, oscillator basis and ability to isolate onsite and pair generators are
premises of the experiment.

## Relation to existing results

The nearest inverse is the [free-field reconstruction](../observer_reconstruction/DERIVATION.md).
It recovers a quadratic stiffness matrix from executed linear-response records;
its altered normal-mode quartic term is a changed-action control, not a local
site interaction. The [canonical scalar construction](../relativistic_scalar/DERIVATION.md)
already executes a supplied local quartic action and bounds particular
disturbed histories. The addition here is identification of the local quartic
coefficient and all pair weights from separate calibration experiments, followed
by prediction of withheld, fully coupled histories. The two existing
constructions remain unchanged.

Hamiltonian identification from controlled input-output observations is an
established problem. Burgarth and Yuasa's
[*Quantum System Identification*](https://arxiv.org/abs/1104.0583) describes how
recoverable information depends on the experimental access and prior knowledge.
Burgarth and Maruyama's
[*Indirect Hamiltonian Identification through a small gateway*](https://arxiv.org/abs/0903.0612)
addresses identification through restricted access with known coupling topology.
Our protocol supplies isolation access to every onsite generator and every pair
service. Its elementary finite-matrix derivation below uses that stronger
access; it does not import a small-gateway theorem or establish identifiability
from passive observations.

## Supplied family and observable boundary

With unit `hbar`, the comparison action is

```text
H = sum_i [p_i^2/2 + A_i q_i^2/2 + g_i q_i^4]
    - sum_{i<j} kappa_ij q_i q_j,
q = (a+a^dagger)/sqrt(2 Omega),
p = i sqrt(Omega/2) (a^dagger-a).
```

For the homogeneous scalar family,
`A_i=m^2+sum_j kappa_ij`, `g=lambda/(24 a_spacing^3)` and
`kappa=c^2/a_spacing^2` on each edge. The supplied basis frequency `Omega=2`
defines preparation and instrument labels. It is not an inferred excitation
gap. Absolute length and physical seconds are not recovered; edge length and
clock gap fix the displayed units.

For occupation dimension `d`, local polynomials mean `P_d q^2 P_d`,
`P_d p^2 P_d` and `P_d q^4 P_d`, formed before compression. Powers of
`P_d q P_d` would give different boundary matrix elements and a different
action. Calibration uses four local levels for each onsite generator and two
levels per endpoint for an isolated bond. Each pair service is queried,
including the pairs with zero coupling. Isolation retains the onsite `A_i`
term while other services are switched off; changing that convention would
change the inverse formulas.

[source.py](source.py) executes number preparation, isolated evolution and a
complete binary pointer instrument through the existing twelve-port source
compiler. The pointer's bit-one effect is the requested projector; its
complement is retained as bit zero. Independent finite targets audit the
executed matrices and never replace them in detector probabilities. The
inverse receives neither those matrices nor their eigenvectors. It processes
floating evaluations of probability laws, not sampled experimental counts.
Its arithmetic is a Python computation on accessible records; that arithmetic
has not itself been compiled into the source observer.

## An off-diagonal generator element from two binary instruments

Let `H` be the real symmetric isolated finite generator. Prepare `|j>` and
read the two projectors onto

```text
|e_±> = (|j> ± i |k>)/sqrt(2),       j != k.
```

Each sign is a fresh preparation. Write

```text
p_±(t) = |<e_±| exp(-i H t) |j>|^2,
s(t) = [p_-(t)-p_+(t)]/2.
```

The state derivative at zero is `-i H|j>`, giving `s'(0)=H_kj`.
Reality of `H` and the initial state gives
`exp(+i H t)|j> = conjugate(exp(-i H t)|j>)`, which interchanges the two
projector probabilities. Thus `s(-t)=-s(t)` exactly. Subtracting the two
projector laws removes the even population terms, including leakage outside
the span of `|j>,|k>`. No two-state invariant subspace is assumed.

Let `B >= ||H-cI||` for some real scalar `c`. The observable whose expectation
is `s` has norm `1/2`. Repeated commutators and
`||[H-cI,O]|| <= 2B||O||` give

```text
|s^(5)(t)| <= (2B)^5/2,
s(t) = H_kj t + b t^3 + r(t),
|r(t)| <= (2B)^5 |t|^5/240.
```

For two decoded durations `0<u<v`, form

```text
h_hat = [v^2 s(u)/u - u^2 s(v)/v] / (v^2-u^2).
```

This cancels the cubic term. Applying the two remainder bounds separately
gives the conservative finite-duration certificate

```text
|h_hat-H_kj| <= (2B)^5 u^2 v^2 (u^2+v^2) / [240(v^2-u^2)].
```

For `u=t/2`, `v=t`, this becomes `B^5 t^4/18`. The implementation uses the
two decoded durations individually, so floating clock evaluations need not
have an exactly dyadic ratio. The declared caps are `B=16` for onsite
calibration and `B=2` for pair calibration. A later backend audit verifies
these caps for the actual source generators. Their validity is a premise
of the interval calculation, not a result inferred from the four laws.

If each of the four probabilities has absolute error at most `delta`, the
additional coefficient error is bounded by

```text
delta [v^2/u + u^2/v] / (v^2-u^2),
```

or `3 delta/t` for dyadic durations. Shorter probes reduce Taylor error but
amplify probability error. [bounds.py](bounds.py) evaluates these certificates
as exact rationals for the represented inputs and exports outward upper
floats. The inferred centers, clock calibration and source matrix execution
still have separate numerical limitations.

For portable replay, the receipt retains these numerical upper bounds but
omits exact rational expansions that depend on the last bits of decoded
floating clock values or propagated radii. The bound helpers regenerate those
exact expansions from their inputs. Fixed rational occupation and splitting
certificates remain in the receipt and are compared exactly.

For `M` independently estimated Bernoulli laws, a union of Hoeffding bounds
gives simultaneous error at most `delta` with failure probability at most
`alpha` when each law has

```text
N >= log(2M/alpha)/(2 delta^2)
```

fresh trials. The receipt reports this prospective resource count for field
calibration and keeps clock laws exact in that sensitivity calculation. It
does not claim those trials were executed or that systematic errors satisfy
this sampling bound.

## Separating the quartic term from the quadratic term

For the onsite generator `p^2/2+Aq^2/2+gq^4`, oscillator algebra gives

```text
H_(n+2,n)/sqrt((n+1)(n+2))
    = (A-Omega^2)/(4 Omega) + g(2n+3)/(2 Omega^2).
```

The calibration separately prepares `|0>` and `|1>`, recovering `H_20` and
`H_31` from the corresponding imaginary-coherence instruments. Set
`u0=H_20/sqrt(2)` and `u1=H_31/sqrt(6)`. Then

```text
g = Omega^2 (u1-u0),
A = Omega^2 + Omega (10u0-6u1).
```

These low generator matrix elements coincide with those of the unrestricted
oscillator polynomial. The finite-duration calibration certificate applies
to the compressed executed evolution; it does not assert an occupation bound
for the excited calibration preparations. Transfer of the inferred
coefficients to other occupations is conditional on the supplied polynomial
family, including the absence of additional operators invisible to these
probes.

Four local levels are important for this identification. On the first two
levels,

```text
P_2 q^4 P_2 = (3/Omega) P_2 q^2 P_2 - 3 I/(4 Omega^2).
```

Consequently a quartic term there can be absorbed into a changed quadratic
coefficient and an unobservable scalar energy shift. A two-level onsite
experiment cannot distinguish those terms. The `0→2` and `1→3` probes avoid
that particular degeneracy. A positive propagated interval for `g` excludes
quadratic onsite generators within this specified canonical family; it does
not exclude arbitrary alternative hidden-state models.

## Pair weights and recovered spatial relations

An isolated pair generator is `-kappa_ij q_i q_j`. Starting from `|00>`,
the same procedure recovers

```text
H_(11,00) = -kappa_ij/(2 Omega),
kappa_ij = -2 Omega H_(11,00).
```

The two-level codes retain this matrix element exactly. All 28 labelled pairs
are probed. Exact topology recovery uses an explicit separation premise:
every weight is zero or at least the supplied positive threshold. An interval
containing zero but lying below that threshold identifies a nonedge only
within that family. Without this premise an arbitrarily weak positive edge
could not be excluded by a finite-width interval. Ambiguous intervals are
rejected.

After this classification, the inverse tests homogeneous positive edge
weights, common `g`, and common positive bare mass squared
`m_i^2=A_i-sum_j kappa_ij`. It then tests whether the resulting graph is the
eight-vertex cube. For this particular graph, shortest-path distance `D_ij`
is the squared Euclidean distance of a unit-edge cube. Therefore

```text
J = I - 11^T/8,
G = -J D J/2
```

is its centered Gram matrix, positive semidefinite of rank three. Geometry
is determined up to isometry and the declared length unit. This is a tested
finite cube realization of a supplied action family. It does not identify
the field sites with the native twelve-port geometry, derive universal
spatial dimension, or produce a spacetime metric or Einstein equation.

## Clock reading and withheld disturbed histories

The separate phase clock uses the supplied unit-gap rotation `exp(+itN)`.
Its two Ramsey laws determine `cos(t)` and `sin(t)`; the chosen phase window
`0<t<pi` selects a reading. Durations differing by `2pi` cannot be
distinguished by these laws alone. Thus the phase window, energy unit and
clock orientation are apparatus premises. The recovered parameter is a
calibrated Hamiltonian model coordinate, not host execution time. The clock
is not a degree of freedom dynamically maintained by the coupled field.

The recovered action predicts histories in two tensor boxes,
`(8,2,2,2,2,2,2,2)` and `(8,4,4,4,4,4,4,4)`. Predictions use sparse
continuous evolution under that inferred action. The withheld source instead
executes onsite and bond words with symmetric splitting. Both start from the
supplied `Omega=2` product vacuum, which is not the interacting vacuum.
At the addressed eight-level site the binary instrument is

```text
K_1 = |e><e|,       |e> = (|0> + sign*i*sqrt(24)|4>)/5,
K_0 = I-K_1.
```

The first record controls a phase on level four, followed by further coupled
evolution and the same binary instrument. If `F_b` denotes that conditional
feedback, each complete record probability is

```text
p(b,c) = ||K_c U(t) F_b K_b U(t) |vacuum>||^2.
```

All four branches remain unnormalized until their probabilities are formed;
there is no postselection. Both instrument signs are predicted before their
corresponding coupled source executions. The native attachment conditions
the sign and consumes the final record using the existing declared interface.
A common stochastic attachment cannot increase total variation. This does
not make the imposed field action or feedback rule intrinsic to the carrier.

## Separate error budgets for the complete history

Let the inferred generator differ from the true finite compressed generator
by coefficient radii `delta_A_i`, `delta_g_i` and `delta_kappa_ij`. For dimension
`d` and the supplied basis frequency,

```text
||P_d q P_d|| <= sqrt(2(d-1)/Omega),
||P_d q^2 P_d|| <= (2d-1)/Omega,
||P_d q^4 P_d|| <= (2d-1)(2d+1)/Omega^2.
```

The first follows from the truncated raising and lowering operator norms.
For the second, use `q^2 <= 2(N+1/2)/Omega` as a quadratic-form inequality.
For the third, `q` sends the first `d` levels into the first `d+1`; applying
the same inequality twice bounds `||q^2 P_d||^2=||P_d q^4 P_d||`.
The triangle inequality then gives

```text
Delta_H = sum_i [delta_A_i ||P q_i^2 P||/2
                 + delta_g_i ||P q_i^4 P||]
          + sum_(i<j) delta_kappa_ij ||P q_i P|| ||P q_j P||
        >= ||H_inferred-H_true||.
```

Pairs declared absent have zero uncertainty only under the separation
premise above. Otherwise their possible small weights must also enter this
sum. Duhamel's inequality bounds evolution error over one interval by
`t Delta_H`. Retaining every instrument outcome makes the full instrument
an isometry into its quantum-record output space. Telescoping the two
intervals through that isometry and the conditional unitary feedback gives

```text
history_TV <= retained_trace_distance <= min(1,2t Delta_H).
```

The existing symmetric-product certificate supplies an additional bound
for the actual default source action. If `S` bounds the sum of centered
local generator norms and `n` is the step count, its two-segment bound is
`epsilon_split=min(1,2 S^3 t^3/(3n^2))`. It is not reapplied to arbitrary
inferred parameters. This suffices: compare the source product with the
true finite continuous action, then compare that action with the inferred
one.

The existing [exact history certificate](../relativistic_scalar/exact_history.py)
compares this true finite action with unrestricted oscillators on the same
eight spatial sites. Its rational Taylor argument starts from the product
vacuum and retains both measurements and feedback. A cube automorphism sends
the actual addressed vertex to its canonical origin, preserving the
homogeneous action, the unique eight-level read site, the other seven equal
cutoffs, the product preparation, instrument sign, and conditional level-four
phase. Consequently it applies after that relabeling; it is not reused for
the isolated excited calibration experiments.

Writing its occupation error as `epsilon_occ` and the coefficient error
above as `epsilon_param`, the three distinct comparisons are

```text
inferred finite continuous law ↔ executed source law:
    TV <= epsilon_param + epsilon_split;
executed source law ↔ original unrestricted-oscillator law:
    TV <= epsilon_split + epsilon_occ;
inferred finite continuous law ↔ original unrestricted-oscillator law:
    TV <= epsilon_param + epsilon_occ.
```

Each bound is capped at one. The last comparison does not pass through the
split execution, so adding splitting error there is unnecessary. Floating
linear algebra, finite clock accuracy and physical service errors remain
outside all three ideal bounds. Receipt numerical guards are not proofs of
those additional error sources.

The matched quartic-deletion control sets `g=0` while retaining the inferred
quadratic and pair coefficients. Its disagreement tests this specified
alternative, while the separate positive-`g` intervals establish the local
identification within the declared family. The finite history discriminator
must not be promoted to an unrestricted-oscillator discriminator when the
available occupation error exceeds its margin. A common controlled spatial
and energy limit, interacting relativistic covariance, an interacting vacuum
and its observable mass gap, and endogenous selection of the action and
clock remain further construction tasks.
