# Comparing twelve clocks on one carrier

The actual native currents give twelve local clocks a common relative rate and
a spatially structured uncertainty. A federation, independent oscillator,
stationary preparation and seconds-per-attempt input are unnecessary for this
result. It uses the original twelve-port, total-two binary carrier, all thirty
seams, the fair placement coin and independent uniform seam choices. There are
no interventions during the comparison. The C and E quantum completions give
the same classical clock law.

Let `n` count completed **moving currents anywhere on this one carrier**. Let
`K_v(n)` count those currents that change occupation at port `v`. The observer
can read these bits from the source current and seam; it need not inspect the
hidden load configuration. Every movement ticks its two endpoints, so

```text
sum_v K_v(n) = 2n                  on every history,
K_v(n)/n -> 1/6                    in probability and mean square.
```

Every pair of local clocks therefore agrees in asymptotic count rate. This is
a stochastic comparison of readings, with finite uncertainty. It does not
assert that their individual ticks coincide or that equal finite counts occur
at the same event. The record interface collects the specified local bits;
it supplies no physical message delivery or simultaneity convention.

## Exact finite uncertainty without stationary preparation

The binary source has 66 configurations. After removing no-movement attempts,
its transition matrix is `J_ij=A_ij/d_i`: `A` is the source move adjacency,
and `d_i` is 8 or 10. Its invariant event law is `nu_i=d_i/600`.
This invariant law is used to solve a finite equation, not imposed as an
initial preparation. The source graph is connected, and every state has a
positive escape probability. Thus each finite movement index is reached
almost surely under the stated scheduler. No-movement attempts remain actual
records and consume resources; they do not increment `n`.

For `g_v(i,j)=abs(x_j(v)-x_i(v))-1/6`, the executable solves and checks

```text
f_v(i) = sum_j J_ij g_v(i,j),
(I-J) h_v = f_v,     sum_i nu_i h_v(i)=0.
```

The solution is exact rational arithmetic on all 66 states and all 12 ports.
Set `D_v(n)=K_v(n)-n/6`. Then

```text
M_v(n) = D_v(n) + h_v(X_n)-h_v(X_0)
```

is a martingale: each increment has zero conditional mean by the checked
Poisson equation. The same calculation gives uniform constants

```text
R = max(h_v)-min(h_v) = 38/25,
V = max_i E_i[(delta M_v)^2] = 30341677/47045000,
abs(delta M_v) <= 241/150.
```

Consequently `E[M_v(n)^2] <= nV` for **every initial configuration**.
Mixtures and every native quantum density have the same bound for this record
family: every complete source Kraus path is a partial injection, whose effect
is diagonal, so clock-record probabilities depend only on the initial
populations. Untouched reference correlations do not alter those marginal
probabilities. This does not replace the conditional quantum output by a
classical state.

The bounded endpoint correction gives an explicit finite error bound:

```text
Pr(abs(K_v(n)/n - 1/6) >= eps)
    <= min(1, n V / (n eps - R)^2),       when n eps > R.
```

When `n eps <= R`, use the bound one. Also
`E[(K_v(n)-n/6)^2] <= 2nV+2R^2`. Both statements prove the convergence above
without a stationary start, a renewal assumption or an independent clock.
Finite-start bias is retained: `abs(E[K_v(n)]-n/6)<=R`, rather than falsely
setting every finite mean to `n/6`. At one million movements the probability
of a local rate error at least `1/60` is below `1/400` under this conservative
bound. The receipt gives its exact rational value.

## Space appears in the clock fluctuations

For the same continuing histories, define

```text
Sigma_vw = lim_(n->infinity) Cov(K_v(n), K_w(n))/n.
```

The exact formula is the invariant mean of the products of the two martingale
increments:

```text
Sigma_vw = sum_i nu_i sum_j J_ij
    (g_v(i,j)+h_v(j)-h_v(i)) (g_w(i,j)+h_w(j)-h_w(i)).
```

This limit holds for every initial native density. For a finite irreducible
chain, Cesaro averages of expected state occupations converge to `nu`; applying
this to the bounded conditional martingale covariance gives the formula.
The covariance contributions of bounded endpoint corrections, divided by `n`,
vanish by Cauchy–Schwarz and the `O(n)` martingale variance bound. This argument
does not need an assumption of convergence at every individual event index.

The executable checks every port pair against its native shortest-path distance:

| Native hop distance | Unordered pairs | `Sigma_vw` | `lim Var(K_v-K_w)/n` |
| --- | ---: | ---: | ---: |
| 0, same port | n/a | `357011/1309500` | `0` |
| 1, a repair seam | 30 | `11281/261900` | `50101/109125` |
| 2 | 30 | `-98989/1309500` | `608/873` |
| 3, antipodal | 6 | `-144091/1309500` | `287/375` |

The three distinct off-diagonal values identify all three distance shells;
the first recovers exactly the thirty source seams. Thus **local clock
fluctuations encode the carrier's native spatial relation**. This is a finite,
unit-independent connection between two actual readout families on one source.
The covariance matrix has zero row sums because the total count is exactly
`2n`; treating the clocks as independent would violate that identity and gives
the wrong single-clock noise as well.

These are theoretical asymptotic covariances, not values that one short history
estimates without error. Recovering the distance classes experimentally requires
sufficient complete records, an estimator with its uncertainty and a controlled
long-block limit. Independent repetitions or a justified mixing analysis must
support such an estimator. This packet supplies the exact target and finite
single-clock error bound, not a finite-sample covariance confidence procedure,
physical lengths, a continuum dimension or Lorentzian spacetime.

## What is removed, and what still matters

Independent identity holds with any fixed positive activity multiply each move
probability and escape probability by the same factor. Their ratio `J_ij` is
unchanged, as is the complete movement-indexed clock law. Attempt costs and
finite attempt caps change. Holding only one part of an otherwise evolving
system is outside this identity.

Relative seam weights still matter. Under the uniform schedule, each port is
incident to `5/30=1/6` of stationary movement events. Doubling the relative
weight of one seam changes this fraction to `6/31` at its endpoints. No global
unit rescaling removes that dimensionless change. Adaptive schedules,
unmodeled feedback and missing current records require their own analysis.

The nearest established source result is the
[native marked current reduction](README.md), with regional continuity in
[`RegionalContinuity.lean`](../../Lean/Screen/RegionalContinuity.lean).
The [clock comparison theorems](../../Lean/Time/ClockReadout.lean) describe
calibrated affine readouts under their premises. The present stochastic
comparison adds a source-specific finite error law; substituting an average
rate for an exact affine premise would be invalid. The proof above is direct
finite matrix and martingale algebra, not a new Lean theorem.

## A laboratory comparison after one calibration

The choice of seconds is a convention and is not an obstruction to these
relations. A lab attachment must identify an apparatus that realizes these
native currents, preparations and scheduling, and collect **paired raw
records**: every movement and local clock bit, together with independently
measured lab timestamps at shared recorded block endpoints. Lost events,
finite recording capacity, timeouts and detector error need bounds. A finite
window stores only a finite tape; no unlimited physical memory is assumed by
the finite tests.

Use one disclosed calibration block to assign a positive lab duration per
count of a chosen clock, with uncertainty. Keep the other clocks and validation
blocks out of that fit. Test their relative count laws and the distance-class
covariances on complete blocks of a declared movement count `n`; independently
check the calibrated count/lab-duration relation on disjoint blocks. One
calibration interval does not establish a constant rate against a lab clock,
and this calculation does not predict how many movements occur in an arbitrary
fixed lab-time interval. The existing source rules fix the count statistics;
the physical correspondence and lab/count relation remain to be established.

The [host-clock attachment](HOST_CLOCK.md) implements this measurement shape
for the software executor: actual native currents and independent host-counter
readings, one calibration block, three held-out checks and a rate-change control.
Its scope is an engineered implementation correspondence. No traceably
calibrated lab apparatus or external physics dataset has been identified with
these native currents, and no empirical confirmation of OPH is claimed.
There is also no identification with the Palm modular phase: native currents
advance the event counts, while the current native maps do not generate that
modular phase action.

## Replay

```sh
python -m model.native_clock.comparison --check model/native_clock/comparison_results.json
python -m pytest -q -W error model/native_clock/test_comparison.py
```

The primitive runtime executes three native attempts for both C and E, starting
from a nonstationary basis preparation. All 125 joint records and full terminal
66-dimensional matrices are compared with an independently enumerated rational
source oracle, including every no-event and zero-probability record. This
execution validates the record interface; exact source equations and the proof
above supply the arbitrary-horizon statement. Floating runtime residuals are
diagnostics, separate from the exact rational clock identities.
