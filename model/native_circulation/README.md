# Gaussian circulation from one continuing observer history

A single native observer can accumulate a circulating current with a Gaussian
long-count limit. Independent copies of the whole history are unnecessary for
this result. The source fixes the limiting covariance. Stationary binary
preparation, the uniform seam scheduler and access to the selected seam/current
record remain declared.

The concrete measurement counts signed circulation around the lexicographically
first triangular face of the source graph. Each native attempt contributes
`g=-1,0,+1`; `S_N` is their retained sum after `N` attempts. The exact prediction is

```text
E[S_N] = 0,               Var(S_N) = N/66.
```

Successive signed increments have zero covariance, but are not independent:
`E[g_n^2 g_(n+1)^2]=1/1980`, while independent increments would give `1/4356`.
The receipt includes the complete finite counter laws for 1,4,16,64,256 attempts,
with every counter value and zero-probability outcome retained. These predictions
need neither hidden initial-load measurement nor a conversion to seconds.

## Source derivation and the existing theorem attachment

Orient each of the thirty seams from its first to second listed port. Let
`B` be the 12-by-30 incidence matrix, negative at the source and positive at
the target. The vector `J_n` records the signed current on the selected seam
and zero on all other seams. Actual native branch enumeration verifies

```text
x_(n+1)-x_n = B J_n,
E[J_n | x_n] = -B^T x_n/60,
E[J_n J_n^T] = I30/198                    (stationary preparation).
```

Put `L=B B^T` and `P=I30-B^T L^+ B`. Exact rational matrix arithmetic checks
that `P` is the orthogonal projection onto `ker B`, of rank 19. Thus
`P J_n` is a bounded martingale difference with respect to the complete native
history, and hence also with respect to the smaller accessible record history.
Its stationary covariance is `P/198`; unequal-attempt covariances vanish by
conditional expectation.

For any seam weights `a`, the exact decomposition is

```text
sum_(n<N) a^T J_n
 = sum_(n<N) (P a)^T J_n + (L^+ B a)^T (x_N-x_0).
```

The second term is bounded on the finite carrier. A gradient-weighted current
therefore has zero asymptotic diffusion. For example, the sum of signed load
changes at one port is just its final minus initial binary load, always between
minus one and one. It cannot have variance proportional to arbitrarily large N.
The triangular weights have `B a=0` and `||a||^2=3`, yielding `1/66` exactly.
The cycle rank 19 counts independent graph currents; it is not a spatial-dimension
prediction.

The nearest existing finite results are
[`GreenKubo.lean`](../../Lean/Thermodynamics/GreenKubo.lean) and
[`NativeRepairFourCycle.lean`](../../Lean/Thermodynamics/NativeRepairFourCycle.lean).
They establish reversible transport/Poisson identities and a different finite
witness. This packet supplies the actual icosahedral binary transition and its
accessible edge-current observable, rather than identifying a supplied four-cycle
with the native carrier. The gradient/circulation decomposition resolves the
finite Poisson correction explicitly. No new Lean central-limit theorem is
claimed.

## The limit and its action

The configuration chain is finite, irreducible and aperiodic: connected seam
exchanges connect all two-token configurations, and each configuration has
positive self-loop probability. Its stationary law is uniform. For the bounded
martingale differences `P J_n`, the normalized sum of conditional covariance
matrices converges to `P/198` by the finite-chain ergodic theorem. Bounded
increments satisfy the martingale Lindeberg condition. The martingale functional
central-limit theorem therefore gives, under stationary preparation,

```text
Y_N(u) = (1/sqrt(N)) sum_(n<floor(Nu)) P J_n
       -> W(u),            Cov(W(u),W(v)) = min(u,v) P/198.
```

Here `u` is a rescaled native-attempt count. It is not calibrated proper time.
The standard martingale method and bounded-endpoint correction are reviewed in
[Kipnis and Varadhan (1985), Sections 3–4](https://www.numdam.org/article/AST_1985__132__65_0.pdf).
Our exact martingale is the simpler finite case; we do not apply their general
state-functional theorem to an unverified edge process. The finite-state
hypotheses and covariance are checked above, and the probabilistic limit theorem
is imported, not implemented by a finite test.

Equivalently, the one-dimensional limit can be obtained from the finite tilted
kernel `K_theta(x,y)=sum_branches P(branch|x) exp(i theta g)`.
`K_0` has a simple eigenvalue1 and all other eigenvalues strictly within the
unit circle. The eigenvalue branch has first derivative zero and second
derivative `-1/66`: its first eigenvector derivative is zero because the
conditional signed drift vanishes. Analytic perturbation at zero then yields
`E exp(i theta S_N/sqrt(N)) -> exp(-theta^2/132)`.

For a fixed partition `0=u_0<...<u_m`, the limiting increment density on the
19-dimensional cycle subspace has quadratic action

```text
S(delta y_1,...,delta y_m) = 99 sum_k ||delta y_k||^2/(u_k-u_(k-1)).
```

This is a finite-partition Brownian increment density, with starting value zero.
It is not a density over differentiable continuous paths. It differs from the
[replica fluctuation action](../native_fluctuations/README.md): the observable,
limit and hypotheses are different. This result supplies no spatial refinement,
relativistic propagation, interacting field, coherent repair selection or
physical duration.

## What executes and what the finite results show

The native C Kraus maps are grouped by their actual triangular current.
The shared primitive runtime executes three consecutive reads from `I66/66`,
retaining all 27 record words and complete 66-dimensional terminal densities.
Every output is compared with independent exact native counts. The checkpoint
after one read is serialized, replayed and continued; a separate test checks
an untouched entangled reference. No native pulse-word compilation is claimed.

For longer histories, an exact integer recursion retains the cumulative counter
and every terminal configuration. It sums over earlier record orders because
the declared readout retains only the counter and performs no order-dependent
feedback. All 60 scheduler/coin alternatives per native attempt are included;
there is no Monte Carlo sampling or probability-tail truncation. The three-step
counter coarsening is checked against the full quantum runtime. The receipt
stores complete counter probabilities and hashes of the computed joint
counter/terminal-state integer arrays.

The limit does not imply a Gaussian finite record law. At 256 attempts the
variance is `128/33` and the excess kurtosis is about **0.6717113**. For the
variance-normalized counter, the characteristic function at1 is about
**0.62319322**, compared with **0.60653066** for a standard Gaussian. At 64 attempts
the corresponding values are **0.66078679** and excess kurtosis **2.30935635**.
The exact finite laws are the predictions at those counts. These comparisons
are numerical diagnostics; no uniform Berry–Esseen bound is asserted.

There is no new laboratory dataset or empirical agreement claim. The result
removes the independent-replica premise for this particular long-count observable,
while retaining the stated finite source and measurement design.

## Reproduce

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.native_circulation --check model/native_circulation/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/native_circulation
```
