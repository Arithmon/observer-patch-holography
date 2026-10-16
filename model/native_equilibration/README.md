# Preparing the stationary native sector, with an error bound

The existing [quantum mixing theorem](../observer/QUANTUM_MIXING.md) already
applies to any finite connected weighted graph, not just one elementary
observer. This packet attaches that result to the actual screen federations.
It replaces **exact stationary preparation** with **discarded native repair
attempts plus a certified finite error**, within each supplied binary
particle-number sector. It does not select that number or an initial mixture
of different sectors.

The architecture supplies the graph; the implemented C or E rule supplies
the Kraus maps. The schedule is independent uniform seam selection and a
fresh fair coin. Forgotten preparation records are not used for conditioning.
These are the hypotheses of this construction, and must not be silently
replaced by the production observer's adaptive schedule.

## The existing theorem and its concrete attachment

For `V` ports, `M` seams and fixed total `k`, write `d=binomial(V,k)`.
The source's matrix-unit difference blocks establish self-adjoint positive
channels. Their diagonal block is symmetric exclusion. The theorem of
[Caputo, Liggett and Richthammer, §4.1.1, equation (4.1)](https://www.stat.berkeley.edu/~aldous/206-RWG/RWGpapers/caputo.pdf)
identifies its gap with the underlying graph's random-walk gap. The existing
OPH quantum proof bounds every coherence block by the same gap. Consequently,

```text
tau = I/d,       gamma = lambda_1(graph Laplacian)/(2M),
distance((Phi^n tensor id_R)(rho), tau tensor rho_R)
    <= min(1, d/2 (1-gamma)^n).
```

This is uniform over input densities and arbitrary untouched references.
It concerns the unread channel, not the accumulated record tape. It also
bounds the discrepancy of **any common subsequent complete instrument**,
including adaptive interventions. The total variation of the resulting
complete record laws cannot exceed the state trace distance. Conditioning
on a rare successful herald can amplify that error; the same unconditional
allowance cannot simply be assigned to the conditional law.

The graph certificate verifies the retained federation graph and gives an exact,
conservative gap certificate without relying on a numerical eigenvalue.
If `D` bounds its diameter, the range/path and variance bounds give
`lambda_1 >= 4/(V D)`. A breadth-first traversal supplies
`D <= min(V-1,2 eccentricity(root))`. Thus `gamma >= 2/(M V D)`.
This bound is deliberately conservative; it is not evidence of rapid
equilibration or a gap uniform in federation size.

The integer preparation budget avoids numerical logarithms and huge rational
powers. For `b=ceil(1/gamma_lower)`, Bernoulli's inequality gives
`(1-gamma_lower)^b <= 1/2`. Enough such blocks make `d/2 * 2^-blocks`
smaller than the requested tolerance. The receipt gives budgets for total
two and half filling at all three retained screen levels. Their large counts
are sufficient bounds, not measured minimum preparation times.

## What executes

On the elementary carrier, the exact gap is `(5-sqrt(5))/60`; the rational
bound `11/240` follows from `sqrt(5)<9/4`. Start with a specified complex
superposition of two configurations, enumerate actual source Kraus records,
and apply their signed-difference blocks repeatedly for the certified attempt
count. Both C and E are executed. The code then measures the complete four-word
law of two successive absolute-current bits and compares it with the tracial
preparation's law. No Hamiltonian or fitted state is used in that propagation.
The conservative budget is **374 discarded attempts**, with a uniform trace
distance allowance **`33/131072 < 1/1024`**. The particular executed initial
state has numerical distance about `1.04e-8`; that diagnostic does not replace
the uniform bound. The large-federation budgets are computed certificates,
not trillions of executed native attempts.

The execution uses the exact source-enumerated native channel with floating
matrix arithmetic; its numerical residuals are diagnostics, separate from
the analytic error certificate. Independent small-graph tests enumerate all
native matrix units, check the exact spectral identification numerically,
and propagate an entangled untouched reference. No fresh Lean theorem is
claimed. The executable attachment supplies preparation/observation error bounds
using the stated general mixing theorem.

## What this resolves, and what it does not

**Within this connected fixed-total unread regime**, the law has a unique
attracting stationary state.
An experiment may prepare it approximately with a declared error budget.
Exact stationarity at a finite attempt count is not asserted. Disconnected
graphs, conserved totals, preparation records and adaptive schedules change
the conclusion and remain explicit.

This result does not select a nontrivial modular clock. The modular flow of
`I/d` on its full matrix algebra is the identity. A Palm reference conditioned
on an event, or a reference on a reconstructed transfer algebra, is a different
operational construction. Native events and changing records continue at
stationarity; a stationary density does not make the history static.
Nor does convergence select a physical vacuum,
the number of spatial dimensions, or a relativistic field action. All counts
here are **discarded native attempts**, not elapsed physical or modular time.

The full unconstrained tensor algebra also contains coherences between the
all-empty and all-full configurations that the native difference records
cannot distinguish. The fixed-sector hypothesis is therefore essential;
there is no claim of a unique state on that larger algebra.

```sh
python -m model.native_equilibration --check model/native_equilibration/results.json
python -m pytest -q -W error model/native_equilibration
```
