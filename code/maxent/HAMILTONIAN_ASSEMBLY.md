# Hamiltonian bookkeeping must not change thermal correlations

This repair under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
starts from main `3d848b7f`. It builds on the observable-coordinate repairs
#1036/#1037 and the sector-energy repair #1046. Those repairs removed
spurious unit and energy-origin dependence. A separate defect remained:
rounding products and partial sums before cancellation changed the actual
Hamiltonian, its Gibbs state, partition function and susceptibility.

## Reproduction

Let `X,Y,Z` be Pauli matrices and `H=X tensor X+Y tensor Y+Z tensor Z`.
The three-term input `[H,H,H]` with coefficients `[1e20,1,-1e20]` sums
exactly to H. All these input entries are exactly represented binary64
numbers. Before this repair:

| Term order | Returned state | XX correlation | Least partial-transpose eigenvalue |
| --- | --- | ---: | ---: |
| `1e20 H + H - 1e20 H` | `I/4` | 0 | 0.25 |
| `1e20 H - 1e20 H + H` | Gibbs state of H | -0.9305533251 | -0.4479149938 |

The first state is separable; the second is entangled. This is an order-one
change in a measurable correlation caused by summation order, not uncertainty
in a tiny reported digit. The first regression commit, `04895b11`, records
**ten failures and two passing controls** before the implementation repair.

Compensating the sum of already rounded products is insufficient. For
`a=1+2^-27`, `b=1-2^-27`, the exact identity `ab-1=-2^-54` is lost when
ab is rounded to one. Multiplying the two coefficients by `-2^54` makes
the omitted operator equal to H itself. Applied to identity terms, the same
defect shifts log Z by one while leaving the normalized state unchanged.
The controls retain both failures, a noncommuting complex remainder, and
all six permutations of the three-term example.

A fixed Hamiltonian may have a redundant decomposition. Constructing its
Gibbs state does not require a unique multiplier representation. The MaxEnt
**inference** API still requires independent constraints modulo identity;
this repair does not remove its rank or uniqueness checks.

## Exact assembly and numerical evaluation

`hamiltonian_assembly.py` accumulates each real and imaginary component of
`sum_a lambda_a S_a` over rational numbers. All finite binary64 components
and their products are rational. No intermediate product or partial sum is
converted back to binary64. The implementation skips structural zeros and
caches repeated products, which keeps the small finite witnesses practical.

For Gibbs evaluation, the scalar `s=Tr(H)/d` is removed **exactly before
conversion**. The returned pair is a rounded traceless variation and a
separately rounded scalar. This preserves diagonal interactions hidden by a
large identity term. The full-matrix API instead rounds the complete sum
once: centering and then numerically adding the origin back can itself erase
a small diagonal entry. The Newton solver uses the same accumulator for its
orthonormal-coordinate Hamiltonian, and the Duhamel response is evaluated at
the repaired state.

Consequently permutations, regroupings and cancellations that preserve the
exact supplied sum preserve this assembly. A regrouping performed outside
the API that already rounds away input information supplies different data;
the accumulator cannot recover it. Products may temporarily exceed the
binary64 range if their final sum is representable. Conversely, a final
nonzero component that would round to zero, or an overflowing final component,
is rejected. No Gibbs eigenvalue floor or extra Hessian term is introduced.

The existing numerical Hermiticity convention is retained by taking the
Hermitian part after validating each input. The new exact-input diagnostic
requires **exactly Hermitian supplied matrices**, so its bounds cannot
silently certify a symmetrized replacement. Masks, Booleans, nonfinite data,
unmatched coefficients and lossy input conversions remain rejected.

## Complete thermal and entanglement law of the witness

Write P for the singlet projector. Direct Pauli multiplication gives
`H=I-4P`, `P^2=P`, and `Tr P=1`. For every finite real beta,

```text
Z(beta) = exp(3 beta)+3 exp(-beta),
p(beta) = 1/[1+3 exp(-4 beta)],
rho(beta) = p P + (1-p)(I-P)/3.
```

Thus the singlet has energy -3 and the three triplet states have energy 1.
Both one-site marginals are `I/2`, while each Pauli correlation is
`<sigma_a tensor sigma_a>=(1-4p)/3`. The complete energy response is
`<H>=1-4p` and `Var(H)=16p(1-p)=-d<H>/d beta`. Identical one-site marginals
therefore do not diagnose whether assembly retained the interaction.

The partial transpose has eigenvalues `1/2-p` once and `(1+2p)/6` three
times. Its negativity is `max(0,p-1/2)`. Negative partial transpose rules out
a separable state, as in [Peres's criterion](https://arxiv.org/abs/quant-ph/9604005).
The converse for this family has an explicit construction, with no numerical
separability solver. For `P_a^+/-=(I+/-sigma_a)/2`, define

```text
R_aligned  = (1/6) sum_{a=x,y,z; s=+,-} P_a^s tensor P_a^s,
R_opposite = (1/6) sum_{a=x,y,z; s=+,-} P_a^s tensor P_a^{-s}.
```

They are mixtures of product states. For `0<=p<=1/2`,
`rho(p)=(1-2p) R_aligned+2p R_opposite`. This proves the entire boundary:
the state is entangled exactly when `beta>log(3)/4`, and separable otherwise.
SymPy checks the symbolic projector, product-state identity and full
partial-transpose characteristic polynomial. Separate matrix exponentials
and numerical thermal evaluations check both sides of the boundary.
This is an established finite Heisenberg thermal family, used here to expose
and quantify a software defect; it is not a new OPH physical model.

## A uniform bound connecting assembly error to observables

For arbitrary finite Hermitian H and K in the dimensionless Gibbs exponent,
put `rho=exp(-H)/Z_H`, `sigma=exp(-K)/Z_K`, and

```text
delta = inf_c ||K-H-c I||_op
      = [lambda_max(K-H)-lambda_min(K-H)]/2.
```

No commutation or lower bound on the Gibbs eigenvalues is needed. Substitution
of the two logarithms gives the exact identity

```text
D(rho||sigma)+D(sigma||rho) = Tr[(rho-sigma)(K-H)].
```

Let `t=||rho-sigma||_1`. Quantum Pinsker in natural-log units gives
`t^2 <= D(rho||sigma)+D(sigma||rho)`. Trace duality, after subtracting cI,
gives the upper bound `delta t`. Hence

```text
||rho-sigma||_1 <= min(2, delta),
D(rho||sigma)+D(sigma||rho) <= delta min(2, delta).
```

The state bound is sharp to first order: `H=-epsilon Z`, `K=epsilon Z`
give `t=2 tanh(epsilon)` and `delta=2 epsilon`. For an observable O,
`|Tr[(rho-sigma)O]| <= inf_c ||O-cI||_op min(2,delta)`.
This converts Hamiltonian assembly error into a bound on any retained read,
independently of the dimension or an arbitrary energy origin. These are
consequences of the existing finite Gibbs and relative-entropy mathematics,
not a selection of the Hamiltonian or its physical temperature.

The diagnostic `hamiltonian_assembly_diagnostics` computes a rational bound
`u=max_i sum_j (|Re E_ij|+|Im E_ij|)` for the error E in the centered matrix,
and an exact scalar error v, then rounds these **upwards**. For Hermitian E,
`||E||_op<=u`, so the returned ideal Gibbs trace-distance bound is
`min(2,u)`. The variational formula `log Z_H=max_tau [S(tau)-Tr(tau H)]`
also gives an ideal log-partition error bound `u+v`. The same formula follows
from nonnegativity of relative entropy to the Gibbs state.

These bounds compare exact Gibbs functions of two finite Hamiltonians: the
exact supplied sum and its returned matrix/scalar representation. They do
**not** cover eigensolver error, exponential/state reconstruction, physical
measurement uncertainty, or input information lost before this API.
Existing faithful-support guards remain necessary. Independent rational
entry errors verify outward rounding; a separate 90-digit matrix exponential
checks the resulting ideal Gibbs and partition bounds below binary64 state
precision. The diagnostic requires no empirical fit or condition-number
estimate.

## Impact and reproduction

The affected consumers are the MaxEnt Gibbs constructor, Duhamel covariance,
Newton solve and independent projection replay. Existing Ising decimation
acceptance results retain their constraint counts, positive generic closure
defect, closed product subfamily and Pinsker conclusions. A replay compares
the existing receipt without rewriting its historical bytes.

The replayed generic closure defects are `0.00018189324934382215` and
`0.0002504907801854086` nats, differing from the recorded values by less
than `7e-16` nats. The product control remains closed to numerical precision
(`6.4e-16` nats). All recorded acceptance decisions are retained.

The sector Gibbs constructor and the three open audits #1048/#1049/#1050
are separate consumers/results; this branch starts directly from main.
The common finite algebra, first-law, conditional-information and protected
record conclusions are not changed. No source model, clock, simulator law,
particle parameter or physical identification is supplied. The model-choice
gates in #1025/#1026 are unchanged. No frozen registration, pinned receipt,
claim payload, paper, book or mandatory runner is regenerated.

```text
PYTHONPATH=code python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent -W error
```

The two new test files include the reproduced failures, the complete
thermal-law controls, rational roundoff checks, independent matrix
exponentials and invalid-input rejection. Numerical assembly and its bounds
are tested separately from the optimizer's convergence status.

All **45 new cases** pass. Eleven isolated implementation mutations are
rejected by those cases:

| Deliberate defect | Failing cases |
| --- | ---: |
| Restore the original assembly implementation | 25 |
| Round each partial sum | 27 |
| Round products before accumulation | 9 |
| Erase imaginary operator components | 10 |
| Ignore the last summand | 33 |
| Round diagonal entries before removing the scalar | 3 |
| Accept a final component that underflows to zero | 1 |
| Report zero rounding bounds | 3 |
| Round bounds to nearest instead of upwards | 2 |
| Omit scalar error from the partition bound | 1 |
| Certify a symmetrized replacement of non-Hermitian input | 1 |

Each mutation runs in its own copied package tree outside the repository
import path. Collection or import errors are not accepted as mutation kills.
The full affected suite passes **827 cases on Linux and 821 on Windows**,
with six existing extended-precision skips on Windows and warnings treated
as errors. Claim, axiom, reader-style, gravity-ladder and whitespace checks
pass; the independent freeze verifier reports `GENERATION_LOCK_VALID`.
