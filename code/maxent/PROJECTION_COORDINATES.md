# MaxEnt inference must not depend on observable coordinates

This repair addresses the constraint-scale finding in
[the maintainer's #1033 audit](https://github.com/FloatingPragma/observer-patch-holography/issues/1033#issuecomment-5986930536).
The defect class is numerical dependence on an arbitrary presentation of a
fixed finite Gibbs family. The implementation previously changed rank,
convergence and even the predicted state when only units or scalar energy
origins changed. The fixes act on the actual observable space modulo identity.

## Reproduced failures

At main `9582fdc3`, the first regression commit reproduces six failing cases
and three passing controls. With `Z = diag(1,-1)`, `X = [[0,1],[1,0]]`:

| Input | Previous result | Correct control |
| --- | --- | --- |
| Project `diag(.7,.3)` onto `[1e-5 Z]` | Rank zero; rejects | Rank one, multiplier `-atanh(.4)/1e-5` |
| Same with units `1e-100` or `-1e-5` | Rank zero; rejects | Same state and inverse-scaled multiplier |
| `[Z + 1e8 I]` | Rank zero; rejects | Same Gibbs family as `[Z]` |
| `[Z, Z + 1e-5 X]` | Rank one; rejects | Rank two, the same family as `[Z,X]` |
| Gibbs state for `[1e20 I + X]`, multiplier `.4` | `I/2` | `(I - tanh(.4) X)/2`; off-diagonal error about `.19` |

A second, distinct failure is premature convergence. For target
`(I + 4e-8 Z)/2` and observable `1e-4 Z`, the old rank test accepts the
family, but its raw gradient is `4e-12`. With tolerance `1e-11`, it returns
multiplier zero and `I/2` without an update. The projected Hilbert--Schmidt
moment error is about `2.83e-8`. Changing units hid an unchanged state error.

The failing controls preceded the fix in git history. They use analytic
Bernoulli/Bloch answers, not comparison with a second call to the producer.

## The object being solved

Let `S_a` be finite Hermitian operators on `C^d` and let the supplied target
`sigma` be a faithful normalized state. Put

\[
 V=\operatorname{span}_{\mathbb R}
       \{S_a-\operatorname{Tr}(S_a)I/d\},\qquad
 \rho_\lambda={e^{-\sum_a\lambda_aS_a}\over
                    \operatorname{Tr}e^{-\sum_a\lambda_aS_a}}.
\]

The solver requires the constraints to be independent modulo identity, so
the original multiplier vector is unique. Counting rank is useful even for
redundant families; solving rejects those families instead of choosing
unadvertised multipliers or deleting constraints. Rank is computed from
normalized real coordinates, not from a Gram matrix whose conditioning is
squared and whose absolute threshold depends on units.

For an invertible real matrix `A` and real vector `b`, set
`S'_a = sum_b A_ab S_b + b_a I`. Then `lambda = A^T lambda'` gives

\[
 \rho'_{\lambda'}=\rho_\lambda,\quad
 \log Z'(\lambda')=\log Z(A^T\lambda')-b^T\lambda',\quad
 m'=Am+b.
\]

Thus the objective `f = log Z + lambda^T m` and the inferred state are
unchanged. The coordinate gradient and Hessian transform as `g'=Ag` and
`K'=AKA^T`; their raw norms and smallest eigenvalues are not invariant.
An absolute Gram cutoff, raw-gradient stopping condition, or fixed Hessian
ridge is consequently not an intrinsic test of this inference problem.

Choose a Hilbert--Schmidt orthonormal Hermitian basis `Q_j` of `V`, and write
the traceless Hamiltonian `H(theta)=sum_j theta_j Q_j`. In these coordinates

\[
 f(\theta)=\log\operatorname{Tr} e^{-H(\theta)}
              +\operatorname{Tr}\sigma H(\theta),\qquad
 g_j=\operatorname{Tr}(\sigma-\rho_\theta)Q_j.
\]

The stopping quantity is

\[
 \eta=\|g\|_2=\|P_V(\sigma-\rho_\theta)\|_{\rm HS}.
\]

Any other orthonormal basis of the same real space is related by an
orthogonal matrix. Therefore `eta` and `||theta|| = ||H||_HS` are invariant,
including under the original affine change of observable coordinates.
Common unitary conjugation preserves them too.

## Existence, uniqueness and an error bound

These are finite mathematical statements about the declared input family,
not physical laws selected by the observer axioms.

**Existence without a fitted Hessian floor.** Let `mu=lambda_min(sigma)>0`
and `d>=2`. If `a=-lambda_min(H)`, then `B=H+aI` is positive semidefinite
and `Tr B=da`. Hence

\[
 \|H\|_{\rm HS}^2=\operatorname{Tr}B^2-da^2
      \le d(d-1)a^2.
\]

Since `log Tr exp(-H) >= a` and `Tr sigma=1`,

\[
 f(\theta)\ge\operatorname{Tr}\sigma(H+aI)
       \ge\mu da
       \ge\mu\sqrt{d/(d-1)}\,\|\theta\|_2.
\]

Thus `f` is coercive, so it attains its minimum. Because `f(0)=log d`, its
minimizer satisfies

\[
 \|\theta_*\|_2\le R={\log d\over\mu}\sqrt{(d-1)/d}.
\]

Here `R` bounds optimizer coordinates; it is not a spatial radius.

**Uniqueness.** For a Hermitian direction `T`, the Hessian quadratic form
is the Kubo--Mori covariance of `T`. In an eigenbasis of the faithful state,
every logarithmic-mean weight is positive. The form vanishes precisely when
`T` is scalar. In the traceless independent basis this means `T=0`, so the
Hessian is positive definite and `f` strictly convex. The minimum is unique
and its gradient is zero, giving moment matching. This proves the finite
existence and uniqueness premise consumed by the closure calculation.

**An error bound for any supplied candidate.** Convexity gives

\[
 0\le f(\theta)-f(\theta_*)
   \le g(\theta)\cdot(\theta-\theta_*)
   \le \eta(\|\theta\|_2+R)=:B(\theta).
\]

Moment matching at the optimum implies the exact Pythagorean identity

\[
 D(\sigma\|\rho_\theta)
 =D(\sigma\|\rho_*)+D(\rho_*\|\rho_\theta),\qquad
 f(\theta)-f(\theta_*)=D(\rho_*\|\rho_\theta).
\]

To verify it directly, insert `log rho_theta = -H(theta)-log Z(theta)`
and cancel `Tr(sigma-rho_*) Q_j=0`. No commutation between constraints is
required. Quantum Pinsker then gives

\[
 \|\rho_*-\rho_\theta\|_1\le\min\{2,\sqrt{2B(\theta)}\}.
\]

This separates **optimization error** `D(rho_* || rho_theta)` from the
**closure defect** `D(sigma || rho_*)`. A correct solver cannot make the
latter disappear. For the controlled qubit example with target
`(I+.3X+.2Y+.4Z)/2` and constraints `[X,Z]`, the optimum is
`(I+.3X+.4Z)/2`: its nonzero closure defect is retained.

The familiar finite exponential-family projection and Pythagorean facts
are established mathematics; see Stephan Weis,
[*Information topologies on non-commutative state spaces*](https://arxiv.org/abs/1003.5671).
His [discussion of boundary Gibbs families](https://arxiv.org/abs/1411.0015)
also explains why the faithful finite setting must not be silently extended
to zero-temperature or singular-target limits. The coercivity estimate above
makes the numerical stopping implication explicit for this repository's
already-declared faithful-target input.

## Implementation and precision boundary

`information_projection.py` separates the generic inference operation from
the Ising acceptance fixture. `maxent_closure_acceptance.py` reexports the
established function names and retains the legacy keyword conventions.
Observable decomposition and thermal diagonalization now live in
`quantum_information/gibbs.py`, shared with the
[sector Gibbs constructor](../collar_alignment/GIBBS_SECTOR_AUDIT.md).
Its input conversion rejects any lost numeric component, including integer
energy gaps rounded away on conversion to binary64.

- Each observable is decomposed into a scalar offset and a bounded
  traceless variation. A diagonal anchor is subtracted before normalization,
  so a large scalar cannot erase an independently stored off-diagonal term.
  These normalized operators supply derivative and inference coordinates.
  The [Hamiltonian assembly repair](HAMILTONIAN_ASSEMBLY.md) accumulates the
  original weighted entries exactly before separating the scalar from the
  thermal matrix. This also prevents multiplication or partial-sum roundoff
  from erasing an interaction. The scalar is restored only to `log Z`.
  Neither operation recovers differences already lost in the caller's input.
- The real traceless Hermitian coordinates use Helmert diagonal components
  and the real/imaginary upper triangle. SVD stays inside that space; it
  cannot create an anti-Hermitian or scalar direction from a small singular
  value. The normalized rank threshold is `64 eps max(shape) s_max`.
  A below-threshold direction is unresolved, not a certified exact dependence.
- Newton uses the actual Kubo--Mori Hessian with no ridge. The logarithmic
  mean uses `expm1` to avoid cancellation. A line search reduces trials that
  leave resolved faithful support; it raises on exhaustion. An objective
  change with no resolved strict Armijo decrease is accepted only with a
  halved gradient norm. Residual norms use scaled `hypot` accumulation so
  a representable error cannot disappear by squaring its components.
- Conversion overflow, erased nonzero components, malformed observables,
  mixed Boolean inputs, unrepresentable output covariance, unresolved
  support and exhausted iteration budgets raise explicit exceptions.
  Nothing silently clips a Gibbs eigenvalue or repairs an invalid state.
- The returned original multipliers are replayed against the original
  observables. Internal-coordinate convergence alone does not pass the gate.

`project_information` returns the state, original multipliers, invariant and
raw residuals, iteration count, `B(theta)` and the trace-distance bound.
`projection_diagnostics` evaluates the same candidate quantities without
running the optimizer or trusting its status. Its inequality is proved
above; its numbers are ordinary floating-point evaluations, **not interval
certificates**. Rank, eigenvalues, moment evaluation and residuals retain
roundoff uncertainty, so a reported zero is not an exact equality proof.
Poorly conditioned reparametrizations or unresolved support can still be
rejected rather than reported as a solution.

`hamiltonian_assembly_diagnostics` separately supplies outward-rounded bounds
on exact-input Hamiltonian assembly and its **ideal** Gibbs consequences.
It requires exactly Hermitian inputs and does not include eigensolver or
state-reconstruction error. Its rational certificate does not promote the
optimizer's floating residual or reported state to an exact certificate.

The compatibility `i_projection` pair remains `(multipliers, raw_residual)`.
Its tolerance now controls the invariant residual; a raw residual in very
large units need not be below that tolerance. Supported tolerance is
`0 < tol <= 1e-6`. This is an intentional convergence-contract correction.
Callers needing both quantities should use `project_information`.

## Independent controls and impact

Tests check analytic qubit and qutrit optima, exact rational constraint rank,
noncommuting Hessians against SciPy Frechet derivatives of the matrix
exponential, high-precision logarithmic means, unitary and affine covariance,
and entropy identities against separate matrix logarithms. Candidate bounds
are tested away from the solution; an incorrect state in the final replay
must fail even after internal Newton convergence. Diagnostics are exercised
with the optimizer disabled. Guards are also run under `python -O`.

| Surface | Effect |
| --- | --- |
| Live MaxEnt closure, #539 | All established count, nonclosure and product-subfamily conclusions survive. Original numeric outputs change at roundoff; the live JSON gains invariant residual and optimizer-error diagnostics. |
| Shared quantum information, collar and Einstein receipts | Uses the established faithful-state validation; no shared entropy/channel definitions change. Their integration suites are rerun. |
| Archived simulator MaxEnt copy | Historical source and receipts retain their bytes; this PR repairs the live repository implementation, not that archived copy. |
| Papers, claim payloads, frozen registrations and custody | Unchanged. These mathematical and numerical repairs supply no new empirical claim or physical identification. |
| #1025/#1026 | No candidate is selected and no positive physics start gate is asserted. |

The live receipt's old `moment_matching_residual` remains the raw norm.
The added fields are `normalized_moment_matching_residual`,
`projection_optimality_gap_bound_nats`, `projection_trace_distance_bound`
and `projection_iterations`. `projection_unique` now uses independence
modulo identity together with the faithful-target premise already checked
by the solver, instead of a units-dependent absolute Hessian threshold.
The raw Hessian minimum remains an informational field in the fixed units.

Reproduce from the repository root with the pinned dependencies:

```text
python -m pytest -q code/maxent
python code/maxent/maxent_closure_acceptance.py
python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent
```

The second command intentionally regenerates only the live #539 receipt.
The existing Windows/Linux quantum-information workflow runs the entire
expanded suite. The mandatory #539 acceptance tests consume the new
diagnostics as well as the original scientific checks.

## Maintainer-style audit follow-up

The audit reproduced a false convergence in the first PR head: for target
`I/2 + 5e-201 X`, the invariant residual is `1e-200/sqrt(2)`, but the direct
sum-of-squares norm underflowed to zero. A requested tolerance of `1e-250`
then accepted the unchanged maximally mixed state and reported a zero error
bound. Both the independent diagnostic and the solve now use `hypot`; the
diagnostic retains the nonzero error, and the solve explicitly rejects the
unresolved tolerance. A rounded-to-zero Armijo decrement no longer passes
the line search as objective progress. The retained tests failed before
these corrections. The existence, uniqueness and error-bound proofs are
unchanged; this repairs evaluation of their numerical premise.

The same audit found that target matrices reached the shared state validator
before the conversion guards used for observables and multipliers. A Linux
extended-precision target with nonzero real or imaginary coherence `1e-400`
was therefore converted into the maximally mixed state and reported as an
exact numerical match. Both projection entry points now guard target
conversion too, rejecting erased components and mixed Boolean values before
the existing faithful-state checks. The extended-range regressions run on
Linux and skip explicitly on platforms without that numeric type.
