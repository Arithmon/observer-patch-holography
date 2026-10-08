# Simulator contract: OPH issue #330 radial lift

## Scope

This contract separates four executable objects:

1. the source screen covariance;
2. the physical radial uniqueness theorem;
3. the branch-typed forward projection;
4. downstream observable transfer.

The first three may support an E4 primordial packet. The fourth is E5. A Planck-shape overlay, CAMB/CLASS run, or likelihood is never an ancestor of an E4 radial receipt.

## Required artifact split

Every run writes separate immutable artifacts:

```text
source/
  geometric_q_receipt.json
  source_release_amplitude.json
  screen_tilt.json
  screen_covariance.json
  physical_mode_basis.json
  source_embedding.json
  radial_dilation_intertwiner.json      # SOURCE_DILATION branch
  radial_cross_covariances.json         # RADIAL_TOMOGRAPHY branch
  radial_null_report.json
  radial_forward_residual.json
  radial_promotion_receipt.json

diagnostic/
  planck_shape_comparison.json
  camb_or_class_outputs/
  plots/
```

No file under `diagnostic/` may occur in the ancestor DAG of a file under `source/`.

## Branch selector

Exactly one of the following is declared:

```json
{"radial_branch": "SOURCE_DILATION"}
```

```json
{"radial_branch": "RADIAL_TOMOGRAPHY"}
```

```json
{"radial_branch": "PRIOR_CONTINUATION"}
```

`PRIOR_CONTINUATION` is never promoted as source-derived E4.

## Common finite forward operator

For a discrete physical spectral measure,

\[
A_{\ell j}=4\pi Z_q^2w_j^{\log k}|\Psi_{\ell}(k_j)|^2,
\qquad C=Ap.
\]

The run records:

- branch type and background-curvature status;
- \(k_j\), weights, units, and density-of-states convention;
- raw normalized radial window and hash;
- \(Z_q\) and theorem parent;
- each \(\Psi_\ell(k_j)\);
- all retained and held-out multipoles;
- numerical precision and quadrature convergence.

## `SOURCE_DILATION` receipt

The theorem-level object is not a fitted slope. The producer first exports the
scale-labelled source-refinement orbit and the commutative square

\[
D_{s,r}U_r=U_rR_{s,r},\qquad
C_{\zeta,r}=U_rC_{{\rm src},r}U_r^*.
\]

This gives the exact residual-transfer identity

\[
D_{s,r}^{-1}C_{\zeta,r}D_{s,r}-u_r(s)C_{\zeta,r}
=U_r\left(R_{s,r}^{-1}C_{{\rm src},r}R_{s,r}
-u_r(s)C_{{\rm src},r}\right)U_r^*.
\]

The finite objects must converge on a common embedded compact-band core,
with local covariance bounds and a vanishing dilation residual on every
core vector, including the error from replacing finite dilation maps by
the physical maps.
Only then does the continuum relation

\[
D_s^{-1}C_\zeta D_s=e^{-\theta s}C_\zeta
\]

hold on the common \(d\ln k\) physical mode basis.

The receipt records:

```json
{
  "receipt": "SCR330_RADIAL_DILATION_INTERTWINER_RECEIPT",
  "physical_mode_basis_id": "...",
  "safe_logk_band": [0.0, 0.0],
  "dilation_maps": ["..."],
  "scale_ratios": [1.0],
  "source_embedding_commutator_norms": [0.0],
  "screen_covariance_naturality_residual_norms": [0.0],
  "physical_operator_residual_norms": [0.0],
  "strong_covariance_cauchy_residuals": [0.0],
  "strong_dilation_cauchy_residuals": [0.0],
  "uniform_covariance_norm_bound": 0.0,
  "finite_to_continuum_passed": false,
  "off_band_leakage": 0.0,
  "max_absolute_log_residual": 0.0,
  "rms_log_residual": 0.0,
  "tolerance": 0.0,
  "refinement_sequence": ["..."],
  "passed": false
}
```

Pass conditions:

- the source mode space is the cofinal scale-labelled refinement orbit, not one angular cut;
- all maps are sourced from the same physical scale bridge;
- the embedding square \(D_{s,r}U_r=U_rR_{s,r}\) passes and its residual converges;
- source and coarse mode projectors commute within the declared residual;
- source covariance naturality and physical covariance naturality have the same residual under \(U_r\);
- covariance convergence and the dilation residual pass on every compact log-wavenumber band of one common operator core; replacement of finite dilation maps by the physical maps has a vanishing error on that core;
- covariance survival equals the screen source cocycle, not a separately chosen exponent;
- safe-band leakage and operator residual converge under refinement;
- finite diagonal checks \(\Delta^2(bk)=b^{-\theta}\Delta^2(k)\) pass;
- source DAG has no observational ancestor.

Negative controls:

- shuffled physical mode labels;
- wrong scale-ratio sign;
- a spectrum with a planted log-periodic wiggle;
- a separately fitted \(n_s\);
- a mode basis from a different source embedding.

Each control must fail.

A nonzero global power law with nonzero tilt has an unbounded multiplication
covariance on `L2(d log k)`. A global uniform covariance norm bound cannot be
required for that limit: a bounded nonzero operator and its unitary conjugate
have equal norm and cannot differ by the scalar `exp(-theta*s) != 1`.
Uniform bounds on each fixed compact band are compatible with this branch;
they do not supply a uniform bound as the band expands.

## Thin-shell Mellin receipt

For `FlatExact` or `FlatAssumed` and a thin shell, compute

\[
I_\ell(\theta)=\frac{\sqrt\pi}{4}
\frac{\Gamma(1+\theta/2)}{\Gamma(3/2+\theta/2)}
\frac{\Gamma(\ell-\theta/2)}{\Gamma(\ell+2+\theta/2)}.
\]

The receipt records the convergence strip, every computed value, precision, and the arithmetic identity

\[
C_\ell^q=A_q\frac{\Gamma(\ell-\theta/2)}{\Gamma(\ell+2+\theta/2)}.
\]

The amplitude is source-derived from

\[
A_\zeta
=\frac{A_q}{\pi^{3/2}Z_q^2(k_\star R_\star)^\theta}
\frac{\Gamma(3/2+\theta/2)}{\Gamma(1+\theta/2)}.
\]

The code must not solve this equation for \(A_\zeta\) using measured CMB \(C_\ell\).

## Finite-window receipt

For a nonnegative normalized radial measure `W`, positive reference radius,
and `0 < theta < 2 ell`, the run evaluates the window transform and analytic bound

\[
\eta_{\ell,W}
=\frac{2\sqrt{J_\ell(\theta)}}{\theta}
\int W(dr)|r^{\theta/2}-R_\star^{\theta/2}|,
\]

\[
J_\ell(\theta)=I_\ell(\theta-2)
-\left[\ell(\ell+1)-\frac{\theta(\theta+1)}2\right]I_\ell(\theta),
\]

\[
|C_{\ell,W}-C_{\ell,R_\star}|
\le4\pi Z_q^2A_\zeta k_\star^\theta\eta
\left(2R_\star^{\theta/2}\sqrt{I_\ell}+\eta\right).
\]

The measure must have a finite `theta/2` radius moment. Without normalization,
`W=2 delta_R` gives `eta=0` while quadrupling the shell spectrum. Signed windows
need a different bound using total variation.

The receipt stores \(I_\ell,J_\ell,\eta\), projected values, the bound, and
the ratio of the quadrature difference to the bound. The Python helper evaluates
the analytic inequality in binary64; it supplies no outward-rounded interval
certificate. Its numerical error and the quadrature error must be controlled
separately before promoting a numerical inequality. Gamma recurrence and
`expm1` preserve the near-scale-invariant derivative norm and radius differences;
finite positive window weights are normalized after rescaling to avoid overflow.

## `RADIAL_TOMOGRAPHY` receipt

The input must contain cross-covariances. Auto-spectra alone fail the receipt:

\[
C_\ell(r_i,r_j).
\]

The fixed-radius uniqueness route requires a strictly positive reference
radius and a cross-covariance section admitting spherical-Hankel inversion in
the declared function or distribution class. At radius zero, all sections
with `ell >= 1` vanish regardless of the spectrum.

The receipt records:

- radial basis and measure \(r^2dr\);
- finite spherical-Hankel matrix and unitarity residual;
- conjugated covariance \(H Q_\ell H^{-1}\);
- off-diagonal leakage in \(k\);
- recovered nonnegative multiplication spectrum;
- refinement convergence;
- reconstruction on held-out radii/windows.

A finite set of auto-windows is not tomography and must fail this receipt.

## Null-space receipt

For every finite operator, publish:

```json
{
  "shape": [2, 2],
  "singular_values": [1.0, 1e-8],
  "rank_threshold": 1e-12,
  "rank": 2,
  "nullity": 0,
  "condition_number_nonzero": 1e8,
  "null_basis_hash": "sha256:...",
  "resolution_kernel_hash": "sha256:..."
}
```

The basis and resolution kernels accompany the evidence bundle, even for a
one-dimensional source branch. Source uniqueness still requires its theorem.

`radial_null_space_report` is a floating-point diagnostic of the supplied raw
matrix \(A\). A scaled SVD uses the relative cutoff
\(\max(\mathtt{rtol},\max(\operatorname{shape}(A))\epsilon_{\rm binary64})\).
Scaling a nonzero entry into the subnormal or zero range is accepted only if
rounding its exact entry-to-scale ratio changes it by at most
\(4\epsilon_{\rm binary64}\) relatively; otherwise the helper refuses that range.
The retained singular values determine the rank, nullity and basis together;
the API's `effective_threshold` corresponds to the receipt's `rank_threshold`.
`relative_cutoff` records the cutoff; `rank_metric` is `raw_operator_euclidean`.
Discarded nonzero directions are numerically unresolved, not proved elements
of the mathematical kernel. The zero matrix has rank zero and a complete
unresolved basis; its `condition_number_nonzero` is `None`. Omit that optional
condition field when packaging a rank-zero report in the receipt schema.

## Prior-selected continuation

For a declared \(p_0,Q\), the exact representative is

\[
p_*=p_0+Q^{-1}A^T(AQ^{-1}A^T)^+(C-Ap_0).
\]

The reference implementation evaluates this theorem numerically. It factors
\(Q=DLL^TD\), where \(D=\operatorname{diag}(\sqrt{Q_{ii}})\) and \(L\) is the
Cholesky factor of the normalized precision. Set \(W=D^{-1}L^{-T}\).
Each nonzero row of \(AW\) and of the corresponding target \(C-Ap_0\) is
divided by that operator row's largest absolute entry; call this diagonal
row scaling \(E\). No Gram matrix \(AQ^{-1}A^T\) is formed. One SVD defines:

\[
B=EAW=U\Sigma V^T,\qquad
p=p_0+WV_r\Sigma_r^{-1}U_r^TE(C-Ap_0),\qquad
R_Q=WV_rV_r^TW^{-1},\quad N_Q=I-R_Q.
\]

`effective_rank` counts retained directions using the same relative cutoff
rule as the null report; `rank_metric` is
`row_equilibrated_prior_whitened_operator`. It can differ from the raw-\(A\)
rank because whitening and row scaling change relative singular values.
The projectors are complementary and
\(Q\)-self-adjoint up to roundoff; truncation alone does not establish
\(AN_Q=0\). Neither the computed rank nor the continuation is an interval
certificate or an exact-arithmetic proof of the constrained minimum.
Before Cholesky, the computed smallest eigenvalue of \(H=D^{-1}QD^{-1}\)
must be positive, with \(\kappa_H=\|H\|_\infty/\lambda_{\min}(H)\) and
\(\eta_Q=n\epsilon_{\rm binary64}\kappa_H\le10^{-7}\), \(n=\dim Q\).
Before inverting the retained spectrum, the solver additionally requires
\[
\eta_Q+\max(\operatorname{shape}(B))\epsilon_{\rm binary64}
\sqrt{\kappa_H}\,\kappa_B\le10^{-7},\qquad
\kappa_B=s_0/s_{\rm last}\quad(0\text{ at rank zero}).
\]
This joint first-order conditioning policy accounts for prior-metric error
and amplification through whitening and inversion; it is neither a certified
error bound nor a relative-accuracy guarantee for every output component.

The inverse and residual interfaces validate original scalar entries before
array coercion and require finite real data; `rtol` is finite with
\(0<\mathtt{rtol}<1\). Zero \(A\) with zero \(C\) is supported.
Nonintegral scalars may round to binary64 within \(4\epsilon_{\rm binary64}\)
relatively, and the computation refers to those accepted binary64 values.
Integral-valued inputs must be exactly representable, regardless of scalar
type or container.
Returned constraints are checked row by row using exact arithmetic on the
accepted binary64 inputs and returned \(p\), without an absolute unit floor.
Unresolved constraints or unrepresentable reported quantities cause an
explicit refusal. Residual norms use scaled hypot evaluation; the objective
is rounded from the exact quadratic expression for the returned vector.
An original-row check of \(AN_Q\) prevents truncation from hiding a resolved
constraint. If rounding \(p_0+\delta\) changes the intended correction by
more than \(10^{-7}\) in relative \(Q\)-norm, the helper refuses that result.

The run publishes \(p_0,Q,R_Q,N_Q\) and sensitivity to declared prior variants.
The helper imposes no positivity constraint: any required positivity active
set needs separate evidence. Its output remains `ConditionalRadialContinuation`,
which cannot by itself promote a source-derived E4 claim.

## Forward residual

Every branch computes a forward residual without re-optimizing source parameters:

\[
r_\ell=C_\ell^q-4\pi Z_q^2\int d\nu(k)\Delta_\zeta^2(k)|\Psi_\ell(k)|^2.
\]

The residual artifact contains raw signed residuals, absolute and relative
norms, numerical error budget, and held-out modes/windows. A failed residual
blocks promotion. `forward_residual` computes \(C-Ap\) from the accepted inputs
without fitting; norms must not vanish merely because squaring underflows.
Its relative norm uses \(\|C\|_2\) without a denominator floor. If \(C=0\),
the relative residual is zero only when the residual is also zero; otherwise
the relative quantity is undefined and the helper refuses explicitly.

## Numerical audit and validation

At main `1ba8a011c4a81fe0680b14dadf4802388a0567a2`, the original inputs
\(A=\operatorname{diag}(1,10^{-8})\), \(C=(10^{-3},10^{-11})\),
\(p_0=0\), \(Q=I\), with default `rtol`, returned \(p=(10^{-3},0)\),
`effective_rank=2` and \(N_Q=\operatorname{diag}(0,1)\).
The Gram cutoff discarded a resolved direction and the absolute residual
floor accepted its missing constraint; the unique solution is
\(p=(10^{-3},10^{-3})\).

The regression and independent-control files retain this failure, exact
Fraction/KKT solutions, metric-dependent ranks, projector identities, zero
maps, original-input validation, and small/large-unit residual controls.
`radial-inverse.yml` executes these two files and `test_oph_radial_lift_330.py`
on Ubuntu and Windows with pinned dependencies and warnings treated as errors.
This repair changes the live numerical helper. It requires no retained
receipt regeneration, frozen-evidence rewrite, paper theorem change, or
physical claim promotion.

## Curved branches

- `FlatExact`: exact flat formulas may be promoted if all other receipts pass.
- `FlatAssumed`: at most conditional physical status.
- `OpenCurved` / `ClosedCurved`: use the declared hyperspherical eigenfunctions and spectral measure; do not use the flat gamma conversion.
- `Unresolved`: fail physical radial promotion.

## E4 promotion

`SCR330_RADIAL_PROMOTION_RECEIPT` passes only if:

- all source/stress/clock/freezeout/phase/scale receipts pass;
- \(A_q,\theta,Z_q,W\) are source outputs;
- the complete null report and forward residual pass;
- exactly one uniqueness branch passes;
- the source DAG is clean;
- `physical_tt_te_ee_claim` is false.

## E5 firewall

A transfer artifact may consume the E4 primordial packet, but the E4 packet may not consume transfer outputs. TT/TE/EE, lensing, recombination, foregrounds, nuisance parameters, covariance, and likelihood live in a separate DAG and separate claim tier.

## Reference implementation

`oph_radial_lift_330.py` implements:

- exact \(I_\ell(\theta)\) and \(J_\ell(\theta)\);
- general-pivot/general-\(Z_q\) amplitude conversion;
- exact thin-shell gamma spectrum;
- finite-window quadrature, analytic bound, and separate numerical error budget;
- finite projection matrix and SVD/null report;
- minimum-prior continuation;
- source-family forward residual;
- finite dilation residual checks;
- approximate-dilation shape bound;
- fail-closed receipts.

`pytest -q` must pass before a receipt bundle is accepted.
