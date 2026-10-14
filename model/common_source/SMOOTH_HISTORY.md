# Complete-output certificates for the continuing smooth observer

[smooth_history.py](smooth_history.py) bounds the complete quantum and
classical output of the existing [smooth continuation](COMPOSITION.md).
It follows both smooth reads, their disturbed states, the later feedback
and Hamiltonian intervals, and the final gauge clock probe. The canonical
program has sixteen retained records. Each native preparation current
`j=-1,0,1` is evaluated separately.

This calculation attaches the polynomial certificate method of
[INFINITE_HISTORY.md](INFINITE_HISTORY.md) to that actual continuation.
The earlier sharp-output comparison concerns a different detector and
history. The general moment bound in [COMPOSITION.md](COMPOSITION.md)
continues to provide its fixed-program occupation convergence theorem;
the present state calculation supplies tighter finite comparisons.
Polynomial states are verification objects, not replacements for the
source-executed finite states.

## Target, preparation and domains

The target is the declared two-site Hamiltonian at `kappa=9/32`, `g=1/4`,
with the blank gauge input and the specified free coupled scalar vacuum.
Hamiltonian times remain model coordinates. The source preparation current
applies `I` for `j=0` and `(3I-4ijX0)/5` for `j=+/-1`. Committed-bit
feedback applies `(3I-4iX0)/5` when its selected bit is one. These exact
rational matrices preserve the norm and scalar occupation.

Use unnormalized Bargmann basis vectors
`|n0,n1,g>_B=(a0*)^n0(a1*)^n1|0,g>`. Their squared norms are
`n0! n1!`. A polynomial stores Gaussian-integer coefficients with one
positive integer denominator. Every inner product, coefficient rounding
residual and discarded-support norm is computed in this factorial metric.

The infinite preparation approximant is the chosen total-occupation
truncation of

\[
 u_0=\exp[-(a_0^*-a_1^*)^2/36]|0,0,0\rangle.
\]

Let `S` be its exactly computed squared norm. It is normalized once by
`sqrt(S)`, shared by all its descendants. The full unnormalized Gaussian
has normalization factor squared `4sqrt(5)/9`. Thus the omitted probability
is bounded above using a directed lower bound on `sqrt(5)` by
`p_tail <= 1-(4sqrt(5)_lower/9) S`. The pure-state trace distance from this
normalized truncation to the true vacuum is at most
`tau=sqrt(p_tail)`, with an outward root bound.

At finite onsite dimension `d`, the target preparation is exactly the
normalized total-occupation projection `N<=d-2`, embedded into the rectangle
`0<=n0,n1<d`. Its normalization `S_d` generally differs from `S`. This
finite preparation has no preparation error relative to its declared target.
Neither initial normalization is replaced by a later branch mass.

For infinite evolution, powers of `H` act on their full expanding support.
Truncating between those powers would change the Taylor polynomial and its
remainder. For finite evolution, the generator is `H_d=Pi_d H Pi_d`, so
compression occurs after **every generator action**. This preserves the
boundary terms inherited from the infinite quadratic Hamiltonian; it does
not replace them by squares of truncated quadratures.

Self-adjointness and the finite-number core are the same as in
[OCCUPATION_REFINEMENT.md](OCCUPATION_REFINEMENT.md). Finite polynomials lie
in the domains of all powers used below. The scalar certificate covers its
declared preparation, not arbitrary input states or a spatial limit.
The executable certificate additionally bounds its work: total Hamiltonian
time at most four, at most two smooth reads and four recorded bits, with
the admitted `Resolution` ranges. These are computation limits, not a
claim that the mathematical induction fails for longer fixed histories.

## Local unitary and pointer bounds

For a self-adjoint generator `A`, a finite polynomial `v`, and the degree
`m` Taylor polynomial `T_m(t,A)`, the unitary integral remainder gives

\[
 \|e^{-itA}v-T_m(t,A)v\|
 \le \frac{|t|^{m+1}}{(m+1)!}\|A^{m+1}v\|.
\]

The right side is evaluated from the next exact polynomial term. There is
no exponential norm factor: the remaining propagator in the integral is
unitary. Each substep then charges the norm of every removed coefficient
and the exact difference introduced by rounding retained coefficients onto
a dyadic grid. The errors use the original normalization `S` or `S_d`.

A smooth read uses `B=a0+a0*` and `alpha=1/sqrt(8)`, with compression
`B_d=Pi_d B Pi_d` in the finite case. Its complete Kraus pair is

\[
 U_-=e^{-i\alpha B},\quad U_+=e^{i\alpha B},\qquad
 K_0=(U_--iU_+)/2,\quad K_1=(U_-+iU_+)/2.
\]

This retains the circuit phases and gives effects `(I+sin(Q0))/2` and
`(I-sin(Q0))/2`, respectively. `Q0=B/sqrt(2)`. The two exponentials use
the same generator, so a rational `alpha0` with `|alpha-alpha0|<=delta`
satisfies

\[
 \|(e^{\pm i\alpha B}-e^{\pm i\alpha_0B})v\|
 \le\delta\|Bv\|.
\]

This follows directly from the spectral theorem, including for unbounded
`B` on these polynomial vectors. The implementation subtracts one grid
unit from an outward upper bound on `sqrt(1/8)` to choose `alpha0`, and
charges that grid unit times an outward bound on `||Bv||`.

For errors `d_-` and `d_+` in the two arms, the parallelogram identity gives

\[
 \|(d_--id_+)/2\|^2+\|(d_-+id_+)/2\|^2
       =(\|d_-\|^2+\|d_+\|^2)/2.
\]

Thus the complete pointer's Taylor error is the root-mean-square of the
two arm bounds, and its coefficient error is at most `delta ||Bv||`.
The two output truncation errors combine in Euclidean norm, as do their
rounding errors. All these contributions are then added. No measured
branch is normalized or removed.

## Complete-history induction

Retain the orthogonal record label with each polynomial vector. At one
stage, let `eta_r` bound the local error on the current approximate input
of old record `r`, including all its new outcomes. The stage increment is
bounded by

\[
 \eta=\left(\sum_r\eta_r^2\right)^{1/2}.
\]

The exact stage is an isometry before classical record dephasing: a
unitary evolution, a complete pointer, a complete gauge read, or a
record-controlled unitary. Hence an incoming complete-vector error `e`
becomes at most `e+eta`. This argument evaluates each local residual on
the approximate vector actually carried; it does not assume that an
approximate Taylor polynomial is itself unitary. Iterating gives the
reported `epsilon`, with zero and rare branches retained.

Let `rho_poly` be the complete CQ state formed from the final polynomials
with their one initial normalization. Its mass may differ from one. If
`rho` is the exact normalized target output, then

\[
 \tfrac12\|\rho-\rho_{\rm poly}\|_1
       \le \tau+\epsilon+\epsilon^2/2.
\]

Indeed the exact dilation vector has norm one and the polynomial vector
has norm at most `1+epsilon`; the rank-one difference is bounded by their
norm sum times their vector distance. Initial preparation contributes its
separate trace distance `tau` by contractivity. For a single bounded
record effect the implementation uses the safe direct probability bound
`tau+2epsilon+epsilon^2`. A nonnormalized polynomial state must not be
treated as though its trace difference were zero.

## Comparing finite and infinite complete outputs

For corresponding polynomial branches `u_r`, `v_r`, with their distinct
initial norms `S`, `S_d`, compute exactly

\[
 p_r=\|u_r\|^2/S,\quad q_r=\|v_r\|^2/S_d,\qquad
 z_r^2=|\langle u_r,v_r\rangle|^2/(SS_d).
\]

The squared overlap is real and nonnegative; no square root of either
normalization is needed in this calculation. Orthogonal record blocks and
the rank-one trace-norm formula give

\[
 D_{\rm poly}=\frac12\sum_r
       \sqrt{(p_r+q_r)^2-4z_r^2}.
\]

The radicands are exact rationals. Cauchy--Schwarz proves nonnegativity;
negative values are errors, not values to clip away. Directed square roots
enclose the sum. Add and subtract the two complete-output error bounds,
then intersect the resulting true-state interval with `[0,1]`. Clipping the
polynomial distance itself to one before accounting for its residuals
would be unjustified.

This comparison retains all quantum output, within-branch phases and
record probabilities. Multiplying an entire record branch by a scalar
phase leaves it invariant, as required for the committed interface in
[SOURCE_SELECTION.md](../SOURCE_SELECTION.md). Equal probabilities alone
do not imply a small distance when relative phases inside a branch differ.

## Native attachments and remaining errors

The finite continuous-work comparison may be widened by the existing sum
of product-formula errors for **all** Hamiltonian intervals. That yields
an ideal source-word comparison. Floating source execution residuals and
physical service errors remain separate obligations; the polynomial
calculation does not certify either by assigning them zero.

For the native attachment, prepare the declared work input independently
of an arbitrary native/reference input. The first complete native record
selects one of the three certified currents. Each block is a product of
its unnormalized native/reference state and the corresponding complete
work state. Its weight multiplies the work trace-distance bound; the
native weights sum to one. Therefore the maximum bound over the three
currents controls the complete joined output. Later record-controlled
native repairs and encoded recoveries are common CPTP continuations, so
contractivity preserves this bound without a branch-count multiplier.

This prepared-work statement does not certify arbitrary initially
entangled native/work inputs or arbitrary gauge preparations. Those inputs
retain the finite instrument identities of [QUANTUM_NATIVE.md](QUANTUM_NATIVE.md)
and [ENCODED_MEMORY.md](ENCODED_MEMORY.md), with their own theorem scopes.
The encoded logical/reference marginal's exact work-independent identity
also remains distinct from this approximation of the full observer output.

An `Enclosure` is an in-process proof object whose digest detects accidental
mutation. The digest is public and is not evidence for externally supplied
error claims. Published evidence is checked by exact recomputation.
Independent tests cover the factorial generator metric, finite compression
between powers, exact vacuum displacement, complete pointer normalization,
support/rounding charges and phase-sensitive block comparisons. Dense finite
history controls separately check the full retained continuation.

## Reproduction

The receipt recomputes all three currents, both executed finite cutoffs,
their complete source/product checks and a coarser infinite-history control.
It reports continuous finite-to-infinite distances separately from the
intervals widened for the finite product formula. The source resource
counts belong to the executed programs; polynomial work is host verification.

```sh
python -m model.common_source.smooth_history --check model/common_source/smooth_history_results.json
python -m pytest -q -W error model/common_source/test_smooth_history_certificate.py model/common_source/test_smooth_history_adversarial.py
```
