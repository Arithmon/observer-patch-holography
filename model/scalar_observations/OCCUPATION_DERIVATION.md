# Certifying occupation error through the continued smooth history

This calculation compares the principal history in [README.md](README.md)
with the same seven **unrestricted** oscillators. It keeps both field reads,
their disturbance, the conditional kick and every surrounding native record.
The torus, seven-mode band, bare action, product blank, detector and
Hamiltonian model times are those supplied in the original packet. Removing
an oscillator occupation cutoff does not remove the field bandwidth or
identify these choices with measurements of nature.

The occupation target uses the exact spatial gradient, `kappa=1`. The
retained grid executions first compare to that same-gradient finite code
using the existing spatial and temporal budgets. Those two comparisons may
then be added. No irrational spatial coefficient is rounded into the exact
occupation calculation.

The retained calculation has two roles. Four-level certificates compare
the existing finite implementation with the complete unrestricted history.
A larger auxiliary polynomial calculation, followed by exact unrestricted
terminal effects, resolves an interaction-dependent **classical record
event**. The auxiliary dimensions are verification resources, not a claim
that a larger source implementation was compiled or executed.

The nearest existing implementation is the
[local residual and rounding certificate](../common_source/INFINITE_HISTORY.md)
for a different coupled-oscillator history. Here the action is quartic in
seven coordinates, and the actual smooth finite-code preparation, readout
and feedback all require comparison. Numerical monitoring of bosonic
occupation error also appears in Woods, Cramer and Plenio,
[*Simulating Bosonic Baths with Error Bars*](https://arxiv.org/abs/1504.01531),
equations (16)--(18). Their quadratic-bath and bounded-interface assumptions
are not imported into this quartic model. The residual argument below is
proved directly for the displayed generators.

## 1. Exact coordinates and the preserved subspace

Use unnormalized oscillator monomials

```text
|n>_B = product_j (a_j*)^n_j |0>,
<n|n>_B = product_j n_j!.
```

Creation increments an occupation, while annihilation decrements it and
multiplies the coefficient by the original occupation. At basis frequency
`Omega=2`,

```text
Q_j = (a_j+a_j*)/2,
P_j = i(a_j*-a_j).
```

The exact-gradient Hamiltonian, with bare mass and coupling `m=g=1`,
therefore has rational coefficients in these coordinates. The free control
`g=0` does also. Its quadratic and quartic monomials are formed on the full
Fock space before any projection. In particular, multiplying already
truncated coordinate matrices is not a way to construct `P H P`.

Let `P` retain the declared rectangular occupations. The unmeasured six
modes start even and remain even: the Hamiltonian contains only even powers
of every coordinate, and all kicks and reads act on `cos_x`. The restriction
to this parity sector is consequently exact for the specified preparation
and history. A uniform cutoff `d` uses
`d * ceil(d/2)^6` active coefficients instead of `d^7`; it does not discard
an additional physical tail. An odd spectator input would invalidate this
reduction and is outside this state-specific certificate.

All finite Hermite polynomials are in the domain of the confining polynomial
Hamiltonian and of `Q_cos_x`. The self-adjoint Hamiltonian is the Friedrichs
operator described in [the original derivation](DERIVATION.md). For `g=0`
it is the positive quadratic oscillator Hamiltonian. These facts are enough
for the finite-support residual argument. No convergence of an infinite
quartic Taylor series, and no analytic-vector hypothesis, is used.

## 2. A local curve certifies its own error

Take either generator `G=H`, optionally with a common scalar energy
subtracted, or a signed `G=Q_cos_x` for an auxiliary kick duration. Put

```text
A = P G P,
L = (I-P) G P.
```

A signed duration in the implementation is equivalently a nonnegative
duration with the sign absorbed into `G`. Centering the Hamiltonian uses
the same scalar on both sides of the comparison and removes the same
overall phase; it does not change the observation law.

The exterior action `L` is computed after the full polynomial action. One
Hamiltonian action has finite support outside the code, including the
quartic transitions. Its coefficients are retained explicitly. Auxiliary
coordinate evolution instead has its one-occupation boundary transition.

For the actual approximate input vector `v`, choose an integer degree `r`
and a nonnegative step duration `h`. Define a finite-code curve

```text
c_j = (-i A)^j v / j!,
p(s) = sum_(j=0..r) c_j s^j,       0 <= s <= h.
```

This is a finite polynomial curve in the operator domain, regardless of
whether `v` is normalized. Its full Schrödinger defect is

```text
i p'(s) - G p(s) = -L p(s) - A c_r s^r.
```

The two terms have orthogonal occupation supports. Hence the exact rational
integral of its squared norm is

```text
R2_full = sum_(j,k=0..r)
            Re <L c_j, L c_k> h^(j+k+1)/(j+k+1)
          + ||A c_r||^2 h^(2r+1)/(2r+1).

R2_finite = ||A c_r||^2 h^(2r+1)/(2r+1).
```

All inner products include the factorial Bargmann weights. Although some
cross terms can be negative, the displayed full sum is a squared-norm
integral and is nonnegative. The finite expression is the defect for the
**same curve** under the compressed generator `A`.

Variation of constants and unitarity give

```text
||exp(-i h G)v - p(h)|| <= integral_0^h ||i p'(s)-G p(s)|| ds
                       <= sqrt(h R2_full),

||exp(-i h A)v - p(h)|| <= sqrt(h R2_finite).
```

The second inequalities are Cauchy--Schwarz on the finite interval. They
are a posteriori residual bounds, not assertions that a Taylor series
converges. Polynomial degree and step duration change the checked residual;
an insufficient choice returns a large bound.

The endpoint is rounded to a declared dyadic coefficient grid. Compute the
norm of the **actual** difference `v_next-p(h)` with its factorial weights
and add an outward upper bound for that norm to both local budgets. This
accounts for every rounded or zeroed coefficient. The approximate state is
not renormalized.

For successive steps, compare each local curve to the unitary acting on
that curve's actual input. The incoming error passes through that unitary
unchanged, so local budgets add. A norm-growth factor for the polynomial
update is neither required nor appropriate. Integer square roots and exact
rational arithmetic enclose the reported budgets; floating trajectories
may be independent checks but are not inputs to these error bounds.

## 3. Preparation and both smooth instruments

The preparation `exp(+/- i Q/4)` and conditional feedback `exp(i Q/4)` are
coordinate evolutions covered by the same theorem. This compares the actual
finite functions of `Q_d=P Q P` with the unrestricted functions of `Q`.
It does not replace the finite operations by projected unrestricted ones.

For a read at `theta=pi/2` and `eta=1`, write

```text
C = cos(Q/2),      S = sin(Q/2),
K0 = (C-S)/sqrt(2),
K1 = (C+S)/sqrt(2).
```

Generate two coherent auxiliary vectors `u_plus=exp(i Q/2)v` and
`u_minus=exp(-i Q/2)v`. Then

```text
C v = (u_plus+u_minus)/2,
S v = (u_plus-u_minus)/(2i).
```

The transformation from `(u_plus,u_minus)` to `(C v-S v,C v+S v)` is
unitary on this direct sum. The physical outgoing vectors have the additional
common factor `1/sqrt(2)`. Thus if the two auxiliary errors are `e_plus`
and `e_minus`, the outgoing flagged-vector error is at most

```text
sqrt((e_plus^2+e_minus^2)/2).
```

The Weyl signs are coherent calculation branches within a Kraus operation;
they are not committed experimental records. Only the final `K0,K1`
outcomes are retained as pointer records. Computing with the raw `C-S,C+S`
vectors and an exact squared-norm weight `1/2` per read avoids introducing
an approximate `sqrt(2)`. After the two reads the raw terminal blocks carry
weight `1/4`.

Earlier records remain orthogonal flags. If raw block `b` has squared-norm
weight `w_b` and a local error budget `e_b`, simultaneous local errors
combine as `sqrt(sum_b w_b e_b^2)`. The conditional kick is applied only
to its designated first-record block, while the other block is carried
through unchanged. No conditional state is normalized and no low-probability
record is removed. The exact complete read is an isometry, and consequently
preserves the previously accumulated error bound. Dephasing its record to
make the actual committed classical output is a channel and cannot enlarge
trace distance.

## 4. What the final budgets mean

Let `Psi_app` be the final flagged approximate output, including both
records and the surviving field. The shared piecewise curves give separate
budgets

```text
||Psi_unrestricted - Psi_app|| <= E_unrestricted,
||Psi_finite       - Psi_app|| <= E_finite.
```

The approximate vector can have a norm different from one. Both ideal
physical histories start with the normalized product blank and use complete
instruments, so their output purifications are normalized. It follows that

```text
TV(ideal finite records, unrestricted records)
  <= trace_distance(ideal finite output, unrestricted output)
  <= min(1, E_unrestricted + E_finite).
```

The same bound holds after the common native continuation and after mixing
the initial native records. Taking the maximum over the two preparation
signs is sufficient. Keeping the exact native probabilities can give the
corresponding weighted bound. This is state-specific; it is not an
operator-norm estimate over arbitrary oscillator preparations or a claim
that a new physical clock has been measured.

Raw approximate probability tables have a different normalization issue.
For a normalized true output and vector error `E`, their half-L1 error is
bounded by `E+E^2/2`. A single event requires the safe bound `2E+E^2` unless
a sharper component calculation is supplied. The ideal-finite versus
unrestricted comparison above instead uses two normalized physical outputs
and needs no such trace-defect correction.

There is also a useful event bound when all approximate probabilities and
their total mass are retained. If the raw table is `q`, its mass is `m`,
and `half_L1(p,q)<=B` for a true normalized law `p`, then for every fixed
event `A`,

```text
|p(A)-q(A)| <= B + |1-m|/2.
```

To see this, split the signed difference `p-q` into its positive and
negative parts. Their sum is `2B` or less and their difference is `1-m`.
Each part is therefore at most `B+|1-m|/2`. This sharper bound uses the
actual complete-table mass; it does not normalize an individual branch.

To compare the already executed finite-grid source with unrestricted
exact-gradient observations, first use the existing temporal and spatial
budgets at that source's occupation code, then add its occupation bound.
Floating compiler error and physical service error remain distinct. A
larger certified occupation reference is not automatically a new source
execution. It can be related to the retained execution only by adding a
separately justified comparison between those codes.

## 5. The retained discriminator uses unrestricted terminal effects

The principal interaction discriminator uses the positive preparation
`exp(i Q/4)`, both Hamiltonian intervals `1/8`, the first complete smooth
read and its conditional `exp(i Q/4)` kick. It retains the same four final
record pairs as the original protocol. The auxiliary interacting code has
dimensions

```text
(6, 12, 6, 6, 6, 6, 6),
```

in mode order `(zero,cos_x,sin_x,cos_y,sin_y,cos_z,sin_z)`. The measured
coordinate receives twelve levels; the other modes receive six. Their even
parities reduce the active field sector to `12*3^6=8748` coefficients before
record flags. The finite-support residual construction runs through the
second evolution and then stops **before** approximating the last read.
Every first-record branch is still present and has its proper weight. Its
budget `E` bounds the entire pre-final flagged vector against the original
unrestricted interacting target.

On that approximate vector evaluate the actual unrestricted final effects,

```text
E0 = (I-sin(Q))/2,
E1 = (I+sin(Q))/2.
```

This avoids a final finite-coordinate function `sin(Q_d)` and its occupation
error. It computes probabilities of the declared detector on the approximate
pre-final state; it does not reset the field or change the target detector.
The surviving infinite-support final quantum vectors are not emitted by
this probability calculation. Its certificate is expressly for the
classical observation law. The separate four-level full-history comparison
continues to cover the actual finite preparation, both finite instruments,
feedback and retained terminal field.

### Exact matrix elements of the final detector

In the same Bargmann coordinates, the canonical commutator gives

```text
exp(i Q) = exp(-1/8) exp(i a*/2) exp(i a/2).
```

For nonnegative occupations `m,n`, its unnormalized matrix element is

```text
<m|exp(i Q)|n>_B = exp(-1/8) sum_(l=max(0,n-m)..n)
    m! n! (i/2)^(m-n+2l) / [l! (n-l)! (m-n+l)!].
```

All coefficients multiplying the common Gaussian factor are rational real
or imaginary numbers. Spectator occupations must match on the two sides;
their factorial norms supply the remaining factors. Therefore each of the
four raw polynomial probabilities has the form

```text
q_(first,last) = A_(first,last) + C_(first,last) exp(-1/8),
```

with exact rational `A,C`. The imaginary part of the Weyl expectation
supplies `sin(Q)`. The degree-32 alternating exponential polynomial is an
upper bound for `exp(-1/8)` and its degree-33 successor is a lower bound.
Their exact rational interval is propagated with the appropriate sign of
`C`. No floating matrix exponential or numerical quadrature supplies these
probability bounds.

Completeness `E0+E1=I` makes the raw table mass exactly the pre-final
polynomial squared norm `m`; the two exponential coefficients cancel in
each first-record block. Applying the complete unrestricted read to the
approximate and true pre-final vectors gives

```text
half_L1(true record law, raw polynomial record law) <= B,
B = E + E^2/2.
```

The fixed event `A={01,10}`, meaning that the two field records differ, is
then enclosed using

```text
|Pr_true(A)-q(A)| <= B + |1-m|/2.
```

The exponential enclosure bounds `q(A)`; the residual and actual mass defect
give the additional radius. Retaining all records is essential to this
mass-aware improvement over the generic single-vector event bound. There
is no normalization of a selected branch, and the raw probabilities remain
raw in the receipt.

### The free comparison traces out spectators exactly

With `g=0`, the original seven-mode Hamiltonian is a tensor sum. The product
blank and the preparation, both reads and conditional kick factor between
the observed `cos_x` mode and the six unobserved modes. Every complete
record amplitude consequently has the same unit-norm spectator factor;
tracing it out leaves its record probability unchanged. This is exact
factorization, not a low-occupation approximation to those spectators.

The observed free Hamiltonian is `P^2/2+Q^2`, with physical oscillator
frequency `sqrt(2)` in the supplied model coordinate. In the fixed
`Omega=2` basis, subtracting its common blank energy `3/4` gives

```text
H_observed - 3/4 = (3/2) N - (a^2+(a*)^2)/4.
```

The same residual method uses twelve occupation levels for this single
mode, stops before the last read, and evaluates the same unrestricted
terminal effects. Its error bound concerns the exact original seven-mode
**record law**, not the seven-mode quantum output state. It does not claim
that the six spectator wavefunctions remain at their initial blanks.
Factorization would fail for a correlated preparation or an apparatus
acting jointly on a spectator, neither of which is used here.

### A fixed event in the full native history

The interacting and free calculations enclose the same event with the same
positive preparation, physical smearing, feedback and model times. Let
`[L1,U1]` and `[L0,U0]` be those enclosures. Then

```text
Pr_g=1(bits differ) - Pr_g=0(bits differ) >= L1-U0.
```

A strictly positive lower bound also bounds the conditional complete-history
TV distance from below. In the full native experiment the event additionally
requires a positive first native current, which selects precisely this
preparation. That current has exact probability `5/66` under the retained
uniform native initial state and repair coin. The unconditional event
difference, and hence the complete native-history TV distance, is at least
`(5/66)*(L1-U0)`. All other native and field outcomes remain part of the
law; this event calculation is not postselection of the experiment.

The event was selected from the existing internal finite-model receipt
before the new certificate. No experimental data selected it or entered
either calculation. The positive result is a discriminator between the
supplied interacting action and its free control in the unrestricted
seven-mode model, not empirical confirmation of that action.

## 6. Optional adjoint refinement for a record event

A complete-output norm bound also resolves high-occupation changes in the
terminal field that a chosen classical record event may not detect. To
certify that event more accurately, use the same forward polynomial and a
separately checked reverse history. This is an additional observable-specific
estimate, not a reclassification of the complete-output error.

Keep the two pointer registers from the beginning, initialized to `00`.
Extend every read to the full unitary on its pointer,

```text
R = [[K0, -K1], [K1, K0]].
```

The conditional kick is controlled by the first pointer. This deferred-read
unitary has the same complete classical law as the committed experiment:
later operations neither mix nor erase the earlier record. Let `U` denote
that unrestricted complete unitary, `v` the approximate terminal vector,
and `M` the projector onto the selected event in the two pointer registers.
The event is evaluated without normalizing or discarding another outcome.
The occupation projection `P` acts only on the field, so it commutes with
`M`. Put

```text
w = M v,
q = ||w||^2,
delta = U psi_initial - v,
||delta|| <= E.
```

The event probability identity gives

```text
|Pr_unrestricted(event)-q|
  <= 2 |<w,delta>| + E^2.
```

Represent `delta` by the forward curve's propagated defects and endpoint
jumps. For a cut at auxiliary time `s`, its contribution pairs the forward
defect with the exact reverse vector `b(s)=U(T,s)* w`. Unitarity gives
`||b(s)||=||w||`. Certify a reverse finite-support polynomial curve
`b_app(s)` with a bound `B` valid at **every** reverse prefix. Since this
curve lies in `P`, it is orthogonal to the exterior forward defect. Write
`D_inside` for the sum of integrated interior defect norms and all actual
endpoint jumps, and `D_exterior` for the integrated exterior defect norms.
Then

```text
|<w,delta>| <= sqrt(q) D_inside + B D_exterior,

|Pr_unrestricted(event)-q|
  <= 2 sqrt(q) D_inside + 2 B D_exterior + E^2.
```

The interior term is bounded against the exact reverse vector, whose norm
is known; no normalization of the approximate reverse curve is assumed.
The exterior term has two residual factors because an in-code reverse
curve pairs to zero with a forward occupation leak. Replacing
`D_exterior` by an upper bound on the integrated **full** forward defect
and jumps is conservative and is what the implementation does.

There are two bookkeeping distinctions from the endpoint proof. First,
the two Weyl curves within a read are simultaneous coherent alternatives.
Their endpoint errors combine in Euclidean norm, but the integral of their
combined defect norm is bounded by the sum of the two individual integral
bounds. The code therefore inflates a read's retained endpoint budget by
`sqrt(2)` when constructing integrated forward defect budgets. Other
stages already have integrated-residual bounds. The same inflation is
applied to the interior budget and its rounding terms.

Second, reversing a read acts on an already present pointer with **both**
inputs. The implementation uses

```text
R* = exp(+i Q tensor Y / 2) exp(+i pi Y / 4).
```

The constant inverse pointer rotation maps `(v0,v1)` to
`(v0+v1,-v0+v1)/sqrt(2)`: exact integer sums and a halved squared state
weight implement it. The subsequent generator `Q tensor Y` also has
rational real and imaginary Bargmann coefficients and is self-adjoint.
The same residual calculation therefore checks the entire inverse read in
one flow, on all existing pointer inputs.

An independent equivalent check recombines the two Weyl evolutions on the
existing pointer. Its raw recombination takes four Weyl/pointer blocks to
two pointer blocks and has operator norm `sqrt(2)`. If the old squared
state weight is `w_old`, a safe error for that alternative is
`sqrt(w_old*(e_plus^2+e_minus^2))`, with separately weighted actual rounding.
The smaller fresh-pointer forward factor cannot be reused for this
inverse unitary. Tests check both this norm and the full inverse circuit
against an independent dense calculation.

Each local residual integral also bounds its restriction to a partial
step. Sums of these nonnegative bounds therefore control all reverse
prefixes, including intermediate points within an evolution or read.
Reverse operations are the inverse second read, inverse second evolution,
inverse first-record kick, inverse first read, inverse first evolution and
inverse preparation. Both pointer registers survive all these operations.

The principal event is `first_bit != last_bit`, containing records `01`
and `10`. It was selected from the already retained finite-model response,
without experimental data. For each supplied action the probability
enclosure is clipped to `[0,1]` and may be intersected with the ordinary
forward bound `q +/- (2 sqrt(q) E+E^2)`. Separation between the `g=1` and
`g=0` intervals establishes a discriminator in the unrestricted
seven-mode model; it is not a fit to or prediction validated by nature.

This adjoint method is implemented and independently tested, but is not
the method producing the main retained discriminator. That discriminator
uses the pre-final curves and exact terminal effects in Section 5. No
production adjoint receipt or larger source execution is implied.

## 7. Interacting observations and remaining physical scope

An interaction discriminator must compare the same event or complete
history with `g=1` and `g=0`. Each action needs its own occupation budget.
The exact-gradient residual calculation covers both; it never divides by
the coupling. A numerical interaction-deletion contrast, less the two
applicable error budgets, establishes a positive unrestricted contrast only
when the resulting lower bound is positive. A cutoff bound larger than the
contrast cannot certify that discriminator.

This calculation addresses the observation error caused by finite
occupation, including the detector and feedback. It does not select the
supplied scalar law from native repair, enlarge the spatial Fourier band,
establish relativistic locality, calibrate an observed interacting mass or
coupling, or compare against an independent experimental dataset. The
original receipt remains frozen with its original finite-code scope; the
additional occupation certificate supplies a new comparison when its
declared settings and bounds apply.
