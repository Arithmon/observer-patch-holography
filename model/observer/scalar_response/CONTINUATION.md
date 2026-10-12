# Two reads of the retained field

The scalar observer now continues the field after its first observation.
One native source and one scalar field are prepared per history. Their
states persist through three native repairs and two field reads. Every
outcome and relative Weyl phase is retained; there is no intermediate reset.

```sh
python3 -m model.observer.scalar_response --continue-field \
  --check model/observer/scalar_response/continuation_results.json
python3 -m pytest -q -W error model/observer/scalar_response
```

## Fixed process and memory

The preparation, action, gamma=4, profiles and four grids are inherited from
the [one-read experiment](README.md). Fix the two model times to `19/20` and
`1`. The settings and 50,000-trial seed `20261012` are declared before the
execution; no outcomes select a grid, time or coefficient.

1. Prepare the native binary total-two shell uniformly and repair `(0,1)`.
   Keep its coin and current j, then use j for the single scalar preparation.
2. At `19/20`, apply the existing centered field instrument and commit y1.
3. Use only y1 to choose the second native seam: `(0,5)` for +1, `(0,1)`
   for -1. Repair the same loads and keep its coin and current j2.
4. Evolve the actual conditional field to `1`, apply the same instrument
   and commit y2. The previous read is not replaced by its unconditioned
   marginal. Use only y2 to choose the third native seam and retain j3.
5. Halt with all five records readable. Halting changes neither source nor
   records. A fresh run requires another explicitly counted preparation.

All 2,112 labelled histories are evaluated. The record-prefix counts at the
six phases are `1,4,8,24,48,124`: 209 reachable prefixes have an implemented
eight-bit encoding, and 124 terminal words need seven fixed-width bits.
This counts the reached controller records, excluding the program, source,
field, preparation apparatus and their service costs. Allocating every
well-typed record sector instead gives 1,099 tape states and 85,722
source-times-tape dimensions before the field. The field Hilbert space is
infinite-dimensional; its four-term representation is not a four-state
physical memory. No cosmological capacity or saturation follows.

## Full-mode state with a small representation

Let `v_a=exp(i Omega t_a) Omega^(-1/2)g/sqrt(2)` and
`alpha=i j gamma Omega^(-1/2)f/sqrt(2)`. Define

```text
G_ab = <i v_a, i v_b>,    m_a = 2 Im <alpha, i v_a>,
X_a = a^dagger v_a + v_a^dagger a,
K_y(a) = [exp(-i X_a/2)-i y exp(i X_a/2)]/2.
```

The kernel accumulates every field mode into G and m using 2D slices.
At q=233 all 12,487,168 modes contribute. A state is stored as
`sum_n c_n D(i sum_a h_na v_a)|alpha>`, where each h is an exact half-integer.
Its Gram matrix is

```text
<h|k> = exp[-(k-h)^T Re(G)(k-h)/2
              + i h^T Im(G) k + i (k-h)^T m].
```

Applying a displacement `D(i s v_a)` changes h to `h+s e_a` and multiplies
its coefficient by `exp[-i s Im(G_a,: h)]`. Both terms of each K are
applied before taking the branch norm. Thus the first read has two terms
and the second has four. This preserves conditional quantum states, not
just their first-read probabilities or Gaussian means. Checkpoints retain
the coefficients, shifts and an execution-kernel binding; restoring them
preserves both subsequent outcomes. A global free vacuum-energy phase
cancels from all states as rays and from all probabilities.

The generic evaluator permits at most seven declared reads and 128 terms;
this protocol uses two and four. It rejects unresolved normalization and
nonfinite arithmetic. Floating evaluation is not an interval certificate.
Independent oscillator-basis tests check complete final vectors, not only
outcome probabilities: reversing all symplectic signs can leave the latter
unchanged. Full-mode coherent-state tests supply another independent check.

## Exact joint law and native continuation

For the same detector at both times put `V=G_11=G_22`, `C=Re G_12`,
`c=2 Im G_12`, so `[X1,X2]=i c`. Write `m_a=j mu_a` and

```text
S_a(j) = exp(-V/2) sin(j mu_a),       D = cos(c/2),
R_j = [exp(-(V-C)) cos(j(mu1-mu2))
       - exp(-(V+C)) cos(j(mu1+mu2))]/2.

p_j(y1,y2) = [1+y1 S_1(j)+y2 D S_2(j)+y1 y2 R_j]/4.
```

Summing the first outcomes conjugates the second effect by the two half
displacements with equal weights, giving the factor D. Weighting those
outcomes gives the centered sandwich
`[cos(X1-X2)-cos(X1+X2)]/2`; the Weyl phases cancel in this expression.
The coherent characteristic function gives R. The implementation evaluates
the conditional states directly; this closed law is an independent check.

Let `A=S_1(1)`, `B=D S_2(1)` and `R0,R1` be the zero/nonzero-current
correlations. Each native history has weight `p_j(y1,y2)/528`. Exact
rational enumeration, independently checked by endpoint averaging, gives

```text
E[j2] = 5 A/264,
E[j3] = 5(B-A)/352 + (R0-R1)/264,
E[N0_final] = 1/6 - 5 A/1056 - 5 B/352 - (R0-R1)/264.
```

At q=233, numerical execution gives `E[j3]` about `-0.00130126` and
`E[y1 y2 | j=0]` about `0.0514838`. The earlier field record improves
prediction of the second bit: deleting it while leaving the experiment
unchanged raises the model's Brier risk by about `0.00064677`. The larger
native-history predictive gain also includes native-source memory and is
not attributed solely to field backaction. Forecasts are calculated from
the declared law before sampling, without hidden-state access or learning.

## Controls distinguish correlation from measurement disturbance

- **Conditional-field donor:** continue an independently drawn first-read
  branch from the same preparation while preserving the actual native
  state and records. This preserves both field marginals and the averaged
  first-read disturbance, but removes the association with the actual y1.
- **Field reset:** explicitly replace the conditional field by the original
  prepared, freely evolved field. The two outcome law then factors given j.
- **Commuting Gaussian control:** retain the same means and real covariance
  but set the commutator to zero. It retains R, including the vacuum record
  correlation, while removing D. Therefore vacuum correlation alone does
  not distinguish quantum backaction from this classical alternative.
- **Independent bit and disconnected actuator:** respectively replace the
  second bit by a fair donor or fix the third seam. Forecast-only record
  deletion and matched record shuffling leave the actual process unchanged.

In particular, on the positive-current preparation,

```text
p_second(+ | first read performed) - p_second(+ | first read omitted)
    = (-r(t2)) [1-cos(c/2)].
```

The rational certificate proves this is greater than `0.000022281162` at
q=233, and greater than `0.0000252366` in the declared compact continuum.
It also proves the distinct vacuum correlation is greater than
`0.045619388905` at q=233. These are properties of the specified instrument
and field. The 50,000-trial software samples test the execution but do not
statistically resolve the small disturbance contrast. No physical data or
universal exclusion of classical hidden-state models is claimed.

## A bound for the complete continued history

The single-read probability certificate cannot be reused as a history
bound: j=0 has a fair single bit but a covariance-dependent two-bit law.
Instead reuse the parent's full one-particle residuals and compact tails.
All quantities below are compared through its same smooth nine-mode
isometry; that comparison is not a truncation of the executed field.

Let b be the smooth detector-vector norm, s the sampling-isometry error,
κ_f and κ_v the parent's full frequency residuals, and θ_h,g, θ_h,f and
θ_infinity the separate graph and continuum tail masses. Set

```text
d_b = (κ_f/2+s)/sqrt(2),    eta = κ_v/sqrt(2),
e_g = d_b+sqrt(θ_h,g/2),    e_f = d_b+sqrt(θ_h,f/2),
e_inf = sqrt(θ_infinity/2).
```

The parent's residual bound depends only on absolute mode coefficients,
so it also applies to the reflected detector. Its negative-power estimate
and finite-band Duhamel formula give detector error `e_g+t eta` and initial
preparation error `gamma e_f`. The corresponding continuum errors are
`e_inf` and `gamma e_inf`. Stationarity lets the cross kernel use only
`Delta=t2-t1=1/20`. Hence graph-to-band errors are

```text
v_h = e_g(2b+e_g),
g_h = b(2e_g+Delta eta)+e_g(e_g+Delta eta),
v_inf = g_inf = e_inf(2b+e_inf),
E_V = v_h+v_inf,           E_G = g_h+g_inf,
E_mu(t) = 2[b gamma e_f + gamma b(e_g+t eta)
             + gamma e_f(e_g+t eta)
             + 2 gamma b e_inf + gamma e_inf^2].
```

These bound variance, complex cross-kernel and mean differences, respectively.
Because `V±C>=0`, the exponential-cosine terms give
`E_R0<=E_V+E_G` and `E_R1<=E_V+E_G+E_mu(t1)+E_mu(t2)`.
The first sine-moment error is bounded by `E_V/2+E_mu(t1)`, sharpened
when useful by twice the parent's probability error. The second has the
same bound plus the attenuation error
`E_D<=max(|c_h|,|c_inf|) E_G/2`. Each moment difference can be capped at 2.

For signed differences a,b,d in the three two-bit Walsh moments, direct
summation over the four atoms gives the exact identity

```text
TV = max(|a|, |b|, |d|, (|a|+|b|+|d|)/2)/2.
```

For j=0 only R differs. Native first-current probabilities are `28/33` for
zero and `5/66` for each sign. Weight the corresponding conditional TV
bounds by these probabilities. Subsequent native updates use the same
stochastic rule, so this bounds the complete labelled history, including
initial source and all coins; forgetting labels cannot increase TV.
The independent directed rational replay yields

| Grid q | Complete-history TV upper bound |
| --- | ---: |
| 13 | 1 |
| 34 | < 0.451806 |
| 89 | 0.052189606637 |
| 233 | 0.025709913643 |

The coarse q=13 and q=34 enclosures do not certify either witness's sign.
Their receipts retain the universal signed lower bound -1; zero would
incorrectly assert nonnegativity. The positive witness bounds above concern
q=233 and the continuum, with q=89 also certified positive.

The same enclosures certify the controls above: `R0=exp(-V)sinh C` is at
least `exp(-V) C` when C>0. For `0<c<1`,
`1-cos(c/2)>=c^2/8-c^4/384`; the parent's strictly negative r(t2) supplies
the positive disturbance bound. No floating execution value enters either
certificate. Replays reject changes to source pins, settings and bounds.

The fixed nine-mode tail bound does not vanish with q. The parent's separate
finite-Weyl-product convergence theorem, allowing enlarged comparison
spaces, supplies the asymptotic limit. These finite error numbers do not
prove a rate or full-state norm convergence. They also exclude numerical
roundoff, apparatus error and physical clock calibration.

The completed addition is an executed retained-field observer and its
complete-history bound. The supplied scalar action, quantization, vacuum,
apparatus and model time keep their original status; no microscopic
derivation or physical identification is inferred from this execution.
