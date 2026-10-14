# Calibrating axial readouts and predicting the remaining observations

This component uses the published zero-simulated-duration characteristic-function
tables from [McGarry et al.](https://arxiv.org/abs/2603.04744) and their
[data deposit](https://doi.org/10.5281/zenodo.18752287). The paper identifies the
real and imaginary entries as spin-measurement means. It reports a displaced,
nearly ground-state preparation, 250 measurements per two-dimensional scan
point, and intermediate-outcome selection. Its zero-duration table includes
values reconstructed by symmetry. These reconstructed values are not additional
trials.

The analysis rule selects a representative quadrant by array indices and
assigns its two axes to calibration and its interior to validation. Neither an
off-axis residual nor a fitted displacement chooses this split. This is a
retrospective test on published observations, not prospective experimental
preregistration. Symmetry copies are identified and excluded from the trial
count. The rule was fixed before this implementation inspected the numeric
table entries; ambiguity in their absolute coordinate convention remains.

## Characteristic function and the model readout

Let `[x,p]=i`, `a=(x+i p)/sqrt(2)`, and `beta=u+i v`. Direct substitution gives

```text
D(beta) = exp(beta a† - beta* a)
        = exp(i sqrt(2) (v x - u p)),
chi(beta) = Tr(rho D(beta)).
```

The scalar packet uses `Q=(a+a†)/sqrt(2 Omega)` and
`P=-i sqrt(Omega/2)(a-a†)`, so `x=sqrt(Omega) Q` and
`p=P/sqrt(Omega)`. Hence

```text
D(u+i v) = exp(i [sqrt(2 Omega) v Q - sqrt(2/Omega) u P]).
```

For `r=sqrt(u^2+v^2)>0`, define
`Q_phi=Q cos(phi)+(P/Omega) sin(phi)`, with
`cos(phi)=v/r`, `sin(phi)=-u/r`, and `eta=sqrt(2 Omega) r`.
Then `D(beta)=exp(i eta Q_phi)`. The smooth two-outcome readout

```text
K0 = cos((eta Q_phi + theta)/2),
K1 = sin((eta Q_phi + theta)/2)
```

has signed output `s=+1` for bit zero and `s=-1` for bit one, with
`E[s]=Re(exp(i theta) chi(beta))`. Phase `theta=0` reads `Re chi`;
phase `theta=-pi/2` reads `Im chi`. The principal scalar history's
`theta=+pi/2` instead reads `-Im chi`; its bit convention must not be silently
reused with the opposite sign. At `beta=0`, `chi=1` exactly.

This is an algebraic identification of the readout family. The principal
seven-mode receipt executes a fixed `Q` read and a different prepared state.
It does not already execute all these rotated quadratures on the experimental
state. No occupation or interaction certificate is transferred to this
component by the identification.

## Axial calibration implies an off-axis Gaussian prediction

Write `mu=(mu_x,mu_p)` and let `Sigma` be the symmetrized covariance matrix:
`Sigma_jk=Tr(rho {R_j-mu_j,R_k-mu_k})/2`, for `R=(x,p)`.
For a Gaussian state, its definition through the Weyl characteristic function
gives

```text
chi(u+i v) = exp(i sqrt(2) (v mu_x-u mu_p)
                -Sigma_xx v^2-Sigma_pp u^2+2 Sigma_xp u v).
```

For a physical single-mode Gaussian, `Sigma` is positive and
`det(Sigma)>=1/4`. Arbitrary displacement and unequal quadrature widths are
allowed. When `Sigma_xp=0`, the exponent separates, proving

```text
chi(u+i v) = chi(u) chi(i v).                         (1)
```

No displacement, temperature, width or detector gain is fitted from validation
entries. The measured complex values on the two axes supply the prediction.
For a Gaussian with nonzero covariance the right side is instead multiplied
by `exp(2 Sigma_xp u v)`. Equation (1) is a state-dependent identity between
expectations; it does not assert that the two Weyl operators commute. Indeed,
`D(u) D(i v)=exp(-i u v) D(u+i v)`.

Agreement with (1) at finitely many points neither proves Gaussianity nor
uniquely identifies a quantum state. The same relation holds for a classical
Gaussian phase-space distribution with independent coordinates. This test
addresses the declared preparation/readout description, not a distinction
between OPH and every competing theory.

The separately fixed nominal prediction assumes a displaced thermal state
with `nbar=0.04`, `mu_x=-1.5`, and `mu_p=0`, using the paper's preparation
values rather than fitting this table:

```text
chi_nominal(u+i v) = exp(-0.54 (u^2+v^2)-i sqrt(2) 1.5 v).
```

This formula additionally requires a verified numerical displacement grid and
its row/column convention. Equation (1) survives an interchange of the two
axes, separate nonzero rescalings, and either axis reversal; its nominal
formula does not. A conflict between source descriptions of the coordinate
range must therefore remain visible and must not be resolved by choosing the
range that fits better. A coordinate-convention comparison is a diagnostic,
not a fitted physical calibration. Both advertised coordinate extents and
both row-sign conventions are evaluated as separate nominal diagnostics;
the best-fitting convention is not selected as an inferred physical fact.

## Propagating complex calibration uncertainty

For each retained, independently measured scalar component let `z_j` be its
reported spin mean and `n_j` its accepted sample count. The outcomes lie in
`[-1,1]`. With `M` distinct measured components and familywise level `delta`,
Hoeffding's bound and a union bound give simultaneous half-widths

```text
r_j = sqrt(2 log(2 M/delta)/n_j),
I_j = [max(-1,z_j-r_j), min(1,z_j+r_j)].               (2)
```

The implementation starts near this analytic width, then certifies the actual
floating value independently of the logarithm and square root. Interpret that
value as its exact rational number `r`, set `x=n r^2/2`, and form
`S=sum_{k=0}^{128} x^k/k!` with rational arithmetic. Since `x>=0`,
`exp(x)>=S`. Round `S` down to 64 fractional binary places to obtain `L`;
then

```text
Pr(any component misses its interval) <= min(1, 2 M/L).
```

The retained certificate includes `r`, `x`, `L` and this final rational bound.
The width is increased until the bound meets the assigned failure allowance.
The finite positive series is only a lower bound on the exponential; no
unproved truncation error or floating transcendental evaluation supplies the
coverage claim. Mean-interval endpoints are rounded outward as well.

The factor in (2) is for a spin mean, not a frequency in `[0,1]`. There are
separate real and imaginary components. Reconstructed mirror entries, reused
axis values and exact identities at the origin do not create new independent
samples. Independence is required within each sample batch. The union bound
does not require independent errors between different components.

For two complex calibration boxes `A=R_A+i I_A` and `B=R_B+i I_B`, ordinary
real interval arithmetic supplies the enclosing prediction box

```text
Re(A B) subset R_A R_B - I_A I_B,
Im(A B) subset R_A I_B + I_A R_B.                     (3)
```

Each real interval product takes the minimum and maximum of its four endpoint
products. Every product and subsequent sum or difference is rounded outward.
The implementation retains the enclosing product box; its components may
also be intersected with `[-1,1]` without losing a characteristic function.
Ignoring dependencies between repeated interval factors enlarges this
enclosure and does not invalidate it. Reuse of one axis in many predictions
needs no new sampling assumption: all predictions follow deterministically on the event
that their original calibration boxes cover their means.

Use (3) to forecast every validation point through (1). Reject the declared
factorization when a predicted component interval and its measured validation
interval are disjoint. If (1) and the sampling premises hold, simultaneous
coverage from (2) makes the probability of any such rejection at most
`delta`. Equivalently, separately allocated calibration and validation levels
give the bound `delta_cal+delta_val`. Failure to reject is compatibility at
this resolution, not acceptance with probability `1-delta`.

The published nominal count of 250 measurements per point does not itself
supply accepted-shot counts after intermediate-outcome selection. A receipt
using `n_j=250` must label its coverage conditional on that denominator and
independent-sampling interpretation. The aggregate tables do not establish
them. Selected-shot means also do not recover the rejected-shot record law.
Unbounded apparatus bias, correlated drift, background correction or extra
normalization is not covered by projection-noise intervals. A quoted
statistical interval is consequently not a complete physical error bound.

## What this attachment identifies

The axial observations calibrate complex readout responses in the apparatus's
declared coordinates. Off-axis values test a withheld relation among those
responses. They do not independently establish an absolute field scale,
prepare the scalar packet's product blank, or identify its seven-mode action.
The `0 ms` label denotes the experiment's simulated evolution duration;
preparation and measurement still have physical durations. Since this
comparison contains no intervening scalar evolution, it cannot infer the
scalar coupling `g`, its mass, a propagation speed, or a conversion from
Hamiltonian model time to seconds. It is component evidence from physical
readout observations, with these premises and limits retained explicitly.
