# A finite executed environment for number dephasing

This construction replaces the preceding benchmark's directly applied analytic
dephasing matrix with a positive finite mixture of executed source words. The
Gaussian phase-noise family, its rate, the independent environment refreshes
and the duration parameter are supplied apparatus. The construction does not
infer those premises from experimental outcomes or derive them from OPH.

## Channel and execution

For oscillator occupation dimension `d`, let

\[
N=\operatorname{diag}(0,1,\ldots,d-1),\qquad
\sigma^2=\gamma t.
\]

Here `t` is the supplied simulated Schrodinger-frame duration in seconds and
`gamma` is the supplied rate per second. The specified noise channel is

\[
\mathcal D_t(\rho)
=\mathbb E_{Z\sim\mathcal N(0,1)}
 [e^{-i\sigma ZN}\rho e^{i\sigma ZN}],
\quad
(\mathcal D_t(\rho))_{nm}
=e^{-\gamma t(n-m)^2/2}\rho_{nm}.
\]

Choose the standard-normal `M`-point Gaussian quadrature rule, with nodes
`z_j` and strictly positive weights `w_j` summing to one. The implemented
finite channel is

\[
\mathcal Q_{M,t}(\rho)
=\sum_{j=1}^{M}w_j U_j\rho U_j^\dagger,
\qquad U_j=e^{-i\sigma z_jN}.
\]

This is CPTP: its Kraus operators are `sqrt(w_j) U_j` and their completeness
sum is the identity. The environment draw is discarded without conditioning;
spin and other reference registers remain coherent. Each invocation requires
an independently refreshed draw. Reusing one draw throughout the whole pulse
sequence would model a different, correlated environment.

For `d=2^b`, little-endian occupation satisfies

\[
N=\sum_{r=0}^{b-1}2^r(I-Z_r)/2.
\]

Consequently each branch is a chronological native word with one identity
phase of angle `-sigma*z_j*(d-1)/2` and `b` native `Z` phases of angles
`sigma*z_j*2^r/2`. Every phase is executed by the existing twelve-port source
compiler. All global phases remain in the target checks. A branch needs
exactly `b+1` native pulses, including zero-angle requests. At `d=64`, this
means 7 pulses per environment realization and 112 pulses to enumerate all
16 branches. The former is the maximum and expected physical branch count;
the latter counts the deterministic all-branch contraction used to calculate
the exact finite ensemble. These are different resource questions.

The host caches the contracted branch matrices and assembles the Schur matrix
from their executed diagonal phases. Applying that matrix is algebraically
identical to summing all the executed Kraus branches. The analytic Gaussian
Schur matrix is used only for independent tests, never as the implemented
prediction channel.

## Uniform finite-rule error

Let `He_M` be the monic probabilists' Hermite polynomial. Its roots are the
quadrature nodes and

\[
\mathbb E[\operatorname{He}_M(Z)^2]=M!.
\]

For `f(z)=exp(-iaz)`, interpolate both its value and derivative at every node
by a degree-at-most-`2M-1` polynomial `H`. The integral divided-difference
remainder, which also applies to this complex-valued function, gives

\[
|f(z)-H(z)|\le
\frac{|a|^{2M}}{(2M)!}\operatorname{He}_M(z)^2.
\]

Gaussian quadrature integrates `H` exactly and `H(z_j)=f(z_j)`. Integration
therefore proves the entrywise bound

\[
\left|\mathbb E e^{-iaZ}-\sum_jw_je^{-iaz_j}\right|
\le\frac{M!}{(2M)!}|a|^{2M}.
\]

Taking `a=sqrt(gamma*t)*(n-m)` and `|n-m|<=d-1` gives

\[
\epsilon_{\rm entry}=
\frac{M!}{(2M)!}[\gamma t(d-1)^2]^M.
\]

This entry bound controls the complete channel including an arbitrarily
large untouched reference. Write a density operator in occupation blocks
`rho_nm` acting on the reference, and let `p_n=tr(rho_nn)`. Positivity gives
`||rho_nm||_1 <= sqrt(p_n*p_m)`. The triangle inequality and Cauchy--Schwarz
then give

\[
\begin{aligned}
\|[(\mathcal Q_{M,t}-\mathcal D_t)\otimes I](\rho)\|_1
&\le\epsilon_{\rm entry}\sum_{n,m}\|\rho_{nm}\|_1\\
&\le\epsilon_{\rm entry}\left(\sum_n\sqrt{p_n}\right)^2\\
&\le d\epsilon_{\rm entry}.
\end{aligned}
\]

Thus the uniform trace-distance, or half-diamond, bound is

\[
\boxed{\epsilon_{\rm env}
=\min\left(1,\frac{dM!}{2(2M)!}[\gamma t(d-1)^2]^M\right).}
\]

The dimension factor is the oscillator dimension, not the oscillator times
spin dimension. Any subsequent CPTP operations and joint measurements cannot
increase this error. Independent environment invocations telescope by summing
their radii; conditional selected detector laws require a separate acceptance
bound. The environment radius alone does not certify splitting error, an
infinite occupation cutoff, experimental calibration or physical service
accuracy.

For the declared `d=64`, `M=16`, `gamma=18/s`, `t=25 microseconds`, the
per-invocation ideal bound is approximately `2.72833e-17`. This very small
quadrature bound is below ordinary floating-point execution uncertainty; it
must not be reported as a total numerical or physical accuracy certificate.

The Gaussian-rule construction and remainder are standard numerical
analysis; see the [NIST DLMF Gaussian quadrature discussion](https://dlmf.nist.gov/3.5#v)
and [Jacobi-matrix characterization](https://dlmf.nist.gov/3.5#vi). The
additional operational step here is compiling each finite environment phase
into this model's source vocabulary and propagating the quadrature error
through reference-entangled states.

## Arithmetic and evidence boundaries

The Jacobi matrix has off-diagonal entries `sqrt(j)`. Its eigenvalues are the
nodes; squared first eigenvector components are the weights. Both are
computed in the existing private 80-decimal-digit arithmetic context,
independent of the platform's binary64 eigensolver. Requested native angles
are rounded to binary64 at the compiler boundary. All source words retain
their hashes, resource counts and independent full-phase matrix checks.

An exactly diagonal native word can acquire off-diagonal roundoff in the
existing polar stabilization. The channel explicitly projects back to its
known diagonal unitary structure, rejects changes above `2e-12`, and retains
the measured change. Classical branch weights receive their disclosed
floating normalization. These adjustments are numerical representations of
the specified finite channel, not discarded measurement outcomes.

`quadrature_bound` evaluates the ideal radius in exact rational arithmetic
from the input binary64 parameter values, stores the exact fraction, and
rounds its displayed positive floating bound upward. This certifies the
ideal quadrature expression. It does not certify errors in the computed
nodes, weights, source matrices or subsequent floating operations. Those
remain separately identified diagnostics. Tests check positive moments,
high-precision characteristic-function errors, arbitrary reference-register
contraction, completeness, entangled-state behavior and source-word counts.
