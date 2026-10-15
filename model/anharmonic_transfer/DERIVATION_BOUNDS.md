# Finite channel execution and selected-readout bounds

These bounds compare two ideal mathematical executions at the **same finite
occupation dimension**, with identical preparation, apparatus and detector
instruments. The reference evolves the supplied finite Lindblad generator
without splitting; the implemented execution uses native unitary words,
symmetric channel splitting and finite environment mixtures. The result bounds
every retained detector probability, including after terminal spin selection.
It does not bound occupation removal, floating arithmetic, preparation fit
uncertainty, laboratory noise-model error or apparatus calibration error.

The earlier anharmonic packet gave subdivision and occupation comparisons as
numerical diagnostics. The present addition is an explicit channel-norm bound
for the finite unsplit comparison, with constants computed for its actual
truncated operators. It requires no fitted decay rate or observed residual.
The time parameter throughout is the supplied simulated Schrödinger-frame
duration in seconds, not host runtime or native word count.

## 1. The complete finite generator

Let the oscillator occupation dimension be `d`, with occupations `0,...,d-1`,
and keep the two-level spin. Define `ad_V(rho)=[V,rho]`, and write

\[
 H_s=\delta N-s\frac{\Omega}{2}X(a+a^\dagger),\qquad
 A=-i\operatorname{ad}_{H_s},\qquad
 D=-\frac\gamma2\operatorname{ad}_N^2.
\]

Each nonnegative-time exponential of `A`, `D` or `A+D` is CPTP, and therefore
has diamond norm one. The number operator includes the identity on spin.
No canonical commutation identity from the infinite oscillator is substituted:

\[
 [a,a^\dagger]=I-d\Pi_{d-1}.
\]

All bounds below remain valid with an arbitrary spectator reference system.

## 2. Symmetric splitting without an exponential growth factor

For `h=tau/subdivisions`, the implemented subpulse is

\[
 S(h)=e^{hD/2}e^{hA}e^{hD/2}.
\]

The standard bounded-generator Strang argument gives

\[
 \|S(h)-e^{h(A+D)}\|_\diamond
 \le\frac{h^3}{24}\|[D,[D,A]]\|_\diamond
     +\frac{h^3}{12}\|[A,[A,D]]\|_\diamond.
\]

This is the contraction case of the double-commutator estimate in
[Iserles and Kropielnicka, Theorem 2.3 and Section 3.3](https://arxiv.org/html/2401.06635v1).
Here its hypotheses hold for bounded finite Lindblad operators in trace norm,
including after adjoining a reference system. We use its integral argument,
not just a formal third-order expansion.

For completeness, let `C=[A,D]`. Differentiating `S(t)` yields the defect

\[
 R(t)=S'(t)-(A+D)S(t)
 =[e^{tD/2},A]e^{tA}e^{tD/2}
  +\tfrac12e^{tD/2}[e^{tA},D]e^{tD/2}.
\]

Writing each commutator as a Duhamel integral and expanding about the shared
middle expression `e^(tD/2) C e^(tA)` cancels the two first-order terms.
The two remaining integrals have total scalar weights `t^2/8` and `t^2/4`,
respectively, multiplying `[D,C]` and `[A,C]`. All intervening exponentials
have **nonnegative** durations. Thus

\[
 \|R(t)\|_\diamond
 \le \tfrac{t^2}{8}\|[D,C]\|_\diamond
      +\tfrac{t^2}{4}\|[A,C]\|_\diamond.
\]

Integrating this defect against `exp((h-t)(A+D))` proves the stated bound.
There is no inverse diffusion, no assumption of unitarity for a dissipative
map and no large factor `exp(h*||D||)`.

## 3. Constants for the actual cutoff Hamiltonian

Put `n=d-1`, `j=Omega*sqrt(n)`, `C_N=ad_N`,
`J=-i[H_s,N]` and `L=-i[H_s,J]`. Centering `N` gives

\[
 \|C_N\|_\diamond\le n,\qquad
 \|J\|\le j,\qquad
 \|\operatorname{ad}_J\|_\diamond\le2j.
\]

The force part of `L` is `Omega^2*(I-d*Pi_(d-1))/2`, and the detuning part
has operator norm at most `delta*j`. Removing the scalar center of the force
part is legitimate for its commutator, giving

\[
 \|\operatorname{ad}_L\|_\diamond
 \le2\delta j+\frac{\Omega^2 d}{2}.
\]

The final term is a real boundary contribution of the finite oscillator.
Dropping it would invalidate the estimate; the tests exhibit a high/low
occupation coherence that violates the resulting incorrect bound.

Since `[A,C_N]=ad_J` and `[A,ad_J]=ad_L`,

\[
 [A,[A,D]]=-\frac\gamma2
 (\operatorname{ad}_L C_N+2\operatorname{ad}_J^2+C_N\operatorname{ad}_L).
\]

Consequently,

\[
 B_{AA}=\gamma\left[
  \left(2\delta j+\frac{\Omega^2d}{2}\right)n+4j^2\right]
 \ge\|[A,[A,D]]\|_\diamond.
\]

For the other double commutator, `F=[C_N,ad_J]=ad_[N,J]` satisfies
`||F||_diamond<=2j`, because `J` changes occupation by precisely one.
Expanding the double commutator with `C_N^2` gives four terms,
`C_N^2 F+2 C_N F C_N+F C_N^2`, up to an irrelevant sign. Therefore

\[
 B_{DD}=2\gamma^2 j n^2
 \ge\|[D,[D,A]]\|_\diamond.
\]

The resulting half diamond distance for one subpulse is at most

\[
 \epsilon_{\rm sub}=\frac{h^3}{2}
 \left(\frac{B_{DD}}{24}+\frac{B_{AA}}{12}\right).
\]

Both constants are independent of the sideband sign and spin rotation angle.
They vanish if either the drive or the dephasing vanishes.

## 4. Whole schedules, environment execution and selection

There are `4*K*subdivisions` sideband subpulses in `K` complete drive steps.
Telescoping CPTP contractions gives

\[
 \epsilon_{\rm split}\le
 \min\{1,\,4K\,\mathrm{subdivisions}\,\epsilon_{\rm sub}\}.
\]

An additional spin rotation with Hamiltonian `delta*N + c*X_spin` commutes
with number dephasing and adds no splitting error when these commuting
channels are composed. The dephasing channel used during its actual finite
duration still needs its environment-implementation budget. Other
noncommuting controls would need additional splitting terms.

For an ideal `M`-point Gaussian-Hermite mixture implementing dephasing for
duration `t`, the environment module supplies

\[
 \epsilon_{\rm env}(t)\le\min\left\{1,
 \frac d2\frac{M!}{(2M)!}
 [\gamma t(d-1)^2]^M\right\}.
\]

Hermite interpolation of `exp(-i*z*x)` has integrated remainder at most
`z^(2M)*M!/(2M)!` under the standard normal measure. For an arbitrary joint
state with a reference system, positivity bounds the sum of occupation-block
trace norms by `d`; this converts the Schur-entry bound to the displayed half
diamond radius. The finite mixture's phase rotations are natively executed;
the mixture weights and Gaussian diffusion law remain supplied.

For each type of environment invocation, its count multiplies its certified
radius. These terms are summed with `epsilon_split` and capped at one. This
composes ideal channels only: numerical contractions and quadrature-node
rounding are separately checked, not included in the mathematical certificate.

After any common terminal instrument, let the retained accepted substates
be `P` and `Q`, with masses `p` and `q`, and retain their rejected branches.
If the complete trace distance is at most `epsilon`, then

\[
 \|P-Q\|_1+|p-q|\le2\epsilon,
 \qquad
 \frac12\left\|\frac Pp-\frac Qq\right\|_1
 \le\frac{\|P-Q\|_1+|p-q|}{2q}
 \le\frac\epsilon q.
\]

Thus every conditional binary probability has absolute error at most
`min(1,epsilon/q)`, and the corresponding spin mean at most twice this value.
The same bound holds using either execution's positive selection mass, or
its positive lower bound. If `q<=epsilon`, the comparison does not guarantee
the other selection mass is positive and reports a vacuous conditional bound.
Treating a floating acceptance as exact does not certify its roundoff error.

## 5. Size and limits of the certificate

With the supplied `delta=2*pi*500`, `Omega=pi/(150 us)`, `gamma=18 /s` and
`tau=50 us`, the unconditioned splitting bounds are:

| Occupation dimension | Subdivisions per 50 us | K=50 | K=78 |
|---|---:|---:|---:|
| 32 | 128 | 0.00033954 | 0.00052968 |
| 64 | 256 | 0.00030686 | 0.00047870 |

These are uniform bounds over input states, detectors and reference systems.
They are selected by the numerical contract, without using later measured
values. They are much smaller than the earlier probability RMSE of about
0.056 before selection; each actual selected-law certificate uses its own
retained selection probability.

The code evaluates the algebra with exact rational inputs and an
outward-rounded rational square root, and rounds exported floats upwards.
A supplied float denotes its exact stored binary value; this does not turn
a laboratory calibration into an exact physical constant. Small subdivision
differences alone would not establish this certificate. Agreement across
occupation cutoffs likewise remains a diagnostic: no unrestricted-oscillator
or complete physical-error bound is claimed.
