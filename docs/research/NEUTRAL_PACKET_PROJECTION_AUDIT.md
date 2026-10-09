# Neutral packet projection: exact circle reduction and numerical repair

This contribution to [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
repairs the evaluation of the already defined neutral Whitney packet. It
derives a closed form for the original circle integral, eliminating angular
quadrature and its aliasing. The supplied action, physical interpretation and
state are unchanged. This is a correction of finite evidence, not selection
of a physical model under #1025/#1026.

## Reproduction before repair

Reviewed baseline: `b8714c986a76ef3b83df88d64d94e67b6e4d8944`.
All examples below have width and hbar one. Indices are zero based: the
first 30 coordinates are radiative; the final 26 are 13 scalar real parts
followed by 13 imaginary parts.

1. Set the center to zero, point coordinates `(0,4,8)` to `(1,1,1)` and
   momentum coordinates `(0,4,8)` to `(2**53,1,-2**53)`. Other entries vanish.
   Every input is exactly representable in binary64. The original dot
   product lost the exact phase one. The returned neutral amplitude was
   `3.1605377469223292e-12` instead of
   `1.7076458324454292e-12 + 2.659500810425262e-12j`, changing interference.
2. Set the center to zero, scalar point coordinate 30 to one and scalar
   momentum coordinate 30 to 256. The default 256-angle quadrature returned
   `+1.3759422606592422e-11`; the original integral gives
   `-4.838248674915097e-12`. Rotating the point by 0.373 radians changed the
   purportedly neutral result as well. The nearby momentum-10 case passed.

The five controls in [test_neutral_packet_projection.py](../../code/electromagnetism/test_neutral_packet_projection.py)
were committed in `614a9e71` before the repair: four scientific failures and
one ordinary-case pass on the baseline. The old function already disclaimed
a certified quadrature error; the repair does not reinterpret that disclaimer
as an accuracy guarantee.

## Exact reduction from the supplied seed

Write a point as `(x_r,y)`, its center as `(c_r,c)` and momentum as `(p_r,p)`.
Here `r` denotes the 30 radiative coordinates, and `y,c,p` have 26 scalar
coordinates. Let `J(u,v)=(-v,u)` and `R_theta=cos(theta) I+sin(theta) J`.
For positive width `s` and hbar `h`, the seed Lebesgue half-density is

\[
 g(x)=(2\pi s^2)^{-14}
 \exp\left[-\frac{|x-c_{\rm full}|^2}{4s^2}
       +\frac{i}{h}p_{\rm full}\mathbin\cdot(x-c_{\rm full})\right].
\]

Rotating scalar center and momentum together gives
`g_theta(x)=(2 pi s^2)^(-14) exp(C+U cos(theta)+V sin(theta))`, with

\[
\begin{aligned}
C&=-\frac{|x_r-c_r|^2+|y|^2+|c|^2}{4s^2}
   +\frac{i}{h}\bigl[p_r\cdot(x_r-c_r)-p\cdot c\bigr],\\
U&=\frac{y\cdot c}{2s^2}+\frac{i}{h}p\cdot y,\\
V&=\frac{y\cdot Jc}{2s^2}+\frac{i}{h}(Jp)\cdot y.
\end{aligned}
\]

In particular the constant scalar phase is `-p dot c/h`. The imaginary
coefficient in `V` is `(Jp) dot y/h`, with the displayed sign.
Expand

\[
 U\cos\theta+V\sin\theta
 =\tfrac12(U-iV)e^{i\theta}+\tfrac12(U+iV)e^{-i\theta}.
\]

The exponential series converges absolutely and uniformly on the compact
circle. Only equal powers of the two Fourier factors survive its average.
Consequently, for all such original real data, including nonzero `B` below
and complex coefficients `U,V`,

\[
 P_0g(x)=(2\pi s^2)^{-14}e^C
   \sum_{k=0}^\infty\frac{[(U^2+V^2)/4]^k}{(k!)^2}
 =(2\pi s^2)^{-14}e^C I_0(\sqrt{U^2+V^2}).
\]

This is an entire function of `U^2+V^2`; there is no square-root branch
choice in the resulting value. The last equality is the defining Bessel
series ([DLMF 10.25](https://dlmf.nist.gov/10.25)); the usual single-coefficient
integral is recorded in [DLMF 10.32](https://dlmf.nist.gov/10.32).
This derivation does not assume successful quadrature or a restricted
momentum orientation.

The [existing norm theorem](../../paper/tex_fragments/WHITNEY_QUANTUM_PACKET.tex)
gives

\[
 A=\frac{|c|^2}{4s^2}+\frac{s^2|p|^2}{h^2},\quad
 B=\frac{p\cdot Jc}{h},\quad Z=\sqrt{A^2-B^2},\quad
 \|P_0g\|^2=e^{-A}I_0(Z)>0.
\]

Cauchy--Schwarz and the arithmetic--geometric mean inequality give
`A >= |B|`. The normalized neutral half-density therefore is

\[
 f_0(x)=(2\pi s^2)^{-14}
        e^{C+A/2}\frac{I_0(\sqrt{U^2+V^2})}{\sqrt{I_0(Z)}}.
\]

This evaluates exactly the state defined in the paper. Neutrality follows
also directly: a scalar rotation of `y` rotates the pair `(U,V)`, preserving
`U^2+V^2` and `|y|^2`. The physical configuration wavefunction still requires
the existing metric factor `rho^(-1/2)`; this function returns the Lebesgue
half-density, not that wavefunction or a propagated state.

## A normalization consequence

For `p=h Jc/(2s^2)`, one has `A=B`, `Z=0` and `U^2+V^2=0`.
Also `p dot c=0`, and the scalar part of `C+A/2` is exactly
`-|y|^2/(4s^2)`. The normalized scalar projection is thus the centered
Gaussian for every center size. This proves why rejecting every state whose
separately reported projection norm underflows would be incorrect. For
`c[0]=1e100`, `p[13]=.5e100`, `s=h=1`, the amplitude at the origin is the
ordinary value `(2 pi)^(-14)`, although the squared projection norm cannot
be returned reliably as binary64. Its standalone reporting API must still
refuse that norm. The normalized-state API evaluates the combined expression.

## Numerical contract and independent controls

Pointwise entry points retain the original finite, exactly binary64-
representable real-scalar input contract. They reject malformed vectors,
Booleans, masked values, nonfinite scalars and inexact coercions before
array conversion. Width and hbar must be positive. Displacements, dot
products, `U^2+V^2`, `A^2-B^2` and `Re(C)+A/2` are computed from the original
inputs as exact fractions. No numerical zero floor is used.

Private high-precision arithmetic evaluates the final normalized expression
before conversion. Its precision budget covers large absolute phases and
cancelling logarithms, with successive higher-precision evaluations required
to agree relatively to `1e-30`. The final complex binary64 result must agree
with that calculation relatively to `1e-12`. This is a normwise numerical
stability policy, **not an interval certificate or proof of an exact node**.
Unresolved zeros and unreportable outputs raise an error. The public seed
logarithm additionally requires its returned unwrapped phase to have at most
`1e-12` radians absolute conversion error; a small relative error in a huge
phase alone would not suffice. The projection evaluates its phase internally
and can remain evaluable when that unwrapped binary64 report is unavailable.

The legacy `nodes` argument remains validated as an integer at least 16, but
has no effect on the closed-form evaluation. Existing callers need no API
change. Raising the angle count is no longer an accuracy intervention.

Controls compare against integration of the original rotated seed and its
positive overlap, and against independently simplified cases. These include
complex scalar centers and momenta, nonunit width and hbar, interference,
near Bessel zeros, rotation, coordinate permutations, normalization
compensation, caller precision contexts and valid/reporting-range pairs.
The controls run with warnings treated as errors on Windows and Linux.
Their numerical integrations are independent diagnostics, not certified
interval enclosures.

## Evidence and claim impact

The live packet receipt stores phase-space, projection-norm and radius
diagnostics; its builder does not call either repaired pointwise function.
It is canonically regenerated to bind the repaired sources and retained
tests. Windows regeneration preserves every numerical and contract value
exactly, including the arbitrary Coulomb frame; only three existing source
pins and one added test pin change. A Linux rebuild of the baseline and
repaired producers also gives identical numerical payloads in that same
environment. Its BLAS-dependent choice of frame is not substituted into the
committed receipt. Existing paper formulas and the `OPH-WHITNEY-QUANTUM-PACKET` claim
remain valid without a statement change. The classical parent trajectory,
frozen receipts, Lean proofs and paper/book sources are unchanged. The new
note is included through canonical active-surface inventory regeneration.

The scientific effect is specific: a public evaluator now preserves the
phase, sign and rotation invariance of its defined state in cases where it
previously failed, and evaluates representable normalized states without
requiring representable intermediate norms. It supplies neither a new
physical clock nor a quantum propagation certificate.

## Initial validation

The twelve affected Whitney modules pass with warnings treated as errors:
719 tests on Linux (Python 3.12), and 718 on Windows (Python 3.13) with one
expected skip because Windows `longdouble` has no additional precision.
The 80 new pointwise controls are included in the Windows/Linux packet
workflow. The refreshed receipt passes the standalone independent verifier;
postdiction ledger parity, observation/premise register checks and the
active-surface inventory check pass. The inventory adds only this note.

Ten isolated, non-equivalent mutations fail their intended controls: restore
the rounded seed dot product or projected dot product, flip the imaginary
`Jp` term, omit the constant scalar phase, omit normalization, use a fixed
80-digit initial budget, discard Bessel phase/sign, disable conversion
precision checks, reject every valid input, or return constant zero. The
unmutated isolated controls pass. In particular, fixed initial precision
can make two numerical evaluations agree on a wrong value after enormous
logarithms cancel; agreement alone cannot replace the scale-derived budget.
These checks cover the specified failure mechanisms, not every possible
implementation error.
