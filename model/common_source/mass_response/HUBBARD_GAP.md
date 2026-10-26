# An alternative architecture with an interaction-generated charge gap

`gap_certificate.py` gives a quadrature-free rational enclosure for a
half-filled repulsive one-dimensional Hubbard chain. This is a separate
coherent architecture from the finite dimerized propagation experiment.
The charge-gap theorem does not turn that primary experiment into a Hubbard
chain, or turn its finite-volume spectra into a measured thermodynamic gap.

## Source values and the architectural map

The certificate calls the existing
`model.native_modular_reference.local_preparation.reference()` and reads all
66 binary source configurations. Core occupation `q=0,1,2` has multiplicities
`45,20,1` and respective weights per configuration

\[
(w_0,w_1,w_2)=(2/225,1/36,2/45).
\]

The two exact statistics are

\[
\Omega=\log(w_2/w_0)=\log5,\qquad
J=\log\frac{w_1^2}{w_0w_2}=\log\frac{125}{64}.
\]

The second statistic removes a normalization and an additive one-body
occupation potential. Its nonzero value is a connected native-reference
interaction statistic. The source modular action alone remains diagonal:
these statistics do not by themselves define propagation or a mass operator.
The certificate consumes the exact ideal reference contract. Finite preparation,
herald, pulse and sampling errors have not been propagated through this
alternative architecture's threshold estimate.

The alternative architecture chooses `t=Omega/2` and `U=J` in

\[
H=-t\sum_{j,\sigma}(c^\dagger_{j\sigma}c_{j+1,\sigma}+\mathrm{h.c.})
 +U\sum_j(n_{j\uparrow}-\tfrac12)(n_{j\downarrow}-\tfrac12).
\]

This model has uniform hopping and no bare staggered or Dirac mass. The
centered interaction sets the particle-hole-symmetric chemical-potential
convention; it does not insert the charge gap. The assignments of the source
statistics to the hopping and quartic channels, the chain architecture,
half filling and ground-state sector are explicit construction choices.
Existing native repair does not uniquely force them. The prepared native
reference is not claimed to be the Hubbard ground state.

## The exact thermodynamic input

[Lieb and Wu, section 7, equations 66--67](https://arxiv.org/html/cond-mat/0207529)
derive the positive jump between addition and removal chemical potentials at
half filling. [Jeckelmann, Gebhard and Essler, equation 2](https://arxiv.org/html/cond-mat/9911281)
state its positive-integral form, with the present hopping convention:

\[
\Delta_\mu=\mu_+-\mu_-
 =\frac{16t^2}{U}\int_1^\infty
 \frac{\sqrt{y^2-1}}{\sinh(2\pi t y/U)}\,dy.
\]

This is the thermodynamic charge and optical threshold for `t>0,U>0`.
It vanishes at `U=0`; positivity for every positive interaction follows
directly from the integral. These are theorem inputs about the stated
Hubbard law, not conclusions about arbitrary repair dynamics.

## A finite exact enclosure without quadrature

Put `a=2*pi*t/U`. Expanding the positive reciprocal hyperbolic sine and
integrating each term gives

\[
\Delta_\mu=\frac{16t}{\pi}\sum_{r\ge0}
 \frac{K_1((2r+1)a)}{2r+1},\qquad
K_1(z)=z\int_1^\infty e^{-zy}\sqrt{y^2-1}\,dy.
\]

The last integral becomes

\[
K_1(z)=z e^{-z}\sqrt2\int_0^\infty
 e^{-zx}x^{1/2}\sqrt{1+x/2}\,dx.
\]

For every `x>=0`, the even-order Taylor polynomial of `sqrt(1+x)` is a
lower bound and the odd-order polynomial is an upper bound. This follows
from the sign of the next derivative on the whole positive axis. Integrating
these polynomials against the positive weight above gives global bounds

\[
\sqrt{\frac\pi{2z}}e^{-z}S_{2n}(z)
 \le K_1(z)\le
\sqrt{\frac\pi{2z}}e^{-z}S_{2n+1}(z),
\]

where `S_n(z)=sum_(k=0)^n c_k/z^k` and

\[
c_0=1,\qquad c_k=c_{k-1}\frac{4-(2k-1)^2}{8k}.
\]

These are finite Taylor-remainder inequalities. Their validity does not
assume convergence of an infinite asymptotic expansion. A lower polynomial
may become negative at other parameters; clamping the final lower endpoint
to zero preserves a valid enclosure.

The certificate retains terms `r=0,1`. Integration by parts gives the
equivalent representation
`K1(z)=integral_0^infinity exp(-z*cosh u)*cosh u du`, hence
`K1(n*a)<=exp(-(n-1)*a)*K1(a)` for `n>=1`. The omitted positive series is
therefore at most

\[
K_1(a)\frac{e^{-4a}}{1-e^{-2a}}.
\]

Every remaining elementary calculation is bounded by exact `Fraction`
arithmetic with outward dyadic rounding:

- `log`: finite `2*atanh((x-1)/(x+1))` sum and a positive geometric tail;
- `pi`: Machin's identity and alternating arctangent remainders;
- square roots: integer `isqrt` endpoint enclosures;
- negative exponentials: reciprocal positive exponential Taylor sums with a
  geometric remainder, range reduction and repeated squaring;
- coefficient, tail and ratio operations: rational interval arithmetic.

The authoritative JSON endpoints are exact rational strings. Its float
endpoint displays also round outward. No floating quadrature enters the
certificate. The default order-8/order-9 enclosure gives approximately

\[
0.0010272435834226310<\Delta_\mu<0.0010272439173696570,
\]

or `0.0006382623246826815 < Delta_mu/Omega < 0.0006382625321756352`.
Higher global Taylor bounds can tighten the interval. Independently
integrating the original Lieb--Wu Bessel expression is a test diagnostic,
not an input or error certificate for these bounds.

The exact zero-interaction control holds the clock and hopping fixed while
disabling the quartic channel. It is a controlled alteration of the coherent
architecture, not another native reference with its full modular flow held
fixed. The exponential upper envelope also tends to zero as `U` decreases
to zero at fixed positive `t`.

## What can be identified

The individual holon/antiholon rest-energy threshold in the symmetric
convention is **`Delta_mu/2`**. The current response creates a neutral pair
and begins at `Delta_mu`. The spin sector remains gapless, so this is not a
positive gap of the full many-body Hamiltonian or an isolated electron pole.

The same primary paper's equations 6--8 describe the **separate 1+1-dimensional
low-energy field-theory limit**: its bare Thirring action has no mass term,
while its charge excitations obey `E(P)^2=P^2+(Delta_mu/2)^2` in the stated
velocity convention. The present certificate supplies no quantitative
fixed-lattice Lorentz-error bound, physical clock calibration, 3+1-dimensional
identification or observed Standard Model mass.

Finite small chains can check operators, instrument execution and interaction
responses. They cannot resolve this tiny thermodynamic threshold by a naive
size extrapolation. Periodic `L=4m` and `L=4m+2` free chains already have
different shell effects. The exact theorem, rather than a small finite gap,
supports the stated thermodynamic conclusion.

Run the certificate or its independent tests from the repository root:

```sh
python -m model.common_source.mass_response.gap_certificate
python -m pytest -q model/common_source/mass_response/test_gap_certificate.py
```
