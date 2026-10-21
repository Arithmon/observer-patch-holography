# One controlled enlargement of the interacting field band

This construction adds two previously omitted spatial harmonics to the
existing scalar field and executes their effect on the same continued
measurement history. It derives the required low-mode mass adjustment,
bounds the complete low-observer output, and certifies a positive
interaction-induced excitation of the added modes. It is a finite,
conditional step toward controlling the field bandwidth.

The [seven-mode scalar observations](../scalar_observations/README.md)
control grid spacing and occupation for a fixed Fourier band. Finer grids
there retain the same seven field degrees of freedom. This packet changes
that band, supplies its interaction terms by exact Fourier integration,
and compares the resulting nine-mode process with a matched seven-mode
description. It preserves the earlier receipts.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.interacting_band_extension --check model/interacting_band_extension/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/interacting_band_extension
```

## Source, band and observation contract

The supplied torus is `[0,2*pi)^3`, with volume `V=(2*pi)^3`. The normalized
real modes are the constant, sine/cosine pairs of each first harmonic,
and the added `cos(3x), sin(3x)` pair. The new pair is a selected resonant
extension; it is not the next complete momentum shell. Pairing preserves
translation closure of this Fourier band before occupation compression.
The rectangular finite occupation code does not preserve all rotations
mixing a sine/cosine pair. No finite-code translation covariance is claimed.

The canonical scalar action has bare `m=c=1`, `g=lambda/V=1` and oscillator
basis frequency `Omega=2`. The occupation dimensions, in the order above,
are `(2,4,4,2,2,2,2,4,4)`: 8,192 field amplitudes. The matched low code has
512 amplitudes. These cutoffs are premises of this comparison; the earlier
unrestricted seven-mode certificate does not extend to the new modes.

All `Q^2`, `P^2`, `Q^3` and `Q^4` matrices compress the full oscillator
polynomial after it is formed. Powers of the compressed `Q` would give a
different action. The finite detector is instead declared as a function of
that compressed `Q`, exactly as in the inherited observation protocol.

The first native repair on the binary total-two carrier selects the sign
of `exp(+-i Q_cos_x/4)` on the product oscillator blank. After model time
`1/16`, the observer retains both Ramsey outcomes
`K0=cos((Q_cos_x+pi/2)/2)`, `K1=sin((Q_cos_x+pi/2)/2)`. Bit one applies
`exp(i Q_cos_x/4)` to the disturbed field. A second interval of `1/16`
and the same detector produce the second bit, which selects the next
native E/F seam. All 1,056 native/field histories remain. The added modes
continue throughout both intervals and are traced only at final readout.

This uses the existing coherent code, initialization, routing and pointer
services. All local matrices are contracted from emitted twelve-port source
pulses, including the three-mode interaction words; independent target
matrices only check them. A 64-step symmetric product executes each
interval. The supplied Hamiltonian model time is shared by the compared
protocols; source pulse counts and host runtime do not calibrate it.
The construction needs no observer federation.

## The new terms and the matching identity

Write `x=Q_cos_x`, `y=Q_sin_x`, `z=Q_cos_3x`, `w=Q_sin_3x`, and let `W7`
be the existing dimensionless quartic polynomial, so the interaction is
`g W7/24`. Exact Fourier integration gives

```text
W9 = W7 + (3/2)(z^2+w^2)^2 + 6(z^2+w^2) sum_old Q_j^2
        + 2[z(x^3-3xy^2) + w(3x^2y-y^3)].
```

The last four terms create odd occupation in the added band. The code
independently obtains every one of the 49 monomials by summing its Fourier
coefficients. A uniform grid with at least 13 points per axis integrates
this quartic exactly. At 12 points `cos(3x)^4` aliases, providing a negative
control. Grid refinement within these nine modes still leaves their
bandwidth fixed. A separate rational `O(N^-2)` complete-history bound covers
the finite-difference kinetic error on that unchanged code.

Let `J psi=psi tensor |0,0>`, `P=JJ*`, and `Q=I-P`, where the two new
oscillator blanks have frequency two. Their exact moments give

```text
J* H9 J = H7 + (g/8) sum_old Q_j^2 + (7/2+g/32) I = H_eff.
delta(m^2) = g/4.
```

Thus keeping the same bare low mass is not the matched comparison. This
identity selects the immediate projected action for the declared blank;
it does not identify a physical pole mass or complete renormalization.
The scalar is retained in the identity. The fixed-duration source record
comparison may omit that common global phase.

There are five orthogonal excited-high output blocks of `Q H9 J`. Put
`S=sum_old P_j Q_j^2 P_j`. Their low operators are

```text
B20 = B02 = sqrt(2) [3/4 + g/32 + (g/16) S],
B22 = g/64,
B10 = (g/24) P_low (x^3-3xy^2) P_low,
B01 = (g/24) P_low (3x^2y-y^3) P_low.
```

The high occupation cutoff four is used here; larger codes add channels.
Exact rational row bounds on the compressed polynomials give
`||Q H9 J||^2 <= R^2 < 7.570328`. The even blocks include motion caused by
the frequency-two blank under the new harmonic's frequency `sqrt(10)`.
Total high excitation therefore is not, by itself, an interaction witness.

## Complete-history comparison theorem

The theorem holds for every initial low density matrix, including an
arbitrarily entangled untouched reference, with high preparation `|0,0>`.
Compare evolution under `H9` against `H_eff`, interleaving identical
low-only instruments and low-only record-controlled operations. Retain
every outcome. Evolution durations are prescribed; the same proof allows
a bounded worst-branch sum `T` when the duration controls remain within
this classical adaptive contract. No operation may read, reset, change or
postselect the added modes.

Purify the input and dilate the low instruments, retaining their records
and ancillary registers. Every intervention commutes with `P` and `Q`.
For the enlarged evolving vector `psi`, write `p=P psi`, `q=Q psi`; let
`phi` be the embedded matched evolution. Block Duhamel and contraction of
the interleaved isometries give

```text
||q|| <= R T,
||p-phi|| <= integral_0^T R ||q(s)|| ds <= R^2 T^2/2.
```

Tracing the added modes removes the `P/Q` cross blocks. The half-trace
norm of the retained output difference is at most

```text
(||p||+||phi||)||p-phi||/2 + ||q||^2/2 <= R^2 T^2.
```

Purification, partial trace and record dephasing preserve this bound. It
covers the complete low quantum output jointly with records and reference,
so it also bounds every complete classical record law without multiplying
by the number of branches. No normalization conditional on staying in the
blank subspace occurs. At `T=1/8` the bound is below `0.118287`.

If the added modes are retained as output, their cross coherences remain.
The direct isometry Duhamel bound is only `R T`. The test Hamiltonian
`X_low tensor X_high` attains reduced distance `sin^2(T)` and full pure-state
distance `|sin(T)|`, verifying this distinction. Tests also cover an
entangled reference, an intervening measurement and feedback, and the
linear error produced by forgetting `P H P` in the comparison action.

This is a direct finite-dimensional projection/Duhamel argument. The new
contribution is its exact mapping to the executed band extension and its
unpostselected observer histories. Hamiltonian truncation and coupling
matching are established methods; see
[Rychkov and Vitale, sections 2–3](https://arxiv.org/html/1412.3460v4).
The separate roles of spatial spacing, state cutoffs and matching are also
explicit in [Jordan, Lee and Preskill, sections 4–5](https://arxiv.org/html/1112.4833v2).
Neither source supplies a continuum theorem for the finite family here.

## Certified resonance and retained numerical controls

In a separate protocol, prepare the full product blank and read whether
the total occupation of the new pair is odd. Two source-compiled CNOTs
write that parity into a blank pointer, and both outcomes remain. The exact first derivative
has squared norm `g^2/768`. This follows from cancellation of the one-quantum
terms in `x^3-3xy^2` and `3x^2y-y^3`. For `g=0`, parity makes the odd event
impossible at every time, despite possible even high excitations.

Center `H9` by its blank expectation. Its exact finite-support second
derivative norm bounds the vector Taylor remainder by
`t^2 ||(H9-<H9>)^2 blank||/2`. At `t=1/1024`, subtracting this remainder
and the 32-step native product bound leaves a strictly positive amplitude.
The ideal executed odd-event probability exceeds `7.79e-10`. This is an
analytic lower bound for the finite declared code, with both read outcomes
retained. This high-mode diagnostic is outside the low-only comparison
theorem's output contract.

The principal replay retains comparisons with the independent sparse
Hamiltonian, the mass-matched low action, the unmatched bare low action,
feedback deletion and a trace-and-reprepare-high control. The last control
retains every high reset Kraus outcome; it does not select the blank outcome.
At 64 steps, the matched low record TV is approximately `3.37e-6`, versus
`1.75e-4` without the mass adjustment; resetting high modes changes the
continued low law by approximately `2.20e-6`. The enlarged source agrees
with its independent sparse reference to `1.24e-9` in record TV. These small
differences are numerical diagnostics, not certified positive lower bounds.
The complete reduced low field state jointly with the two records is also
compared, giving trace distance approximately `0.00033317` between the
enlarged and matched source outputs. Equality of bit probabilities alone
would not certify that stronger comparison.

The palindrome agrees with the exact exponential through degree two.
Centered local norm sum `S` gives operator error `S^3 t^3/(3 n^2)` per
interval by integral Taylor remainders and unitary telescoping. Adding the
two source errors to the projection bound compares the complete enlarged
source and matched low source histories; the combined bound is below
`0.129474`. The bound covers exact ideal
compilation. Floating pulse contractions, guarded polar corrections and
physical service errors are outside it. Raw full-phase residuals and
source/reference errors are reported separately, without branch rescaling.
The principal worst branch requests 2,950,900 source pulses and 5,901,800
code load/export services, plus the two native repairs. The resource ledger
counts requested operations rather than cached host computations;
preparation, hardware fidelity and classical control costs remain supplied.

## Interacting readouts on the same radial product

`same_product_interaction` in the retained receipt joins an interacting
finite field to the [radial-product geometry](../relative_observer_geometry/README.md#lorentzian-density-of-the-product-wave).
Take `(0,pi)` times the same angular screen, measure
`dr/pi dmu_screen`, and the two orthonormal modes `sqrt(2)sin(r)` and
`sqrt(2)sin(2r)` times the angular constant. Their squared free frequencies
are `1,4`. For a declared positive effective coupling `g`, the projected
Hamiltonian is

```text
H = (P1^2+Q1^2+P2^2+4 Q2^2)/2
    + g/24 [3(Q1^4+Q2^4)/2 + 6 Q1^2 Q2^2].
```

The quartic coefficients follow by integrating the four radial sine
factors; the executable tests independently integrate them. All oscillator
polynomials are formed before compression to four occupations per mode.
The frequency-two oscillator basis is only a matrix representation. An
independent padded-oscillator construction checks the complete Hamiltonian.
The interaction is positive and couples the two modes. Nonlinear evolution
also excites omitted modes: even `psi1^3` has nonzero projection onto
`psi3`. This is a Galerkin effective theory, with explicit spatial and
occupation cutoffs, and does not claim an invariant two-mode sector.

The geometric map is explicit. Write the product wave coefficient density
as `C=mu diag(-1,a)` and its Lorentzian metric as
`g_metric=sqrt(-det C) C^-1`. The lapse satisfies
`N^2=mu sqrt(det a)` and `sqrt(abs(det g_metric))=mu N^2`.
Thus a constant covariant coupling `lambda phi^4/24` projects to

```text
g_eff = lambda/pi integral mu_screen^2 sqrt(det a_screen) dx.
```

This invariant integral is positive. Setting `g_eff=1` chooses a coupling;
repair does not select it. The angular constant is likewise only a retained
mode when `N` varies. Smoothness of the **whole** coefficient density gives
a smooth Lorentzian background as specified in the linked geometry proof;
the native measurable coefficient theorem alone does not give smoothness.

The run uses `g=1`, `beta=1/4`, two model-time intervals `1/16`, and 64
palindrome steps per interval. It actually compiles and executes a
purification of `rho=exp(-beta H)/Tr exp(-beta H)` on four work and four
reference qubits, keeping every reference amplitude. The source output is
not normalized to repair drift. Its raw norm, preparation error and complete
resources are retained. The ideal observer algebra has
`sigma_s(A)=rho^(is) A rho^(-is)=alpha_(-beta s)(A)`, so the interacting
state and interacting clock use the same `H`. The raw prepared marginal's
modular action is checked separately. A free Gibbs state has nonzero
commutator with this interacting Hamiltonian and cannot supply this clock.
The parameter is supplied Hamiltonian model time; this protocol does not
calibrate its intervals against incident counts or laboratory seconds.

Both source and independent Hamiltonian executions retain measurement
disturbance, conditional feedback and all four two-bit outcomes, jointly
with the final work and reference. Native sign control yields 1,056
complete histories. The numerical complete-output trace error is about
`1.97e-8`; deleting the interaction changes the record law by about
`7.31e-4`, and deleting feedback by `1.95e-3`. The ideal splitting bound is
retained separately from the telescoped numerical source allowance.
The latter includes preparation, every unitary step, both measurement
dilations and both kicks, using the inherited raw-target and guarded polar
correction limits. Neither allowance certifies physical service fidelity.
The worst branch requests 217,110 source pulses, 434,220 proper-code
load/export events and two native repairs. These are a finite executable
interaction and clock consistency check, not a public-data comparison.

The inverse check fits `g`, `beta` and a radial frequency scale `c` from
six independent probabilities at model interval `1/16`, using an independent
padded-oscillator Hamiltonian. Three distinct initial guesses recover
`(1,1/4,1)` within `1e-4`; their fitted model predicts the withheld source
interval `1/8` within `1e-8` in every outcome probability. The local
sensitivity singular values are approximately `0.264713,0.001489,0.000110`:
the menu distinguishes these three parameters locally, but small probability
errors can amplify substantially. These are noiseless synthetic probabilities,
not finite-sample confidence bounds or a proof of global uniqueness. The
declared comparison family changes radial frequencies to `ck`; it corresponds
to length `pi/c` with normalized radial measure and independent effective `g`.
Holding the covariant coupling fixed would additionally change `g_eff` by `c`.
This fit identifies neither radial topology nor the full angular tensor.

An absolute energy scale is still needed: if `H_phys=E_* H_model`, then
`beta_model=E_*/(k_B T)` and `t_model=E_* t_phys/hbar`. Simultaneously scaling
the Hamiltonian by a constant and dividing inverse temperature and duration
by it preserves all these records. A legitimate calibration can fix that
unit; the remaining ratios and readout predictions must then be tested
without further adjustment.

## What remains open

One interacting bandwidth enlargement and its matching error are established
on a finite occupation domain. Neither repeated band enlargement nor an
unrestricted occupation limit for the enlarged process is proved. The
matched mass shift is blank- and cutoff-dependent, and later enlargement
can induce further operators and memory. A finite torus, finite anisotropic
band, finite code and prescribed model clock do not establish Lorentz
symmetry, spacelike commutativity, scattering, an interacting vacuum, or
observed spacetime. Those need a compatible larger-domain construction,
controlled omitted-mode effects and an independently fixed observational
meaning. The executable source compilation supplies none of these physical
identifications by itself. The same-product example supplies an explicit
conditional interacting effective description and its Gibbs modular clock;
it does not establish a nontrivial interacting ultraviolet limit, exact
finite-cutoff microcausality, a CFT, or the observed interacting QFT.
