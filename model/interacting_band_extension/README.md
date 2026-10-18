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
identifications by itself.
