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
`psi3`. This is a Galerkin effective theory with a fixed spatial band. Its
executed code has four occupations per mode; the occupation limit below
keeps the same two spatial modes and does not make their sector invariant.

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
The parameter in this fixed-duration run is supplied Hamiltonian model
time. The incident controller below implements that same evolution on a
native event clock; no laboratory-second calibration is supplied.

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

## Native incidents drive the same interacting field

`product_incident_clock_control` executes the interacting Gibbs preparation,
Hamiltonian, disturbed reads and feedback under a native timing controller.
Every incident on its selected port triggers one compiled interval `1/16`,
including incidents with zero repair current. The timing carrier pauses
during field services, and the work/reference register is protected between
them. The finite fixture retains all 3,600 two-attempt seam/coin tapes and
4,900 tape/field-prefix outcomes. Zero or one incident leaves an explicit
timeout with its full surviving work/reference state; the timeout probability
is `35/36`. Grouping identical source words is computational reuse, not
postselection. Counting only nonzero currents changes the clock law.

The separate refinement contract uses the full screen's IID seam schedule.
For port degree `d`, seam count `M_l=630*4^l` and normalized attempt rate
`a_l=1260*16^l`, incidence probability is `p=d/M_l` and rate is
`R_l=2d*4^l`. An incident executes the same interacting generator for
`delta=1/R_l`; reads occur at counts `n=ceil(R_l t)` and `2n`. Doob's
inequality for centered incident counts through
`K=ceil[a_l(2n/R_l+eta)]` bounds timeout or timing failure by
`K p(1-p)/(R_l^2 eta^2)`. On the complementary event, both arrival
times are within `eta+1/R_l` of the rounded model read times. This is a
bound in the normalized native-attempt coordinate; paused service durations
are outside it. The 12-port fixture is not a full refined-screen execution.

Let `S` bound the centered local generator norms at the chosen occupation,
`m=3` be the term count, `s` the palindrome subdivisions per incident and
`J=2n`. On the **same retained native tape**, the complete source output
differs from the exact interacting history by at most

```text
E_preparation + (2 s m J + 4) epsilon + J S^3/(3 R_l^3 s^2).
```

Here `epsilon` is a complete-instrument half-diamond error including
failures, not a successful-matrix residual. Comparing work outputs at
nominal read intervals adds `2S(n/R_l-t)` and the timeout bound. The clock
tape stays in both compared outputs; stochastic timestamps are not claimed
to converge in total variation to deterministic point masses. The default
requirements `E_preparation=R_l^-2`,
`epsilon=[R_l^2(2smJ+4)]^-1` and
`eta=2^(-floor(l/2)-1)` make these errors vanish at fixed occupation.
They are service-accuracy requirements, not hardware measurements.
Since the Gibbs state uses this same `H`, an ideal incident advances modular
parameter by `-delta/beta`; at `beta=1/4` the executed fixture advances it
by `-1/4`. Source splitting and precision errors are charged above. The
active coupling rule, independent timing carrier and protected service
schedule are supplied, and do not establish autonomous repair dynamics or
SI time units.

## Unrestricted occupation for the same two-mode history

Fix `g>=0` and `beta>0`, keeping the radial modes, effective coupling,
detector and Hamiltonian model times above. On `L2(R^2)` the displayed
Hamiltonian is the Friedrichs operator of the closed positive form
`h[psi]=||grad psi||^2/2+integral V|psi|^2`. The quartic part of `V` is
nonnegative, so `H>=H_harm`, whose energies are `n1+2n2+3/2`. Let `P_d`
retain the first `d` frequency-two Hermite functions in each coordinate,
and `H_d=P_d H P_d` on that finite subspace. Forming each polynomial before
compression is necessary for this identity.

These subspaces form a common form core. Spatial cutoffs and mollification
approximate the form domain by smooth compactly supported functions;
their Hermite coefficients decay faster than any power. Hermite truncation
therefore converges with every polynomial weight and derivative used in
`h`. Galerkin orthogonality makes `(H_d+1)^-1 P_d f` the form-norm
projection of `(H+1)^-1 f`, proving strong resolvent convergence. Equivalently,
these are decreasing closed forms with increasing finite domains in
[Simon's nondense-form convergence theorem, Sections 3–4](https://math.caltech.edu/SimonPapers/81.pdf).
Consequently `P_d exp(-it H_d) P_d` converges strongly to `exp(-it H)`.
The zeros outside `P_d` are embedding conventions: exponentiating a
zero-extended Hamiltonian would incorrectly add an infinite identity to
the Gibbs operator.

Min-max gives `lambda_j(H_d)>=lambda_j(H_harm)` and, for each fixed index,
`lambda_j(H_d)` decreasing to `lambda_j(H)`. The harmonic form has compact
embedding, so `H` has compact resolvent and

```text
0 < exp[-beta(13/8+5g/128)] <= Z_d,Z
  <= exp(-3beta/2)/[(1-exp(-beta))(1-exp(-2beta))] < infinity.
```

The lower bound is Jensen's inequality on the frequency-two product blank,
which belongs to every code. Extend the finite eigenvalue list by infinity.
Harmonic eigenvalue domination and dominated summation give `Z_d -> Z`.
Strong convergence of the positive Gibbs operators plus convergence of
their traces gives trace-norm convergence. Explicitly, for any positive
trace-class `T` and finite-rank `P`, with `Q=I-P`,
`||T-PTP||_1 <= Tr(QT)+2 sqrt(Tr(T) Tr(QT))`; apply this to the Gibbs
operators after making the limiting tail small. Thus the embedded
`rho_d=exp(-beta H_d)/Z_d` converges to the faithful `rho` in trace norm.
The exact rational `product_gibbs_tail` bounds the normalized ordered
spectral tail beyond rank `k(k+1)`: putting `q=1/(1+beta)`, its
unnormalized upper bound is
`q^(2k)[k+1-k q^2]/[(1-q)(1-q^2)]`. This controls Gibbs tightness; it is
not an error bound for an occupation box or for a four-level observation.

The untouched reference also has a specified comparison. If the compiled
purification is `V_d diag(sqrt(w_d))`, its reference labels denote energy
eigenvectors. Map that reference by `conjugate(V_d)` into the common
Hermite reference; the resulting amplitude matrix is `sqrt(rho_d)`.
These canonical purifications converge in Hilbert–Schmidt norm: their
square roots converge strongly, their squared Hilbert–Schmidt norms are
one, and finite-rank tests give weak Hilbert–Schmidt convergence. This
reference isometry changes no work measurement or native record and
retains the entire reference. Arbitrary eigenvector signs or degenerate
reference labels cannot be compared across cutoffs without it.

The actual detector uses functions of `Q_d=P_d Q P_d`, rather than
`P_d f(Q) P_d`. Finite Hermite vectors are an operator core for `Q`, and
`Q_d` agrees with `Q` on each such vector for sufficiently large `d`.
Hence its bounded continuous functions converge strongly, including the
signed preparation and feedback kicks and both cosine/sine Kraus
operators. Their adjoints converge also; the two Kraus operators obey
exact completeness at every cutoff. Each fixed finite adaptive history
is a finite product of these contractions and propagators. Acting on
the convergent canonical purification gives convergence of every
unnormalized branch vector, and therefore trace-norm convergence of the
complete output with all classical records, surviving work and reference.
Appending the common native continuation is a channel and preserves this
convergence. This extends the existing
[smooth-history occupation argument](../scalar_observations/OCCUPATION_DERIVATION.md)
to this two-mode Gibbs preparation. It gives neither uniform diamond-norm
convergence over all input energies nor a numerical error for the existing
four-level run.

The same limiting state and dynamics satisfy, for bounded observables,
`sigma_s(A)=rho^(is) A rho^(-is)=alpha_(-beta s)(A)`. The finite modular
flows, on `P_d A P_d`, converge through the strong propagator limit; no boundedness of
`log(rho)` is asserted. Occupation refinement therefore preserves the
interacting Gibbs clock and complete histories. Coupling, temperature,
radial boundary conditions and the two-mode spatial projection remain
declared choices.

`product_occupation_refinement` retains ideal comparisons at `d=4,8,12`.
The partition functions are approximately `3.64485,5.04532,5.31594`; the
successive complete work/reference/record trace distances are `0.56045`
and `0.24501`. Thus the executed four-level thermal history has substantial
occupation dependence. These finite diagnostics neither certify its
infinite-occupation error nor replace the convergence proof. Independent
controls check nested polynomial compression, harmonic eigenvalue
domination, reference-basis alignment and the growing splitting budget.

To transfer this ideal limit to source executions, choose a diagonal
schedule with preparation error and summed complete-operation errors
tending to zero. At each finite `d`, the centered norm sum `S_d` gives
the retained palindrome error `S_d^3 t^3/(3 n_d^2)` per interval; choose
`n_d` accordingly. The incident-clock comparison likewise needs its
finite-`d` norm and waiting-time terms to vanish along the same schedule.
Fixed step counts, fixed source precision and a fixed native clock
refinement do not suffice uniformly in `d`. Protected services and
arbitrarily accurate growing-code preparation/control are implementation
premises; the source execution retained here uses `d=4`.

## What remains open

One interacting bandwidth enlargement and its matching error are established
on a finite occupation domain. Repeated band enlargement and unrestricted
occupation for that nine-mode comparison are not proved. The two-mode
occupation theorem above keeps its spatial band fixed. The
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
