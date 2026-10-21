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
active coupling rule, independent timing carrier and protected source
services are supplied. The finite scheduling logic has the homogeneous
realization below; this does not derive the apparatus from passive repair
or calibrate SI time units.

<a id="one-fixed-instrument-controls-the-clock-and-field"></a>

### One fixed instrument controls the clock and field

`ProductControllerInstrument` realizes this same finite history with one
time-independent instrument. Its classical registers store an instruction
pointer, local-word cursor, attempt count, committed field bits, native
loads, the complete native tape and a failure flag. The admitted domain is
the consistent checkpoints: every native tape reverses to a binary
two-particle state and replays to its recorded currents and final loads;
the instruction and field prefix agree with the number of pending
incidents. Malformed register patterns are outside this contract.

The fixed instruction table first kicks the prepared field. A `ready`
instruction enables the next uniformly selected seam and placement coin.
Every port incident, including a null-current incident, enters `evolve`.
That instruction advances through the source-executed local palindrome;
`read` commits its pointer bit, and the first bit alone selects `feedback`.
No native repair is enabled while these services are busy. After two native
attempts the controller halts, retaining zero-, one- and two-read outputs.
The selector receives no host iteration, model-time variable or hidden
quantum amplitude. The Gibbs preparation is the same source-executed prefix.

For a consistent classical checkpoint c and emitted outcome a, write its
next checkpoint as `delta(c,a)`, native weight as `p_(c,a)` and source work
map as `K_(c,a)`. The instrument has Kraus operators

```text
J_(c,a) = |delta(c,a)><c| tensor sqrt(p_(c,a)) K_(c,a).
```

The incoming controller is classical: these operators dephase its different
checkpoint labels. At `ready`, the sixty weights sum to one; at `read`, the
two source pointer Kraus maps are complete; other active services are the
existing exact-code unitary words, and halt is absorbing. Consequently
`sum J*J=I` on the admitted control and valid work-code domain. Every
finite branch has exactly the native and source word of the original
history. Induction on instructions therefore proves equality of the
complete instruments for every work/reference state, including correlations
with a classical native-load register. It asserts no coherent
superposition of classical timing-carrier loads. The executable audits all
sixteen work-input columns, 3,600 native tapes and 4,900 terminal outcomes;
common busy subroutines are reused only when their entire source words
coincide. Removing feedback or ignoring null incidents breaks the audit.

The same controller also realizes the
[postcommit persistence-observer interface](../common_source/DERIVATION.md#a-native-postcommit-operational-observer).
Restrict to its invariant domain after the first native record has been
written. Each first-record projector `Q_p` is fixed by `T*`, because every
transition copies the old native prefix. After a field bit is committed,
its prefix has the same property on that later invariant domain. The
executable checks every emitted native transition and each reused busy
instruction class, including absorbing halt/failure and idle transitions.
Thus `G=T*-I` annihilates these public projectors. Use the actual first-write
state, public/access algebra `span{Q_p}`, and
`E(X)=sum_p Tr(Q_p X Q_p)/rank(Q_p) Q_p`; ranks may differ. Continuing the
same instrument or idling are available interventions and both fix this
algebra. One may choose between them by the retained label. With the same
issued/retained evidence labels, prediction and unused-buffer embeddings as
the cited construction, every finite observer-interface clause follows on
this product program. These buffers and interface choices remain declared;
no new Lean inhabitant is asserted. The first write itself changes its
projector from zero to probability `1/60` for a selected seam/coin record,
so global blank-to-written durability is false. The proof preserves public
prefixes, not arbitrary quantum work observables, and does not identify
the product registers with the separate native checkpoint carrier.

The existing six-level code-export Kraus maps separately give all genuine
complement outcomes. `product_controller_export` routes each to an
absorbing failure sector and retains its reference state and complement
label. These sectors form a direct sum with the successful code and tensor
with untouched spectators. Their probability is exactly zero on valid-code
source words. The 4,900-outcome run is the complete valid-code history;
the independent export control includes a nonzero off-code component and
its entangled reference. No failed branch is normalized or discarded.

There are at most 775 controller instructions after preparation, with the
same 217,110 source pulses including preparation, 434,220 load/export
events and two native timing repairs. One instruction may call an entire
compiled source word; its primitive cost remains counted. Controller
instruction count is neither incident count nor physical duration.
For fixed `0<=q<1`, the single lazy instrument `q Id+(1-q) T` merely
inserts recorded idle outcomes. If L bounds active instructions and
`mu=(1-q)J>=L`, then after J invocations

```text
Pr(unfinished) <= exp(-(mu-L)^2/(2 mu))
               <= 1 / (1+(mu-L)^2/(2 mu)),
E[invocations until halt] <= L/(1-q).
```

All unfinished work/reference states are retained at a finite invocation
budget. Erasing only idle labels from completed histories gives the same
native/field instrument; conditioning a finite-budget run on finishing
would generally bias it and is not used. An implementation error epsilon
per complete controller invocation adds at most `J epsilon` to the
existing source and splitting budget. This is a construction of the finite
scheduler, not a certificate of controller hardware, storage or source
service accuracy. The instruction table and architecture remain declared.
The passive mixing channel is still distinct from the nontrivial modular
unitary of the interacting Gibbs state.

<a id="one-reference-clock-for-complete-intervened-histories"></a>
### One reference clock for complete intervened histories

`product_modular_history_control` connects the same fixed controller to
the interacting Gibbs reference on its complete recorded channel. Write
`rho=exp(-beta H)/Z`, `sigma_s(A)=rho^(is) A rho^(-is)`, and
`s=-t/beta`. The Schrödinger modular segment is
`rho^(-is)=Z^(is) exp(-itH)`. Substituting this segment into every
elapsed-time interval changes a history Kraus operator only by the scalar
phase `Z^(i n s)`, where n is that history's zero-inclusive incident count.
Every intervening read and record-conditioned feedback stays in its original
order. The record labels are committed classical outputs. Thus the ideal
complete channels agree, including timeouts, all
native tapes, field records, surviving work and arbitrary references. This
is a statement about the entire instrument, not only a two-point function.
Each tick in the executed fixture has `t=1/16`, `s=-1/4` and `beta=1/4`;
the refining clock contract above replaces these by `t=1/R_l` and
`s=-1/(beta R_l)` with its existing timing and source errors.

The actual compiled word approximates this channel. For its maps A_j and
phase-aligned ideal modular maps B_j, including square roots of the native
branch weights, the direct-sum Stinespring estimate is

```text
e = || sum_j (A_j-B_j)^* (A_j-B_j) ||^(1/2),
half_diamond(Phi_A-Phi_B) <= e (1+||sum_j A_j^* A_j||^(1/2))/2.
```

The 4,900-output, all-sixteen-input execution gives a numerical bound
`1.69e-8`; its actual prepared work/reference output differs by `3.44e-9`.
Both fit the existing splitting and complete service budget, without
normalizing source amplitudes or discarding records. The all-input check
does not establish access to every work observable. The Hamiltonian,
reference, reads and relative cocycle preserve first-mode parity and
therefore restrict to its block algebra. Conjugating the full
ideal intervention menu, including conditional feedback, by a common
`V=exp(-iuH)` gives `K'_j=V K_j V*` on every input. Leaving only feedback
untranslated fails this identity. The translated menu is an algebraic
control; no extra compiled implementation is claimed.

The fixed reference must also be distinguished from a conditioned work
state. The faithful first-read posterior has a different modular action:
the retained probe mismatch is about `0.0858`. For that posterior eta,
`u_s=eta^(is) rho^(-is)` obeys
`u_(s+t)=u_s sigma_s(u_t)` and
`sigma_s(A)=u_s* sigma_s^eta(A) u_s`, verified independently in the control.
This relative cocycle describes the change of reference; it is not an
additional operation executed by the controller. The same-H continuation
therefore uses the fixed Gibbs reference throughout the disturbed history.
The committed classical record algebra alone has trivial modular flow on
its state support.
Neither choosing this reference nor realizing the channel identifies
seconds, a geometric modular boost or an observed physical field.

<a id="the-interacting-detector-and-kick-obey-the-same-gibbs-law"></a>

### The interacting detector and kick obey the same Gibbs law

`product_detector_kms_control` applies the existing
[finite Gibbs/KMS theorem](../../Lean/QFT/StructuralInheritance.lean)
to this interacting preparation and its actual detector. Its bit effect is
`A=K0* K0-K1* K1=-sin(Q)`; its preparation kick is `exp(i epsilon Q)`.
For energy eigenvalues E_i and Gibbs weights p_i, put
`z_ij=A_ij Q_ji` and `omega_ij=E_j-E_i`. These signed cross weights obey

```text
W_AQ(t) = sum_ij p_i z_ij exp(-i omega_ij t),
p_j z_ij = exp(-beta omega_ij) p_i z_ij,
W_AQ(t-i beta) = Tr(rho Q A(t)),
d_epsilon <A(t)>|0 = i Tr(rho [A(t),Q]).
```

Thus the symmetrized cross-noise coefficient `(p_i+p_j)z_ij/2`
and commutator coefficient `(p_i-p_j)z_ij` have ratio
`coth(beta omega_ij/2)/2` when the denominator is nonzero. At zero
frequency the undivided detailed-balance identity applies. No independent
noise parameter enters; this is a cross-response to Q, not an A-generated
kick or a probability distribution of cross weights.

The comparison reuses the existing fixed-interval native-sign histories;
the homogeneous incident controller has a fixed positive kick. These
complete source histories with kicks `epsilon=+-1/4` give first-read
central slope `-0.0345168856` at model time `1/16`; the same-H ideal slope
is `-0.0345168861` and its infinitesimal response is `-0.0347883501`.
The finite-kick remainder is bounded by
`4 epsilon^2 ||Q||^3/3`, and complete history half-trace error delta adds
`2 delta/epsilon` to the central slope. If the raw prepared purification
has vector error e, its unnormalized marginal has trace error at most
`e(2+e)`; ordered cross correlations and response acquire allowances
`e(2+e)||Q||` and `2e(2+e)||Q||`. The retained source amplitudes are not
renormalized. Equilibrium cross-correlation measurements are algebraic
predictions here, not extra executed pointer protocols. A free Gibbs
preparation under the interacting clock violates the displayed KMS boundary.

The accessible algebra also preserves first-mode parity. Reweighting its
two central sectors by positive constants preserves restricted modular
flow and detailed balance, while changing accessible populations and
cross noise. The control changes the even-sector mass by `0.15320` and
the equal-time cross correlation by `0.001384`. KMS therefore does not
select these central weights. Nor does controller stationarity select the
Gibbs state: on the absorbing halt sector, every work/reference density is
fixed by the same instrument. These finite counterexamples leave Hamiltonian,
temperature and preparation selection as source-law choices.

<a id="a-held-out-record-law-after-central-weight-calibration"></a>

### A held-out record law after central-weight calibration

The [complete-record inverse](../common_source/DERIVATION.md#central-weight-from-complete-records) identifies the
even-sector weight w within `rho_w=w rho_even+(1-w)rho_odd`, with the same
`H`, `beta=1/4`, kicks, pointer and feedback fixed. Its calibration event A
contains terminal bit patterns `(0)`, `(0,0)`, `(1,1)` at model interval
`1/16`; all seven patterns, including native-clock timeouts, enter the trial
count. A distinct validation event B is equal bits in the existing complete
two-read program at interval `1/8`. That validation has four unconditional
outputs and no native timer. It does not reuse a selected subset of the
calibration trials or identify timestamps between these two protocols.

For sector probabilities c and v, linearity gives

```text
p_A = c_odd + w (c_even-c_odd),
p_B = v_odd + (v_even-v_odd)/(c_even-c_odd) (p_A-c_odd).
```

The certified validation contrast exceeds `0.00363067148365`; the affine
slope is about `23.99193`. `product_central_holdout_interval` propagates the
exact sector enclosures through the inverse interval, then adds validation
service and sampling allowances. For exact coefficients, calibration error
e_A and validation error e_B give the simpler bound
`|p_B_observed-p_B_predicted| <= e_B+|slope| e_A`. Each e includes its
service, preparation and sampling error. Fresh independent IID samples of
sizes N_A,N_B with radii r_A,r_B have joint failure probability at most
`2 exp(-2 N_A r_A^2)+2 exp(-2 N_B r_B^2)`.

The same source preparation gives `p_A=0.1569589724` and held-out
`p_B=0.6900830597`, retaining every work/reference output. The synthetic
count example uses `10^11` calibration and `10^9` validation trials as a
requirement, not collected data; it predicts `[0.68877855,0.69138758]` with
the declared service bounds and conservative sampling failure below
`3.82e-6`. This conservative sufficient budget concerns this calibration
menu; it is not a necessary sample lower bound for other admitted readouts.
Changing only the validation interaction, temperature, feedback or interval
violates that interval in the retained controls. No validation parameter
is fitted. The result predicts classical record probabilities within this
fixed family; it reconstructs neither unknown reference correlations nor a
physical state, clock dictionary or observed Universe from modular flow.

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
