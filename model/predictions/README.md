# Fixed observer predictions and complete-count comparison

This package distinguishes two specified finite observer instruments using
an interference measurement. Both apply the same native repair to classical
loads. The current-record instrument C and contrast-record instrument E give
different probabilities when the source preparation is coherent. The settings,
sample sizes, error allowances and decision rule are fixed before comparison.

The package also audits the inherited P-dependent electromagnetic endpoint
against published data and replays a [public Bell-record comparison](BELL_RECORDS.md).
The latter constrains a local-factorized extension of the diagonal recorded
branch, under the experiment's stated setting and locality contract. These
retrospective calculations and the prospective finite-instrument protocol
have different evidential roles.

The [complete repair identification](../repair_identification/README.md)
extends the two-coordinate phase probe to all nine real instrument
coordinates. The [compiled-history forecast](../compiled_clock_history/README.md#calibrated-repair-law-and-complete-readout-predictions)
then carries calibration uncertainty through a complete two-tick experiment
with feedback and source-emitted apparatus. Under the stated synthetic
probability tubes, its two intervention settings retain a resolved joint
clock/detector contrast. The [resource-clock forecast](../laboratory_clock/README.md)
additionally predicts deadline-conditioned detector events under a joint
repair/pulse-duration calibration region. These are conditional source
predictions, not new discrepancies with SM + GR, externally registered
predictions or physical confirmations. The
[public-data audit](../../docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md)
tests a direct event-count/optical-phase identification against public
Ar13+/Yb+ records. Its complete-block range test charges missing data and
calendar-window selection. It rejects that simple attachment under its
explicit bounded-error and readout assumptions; collective clocks and
interacting physical identification require their own dictionary. The
[regional extension](#collective-regional-optical-attachment) derives all
shared-event collective clocks on one carrier and tests every distinct pair
under the same direct readout hypothesis. No clock definition needs
external approval, and no atomic transition is derived by choosing a unit.

The [simultaneous source-clock limit](../relative_observer_geometry/README.md#incident-clocks-on-the-simultaneous-density-continuum-scale)
also predicts joint clock jitter and density observations. For fixed bounded
incidence patterns, centered incident counts have Gaussian finite-dimensional
covariance `min(t,u)*c_ij`, where `c_ij` counts shared seams. Adjacent degree-five
ports have limiting correlation `1/5`; the degree-six endpoints of an external
seam have `1/6`. Conditional native stirring proves asymptotic independence
from density snapshot vectors, while exact finite nonlinear clock–field
dependence is nonzero. Uniform inverse-clock bounds let the same local clocks
read the proved continuous density paths. These predictions use the declared
stationary preparation, IID seams and observation access. Public optical ratios
do not supply the joint density readouts or identify this source with atomic
physics; no additional empirical verdict follows.

## Run

From the repository root:

```sh
# Replay the retained exact predictions, statistical checks and endpoint audit.
python3 -m model.predictions --check model/predictions/results.json
python3 -m model.predictions.bell_records --check model/predictions/bell_results.json
python3 -m model.predictions.optical_clock --check model/predictions/optical_clock_results.json
python3 -m model.predictions.regional_optical --check model/predictions/regional_optical_results.json
python3 -m model.predictions.modular_optical --check model/predictions/modular_optical_results.json
python3 -m model.predictions.repair_response --check model/predictions/repair_response_results.json

# Produce an unfilled count template; the destination must not exist.
python3 -m model.predictions --template /tmp/observer-counts.json

# After independent acquisition supplies every count, evaluate the file.
python3 -m model.predictions --evaluate /tmp/observer-counts.json

# Exercise the decision rule with explicitly labeled synthetic data.
python3 -m model.predictions --demo C
python3 -m model.predictions --demo E

python3 -m pytest -q -W error model/predictions
```

The template contains `null` counts and cannot pass validation. Every setting
requires exactly 8,192 attempted trials, with all twelve common
`(coin,current,outcome)` bins and an abort count supplied explicitly. Counts
must be nonnegative integers, including explicit zeros. The evaluator rejects
missing or duplicated bins, missing trials, a different protocol hash and
undeclared fields. It never divides by the number of successful trials.

## Fixed apparatus and predictions

Use the [native observer's](../observer/README.md) 78-dimensional load carrier:
twelve nonnegative integer port loads with total load two. Ports use indices
0 through 11. Let `x={2,3}` mean one load at port 2 and one at port 3, and
let `z={1,2}` have the analogous meaning. Prepare

```text
|phi_theta> = (|x> + exp(i theta)|z>)/sqrt(2).
```

Apply one repair on seam `(0,1)`, retaining its record. Measure the binary
source projector `P_plus=|phi_0><phi_0|` and its complement. Readout acts
inside each committed-record block; it does not erase or interfere different
classical record values. Fresh preparations supply independent trials.

| Fixed preparation | C: probability of plus | E: probability of plus |
| --- | --- | --- |
| Phase 0 | `5/8` | `3/8` |
| Phase pi | `1/8` | `3/8` |
| Phase pi/2 | `3/8` | `3/8` |
| Equal incoherent mixture of x and z | `3/8` | `3/8` |

The phase-reversal contrast is `1/2` for C and zero for E. All four
preparations have identical classical load probabilities. Consequently a
classical endpoint comparison alone cannot distinguish these instruments.

The two columns are special cases of one formula. For every quantum record
law compatible with the native seam rule, the plus probability is
`3/8 + Re(g e^{-iθ})/2`, where `g` is the law's coherence between the empty
seam of `x` and the waiting single load of `z`; C has `g=1/2` and E has
`g=0`. The [record-law classification](../record_laws/README.md) proves this
and gives the other distinguished laws: recording only the current gives
`g=1/√2`, and recording the full seam occupation gives `g=0`. The aggregate
plus frequencies of the four settings therefore estimate `Re g` and `Im g`
for any admissible law, while the locked decision rule below compares the two
executed laws.
[DERIVATION.md](DERIVATION.md) gives the complete joint law and the exact
channel calculation.

The [instrument code](observer_instrument.py) derives these probabilities
from the native endpoint map. It exposes every record, outcome and
unnormalized terminal source density, including zero-probability events.
E retains its additional contrast label internally. The shared comparison
interface sums probabilities over that label; this does not recover C's
quantum coherence.

## Error allowances and verdicts

[protocol.json](protocol.json) binds the predictions to the implementation
hashes. For every setting the plan fixes:

- 8,192 independent attempted trials, with no optional stopping;
- a sampling half-width of `1/32`;
- an independently established apparatus probability-error budget of `1/64`;
- a familywise sampling-error requirement below `0.001`.

Four settings each contribute twelve joint-bin frequencies, the aggregate
plus frequency and the abort frequency: 56 statistics. An exact rational
certificate bounds the simultaneous sampling failure probability by less
than `0.000012606`. Each observed interval is compared with the candidate's
ideal probability widened by the apparatus budget. One disjoint interval
rejects that candidate under the declared contract. Testing C and E against
the same simultaneous intervals needs no additional multiplicity factor.

The apparatus budget covers preparation, repair, storage and readout errors
in each scored probability. It is a physical premise requiring evidence
independent of the comparison counts. If its calibration guarantee has
failure probability beta, add beta to the stated sampling bound for an
unconditional error guarantee. A count file cannot establish that calibration,
trial independence or its own provenance.

An abort is an attempted trial without a complete valid record and readout.
It counts in the denominator and has ideal probability zero. A statistically
resolved excess over its allowance returns `apparatus_budget_exceeded`, so
the evaluator cannot attribute the discrepancy to the ideal repair law.

Without supplied counts the verdict is `no_data`. Software demonstrations
are marked `synthetic`; external counts are marked `external_unverified`.
Both report `empirical_confirmation=false` and `provenance_verified=false`.
If one candidate survives and the other is rejected, the counts discriminate
the conditional instruments. If both survive, this test does not separate
them. If both are rejected, neither specified apparatus model explains the
counts within the declared premises. Survival means `not_rejected`, not proof.

## Physical scope

The native load interpretation, coherent preparation and source-projector
measurement are explicit apparatus choices. Both C and E are ordinary
quantum instruments. Their separation is a test of the chosen record law,
not a novel prediction distinguishing OPH from quantum mechanics.

Neither finite instrument result has a valid comparison with publicly
available physical data. The C/E protocol requires a specific coherent
12-port preparation, native fair-coin repair, committed coin/current record,
source-projector readout and complete attempted-trial counts; no matching
dataset or calibrated implementation is supplied. The
[incident-clock field construction](../observer_time_fields/README.md) predicts
local laws in dimensionless incident-repair counts, while the
[screen comparison](../native_action/INCIDENT_REFINEMENT.md) is conditional on
source port maps, matched initial local/reference states and a particular
incident-event record interface. There is no physical calibration connecting
those quantities to published laboratory time, position or field readouts.
The ideal modular-phase target is a supplied reference in a separate finite
model, and the CODATA endpoint comparison below is retrospective; neither is
an empirical test of these incident-clock results or the C/E protocol.

## Collective regional optical attachment

[regional_optical.py](regional_optical.py) extends the direct optical test to
every pair of distinct event sets generated by nonempty port regions on one
native carrier. The [collective clock theorem](../laboratory_clock/README.md#collective-clocks-on-the-same-carrier)
derives the full renewal law, its mean ratio `b/a`, and its residual fluctuations.
No independent copies, noise parameter or atomic model is introduced. The
tested attachment is

```
reported 300-second block ratio = gain*a*B(N)/(b*N) + error,
gain >= 1,    |error| <= 10^-12.
```

Here `N` is one fixed integer number of A ticks per block in the declared
`300*[642121496772644,642121496772646]` interval. The frequency bracket is a
hypothesis, not a hard bound inferred from a quoted standard uncertainty.
Each consecutive candidate calendar block, including every unobserved block,
must contain exactly N consecutive A intervals. This makes the block count
laws independent under the declared IID scheduler.
The normalization `a/b` follows from the count law and permits unequal native
clock rates. The clocks require ordered regional tape unions and coincidences.
The hypothesis includes no additional temporal filtering.

Rejection uses a point-probability bound, not variance alone. If the geometric
part has failure probability `q`, it is the sum of an independent
`Bernoulli(q/(1+q))` and twice a geometric variable with failure `q^2`.
The larger variance `v` of this Bernoulli and the terminal `Bernoulli(c/a)`
is strictly positive for every distinct event pair. Convolution and Fourier
inversion give `sup_k Pr(B(N)=k) <= sqrt(pi/(8*N*v))`: use
`|1-r+r*exp(it)|^N <= exp(-2*N*r*(1-r)*sin(t/2)^2)` and
`sin(|t|/2)>=|t|/pi` on `[-pi,pi]`, then enlarge the integral to the real line.
The code rounds outward using `pi<22/7` and rational arithmetic.

For an observed range `w`, the admissible B-count interval contains at most
`ceil((b/a)*N_upper*(w+2*10^-12)/gain_lower)+1` integers. All 2,848 incidence
classes are checked; the weakest bound occurs at `(a,b,c)=(5,30,5)`.
If any K of T candidate blocks lie in an interval of mass at most `p`, anchor
its smallest value and union over the other K−1 values. This gives
`T*choose(T-1,K-1)*p^(K-1)`. Arbitrary block removal and a fitted constant gain
or interval centre are allowed. The calculation additionally charges
`105120^3` for selecting calendar endpoints and K inside 2021, and all
11,782,056 ordered distinct event-set pairs. It covers neither searches in
other years nor other block durations.

Using complete UTC-grid summaries from the [public Ar13+/Yb+ data](https://zenodo.org/records/6901524),
the direct regional bounds are `10^-121` for 36Ar and `10^-22` for 40Ar.
A separate conditional extension averages eight **independent native replicas**
through a preselected causal, nonnegative, unit-integral linear filter with
**compact support `[0,30 s]`**. A ticks are equally spaced in the declared
calendar coordinate. Interior `[0,270 s]` ticks contribute their full filter
mass to one 300-second block. Conditioning on every noninterior increment
leaves independent block innovations with arbitrary fixed outside shifts;
the same anchor-and-subset bound applies. Binary geometric digits and a
split Fourier integral give primitive-port bounds `10^-123`, `10^-38` and
regional bounds `10^-123`, `10^-25` for 36Ar, 40Ar respectively. The regional
worst class is `(29,30,29)`. Filter selection from the measured values is not covered.

The [experiment](https://arxiv.org/html/2205.13053v2) reports eight servo
outputs and a response time constant around 30 seconds. Those facts do **not**
establish independent native replicas, equally spaced native ticks, or a
compact impulse response: a decaying filter can retain an infinite tail.
The extension rejects its stated conjunction, not the actual atomic servo.

A separate **moment-matching diagnostic** equates the unfiltered, zero-error
population variance to the sample Allan value. It requires at least
`2544097785171`, `2407037601663` independent primitive replicas, or
`53002037192`, `50146616702` at the finite carrier's regional floor, for
36Ar, 40Ar respectively. These are conditional algebraic targets, not
confidence bounds on physical resources; the sample Allan value is not a
certified population upper bound. No continuum refinement level follows.
The refining source has `M_l=630*4^l` seams and `a_l=1260*16^l` attempts per
density-coordinate unit, hence a fixed degree-d port has **expected**
`a_l*d/M_l=2*d*4^l` incidents. Its seconds calibration remains free, and the
thirty-seam regional floor does not transfer to arbitrary refined regions.
Identical regional event sets satisfy `B(N)=N` and remain outside the test.

The clock and field comparison also reuses the existing
[pair-incidence reconstruction](../relative_observer_geometry/README.md#inferring-nonuniform-scheduling-from-clock-comparisons)
in a directly testable moment form. Let `I` mark all ports incident to one
source event, `p=E[I]`, and `Sigma=Cov(I)`. Fair repair gives
`P=Id-diag(p)+(Sigma+p*p^T)/2` for the mean density transfer. Every event is
present, including zero current, and ordered coincident tapes supply the
index. The optical module checks this identity on all 3,960 native branches;
no stationary load preparation or second field rate is needed. For a fixed
window of N IID events, estimate `Q=E[I I^T]`, reusing its diagonal for p.
Chebyshev on its 78 distinct entries gives operator error at most `6*delta`
with failure at most `78/(4*N*delta^2)`. This is a finite sampling certificate
for the existing clock-to-field correspondence, not an independent physical
derivation. Public ratio files lack these joint clock/occupation records.

The protocol uses event order: preparation, one repair, then readout. It
requires no conversion of a solver counter to physical seconds. A laboratory
implementation must nevertheless implement and calibrate those operations.
No spatial continuum, gravity or universal physical clock is inferred from
this finite test.

## Recovering the modular clock from visible current records

[inverse_clock.py](inverse_clock.py) ties a current record to the coherent
readouts of the [same three-mode transfer](../transfer_modular_probe/README.md).
Its input is the complete movement-bit tape `B_t=|j_t|`, **including every
zero-current attempt**, under the uniform IID seam schedule and fair coin on
the binary total-two carrier. This is an admitted complete-attempt record
service, stronger than an isolated port's incident tape. Record position
defines its clock; no physical seconds conversion enters the construction.

Write `p=E[B_t]` and `c_n=E[B_t B_(t+n+1)]-p^2`. Thus `c_0` is the
covariance of consecutive records, not the same-time Bernoulli variance.
The actual marked strong lumping gives
`c_n=w_+ lambda_+^n+w_- lambda_-^n` with positive weights and distinct
positive roots. Four consecutive covariances invert this relation:

```
D = c_0*c_2-c_1^2 = w_+*w_-*(lambda_+-lambda_-)^2
S = (c_0*c_3-c_1*c_2)/D = lambda_++lambda_-
P = (c_1*c_3-c_2^2)/D  = lambda_+*lambda_-
lambda_± = (S ± sqrt(S^2-4P))/2
w_+ = (c_1-lambda_-*c_0)/(lambda_+-lambda_-)
w_- = (lambda_+*c_0-c_1)/(lambda_+-lambda_-)
```

Full 66-state marked path sums, computed independently of the forward
spectral formula, give `p=5/33`, `D=1/9702990000`, `S=5/3` and `P=311/450`.
Consequently the observed two-mode domain determines
`lambda_±=(25±sqrt(3))/30` and `E_±=-log(lambda_±)`. Their unit-free ratio
is `E_-/E_+≈2.2033420301466635`. The one-exponential approximation matching
`c_0,c_1` misses `c_2` by exactly `1/2673000`; matching an average rate or
one decay does not establish the joint prediction. This approximation is
a moment control, not a claimed alternative native instrument.

`identify` accepts rational probability intervals and encloses every
compatible model in this positive two-mode domain. It verifies positive
Hankel determinant, separated real roots in `[1/2,1)`, and positive residues;
failure produces an inconclusive result. Repeated uncertain parameters are
shared in the underlying model; outward interval arithmetic enlarges their
image. A successful enclosure does not certify that an arbitrary input box
contains a realizable model, or rule out additional hidden modes outside the
declared domain. It also does not distinguish C from E, whose classical
current tapes agree.

`predict_readout` uses those same root intervals to predict the existing
coherent code's X and Y quadratures: `q_X=(1+cos(s E_j))/2` and
`q_Y=(1+sin(s E_j))/2`. Preparation in a specified spectral mode, the
coherent phase operation and the detector basis are declared services;
passive repair is not identified with this unitary. The retained example
uses `s=4`, both modes, and all four outcomes including padding. Synthetic
source-centred input tubes of radius `1e-12` predict the four zero-address
probabilities within approximately `[.946797,.948659]`, `[.721645,.723507]`,
`[.762113,.764252]` and `[.924060,.926199]`; the compiled source outputs lie
inside the corresponding outward bounds. The same inference predicts the
unused lag `c_4`. No separate field frequency or clock scale is fitted.
A common unit convention rescales both frequencies together. The raw arrays
retain padding without normalization and receive twice their complete
half-diamond error as an event-weight allowance.

The probability tubes have a dependent-tape coverage contract. Use the same
`M` starting positions from `M+5` retained bits to estimate the mean and five
lag products, reserving the fifth product for validation. A bounded function
of `L` consecutive marked events has variance at most `1/4`. For separated
windows whose starting positions differ by `k>=L`, conditioning on their
endpoint states and the reversible marked lumping gives absolute covariance
at most `(9/10)^(k-L)/4`, since its nonconstant eigenvalues are below `9/10`.
Overlapping windows retain the `1/4` bound. Summing covariances gives
`Var(mean_L) <= (2L-1+2/(1-9/10))/(4M)`. For the six statistics (`L=1..6`),
Chebyshev and a union bound give simultaneous tube failure at most
`39/(M*delta^2)`, without independent-window or independent-validation data
assumptions. For an arbitrary initial binary density, discarding a prefix
of 221 native attempts adds at most `1e-3` by the existing
[unread warm-up theorem](../record_laws/README.md); no retained controller
record is used to select the subsequent tape.

At `delta=1e-12`, this conservative proof gives total failure at most `.01`
for `M=4333333333333333333333333334` windows. This is a sufficient bound,
not a minimal data requirement, a collected dataset, or a feasible
experimental protocol. Broad tubes can be inconclusive. The retained
receipt exercises algebraic inference, interval coverage and complete
source readouts; the public optical data below do not supply this marked
current tape or the coherent-mode measurement dictionary. An atomic or
continuum-field identification must preserve these joint constraints before
those observations can test this correspondence.

Replay with `python -m model.predictions.inverse_clock --check
model/predictions/inverse_clock_results.json`.

### Preparing the finite edge reference from source pulses

The same replay implements an active preparation for the
[native edge generator](../observer/COHERENT_CONTINUATION.md#a-finite-edge-generator-representation-from-native-repair).
Its density transfer is `K12=I-L/60`; `L` has eigenvalues
`0, 5-sqrt(5), 6, 5+sqrt(5)` and sector ranks `1,3,5,3`. A closed-history
trace represents a Gibbs functional but does not prepare a coherent density.
Here a supplied protected register in `I12/12`, a declared event preparation,
and an explicit herald produce the density itself.

Prepare each of the thirty even event addresses with probability `(2-p)/60`
and each odd address with probability `p/60`. The existing compiled reversible
shift fixes an even address and swaps the corresponding seam endpoints at an
odd address. Unprepare the event register and herald a valid site address,
event zero and clean scratch.
The success operator is exactly
`F_p=(1-p)I+pK12=I-pL/60`. The retained instance uses `p=1/2` and two
successful heralds, producing `rho=F_p^4/Tr(F_p^4)`. The source routines are
executed pulse by pulse on all columns, then reused on the fifteen-qubit
register for all twelve logical input columns. Every first-failure output,
including address padding and dirty scratch, contributes to the total
probability. The complete raw filter is reused for the second herald.
Fixed numerical error guards test those arrays; they are not interval or
physical precision certificates. An executed source rotation
`exp(i*pi*Y_01/4)` followed by all sixteen address outcomes reads the
normalized zero-minus-one contrast `260343/43825510`; dephasing the prepared
density reduces this contrast to zero. Sector projectors below are
mathematical readout predictions, not a compiled sector detector.

For a fixed preparation parameter `t>=0`, set `p=30t/n<=1`. After `n`
successful heralds,

```
rho_n = (I-tL/(2n))^(2n) / Z_n  ->  sigma_t = exp(-tL)/Z
p_R(n,t) = d_R (1-t lambda_R/(2n))^(2n) / Z_n
Pr(all n heralds succeed) = Z_n/12 >= 1/12.
```

These relations tie preparation and sector readouts to the same generator,
with no independently fitted sector temperatures. Since `lambda_max<8`,
`0 <= -2n log(1-t lambda/(2n))-t lambda <= e_n`,
where `e_n=16t^2/[n(1-4t/n)]`. The finite weights are the Gibbs weights
multiplied by factors in `[exp(-e_n),1]`. Rejection coupling gives trace
distance at most `min(1,e_n)`; the normalized log ratios give
`D(rho_n||sigma_t)<=e_n`. The full-register modular gap inferred from
`log[(p_S/d_S)/(p_R/d_R)]` differs from `t(lambda_R-lambda_S)` by at most
`e_n`; it is fixed by the same preparation and sector law. Thus the Gibbs variational identity
`S(rho)-t Tr(rho L)=log Z-D(rho||sigma_t)` supplies an explicit
free-entropy deficit and stationarity at the limit on the trace-one fibre.
This is not a derivation of bulk-plus-area entropy stationarity.

The density has nonzero off-diagonal entries: dephasing preserves site
probabilities but changes the prepared reference. Its sector weights use
`d_R`, not the regular-group multiplicities `d_R^2`; the represented algebra
has dimension 44, with no four-dimensional irrep. Full group sewing and a
compact-group limit do not follow. The edge modular generator is `tL` up to
a scalar on the full logical algebra `B(C^12)`, which includes operators
between sectors. On the represented algebra `M1+M3+M5+M3`, the state is
central and its modular flow is the identity. The full-register generator
is distinct from the wave generator `sqrt(L)`; one clock conversion
cannot identify their gaps across all bands. The ideal growing filter has
the proved limit, while source precision and service budgets for that entire
family require a separate bound. Biased coherent preparation, protected
storage, initial density and herald readout are admitted observer services,
not a state selected by passive repair. Public clock-ratio datasets do not
contain this joint preparation and sector measurement, so these are
conditional model predictions, not an empirical discriminator against SM
or GR.

## Native repair response from passive clock and field records

[repair_response.py](repair_response.py) predicts a finite response using the
same repair law, fields and incident clocks that define the passive process.
There is no fitted force, response coefficient or conversion to seconds.
Prepare the uniform fixed-number binary sector. Disclose the occupations at
ports 0, 1 and 2. If the port-0 bit is one, apply one fair native repair on
seam `(0,1)`; otherwise apply no repair. The trigger, coin and signed current
are retained. Postpulse fields are reconstructed from the disclosed bits and
that current, without another state query. The pulse contributes one tick to
each incident observer clock when it occurs, including a zero-current pulse.

Write `a=n_0`, `b=n_1` for the initial occupations in a *passive continuation*
experiment. The corresponding post-intervention law has density

```text
w = 1 + (b-a)/2,
Pr_feedback(F) - Pr_passive(F) = E_passive[(b-a) 1_F]/2.
```

Here `F` is any fixed subsequent record event, with the same continuation
rule in both experiments, independent of the earlier pulse mark. It may
include fields, incident-clock comparisons, feedback on subsequent records,
stopping and timeouts. A continuation that consults the pulse mark must use
the separate marked weights below, not the summed density. The identity follows from detailed
balance of the local fair repair: for its self-adjoint operator `R`,
`R a=(a+b)/2`, and the pulse policy has operator `I + a(R-I)`. This is the
finite Markov adjoint identity applied to the actual native intervention,
not a claim of a new general fluctuation-response theorem. It is exact at
the specified pulse amplitude; it is not an infinitesimal-response fit.

The complete marked prediction uses these passive-record weights:

| Retained trigger, coin, current | Weight |
| --- | --- |
| No pulse (`0,-1,0`, with `-1` the unused-coin label) | `1-a` |
| Pulse (`1,0,1`) | `(1-a)*b/2` |
| Pulse (`1,1,0`) | `a/2` |
| Pulse (`1,0,0`) | `a*b/2` |
| Each other declared pulse mark | `0` |

These are importance densities for unnormalized joint events, not
conditional probabilities of drawing a mark from each passive sample.
Their sum need not equal one at each sample. Their equilibrium averages
produce a normalized complete law. No division by that sum is permitted.
Pre-pulse disclosed fields can be recovered from the postpulse fields and
retained current, so the prediction preserves their record as well.

The executed forecast stops at the first subsequent port-0 incident event,
with a supplied watchdog cap of one or two uniform native attempts. The
second observer counts incidents at port 1 or port 2. Every signed-current
record meeting the three read ports, zero-current events, pulse marks and
timeouts is included. The exact timeout probability is `(5/6)^cap`.
For the event "the stop succeeds, the final port-1 field is one, and the
second clock has at least one *subsequent* tick", the cap-two predictions are:

| Second clock | Passive probability | Feedback probability | Response from passive covariance |
| --- | --- | --- | --- |
| Port 1 | `71/5400` | `34/2475` | `7/11880` |
| Port 2 | `1/216` | `5/792` | `1/594` |

The cap-one response is exactly zero for both clocks. Total clock readouts
also retain the pulse tick: subtracting that known contribution in the
specified event merely selects equal continuation menus for comparison.
The watchdog attempt count supplies no physical duration. Independent
primitive-path expansion checks all one-attempt record labels; exact native
propagation checks the capped laws. An unconditional extra repair has zero
marginal response at equilibrium. A nonstationary preparation violates the
stationary weight formula; the retained control has predicted total mass
`1/2` against actual mass one and is rejected without renormalization.

For an arbitrary finite graph with `V` ports, `M` seams, fixed total `k`,
`chi=k*(V-k)/(V*(V-1))`, and native density transfer `K=I-L/(2M)`, the same
identity gives the complete mean occupation response to a pulse on `(a,b)`:

```text
Delta <n_c(n)> = chi/2 * [(K^n)_(b,c) - (K^n)_(a,c)].
```

For any chosen spatial normalization `s>0`, let `r=2*M*s/(2*M*s+z)`, `z>0`.
Then `(1/(2*M*s))*sum_(n>=0) r^(n+1)*Delta<n(n)>` equals
`chi/2*(e_b-e_a)^T*(s*L+z*I)^-1`. Full 66-configuration inverses independently
verify this response on the native carrier. This makes the
[spatial Green operator](../observer_time_fields/REFINEMENT.md#native-temporal-correlations-and-the-spatial-green-field)
a response to an observer's explicit repair operation, as well as a passive
correlation. One microscopic pulse does not by itself have a nonzero
continuum normalization. The following controlled source supplies one.

### A controlled continuum mean with an executed-event clock

Fix a real port profile `h`, disclose the initial binary load and reconstruct
subsequent occupations from retained currents. In each selection slot choose
a uniform seam `(a,b)` and execute its unchanged fair native repair with
probability

```text
p_execute = [1 + epsilon*(h_b-h_a)*(n_a-n_b)]/2.
```

Every selection, idle/execution flag, repair coin and current is retained.
The source certifies `|epsilon*(h_b-h_a)|<=1/2`. The designated dynamical
clock counts executed repairs, including zero-current repairs. It is an
aggregate event clock; it is not a primitive port-incident clock. Selection
slots and controller service have separate counts. Native fair bits can be
read by repairing a disclosed equal-occupation seam: its load stays fixed,
but every service write remains a charged native event. An equal pair exists
in every binary carrier triangle. These resource counts do not identify
seconds, an atomic transition or an autonomous force law.

For the uniform fixed-number preparation the execution subkernel has
stationary incoming mass `pi/2`, even at nonzero epsilon; its row activity
depends on the load. Dividing by that activity gives the actual kernel at
successive executed repairs. Differentiating this normalization at zero
produces the extra `K` in the exact response:

```text
executed repair count r: d< n >/d epsilon|_0 = 2*chi*K*(I-K^r)*h,
selection slot count q: d< n >/d epsilon|_0 = 2*chi*(I-Q^q)*h,
Q = (I+K)/2.
```

Full 66-configuration propagation checks both identities. Replacing `q`
by `2r` is not exact: the retained finite controls have different responses.
This correction retains the correlation between clock stopping and load.

On the actual level-`l` carrier set `s=4^l`, `N=240s`, `M=630s` and density
one third, with its canonical susceptibility `chi`. For a fixed mean-zero
flat coordinate `f`, put `h_l=cell_average(f)/sqrt(N*chi)`, copied to the
ports of each cell, and read
`eta_l(g)=sum_p cell_average(g)_p*(n_p-1/3)/sqrt(N*chi)`.
The [common native generator](../observer_time_fields/REFINEMENT.md#the-full-native-dynamics-have-the-same-limit)
then gives, along the same spatial subsequence, at
`r=floor(2*M*s*tau)=floor(1260*16^l*tau)` executed repairs,

```text
d E[eta_l(g)]/d epsilon|_0 --> 2 <g,(I-exp(-tau*A))f>.
```

The small-`tau` source is `2Af` in the dual of the closed form domain when
`f` lies in that domain. This derives a controlled mean response through
the same `A` as the passive field and heat limit. It does not identify a
physical force or quantum field interaction.

The statement also has a finite, level-independent amplitude. At `tau=1`
stop at `r` executed repairs or `q=4r` selection slots, retaining the actual
terminal field on timeout. Against the zero-bias law the stopped likelihood
is `L=product(1+epsilon*s_j)`, where `s_j` is the signed activation score
and is zero after stopping. Its conditional mean is zero and variance is
at most `v=h_l^T L_graph h_l/M`. The
[coordinate energy bound](../physical_space_time/OCCUPATION_FIELD.md#a-common-positive-spatial-generator-from-the-native-green-field)
gives `q*v<=8/chi<36`. Writing `S=sum s_j`, martingale conditioning yields

```text
E L^2 <= (1+epsilon^2*v)^q <= exp(epsilon^2*q*v),
E (L-1-epsilon*S)^2 <= (1+epsilon^2*v)^q-1-q*epsilon^2*v.
```

Indeed the remainder increments are `epsilon*s_j*(L_previous-1)`;
conditional centering removes their cross terms. When `epsilon^2*q*v<=1/2`,
Cauchy--Schwarz bounds the nonlinear mean error by
`36*epsilon^2*||g||_L2(mu)`, uniformly in refinement.

The cap comparison retains the stopping correlations. At zero bias the
uncapped slot count `T_r` has mean `2r` and
`P(T_r>4r)<=2^(-r/2)`. With `B` the maximum absolute terminal field, the
squared difference of capped and uncapped derivatives is at most
`8*B^2*r*v*P(T_r>4r)`: use the complete stopped score, whose second moment
is at most `2r*v`, and the fact that the fields differ only on timeout.
The biased timeout obeys
`P_epsilon(timeout)^2<=E_0[L^2]*P_0(timeout)`; it is not assigned the
unbiased watchdog probability. The executable propagates all likelihood,
score and field moments on the complete finite state space and verifies
these correlations, including a nonzero nonlinear mean correction.

For a flat coordinate of the unit carrier, `||f||<=1` and
`||f||^2>=3/16`. Since `A>=3/4096`,
`2<f,(I-exp(-A))f> >= 9/32792`. Consequently the fixed choice
`epsilon=10^-6` has the strictly positive bound

```text
liminf_l E_epsilon[eta_l(f)] >= epsilon*9/32792 - 36*epsilon^2
                            = 244359/1024750000000000 > 0.
```

This is a nonvanishing controlled continuum signal. It does not prove a
unique nonlinear law at fixed epsilon. A finite fair-bit implementation
rounds each activation probability down using
`b=ceil(log2(16*N*q^2))` bits. Coupling bounds the path discrepancy by
`q*2^-b`; the squared terminal-mean error is at most
`1/(64*q^2*N*chi)`, which vanishes. At most `b*q` native coin-service
writes are charged, plus the separately supplied uniform-seam selection
service. This implementation is a declared observer protocol, not a
claim that passive repair spontaneously chooses this drive or clock.

The [published Ar13+/Yb+ files](https://zenodo.org/records/6901524) contain
corrected optical frequency ratios. They do not contain these disclosed
occupations, trigger/current records and corresponding intervention arms.
Their [experimental description](https://arxiv.org/html/2205.13053v2) supplies
an atomic apparatus; this packet supplies no identification of its readouts
with the native pulse. These exact native forecasts therefore add no public
measurement agreement or discriminator against SM + GR. Choosing event-count
units is sufficient to define the model clocks; identifying an observed
physical response must preserve the joint relations above.

Replay with `python -m model.predictions.repair_response --check
model/predictions/repair_response_results.json`.

## Minimal modular-transition optical attachment

The source-derived [three-mode transfer generator](../transfer_modular_probe/README.md)
has energies `0`, `E_+=-log((25+sqrt(3))/30)` and
`E_-=-log((25-sqrt(3))/30)`. Its positive gaps are `E_+`, `E_-` and
`E_--E_+`. [modular_optical.py](modular_optical.py) tests the explicit
hypothesis that the measured Ar13+ and Yb+ optical transitions correspond
directly to two such gaps with **one common positive frequency scale**.
That scale cancels in their ratio. A common Gibbs inverse temperature also
cancels, since `-log(T3^beta/Tr(T3^beta))=beta H+constant`. The test uses the
real logarithm of the positive transfer; it does not add integer harmonics,
alternative logarithm branches, species-dependent gains or offsets.

All nine ordered assignments, including the three same-gap assignments,
are evaluated by exact outward intervals. Their distinct ratios are
approximately `0.453856`, `0.546144`, `0.831019`, `1`, `1.203342`, `1.831019`
and `2.203342`. King et al. report the dimensionless ratios
`R40=1.05776938758748094(11)` and `R36=1.05776646273518748(13)` in
[Extended Data Table 1 of the accepted manuscript](https://arxiv.org/pdf/2205.13053v2#page=20),
published in [Nature 611, 43–47 (2022)](https://doi.org/10.1038/s41586-022-05245-4).
These table values and their quoted standard uncertainties are retained in
[the small data record](data/modular_optical_ratios.json); this comparison
does not reanalyse the measured time traces.

Under the declared **total absolute ratio error at most `1e-12`**, every
assignment is disjoint from either observation. The closest ratio is one;
the certified absolute separation exceeds `0.057` for each isotope. The
error budget combines measurement and attachment errors and is part of
the tested hypothesis. Quoted standard uncertainties are not converted
into hard bounds or a significance claim, and the two observations are
compared separately without an independence assumption. A `0.06` total
ratio-error control does permit an assignment. Allowing separate arbitrary
frequency scales would also remove the test's constraint.

This is a retrospective failure of the minimal direct-transition assignment.
The source does not derive atomic structure, identify these ions with the
three-state code, or supply an optical detector mapping. It therefore neither
excludes all modular clocks nor supplies a new discriminator against SM + GR.
It records a concrete restriction on a proposed physical dictionary while
keeping the source spectrum and its observer readout independent of the
measured ratios.

## Inherited coupling comparison

[endpoint_audit.py](endpoint_audit.py) compares the pinned
`thomson_structured_running` output used by the
[P-dependent observer](../observer/p_response/README.md) with the
[CODATA 2022 value published by NIST](https://physics.nist.gov/cuu/pdf/wall_2022.pdf),
`alpha^-1=137.035999177(21)`. The inherited output is approximately
`136.994835177413`, a deficit of **300.388218 parts per million** relative
to that central value. Exact rational arithmetic makes its narrow numerical
interval disjoint from the central value plus or minus five reported
standard uncertainties.

The [P-derivation guide](../../code/P_derivation/README.md) states the
incomplete physical transport obligations. The numerical interval does not
bound those missing physical contributions. The audit rejects the additional
hypothesis that this output is the complete physical endpoint with zero
missing correction. It assigns no theory sigma significance, performs no fit,
and does not recalibrate P from the target. This is a retrospective endpoint
gap, not a refutation of OPH as a whole or a prospective prediction success.
