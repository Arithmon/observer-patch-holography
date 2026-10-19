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

## Run

From the repository root:

```sh
# Replay the retained exact predictions, statistical checks and endpoint audit.
python3 -m model.predictions --check model/predictions/results.json
python3 -m model.predictions.bell_records --check model/predictions/bell_results.json
python3 -m model.predictions.optical_clock --check model/predictions/optical_clock_results.json
python3 -m model.predictions.regional_optical --check model/predictions/regional_optical_results.json
python3 -m model.predictions.modular_optical --check model/predictions/modular_optical_results.json

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

Using the unchanged complete UTC-grid block summaries from the
[public Ar13+/Yb+ data](https://zenodo.org/records/6901524), the conservative
conditional bounds are `10^-121` for 36Ar and `10^-22` for 40Ar. These are
retrospective tests of the entire attachment conjunction. The
[experiment](https://arxiv.org/html/2205.13053v2) uses averaged servo outputs
and a roughly thirty-second servo response; this packet does not derive that
readout from native repairs. Consequently it excludes neither the experiment's
atomic mechanism nor every OPH clock. Identical regional event sets satisfy
`B(N)=N` exactly and are outside the distinct-pair test; calling one event set
by two region names supplies no independent atomic stability explanation.
The source theorem and this conditional comparison derive no minimum number
of atoms, independent replicas or units of quantum coherence.

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
A common unit convention rescales both frequencies together. The raw
compiled arrays receive an event-weight allowance of twice their complete
half-diamond error, retain padding and are not renormalized.

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
