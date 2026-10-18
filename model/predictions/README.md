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
interacting physical identification remain open. No clock definition needs
external approval, and no atomic transition is derived by choosing a unit.

## Run

From the repository root:

```sh
# Replay the retained exact predictions, statistical checks and endpoint audit.
python3 -m model.predictions --check model/predictions/results.json
python3 -m model.predictions.bell_records --check model/predictions/bell_results.json
python3 -m model.predictions.optical_clock --check model/predictions/optical_clock_results.json

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

Neither new finite result currently has a valid comparison with publicly
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

The protocol uses event order: preparation, one repair, then readout. It
requires no conversion of a solver counter to physical seconds. A laboratory
implementation must nevertheless implement and calibrate those operations.
No spatial continuum, gravity or universal physical clock is inferred from
this finite test.

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
