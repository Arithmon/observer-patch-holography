# A shared-event count read from the compiled clock and field

This successor of the [source Ramsey join](RAMSEY_SERVICE.md) measures a
finite, operational clock comparison on the **unchanged C and E laws**. One
committed port-0 incident asks for a transported modular word at parameter
`s=4` and advances the protected native density wave by one step. The native
history is hidden from the readout. An observer receives only the terminal
`(1,5)` detector current and the compiled field Ramsey sign, and estimates
`N`, the number of port-0 incidents in the two-attempt policy. Every complete
native/clock/field record remains in the source instrument; the statistical
claim below is conditional on successful finite preparation and a terminal
joint read. The [campaign contract](clock_count_inference_campaign.json)
fixes that selection, the matched-marginal control, C/E, and predecessor
result hashes before execution.

## Exact uncertainty reduction

The exact complete target law gives `P(N=1)=307/360` and `P(N=2)=53/360` for
both laws. For binary `N`, the optimal squared-error predictor after an
outcome `y` is `1+P(N=2|y)`. Its remaining risk is
`sum_y P(N=1,y) P(N=2,y) / P(y)` in squared incident-count units. Rational
cosine enclosures from the existing parameter-four source certificate give:

| Readout | C risk | E risk |
| --- | ---: | ---: |
| No readout | 0.125547839506 | 0.125547839506 |
| Field Ramsey sign | 0.121297515512 | 0.121297515512 |
| Modular detector current | 0.096165175202 | 0.096161912593 |
| Both outputs | 0.093432667201 | 0.093428964635 |

These are display values. The [retained receipt](clock_count_inference_results.json)
contains outward rational intervals for each probability, posterior, risk and
gain. The certified joint gain over the detector alone exceeds `0.0027325`
for C and `0.0027329` for E. For the outcome `(current=-1, field=-)`, the
posterior probability of two incidents exceeds `0.87`; for
`(current=-1, field=+)`, it is below `0.62`. This is information about the
shared count even when its native tape is concealed from that readout.

The countercontrol independently resamples the field bit from its exact
unconditional marginal while retaining the complete clock marginal and
field marginal. It removes the extra count information: the combined risk
then equals the detector-only risk. This tests the common-event coupling, not
merely two separate changing marginals. The C/E difference is retained rather
than averaged away.

## Source execution and theorem attachment

The producer freshly runs the existing source clock frames on all admitted
input columns and the source-pulse Ramsey field factors for both counts. It
retains each source export-abort position and the three finite-preparation
failure flags before conditioning on a read. It checks complete finite mass,
the source-read table and its Bayes risk. The three predecessor receipts are
hash-verified against their full source closures, including the native-return
wave-reference construction; their result bytes and the new code, test,
contract and this document are bound into the new receipt. The exact risk
certificate uses the predecessor's full rational coefficient law, not a fit
to the newly executed floating outputs.

The faithful transported target on `H66` has a finite self-adjoint
`H_A=-log(sigma_A)`. [PrivateInner](../../Lean/Dynamics/PrivateInner.lean)
therefore supplies the specified inner modular action, while the inherited
source-word certificate establishes its implementation. This theorem does not
select that reference or create a physical time unit.

The existing [complete Kraus theorem](../../Lean/Dynamics/ChoiCPTP.lean)
applies because the predecessor has checked the native partial-injection
maps, proper-code complement, both Ramsey outcomes, every source abort, and
finite-preparation flags. It licenses composition and untouched-reference
extension under those supplied instrument premises. Here the source-clock
and compiled-field half-trace allowances sum to less than `4.41e-7` for the
full flagged output. With ideal read probability `5/33`, elementary
conditioning gives a conservative read-law total-variation allowance below
`5.83e-6`. Optimal binary squared-error risk is 1-Lipschitz under total
variation; the difference of the two risks therefore retains a strict gain
above `0.00272` for the executed finite-preparation readout. The receipt
keeps the exact fractions and source float checks separate. Its routine
floating guards are not interval proofs of hardware or of the numerical
backend. Source diagnostic values are displayed as nearest `10^-8` grid
centers with one grid unit of margin; the source guards act on unrounded
values.

The [native-return wave reference](../observer_time_fields/README.md#a-native-return-reference-for-the-coherent-wave-generator)
reconstructs a reference for the same `arccos(K)` density-wave generator.
That reference is not prepared or enacted as the clock in this instrument;
the transported clock uses a different faithful three-charge reference.
The native-return cap-512 generator allowance is about `0.03588`;
a direct two-step Duhamel bound is about `0.07176`, wider than this small
uncertainty gain. Thus that cap does not certify the gain if one *substituted*
its approximate generator for the active wave. This packet uses the actual compiled wave service and its
own bound instead. Its mean-zero 720-arc representation is intertwined with
the native 3960-arc density sector for the selected probe, not with every
arc observable or a regional full-field algebra.

[ClockReadout](../../Lean/Time/ClockReadout.lean) requires an individual
strictly monotone record reading; the stochastic detector outcome does not
inhabit that interface. [ClockComparison](../../Lean/Time/ClockComparison.lean)
requires two ordered shared events and a distinct third event to check an
affine calibration. Here `s=4N` is a programmed dimensionless request per
successful incident, and the displayed gain is a readout of `N`; it is not
a three-event affine-rate test. Host execution time, native attempts, modular
parameter, event ticks, and calibrated physical seconds remain separate.

## Complete resource boundary

The reused clock source budget permits at most two parameter-four words and
`6,180` clock/native source pulses in a successful two-attempt history. The
newly reused field compiler uses at most `285,030` pulses for two steps,
including preparation and Ramsey analysis. Thus the combined successful
path has at most `291,210` source pulses and `582,420` coherent processor
load/export calls, with the six-level processor shared. The quantum register
inventory is nine clock-address qubits plus sixteen field/Ramsey/scratch
qubits. The inherited compact field-program representation is `6,311,066`
UTF-8 bytes, and the core/repair controller has 132 policy states; program
storage, the FIFO, exported record tape and host workspace are accounted
separately. The native finite preparation retains its separate upper limit
of `659,460,352` attempted steps. This packet adds six classical readout
cells and two candidate count labels, not a new quantum service. Physical
pulse, storage and transfer costs remain supplied apparatus parameters.

The result identifies a finite event-count comparison. It supplies no
observer-clock unit in seconds, no region-to-region physical rate, and no
CP02 integrated-gate acceptance.

Run the source-bound replay and independent tests with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.clock_defect_energy.clock_count_inference --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/clock_defect_energy/test_clock_count_inference.py
```
