# A complete joint event-clock and field protocol

This packet executes a finite join between the transported modular clock and
the native field's coherent walk. One committed native incident count advances
both, and their terminal outcomes are retained as one distribution. The clock
and field generators need not coincide. The join uses their existing coherent
service contracts; it does not derive a single native clock/field primitive or
a calibration in seconds.

The source is the pending PR1149 head recorded by
[the successor contract](joins_campaign.json). The earlier
[clock/wave diagnostic](CLOCK_WAVE.md) is preserved. This construction advances
the CP02 interface beyond propagating a clock interval through a separate field
experiment: the clock phase and wave are advanced within the same history.

## Protocol and complete records

The clock source is binary `H66`, with the existing local prepared reference
whose three ideal weights are `2/225`, `1/36`, and `2/45`. Those weights define
the compiled target; the unknown finite prepared state is compared using its
existing bounds. The field is a **separate protected register** of 3960 native
arc coordinates. Its supplied preparation is `A(6*n_0-1)/sqrt(330)`, with a
Ramsey reference arm and a live arm. Passive source repairs do not act on this
register, and source time is paused during coherent services.

1. Run the existing bounded native clock preparation. Keep success, herald
   failure, warmup timeout, and completion timeout. Only success continues.
   Its earlier fine tape has exactly the predecessor's coarse-grained interface
   and is unavailable to the later controller; success is never assumed certain.
2. Apply the existing two-address probe rotation. Attempt seam `(0,5)`.
   Record the native coin and the full C current or E contrast.
3. Update the reference core from that committed record. This port-0 incident,
   including zero current, requests the existing transported parameter-4 word.
   On successful export advance the live wave arm by one `W` step.
4. If the first current was positive, attempt `(0,7)`; otherwise attempt `(1,8)`.
   Retain the second native record. Only the former seam is a port-0 incident,
   so it requests a second transported parameter-4 word and a second `W` step.
5. Apply the predecessor's conditioned probe analyzer, read the complete native
   `(1,5)` detector instrument, and read both field Ramsey outcomes.

Thus the accepted complete histories have `N=1` or `N=2` incident ticks,
clock parameter `s=4N`, and wave step count `N`. Every source export abort is a
separate terminal output, with the already accumulated wave state and the
untouched external reference retained. No abort is normalized away. No source
or field is reset between the two attempts. Unattained native raw-register
slots are zero maps on all inputs, using the predecessor's explicit raw
decoder; zero probability is distinct from dropping a supported record.

Only already certified parameter-4 words are needed here. The broader all-core
inventory in the same campaign is independent progress, not a premise smuggled
into this two-attempt join. The finite alphabet and event count are known to
the controller. The terminal phase detector is statistically different for
one and two words; it is not a perfect one-shot tick decoder or an unrestricted
phase-unwrapping scheme.

## Joint law and matched-marginal control

The exact ideal-target clock probabilities have coefficients in
`{1, cos(4*log(8/5)), cos(8*log(8/5))}`. They are derived from the actual C/E
partial-injection branches and the two-address rotations. Rational phase
enclosures come from the existing source-word certificate, independently of
floating execution. The complete executed source arrays are checked against
that law; they include all 512 source addresses before export.

For either C or E, the ideal target gives
`P(N=1)=307/360` and `P(N=2)=53/360`. The existing exact walk gives
`P(field+|N=1)=21/22` and `P(field+|N=2)=91/110`.
The receipt retains the joint law with native detector current, both field
outcomes, and a digest of every attained complete clock-record coefficient.

The countercontrol independently resamples the field's Ramsey outcome from
its unconditional marginal. It preserves **the entire clock-record marginal**
and the field-outcome marginal exactly, while changing their joint law. Its
total variation distance is computed exactly. The negative native detector
current and field-plus outcome also have a nonzero covariance. A separate
control applies only the first clock word while retaining both native attempts
and the same wave triggers; its changed detector law verifies that the clock
phase itself responds to the shared incident count.

Including ideal preparation failures multiplies the first control distance by
the herald mass `5/33`; the failure outputs are unchanged. The actual bounded
preparation is executed as well and its complete flagged joint distribution
is checked against the retained preparation-plus-source allowance. The
conditional displayed law therefore does not hide the cost of preparation.

These are finite apparatus diagnostics, not new empirical comparisons.
Classical record-mediated correlation is the admitted coupling; no direct
clock-field Hamiltonian interaction or spacetime metric is inferred.

## Full quantum and error statement

Let `J_h` denote an executed clock-source operator for a complete successful
native history and detector read. Its joint terminal Kraus operators are

`J_h tensor (I +/- W^N(h))/2 tensor I_R`.

An abort uses its actual source export operator, the accumulated controlled
walk, and identity on `R`; padding is traced, retaining the reference density.
Since `W` is unitary, the sum of the two Ramsey effects is identity. Combining
this with the complete native/source instrument proves completeness for
arbitrary clock-field entanglement and arbitrary untouched references, not
only the product input used for the displayed law. The tests execute an
entangled source/field/reference input, retain an off-diagonal reference, and
inject a nonzero export failure to test the abort convention.

The word certificates are operator-norm bounds on the actually executed source
arrays. The existing finite product/channel theorem supplies their complete
half-diamond error, including references and branch-conditioned words. At most
two clock words and two analyzer words occur. Native frame errors and complete
preparation errors are charged as well. The separate parameter-8 preparation
stability certificate concerns the initial reference's modular channel; it
does not assert exact covariance of an unknown prepared reference through all
record projections.

The exact arc walk is an admitted coherent operation. This packet does not
compile the whole 3960-dimensional controlled walk from source pulses. If
field preparation, each controlled walk, and the Ramsey read have CPTP
half-diamond allowances `epsilon_p`, `epsilon_W`, and `epsilon_r`, their
additional complete allowance is bounded by
`source_gain*(epsilon_p+2*epsilon_W+epsilon_r)`. These errors and any physical
hardware error remain explicit service premises. The receipt gives both
source-only and complete preparation-plus-source bounds, along with pulses,
loads/exports, source and wave registers, record bits, and finite attempt caps.
The host interpreter's memory is not claimed to have been physically bounded.

## Finite theorem reuse

| Existing result | Predicates checked here | Conclusion used and boundary |
| --- | --- | --- |
| [Native history transport](../native_modular_reference/history_transport.py) | Actual complete C/E partial-injection records; faithful diagonal reference; core update from current; no discarded supported outcome | `U_h(s) K_h = K_h U_0(s)` on the admitted branch, so interleaved incident words equal the count-dependent ideal action. No unknown-reference covariance or phase gate inferred from records alone. |
| [Transported source clock](../native_modular_reference/history_clock.py) | Same 512-address encoding, already certified parameter-4 diagonal words, actual analyzer pairs and native detector | Executed source operations and complete export/reference error bounds. Coherent load/export and phase services remain supplied. |
| [Coherent native continuation](../observer/coherent_continuation.py) | The actual 3960-arc reversal, coherent event embedding, phase-minus reflection, and unchanged register between steps | Unitary `W`, native density recurrence, and exact two-step Ramsey probabilities. The walk is an added coherent continuation, not passive native repair. |
| [Existing active wave clock](../relative_observer_geometry/temporal_records.py) | Separate protected field register; zero-inclusive incidence; timer paused during services | Reuses the operational event-trigger architecture. The new addition is its **joint transported modular detector**, with C and E retained. |
| [Finite modular clock comparison](../common_source/THERMAL_CLOCK.md) | Clock and field act on distinct registers; record-conditioned services declared | No requirement that clock and field generators agree. Full-group equality is a stronger, optional Gibbs/reference statement. |
| [Clock readout](../../Lean/Time/ClockReadout.lean) and [clock comparison](../../Lean/Time/ClockComparison.lean) | Parameter and event dictionary declared; no unmeasured physical calibration | The finite instrument gives an explicit operational relation. It does not use affine regrading theorems to create missing anchors or proper time. |

Run the source-bound replay with
`python -m model.clock_defect_energy.joint_clock --check` and the targeted
tests with `pytest model/clock_defect_energy/test_joint_clock.py`.
