# Source-pulse Ramsey control of the native density wave

This packet compiles and executes the missing Ramsey/controller service for
the existing finite density walk, then uses it in the complete joint clock
and field instrument. The [versioned boundary](ramsey_campaign.json) keeps
C and E, the transported parameter-four clock menu, protected field storage,
blank workspace and terminal computational-basis read explicit. It does not
identify a physical duration or accept the integrated CP02 gate.

The starting [density compiler](../observer/COHERENT_CONTINUATION.md#source-compilation-of-the-finite-density-walk)
has 720 valid port/event arcs, 304 address-padding values and five clean
scratch bits. Its 15-qubit walk uses 114,910 source pulses. Its exact
mean-zero density-sector intertwiner identifies the cyclic density
amplitudes with those of the 3960-arc binary source. This packet compiles
that density representation; it makes no claim to the full 3960-arc unitary
or its higher-configuration observables.

## Controlled source program

One additional qubit stores the Ramsey arm. The reference arm is the
identity. The live arm executes the existing reflection and reversal.
The source compiler emits the arm-controlled event reflection as a phase
word between the inherited event-preparation word and its inverse. Global
phases are kept: replacing the live walk by its negative preserves a
standalone unitary channel and changes the Ramsey measurement.

The reversal uses the original reversible event decoder. Adding the Ramsey
condition to each endpoint-target operation uses an available scratch bit,
with computation and uncomputation on both arms. No extra scratch is added.
Exact integer execution checks all 65,536 basis inputs for bijectivity, the
entire reference arm for identity, and all 2048 clean-scratch address/arm
inputs for the intended permutation and clean return. Address padding has
its declared unitary extension. Dirty-scratch live-arm behavior is retained
but is not identified with the predecessor's extension.

The resulting controlled walk uses **123,646 native pulses**, each with a
coherent proper-code load and export. Source words also prepare the Ramsey
arm, analyze it, prepare the uniform event state and prepare the fixed
normalized port probe `(11,-1,...,-1)/sqrt(132)`. The fixed probe corresponds
under the exact density map to the original centered binary load probe.
Its Gray-elimination angles use the integer amplitude ratios at eighty
decimal digits before the declared dyadic rounding, as in the inherited
uniform-event compiler.

The field and encoded clock use 25 logical register qubits in total, with
the shared six-level source processor counted separately. The receipt
measures compact program bytes, expanded-word bytes and the 17-bit walk
position counter. This interpreter storage is separate from the 132-state
core/repair policy, the FIFO and its exported records. The accounting does
not derive physical memory or bound the host interpreter's workspace.

Every source routine is executed pulse by pulse on all its input columns.
The 16-qubit calculation reuses those original arrays as tensor factors.
It retains scratch, address padding and arbitrary spectator columns. It
does not replay the expanded 123,646-pulse word literally on every full
register column, insert ideal target arrays, normalize the output, or apply
a polar correction. A separate small-register test compares literal source
pulses with factor reuse through overlapping controls.

## Complete joint law and continuation

The existing two-attempt C/E policy chooses one or two incident ticks.
For each terminal source frame `J_h`, this packet executes the complete
field service and retains both Ramsey outputs. All clock export aborts
retain the unmeasured field at the number of successful ticks already
accumulated, along with the untouched reference. The predecessor's four
native preparation flags remain present; failed preparations bypass the
field service. Computational-basis arm reading is a declared terminal
instrument; coherent preparation and analysis are the new source words.

The source arrays reproduce the original density-sector plus probabilities
`21/22` and `91/110` for one and two ticks within the declared numerical
allowance. Tests also use complex source/field/reference entanglement,
off-diagonal reference entries, explicit clock-export leakage and invalid
processor input with nonzero proper-code failure. Complete branch maps
compose by finite Kraus amplification; no branch is normalized away.

A separate four-attempt retained history moves the reference core through
`{1,7}`, `{1,6}`, `{1,2}` and `{1,3}`. Each decoded native record selects an
executed word from the all-66 clock inventory and advances the same live
wave register. Three records leave the one-slot FIFO. Neither reference
nor field is prepared again. The four clock-abort positions preserve field
prefixes with zero, one, two and three successful ticks. This selected
history is unnormalized; its surrounding complete instrument is supplied
by the inherited full C/E branch kernel and every source-export outcome.
It is not reported as an exhaustive four-step probability distribution.

## Bounds and remaining premises

For an operation routine, the checked full-array Frobenius error bounds its
operator error on arbitrary input and reference. State-preparation checks
use the original blank-input column, with arbitrary spectator references.
If the sum of the counted routine allowances is `s<1`, the complete product
operator allowance is `delta=s/(1-s)` and its half-trace allowance is
`delta*(2+delta)/2`. Classical branch selection uses the maximum service
allowance, without multiplying by the number of records. The complete
clock/field allowance composes this with the predecessor's source-clock
gain and preparation bound.
The replay executes the source clock arrays and consumes their unchanged
complete-error allowance from the hash-verified joint-clock receipt. It
does not repeat the predecessor's exact native/probe certificate solve.
Clock words selected along the four-attempt exported history are freshly
executed and certified by the all-core service interface.

These routine guards are floating execution checks, not interval proofs of
the numerical backend or physical precision certificates. Physical pulse
errors and coherent transfer errors retain separate allowances; `P` pulses
with operator error at most `epsilon` add at most
`P*epsilon/(1-P*epsilon)` when `P*epsilon<1`. The proper-code source
generators have exact invariant codes, so their complement outcomes have
zero probability on the admitted inputs. The complete exporter retains
those outcomes and is tested on invalid processor inputs where they are
nonzero. Processor failure and application-address padding are distinct.

The finite theorem attachment uses the source density-sector Gram identity
and walk intertwiner, complete Kraus composition from
[ChoiCPTP](../../Lean/Dynamics/ChoiCPTP.lean), and the finite parameter/readout
interface in [ClockReadout](../../Lean/Time/ClockReadout.lean). No new Lean
formalization is claimed. The coherent continuation law remains supplied:
the existing identity/reflection and record-phase countercontrols are
preserved. Physical clock calibration, regional reference compatibility,
coherent storage and the physical implementation of transfer/read services
remain additional tasks.

Replay with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.clock_defect_energy.ramsey_service --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/clock_defect_energy/test_ramsey_service.py
```
