# A common interacting observer process

**One finite source program couples quantum electric flux, charged
fermions and a neutral scalar, measures that shared state, and uses its
records for feedback and native repair.** An exact certificate proves a
nonzero scalar response to the recorded matter feedback. A separate exact
certificate proves that this response survives infinite oscillator
occupation and the actual coordinate-sign readout. This is a construction
toward OPH's [physical completion target](../CONSTRUCTION.md#completion-target-and-acceptance-criteria).
That target remains open.

The canonical finite [executable specification](COMPOSITION.md) and
[program runtime](protocol.py) now compose these operations on retained
states. Both original sharp and smooth histories use this engine. The new
[four-read continuation](composition.py) adds a second smooth field read,
record-controlled feedback and a later gauge clock probe to the same
preparation. Its complete native attachment retains all 4,224 joint outcomes
at d=2 and d=4. The smooth occupation theorem now covers those later reads
and evolution, including their quantum output rather than only probabilities.

```sh
python -m model.common_source.composition --check model/common_source/composition_results.json
```

The [smooth-history certificate](SMOOTH_HISTORY.md) compares all 16 ideal
finite work outputs at d=2,4 with the infinite history for each of the three
native currents, including later feedback and the gauge clock probe. Its
[producer](smooth_history.py) and [receipt](smooth_history_results.json)
fix the blank-gauge and scalar-vacuum preparation. The broader uniform
moment bound and the missing numerical and physical service bounds remain
separate.

The [coherent native attachment](QUANTUM_NATIVE.md) executes the native
repairs around that same continuation with the existing source controls.
It accepts the full 78-state native carrier with arbitrary reference
coherence, retains every raw outcome, and uses an executed port permutation
for the final record-selected seam. C and E remain distinct quantum record
laws. Its [runtime](quantum_native.py) and
[receipt](quantum_native_results.json) keep initialization, code services,
storage and physical calibration as explicit inputs. The earlier classical
native benchmark remains available with its original scope and counts.

The [encoded memory variant](ENCODED_MEMORY.md) uses
`native_mode='protected_qubit'` to recover the existing four-configuration
logical qubit after each native event. The same field continues through its
four reads; the `passive_qubit` control uses the same encoded input and seam
policy without recovery. This is active correction under the existing
coherent services, with a separate contract from full-carrier recovery or
physical memory calibration. The [producer](encoded_observer.py) and
[receipt](encoded_observer_results.json) compare the active and passive
complete instruments, including their source costs.

The [current source synthesis](QUANTUM_NATIVE.md#exact-gray-parity-synthesis)
uses the established Walsh/Gray parity method to reduce each recurring
nine-bit controlled flip from 11,782 to 2,044 native pulses. It preserves
the complete 512-dimensional endpoint operator, including padding and
reference coherence, with no extra qubits. The encoder now uses 28,635
pulses; the `d=2` protected C history uses at most 1,158,423. These counts
include every repeated source request. Physical services remain supplied,
and the shorter words need not preserve intermediate source states.

The [observer readout bridge](OBSERVER_READOUTS.md) reconstructs Hilbert
representations and modular flow from this same executed state, extends
its complete records with a source-compiled binary modular clock, and
extracts the interactions between its declared regions. The
[occupation analysis](OCCUPATION_REFINEMENT.md) adds an interacting
two-oscillator limit. The [infinite-history certificate](INFINITE_HISTORY.md)
encloses the limiting sharp-readout response itself. The
[complete-output comparison](READOUT.md#complete-comparisons-with-the-actual-finite-pointer)
also bounds the finite detector's retained quantum state, for the same
zero-current preparation with feedback enabled and disabled. At d=2 its
distance from the infinite sharp history is about 0.455 with feedback;
even at d=256 every output confined to that occupation code has error
above 0.122. Accurate sign probabilities alone do not certify an accurate
postmeasurement state.

A separate [finite-resolution detector](SMOOTH_REFINEMENT.md) replaces that
sharp sign measurement with a source-executed Ramsey readout. Its complete
interacting record/quantum output has a quantitative occupation error bound,
including measurement backreaction. The infinite-model unconditional
measurement energy increase is exactly `1/8` in the declared model units.
The executed finite cutoffs retain their actual energy changes and boundary
corrections. This is a different detector; the original sharp-record clock
is not transferred to it.

The [history repair construction](HISTORY_REPAIR.md) solves the original
sharp-readout event history as local compatibility constraints. Every fair update schedule
reaches its unique completed assignment for the fixed apparatus inputs.
Its observer states and modular readouts are consequently independent of
the auxiliary repair counter. The supplied record meanings still select
the actual feedback operations inside that history.

The [native slow-geometry readout](SLOW_GEOMETRY.md) attaches the existing
rank-three response carrier to all 1,056 joint histories. The original
two-seam controller spans two displacement directions. An additional
binary continuation uses three independent seam directions across its two
events, with every new native outcome and its quantum reference retained.
These are dimensionless carrier displacements; physical position, momentum
and spacetime remain separate attachments.

The [calibration checks](CALIBRATION.md) compare carrier distances, sampled
modular-clock readings and finite detector energy with their declared
mathematical targets. They also show why choosing the physical speed of
light would relate length and time units without determining an energy
scale. A propagation law connecting these readouts remains to be built.

The [fixed-reference clock bridge](THERMAL_CLOCK.md) distinguishes a
stationary equilibrium reference from the actual changing observation state.
Choosing a faithful Gibbs reference for the existing full Hamiltonian makes
its continuous between-event dynamics modular, with `t=-beta s`.
Source-executed clock probes retain the earlier records. The reference is
constructed from the supplied action; its selection by repair and a physical
time calibration remain unproved. The earlier local-clock construction and
the failed identification using the local post-measurement state are retained.
The earlier clock probes appended to the sharp history are outside its
infinite-readout certificate. The [smooth continuation](COMPOSITION.md)
separately controls the actual gauge-probe instrument after two smooth reads;
it supplies no Gibbs-reference limit or physical clock calibration.

The [native propagation check](CAUSALITY.md) verifies exact dependency cones
for fixed seam schedules on all 78 admitted load states. Concrete paired
histories show information arriving only after the intervening local updates.
Choosing one edge and one causal layer as units gives a graph speed bound of
one. Relating that causal depth to modular readings is an operational
consistency requirement, not a demand for an external fundamental ruler.

The [native field calculation](NATIVE_FIELDS.md) derives local continuity
and the conditional mean current from the existing integer repair. Its
full 78-state kernel checks graph diffusion for the existing uniform-seam
calibration. This mean relaxes while individual recorded histories continue;
its real decay modes cannot be identified with undamped field oscillations.
The result adds no field Hamiltonian or clock convention.

## Run and inspect

From the repository root, after installing `model/requirements.txt`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.common_source --check model/common_source/results.json
python3 -m model.common_source.infinite_history --check model/common_source/infinite_results.json
python3 -m model.common_source.sharp_output --check model/common_source/sharp_output_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.common_source.bridge --check model/common_source/bridge_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.common_source.smooth_bridge --check model/common_source/smooth_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.common_source.thermal_bridge --check model/common_source/thermal_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.common_source.causality --check model/common_source/causality_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.common_source.native_fields --check model/common_source/native_fields_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/common_source
```

`--write PATH` retains a new receipt and refuses to overwrite an existing
file. [results.json](results.json) contains the complete quantum branches,
all 1,056 native/quantum histories, controls, source-service counts and
rational error certificates. [bridge_results.json](bridge_results.json)
retains the Hilbert/modular/regional readouts and every appended clock branch.
[smooth_results.json](smooth_results.json) retains the separate finite-resolution
histories, feedback and disconnected controls, energy ledger and error bounds.
[infinite_results.json](infinite_results.json) retains the exact infinite
sharp-readout signal, state-error budgets, polynomial hashes and controls.
[sharp_output_results.json](sharp_output_results.json) retains the complete
finite/infinite output comparisons, unavoidable code-support errors and
independently checked finite-projector residuals.
[thermal_results.json](thermal_results.json) retains the separate fixed-reference
clock continuations; [causality_results.json](causality_results.json) retains
the fixed-schedule native dependency and delayed-influence checks.
The [derivation](DERIVATION.md) supplies the
carrier maps, proof and theorem boundaries. `python3 -m model.verify`
checks this package's hashes and replays it with the other model packets.

## The same process from preparation to feedback

The quantum working register contains two opposite-charge fermion
channels, a quantum electric link tied to their occupation by operator
Gauss constraints, and two neutral scalar sites. Scalar occupation
dimensions 2 and 4 are executed. The flux-dressed four-state code represents
the declared gauge-invariant matter/link operations on an exact invariant
sector. Its electric energy
does not commute with charged hopping.

The dimensionless finite action is

```text
H = H_hop + E^2 + H_scalar + g sqrt(2) Q0 (n_L0-n_R0),
g = 1/4,     kappa_scalar = 9/32.
```

Each scalar onsite quadratic is compressed from the infinite oscillator
before taking the cutoff; squaring the finite quadratures would change
the boundary term. The scalar/matter interaction acts on the same charged
register as the electric and hopping terms. Its state becomes entangled
across sectors. g, the action, preparation and schedule are declared
construction choices, without an identification of g with alpha or P.

The complete experiment proceeds as follows:

1. Native repair records a current. It controls an emitted gauge/matter
   preparation pulse on the common quantum register.
2. That register evolves for model-coordinate interval 1. A native source
   pointer measures the right occupation of the first charged species and
   commits bit a.
3. Only a selects a further charged-hopping pulse. The register evolves for
   interval 1/2 under the same interacting action.
4. A second source pointer measures the sign of the first scalar coordinate
   and commits bit b. Only b selects the next native repair seam.

All quantum preparations, evolution factors, pointer gates and controlled
pulses execute through the inherited twelve-port source compiler. In this
earlier sharp-history benchmark the native integer-repair channel remains
the declared primitive. The separate [coherent attachment](QUANTUM_NATIVE.md)
executes the full native interface around the four-read smooth continuation.
Coherent code transfer, blank preparation, storage and readout remain
supplied services in both constructions.

At 512 steps per interval, the quantum protocol uses the following emitted
resources per realized history, including preparation, both pointers,
feedback and zero-angle requests:

| Scalar dimension | Working qubits | Native pulses | Coherent load/export events |
| --- | --- | --- | --- |
| 2 | 4 | 176,182 | 352,364 |
| 4 | 6 | 1,879,628 | 3,759,256 |

Two fresh pointer blanks and two committed quantum bits are additional.
The native source and native records are accounted separately. These are
expanded word counts; the numerical replay reuses executed step matrices.
Physical service duration and initialization/storage/readout costs remain
unspecified. No efficiency or physical cost bound follows from host runtime.

Every complete path retains both quantum bits, both native records, final
native loads and its unnormalized terminal quantum state. Probabilities
are their squared norms, with no successful-branch selection. Committed
records keep their values while the working state changes. Complete
quantum histories, rather than independently supplied scalar/matter
probabilities, feed the native continuation. The common record law gives
compatible marginal statistics for these observers; it does not prove
general inter-observer consensus dynamics.

## A rigorously nonzero finite response

For scalar dimension 2, zero native preparation current and the declared
intervals, exact rational Taylor arithmetic encloses the unconditional
probability of a positive scalar coordinate:

| Finite continuous-Hamiltonian protocol | Exact outward interval |
| --- | --- |
| Matter-record feedback enabled | [0.370697726608684, 0.370697726608685] |
| Feedback disabled | [0.402359206973491, 0.402359206973492] |
| Enabled minus disabled | [−0.031661480364808, −0.031661480364807] |

Including both protocols' ideal temporal splitting errors at 512 steps per
interval gives **[−0.033297300667839, −0.030025660061775]** for the
ideal compiled-word contrast. It excludes zero. The bound is derived from
the exact action and complete instrument; numerical predictions are not
inputs to the certificate. With g set to zero, gauge-only feedback cannot
change the disconnected scalar marginal. This control is executed with
the same requests and counted services, including zero-angle pulses.

A separate exact electric echo changes a charged occupation probability
from zero to `576/3125`. It checks electric backreaction on the same matter
carrier. The retained joined experiment also compares actual native
feedback with a fixed seam and an independent actuator bit having the same
unconditional frequency.

The temporal certificate compares exact finite-Hamiltonian evolution with
ideal exact-factor source words, for the entire recorded instrument and an
arbitrary spectator reference. It excludes binary64 coefficient and
execution errors and physical service errors. Actual floating source
execution is separately checked against independent matrices. The two
executed cutoffs provide a numerical diagnostic. The following independent
certificate supplies the interacting occupation and infinite-readout result.

## Removing the occupation and readout approximations

For the same interacting action, infinite coupled free vacuum, zero native
preparation current and feedback schedule, the actual coordinate half-line
readout has unconditional enabled-minus-disabled probability in
**[−0.024045670143, −0.024043450955]**. This interval excludes zero after
including initial-state projection, full-Hamiltonian Taylor remainders,
discarded occupation support, coefficient rounding and the final readout.
Removing the interaction gives exactly zero contrast.

This is an exact rational certificate for two infinite oscillators, with
all records and quantum outcome states retained. The
[derivation](INFINITE_HISTORY.md) and compact [receipt](infinite_results.json)
are reproducible without floating simulation as an input. The terminal
polynomials can be regenerated with `--polynomials --write PATH.json`.

The [readout audit](READOUT.md) also identifies a necessary distinction:
the sign projector of a truncated coordinate is not the compression of the
infinite sign projector. Their operator difference already has norm about
0.101058 at d=2. The previous finite signal remains valid for its detector;
the new certificate uses the actual infinite detector, including its
postmeasurement tails.

A separate [occupation theorem](OCCUPATION.md) bounds the rectangular
finite-Hamiltonian family through the first retained matter measurement
and feedback, uniformly over arbitrary initial matter states and untouched
references with the fixed scalar vacuum. At onsite d=64 the pre-readout
trace-distance error is below 0.000002498. The final scalar read requires
the separate readout treatment above. The independent
[fourth-moment analysis](OCCUPATION_REFINEMENT.md) supplies a further
pre-readout bound and proves strong convergence of the complete sharp
history for all three native currents, including its native continuation.
Neither pre-readout estimate alone bounds the final sharp detector. The
separate complete-output calculation now adds the missing projector
comparison for the zero-current preparation, with every quantum branch
retained. It gives finite error enclosures at d=32,64,128 and a direct d=2
comparison including the ideal source splitting budget. These mathematical
cutoff comparisons supply no additional source execution, spatial continuum
or physical service error budget.

For zero native preparation current, the certified complete continuous-output
distances are:

| Onsite occupation dimension | Feedback enabled | Feedback disabled |
| --- | --- | --- |
| 2 | [0.455154, 0.455155] | [0.457642, 0.457643] |
| 32 | [0.203713, 0.226530] | [0.203935, 0.226752] |
| 64 | [0.180917, 0.180924] | [0.181090, 0.181097] |
| 128 | [0.152201, 0.152204] | [0.152341, 0.152344] |

These intervals are rounded outwards from the exact receipt. At d=2 the
additional ideal 512-step source budget gives [0.454336, 0.455973] with
feedback and [0.456824, 0.458461] without it. At d=64 the pre-readout error
is below 0.000002498, so it cannot explain the much larger final error:
the sharp detector dominates. For the vacuum control, optimal approximation
within these occupation codes requires dimension `Omega(epsilon^-4)` for
complete-output error epsilon. That resource bound concerns the specified
Fock code and that control, not every abstract finite encoding or the
interacting histories. The separate smooth detector avoids this discontinuous readout
while changing the observable being measured.

The infinite-history derivation also proves that the finite sign instruments
converge to this same infinite instrument, including terminal states. With
`d=2^k` and `d^2` exact palindromic steps per interval, the mathematical
finite-family contrast is eventually strictly negative. This is a limit
theorem, without an explicit sufficient d or a larger native compiler run.
It concerns zero native preparation current. The sharp postread quantum
states retain their full half-line wavefunctions; their kinetic energy
need not be finite. The pre-readout moment bounds therefore cannot be
reused automatically through further field evolution and sharp reads.

For the separate smooth detector, the full ideal history error decreases as
`O(1/d)` at fixed model time. Its conservative bound is about `0.2063` at
occupation dimension 64, `0.04096` at 256 and `0.008840` at 1024. These are
analytic bounds, not executions of those larger source programs. At the
executed dimensions 2 and 4 the bound is vacuous, so agreement between those
small cutoffs does not certify continuum precision. The joint ideal limit
uses `n(d)=d^2` temporal steps; floating execution and physical service errors
remain separate.

The earlier local-reference and post-sharp fixed-Gibbs-reference probes have
no proved occupation limit in those detector analyses. The separate
[continued smooth program](COMPOSITION.md) bounds its later gauge probe on
the retained state, without asserting convergence of Gibbs modular data.

## What this establishes and what remains

This packet binds source control, CAR, operator Gauss, record readback and
feedback to one explicit interacting finite process. It reuses the
existing abstract common source and theorems, and adds executed carrier
maps and a certified operational response. [Theorem attachments](DERIVATION.md#compatible-theorem-use-and-open-physical-connections)
distinguish committed-record persistence from writing, block readback
from commutative record readback, and exact hardware refinement from
approximate occupation cutoffs. No full Lean observer or consensus-tower
instantiation is asserted.

The action, initial projected free vacuum, couplings and model-coordinate
schedule are supplied. The local-reference clock reads an ensemble
probability on a bounded phase window. The fixed-Gibbs-reference bridge
relates modular and Hamiltonian parameters on continuous evolution segments;
it does not turn the intervening measurements and feedback into one modular
unitary flow. Both protocols use declared references and pulse schedules,
without physical duration calibration. The actual preparation remains the
projected free vacuum, not the Gibbs reference or an established ground state
of the coupled Hamiltonian. P and N retain their existing
closure and physical-readback contracts; they do not fix these additional
choices. This construction has no spatial refinement, established
interacting 3+1-dimensional effective field theory, responsive gravity,
physical clock calibration or independent empirical validation. Those
must refer to observations of this same construction before they can
complete OPH. Fundamental fields or spacetime are not required; their
observational content is.
