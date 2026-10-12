# A common interacting observer process

**One finite source program couples quantum electric flux, charged
fermions and a neutral scalar, measures that shared state, and uses its
records for feedback and native repair.** An exact certificate proves a
nonzero scalar response to the recorded matter feedback. A separate exact
certificate proves that this response survives infinite oscillator
occupation and the actual coordinate-sign readout. This is a construction
toward OPH's [physical completion target](../CONSTRUCTION.md#completion-target-and-acceptance-criteria).
That target remains open.

## Run and inspect

From the repository root, after installing `model/requirements.txt`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.common_source --check model/common_source/results.json
python3 -m model.common_source.infinite_history --check model/common_source/infinite_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/common_source
```

`--write PATH` retains a new receipt and refuses to overwrite an existing
file. [results.json](results.json) contains the complete quantum branches,
all 1,056 native/quantum histories, controls, source-service counts and
rational error certificates. The [derivation](DERIVATION.md) supplies the
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
pulses execute through the inherited twelve-port source compiler. The
native integer-repair channel remains the existing declared primitive in
the joint instrument. Its entire 78-dimensional channel is not compiled
into the six-dimensional processor here. Coherent code transfer, blank
preparation, storage and readout remain supplied services.

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
the separate readout treatment above. Neither proof supplies a spatial
continuum or a physical service error budget.

The infinite-history derivation also proves that the finite sign instruments
converge to this same infinite instrument, including terminal states. With
`d=2^k` and `d^2` exact palindromic steps per interval, the mathematical
finite-family contrast is eventually strictly negative. This is a limit
theorem, without an explicit sufficient d or a larger native compiler run.

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
clock are supplied. In particular that preparation is not established as
the ground state of the coupled Hamiltonian. P and N retain their existing
closure and physical-readback contracts; they do not fix these additional
choices. This construction has no spatial refinement, established
interacting 3+1-dimensional effective field theory, responsive gravity,
physical clock calibration or independent empirical validation. Those
must refer to observations of this same construction before they can
complete OPH. Fundamental fields or spacetime are not required; their
observational content is.
