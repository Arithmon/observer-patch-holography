# The runnable OPH model

This directory contains the executable observer construction: native repair,
quantum instruments, recorded feedback, observer-accessible algebras and
checkpoint continuation. Its finite predictions are computed from declared
source operations and compared with exact mathematical targets or public
measurement data under explicit apparatus assumptions.

The construction supports interacting quantum dynamics, local event clocks,
effective field descriptions with memory, and quantitative observation bounds
across federations of elementary observers. **A physical identification with
our Universe, including interacting relativistic spacetime, remains open.**
The [recovery assessment](OBSERVER_RELATIVE_RECOVERY.md) states the current
results and remaining obligations. The
[construction specification](CONSTRUCTION.md#completion-target-and-acceptance-criteria)
defines the acceptance criteria.

## Architecture and execution

The [primitive observer runtime](primitive_observer/README.md) executes
complete classical/quantum instruments. Each operation acts on a continuing
quantum state, retains its declared classical outcomes, and may select an
intervention using those records. Hidden alternatives are summed as CP maps.
Failed preparations, aborts and finite-cap timeouts remain in the complete
observation law. Serialized checkpoints retain the state and records and are
validated by replay of their declared prefix.

The native construction uses twelve-port elementary carriers and their
seam wiring. Federations join those same elementary factors; they do not
introduce a different kind of larger elementary observer. C and E are two
specified quantum continuations of the native repair: their configuration
laws agree, while the records that distinguish their Kraus alternatives
produce different coherence laws. The
[candidate register](CANDIDATES.md) describes their domains and interfaces.

A declared instance supplies its carrier, accessible algebra, instruments,
repair law, preparation domain and scheduling protocol. Consequences of
those inputs can be derived without proving that the instance is uniquely
selected by the OPH axioms. The [source-selection assessment](SOURCE_SELECTION.md)
tracks that distinction. In particular, committed native repairs alone
cannot create coherence from a diagonal input. The
[compiled native interface](compiled_native_interface/README.md) implements
phase preparation and detection through coherent proper-code services:
a balanced mixer uses 34 source pulses, and its imaginary quadrature adds
one pulse. Source-array comparison bounds retain all encoding and padding
coordinates; they do not certify physical gate fidelity.

Every rate or time statement names its parameter. Native attempt counts,
local recorded events, modular parameters, calibrated observer clocks and
auxiliary solver counters are distinct. A preferred universal physical
clock is not required. Their operational correspondences are tested through
joint observations rather than asserted from a shared name. The source record
laws are invariant under any positive rescaling of attempt count as duration,
so seconds and hertz require paired event-time records and independent
calibration; no such native-to-lab clock attachment is established here.

## Native action, local clocks and effective fields

The [local native action](native_action/README.md) is derived directly from
the endpoint instruments. On operator space,
`H_native = sum_e (I-Phi_e)` is positive and local, with exclusion exchange
and C/E coherence-killing coefficients fixed by the repair law. It is a
Liouville operator, not a Hamiltonian on the carrier's state vectors. The
exact uniform discrete transfer is `I-H_native/M`. Independent rate-one
Poisson scheduling per seam gives `exp(-s H_native)` as a separate declared
scheduler construction; its parameter is not identified with physical or
modular time.

The integer-load rule also has an exact local cost-minimization law: at fixed
selected-pair sum, each repair minimizes the pair's quadratic load cost. The
global potential `V=||x||^2/2` decreases by `(d^2-r)/4`, where
`r=|d| mod 2`. Under the separately declared
uniform all-seam schedule, its mean is discrete gradient descent of graph
disagreement. This is a source-derived dissipative optimization law, not a
spacetime path action. The executed quantum channel strictly attenuates a
coherent mode, ruling out unitary conjugation on this carrier; the principal
logarithm's half-step map is not positive, so it does not generate a
completely positive semigroup. Other dilations or generators are not thereby
excluded. The same finite kernel yields an exact additive conditional path
surprisal in nats, with the initial-state prior kept separate. It is not the
quadratic repair cost or mean-field generator and supplies no physical action
units.

The [native quantum transfer](native_quantum_history/README.md) and
[effective fields](native_effective_fields/README.md) provide finite
operator reconstruction, non-Gaussian current correlations and controlled
logarithm comparisons. The [regional field construction](relative_observer_geometry/README.md)
recovers local occupations from disclosed preparations and incident-current
records, with compatible overlaps and explicit clock-comparison tests. Those
records also reconstruct an exact, unitless rank-three metric on the twelve
ports: all thirty native exchange seams are exactly its shortest shell, and
the native field generator is the sum of exchange terms on those edges. This
is a finite spatial interpretation of the declared carrier, not a calibrated
or extended space. On levels 0→1 and 1→2, the existing port-ancestry map fails
to preserve seam adjacency: respectively, 74 of 2,520 and 280 of 10,080 fine
seams map to non-seams, while 2 and 4 fine seams collapse to one coarse port.
The cell-ancestry map is not record-derived, and cross-carrier routes use a
supplied cell-centre/local-tangent-frame convention (`source_derived: false`);
the local carrier Gram itself is record-derived. Exact native-hop balls give
finite dimensionless distance data, with roughly quadratic growth over part of
the tested range, but this does not establish a three-dimensional physical
bulk or path-independent global frames; continuum spacetime remains open.

The [native causal-cone check](common_source/CAUSALITY.md) adds one exact
space/time relation at this finite level: for a fixed seam schedule and shared
local coins, dependence can cross at most one seam per ordered repair event,
and retained signal witnesses attain that bound. This is a speed bound in graph
hops per native event. Native event counts, local event-clock readouts and the
modular parameter have not been identified with one another; the result gives
no laboratory speed or continuum light cone.

[Local event clocks](relational_event_clock/README.md) tick on actual
occupation changes at an observer port while remote repairs continue. Their
completed quantum instruments sum the unresolved waiting process, with
explicit stopping tails. The [joint clock/field construction](joint_clock_process/README.md)
executes two such clocks across elementary observers and preserves their
complete quantum field, intervention and reference outputs. Exact block
elimination supplies the field description's memory terms. At completed
events this experiment admits a four-mark predictor; unresolved waits
retain coherence that the completed event instrument erases. C and E have
the same completed clock law here, yet predict field-plus probabilities
`35/54` and `14/27` after three native attempts under the stated probe.

Clock conversion preserves the full distribution of one observer's ticks
between another's ticks, including conditional field states and references.
Without feedback, each two-A-tick cycle has mean two B ticks and variance
16. Equal mean rates therefore do not imply identical readings. A specified
B-triggered field intervention changes that law; direct native elimination
and conversion through the joint event instrument agree exactly.

The [incident-clock field connection](observer_time_fields/README.md) counts
every attempted repair meeting the observer region, including zero currents.
Its remote/incident split derives a positive effective transfer and a
controlled logarithm in those clock units. Every stationary local field word
and compatible recorded intervention has the same probability in the
reconstructed description. Complete four-read protocols execute on all
28,680 two-particle states of the 240-port federation, with initial labels,
remote evolution, feedback ticks and timeouts retained. Arithmetic and
stopping errors are bounded. Movement-only clocks have a negative mode and
cannot use that positive-transfer logarithm. The
[finite-density extension](observer_time_fields/REFINEMENT.md) connects the
positive transfer to the actual levels-three/five stopped-clock comparison:
its complete four-bit law retains a nonzero connected statistic at density
`1/3`, and native feedback reverses that statistic's sign. Exact exclusion
duals, complete timeout records and source executions with serialized field
counters supply the joint field/refinement certificate. This does not
identify a physical native unitary or an all-level continuum limit.

[Native modular dynamics](relative_modular_dynamics/README.md) identify a
protocol algebra preserved by both event operations and the chosen modular
flow. Their finite covariance defect shows that these operations cannot
simply be equated. The [transfer modular probe](transfer_modular_probe/README.md)
executes a reference-derived modular flow on a reconstructed effective
algebra. The source-derived first-current Palm reference has a faithful
three-state conditional restriction on each existing cross-degree coherent
mode, with weights `(5/13,4/13,4/13)`, Palm mass `13/300`, and modular
frequency `log(5/4)`. These are restrictions to the binary two-load subspace
of a three-port support, not autonomous port factors. An existing coherent
preparation and Y readout give an exact, periodic response to the declared
modular-phase target: `1/2 - sin(s log(5/4))/140` for C and
`1/2 - sin(s log(5/4))/220` for E. This is a candidate subjective modular
time readout for target parameter `s`, not a phase advance generated by native
repair. The successful coherent state is real and every current native Kraus
branch is a real partial injection, so every native-attempt iterate remains
real and leaves the Y outcome unbiased. The periodic response is not a unique
timestamp, and no conversion to native attempts, laboratory duration or proper
time is derived.

The [observer federation](federated_observer/README.md) uses screen levels
with 20, 80 and 320 elementary carriers. Its
[quantum extension](federated_quantum_sector/README.md) preserves complete
source records and arbitrary untouched references, including executed
cross-carrier coherence probes. Local dependency bounds in the
[native action construction](native_action/README.md) control the whole
local quantum and record output under local feedback. They give a
conditional infinite-graph process when neighborhoods, instruments and
initial restrictions are compatible. The
[screen comparison](native_action/README.md#compatible-neighborhoods-of-the-actual-screen)
identifies source-port neighborhoods at levels 1/3, 2/4 and 3/5. For the
15,360- and 245,760-port screens, matched radius-nine initial restrictions
give complete root-observer output distance below `0.000382` at auxiliary
per-seam parameter `s=1/4`. A fixed-total preparation can add a separate
state mismatch. The [incident-clock comparison](native_action/INCIDENT_REFINEMENT.md)
also gives a complete-output distance below `6.656e-9` after four local
incident attempts, with no auxiliary time coordinate. At fixed density
`1/3`, a certified two-time field covariance remains above `0.15332581621`;
full-source histories execute native feedback and checkpoint continuation.
All-level routing compatibility, Lorentz symmetry and an
interacting spacetime scaling limit remain open.

## Native measurement predictions

[Finite native equilibration](native_equilibration/README.md) derives
preparation of the stationary state within a connected fixed-number sector,
with an error bound uniform over untouched references. It does not select
the conserved total or its mixture across sectors. Preparation records
must be discarded for the stated unread convergence theorem.

The [continuing-carrier discrimination protocol](closed_prediction_protocol/README.md)
uses 379 discarded attempts before each of 3,927 trials. Its locked design
retains every failed herald and requires no exact reset state or IID premise
for the actual trials. Conditional-law bounds give a wrong-candidate
probability at most `3/400` under the exact-source model. A further uniform
physical service error allowance of `1/1570800` per complete trial would
suffice for a 1% bound; that physical allowance is not authenticated.
The producer checks consecutive trials and full reference outputs, rather
than presenting computed probabilities as measured trial counts.

Additional native predictions include
[event-clock readouts](native_relational_readout/README.md),
[one-history circulation](native_circulation/README.md),
[event-count circulation](event_circulation/README.md), and
[native-reference modular probes](native_modular_reference/README.md).
Their stationary preparations, record access, coherent services and attempt
or event parameters are stated individually. No public laboratory dataset
is identified with these native protocols.

## Empirical comparisons

The public oscillator comparisons execute a **supplied experimental model**
through the source interface. Their agreement tests that conditional model;
it does not derive its oscillator action or apparatus from native repair.
Measurement coordinates, preparation assumptions, accepted-trial information
and physical service errors remain explicit limitations.

The [executable Gaussian readout comparison](scalar_observations/experimental_readout/EXECUTABLE_READOUT.md)
uses two calibrated displaced-thermal preparation parameters and executes
100 reserved complex readouts. Its primary per-spin-mean RMSE is `0.06647`.
The reported ideal occupation probability bound is below `3.82e-7` at that
preparation. This is a retrospective calibration/evaluation split, with
unresolved coordinate conventions and experimental postselection.

The [driven anharmonic comparison](anharmonic_observations/README.md)
predicts 707 positive-duration detector means from a declared pulse sequence,
preparation calibration and dephasing channel. The
[changed-potential comparison](anharmonic_transfer/README.md) evaluates five
preparations against 125 positive-duration derived observations. Its primary
RMSE is `0.2876427868` in dimensionless two-point position-estimator units;
material trajectory discrepancies remain. These units differ from binary
probability and spin-mean errors and must not be compared directly.

[Readout identity controls](readout_controls/README.md) isolate
residuals that a changed oscillator action cannot remove under the stated
measurement interpretation. The [two-order characteristic constraints](characteristic_constraints/README.md)
require neither a Gaussian state nor an evolution law. They apply quantum
Gram positivity to both displacement orders at the unresolved coordinate extents 4 and 5. Their strongest rational witnesses require a deterministic
correction to at least one complex mean of at least `0.02439` or `0.03129`,
respectively. These are descriptive consistency constraints on noisy means,
not estimates of detector bias, statistical rejections or OPH-specific
empirical confirmation. Accepted counts, coordinate calibration and
apparatus error bounds are needed for a statistical assessment.

A randomized Gram-witness protocol permits preparation drift, feedback and
quantum memory between trials. Fresh settings, complete outcomes and
conditional detector-error bounds replace an IID preparation premise. Its
source demonstration retains 100 two-trial outcomes. The deposited means
select witnesses for a future test; they do not supply the random-setting
trial records required to apply its statistical bound.

The [calibrated measurement protocol](calibrated_measurement/README.md)
provides explicit calibration and sampling bounds for future counts.
[Preparation and detector probes](empirical_diagnostics/README.md) identify
measurements that distinguish otherwise unresolved preparations.
[Identifiability checks](identifiability_checks/README.md) retain calibration
and geometry counterexamples and empirical controls. Software fixtures are
labelled separately from public observations.

## Supplied field models and reconstruction tests

The [common interacting source](common_source/README.md) places electric
flux, charged matter and a scalar field on one continuing register. Its
[finite executable specification](common_source/COMPOSITION.md) joins
preparation, readouts, retained feedback and native actuation. Certificates
cover [mixed-state continuation](common_source/OPEN_HISTORY.md),
[complete numerical outputs](primitive_observer/COMPLETE_CERTIFICATE.md),
[occupation refinement](common_source/OCCUPATION_REFINEMENT.md), and
[finite-resolution measurements](common_source/SMOOTH_REFINEMENT.md) under
their stated domains. Numerical, ideal-target and physical service errors
are separate quantities.

[Shared scalar observations](scalar_observations/README.md) compare fixed
field smearings across interacting grids with complete measurement and
occupation bounds. The [spatial interacting family](spatial_interaction/README.md)
executes a supplied repulsive Bose action with controlled infinite-volume
observation laws at fixed lattice spacing. The
[canonical scalar regulator](relativistic_scalar/README.md) supports
number-changing quartic interactions, an untruncated finite-volume gap
bound and a separate free Klein–Gordon refinement. None selects its mass,
geometry or interacting relativistic action from native repair.

[Record-only reconstruction](observer_reconstruction/README.md),
[interacting reconstruction](interacting_reconstruction/README.md), and
[coupled calibration](coupled_reconstruction/README.md) recover parameters
and geometry within supplied action families, then predict complete
measurement histories with disturbance and feedback. Basis, preparation,
instrument access and calibrated clock hypotheses remain explicit.
Identification of the reconstructed field geometry with native observer
geometry is unresolved.

## Run and verify

From the repository root, use Python 3.11 or newer. The exact native command
uses only the standard library:

```sh
python3 -m model.observer native --check model/observer/results/native.json
python3 -m model.observer run --trials 10000 --seed 42
```

The sequential sample run is a software model with its declared preparation
and controller. A Python invocation is not assigned a physical duration.
For quantum constructions and the full verification suite:

```sh
python3 -m venv .venv-model
. .venv-model/bin/activate
python3 -m pip install -r model/requirements.txt
python3 -m model.verify
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model
```

`model.verify` checks artifact hashes and replays retained results. Individual
current constructions can also be checked directly:

```sh
python3 -m model.native_action --check model/native_action/results.json
python3 -m model.joint_clock_process --check model/joint_clock_process/results.json
python3 -m model.closed_prediction_protocol --check model/closed_prediction_protocol/results.json
python3 -m model.characteristic_constraints --check model/characteristic_constraints/results.json
```

Use `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` for reproducible resource usage
on numerical runs. Package READMEs document their hypotheses, tests and
output options. Result writers refuse to overwrite existing evidence files.

## Research obligations and navigation

The central obligation is an observational correspondence: compatible field,
clock and intervention laws must recover the physical regimes observers
measure. A finite algebra, a clock comparison or a graph dimension alone
does not supply that correspondence. The remaining work includes compatible
federation/refinement maps, controlled interacting scaling limits, a physical
apparatus mapping and discriminating data with justified errors. Supplied
benchmark structure remains a premise until its connection is demonstrated.

- [Observer-relative recovery](OBSERVER_RELATIVE_RECOVERY.md): current native
  action, clock, field and federation results and their limits.
- [Source selection](SOURCE_SELECTION.md): implemented inputs, derived
  consequences and physical selection obligations.
- [Construction guide](CONSTRUCTION.md) and
  [construction paper](../paper/executable_observer_construction.pdf):
  contracts, acceptance criteria and finite-to-physical theorem routes.
- [Model–theorem premise map](PREMISES.md): hypotheses and mathematical
  attachments; [candidate register](CANDIDATES.md): declared alternatives.
- [Mathematical foundations](../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md) and
  [working terms](../docs/INTERPRETATION_AND_TERMS.md): compatibility,
  scheduling, observer access and the distinction between history and state.
- [Supporting executable evidence](../code/README.md),
  [reproduction guide](../REPRODUCE.md), and
  [scientific registers](../tracking/README.md): source calculations,
  reproduction contracts and physical observation status.

Software is licensed under [Apache-2.0](../code/LICENSE). Documentation and
data follow the repository's [licensing map](../LICENSE).
