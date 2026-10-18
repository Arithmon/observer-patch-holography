# OPH model construction from local repair

This is the construction guide for the [model implementation](README.md).
The [repair-law guide](CANDIDATES.md) explains which parts of the repair law
are forced and which are free; the generated [premise map](PREMISES.md)
connects the constructions to the theorem library.

**Build the simplest shared repair process that satisfies the required
constraints.**

## What is being built

OPH proposes that the physical world is what overlapping observers agree on
once they have repaired their disagreements. The construction makes that
proposal executable. It fixes a finite carrier, a local repair law, observers
that read records and act on them, and the services those observers use. It
then runs the specification and derives what its observers can measure. The
aim is a specification whose observers' joint statistics admit standard
physics as an effective description.

Mathematically, an instance is a tuple

```text
(carrier H, accessible algebra A, repair instruments {Φ_{e,r}}, scheduler σ,
 observer instruments and feedback maps, preparation domain, declared services).
```

Each seam `e` carries an instrument whose outcomes `r` are committed records.
An observer is a finite tape of such records together with a controller that
chooses the next operation from the tape alone. Everything an observer can
report is a probability computed from these maps.

The working method has four steps. Build the simplest instance that meets the
stated requirements. Derive its consequences exactly. Keep every supplied
input visible. Find the observations that distinguish it from the
alternatives the requirements leave open. A finite implementation is complete
when its declared observation, record and feedback contract runs and passes
its checks. Proving that self-consistency uniquely selects that
implementation is a separate question and is not a prerequisite for
constructing it. A complete physical model additionally needs compatible
physical identifications and tests on the same process.

## Completion target and acceptance criteria

**Build an executable realization of the OPH observer algebra whose
compatible inter-observer measurements reproduce the observed Universe.**
Their joint observations must admit a nontrivial 3+1-dimensional QFT as a
controlled effective description, including the observed matter and
interactions, and classical physics and dynamical gravity under appropriate
coarse-graining. Four-dimensional QFT and 3+1-dimensional spacetime refer
to the same dimension count of that description.

**No fundamental QFT or spacetime recovery is required.** The target concerns
what observers can measure and jointly interpret. Fields, particles and
spacetime geometry can be effective descriptions of those observations,
without being primitive objects or an underlying observer-independent
substrate. Implementation data structures do not change that ontology.
Success requires agreement of the accessible joint and conditional
statistics and intervention responses in each description's validity
regime, with controlled errors and corrections. It does not require every
effective description to remain valid at arbitrarily fine scales.

Meaning assignment connects actual preparations, instruments, interventions
and records to shared physical quantities. Consensus and repair must preserve
compatible joint and conditional statistics on overlapping observer domains.
This is the proposed operational sense in which an objective reality emerges.
It does not follow merely from observers agreeing on an interpretation.
See the [interpretation contract](../docs/INTERPRETATION_AND_TERMS.md#observer-measurements-and-operational-objectivity).

| Completion obligation | Acceptance evidence |
| --- | --- |
| Executable observer algebra | Explicit state, operation, preparation, measurement, retained-record and feedback maps; their algebraic identities and observer-access limits; accounted memory and service assumptions; actual execution of the claimed protocols. A finite two-event implementation discharges its own contract, not every observer or federation interface. |
| Shared operational meaning | Interpretation maps specified independently of the validation outcomes, compatible across observers and refinement, and preserving the declared intervention statistics and quantum correlations. Shared labels alone are insufficient; incompatible unperformed measurements need no common table of classical outcomes. |
| Effective quantum-field description | A common refining family of observer statistics and interpretation maps whose errors relative to a nontrivial 3+1-dimensional QFT tend to zero in the declared continuum domain, including the observed matter and interactions. The comparison concerns measurable correlations and intervention responses, not recovery of fundamental field objects. Where an effective vacuum-sector mass gap is defined, it is nonnegative and may be zero. |
| Effective classical and gravitational descriptions | The same observer construction yields observations that admit compatible classical, geometric, stress and clock descriptions in their stated regimes, reproducing the gravitational predictions with controlled corrections. A fundamental metric or manifold is not a prerequisite. Separate calculations on a supplied geometry need a map to these observer statistics before discharging this obligation. |
| Physical predictions | Explicit parameter/input ancestry, calibration and uncertainty; predictions and rejection rules fixed before independent comparison data; agreement with the observed statistics, including alternatives and distinguishing tests. Reusing a calibration datum is not an independent prediction. |

Exact equality of limiting observable laws is the mathematical target.
Finite simulations and experiments are assessed with declared approximation
and measurement errors. They do not prove arbitrary empirical precision or
force a particular random outcome in an independently repeated experiment.
Agreement must be constrained by measurements; a consistent world model can
describe different physics from nature.

### The observable correspondence to establish

Fix an observer region, preparation family, accessible instruments, feedback
rules and finite resource domain before comparing descriptions. A physical
dictionary maps these operations and their complete outputs to the proposed
effective field description. It must include local clock comparisons, signal
records and field readouts together. Different observers' dictionaries must
agree on their shared records and allowed joint interventions. Choosing
seconds and metres sets units; it does not determine these response laws.

For a protocol `p` in that fixed menu, let `rho_native(p)` and `rho_eff(p)`
be its complete record/quantum outputs, including failures, unresolved waits
and any admitted untouched reference. Map both outputs to one declared
comparison algebra with decoders `D_native` and `D_eff`. The finite target is

```text
sup_p (1/2) ||D_native(rho_native(p)) - D_eff(rho_eff(p))||_1 <= epsilon.
```

The supremum also covers the related initial states and references declared
by the contract. For classical records this is total variation of the full
joint law. Matching individual expectations, two-point functions or clock
means is weaker. Conditioning on a rare successful outcome needs its own
probability-dependent error bound; it cannot erase failed branches from the
comparison. Decoders must retain the observations being claimed: mapping
both outputs to a constant is not an observation-faithful correspondence.

The finite source results provide parts of this contract. A physical QFT
attachment must additionally supply an operational domain with three spatial
and one temporal dimension, the appropriate covariance/locality and state
conditions, and measured interaction responses. Its approximation errors
must be controlled across the claimed domains or refinement. A finite
positive transfer and its spectral logarithm alone do not provide those
hypotheses. The [premise map](PREMISES.md) records which theorem interfaces
are instantiated and which require extra inputs.

The empirical step uses the *same* dictionary and calibrated parameters on
independent observation protocols. Existing public data can supply that test;
a physical emulator or new laboratory is unnecessary. An assumed
experimental outcome is a conditional forecast and cannot validate the
dictionary. The [public-data audit](../docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md)
keeps the ACT comparison separate from native-observer derivation.

Hardware, protocol and implementation choices are open design variables.
No uniqueness theorem for the chosen rule is a prerequisite. For the current
P/N construction, closure-derived P and N are accepted inputs under their
declared contracts; that assumption does not discharge the other physical
interfaces. Supplied background geometry, action, state, interpretation or
clock data must stay visible until the same construction accounts for them.

These obligations can advance in parallel. Passing CI verifies the retained
software contracts; it does not mark the physical target complete. The
current native observer, P detector and source-compiled two-site scalar
apparatus are completed finite components under their stated premises.
The controlled effective interacting-field and gravitational descriptions
of the observer statistics, and independent empirical completion, are open.

## Why this specification

Every choice below is either forced by a stated requirement, chosen as the
simplest member of an explicitly classified family, or declared as an input
with its alternatives kept. The [repair-law guide](CANDIDATES.md) gives the
law-level arguments in full.

| Choice | What the model uses | Reason | What stays open |
| --- | --- | --- | --- |
| Carrier | Twelve ports and thirty seams of the source icosahedron; integer loads with total two | Declared by the architecture premise; total two is the smallest sector where a single load waits or moves and two loads collide | Selection of the architecture from an unrestricted menu; extended carriers and refinement |
| Classical seam rule | Endpoints become floor and ceiling of their mean; a fair coin breaks ties | Forced by conservation, one-step balance, a tie-break that reads no load, and endpoint symmetry | Nothing on its domain; the four requirements are premises |
| Record law | C executes; E is the comparison law | C records the classical rule's own random variables, the coin and the current, and nothing more | Which coherence `g` nature keeps; the admissible laws form one convex family |
| Scheduler | Uniform independent seam attempts | The symmetric baseline; prescribed seams serve as controls | Identification of attempts with a physical clock |
| Observer | A finite record tape and a controller that reads committed records and chooses the next seam | The smallest object that records, predicts and acts on its own records | Larger observers, record delivery and refinement |
| Services | Coherent preparation, readout, mixers, record storage and compiled source words | Declared and costed in every error budget | Physical realization and authenticated service errors |
| Closure inputs | P and N as derived by the closure equations | Accepted under their declared contracts | Physical endpoint corrections |
| Field benchmarks | Supplied scalar, gauge and matter actions on a common source | Conditional comparisons with explicit errors | Derivation from native repair |

**The record law in one paragraph.** Requiring seam locality, the exact
classical law, local number-phase covariance and endpoint symmetry leaves a
convex family of quantum record laws. Each member is a table of five
coherence multipliers, one positive matrix per transported current. Four
members stand out: C0 records only the current, C adds the coin, E adds the
seam parity and T adds the seam total. Exactly E and T are idempotent, and
every law that keeps an empty seam coherent with a full seam completes to E
under repetition. All members share every classical observation. One
coherent interference experiment reads the multiplier `g` directly, with
predictions `3/8 ± Re(g)/2` at phases `0` and `π`. The model executes C because
its records are exactly the classical event data. The
[record-law classification](record_laws/README.md) proves these statements
on the 78-state carrier and checks them against the executed C and E
instruments.

**The observer in one paragraph.** An observer here is the smallest system
that keeps records and uses them. The [native observer](observer/README.md)
commits the coin and current of one repair, forecasts the next current from
that record and chooses its next seam from the current alone. Its controller
sees no hidden load. Because every admissible record law records the current,
the same controller runs under each of them. Larger observers compose such
tapes, federations join elementary carriers, and detectors attach to the
record through declared couplings.

## Current executable result

The [primitive observer runtime](primitive_observer/README.md) executes
complete instruments on mixed states, committed records, feedback and
untouched references. Source-derived native dynamics and supplied physical
benchmark programs use this same execution interface. Their preparations,
clocks, interpretation maps and error budgets remain explicit.

### Local construction and its extension

The [coherent screen observer](coherent_screen_observer/README.md) executes
one quantum extension on the actual matched level-three and level-five
screens. The input is a logical qubit on two root ports, with any untouched
reference, and a separate Bernoulli environment. Its two native incident
ticks, seam-category records, local phase feedback and complete reduced
quantum output compose without a state reset. A bounded operator dual
evaluates the common 222-port channel; the existing complete-output locality
theorem controls its comparison with both full screens. This discharges a
finite coherent screen comparison under that preparation contract, while
retaining the all-level and physical interpretation obligations.

The local proof obligations start on one twelve-port carrier: specify its
accessible algebra, derive repair and record laws, compare its local clocks,
and reconstruct spatial relations from those observations. Observer regions
may overlap inside the same carrier. A federation is useful for tests involving
more carriers, larger distances or refinement; none is required merely to
define an observer or compare two local record views.

An extended spacetime claim requires a compatible family of operational
domains and an effective regime with controlled errors. The existing
spherical federation is a candidate for that extension. Its geometry and
refinement must earn their physical interpretation through observations;
they are not prerequisites imposed on every finite construction. A single
carrier can establish a finite contract without establishing that extension.

Action representations also have different domains. The endpoint Kraus
instrument determines an exact reversible representation of its finite
histories. Choosing a unitary completion and a Hermitian logarithm is an
encoding choice. It supplies no prediction for an extra intermediate-time
measurement until such access and evolution are included in the observer
contract. A valid construction need not be uniquely selected; its advertised
observations must follow from its declared operations.

The binary current calibration and quantum clocks share the same 66-state
sector. Native propagation, occupation covariance and the slow spatial Gram
can therefore be compared without using a doubled-load preparation. A
separate current-fluctuation action is derived in one port's tick count, with
an exact Palm-conditioned triangle variance `3m/5`. Actual host-counter
records test a software clock correspondence with held-out blocks and a
workload control. The [joint clock bridge](native_clock/CLOCK_BRIDGE.md)
derives incident-hit/current-count laws and finite conversion errors from any
initial configuration, including zero hits and a separate feedback-count law.
Its reanalysis preserves the host calibration and workload rejection. These
results retain their total-two preparation domain; neither a traceable lab
clock nor the full-screen field protocol is calibrated by that comparison.

### Native action, federations and observer-relative clocks

The [quantum federation construction](federated_quantum_sector/README.md)
joins identical elementary carriers using the declared source wiring. Native
Kraus maps preserve signed-difference operator blocks, including individual
outcomes and reference correlations. Its complete cross-carrier experiments
retain coherent interventions, all records and conditional quantum states.
The [local action](native_action/README.md) is fixed by these instruments:
`H_native=sum_e(I-Phi_e)`. Its positive local terms reduce to exclusion plus
coherence killing. The original discrete transfer is `I-H_native/M`; the
optional per-edge Poisson extension is `exp(-s H_native)`. This is a local
Liouville action, not a selected relativistic Hamiltonian or physical clock.

The same construction controls truncation of distant ports for complete
local record, quantum and reference outputs. For an observed region A,
maximum degree Delta and retained radius R, the Poisson-parameter error is
at most `min(1, |A|(Delta s)^(R+1)/(R+1)!)`; the fixed-count bound is given
separately. A conditional infinite-federation process follows on a fixed
bounded-degree graph with compatible local states and operations. The
[declared screen refinements](federated_observer/REFINEMENT.md) require
compatible access and dynamics; graph rank or finite positivity alone does
not establish a relativistic spacetime limit.

Actual source neighborhoods are identified at levels 1/3, 2/4 and 3/5,
with preserved local port labels and explicit complete-output comparison
bounds. This supplies finite cross-size observation maps rather than assuming
that parent/child coarse averages intertwine dynamics. Matched initial
restrictions are necessary; a uniform global total-two preparation dilutes
locally as the screen grows. An all-level routing argument and a controlled
nontrivial physical scaling regime are open.

The [incident-clock attachment](native_action/INCIDENT_REFINEMENT.md) compares
these actual screens after the same number of local repair opportunities,
without supplying a duration per operation. Its stopped-word identity and
source boundary counts bound the complete quantum and record output, including
local feedback and untouched references. Four incident ticks give a pair
allowance below `6.656e-9`. A compatible density-`1/3` preparation supplies a
certified nonzero two-time field covariance, while full-screen executions
retain the actual seam scheduler, current records, feedback repair costs,
timeouts and serialized continuation. This closes the finite clock/refinement
connection; the all-level physical interpretation obligations remain separate.

The [finite-density field attachment](observer_time_fields/REFINEMENT.md)
uses that same preparation, root clock, current support and native feedback
to certify a complete four-bit field law, including every finite-cap timeout
prefix. The positive incident transfer preserves this law with the Bernoulli
stationary boundary or with each screen's own canonical fixed-number
restriction. The latter removes a product-preparation premise: native
connected exchange derives the unique diagonal equilibrium conditional on
the conserved total. Its exact degree-two exclusion dual retains conservation
correlations without enumerating the full binary space. At density `1/3`,
feedback reverses the connected four-field statistic on both screens, while
their complete field laws have a resolved finite-size difference. Two field
counters and controller memory survive serialized source continuation;
arithmetic, cutoff, spatial and preparation differences are accounted
separately. Selecting the conserved density and realizing equilibrium
preparation remain inputs, with no large-screen mixing-time claim.

The [temporal-record representation](relative_observer_geometry/README.md#temporal-hilbert-space-reconstructed-from-field-records)
reconstructs a seven-dimensional cyclic positive transfer from exact field
correlations in the elementary total-two, two-port-clock domain. Exact
flatness proves its recurrence under the established positive-transfer
class; ordinary samples do not certify that rank. On the distinct full-screen
product preparation and twelve-port root clock, bounded correlations certify
a two-vector temporal subspace and positive transfer compression. Its Ritz
values are not closed modes or a physical energy/clock identification.

The [local event clock](relational_event_clock/README.md) counts actual
occupation changes at one observer port. Its stopped instrument retains
remote repairs during the wait, finite timeouts and conditional quantum
outputs. The [joint clock process](joint_clock_process/README.md) compares
two such clocks and an observer-chosen field algebra on the same continuing
cross-carrier history. Exact elimination of the hidden factor gives a
record-conditioned quantum memory recurrence, valid with arbitrary untouched
references. Completed events in this specified experiment admit a four-mark
measure-and-prepare description; unresolved waits retain quantum information.
Instrument-valued clock conversion also preserves the count of B events
between A events. Direct native elimination and elimination through the
joint event instrument agree, including B-triggered field feedback. The
passive two-A-tick renewal cycle has B-count mean two and variance 16;
equating mean rates would discard observable clock fluctuations.
Neither reduction identifies native attempts or event counts with modular
flow. The [recovery contract](OBSERVER_RELATIVE_RECOVERY.md) specifies the
remaining observational and clock correspondences.

The [native field identities](common_source/NATIVE_FIELDS.md) derive exact
continuity, current noise and mean diffusion under the uniform scheduler.
[Quantum transfer](native_quantum_history/README.md),
[record-clock laws](native_clock/README.md),
[joint current/phase readouts](native_relational_readout/README.md),
[replicated fluctuations](native_fluctuations/README.md) and
[one-history circulation](native_circulation/README.md) provide additional
controlled finite or fluctuation descriptions of that source. Their
scheduler, accessible records and limit hypotheses are distinct. Replication
is not a spatial refinement, and a Gaussian limit does not select a
fundamental Gaussian field action.

### Continuing measurement predictions

[Native equilibration](native_equilibration/README.md) replaces an exact
stationary preparation by a finite sequence of discarded repairs within a
declared fixed-number sector, uniformly over untouched references. The
[sequential prediction protocol](closed_prediction_protocol/README.md)
applies that preparation before each trial on the same continuing carrier.
It bounds the dependent complete transcript without an IID or exact reset
premise. The locked design uses 379 discarded attempts before each of 3,927
trials and three unread attempts in each useful phase probe. Its exact-source
wrong-candidate probability is at most 0.75%; a separate, unauthenticated
physical service allowance would suffice for a 1% target. Failed heralds
remain trials. The producer executes two consecutive trials per candidate
and reference controls, not the entire proposed experimental campaign.

The [native interface](native_interface/README.md) obtains the coherent
readout from a source mixer, native herald and native pointer. The
[compiled interface](compiled_native_interface/README.md) verifies the
source control against its ideal target, including padding. Source matrix
comparisons, floating execution diagnostics and physical device errors are
separate. The [finite-law discriminator](predictions/README.md) and its
[encoded realization](observer/source_discrimination/README.md) supply
additional complete protocols under their own preparation and service
contracts. Mean loads and native classical histories assign the same
probabilities under every admissible record law; the discriminator reads the
coherence `g` that separates them. No public laboratory dataset is
identified with these native-law experiments.

The [public Bell-record comparison](predictions/BELL_RECORDS.md) tests the
scope of a classical interpretation. The diagonal native branch obeys
`S<=2` under independent settings and local factorized responses. Replaying
the deposited Hensen records gives 196 wins in 245 eligible trials and
`S=2.42250`; the published fixed-count calculation has a tail bound near
`0.03917`, conditional on its setting calibration and experimental access
contract. This retrospective null rejection excludes that extension, not
the full quantum C/E construction, and supplies no lab-clock identification
or a prospective validation of OPH.

The [measurement bounds](measurement_bounds/README.md),
[empirical diagnostics](empirical_diagnostics/README.md) and
[action-independent readout identities](readout_controls/README.md)
constrain the relation between a declared apparatus and data. They do not
supply an OPH-specific empirical confirmation. Missing raw counts are not
inferred, and calibration is not refitted from held-out outcomes.

### Supplied interacting work and quantum continuation

The [common interacting source](common_source/README.md) carries electric
flux, a charged CAR pair and two neutral scalar sites on one work register.
Its matter readout controls a feedback operation; the scalar readout controls
a native repair. The complete history retains 1,056 native/quantum records
and their terminal states. Nonzero feedback contrasts have finite and
[infinite-occupation certificates](common_source/INFINITE_HISTORY.md).
The action, initial state, coherent services and Hamiltonian model time are
declared inputs; compiling them does not identify native repair with that
Hamiltonian. The [theorem attachments](common_source/DERIVATION.md#compatible-theorem-use-and-open-physical-connections)
state the action and physical identification hypotheses.

The [mixed-state continuation](common_source/OPEN_HISTORY.md) implements
an explicitly supplied unobserved environment. Discarded phase draws occur
at declared boundaries; their labels cannot enter feedback or postselection.
It preserves visible histories, disturbed states, reference correlations and
checkpoint continuation. Environment variance, scheduling and physical
service errors remain inputs. The [typed protocol](common_source/COMPOSITION.md) and its
[runtime](common_source/protocol.py) execute smooth reads and feedback without
resetting the work register.
[Quantum native continuation](common_source/QUANTUM_NATIVE.md) and
[encoded recovery](common_source/ENCODED_MEMORY.md) attach full native
instruments and active logical memory under explicit code-service contracts.

The [native-clocked common observer](clocked_common_observer/README.md)
uses an actual current event to trigger that work's read-and-feedback suffix.
Work evolves during the wait. Its 128-slot execution retains 4,097 stopped
and timeout outputs, with complete finite-target error at most 0.000019,
including the specified binary arithmetic. All-input channel and tail bounds
control the finite composition. Work duration per attempt, the supplied
action and physical calibration remain separate inputs.

The [causality identities](common_source/CAUSALITY.md) control committed
records and permitted intervention order. Error bounds apply to their stated
instruments and preparation domains.
The [complete-output comparison](common_source/READOUT.md#complete-comparisons-with-the-actual-finite-pointer)
shows why agreement of one sign probability need not imply a close quantum
output: the sharp detector has a long occupation tail. The [occupation refinement](common_source/OCCUPATION_REFINEMENT.md) and
[smooth-detector refinement](common_source/SMOOTH_REFINEMENT.md) distinguish
their pointer domains. The
[smooth-history bound](common_source/SMOOTH_HISTORY.md) has a different
instrument and a specified preparation contract. The
[direct numerical certificate](common_source/OUTPUT_CERTIFICATE.md) and
[complete observer certificate](primitive_observer/COMPLETE_CERTIFICATE.md)
account for their actual executed arrays; neither is a physical service
certificate. [Observer readouts](common_source/OBSERVER_READOUTS.md) and
[thermal-clock comparison](common_source/THERMAL_CLOCK.md) retain the
declared reference state, access and clock assumptions.

### Effective field benchmarks and empirical readout

The [two-site source compiler](observer/source_scalar/README.md) executes
preparation, evolution and both detector outcomes from source words, with
complete occupation and splitting bounds. Its separate action and apparatus
do not inherit another spatial family’s continuum result.

The [shared scalar observations](scalar_observations/README.md) compare
complete field detector, preparation and feedback laws across a fixed-band
spatial family. Their [occupation proof](scalar_observations/OCCUPATION_DERIVATION.md)
and [receipt](scalar_observations/occupation_results.json) control projection
onto finite source codes. The exact-gradient endpoint is
`GalerkinParameters(grid=None)`: increasing quadrature resolution does not
add modes. Extending the field needs additional modes, matched observation
maps and controlled omitted-mode effects.

The [interacting band extension](interacting_band_extension/README.md)
performs one such extension, adding the cosine and sine at momentum
`(3,0,0)`. The common quartic action fixes their new couplings. Averaging
their specified oscillator blank induces a low-mode mass-squared shift
`g/4`; the comparison includes that shift instead of equating bare and
matched parameters. A block-Duhamel bound controls complete low-mode
observations through measurements and feedback, with the added modes
retained between reads and traced only at the comparison output. Its
finite occupation domain and selected anisotropic band remain explicit.
This does not remove all omitted modes or prove relativistic locality.

The [relativistic scalar regulator](relativistic_scalar/README.md) executes
an eight-site three-dimensional quartic action with an unrestricted-
oscillator comparison. Its ideal record TV error is below 0.000928 at the
stated finite parameters. The supplied positive mass gives a finite-box
gap bound, and the free spatial family has a controlled Klein--Gordon limit.
An interacting relativistic limit requires additional spatial/energy control
and parameter matching. The [spatial interaction model](spatial_interaction/README.md)
has controlled growing-volume observation laws and a nonzero interaction
witness for a supplied repulsive Bose action. Growing a box at fixed spacing
does not remove its spatial regulator.

The [observer reconstruction](observer_reconstruction/README.md) recovers
free propagation and finite geometry from executed detector laws and a
specified reference-clock readout, with withheld predictions and controlled
occupation errors. Its quadratic action and mode apparatus are supplied.
A cube adjacency graph does not by itself determine physical dimension;
the [identifiability checks](identifiability_checks/README.md) include a
planar control and alternative actions invisible to a restricted apparatus.

The [published readout comparison](scalar_observations/experimental_readout/README.md)
calibrates a Gaussian response relation on zero-duration trapped-ion data;
100 off-axis means test that relation with spin-mean RMSE 0.07639. The
[executable detector](scalar_observations/experimental_readout/EXECUTABLE_READOUT.md)
contracts source preparation, rotation and detector words under the declared
laboratory model. Published coordinate conflicts, postselection and
normalization limitations remain explicit. This component postdiction
establishes neither the seven-mode dynamics nor a physical clock.

The [occupation-control construction](occupation_control/README.md) bounds
the continuous polynomial-path defect, occupation boundary, endpoint jumps,
thermal initialization, displacement and pointer for the supplied pulse
apparatus. It includes arithmetic under its stated IEEE-754 contract, not
the separate native compiler or laboratory error. At 64 levels the case-a
estimator radius is about 9.91e-8 at 0.4 ms, 0.0163 at 2 ms, 0.211 at 5 ms
and 1.68 at 10 ms; the 32-level conditional bound is unavailable at the last
two durations. Agreement between cutoffs is not an error certificate.
The independently controlled unrestricted-versus-exact-finite radius can
compose with a valid native execution bound.

The [identification design](identification_design/README.md) computes the
response rank and nullspace for a declared correction family. Four-level
access cannot detect its `binom(N,4)` term; coherent levels four and five
expose those directions, apart from the global-energy gauge. A 41-pulse
source preparation and phase read distinguish two supplied actions, while
occupation commitment before the read erases the distinction. The efficient
sequence is to identify blind directions, choose a separating intervention,
then compile its locked measurement. A component sampling bound does not
transfer to a noncommuting field model or physical apparatus.

### Source execution and numerical contracts

Exact finite-code identities give controlled occupation parity in two native
pulses and dimension-64 number drift in seven. They do not import unrestricted
oscillator commutators at a cutoff. Parameter-keyed compilation caches share
source words, never observed trajectories, fitted parameters or prediction
outcomes. Canonical verifiers execute every declared calculation and
fail on any violated contract.

Pulse matrices use the fixed-guard polar representation of the source
library. Raw, incremental and corrected-target checks bound the representation
change; no state or measurement outcome chooses that correction. Density
matrices are not renormalized to repair execution. These numerical controls
are distinct from physical gate accuracy. Frozen evidence retains its
original bytes and scope; live receipts reproduce from canonical producers.

## Native repair and its observer contract

Use the twelve-port source, with 78 nonnegative integer configurations of
total load two. Its 66 binary configurations form an invariant shell. For
an oriented seam `(a,b)` and fair placement coin `c`, native C uses

```text
d = x_a − x_b,
j = [d − (2c−1)(|d| mod 2)]/2,
y = x + j(unit_b − unit_a).
```

The load is conserved. Its quantum current instrument has Kraus operators
`L_(c,j)=(1/sqrt(2)) Σ_(x:j(x,c)=j) |y(x,c)><x|`.
The record is `(c,j)`. C is one member of the admissible family in the
[record-law classification](record_laws/README.md). Keeping the pre-repair
contrast instead gives E, dropping the coin gives C0, and recording both
endpoint loads gives T. All of them have identical classical transitions.
The classical rule fixes the diagonal of every member and none of its
coherences.

The [native observer](observer/README.md) first applies C on `(0,1)`. Its committed first record alone
selects `(0,5)` when `j1>0`, otherwise `(0,1)`, and supplies the preissued
forecast of a positive second current. No source loads or future outcomes
are accessible to that controller. Full simulator checkpoints and external
occupation audits are distinct from the observer's record access.

The [derivation](observer/DERIVATION.md) gives one fixed CPTP map on the
phase/tape/source carrier and shows classical-prefix preservation. Genuine
feedback can reduce coherence between different committed record values while
preserving their classical probabilities. Halting is an identity operation
on the terminal sector, not a claim of unlimited continuing memory.

Exact controls separate record information from knowledge of the chosen
action. Native-record prediction improves on the action-known prior and a
within-action shuffled record. Native feedback gives final port-0 occupation
`35/264`, versus `1/6` when the selector is disconnected or independently
shuffled with the same action frequencies. Occupation here is an external
audit; the native current remains the implemented observation.

## P-dependent response and finite capacity

For this design, **P and N are accepted as derived from OPH's self-reference
equations**. Apparatus, protocol and simulation hardware remain construction
choices. The [detector implementation](observer/p_response/README.md)
evaluates the existing forward-closure P input and source-only
structured-running map `α(P)`, then sets `g(P)=sqrt(4π α(P))` in natural units.
The observer does not fit P or N to its desired output.
Its chosen charge assignment is `Q_v=g(P)N_v`, so native current j transfers
charge `g(P)j` while conserving the total.

A prepared path qubit uses the first current to orient a calibrated loop
exposure χ. A Y measurement writes `y=±1` with

```text
Pr(y | j1,P,χ) = [1 + y sin(g(P)j1χ)]/2,
E[g(P)j2] = (5g(P)/264) sin(g(P)χ).
```

Only y chooses the second seam. At χ=1 the mean native second current is
approximately `0.00564883566843`. Disconnecting the detector or independently
shuffling y removes the mean response while preserving its half/half outcome
and seam frequencies. Weak exposure gives
`E[g j2]=(5π/66) α(P)χ + O(χ³)`: P affects an actual response magnitude.
The loop geometry can separately use `A_loop=K_loop P ell_star²`.

The chosen charge, prepared probe, calibrated flux, readout and controller
are specified apparatus. This is a software execution of that conditional
design, not a report of fabricated electromagnetic hardware. The inherited
source-only α branch retains its documented physical endpoint status;
executing the adapter does not turn its output into an independent precision
measurement. The [endpoint check](predictions/endpoint_audit.py) retains
the approximately 300 ppm inverse-alpha deficit of that incomplete output;
a zero-correction precision endpoint is therefore unsupported. It does not
reject other OPH constructions or fit P to the target. P does not alone
determine a whole observer's Hilbert dimension,
local exposure, temperature or event rate.

The [capacity calculation](observer/CAPACITY.md) reuses the existing
correctable-public-record theorem. Under the binary preparation, 11 terminal
native words and 16 reachable prefixes allow a four-bit classical prefix
encoding. On the full 78-state source those counts are 17 and 24. These are
reachable record counts, not the dimension of the whole machine. The
P-detector's richer audit has 24 terminal words and 37 prefixes; its actuator
receives a single bit. Program, source, probe and service capacity are
additional.

For a declared capacity carrier of dimension D, `N=ln D` bounds a correctable
classical record code by `ln M0≤N`. The cosmological relation
`N=3π/(Λ ell_star²)` needs the stated physical capacity/area identification.
A finite controller does not establish global saturation or cosmic closure.
See the existing [capacity implementation](../code/capacity_readback/README.md)
and [N-closure theorem](../paper/tex_fragments/OBSERVERS_N_CLOSURE_PRO5.tex).

## Geometry, response and local event clocks

The native current process supplies a finite port geometry and dimensionless
clock relation under its declared uniform binary reference and IID seam
schedule. These are readouts of the same recorded events, with a finite port
chart and attempt parameter whose physical identification requires an
independent correspondence.
The [repair-law guide](CANDIDATES.md) retains the exact baseline relations.

N2 lets the native current drive the existing local six-level response;
N3 is its simpler flat alternative. They preserve the full native load
instrument, so internal preparations cannot feed back on that marginal.
N4 instead conjugates binary C/E by a coherent load–response exchange.
It preserves locality and total occupation and supplies finite reciprocal
feedback, with an added coupling θ and specified internal preparation.
N4's record label j is not a bare-load current.

The [central-address extension](../code/native_response/README.md#a-central-port-carrier-that-retains-the-quantum-repair)
places the full quantum algebra in twelve identical blocks. Refreshing the
central address and repairing a uniformly selected incident seam reproduces
the entire original channel after decoding. With its declared collective
response, chart implementation, scalar seams and all-level refinement, this
gives a complete abstract A1–A3 source on the binary shell. Reversible
response words and chart paths have control sequences fixed independently
of measured records; feedback uses the declared repair/reset templates.
Allowing records to select reversible words would add 144 independent
address-controlled directions and change this source grammar. The minimum
twelve-block dimension assumes a unital embedding of the full working
matrix algebra and exactly twelve primitive central ports; it is not a
universal physical encoding bound. It does not extend N4 to all 78 loads.

Native and N4 versions share the full tracial quantum reference, collective
response and normalized occupation geometry while their finite kinetic
invariants differ. These requirements therefore allow both implementations;
uniqueness is unnecessary for constructing either. A central control address
is not an observed occupied port. For the prescribed occupied-port effects,
every instrument erases coherences between disjoint load supports, with
recovery error at least `1/2` on an explicit orthogonal pair. Endpoint
comparison consumes a separate record service; scalar seams carry no remote
record. The explicit occupied-port endpoint instrument recovers the
certified global-attempt statistic under the tracial preparation; adding
intermediate occupation reads changes that experiment. Physical realization and service costs remain additional contracts.

The [native-response calculation](../code/native_response/README.md) uses
one common tracial state and fixed bare occupations to recover the same
port Gram geometry and a certified two-attempt memory excess for N4.
For one local occupation under global IID attempts, the corrected statistic
is `b_local=c2−c1²−(1−c1)²/10`. The subtraction removes ordinary spatial-mode
mixing and gives a zero native null. All outcomes are retained and there is
no intermediate measurement or reset.

The [local-tick extension](local_tick/README.md) counts every attempted seam
incident to the chosen port, including unchanged loads. All intervening
remote attempts remain in the evolution. Under uniform IID seam selection
the incident probability is `p=1/6`; the first-hit Schrödinger channel is
`K=B(I−A)^−1` for remote/incident channel sums A and B. It requires a different
correction:

```text
r = 1−c1,
b_tick = c2−c1²−r²(209r+180)/(1331r²+3300r+1800).
```

Native repair has `c1=5/11`, `c2=93/407`, `b_tick=0`. Reusing the existing
edge-Gram certificate `G≥mI` proves `b_tick≥3m/50` for N4. Certified lower
bounds are `0.00001119970056` at θ=1/16 and `0.00022670580642` at θ=1/8.
These are lower bounds, not numerical evaluations of the full N4 correlation.
C and E share c1 and these bounds; the local-tick result also gives
`c2_C≥c2_E`, without asserting strictness.

Preparation, occupation readout, incident-event recognition and calibration
remain specified services. The N1 current detector does not implement this
N4 occupation apparatus. Local attempt count, model time, host runtime,
modular parameter and calibrated observer duration are distinct. The
[working terms](../docs/INTERPRETATION_AND_TERMS.md) impose no preferred
universal physical clock.

<a id="retained-scalar-and-entropy-constructions"></a>

## Scalar and entropy constructions

### Executed scalar detector and native continuation

The [scalar adapter](observer/scalar_response/README.md) prepares the original
finite-action vacuum with signed momentum mean `4j f`, where the first native
current is `j∈{−1,0,1}` and f is the inherited left-slab profile. The supplied
positive scalar Hamiltonian evolves this field. A binary Weyl instrument
reads the inherited detector profile at model time `t∈[19/20,1]`, retains
its conditional output state and commits `y=±1`. Only y chooses the second
native seam. The first current and the detector bit remain readable after
feedback. Conditional states use finite superpositions of Weyl-displaced
coherent states, with their coefficients and Gram overlaps retained; a full
Fock-space vector need not be allocated. The instrument acts on the field;
it is not a draw from a stored theoretical receipt or execution of a compiled
native quantum-control word.

Write `r_q(t)=Pr(y=+1|j=1)−1/2` at refinement q. Symmetric signed preparation
gives the native continuation statistic

```text
E[j2] = 5 r_q(t)/132.
```

The comparison uses the existing golden refinements `q=13,34,89,233`, the same
preparation and detector profiles, and the same model-time window. If the
parent scalar bound is `|r_q−r_cont|≤ε_q`, the complete labelled outcome-history
law differs by total variation at most `5ε_q/33`, and the second-current mean
differs by at most `5ε_q/132`. These transfer bounds concern this specified
observable and protocol; they do not assert convergence of the full Fock
state or an interacting theory. The first two grids retain the parent's
trivial error bound; q=89 and q=233 give nontrivial bounds. The analytic parent
establishes the limit independently of these four numerical executions.

The [local evolution](observer/scalar_response/local_evolution.py) executes
the canonical kick–drift–kick recurrence with the local tensor stencil at
q=13 and q=34. It advances the mean and detector pullback columns; the
original-vacuum covariance then supplies the squeezed detector variance.
This local receipt performs one field read and two native repairs.
The q=13 step-halving control compares at the same actual model time. The
separate [receipt](observer/scalar_response/local_results.json) retains
checkpoint summaries and sampled native records. It is not an authenticated
log of every field update or a compiled native quantum-control word.

The existing time-refinement certificate adds a temporal detector error to
the spatial one on all four grids, including the change in the original
vacuum under split evolution. The same `5/33` and `5/132` factors transfer
that combined bound to complete histories and the mean current. A nearest
grid-time comparison also includes the parent's timing bound. Only q=13
and q=34 execute the local loop here; their inherited spatial bound is
trivial. The tighter q=89 and q=233 spatial/temporal bounds remain analytic
comparisons for this extension. Numerical agreement with the spectral
oracle is a software check, not a certified floating-point error bound.

The [prediction protocol](observer/scalar_response/PREDICTION.md) fixes the
feedback statistic with disconnected, zero-preparation and sign-reversal
controls. Independent mathematical reconstruction and seeded execution test
the software prediction. A physical comparison additionally needs a realized
field, apparatus, calibration and independently acquired records. The scalar
action, quantization, vacuum, preparation, Born readout and model clock remain
explicit inputs; the finite native observer does not select them.

### Field and record continuation

The [continued scalar observer](observer/scalar_response/CONTINUATION.md)
retains one field through two noncommuting reads at model times `19/20` and
`1`, with three repairs of the same native loads. All 2,112 labelled histories
are kept. Its 209 reachable record prefixes have an explicit eight-bit
encoding; source, field and service memory are separate. Four coherent Weyl
terms retain the terminal quantum state, while every one of the 12,487,168
q=233 modes contributes to their Gram kernel.

This extension uses continuous full-mode spectral evolution. The local
split-update engine above and this two-read conditional-state continuation
are separate implementations with separate error bounds.
The two-site source compiler executes a separate single-read apparatus;
compiling this two-read history requires its own implementation and bounds.

The directed rational bound puts the q=233 complete-history TV error
below 0.02571. It includes the zero-current histories, whose two-bit
correlation is not controlled by the single-read bias bound. It also
certifies a positive second-read disturbance above 0.000022281162 relative
to omitting the first read on the positive-current preparation. A commuting
Gaussian control retains ordinary temporal correlation but removes that
disturbance. Independent state-vector, conditional-field donor, reset and
record-deletion controls distinguish these mechanisms. The existing supplied
action, vacuum, instrument and time contract is preserved.

### Scalar action and entropy premises

The separate scalar benchmark factors weighted agreement into its supplied
stiffness. If repair capacity also measures inertia, linear dynamics, energy
conservation, time reversal and the declared velocity reading give
`M f_ddot=−(K_D+m²M)f` and the canonical Hamiltonian. These hypotheses determine
the equation within that class. The `I+εA` kinetic counterfamily has the same
normalized classical configuration Gibbs marginal with different frequencies;
the full phase-space reference and quantum vacuum change. Three-mode curvature
and a distance-two response distinguish it after common clock/mass/speed
refits. Nonzero deformation violates the same-capacity inertia premise;
configuration equilibrium alone does not determine inertia.

Local coordinates `χ=(m q+i p,Dq)` represent the original canonical fields
without independent edge oscillators. The existing free-scalar continuum,
local force and retained-pointer instruments remain conditional results;
interactions and physical matter identification are further work. See the
[source instruments](../code/source_scalar_instruments/README.md) and
[finite implementation](../code/source_scalar_finite_instrument/README.md).

The retained q=3 eight-mode preparation/history uses the same 8,192-level
onsite code, with combined mathematical trace-distance error below `0.001237`
before physical service errors. Its ramp error is below `0.000433450` at
model duration 512; `2^24` local time steps add below `0.000695044`.
Existing source controls give a finite synthesis recipe, and exact commuting
compilation removes an unnecessary inner approximation loop. These are
construction/error results; no enormous native pulse word or 106-qubit
field/pointer state was executed. The [premise map](PREMISES.md#re-10)
retains the detailed preparation, history and resource assumptions.

Regional maximum-entropy completion, nested inclusion towers, protected
records and charge-resolved collars remain available from the
[finite algebra results](../code/quantum_information/ALGEBRAS_AND_REPAIR.md).
A charge center is not by itself a physical gauge-edge factor, area law or
stress coupling. Nonnested completions need not commute. The
[gravity premise ladder](../docs/registers/PREMISE_REGISTER_V3.md#gravity-rungs) retains
the physical assumptions needed for the existing gravitational conclusions.

## Connecting finite objects to the physical theorem routes

The same word can name different objects in these constructions. Use the
actual carrier, operation and accessible record when applying an existing
theorem. The [generated map](PREMISES.md) links their executable sources and
proofs; the [observation ledger](../docs/registers/OBSERVATION_LEDGER_V3.md)
records the broader physical targets without promoting their status from a
finite software run.

| Constructed object | Result that can be reused | Required connection for a broader claim |
| --- | --- | --- |
| Common flux-dressed CAR/scalar register and two recorded instruments | Source execution, operator Gauss, committed records, finite temporal and occupation bounds, exact sharp-readout contrasts, state-specific complete-output cutoff enclosures and unavoidable occupation-code output errors | Spatial refinement, vacuum/coupling identification, numerical/physical source services and responsive stress/geometry on this same process; post-sharp clock continuation needs its own error analysis |
| Observer representations and appended clocks on the finite sharp history | GNS and local-reference modular readouts, fixed full-work Gibbs/modular correspondence, source-compiled clock probes and CP regional restrictions | Declared observer access, factorization and references; reference selection by repair, physical clock calibration, growing spatial family and occupation limits of the appended clocks |
| Separate smooth Ramsey detector on the common action | Quantitative complete-output occupation bound, finite-energy backreaction and a separate prepared four-read history certificate | Different instrument and bit dictionary; the uniform moment bound is vacuous at small executed cutoffs, and the sharp-history clock does not transfer |
| N1 source, phase and finite record tape | All-input CPTP law, repeatable committed classical records, exact forecast and feedback controls | Extension of the specified finite observer/federation contracts to the intended access, delivery and refinement family; no full physical `OperationalObserver` realization is asserted |
| P-dependent path-qubit detector | Complementary binary effects, branch charge covariance and bit-controlled native feedback under accepted closure inputs | Physical realization of the declared apparatus; a same-quantity screen/Thomson readback and precision endpoint claim keep their existing contracts |
| Native-current scalar preparation, binary field instrument and native continuation | Executed conditional instrument, feedback mean `5r_q/132`, inherited history and mean refinement bounds | Source identification of the scalar action, vacuum, preparation/readout and clock; a physically calibrated independent detector comparison |
| Reachable native record words on the declared carrier | Existing correctable-record evaluator and explicit classical prefix encoding | A compatible physical capacity cut and global closure family before identifying this count with cosmological N |
| Twelve central addresses with binary native/N4 quantum blocks | Full-channel decoding, complete abstract scalar-seam grammar, shared full reference and normalized geometry with distinct finite kinetics | Physical realization, retained-record delivery and compatible apparatus; the address is not an occupation measurement |
| N4 tracial bare occupation at incident attempts | Corrected first-hit null and inherited strict finite memory bounds | N4 occupation preparation, detector and incident-event recognition; the N1 current detector does not measure this occupation |
| Supplied canonical scalar action, vacuum and Weyl observables | Existing finite instrument and controlled free-field continuum theorems | Explicit maps from the chosen native implementation with controlled state, observable and service errors along the same refinement family |
| Native C/E endpoint maps and their federation extensions | Local Liouville action, exact instrument blocks and complete quantum observation locality | Compatible screen refinement, physical clock/field correspondence and controlled relativistic scaling |
| Joint local clock and field instruments | Exact record-conditioned quantum memory; completed-event work/reference measure-and-prepare closure | Declared support, schedule, preparation and joint apparatus; unresolved waits must remain in the quantum history |
| Coherent overlapping observers on the native total-two carrier | Exact complete-output correspondence for a transported matrix cone, including arbitrary untouched references and number-preserving union-supported quantum interventions | Declared coherent preparation and apparatus; the inverse encoding is not positive on arbitrary effective matrices, and the correspondence supplies no relativistic field dynamics |
| Seven-to-nine-mode interacting scalar extension | Matched low-mode action and a complete reduced-output error bound through continuing observation and feedback | Supplied scalar action, oscillator blanks, finite occupation codes and model-time intervals; further bands, physical parameters and Lorentzian observation laws need additional control |
| Continuing native preparation and likelihood test | Uniform finite preparation and dependent-transcript decision bounds | Physical apparatus mapping and authenticated conditional service errors; the finite checks are not an executed laboratory campaign |
| Declared gauge, edge, stress and geometry structures | Existing conditional gauge/gravity implications | Joint physical carrier, interactions, conserved stress and compatible geometric variations on that family |

No observation rung or physical premise is discharged merely by putting these
sources in the same repository. Conversely, the remaining physical joins do
not reopen the completed finite observer contract. Model maintenance belongs
here; changes to the broader papers follow their separate review policy.

<a id="retained-audit-controls"></a>

## Counterexamples and identifiability

The [fixed audit package](../code/programme_obstruction_audit/README.md)
retains exact author and independent controls for the following boundaries.
Run `python3 code/programme_obstruction_audit/run_all.py` from the repository
root with assertions enabled. Its seven scripts validate the frozen inputs
and write fresh receipts outside the package. These controls reject the
named identifications or overly broad inferences; none rejects A1–A3 or
requires a unique designed implementation.

| Audit issue | Exact retained boundary | Construction consequence |
| --- | --- | --- |
| [#1087: complete carrier](https://github.com/FloatingPragma/observer-patch-holography/issues/1087) | `Z(M66)=Z(M78)=C I` cannot supply twelve nonzero primitive central ports. Binary occupations are noncentral, overlap and sum to `2I`; total-two occupations also admit eigenvalue two. | Keep the native matrix algebra as a sector. The twelve-block central-address extension supplies the missing central algebra and real port-controlled operations; occupation measurement remains a separate contract. |
| [#1088: marked triplets](https://github.com/FloatingPragma/observer-patch-holography/issues/1088) | The slow Gram satisfies `G²=4G`, `tr G=12`, `(5I−A)G=(5−sqrt(5))G`. One fixed order-five action has slow-frame trace `(1−sqrt(5))/2` and lower response-triplet trace `(1+sqrt(5))/2`, excluding an invertible same-action intertwiner. | Keep spatial and internal response domains distinct. Whole-law Galois substitution or an explicit action twist needs its own attachment; incidence recharting alone does not identify the marked triplets. |
| [#1089: feedback class](https://github.com/FloatingPragma/observer-patch-holography/issues/1089) | At fixed load `I_odd/2`, the rational number-preserving exchange gives endpoint probabilities `769/1250` and `481/1250`, difference `144/625`. Four Kraus operators certify CPTP; its private commutant has dimension one, excluding every single-product exponential representation, including alternate logarithm branches. | Preserve the all-input reduced-load and single-product no-feedback theorems within their stated classes. The escape is a finite channel control; importing geometry, current or clock semantics requires a common observable map. |
| [#1090: configuration equilibrium](https://github.com/FloatingPragma/observer-patch-holography/issues/1090) | For `A=I+L_path3`, the positive `ε=±1/16` controls retain the normalized configuration Gibbs marginal but change the full reference and vacuum. For stiffness response `Y=γ(A+εA²)`, `Y13=γε` and the third-mode fit residual is `6γε`; actual acceleration is `−Yq`. | The same-capacity inertia theorem retains its premise. Separate this configuration-only control from the native/N4 counterpair, which fixes the full quantum reference and normalized geometry but changes finite kinetics. Neither supplies two physical realizations. |

The [independent joint reconstruction](../code/programme_obstruction_audit/joint_independent.py)
checks the marked path and occupation identities; the
[independent channel reconstruction](../code/programme_obstruction_audit/channel_independent.py)
checks all 64 operator matrix units and solves the private commutant.
The [kinetic control](../code/programme_obstruction_audit/kinetic_control.py)
retains all eight prescribed rational cases, including rejection of singular
`ε=−1/4`. The source extension does not erase these negative controls.

## Remaining construction work

The missing joins concern one compatible family of observer observations,
not an inventory of additional finite theorems:

- Establish compatible federation neighborhoods, states, accessible algebras
  and instruments across the intended scaling family. The conditional
  quasilocal quantum construction has these hypotheses; the existing screen
  refinement does not discharge them.
- Connect the derived native local action and event-clock histories to an
  effective interacting field/spacetime description through complete joint
  observation and intervention laws. Preserve quantum memory where a reduced
  field is not Markovian. Modular parameters and event records need a checked
  correspondence on that process, not a preferred universal clock.
- Identify the physical sector, apparatus and parameters with independently
  calibrated measurements. Specify data ancestry, permitted fits and units,
  then lock predictions and rejection rules before comparison. Source control
  words implement a supplied action without selecting its physical meaning.
- Realize the one-event interference apparatus that reads the record law's
  coherence `g`, authenticate its full service budget and collect complete
  records. Finite native preparation and a dependence bound replace exact
  reset and IID premises; physical controls and a laboratory mapping are
  open.
- Connect interacting matter, conserved stress and responsive geometry to the
  same observations and observer-relative clock comparisons. Fixed-volume
  occupation control and growing volume at fixed spacing address different
  limits; neither alone removes a spatial regulator.
- Account for long-running record storage, coherent services and fault
  budgets while preserving finite protocol behavior and references. Checked
  simulator restoration is not a physical reset or autonomous memory service.
  N4 additionally requires its stated occupation preparation/readout and
  event-recognition apparatus.

Useful tests target actual alternatives. The
[binary coherent continuation](observer/COHERENT_CONTINUATION.md) implements
a reversible 3,960-arc lift and a compressed-amplitude recurrence. Different
continuations or cross-record phases give distinct Ramsey predictions while
agreeing on the committed native instrument. A selecting restriction must
constrain that additional interface freedom. The
[identifiability checks](identifiability_checks/README.md) likewise separate
apparatus restrictions, unseen action directions and geometry ambiguities.
The stronger effective-Hamiltonian comparator remains relevant to empirical
protocol design; missing timing, normalization or raw counts cannot be
inferred by fitting away residuals.

The [mathematical reuse guide](../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md)
requires compatibility defined independently of the update and a mapping of
the actual domain, records and scheduler to the nearest theorem. Complete
history, instantaneous fixed state, stationary distribution and continuing
operation are distinct claims. A self-consistent history can contain changing
records and interventions.

## Source selection: the three primitive obligations

The [source-selection guide](SOURCE_SELECTION.md) organizes source selection
around carrier architecture, observer interface and implementation, and the
complete repair rule. For the repair rule the selection problem is sharp:
the classical seam rule is forced by four stated requirements, the record law
ranges over an exactly classified convex family, and one interference
experiment reads the coherence that separates its members. The
[repair-law guide](CANDIDATES.md) states the principle that singles out each
distinguished member. Architecture, observer implementation and repair law
may supply the instance inputs. Action, state, scheduling and control consequences require
explicit derivations from those inputs or a separate declared status; renaming a
supplied input does not discharge it. The existing finite field programs
retain their comparison scope.

Selection is up to the observations admitted by the interface. Scalar phases
on committed Kraus branches cancel through every finite adaptive recorded
history. They become observable when uncommitted record interference is
admitted. The existing source words execute a fixed read with `3/4` versus
`1/4` plus probability for two phase implementations of the same committed
C/E instrument; the commit-first control gives `1/2` for both. The distinct
record dictionary is explicit. Coherent continuation alternatives also
survive all proper carrier rotations, real orthogonality, native support and
a common algebraic time reversal. A source-selection criterion must exclude an
actual alternative under the relevant interface, rather than infer a
preferred wave from these shared structural properties.

## Work ownership and claim accounting

`model/` is the live home of observer and repair-law construction. Its
implementation, tests, candidate comparison and this guide evolve together.
GitHub issues own tasks; the [claim registry](../claims/claim_registry.yaml)
owns established claim scope; [PREMISES.md](PREMISES.md) is generated
navigation, not a proof that every theorem antecedent has been discharged.

The connected responsibilities are [common construction #1026](https://github.com/FloatingPragma/observer-patch-holography/issues/1026),
[shared repair law #1076](https://github.com/FloatingPragma/observer-patch-holography/issues/1076),
[complete source history #1077](https://github.com/FloatingPragma/observer-patch-holography/issues/1077),
[physical tests #1025](https://github.com/FloatingPragma/observer-patch-holography/issues/1025),
[correctness #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
and [obstructions #1078](https://github.com/FloatingPragma/observer-patch-holography/issues/1078).
They are parallel responsibilities, not twelve sequential completion hurdles.
Finite results retain their stated validity unless evidence invalidates
them; larger claims add declared requirements. A negative result rejects its stated family.

## Research programme and paper-release hold

Publish tested model code and reproducible results in this directory. The
[paper set](../docs/PAPER_SET.md) records the model milestone around which
the prose papers are organized and the role of each paper. **A paper
release follows the owner's approval of a recorded model milestone.**
Audits, necessary correctness fixes, errata, preview builds and
reproducibility checks continue between releases. Frozen evidence and
deposited releases retain their original bytes. Passing CI or publishing the
model code does not itself release a paper.
