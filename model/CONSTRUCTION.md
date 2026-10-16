# OPH model construction from local repair

This is the construction guide for the [live model implementation](README.md).
The [candidate register](CANDIDATES.md) gives law-level comparisons; the
generated [premise map](PREMISES.md) connects those constructions to the
existing theorems.

**Build the simplest shared repair process that satisfies the required
constraints.** A finite implementation is complete when its declared
observation, record and feedback contract runs and passes its checks.
Proving that self-consistency uniquely selects that implementation is a
separate question, not a prerequisite for constructing it. A complete
physical model additionally needs compatible physical identifications and
tests on the same process.

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

Hardware, protocol and implementation choices remain open design variables.
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
of the observer statistics, and independent empirical completion, remain open.

## Current executable result

The [primitive observer runtime](primitive_observer/README.md) executes
complete instruments on mixed states, committed records, feedback and
untouched references. Source-derived native dynamics and supplied physical
benchmark programs use this same execution interface. Their preparations,
clocks, interpretation maps and error budgets remain explicit.

### Local construction and its extension

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
nontrivial physical scaling regime remain open.

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
stationary boundary. Its exact degree-two exclusion dual avoids enumerating
the full binary configuration space without factorizing correlations.
Feedback reverses the sign of the certified connected four-field statistic
on both actual screens. Two field counters and the controller memory survive
serialized source continuation; arithmetic, cutoff and spatial errors are
accounted separately. This completes their finite joint time/field domain.

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
contracts. Mean loads and native classical histories alone cannot select C
from E. No public laboratory dataset is identified with these native-law
experiments.

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
The record is `(c,j)`. Keeping the pre-repair contrast instead gives the
alternative instrument E, with identical classical transitions but different
quantum coherences. The classical repair rule alone does not select C or E.

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
The [candidate register](CANDIDATES.md) retains the exact baseline relations.

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
- Realize the native C/E discrimination apparatus, authenticate its full
  service budget and collect complete records. Finite native preparation and
  a dependence bound replace exact reset and IID premises; physical
  controls and a laboratory mapping remain open.
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
complete repair rule. Architecture, observer implementation and repair law
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
