# OPH model construction from local repair

This is the construction guide for the [live model implementation](README.md).
It replaces the former model-search ledger with the current contracts,
results and next connections. The [candidate register](CANDIDATES.md) gives
law-level comparisons; the generated [premise map](PREMISES.md) connects
those constructions to the existing theorems.

**Build the simplest shared repair process that satisfies the required
constraints.** A finite implementation is complete when its declared
observation, record and feedback contract runs and passes its checks.
Proving that self-consistency uniquely selects that implementation is a
separate question, not a prerequisite for constructing it. A complete
physical model additionally needs compatible physical identifications and
tests on the same process.

## Current executable result

The [native observer](observer/README.md) executes two repair events with a
finite tape. It records the first current, issues a supplied forecast from
that record, chooses the next seam, records the result and halts. Exhaustive
checks retain all 264 binary-source paths. A 20,000-trial seeded execution
reproduces the exact predictions within the reported sample differences.

The [P-dependent detector](observer/p_response/README.md) computes the
existing P-to-electromagnetic-coupling map, writes a binary detector result
and lets that bit choose the second repair. All 528 weighted histories are
retained. A separate 200,000-trial execution agrees with the analytic
response within its reported sampling uncertainty. Only the detector bit
reaches the actuator; feedback preserves its committed value.

These are executable finite observation–record–feedback constructions.
The initial ensemble, tape, forecast table, detector and control protocol
are explicit design inputs. Two events means two successive repair
invocations; it does not mean that the source or complete history has only
two alternatives. The binary detector is an access interface to a richer
source and record history.

The [scalar-response construction](observer/scalar_response/README.md)
extends this protocol through the existing scalar field and its bounded
continuum comparison. The native current sets a local preparation sign;
the field detector writes the bit that controls the next native event.
This joins an executed observer to a supplied field model without identifying
native integer averaging with the field Hamiltonian.
Its local evolution extension advances that field by neighboring-site
canonical updates at q=13 and q=34, propagates the detector's vacuum variance,
and uses the resulting probability in the same native continuation. It
tests temporal refinement while retaining the existing continuum theorem.

The [source-compiled two-site construction](observer/source_scalar/README.md)
emits and executes the field apparatus through the existing twelve-port
controls. Blank registers are transformed into an entangled projected vacuum;
the first native current supplies an impulse; finite quantum evolution and
a quantum pointer produce both recorded outcome branches. Only the pointer
bit controls the next native repair. The four- and eight-level site codes
run at 64 and 128 steps with exact initial preparation-tail and finite
splitting certificates, independent branch-density checks and retained
200,000-trial feedback runs per code size. The finest protocol uses seven
work qubits, 1,774,010 native pulses and 3,548,020 declared proper-code
service events, evaluated by caching executed subroutine contractions.

This is a separate two-site action, not another member of the earlier
spatial refinement family. Its physical occupation projection retains the
boundary terms in both quadratic operators; squaring truncated coordinates
would implement a different Hamiltonian. The source controls implement
preparation, evolution and readout under supplied coherent transport,
reset and storage services. Full infinite-occupation propagation/readout
bounds, spatial refinement, physical service errors and empirical
identification require separate bounds and evidence. The earlier continuum bounds cannot be
transferred merely because both experiments use a scalar field.

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

The observer first applies C on `(0,1)`. Its committed first record alone
selects `(0,5)` when `j1>0`, otherwise `(0,1)`, and supplies the preissued
forecast of a positive second current. No source loads or future outcomes
are accessible to that controller. Full simulator checkpoints and external
occupation audits are distinct from the observer's record access.

The [derivation](observer/DERIVATION.md) gives one fixed CPTP map on the
phase/tape/source carrier and shows classical-prefix preservation. Genuine
feedback can reduce coherence between different old record values while
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
measurement. P does not alone determine a whole observer's Hilbert dimension,
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
schedule. These are readouts of the same recorded events, without identifying
a finite port chart with physical space or attempt count with proper time.
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

## Retained scalar and entropy constructions

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

### Retained field and record continuation

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

The new directed rational bound puts the q=233 complete-history TV error
below 0.02571. It includes the zero-current histories, whose two-bit
correlation is not controlled by the old single-read bias bound. It also
certifies a positive second-read disturbance above 0.000022281162 relative
to omitting the first read on the positive-current preparation. A commuting
Gaussian control retains ordinary temporal correlation but removes that
disturbance. Independent state-vector, conditional-field donor, reset and
record-deletion controls distinguish these mechanisms. The existing supplied
action, vacuum, instrument and time contract is preserved.

### Scalar action and entropy results retained as inputs

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
[gravity premise ladder](../docs/registers/GRAVITY_PREMISE_LADDER.md) retains
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
| N1 source, phase and finite record tape | All-input CPTP law, repeatable committed classical records, exact forecast and feedback controls | Composition into the intended observer menu, delivery/refinement and continuation contract; no full `OperationalObserver` or federation instantiation is asserted |
| P-dependent path-qubit detector | Complementary binary effects, branch charge covariance and bit-controlled native feedback under accepted closure inputs | Physical realization of the declared apparatus; a same-quantity screen/Thomson readback and precision endpoint claim keep their existing contracts |
| Native-current scalar preparation, binary field instrument and native continuation | Executed conditional instrument, feedback mean `5r_q/132`, inherited history and mean refinement bounds | Source identification of the scalar action, vacuum, preparation/readout and clock; a physically calibrated independent detector comparison |
| Reachable native record words on the declared carrier | Existing correctable-record evaluator and explicit classical prefix encoding | A compatible physical capacity cut and global closure family before identifying this count with cosmological N |
| Twelve central addresses with binary native/N4 quantum blocks | Full-channel decoding, complete abstract scalar-seam grammar, shared full reference and normalized geometry with distinct finite kinetics | Physical realization, retained-record delivery and compatible apparatus; the address is not an occupation measurement |
| N4 tracial bare occupation at incident attempts | Corrected first-hit null and inherited strict finite memory bounds | N4 occupation preparation, detector and incident-event recognition; the N1 current detector does not measure this occupation |
| Supplied canonical scalar action, vacuum and Weyl observables | Existing finite instrument and controlled free-field continuum theorems | Explicit maps from the chosen native implementation with controlled state, observable and service errors along the same refinement family |
| Declared gauge, edge, stress and geometry structures | Existing conditional gauge/gravity implications | Joint physical carrier, interactions, conserved stress and compatible geometric variations on that family |

No observation rung or physical premise is discharged merely by putting these
sources in the same repository. Conversely, the remaining physical joins do
not reopen the completed finite observer contract. Model maintenance belongs
here; changes to the broader papers follow their separate review policy.

## Retained audit controls

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

The finite observer is a completed component to extend and preserve. The
next concrete target is to implement this scalar preparation, local update
and binary detector on the chosen microscopic source, retaining their
records and bounding their joint error on the same refinement family.
That tests the missing source-to-field connection directly. A physical
completion additionally needs the compatible interacting and geometric
interfaces and an independent comparison below:

- Implement longer or composed observers with explicit memory, reset,
  continuation and resource contracts while preserving current behavior.
- Join the N4 preparation, occupation detector and local event counter if
  claiming its certified memory test as an implemented apparatus.
- Attach interacting matter, responsive geometry and physical clock
  calibration to the same state, event law and retained records.
- Specify a physically realizable preparation and observable, compare
  alternatives, and freeze a prediction and rejection rule before independent
  comparison data. Existing seen-data comparisons remain postdictions.

The [established-mathematics guide](../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md)
governs reuse: define compatibility independently of the update, map the
actual domain and scheduler to the nearest theorem, and identify what is new.
Complete history, instantaneous fixed state, stationary distribution and
continuing operation remain different claims. A history-level self-consistent
structure can contain changing records and interventions.

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
Finite results stay completed unless new evidence invalidates them; larger
claims add declared requirements. A negative result rejects its stated family.

## Research programme and paper-release hold

Publish tested model code and reproducible results in this directory.
**Paper releases and broad restructuring of the prose papers remain on hold**
pending a recorded review of a coherent model milestone. Audits, necessary
correctness fixes, errata, preview builds and reproducibility checks continue.
Frozen evidence and deposited releases retain their original bytes. Passing
CI or publishing the model code does not itself release a new paper.
