# OPH model search: iterative construction from local repair

**Governing principle: find the simplest shared repair rule that satisfies
the known constraints.**

**Programme completion requires finding that shared repair rule and showing
that it satisfies the required constraints on one common model.** Component
proofs, simulations, candidate implementations and issue closures are bounded
milestones; they do not by themselves complete the programme.

**Current status.** [N4's finite construction](../code/native_response/README.md)
supplies bare-occupation geometry and a certified two-attempt memory signal.
The corrected single-port statistic removes the native spatial-mode-mixing
term under the stated tracial preparation and IID global-attempt protocol. It
distinguishes N4 from the matched native heat description at two couplings; it
does not select C versus E or uniquely identify N4. The remaining model-search
gap is a single-source realization of preparation, local readout, record
delivery and clock comparison, together with a separate null for
observer-local stopping times.

The [candidate register](MODEL_CANDIDATES.md) compares the concrete active
laws, their tested properties, added assumptions and next decisive tests.

The [agreement research foundations](AGREEMENT_RESEARCH_FOUNDATIONS.md)
govern reuse of distributed self-stabilization, sheaf diffusion and extension
results. The search must map their premises to the candidate before adding
mechanisms or proving another abstract convergence result. Recovery, exact
agreement, boundary uniqueness and physical law selection remain separate
obligations; the existing source and scalar evidence is retained.

The [dual description](AGREEMENT_RESEARCH_FOUNDATIONS.md#one-structure-complete-history-and-internal-computation)
concerns one structure: a complete consistent history whose internal local
dependencies constitute its computation. Construct the same-history record
and intervention maps; distinguish solver scheduling from observer-local time.
This interpretation needs no external host or global clock and makes no new
endpoint-convergence claim.

**Programme organization:** research is organized around what one repair
process must supply. The [work ownership](#work-ownership-and-claim-accounting)
below assigns construction, physical tests and audits in parallel. Paper
releases and broad prose restructuring are subject to the
[model milestone](#research-programme-and-paper-release-hold).

**Current microscopic tests:** retain N2's local current-driven response and
the simpler flat N3 as preservation benchmarks. Test N4, a coherent exchange
using the same private factors and response matrices, for internal-state
influence on load repair. N4 supplies this finite feedback while inheriting
completion, at the cost of changing the load marginal and the meaning of its
current labels. A shared bare-occupation protocol recovers the same
finite port geometry and distinguishes its kinetics from a common constant
rescaling of native attempt time. The
[corrected single-port test](#a-single-port-test-that-removes-spatial-mode-mixing)
needs only local endpoint outcomes under the specified tracial preparation
and one-/two-attempt protocol. Full source grammar and physical selection require
justification. The
[general search findings](#general-findings-that-guide-the-search) explain
which tests discriminate laws and how added assumptions are counted.

**Current scalar benchmark:** the separate minimal canonical inertial
completion supplies controlled preparation, free evolution and retained
observations, with a same-code mathematical error below 0.001237 before
physical service errors. It is a component benchmark, not a selection of the
microscopic law. Its construction and compiler results remain available below.

**Model-selection priority:** one microscopic repair process must supply
compatible geometry, reversible response and clock readouts. The scalar
benchmark tests one component of that requirement. Simplifying its simulation
or compiler does not select the full source law. Use the existing theorem
interfaces below to constrain the next candidate before adding mechanisms.

## 1. Headline results and missing-link completion status

- **The reciprocal candidate keeps port geometry but changes its kinetics.**
  At the supplied tracial reference, N4's bare-occupation Dirichlet matrix is
  exactly `α(θ)L_ico/396`, with `α≥1/9` for every real coupling. The existing
  rank-three port Gram therefore survives. Two endpoint measurements, with
  no intermediate dephasing, also expose a strictly positive two-attempt
  memory excess at θ=1/8 and 1/16. This rejects an autonomous native heat
  equation with a fitted constant attempt rate on the same bare observables. The
  [derivation and finite certificate](#bare-occupation-geometry-and-a-clock-independent-memory-test)
  use the same N4 states and events; they add no interaction or state factor.
- **One coherent exchange supplies finite load feedback.** N4 uses the
  existing binary repair, response matrices, local factors and one supplied
  coupling. Changing only the private preparation changes the endpoint
  occupation probability by approximately 0.035959 at θ=1/8. Both simpler
  single-product controls have no feedback, as predicted analytically.
  All 187 fixed checks and an identical isolated replay passed. Locality,
  conservation and completion transfer;
  original current-label semantics and geometry/clock identities do not.
  This is a candidate interaction on the binary shell, not a physical force
  derivation or an extension to all integer-load states.
- **Candidate discrimination needs fewer resources than the retained-history test.**
  N2 and its flat alternative N3 differ after one seam attempt and one local
  six-level preparation/effect, averaging both coin outcomes. At supplied
  θ=1/8 the probability difference is approximately 0.049095. All 67 fixed
  checks and an identical isolated replay passed. No entangled reference,
  coherent history or postselection is needed; physical preparation and
  readout are separate inputs.
- **One native event can join repair and local reversible response.**
  Attach the existing six-level response to the recorded current at both seam
  endpoints. The original load channel, classical records and disjoint-seam
  locality are preserved. A unitary change of coordinates reuses the existing
  completion proof; positive seam mixtures converge on the finite binary
  shell while retaining protected local block populations. All 268 checks and
  an identical isolated replay passed. This is a one-way finite coupling,
  with extra local response factors and a supplied coupling; no physical
  clock, full A1 grammar or reciprocal force is established.
- **The native repair connects port geometry and a dimensionless
  clock relation.** In the declared total-two nearest-integer model at uniform
  binary equilibrium, the same
  recorded events reconstruct `L_ico=3B/ν`, its rank-three slow Gram
  `G=4P_slow`, and `κ γ_s/ν_t=(5−sqrt(5))/6`. Independent pile interventions
  recover the same mean response; state-dependent activity controls delimit
  its interpretation. Preserve this established common-process result as a
  primary law-level benchmark. Its event-to-proper-time map and reversible
  gauge-current attachment are separate missing connections.
- **One conditional source construction joins preparation and retained
  observations on the same code.** The original eight-site q=3 blank,
  coupled-vacuum preparation, active free field, force and two coherent
  retained pointers use 8,192 onsite levels throughout. Existing source
  controls give a finite synthesis recipe. The combined mathematical
  trace-distance error is **<0.001237**, including preparation, time steps,
  occupation cutoffs and both synthesis stages. Complete readout preserves
  the bound without postselection. Four reviewed packets passed **535
  checks and four identical isolated replays**; these include analytic-bound,
  small-matrix and bookkeeping checks, not a simulation of the 106-qubit
  field/pointer state. Physical initialization and service errors are additional.
- **A simpler exact compiler removes an unnecessary approximation loop.**
  Diagonalize the chosen finite quadrature and compile that basis with existing
  source controls; the remaining diagonal phases commute. All 45 fixed cases
  at d=2,4,8 passed full-phase native-word comparisons and an identical replay,
  with largest numerical operator residual below 2.2×10^-13. This removes
  the inner Strang repetition under ideal controls without changing the law.
  The older full-history ceiling <3.44×10^34 pulses is a conservative baseline,
  not an intrinsic cost. Large-cutoff precision and complete source execution
  retain their own obligations; RE-10 is conditional and #1077 stays open.
- **The scalar benchmark can be evaluated directly at modest classical cost.**
  One onsite Gaussian calculation prepares the original eight-mode state and
  follows all five retained-pointer histories. With fixed steps chosen before
  execution, every prepared joint-history distance to the continuous reference
  was numerically below 0.000493. All 90 checks and the identical isolated
  replay passed. This is a practical simulator of the supplied quadratic
  candidate; its numerical errors do not replace the rigorous source-code
  budget or establish native execution.
- **The repair force has a minimal conditional Hamiltonian completion.**
  Weighted agreement factors the supplied scalar stiffness exactly. If the
  same capacity measures inertia, linear autonomous dynamics, the declared
  velocity reading, energy conservation and time reversal uniquely give
  `M f_ddot=−(K_D+m²M)f`. The canonical Hamiltonian and action follow.
  Choosing that inertia and those premises remains substantive.
- **Local field access is compatible with that completion.** The constrained
  variables `χ=(m q+i p,Dq)` obey a local self-adjoint first-order equation
  with `DᵀD=A−m²I`. Recover q and p locally from its vertex quadratures.
  No independent edge oscillator or nonlocal positive-branch preparation is
  needed. χ is a field-operator vector, not a one-particle wavefunction; its
  commutators and quantum ordering are explicitly retained.
- **The existing scalar continuum and quantum instruments join exactly.**
  The original positive canonical Hamiltonian, fields, vacuum, local force
  and pointer couplings are unchanged. Their massive three-dimensional
  free-field continuum and conditional operational results apply directly.
  This route needs neither a doubled sheet nor the arc walk's phase-error
  correction or population-dependent event normalization. Those are valid
  alternative constructions, not requirements of this canonical completion.
- **Finite source controls implement the scalar operation components.**
  The existing twelve-port response and coherent code service compile local
  force and pointer factors into proper qubit codes. Direct infinite-oscillator
  bounds include escaped amplitude and arbitrary reference entanglement.
  All 18 fixed cases passed. No extra primitive oscillator gate is assumed;
  complete source execution and physical service budgets remain separate.
- **A local schedule prepares the original coupled vacuum.** Starting from
  independent frequency-one oscillator blanks, smoothly turn on the existing
  stiffness. The path stays gapped and has the original eight-mode q=3 vacuum
  as its target. A sharper endpoint argument certifies full-state error
  **<0.000433450 at the same duration 512**, improving the previous
  <0.070155 bound by more than 160 times without changing the law. The fixed
  duration 1024 gives <0.000108323; 37 new checks and an identical isolated
  replay passed. The numerical distance at 512 is about 0.00016024 and stays
  diagnostic. Blanks, scheduled Hamiltonian control and model time are supplied.
- **The same preparation has an explicit local time-step construction.**
  Midpoint kick–drift–kick includes the changing-stiffness error and has a
  full Gaussian state bound of order N^-2. At T=512 and N=2^24 the bound is
  <0.000695044 relative to the continuous ramp. Together with its endpoint
  bound, this gives **<0.001129 to the original coupled vacuum**. All 62
  checks and an identical isolated replay passed. This is an oscillator
  construction; native pulse execution and service costs remain open. Computing its small symplectic matrices is not executing a native word.
- **The new interfaces passed bounded tests and countercontrols.** Exact
  source/action identities and full-matrix checks passed at 8 and 64 sites,
  including all frequency signs, local constraints, quantum commutators,
  preparation and detector comparison. The counterfamily with kinetic
  mobility `I+εA` has the same classical configuration equilibrium but
  different frequencies: equilibrium alone does not select the law.
- **The `I+εA` kinetic alternatives have a calibration-resistant discriminator.**
  Three squared frequencies detect a curvature that no common clock, mass
  and speed refit removes. An exact distance-two acceleration coefficient
  independently proves the distinction. All 16 fixed cases at 8 and 64 sites
  passed 258 checks, including unstable and affine controls. Vanishing kinetic
  perturbations can nevertheless share the position-Weyl continuum limit:
  that limit alone cannot select a unique microscopic inertia.
- **The remaining physical choices are sharper, not erased.** Justify the
  same-capacity inertia, source coefficients, quantum representation/state,
  controls and clock on one microscopic process; then attach interactions,
  physical matter and responsive geometry. The finite native-control,
  periodic continuum, entropy, geometry and protected-memory results remain
  available below with their original premises. No complete interacting QFT
  or empirical identification with our Universe is claimed.

**Status convention:** “finite complete” closes the stated mathematical
construction, under its declared inputs. It does not close the entire physical
entry. A construction, selection of that construction by OPH, and agreement
with observations are separate claims. The experiments retain analytical
proofs, exact controls, separate agent reviews and recorded replays locally;
this document is not an empirical verification bundle or an axiom amendment.

| Link | Completion status and constructed result | What completes the remaining physical connection |
| --- | --- | --- |
| **M01: common source** | **Finite source constructions and canonical field join complete.** Native C/E and coherent extensions; exact weighted factorization and local canonical scalar fields; N2's current-driven response preserves the load marginal, while N4's coherent exchange supplies private-state influence on bare occupations. | Select the source coefficients, inertia, quantum representation and controls. Conservative scalar evolution and dissipative repair are distinct laws. The added response factors and coupling need a complete source grammar; N4's shared bare-occupation protocol needs source-realized preparation/readout and a complete record/service interpretation. Finite feedback does not establish a physical force. |
| **M02: observers and agreement** | **Finite observer contracts, candidate equalizers and nested completions complete.** Proper overlaps, protected marginals and actual interventions. E completes the supplied exchange/phase-invariant algebra. Exterior repair realizes regional maximum-entropy completions with an exact inclusion tower. The current/response extension inherits a finite binary-network fixed-algebra limit. | Justify those interpretations from the required observer menu; translated-chart naturality alone does not select them. Supply the used completeness, acceptance, confluence, authentication and broader history contracts. Nonnested completions can fail to commute; the supplied targets are not the complete A2 diagram. |
| **M03: geometry** | **Finite geometry and controlled coarse energy constructed.** Native rank-three response; uniform memory/gap bounds and local coarse energies; periodic acoustic tensor; existing golden coordinates and uniquely reconstructible relative action capacities. | Identify physical regions, spacetime metric, area/volume and bulk refinement; select the response being reconstructed. Event-count weights and nonuniform action capacities are distinct finite references. A positive Schur energy need not be a Markov generator. |
| **M04: physical time** | **Finite clocks and explicit model-time laws constructed; physical attachment open.** Canonical flow reuses the original scalar time. The weighted walk and periodic exchange/wave scalings remain alternatives. | Justify and calibrate the physical law and clock with complete execution costs. Avoiding the walk normalization does not identify canonical time with accepted source-event counts or establish an operational physical cone. |
| **M05: fields and current** | **Finite currents, full scalar stiffness and local canonical fields joined.** Exact energy/symplectic preservation and local position/momentum decoder; no independent edge field is added. | Select physical charge, energy, stress and units; attach interactions and reciprocal geometry response. Local field coordinates do not supply stress coupling or native many-body closure. |
| **M06: quantum dynamics and regime** | **Conditional scalar continuum and one finite source construction joined.** The unchanged Hamiltonian and instruments share the same q=3 onsite code. Preparation, active evolution, force and two retained pointers have a combined mathematical error <0.001237, including time, cutoff and source synthesis bounds; the same complete readout contracts it. | Emit and execute an efficient native word; supply physical initialization, drift, pulse/service precision, storage and readback, and include their errors and resources. Justify representation, controls and physical time. The earlier k=3 local-code vacuum bound remains inconclusive. This fixed construction is not autonomous state selection, a refinement-uniform resource bound or an interacting QFT. |
| **M07: entropy/area stationarity** | **Finite counting stationarity, exterior completion and bulk collar splits complete.** `S_A+〈Z_B〉=log66−D(ρ_A‖τ_A)`. Exact recovery and a nontrivial 3×4 bulk block are constructed with declared charge centers. | Identify the formal split with the required physical edge factors and derive area, stress/flow, fixed-volume variations, normalization and common refinement. Z_B contains remote bulk entropy. Charge retention alone does not give all-state screening, and unrestricted boundary repair can break alignment. |
| **M08: gravity** | **Conditional OPH reconstruction exists; this candidate's physical attachment is open.** Finite entropy and controlled charge-collar constructions supply explicit candidate inputs. | Realize the flagship's G1–G6 premises on this same family, including its physical EC factors, conserved stress, universal coupling and gravitational constraints, or justify a different route. The nontrivial bulk collar factors do not discharge the physical gauge-edge premises; the tested spatial-symmetry identification fails. |
| **M09: parameters and P** | **Relative-capacity identification and conditional kinetic-law uniqueness constructed.** Reciprocity reconstructs normalized capacities for an admissible connected response. The same-capacity linear reversible completion uniquely gives its canonical equation; exact kinetic discriminators survive clock/mass/speed refits. | Select the response, absolute action scale, inertia identification, mass and physical clock. The counterfamily shows why configuration equilibrium alone cannot do this. No P-to-amplitude, temperature or rate map follows. |
| **M10: physical matter** | **Declared free scalar benchmarks joined to repair stiffness.** The canonical massive 3D scalar preserves the existing full field/instrument theory; periodic massless 2D and weighted arc versions remain alternatives. | Attach statistics, physical sectors, interactions, stable poles and charges to the same process. Supplied Bose CCR and a free scalar limit do not select particle content or establish the fermionic compiler premises. |
| **M11: cosmology/dark sector** | **Outside the active minimal target; open if claimed.** Existing conditional cosmology constructions are retained. | Attach common expansion/source dynamics, stress, physical scales and observational transfer; add abundance, pressure/slip and lensing where used. |
| **M12: eligible physical test** | **Internal discriminators complete; empirical test open.** Decay versus oscillation, retained-record menus and alternative inertias distinguish supplied laws. One local preparation/effect distinguishes N2/N3; private preparations distinguish N4's load response from its single-product controls. The corrected single-port tracial statistic distinguishes N4 from the matched native heat description at the certified couplings. | Select a physical observable and rejection rule; enumerate freedom, calibration, uncertainties and baseline, then freeze an unavoidable relation before independent comparison. The finite probes are diagnostic resources, not independently realized physical instruments. |

M01–M12 are a dependency map. Develop the connected cluster needed by a
specific claim; they are neither twelve sequential hurdles nor a requirement
to solve all physics at once. Finite results stay completed unless a new test
invalidates them. Enlarging the physical scope adds an explicit obligation.

### One process, with geometry, gauge response and clocks read from it

The existing arguments constrain the shape of the repair law, not a minimum
number of mechanisms. The scalar-equalizer counterexample below does not
prove that a simple common law is impossible. The target is
**one microscopic process with several verified readouts**, rather than one
scalar matrix relabelled as geometry, gauge current and time. Its shared
states, operations, retained records and reduction maps must be explicit.
Independent fixtures placed beside one another do not establish this connection.

**Clock convention.** The target is observer-local clock readouts and compatible
comparisons tied to the same event history, without a preferred universal clock
or common instantaneous slice. A shared time parameter in a finite benchmark
is a supplied model choice. For example, the
[scalar action clock](../paper/tex_fragments/SOURCE_SCALAR_CLOCK.tex) assigns
one duration to a whole configuration layer; extending it to independent local
clocks requires further work. The
[count clock](../paper/tex_fragments/SOURCE_COUNT_CLOCK.tex) instead reads
proper-time ratios under its stated order/count-limit and ancestry assumptions.

| Existing argument | What the same repair process must supply | Present connection and limitation |
| --- | --- | --- |
| Rank-three carrier and slow band: [microphysics](../paper/screen_microphysics_and_observer_synchronization.tex), `thm:echosahedral-source-selector-microphysics` and `thm:a5-register-commutant` | Its own response/current records recover the declared carrier operator, the selected slow projector and the Gram geometry. State the regime and normalization of that readout. | At uniform binary equilibrium, the native recorded model reconstructs `L_ico`, with `4 P_slow=G` of rank three, and the rate-independent relation `κ γ_s/ν_t=(5−sqrt(5))/6`. The golden scalar is a separate supplied field construction; its eight-site generator is a different object. Their microscopic-to-field reduction and compatible local clock readouts need explicit maps. Rank three alone is not a bulk dimension or a physical generation count. |
| Finite gauge Lie algebra: [gauge paper](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex), `thm:current-algebra` | A source-bound injective twelve-port reversible response, compact commutator closure, and A5-equivariance with the twelve-port permutation module, realized by endogenous inner proper-carrier transport on that same response. Reuse the classification once these inputs are realized. | The existing six-level response supplies the compiler's controls. Choosing those controls does not derive a gauge-field/matter action from scalar repair or show that autonomous events select the response. Global form, physical couplings and refinement retain their separate premises. |
| Operational clocks: [count clock](../paper/tex_fragments/SOURCE_COUNT_CLOCK.tex), [scalar action clock](../paper/tex_fragments/SOURCE_SCALAR_CLOCK.tex), and [curl/Stone comparison](../Lean/Screen/CurlStoneClockBridge.lean) | The same retained history and accessible records instantiate a clock theorem's hypotheses; the resulting time must intertwine the candidate's observable evolution. Charge readout, transport and service costs. | Count clocks require controlled order/count limits and complete ancestry. The scalar clock uses a supplied action and nondegenerate layer data. The curl comparison has explicit duration/frequency calibration premises. None identifies canonical model time with processor pulse count by notation. |

The [existing record-derived geometry and clock relation](#current-derived-port-geometry-and-clock-relations) is the positive starting point. Keep its
nearest-integer event rule, declared schedule, same-record response tests and
preparation controls when testing an extension; do not replace that evidence
with a separate scalar spectrum match.

A useful restriction follows before any simulation: a self-adjoint linear
response on the oriented 24-port register that commutes with `A5 × Z2` has
**eight sector coefficients**, by the existing register-commutant theorem.
This can replace an unrestricted matrix search for that response class.
Locality, admissibility and measured response can constrain it further.
It does not fix the eight coefficients, nor restrict the reversible gauge
current to an eight-dimensional commuting algebra: those are different objects.

There is also an existing decisive negative control. The scalar seam
equalizers have commutator span of rank 49 and Lie closure `so(11)` of
dimension 55 on the zero-sum port space, as documented in the gauge paper.
They are not the required twelve-dimensional compact current. A simple
diagonal phase lift is abelian; adding the tested adjacency tangent gives
`u(12)`, which also misses the target. These results guide candidate rejection;
they do not reject every non-diagonal reversible response on the carrier.

Match carrier actions as well as dimensions. In the six-level response of
[source record gluing](../paper/tex_fragments/SOURCE_RECORD_GLUING.tex), the
simple su(2) ideal uses the Galois-conjugate triplet; native slow geometry uses
the `5−sqrt(5)` band. A proposed identification must justify the intertwining
map. Keeping internal response and spatial geometry distinct is also allowed;
their equal dimension alone does not identify them.

**Working hypothesis under test:** a local operation retains a reversible,
phase-sensitive private response and produces repairable public records.
Forgetting or retaining different records gives distinct, explicitly related
views of that operation. A quantum instrument with retained outcomes is one
available candidate description, used by the native C/E constructions. The
[endpoint-current construction](#a-local-response-driven-by-the-actual-repair-record)
tests a finite attachment, with a flat local alternative. Neither is selected
uniquely by these arguments. The noncommuting response,
its coupling to repaired records and its event clock must be derived or tested
together. Simply appending an independent gauge sector or oscillator clock
would leave the shared-process obligation unfilled.

For each candidate, use this order of work:

1. Pin the microscopic states, local event and retained-record contract.
2. Derive its response, reversible perturbations and record-clock candidates
   before importing a desired continuum action.
3. Match the actual hypotheses of the three theorem interfaces above on that
   same family, including cross-response under the same interventions.
4. Use the scalar benchmark and compiler as component checks where applicable;
   reject or revise a candidate that fails a required source interface.

This directs the simplification search toward a common source law. It does
not require all physical limits to be solved at once, and it does not turn
conditional scalar compatibility into a complete physical model.

### The reverse-engineering table: hypotheses as well as results

The [model–theorem premise map](MODEL_PREMISE_MAP.md) connects this M01–M12
summary to the existing canonical premise register, claim registry and gravity
ladder. It records what each construction supplies, what it assumes,
which existing result can be reused, and what evidence would discharge the
remaining obligation **on the same model and at the claimed scale**.

A theorem proves an implication; a model can construct its inputs. That is
real progress even before physical selection. However, moving a premise into
the candidate's definition does not derive it. When a construction uses new
hypotheses to fill an old interface, keep those hypotheses in the total input
budget. Shared requirements are counted once by identity, with every use
retained; competing branches are not combined into one mandatory checklist.

The map covers every registered canonical PR identifier, including
inactive branches. The larger named-assumption/import inventory is reported
separately. This is a checked cross-reference and research ledger, **not a Lean
proof that every antecedent in every paper has been extracted or discharged**.
The original claim classifications and premise dispositions remain authoritative.
In particular, an excluded physical conclusion is not silently turned into an
unproved hypothesis of a completed mathematical theorem.

| Existing result to reuse | Connection supplied by the current construction | What has to be supplied |
| --- | --- | --- |
| [Action-measure identification](../Lean/Geometry/SourceActionMeasure.lean), `normalized_symmetrizer_unique` | A connected response admitting a positive diagonal symmetrizer fixes normalized relative capacities; the weighted repair factorization uses that same capacity/stiffness interface. | Selection of the response, absolute scale and equality of repair capacity with inertia. |
| [Source-to-Hamiltonian composition](../Lean/Variational/SourceToHamiltonianComposed.lean), `source_to_hamiltonian_composed` | Reuse the action/Legendre machinery once its real enrichment and variational hypotheses hold. The present quadratic completion gives its own explicit conditional action. | The theorem does not choose the real enrichment, kinetic law, quantum representation or physical clock. |
| [Common scalar packet](../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex), `prop:common-scalar-continuum-limit` | The canonical completion and constrained local field coordinates retain the original action, fields and vacuum, so the same free scalar limit applies. | Native implementation and physical identification; no interacting limit follows. |
| [Local preparation](../paper/tex_fragments/SOURCE_SCALAR_LOCAL_PREPARATION.tex), “Preparation with an active free field” | The local ramp supplies a bounded initial-state error for the same original **q=3** coupled vacuum, which propagates by channel contractivity. | Supplied blank/control/time; it does not supply the separately used q=5 or q=233 numerical preparation budgets. |
| [Finite instruments](../paper/tex_fragments/SOURCE_SCALAR_FINITE_INSTRUMENT.tex), `prop:source-scalar-finite-pulse` | Exact force and retained-pointer histories reuse the solvable Gaussian pulse and keep its relative pointer phases. | Native pulse/service accuracy and complete readback. The fixed commuting histories have a common onsite code with controlled intermediate tails. |
| [Sequential instruments](../paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex), `thm:source-scalar-sequential-instrument` | Reuse the joint field/record channel and trace-distance contraction; do not rederive a separate probability law for each history. | Source production of preparations, instruments and readout; adaptive feedback requires its own explicit history. |
| [Time refinement](../paper/tex_fragments/SOURCE_SCALAR_TIME_REFINEMENT.tex), `prop:source-scalar-time-refinement` | Reuse the fixed-action split-step structure. The local ramp has its own midpoint variation and stable splitting estimate, giving a full Gaussian state bound of order N^-2. | The fixed onsite preparation/history cutoffs have separate state-specific bounds. Native pulse/service error and a physical clock require their own inputs; the Gaussian estimate is not an arbitrary-input channel norm. |
| [Gravity premise ladder](registers/GRAVITY_PREMISE_LADDER.md) | Reuse the conditional entropy-to-gravity implications once the same candidate supplies their typed inputs. | Physical edge alignment, joint state/geometry stationarity, stress, common continuum and calibration remain distinct; scalar preparation does not close them. |

## 2. The iterative minimalist method

**Compatibility with OPH's mathematical requirements matters; reproducing
every detail of the current implementation is not the goal.** We are
reverse-engineering hardware and software. Protected registers, a quantum
representation or a changed operation are permitted candidate hypotheses
when their purpose, cost and compatibility are explicit. The simulator is a
tested starting realization, not an immutable specification of nature.

**Find the simplest shared repair rule that satisfies the known constraints;
let a decisive test determine the next change.** Simplicity is a research preference, not a proof that nature
selects the shortest formula. Preserve the actual carrier and load-bearing
premises of every imported result.

1. **Name one missing interface.** State which objects must be related, which
   result needs that relation, and what would falsify the proposed attachment.
2. **Map an established theorem and existing construction first.** Use the
   [reuse contract](AGREEMENT_RESEARCH_FOUNDATIONS.md) to identify its state,
   restriction, boundary, fault and scheduling assumptions. Record the
   candidate-specific missing premise. Current, recorded order,
   first passage and regional entropy have supplied useful objects without
   adding independent fields. Test those before introducing new state variables.
3. **Add only what a demonstrated gap requires.** If the existing source fails,
   specify the smallest compatible hardware or law change. Explain why it is
   needed and what demonstrated capabilities it preserves.
4. **Test one common process.** Use the same states, operations, records and
   maps. A matching mean, symmetry, stationary histogram or parameter name
   does not join otherwise separate constructions.
5. **Retain countercontrols and bounded decisions.** Pin sources and inputs,
   declare the computation cap and stopping rule, compare an independent
   calculation, review the derivation, and replay. Do not retune after failure
   and call the original prediction successful.
6. **Update the completion ledger.** Close the precise finite subproblem,
   retain the remaining parent criterion, and move a rejected route to the
   historical section. Do not repeatedly reopen a completed contract merely
   because physical interpretation remains to be supplied.

Each extension declares state, law, preparation, schedule, reference, units,
parameters, quotient, readout and memory budget. Label each as an axiom input,
derived object, candidate hypothesis, calibration or target-exposed choice.
Adaptation is allowed if its state and update belong to the same process.

### General findings that guide the search

**Prefer tests that eliminate laws, and count assumptions rather than symbols.**
The following lessons govern candidate selection and the
[law-search issue](https://github.com/FloatingPragma/observer-patch-holography/issues/1076).
They are methodological findings, not additional microscopic axioms.

| Finding | Consequence for the search |
| --- | --- |
| Locality, convergence and protected memory admit many different internal responses. | Before adding a variant, name the requirement that could reject it. Passing a shared convergence test cannot select the response algebra. |
| Exact preservation of the old load channel for every internal input excludes response-to-load feedback. | Preserve demonstrated reference behavior where needed, but do not freeze an implementation choice that prevents the missing capability. Specify which preparations and readouts must be preserved. |
| Reusing source matrices does not transfer the source's entire algebra or operation grammar. | Bind states, overlap restrictions, admissible interventions and records before importing a theorem's conclusion. Equal dimensions or matching spectra are insufficient. |
| Algebraic regrouping can shorten a formula without reducing its choices. | Count state variables, independent coefficients, preparations, clocks and services. Declare fixed coefficient ratios as choices even when they are not fitted. |
| A demanding diagnostic can obscure a simple distinction. | Try one local state, one event and one local effect before entangled references, long coherent histories or large simulations. N2/N3 admit this reduction. |
| Mixing closed spatial modes produces apparent local memory. | Derive the null for the actual observable and sampling protocol. N4's single-port statistic subtracts the native mode-mixing term; changing from global attempts to local stopping times needs a new null. |
| A successful compiler implements a supplied law; it does not select that law. | Prioritize discriminating microscopic tests. Use the scalar compiler and continuum benchmark where a concrete shared-source question needs them. |
| A conjugated channel may only change coordinates if its observations are conjugated too. | Fix operational preparations and observables before comparing laws; track changes in record meaning explicitly. |
| Conditional mathematical consistency does not establish physical correctness or uniqueness. | Retain viable alternatives and allow a failed attachment to revise the candidate. Physical selection requires an identifiable prediction and independent evidence. |

The criterion is the **simplest shared rule that satisfies the known
constraints**, rather than the simplest formula that passes a permissive test.
An additional term earns a place only by resolving a demonstrated obstruction
and surviving its distinguishing test. A necessity argument restricted to one
candidate family must not be promoted to a theorem about every possible law.
If all retained candidates make the same accessible predictions, record that
non-identifiability and seek a stronger justified interface or observation;
adding indistinguishable variants does not solve selection.

#### Search for programme-wide obstructions

The [active obstruction search (#1078)](https://github.com/FloatingPragma/observer-patch-holography/issues/1078)
asks whether the adopted OPH programme can be ruled out. A decisive negative
result is useful. Pin the axiom revision, target conclusion, mandatory physical
premises, admissible interpretations and calibration freedom before testing.
Distinguish those commitments from optional implementation choices.

Search for a contradiction in the adopted basis, impossibility of jointly
realizing necessary structures on one source, independence of allegedly
derived physical inputs, or an unavoidable failed physical prediction. Reuse
the existing claim graph and premise map to locate the actual theorem
statements; do not conjoin alternative branches as mandatory premises.
Count the content of assumptions, including the response architecture in A1;
the number of axiom headings does not measure explanatory economy.

A candidate or family no-go excludes the full formulation only if it covers
every admissible realization. Premise-checked countermodels with different
outputs can refute a claimed derivation; they do not automatically eliminate
all shared predictions. An unfinished construction proves neither impossibility
nor independence. Conversely, failing to find an obstruction supplies no
physical confirmation.

Record the strongest surviving alternatives and obtain independent scope
review before promoting an exclusion. Changing a core commitment creates a
revised target and preserves the original negative result. A maintainer may
also stop for insufficient predictive or explanatory return under explicit
milestones and costs; that is a research decision, distinct from a proof of
impossibility. Opening this search makes no stop decision. Existing
[claim- and branch-specific falsifiers](OPH_FALSIFICATION_PROGRAM.md) retain
their exact scope.

#### A deletion test exposes what convergence cannot select

For N2, the load-preservation, locality and completion arguments work with
any finite skew-Hermitian family A_p in place of D_p. Covariance also transfers
when that family has the declared carrier action. In particular, retaining
only the source's lower block,

```text
A_p = diag(0_3, hat(σ(v_p)/2)),
```

preserves those properties and the protected 3+3 populations. It is a
noncommuting, carrier-covariant response, yet its port-response map has rank
three and a nine-dimensional kernel. It fails the injective twelve-direction
premise of the [gauge theorem](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex).
Thus convergence and covariance do not replace that premise. This is an
analytic deletion control, not an admissible complete OPH source.

The [source witness](../paper/tex_fragments/SOURCE_RECORD_GLUING.tex) uses
`direct_sum_(12 ports) M_6`, with center `C^12`. The candidate adds a quantum
load space and twelve simultaneous tensor factors. Reusing D preserves its
matrix relations, not that direct-sum center or complete operation grammar.
Establishing an admissible common source is a selection requirement.

Reuse the [finite entropy-completion result](../code/quantum_information/ALGEBRAS_AND_REPAIR.md)
for the binary candidate's trace-preserving conditional expectations and
their specified full retained-algebra fibers. Its tracial reference makes
modular invariance automatic; this does not select a coupling or agreement
algebra. There is no need to add a separate entropy postulate to obtain this
finite conditional result, or to adopt the repair RFC to use it.

### Work ownership and claim accounting

**The simplest shared rule is the research target.** Added mechanisms must
answer a demonstrated missing connection; a simpler rule is preferred when
it satisfies the same constraints with fewer independent physical inputs.
Compiler cost is an implementation cost, not evidence about how complicated
the law of nature must be.

The programme targets one common microscopic process whose constructed
objects satisfy the premises of the physical conclusions being claimed.
Separate conditional theorems can hold on different states, scales or laws;
their validity alone does not construct that common realization. The research
question is which inputs the repair law supplies, which it leaves free, and
which observable can distinguish the resulting candidate.

The first outputs under this method are the
[recovery contract and linear extension calculation](AGREEMENT_RESEARCH_FOUNDATIONS.md#7-work-products-and-scientific-contribution):
define the native legitimate set and fault span, identify its actual scheduler,
and factor the scalar comparison into its fixed overlap maps. Compute the
protected nullspace, boundary extension residual and uniqueness condition.
Reuse established conclusions when these objects meet their hypotheses;
then spend model-selection effort on the missing source, instrument and clock
connections and a transient observable that distinguishes viable laws.

| Owner | Responsibility and bounded result |
| --- | --- |
| [#1026](https://github.com/FloatingPragma/observer-patch-holography/issues/1026) | Common-model construction umbrella: one law, compatible theorem attachments, explicit remaining assumptions and tested histories. |
| [#1076](https://github.com/FloatingPragma/observer-patch-holography/issues/1076) | Find the simplest shared rule compatible with carrier, reversible-current and clock interfaces; pin its inputs and test dynamical alternatives. |
| [#1077](https://github.com/FloatingPragma/observer-patch-holography/issues/1077) | Complete one source-coded preparation, evolution and retained-record history with a shared error and resource budget. |
| [#1025](https://github.com/FloatingPragma/observer-patch-holography/issues/1025) | Develop physical tests alongside construction: observable, interpretation, calibration, baseline, uncertainty, exposure and rejection rule. |
| [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033) | Audit and repair existing proofs, code, evidence and scientific wording, including effects on the candidate. |
| [#1078](https://github.com/FloatingPragma/observer-patch-holography/issues/1078) | Search for programme-wide inconsistency, joint-realization no-go results, failure of explanatory unification and empirical exclusion; distinguish scoped rejection from a practical stop decision. |

These are parallel responsibilities. A link to another issue is not a hard
dependency; declare a dependency only when a concrete required output is
absent. The two construction tasks use a pinned conditional candidate while
its physical selection is examined. Extend that family towards continuum
control, interactions and joint matter/geometry only through a specific
missing interface and a bounded test. The M01–M12 and RE groups are a
scientific map, not a mandate to create an issue for every row.

No standalone finite-theorem campaign, theorem-count target or generic
formalization quota belongs to this programme. A proof is useful when it
repairs a defect or discharges a necessary model premise. Try existing
theorems and simpler constructions first. Any added mathematical structure
must resolve a demonstrated obstruction and enter the candidate's total
assumption and resource budget; renaming an input does not derive it.

Each investigation declares its calculation, controls, acceptance/rejection
rule and local compute budget. Close a bounded task with a supported result,
a candidate rejection, or a precise obstruction. An obstruction changes the
next experiment; it does not establish or refute every OPH model. There is no
calendar deadline for declaring nature solved and no automatic promotion
from an issue closure to physical completion. Large runs and cloud spending
require their own authorization.

Use the existing records together:

- **GitHub issues** own live work, dependencies and concise decisions. Keep
  exploratory logs and unapproved private contributions local.
- **The [candidate register](MODEL_CANDIDATES.md)** owns the concise comparison
  of concrete laws and their properties. This document supplies their
  derivations and test summaries; the [premise map](MODEL_PREMISE_MAP.md)
  owns model-to-theorem connections. Edit
  [`tracking/model_theorem_map.json`](../tracking/model_theorem_map.json) and
  regenerate the map; do not create a competing claim registry.
- **The [claim registry](../claims/claim_registry.yaml) and scientific
  registers** own established claim identities, scope and evidence. Preserve
  frozen prediction contracts and historical classifications.

Every claimed connection records the candidate revision, exact theorem/claim
and premise references, matching state/domain/refinement, what was
constructed, all retained or added assumptions, evidence and countercontrols,
and the condition for closure or rejection. State mathematical construction,
physical selection and empirical support separately. The reference checker
does not certify semantic premise discharge or exhaustive antecedent coverage.

### Research programme and paper-release hold

**Paper releases and broad restructuring of the prose papers are on hold.**
Keep model work in its issues, this status document and the premise map while
the common construction takes shape. Repeatedly reorganizing the papers
around component results would obscure which physical premises the same
process actually supplies.

Audits, necessary correctness fixes and dated errata continue, including
affected paper passages and rebuilt review PDFs where required. Existing
releases and frozen evidence retain their history. Preview builds and CI
continue; passing them does not lift the hold. A necessary corrective release
must identify its correction scope and receive an explicit maintainer
publication decision, without adding a broad model narrative.

Reconsider the hold at a recorded maintainer review when one explicit
candidate has a complete tested finite history, supplies a substantial
connected set of existing theorem premises on the same family, and exposes
its remaining physical assumptions with an error and resource account. A
meaningful model milestone is sufficient; completing every physical sector
is not required. The review states which claims warrant a paper revision and
whether release work may resume. The resulting exposition should follow
state and repair law → constructed objects → applicable theorems → physical
consequences → tests.

The connected interfaces below guide bounded tasks; only the issue owners
above assign active work:

| Next question | Smallest proposed move | Decisive completion or rejection condition |
| --- | --- | --- |
| **Can one simple local rule supply the carrier, current and clock?** | Reuse the native response and gauge/clock interfaces above, with one shared state and event record. | Exhibit the response reduction, reversible current and same-history clock map; independent attached fixtures or the rejected scalar-current identification do not qualify. |
| **Can N4's memory distinction be read locally?** | Use the corrected single-port endpoint statistic below, retaining the tracial preparation and one-/two-global-attempt schedule. | Realize preparation and readout in the source, bound the complete measured statistic's error, and derive a new null before substituting observer-local ticks. Positive uncorrected correlation excess alone fails this gate. |
| **Why should repair capacity also be inertia?** | Attach the tested three-mode curvature or distance-two acceleration discriminator to an admissible preparation/readout. | Derive the identification from microscopic operations or constrain it independently. Equilibrium and the shared position-observable continuum limit do not select it; retain clock, velocity and linearity premises. |
| **Which weighted source law is selected?** | Reuse the exact stiffness factorization and action-measure theorem. | Justify conductances, reference, mass and boundaries independently of the target action; reciprocity identifies normalized capacities given an admissible response, not that response's physical selection. |
| **Can the finite quantum components become a full source history?** | Join the local vacuum-preparation schedule and exact retained-pointer history through the existing state-error contract; reuse source controls for implementation. | Use the constructed common onsite preparation/history and finite source compiler to specify an executable word with complete outcomes and a shared resource/error budget. Fixed-path temporal and occupation errors are controlled; source precision, storage, drift and service costs require implementation evidence. |
| **Which model time is physical?** | Attach the canonical flow and existing local split implementation to a declared source clock. | Establish preparation, readout, event ancestry and all execution costs; a model-time parameter or a local matrix is insufficient. |
| **Does the archived spherical family converge?** | Use its uniform internal gap, exact memory and positive local coarse approximants. | Prove compatible spatial operator/state convergence on that different routing; neither of the two completed continuum examples substitutes for this. |
| **Can the field assumptions be supplied by the microscopic process?** | Test the canonical phase reference, state preparation and statistics against retained records and resources. | Justify or replace the supplied Bose CCR and vacuum. The constrained field decoder is local; the free Fock representation requires justification. |
| **Which interaction is the smallest justified extension?** | Start with actual multi-occupation repair and its retained correlations before adding a fitted continuum interaction. | Derive a controlled many-body limit or a surviving interaction and a distinguishing observable. A free field does not determine the Standard Model. |
| **Can the collar become a physical edge?** | Seek an internal paired-representation map compatible with the update, observables and state. | Supply Gauss/sector alignment, area/stress response and common refinement. Spatial covariance and bulk dimensions alone do not supply them. |
| **What is the first eligible physical test?** | Attach an unavoidable relation to a calibrated physical preparation and readout. | Enumerate remaining choices, uncertainties and baseline; freeze the prediction before independent comparison. |

## 3. The candidate repair law being tested

### A local response driven by the actual repair record

The current microscopic test joins native repair and the existing reversible
response through the same recorded event. For a seam e=(a,b), the current j
fixes the native load change. Keep its current-record Kraus operator L_r and
apply U_r=exp[θj(D_b−D_a)] at both endpoint response factors. The D_p are the
existing six-level [source response](../paper/tex_fragments/SOURCE_RECORD_GLUING.tex),
with an explicit incidence-preserving port identification. The working load
space has 78 total-two configurations; the declared response space is
`(C^6)^⊗12`. With identity on the other response factors, the event is

```text
K_(e,c,j) = L_(e,c,j) ⊗ U_(e,j)^(a) ⊗ U_(e,j)^(b),
U_(e,j) = exp[θ j (D_b−D_a)].
```

The local factors and θ are declared hypotheses.

#### What is supplied by that one event

For any correlated load/response input, discarding the response gives exactly
the original native repair channel, branch by branch for classical records.
Its classical geometry and dimensionless clock relations therefore retain
their original preparation and scheduling conditions. The coherent record
off-diagonal blocks can change. Disjoint seams commute because they operate
on disjoint endpoint factors; proper-carrier covariance follows from the
existing response transport. A single shared response factor is a useful
rejected control: it permits disjoint seams to interfere.

The main simplification is exact reuse of the completion argument. Set
s_e(x)=(x_b−x_a)/2 and let W_e apply exp[θs_e(x)(D_b−D_a)] at the two endpoints,
controlled by the load configuration. A repair changes s_e by j, hence
K_r=W_e(L_r⊗I)W_e†. The coupled C_e and E_e are conjugates of their original
channels and inherit C_e^n=E_e+2^(1−n)(C_e−E_e), for n≥1. No additional
relaxation law is needed for this single-seam result.

On the 66-state binary shell, the existing Hilbert–Schmidt positive-contraction
and quadratic-form argument transfers too. Any fixed normalized positive seam mixture
of the coupled atomic channels converges to its common fixed algebra after
fresh records are discarded. This is convergence of reduced states, not of
the growing coherently retained history. The completed-channel mixture has
the same fixed algebra, and their gaps satisfy
γ_E/2≤γ_C≤γ_E on its orthogonal complement. No absolute gap or continuum-uniform
rate has been computed for this enlarged model; the old uncoupled numerical
gap cannot be imported through different seam-dependent W_e.

The local 3+3 block populations commute with every D_p and survive every
event. Thus this finite model can settle while retaining protected local
distinctions. This does not establish arbitrary private coherence, a unique
equilibrium or the full physical observer algebra. Sustained recorded
operation needs the same explicit storage/disposal policy as the native law.

#### A simpler alternative keeps the selection question testable

Let J=Σ_pD_p and Dbar_p=D_p−J/12. The endpoint operations
exp(−θjDbar_a) and exp(+θjDbar_b) define another local current-record lift.
They are induced by the single global change of coordinates
W_flat(x)=⊗_p exp(θx_pDbar_p), so every closed load history has identity
response holonomy. This alternative preserves the same native repair and
locality. Distinguishing it from the main lift requires response-history
measurements, not another match of the unchanged load spectrum.

In the main lift, a closed load loop can leave local path memory even with
commuting response matrices. Loop memory alone must not be called physical
gauge curvature. The tests separate order-sensitive noncommutation from this
ordinary path memory with a commuting-generator control and retain all
outcomes of the fixed short history.

#### Fixed tests and what they discriminate

The frozen cases use θ=0, 1/16 and 1/8; all 78 total-two load states, all
30 native seams and all 60 proper carrier transports are checked. The local
support argument covers all 315 disjoint seam pairs. Exact current grouping
and C/E identities are combined with correlated-input matrix checks and a
nonzero C−E coherence test. All **268 checks** passed and the sole independent
fresh-directory replay reproduced the scientific output byte for byte.

At θ=1/8 the shared-processor negative control has a disjoint-order channel
infidelity of approximately 0.000844, whereas endpoint support preserves
commutation. On a fixed three-seam schedule the two local alternatives have
complete record-plus-response trace distance **0.137390**, numerically, for
the declared product of local maximally entangled inputs and references.
All eight equiprobable coin records are included. This is a source-model
comparison with supplied probe resources, not an experimental observation or
a physical coupling estimate. The commuting-generator control also retains
loop memory, so that signal alone is not evidence for non-Abelian curvature.

The all-seam convergence argument reuses the existing positive-contraction
proof; the enlarged network matrix and its absolute gap were not computed.

#### A simple local readout distinguishes the alternatives

Prepare load at ports {0,2}, attempt seam (0,1), and prepare endpoint 0's
private factor in P=|0><0|. With the same source labeling, put
`A=[D_1,P]`, `a=||A||_F` and `F=(I+A/a)/2`. A is Hermitian and `||A||_op≤a`,
so F is a valid two-outcome effect. Measure F only at endpoint 0. Averaging the wait and move
outcomes gives

```text
p_k(θ) = Tr(F [P+exp(θk) P exp(−θk)])/2,
k_N2 = D_1−D_0,       k_N3 = −(D_0−J/12),       J=Σ_p D_p.
```

Since `[J,P]=0`, the difference has derivative `a/4>0` at θ=0. The source
formula gives `a=(1+sqrt(5))/2` for this chosen probe. At θ=1/8, the numerical
probabilities are 0.549068836 for N2 and 0.499973968 for N3, a difference
of 0.049094868. At θ=1/16 the difference is 0.025096764; at zero it vanishes.
Both probabilities were checked separately against independent source
reconstruction and exponential evaluation. All **67 fixed checks** passed,
with an identical isolated replay.

This diagnostic uses one supplied local state/effect and both classical
outcomes. It removes the earlier comparison's entangled reference and
retained history. It does not construct an autonomous source preparation or
readout, give interval-certified probabilities, or compare with physical data.

#### Remaining connection

This is a finite one-way coupling: private interventions cannot change the
unconditional native repair probabilities. Reciprocal backreaction is absent.
For an input that factors into load and product private states, the
record-conditioned response factors across ports and creates no
private-to-private entanglement.
The added tensor factors and conditional response operations have not been
shown to satisfy the full A1 response grammar on the enlarged source.
Conserved difference drives omit the central relative-phase direction of D;
the complete twelve-direction response remains a supplied interface.

The candidates leave θ, physical clock, microscopic state selection, responsive
geometry and interacting continuum dynamics open. Their common finite
convergence supplies a useful requirement, but does not select the physical
law or close an M01–M12 parent entry. The next useful test must discriminate
the response histories or supply reciprocal coupling while preserving locality
and a justified completion argument. Reciprocal coupling requires relaxing
exact preservation of C on specified internal preparations: retaining that
marginal for every input mathematically forbids feedback. The archived kernel
is a reference benchmark, not an immutable microscopic implementation.

### A reciprocal coherent-exchange candidate

N4 asks whether the same local repair can respond to its private sector. It
changes the coupling, using the binary shell and the same twelve six-level
factors as N2. No additional state variable is introduced. One coupling θ and
a fixed relative choice of response directions are supplied.

#### The simpler attempt has an exact obstruction

On a seam's odd-occupation block {|10>,|01>}, C depolarizes the load qubit.
It preserves the even block and multiplies even/odd coherences by 1/2. It is
therefore covariant under every local number-preserving load unitary. For
any Hermitian such A and private skew-Hermitian G, conjugating C by
`W=exp(θ A⊗G)` leaves its reduced load channel unchanged, even on correlated
inputs: diagonalize G and apply load covariance on each surviving block.
Thus a single product of this form cannot supply feedback. This is an
obstruction within that family, not a lower bound on every possible law.

#### One exchange, with its choices exposed

Let `Z=N_a−N_b`, `X=S_e Z²`, `Y=i S_e Z`; X and Y are Pauli operators on
the odd block and vanish on the even block. Using the existing source center
`J=Σ_p D_p`, set

```text
g=D_b−D_a,                  h=D_a+D_b−J/6,
G=g^(a)+g^(b),              H=h^(a)+h^(b),
T_e=Y⊗G+X⊗H,               W_e=exp(θ T_e),
K_r=W_e (L_r⊗I) W_e†.
```

Equivalently, with `s_+=|10><01|`, `s_−=s_+†` and `B=H−iG`,
`T_e=s_+⊗B−s_−⊗B†`: one coherent exchange and its adjoint. This compact
form does not eliminate the supplied equal weighting of its two real
directions. Both directions have zero-sum port coefficients, so this coupling
does not restore the missing central response. The amplitude obeys
`[B,B†]=2i[G,H]`; a normal amplitude reduces
to private-controlled load rotations and cannot evade the preceding
obstruction. The chosen source amplitude has a nonzero commutator.

The event is local, conserves total occupation, respects proper carrier
covariance and protects both endpoint 3+3 block projectors. Conjugating E
with the same W transfers the exact C/E completion identities. The existing
binary-network positive-contraction argument gives convergence of fixed
positive seam mixtures to their common conditional-expectation algebra.
Its full algebra and absolute gap are uncomputed. This concerns reduced
channels with records discarded, not unlimited retained histories.

#### Feedback test and its operational cost

On a neighboring seam, keep the odd load state fixed at I_2/2. The original
diagnostic labels it (0,1); in the authoritative response chart used by the
reproduction below, a representative is (0,2). Put
`q=i[g,h]`, `ell=sqrt(49+19 sqrt(5))` and prepare each private endpoint in
`σ_±=(I±q/(2 ell))/6`. The source formula gives `||q||_F=ell`, so these
are normalized full-rank states with eigenvalues at least 1/12. Measure the
original endpoint occupation N_a, without changing its definition between
preparations. Independently expanding the conjugated channel gives

```text
P(a occupied | σ_+⊗σ_+) − P(a occupied | σ_−⊗σ_−)
    = −θ² ell/3 + O(θ³).
```

At θ=1/8 the numerical probabilities are **0.482020673** and **0.517979327**;
at θ=1/16 they are 0.494257288 and 0.505742712. Both one-term controls give
1/2 and preserve the load marginal on a complex correlated input, as the
analytic obstruction requires. All **187 fixed checks** and an identical
isolated replay passed, including
independent source/derivative controls, nonzero C−E coherence, locality,
covariance and record semantics. Finite-angle probabilities are numerical
witnesses; the Taylor expression is not a certified finite-angle error bound.

The branch index j denotes a transformed occupation transfer. It obeys
`[W_e(N_a⊗I)W_e†,K_r]=−j K_r`, while the identity with N_a fails in the
test. These transformed observables depend on the seam. A shared observable
and record interpretation must therefore be established before using a clock
or geometry theorem. Measuring loads before and after each event would
define a different, additionally dephased instrument.

N4 closes a finite feedback construction under supplied quantum inputs. The
following readout recovers its port geometry without reinterpreting j as a
bare current. Full-78 transient repair, complete source-response grammar and
physical clock selection are not supplied. On one seam N4 is unitarily
equivalent to C; its distinction is relative to fixed preparations and bare
occupation readouts. No M01–M12 parent is closed by these finite tests.

#### Bare-occupation geometry and a clock-independent memory test

Fix the existing N4 law on the 66 binary total-two load states and twelve
six-level factors. Let τ be normalized trace on that entire space and Φ the
uniform, independent 30-seam mixture, with records discarded. This is a
supplied stationary preparation and schedule, not a uniqueness claim for the
stationary state. Use the **same bare occupations N_p on every seam**.
Prepare τ, measure and retain the initial binary loads, apply one or two
unmodified repair events, and measure final loads. The initial measurement
leaves τ unchanged when its outcomes are averaged; there is no intervening
measurement. For `f_v=Σv_p N_p`, the measured correlation is exactly
`τ[f_v Φⁿ(f_v)]`. These are independently prepared experiments at each n,
not a repeatedly dephased trajectory instrument.

On one odd load block define `R=Tr_load(W† Z W)/2` and
`α=1−Tr_36(R²)/36`. Locality, number conservation and source covariance give

```text
K_pq = τ[N_p(I−Φ)N_q] = α L_ico,pq/396,
E[(X_1−X_0)(X_1−X_0)ᵀ] = 2K,
60 K/Tr(K) = L_ico,
max{1/9, 1−(4/3)(49+19√5)θ⁴} ≤ α(θ) ≤ 1.
```

The lower bound holds for **every real θ**.
Indeed R is a traceless Hermitian contraction on each of four protected
nine-dimensional private blocks, so each has squared trace at most eight.
The quartic bound follows from the Frobenius Lipschitz bound for
`cos(2θ√x)` and the exact source commutator norm. The smallest positive
eigenspace of this Dirichlet matrix recovers `G=4P_slow`, the existing
unit-diagonal rank-three port Gram. This is a finite equilibrium response
geometry; it does not identify a spacetime metric or a long-time density mode.

Does preserving that geometry also preserve the native diffusion clock?
Let P project onto centered linear occupations in the τ inner product and
`A=I−Φ`. For a centered carrier eigenvector `L_ico v=λv`, put
`c_n=τ[f_v Φⁿ(f_v)]/τ[f_v²]`. Orthogonal projection gives the exact test

```text
r = αλ/60,       c_1 = 1−r,
c_2−c_1² = ||(I−P)A f_v||²_τ / ||f_v||²_τ = b_v ≥ 0.
```

N1 and either one-term exchange control have `b_v=0`. For N4, a local
overlap calculation and interval certificate give **b_v>0 for every
nonconstant carrier vector at θ=1/8 and 1/16**. The calculation uses only
local matrices; it does not instantiate the enormous twelve-factor network.
For the rank-three carrier band, floating-point evaluations are:

| θ | α | b_v | b_v/r² |
| --- | --- | --- | --- |
| 1/8 | 0.975986704560 | 0.00004178223344 | 0.0206705310 |
| 1/16 | 0.998236412623 | 0.00000329877975 | 0.0015600348 |

The decimals describe evaluations, while the certificate proves positivity
at those exact rational couplings. With a supplied Poisson attempt rate ν,
the continuous correlation has initial decay νr and curvature `ν²(r²+b_v)`.
The ratio `b_v/r²` is invariant under a common change of time units.
Consequently no common constant clock rescaling makes these correlations a native
single-exponential heat mode. Agreement of Dirichlet matrices alone would
have missed this failure of closure of the linear-occupation subspace.

The [reproduction and proof](../code/native_response/README.md) bind the
actual source chart, protected blocks, trace normalization, local overlap
weights and finite error enclosures. This resolves the finite shared
geometry/kinetics readout for the declared N4 preparation, with a genuine
discriminator against simpler preserved-marginal laws. It neither selects
N4 over those laws physically nor restores j's original bare-current
meaning. The coupling, tracial preparation, schedule, endpoint measurement
access and physical conversion of event count remain explicit inputs.

#### A single-port test that removes spatial mode mixing

The scalar eigenmode test must not be applied unchanged to an arbitrary local
occupation. For `f=N_p−1/6`, define `c_n=τ[f Φⁿ(f)]/τ[f²]` under the same
tracial preparation and IID seam attempts. Here `τ[f²]=5/36`,
`(L_ico)_pp=5` and `(L_ico²)_pp=30`. The known compression therefore gives

```text
r = 1−c_1 = α/11,
||P A f||²_τ / ||f||²_τ = α²/110,
b_local = ||(I−P)A f||²_τ / ||f||²_τ
        = c_2−c_1²−(1−c_1)²/10.
```

Self-adjointness of Φ makes the last line an orthogonal-decomposition
identity. Native repair has `c_2−c_1²=1/1210` despite exact mean closure:
that positive value comes from mixing the spatial bands. The corrected
`b_local` is zero for native repair and the one-term controls. If the existing
edge-Gram certificate is `B≥mI`, five seams incident to a port give
`b_local≥m/100`. Thus #1085's certificate proves strict positivity
at θ=1/16 and 1/8. No new positivity certificate or interaction is required.

Measure only `N_p` at the two endpoints of each independently prepared trial.
For its spectral projector `E_a` and full load-pattern projectors `P_x`,
`Σ_{x:x_p=a}P_x τ P_x=E_a τ E_a=τ E_a`. Linearity of the intervening channel
then gives the same local joint law as marginalizing the all-load experiment.
There is no intermediate measurement, private reset or postselection. The
[derivation and controls](../code/native_response/README.md#a-corrected-single-port-readout)
therefore reduce the readout requirement without supplying its physical
implementation.

The stopping rule counts **one or two global IID attempts**. Local
receiver ticks, state-dependent schedules and timeouts define different
conditional channels; the formula cannot be transferred to them unchanged.
Here “global” counts all seams in this supplied benchmark; it does not posit
a preferred universal physical clock.
Under a separately supplied constant Poisson rate ν>0, the initial-derivative
ratio is `c''(0)/c'(0)²−11/10=b_local/r²`, invariant under one common constant
time-unit change. It supplies no calibrated observer duration.

Before reporting a measured positive value, fix preparation/readout/schedule
errors, trial allocation and a joint confidence or error rule for the whole
corrected statistic. Quadratic plug-in estimates have finite-sample bias;
a positive point estimate alone is insufficient. C and E share these first
two correlations, and other enlarged classical or quantum laws can have
positive hidden-component norms. This excludes the matched closed native
heat description under the stated contract, not every alternative to N4.

### Weighted agreement and its minimal inertial completion

For a chosen overlap e=(i,j), reconcile the two readings to their common
capacity-weighted mean, leaving every other reading untouched:

```text
E_e^M(f_i,f_j)=(a,a),    a=(M_i f_i+M_j f_j)/(M_i+M_j),
lambda_e=c_e(1/M_i+1/M_j),
Σ_e lambda_e(I−E_e^M)=M^−1 K_internal.
```

M_i is a positive local capacity and c_e a positive seam conductance. Each
map preserves `Σ_i M_i f_i`. The rates and reference are declared inputs.
This real-valued extension is distinct from the native integer rule below.
Grounded boundary seams and the supplied constant mass give the full
configuration potential `V=fᵀ(K_D+m²M)f/2`.

The scalar benchmark adds a local velocity ν and uses **the same M for
inertia**, retaining the repair restoring force:

```text
f_dot=ν,    M ν_dot=−(K_D+m²M)f,
Π=Mν,      H_can=ΠᵀM^−1Π/2+V(f).
```

This is a conservative reversible completion, not dissipative averaging:
it generally oscillates and does not converge asymptotically to agreement.
Any claimed settling, record formation or environment coupling needs its
own compatible mechanism and accounting. The conditional uniqueness result
below explains exactly which premises fix this completion.

Its simplest local first-order coordinates are
`χ=(m q+i p,Dq)`, where q and p are the mass-whitened position and momentum
and `DᵀD=M^−1/2K_D M^−1/2`. Edge entries are dependent differences; both
frequency signs describe the original real field. Quantum CCR, state,
controls and clock are explicit further hypotheses. The positive canonical
state Hamiltonian stays unchanged. This avoids adding an arc workspace or
using a spectral positive-energy decoder merely to represent that field.
The doubled-sheet coherent walk remains an alternative under different
implementation and clock premises, detailed in section 4.

### Exact local field coordinates for the canonical completion

Factor the full grounded stiffness as `K_D=BᵀCB`, where B has one signed
incidence row for each internal edge and one single-endpoint row for every
Dirichlet boundary stub. C is the diagonal positive conductance matrix. For
constant mass m>0, whiten the original scalar variables:

```text
q=M^(1/2)f,    p=M^(−1/2)Π,    D=C^(1/2)BM^(−1/2),
A=m²I+DᵀD,    H_can=(pᵀp+qᵀAq)/2,    Ω=sqrt(A).
```

The following local first-order coordinates give **exactly the existing
canonical scalar dynamics**:

```text
χ=(u,v)=(m q+i p,Dq),
i χ_dot=H_inc χ,    H_inc=[[mI,Dᵀ],[D,−mI]],
H_inc²=diag(A,m²I+DDᵀ).
```

The edge coordinates are dependent local differences. The invariant physical
image is real-linear, defined by `Im(v)=0` and `v=D Re(u)/m`. Keep both
frequency signs. The local inverse `q=Re(u)/m`, `p=Im(u)` uses no spectral
projection or inverse square root. Classically `||χ||²=2H_can` and
`Im〈χ,χ'〉/m=qᵀp'−pᵀq'`, so the representation preserves the energy and
canonical symplectic form. It does not introduce independent edge oscillators.

With supplied canonical quantum commutators, χ denotes **field operators**,
not a one-particle state or independent annihilators. Entrywise adjoints give
`q=(u+u†)/(2m)` and `p=(u−u†)/(2i)`. On the oscillator core,

```text
[χ_a,χ_b†]=ℏ(H_inc+mI)_ab,
[χ_a,χ_b]=ℏ[[0,Dᵀ],[−D,0]]_ab,
Σ_a χ_a†χ_a=2H_can−m n ℏI.
```

Here n is the number of mutable sites. These identities retain the quantum
ordering correction. The stable state
Hamiltonian is the **original** `H_can`; after subtracting its vacuum energy
it is `ℏ dΓ(Ω)>=0`. Unrestricted bosonic `dΓ(H_inc)` instead has arbitrarily
negative energy and is not this construction. The extra `−m` edge-kernel
modes are excluded by the local field constraints, not populated as additional
particles. The earlier positive-branch map
`[Ω+m;D][2Ω(Ω+m)]^(−1/2)` remains spectral and generally nonlocal; changing
the represented object to constrained field quadratures avoids needing that
map for the canonical field interface.

At q=3 and q=5, exact Q(phi) checks reconstruct every stiffness entry from
incidence and verify weighted symmetry, both coefficient intertwiners,
energy/symplectic forms, local constraints and both commutator matrices.
The q=5 capacities, addresses, preparation and detector agree entry by entry
with the existing scalar execution receipt. Full-matrix checks retain every
frequency sign and edge mode, rather than selecting a comparison band.
They also compare the original force pulse and field detector through this
local representation. Countercontrols retain the extra unconstrained modes,
incorrect independent-boson interpretation, varying-mass obstruction and
alternative inertia. These finite checks support the algebraic join; the
continuum conclusion uses the existing proof.

### Reuse the scalar continuum and instruments without changing their law

The [common scalar packet](../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex)
proves the massive three-dimensional Dirichlet free-field
Weyl-product limit for its stated bounded piecewise-continuous L² smearings
with Riemann-integrable products. The canonical completion above retains
**the same** action, field normalization, vacuum, evolution and preparation.
Consequently its theorem applies directly, with no walk-phase approximation,
population-dependent event normalization or new continuum fit.

The same local f_i couples to the original
[force preparation](../paper/tex_fragments/SOURCE_SCALAR_LOCAL_PREPARATION.tex)
and [finite pointer instrument](../paper/tex_fragments/SOURCE_SCALAR_FINITE_INSTRUMENT.tex).
Substitute `f_i=(u_raw,i+u_raw,i†)/(2m)` with
`u_raw=m f+iM^−1Π`: all coupling coefficients and their local supports remain
unchanged. Existing [sequential-read results](../paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex)
therefore transfer with their original control, entanglement, reset, noise
and Born-readout premises. This is a compatible quantum instrument
construction, not an executed quantum observation history or a
realization of these controls in native C.

Two distinctions matter. The
[regional reconstruction theorem](../paper/tex_fragments/SOURCE_SCALAR_REGIONAL_TIME_SLICE.tex)
requires `row(A_RE) subset row(C_read)` when access consists of two regional
position layers and exterior linear reads; its minimum exterior rank remains
`rank(A_RE)`. Writing a local momentum quadrature does not provide access to
that quadrature from the narrower menu. Also the
[local split-step theorem](../paper/tex_fragments/SOURCE_SCALAR_TIME_REFINEMENT.tex)
retains its step restriction and full-mode vacuum-squeezing bound. Its finite
split flow is not the continuous flow above. A local continuous generator
has finite-cutoff tails; it does not supply a strict microscopic graph cone.

The precise closure is **repair stiffness → conditional canonical law →
local field representation → existing scalar continuum and instruments**.
Quantum CCR, state, controls, the kinetic choice and physical time
need a common microscopic justification or implementation. No new defect
in the cited core theorems was established by this targeted review.

### A finite source code for the scalar force and pointer operations

The existing six-level source controls and coherent proper-code service give
an explicit finite implementation of two local scalar interaction factors.
This reuses the [native source realization](../code/m1_source_realization/README.md)
and [coherent transfer contract](research/COHERENT_SOURCE_CLOCKS.md), rather
than adding a primitive oscillator gate. With supplied canonical fields and
local occupation basis of reference frequency one, write

```text
Q_i=√M_i f_i,       X_d=P_d(a+a†)P_d/√2,       d=2^b,
U_F=exp(+iκ_F X_d), κ_F=√M_i v_i,
U_P=exp(−iκ_P Z⊗X_d), κ_P=√M_i g_i/2.
```

The occupation basis is an encoding choice; its product vacuum does not
replace the coupled scalar vacuum. Encode one oscillator in b proper qubit
codes and retain the pointer. Expand X_d in Pauli strings without dropping
nonzero coefficients. The actual twelve-port controls implement each string
through single-qubit rotations and a phase-exact three-pulse CNOT. A fixed
palindromic product with r repetitions has operator error at most
`|κ|³L³/(3r²)`, where `L=Σ|c_P|`. The bound holds with an arbitrary reference.
Packing, transport, resets, analogue precision and complementary failure
outcomes remain part of the source-service contract.

The comparison reaches the **infinite oscillator**, including amplitude
outside the code. On entrance levels 0,…,k−1, put

```text
t_j=|κ|^j/j! sqrt[2^j(k+j−1)!/(k−1)!],
ε_cut=min(2,2Σ_(j≥d−k+1)t_j).
```

The exact and truncated exponential powers agree below that order, giving
an isometry error at most ε_cut. This also bounds the half diamond distance
restricted to that entrance code, after capping at one. Add the compilation
and service errors. A trace-preserving overflow encoder retains its failure
flag; for occupation tail p_k the general-input error is bounded by
`min(1,2√p_k+ε_cut+ε_compile+ε_service)`.
This mathematical encoder is not a demonstrated vacuum-preparation device.

The bounded checks retain all original 64-site coefficients and select the
largest local force and pointer strengths by exact comparison. All 18 fixed
cases use d=8,16,32 and r=8,16,32, with entrance k=3. They compare complete
operators, entangled-reference outputs and escaped amplitude; phase-erasing
controls fail. At d=8,r=32, the ideal-service restricted half-diamond bounds
are below 0.0022 for the force and 0.0000027 for the pointer, including source
coefficient rounding. The uncancelled constructions require 10,272 and
14,688 native pulses respectively, plus packing/export, routing, reset and
preparation costs. Increasing d at fixed r worsens the synthesis bound;
cutoff and compilation effort must scale together.

The full coupled-vacuum occupations are approximately 2.54 and 2.75 in these
local code coordinates. Conservative analytical occupation bounds give a
**vacuous combined vacuum-output guarantee of one** at k=3. This does not
invalidate the low-input operation bound, but it prevents claiming accurate
control of the original vacuum from these tests. Tiny cutoff errors are not
a replacement for a state-preparation and energy budget.

The completed part is a conditional finite **operation** interface. These
are individual force and pointer factors from disjoint original supports,
not a transported response. Those earlier q=5 component tests do not supply
concurrent free evolution, the complete finite-duration pulse, GHZ pointer
preparation, source vacuum encoding or an executed readout history. The q=3
construction below supplies a separate conditional common-code history; it
does not substitute a q=3 preparation budget for this q=5 target.
Exact CCR cannot hold in a finite code: the boundary defect is
`[X_d,Y_d]=i(I−d|d−1><d−1|)`. Controlled approximation replaces that impossible
identity. Source controls remain supplied; compiling a chosen operation does
not show that autonomous repair selects it or calibrates physical time.

### Local vacuum preparation and one retained-pointer history

The simplest new preparation keeps the original oscillator variables and
turns on the existing stiffness locally:

```text
H(t) = [pᵀp + qᵀ(I + s(t/T)(A−I))q]/2,
s(u) = 3u²−2u³,   0 ≤ u ≤ 1.
```

For the original q=3 model, start from independent frequency-one oscillator
blanks and use the original eight-mode matrix A. The path has frequency gap
at least one. The original integration-by-parts estimate certified full
infinite-oscillator trace distance <0.070155 at T=512. A second endpoint
integration by parts gives **O(T^-2)** for the same path. The argument uses
an independently established uniform intermediate-time bound and exact
Bogoliubov normalization; it does not assume that a small final excitation
implies small excitation throughout the path.

Outward interval enclosures of the exact source eigenvalues and three
nonnegative integrals give:

| Fixed model duration | Certified full-state preparation error |
| --- | --- |
| 512 | <0.000433450 |
| 1024 | <0.000108323 |

All 37 checks passed, with independent review and a byte-identical isolated
replay. The earlier numerical distance 0.00016024 at T=512 remains a
diagnostic. This improves the error certificate without changing the source,
state variables or ramp. The schedule, blanks and model time are supplied;
autonomous vacuum selection and a refinement-uniform resource bound remain
separate questions.

A second constructive preparation encodes squeezed finite states followed by
passive mixing into the existing source-qubit controls. For the same original
q=3 vacuum, its exact retained-mass formula gives trace distance
≤0.059059008; 196 checks passed. Its conservative compiler has an upper bound of about
396 billion native pulses before the declared service costs. The complete
word was not executed or precision-certified. This existence construction
supports the interface but makes the local ramp the preferred next route.

For free evolution, classical force and commuting retained-pointer controls,
the existing finite-pulse solution extends directly to piecewise histories.
Write `H(t)=H_0−j(t)ᵀq`, where j includes the signed pointer Z controls. Then

```text
U(T)=exp(−iH_0T) exp(iΦ) exp[i(aᵀq+bᵀp)],
a=∫cos(Ωt)j(t)dt,   b=∫Ω^−1 sin(Ωt)j(t)dt,
Φ=(1/2)∫_(t>s) j(t)ᵀΩ^−1 sin[Ω(t−s)]j(s) ds dt.
```

The Magnus expansion terminates because those controls commute. Φ includes
identity, Z and ZZ terms; relative pointer phases cannot be discarded.
Five fixed eight-mode histories with two retained pointers and three code
cutoffs passed 100 checks against an independent coherent-state evolution.
Dropping the relative phases changes the two-pointer output, by trace distance
about 0.00298 in the force-plus-two-pointer case. Noncommuting controls lie
outside this exact formula's hypotheses.

The normal-mode finite-code comparison includes escaped output and arbitrary
reference entanglement. It is an operation target, not a free native basis
change. Joining the results uses the existing contract
`ε_total ≤ ε_vac + ε_history + ε_service`, with a compatible encoding,
decoder and error domain. The sharper continuous-preparation bound is smaller
than the approximately 0.003 phase diagnostic, but that alone is not a certified discrimination:
the two complete error balls must be disjoint, and the diagnostic signal
itself needs an error enclosure. Site-to-mode access belongs to this older
normal-mode route; the onsite construction below avoids it. Source control
precision, drift during service, reset and all retained outcomes require
one implementation budget.

### Local time steps for the same preparation

Let `A_k=I+s((k+1/2)/N)(A−I)` and `h=T/N`. Use the local schedule

```text
p ← p − (h/2) A_k q
q ← q + h p
p ← p − (h/2) A_k q.
```

The potential kicks factor into commuting onsite and original neighboring
quadratic terms; the drift is onsite. The proof separates the changing
midpoint stiffness from the frozen local splitting. Instantaneous energy
coordinates control the first term, while the kick–drift–kick modified
energy controls long products without an exponential estimate. For each
mode with final frequency W and b=W²−1, put

```text
r = sqrt(1+3b/(4N)),    c = sqrt(1−(hW)²/4),    hW ≤ 1,
δ = W { b[T/4+rT²/8] + T³W³/(5c) }/N².
```

The antilinear relative Bogoliubov coefficient is at most δ. Summing
`δ_j²/2` over all eight modes, including degeneracies, bounds the squared
trace distance between the two pure Gaussian outputs from the same blank.
This state-specific estimate avoids an operator norm of an unbounded
Hamiltonian. It is not a channel-norm bound on arbitrary unknown inputs.

At the fixed T=512 and N=16,777,216, outward intervals give error
**<0.000695044 relative to the continuous ramp**. Adding the independently
certified endpoint error gives **<0.001129 to the original coupled vacuum**.
All 62 checks and the sole isolated replay passed, including an independent
onsite assembly, constant-frequency controls and the original weaker
splitting bounds. The large step count was evaluated as small symplectic
matrix products, not as a finite source pulse word. Floating ODE/trajectory
comparisons remain diagnostics; the analytical interval estimate supplies
the certificate.

This supplies the changing-stiffness estimate for the scheduled oscillator
preparation. Its finite onsite-code comparison is given below. Source pulse
precision, transport, reset, drift during service and all retained history
outcomes require their own contributions to the total error and resource
budget. Model-time step size does not select physical time.

### Direct Gaussian simulation of the fixed scalar benchmark

Use onsite coordinates with one complex symmetric 8×8 width matrix W and
branch coefficients l,c. Potential kicks update W and l locally; free kinetic
drifts solve with `I+i h W`. All twenty branches of the five force/pointer
histories share W. The common Gaussian normalization phase can be removed
uniformly while retaining every relative branch phase and readout outcome.
This compact representation simulates the infinite-oscillator quadratic law;
it does not assume false finite-cutoff canonical commutators.

The fixed test uses the original q=3 source, frequency-one product blank,
cubic preparation ramp with T=512 and 65,536 steps, then 1,280 steps per
half-unit history. Independent real covariance/mean/Weyl propagation checks
the onsite calculation, and the original continuous histories provide the
comparison. Full joint fidelity uses one coherent sum over pointer branches,
without independent phase alignment.

The preparation distance was numerically 0.000492825, and the largest
prepared joint-history distance was 0.000492026, below the prospective
0.0011 criterion. The largest pointer-state distance was below 8.97×10^-7.
All 90 checks passed; the sole isolated replay was byte-identical. Production
and replay together used 2.895 CPU seconds rounded up on the test machine.
These are floating-point comparisons, not interval-certified errors. The
phase-omission controls remain detectable and all four outcomes are retained.

This makes the component benchmark inexpensive to investigate without changing
the candidate law. It supplies neither a microscopic selection argument nor a
physical clock, and does not execute the finite source pulse word. Keep its
simulation results separate from the conservative construction bounds below.

### Finite onsite preparation with intermediate occupation control

Keep the same action, blank, T=512 and N=2^24. At each site use the number
levels `0,...,d−1`, with **d=8192**, and set `Q_d=Π_d q Π_d`,
`P_d=Π_d p Π_d`. The finite factors are `Q_d²`, `P_d²` and `Q_di Q_dj`.
All potential factors commute, so the existing local schedule is retained.
The frequency-one boundary identity is essential:

```text
Π_d q² Π_d = Q_d² + (d/2)|d−1><d−1|.
```

Replacing the chosen polynomial factor by compressed q² changes the word
and can break that commutation. Tests explicitly detect this substitution.
Each finite generator is extended by identity on spectator spaces for the
comparison; escaped amplitude on the active sites is included.

The proof supplies a covariance bound `V<65I` through **every partial ideal
quadratic factor**, including its ordering within each kick and drift.
For a one- or two-site marginal, the Gaussian occupation generating function
is `M(z)=det[I−(z−1)(V−I/2)]^(-1/2)`. The standard formula is taken from
[Cardin and Quesada, equations (28)–(31)](https://arxiv.org/html/2212.06067v2#S3.E30);
its positive-z convergence domain and the required weighted tails are
checked explicitly. The marginal may be mixed and correlated.

With z=130/129 and threshold d−2, weighted onsite and pair tail norms are
less than `280 exp(−30)` and `700 exp(−30)`. The quadratic defect constants
are 2 and 1. Duhamel's formula and unitary composition along the ideal
prefixes give

```text
D_cutoff ≤ 2 u_d W_onsite + v_d W_pair,
W_onsite ≤ 264192,    W_pair ≤ 262144,
D_cutoff < 0.000033144832.
```

The weights sum integrated coefficients over the full duration; there is
no extra step-count multiplier. Adding the independently checked endpoint
and time-discretization bounds gives **<0.001162** to the original coupled
vacuum. This is an input-specific full-state bound from the supplied blank,
not an arbitrary-input oscillator channel bound.

All **169 checks** passed, followed by an identical isolated replay.
Exact rational inequalities and the reviewed analytic argument supply the
certificate. Independent small-code controls check complete escaped outputs,
number rotations, Pauli reconstruction, identity phases, correlated Gaussian
normalization and the cutoff boundary correction. No 2^104-dimensional
state or native preparation word was executed.

This completes the preparation cutoff comparison on a fixed onsite code.
Native Pauli synthesis uses its separate finite-code operator-norm bound;
its intermediate states need not be Gaussian. The force and retained-pointer
history below uses this same code with separately proved time, displacement
and cutoff bounds. Complete source services have their own error budget.
Initial preparation error propagates by channel
contraction; small trace distance alone does not control a noisy state's
occupation moments. No autonomous vacuum selection, refinement-uniform
resource bound or physical clock is supplied by this finite construction.

### The same onsite code through force, evolution and retained pointers

The original five q=3 histories admit a local construction on the preparation
code, with **d=8192 throughout**. Keep active free evolution, the original
force and two nonoverlapping pointer windows, all relative phases and both
retained readouts. Use history step `h=1/(10·2^16)` for 327,680 steps over
model duration 1/2. No physical normal-mode transform is inserted.

For a constant pointer branch j, the local kick–drift–kick unitary equals
an exact shifted quadratic evolution. With `a=A^-1 j`, each mode obeys

```text
θ=2 asin(hω/2),  κ=ω sqrt(1−h²ω²/4),
B=θ/(hκ),       C=θκ/h,
H_modified = [pᵀBp+(q−a)ᵀC(q−a)]/2 − jᵀA^-1j/2.
```

The Mehler kernel fixes the unitary phase as well as the canonical map.
The scalar term and ordered cross-window Weyl phases are retained. On
these unchanged profiles, a full-state Duhamel estimate gives temporal
error **<0.000005**. Partial-factor covariances and means are bounded by
`V<7I` and `|mean|<6`. Displaced Gaussian weighted tails then give history
cutoff error **<10^-81** on the same large preparation code.

The retained pointer Z sectors are orthogonal. Bounding their field defects
branch by branch and taking the maximum preserves the coherent output;
it does not dephase or postselect the pointers. Arbitrary reference
entanglement with those pointers is retained. The final complete readout
contracts the same trace-distance bound. The initial field is the declared
coupled vacuum; approximation of that field state enters separately.

All **97 checks** and a byte-identical isolated replay passed. Exact source
and rational checks support the analytic bounds. Independent complex
Schrödinger wavefunction evolution agrees with the complete Weyl-history
pointer matrices, including their phases. Small-code controls retain escaped
outputs and the top-occupation correction. Preparation plus this ideal
finite history has error **<0.001167**, before native synthesis and services.
The phase-omission separation near 0.002979 remains a numerical diagnostic;
its rigorous separation and complete implementation error require their own
comparison.

### Compiling the finite preparation with existing source controls

The existing six-level response, proper qubit codes, phase-exact CNOT and
Pauli gadgets implement the finite quadratic preparation to a prescribed
accuracy. A binary-carry decomposition bounds its Pauli coefficient sums
without generating a huge list: `L(Q_d)<(1+sqrt(2))d/2`. Products and the
exact number rotation give `L≤3d²/2` for all three quadratic generators.
The unchanged preparation coefficients satisfy `Σ|θ|³≤129/128`.

Using the existing symmetric product formula with fixed `r=2^46` gives

```text
ε_synthesis ≤ L³ Σ|θ|³/(3r²) ≤ 1161/2^24 < 0.0001.
```

This bound is uniform on the finite code, including references. Identity
Pauli phases are supplied by the stored `phase_0+phase_1` control and charged
once per onsite factor. The native intermediate states need not be Gaussian;
the occupation comparison follows ideal quadratic substeps, and the compiler
comparison covers the whole finite code.

All **177 checks** passed with an identical isolated replay. They include
independent binary-carry and Pauli reconstructions, actual full-six-level
control identities, complete proper-code failure channels and phase controls.
The large construction was neither emitted nor executed. Its sufficient
uncancelled upper count is about **3.27×10^34 native pulses**, plus twice as
many code-service events. This enormous upper bound is not a lower bound on
other compilers or evidence of efficient implementation.

Analogue precision, transport, storage, reset, physical processor drift,
source clocks and complete preparation/readout services remain charged
inputs. A matrix synthesis theorem does not establish their physical
availability or identify source execution time with the model-time ramp.

### Exact commuting compilation: simplify the implementation first

The same finite targets admit an algebraically exact factorization. Write
`Q_d=O Λ Oᵀ` and retain the existing `P_d=F Q_d F†`, with `F=diag(i^n)`.
Then Q, Q², QQ and pointer-ZQ operations are diagonal phase operations in a
common onsite quadrature basis; P² uses the same basis with the number phase.
An entire potential layer can share its charged entrance and exit basis
changes. No inner noncommuting product-formula repetition is required.

The basis O is itself compiled. Gray-neighbor Givens rotations reduce to
controlled-Y rotations; expanding the spectator projector gives mutually
commuting Pauli terms. Diagonal phases use commuting Walsh-Z terms, including
the identity phase. These words use the same twelve-port source controls and
proper single/pair codes. The eigensolver is a classical compiler step, not a
new native quantum gate or a free change of chart.

All **45 fixed cases**, with d=2,4,8 and angles 0, 1/7, −1/5, passed against
independent full matrix exponentials. The source pulse words were emitted and
numerically evaluated from their full six-level generators. The largest
full-phase operator residual was **2.181×10^-13**; the sole isolated replay
reproduced results, checks and emitted words byte for byte. At d=8 and angle
1/7, the QQ operation uses 4,870 native pulses, including its basis changes,
plus 9,740 pack/export events under the declared service. This is a component
cost at that cutoff, not the complete preparation/history cost.

The controls detect wrong Gray orientation, omitted spectator restrictions,
lost phases and the wrong top-level quadratic compression. In particular,
finite cutoff dynamics are not exact oscillator shears: Q₂²=I/2 makes that
purported identity fail. Cheap Gaussian simulation and finite source-word
execution therefore have separate verification contracts.

This removes the generic compiler's inner Strang error in ideal exact
arithmetic. Floating angle/control errors, full service costs and practical
large-d execution need budgets; the elementary basis construction has
an O((log d)d³) upper count. No large-cutoff word or physical hardware was run.
The following full-history bound retains its original compiler and error
budget until a replacement has been checked at that scope. This simplification
is evidence about implementation, not selection of a microscopic source law.

### One composed preparation–history error and resource budget

Use the same **106 work qubits**: 104 for the eight onsite field codes and
two retained pointers. The original blank, ramp, force, pointer profiles and
readout remain fixed. No physical normal-mode conversion or new primitive
matrix gate is added. For the history, the joint quadratic angle sum is at
most 258 and the linear sum below 15. Thus `Σ|θ|≤273` and
`max |θ|≤1/10240`. The same compiler repetition count `r=2^46` gives

```text
ε_history synthesis ≤ L³ Σ|θ|³/(3r²)
                    ≤ 2457/(25·2^39) < 2×10^-10.
```

Multiplying a local Q by a retained-pointer Z preserves its Pauli
coefficient sum. This estimate is uniform over the finite code and untouched
references; it requires no Gaussian assumption about intermediate native gates.
The resulting mathematical comparison with the original continuous history
has the following budget. Displayed component decimals are rounded upward;
the exact rational sum is below 0.001235681.

| Contribution | Sufficient upper bound |
| --- | ---: |
| Preparation endpoint, local time steps and onsite cutoff | 0.001161638 |
| Preparation source synthesis | 0.000069201 |
| History time steps, including pointer phases | 0.000004843 |
| History onsite cutoff | 10^-81 |
| History source synthesis | 2×10^-10 |
| **Combined mathematical error, before physical service errors** | **<0.001237** |

The order of comparison matters. First apply the uniform finite-history
compiler bound to the actual encoded prepared state. Under the same extended
finite history, replace that state by the ideal coupled vacuum and charge
its preparation and synthesis errors. Only then use Gaussian cutoff and time
bounds on the ideal history branches. Finally apply the same complete CPTP
readout. This neither assumes that trace-distance-close states have small
occupation moments nor discards pointer coherence or failed outcomes.

The join passed **92 exact arithmetic, input-custody and accounting checks**
and an identical isolated replay. These checks do not formally prove the
analytic composition; that argument was reviewed separately. Together with
the three component packets, the total is 535 checks and four replays.

The combined evolution recipe has a conservative upper inventory of
`34,316,088,908,833,975,917,332,078,914,437,120` native pulses, below 2^116,
and at most twice as many pack/export events, below 2^117. The 106 work
qubits exclude the processor, fresh service buffers, references and records.
Initialization, routing, reset, storage and readout need additional resources.
Neither the Pauli lists nor this native word were emitted or executed.

For an actual source execution, the error budget needs

```text
ε_actual ≤ ε_mathematical + ε_field blank + ε_pointer initialization
           + Σ ε_pulse + Σ ε_code service + ε_storage/drift
           + ε_readout + other admitted service errors.
```

Service bounds must retain references and all outcomes, for example through
half diamond-norm estimates. No values or attainability for these terms are
supplied here. **RE-10 is conditional, not closed**: its common-code
mathematical composition is constructed. Native implementation and
physical attachment require the additional inputs listed above. No separating physical measurement,
continuum-uniform implementation or selected physical clock follows.

### A short selection theorem isolates the substantive extra assumption

Let `K=K_D+m²M>0`. Suppose `(f,ν)` evolves linearly and autonomously,
`f_dot=ν`, its conserved energy is
`E=(fᵀKf+νᵀMν)/2`, and time reversal sends `(f,ν)` to `(f,−ν)`.
Writing `ν_dot=Ff+Gν`, conservation for every state forces
`MF=−K` and `MG+(MG)ᵀ=0`. Time reversal forces G=0. Thus

```text
f_dot=ν,    M ν_dot=−Kf,
Π=Mν,      H_can=(ΠᵀM^−1Π+fᵀKf)/2.
```

This proves the unique dynamics **within those hypotheses** and exhibits
its local canonical form and Legendre action. It does not select the
hypotheses, an arbitrary nonlinear Poisson structure, quantum CCR or physical
time. In particular, identifying repair capacity with inertia is a candidate
principle, not a consequence of the equilibrium alone. Linearity also
matters: reversible nonlinear forces perpendicular to velocity can conserve
the same energy.

A concrete counterfamily shows the missing information. In whitened
coordinates, `H_ε=(pᵀ(I+εA)p+qᵀAq)/2`, ε≥0, is positive, local and
reversible and preserves every orthogonal source symmetry commuting with A.
Its normalized **classical configuration Gibbs marginal** is unchanged,
but its frequencies are `sqrt(λ(1+ελ))`. For distinct eigenvalues these changes
cannot be absorbed into one clock rescaling. Its full phase-space Gibbs
reference and quantum vacuum change. Configuration agreement and the
restoring force therefore do not fix the inertia; dynamical evidence can
distinguish these laws.

### Distinguishing kinetic laws after calibration

Put `A=m²I+L`, where L is the same source spatial stiffness. The kinetic
alternative has the exact acceleration operator

```text
Y_ε=A+εA²
   =m²(1+εm²)I+(1+2εm²)L+εL².
```

A common clock rescaling within this mode comparison multiplies Y by γ>0.
Refitting a canonical mass and speed gives an affine function of a spatial
eigenvalue μ; the alternative
has curvature γε. For three equally spaced source modes with spacing d,

```text
C = (Y₂−2Y₁+Y₀)/(Y₂−Y₀)
  = εd/[1+2ε(m²+μ₁)],
```

provided the denominator is nonzero. C is independent of the common clock.
At graph distance two, `I_ij=L_ij=0` but `(L²)_ij>0`, so
`Y_ij=ε(L²)_ij` independently proves the same obstruction to an affine refit.
The exact raw golden-source coefficients are `40+64φ` at q=3 and
`178+288φ` at q=5. Positive mass whitening preserves their signs and zeros.

The fixed tests used ε=0, ±1/65536 and −1/16 with two clock multipliers.
All 16 cases passed 258 checks and an identical isolated replay. The two
small perturbations are positive-kinetic; −1/16 is rejected by exact
instability bounds. With unit clock and ε=1/65536, fitting the first two
modes leaves third squared-frequency residuals about 0.00928 and 0.02126 at
q=3 and q=5, respectively. Both fitted masses remain positive. These are
internal source comparisons with fixed spatial/readout identities, not
laboratory measurements or a selected physical ε.

Two useful simplifications prevent redundant searches. First,
`I+εA=(1+εm²)(I+ε'L)`, with `ε'=ε/(1+εm²)`. The corresponding canonical
field rescaling and clock change identify these parameterizations only when
force and readout calibrations are changed consistently. Second, a fixed
positive ε retains a `k⁴` dispersion term, and a fixed negative ε eventually
loses kinetic positivity under refinement. Neither disappears by taking a
finer mesh.

In contrast, ε_q≥0 tending to zero can share the existing **position-field
Weyl and force-to-position limit**, even without ε_q||A_q|| tending to zero.
Reuse the scalar theorem's full fixed-band stiffness residual and Riemann-tail
argument: for `g_ε(x)=sqrt(x+εx²)`,
`g_ε'(x)²=ε+1/[4x(1+εx)]`; the vacuum field coefficient is
`(A_q^-1+ε_q I)^(1/4)`, uniformly close to `A_q^-1/4`.
Negative vanishing ε_q additionally needs a uniform positive kinetic margin.
This is an analytical, observable-specific result. It does not transfer the
full velocity menu, equate finite vacua or provide a native continuum
implementation. Its consequence for the search is precise: a shared free
continuum is insufficient to identify unique microscopic inertia.

This is consistent with the existing
[source-to-Hamiltonian theorem](../Lean/Variational/SourceToHamiltonianComposed.lean),
which uses real enrichment, regularity and real fixed-endpoint minimality,
and the [realized-history nonidentifiability result](../Lean/Variational/RealizedHistoryLegendreNoGo.lean).
The [field energy-inner-product](../Lean/Screen/FieldSectorEnergyInnerProduct.lean)
and [clock conversion](../Lean/Screen/CurlStoneClockBridge.lean) results
supply conditional Hamiltonian readings of their committed field
step. These conversions need not be reproved. The next question is why the
source chooses that kinetic step and clock, together with the original
capacities, mass and conductances.

## 4. Retained native results and coherent alternatives

These constructions remain valid under their stated inputs. They are
alternatives or native baselines, not additional steps required by the
canonical scalar completion above.

### Native working state and atomic update retained as the baseline

Use the documented twelve-port icosahedral carrier with thirty seams. The
current exact test sector is

```text
x_i ∈ nonnegative integers,   Σ_i x_i = 2,
H_work = span{|x〉}: dimension 78.
```

These are configurations of the entire carrier, not the twelve primitive A1
atoms. Total two is a declared test sector, not a selected cosmic state. The
binary shell `x_i∈{0,1}` has 66 configurations and is invariant under repair.

Orient a selected seam e=(i,j), set `d=x_i−x_j`, and take a fair placement
coin `c∈{0,1}`. With `b_e=unit_i−unit_j`, define

```text
j(x,c) = [d − (2c−1)(abs(d) mod 2)]/2,
x' = x − b_e j(x,c).
```

Equivalently, replace the two endpoint loads by the floor and ceiling of
half their sum; c=1 places the ceiling at i. This is the native nearest-integer
law. It conserves total load. Large imbalances lower `V=Σ_i x_i²`; difference
one gives a swap or wait, and difference zero gives a wait. Some neighboring
seam disagreements can increase even when V decreases.

**Atomic instruction:** apply that local reconciliation and retain its
transported current with the coin and event identity. The quantum candidate
coherently groups configurations with the same `(c,j)`:

```text
L_(e,c,j) = (1/√2) Σ_{x:j(x,c)=j} |x−b_e j〉〈x|,
C_e(ρ) = Σ_(c,j) L_(e,c,j) ρ L_(e,c,j)*,
Σ_(c,j) L* L = I.
```

Each branch is an injective translation on its domain. The completeness
identity makes C_e completely positive and trace preserving (CPTP), with
exactly the native classical transition probabilities. The Hilbert-space
representation, quantum probabilities and record registers are explicit
candidate structure. The classical law alone does not uniquely select them.

### Recorded events, observers and finite hardware

Keeping the coherent record implements the isometry
`W_e=Σ_(c,j) |c,j〉⊗L_(e,c,j)`. Knowing the output and current recovers the
classical predecessor, `x=x'+b_e j`. Keeping the full coherent joint state
retains more information than measuring those labels.

For total two there are six `(c,j)` letters. The fixed-capacity recorder is

```text
H_cap = ⊕_(p=0)^L [(C^6)^⊗p ⊗ H_work],
T_C(ρ) = Σ_(p<L) V_p P_p ρ P_p V_p* + P_L ρ P_L.
```

The pointer p selects the programmed seam; V_p appends its record and applies
C coherently. All sectors are allocated initially. Cross-pointer coherence
is dephased, and the full-capacity sector halts. At L=2, dimension is 3354;
the tested word uses (0,1) then (0,5). A separately declared seven-state,
single-use meter supplies actual working-number and record measurements.
The total space including that meter has dimension 23478.

Proper patch restrictions trace unowned ports and records whose entire seam
is not owned, retaining event-slot addresses. Their overlap restrictions
agree. Later repair and an outcome-averaged working measurement preserve the
whole old record marginal. Selecting a measurement outcome can change its
conditional state. Reading the record itself can destroy its quantum
coherence while preserving repeatable classical values.

A proof for arbitrary finite words requires the corresponding capacity and
schedule. The two-write recorder does not provide indefinitely active repair
or a free reusable memory. Its explicit halt is a completed finite resource
policy; sustained operation requires a further justified policy.

A constructed alternative uses a fixed 180-state temporary register Z for
`(seam,coin,current)`. With uniform independent seam choices, let

```text
J(ρ) = (1/30) Σ_e |e〉〈e| ⊗ W_e ρ W_e*,
F(σ_ZW) = J(Tr_Z σ_ZW),
F^n = J Q_C^(n−1) Tr_Z,   n≥1.
```

This CPTP law runs on the same 14,040-dimensional working/diagnostic space
at every tick. Its unique fixed state is `Ω=J(τ_B)`; the working reduction is
τ_B. A separate protected register retains its marginal, with joint limit
`ρ_R⊗Ω`. The prior diagnostic is explicitly discarded. There are `60^n`
distinct seam/coin histories after n active attempts, so indefinite exact
readability of all history from the record alone requires that many
distinguishable record states. Fixed finite memory cannot provide it. This is an open
channel implementation contract, not a closed-unitary Universe or a supplied
physical reset environment. Neither policy changes native C itself.

### Reading correlations from two repair events

For each occupied pair {i,k}, choose native seams (i,j) and (k,l) with four
distinct endpoints. Such seams exist for every pair. The two recorded
outward transfers, with their specified coins, can occur only from {i,k}.
Their joint effect on the initial shell state is

```text
effect = |ik〉〈ik|/4       conditioned on the ordered seam choices,
effect = |ik〉〈ik|/3600    including the independent uniform seam schedule.
```

Consequently the 66 query probabilities determine all diagonal pair
probabilities. They determine regional number distributions, exterior
completion counts and the entropy differential at τ_B. No added field or
different microscopic operation is needed. This exact identification uses
an ensemble of the same preparation; finite-sample accuracy and sample cost
are additional questions. Both event labels must be retained. A single-event
scratch register alone does not implement this observation menu.

One event provides less information: its full label menu spans the constant,
port loads and native-edge pair occupations, with rank 42 and a
24-dimensional invisible diagonal subspace. The original entropy/activity
counterexample below is visible to that full one-event menu. A
different tangent, `(δ_02+δ_13−δ_03−δ_12)/66`, is invisible to every
one-event label but changes entropy on the connected region {0,2,10} by
`log(9/4)/66`; two events detect it.

Classical record labels do not determine quantum coherences. Every finite
record word is an injective partial translation, so its probability effect
is diagonal. Coherent and dephased preparations can therefore give identical
label statistics and different regional entropies. The constructed
working interventions remain necessary for quantum claims. At τ_B, the
regional modular observable is diagonal, so the first entropy differential
is nevertheless accessible through this classical reconstruction.

### Coherent current histories give full ensemble tomography

Keeping the histories coherent extends what can be read without changing C.
For any two shell configurations x and y, choose a native exclusion path from
y to x. All pairs have such a path of length at most four. One record history
holds x fixed; the other moves y to x along that path. The histories are
orthogonal, have unique input preimages and end at the same configuration.

Measure that final occupation and interfere the two record histories with
relative phases `1, −1, −i, +i`. For a path of m events the effects are

```text
F_t = 2^(−m−1) [|x〉〈x|+|y〉〈y|+conj(t)|x〉〈y|+t|y〉〈x|],
Re ρ_xy = 2^(m−1)(p_1−p_−1),
Im ρ_xy = 2^(m−1)(p_−i−p_i).
```

Together with diagonal readout these span all 4,356 Hermitian directions
on the 66-state shell. The construction uses 4,291 unoptimized measurement
settings. Four events suffice; no minimal-depth claim is made. With only
real interference the imaginary coherences remain invisible; dephasing the
current record before interference destroys the relevant signal.

The additional resources are explicit: repeated identical preparations,
full working occupation readout, a collective coherent-history measurement
and a complex phase reference. This is ensemble tomography, not a way to
copy or fully read an unknown single quantum state. At the uniform shell
reference a length-four two-outcome fringe occurs with probability 1/528
under prescribed seams, or 1/427,680,000 including its passive uniform seam
word. The four-event six-letter recorder has dimension 102,630 on this
shell. Informational completeness does not imply economical measurement.

### An independently defined candidate agreement target

On the binary shell, define S_e by exchanging the seam's endpoint occupancies
and `U_e=exp[iπ(N_i−N_j)/2]`. They are local operators lifted to the common
total-two space, independently of the repair channel. The proposed public
algebra is the equalizer of `Id`, `Ad(S_e)` and `Ad(U_e)`. The eight-element
group generated by S_e and U_e has an exact twirl equal to E_e. Its algebra is
`M_46 ⊕ (I_2⊗M_10)`, of dimension 2216. The corresponding common equalizers
have dimension 1378 for the tested two seams and one for all thirty seams.
With protected R, the common target includes the full protected algebra
`B(R)⊗I_binary`.

This closes the candidate target-to-completion calculation. The public
interpretation explicitly ignores those representation changes. A correctly
translated covariant chart can instead retain the full working algebra, so
A2 naturality alone does not select this target. Endpoint number remains a
changing working probe, not an invariant committed meaning by definition.
Completion preserves the chosen retained data; observers with different
retained data are not made equal merely by applying this twirl.
Overlapping primitive completions need not commute; a common ensemble limit
does not prove transaction-level confluence of individual histories.

### Atomic steps, completed repair and network equilibrium

The completed comparison law E_e records `(coin, pre-repair contrast)` rather
than `(coin,current)`. It has the same classical kernel but different quantum
coherences. On one seam,

```text
E_e²=E_e,   E_e C_e=C_e E_e=E_e,
C_e^n=E_e+2^(1−n)(C_e−E_e),   n≥1.
```

Thus C is a partial reconciliation whose repeated use reaches E. On the
binary shell E is a faithful trace-preserving conditional expectation;
its fixed algebra equals the candidate invariant target above, while
identification with the complete source A2 overlap diagram is open.
The corresponding contrast recorder remains a valid comparison construction,
not an abandoned model or the same quantum instrument as C.

For independent uniform choices among all thirty seams, let
`Q_C=Σ_e C_e/30`, `Q_E=Σ_e E_e/30` and `τ_B=I_binary/66`. Both channels obey

```text
Q_C^n(ρ) → Tr(ρ) τ_B,
Q_E^n(ρ) → Tr(ρ) τ_B.
```

The result includes all 78 input states and transient piles. On the 66-state binary shell,
Hilbert–Schmidt contraction is bounded by `(3239/3250)^n`; the bound is
certified but nonsharp. The uniform two-seam mixture retains a larger
outside algebra and has recurrent contraction factor 3/4. The single-seam
factor 1/2 is not a network rate. Positive fixed weights and bounded-coverage
reduced schedules have their stated extensions; unrestricted adaptive,
outcome-conditioned trajectories are not covered by the ensemble claim.

An untouched protected register R gives the all-thirty reduced limit
`ρ_RW → ρ_R⊗τ_B`, including initially entangled states. Protected distinctions
survive while working correlations relax. The written coherent records
retain information outside that reduced state. This is an explicit protected
hardware construction, not a derivation of protection from bare working loads.

The reference and continuation matter: uniform binary equilibrium can coexist
with moving swap/wait realizations. An exactly all-equal native state is
absorbing. This law does not generate the nonzero-remainder sector or the
first physical disturbance from perfect symmetry. State selection and the
source of incoming interactions remain separate questions.

### The operator follows from the native binary kernel

Let S_e swap the endpoint occupations. The native fair coin gives
`P_e=(I+S_e)/2`. For uniform seam attempts on one carrier, define P as the
classical transition matrix on the 66 configurations. A unit-rate Poisson
embedding gives

```text
P = (1/60) Σ_e(I+S_e),       Q = P−I,
H = −Q = (1/60) Σ_e(I−S_e) ≥ 0,
H = (1/120) Σ_e[I−σ_i·σ_j]   on the full binary-qubit space.
```

The last expression is a ferromagnetic exchange Hamiltonian; its qubits are
occupation coordinates, not identified physical spins. The connected
fixed-total shell has unique uniform ground amplitude `g_x=1/√66`. The
full qubit space has multiple ground states across occupation sectors.
No independent exchange strength was fitted, once the sector, fair coin,
uniform schedule and Poisson normalization were declared.

For diagonal observables F and G, stationary repair has the exact relation

```text
〈g|F exp(−sH) G|g〉 = E_τ[f(X_s)g(X_0)],       s≥0.
```

The analogous ordered multi-time identities also hold. This is the established
[Rokhsar–Kivelson correspondence](https://arxiv.org/abs/cond-mat/0311345),
with [stochastic-matrix-form generalizations](https://arxiv.org/abs/cond-mat/0502068),
instantiated from the actual source events. It does not equate repair
trajectories to `exp(−itH)` evolution. The density-channel generator
`Q_C−Id` acts on 4,356 operator directions; H acts on 66 amplitudes.
The native mixed reference `I/66` is not the pure vector state `|g〉〈g|`.

The distinction is testable at finite size. A neighboring transition
has native probability slope 1/60, whereas the unitary transition probability
starts quadratically, with second derivative 1/1800. No common positive
constant time rescaling equates them. The specified reduced C semigroup has
no nonconstant undamped sector. This does not exclude coherent environments,
controlled records or a different limiting physical interpretation.

### The same C can reveal reversible conditional gates

On local basis `(00,01,10,11)`, write its four active Kraus operators as

```text
a = diag(1,1,0,1)/√2,    b = |01〉〈10|/√2,
c = diag(1,0,1,1)/√2,    d = |10〉〈01|/√2.
```

For the four sign triples `(+++)`, `(+−−)`, `(−+−)`, `(−−+)`, let U_s be
identity on 00 and 11 and, on 01/10, let

```text
U_s = [I+i(s_x σ_x+s_y σ_y+s_z σ_z)]/2,
C_e(ρ) = (1/4) Σ_s U_s ρ U_s*.
```

A specified unitary change of basis on the existing coherent current record
gives these four Kraus operators `U_s/2`. Every outcome has probability 1/4
independently of the input, and its conditional working operation is unitary.
Known records plus supplied inverse control reverse even overlapping events.
This is an exact unravelling of the same channel, not a changed classical
law or a claim that unobserved repair is reversible.

Sharp current readout and this tetrahedral basis are incompatible measurement
menus. Measuring current first destroys the unitary conditional property.
Their unchanged discarded channel does not make them simultaneously available
sharp records. The tetrahedral gate logarithms average to zero, so this
unravelling does not select H. It is valid on binary occupations; the full
78-state C is nonunital and cannot be an average of unitary channels.

### Adaptive record readout controls the full finite shell

The four-gate menu need not stay fixed. Rotating its tetrahedral axes by any
supplied `R∈SO(3)` is another orthonormal measurement basis of the same four
active current records. Each outcome applies identity on even occupations
and `exp[i(π/3)n·σ]` on odd occupations, with probability 1/4. Every branch
is unitary, including outcomes other than the designated one. At each fixed
scheduled step, discarding its outcome gives C_e. For a prescribed seam word
and count, averaging all outcomes gives the original channel composition.
Adaptive stopping and seam selection define a different controlled protocol.

Any odd-sector SU(2) residual target is a product of three designated gates
of that form. After each three-event block, track the actual accumulated
gate and choose a new residual block. A designated triple has conditional
probability 1/64, independently of the unknown input and previous outcomes:

```text
Pr(not completed after B blocks) ≤ (63/64)^B,
E[native events per local odd-SU(2) target] ≤ 192.
```

No working trial is discarded. These are ideal exact-control bounds, using
continuous record-basis settings, phase references, gate tracking, fresh
records and adaptive stopping. They do not establish finite-precision
hardware efficiency; a finite recorder has a failure tail. Retaining the
accumulated gate is not retaining every historical record.

One seam alone cannot implement arbitrary exchange phases while preserving
all even/odd coherences. Native overlapping seams overcome that restricted
obstruction. Let `X_ij=σ_i^+σ_j^-+σ_j^+σ_i^-` and let N_i be occupation.
The native triangle `(0,1,5)` supplies

```text
K = i[i[X_ij,X_jk],X_ik] = −2(I−2N_j)(N_i−N_k),
[K+2(N_i−N_k)]/4 = N_j(N_i−N_k).
```

Connected native iSWAP controls permute all number projectors; their phases
cancel under this conjugation. They therefore supply every pair-projector
difference `N_a N_b−N_b N_c`. On the total-two shell, these span all 65
traceless diagonal directions. Double commutators with the connected native
exclusion transitions isolate the off-diagonal matrix directions. The
resulting control Lie algebra is **su(66)**, established structurally rather
than by a numerical rank estimate.

Since `Tr H=10`, the traceless target is `H−(5/33)I`. Its unitary evolution
is controllable on this fixed shell; the scalar difference is a global phase
and gives the same density evolution as H. This is a full finite join of
the actual C records to the particular source-derived Hamiltonian, under
the stated controller resources. It does not control relative phases between
charge sectors or establish the same result on the nonunital 78-state channel.

The 192-event expectation is **per local target gate**, not a bound for a
complete H synthesis. The Lie-algebra result supplies no practical global
word length, finite-precision error budget or continuum resource estimate.
It removes the trial-discard cost of the simpler protocol below; it does not
select a preferred controller, Hamiltonian program or autonomous physical
clock. The requested evolution parameter t is a control input, and random
completion time is a separate observable.

### Two native events implement the particular exchange gate conditionally

The join can be made to H itself. Write xy for matrix product K_x K_y,
whose chronological record is `(y,x)`. The actual two-event words obey

```text
I = 2(aa+cc−ac),        S_e = 2(ac+bc+da).
```

Consequently, a specified coherent projection on the two records implements
`U_θ=exp[−iθ(I−S_e)]`. One five-word projection has success probability
1/12 for every input. The least-norm projection using all ten nonzero words
has

```text
p_θ = 3/[4(4−cos²θ)],       3/16 ≤ p_θ ≤ 1/4.
```

This optimum is for one rank-one herald outcome in exactly two native events,
with no extra working control. It does not bound feedback, longer records or
multiple accepted outcomes. The success/failure instrument is trace preserving
and sums back to the unchanged C_e². Failure generally changes the state;
an unknown input cannot simply be retried unchanged. Dephasing records can
leave a plausible herald frequency while destroying the target coherence.

A prescribed sweep of the thirty seam gates at `θ=t/(60r)`, repeated r times,
conditionally approximates the full exchange evolution with

```text
‖U_r(t)−exp(−itH)‖ ≤ t²/(30r).
```

The bound uses the native 120 overlapping seam pairs and commuting disjoint
pairs. There are 30r heralds and 60r primitive events. The all-success
probability is `p_(t/(60r))^(30r)`, at most `4^(−30r)`, or `12^(−30r)` for
the simpler construction. Thus this is a rigorous finite operational join
with exponential cost, not an efficient deterministic implementation or an
autonomous physical law. The phase angle, schedule and coherent measurement
are supplied. Infinite Trotter refinement is not carrier/UV refinement.

The exponential cost belongs to this protocol, not to every possible record
control. An analytical follow-up uses the tetrahedral gates' finite
24-element group: tracking the gate product and stopping at an attainable
target gives almost-sure arrival with finite expected waiting, provided
fresh coherent records and adaptive stopping are available. A finite recorder
has a failure tail. On the complete local binary space these fixed-menu words
are identity on the even sector and have odd-sector determinant one; they
cannot equal a nontrivial U_θ above. On an odd-only sector global phase is
unobservable and some projective swaps are reachable. Changing readout bases,
overlapping seams or feedback may offer broader control and is not excluded
by that restricted observation.

A separate supplied coherent seam-choice register permits a one-event
postselected filter `P/√2`, with success `Tr(P²ρ)/2`. A classical mixture of
seam choices does not supply the required inter-seam interference. This
provides an operational transfer filter, distinct from both H and a unitary.

### A common preparation and observable distinguish the time laws

The port Laplacian has density bands `λ_±=5±√5`. Choose
`v=P_λ unit_0` from its exact spectral projector and define
`F(x)=Σ_i v_i x_i`. The same shell gives

```text
〈F〉_τ = 〈F³〉_τ = 0,     σ²=〈F²〉_τ=5/132,
χ_x = F(x)/(√66 σ),       Hχ=(λ/60)χ,
ψ(0) = (g+χ)/√2.
```

This explicit coherent preparation and the same working-number observable
then predict

```text
unitary H:    〈F〉(t)=σ cos(ωt),
native C:     〈F〉(s)=σ exp(−ωs),       ω=λ/60.
```

A negative unitary half-cycle cannot equal positive native decay under any
forward clock relabeling. The fast/slow frequency ratio is `(3+√5)/2`,
independent of a common rate. Period divided by the matching relaxation time
is `2π` under the declared common normalization; independently calibrated
clocks cannot silently be equated.

Coherent quadratures give phase modulo `2π` through the constructed ensemble
readout. Cycle counting, persistent records and measurement back-action
matter for an individual clock. Neither law creates this preparation from
its stationary mixed reference. This closes a finite operational discriminator,
not physical seconds, energy units or agreement with an observed system.

### Other routes remain explicit alternatives

| Route from repair | What is constructed or available | Additional selection and unresolved cost |
| --- | --- | --- |
| **Poisson exchange: H=I−P** | Exact finite positive local H, Euclidean correlations, adaptive native-record controllability and the number-mode clock test above. | Physical real time, preparation and units. Adaptive control lacks a practical global synthesis bound; the simpler herald protocol has exponential all-success cost. |
| **Discrete transfer: H_T=−log P** | The declared discrete transfer defines this spectral operator on the positive spectral support; coherent P filtering is constructed. | A zero eigenvalue requires support restriction. Logarithms generally introduce longer-range terms and different frequencies; no local physical-time identification is proved. |
| **Retained-record unitaries** | Same-C rotated readouts, ideal local feedback and overlapping-seam su(66) control; two-event heralds give a direct exchange benchmark. | Continuous basis control, phase reference, adaptive schedule and resources. Controllability or an average channel does not select the physical Hamiltonian. |
| **Local coherent walk** | Vertex-local weighted reflection and arc reversal give phases `±arccos(λ_P)`, a strict graph-distance cone, a periodic native-carrier free-field limit and the weighted golden-source scalar connection. | It changes the quantum instrument from C. Arc states, clock/space scaling, positive branch and field interpretation are supplied; the actual physical refinement remains open. |
| **Reference/modular generator** | The selected state defines `K=−log τ`; regional modular operators are derived below. | The global shell reference is `I/66`, so its K is constant and generates no nontrivial global motion. Another state or a thermal-time identification requires an additional argument. |
| **Unitary extension of retained records** | The native W_e is an isometry; extending its action to a unitary on a larger space permits a Hermitian logarithm. | The extension, logarithm branch, clock, record supply and locality must be specified. A chosen dilation is not a uniquely selected physical Hamiltonian. |
| **Incidence/first-order lift** | The actual port-graph incidence gives a local Hermitian block operator squaring to vertex/edge Laplacians; cycle and constant zero modes are counted. | Edge-space and energy-sign interpretation, physical clock and compatible polar/edge maps remain supplied. Continuous-time propagation has tails; a local matrix does not imply a strict discrete light cone. |

In the earlier carrier-wide construction, a fair swap/wait coin supplies
the same P as native placement. Let A append a uniform 60-state seam/coin register, T perform its
controlled swap and `R=2AA*−I`. Then

```text
U_walk=TR,        A*TA=P,
U_walk(Av)=TAv,   U_walk(TAv)=2λ_P TAv−Av.
```

For `|λ_P|<1`, this two-dimensional span has characteristic polynomial
`z²−2λ_P z+1`; at `λ_P=1` the vectors coincide. Thus
`θ=arccos(1−ε)~√(2ε)` near a vanishing relaxation gap ε. This instantiates
[quantized-walk spectral theory](https://arxiv.org/abs/quant-ph/0401053).
If a controlled spatial limit later gave `ε(k)∝|k|²`, this branch could have
linear small-momentum phase. That conditional observation does not establish
Lorentz symmetry or a continuum field.

Discarding that swap register gives `(ρ+S_eρS_e)/2`, which differs from C_e
on coherences. Repreparing it at every step loses the closed walk. Its global
reflection was an implementation gap. The following construction supplies a
local realization on the one-particle port graph using vertex-local arc
states; it does not establish locality of the earlier binary66 walk.

### A fixed local coherent update on the actual graphs

On the archived port-pair graph let L be the port Laplacian, with maximum
degree six. Define `P=I−L/63`. Append an outgoing directed arc at each vertex,
including a self-loop, with amplitude `sqrt(P_yx)`. Let A be that vertex-local
isometry and S reverse each arc. Then

```text
U=S(2AA*−I),       A*SA=P.
```

The reflection acts inside one vertex; reversal acts along one edge. Their
fixed composition is a two-layer local unitary with a strict graph-distance
light cone: after n steps, support moves at most n graph edges. The actual
240- and 960-port examples use 1,500 and 6,000 arc states. No global seam
reflection or adaptive gate controller is needed for this candidate update.

The kernel matched here is the source's **synchronous one-particle Euler
kernel**, with denominator 63. One globally sampled seam attempt has
denominator `63N` for N carriers; N such attempts and a Poisson semigroup
are different maps. The local walk's coherent instrument differs
from C, and its measured positions after two steps need not follow P².
Coherent arc retention, the two-layer timing and state preparation remain
candidate inputs. A local fixed rule does not by itself select physical time.

On each nonconstant vertex mode the inherited walk eigenphases are
`±acos(λ_P)`. A normalized positive-frequency embedding is explicit, but
uses an inverse square root of `I−P²`: locality of U does not imply local
spectral preparation or readout. Extra cycle and symmetric arc modes are
present and counted, rather than silently removed or labeled physical gauge
particles. The general spectral mapping is established
[quantum-walk theory](https://arxiv.org/abs/1506.06457).

The actual graph incidence also gives a local Hermitian vertex/edge block
operator whose square contains the two Laplacians. Its added cycle zero
modes and constant mode are explicit. Continuous-time evolution has nonzero
propagation tails, so this alternative does not inherit the walk's strict
stepwise light cone. Convergence of its vertex Laplacian alone would not
supply compatible edge/polar maps or select its energy-sign interpretation.

### Controlled coarse memory and positive local energy approximations

For the actual matching-glued carrier family write `L=L0+E`, where E contains
external seams. Let V repeat a cell value with normalization `1/sqrt12`,
`Q=I−VV*`, and write L in retained/hidden blocks `[[A,B],[B*,D]]`. D is
restricted to Q, not inverted on its full-space zero modes. The native
internal spectrum and the external matching imply

```text
g=5−sqrt5,       D≥gI,
A=L_dual/12,     BB*=2A−A²≤(3/4)I.
```

These are bounds independent of cell count under the stated routing
hypotheses. Exactly N one-particle eigenvalues lie in `[0,1/2]`, and the
other 11N are at least g. This identifies a separated coarse band; it does
not itself identify a physical continuum or many-body field.

In the declared per-carrier Poisson clock, eliminating hidden v gives the
exact equation

```text
u_dot(t)=−Au(t)/63
  +(1/63²) integral_0^t B exp[−D(t−s)/63] B* u(s) ds
  −B exp(−Dt/63)v(0)/63.
```

The kernel norm is at most `(3/4)exp(−gt/63)`. Hidden preparation v(0) remains
an input; equal initial cell averages do not determine it. The unitary
exchange counterpart has an oscillatory kernel, with no corresponding
finite-dimensional decay bound.

The static coarse energy `S_c=A−BD^−1B*` obeys

```text
[g/(g+2)] A ≤ S_c ≤ A.
```

For real Laplace parameter z>0 in L units, the exact compressed resolvent lies
between `(zI+S_c)^−1` and that operator divided by `1+δ(z)`, where
`δ(z)=2(g+2)z/[g²(g+z)]`. This gives a controlled low-frequency response,
uniform in cell count. It is not an all-time autonomous coarse evolution.

There is also a uniform locality result. Let `D0=Q L0 Q` on Q and
`K=D0^−1/2 Q E Q D0^−1/2`. Then `0≤K≤rI`, `r=2/g<1`. Truncating
`(I+K)^−1` after power m gives an effective energy S_m with

```text
0 ≤ (1−2/g)A+A²/g ≤ S_m ≤ A,
||S_c−S_m|| ≤ 3 r^(m+1)/(4g).
```

S_m couples only cells at graph distance at most m+2. Thus positive
finite-range energies approximate the true coarse energy exponentially well.
This analytical extension is independently reviewed; it does not assume a
new routing or fit a long-range cutoff.

**The probability interpretation fails a separate exact test.** On the
actual level-0 graph, `S_c[16,17]≈+0.0019394`, with its positive sign certified
by rational arithmetic. Therefore `−S_c` is not a Markov generator. Positive
energy and a controlled response do not imply nonnegative transition rates.
The exact projected diffusion remains stochastic but has memory. These
results close a controlled coarse-response interface without hiding the
previous failure of cell-only autonomy.

### A periodic native-carrier continuum benchmark

To establish an actual limit, use a declared periodic screen rather than
infer convergence from the two archived sphere levels. Each Z² period has
two native twelve-port carriers, retaining all thirty seams in each. Match
face ports `(0,1,5)` across the two carriers with translations `(0,0)`,
`(1,0)`, `(0,1)`. There are 24 ports and 63 seams per period. The gluing and
coordinate interpretation are explicit additional choices.

The 24×24 Bloch Laplacian has a simple acoustic band. Exact rational cell
correctors, independently checked against the native three-terminal network
reduction, give

```text
G=(5/246) [[1,−1/2],[−1/2,1]],
lambda_0(L(k))=kᵀGk+O(|k|⁴).
```

Native adjacent-port resistance `11/30` fixes the reduced face conductance
`20/11`; the coefficient is derived from the source seams, not fitted to
sampled eigenvalues. The general Bloch/effective-mass method is established
[periodic-graph analysis](https://arxiv.org/abs/1502.02989).

Take an N×N periodic array, spacing `h=1/N`, and declare one local walk step
to take time h. On the positive-frequency vertex sector,

```text
Omega_h=acos(I−L_h/63)/h,
omega(p)²=(5/7749)[(2*pi*p1)²+(2*pi*p2)²−(2*pi*p1)(2*pi*p2)].
```

The limiting wave operator is `−div[(2G/63)grad]` on the declared flat
unit two-torus. First-order intracellular correctors produce a **full
frequency residual O(h)** on every fixed Fourier band, including optical
modes outside the acoustic subspace. Its weighted-resolvent bound is checked
on five prescribed momenta at N=16,32,64. Maximum residuals are approximately
`0.03960, 0.02017, 0.01019`; the uncorrected constant-vector residual does
not vanish in this comparison. The infinite-refinement assertion follows
from the analytical proof, not from fitting those three numbers.

A graph-path estimate gives the uniform bound `Omega_h²≥1/3654` after
removing the single constant mode. Combining finite-band evolution control
with that infrared bound and Fourier-tail approximation gives the following
scoped quantum result:

**With Bose CCR, Gaussian vacuum and positive-frequency interpretation
supplied, bounded Weyl-product expectations converge for arbitrary real
mean-zero L² field smearings, uniformly on bounded time intervals.** The
same holds in states prepared by finitely many such Weyl displacements.
All optical lattice modes remain in the operator and vacuum. The limiting
one-particle space contains all nonconstant Fourier modes, not a fixed
finite set. Discrete walk powers give the corresponding evolution at
stroboscopic times, with the rounding error controlled.

This is a constructed massless free-field observable limit in two spatial
dimensions. Its precise boundaries are part of the result:

- Periodic gluing, coordinates, time scaling, arc states, Bose statistics,
  vacuum and spectral preparation/readout are declared inputs. They are not
  derived from native integer repair or selected by matching its kernel.
- The Fock extension is not the native hard-core many-body C/H theory.
  Interactions, particle statistics selection and physical gauge/matter
  sectors have not been derived by this free construction.
- The massless constant mode has no normalizable oscillator vacuum and is
  excluded explicitly. The claim covers the zero-mean sector and suitable
  derivative-field observables, not an unrestricted point-field or scalar
  microcausality statement including that mode.
- Arbitrary canonical momentum smearings requiring positive frequency
  powers, ultraviolet composite observables and a physical 3D bulk are not
  supplied by the stated Weyl limit.

The existing [golden-address scalar detector theorem](../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex)
uses a supplied scalar action on different source coordinates. The new
benchmark computes its acoustic coefficient and frequency from an explicit
native repair graph and its declared coherent walk. Their premises differ;
neither may silently substitute for the other or for the missing physical
attachment of C.

### Alternative: the doubled-sheet coherent walk to the same scalar limit

The [golden-address scalar packet](../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex)
supplies the actual sorted coordinates, positive tensor capacities,
Dirichlet action and a massive free-field detector limit. The earlier
coherent-walk connection factors its full operator through weighted agreement:

```text
Σ_e lambda_e(I−E_e^M) = M^−1 K_internal,
A_q = I + M^−1/2 K_D M^−1/2.
```

For this alternative, define `P=Σ_e(lambda_e/Λ)E_e^M`,
`Λ=Σ_e lambda_e`: one completed weighted agreement event. Keep the original
spatial conductances on each of two copies of every site and join the
copies at site i with conductance `(d_i+m²M_i)/2`, where d_i is the summed
boundary grounding. The full doubled process is conservative. The binary
copy and its odd sector are supplied additional structure. In P and Λ_q,
the sum is over every edge of this full doubled graph, whose capacities are
(M,M); it is not the original interior graph normalization.

The interior massive operator is not a conservative Markov generator:
`(M^−1 K_D+I)1=M^−1 d+1`, which is nonzero. In the doubled conservative
extension, the normalized odd embedding `J_-v=(v,−v)/sqrt(2)` instead obeys
`B_extended J_-=J_- A_q` exactly. This keeps every mode and supplies the
missing mass/boundary terms without destroying the full graph's equilibrium.
The signed odd reading is not a nonnegative probability density.

Because the whitened endpoint expectations are orthogonal projections,
P's symmetric discriminant lies between zero and identity. The local arc
walk has the odd-sector frequency and declared step

```text
delta_q=1/Λ_q,    a_q=sqrt(2 delta_q),
Omega_q=acos(I−delta_q A_q)/a_q.
```

This P is exactly one weighted agreement event, not an Euler proxy for one
native attempt. The coherent walk itself uses a whole local coin layer and
reversal; a shared one-step kernel does not identify their execution costs.
For golden refinement, `Λ_q>=3q²(q−2)(q−1)²`, hence `delta_q→0`.
The clock is population-dependent. Graph locality does not establish a
uniform operational light cone in the supplied coordinates or identify the
existing count clock.

For any fixed target sine-mode frequency w and normalized sampled column j,
the scalar spectral theorem gives, eventually along refinement,

```text
||(Omega_q−w)j||
  <= sqrt(2) ||(sqrt(A_q)−w)j|| + delta_q w³/6.
```

This controls the **full** residual, including leakage outside the comparison
band. It follows from the derivative bound `g_delta'(s)<=sqrt(2)` for
`g_delta(s)=acos(1−delta s²)/sqrt(2 delta)`, where `delta s²<=1`, and the
scalar phase-error bound when `sqrt(2 delta)w<=1`. It requires no unproved
operator-Lipschitz assertion for arbitrary matrices.

Retain the original scalar field map, Fock vacuum and displaced preparations:
use `sqrt(A_q)^−1/2` to normalize the fields, while replacing only the
evolution `dGamma(sqrt(A_q))` by `dGamma(Omega_q)`. Both are number
conserving and preserve the same vacuum. The transferred residual, old
negative-power bounds and Riemann-sampling tail argument yield the **same
three-dimensional massive Dirichlet free-field Weyl-product limit** for the
old theorem's bounded piecewise-continuous L² smearings with
Riemann-integrable products. Integer-step rounding vanishes as well.
Finite-q phases change; the old finite detector error bars need an added
comparison term and cannot be reused unchanged.

Exact arithmetic checked all 64 site capacities and every action entry
against the existing q=5 execution receipt, including its 144 spatial edges.
The conservative extension has 128 sites and 352 seams. Whole-operator arc
and residual checks also passed at q=3 and q=5. These small checks validate
the join; the continuum statement follows from the analytical transfer,
not extrapolation from two grids.

The positive-frequency embedding is spectral and need not be spatially
local. Local scalar preparation/readout on the arc workspace is therefore
an attachment to construct. The declared source action, rates, mass,
boundary, sheet sector, coherent lift, Bose representation and clock are
not selected by this identity. It is a concrete common construction under
those inputs, not a proof of the physical law or an interacting QFT.

### What the reread simplifies, and what remains

The [action-measure theorem](../paper/tex_fragments/SOURCE_ACTION_DENSITY_IDENTIFICATION.tex)
reconstructs normalized positive capacities from a connected directed
response that admits a positive diagonal symmetrizer, using reciprocity. Thus capacities need not be fitted independently
once that response is given. Its selection, overall action scale and
identification with physical volume remain open. Equal source event-count
weights are not these nonuniform finite action capacities.

Existing [time-refinement](../paper/tex_fragments/SOURCE_SCALAR_TIME_REFINEMENT.tex),
[regional reconstruction](../paper/tex_fragments/SOURCE_SCALAR_REGIONAL_TIME_SLICE.tex)
and [quantum instrument](../paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex)
results supply conditional interfaces. The canonical completion keeps their
law and ports unchanged. For the alternative changed walk, finite execution
certificates require comparison; sharing addresses and allowed edges
alone does not preserve event ancestry or a clock certificate. The regional
collar theorem's minimum exterior linear access is `rank(A_RE)`; access to
those quantum observables remains an explicit hypothesis.

The [coherent source criteria](research/COHERENT_SOURCE_CLOCKS.md) give a
specific next test: preserve both complete computational and Fourier read
tables, including inputs entangled with a reference. Under the stated channel
contract this forces coherent identity transfer, whereas one-basis agreement
permits dephasing. This can attach an existing quantum source construction;
a spectral intertwiner alone does not supply that transport/readout service.
The consensus Hamiltonian theorem also remains available once an actual
continuous automorphism group is established; it does not choose that group's
operator or clock from a dissipative repair channel.

The archived spherical graphs retain their own convergence question. Their
cell-constant embeddings fail dynamical intertwining: a fine parent indicator
gives child values `[2,2,2,0]`, which no overall time rescaling fixes. The
uniform memory and locality bounds give tools for that distinct routing.
Neither the golden scalar join nor the periodic benchmark proves its limit.

The remaining connected targets are:

1. **Select the source law.** Justify the weights, reference, mass, inertia
   identification and quantum operations, preserving the native results used. Fitting the desired action backwards is not
   physical selection.
2. **Attach quantum access and time.** Construct the actual preparations,
   readouts, code transport and clock with complete resource and event-order
   accounting. Reuse the existing conditional interfaces where their
   hypotheses hold.
3. **Join interacting matter and responsive geometry.** Attach statistics,
   charges, stable excitations, conserved stress and universal coupling to
   one common process. The existing conditional matter and gravity theorems
   retain their premises; a free scalar does not discharge them.
4. **Freeze an eligible physical test.** Identify an unavoidable calibrated
   relation, enumerate remaining choices and uncertainties, and compare it
   with independent observations under a stated rejection rule.

No OPH core theorem defect was established by this targeted reread. The
gain is an exact local field/instrument join, a conditional selection
lemma and a finite source-control compiler for its local operation factors. The remaining task is physical
selection and joint operational realization, rather than rebuilding its
free-field theorem or treating its field coordinates as independent particles.


## 5. Entropy, geometry and memory connections retained

| Constructed connection | Exact result | Physical boundary |
| --- | --- | --- |
| **Transport and phase** | `Δx=−b_e j`; every region's change equals its boundary flux. A branch gains phase `exp[−i(θ_i−θ_j)j]` under local number phases, compensated by its record. | Covariance constrains quantum lifts; it does not uniquely select C or identify electromagnetism. |
| **Response and fluctuations** | For uniform active seam attempts, `E[Δx∣x]=−L_ico x/60`. At uniform binary equilibrium, susceptibility is `χ=(5/33)(I−11ᵀ/12)`. The exact discrete balance uses conditional innovations rather than raw increments. | This is preparation response and finite fluctuation balance, not a supplied temperature, external field or physical stress law. |
| **Observer intervention** | A working read between overlapping C events changes a later coherence expectation from 1/4 to 0 while preserving the averaged old record state. All 312 two-event classical paths reverse; all 6084 working matrix units match C_f C_e. | Compatible common-state restrictions also hold before repair, so they do not independently identify an A2 repair target. Hidden-neighbor dependence remains. |
| **Event order** | Disjoint recorded operations commute after occurrence alignment. Linear extensions of a fixed dependency order give the same completed operation. | Intermediate probes and timestamps must be included in the dependency order; the theorem does not identify all observed histories. |
| **Entropy completion** | On the shell, `D(ρ‖Eρ)=S(Eρ)−S(ρ)`; Eρ uniquely maximizes entropy at fixed retained-algebra data. C approaches that completion. | The entropy deficit need not follow the channel's exact halving rate. The full 78-state space has no faithful stationary reference for this repair. |

### Current-derived port geometry and clock relations

Let D be oriented incidence, `G=L_ico⁺`, `P_c=I−DᵀGD`, and
`A_n=Σ_(k≤n)J_k` the integrated recorded current. Then `P_c A_n` is a
martingale, and at uniform binary equilibrium,

```text
Cov(P_c A_n) = n(5/33)P_c/30,
Var(Y_n) = n/66                  for circulation around a native triangle,
(I−P_c) A_n = DᵀG(X_0−X_n).
```

The gradient component stays bounded; cycle records retain information absent
from endpoint densities. The circulation variance measures ensemble attempt
count. A single squared circulation is noisy and need not grow monotonically.
Uniform pile preparations with the same mean loads give first-attempt variance
1/60 instead, so the calibration depends on the declared state.

There is a direct connection to the existing native **port** geometry. Let
`B=E[ΔX ΔXᵀ]` and `ν=E[(cᵀJ)²]` for one triangular circulation attempt.
At uniform binary equilibrium, branch counting gives

```text
B=L_ico/198,   ν=1/66,   L_record=3B/ν=L_ico.
```

The factor three is the loop length. A common event rate and load-unit scale
cancel. The same operator is independently recovered from the mean response
to the twelve pile preparations. Its Green matrix gives port resistances
`11/30`, `7/15` and `1/2` for the three nontrivial pair classes. Its slow
eigenspace has rank three, and four times that projector is exactly the
documented unit-diagonal native Gram. The finite centered response has rank
eleven; rank three is the specified spectral/long-time readout. This recovers
the existing source geometry from records, within that declared regime.

With `κ=5/33`, continuous Poisson attempt rate λ gives slow density decay
`γ_s=λ(5−√5)/60` and circulation variance rate `ν_t=λ/66`. Thus

```text
κ γ_s/ν_t = (5−√5)/6.
```

This is a common-process clock/response relation with no fitted rate.
An exactly localized pile instead yields a rank-five activity operator;
a half-equilibrium/half-localized mixture yields a connected but differently
weighted operator, while measured mean propagation remains unchanged. These
controls reject identifying arbitrary state-dependent activity with the
metric controlling that mean response. The physical source-clock and
responsive-geometry attachments need their explicit maps.

### An entropy-driven metric without changing the dynamics

Let P be the shell transition matrix and use continuous Poisson attempt time
at declared rate one. Configurations connected by one native exclusion move
have generator rate `q_xy=1/60`. For positive probabilities p, define

```text
θ(a,b) = (a−b)/(log a−log b),   θ(a,a)=a,
K_p = Σ_{unordered x~y} q_xy θ(p_x,p_y) b_xy b_xyᵀ,
p_dot = (P−I)p = −K_p log(p/τ),   τ_x=1/66.
```

Here b_xy is configuration incidence. K_p is positive on zero-sum forces;
its inverse there gives a tangent metric. Relative entropy dissipates exactly:
`dD(p‖τ)/dt=−log(p/τ)ᵀK_p log(p/τ)`. The general reversible-chain
construction is established mathematics; this is its normalized application
to native repair ([Maas, 2011](https://arxiv.org/abs/1102.5238)).

For occupancy matrix X, the projected mobility is `M_p=XᵀK_pX`.
At the specified full state, the least configuration tangent cost inducing
a zero-sum port velocity v is `vᵀM_p⁺v`. Thus it has a precise transport
meaning. It depends on correlations beyond the mean `μ=Xᵀp`.
At equilibrium, `M_τ=L_ico/396=B_τ/2`; in general `M_p≤B_p/2`, where
B_p is the raw increment second-moment rate. The mean law satisfies
`μ_dot=−L_ico μ/60`: changing mobility and entropy force combine to give
that same response. Fixed mean propagation does not rule out this metric.

On the entropy counterexample `p_ε=τ+ε(δ_01+δ_23−δ_02−δ_13)/66`, the
full port mobility changes at first order while every mean stays fixed.
For A={0,1,7}, its regional contraction is
`1_AᵀM_pε 1_A=1/44−ε²/3960+O(ε³)`. The regional first derivative
is zero, so this metric does not itself cancel the first-order entropy
witness. It is a classical configuration transport metric and projected
tangent cost, with a supplied count-time scale. A spacetime metric,
physical action and quantum refinement have not been identified with it.

The binary configuration transition matrix P gives a genuine resistance
metric `R_ab` through `(I−P)⁺`. Mean first-passage times satisfy
`H_ab+H_ba=66R_ab`. A common supplied Poisson rate λ cancels in

```text
mean commute duration × triangle-circulation variance rate = R_ab.
```

This joins clock and response without fitting a rate. It does not select
physical units or make configurations into spatial positions. Hitting probes
use continuing active repair: a capacity-two halt can leave a target unhit
forever with positive probability. The supplied port chart and exact linear
response remain fixed; these new probes require a physical geometric map.

### Retained-port costs and exact boundary memory

For A={0,1,7}, the native one-hop exterior collar is
C={5,6,8,9,10,11}; the remote exterior is D={2,3,4}. There is no A–D seam.
Let U=A∪C. At a specified positive full probability p, prescribing a retained
velocity v_U while allowing hidden velocities to vary has minimum cost
`v_Uᵀ(M_UU)⁻¹v_U`. Requiring hidden velocities to vanish instead gives

```text
S_p = M_UU − M_UD (M_DD)⁻¹ M_DU,
cost = v_Uᵀ S_p⁺ v_U,   Σ v_U=0.
```

Both costs follow by minimizing the same full configuration tangent metric.
They respond to the tested pair-correlation perturbation despite unchanged
mean loads. This distinction prevents using a Schur complement for the wrong
constraint. Neither minimization makes the hidden dynamics instantaneous.

The actual mean equation is `μ_dot=−Hμ`, with `H=L_ico/60`. Eliminating D,
writing u=μ_U and w=μ_D, gives the exact retained evolution

```text
u_dot(t) = −H_UU u(t) − H_UD exp(−H_DD t) w(0)
           + ∫₀ᵗ H_UD exp[−H_DD(t−s)] H_DU u(s) ds.
```

This specifies both the memory kernel and its initial-state input. The
static response uses `H_UU−H_UD(H_DD)⁻¹H_DU`; finite-frequency response uses
`(zI+H_DD)⁻¹` instead. The existing native-transport paper establishes that
general distinction. Here it is joined to the entropy transport metric and
tested on this exact collar, with explicit time normalization.

The preparations `(I₆₆/66+|23〉〈23|)/2` and `(I₆₆/66+|24〉〈24|)/2` have identical
complete U marginals and regional/collar charge distributions. Their collar
mean derivatives differ; their A means differ at second derivative. Preparing
the uniform exterior completion Φ_U once also fails to keep that family
invariant under all-thirty repair. Thus unrestricted reduced dynamics needs
the retained memory or an independently justified approximation. The local
microscopic law itself is unchanged.

### Two source contrasts close the retained mean dynamics

The identified memory has a finite construction. For the remote triangle,
`H_DD=(6I−J)/60`, where J is the three-by-three all-ones matrix. Its decay
rates are 1/20 on the uniform direction and 1/10
on its two-dimensional contrast plane, in the declared attempt-time units.
Conservation fixes the uniform remote load from the retained sum. Define

```text
z = (μ₂−μ₄, μ₃−μ₄),
B = [[1,0,−1],[0,1,−1]],
T = [[2,−1],[−1,2],[−1,−1]]/3,
w = (2−Σu)1_D/3 + Tz.
```

For centered retained means `ū=u−1_U/6`, the exact autonomous system is

```text
ū_dot = −[H_UU−H_UD 1_D 1_Uᵀ/3] ū − H_UD Tz,
z_dot = −B H_DU ū − z/10.
```

The reconstruction and generator agree with native updates on every shell
basis mean. The hidden-to-retained matrix `H_UD T` has rank two. Positive
preparations with the same entire retained marginal realize both independent
hidden directions, proving that fewer than two linear hidden summaries cannot
close these retained derivatives over the declared preparation family.

This is nine retained means plus two contrasts: all twelve first moments
modulo total conservation. Eliminating z gives a single-exponential memory
kernel with two independent amplitudes. Initial contrasts and access to their
current increments must be supplied; a mean is an ensemble quantity, rather
than a deterministic individual trajectory. The original four-pair entropy
tangent changes none of these eleven coordinates while changing entropy and
mobility. Thus mean closure is complete, with full-state closure a separate
question. An open patch may use the exact memory description; autonomous
local evolution is required only for a claim that uses it.

### Full diagonal closure costs twenty additional coordinates

The retained region U=A∪C has 46 distinct occupation patterns. Ten patterns
have three compatible remote states; the other 36 have one. In each
three-state fiber, retain two probability contrasts alongside its total y:

```text
z1=p1−p3,  z2=p2−p3,
p3=(y−z1−z2)/3,  p1=(y+2z1−z2)/3,  p2=(y−z1+2z2)/3.
```

These twenty extra coordinates give an invertible coordinate map R of the
66 probabilities and hence an exact transformed generator `R Q R⁻¹`.
Observability ranks of the regional readout J through two derivatives are
`rank(J)=46`, `rank(J,JQ)=65`, `rank(J,JQ,JQ²)=66`. Therefore twenty is
necessary and sufficient among linear autonomous extensions valid for all
positive diagonal preparations. This restores every normalized diagonal
degree of freedom; it is not compression or a nonlinear minimality theorem.

Positive states with identical regional data and all port means can have
unequal regional derivatives. Another positive pair agrees even in its first
regional derivative but differs at the second. The two-contrast mean closure
above cannot be promoted to full-state closure. Classical current histories
can initialize the diagonal coordinates through repeated preparations, while
quantum coherences require the distinct coherent readout constructed above.
Tomography alone does not make an open regional quantum channel autonomous.

### Regional entropy and recorded modular balance

For A={0,1,7}, tracing τ_B gives probability 6/11 for the empty pattern,
3/22 for each one-occupied pattern, and 1/66 for each two-occupied pattern.
On the seven-dimensional feasible support,

```text
K_A = −log ρ_A = log(11/6) I + log4 Q_A + log(9/4) choose(Q_A,2),
δS_A = Tr(δρ_A K_A).
```

The last term follows from fixed-total counting, not an added microscopic
interaction. For outward recorded current J_A, `ΔQ_A=−J_A` and

```text
ΔK_A = −log4 J_A + log(9/4)[−Q_A J_A + J_A(J_A+1)/2].
```

A preparation tilt at port 0 gives regional entropy derivative
`(4log2+log3)/22`. Applying the crossing seam (0,5) changes it to
`(4log2+log3)/66`; the difference is the same recorded modular flux.
This is an actual regional-state calculation, not an inference from global
stationarity. Physical area, stress and boost/modular-flow correspondence
have not been derived from it.

The complete first differential for a region with b exterior ports is,
on its feasible support and for b≥2,

```text
δS_A = log[(b−1)/2] δ〈Q_A〉 + log[2b/(b−1)] δ〈choose(Q_A,2)〉.
```

This covers arbitrary shell-density variations, including directions outside
the port-tilt family. The
pair term exposes an omitted response. For A={0,1,7}, take the positive
nearby preparation curve

```text
p_ε = τ_B + ε(δ_01 + δ_23 − δ_02 − δ_13)/66,   |ε|<1,
```

where ij denotes the occupied pair. Every individual mean load is constant,
but `dS_A/dε=log(9/4)/66`. Fixed cut capacity and boundary current activity
both have zero derivative. They cannot supply a compensating area term on
this tangent. The tested one-crossing-event mutual information between
`(seam,coin,current)` and the post-event regional readout also fails a
constant coefficient across independent port tilts. Its failure is distinct
from the capacity/activity witness; not every information readout is ruled out.

These tests use fixed region cardinality as a candidate volume and do not
construct the physical fixed-volume family. They reject those attachments;
the following counting construction supplies a response they missed.

### Positive entropy balance from exterior completions

Use the existing fixed-total shell split, with a interior and b exterior
ports. In regional charge sector q there are `d_A(q)=choose(a,q)` interior
patterns and `d_B(q)=choose(b,2−q)` compatible exterior patterns. The reference
regional state has eigenvalue `d_B(q)/66` on each feasible interior pattern.
Define a central observable by counting the exterior completions:

```text
Z_B = Σ_q log d_B(q) Π_q,
−log τ_A = log66 I − Z_B,
S(ρ_A)+〈Z_B〉 = log66 − D(ρ_A‖τ_A).
```

The last identity holds for every allowed shell density, including quantum
coherences within a regional charge sector. Its maximum is log66, attained
exactly when `ρ_A=τ_A`. Hence every normalized first variation there obeys
`δS_A+δ〈Z_B〉=0`, with a positive coefficient of one and an exact
relative-entropy deficit away from the maximum. Z_B was counted from source
multiplicities, not fitted to cancel a selected entropy perturbation.

On the original four-pair tangent, `δ〈Z_B〉=log(4/9)/66`, precisely
balancing the missing entropy response. Two-event diagonal reconstruction
reads this observable and the first entropy variation, even for quantum
perturbations at the reference. It does not recover arbitrary off-reference
von Neumann entropy when local coherences are unknown.

This closes a finite positive counting-stationarity subproblem. It also
sharpens the next construction: **realize these multiplicities in a physical
collar and derive its area map**. Exterior bulk completion counts are not
automatically gauge-representation dimensions. In particular, ordinary U(1)
charge irreducible representations have dimension one; their representation
edge term alone cannot be renamed `log choose(b,2−q)`. The flagship's edge
entropy theorem requires its actual collar decomposition. Locality,
refinement, area normalization, physical fixed volume and stress/modular-flow
identification need that common construction.

Fixed total charge makes ρ_A block diagonal in regional q, even when the
global state has coherence between different q sectors. Global sector-entropy
bookkeeping must retain the pinching correction
`C_Q=S(Δ_Qρ)−S(ρ)`. Exact coherent controls distinguish that correction
from the regional identity above. The counted maximum does not by itself
prove that every region's functional grows monotonically under native repair.

### Exterior repair implements the counted completion

The completion has an implementation with the same primitives. Fix A and
schedule only native seams wholly in its exterior B. If that exterior graph
is connected and has at least three ports, positive fixed seam weights give
the reduced C and E limits

```text
Φ_A(ρ) = ⊕_q ρ_A,q ⊗ I_B,2−q / choose(b,2−q),
fixed algebra = ⊕_q B(H_A,q) ⊗ I_B,2−q.
```

Here ρ_A,q is the unnormalized regional charge block. Relative-number phases
force surviving exterior bra and ket patterns to coincide; connected endpoint
swaps then make their coefficients uniform at fixed exterior charge. This
independently defines the fixed algebra. Each native C_e is self-adjoint with
spectrum in {0,1/2,1} and has the same fixed space as E_e; their positive
mixture converges to the stated conditional expectation.

For A={0,1,7}, all 4,356 matrix units verify the result on 18 exterior seams.
The algebra has dimension `1²+3²+3²=19`. Φ_A preserves the entire regional
density, including within-charge quantum coherence, and uniquely maximizes
global entropy at that fixed regional state:

```text
S(Φ_Aρ) = S(ρ_A)+〈Z_B〉.
```

The completed state retains the shared charge-sector correlation across the
boundary; it removes within-sector correlations with the exterior.

This joins the counting identity to an actual scheduled repair process.
The restrictions matter: one- or two-port exteriors can retain additional
vacuum/full coherences; disconnected exteriors retain component charges.
The all-thirty schedule includes boundary crossings and changes ρ_A, so it
does not implement this regional-preserving limit. A finite recorder that
halts after two writes cannot implement an infinite-time completion either;
continuing reduced evolution uses its declared memory/reset policy. Physical
collar locality and area identification require further work on this same
construction.

Coherent record retention preserves total joint entropy. Forgetting records
can increase reduced entropy. The full 78-state uniform state instead loses entropy
under repair, so shell assumptions cannot be dropped. Likewise, full-space
norm convergence to τ_B need not make relative entropy to that singular
reference finite while transient population remains.

### Nested regions compose exactly

For A contained in R, their retained algebras satisfy `A_A⊂A_R`. The
trace-preserving conditional expectations therefore obey

```text
Φ_A Φ_R = Φ_R Φ_A = Φ_A,
D(ρ‖Φ_Aρ) = D(ρ‖Φ_Rρ)+D(Φ_Rρ‖Φ_Aρ).
```

Applying Φ_A last uses preservation of the A marginal; the other order uses
the exterior counting convolution. The supports satisfy
`supp(ρ)⊆supp(Φ_Rρ)⊆supp(Φ_Aρ)`. Each relative entropy equals the finite
entropy increase of its completion, so those differences telescope, including
singular inputs. Logarithms are restricted to the relevant supports. The
same chain rule holds along longer inclusion chains.

The native chain `{0,1,7}⊂{0,1,2,7}⊂{0,1,2,5,7}` has connected exteriors
with 18, 13 and 10 native seams. Its fixed algebras have dimensions 19, 53
and 126. All 4,356 shell matrix units verify both composition orders, with
independent phase/permutation commutants and exact coherent entropy controls.
This realizes the established expectation-tower identity on the same source.

Inclusion is load-bearing. For A={0,1,7} and N={0,1,5}, completing the
occupied pair {2,7} in the two orders gives port-7 means `2/9` and `17/81`.
The nested tower therefore does not prove general confluence. Limits under a
repeated schedule require their own analysis; a two-operation order effect
does not refute such convergence. The positive result concerns nested regions
and their declared exterior-only limits on one finite carrier, rather than
a spatial continuum.

### Charge-resolved collar and a controlled alignment bound

Use the same A, one-hop collar C and remote D defined above. At the uniform
fixed-total reference, the bare collar has `I(A:D∣C)>0`: conditioning on C
leaves a charge anticorrelation between A and D. Retain Q_A in a declared
three-value register, using the isometry `V|x〉=|x〉|Q_A(x)〉`. For
`σ=Φ_Aρ`, the registered state has the explicit form

```text
VσV* = ⊕_q p_q ρ_A,q^normalized ⊗ ω_CD,2−q ⊗ |q〉〈q|,
ω_CD,2−q = I_CD,2−q / choose(9,2−q),
I(A:D∣C,Q_A) = 0,   I(A:CD∣Q_A) = 0.
```

The preselected collar split has center q, a one-dimensional left factor,
and all C degrees of freedom on the right. Reading `(Q_A,Q_C)` and appending
the uniform compatible D charge sector gives an explicit CPTP recovery;
all 4,356 matrix units verify its action on the completed family. Replacing
the normalized right environment within each q sector preserves the entire
A density, including its within-charge coherence.

This is a concrete finite factorization, with an explicit resource cost.
Current records update Q_A from an initial charge offset; they do not supply
that offset automatically. Coherent copying preserves global coherence,
whereas reading or discarding the charge register pinches cross-charge
coherences. The completed family is charge block diagonal, so that distinction
does not disturb it. Physical implementation and timing of the register need
their own map.

For arbitrary shell inputs, the completion supplies a controlled comparison:

```text
δ = S(Φ_Aρ)−S(ρ) = D(ρ‖Φ_Aρ),
‖VρV*−V(Φ_Aρ)V*‖₁ ≤ √(2δ).
```

Every exterior-only native channel fixes the comparison and preserves the
regional marginal, so data processing makes δ nonincreasing. The connected
exterior limit drives δ to zero. Thus the unchanged repair law approaches
this declared aligned family with an explicit error certificate. This bound
uses completion entropy deficit; it is neither a bound from raw conditional
mutual information alone nor a uniform continuum estimate.

The controls enforce the distinction between Markovity and alignment. A
single crossing and one all-thirty average preserve registered Markovity on
the tested completed preparation while breaking the preselected q product
split. The two-event path `C_(0,5)` then `C_(4,5)` breaks screening too.
These controls refresh the register to the resulting Q_A; the failure
persists with that readout resource granted.
The [Markov-state structure theorem](https://arxiv.org/abs/quant-ph/0304007)
provides some split at zero conditional mutual information; the core's MSA
condition requires alignment with the specified physical factors. The
first trivial-left construction is not identified with those factors. A
stronger native bulk split is available: incidence divides C into
L={5,8,10} and R={6,9,11}. Adding Q_L to the declared center gives a 3×4
block with both factors nontrivial, and exact alignment of completed states.
That block has empty A; it is not evidence of a nontrivial matter region
coupled to a physical gauge edge. Boundary repair can change Q_L, so this
center is not protected for arbitrary inputs.

The actual face stabilizers are spatial C3 or S3 actions. Their pointwise-A
stabilizer is trivial; they do not provide the required nontrivial
collar-supported gauge action while preserving the full A observables.
Covariance also does not imply a Gauss invariant-vector constraint: a tested
spatial invariant vector acquires leakage 4/27 under a native exterior event.
An appended invariant S3 pair supplies a constant log2 entropy only by adding
a prepared spectator resource; it has zero area response. These tests reject
those particular physical-edge identifications, not every internal symmetry
or encoding. The physical edge/area connection remains open.

Finally, a collar does not turn all exterior entropy into area. For
`w(s∣q)=choose(6,s)choose(3,2−q−s)/choose(9,2−q)`, counting gives

```text
log choose(9,2−q) = H(w(·∣q)) + Σ_s w(s∣q) log choose(6,s)
                                + Σ_s w(s∣q) log choose(3,2−q−s).
```

The last term is remote bulk entropy; at the reference its mean is
`(5/11)log3`. Ordinary U(1) charge irreps have dimension one. A physical
edge-dimension/area map must therefore do more than retain Q_A or rename Z_B.

## 6. OPH results retained as inputs and attachment targets

The current work uses the documented source and preserves its existing
results. Their premises must be realized together, not replaced by simpler
analogies or treated as if the results did not exist.

| Relevant foundation | What is retained and what needs joining |
| --- | --- |
| [A1–A3 and carrier](AXIOM_REFERENCE.md), [microphysics](../paper/screen_microphysics_and_observer_synchronization.tex) | Oriented twelve-port incidence, inverse pairing, A5 symmetry, exact rank-three Gram geometry and native scalar repair. The 78-state configuration lift is an added representation of a declared sector, not the entire A1 algebra. |
| [Consensus and observers](../paper/reality_as_consensus_protocol.tex), [shared observer/support construction](../paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex) | Concrete protected records, proper overlaps, repair/normalizer and support/refinement constructions exist. Their acceptance, confluence, quotient and naturality hypotheses must be attached to the candidate's operational menu. |
| [Canonical repair RFC](CANONICAL_REPAIR_LAW_RFC.md), [finite completion theorem](../code/quantum_information/ALGEBRAS_AND_REPAIR.md) | Faithful-reference, modular-invariant retained algebras admit canonical conditional expectations. The RFC is an unadopted strengthening: agreement algebras, complete grammar, reference and path-law premises are substantive inputs. This research does not adopt it or select an A2 equalizer by defining it as Fix(C). |
| [Spacetime construction](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex), [operational clocks](research/MASSIVE_OPERATIONAL_CLOCKS.md) | Flat source causal/count limits and moving/matter-clock constructions exist. Their event populations, probes, refinement and rate/scale maps must join the same repair instrument. Ordinal order and Poisson attempt time do not discharge those contracts. |
| [Flagship fields and gravity](../flagship/from_observer_consensus_to_standard_physics.tex), [gravity premise ladder](registers/GRAVITY_PREMISE_LADDER.md) | Scalar/Maxwell/interacting constructions and conditional Einstein reconstruction are available. Their action, state, complete response, stress, entropy/area and common-refinement premises remain load-bearing. Gauge type and ambient dimension do not alone supply the physical model. |
| [Geometric identifiability](../paper/tex_fragments/NATIVE_GEOMETRIC_SOURCE_IDENTIFIABILITY.tex), [physical-identification registry](../claims/physical_identification_registry.json) | Fixed chart volumes cannot become fluctuating geometry by relabeling load noise. New operational constitutive maps are permitted but need tests. P is not a repair coin or physical rate by naming convention. |

The [claim registry](../claims/claim_registry.yaml) and linked proofs remain
authoritative. This is a research ledger, not a new theorem registration.
[#1025](https://github.com/FloatingPragma/observer-patch-holography/issues/1025)
tracks physical test design alongside construction;
[#1026](https://github.com/FloatingPragma/observer-patch-holography/issues/1026)
tracks the common model and its theorem attachments. These finite
results do not establish their physical completion. Existing seen-data
comparisons retain the [postdiction ledger](POSTDICTION_LEDGER.md) and
[frozen-prediction rules](FROZEN_PREDICTION_LADDER.md).

### Dongyang Stephen Chen's contribution to the method

Both **Genesis / The Ladder of Symmetry** and **The Balance Bridge** were
supplied by Dongyang Stephen Chen. Genesis motivated minimal relational
records and exact reconciliation; the Balance Bridge motivated loop-sensitive
readouts, symmetry comparisons and tests against null controls. These ideas
support the current “local event, retained distinction, test the consequence”
method. The particular native quantum instrument, network proof, current
clock and modular-response calculations are constructions tested in this
campaign, not outputs attributed to unseen Genesis or Prism software.

## 7. Retired approaches and lessons retained

The following routes are **obsolete as the active construction strategy**.
A successful diagnostic remains a valid result in its own stated domain;
retiring a route does not declare all its mathematics false. The active
inertial completion is explicitly conditional: its source stiffness is
factored through repair, its kinetic assumption is isolated and tested, and
its operational and continuum interfaces reuse the original action. This
does not reinstate importing a fitted physical action as a derivation. The comparison
law E and the established OPH constructions above remain active inputs.

| Earlier route | What it established | Why it was retired from the main search |
| --- | --- | --- |
| **Import a continuum action or Hamiltonian first** | Supplied Einstein–scalar–Maxwell and quantum models supplied consistency, propagation and binding benchmarks. | Importing the desired physics did not derive it from repair. Too many action, coupling and clock choices preceded the microscopic connection. |
| **Add geometry variables, kinetic terms and counterterms successively** | Positive coupled matter/geometry models and restricted shared principal metrics; some fixed-background free continuum limits. | The explanatory burden grew with each supplied mechanism. These remain benchmarks, while the active approach first tests source readouts and only adds structure for a named failed interface. |
| **Adaptive equilibrium/heat-bath repair as a universal solution** | Controlled finite equilibria and induced interactions in selected models. | Component factorization, volume pinning, positivity failures in some constraints, and cold-limit traps prevented a common physical construction. A stationary distribution alone did not establish the required dynamics. |
| **Direct weak-field or entropy/area attachments** | Explicit comparisons and regional-response calculations. | Tested attachments failed atomic-clock or weak-field checks; a one-sided area proxy failed a common stationarity coefficient across independent perturbations. This rejects those attachments, not every possible geometry or area map. |
| **Generic phase loops, cubic grids and scalar toy worlds** | Agreement with nonzero holonomy, finite positive equilibria, greedy traps and residual motion. | They lacked a demonstrated map preserving the full native twelve-port contract. Loop frustration can be static; phase covariance alone did not produce physical gauge dynamics. |
| **Exploratory uphill repair as the missing microscopic law** | On restricted odd Z3 cycles, any positive exploration restores almost-sure agreement, with an inverse-probability waiting cost. | It resolves a trap of that move menu. Native repair has a different admissible law and can increase some local mismatches while decreasing its actual potential. No P-to-exploration relation was derived. |
| **Separate record repair plus independently evolving residuals** | Agreement could coexist with persistent distinctions and activity; conditional diffusive limits were constructed. | Autonomous component marginals did not establish reciprocal coupling, compatible local clock/response maps or relativistic propagation. The current candidate tests full instruments and records together. |
| **Identify dynamics by matching a mean, channel or stationary histogram** | Exact controls produced identical means with different histories, and identical forgotten channels with different retained instruments. | These summaries cannot select the physical law. Ordered events, interventions, record preservation and same-source maps are required. |
| **Treat an analogy or reported measurement as the missing proof** | Chen's bounded Genesis moment formula separates support at most k from exactly k+1 under distinct nonzero nodes and weights; loop tests motivated useful calculations. | Universal overflow rejection failed finite-field controls. Shared algebraic form, a unit mode or a reported power law does not select physical dynamics. Supplied Prism measurements and broad software claims were not independently reproduced. |

The retained lesson is operational: a simple local rule can support conserved
transport, protected distinctions, entropy completion and nontrivial histories.
Each physical interpretation must earn its place through an explicit map and
a test on that same rule. Historical details and experiment folders remain
local rather than becoming the structure of the active research plan.
