# OPH reverse-engineering table: model and theorem premises

As of 2026-10-11. Generated from [the navigation map](../tracking/model_theorem_map.json)
by `python tools/build_model_theorem_map.py`. Edit the input, then regenerate.

Complete navigation coverage of the current canonical PR identifiers, plus explicit active-model obligations. Not an exhaustive extraction or proof of discharge of all paper/Lean antecedents. Branches are alternatives where stated; links do not alter claim-registry dependency roles.

**Governing search principle: find the simplest shared repair rule that satisfies the known constraints.** Reuse existing results before adding mechanisms. Geometry, reversible response and operational clocks require compatible maps from the same source process; separate component constructions do not establish that common realization.

These groups are workgroup navigation crosswalks, not claim-consumer edges or changes
to proof status. A group links related construction, branch, downstream or metrology
work; that association does not make every listed premise an antecedent of every
listed theorem. Existing claim-registry edges and `explicit_non_consumer` declarations
remain authoritative. Physical attachments excluded by a theorem are not hidden premises
of that theorem. Group statuses describe the remaining connection, not whether a theorem
is proved; the constructed results are recorded separately below. This schema contains
no proof-discharge edges.

## Inventory and limits

Counts below come from the [canonical premise register](../tracking/premise_register.json)
and [claim registry](../claims/claim_registry.yaml). Ordinary assumption names are
explained in the [assumption dictionary](../claims/assumption_dictionary.md).

| Inventory | Count |
| --- | ---: |
| Canonical premise IDs represented in this crosswalk | 84/84 |
| Model criteria represented | 12/12 |
| Registered claims | 326 |
| Unique ordinary `assumptions` names | 771 |
| Unique `imported_results` names | 829 |
| Workgroups | 26 |
| Assumption names explicitly linked by optional group fields | 18 |
| Claims explicitly linked by optional group fields | 9 |

Canonical-ID coverage is a reference check. Ordinary assumptions and imported names
are inventoried only; they have not all been individually, semantically bound to the
model by this page. Names can overlap or encode several obligations. These inventories
cannot be added to obtain a certified total of unproved premises, and coverage is not
a global completeness or proof claim.

| Existing claim-registry role | Unique canonical IDs in that role |
| --- | ---: |
| `consumed` | 58 |
| `open` | 34 |
| `boundary` | 15 |

These role sets overlap. They record metadata on existing claims; `consumed`
does not mean discharged, and `boundary` does not mean an unproved antecedent.

Additional assumption fields are counted separately, not silently folded into
the ordinary-assumption inventory:

| Field | Claims containing field | Occurrences | Unique names | Also ordinary |
| --- | ---: | ---: | ---: | ---: |
| `null_decision_assumptions` | 1 | 2 | 2 | 0 |
| `threshold_kinematic_assumptions` | 1 | 6 | 6 | 0 |

Existing claim classifications: `explicit_edges`: 170; `explicit_non_consumer`: 156.

## Remaining-input and dependency rule

For a selected scope G, let I(g) denote each group's required inputs on its
declared family, state, observable algebra, clock, regime and quantitative error domain.
The remaining-input set is U(G) = union over g in G of I(g), with genuinely identical
inputs counted once. This is a bookkeeping rule, not a machine-computed logical set:
the prose, ordinary assumptions and imports require semantic review. Active
construction, optional branches, downstream physical attachments and metrology are
separate scopes; unioning them does not turn all of them into theorem antecedents.

A result H implies C discharges an input P only when C establishes P on that same
domain, with the required errors, and H has independent supplied or established
support. Substitution retains any outstanding inputs of H. Circular references,
a change of regulator/state/code, or an imported calibration relabeled as a derivation
do not close an obligation. Establishing a conditional theorem and realizing its
hypotheses in the candidate are distinct tasks. This tool checks neither entailment
nor acyclicity of scientific arguments and changes no existing dependency edges.

## Workgroup crosswalk

For the active path, RE-05–RE-12 form the connected action, inertia, quantum
representation, state, history, refinement and clock cluster. Downstream rows
are not additional prerequisites for the supplied finite scalar theorem.

| Group | Model criteria | Scope | Work status | Canonical navigation IDs |
| --- | --- | --- | --- | --- |
| [RE-01: Common repair reference and law](#re-01) | M01, M02 | active | conditional | PR-07 |
| [RE-02: Confluence, acceptance and protected observations](#re-02) | M01, M02 | active | conditional | PR-01 |
| [RE-03: Scheduler, reference and resources](#re-03) | M01, M02, M04 | branch | conditional | PR-05, PR-20, PR-60, PR-61, PR-62, PR-63, PR-77, PR-78, PR-79, PR-80, PR-84, PR-85 |
| [RE-04: Observer regions, causal provenance and time-slice access](#re-04) | M02, M03, M04, M06 | active | conditional | PR-32, PR-33, PR-34, PR-58, PR-81, PR-82 |
| [RE-05: Selected scalar action and relative capacities](#re-05) | M01, M03, M05, M09 | active | conditional | None |
| [RE-06: Inertia and the real dynamical completion](#re-06) | M01, M04, M05, M06, M09 | active | conditional | PR-06, PR-21, PR-22, PR-45, PR-66 |
| [RE-07: Quantum representation and operational probabilities](#re-07) | M02, M06 | active | conditional | PR-02, PR-03, PR-04, PR-43, PR-44 |
| [RE-08: Native coherent controls, phase and reset](#re-08) | M01, M02, M06 | active | conditional | PR-64, PR-65 |
| [RE-09: Preparation of the same coupled vacuum](#re-09) | M06 | active | conditional | None |
| [RE-10: One encoded preparation, evolution and record history](#re-10) | M02, M04, M06 | active | conditional | None |
| [RE-11: Controlled refinement of the selected quantum family](#re-11) | M03, M06, M10 | active | conditional | None |
| [RE-12: Source clock, physical time and execution cost](#re-12) | M04, M06, M09 | active | open | PR-15, PR-22 |
| [RE-13: Physical geometry and responsive metric](#re-13) | M03, M07, M08 | downstream | conditional | PR-09, PR-10, PR-16, PR-52 |
| [RE-14: Gauge algebra to physical gauge fields](#re-14) | M05, M10 | downstream | conditional | PR-11, PR-35, PR-46, PR-54 |
| [RE-15: Physical edge factors and collar recovery](#re-15) | M02, M07, M08 | downstream | conditional | None |
| [RE-16: Modular flow, null translations and stress](#re-16) | M04, M05, M07, M08 | downstream | conditional | PR-23, PR-24, PR-25 |
| [RE-17: Joint state–geometry entropy stationarity](#re-17) | M07, M08 | downstream | open | PR-26 |
| [RE-18: One gravitational branch and its normalization](#re-18) | M05, M08, M09 | downstream | open | PR-27, PR-28 |
| [RE-19: Matter, statistics and interactions](#re-19) | M05, M06, M10 | downstream | open | PR-12, PR-36, PR-41, PR-47, PR-48, PR-49, PR-50, PR-59, PR-67, PR-68, PR-69 |
| [RE-20: Parameter selection and P](#re-20) | M09, M10, M12 | downstream | open | PR-13, PR-17, PR-37, PR-38, PR-39, PR-40, PR-51 |
| [RE-21: Cosmology and dark-sector attachment](#re-21) | M08, M11, M12 | downstream | open | PR-18, PR-19, PR-29, PR-30, PR-31, PR-42 |
| [RE-22: Physical propagation and an eligible test](#re-22) | M03, M04, M05, M12 | downstream | open | PR-53, PR-70, PR-71, PR-72, PR-73, PR-74, PR-75, PR-76 |
| [RE-23: External calibrations and empirical imports](#re-23) | M04, M09, M11, M12 | metrology | imported | PR-14, PR-15 |
| [RE-24: Baryon-specific branch](#re-24) | M10, M12 | branch | open | PR-55, PR-56, PR-57 |
| [RE-25: Low-temperature and common-reference branch](#re-25) | M06, M07 | branch | conditional | PR-08 |
| [RE-26: Common-model compatibility and remaining named inputs](#re-26) | M01, M02, M03, M04, M05, M06, M07, M08, M09, M10, M11, M12 | active | open | None |

<a id="re-01"></a>

### RE-01: Common repair reference and law

**Constructed:** Finite conditional repair and candidate equalizers exist. Weighted averaging factors the supplied scalar stiffness exactly. A current-record lift couples native C to the existing six-level response at both seam endpoints, preserving its full load marginal and disjoint-seam locality. N4 instead conjugates binary C/E by one coherent load-response exchange and supplies private-state influence on fixed bare occupation probabilities. At the declared full tracial reference, N4 now supplies a shared bare-occupation geometry and two-attempt kinetics protocol; its finite memory excess rules out native heat closure at the tested couplings.

**Remaining input:** For the conditional-resampling branch, choose one state space, faithful reference and visible fibres. The conservative canonical extension and dissipative repair are different dynamics; a faithful thermal reference is not mandatory for every pure-state candidate. The full-model candidate must link its repair marginal, reversible port response and record clock through shared states and interventions; the scalar benchmark alone does not supply this. This lift adds local response factors and a supplied coupling, is one-way at the repair marginal, and does not establish the complete response grammar on the enlarged source. N4 changes the load marginal and current-label meaning; its supplied coupling and the physical admissibility of the shared bare-occupation preparation/readout require independent justification.

**Closure condition:** Exhibit one local process and explicit restriction/intervention maps supplying the intended repair, carrier-response, reversible-current and clock contracts. Reuse the existing theorem implications once these interfaces hold; derive each extra choice or retain it as an explicit model hypothesis. Separate matched component fixtures do not establish a common source.

Explicit assumption-name navigation: ` faithful_common_reference_law `, ` complete_repaired_visible_datum_fibre `.

Sources:

- [Finite conditional repair](../Lean/Thermodynamics/FiniteConditionalRepair.lean).
- [Unadopted canonical repair-law proposal](../docs/CANONICAL_REPAIR_LAW_RFC.md).
- [Current candidate and native baselines](../docs/MODEL_SEARCH.md).
- [N4 shared bare-occupation response and finite memory certificate](../code/native_response/README.md).

<a id="re-02"></a>

### RE-02: Confluence, acceptance and protected observations

**Constructed:** The normal-form and consensus implications are proved for their explicit repair and transactional contracts; finite protected-memory examples are constructed. The endpoint-current lift inherits C/E completion by seamwise unitary conjugation; positive seam mixtures on its finite binary shell converge to their common fixed algebra after records are discarded, with protected local block populations. The same finite completion and network-positivity arguments transfer to N4 on the binary shell; no new damping mechanism is needed.

**Remaining input:** The required observer menu must supply completeness, local diamonds, well-founded aggregate descent, semantic read closure and valid acceptance. Mere agreement or oscillatory Hamiltonian evolution does not supply termination. This candidate fixed algebra is not identified with the full A2 observer equalizer; its absolute gap and refinement control have not been established. Neither conjugated fixed algebra is identified with the complete physical A2 data; full integer-pile repair is not established for N4.

**Closure condition:** Verify the complete contract on the actual operations and retained observations, including conflicting updates; do not infer it from translated-chart naturality.

Explicit assumption-name navigation: ` repair_completeness `, ` local_diamond `, ` strict_well_founded_transactional_descent `.

Sources:

- [Observer confluence](../Lean/ObservableNormalForms/ObservableNormalForms/ObserverConfluence.lean).
- [Consensus paper](../paper/reality_as_consensus_protocol.tex).
- [Transactional and acceptance assumptions](../claims/assumption_dictionary.md).

<a id="re-03"></a>

### RE-03: Scheduler, reference and resources

**Constructed:** Native counting and bounded scheduling results retain their declared reference and fairness/resource classes.

**Remaining input:** Select the scheduler and complete path reference actually used. Adaptive repair, fixed-node attempts and the canonical Hamiltonian are different branches; their resource hypotheses are not automatically interchangeable.

**Closure condition:** Declare one branch and prove its fairness, locality and attempt/waste bounds on complete histories, charging all operations.

Sources:

- [Canonical scheduler/reference entries](../docs/registers/PREMISE_REGISTER_V3.md).
- [Consensus scheduler contracts](../paper/reality_as_consensus_protocol.tex).

<a id="re-04"></a>

### RE-04: Observer regions, causal provenance and time-slice access

**Constructed:** Finite carriers, authenticated source order and regional reconstruction contracts exist. The canonical executed-history constructor supplies an ancestry rank; the abstract semantic-log interface requires one.

**Remaining input:** Identify the actual joint carrier, regions, correlations, access menu and semantic read support. Local momentum coordinates do not supply momentum access from a two-position-layer menu.

**Closure condition:** Construct restriction-compatible regions and retained histories on the same source. For the scalar regional theorem verify its exterior row-space criterion; local coordinates alone do not establish that criterion.

Existing claim-row navigation (no new edges): ` OPH-SCALAR-REGIONAL-COLLAR-RECONSTRUCTION `.

Sources:

- [Regional scalar reconstruction](../paper/tex_fragments/SOURCE_SCALAR_REGIONAL_TIME_SLICE.tex).
- [Observer and provenance entries](../docs/registers/PREMISE_REGISTER_V3.md).

<a id="re-05"></a>

### RE-05: Selected scalar action and relative capacities

**Constructed:** The golden scalar action is explicit. A connected response admitting a positive diagonal symmetrizer identifies normalized capacities; incidence factorization and the derived rates then need no independent fit.

**Remaining input:** The response, conductances, boundary conditions, mass and absolute action scale remain supplied. Algebraic reconstruction does not select the physical response or turn capacities into physical volumes.

**Closure condition:** Supply an independent source-selection or comparison argument for the action family and its scale; retain symmetrizer positivity and compatibility hypotheses.

Explicit assumption-name navigation: ` declared_golden_tensor_scalar_action_and_compact_instruments `, ` connected_negative_response_action_measure `.

Existing claim-row navigation (no new edges): ` OPH-SOURCE-ACTION-MEASURE-IDENTIFICATION `.

Sources:

- [Declared common scalar action](../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex); exact source anchor: ` eq:common-scalar-action `.
- [Unique normalized symmetrizer](../Lean/Geometry/SourceActionMeasure.lean); exact source anchor: ` normalized_symmetrizer_unique `.
- [Action-measure identification](../paper/tex_fragments/SOURCE_ACTION_DENSITY_IDENTIFICATION.tex); exact source anchor: ` prop:source-action-measure-identification `.

<a id="re-06"></a>

### RE-06: Inertia and the real dynamical completion

**Constructed:** Same-capacity linear reversible dynamics gives the displayed canonical Hamiltonian. Existing action-to-Hamiltonian and field-clock conversions remain proved under their own hypotheses. Three-mode curvature and an exact distance-two acceleration coefficient distinguish the I+epsilon A counterfamily after common clock/mass/speed calibration.

**Remaining input:** Repair capacity equals inertia, velocity reading, linearity, energy and time reversal are candidate premises. Auxiliary oscillator and temporal Maxwell lifts are related alternatives, not consequences of the scalar construction.

**Closure condition:** Derive the kinetic identification from admissible operations or independently constrain dynamical response; exclude the same-equilibrium alternative inertias. Conversion of a supplied action is not its selection.

Sources:

- [Conditional kinetic uniqueness and counterfamily](../docs/MODEL_SEARCH.md).
- [Source-to-Hamiltonian composition](../Lean/Variational/SourceToHamiltonianComposed.lean); exact source anchor: ` source_to_hamiltonian_composed `.
- [Field clock conversion](../Lean/Screen/CurlStoneClockBridge.lean).

<a id="re-07"></a>

### RE-07: Quantum representation and operational probabilities

**Constructed:** The scalar Fock/Weyl construction and constrained local field representation are mathematically compatible. The Hamiltonian is the original positive canonical one.

**Remaining input:** CCR, Hilbert representation, subsystem split, accessible phase effects, channel semantics and Born readout remain declared. Finite-code CCR are approximate with an explicit boundary defect.

**Closure condition:** Realize the specified state/observable/instrument relations on the same source, with controlled approximation. Shared spectra or a classical stationary distribution do not supply this interface.

Sources:

- [Scalar quantum representation](../paper/tex_fragments/SOURCE_SCALAR_QUANTUM.tex).
- [Complete quantum instrument](../paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex); exact source anchor: ` thm:source-scalar-sequential-instrument `.
- [Coherent source assumptions](../docs/research/COHERENT_SOURCE_CLOCKS.md).

<a id="re-08"></a>

### RE-08: Native coherent controls, phase and reset

**Constructed:** Existing proper-code source controls compile individual force/pointer factors. Complete retained-pointer phases are checked in the exact history target. N2/N3 admit a one-event, one-port state/effect discriminator without an external quantum reference or coherent record history.

**Remaining input:** Coherent code gluing, pure blanks, preparation, reset, precision, transport, drift and complete readback remain supplied service requirements. Event-conditioned raw responses in the endpoint-current candidate require a grammar/implementation check beyond the existing proper-code controller contract. Those local preparations and effects are supplied diagnostic resources, not an autonomous hardware execution. N4 uses existing D matrices in a supplied load-response interaction; this does not establish admissibility in the old primitive-response grammar.

**Closure condition:** Execute or bound the full source instrument with its complementary outcomes and all service errors; a requested unitary or accurate coefficient is not a complete channel certificate.

Sources:

- [Coherent code and service contract](../docs/research/COHERENT_SOURCE_CLOCKS.md).
- [Finite pulse and noise theorems](../paper/tex_fragments/SOURCE_SCALAR_FINITE_INSTRUMENT.tex); exact source anchor: ` prop:source-scalar-finite-noise `.

<a id="re-09"></a>

### RE-09: Preparation of the same coupled vacuum

**Constructed:** For the original eight-mode q=3 action, the unchanged local gapped ramp has rigorous full-state error below 0.000433450 at model time 512 and 0.000108323 at 1024. At T=512, local time discretization with N=2^24 adds less than 0.000695044; a fixed d=8192 onsite code adds less than 0.000033144832 through all partial quadratic factors. The ideal finite preparation error is below 0.001162. Existing source controls compile that word with uniform synthesis error at most 1161/2^24; its preparation and subsequent history share one conditional mathematical budget in RE-10. The older finite-code preparation and its 0.059059008 bound remain a separate construction.

**Remaining input:** Independent oscillator blanks, control of the ramp and model time are hypotheses. The native word has not been emitted or executed; physical pulse and code-service precision, storage, initialization and readout costs remain additional. Preparation is not autonomous vacuum selection; no uniform refinement resource theorem is supplied.

**Closure condition:** Implement the native preparation within the same-code history of RE-10 and charge actual initialization, pulse and service errors. Reuse the uniform synthesis estimate and transport the state-specific preparation error by channel contraction. Do not substitute q=3 evidence for q=5/q=233 budgets or a numerical trajectory for a certified bound.

Sources:

- [Bounded local preparation result](../docs/MODEL_SEARCH.md).
- [Preparation error contract](../paper/tex_fragments/SOURCE_SCALAR_LOCAL_PREPARATION.tex); exact source anchor: ` Preparation with an active free field `.

<a id="re-10"></a>

### RE-10: One encoded preparation, evolution and record history

**Constructed:** The original five q=3 force/retained-pointer histories and coupled-vacuum preparation compose on one d=8192 onsite code with 104 field and two pointer work qubits. Time, displaced occupation cutoffs and existing-source synthesis have a combined mathematical trace-distance bound below 0.001237, preserved by the same complete readout. All relative phases and outcomes are retained. The finite recipe has a conservative evolution pulse upper bound below 2^116, with fewer than 2^117 pack/export events; neither word emission nor native execution was performed.

**Remaining input:** Physical blank/pointer initialization, pulse/control and coherent code-service accuracy, storage/drift, routing/reset, complete readback and their resource budgets remain supplied or unimplemented. The enormous constructive upper count does not establish efficient execution. Physical clock selection and refinement-uniform resources remain separate obligations.

**Closure condition:** Emit and execute an adequate native implementation and add reference-stable physical initialization, pulse, service, storage and readout errors to the constructed mathematical budget. Preserve complete records and the comparison order: uniform finite compiler on actual code, preparation error under a common extended finite history, state-specific ideal Gaussian history estimates, then complete readout. A common-code mathematical recipe does not close the physical experiment.

Existing claim-row navigation (no new edges): ` OPH-SCALAR-SEQUENTIAL-INSTRUMENT `.

Sources:

- [One composed preparation–history error and resource budget](../docs/MODEL_SEARCH.md).
- [Exact finite pulse](../paper/tex_fragments/SOURCE_SCALAR_FINITE_INSTRUMENT.tex); exact source anchor: ` prop:source-scalar-finite-pulse `.
- [Sequential composition](../paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex); exact source anchor: ` thm:source-scalar-sequential-instrument `.
- [Fixed-action local splitting](../paper/tex_fragments/SOURCE_SCALAR_TIME_REFINEMENT.tex); exact source anchor: ` prop:source-scalar-time-refinement `.

<a id="re-11"></a>

### RE-11: Controlled refinement of the selected quantum family

**Constructed:** The existing common scalar theorem gives the stated massive 3D free-field Weyl-product limit. The canonical field join preserves the same action, vacuum and observables.

**Remaining input:** Native preparation/control/record errors must remain controlled along that family. The archived spherical routing and interacting many-body limits are different targets. Vanishing kinetic perturbations can share the narrower position-Weyl limit, so that limit alone does not select unique microscopic inertia.

**Closure condition:** Prove uniform or summable errors and compatible states/observables for the actual implementation; reuse the existing scalar limit rather than reprove it. Interactions require their own continuum theorem.

Existing claim-row navigation (no new edges): ` OPH-SCALAR-SPACE-TIME-DETECTOR-LIMIT `.

Sources:

- [Existing scalar continuum](../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex); exact source anchor: ` prop:common-scalar-continuum-limit `.
- [Joint observable refinement](../paper/tex_fragments/SOURCE_SCALAR_TIME_REFINEMENT.tex).
- [Uniform refinement hypotheses](../claims/assumption_dictionary.md).

<a id="re-12"></a>

### RE-12: Source clock, physical time and execution cost

**Constructed:** Canonical model time, supplied event-count clock comparisons and conditional scalar action-clock reconstruction exist. The endpoint-current lift preserves the native load marginal and classical record blocks and therefore its scoped geometry/count relations; this adds no physical clock calibration. For N4, bare-occupation measurements at the endpoints of an unmodified one- or two-event experiment recover the native port Dirichlet geometry up to a positive scalar. Certified two-attempt memory at theta=1/8 and 1/16 rejects native heat dynamics even after a clock-rate fit.

**Remaining input:** Choose a physical clock on the same process and account for preparation, service, readout and ancestry costs. A source scale identification and conversion to laboratory units are distinct. The record clock must be read from the histories that also supply the carrier response and reversible source, with an explicit relation to scalar model time. N4 does not preserve C: its record j is conjugate to a seam-dependent transformed occupation. The new finite occupation protocol does not restore the native current-label identities or calibrate a physical clock; those require compatible physical current and time readouts.

**Closure condition:** Prove the same-history clock/observable intertwiner and its domain, then calibrate independently. A local generator alone does not establish a microscopic causal cone.

Sources:

- [Source count-clock theorem](../paper/tex_fragments/SOURCE_COUNT_CLOCK.tex); exact source anchor: ` thm:source-count-clock `.
- [Scalar action clock](../paper/tex_fragments/SOURCE_SCALAR_CLOCK.tex).
- [Conditional clock conversion](../Lean/Screen/CurlStoneClockBridge.lean).
- [N4 shared bare-occupation response and finite memory certificate](../code/native_response/README.md).

<a id="re-13"></a>

### RE-13: Physical geometry and responsive metric

**Constructed:** Finite source geometry, golden coordinates and conditional Lorentz/spacetime constructions exist. The canonical source includes its adopted oriented boundary and complete-response clauses.

**Remaining input:** Choose the physical region/metric/area/volume attachment and a compatible source-causal continuum. A frozen chart cannot itself provide changing physical curvature.

**Closure condition:** Construct operational metric and geometric response on the same source and refinement family. Preserve the distinction between action capacity, event-count measure and physical volume.

Sources:

- [Native geometric-source boundary](../paper/tex_fragments/NATIVE_GEOMETRIC_SOURCE_IDENTIFIABILITY.tex).
- [Adopted source clauses](../docs/AXIOM_REFERENCE.md).
- [Continuum and metric interfaces](../docs/registers/GRAVITY_PREMISE_LADDER.md).

<a id="re-14"></a>

### RE-14: Gauge algebra to physical gauge fields

**Constructed:** The finite local Standard Model Lie-algebra result remains proved under the adopted A1 complete-response/closure clauses and A2 endogenous transport. Compact-Lie classification is an explicit mathematics import.

**Remaining input:** The physical current/action, same-source category, global form and coupling attachment remain separate. The scalar field construction does not instantiate a physical gauge theory. Compatibility with the twelve-port reversible response is an immediate law-selection constraint, even though constructing physical gauge fields is downstream. Scalar seam equalizers cannot serve as that current: their tested Lie closure is so(11), not the required twelve-dimensional algebra. Conserved difference drives in the endpoint-current lift omit the central relative-phase generator; the full twelve-direction response is a supplied interface, not derived physical gauge dynamics. A carrier-covariant lower-block deletion retains N2 locality, completion and protected populations while failing injective twelve-direction response; those successful properties cannot replace the gauge theorem premise. The original direct-sum source algebra and the candidate tensor factors require an explicit interface.

**Closure condition:** Bind the proved algebra and representation to actual source observables, dynamics and refinement; supply the categorical and physical gauge interfaces used by the selected branch.

Sources:

- [Finite gauge structure and physical attachments](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex); exact source anchor: ` thm:current-algebra `.
- [Gauge source and global-form entries](../docs/registers/PREMISE_REGISTER_V3.md).

<a id="re-15"></a>

### RE-15: Physical edge factors and collar recovery

**Constructed:** Finite charge centers, regional completion and nontrivial bulk collar splits are constructed. Existing exact Markov and conditional recovery theorems remain available.

**Remaining input:** The physical edge-center split and state alignment, edge-reference factorization and log-d normalization must agree. Approximate recovery needs strong conditional mixing and the sharp collar rate on one family. These are alternative recovery routes. Strong matrix-mixing gives direct log-defect control without a universal spectral floor; the physical EC-replacement modulus remains separate.

**Closure condition:** Realize the selected physical split and prove either its exact central-interface condition or its complete uniform approximate-recovery hypotheses. Charge retention alone is insufficient.

Explicit assumption-name navigation: ` central_interface_clause `, ` EC_aligned_Markov_replacement_control `, ` strong_conditional_gibbs_mixing `, ` sharp_collar_rate_receipt `, ` edge_normalization_log_d `.

Existing claim-row navigation (no new edges): ` OPH-GR-D4D-ENTROPY-BRIDGE `.

Sources:

- [Gravity rungs3–5](../docs/registers/GRAVITY_PREMISE_LADDER.md).
- [EC alignment, mixing and edge normalization](../claims/assumption_dictionary.md).

<a id="re-16"></a>

### RE-16: Modular flow, null translations and stress

**Constructed:** Typed modular/half-sided inclusion and conditional local-stress implications are available; finite cap diagnostics are constructed.

**Remaining input:** Support-visible BW identification, common null translation, effective locality, tomography and the same-source physical stress attachment remain required.

**Closure condition:** Instantiate the modular and stress interfaces on the candidate with the required common domains, locality and normalization; algebraic modular flow alone is not physical time or stress.

Explicit assumption-name navigation: ` null_stress_bridge `.

Sources:

- [Gravity rungs7–8](../docs/registers/GRAVITY_PREMISE_LADDER.md).
- [Modular and stress construction](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex).

<a id="re-17"></a>

### RE-17: Joint state–geometry entropy stationarity

**Constructed:** Finite counting stationarity and fixed-cap MaxEnt/envelope identities are proved under their declared variations.

**Remaining input:** The candidate must supply constraint-tangent generalized-entropy stationarity for joint state and geometry variations, with the physical area/edge identification.

**Closure condition:** Construct the permitted variation family and show the full generalized entropy has zero first variation at fixed physical volume. A fixed-reference entropy optimum alone does not do this.

Explicit assumption-name navigation: ` fixed_cap_MaxEnt_variation_family `.

Sources:

- [Gravity rung6](../docs/registers/GRAVITY_PREMISE_LADDER.md).
- [Einstein derivation scope](../paper/tex_fragments/EINSTEIN_DERIVATION_SCOPE.tex).
- [Stationarity and fixed-cap assumptions](../claims/assumption_dictionary.md).

<a id="re-18"></a>

### RE-18: One gravitational branch and its normalization

**Constructed:** The conditional finite-to-Einstein composition and its countermodels identify the required interfaces.

**Remaining input:** One inhabited physical tower must supply small-ball geometry, conserved stress, Ward/Bianchi transport, universal coupling, absolute vacuum baseline and common normalization. All required small-ball remainders must vanish at the prescribed rate on that one family.

**Closure condition:** Exhibit all consumed receipts on one family, including compatible controls and limits. Separate successful examples of individual premises do not establish joint satisfiability.

Explicit assumption-name navigation: ` einstein_typed_common_domain_tower `, ` universal_coupling_receipt `, ` independent_physical_scale_receipt `, ` vacuum_reference_receipt `.

Existing claim-row navigation (no new edges): ` OPH-GR-E2E-BRANCH-ENTRY `.

Sources:

- [Gravity rungs8–11](../docs/registers/GRAVITY_PREMISE_LADDER.md).
- [Flagship gravitational interface](../flagship/from_observer_consensus_to_standard_physics.tex).
- [Common tower and universal coupling](../claims/assumption_dictionary.md).

<a id="re-19"></a>

### RE-19: Matter, statistics and interactions

**Constructed:** Finite carrier/representation and declared Higgs/Yukawa packets exist; the active scalar candidate has a conditional free Bose limit.

**Remaining input:** Select physical sectors, spin/chirality/statistics, kinetic attachments, interactions, coupling coefficients and stable poles on the same process. A declared finite grammar is not its physical selection.

**Closure condition:** Construct the selected interacting source family and its controlled observables/continuum, reusing the algebraic packets without treating supplied matter content as derived.

Sources:

- [Matter and QFT interfaces](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex).
- [Matter selection and action hypotheses](../claims/assumption_dictionary.md).

<a id="re-20"></a>

### RE-20: Parameter selection and P

**Constructed:** Exact conditional identities, normalized-capacity reconstruction and selection boundaries are recorded in the existing ledgers.

**Remaining input:** P and N remain proposed fundamental numerical inputs. Closure maps, same-quantity readback, Koide balance and screen/electroweak attachment are not supplied by the minimal repair experiment.

**Closure condition:** Give a target-independent source selection and same-quantity observable bridge for each claimed value. Do not identify P with a detuning amplitude, temperature or rate without that map.

Sources:

- [Selection ledger](../docs/registers/SELECTION_LEDGER.md).
- [Parameter and comparison premises](../docs/registers/PREMISE_REGISTER_V3.md).
- [Frozen comparison boundaries](../docs/FROZEN_PREDICTION_LADDER.md).

<a id="re-21"></a>

### RE-21: Cosmology and dark-sector attachment

**Constructed:** Conditional radial, horizon, cosmological and dark-source calculations remain available outside the active free-scalar target.

**Remaining input:** Physical expansion/source dynamics, horizon charges, EOS, stress and transfer to observables remain branch-specific. Euclidean vacuum selection needs target-independent source slab data, reflection-positive reconstruction and compatible transfer/refinement; positive transfer alone does not select a vacuum; cosmological transport has additional field and collision inputs.

**Closure condition:** Construct the common source/stress/scale history and each claimed observational transfer. Fixed-history selection alone cannot determine volume response or the equation of state.

Existing claim-row navigation (no new edges): ` OPH-COSMO-VACUUM-ENSEMBLE-CLASSIFICATION `, ` OPH-COSMO-BOLTZMANN-TRANSPORT-REDUCTION `, ` OPH-DARK-SOURCE-DYNAMICS-LENSING `.

Sources:

- [Equation-of-state selection boundary](../Lean/Thermodynamics/EquationOfState.lean).
- [Cosmological vacuum hypotheses](../cosmology/oph_cosmological_vacuum_and_structure_formation.tex).
- [Boltzmann reduction inputs](../cosmology/oph_boltzmann_transport_derivation.tex).
- [Dark source and lensing](../docs/research/DARK_SOURCE_DYNAMICS_LENSING.md).

<a id="re-22"></a>

### RE-22: Physical propagation and an eligible test

**Constructed:** Finite seam symbols, exact coefficients and conditional frozen comparison contracts exist. Internal scalar discriminators are not physical measurements. A local N2/N3 probe and N4 private-preparation/load-response test give finite internal discriminators under supplied quantum preparations and coupling, without empirical data.

**Remaining input:** Attach the same source operator, sector, physical scale, frame/boost, refinement, nuisance model and calibrated readout. Exclusivity and a positive lower scale bound are needed where the registered test requires them.

**Closure condition:** Freeze an unavoidable relation and rejection rule with full remaining freedom, baseline and uncertainties before an independent eligible observation; preserve the frozen ladder unchanged.

Sources:

- [Existing frozen targets and eligibility](../docs/FROZEN_PREDICTION_LADDER.md).
- [Propagation attachment premises](../docs/registers/PREMISE_REGISTER_V3.md).

<a id="re-23"></a>

### RE-23: External calibrations and empirical imports

**Constructed:** The canonical register explicitly distinguishes hadronic empirical payloads and laboratory clock/energy anchors from source deductions.

**Remaining input:** Every result using these inputs retains their provenance, uncertainty and transfer contract. Fundamental-input counts do not absorb or erase them.

**Closure condition:** Supply admissible independent calibration/import data, or replace a payload with an independently justified source computation; do not count calibration as a theorem proof.

Sources:

- [Empirical imports PR-14/15](../docs/registers/PREMISE_REGISTER_V3.md).
- [Hadronic import policy](../docs/policies/HADRON.md).

<a id="re-24"></a>

### RE-24: Baryon-specific branch

**Constructed:** Declared charge labels and finite operator-eligibility/contraction calculations exist.

**Remaining input:** Physical baryon/proton states, baryon-violating effective action and decay readout are separate inputs, outside the active scalar target.

**Closure condition:** Instantiate the same physical state, dynamics and observable contract before using the finite operator result as a physical prediction.

Sources:

- [Baryon entries](../docs/registers/PREMISE_REGISTER_V3.md).
- [Baryon and proton attachments](../claims/assumption_dictionary.md).

<a id="re-25"></a>

### RE-25: Low-temperature and common-reference branch

**Constructed:** Finite Gibbs/repair implications and explicit common-reference obstructions are available. The active pure-vacuum preparation is a different construction.

**Remaining input:** A continuum low-temperature claim needs its common reference, physical energy/temperature identification and refinement-uniform spectral-tail bound; collar recovery uses its separately stated exact or mixing hypotheses. A fitted finite stationary state does not supply these inputs.

**Closure condition:** Exhibit the required common objects and uniform bounds or explicitly choose another branch; retain negative controls rather than silently changing the reference.

Sources:

- [Common object interfaces](../Lean/Thermodynamics/CommonObjectBinding.lean).
- [Physical Gibbs energy identification](../Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean).
- [Refinement-uniform thermal premise](../docs/registers/PREMISE_REGISTER_V3.md).

<a id="re-26"></a>

### RE-26: Common-model compatibility and remaining named inputs

**Constructed:** This crosswalk covers the canonical PR inventory and identifies actual scalar theorem reuse. The authoritative claim registry separately records ordinary assumptions, imports and explicit non-consumers. N4 now joins equilibrium port geometry and a nontrivial kinetics discriminator on the same bare occupations, preparation and events, without adding state factors or a new interaction.

**Remaining input:** The whole paper/theorem body has not been exhaustively normalized into one semantic dependency graph. Ordinary assumptions and imported-result names are inventoried, not individually discharged by this crosswalk. The endpoint-current candidate declares twelve local six-level factors, conditional response operations, retained records and a coupling theta; its flat alternative shows that preservation and convergence alone do not select that coupling law. N4 reuses the same factors but adds a fixed choice of response quadratures; writing it as one exchange does not derive that choice. OPH physical correctness is a hypothesis under investigation, not a consequence of these finite compatibility checks.

**Closure condition:** For each promoted claim bind its actual assumptions and imports, all transitive antecedents and model-added hypotheses to compatible objects with scoped evidence. Keep unreviewed bindings visible and reject circular discharge.

Sources:

- [Per-claim assumptions and edge classifications](../claims/claim_registry.yaml).
- [Named assumption meanings](../claims/assumption_dictionary.md).
- [Canonical register scope](../docs/registers/PREMISE_REGISTER_V3.md).
- [N4 shared bare-occupation response and finite memory certificate](../code/native_response/README.md).
