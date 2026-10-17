# OPH reverse-engineering table: model and theorem premises

As of 2026-10-17. Generated from [the navigation map](../tracking/premise_register.json)
by `python tools/build_model_theorem_map.py`. Edit the input, then regenerate.

This map connects the canonical premise identifiers and active model obligations to their current construction evidence. It is navigation coverage, not an extraction of every paper or Lean antecedent and not a proof of discharge. Alternative branches retain their own hypotheses; cross-references do not change claim-registry dependency roles.

**Construct the simplest shared repair process satisfying the declared observational contract.** Architecture, observer services and the complete repair law can supply inputs. Their consequences need not be postulated again. A theorem attachment must identify the same states, instruments, retained records, access rules and scheduling domain; matching component formulas alone is insufficient. Unique selection of an implementation is separate from executing that implementation and testing its predictions.

The live construction is `model/`. Its native C/E instruments determine local exchange and coherence laws, finite stationary preparation, field and event-clock statistics, and complete quantum observation bounds. The finite transfer logarithm, local Liouville generator and modular action on a chosen algebra/reference are distinct constructions. An observer-relative clock is an accessible record process; no preferred universal physical clock is required. Supplied scalar, gauge and matter actions remain conditional benchmarks with their own occupation, temporal and spatial error contracts.

**Completion concerns observations.** The intended effective quantum-field, classical and gravitational descriptions must reproduce compatible joint statistics and intervention responses within controlled regimes. Fields, metrics and manifolds need not be fundamental objects. Native source execution does not by itself establish physical spacetime or identify a laboratory apparatus. Public oscillator comparisons retain their calibration, postselection, coordinate and sampling limitations; numerical agreement, mathematical consistency and code fidelity are separate from physical identification.

The entries distinguish finite executions, analytic bounds, imported theorem hypotheses and remaining physical attachments. Separate actions, preparations, detectors and approximation domains keep separate certificates. Reusing a theorem requires a valid map of its actual hypotheses, not a count of available results or the simultaneous imposition of incompatible branches.

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
| Registered claims | 339 |
| Unique ordinary `assumptions` names | 818 |
| Unique `imported_results` names | 896 |
| Workgroups | 27 |
| Assumption names explicitly linked by optional group fields | 33 |
| Claims explicitly linked by optional group fields | 18 |

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

Existing claim classifications: `explicit_edges`: 171; `explicit_non_consumer`: 168.

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
| [RE-27: Executed common scalar, gauge, matter and record process](#re-27) | M01, M02, M04, M05, M06, M10 | active | conditional | None |

<a id="re-01"></a>

### RE-01: Common repair reference and law

**Constructed:** Finite conditional repair and candidate equalizers exist. Weighted averaging factors the supplied scalar stiffness exactly. A current-record lift couples native C to the existing six-level response at both seam endpoints, preserving its full load marginal and disjoint-seam locality. N4 instead conjugates binary C/E by one coherent load-response exchange and supplies private-state influence on fixed bare occupation probabilities. At the declared full tracial reference, N4 supplies a shared bare-occupation geometry and two-attempt kinetics protocol; its finite memory excess rules out native heat closure at the tested couplings. The executable N1 observer composes current observation, retained classical history, record-only prediction and feedback for two events, with all benchmark paths retained. Its detector extension maps the accepted closure P through the existing electromagnetic coupling to a binary record and native feedback. A minimal twelve-block address carrier preserves the complete quantum channel and joins binary native/N4 repair with collective response in an abstract scalar-seam A1–A3 source. The two laws share the full quantum reference and normalized geometry but differ in certified finite kinetics. A separate coherent phase/readout protocol distinguishes C and E exactly while preserving their common classical load transitions; it supplies a candidate comparison without requiring a uniqueness theorem. The common-source joint instrument lets the actual first native current prepare the gauge/scalar experiment and lets its final scalar record control the second native repair, retaining all 1056 histories and their quantum terminal states. In that earlier classical attachment, native M78 integer repair remains a supplied primitive while the quantum apparatus uses the existing M6 source controls. The independent uniform-seam forgotten-record C/E channels have exact quantum gap (5-sqrt(5))/60 on both the binary shell and the full total-two carrier, by an exact signed-difference reduction. Their sharp binary coherence rates differ despite identical classical geometry and relaxation. The coherent native attachment places source-compiled C or E repairs before and after the same retained four-read work program. Its first decoded current controls the gauge kick and its second smooth bit selects the final native seam. The complete joint channel admits arbitrary native/work/reference inputs; default executions use the declared scalar vacuum. This supplies a common finite interface without selecting a unique native record law or a physical encoding. A separate coherent-continuation experiment checks all 3960 native binary arcs and a bijective unread-record reversal. With supplied uniform coherent event embedding A and reflection, W=S(2AA*-I) has compression A*W^n A=T_n(K), and the actual centered load columns obey KX=X(I-L/60). This gives a compressed-amplitude wave recurrence, not a mean-load or probability law. Two continuation phases share the same first coherent event but have two-step Ramsey plus probabilities 1 and 113/150 for the specified basis input. Record-label phases invisible to complete committed C/E CP branches change the compressed spectral trace by -1/900. Thus the native instrument does not select the coherent extension; collision states obstruct this same-domain arc reversal on the full carrier. Coherent preparation, unread records, reflection and interferometric readout are additional controls; no source pulse realization, spatial refinement or physical clock is claimed. Proper-carrier covariance and the declared real coherent-continuation criteria leave the two executed alternatives. A source-executed record-phase experiment has the same complete committed C/E instrument but a fixed coherent read with probabilities 3/4 versus 1/4; commitment before read gives 1/2 for both.

**Remaining input:** For the conditional-resampling branch, choose one state space, faithful reference and visible fibres. The conservative canonical extension and dissipative repair are different dynamics; a faithful thermal reference is not mandatory for every pure-state candidate. The full-model candidate must link its repair marginal, reversible port response and record clock through shared states and interventions; the scalar benchmark alone does not supply this. The one-way current-record lift adds local response factors and a supplied coupling. For binary N4, collective controls supply a complete abstract response grammar on the central-address carrier; independent-factor or address-controlled reversible probes are excluded from that grammar. N4 changes the load marginal and current-label meaning; its supplied coupling and the physical admissibility of the shared bare-occupation preparation/readout require independent justification. Exact all-input preservation of the original reduced load channel forbids unconditional internal-to-load feedback. The retained rational CPTP control has feedback difference 144/625 and a scalar-only private commutant, so the restricted single-product exclusion cannot be extended to all number-preserving conjugations (#1089). Reversible words and chart paths are fixed independently of measured records; admitted feedback uses the repair/reset templates. Unrestricted record-conditioned reversible words add 144 independent address-controlled response directions and are excluded from this abstract grammar. The coherent native attachment implements its declared finite encoding; source compilation alone supplies neither a physical native-load interface nor independently calibrated preparation, code-transfer, storage and readout errors. A source-selection argument must constrain cross-record phases and coherent continuation, or identify observables independent of those choices. Selection is relative to the admitted observer interface: scalar branch phases invisible throughout committed histories need no separate operational selection there, while coherent access makes their full implementation part of the repair-rule obligation.

**Closure condition:** Exhibit one local process and explicit restriction/intervention maps supplying the intended repair, carrier-response, reversible-current and clock contracts. Reuse the existing theorem implications once these interfaces hold; derive each extra choice or retain it as an explicit model hypothesis. Separate matched component fixtures do not establish a common source.

Explicit assumption-name navigation: ` faithful_common_reference_law `, ` complete_repaired_visible_datum_fibre `, ` declared_uniform_independent_seam_schedule_and_forgotten_records `, ` declared_common_source_coherent_native_attachment `.

Existing claim-row navigation (no new edges): ` OPH-MODEL-NATIVE-OBSERVER-EXECUTION `, ` OPH-MODEL-P-RESPONSE-EXECUTION `, ` OPH-MODEL-COMMON-SOURCE-EXECUTION `.

Sources:

- [Finite conditional repair](../Lean/Thermodynamics/FiniteConditionalRepair.lean).
- [Unadopted canonical repair-law proposal](../docs/CANONICAL_REPAIR_LAW_RFC.md).
- [Current candidate and native baselines](../model/CONSTRUCTION.md).
- [N4 shared bare-occupation response and finite memory certificate](../code/native_response/README.md).
- [Executable native observer and finite record/feedback contract](../model/observer/README.md).
- [Native C and exact observation-record-feedback execution](../model/observer/native.py); exact source anchor: ` def run() `.
- [P-dependent detector and bit-controlled continuation](../model/observer/p_response/experiment.py); exact source anchor: ` def controlled_native_step( `.
- [Exact CPTP feedback escape and single-product boundary](../code/programme_obstruction_audit/channel_independent.py).
- [Central address refresh, incident-star repair and decoding](../code/native_response/carrier.py).
- [Exact recorded C/E instruments and common phase/readout probabilities](../model/predictions/observer_instrument.py).
- [Actual native-current preparation and final-bit continuation](../model/common_source/native_join.py).
- [Exact quantum gap, coherence rates and complete-carrier memory boundary](../model/observer/QUANTUM_MIXING.md).
- [Complete coherent native/work composition and retained physical boundaries](../model/common_source/QUANTUM_NATIVE.md).
- [Source-compiled native attachment and complete-output replay](../model/common_source/quantum_native_results.json).
- [Native coherent attachment and continuation nonselection](../model/observer/COHERENT_CONTINUATION.md).
- [Exact binary arc and phase countercontrol receipt](../model/observer/results/coherent_continuation.json).
- [Three primitive selection obligations and interface-equivalence tests](../model/SOURCE_SELECTION.md).

The [work-observer identifiability control](../model/identifiability_checks/README.md#what-a-current-controlled-work-observer-can-distinguish)
retains exact native diagonal-block equivalence for C/E, including initial
native/work/reference correlations and common-record feedback in both
directions. Work-only instruments cannot expose the discarded native
coherences. This narrows the physical comparison interface; it does not
discharge repair-law selection or the native-to-laboratory identification.

<a id="re-02"></a>

### RE-02: Confluence, acceptance and protected observations

**Constructed:** The normal-form and consensus implications are proved for their explicit repair and transactional contracts; finite protected-memory examples are constructed. The endpoint-current lift inherits C/E completion by seamwise unitary conjugation; positive seam mixtures on its finite binary shell converge to their common fixed algebra after records are discarded, with protected local block populations. The same finite completion and network-positivity arguments transfer to N4 on the binary shell; no new damping mechanism is needed. The finite observer explicitly preserves committed classical record values through feedback while allowing disturbance of coherence between record values. The common-source record adapter instantiates exact postcommit checkpoint projector preservation, block and classical public conditional expectations, and all-input controlled Kraus intertwiners. The distinction between a new write on a blank tape and subsequent preservation is retained; no universal confluence or full observer-federation interface follows. The same finite source apparatus has a twelve-node acyclic compatibility system defined before its local update rule. Parent-to-child repair reaches its unique complete history under every node-fair schedule from arbitrary initial assignments; five covering rounds suffice. Ideal contractivity gives a strictly decreasing real residual potential, while finite recovery uses the separate depth/fairness induction. All complete terminal vectors and defined deterministic observer extractions are independent of the solver schedule. Complete classical native histories have an exact partial-permutation/projector form: every initial load population is recoverable from final source and tape, while classical commitment removes cross-history coherences. Record-only E forests identify the full or binary input in eight or seven attempts; a C stopping controller does so almost surely with mean 24 or 21. Every unfinished prefix retains its raw probability. Bounded record-only C policies retain an entangled branch on this carrier; this is compatible with the doubled-forest exclusion of every passive exact code. The coherent attachment preserves committed record-prefix probabilities while retaining the native and work quantum output within each branch. It composes the complete C/E source instruments, including every raw padding and abort outcome; a first-event abort halts with the prepared field retained. This is an instrument-composition result, not new confluence, full-state recovery after classical commitment or a native global-history solver. A separate active encoded variant compiles recovery after each native event while preserving the committed tape and the continuing work state. Its logical-code branch identity applies with arbitrary references; a matched passive schedule demonstrates that retaining records alone does not preserve the logical phase.

**Remaining input:** The required observer menu must supply completeness, local diamonds, well-founded aggregate descent, semantic read closure and valid acceptance. Mere agreement or oscillatory Hamiltonian evolution does not supply termination. This candidate fixed algebra is not identified with the full A2 observer equalizer; its absolute gap and refinement control have not been established. Neither conjugated fixed algebra is identified with the complete physical A2 data; full integer-pile repair is not established for N4. This fixed-input acyclic solver neither implements global history repair on native M78 hardware nor discharges the integer BareConsensusTower, unrestricted observer menu or arbitrary conflicting-update contracts. It solves the given event graph and action; it does not derive them. The final-only native recovery results require classical commitment and the declared seam policies, and do not apply to an unread coherent tape. The active encoded variant adds intervening quantum controls on a selected two-dimensional code. It establishes neither correction of arbitrary full-carrier states nor compatibility with the passive geometry/clock calibration, and it supplies no physical noise or memory-service certificate.

**Closure condition:** Verify the complete contract on the actual operations and retained observations, including conflicting updates; do not infer it from translated-chart naturality.

Explicit assumption-name navigation: ` repair_completeness `, ` local_diamond `, ` strict_well_founded_transactional_descent `, ` declared_common_history_constraint_repair `, ` declared_native_classical_history_and_final_recovery `, ` declared_common_source_coherent_native_attachment `, ` declared_common_source_active_encoded_native_memory `.

Existing claim-row navigation (no new edges): ` OPH-MODEL-NATIVE-OBSERVER-EXECUTION `, ` OPH-MODEL-COMMON-SOURCE-EXECUTION `.

Sources:

- [Observer confluence](../Lean/ObservableNormalForms/ObservableNormalForms/ObserverConfluence.lean).
- [Consensus paper](../paper/reality_as_consensus_protocol.tex).
- [Transactional and acceptance assumptions](../claims/assumption_dictionary.md).
- [Fixed CPTP observer and classical-record preservation](../model/observer/DERIVATION.md).
- [Independent native observer tests](../model/observer/test_native_observer.py).
- [Exact checkpoint and record-controlled instrument adapters](../model/common_source/records.py).
- [Fixed-history compatibility, local residual descent and fair-schedule proof](../model/common_source/HISTORY_REPAIR.md).
- [Local constraints using the actual source operators and complete quantum outputs](../model/common_source/history_repair.py).
- [Local density-history completion of the supplied anharmonic apparatus and its action control](../model/identifiability_checks/consensus_observability.py).
- [Native/work observational equivalence and the separate number-dephasing obstruction](../model/identifiability_checks/README.md#repair-and-consensus-in-an-empirical-comparison).
- [Complete-record recovery, encoded-qubit countercontrol and exact identifying schedules](../model/observer/RECORD_RECOVERY.md).
- [Source-bound exact history, recovery and stopping replay](../model/observer/results/record_recovery.json).
- [Complete coherent native/work composition and retained physical boundaries](../model/common_source/QUANTUM_NATIVE.md).
- [Source-compiled native attachment and complete-output replay](../model/common_source/quantum_native_results.json).
- [Active encoded recovery, continuing quantum memory and passive loss control](../model/common_source/ENCODED_MEMORY.md).

<a id="re-03"></a>

### RE-03: Scheduler, reference and resources

**Constructed:** Native counting and bounded scheduling results retain their declared reference and fairness/resource classes. The native observer declares its two-event tape and halting rule. The correctable-record evaluator gives 11 reachable terminal words and 16 prefixes under its binary preparation, with a lossless four-bit classical prefix encoding. On one carrier, complete even-load calibration laws reconstruct the unweighted structural response after dividing out each positive seam exposure. This removes uniform scheduling and fair-coin requirements for that structural readout; complete known preparations and full seam exposure remain. Deterministic even-load currents also permit exact reconstruction from 360 executed native shots, one per preparation/seam, without a probability-law oracle. This removes native sampling uncertainty for that protocol, not laboratory preparation or detector errors. A separate current-only binary calibration uses twelve spanning preparations and addressed seams. Its 110 stochastic settings have cap 24 and total timeout probability at most 55/4194304 under independent fair placement draws; 250 constant-current settings need one attempt each. The retained seeded execution uses 578 attempts. Both outcome values, rather than hidden coin records, support exact reconstruction under the declared fair law. Preparation/reset and record access remain inputs; the samples do not establish fairness. The exact marked incident/movement-clock law retains zero-current local hits, remote waits, field-conditioned continuations and feedback hits. For the uniform total-two single-port source, I/K tends to 33/5; an exact Poisson corrector supplies a finite conversion-error bound indexed by fixed global movement count, uniformly over initial states. The retained host record comparison preserves baseline acceptance and workload-control rejection under that conversion uncertainty without refitting. It does not transfer the ratio or calibration to the different root-clock product/canonical-sector field protocol, identify a physical clock, or authenticate the pseudorandom entropy.

**Remaining input:** Select the scheduler and complete path reference actually used. Adaptive repair, fixed-node attempts and the canonical Hamiltonian are different branches; their resource hypotheses are not automatically interchangeable.

**Closure condition:** Declare one branch and prove its fairness, locality and attempt/waste bounds on complete histories, charging all operations.

Existing claim-row navigation (no new edges): ` OPH-MODEL-OBSERVER-RECORD-CAPACITY `, ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Sources:

- [Canonical scheduler/reference entries](../docs/registers/PREMISE_REGISTER_V3.md).
- [Consensus scheduler contracts](../paper/reality_as_consensus_protocol.tex).
- [Finite source witnesses, record capacity and controller encoding](../model/observer/CAPACITY.md).
- [Reachable-word capacity and compact controller execution](../model/observer/capacity.py); exact source anchor: ` def compile_controller( `.
- [Independent capacity checks](../model/observer/test_observer_capacity.py).
- [Single-carrier exposure-normalized geometry and scheduling control](../model/relative_observer_geometry/intrinsic_geometry.py).
- [Current-only binary calibration and exact occupation propagation](../model/relative_observer_geometry/binary_dynamics.py).
- [Finite binary calibration, stopping bound and quantum response receipt](../model/relative_observer_geometry/binary_results.json).
- [Joint incident and movement clocks with finite conversion uncertainty](../model/native_clock/CLOCK_BRIDGE.md).
- [Retained joint clock and host comparison](../model/native_clock/clock_bridge_results.json).

<a id="re-04"></a>

### RE-04: Observer regions, causal provenance and time-slice access

**Constructed:** Finite carriers, authenticated source order and regional reconstruction contracts exist. The canonical executed-history constructor supplies an ancestry rank; the abstract semantic-log interface requires one. The executable native observer specifies a finite source/tape carrier, a record-only controller menu and repeatable committed-prefix reads. The electromagnetic detector distinguishes the one-bit actuator algebra from its complete recorded trajectory. A separate binary native/N4 source admits twelve primitive central address blocks with a unital inclusion of its full quantum working algebra. Refresh followed by incident-star repair decodes to the original channel on arbitrary correlated inputs and repeated histories.

**Remaining input:** Identify the actual joint carrier, regions, correlations, access menu and semantic read support. Local momentum coordinates do not supply momentum access from a two-position-layer menu. This single local protocol has no federation, transported record service or continuum regional net; it does not instantiate the complete abstract OperationalObserver structure. The unextended full M66/M78 working algebra has scalar centre and cannot itself be the twelve-central-port A1 carrier (#1087). The addressed abstract source has scalar inter-carrier seams; its endpoint records require a separate retention service, and it supplies no physical field net.

**Closure condition:** Construct restriction-compatible regions and retained histories on the same source. For the scalar regional theorem verify its exterior row-space criterion; local coordinates alone do not establish that criterion.

Existing claim-row navigation (no new edges): ` OPH-SCALAR-REGIONAL-COLLAR-RECONSTRUCTION `, ` OPH-MODEL-NATIVE-OBSERVER-EXECUTION `.

Sources:

- [Regional scalar reconstruction](../paper/tex_fragments/SOURCE_SCALAR_REGIONAL_TIME_SLICE.tex).
- [Observer and provenance entries](../docs/registers/PREMISE_REGISTER_V3.md).
- [Finite source/tape access contract and record algebra](../model/observer/DERIVATION.md).
- [One-bit actuator versus full audit algebra](../model/observer/p_response/ACCESS_CAPACITY.md).
- [Central-port carrier, abstract source and recording boundary](../code/native_response/README.md); exact source anchor: ` ## A central-port carrier that retains the quantum repair `.
- [Correlated quantum and complete-history decoding tests](../code/native_response/test_carrier.py).
- [Whole-simple-algebra and occupation controls](../code/programme_obstruction_audit/joint_independent.py).

<a id="re-05"></a>

### RE-05: Selected scalar action and relative capacities

**Constructed:** The declared golden scalar action admits the stated capacity and conductance reconstruction under its positive diagonal symmetrizer hypotheses. The two-site oscillator apparatus supplies mass one and spring kappa=9/32 and is a different action family. For native binary repair, the actual mean response fixes the relative quadratic form without a fitted oscillator coefficient. More fully, the C/E endpoint instruments derive H_native=sum_e(I-Phi_e) on operator space, with positive local terms, exclusion exchange and candidate-dependent coherence killing. The exact discrete transfer is I-H_native/M. Independent rate-one Poisson seam attempts give exp(-s H_native), where s is auxiliary expected attempts per seam. The finite transfer logarithm, this local Liouville action and a state-vector Hamiltonian are distinct objects; each has its own observation contract.

**Remaining input:** The supplied scalar benchmarks retain their masses, boundary conditions, response and absolute action scale. Native exchange and coherence coefficients follow once the complete repair instrument and scheduler are specified; they need no separate field-action fit. Selecting the physically applicable instrument, identifying its accessible fields and clock records with measurements, and obtaining relativistic observable laws remain open. Algebraic reconstruction does not turn capacities into physical volumes or the dissipative native action into a unitary Hamiltonian.

**Closure condition:** Supply an independent source-selection or comparison argument for the action family and its scale; retain symmetrizer positivity and compatibility hypotheses.

Explicit assumption-name navigation: ` declared_golden_tensor_scalar_action_and_compact_instruments `, ` connected_negative_response_action_measure `.

Existing claim-row navigation (no new edges): ` OPH-SOURCE-ACTION-MEASURE-IDENTIFICATION `, ` OPH-MODEL-SCALAR-RESPONSE-EXECUTION `, ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Sources:

- [Declared common scalar action](../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex); exact source anchor: ` eq:common-scalar-action `.
- [Unique normalized symmetrizer](../Lean/Geometry/SourceActionMeasure.lean); exact source anchor: ` normalized_symmetrizer_unique `.
- [Action-measure identification](../paper/tex_fragments/SOURCE_ACTION_DENSITY_IDENTIFICATION.tex); exact source anchor: ` prop:source-action-measure-identification `.
- [Executed scalar detector and native feedback contract](../model/observer/scalar_response/README.md).
- [Binary instrument, history law and inherited refinement bound](../model/observer/scalar_response/DERIVATION.md).
- [Retained scalar-response execution and controls](../model/observer/scalar_response/results.json).
- [Declared two-site action with the physical occupation-boundary terms](../model/observer/source_scalar/field.py); exact source anchor: ` def physical_hamiltonian( `.
- [Local quantum action derived from complete native instruments](../model/native_action/README.md).
- [Native hard-core exchange and controlled transfer logarithm](../model/native_effective_fields/README.md).
- [Source attachment of relative capacities and mean action](../model/primitive_observer/README.md).

<a id="re-06"></a>

### RE-06: Inertia and the real dynamical completion

**Constructed:** Same-capacity linear reversible dynamics gives the displayed canonical Hamiltonian. Existing action-to-Hamiltonian and field-clock conversions remain proved under their own hypotheses. Three-mode curvature and an exact distance-two stiffness-response coefficient distinguish the I+epsilon A counterfamily after common clock/mass/speed calibration. For Y=gamma(A+epsilon A^2), Y13=gamma epsilon while actual acceleration is -Yq. The configuration-Gibbs control preserves only the normalized q marginal; its full phase-space reference and quantum vacuum change. A distinct addressed native/N4 counterpair fixes the full quantum reference and collective response while its finite kinetic invariant differs (#1090). The scalar-response adapter propagates the supplied Hamiltonian and original vacuum. Its feedback statistic does not select inertia from native repair or configuration equilibrium. Each native C/E endpoint instrument has an explicit unitary completion and dimensionless Hermitian history action. Keeping hidden branch records reproduces arbitrary finite record-conditioned histories and untouched-reference outputs. The unused completion and logarithm are representation choices: equal native endpoints can have different auxiliary intermediate readouts. This supplies a finite action representation without selecting an autonomous carrier Hamiltonian or physical inertia. The binary calibration recovers the same Gram within the 66-state quantum clock sector. Both C/E native channels give the occupation response T=I-L/60 on every initial density. Under the uniform stationary reference, occupation covariance at attempt lag k is (5/33)T^k Q, with Q=I-11^T/12 the centered port projection. Exact two-hop response and the three contracting nonconstant mode multipliers connect this geometry to native propagation. This direct readout has no phase rotation under positive scalar clock rescaling; other quantum reconstructions are outside that restriction.

**Remaining input:** Repair capacity equals inertia, velocity reading, linearity, energy and time reversal are candidate premises. Auxiliary oscillator and temporal Maxwell lifts are related alternatives, not consequences of the scalar construction. Neither counterpair violates the canonical same-capacity inertia theorem. The fixed-reference pair is an abstract binary source construction, not two physically identified realizations.

**Closure condition:** Derive the kinetic identification from admissible operations or independently constrain dynamical response; exclude the same-equilibrium alternative inertias. Conversion of a supplied action is not its selection.

Existing claim-row navigation (no new edges): ` OPH-MODEL-SCALAR-RESPONSE-EXECUTION `, ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Sources:

- [Conditional kinetic uniqueness and counterfamily](../model/CONSTRUCTION.md).
- [Source-to-Hamiltonian composition](../Lean/Variational/SourceToHamiltonianComposed.lean); exact source anchor: ` source_to_hamiltonian_composed `.
- [Field clock conversion](../Lean/Screen/CurlStoneClockBridge.lean).
- [Configuration-only Gibbs and stiffness-response control](../code/programme_obstruction_audit/kinetic_control.py).
- [Separate fixed-full-reference abstract kinetic counterpair](../code/native_response/README.md); exact source anchor: ` ### A fixed-full-reference counterpair `.
- [Executed scalar detector and native feedback contract](../model/observer/scalar_response/README.md).
- [Binary instrument, history law and inherited refinement bound](../model/observer/scalar_response/DERIVATION.md).
- [Retained scalar-response execution and controls](../model/observer/scalar_response/results.json).
- [Exact native instrument history action and endpoint-equivalent interpolations](../model/native_action/history_action.py).
- [Retained finite unitary history certificate](../model/native_action/history_action_results.json).
- [Current-only binary calibration and exact occupation propagation](../model/relative_observer_geometry/binary_dynamics.py).
- [Finite binary calibration, stopping bound and quantum response receipt](../model/relative_observer_geometry/binary_results.json).

<a id="re-07"></a>

### RE-07: Quantum representation and operational probabilities

**Constructed:** The scalar Fock/Weyl construction and constrained local field representation are mathematically compatible. The Hamiltonian is the original positive canonical one. The finite native observer has an all-input CPTP map with coherent write isometries and an explicit classical readout. The binary electromagnetic detector has complementary effects on arbitrary coherent record inputs; its outcome probabilities and complete continuation are implemented. A separate binary native/N4 carrier has twelve primitive central address blocks and an exactly recoverable full quantum sector. Every instrument reporting effects N_p/2 erases disjoint-load coherences; an explicit pair has unavoidable recovery error at least 1/2. The scalar-response adapter executes a binary Weyl instrument on the evolved scalar field and retains conditional branch states before native continuation. The pointer, quantization and Born readout remain supplied; drawing independent unconditional receipt marginals is not this instrument execution. The two-site source implementation compiles an entangled finite vacuum preparation, signed impulse, scalar propagation and pointer coupling from the existing twelve-port M6 controls. Both pointer outcomes and their complete conditional finite field states are retained before the native bit-controlled continuation. Its physical-square projection includes the d/2 top-level correction; an explicit finite-CCR defect prevents identifying these gates with exact Gaussian shears. A full 78-state exact instrument comparison gives plus probabilities 5/8 and 1/8 for C at phases 0 and pi, versus 3/8 for E; phase pi/2 and diagonal controls agree at 3/8. Both readout outcomes and committed old records are retained. The separate two-site Gaussian-history cutoff bounds complete output states for the declared prepared histories, including unnormalized terminal field branches. The common interacting terminal vectors supply explicit matrix GNS representations for record/work, classical-record and scalar-site access. The exact rank-r formula gives dimension d*r on M_d; numerical spectra and rank thresholds are reported separately. Forgetting records is an explicit access restriction, with no postselection or added noise floor. A separate faithful full-work Gibbs reference is attached to the same finite Hamiltonian. The actual preparation and all recorded quantum outputs retain their original states; faithfulness is obtained by choosing the reference, not by adding noise to those outputs. Positive weights on the record/work blocks leave the classical centers fixed under modular flow. The continued adapter retains conditional Weyl branches through two noncommuting reads of the same field; all modes enter a small Gram kernel and independent final Fock states check the relative phases. A commuting Gaussian control preserves temporal correlation while removing the separately certified disturbance. The source-discrimination construction embeds the full 78-state carrier in nine qubits and intertwines both recorded C/E channels exactly on all 6084 matrix units, including arbitrary reference systems. Its binary pointer realizes the complete Lueders readout, preserving the rank-77 native minus subspace and every old record block. The contrast instrument E equals C followed by one recorded seam-parity bit, implemented by a three-pulse source CNOT. This is a minimal classical refinement of the fixed coin/current labels. Discarding its label does not undo the resulting dephasing; the shared source controls and classical endpoint map leave an operationally distinguishable record-law choice. The same uniform forgotten-record laws are entanglement-breaking on the binary shell by attempt 134. An explicit collision/binary reference-entangled input proves that neither full 78-state channel is entanglement-breaking at any finite attempt count, even though all inputs converge to the rank-66 uniform binary state. This separates the invariant-shell result from the complete carrier and from retained-record or controlled quantum memory. The complete classical-record recovery theorem gives optimal maximally-mixed entanglement fidelity D^(-2) sum_h w_h rank(P_h)^2 and arbitrary-reference iid dephasing bounds. Existing source inverse words restore each history-projected input while retaining its archive, with all raw outcomes and complex phases checked. A concrete four-configuration qubit satisfies exact Knill-Laflamme identities after one C/E event; that component supplies an algebraic final-recovery certificate. The separate active encoded attachment executes its source encoder and record-conditioned recovery, and composes logical identity with the retained work instrument on arbitrary logical/work/reference inputs. A fixed passive (0,1) then (2,3) history loses logical phase and has optimal Bell entanglement fidelity 1/2; active correction after each event preserves the ideal logical qubit. Exact incident-clock field correlations reconstruct a seven-dimensional cyclic detector Hilbert space and its transfer spectrum on the elementary total-two source; independent full 66-state propagation checks withheld lags. The positive-transfer moment class is supplied by the existing theorem. Exact rank is not inferred from sampled records, and the seven-mode Gram is severely ill-conditioned. On each actual level-3/5 screen with product preparation, bounded record correlations also certify a two-vector temporal subspace and compressed transfer values near 0.88702 and 0.97079 with error below 2e-5. These compressed values are not a closed two-mode spectrum or laboratory energy calibration.

**Remaining input:** CCR, Hilbert representation, subsystem split, accessible phase effects, channel semantics and Born readout remain declared. Finite-code CCR are approximate with an explicit boundary defect. These finite protocols use supplied quantum instrumentation and do not derive the Born rule or the scalar Fock representation from the native load law. For this source implementation, coherent proper-code service, reset blanks and Born readout retain their declared status. The two-site ideal occupation limit uses the declared Gaussian action and preparation; source-word execution does not derive canonical quantization or supply a spatial field limit. The C/E discriminator assumes calibrated coherent source preparation and readout; it is not compiled by the separate two-site source word. The separate encoded C/E source construction retains the supplied physical load-to-code identification, coherent services, reset and Born readout; its exact channel intertwining does not calibrate those interfaces. The final-only native recovery results require classical commitment and the declared seam policies, and do not apply to an unread coherent tape. The active encoded variant adds intervening quantum controls on a selected two-dimensional code. It establishes neither correction of arbitrary full-carrier states nor compatibility with the passive geometry/clock calibration, and it supplies no physical noise or memory-service certificate.

**Closure condition:** Realize the specified state/observable/instrument relations on the same source, with controlled approximation. Shared spectra or a classical stationary distribution do not supply this interface.

Explicit assumption-name navigation: ` declared_uniform_independent_seam_schedule_and_forgotten_records `, ` declared_native_classical_history_and_final_recovery `, ` declared_common_source_active_encoded_native_memory `.

Existing claim-row navigation (no new edges): ` OPH-MODEL-NATIVE-OBSERVER-EXECUTION `, ` OPH-MODEL-P-RESPONSE-EXECUTION `, ` OPH-MODEL-SCALAR-RESPONSE-EXECUTION `, ` OPH-MODEL-COMMON-SOURCE-EXECUTION `.

Sources:

- [Scalar quantum representation](../paper/tex_fragments/SOURCE_SCALAR_QUANTUM.tex).
- [Complete quantum instrument](../paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex); exact source anchor: ` thm:source-scalar-sequential-instrument `.
- [Coherent source assumptions](../docs/research/COHERENT_SOURCE_CLOCKS.md).
- [All-input native CPTP construction](../model/observer/DERIVATION.md).
- [Binary detector completeness on coherent record inputs](../model/observer/p_response/DERIVATION.md).
- [Central addresses and the occupied-port information bound](../code/native_response/README.md).
- [Executed scalar detector and native feedback contract](../model/observer/scalar_response/README.md).
- [Binary instrument, history law and inherited refinement bound](../model/observer/scalar_response/DERIVATION.md).
- [Retained scalar-response execution and controls](../model/observer/scalar_response/results.json).
- [Finite occupation representation, binary instrument and complete outcome contract](../model/observer/source_scalar/DERIVATION.md).
- [Independent projected-action, vacuum and conditional-field-state checks](../model/observer/source_scalar/test_experiment.py).
- [Exact C/E discrimination with complete retained readout branches](../model/predictions/observer_instrument.py).
- [Prepared-history instrument trace-distance bound, not an all-input channel norm](../model/observer/source_scalar/cutoff.py).
- [Same-history observer GNS and state pairing](../model/common_source/operational.py).
- [Actual states, full-work thermal reference and retained record algebra](../model/common_source/THERMAL_CLOCK.md).
- [Retained field, record capacity and complete-history proof](../model/observer/scalar_response/CONTINUATION.md).
- [Executed two-read history and outward rational refinement bounds](../model/observer/scalar_response/continuation_results.json).
- [Exact full-carrier C/E source intertwining and arbitrary-reference extension](../model/observer/source_discrimination/certificate.py); exact source anchor: ` def channel_certificate( `.
- [Exact quantum gap, coherence rates and complete-carrier memory boundary](../model/observer/QUANTUM_MIXING.md).
- [Complete-record recovery, encoded-qubit countercontrol and exact identifying schedules](../model/observer/RECORD_RECOVERY.md).
- [Source-bound exact history, recovery and stopping replay](../model/observer/results/record_recovery.json).
- [Active encoded recovery, continuing quantum memory and passive loss control](../model/common_source/ENCODED_MEMORY.md).
- [Detector temporal representation recovered from field correlations](../model/relative_observer_geometry/temporal_records.py).
- [Exact visible modes and certified screen transfer compression](../model/relative_observer_geometry/temporal_results.json).

<a id="re-08"></a>

### RE-08: Native coherent controls, phase and reset

**Constructed:** Existing proper-code source controls compile individual force/pointer factors. Complete retained-pointer phases are checked in the exact history target. N2/N3 admit a one-event, one-port state/effect discriminator without an external quantum reference or coherent record history. The N1-based detector software executes all native and binary-probe outcomes, retains its committed bit and uses only that bit for the next seam. The two-site apparatus emits and executes a composed preparation, field and pointer word. Every one/pair primitive is a full M6 exponential reconstructed from the existing twelve source-control coefficients. Gray rotations and commuting Walsh phases charge all basis changes and scalar phases; zero-angle requests retain their load/export services. The complete d=2,4 export channels retain complement failures and are checked on every matrix unit and reference-entangled inputs. Executed valid-code histories have zero service-failure probability in exact arithmetic. The separate C/E comparison executes its finite quantum instrument in exact Gaussian rationals under declared coherent phase preparation and source-projector readout. Preparation above six qubits uses a binary amplitude tree with commuting Walsh-Y factors through the same source primitives; the smaller retained preparation words are preserved. The eight-field-qubit d=16 preparation emits 5379 native pulses. The existing coherent M6 proper-code controls execute a joint scalar, CAR and electric-flux interaction and both pointer instruments. The abstract source and signed fermionic control constructions predate this execution; this is a composed finite execution using those capabilities. A d=2 continuation executes a source Z probe, modular-gap-dependent X pulse and fresh number pointer on every terminal branch. Its eight retained outcomes preserve the old two-bit marginal and obey the independently checked modular sinusoid. A separate source-compiled smooth Ramsey pointer acts on the same interacting history and retains all four gauge/smooth-detector outcomes. Its exact finite target has effects (I+sin Q0,d)/2 and (I-sin Q0,d)/2; numerical source execution is compared with the full finite Kraus maps, including their phases. Its declared bit dictionary differs from the sharp-sign protocol. Source continuations of every actual post-sharp branch at model times 0, 1/2 and 1 append a gauge-number pointer and retain eight outcomes. Their fixed full-work Gibbs-reference target is the same finite Hamiltonian evolution in modular parameter s=-t/beta; exact target alignment and source-word approximation errors are distinguished. A separate emitted and executed proper-code construction supplies the C/E local dilations, fixed phase preparations, retained fair-coin and diagonal-randomizer branches, and complete binary source-projector readout. All primitive source phases and load/export requests are charged; full-column and reference-entangled tests compare the executed words with independent targets. A phase-exact cyclic Gray parity compiler reduces the full nine-bit controlled-flip macro from 11782 to 2044 source pulses. Exact integer phase accounting and literal full-register executions preserve the endpoint operator, scalar phase, padding and arbitrary-reference action, reducing routing and active-memory resources without adding ancillas or removing source services. The intermediate pulse history changes; this is an endpoint-preserving compiler simplification. Existing C/E source words execute the fixed coherent record-phase read with every outcome retained; all 512 committed branch columns are preserved up to scalar phase. The coherent-read dictionary is distinct from a native current record.

**Remaining input:** Coherent code gluing, pure blanks, preparation, reset, precision, transport, drift and complete readback remain supplied service requirements. Event-conditioned raw responses in the endpoint-current candidate require a grammar/implementation check beyond the existing proper-code controller contract. Those local preparations and effects are supplied diagnostic resources, not an autonomous hardware execution. N4 uses existing D matrices in a supplied load-response interaction; this does not establish admissibility in the old primitive-response grammar. The P-dependent detector's chosen probe preparation, flux exposure, routing and readout remain specified apparatus; that separate detector does not execute the scalar source word or realize the N4 occupation apparatus. The two-site compiler uses the declared coherent proper-code completion; its full export-channel definition is not a hardware realization of transport, reset or error control. The C/E phase/readout instrument is not an emitted proper-code hardware word; its independently bounded apparatus error is a test premise. The programmed modular clock uses supplied coherent controls and a state-derived analogue angle; neither its schedule nor its state estimate is an autonomous clock mechanism or physical service certificate. A compositional trace/channel error contract gives sufficient pulse, service, initialization, record and readout allowances, while compilation, numerical execution and independently established physical error bounds remain separate premises.

**Closure condition:** Execute or bound the full source instrument with its complementary outcomes and all service errors; a requested unitary or accurate coefficient is not a complete channel certificate.

Existing claim-row navigation (no new edges): ` OPH-MODEL-P-RESPONSE-EXECUTION `, ` OPH-MODEL-SCALAR-RESPONSE-EXECUTION `, ` OPH-MODEL-COMMON-SOURCE-EXECUTION `.

Sources:

- [Coherent code and service contract](../docs/research/COHERENT_SOURCE_CLOCKS.md).
- [Finite pulse and noise theorems](../paper/tex_fragments/SOURCE_SCALAR_FINITE_INSTRUMENT.tex); exact source anchor: ` prop:source-scalar-finite-noise `.
- [Complete detector and feedback branches](../model/observer/p_response/experiment.py); exact source anchor: ` def histories( `.
- [Source-response reconstruction and complete proper-code export](../model/observer/source_scalar/compiler.py); exact source anchor: ` def native_pulse( `.
- [Independent source, full-phase gate and reference-stable service tests](../model/observer/source_scalar/test_compiler.py).
- [Declared phase/readout services and independent apparatus-error budget](../model/predictions/protocol.py).
- [Shared source words for coupled finite scalar and charged-matter evolution](../model/common_source/scalar_coupling.py).
- [Actual retained-branch source word for the modular clock](../model/common_source/modular_clock.py).
- [Complete source Ramsey detector on the common interacting carrier](../model/common_source/smooth_protocol.py).
- [Complete source probe continuation for the fixed thermal reference](../model/common_source/thermal_clock.py).
- [Emitted C/E source permutations, fixed preparations and complete pointer word](../model/observer/source_discrimination/circuit.py).
- [Independent full-phase, complete-column and arbitrary-reference source execution tests](../model/observer/source_discrimination/test_circuit.py).
- [Explicit compositional apparatus-error contract with retained failure branches](../model/observer/source_discrimination/DERIVATION.md); exact source anchor: ` ## Compositional apparatus contract `.
- [Full-register phase and source-resource preservation tests](../model/common_source/test_native_gray_synthesis.py).
- [Source-executed coherent-versus-committed record interface](../model/observer/results/source_selection.json).

<a id="re-09"></a>

### RE-09: Preparation of the same coupled vacuum

**Constructed:** For the original eight-mode q=3 action, the unchanged local gapped ramp has rigorous full-state error below 0.000433450 at model time 512 and 0.000108323 at 1024. At T=512, local time discretization with N=2^24 adds less than 0.000695044; a fixed d=8192 onsite code adds less than 0.000033144832 through all partial quadratic factors. The ideal finite preparation error is below 0.001162. Existing source controls compile that word with uniform synthesis error at most 1161/2^24; its preparation and subsequent history share one conditional mathematical budget in RE-10. The older finite-code preparation and its 0.059059008 bound remain a separate construction. For the scalar-response adapter, the first native current sets the signed momentum preparation 4j f on the original finite-action vacuum. This is an executed conditional preparation, with vacuum and admissible controls supplied rather than source selected. For the separate two-site action, the source word prepares an entangled site state from blank registers at onsite dimensions d=4, d=8 and d=16. It retains the antisymmetric normal-mode squeezed vacuum through total occupation d-2 and compiles the resulting finite vector. Directed rational bounds certify its initial trace distance from that action's infinite coupled vacuum. The target vector supplies compiler coefficients, not the executed initial state. The d=8 initial vacuum trace-distance upper bound is 0.000079898200718; it does not include propagation or readout cutoff error. A separate uniform Gaussian number-tail and projected Duhamel calculation adds the impulse, evolution, both pointer arms and final outside component on the same embedded two-site history. The actual-screen field protocol also admits uniform fixed-number preparations: connected native exchanges derive the unique diagonal stationary law within each conserved sector. Exact hypergeometric restrictions, rather than independent product occupations, certify all 44 field/timeout outcomes at totals 5120/15360 and 81920/245760. Positive baseline and negative feedback fourth cumulants survive. Each screen uses its own initial restriction; the complete field-law distance is nonzero and certified near 3.1e-5. Eight paired full-source histories preserve exact totals, native feedback and serialized field-counter continuation. The conserved total and exact sampling service remain inputs; stationarity supplies no finite physical preparation certificate.

**Remaining input:** Independent oscillator blanks, control of the ramp and model time are hypotheses. The original q=3, d=8192 preparation word has not been emitted or executed; physical pulse and code-service precision, storage, initialization and readout costs remain additional. Preparation is not autonomous vacuum selection; no uniform refinement resource theorem is supplied. The initial two-site preparation bound alone controls only a common subsequent channel; its complete-history cutoff certificate supplies the separate finite-versus-infinite propagation and readout comparison. This action and state differ from the original eight-mode ramp, and physical initialization/control errors require additional bounds. The native canonical-sector state is a distinct stationary classical source family, not the supplied interacting oscillator vacuum. Particle density is conserved boundary information, and neither finite waiting resources nor quantum vacuum selection follows from its stationary characterization.

**Closure condition:** Implement the native preparation within the same-code history of RE-10 and charge actual initialization, pulse and service errors. Reuse the uniform synthesis estimate and transport the state-specific preparation error by channel contraction. Do not substitute q=3 evidence for q=5/q=233 budgets or a numerical trajectory for a certified bound. For the two-site branch, use its complete ideal history certificate and add independently justified numerical and physical service errors before claiming a realized infinite-oscillator comparison.

Existing claim-row navigation (no new edges): ` OPH-MODEL-SCALAR-RESPONSE-EXECUTION `.

Sources:

- [Bounded local preparation result](../model/CONSTRUCTION.md).
- [Preparation error contract](../paper/tex_fragments/SOURCE_SCALAR_LOCAL_PREPARATION.tex); exact source anchor: ` Preparation with an active free field `.
- [Executed scalar detector and native feedback contract](../model/observer/scalar_response/README.md).
- [Binary instrument, history law and inherited refinement bound](../model/observer/scalar_response/DERIVATION.md).
- [Retained scalar-response execution and controls](../model/observer/scalar_response/results.json).
- [Two-site coupled-vacuum target and directed rational tail bound](../model/observer/source_scalar/field.py); exact source anchor: ` def vacuum_vector( `.
- [Deterministic source preparation from blank registers](../model/observer/source_scalar/compiler.py); exact source anchor: ` def compile_real_state( `.
- [Preparation overlap and boundary of trace-distance transport](../model/observer/source_scalar/DERIVATION.md).
- [Retained two-site source execution, states, counts and separate error certificates](../model/observer/source_scalar/results.json).
- [Normalized preparation and full projected-history cutoff bound](../model/observer/source_scalar/cutoff.py).
- [Exact canonical-sector field law and source execution](../model/observer_time_fields/canonical_sector.py).
- [Fixed-number full-screen field predictions and finite-size comparison](../model/observer_time_fields/canonical_sector_results.json).
- [Common clock/field proof and conserved-sector extension](../model/observer_time_fields/REFINEMENT.md).

<a id="re-10"></a>

### RE-10: One encoded preparation, evolution and record history

**Constructed:** The original five q=3 force/retained-pointer histories and coupled-vacuum preparation compose on one d=8192 onsite code with 104 field and two pointer work qubits. Time, displaced occupation cutoffs and existing-source synthesis have a combined mathematical trace-distance bound below 0.001237, preserved by the same complete readout. All relative phases and outcomes are retained. The finite recipe has a conservative evolution pulse upper bound below 2^116, with fewer than 2^117 pack/export events; neither word emission nor native execution was performed. Separately, model/observer executes a finite native two-event history and a P-dependent bit-controlled continuation; those are smaller completed protocols with different states and operations. The scalar-response adapter joins a first native repair, signed scalar preparation/evolution, binary field instrument and retained-bit-controlled second native repair. Its conditional field branches and complete native histories execute in software. This does not emit the separate large common-code scalar/pointer word or instantiate its native quantum service schedule. A separate local scalar execution advances the mean and two detector pullback columns by neighboring-site canonical kick-drift-kick recurrences at q=13 and q=34, includes original-vacuum squeezing in its terminal detector law, and samples the same bit-controlled native continuation. A q=13 step-halving control fixes the attained model time. Its retained traces contain checkpoint summaries and native records, not a full authenticated field-event history or a compiled native quantum-control word. A distinct two-site construction emits and executes one source-controlled preparation, signed impulse, 64- or 128-step finite evolution and terminal pointer word at d=4, d=8 and d=16, using five, seven or nine work qubits including the pointer. Repeated routine matrices are contractions of executed full M6 pulses, not substituted target gates. Both conditional field outputs and native feedback histories are retained. Exact rational nested-commutator bounds control ideal split evolution against continuous evolution of the same compressed finite Hamiltonian; the d=8 bound is below 0.000174 at 128 steps. Preparation, finite splitting and numerical compiler checks have distinct error scopes. Uniform Gaussian number moments and projected Duhamel residuals bound the complete ideal two-oscillator instrument, including terminal field branches, by less than 0.077171 at d=8 and 0.000196 at d=16. The d=16 source word executes at both step counts; d=32 and d=64 occupation certificates are analytic comparisons. Adding the finite splitting bound gives the ideal split-versus-infinite comparison for the same prepared history, with combined error below 0.001465 at d=16 and 128 steps. Numerical execution and physical service errors require separate terms. The common-source extension prepares one gauge/scalar register, evolves its noncommuting coupled Hamiltonian, records charged occupation, applies feedback controlled by that committed bit, then records scalar-coordinate sign. It retains all four quantum branches and joins them to the first and second native repair records. No probability is supplied by a separately evolved field process. For zero native preparation current, the same declared interacting action has an exact infinite-oscillator feedback contrast [-0.024045670143, -0.024043450955]. Full-Hamiltonian Taylor powers on expanding support charge initial projection, every support discard and exact coefficient rounding; both first records and both infinite half-line readout outcomes retain their unnormalized quantum states. The final detector uses Hermite half-line integrals rather than the sign of a compressed coordinate matrix. A separate non-Gaussian exponential-number-moment and projected-Duhamel proof bounds the pre-readout occupation refinement uniformly over arbitrary matter/reference inputs with the fixed free coupled scalar vacuum. These are mathematical infinite-history and cutoff certificates, not emitted infinite source words. A completed two-read extension retains one field and the same native loads through three repairs, all 2,112 labelled histories and 209 reached record prefixes, with an explicit eight-bit prefix code, checkpoint restoration and final halt. It uses no intermediate source or field reset. The canonical finite program executes a retained four-read smooth continuation and its native attachment at d=2,4, with all 4224 joint outcomes. Its prefix Kraus identities, resumable simulator state and stagewise moment/Duhamel bound connect repeated smooth measurements to later model-time evolution and a gauge clock probe. Finite temporal, occupation and missing physical-service budgets remain distinct. The coherent native extension executes the complete C/E repairs around that same four-read smooth program. A source-executed exchange of port labels 1 and 5 routes the existing seam (0,1) circuit to (0,5), with its inverse included. The 78 admitted native inputs occupy a 512-dimensional data code; the product-channel comparison retains native/reference coherence and all raw outcomes, including padding and aborts. Resource accounting includes the actual routing and complete attempted repairs. This is a separate finite common-source history, not execution of the original q=3 large-code recipe. The protected-qubit mode inserts source-compiled recovery after the first (0,1) repair and the final repair, selected between (2,3) and (0,1) by the same second smooth bit. Source encoding, recovery and terminal decoding belong to its attempted-path costs; the same encoded-input passive control omits recovery.

**Remaining input:** Physical blank/pointer initialization, pulse/control and coherent code-service accuracy, storage/drift, routing/reset, complete readback and their resource budgets remain supplied or unimplemented. The enormous constructive upper count does not establish efficient execution. Physical clock selection and refinement-uniform resources remain separate obligations. The executed two-site word does not emit or execute the original q=3, d=8192 scalar field/pointer word. The ideal occupation certificate supplies no golden-grid spatial limit, certified floating execution error or physical service-error budget. The earlier free two-site d=8 bound is larger than its detector bias and cannot certify that bias against its infinite reference. The new interacting infinite-history certificate has its own half-line readout and does not certify floating or physical source services, arbitrary scalar-input diamond error, or an infinite apparatus. The coherent native attachment implements its declared finite encoding; source compilation alone supplies neither a physical native-load interface nor independently calibrated preparation, code-transfer, storage and readout errors. Active logical correction is restricted to the supplied code and specified native instrument; external physical noise, raw unknown-state loading and storage calibration retain their separate contracts.

**Closure condition:** Emit and execute an adequate native implementation and add reference-stable physical initialization, pulse, service, storage and readout errors to the mathematical budget. Preserve complete records and the comparison order: actual code implementation, preparation and evolution errors for that action, then its complete readout. The executed finite two-site programs and the interacting infinite-history certificate provide their respective occupation evidence; neither emits the original large common-code recipe or an infinite apparatus. A laboratory claim requires independently established physical-service and clock attachments on the same process.

Explicit assumption-name navigation: ` declared_common_source_coherent_native_attachment `, ` declared_common_source_active_encoded_native_memory `.

Existing claim-row navigation (no new edges): ` OPH-SCALAR-SEQUENTIAL-INSTRUMENT `, ` OPH-MODEL-SCALAR-RESPONSE-EXECUTION `, ` OPH-MODEL-COMMON-SOURCE-EXECUTION `.

Sources:

- [One composed preparation–history error and resource budget](../model/CONSTRUCTION.md).
- [Exact finite pulse](../paper/tex_fragments/SOURCE_SCALAR_FINITE_INSTRUMENT.tex); exact source anchor: ` prop:source-scalar-finite-pulse `.
- [Sequential composition](../paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex); exact source anchor: ` thm:source-scalar-sequential-instrument `.
- [Fixed-action local splitting](../paper/tex_fragments/SOURCE_SCALAR_TIME_REFINEMENT.tex); exact source anchor: ` prop:source-scalar-time-refinement `.
- [Separate implemented two-event observer, not scalar-word execution](../model/observer/README.md).
- [Executed scalar detector and native feedback contract](../model/observer/scalar_response/README.md).
- [Binary instrument, history law and inherited refinement bound](../model/observer/scalar_response/DERIVATION.md).
- [Retained scalar-response execution and controls](../model/observer/scalar_response/results.json).
- [Local scalar updates and terminal detector probability](../model/observer/scalar_response/local_evolution.py).
- [Local q=13,34 execution and temporal refinement control](../model/observer/scalar_response/local_results.json).
- [Independent inherited temporal detector certificate](../code/source_scalar_time_refinement/verify_source_scalar_time_refinement.py).
- [Emitted source routines and their composed field/pointer execution](../model/observer/source_scalar/experiment.py); exact source anchor: ` def source_history( `.
- [Exact rational bound for splitting the same finite action](../model/observer/source_scalar/bounds.py); exact source anchor: ` def split_error_bound( `.
- [Source-word composition, physical cutoff and remaining error contract](../model/observer/source_scalar/DERIVATION.md).
- [Retained two-site source execution, states, counts and separate error certificates](../model/observer/source_scalar/results.json).
- [Complete ideal occupation-history and output-state error budget](../model/observer/source_scalar/cutoff.py).
- [One coupled quantum history with two source-executed pointer instruments](../model/common_source/experiment.py).
- [Infinite interacting feedback signal and complete error proof](../model/common_source/INFINITE_HISTORY.md).
- [Exact full-H expanding-support history and retained polynomial replay](../model/common_source/infinite_history.py).
- [Retained infinite response, cutoff comparisons and polynomial hashes](../model/common_source/infinite_results.json).
- [Independent Hamiltonian, normalization, support-loss and history checks](../model/common_source/test_infinite_history.py).
- [Uniform interacting pre-readout occupation-refinement proof](../model/common_source/OCCUPATION.md).
- [Rational arbitrary matter/reference cutoff bounds](../model/common_source/occupation.py).
- [Independent interacting moment and boundary-generator checks](../model/common_source/test_occupation.py).
- [Infinite half-line instrument and finite-sign correction](../model/common_source/READOUT.md).
- [Exact Hermite half-line weights and directed pi enclosures](../model/common_source/halfline.py).
- [Independent half-line and complete-readout checks](../model/common_source/test_halfline.py).
- [Retained field, record capacity and complete-history proof](../model/observer/scalar_response/CONTINUATION.md).
- [Executed two-read history and outward rational refinement bounds](../model/observer/scalar_response/continuation_results.json).
- [Canonical finite observer and smooth continuation theorem](../model/common_source/COMPOSITION.md).
- [Complete four-read source/native history and budgets](../model/common_source/composition_results.json).
- [Complete coherent native/work composition and retained physical boundaries](../model/common_source/QUANTUM_NATIVE.md).
- [Source-compiled native attachment and complete-output replay](../model/common_source/quantum_native_results.json).
- [Executed native port routing on the complete data code](../model/common_source/native_routing.py).
- [Complete source-native and retained-work instrument](../model/common_source/quantum_native.py).
- [Active encoded recovery, continuing quantum memory and passive loss control](../model/common_source/ENCODED_MEMORY.md).
- [Source encoder, complete conditional recovery and inverse decoding](../model/common_source/encoded_memory.py).
- [Complete active and passive encoded observer producer](../model/common_source/encoded_observer.py).
- [Encoded observer complete-output comparisons and source costs](../model/common_source/encoded_observer_results.json).

<a id="re-11"></a>

### RE-11: Controlled refinement of the selected quantum family

**Constructed:** The supplied scalar families have separate spatial, occupation and temporal refinement contracts. The two-site source apparatus has an explicit occupation certificate; the shared fixed-band scalar construction increases quadrature resolution without adding modes. The repulsive Bose family controls growing volume at fixed spacing, while the canonical scalar regulator has a free Klein-Gordon refinement and finite-volume interacting bounds. These results do not transfer an interacting relativistic limit between different actions.

On the common interacting scalar/gauge/matter action, moment and projected-Duhamel estimates control pre-readout evolution. The sharp half-line detector has a complete strong-limit construction and a state-specific certified infinite-output comparison; its long occupation tail prevents small probability error from certifying small complete-output error. Finite sign-projector and rectangular-code support certificates retain their preparation domains. Sharp measurement may destroy finite kinetic-energy moments, so a later evolution needs its own estimate. The distinct smooth Ramsey instrument has quantitative occupation/temporal bounds and finite-energy control through repeated reads, feedback and a later gauge-clock probe. Its retained sixteen-output certificate is for the fixed blank-gauge/coupled-vacuum preparation, with independent native/reference inputs; it is not an arbitrary-gauge or correlated-work certificate.

For native C/E federations, the complete local-output bound controls retained quantum states, local records, local feedback and untouched references. Compatible restrictions on a fixed bounded-degree graph yield a quasilocal quantum process at bounded auxiliary per-seam Poisson parameter. Actual screen parent/child mean restrictions require memory and do not automatically intertwine dynamics. Each locality or refinement argument must use compatible initial states and access maps on the same source family. Actual source-screen pairs 1/3, 2/4 and 3/5 have native-port-labelled neighborhoods agreeing through radii 3,6 and 9. Matched radius-nine initial restrictions give complete root-output distance at most 2187/5734400 at per-seam s=1/4 between 15360 and 245760 ports. The explicit product-state family is compatible; fixed global total-two preparations have a separately charged local mismatch. Exact mean-density return kernels agree through uniformization power 22 and first differ at 23; this stronger first-moment check is not a complete-history bound. The actual-screen field protocol also admits uniform fixed-number preparations: connected native exchanges derive the unique diagonal stationary law within each conserved sector. Exact hypergeometric restrictions, rather than independent product occupations, certify all 44 field/timeout outcomes at totals 5120/15360 and 81920/245760. Positive baseline and negative feedback fourth cumulants survive. Each screen uses its own initial restriction; the complete field-law distance is nonzero and certified near 3.1e-5. Eight paired full-source histories preserve exact totals, native feedback and serialized field-counter continuation. The conserved total and exact sampling service remain inputs; stationarity supplies no finite physical preparation certificate.

**Remaining input:** A common controlled interacting relativistic observation limit remains open. Growing volume, increasing occupation and reducing a temporal splitting step are different operations. All-level native source-screen compatibility, nontrivial operational scaling and laboratory identification must be established on the actual family. Supplied scalar actions retain their own geometry and preparation inputs; their bounds do not select the native law. Numerical source-array certificates and mathematical ideal limits do not authenticate physical service errors. Postselected output comparisons need their own success-probability control, and infinite stopping needs a tail bound.

**Closure condition:** Use each proved occupation and readout certificate with its actual initial-state, record and observable contract. Keep the sharp-sign and smooth Ramsey protocols distinct, and preserve the separate finite modular-clock scope. Control source compilation and physical-service errors before an apparatus claim. Other interactions or spatially refining families require their own uniform errors and compatible states/observables; the existing free-scalar spatial limit remains a different theorem.

Explicit assumption-name navigation: ` declared_common_source_smooth_ramsey_instrument `.

Existing claim-row navigation (no new edges): ` OPH-SCALAR-SPACE-TIME-DETECTOR-LIMIT `, ` OPH-MODEL-SCALAR-RESPONSE-EXECUTION `, ` OPH-MODEL-COMMON-SOURCE-EXECUTION `, ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Sources:

- [Existing scalar continuum](../paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex); exact source anchor: ` prop:common-scalar-continuum-limit `.
- [Joint observable refinement](../paper/tex_fragments/SOURCE_SCALAR_TIME_REFINEMENT.tex).
- [Uniform refinement hypotheses](../claims/assumption_dictionary.md).
- [Finite observer and scalar continuum scope distinction](../model/CONSTRUCTION.md).
- [Executed scalar detector and native feedback contract](../model/observer/scalar_response/README.md).
- [Binary instrument, history law and inherited refinement bound](../model/observer/scalar_response/DERIVATION.md).
- [Retained scalar-response execution and controls](../model/observer/scalar_response/results.json).
- [Local scalar updates and terminal detector probability](../model/observer/scalar_response/local_evolution.py).
- [Local q=13,34 execution and temporal refinement control](../model/observer/scalar_response/local_results.json).
- [Independent inherited temporal detector certificate](../code/source_scalar_time_refinement/verify_source_scalar_time_refinement.py).
- [Reference-stable finite temporal splitting certificate](../model/observer/source_scalar/bounds.py).
- [Complete ideal occupation bound and separation from spatial continuum evidence](../model/observer/source_scalar/DERIVATION.md).
- [Rational Gaussian-moment and projected Duhamel occupation-refinement certificate](../model/observer/source_scalar/cutoff.py).
- [Independent Bargmann, boundary-generator and complete-output checks](../model/observer/source_scalar/test_cutoff.py).
- [Independent finite interacting signal enclosure and temporal error composition](../model/common_source/finite_signal.py).
- [Regional CP maps and the exact occupation-compression obstruction](../model/common_source/spatial.py).
- [Interacting occupation proof and separate sharp-readout convergence scope](../model/common_source/OCCUPATION_REFINEMENT.md).
- [Exact rational pre-readout moment and Duhamel bound](../model/common_source/refinement.py).
- [Complete smooth-instrument occupation bound and infinite ensemble energy identity](../model/common_source/SMOOTH_REFINEMENT.md).
- [Directed full-instrument error budget with explicit small-register vacuity](../model/common_source/smooth_refinement.py).
- [Infinite interacting feedback signal and complete error proof](../model/common_source/INFINITE_HISTORY.md).
- [Exact full-H expanding-support history and retained polynomial replay](../model/common_source/infinite_history.py).
- [Retained infinite response, cutoff comparisons and polynomial hashes](../model/common_source/infinite_results.json).
- [Independent Hamiltonian, normalization, support-loss and history checks](../model/common_source/test_infinite_history.py).
- [Uniform interacting pre-readout occupation-refinement proof](../model/common_source/OCCUPATION.md).
- [Rational arbitrary matter/reference cutoff bounds](../model/common_source/occupation.py).
- [Independent interacting moment and boundary-generator checks](../model/common_source/test_occupation.py).
- [Infinite half-line instrument and finite-sign correction](../model/common_source/READOUT.md).
- [Exact Hermite half-line weights and directed pi enclosures](../model/common_source/halfline.py).
- [Independent half-line and complete-readout checks](../model/common_source/test_halfline.py).
- [Retained field, record capacity and complete-history proof](../model/observer/scalar_response/CONTINUATION.md).
- [Executed two-read history and outward rational refinement bounds](../model/observer/scalar_response/continuation_results.json).
- [Complete sharp-output certificate](../model/common_source/sharp_output.py).
- [Complete sharp-output certificate](../model/common_source/sharp_pointer.py).
- [Complete sharp-output certificate](../model/common_source/sharp_output_results.json).
- [Independent sharp-output and adversarial controls](../model/common_source/test_sharp_output.py).
- [Independent sharp-output and adversarial controls](../model/common_source/test_sharp_pointer.py).
- [Independent sharp-output and adversarial controls](../model/common_source/test_sharp_rate.py).
- [Canonical finite observer and smooth continuation theorem](../model/common_source/COMPOSITION.md).
- [Complete four-read source/native history and budgets](../model/common_source/composition_results.json).
- [Prepared complete smooth-history certificate and native/encoded lifting scope](../model/common_source/SMOOTH_HISTORY.md).
- [Exact sixteen-output continuation enclosure and canonical replay](../model/common_source/smooth_history.py).
- [All-current prepared finite/infinite comparisons with separate temporal errors](../model/common_source/smooth_history_results.json).
- [Complete quantum and record locality and conditional local limit](../model/native_action/README.md).
- [Actual source-screen coarse-density memory](../model/federated_observer/MEMORY.md).
- [Fixed-band shared observation family](../model/scalar_observations/README.md).
- [Supplied interacting fixed-spacing volume limit](../model/spatial_interaction/README.md).
- [Canonical regulator and separate free relativistic refinement](../model/relativistic_scalar/README.md).
- [Exact actual-screen local instrument maps and state mismatch](../model/native_action/screen_compatibility.py).
- [Source-screen neighborhood, runtime and density-response certificates](../model/native_action/screen_results.json).
- [Exact canonical-sector field law and source execution](../model/observer_time_fields/canonical_sector.py).
- [Fixed-number full-screen field predictions and finite-size comparison](../model/observer_time_fields/canonical_sector_results.json).
- [Common clock/field proof and conserved-sector extension](../model/observer_time_fields/REFINEMENT.md).

<a id="re-12"></a>

### RE-12: Source clock, physical time and execution cost

**Constructed:** Source count and finite execution-cost results retain their declared clock and service domains. Native local event instruments sum actual unresolved repairs, retain remote evolution and include finite-cap timeouts. The joint clock/field process supplies exact reduced-memory identities and a four-mark completed-event closure on its specified one-particle support; conditional reference states and controller memory remain. Equal completed clock statistics need not imply equal phase observations. Independent identity holds cancel from completed stopped laws without fixing a laboratory duration.

A chosen faithful algebra/reference defines its modular action. The native Palm reference, reconstructed transfer reference and supplied full-work Gibbs reference are separate pairs. Source probes execute the relevant finite modular automorphisms on stated accessible codes. For sigma_beta proportional to exp(-beta H), Hamiltonian segments have t=-beta*s, with reference choice, beta and full-work access supplied. Exact native covariance defects prevent equating particular event operations with those modular flows. Actual interventions and feedback are retained as instruments, rather than treated as one autonomous modular group. Pulse, native-attempt, event and load/export counts are separately accounted resources. Instrument-valued clock conversion equals direct native stopping, with every remote-count coefficient and conditional quantum/reference output retained. The passive two-A-tick renewal cycle has B-count mean 2 and variance 16; specified B-entry feedback changes the law. A uniform two-event survival bound certifies completion, while capped runtime execution retains all timeouts. This derives a stochastic relative-clock law, not a deterministic conversion from equal mean rates. On the 66-state binary total-two sector of one twelve-port carrier, actual local movement clocks have asymptotic count rate 1/6 per carrier movement from any initial density in that sector. An exact Poisson corrector supplies finite count-error bounds without stationary preparation or an independent clock. Their long-run covariance distinguishes all three graph-distance shells. A positive seconds calibration changes no such dimensionless law; empirical attachment requires a matched source apparatus and held-out paired records. A recorded software attachment executes the same 264532 native attempts in each of two runs against a host monotonic counter. One baseline block calibrates port 2; three held-out blocks measured with port 3 have timing errors 2.3% to 7.8%, within the fixed 25% engineering tolerance. Host-only workload preserves all native records and rejects that calibration in every control block. Replay checks the seeded choices and source currents, not timestamp authenticity or clock traceability; physical proper time and modular flow are not inferred. Native circulation has an action in one port's own tick count. A local-tick Palm herald succeeds with probability 5/198; failure remains recorded. The stopped occupation identity gives six carrier movements per local tick, exact cycle covariance P/5 and triangle variance 3m/5 after m local ticks. Finite fourth moments and adjacent squared-current dependence distinguish the law from independent Gaussian blocks. The long-count Gaussian increment action has coefficient 5/2 on the cycle subspace, without replicas or an independent clock. Carrier-wide current access and the uniform scheduler are declared; this is not a relativistic action. Fixed-window quantum execution preserves herald failures and zero-tick outcomes, while record-controlled seam selection changes the finite law. Exact incident-clock field correlations reconstruct a seven-dimensional cyclic detector Hilbert space and its transfer spectrum on the elementary total-two source; independent full 66-state propagation checks withheld lags. The positive-transfer moment class is supplied by the existing theorem. Exact rank is not inferred from sampled records, and the seven-mode Gram is severely ill-conditioned. On each actual level-3/5 screen with product preparation, bounded record correlations also certify a two-vector temporal subspace and compressed transfer values near 0.88702 and 0.97079 with error below 2e-5. These compressed values are not a closed two-mode spectrum or laboratory energy calibration. The exact marked incident/movement-clock law retains zero-current local hits, remote waits, field-conditioned continuations and feedback hits. For the uniform total-two single-port source, I/K tends to 33/5; an exact Poisson corrector supplies a finite conversion-error bound indexed by fixed global movement count, uniformly over initial states. The retained host record comparison preserves baseline acceptance and workload-control rejection under that conversion uncertainty without refitting. It does not transfer the ratio or calibration to the different root-clock product/canonical-sector field protocol, identify a physical clock, or authenticate the pseudorandom entropy.

**Remaining input:** No preferred universal clock is required. Relating an observer's records to modular evolution, laboratory durations or effective proper time requires an operational correspondence preserving joint and conditional observations. A common mean rate or affine label comparison does not supply that correspondence. Physical service accuracy, reference preparation and clock calibration remain explicit. The scalar benchmark time and native event parameter cannot be identified merely because they occur in one program; unbounded event protocols need well-posed stopping and complete tail control. Choosing seconds is a unit convention. The dynamical obligation is a compatible joint record law and an independently validated apparatus correspondence, rather than an absolute duration derived without calibration.

**Closure condition:** Prove the same-history clock/observable intertwiner and its domain, then calibrate independently. A local generator alone does not establish a microscopic causal cone.

Explicit assumption-name navigation: ` declared_common_source_thermal_reference `.

Existing claim-row navigation (no new edges): ` OPH-MODEL-N4-LOCAL-TICK-CERTIFICATE `, ` OPH-MODEL-SCALAR-RESPONSE-EXECUTION `, ` OPH-MODEL-COMMON-SOURCE-EXECUTION `, ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Sources:

- [Source count-clock theorem](../paper/tex_fragments/SOURCE_COUNT_CLOCK.tex); exact source anchor: ` thm:source-count-clock `.
- [Scalar action clock](../paper/tex_fragments/SOURCE_SCALAR_CLOCK.tex).
- [Conditional clock conversion](../Lean/Screen/CurlStoneClockBridge.lean).
- [N4 shared bare-occupation response and finite memory certificate](../code/native_response/README.md).
- [Incident-attempt stopping channel and certified N4 bounds](../model/local_tick/README.md).
- [Exact local stopping null and certificate replay](../model/local_tick/experiment.py); exact source anchor: ` def calculate() `.
- [Independent finite-chain and quantum stopping checks](../model/local_tick/test_local_tick.py).
- [Occupied-port endpoint instrument and corrected record statistic](../code/native_response/port_readout.py).
- [Fixed model-duration recipes and actual pulse/service counts](../model/observer/source_scalar/experiment.py); exact source anchor: ` def program( `.
- [Counts for native pulses and proper-code load/export requests](../model/observer/source_scalar/compiler.py); exact source anchor: ` def resources( `.
- [Retained two-site source execution, states, counts and separate error certificates](../model/observer/source_scalar/results.json).
- [Observer Hilbert space and the distinction between modular and physical time](../model/common_source/OBSERVER_READOUTS.md).
- [Source-executed modular clock with retained outcomes and resources](../model/common_source/modular_clock.py).
- [Fixed-reference thermal alignment and physical clock boundary](../model/common_source/THERMAL_CLOCK.md).
- [Record-preserving source probes and independent thermal targets](../model/common_source/thermal_clock.py).
- [Retained field, record capacity and complete-history proof](../model/observer/scalar_response/CONTINUATION.md).
- [Executed two-read history and outward rational refinement bounds](../model/observer/scalar_response/continuation_results.json).
- [Canonical finite observer and smooth continuation theorem](../model/common_source/COMPOSITION.md).
- [Complete four-read source/native history and budgets](../model/common_source/composition_results.json).
- [Complete field memory and two source event clocks](../model/joint_clock_process/README.md).
- [Stopped source instrument and corner modular obstruction](../model/relational_event_clock/README.md).
- [Invariant protocol algebra and native covariance defect](../model/relative_modular_dynamics/README.md).
- [Source-executed transfer-reference modular observations](../model/transfer_modular_probe/README.md).
- [Finite recurrence-clock contracts and scope](../Lean/QFT/SourceRecurrenceClock.lean).
- [Complete count-generating clock conversion and native stopping oracle](../model/joint_clock_process/clock_change.py).
- [Single-carrier clock comparison, spatial covariance and laboratory attachment contract](../model/native_clock/COMPARISON.md).
- [Measured host-clock correspondence and timing-control scope](../model/native_clock/HOST_CLOCK.md).
- [Raw native record replay and held-out timing comparison](../model/native_clock/host_clock.py).
- [Recorded native choices, currents and host-counter block endpoints](../model/native_clock/host_clock_records.json).
- [Calibration, held-out blocks and workload-control verdicts](../model/native_clock/host_clock_results.json).
- [Exact stopped current moments and action in one local clock](../model/native_action/current_clock_action.py).
- [Local-Palm action, finite higher moments and complete observation-window receipt](../model/native_action/current_clock_action_results.json).
- [Detector temporal representation recovered from field correlations](../model/relative_observer_geometry/temporal_records.py).
- [Exact visible modes and certified screen transfer compression](../model/relative_observer_geometry/temporal_results.json).
- [Joint incident and movement clocks with finite conversion uncertainty](../model/native_clock/CLOCK_BRIDGE.md).
- [Retained joint clock and host comparison](../model/native_clock/clock_bridge_results.json).

<a id="re-13"></a>

### RE-13: Physical geometry and responsive metric

**Constructed:** Finite source geometry, golden coordinates and conditional Lorentz/spacetime constructions exist. The canonical source includes its adopted oriented boundary and complete-response clauses. The common-source attachment reuses the existing native twelve-port carrier and proves I-L_ico/60 by exact enumeration of every total-two state, seam and placement coin under a uniform all-seam calibration. The lowest positive Laplacian band has rank three and Gram G=4P3, with an exact icosahedral coordinate witness. Every actual joint history receives a dimensionless displacement readout from its committed currents and declared seams, checked against working-load changes without giving the controller extra state access. The original E/F controller has a rank-two displacement span. A declared binary F/G continuation after the first E event reexecutes every native branch while retaining the same quantum history, and its three admitted seam displacements span the rank-three carrier. The existing native integer rule also gives pathwise load continuity and fair-coin conditional current mean (x_a-x_b)/2 and variance parity(x_a-x_b)/4. Exact enumeration of all 4680 state/seam/coin events builds the full 78-state transition kernel for the existing uniform geometry calibration and proves its closed mean law I-L_ico/60. The binary 66-state uniform distribution is stationary with continuing load fluctuations. All nonconstant mean multipliers are positive contractions, so this fixed linear mean readout cannot be an undamped wave under a positive scalar clock rescaling. The actual record-controlled E/F protocol is executed as a distinct control; no new clock, Hamiltonian, physical field identification or continuum premise is supplied. Exposure-normalized even-load current laws recover the same unweighted rank-three structural Gram without a uniform schedule or fair coins. The actual slowest response under a doubled seam weight has multiplicity two. The full inverse-response Gram instead has rank eleven, so choosing the slow-band metric is an operational readout choice, not a proof of physical bulk dimension. These exact single-carrier constructions require no federation. On the same uniform single carrier, asymptotic local-tick covariance independently recovers the native distance shells from a continuing movement-count history. Deterministic even-load currents also permit exact reconstruction from 360 executed native shots, one per preparation/seam, without a probability-law oracle. This removes native sampling uncertainty for that protocol, not laboratory preparation or detector errors. Native circulation has an action in one port's own tick count. A local-tick Palm herald succeeds with probability 5/198; failure remains recorded. The stopped occupation identity gives six carrier movements per local tick, exact cycle covariance P/5 and triangle variance 3m/5 after m local ticks. Finite fourth moments and adjacent squared-current dependence distinguish the law from independent Gaussian blocks. The long-count Gaussian increment action has coefficient 5/2 on the cycle subspace, without replicas or an independent clock. Carrier-wide current access and the uniform scheduler are declared; this is not a relativistic action. Fixed-window quantum execution preserves herald failures and zero-tick outcomes, while record-controlled seam selection changes the finite law.

**Remaining input:** Choose the physical region/metric/area/volume attachment and a compatible source-causal continuum. A frozen chart cannot itself provide changing physical curvature. Uniform calibration and the two-event binary controllers are separate protocols. Finite calibration kernels have eleven nonconstant modes; rank three refers to their slow response quotient. The added controller spans that carrier but establishes no isotropic or autonomous motion. Neither this readout nor the two-scalar interaction graph identifies physical position, momentum, length scale or a spatial continuum. The uniform-scheduler premise belongs to the actual slow dynamics; dividing exposure out of exact laws does not remove that premise for motion. A compatible extended family is required for a continuum claim, but the specific spherical federation is optional.

**Closure condition:** Construct operational metric and geometric response on the same source and refinement family. Preserve the distinction between action capacity, event-count measure and physical volume.

Explicit assumption-name navigation: ` declared_native_slow_geometry_readout `.

Existing claim-row navigation (no new edges): ` OPH-MODEL-COMMON-SOURCE-EXECUTION `, ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Sources:

- [Native geometric-source boundary](../paper/tex_fragments/NATIVE_GEOMETRIC_SOURCE_IDENTIFIABILITY.tex).
- [Adopted source clauses](../docs/AXIOM_REFERENCE.md).
- [Continuum and metric interfaces](../docs/registers/PREMISE_REGISTER_V3.md#gravity-rungs).
- [Native response calibration and original and extended binary-controller scopes](../model/common_source/SLOW_GEOMETRY.md).
- [Exact slow projector and committed-current displacements on all joint paths](../model/common_source/slow_geometry.py).
- [Exact native continuity, mean transport and their restricted field interpretation](../model/common_source/NATIVE_FIELDS.md).
- [Exact structural recovery and distinct response-metric domains](../model/relative_observer_geometry/intrinsic_geometry.py).
- [Single-carrier geometry and scheduler certificates](../model/relative_observer_geometry/intrinsic_results.json).
- [Single-carrier clock comparison, spatial covariance and laboratory attachment contract](../model/native_clock/COMPARISON.md).
- [Exact stopped current moments and action in one local clock](../model/native_action/current_clock_action.py).
- [Local-Palm action, finite higher moments and complete observation-window receipt](../model/native_action/current_clock_action_results.json).

<a id="re-14"></a>

### RE-14: Gauge algebra to physical gauge fields

**Constructed:** The finite local Standard Model Lie-algebra result remains proved under the adopted A1 complete-response/closure clauses and A2 endogenous transport. Compact-Lie classification is an explicit mathematics import. The designed N1 charge assignment Q_v=g(P)N_v gives exact branch charge transfer and local U(1) covariance; its binary loop detector makes g(P) affect the native feedback response. The common-source gauge module embeds two oppositely charged CAR channels and an integer electric link into one invariant four-dimensional code. Dressed hopping, electric energy and neutral scalar coupling preserve operator Gauss and even parity on that code, with native source execution. This finite U(1) example supplies neither non-Abelian dynamics nor a chiral gauge continuum or physical electromagnetic normalization.

**Remaining input:** The physical current/action, same-source category, global form and coupling attachment remain separate. The scalar field construction does not instantiate a physical gauge theory. Compatibility with the twelve-port reversible response is an immediate law-selection constraint, even though constructing physical gauge fields is downstream. Scalar seam equalizers cannot serve as that current: their tested Lie closure is so(11), not the required twelve-dimensional algebra. Conserved difference drives in the endpoint-current lift omit the central relative-phase generator; the full twelve-direction response is a supplied interface, not derived physical gauge dynamics. A carrier-covariant lower-block deletion retains N2 locality, completion and protected populations while failing injective twelve-direction response; those successful properties cannot replace the gauge theorem premise. The original direct-sum source algebra and the candidate tensor factors require an explicit interface. This chosen U(1) response does not derive the Standard Model representation, Maxwell field-source dynamics or the omitted nonabelian interactions. For the original marked A5 action, the slow geometry triplet and lower response triplet have unequal order-five traces (1-sqrt(5))/2 and (1+sqrt(5))/2; no invertible same-action intertwiner exists (#1088). The addressed construction keeps these spatial and internal domains distinct. A whole-law mirror or action twist needs an explicit attachment.

**Closure condition:** Bind the proved algebra and representation to actual source observables, dynamics and refinement; supply the categorical and physical gauge interfaces used by the selected branch.

Existing claim-row navigation (no new edges): ` OPH-MODEL-COMMON-SOURCE-EXECUTION `.

Sources:

- [Finite gauge structure and physical attachments](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex); exact source anchor: ` thm:current-algebra `.
- [Gauge source and global-form entries](../docs/registers/PREMISE_REGISTER_V3.md).
- [Chosen charge transfer and branch gauge covariance](../model/observer/p_response/DERIVATION.md).
- [Exact fixed-marked-action triplet obstruction](../code/programme_obstruction_audit/joint_independent.py).
- [Collective response, actual chart implementers and marked-triplet tests](../code/native_response/test_carrier_response.py).
- [CAR signs, quantum electric flux and invariant Gauss-sector source operations](../model/common_source/gauge_matter.py).

<a id="re-15"></a>

### RE-15: Physical edge factors and collar recovery

**Constructed:** Finite charge centers, regional completion and nontrivial bulk collar splits are constructed. Existing exact Markov and conditional recovery theorems remain available.

**Remaining input:** The physical edge-center split and state alignment, edge-reference factorization and log-d normalization must agree. Approximate recovery needs strong conditional mixing and the sharp collar rate on one family. These are alternative recovery routes. Strong matrix-mixing gives direct log-defect control without a universal spectral floor; the physical EC-replacement modulus remains separate.

**Closure condition:** Realize the selected physical split and prove either its exact central-interface condition or its complete uniform approximate-recovery hypotheses. Charge retention alone is insufficient.

Explicit assumption-name navigation: ` central_interface_clause `, ` EC_aligned_Markov_replacement_control `, ` strong_conditional_gibbs_mixing `, ` sharp_collar_rate_receipt `, ` edge_normalization_log_d `.

Existing claim-row navigation (no new edges): ` OPH-GR-D4D-ENTROPY-BRIDGE `.

Sources:

- [Gravity rungs3–5](../docs/registers/PREMISE_REGISTER_V3.md#gravity-rungs).
- [EC alignment, mixing and edge normalization](../claims/assumption_dictionary.md).

<a id="re-16"></a>

### RE-16: Modular flow, null translations and stress

**Constructed:** Typed modular/half-sided inclusion and conditional local-stress implications are available; finite cap diagnostics are constructed. The common finite history supplies an explicit faithful scalar marginal and its modular flow, checked against GNS pairing and finite KMS identities. The classical record algebra has identity flow; sharp-record conditional states are nonfaithful on the full scalar algebra. The fixed full-work Gibbs reference provides an exact finite modular/Hamiltonian identity with t=-beta*s and a distinct evolving actual state. This uses the specified finite Hamiltonian and is not a support-visible BW or null-translation attachment.

**Remaining input:** Support-visible BW identification, common null translation, effective locality, tomography and the same-source physical stress attachment remain required.

**Closure condition:** Instantiate the modular and stress interfaces on the candidate with the required common domains, locality and normalization; algebraic modular flow alone is not physical time or stress.

Explicit assumption-name navigation: ` null_stress_bridge `, ` declared_common_source_thermal_reference `.

Existing claim-row navigation (no new edges): ` OPH-MODEL-COMMON-SOURCE-EXECUTION `.

Sources:

- [Gravity rungs7–8](../docs/registers/PREMISE_REGISTER_V3.md#gravity-rungs).
- [Modular and stress construction](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex).
- [Finite modular readout on the actual interacting observation state](../model/common_source/operational.py).
- [Finite thermal modular identity with BW and stress obligations retained](../model/common_source/THERMAL_CLOCK.md).

<a id="re-17"></a>

### RE-17: Joint state–geometry entropy stationarity

**Constructed:** Finite counting stationarity and fixed-cap MaxEnt/envelope identities are proved under their declared variations.

**Remaining input:** The candidate must supply constraint-tangent generalized-entropy stationarity for joint state and geometry variations, with the physical area/edge identification.

**Closure condition:** Construct the permitted variation family and show the full generalized entropy has zero first variation at fixed physical volume. A fixed-reference entropy optimum alone does not do this.

Explicit assumption-name navigation: ` fixed_cap_MaxEnt_variation_family `.

Sources:

- [Gravity rung6](../docs/registers/PREMISE_REGISTER_V3.md#gravity-rungs).
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

- [Gravity rungs8–11](../docs/registers/PREMISE_REGISTER_V3.md#gravity-rungs).
- [Flagship gravitational interface](../flagship/from_observer_consensus_to_standard_physics.tex).
- [Common tower and universal coupling](../claims/assumption_dictionary.md).

<a id="re-19"></a>

### RE-19: Matter, statistics and interactions

**Constructed:** Finite carrier/representation and declared Higgs/Yukawa packets exist; the active scalar candidate has a conditional free Bose limit. The common finite register has an actual neutral-scalar/density interaction g sqrt(2) Q0(nL0-nR0), with g=1/4 supplied. The source-executed word couples the field to charged hopping rather than leaving tensor factors independent. Its finite Gauss sector and detector response do not select matter content, empirical parameters, scattering amplitudes or an interacting 3+1-dimensional limit.

**Remaining input:** Select physical sectors, spin/chirality/statistics, kinetic attachments, interactions, coupling coefficients and stable poles on the same process. A declared finite grammar is not its physical selection.

**Closure condition:** Construct the selected interacting source family and its controlled observables/continuum, reusing the algebraic packets without treating supplied matter content as derived.

Existing claim-row navigation (no new edges): ` OPH-MODEL-COMMON-SOURCE-EXECUTION `.

Sources:

- [Matter and QFT interfaces](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex).
- [Matter selection and action hypotheses](../claims/assumption_dictionary.md).
- [Declared finite scalar/matter interaction on the shared source](../model/common_source/scalar_coupling.py).

<a id="re-20"></a>

### RE-20: Parameter selection and P

**Constructed:** Exact conditional identities, normalized-capacity reconstruction and selection boundaries are recorded in the existing ledgers. The finite detector accepts P and N as closure-derived premises, executes the existing P-to-alpha map, and uses its chosen charge assignment and calibrated loop exposure to produce a retained detector bit and P-dependent native feedback. An exact retrospective endpoint audit holds P and the declared branch fixed and reports an inverse-alpha deficit of about 300.388218 ppm against CODATA 2022. It separates the pending physical transport contract from the inconsistent additional hypothesis of a complete endpoint with zero correction.

**Remaining input:** The detector is a conditional design under accepted closure inputs; its charge assignment, probe, flux, readout and protocol are specified apparatus. Executing the existing source-only coupling branch does not upgrade its physical endpoint status. Koide balance, the complete screen/electroweak attachment, physical realization and independent comparison retain their own scopes. The five-standard-uncertainty comparison allowance is not a rigorous confidence interval; no theory-error significance or derived missing transport term is supplied.

**Closure condition:** For independent physical parameter claims, give the claimed source selection and same-quantity observable bridge. The finite detector supplies a conditional P-to-response map, not a P-only temperature or event-rate law; design uniqueness is not required for its implementation.

Existing claim-row navigation (no new edges): ` OPH-MODEL-P-RESPONSE-EXECUTION `.

Sources:

- [Selection ledger](../docs/registers/SELECTION_LEDGER.md).
- [Parameter and comparison premises](../docs/registers/PREMISE_REGISTER_V3.md).
- [Frozen comparison boundaries](../docs/FROZEN_PREDICTION_LADDER.md).
- [P-dependent electromagnetic observation and native feedback](../model/observer/p_response/README.md).
- [Executed declared P-to-coupling map](../model/observer/p_response/experiment.py); exact source anchor: ` def coupling( `.
- [Pinned retrospective endpoint residual and separate incomplete-map verdict](../model/predictions/endpoint_audit.py).

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

**Constructed:** Native C/E instruments give different complete coherent-readout laws under shared classical endpoint dynamics. Source compilation executes both alternatives on the full 78-state encoding, with all outcomes, padding and arbitrary references. The continuing-carrier discriminator uses finite native mixing before each trial and a uniform conditional-law coupling bound, so actual trial independence and exact stationary resets are unnecessary. Its locked design uses 379 discarded attempts before each of 3927 trials and bounds wrong-candidate probability by 3/400 in the exact-source model. The retained execution verifies two consecutive trials, not a sampled full campaign. A separately authenticated physical per-trial conditional error at most 1/1570800 would suffice for a one-percent bound.

Public oscillator comparisons execute a supplied experimental action and apparatus model, with retrospective calibration/evaluation splits. Zero-displacement, opposite-probe and two-order Gram identities constrain observations independently of the oscillator action; their aggregate-mean residuals are deterministic diagnostics. The calibrated measurement evaluator retains null settings, complete counts, apparatus budgets and simultaneous sampling error. These controls provide explicit measurement requirements without treating software trials as physical data. A randomized Gram witness has a fixed-horizon conditional concentration bound despite preparation drift and quantum memory. Fresh settings must remain independent of the joint prepared work/reference state conditional on the past. Its two-trial source program retains 100 outcomes; physical detector-error bounds and fresh complete records are needed. Public means determine future threshold-crossing budgets conditional on repeated negative scores, not power or a present-data rejection. A retrospective replay of the deposited Hensen Bell records reproduces 196 wins in 245 eligible trials, S=2.42249987046 and a fixed-horizon local-response null tail 0.0391746435588 using the published setting allowance. Exact native local readbacks instantiate the classical branch and all sixteen deterministic local response tables establish its CHSH bound under setting independence and local factorization. Eligibility, locality, setting predictability and fixed-horizon validity remain external premises. This constrains a local-factorized extension of the diagonal field branch; it excludes neither full quantum OPH nor either C/E candidate and supplies no native lab calibration.

**Remaining input:** The native discriminator needs a justified laboratory map and physical service error bounds; no deposited dataset is identified with that native protocol. The public oscillator means do not authenticate accepted shot counts, common preparation or all coordinate conventions. A statistical rejection needs a valid sampling/selection model and calibrated instruments. Mathematical or source-array error is not hardware fidelity. Retrospective witnesses and comparisons cannot be relabeled as preregistered observations, and a general quantum consistency test is not an OPH-specific prediction.

**Closure condition:** Freeze an unavoidable relation and rejection rule with full remaining freedom, baseline and uncertainties before an independent eligible observation; preserve the frozen ladder unchanged.

Explicit assumption-name navigation: ` declared_randomized_Gram_witness_protocol `.

Existing claim-row navigation (no new edges): ` OPH-MODEL-NATIVE-OBSERVER-EXECUTION `, ` OPH-MODEL-P-RESPONSE-EXECUTION `, ` OPH-MODEL-N4-LOCAL-TICK-CERTIFICATE `, ` OPH-MODEL-SCALAR-RESPONSE-EXECUTION `, ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Sources:

- [Existing frozen targets and eligibility](../docs/FROZEN_PREDICTION_LADDER.md).
- [Propagation attachment premises](../docs/registers/PREMISE_REGISTER_V3.md).
- [Native seeded execution and paired controls](../model/observer/results/run.json).
- [Detector exact histories and seeded execution](../model/observer/p_response/results.json).
- [Local-tick null and inherited strict lower bounds](../model/local_tick/results.json).
- [Executed scalar detector and native feedback contract](../model/observer/scalar_response/README.md).
- [Binary instrument, history law and inherited refinement bound](../model/observer/scalar_response/DERIVATION.md).
- [Retained scalar-response execution and controls](../model/observer/scalar_response/results.json).
- [Source-generated detector law, complete continuation and sampling controls](../model/observer/source_scalar/experiment.py).
- [Independent finite-target and complete outcome checks](../model/observer/source_scalar/test_experiment.py).
- [Retained two-site source execution, states, counts and separate error certificates](../model/observer/source_scalar/results.json).
- [Locked candidate predictions, complete counts and simultaneous decision rule](../model/predictions/protocol.py).
- [Mutation, complete-count, statistical-budget and synthetic-discrimination controls](../model/predictions/test_protocol.py).
- [Retained field, record capacity and complete-history proof](../model/observer/scalar_response/CONTINUATION.md).
- [Executed two-read history and outward rational refinement bounds](../model/observer/scalar_response/continuation_results.json).
- [Executed source-word C/E discriminator against the unchanged four-setting predictions](../model/observer/source_discrimination/experiment.py).
- [Sufficient complete-instrument error budget for the unchanged 1/64 apparatus allowance](../model/observer/source_discrimination/certificate.py); exact source anchor: ` def apparatus_error_contract( `.
- [Retained encoded C/E source execution and explicit conditional error allocation](../model/observer/source_discrimination/results.json).
- [Continuing native-carrier discrimination without an IID premise](../model/closed_prediction_protocol/README.md).
- [Two-order positivity constraints and empirical scope](../model/characteristic_constraints/README.md).
- [Complete-count calibration and sampling contract](../model/calibrated_measurement/README.md).
- [Action-independent empirical readout controls](../model/readout_controls/README.md).
- [Certified randomized witness and conditional fixed-horizon test](../model/characteristic_constraints/randomized.py).
- [Continuing source execution and conditional-score comparison](../model/characteristic_constraints/randomized_source.py).
- [Classical branch bound and retrospective public Bell comparison](../model/predictions/BELL_RECORDS.md).
- [Complete published-data selection and exact null-tail replay](../model/predictions/bell_results.json).

<a id="re-23"></a>

### RE-23: External calibrations and empirical imports

**Constructed:** The canonical register explicitly distinguishes hadronic empirical payloads and laboratory clock/energy anchors from source deductions. The endpoint audit imports only the cited CODATA 2022 inverse-alpha comparison value and standard uncertainty. It is a known-data retrospective check, with P and the inherited incomplete source branch unchanged, no fitted correction and no empirical validation verdict. A separate C/E apparatus comparison declares independently established coherent preparation/readout errors and fresh-trial independence; its software evaluator does not authenticate those premises. The native clock comparison supplies dimensionless count-ratio and covariance predictions plus an explicit paired-record calibration contract. No existing laboratory apparatus or dataset is identified with that native source. Actual paired source-count/host-timestamp records are retained for an engineered software implementation. Their held-out timing check and rate-change control test that execution correspondence; they supply no independently matched natural apparatus, traceable clock-accuracy bound, or empirical OPH confirmation. A retrospective replay of the deposited Hensen Bell records reproduces 196 wins in 245 eligible trials, S=2.42249987046 and a fixed-horizon local-response null tail 0.0391746435588 using the published setting allowance. Exact native local readbacks instantiate the classical branch and all sixteen deterministic local response tables establish its CHSH bound under setting independence and local factorization. Eligibility, locality, setting predictability and fixed-horizon validity remain external premises. This constrains a local-factorized extension of the diagonal field branch; it excludes neither full quantum OPH nor either C/E candidate and supplies no native lab calibration.

**Remaining input:** Every result using these inputs retains their provenance, uncertainty and transfer contract. Fundamental-input counts do not absorb or erase them. The endpoint residual does not provide the missing physical transport or its theory-error budget; five reported standard uncertainties are a comparison allowance, not a rigorous confidence interval.

**Closure condition:** Supply admissible independent calibration/import data, or replace a payload with an independently justified source computation; do not count calibration as a theorem proof.

Existing claim-row navigation (no new edges): ` OPH-MODEL-P-RESPONSE-EXECUTION `, ` OPH-MODEL-NATIVE-OBSERVER-EXECUTION `, ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Sources:

- [Empirical imports PR-14/15](../docs/registers/PREMISE_REGISTER_V3.md).
- [Hadronic import policy](../docs/policies/HADRON.md).
- [Cited CODATA target, inherited source pins and retrospective comparison scope](../model/predictions/endpoint_audit.py).
- [Independent calibration and sampling premises of the finite instrument test](../model/predictions/protocol.py).
- [Single-carrier clock comparison, spatial covariance and laboratory attachment contract](../model/native_clock/COMPARISON.md).
- [Measured host-clock correspondence and timing-control scope](../model/native_clock/HOST_CLOCK.md).
- [Raw native record replay and held-out timing comparison](../model/native_clock/host_clock.py).
- [Recorded native choices, currents and host-counter block endpoints](../model/native_clock/host_clock_records.json).
- [Calibration, held-out blocks and workload-control verdicts](../model/native_clock/host_clock_results.json).
- [Classical branch bound and retrospective public Bell comparison](../model/predictions/BELL_RECORDS.md).
- [Complete published-data selection and exact null-tail replay](../model/predictions/bell_results.json).

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

**Constructed:** The live model connects native repair, recorded observations, feedback and checkpoint continuation. Full native C/E source compilation, common interacting work, coherent native/work composition and active encoded memory have explicit but distinct interface contracts. The native action, finite equilibration, local event clocks and field memory are consequences of their specified common source operations. The supplied scalar/gauge/matter benchmark executes a complete joint history with occupation and temporal bounds on its declared preparations and detectors.

Theorem attachment is checked on the actual objects: native signed-difference blocks supply exchange/coherence rates; local endpoint supports and tensor splits supply dependency-cone and no-signalling hypotheses; stopped instruments retain remote evolution and conditional references; positive state expectations supply measurement Gram constraints. These are compatible finite and analytic attachments. They do not assert that every theorem branch uses the same reference, observer algebra or action, or that an all-theorem physical model has been completed. The field/refinement attachment uses a common finite-density preparation, clock, record dictionary and feedback. Its canonical extension removes mandatory product independence while preserving exact fixed-number correlations. Record-derived temporal representations and joint incident/movement clocks have explicit domain maps; their distinct sectors, observer regions and empirical contracts remain separate.

**Remaining input:** Architecture and admitted observer services may legitimately be design inputs. Unique selection of an implementation is a separate question. Physical interpretation and quantitative service fidelity require independent evidence. Complete native record laws, benchmark field laws and laboratory apparatus must be connected through compatible observations and interventions before their physical claims can be combined. Source-screen scaling, relativistic field behavior, responsive gravity and observed matter remain open. A programmed modular probe is not automatically the evolution of the continuing observer. No cross-reference transfers a certificate to a different state, detector, action, clock or approximation domain.

**Closure condition:** For each promoted claim bind its actual assumptions and imports, all transitive antecedents and model-added hypotheses to compatible objects with scoped evidence. Keep unreviewed bindings visible and reject circular discharge.

Existing claim-row navigation (no new edges): ` OPH-MODEL-NATIVE-OBSERVER-EXECUTION `, ` OPH-MODEL-OBSERVER-RECORD-CAPACITY `, ` OPH-MODEL-P-RESPONSE-EXECUTION `, ` OPH-MODEL-N4-LOCAL-TICK-CERTIFICATE `, ` OPH-MODEL-SCALAR-RESPONSE-EXECUTION `, ` OPH-MODEL-COMMON-SOURCE-EXECUTION `, ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Sources:

- [Per-claim assumptions and edge classifications](../claims/claim_registry.yaml).
- [Named assumption meanings](../claims/assumption_dictionary.md).
- [Canonical register scope](../docs/registers/PREMISE_REGISTER_V3.md).
- [N4 shared bare-occupation response and finite memory certificate](../code/native_response/README.md).
- [Live observer and repair-law implementation](../model/README.md).
- [Reproducible package receipt and source verification](../model/verify.py).
- [Counterexamples, identifiability and compatible constructions](../model/CONSTRUCTION.md); exact source anchor: ` ## Counterexamples and identifiability `.
- [Executed scalar detector and native feedback contract](../model/observer/scalar_response/README.md).
- [Binary instrument, history law and inherited refinement bound](../model/observer/scalar_response/DERIVATION.md).
- [Retained scalar-response execution and controls](../model/observer/scalar_response/results.json).
- [Composed finite source apparatus and retained physical obligations](../model/observer/source_scalar/DERIVATION.md).
- [Current executable source-compiled two-site construction](../model/observer/source_scalar/experiment.py).
- [Operational completion criteria and effective-description scope](../model/CONSTRUCTION.md); exact source anchor: ` ## Completion target and acceptance criteria `.
- [Complete ideal two-oscillator output bound on the executed construction family](../model/observer/source_scalar/cutoff.py).
- [Specified finite candidate comparison and complete-data decision rule](../model/predictions/protocol.py).
- [Retrospective source-output mismatch with inherited physical-completion boundary](../model/predictions/endpoint_audit.py).
- [Finite interacting common-source model and explicit completion boundary](../model/common_source/README.md).
- [Finite schedule-independent recovery and solver/hardware distinction](../model/common_source/HISTORY_REPAIR.md).
- [Native geometry readout without physical position or momentum identification](../model/common_source/SLOW_GEOMETRY.md).
- [Conditional mathematical limit of the chosen smooth observation protocol](../model/common_source/SMOOTH_REFINEMENT.md).
- [Infinite interacting feedback signal and complete error proof](../model/common_source/INFINITE_HISTORY.md).
- [Current native clock, field and federation attachment](../model/OBSERVER_RELATIVE_RECOVERY.md).
- [Implemented inputs and separate physical selection obligations](../model/SOURCE_SELECTION.md).
- [Source-derived action and local quantum observation contract](../model/native_action/README.md).
- [Exact complete-observation memory and local clocks](../model/joint_clock_process/README.md).

<a id="re-27"></a>

### RE-27: Executed common scalar, gauge, matter and record process

**Constructed:** The existing coherent proper-code source executes a common finite Hamiltonian H_hop+E^2+H_scalar+(sqrt(2)/4)Q0(nL0-nR0). Its four-state ordered-CAR/integer-flux sector obeys both operator Gauss constraints and even parity. Source words prepare the scalar projection, evolve all sectors jointly and perform two complete quantum pointer instruments; the first record controls a dressed hopping operation and the final scalar-sign record controls native repair. Exact postcommit record adapters distinguish classical public readback from block-algebra retention and satisfy the dual fixed-observable criterion. All four quantum outcomes and all 1056 native/quantum histories are retained, with probabilities attached to their terminal states. Ideal temporal refinement has a rational all-input O(n^-2) channel bound and an independent exact d=2 signal enclosure. Disconnected and disabled-feedback controls test the cross-sector causal response. GNS and modular readouts restrict the same complete terminal state to declared observer algebras. A source-executed d=2 modular clock extends its actual branches without changing committed records. Supplied scalar factors admit CP regional restriction and Hamiltonian-derived interaction support, with the charged matter/link sector kept shared. The same ideal interacting two-oscillator family has a quantitative pre-sharp-readout occupation bound and a strong-limit proof for its complete retained sharp-sign instrument and native continuation. These are separate from the finite temporal budget and from physical service error. Its twelve-node event constraints admit local parent-to-child repair with a unique schedule-independent complete history under node fairness, and an exact ideal real-residual descent inequality. The existing native slow projector is connected to every retained quantum/native path through exact committed-current displacement identities. An independent uniform thirty-seam calibration computes the rank-three response quotient. The original E/F controller retains its rank-two displacement span, while a declared final-bit F/G selector after E reexecutes all native continuations and spans the rank-three carrier. A separate emitted smooth Ramsey instrument replaces the final sharp-sign detector and retains its own bit dictionary and all quantum outputs. The full chosen initial-state instrument has an O(1/d) occupation bound, and n(d)=d^2 product steps give a joint ideal temporal/occupation limit. Its smooth ensemble preserves a finite fourth number moment and adds exactly 1/8 to infinite-model energy in the declared dimensionless convention. Finite source energy changes are checked against finite Kraus operators, not the infinite CCR identity. A supplied fixed faithful full-work Gibbs reference aligns its modular automorphisms with continuous propagation of this same Hamiltonian at t=-beta*s. The actual prepared/recorded state remains distinct. Complete source gauge-number probe continuations preserve both old record marginals and retain all eight outputs, with the existing temporal approximation budget against the continuous target. For zero native preparation current, the same declared interacting action has an exact infinite-oscillator feedback contrast [-0.024045670143, -0.024043450955]. Full-Hamiltonian Taylor powers on expanding support charge initial projection, every support discard and exact coefficient rounding; both first records and both infinite half-line readout outcomes retain their unnormalized quantum states. The final detector uses Hermite half-line integrals rather than the sign of a compressed coordinate matrix. A separate non-Gaussian exponential-number-moment and projected-Duhamel proof bounds the pre-readout occupation refinement uniformly over arbitrary matter/reference inputs with the fixed free coupled scalar vacuum. These are mathematical infinite-history and cutoff certificates, not emitted infinite source words. Strong resolvent convergence of the compressed coordinates and absence of an infinite-coordinate atom at zero also prove convergence of the complete finite sign instruments. For the mathematically defined exact palindromic family with d=2^k and n=d^2 steps per interval, the two-interval temporal error is O(1/d); the contrast therefore tends to a value in the certified negative interval and is eventually negative. This gives no explicit sufficient cutoff or additional native compiler execution beyond d=2,4. A separate complete-output certificate compares the actual zero-current sharp histories, retaining both matter records and both sign outcomes. Exact half-line overlaps enclose the d=2 ideal 512-step source-word distance from the infinite output in [0.454336,0.455973] with feedback and [0.456824,0.458461] without it. Source-bound factorial-metric eigenbasis residuals certify finite sign projectors and complete continuous-output enclosures at d=32,64,128. A separate witness proves every output in the named rectangular d=256 occupation code has distance above 0.122 for either prepared history. For the vacuum control the optimal such code approximation scales as Theta(d^(-1/4)), even though both sign probabilities can be exactly one half. These are fixed-state occupation-code results, not bounds for every abstract encoding, arbitrary-matter instruments or physical apparatus; post-sharp clock continuation is not included. The canonical finite program executes a retained four-read smooth continuation and its native attachment at d=2,4, with all 4224 joint outcomes. Its prefix Kraus identities, resumable simulator state and stagewise moment/Duhamel bound connect repeated smooth measurements to later model-time evolution and a gauge clock probe. Finite temporal, occupation and missing physical-service budgets remain distinct. A separate coherent native attachment composes source-compiled full-carrier C/E repairs with the retained four-read smooth program. The first native current selects the gauge kick and one minus the second smooth bit selects the final seam, using source-executed port routing. Complete raw outputs include padding and aborts. The product-channel identity holds for arbitrary native/work/reference inputs; the default field vacuum retains its own occupation hypotheses. Earlier classical native attachments and their counts remain separate evidence. The active encoded variant uses the existing four-configuration qubit, a source encoder and a committed-record-controlled recovery after both native events. Its first seam is (0,1); the final seam is (2,3) for second smooth bit zero and (0,1) otherwise. Before the terminal decoded read, each ideal corrected branch acts as logical identity tensored with its intended work map, up to the explicit record amplitude. The passive mode uses the same encoded input and work-controlled seam policy without recovery. Full output padding and aborts remain retained, and the vacuum occupation premise is unchanged. A separate exact prepared-history certificate retains all sixteen smooth-continuation outputs for every native current, through later feedback and the gauge probe. Its maximum-current ideal-product output distance is below 0.264092 at d=2 and 0.053956 at d=4 with 512 splits per interval. The blank-gauge and coupled scalar-vacuum input is fixed; complete native and protected-memory continuations preserve this bound with arbitrary independent native/reference input. It does not cover arbitrary initial gauge states or correlated native/work inputs, and supplies no floating-execution or physical-service bound.

**Remaining input:** This finite action, coupling, preparation, CAR/link representation, coherent code service and model-time schedule are declared construction choices. The earlier classical attachment retains the native primitive. The coherent attachment executes the encoded C/E interface, but its physical load-to-code identification, initialization, coherent transfer, storage and readout accuracy remain supplied. The prepared pure state is not the unconstrained A3 tracial optimizer. The programmed local modular instrument and the separate full-work Gibbs-reference alignment have different reference states and scopes. The latter supplies model-time conversion for the chosen Hamiltonian. The reference and beta are declared choices, with no Gibbs-preparation service claimed; physical clock and temperature calibration and BW identification require separate evidence. Measurement and feedback interventions are not generated by one autonomous modular group. Numerical and physical service errors, spatial refinement for this interacting family, non-Abelian/chiral matter, physical parameter and metrological attachment, responsive gravity and independent empirical tests remain separate. The interacting family uses its own moment proof; no free Gaussian evolution bound is transferred, and the limit leaves the number of scalar sites fixed. No full Lean OperationalObserver or federation instantiation and no simultaneous realization of incompatible finite hypotheses is asserted. The auxiliary history solver is not a source-compiled global repair, its real potential is not an integer physical-tower rank, and the fixed action and boundary remain supplied. Native response coordinates are dimensionless; rank-three calibration alone identifies neither physical momentum nor an extended three-dimensional field. The complete smooth-instrument bound does not repair the old sharp-sign rate or inherit that protocol's modular-clock attachment. Its uniform moment-based occupation bounds are vacuous at the retained small registers, and its infinite two-oscillator/gauge limit supplies no observed-matter, spatial-field or gravitational identification. The exact infinite sharp-sign response and its complete finite-output comparison are state-specific to zero native preparation current; the latter does not certify arbitrary matter preparations or every abstract finite encoding. The exponential-number theorem is uniform over arbitrary matter/reference inputs only with the fixed scalar vacuum and before final readout. Neither that theorem nor the exact sharp-history arithmetic extends the appended local modular clock or full-work Gibbs-reference probes, or emits an infinite source program. The active encoded variant does not restore arbitrary full-carrier coherence or certify fault-tolerant physical memory, and its interventions do not inherit the passive native geometry/clock calibration.

**Closure condition:** Preserve this actual common state, operation and record dictionary when using its occupation proofs, exact infinite sharp response, smooth instrument, finite local modular clock or fixed full-work Gibbs reference. Extend spatial families and attach physical clocks, matter and gravitational observation interfaces with their own controlled errors. Reuse each theorem only after exhibiting its hypotheses on these same objects. Different interactions require new approximation bounds; compiling or rigorously approximating a supplied action does not derive or empirically validate that action.

Explicit assumption-name navigation: ` declared_common_gauge_CAR_scalar_action_and_preparation `, ` declared_common_source_two_record_intervention_protocol `, ` declared_native_current_to_common_quantum_history_attachment `, ` declared_coherent_proper_code_services_and_analogue_source_controls `, ` declared_common_source_observer_algebra_and_access `, ` declared_common_source_programmed_modular_instrument `, ` declared_common_history_constraint_repair `, ` declared_native_slow_geometry_readout `, ` declared_common_source_smooth_ramsey_instrument `, ` declared_common_source_thermal_reference `, ` declared_common_source_coherent_native_attachment `, ` declared_common_source_active_encoded_native_memory `.

Existing claim-row navigation (no new edges): ` OPH-MODEL-COMMON-SOURCE-EXECUTION `, ` OPH-MODEL-NATIVE-OBSERVER-EXECUTION `, ` OPH-SOURCE-COHERENT-CODE-COMPLETION `, ` OPH-SOURCE-LOCAL-FERMIONIC-REALIZATION `.

Sources:

- [Scope of the live coupled finite construction](../model/common_source/README.md).
- [Common carrier, ideal temporal bound and record composition](../model/common_source/DERIVATION.md).
- [Exact CAR/flux sector and native gauge operations](../model/common_source/gauge_matter.py).
- [Same-register scalar coupling and source execution](../model/common_source/scalar_coupling.py).
- [Postcommit projector and conditional-expectation adapters](../model/common_source/records.py).
- [All native and quantum branches in one complete history](../model/common_source/native_join.py).
- [Exact rational finite detector-response enclosure](../model/common_source/finite_signal.py).
- [Dual fixed-observable criterion reused by the checkpoint adapter](../Lean/Dynamics/ProtectedCharge.lean).
- [Existing coherent M6 source](../code/m1_source_realization/README.md).
- [Existing signed fermionic source controls and separate clock scope](../code/m1_fermionic_source/README.md).
- [Same-history observer representation, modular clock and regional interpretation](../model/common_source/OBSERVER_READOUTS.md).
- [Finite GNS and modular diagnostics on actual observation states](../model/common_source/operational.py).
- [Complete source-executed modular clock instrument](../model/common_source/modular_clock.py).
- [CP regional restriction and common-Hamiltonian interaction support](../model/common_source/spatial.py).
- [Interacting occupation proof and separate sharp-readout convergence scope](../model/common_source/OCCUPATION_REFINEMENT.md).
- [Exact rational pre-readout moment and Duhamel bound](../model/common_source/refinement.py).
- [Fixed-history compatibility, local residual descent and fair-schedule proof](../model/common_source/HISTORY_REPAIR.md).
- [Local constraints using the actual source operators and complete quantum outputs](../model/common_source/history_repair.py).
- [Native response calibration and original and extended binary-controller scopes](../model/common_source/SLOW_GEOMETRY.md).
- [Exact slow projector and committed-current displacements on all joint paths](../model/common_source/slow_geometry.py).
- [Distinct smooth detector, complete limit and energy scope](../model/common_source/SMOOTH_REFINEMENT.md).
- [Actual source pointer and complete finite energy ledger](../model/common_source/smooth_protocol.py).
- [Full recorded-instrument occupation and temporal certificates](../model/common_source/smooth_refinement.py).
- [Infinite interacting feedback signal and complete error proof](../model/common_source/INFINITE_HISTORY.md).
- [Exact full-H expanding-support history and retained polynomial replay](../model/common_source/infinite_history.py).
- [Retained infinite response, cutoff comparisons and polynomial hashes](../model/common_source/infinite_results.json).
- [Independent Hamiltonian, normalization, support-loss and history checks](../model/common_source/test_infinite_history.py).
- [Uniform interacting pre-readout occupation-refinement proof](../model/common_source/OCCUPATION.md).
- [Rational arbitrary matter/reference cutoff bounds](../model/common_source/occupation.py).
- [Independent interacting moment and boundary-generator checks](../model/common_source/test_occupation.py).
- [Infinite half-line instrument and finite-sign correction](../model/common_source/READOUT.md).
- [Exact Hermite half-line weights and directed pi enclosures](../model/common_source/halfline.py).
- [Independent half-line and complete-readout checks](../model/common_source/test_halfline.py).
- [Adversarial exact receipt, polynomial export and CLI checks](../model/common_source/test_infinite_receipt.py).
- [Fixed thermal reference and complete actual observation-history scope](../model/common_source/THERMAL_CLOCK.md).
- [Actual source probes under the full-work thermal reference](../model/common_source/thermal_clock.py).
- [Retained thermal-reference execution receipt and replay](../model/common_source/thermal_bridge.py).
- [Complete sharp-output certificate](../model/common_source/sharp_output.py).
- [Complete sharp-output certificate](../model/common_source/sharp_pointer.py).
- [Complete sharp-output certificate](../model/common_source/sharp_output_results.json).
- [Independent sharp-output and adversarial controls](../model/common_source/test_sharp_output.py).
- [Independent sharp-output and adversarial controls](../model/common_source/test_sharp_pointer.py).
- [Independent sharp-output and adversarial controls](../model/common_source/test_sharp_rate.py).
- [Canonical finite observer and smooth continuation theorem](../model/common_source/COMPOSITION.md).
- [Complete four-read source/native history and budgets](../model/common_source/composition_results.json).
- [Complete coherent native/work composition and retained physical boundaries](../model/common_source/QUANTUM_NATIVE.md).
- [Source-compiled native attachment and complete-output replay](../model/common_source/quantum_native_results.json).
- [Executed native port routing on the complete data code](../model/common_source/native_routing.py).
- [Complete source-native and retained-work instrument](../model/common_source/quantum_native.py).
- [Active encoded recovery, continuing quantum memory and passive loss control](../model/common_source/ENCODED_MEMORY.md).
- [Source encoder, complete conditional recovery and inverse decoding](../model/common_source/encoded_memory.py).
- [Complete active and passive encoded observer producer](../model/common_source/encoded_observer.py).
- [Encoded observer complete-output comparisons and source costs](../model/common_source/encoded_observer_results.json).
- [Prepared complete smooth-history certificate and native/encoded lifting scope](../model/common_source/SMOOTH_HISTORY.md).
- [Exact sixteen-output continuation enclosure and canonical replay](../model/common_source/smooth_history.py).
- [All-current prepared finite/infinite comparisons with separate temporal errors](../model/common_source/smooth_history_results.json).
