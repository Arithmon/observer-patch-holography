# OPH reverse-engineering table: model and theorem premises

As of 2026-10-18. Generated from [the navigation map](../tracking/premise_register.json)
by `python tools/build_model_theorem_map.py`. Edit the input, then regenerate.

This map connects the canonical premise identifiers and active model obligations to their current construction evidence. It is navigation coverage, not an extraction of every paper or Lean antecedent and not a proof of discharge. Alternative branches retain their own hypotheses; cross-references do not change claim-registry dependency roles.

**Construct the simplest shared repair process satisfying the declared observational contract.** Architecture, observer services and the complete repair law can supply inputs. Their consequences need not be postulated again. A theorem attachment must identify the same states, instruments, retained records, access rules and scheduling domain; matching component formulas alone is insufficient. Unique selection of an implementation is separate from executing that implementation and testing its predictions.

The live construction is `model/`. Its native C/E instruments determine local exchange and coherence laws, finite stationary preparation, field and event-clock statistics, and complete quantum observation bounds. The finite transfer logarithm, local Liouville generator and modular action on a chosen algebra/reference are distinct constructions. An observer-relative clock is an accessible record process; no preferred universal physical clock is required. Supplied scalar, gauge and matter actions remain conditional benchmarks with their own occupation, temporal and spatial error contracts.

**Completion concerns observations.** The intended effective quantum-field, classical and gravitational descriptions must reproduce compatible joint statistics and intervention responses within controlled regimes. Fields, metrics and manifolds need not be fundamental objects. Native source execution does not by itself establish physical spacetime or identify a laboratory apparatus. Public oscillator comparisons retain their calibration, postselection, coordinate and sampling limitations; numerical agreement, mathematical consistency and code fidelity are separate from physical identification.

The entries distinguish finite executions, analytic bounds, imported theorem hypotheses and remaining physical attachments. Separate actions, preparations, detectors and approximation domains keep separate certificates. Reusing a theorem requires a valid map of its actual hypotheses, not a count of available results or the simultaneous imposition of incompatible branches.

The [current physical bridge](SOURCE_SELECTION.md#current-bridge-to-physical-observations) distinguishes finite observer field/clock results, the remaining effective QFT and spacetime identification, and actual public-data tests. No federation is required for the completed single-carrier calculations. Issue #1077 retains its original q=3 scalar-word deliverable; smaller executed histories do not replace it.

The exhaustive contract inventory below covers every registered claim and tracked Lean source module. Restricted same-model bindings retain all canonical inputs and require explicit domain maps before composition. Missing model realizations remain visible; coverage and typed compatibility do not certify scientific entailment or physical completion.

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
| Unique ordinary `assumptions` names | 819 |
| Unique `imported_results` names | 898 |
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

**Constructed:** Finite conditional repair and candidate equalizers exist. Weighted averaging factors the supplied scalar stiffness exactly. A current-record lift couples native C to the existing six-level response at both seam endpoints, preserving its full load marginal and disjoint-seam locality. N4 instead conjugates binary C/E by one coherent load-response exchange and supplies private-state influence on fixed bare occupation probabilities. At the declared full tracial reference, N4 supplies a shared bare-occupation geometry and two-attempt kinetics protocol; its finite memory excess rules out native heat closure at the tested couplings. The executable N1 observer composes current observation, retained classical history, record-only prediction and feedback for two events, with all benchmark paths retained. Its detector extension maps the accepted closure P through the existing electromagnetic coupling to a binary record and native feedback. A minimal twelve-block address carrier preserves the complete quantum channel and joins binary native/N4 repair with collective response in an abstract scalar-seam A1–A3 source. The two laws share the full quantum reference and normalized geometry but differ in certified finite kinetics. A separate coherent phase/readout protocol distinguishes C and E exactly while preserving their common classical load transitions; it supplies a candidate comparison without requiring a uniqueness theorem. The common-source joint instrument lets the actual first native current prepare the gauge/scalar experiment and lets its final scalar record control the second native repair, retaining all 1056 histories and their quantum terminal states. In that earlier classical attachment, native M78 integer repair remains a supplied primitive while the quantum apparatus uses the existing M6 source controls. The independent uniform-seam forgotten-record C/E channels have exact quantum gap (5-sqrt(5))/60 on both the binary shell and the full total-two carrier, by an exact signed-difference reduction. Their sharp binary coherence rates differ despite identical classical geometry and relaxation. The coherent native attachment places source-compiled C or E repairs before and after the same retained four-read work program. Its first decoded current controls the gauge kick and its second smooth bit selects the final native seam. The complete joint channel admits arbitrary native/work/reference inputs; default executions use the declared scalar vacuum. This supplies a common finite interface without selecting a unique native record law or a physical encoding. A separate coherent-continuation experiment checks all 3960 native binary arcs and a bijective unread-record reversal. With supplied uniform coherent event embedding A and reflection, W=S(2AA*-I) has compression A*W^n A=T_n(K), and the actual centered load columns obey KX=X(I-L/60). This gives a compressed-amplitude wave recurrence, not a mean-load or probability law. Two continuation phases share the same first coherent event but have two-step Ramsey plus probabilities 1 and 113/150 for the specified basis input. Record-label phases invisible to complete committed C/E CP branches change the compressed spectral trace by -1/900. Thus the native instrument does not select the coherent extension; collision states obstruct this same-domain arc reversal on the full carrier. Coherent preparation, unread records, reflection and interferometric readout are additional controls; no source pulse realization, spatial refinement or physical clock is claimed. Proper-carrier covariance and the declared real coherent-continuation criteria leave the two executed alternatives. A source-executed record-phase experiment has the same complete committed C/E instrument but a fixed coherent read with probabilities 3/4 versus 1/4; commitment before read gives 1/2 for both. Exact native diagonal-block equivalence makes C/E indistinguishable through work-only instruments and common-record feedback, including arbitrary initial native/work/reference correlations. Coherent native access can distinguish them. A separate cut-metric obstruction rules out identifying fixed occupation-diagonal binary-pinch mixtures directly with the supplied quadratic number-dephasing rates; it is not an empirical rejection of C or E. Complete current-law identification executes nine independent scalar quadratures and reconstructs all nine real coordinates of the classified 78-state coarse definite-current instrument. It includes collision coherence invisible on the binary sector, source-compiled probes and continuing reference histories. This determines the model instrument only within the stated family and probe interface.

**Remaining input:** For the conditional-resampling branch, choose one state space, faithful reference and visible fibres. The conservative canonical extension and dissipative repair are different dynamics; a faithful thermal reference is not mandatory for every pure-state candidate. The full-model candidate must link its repair marginal, reversible port response and record clock through shared states and interventions; the scalar benchmark alone does not supply this. The one-way current-record lift adds local response factors and a supplied coupling. For binary N4, collective controls supply a complete abstract response grammar on the central-address carrier; independent-factor or address-controlled reversible probes are excluded from that grammar. N4 changes the load marginal and current-label meaning; its supplied coupling and the physical admissibility of the shared bare-occupation preparation/readout require independent justification. Exact all-input preservation of the original reduced load channel forbids unconditional internal-to-load feedback. The retained rational CPTP control has feedback difference 144/625 and a scalar-only private commutant, so the restricted single-product exclusion cannot be extended to all number-preserving conjugations (#1089). Reversible words and chart paths are fixed independently of measured records; admitted feedback uses the repair/reset templates. Unrestricted record-conditioned reversible words add 144 independent address-controlled response directions and are excluded from this abstract grammar. The coherent native attachment implements its declared finite encoding; source compilation alone supplies neither a physical native-load interface nor independently calibrated preparation, code-transfer, storage and readout errors. A source-selection argument must constrain cross-record phases and coherent continuation, or identify observables independent of those choices. Selection is relative to the admitted observer interface: scalar branch phases invisible throughout committed histories need no separate operational selection there, while coherent access makes their full implementation part of the repair-rule obligation. Self-read cannot distinguish observationally equivalent fine refinements or uncommitted dilations. Calibration requires the declared coherent apparatus, family membership and independently chosen sampling design; instrument identification is not a theorem that selects nature's source or an observer-independent eternal truth.

**Closure condition:** Exhibit one local process and explicit restriction/intervention maps supplying the intended repair, carrier-response, reversible-current and clock contracts. Reuse the existing theorem implications once these interfaces hold; derive each extra choice or retain it as an explicit model hypothesis. Separate matched component fixtures do not establish a common source.

Explicit assumption-name navigation: ` faithful_common_reference_law `, ` complete_repaired_visible_datum_fibre `, ` declared_uniform_independent_seam_schedule_and_forgotten_records `, ` declared_common_source_coherent_native_attachment `.

Existing claim-row navigation (no new edges): ` OPH-MODEL-NATIVE-OBSERVER-EXECUTION `, ` OPH-MODEL-P-RESPONSE-EXECUTION `, ` OPH-MODEL-COMMON-SOURCE-EXECUTION `.

Sources:

- [Finite conditional repair](../Lean/Thermodynamics/FiniteConditionalRepair.lean).
- [Canonical premises and conditional repair-law selection](../docs/AXIOM_REFERENCE.md).
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
- [Exact C/E diagonal-block intertwining and correlated adaptive work-observer replay](../model/identifiability_checks/rule_observability.py).
- [Binary-pinch rate obstruction to direct number-dephasing identification](../model/identifiability_checks/dephasing_observability.py).
- [Complete coarse-current identification and operational equivalence boundary](../model/repair_identification/README.md).

<a id="re-02"></a>

### RE-02: Confluence, acceptance and protected observations

**Constructed:** The normal-form and consensus implications are proved for their explicit repair and transactional contracts; finite protected-memory examples are constructed. The endpoint-current lift inherits C/E completion by seamwise unitary conjugation; positive seam mixtures on its finite binary shell converge to their common fixed algebra after records are discarded, with protected local block populations. The same finite completion and network-positivity arguments transfer to N4 on the binary shell; no new damping mechanism is needed. The finite observer explicitly preserves committed classical record values through feedback while allowing disturbance of coherence between record values. The common-source record adapter instantiates exact postcommit checkpoint projector preservation, block and classical public conditional expectations, and all-input controlled Kraus intertwiners. The distinction between a new write on a blank tape and subsequent preservation is retained; no universal confluence or full observer-federation interface follows. The same finite source apparatus has a twelve-node acyclic compatibility system defined before its local update rule. Parent-to-child repair reaches its unique complete history under every node-fair schedule from arbitrary initial assignments; five covering rounds suffice. Ideal contractivity gives a strictly decreasing real residual potential, while finite recovery uses the separate depth/fairness induction. All complete terminal vectors and defined deterministic observer extractions are independent of the solver schedule. Complete classical native histories have an exact partial-permutation/projector form: every initial load population is recoverable from final source and tape, while classical commitment removes cross-history coherences. Record-only E forests identify the full or binary input in eight or seven attempts; a C stopping controller does so almost surely with mean 24 or 21. Every unfinished prefix retains its raw probability. Bounded record-only C policies retain an entangled branch on this carrier; this is compatible with the doubled-forest exclusion of every passive exact code. The coherent attachment preserves committed record-prefix probabilities while retaining the native and work quantum output within each branch. It composes the complete C/E source instruments, including every raw padding and abort outcome; a first-event abort halts with the prepared field retained. This is an instrument-composition result, not new confluence, full-state recovery after classical commitment or a native global-history solver. A separate active encoded variant compiles recovery after each native event while preserving the committed tape and the continuing work state. Its logical-code branch identity applies with arbitrary references; a matched passive schedule demonstrates that retaining records alone does not preserve the logical phase. The same conditional finite-DAG argument has an executed density-matrix realization for all 130 endpoints of the supplied anharmonic apparatus, preserving all four detector outcomes and reproducing the existing 125 comparison scores without a new fit. Fair schedules agree, unfinished schedules fail closure, and a different supplied action also closes with different predictions. This is the same empirical comparison through local history repair, not an independent confirmation or a native C/E implementation of the history solver.

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
- [Complete-record recovery, encoded-qubit countercontrol and exact identifying schedules](../model/observer/RECORD_RECOVERY.md).
- [Source-bound exact history, recovery and stopping replay](../model/observer/results/record_recovery.json).
- [Complete coherent native/work composition and retained physical boundaries](../model/common_source/QUANTUM_NATIVE.md).
- [Source-compiled native attachment and complete-output replay](../model/common_source/quantum_native_results.json).
- [Active encoded recovery, continuing quantum memory and passive loss control](../model/common_source/ENCODED_MEMORY.md).
- [Local density-history completion of the supplied anharmonic apparatus and its action control](../model/identifiability_checks/consensus_observability.py).
- [Repair, consensus and their empirical identifiability boundaries](../model/identifiability_checks/README.md).

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

**Constructed:** Finite carriers, authenticated source order and regional reconstruction contracts exist. The canonical executed-history constructor supplies an ancestry rank; the abstract semantic-log interface requires one. The executable native observer specifies a finite source/tape carrier, a record-only controller menu and repeatable committed-prefix reads. The electromagnetic detector distinguishes the one-bit actuator algebra from its complete recorded trajectory. A separate binary native/N4 source admits twelve primitive central address blocks with a unital inclusion of its full quantum working algebra. Refresh followed by incident-star repair decodes to the original channel on arbitrary correlated inputs and repeated histories. The current-record regional protocol reconstructs known diagonal fields and their common-prefix overlaps on the actual native carrier. Its pullback is a finite commutative GNS isometry. This does not supply the homomorphic restriction maps or other quantum-net hypotheses of FiniteCausalObserverNet; shared marginals alone do not reconstruct joint laws.

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
- [Actual regional record gluing and quantum-net attachment boundary](../model/relative_observer_geometry/README.md).

<a id="re-05"></a>

### RE-05: Selected scalar action and relative capacities

**Constructed:** The declared golden scalar action admits the stated capacity and conductance reconstruction under its positive diagonal symmetrizer hypotheses. The two-site oscillator apparatus supplies mass one and spring kappa=9/32 and is a different action family. For native binary repair, the actual mean response fixes the relative quadratic form without a fitted oscillator coefficient. More fully, the C/E endpoint instruments derive H_native=sum_e(I-Phi_e) on operator space, with positive local terms, exclusion exchange and candidate-dependent coherence killing. The exact discrete transfer is I-H_native/M. Independent rate-one Poisson seam attempts give exp(-s H_native), where s is auxiliary expected attempts per seam. The finite transfer logarithm, this local Liouville action and a state-vector Hamiltonian are distinct objects; each has its own observation contract. The binary record-table family gives a self-adjoint positive I-Phi exactly when k, g and h are real. Its Hilbert-Schmidt skew part is a Hamiltonian commutator exactly when Im(k)+Im(h)-Im(g)=0. A phase-chosen C instrument has the explicit interacting component K_theta=-(sin(theta)/2) sum_e P_odd,e under the declared Poisson scheduler, with the same classical selective histories and compatible fixed space. The phase, physical time and Lorentzian field identification remain unselected.

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
- [Action classification and phase-dependent interaction](../model/repair_action_audit/README.md).

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

**Closure condition:** Emit and execute an adequate native implementation and add reference-stable physical initialization, pulse, service, storage and readout errors to the mathematical budget. Preserve complete records and the comparison order: actual code implementation, preparation and evolution errors for that action, then its complete readout. The executed finite two-site programs and the interacting infinite-history certificate provide their respective occupation evidence; neither emits the original large common-code recipe or an infinite apparatus. A laboratory claim requires independently established physical-service and clock attachments on the same process. Issue #1077 retains the original q=3, d=8192 onsite scalar-word deliverable. The later executed finite native and two-site protocols are separate milestones; accepting one as a replacement would require an explicit scope revision, not an inferred closure.

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

The coherent screen observer executes a complete logical-qubit channel at native incident-clock count two on the actual matched level-three/five screens. A degree-two exclusion dual and fixed signed-difference operator compute all pair/reference populations and coherences, including mask-controlled phase feedback and all 21 capped records. Spatial, cutoff and directed-arithmetic bounds put the computed table within 7.779e-10 of either completed full-screen output, uniformly over logical inputs and untouched references under the independent Bernoulli-environment preparation. The two completed source channels differ by at most 2.731e-10. C and E share the classical record law but have completed Bell-output distance 101/968; an unnormalized retained C block certifies entanglement. Preparation and phase services remain declared, and this finite comparison does not establish an all-level or physical-time limit.

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
- [Complete coherent screen-observer channel, proof and resource budgets](../model/coherent_screen_observer/README.md).
- [Native coherent screen comparison and independent full-state controls](../model/coherent_screen_observer/results.json).

<a id="re-12"></a>

### RE-12: Source clock, physical time and execution cost

**Constructed:** Source count and finite execution-cost results retain their declared clock and service domains. Native local event instruments sum actual unresolved repairs, retain remote evolution and include finite-cap timeouts. The joint clock/field process supplies exact reduced-memory identities and a four-mark completed-event closure on its specified one-particle support; conditional reference states and controller memory remain. Equal completed clock statistics need not imply equal phase observations. Independent identity holds cancel from completed stopped laws without fixing a laboratory duration.

A chosen faithful algebra/reference defines its modular action. The native Palm reference, reconstructed transfer reference and supplied full-work Gibbs reference are separate pairs. Source probes execute the relevant finite modular automorphisms on stated accessible codes. For sigma_beta proportional to exp(-beta H), Hamiltonian segments have t=-beta*s, with reference choice, beta and full-work access supplied. Exact native covariance defects prevent equating particular event operations with those modular flows. Actual interventions and feedback are retained as instruments, rather than treated as one autonomous modular group. Pulse, native-attempt, event and load/export counts are separately accounted resources. Instrument-valued clock conversion equals direct native stopping, with every remote-count coefficient and conditional quantum/reference output retained. The passive two-A-tick renewal cycle has B-count mean 2 and variance 16; specified B-entry feedback changes the law. A uniform two-event survival bound certifies completion, while capped runtime execution retains all timeouts. This derives a stochastic relative-clock law, not a deterministic conversion from equal mean rates. On the 66-state binary total-two sector of one twelve-port carrier, actual local movement clocks have asymptotic count rate 1/6 per carrier movement from any initial density in that sector. An exact Poisson corrector supplies finite count-error bounds without stationary preparation or an independent clock. Their long-run covariance distinguishes all three graph-distance shells. A positive seconds calibration changes no such dimensionless law; empirical attachment requires a matched source apparatus and held-out paired records. A recorded software attachment executes the same 264532 native attempts in each of two runs against a host monotonic counter. One baseline block calibrates port 2; three held-out blocks measured with port 3 have timing errors 2.3% to 7.8%, within the fixed 25% engineering tolerance. Host-only workload preserves all native records and rejects that calibration in every control block. Replay checks the seeded choices and source currents, not timestamp authenticity or clock traceability; physical proper time and modular flow are not inferred. Native circulation has an action in one port's own tick count. A local-tick Palm herald succeeds with probability 5/198; failure remains recorded. The stopped occupation identity gives six carrier movements per local tick, exact cycle covariance P/5 and triangle variance 3m/5 after m local ticks. Finite fourth moments and adjacent squared-current dependence distinguish the law from independent Gaussian blocks. The long-count Gaussian increment action has coefficient 5/2 on the cycle subspace, without replicas or an independent clock. Carrier-wide current access and the uniform scheduler are declared; this is not a relativistic action. Fixed-window quantum execution preserves herald failures and zero-tick outcomes, while record-controlled seam selection changes the finite law. Exact incident-clock field correlations reconstruct a seven-dimensional cyclic detector Hilbert space and its transfer spectrum on the elementary total-two source; independent full 66-state propagation checks withheld lags. The positive-transfer moment class is supplied by the existing theorem. Exact rank is not inferred from sampled records, and the seven-mode Gram is severely ill-conditioned. On each actual level-3/5 screen with product preparation, bounded record correlations also certify a two-vector temporal subspace and compressed transfer values near 0.88702 and 0.97079 with error below 2e-5. These compressed values are not a closed two-mode spectrum or laboratory energy calibration. The exact marked incident/movement-clock law retains zero-current local hits, remote waits, field-conditioned continuations and feedback hits. For the uniform total-two single-port source, I/K tends to 33/5; an exact Poisson corrector supplies a finite conversion-error bound indexed by fixed global movement count, uniformly over initial states. The retained host record comparison preserves baseline acceptance and workload-control rejection under that conversion uncertainty without refitting. It does not transfer the ratio or calibration to the different root-clock product/canonical-sector field protocol, identify a physical clock, or authenticate the pseudorandom entropy. A separate exact affine event/timestamp contract preserves the correlated rate/offset calibration set and tests held-out endpoints jointly without retuning. Missing timing or count-error bounds remain unresolved; uncertain integer counts use an explicitly disclosed continuous outer relaxation. Under a declared fixed attempt cadence, the native waiting distribution yields deadline-probability intervals. Synthetic controls and retained host measurements are distinguished; no native-event laboratory correspondence or traceability is established. A resource-clock refinement retains total native attempts alongside the ideal incident-clock history and attaches the source pulse invocation costs. Forgetting the added count reproduces every original conditional source matrix. A joint duration polytope provides conditional deadline bounds and a disjoint held-out test under an explicit affine timestamp hypothesis.

**Remaining input:** No preferred universal clock is required. Relating an observer's records to modular evolution, laboratory durations or effective proper time requires an operational correspondence preserving joint and conditional observations. A common mean rate or affine label comparison does not supply that correspondence. Physical service accuracy, reference preparation and clock calibration remain explicit. The scalar benchmark time and native event parameter cannot be identified merely because they occur in one program; unbounded event protocols need well-posed stopping and complete tail control. Choosing seconds is a unit convention. The dynamical obligation is a compatible joint record law and an independently validated apparatus correspondence, rather than an absolute duration derived without calibration. The controller attempt-count record is an extra read service. This refinement neither propagates the raw compiled arrays nor inherits coarse-current calibration bounds for its finer record. Native attempts, pulse costs and physical seconds remain distinct; no measured timing evidence is supplied by the fixture.

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
- [Paired event/clock calibration and native waiting-time push-forward](../model/laboratory_clock/README.md).
- [Ideal-history resource refinement and conditional two-counter duration attachment](../model/laboratory_clock/resource_clock.py); exact source anchor: ` def resource_law( `.

<a id="re-13"></a>

### RE-13: Physical geometry and responsive metric

**Constructed:** Finite source geometry, golden coordinates and conditional Lorentz/spacetime constructions exist. The canonical source includes its adopted oriented boundary and complete-response clauses. The common-source attachment reuses the existing native twelve-port carrier and proves I-L_ico/60 by exact enumeration of every total-two state, seam and placement coin under a uniform all-seam calibration. The lowest positive Laplacian band has rank three and Gram G=4P3, with an exact icosahedral coordinate witness. Every actual joint history receives a dimensionless displacement readout from its committed currents and declared seams, checked against working-load changes without giving the controller extra state access. The original E/F controller has a rank-two displacement span. A declared binary F/G continuation after the first E event reexecutes every native branch while retaining the same quantum history, and its three admitted seam displacements span the rank-three carrier. The existing native integer rule also gives pathwise load continuity and fair-coin conditional current mean (x_a-x_b)/2 and variance parity(x_a-x_b)/4. Exact enumeration of all 4680 state/seam/coin events builds the full 78-state transition kernel for the existing uniform geometry calibration and proves its closed mean law I-L_ico/60. The binary 66-state uniform distribution is stationary with continuing load fluctuations. All nonconstant mean multipliers are positive contractions, so this fixed linear mean readout cannot be an undamped wave under a positive scalar clock rescaling. The actual record-controlled E/F protocol is executed as a distinct control; no new clock, Hamiltonian, physical field identification or continuum premise is supplied. Exposure-normalized even-load current laws recover the same unweighted rank-three structural Gram without a uniform schedule or fair coins. The actual slowest response under a doubled seam weight has multiplicity two. The full inverse-response Gram instead has rank eleven, so choosing the slow-band metric is an operational readout choice, not a proof of physical bulk dimension. These exact single-carrier constructions require no federation. On the same uniform single carrier, asymptotic local-tick covariance independently recovers the native distance shells from a continuing movement-count history. Deterministic even-load currents also permit exact reconstruction from 360 executed native shots, one per preparation/seam, without a probability-law oracle. This removes native sampling uncertainty for that protocol, not laboratory preparation or detector errors. Native circulation has an action in one port's own tick count. A local-tick Palm herald succeeds with probability 5/198; failure remains recorded. The stopped occupation identity gives six carrier movements per local tick, exact cycle covariance P/5 and triangle variance 3m/5 after m local ticks. Finite fourth moments and adjacent squared-current dependence distinguish the law from independent Gaussian blocks. The long-count Gaussian increment action has coefficient 5/2 on the cycle subspace, without replicas or an independent clock. Carrier-wide current access and the uniform scheduler are declared; this is not a relativistic action. Fixed-window quantum execution preserves herald failures and zero-tick outcomes, while record-controlled seam selection changes the finite law. A disclosed total-one preparation supplies exact native first-arrival records in the receiver incident-attempt clock. Mean arrival counts are 22, 28 and 30 for graph-distance shells one through three; the resulting passage metric has centered Gram rank eleven. Complete arrival and timeout laws retain remote repairs. A separate total-two source witness shows a two-seam intervention response at the first receiver incident tick, disproving a one-seam-per-local-tick identification while preserving native-attempt locality.

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
- [Single-carrier passage geometry and local-clock causal control](../model/physical_space_time/README.md).

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

A retrospective comparison consumes the actual four ordered P-ACT-LB LCDM CAMB posterior chains: Planck plus ACT DR6 primary spectra, ACT/Planck lensing and DESI 2024 BAO. Official half-chain burn-in leaves multiplicity 508985; the conditional fixed tilt 1-P/48 lies at empirical lower-tail mass 3114/508985. Exact decimal/rational replay retains per-chain and omitted-chain sensitivity. This is an exposed posterior-position diagnostic under the released cosmological model and priors, not a p-value, likelihood ratio, independent-row estimate or FZ-15 prospective verdict.

A finite seven-to-nine-mode scalar extension derives added couplings from the same quartic integral, accounts for the high-blank-induced mass shift, executes native source words and bounds the complete low-mode output after tracing the added modes. Their state continues between measurements, with an explicit reset control and an interaction-specific odd-excitation witness. This supplies one bounded bandwidth extension under a declared action and occupation code, not relativistic locality, an ultraviolet limit or a physical timing map. One continuing native occupation history has a fixed-graph Gaussian fluctuation limit whose covariance and precision follow from the actual binary generator. The incident-clock correction is explicit and independently checked against the 66-state kernel. A separately supplied periodic square family has fixed-band Gaussian convergence in the iterated observation-count then spatial-mesh limits.

**Remaining input:** The native discriminator needs a justified laboratory map and physical service error bounds; no deposited dataset is identified with that native protocol. The public oscillator means do not authenticate accepted shot counts, common preparation or all coordinate conventions. A statistical rejection needs a valid sampling/selection model and calibrated instruments. Mathematical or source-array error is not hardware fidelity. Retrospective witnesses and comparisons cannot be relabeled as preregistered observations, and a general quantum consistency test is not an OPH-specific prediction. The tilt relation uses an alpha-calibrated P and a separate conditional primordial embedding. No theorem here identifies it with the runnable observer field or clock. A physical parameter restriction can be tested without changing SM/GR dynamics, but its source ancestry and comparison uncertainty must be retained. Full occupation-vector readback needs initial disclosure and every seam/current record, including remote records. A port-only observer reads only its own occupation. The fluctuation action is not a Lorentzian dynamical action; the torus family is not an established continuum of the source screens, and no uniform joint limit or physical identification follows.

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
- [Public-data candidate eligibility and conditional comparison scope](../docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md).
- [Actual ACT combined posterior projection and exact weighted comparison](../code/cosmology/tilt_posterior/README.md).
- [Retained actual-data tilt comparison without a prospective verdict](../code/cosmology/tilt_posterior/results.json).
- [Controlled finite bandwidth extension and matched low-mode observations](../model/interacting_band_extension/README.md).
- [One-history occupation field, clock correction and separate torus limit](../model/physical_space_time/OCCUPATION_FIELD.md).

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

An executed single-carrier protocol joins current-reconstructed fields, incident clocks, record-controlled seam repair and terminal populations for two overlapping observers. Direct native branch enumeration and a source-derived capped transfer agree on all 5768 output words, including 5079 zero-probability words and timeout mass 369/625. The unread completed transfer is positive in its specified weighted inner product; individual capped maps use a transported positive cone. This is an exact diagonal observation correspondence, without a certificate for arbitrary quantum inputs, a relativistic QFT or physical identification.

The separate coherent overlapping-observer protocol omits initial occupation disclosure and establishes the complete native source/reference correspondence using exact remote maps, effects and inverse on all 4356 operator matrix units, followed by composition with the admitted local services. A CPTP remote encoding and its restricted-cone inverse preserve the declared incident-clock, feedback and timeout experiment, including quantum references. An independent native endpoint oracle rejects separately matched clock/detector marginals. The separate compiled clock history emits preparation, phase/mixing and terminal-pointer words, retains all 512 code addresses and 1656 clock/timeout/complement-abort outputs, and propagates original binary arrays exactly. Four complete histories, including a compiled Bell preparation, agree with the ideal correspondence within 1e-11; serialized continuation preserves every quantum entry. The maximum apparatus budgets are 51 or 54 pulses and twice as many load/export services, plus five native attempts. Native repair, encoding and physical service fidelity remain declared interfaces. The supplied scalar branch separately adds an omitted Fourier pair with induced low-mode mass matching and a finite-code reduced-output bound through continuing measurement and feedback. These results have different source/action/clock domains and do not identify them physically. The complete claim and Lean-module inventory retains every canonical hypothesis and marks absent same-model mappings open. Restricted model connections check all seven declared domain coordinates, require sourced maps for changes, retain the whole canonical input inventory and reject circular fact dependencies. The calibration-to-incident-history connection restricts the identified 78-state current family to its invariant binary sector and derives the cap-two expected exposure 408/125, including feedback, before adding the independent source compilation error.

**Remaining input:** Architecture and admitted observer services may legitimately be design inputs. Unique selection of an implementation is a separate question. Physical interpretation and quantitative service fidelity require independent evidence. Complete native record laws, benchmark field laws and laboratory apparatus must be connected through compatible observations and interventions before their physical claims can be combined. Source-screen scaling, relativistic field behavior, responsive gravity and observed matter remain open. A programmed modular probe is not automatically the evolution of the continuing observer. No cross-reference transfers a certificate to a different state, detector, action, clock or approximation domain. The single-carrier passage, action and clock packets provide explicit domain maps and countercontrols without federation. An integrated replay can check those maps without removing their preparation, scheduler, interface or physical-attachment premises. Existing public-data comparisons must stay attached to their own observable dictionaries until a derivation connects them to that same process. Typed contract checks validate the stated composition ledger, not its mathematical entailment. The finite A4 full-premise witness is a synthetic branch; it is not a source-derived realization of the common physical tower. Unbound theorem contracts retain their own state, action, clock and refinement premises.

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
- [Current finite-to-physical obligations and original scalar milestone scope](../model/SOURCE_SELECTION.md#current-bridge-to-physical-observations).
- [Single-carrier space, action, clock and federation audit](../model/PHYSICAL_AUDIT.md).
- [Executed common field/clock/current history and exact effective correspondence](../model/observer_correspondence/README.md).
- [Controlled finite bandwidth extension and matched low-mode observations](../model/interacting_band_extension/README.md).
- [Coherent native correspondence with arbitrary references and joint clock controls](../model/quantum_observer_correspondence/README.md).
- [Source-compiled preparation and intervention in the continuing clock/field history](../model/compiled_clock_history/README.md).
- [Exact original-array execution, complete-output and resource certificates](../model/compiled_clock_history/results.json).
- [Common-current calibration propagated through the same complete compiled history](../model/compiled_clock_history/law_calibration.py); exact source anchor: ` def complete_bound( `.
- [Explicit synthetic A4 premise inhabitant and its source-selection boundary](../Lean/Tower/FullPremiseInhabitant.lean); exact source anchor: ` theorem branch_endpoint_unconditional `.

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

<a id="finite-contract-coverage"></a>

## Finite-contract coverage and compatible model connections

This inventory covers every canonical claim, including non-theorem rows, and
every tracked Lean source module. This is a deliberate superset of the finite
contracts; a lexical declaration count does not decide finiteness or proof validity.
A registry addition or source module addition requires an explicit inventory update.
The existing Lean builds remain responsible for kernel checks.

A restricted binding below records the named model clause and its source argument.
It does not instantiate the entire registry row. All that row's assumptions and
imports remain attached, including inputs outside the restricted clause. Claim
classifications and consumed/open/boundary roles are preserved without conversion.
A row with no binding has an open model mapping: realization of its own hypotheses
and conclusion on this common observer process is not established here. A theorem
on a different supplied fixture remains a valid theorem on that fixture.

The checker requires exact produced/consumed fact names, a noncircular connection
graph, and an explicit sourced map for every changed carrier, law, preparation,
record access, control, clock or regime coordinate. These are machine-checked
contract and coverage constraints. The source arguments, their mathematical
entailment and physical relevance require separate review; strings and source
anchors are not proofs. No entry silently discharges a canonical premise.

Coverage: **339 claims**, **659 Lean source modules**,
**12679 public theorem/lemma declarations**;
**8 restricted model bindings** and **5 explicit connections**.

### Named model domains

**native_coherent.** carrier: 66 binary total-two states with arbitrary untouched finite references; law: native C with uniform independent seam and fair coin; preparation: arbitrary compatible quantum input; executed coherent and entangled preparations; records: committed regional currents and incident ticks, complete conditional quantum output; controls: ideal number-preserving union operations on ports 0,1,5; clock: two observer incident counters; scheduled attempts for cap and error bookkeeping; regime: finite two-interval protocol, cap two; no continuum or laboratory identification.

**identified_current_family.** carrier: 78 total-two states with arbitrary spectator and reference coherences; law: nine-real-coordinate admissible definite-current family; preparation: nine coherent-pair quadratures on declared source code; records: coarse definite current and complete plus/complement probe outcome including abort; controls: source-compiled coherent preparations and reads; common public-record feedback; clock: addressed repair attempts, including zero current; regime: finite complete-current instrument identification with fixed-design probability tubes.

**compiled_incident.** carrier: 66 binary states isometrically embedded in 512 chart states; explicit complement abort; arbitrary untouched references; law: native C with uniform independent seam and fair coin; retained complement-abort extension; preparation: declared compiled coherent or entangled preparation, no initial occupation disclosure; records: committed regional currents, incident ticks, feedback records, terminal bit, timeouts and aborts with full quantum states; controls: same union operations emitted by source compiler; exact rational propagation of raw binary arrays; clock: two observer incident counters; scheduled attempts for cap and error bookkeeping; regime: finite two-interval protocol, cap two; source-array error certificate.

**calibrated_incident.** carrier: 66 binary states isometrically embedded in 512 chart states; explicit complement abort; arbitrary untouched references; law: comparison of executed raw native C with each admissible ideal current table in common tubes containing C; same independent uniform seam scheduler; preparation: executed nominal-C compiled preparation compared with its matched ideal coherent or entangled input; no initial occupation disclosure; records: committed regional currents, incident ticks, feedback records, terminal bit, timeouts and aborts with full quantum states; controls: executed nominal-C raw apparatus arrays compared through ideal C with identical ideal controls under each admissible alternative law; clock: two observer incident counters; scheduled attempts for cap and error bookkeeping; regime: finite two-interval protocol, cap two; law uncertainty plus separate source-array error.

**supplied_scalar_band.** carrier: 512 low occupation states embedded in 8192 nine-mode states; law: supplied quartic torus Hamiltonian, mass-matched low projection; preparation: declared low preparation and high-mode blank; records: complete low-record histories and reduced terminal states; high modes traced only at end; controls: compiled low instruments; high dynamics continue between reads; clock: supplied scalar model time; regime: one seven-to-nine-mode enlargement on declared finite occupation codes.

**synthetic_A4.** carrier: two public Boolean flags and one hidden Boolean bit; law: set-left and set-right accepted repairs; preparation: arbitrary raw configuration; records: public flag pair, erasing hidden bit; controls: two flag repairs with an actual commuting choice peak; clock: finite accepted rewrite order; regime: finite synthetic A4 premise inhabitant, no source selection or cross-regulator map.

**native_occupation_average.** carrier: 66 binary total-two configurations on the twelve-port source; law: common diagonal native C/E law, uniform independent seams and fair placement coins; preparation: uniform diagonal stationary fixed-total-two law with disclosed initial occupations; records: full seam/current tape for vector readback; local port-only tape supports only its scalar occupation; controls: no feedback or interventions during one continuing occupation average; clock: all native attempts or port-zero incident attempts, including zero current; regime: exact finite laws and finite-count covariance; count-to-infinity CLT on the fixed graph; separate supplied torus has a subsequent mesh limit.

**ideal_resource_clock.** carrier: 66 binary total-two source states; executed source density matrices without an added reference; law: native C with uniform independent seam and fair coin; preparation: declared coherent preparation of the ideal two-interval experiment; records: complete native observer output refined by a disclosed controller attempt count and corresponding pulse cost; controls: ideal union operations with source pulse invocation costs; an additional controller-count record service is admitted; clock: native attempt count and source pulse count; a separate common affine duration/offset hypothesis may attach timestamps; regime: finite cap-two ideal history with all timeouts; calibration and disjoint held-out timestamp intervals.

### native-coherent-instrument

Domain: **native_coherent**. The fixed native C protocol preserves every finite record and conditional source/reference state through the allowed local operations. All 4356 operator units support the structural correspondence proof.

Established restricted facts: ` complete_native_coherent_history `.

Required connected facts: None.

Canonical claims: ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Retained ordinary assumptions: ` declared_coherent_proper_code_services_and_analogue_source_controls `, ` declared_incident_clock_field_protocol `, ` declared_incident_screen_refinement_protocol `, ` declared_native_binary_federation_instruments `, ` declared_native_clock_field_and_initial_state_contract `, ` declared_native_local_access_and_scheduler `.

Retained imported results: ` Caputo_Liggett_Richthammer_exclusion_spectral_gap `, ` Gaussian_Green_covariance_and_generalized_free_field_construction `, ` Poisson_factorial_moments_and_finite_dependency_cones `, ` classical_CHSH_local_response_bound_and_published_Bell_null_contract `, ` conditional_law_sequential_coupling `, ` directed_integer_probability_enclosure `, ` finite_Kraus_Stinespring_construction `, ` finite_matrix_block_elimination_and_stopped_instrument_resolvents `, ` finite_state_martingale_functional_central_limit_theorem `, ` ordered_subsequence_count_at_incident_stopping_time `, ` positive_transfer_spectral_logarithm `, ` stationary_local_cylinder_similarity_under_commuting_interventions `, ` trace_distance_contraction_and_unitary_channel_composition `.

Boundary: The effective inverse is positive only on the transported state cone. Ideal union controls are supplied. This finite clause neither instantiates the full claim nor selects physical space, time or action.

Exact source anchors:

- [model/quantum_observer_correspondence/construction.py](../model/quantum_observer_correspondence/construction.py); ` def structural_certificate(): `.
- [model/quantum_observer_correspondence/construction.py](../model/quantum_observer_correspondence/construction.py); ` def execute( `.

### complete-current-identification

Domain: **identified_current_family**. Nine scalar quadratures determine the classified coarse current instrument on all 78 states. The exact CP tests and reference-uniform bound epsilon <= min(1,7 delta) concern this admitted family and interface.

Established restricted facts: ` full_current_instrument_error_bound `.

Required connected facts: None.

Canonical claims: ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Retained ordinary assumptions: ` declared_coherent_proper_code_services_and_analogue_source_controls `, ` declared_incident_clock_field_protocol `, ` declared_incident_screen_refinement_protocol `, ` declared_native_binary_federation_instruments `, ` declared_native_clock_field_and_initial_state_contract `, ` declared_native_local_access_and_scheduler `.

Retained imported results: ` Caputo_Liggett_Richthammer_exclusion_spectral_gap `, ` Gaussian_Green_covariance_and_generalized_free_field_construction `, ` Poisson_factorial_moments_and_finite_dependency_cones `, ` classical_CHSH_local_response_bound_and_published_Bell_null_contract `, ` conditional_law_sequential_coupling `, ` directed_integer_probability_enclosure `, ` finite_Kraus_Stinespring_construction `, ` finite_matrix_block_elimination_and_stopped_instrument_resolvents `, ` finite_state_martingale_functional_central_limit_theorem `, ` ordered_subsequence_count_at_incident_stopping_time `, ` positive_transfer_spectral_logarithm `, ` stationary_local_cylinder_similarity_under_commuting_interventions `, ` trace_distance_contraction_and_unitary_channel_composition `.

Boundary: Coarse observations do not identify fine record refinements or coherent dilations. Binary-only probes cannot identify the collision coordinate n. Family membership, apparatus errors and fixed sampling design remain inputs; synthetic records do not select nature's law.

Exact source anchors:

- [model/repair_identification/core.py](../model/repair_identification/core.py); ` def reconstruct( `.
- [model/repair_identification/core.py](../model/repair_identification/core.py); ` def instrument_bound( `.
- [model/repair_identification/README.md](../model/repair_identification/README.md); ` **Theorem 2.** `.

### compiled-incident-history

Domain: **compiled_incident**. The source compiler supplies the preparation, union intervention and detector in the same complete incident-clock experiment. Every 512-state output, timeout and complement abort is retained. The cap-two expected native repair count is 408/125, including feedback.

Established restricted facts: ` complete_compiled_output_error_bound `; ` incident_attempt_expectation `.

Required connected facts: ` complete_native_coherent_history `.

Canonical claims: ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Retained ordinary assumptions: ` declared_coherent_proper_code_services_and_analogue_source_controls `, ` declared_incident_clock_field_protocol `, ` declared_incident_screen_refinement_protocol `, ` declared_native_binary_federation_instruments `, ` declared_native_clock_field_and_initial_state_contract `, ` declared_native_local_access_and_scheduler `.

Retained imported results: ` Caputo_Liggett_Richthammer_exclusion_spectral_gap `, ` Gaussian_Green_covariance_and_generalized_free_field_construction `, ` Poisson_factorial_moments_and_finite_dependency_cones `, ` classical_CHSH_local_response_bound_and_published_Bell_null_contract `, ` conditional_law_sequential_coupling `, ` directed_integer_probability_enclosure `, ` finite_Kraus_Stinespring_construction `, ` finite_matrix_block_elimination_and_stopped_instrument_resolvents `, ` finite_state_martingale_functional_central_limit_theorem `, ` ordered_subsequence_count_at_incident_stopping_time `, ` positive_transfer_spectral_logarithm `, ` stationary_local_cylinder_similarity_under_commuting_interventions `, ` trace_distance_contraction_and_unitary_channel_composition `.

Boundary: The compiled raw CP arrays need not be exactly trace preserving; their mass defect is charged. Source numerical error is separate from physical hardware fidelity. Pulse cost does not advance the native incident clock.

Exact source anchors:

- [model/compiled_clock_history/history.py](../model/compiled_clock_history/history.py); ` def execute( `.
- [model/compiled_clock_history/receipt.py](../model/compiled_clock_history/receipt.py); ` def bound( `.
- [model/compiled_clock_history/law_calibration.py](../model/compiled_clock_history/law_calibration.py); ` def exposure( `.

### calibrated-complete-incident-law

Domain: **calibrated_incident**. Restrict the identified full-carrier current map to the invariant binary sector and propagate uncertainty through the ideal stopped experiment. Nominal C must belong to the tubes. Expected exposure gives ideal half-trace error min(1,(408/125) epsilon). Adding the executed nominal-C source-array discrepancy bounds that raw output against each alternative ideal law.

Established restricted facts: ` calibrated_complete_incident_history_bound `.

Required connected facts: ` full_current_instrument_error_bound `; ` complete_compiled_output_error_bound `; ` incident_attempt_expectation `.

Canonical claims: ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Retained ordinary assumptions: ` declared_coherent_proper_code_services_and_analogue_source_controls `, ` declared_incident_clock_field_protocol `, ` declared_incident_screen_refinement_protocol `, ` declared_native_binary_federation_instruments `, ` declared_native_clock_field_and_initial_state_contract `, ` declared_native_local_access_and_scheduler `.

Retained imported results: ` Caputo_Liggett_Richthammer_exclusion_spectral_gap `, ` Gaussian_Green_covariance_and_generalized_free_field_construction `, ` Poisson_factorial_moments_and_finite_dependency_cones `, ` classical_CHSH_local_response_bound_and_published_Bell_null_contract `, ` conditional_law_sequential_coupling `, ` directed_integer_probability_enclosure `, ` finite_Kraus_Stinespring_construction `, ` finite_matrix_block_elimination_and_stopped_instrument_resolvents `, ` finite_state_martingale_functional_central_limit_theorem `, ` ordered_subsequence_count_at_incident_stopping_time `, ` positive_transfer_spectral_logarithm `, ` stationary_local_cylinder_similarity_under_commuting_interventions `, ` trace_distance_contraction_and_unitary_channel_composition `.

Boundary: Raw arrays under alternative repair laws are not executed. The homogeneous relabelled seam law is a premise, not inferred by probing one seam. Raw-array maps need not preserve trace, so their event-weight allowance is twice their half-trace bound. No physical law or duration is selected; admitted services and sampling assumptions remain explicit.

Exact source anchors:

- [model/compiled_clock_history/law_calibration.py](../model/compiled_clock_history/law_calibration.py); ` def complete_bound( `.
- [model/compiled_clock_history/law_calibration.py](../model/compiled_clock_history/law_calibration.py); ` def from_tubes( `.
- [model/compiled_clock_history/law_calibration.py](../model/compiled_clock_history/law_calibration.py); ` def from_counts( `.

### supplied-first-band-control

Domain: **supplied_scalar_band**. Projecting the supplied nine-mode quartic action onto its declared high blank derives the low mass adjustment. The finite off-code coupling controls complete low output with untouched references and continuing high dynamics.

Established restricted facts: ` first_band_mass_matching_and_complete_low_output_bound `.

Required connected facts: None.

Canonical claims: ` OPH-MODEL-SHARED-SCALAR-OBSERVATIONS `.

Retained ordinary assumptions: ` declared_coherent_proper_code_services_and_analogue_source_controls `, ` declared_fixed_band_scalar_observation_protocol `, ` declared_scalar_band_extension_protocol `.

Retained imported results: ` canonical_lattice_scalar_quantum_simulation `, ` exact_rational_arithmetic_with_directed_root_and_pi_enclosures `, ` quadratic_form_Leibniz_identity_for_smooth_measurement `, ` symmetric_product_formula_integral_Taylor_remainder `, ` weighted_number_moment_and_projected_Duhamel_cutoff_bounds `.

Boundary: This action, occupation code, blank and model time are supplied. The one-band R-squared T-squared estimate is a reduced-output bound, not a full joint-state estimate, a native-C Hamiltonian derivation or an all-band continuum result.

Exact source anchors:

- [model/interacting_band_extension/README.md](../model/interacting_band_extension/README.md); ` R^2 T^2 `.
- [model/interacting_band_extension/construction.py](../model/interacting_band_extension/construction.py); ` def reduced_output_distance( `.

### finite-A4-existence-witness

Domain: **synthetic_A4**. An explicit finite branch combines all five A4 endpoint premises, including a nontrivial choice peak and hidden representatives. Its unconditional corollary quantifies over this constructed branch.

Established restricted facts: ` synthetic_branch_public_endpoint_unique `.

Required connected facts: None.

Canonical claims: ` OPH-PUBLIC-WORLD-ENDPOINT-FINITE `.

Retained ordinary assumptions: ` finite_inhabited_public_world_presentation `, ` injective_OPH_to_tower_public_readback_adapter `, ` repair_output_and_enabledness_quotient_congruence `, ` source_and_cross_regulator_endpoint_realization_deferred `, ` terminating_confluent_semantically_complete_labelled_repair `.

Retained imported results: ` existing_OPH_repair_and_gauge_theorems `, ` finite_kernel_quotients `, ` well_founded_confluent_rewriting `.

Boundary: The source does not select this synthetic branch. Its word unconditional removes the branch theorem's hypotheses; it does not establish a physical common tower, a unique Universe or a native C confluence theorem.

Exact source anchors:

- [Lean/Tower/FullPremiseInhabitant.lean](../Lean/Tower/FullPremiseInhabitant.lean); ` theorem branch_endpoint_unconditional `.
- [Lean/Tower/FullPremiseInhabitant.lean](../Lean/Tower/FullPremiseInhabitant.lean); ` theorem branchQuotientCompatible `.

### native-occupation-fluctuation-field

Domain: **native_occupation_average**. The unchanged native binary process yields the covariance and precision of its one-history occupation-average Gaussian limit. Exact 66-state transitions check the field closure, stationary lag law and clock correction. A separate supplied torus family admits fixed-band Gaussian convergence after the time limit.

Established restricted facts: ` native_occupation_covariance_and_clock_corrected_precision `.

Required connected facts: None.

Canonical claims: ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Retained ordinary assumptions: ` declared_coherent_proper_code_services_and_analogue_source_controls `, ` declared_incident_clock_field_protocol `, ` declared_incident_screen_refinement_protocol `, ` declared_native_binary_federation_instruments `, ` declared_native_clock_field_and_initial_state_contract `, ` declared_native_local_access_and_scheduler `.

Retained imported results: ` Caputo_Liggett_Richthammer_exclusion_spectral_gap `, ` Gaussian_Green_covariance_and_generalized_free_field_construction `, ` Poisson_factorial_moments_and_finite_dependency_cones `, ` classical_CHSH_local_response_bound_and_published_Bell_null_contract `, ` conditional_law_sequential_coupling `, ` directed_integer_probability_enclosure `, ` finite_Kraus_Stinespring_construction `, ` finite_matrix_block_elimination_and_stopped_instrument_resolvents `, ` finite_state_martingale_functional_central_limit_theorem `, ` ordered_subsequence_count_at_incident_stopping_time `, ` positive_transfer_spectral_logarithm `, ` stationary_local_cylinder_similarity_under_commuting_interventions `, ` trace_distance_contraction_and_unitary_channel_composition `.

Boundary: The full vector needs the disclosed initial occupations and every seam/current record; it is not readable from a local incident tape. The field action describes fluctuation statistics, not a Lorentzian evolution. The torus is a different supplied architecture with an iterated time-then-mesh limit, no uniform joint limit and no continuum identification for the source screens.

Exact source anchors:

- [model/physical_space_time/occupation_field.py](../model/physical_space_time/occupation_field.py); ` def configuration_kernel( `.
- [model/physical_space_time/occupation_field.py](../model/physical_space_time/occupation_field.py); ` def covariance_at_count( `.
- [model/physical_space_time/occupation_field.py](../model/physical_space_time/occupation_field.py); ` def incomplete_tape_control( `.
- [model/physical_space_time/OCCUPATION_FIELD.md](../model/physical_space_time/OCCUPATION_FIELD.md); ` ## Finite-graph theorem and the source-derived action `.

### ideal-resource-clock-attachment

Domain: **ideal_resource_clock**. An independent endpoint replay refines each ideal native history by its total scheduled-plus-feedback attempt count. Forgetting that count reproduces every ideal conditional source matrix. Source pulse costs supply a second operation counter; a common feasible duration polytope gives conditional deadline bounds and held-out tests.

Established restricted facts: ` ideal_history_resource_refinement_and_conditional_duration_bounds `.

Required connected facts: ` complete_native_coherent_history `.

Canonical claims: ` OPH-MODEL-NATIVE-PROCESS-RECOVERY `.

Retained ordinary assumptions: ` declared_coherent_proper_code_services_and_analogue_source_controls `, ` declared_incident_clock_field_protocol `, ` declared_incident_screen_refinement_protocol `, ` declared_native_binary_federation_instruments `, ` declared_native_clock_field_and_initial_state_contract `, ` declared_native_local_access_and_scheduler `.

Retained imported results: ` Caputo_Liggett_Richthammer_exclusion_spectral_gap `, ` Gaussian_Green_covariance_and_generalized_free_field_construction `, ` Poisson_factorial_moments_and_finite_dependency_cones `, ` classical_CHSH_local_response_bound_and_published_Bell_null_contract `, ` conditional_law_sequential_coupling `, ` directed_integer_probability_enclosure `, ` finite_Kraus_Stinespring_construction `, ` finite_matrix_block_elimination_and_stopped_instrument_resolvents `, ` finite_state_martingale_functional_central_limit_theorem `, ` ordered_subsequence_count_at_incident_stopping_time `, ` positive_transfer_spectral_logarithm `, ` stationary_local_cylinder_similarity_under_commuting_interventions `, ` trace_distance_contraction_and_unitary_channel_composition `.

Boundary: The extra controller count is an admitted read service. This is a refinement of the ideal native history, not an execution of the raw compiled-array history. PR1117 coarse-current calibration bounds do not automatically control the finer resource record. The affine timing law and timestamp errors are external hypotheses; fixture timestamps are not physical measurements.

Exact source anchors:

- [model/laboratory_clock/resource_clock.py](../model/laboratory_clock/resource_clock.py); ` def counted_history( `.
- [model/laboratory_clock/resource_clock.py](../model/laboratory_clock/resource_clock.py); ` def resource_law( `.
- [model/laboratory_clock/resource_clock.py](../model/laboratory_clock/resource_clock.py); ` def deadline( `.

### Sourced domain connections

**native-coherent-instrument → compiled-incident-history** consumes ` complete_native_coherent_history `.

- carrier: Use the checked 66-to-512 isometry and retain the complement as an explicit abort branch.
- law: Embed every native code Kraus map and add the common complement-abort map; do not reinterpret padding as native states.
- preparation: Specialize the arbitrary-input identity to the declared coherent/reference preparation and charge its emitted-array error.
- records: Append complement-abort labels and the terminal detector while retaining all old words, including zero branches.
- controls: Compare the emitted raw arrays with the same ideal union operations in complete output norm.
- regime: Keep cap two and add the independently certified source-array error, including mass defect.

Source: [model/compiled_clock_history/README.md](../model/compiled_clock_history/README.md) at `` Write `P=VV*` and `Q=I-P`. ``.

**complete-current-identification → calibrated-complete-incident-law** consumes ` full_current_instrument_error_bound `.

- carrier: Restrict the same complete 78-state current instrument to the invariant 66-state binary sector; isometric embedding and common abort extension cannot increase its diamond discrepancy.
- law: Require nominal C as a feasible tube witness and use any admissible table in those same tubes, with a shared seam law at every call.
- preparation: The all-input/reference norm bound covers the incident protocol preparation without a new calibration fit.
- records: Forget inaccessible probe/private labels and retain the common current/incident-clock history; do not compare fine coin records or coherent environments.
- controls: Use identical ideal intervening controls for law comparison and account for compiled-array error separately.
- clock: Sum call error over the state-independent live-branch weights of the incident scheduler, including the feedback repair.
- regime: Apply the bound to the complete cap-two two-interval output with every timeout retained.

Source: [model/compiled_clock_history/law_calibration.py](../model/compiled_clock_history/law_calibration.py) at ` def from_tubes( `.

**compiled-incident-history → calibrated-complete-incident-law** consumes ` complete_compiled_output_error_bound `.

- law: Compare the identical compiled apparatus first with ideal C, then with each admissible calibrated ideal current table; no raw-map trace-preserving assumption is inserted.
- regime: Add the reference-uniform calibration discrepancy and compilation discrepancy on the common complete output alphabet.
- preparation: Compare the nominal-C raw preparation to the matched ideal preparation before replacing ideal C by another admissible ideal instrument.
- controls: The native-law comparison keeps ideal controls identical; the separate source-array discrepancy is for the executed nominal-C apparatus only.

Source: [model/compiled_clock_history/law_calibration.py](../model/compiled_clock_history/law_calibration.py) at ` def complete_bound( `.

**compiled-incident-history → calibrated-complete-incident-law** consumes ` incident_attempt_expectation `.

- law: Compare the identical compiled apparatus first with ideal C, then with each admissible calibrated ideal current table; no raw-map trace-preserving assumption is inserted.
- regime: Add the reference-uniform calibration discrepancy and compilation discrepancy on the common complete output alphabet.
- preparation: Compare the nominal-C raw preparation to the matched ideal preparation before replacing ideal C by another admissible ideal instrument.
- controls: The native-law comparison keeps ideal controls identical; the separate source-array discrepancy is for the executed nominal-C apparatus only.

Source: [model/compiled_clock_history/law_calibration.py](../model/compiled_clock_history/law_calibration.py) at ` def complete_bound( `.

**native-coherent-instrument → ideal-resource-clock-attachment** consumes ` complete_native_coherent_history `.

- carrier: Specialize the reference-complete native correspondence to the executed ideal source density matrices.
- preparation: Use the same coherent preparation and omit the separate entangled-reference setting.
- records: Admit a classical controller attempt-count record in the endpoint replay; equality to the earlier output holds after forgetting it. Coarse law calibration does not identify this finer instrument.
- controls: Keep the ideal operations and attach their pulse invocation costs; admit the count-record service explicitly rather than inferring it from incident ticks.
- clock: Keep attempt and pulse counters distinct; timestamp conversion additionally assumes a shared affine duration law, with calibration and held-out endpoints.
- regime: Refine the complete cap-two ideal output; do not normalize or reinterpret raw compiled-array weights as probabilities.

Source: [model/laboratory_clock/resource_clock.py](../model/laboratory_clock/resource_clock.py) at ` def resource_law( `.

### Every canonical contract

A/I count ordinary assumption/import names, not independent or unproved premises.
Each identifier points to its owning evidence; the complete names and statements
remain in the [canonical registry](../claims/claim_registry.yaml). An open mapping
does not add an antecedent to an explicit non-consumer theorem.

| Canonical claim | Restricted model mapping | A / I | Original premise roles |
| --- | --- | ---: | --- |
| [OPH-SYN-FORMAL-BASIS](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-02 |
| [OPH-OBS-PATCH-DEF](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 3 / 2 | explicit_non_consumer |
| [OPH-CARRIER-PRESENTATION-INVARIANCE](../paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-UNIFIED-TYPED-SPINE](../paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex) | Open: no same-model instantiation recorded | 7 / 5 | open: PR-03, PR-04, PR-51, PR-52, PR-53, PR-54, PR-64, PR-65 |
| [OPH-CONS-D1](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 11 / 3 | consumed: PR-01 |
| [OPH-CONS-ADAPTIVE-REPAIR-STRATIFICATION](../extra/observable_normal_forms.tex) | Open: no same-model instantiation recorded | 2 / 2 | consumed: PR-77, PR-78 |
| [OPH-CONS-CUMULATIVE-ATTEMPT-CAPACITY](../extra/observable_normal_forms.tex) | Open: no same-model instantiation recorded | 2 / 3 | consumed: PR-79, PR-80; boundary: PR-78 |
| [OPH-CONS-CAPACITY-BOUNDED-PUBLIC-ENDPOINT](../extra/observable_normal_forms.tex) | Open: no same-model instantiation recorded | 4 / 3 | consumed: PR-01, PR-79, PR-80 |
| [OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 2 / 3 | consumed: PR-81; boundary: PR-82 |
| [OPH-CONS-ADAPTIVE-PUBLIC-ENDPOINT](../extra/observable_normal_forms.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-01, PR-77 |
| [OPH-CONS-BFT-PREPARED-LOCK](../paper/appendix_B_bft_qecc_extensions.tex) | Open: no same-model instantiation recorded | 15 / 1 | explicit_non_consumer |
| [OPH-CONS-REFINEMENT-MODULI](../extra/observable_normal_forms.tex) | Open: no same-model instantiation recorded | 8 / 2 | explicit_non_consumer |
| [OPH-CONS-PROTECTED-BEHAVIOR-PROFILE](../extra/observable_normal_forms.tex) | Open: no same-model instantiation recorded | 6 / 4 | explicit_non_consumer |
| [OPH-CONS-CROSS-SOURCE-ENDPOINT](../extra/observable_normal_forms.tex) | Open: no same-model instantiation recorded | 3 / 1 | explicit_non_consumer |
| [OPH-GR-PORT-RESPONSE-COMPLETION](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 5 / 3 | explicit_non_consumer |
| [OPH-GR-D2-COLLAR-RECOVERY](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 4 / 3 | explicit_non_consumer |
| [OPH-GR-D5](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 8 / 5 | consumed: PR-23, PR-24, PR-25, PR-26, PR-27, PR-28 |
| [OPH-GR-D3-CAP-H3](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 5 / 4 | explicit_non_consumer |
| [OPH-GR-D3-H3-LOCALIZATION](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 7 / 4 | explicit_non_consumer |
| [OPH-GR-D3H-GEOMETRY-PRODUCER](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 6 / 3 | explicit_non_consumer |
| [OPH-GR-D4-NULLNET-STANDARDNESS](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 8 / 5 | explicit_non_consumer |
| [OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 7 / 10 | consumed: PR-16; open: PR-15, PR-52 |
| [OPH-GR-FINITE-LOCAL-DOMAIN](../paper/screen_microphysics_and_observer_synchronization.tex) | Open: no same-model instantiation recorded | 4 / 4 | explicit_non_consumer |
| [OPH-GR-D4C-LOCAL-STRESS](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 6 / 2 | consumed: PR-24, PR-25 |
| [OPH-GR-D4D-ENTROPY-BRIDGE](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 3 / 2 | explicit_non_consumer |
| [OPH-GR-D5A-ABSOLUTE-EINSTEIN](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 17 / 2 | consumed: PR-23, PR-24, PR-25, PR-26, PR-27, PR-28 |
| [OPH-GR-E2E-BRANCH-ENTRY](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 16 / 1 | consumed: PR-23, PR-24, PR-25, PR-26, PR-27, PR-28 |
| [OPH-GR-SOURCE-DIRECTION-EINSTEIN-COMPOSITION](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 5 / 3 | explicit_non_consumer |
| [OPH-GR-D6-CAPACITY](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 9 / 4 | explicit_non_consumer |
| [OPH-GR-D6-READBACK-NOGO](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 5 / 5 | boundary: PR-51 |
| [OPH-GR-D6-A5-FINITE-CONTROL](../paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex) | Open: no same-model instantiation recorded | 4 / 3 | explicit_non_consumer |
| [OPH-GR-D6-HORIZON-RECORD](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 3 / 2 | open: PR-42 |
| [OPH-GR-DS-MU2-IDENTITY](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex) | Open: no same-model instantiation recorded | 1 / 2 | explicit_non_consumer |
| [OPH-GR-DS-CAPACITY-TRANSFER](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex) | Open: no same-model instantiation recorded | 2 / 2 | explicit_non_consumer |
| [OPH-GR-DS-SHOCK-SIGN-ATTACHMENT](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex) | Open: no same-model instantiation recorded | 5 / 2 | open: PR-42, PR-52 |
| [OPH-GR-DS-DISCRETE-SHOCK-SPECTRUM](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex) | Open: no same-model instantiation recorded | 2 / 2 | explicit_non_consumer |
| [OPH-REPAIR-SHOCK-DOMAIN-BOUNDARY](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex) | Open: no same-model instantiation recorded | 1 / 1 | explicit_non_consumer |
| [OPH-COSMO-SCREEN-SPECTRUM](../paper/tex_fragments/SCREEN_SPECTRUM_THEOREMS.tex) | Open: no same-model instantiation recorded | 10 / 5 | explicit_non_consumer |
| [OPH-WZ-STRICT-1L-POLE-MAP](../paper/tex_fragments/STRICT_ONE_LOOP_WZ_POLE_MAP_PORT.tex) | Open: no same-model instantiation recorded | 4 / 3 | open: PR-47, PR-48, PR-49, PR-50, PR-54 |
| [OPH-SM-Q1-LOCAL-G6](../paper/tex_fragments/QFT_Q1_Q4_WZ_THEOREM_PORT.tex) | Open: no same-model instantiation recorded | 3 / 2 | open: PR-47, PR-48, PR-49, PR-50, PR-54 |
| [OPH-SM-Q2E-CHIRAL-MEASURE-CRITERION](../paper/tex_fragments/QFT_Q1_Q4_WZ_THEOREM_PORT.tex) | Open: no same-model instantiation recorded | 1 / 2 | open: PR-47, PR-54 |
| [OPH-SM-Q2H-POSITIVE-HAMILTONIAN-SOUNDNESS](../paper/tex_fragments/QFT_Q1_Q4_WZ_THEOREM_PORT.tex) | Open: no same-model instantiation recorded | 1 / 2 | open: PR-54 |
| [OPH-SM-Q3-BV-RESTORATION](../paper/tex_fragments/QFT_Q1_Q4_WZ_THEOREM_PORT.tex) | Open: no same-model instantiation recorded | 1 / 2 | open: PR-54 |
| [OPH-SM-Q4-OS-OBSERVABLE-SECTOR](../paper/tex_fragments/QFT_Q1_Q4_WZ_THEOREM_PORT.tex) | Open: no same-model instantiation recorded | 1 / 1 | open: PR-54 |
| [OPH-SM-Q4-RESONANCE-CONTINUATION](../paper/tex_fragments/QFT_Q1_Q4_WZ_THEOREM_PORT.tex) | Open: no same-model instantiation recorded | 1 / 2 | open: PR-54 |
| [OPH-SCREEN-ECHOSAHEDRAL-SELECTOR](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex) | Open: no same-model instantiation recorded | 4 / 3 | explicit_non_consumer |
| [OPH-SCREEN-GALOIS-FRAME-BOUNDARY](../code/a5_closure/galois_port_frame_certificate.md) | Open: no same-model instantiation recorded | 4 / 3 | boundary: PR-53 |
| [OPH-A5-PRIMITIVE-PORT-SPIN6](../flagship/from_observer_consensus_to_standard_physics.tex) | Open: no same-model instantiation recorded | 6 / 3 | open: PR-15, PR-53 |
| [OPH-A5-CARRIER-CLASS-DISPERSION-BAND](../flagship/from_observer_consensus_to_standard_physics.tex) | Open: no same-model instantiation recorded | 6 / 4 | explicit_non_consumer |
| [OPH-POSITIVE-COSINE-FREQUENCY-CONTRACTION](../paper/screen_microphysics_and_observer_synchronization.tex) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-TIME-ORDER-TYPE-LEDGER](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 3 / 2 | explicit_non_consumer |
| [OPH-BOUNDED-OBSERVER-TIME-CALIBRATION](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 4 / 3 | explicit_non_consumer |
| [OPH-PUBLIC-RECORD-ALGEBRA-NOCLONING](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-02 |
| [OPH-CONSENSUS-TOWER-INTERFACE](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-02 |
| [OPH-PUBLICIZATION-DYNAMICS-FINITE](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 7 / 4 | consumed: PR-02; open: PR-15 |
| [OPH-PUBLIC-PRIVATE-DYNAMICS-FINITE](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 5 / 5 | consumed: PR-43 |
| [OPH-FINITE-QUANTUM-LIMITATION-SUITE](../extra/machine_checked_finite_event_algebras.tex) | Open: no same-model instantiation recorded | 3 / 4 | consumed: PR-02 |
| [OPH-FINITE-EXTERIOR-COMPONENT-BRIDGE](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 4 / 2 | consumed: PR-59; open: PR-47 |
| [OPH-FINITE-BORN-FRAME-RANK-GAP](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 6 / 5 | consumed: PR-02, PR-03; open: PR-04, PR-64, PR-65 |
| [OPH-THERMO-FOUR-LAW-PACKAGE](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 8 / 4 | consumed: PR-07, PR-08, PR-15 |
| [OPH-THERMO-SOURCE-REALIZATION-AUDIT](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 4 / 2 | boundary: PR-07 |
| [OPH-FINITE-TRANSPORT-GREEN-KUBO-DIFFUSION](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 6 / 3 | open: PR-07, PR-15 |
| [OPH-FINITE-LOCALITY-NOSIGNALLING](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 4 / 3 | consumed: PR-44; open: PR-52 |
| [OPH-ADAPTIVE-SCHEDULER-LOCALITY-HELPER](../Lean/ObserverPatchHolography/Locality/AdaptiveScheduler.lean) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-44; open: PR-52 |
| [OPH-FINITE-HISTORY-VARIATIONAL-HELPERS](../Lean/Variational/MechanicsAdequacySurface.lean) | Open: no same-model instantiation recorded | 7 / 5 | consumed: PR-05, PR-06, PR-45; open: PR-15 |
| [OPH-A5-SEAM-CURRENT-EDGE30](../code/a5_fingerprint/seam_current_edge_prediction.py) | Open: no same-model instantiation recorded | 14 / 3 | consumed: PR-20, PR-21, PR-22; open: PR-15, PR-53 |
| [OPH-MODAL-MAXWELL-FACTORIZATION-BOUNDARY](../paper/screen_microphysics_and_observer_synchronization.tex) | Open: no same-model instantiation recorded | 2 / 3 | open: PR-53, PR-54 |
| [OPH-SCREEN-PORT-CURRENT-INNER](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex) | Open: no same-model instantiation recorded | 8 / 4 | consumed: PR-09, PR-10, PR-11; open: PR-54 |
| [OPH-SM-SIMPLE-GUT-XY-CHANNEL-ABSENCE](../code/a5_closure/receipts/port_current_inner_reference.receipt.json) | Open: no same-model instantiation recorded | 1 / 2 | explicit_non_consumer |
| [OPH-SCREEN-SUPER-TANNAKIAN-MATTER-LIFT](../code/a5_closure/super_tannakian_matter_lift_certificate.md) | Open: no same-model instantiation recorded | 7 / 5 | consumed: PR-59; open: PR-46, PR-47, PR-48, PR-54 |
| [OPH-SCREEN-CARRIER-SUPPORT-BRIDGE](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 4 / 3 | explicit_non_consumer |
| [OPH-FINITE-CONSERVATION-WARD-PRECURSOR](../paper/screen_microphysics_and_observer_synchronization.tex) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-GAUGE-D7-D9](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex) | Open: no same-model instantiation recorded | 6 / 4 | consumed: PR-59; open: PR-35, PR-46, PR-47, PR-54 |
| [OPH-QUARK-REGISTER-CLEBSCH](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 8 / 4 | consumed: PR-59; open: PR-50 |
| [OPH-QUARK-YUKAWA-TWO-MODULUS-NONIDENTIFIABILITY](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 6 / 3 | open: PR-48, PR-50 |
| [OPH-KOIDE-CIRCULANT-IDENTITY](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 3 / 4 | consumed: PR-13; open: PR-47, PR-50 |
| [OPH-A5-REAL-AXIS-CABIBBO-NOGO](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 3 / 3 | boundary: PR-50 |
| [OPH-W5-STABILISER-POTENTIAL-BOUNDARY](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 4 / 4 | boundary: PR-50 |
| [OPH-SCREEN-A5-SM-RECOGNITION](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex) | Open: no same-model instantiation recorded | 12 / 4 | consumed: PR-09, PR-10, PR-11, PR-36, PR-59; open: PR-35, PR-46, PR-47, PR-48, PR-54 |
| [OPH-GAUGE-KINETIC-HISTORY-BINDING](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex) | Open: no same-model instantiation recorded | 4 / 2 | consumed: PR-05; open: PR-54 |
| [OPH-SM-ROUTE-IDENTIFICATION](../paper/tex_fragments/UNIFIED_OBSERVER_PHYSICS_SPINE.tex) | Open: no same-model instantiation recorded | 4 / 3 | open: PR-35, PR-46, PR-47, PR-54 |
| [OPH-GAUGE-EXTERIOR-SM-WITNESS](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-59; open: PR-47, PR-48 |
| [OPH-Q-PN-CLOSURE](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 7 / 2 | consumed: PR-17, PR-37; open: PR-38, PR-39, PR-51 |
| [OPH-Q-N-RESERVE-CANDIDATES](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 13 / 3 | consumed: PR-17, PR-51 |
| [OPH-HIER-EW](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 9 / 1 | consumed: PR-17, PR-51, PR-59 |
| [OPH-SCREEN-24-CLOCK-DETERMINANT](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 2 / 3 | open: PR-15, PR-50 |
| [OPH-SOURCE-DIMENSIONLESS-GAP](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 4 / 3 | explicit_non_consumer |
| [OPH-SCREEN-MODULAR-GEARING](../paper/screen_microphysics_and_observer_synchronization.tex) | Open: no same-model instantiation recorded | 5 / 4 | open: PR-07, PR-15, PR-42 |
| [OPH-ALPHA-PIXEL](../extra/fine_structure_constant_derivation.tex) | Open: no same-model instantiation recorded | 5 / 2 | consumed: PR-17, PR-37; open: PR-14, PR-38, PR-39 |
| [OPH-DM-CONT](../cosmology/oph_dark_matter_paper.tex) | Open: no same-model instantiation recorded | 5 / 6 | explicit_non_consumer |
| [OPH-EDGE-HEAT-KERNEL-FINITE-DIAGNOSTICS](../paper/tex_fragments/PAPER.tex) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-YM-GAP](../extra/yang_mills_gap_clay_problem.tex) | Open: no same-model instantiation recorded | 7 / 3 | open: PR-54 |
| [OPH-STRING-SELECTOR](../paper/archive/observer_patch_holography_as_string_vacuum_selector.tex) | Open: no same-model instantiation recorded | 8 / 2 | explicit_non_consumer |
| [OPH-IBM-GENERATIVE-REPAIR](../code/ibm_quantum_cloud/protocols/BLINDED_GENERATIVE_REPAIR_KERNEL.md) | Open: no same-model instantiation recorded | 7 / 4 | explicit_non_consumer |
| [OPH-THINKING-PATCH](../paper/archive/thinking_as_patch_net_fixed_point_search.tex) | Open: no same-model instantiation recorded | 5 / 3 | explicit_non_consumer |
| [OPH-PARADISE](../paper/archive/paradise_as_fixed_point_consensus.tex) | Open: no same-model instantiation recorded | 3 / 1 | explicit_non_consumer |
| [OPH-PUBLIC-WORLD-ENDPOINT-FINITE](../paper/observers_are_all_you_need.tex) | finite-A4-existence-witness | 5 / 3 | consumed: PR-01 |
| [OPH-CANONICAL-LORENTZ-MODULE](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 4 / 4 | explicit_non_consumer |
| [OPH-FINITE-CAUSAL-OBSERVER-NET](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 4 / 4 | consumed: PR-32, PR-33, PR-34, PR-44; open: PR-52 |
| [OPH-FINITE-HOLONOMY-CHARACTER-PHASE](../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex) | Open: no same-model instantiation recorded | 5 / 4 | open: PR-46, PR-52, PR-54 |
| [OPH-ALGEBRAIC-EVENT-FRAME-SOLDERING](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 5 / 4 | explicit_non_consumer |
| [OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 3 / 3 | open: PR-52, PR-54, PR-58 |
| [OPH-FINITE-UNITARY-SCATTERING-LIMIT-NO-GO](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 1 / 4 | explicit_non_consumer |
| [OPH-GEOMETRY-INVERSE-SQUARE-SHELL](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 3 / 2 | consumed: PR-18, PR-19; open: PR-52 |
| [OPH-GEOMETRY-KINEMATIC-MASS-SHELL](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 3 / 2 | open: PR-15, PR-52 |
| [OPH-GEOMETRY-LAYERED-DISCRETE-GAUSS](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 5 / 3 | consumed: PR-29, PR-30, PR-31 |
| [OPH-EM-POSITION-SPACE-MAXWELL-COULOMB](../paper/screen_microphysics_and_observer_synchronization.tex) | Open: no same-model instantiation recorded | 1 / 3 | open: PR-53, PR-54 |
| [OPH-QUANTUM-TSIRELSON-SATURATION](../extra/machine_checked_finite_event_algebras.tex) | Open: no same-model instantiation recorded | 1 / 2 | consumed: PR-02, PR-44; open: PR-52 |
| [OPH-QUANTUM-OPERATIONAL-PHASE-INTERFACE](../extra/machine_checked_finite_event_algebras.tex) | Open: no same-model instantiation recorded | 1 / 2 | consumed: PR-02; open: PR-03, PR-04, PR-64, PR-65 |
| [OPH-QUANTUM-TIME-SLICE-INTERFACE-DELIMITATION](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 1 / 2 | consumed: PR-32, PR-33, PR-34; open: PR-52, PR-58 |
| [OPH-MECH-ENRICHMENT-CHARACTERIZATION](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 2 / 2 | boundary: PR-06 |
| [OPH-THERMO-COHERENT-REFINEMENT-ENVELOPE](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 1 / 2 | consumed: PR-07, PR-08, PR-15 |
| [OPH-EM-LOCAL-FACE-MAXWELL-STATIC-COMPOSITION](../paper/screen_microphysics_and_observer_synchronization.tex) | Open: no same-model instantiation recorded | 1 / 2 | open: PR-53, PR-54 |
| [OPH-QUANTUM-TRIPLE-CARRIER-JOIN](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 1 / 2 | consumed: PR-32, PR-34; open: PR-52, PR-58 |
| [OPH-THERMO-COFINAL-SPECTRAL-TAIL-COMPOSITION](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 1 / 2 | consumed: PR-07, PR-08, PR-15 |
| [OPH-EM-TEMPORAL-MAXWELL-EVOLUTION-CONTINUITY](../paper/screen_microphysics_and_observer_synchronization.tex) | Open: no same-model instantiation recorded | 4 / 2 | consumed: PR-66; open: PR-15, PR-53, PR-54 |
| [OPH-QUANTUM-HINGE-FIBERED-JOIN-AND-SLICE-REGION-DICHOTOMY](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 3 / 2 | consumed: PR-32, PR-33, PR-34; open: PR-52, PR-58 |
| [OPH-QUANTUM-PHASE-MODEL-DETERMINATION-BOUNDARY](../extra/machine_checked_finite_event_algebras.tex) | Open: no same-model instantiation recorded | 1 / 2 | consumed: PR-02; open: PR-03, PR-04, PR-64, PR-65 |
| [OPH-GEOMETRY-INTERNAL-CLOCK-REST-FREQUENCY](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 4 / 2 | open: PR-15, PR-52 |
| [OPH-SCREEN-ELECTROWEAK-BREAKING-COMPOSITION](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 5 / 2 | consumed: PR-59, PR-67, PR-68, PR-69; open: PR-47, PR-48, PR-49, PR-50, PR-54 |
| [OPH-MECH-POLYNOMIAL-BOOST-SELECTION](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-45; boundary: PR-06 |
| [OPH-SCREEN-ASSEMBLED-STRUCTURAL-ACTION](../Lean/Screen/AssembledActionComposition.lean) | Open: no same-model instantiation recorded | 4 / 2 | consumed: PR-59, PR-67, PR-68, PR-69; open: PR-35, PR-36, PR-46, PR-47, PR-48, PR-49, PR-50, PR-54 |
| [OPH-QFT-CAUSAL-ORDER-COMPOSITION](../Lean/Geometry/CausalOrderComposition.lean) | Open: no same-model instantiation recorded | 4 / 3 | open: PR-16, PR-44, PR-52 |
| [OPH-GEOMETRY-FREE-EVOLUTION-PERSISTENCE](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 4 / 3 | open: PR-15, PR-52, PR-54 |
| [OPH-EVENTALGEBRA-SOURCE-REACHABILITY-DIAGONAL-CHSH](../Lean/EventAlgebra/SourceReachabilityDelimitation.lean) | Open: no same-model instantiation recorded | 3 / 2 | consumed: PR-02, PR-44; open: PR-52 |
| [OPH-EVENTALGEBRA-PAULI-MINIMAL-COUPLING-SQUARE](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 3 / 2 | open: PR-15, PR-47, PR-52, PR-54 |
| [OPH-SCREEN-GLOBAL-FORM-CHARACTER-DESCENT](../Lean/Screen/GlobalFormCharacterDescent.lean) | Open: no same-model instantiation recorded | 4 / 2 | consumed: PR-35, PR-59; open: PR-46 |
| [OPH-GEOMETRY-EVENT-POPULATION-CHART-INTERFACE](../Lean/Geometry/EventPopulationChartInterface.lean) | Open: no same-model instantiation recorded | 4 / 3 | open: PR-16, PR-52 |
| [OPH-SCREEN-FERMION-SECTOR-ASSEMBLY](../Lean/Screen/FermionSectorAssembly.lean) | Open: no same-model instantiation recorded | 4 / 3 | open: PR-35, PR-36, PR-46, PR-47, PR-48, PR-49, PR-50, PR-54, PR-59 |
| [OPH-GEOMETRY-COMMON-WORLD-TWO-ISLAND-WITNESS](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 5 / 3 | consumed: PR-59, PR-66, PR-67, PR-68, PR-69; open: PR-15, PR-16, PR-44, PR-47, PR-48, PR-49, PR-50, PR-52, PR-53, PR-54 |
| [OPH-QUANTUM-PHASE-OPERATION-RECORDED-DECISION-ATTAINMENT](../extra/machine_checked_finite_event_algebras.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-02, PR-04; open: PR-03, PR-52, PR-64, PR-65 |
| [OPH-GEOMETRY-REST-FIBER-SHELL-TRANSPORT](../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) | Open: no same-model instantiation recorded | 5 / 3 | consumed: PR-18, PR-19; open: PR-52 |
| [OPH-SCREEN-NEUTRAL-CURRENT-DICTIONARY](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 5 / 2 | consumed: PR-59, PR-67, PR-68; open: PR-47, PR-48, PR-49, PR-54 |
| [OPH-QFT-JOIN-RECORD-NET-MORPHISM](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 5 / 2 | consumed: PR-32, PR-33, PR-34; open: PR-52, PR-58 |
| [OPH-CONS-SCHEDULER-CLASS-CAPACITY-SPLIT](../extra/observable_normal_forms.tex) | Open: no same-model instantiation recorded | 3 / 2 | consumed: PR-60, PR-61, PR-62, PR-63 |
| [OPH-QUANTUM-OPERATIONAL-ADDITIVITY-BOUNDARY](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-02, PR-04; open: PR-52, PR-64, PR-65; boundary: PR-03 |
| [OPH-QFT-OPERATOR-LAYER-MORPHISM-BOUNDARY](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 4 / 2 | consumed: PR-32, PR-33, PR-34; open: PR-52, PR-58 |
| [OPH-GEOMETRY-COMMON-WORLD-ISLAND-BRIDGE](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 4 / 3 | consumed: PR-18, PR-19; open: PR-15, PR-52, PR-53, PR-54 |
| [OPH-SCREEN-GAUGE-SECTOR-BRACKET-COMPLETION](../Lean/Screen/GaugeSectorBracketCompletion.lean) | Open: no same-model instantiation recorded | 4 / 2 | open: PR-35, PR-36, PR-46, PR-47, PR-48, PR-49, PR-50, PR-54, PR-59 |
| [OPH-SCREEN-CHARGED-CURRENT-DICTIONARY](../paper/deriving_the_particle_zoo_from_observer_consistency.tex) | Open: no same-model instantiation recorded | 5 / 2 | consumed: PR-59, PR-67, PR-68; open: PR-47, PR-48, PR-49, PR-54 |
| [OPH-COMPUTATION-REPAIR-PROGRAM-BOUNDARY](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 3 / 2 | explicit_non_consumer |
| [OPH-COMPUTATION-FIXED-FEDERATION-PROGRESS](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 1 / 2 | consumed: PR-84 |
| [OPH-COMPUTATION-FIXED-FEDERATION-ENDPOINT](../extra/observable_normal_forms.tex) | Open: no same-model instantiation recorded | 1 / 1 | consumed: PR-84 |
| [OPH-COMPUTATION-FIXED-FEDERATION-EXECUTION-CLASSIFICATION](../paper/reality_as_consensus_protocol.tex) | Open: no same-model instantiation recorded | 2 / 3 | consumed: PR-84, PR-85 |
| [OPH-QUANTUM-LUEDERS-PHASE-INSTRUMENT](../extra/machine_checked_finite_event_algebras.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-02, PR-04; open: PR-03, PR-52, PR-64, PR-65 |
| [OPH-EVENTALGEBRA-SOURCE-PHASE-EFFECT-SELECTION](../Lean/EventAlgebra/SourcePhaseSelection.lean) | Open: no same-model instantiation recorded | 4 / 3 | consumed: PR-02, PR-04; open: PR-03, PR-52, PR-64, PR-65 |
| [OPH-EVENTALGEBRA-SOURCE-PHASE-INSTRUMENT-OUTCOME-BRIDGE](../extra/machine_checked_finite_event_algebras.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-02, PR-04; open: PR-03, PR-52, PR-64, PR-65 |
| [OPH-EVENTALGEBRA-SOURCE-PHASE-COMMON-PREPARATION-HULL](../extra/machine_checked_finite_event_algebras.tex) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-02, PR-04; open: PR-03, PR-52, PR-64, PR-65 |
| [OPH-EVENTALGEBRA-SOURCE-PHASE-BORN-WEIGHT-BOUNDARY](../extra/machine_checked_finite_event_algebras.tex) | Open: no same-model instantiation recorded | 4 / 3 | consumed: PR-02, PR-04; open: PR-03, PR-52, PR-64, PR-65 |
| [OPH-EVENTALGEBRA-PRODUCT-SPLIT-SLOT-LOCAL-CHSH](../Lean/EventAlgebra/ProductSplitSeparability.lean) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-02, PR-44; open: PR-52 |
| [OPH-GR-EINSTEIN-BRANCH-INHABITANT](../Lean/Geometry/EinsteinBranchInhabitant.lean) | Open: no same-model instantiation recorded | 3 / 4 | consumed: PR-07, PR-42; open: PR-52; boundary: PR-15 |
| [OPH-EM-SCALED-MAXWELL-STABILITY-ACTION](../paper/screen_microphysics_and_observer_synchronization.tex) | Open: no same-model instantiation recorded | 5 / 3 | consumed: PR-66; open: PR-15, PR-53, PR-54 |
| [OPH-GEOMETRY-COMMON-WORLD-INSTRUMENT-JOIN](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 6 / 4 | consumed: PR-02, PR-04, PR-07, PR-18, PR-19, PR-42, PR-66; open: PR-03, PR-15, PR-52, PR-53, PR-54, PR-64, PR-65 |
| [OPH-QFT-COLIMIT-NORMED-COMPLETION](../Lean/QFT/ColimitNormedCompletion.lean) | Open: no same-model instantiation recorded | 3 / 3 | open: PR-52, PR-58 |
| [OPH-QFT-COMPLETED-COLIMIT-GNS-INTERFACE](../Lean/QFT/ColimitGNSRepresentation.lean) | Open: no same-model instantiation recorded | 1 / 2 | open: PR-52, PR-58 |
| [OPH-QFT-SELECTED-STATE-GNS-ATTACHMENT](../Lean/QFT/ColimitSelectedStateGNS.lean) | Open: no same-model instantiation recorded | 1 / 4 | open: PR-52, PR-58 |
| [OPH-QFT-FINITE-TWO-SITE-ISING-HILBERT-DYNAMICS](../Lean/QFT/FiniteTwoSiteIsingField.lean) | Open: no same-model instantiation recorded | 1 / 3 | open: PR-15, PR-52, PR-54, PR-58 |
| [OPH-QFT-SELECTED-GNS-ISING-HAMILTONIAN-COMPOSITION](../Lean/QFT/FiniteTwoSiteIsingColimitGNS.lean) | Open: no same-model instantiation recorded | 1 / 3 | open: PR-15, PR-52, PR-54, PR-58 |
| [OPH-QFT-SOURCE-HISTORY-GNS-HAMILTONIAN-ATTACHMENT](../Lean/QFT/SourceHistoryGNSDynamics.lean) | Open: no same-model instantiation recorded | 1 / 4 | open: PR-15, PR-52, PR-54, PR-58 |
| [OPH-QFT-SOURCE-HISTORY-THREE-SLOT-LOCAL-GNS-DYNAMICS](../Lean/QFT/SourceHistoryThreeSlotLocalGNS.lean) | Open: no same-model instantiation recorded | 1 / 3 | open: PR-15, PR-52, PR-54, PR-58 |
| [OPH-DM-COLLAR-SATURATION](../Lean/ObserverPatchHolography/EinsteinBranch/CollarScaleSaturation.lean) | Open: no same-model instantiation recorded | 2 / 1 | explicit_non_consumer |
| [OPH-DM-DEEP-PROFILE-CHARACTERIZATION](../Lean/ObserverPatchHolography/EinsteinBranch/DeepProfileClosure.lean) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-QFT-SOURCE-HISTORY-STATE-PRESERVING-RECORD-NET](../Lean/QFT/SourceHistoryExpectationNet.lean) | Open: no same-model instantiation recorded | 1 / 3 | open: PR-15, PR-52, PR-54, PR-58 |
| [OPH-EM-CERTIFIED-SCALED-STEP-INSTRUMENT](../Lean/Screen/CertifiedScaledStepInstrument.lean) | Open: no same-model instantiation recorded | 3 / 2 | consumed: PR-66; open: PR-15, PR-53, PR-54 |
| [OPH-DM-PER-CUT-COMPOSITION](../Lean/ObserverPatchHolography/EinsteinBranch/PerCutCollarComposition.lean) | Open: no same-model instantiation recorded | 2 / 1 | explicit_non_consumer |
| [OPH-DM-SPARSE-DEFECT-MECHANISM](../Lean/ObserverPatchHolography/EinsteinBranch/SparseRecordDefectWitness.lean) | Open: no same-model instantiation recorded | 2 / 1 | explicit_non_consumer |
| [OPH-DM-COLLAR-PREMISE-DERIVATION](../Lean/ObserverPatchHolography/EinsteinBranch/CollarPremiseDerivation.lean) | Open: no same-model instantiation recorded | 2 / 2 | explicit_non_consumer |
| [OPH-DM-FINITE-CAP-GENERATOR-SPLIT](../Lean/ObserverPatchHolography/EinsteinBranch/FiniteCapGeneratorSplit.lean) | Open: no same-model instantiation recorded | 2 / 1 | explicit_non_consumer |
| [OPH-QUANTUM-SOURCE-BOUND-INSTRUMENT-INTERFACE](../Lean/EventAlgebra/SourceBoundInstrumentInterface.lean) | Open: no same-model instantiation recorded | 2 / 1 | consumed: PR-02, PR-04; open: PR-64, PR-65; boundary: PR-03 |
| [OPH-GEOMETRY-COMMON-WORLD-MAXWELL-CLOCK-JOIN](../Lean/Geometry/CommonWorldMaxwellClockJoin.lean) | Open: no same-model instantiation recorded | 2 / 1 | consumed: PR-02, PR-04, PR-07, PR-18, PR-19, PR-42, PR-66; open: PR-03, PR-15, PR-52, PR-53, PR-54, PR-64, PR-65 |
| [OPH-DE-FIXED-CAPACITY-WLAW](../Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean) | Open: no same-model instantiation recorded | 3 / 2 | explicit_non_consumer |
| [OPH-GW-INTEGER-K-COMB-TEMPLATE](../code/gravitation/ringdown_comb/integer_k_comb_template.py) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-CAL-SI-ANCHOR-IMPORT](../Lean/Thermodynamics/PhysicalCalibrationImport.lean) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-15 |
| [OPH-GEOMETRY-COMMON-WORLD-JOINT-ACTION](../Lean/Geometry/CommonWorldJointAction.lean) | Open: no same-model instantiation recorded | 3 / 1 | consumed: PR-02, PR-04, PR-07, PR-18, PR-19, PR-42, PR-66; open: PR-03, PR-15, PR-52, PR-53, PR-54, PR-64, PR-65 |
| [OPH-EM-DISPERSION-ARMING-INTERFACE](../Lean/Screen/DispersionArmingInterface.lean) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-70, PR-71, PR-72, PR-73, PR-74, PR-75, PR-76 |
| [OPH-QFT-SOURCE-CLOCK-CANDIDATE](../Lean/QFT/SourceClockCandidate.lean) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-15 |
| [OPH-QUANTUM-INS03-EXPORT-VALIDATOR](../code/phase_instrument_export/ins03_export_validator.py) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-64, PR-65 |
| [OPH-GEOMETRY-SCREEN-CARRIER-MAP-CANDIDATE](../Lean/Geometry/ScreenCarrierMapCandidate.lean) | Open: no same-model instantiation recorded | 2 / 2 | open: PR-15, PR-52, PR-53 |
| [OPH-DM-MATCHED-OBSERVABLE-DIAGNOSTIC](../code/cosmology/rar_deep_regime/matched_observable_diagnostic.py) | Open: no same-model instantiation recorded | 2 / 2 | explicit_non_consumer |
| [OPH-QUANTUM-INS03-V2-REPLAY-GATE](../code/phase_instrument_export/v2/PRIMITIVE_GENERATION_SEMANTICS.md) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-64, PR-65 |
| [OPH-YM-FIBER-RATE-COMPARISON](../Lean/Screen/KogutSusskindFiberRateComparison.lean) | Open: no same-model instantiation recorded | 2 / 1 | explicit_non_consumer |
| [OPH-GEOMETRY-CARRIER-DYNAMICS-COMPATIBILITY](../Lean/Geometry/CarrierDynamicsCompatibility.lean) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-15, PR-52, PR-53 |
| [OPH-DM-ROTATION-CURVE-PROFILE-LIKELIHOOD](../code/cosmology/rar_deep_regime/joint_likelihood/joint_rar_likelihood.py) | Open: no same-model instantiation recorded | 2 / 2 | explicit_non_consumer |
| [OPH-QFT-ENLARGED-TARGET-RECORD-EMBEDDING](../Lean/QFT/EnlargedTargetRecordEmbedding.lean) | Open: no same-model instantiation recorded | 1 / 1 | open: PR-15, PR-52 |
| [OPH-YM-GAUGE-ORBIT-QUOTIENT-DESCENT](../Lean/Screen/GaugeOrbitQuotientGap.lean) | Open: no same-model instantiation recorded | 2 / 1 | explicit_non_consumer |
| [OPH-QFT-SOURCE-RECURRENCE-CLOCK-PRECURSOR](../Lean/QFT/SourceRecurrenceClock.lean) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-15 |
| [OPH-GEOMETRY-CHARGE-FIXED-INTERACTION](../Lean/Geometry/ChargeFixedInteraction.lean) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-15, PR-52, PR-54 |
| [OPH-GEOMETRY-INTERNAL-ENERGY-INERTIA-PRECURSOR](../Lean/Geometry/InternalEnergyInertia.lean) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-15, PR-52, PR-54 |
| [OPH-EM-LOCAL-ENERGY-BALANCE](../Lean/Screen/LocalEnergyBalance.lean) | Open: no same-model instantiation recorded | 1 / 1 | open: PR-52, PR-53 |
| [OPH-GEOMETRY-PORT-CHARGE-MINIMAL-COUPLING](../Lean/Geometry/PortChargeMinimalCoupling.lean) | Open: no same-model instantiation recorded | 2 / 2 | open: PR-15, PR-52, PR-54 |
| [OPH-GEOMETRY-COMPOSITE-MOMENTUM-COVARIANCE](../Lean/Geometry/CompositeMomentumCovariance.lean) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-15, PR-52, PR-54 |
| [OPH-GEOMETRY-PROPER-TIME-INTERNAL-ACTION](../Lean/Geometry/ProperTimeInternalAction.lean) | Open: no same-model instantiation recorded | 3 / 2 | open: PR-15, PR-52, PR-54 |
| [OPH-EM-CARRIER-MODE-OSCILLATORS](../Lean/Screen/CarrierModeOscillators.lean) | Open: no same-model instantiation recorded | 2 / 2 | open: PR-15, PR-52, PR-53 |
| [OPH-THERMO-GIBBS-REFERENCE-ENERGY-IDENTIFICATION](../Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean) | Open: no same-model instantiation recorded | 2 / 2 | open: PR-15, PR-52 |
| [OPH-THERMO-MODULAR-ENERGY-ADDITIVITY](../Lean/Thermodynamics/ModularEnergyAdditivity.lean) | Open: no same-model instantiation recorded | 3 / 1 | open: PR-15, PR-54 |
| [OPH-GEOMETRY-WORLDLINE-HOP-TRANSPORT](../Lean/Geometry/WorldlineHopTransport.lean) | Open: no same-model instantiation recorded | 3 / 1 | open: PR-15, PR-52, PR-53, PR-54 |
| [OPH-THERMO-COLLAR-TEMPERATURE-READING](../Lean/Thermodynamics/CollarTemperatureReading.lean) | Open: no same-model instantiation recorded | 2 / 2 | open: PR-15, PR-52 |
| [OPH-EM-CARRIER-MODE-EQUIVARIANCE](../Lean/Screen/CarrierModeEquivariance.lean) | Open: no same-model instantiation recorded | 1 / 2 | open: PR-15, PR-52, PR-53 |
| [OPH-GEOMETRY-TRANSPORTED-CHARGE-FORCE-LAW](../Lean/Geometry/TransportedChargeForceLaw.lean) | Open: no same-model instantiation recorded | 3 / 1 | open: PR-15, PR-52, PR-53, PR-54 |
| [OPH-GEOMETRY-SEAM-STEP-SPEED-LIMIT](../Lean/Geometry/SeamStepSpeedLimit.lean) | Open: no same-model instantiation recorded | 2 / 1 | open: PR-15, PR-52, PR-53 |
| [OPH-GEOMETRY-SOURCE-CLOCK-RATE-ALONG-WORLDLINES](../Lean/Geometry/SourceClockRateAlongWorldlines.lean) | Open: no same-model instantiation recorded | 3 / 1 | open: PR-15, PR-52, PR-54 |
| [OPH-EM-GOLDEN-SECTOR-CHARACTERS](../Lean/Screen/GoldenSectorCharacters.lean) | Open: no same-model instantiation recorded | 1 / 2 | open: PR-52, PR-53 |
| [OPH-GEOMETRY-PROPER-LENGTH-CLOCKED-CHAIN](../Lean/Geometry/ProperLengthClockedChain.lean) | Open: no same-model instantiation recorded | 3 / 1 | open: PR-15, PR-52, PR-54 |
| [OPH-EM-GOLDEN-SECTOR-IRREDUCIBILITY](../Lean/Screen/GoldenSectorIrreducibility.lean) | Open: no same-model instantiation recorded | 1 / 1 | open: PR-52, PR-53 |
| [OPH-GEOMETRY-TIMELIKE-CLASS-FORCE-LAW](../Lean/Geometry/TimelikeClassForceLaw.lean) | Open: no same-model instantiation recorded | 3 / 1 | open: PR-15, PR-52, PR-53, PR-54 |
| [OPH-EM-CARRIER-EVOLUTION-FLOW](../Lean/Screen/CarrierEvolutionFlow.lean) | Open: no same-model instantiation recorded | 1 / 1 | open: PR-15, PR-52, PR-53 |
| [OPH-QUANTUM-INSTRUMENT-SELECTION-BY-CERTAIN-STATES](../Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean) | Open: no same-model instantiation recorded | 3 / 2 | open: PR-02, PR-03, PR-04, PR-52, PR-64, PR-65 |
| [OPH-EM-SOURCED-FIELD-CONSISTENCY](../Lean/Screen/SeamChargeContinuity.lean) | Open: no same-model instantiation recorded | 5 / 2 | open: PR-15, PR-53, PR-54 |
| [OPH-EM-GOLDEN-SECTOR-COMPLEX-IRREDUCIBILITY](../Lean/Screen/GoldenSectorComplexIrreducibility.lean) | Open: no same-model instantiation recorded | 1 / 2 | open: PR-52, PR-53 |
| [OPH-EM-FIELD-SECTOR-ENERGY-INNER-PRODUCT](../Lean/Screen/FieldSectorEnergyInnerProduct.lean) | Open: no same-model instantiation recorded | 2 / 2 | open: PR-15, PR-52, PR-53 |
| [OPH-EM-CURL-SECTOR-EIGENBASIS](../Lean/Screen/CurlSectorEigenbasis.lean) | Open: no same-model instantiation recorded | 2 / 2 | open: PR-15, PR-52, PR-53 |
| [OPH-SEAM-MAXWELL-CONTINUUM](../paper/screen_microphysics_and_observer_synchronization.tex) | Open: no same-model instantiation recorded | 5 / 4 | explicit_non_consumer |
| [OPH-SERIAL-MAXWELL-READOUT](../Lean/Screen/SerialMaxwellReadout.lean) | Open: no same-model instantiation recorded | 4 / 2 | explicit_non_consumer |
| [OPH-CONE-WHITNEY-BRIDGE](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-WHITNEY-MAXWELL-DYNAMICS](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 5 / 4 | explicit_non_consumer |
| [OPH-WHITNEY-RADIATIVE-QUANTUM](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 4 / 5 | explicit_non_consumer |
| [OPH-WHITNEY-CHARGED-MATTER](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 5 / 4 | explicit_non_consumer |
| [OPH-WHITNEY-SPATIAL-CONSISTENCY](../paper/tex_fragments/WHITNEY_SPATIAL_CONSISTENCY.tex) | Open: no same-model instantiation recorded | 6 / 3 | explicit_non_consumer |
| [OPH-WHITNEY-INTERACTING-QUANTUM](../paper/tex_fragments/WHITNEY_INTERACTING_QUANTUM.tex) | Open: no same-model instantiation recorded | 8 / 7 | explicit_non_consumer |
| [OPH-WHITNEY-CHARGED-EXECUTION](../paper/tex_fragments/WHITNEY_CHARGED_EXECUTION.tex) | Open: no same-model instantiation recorded | 6 / 3 | explicit_non_consumer |
| [OPH-WHITNEY-REAL-CONTINUUM](../paper/tex_fragments/WHITNEY_REAL_CONTINUUM.tex) | Open: no same-model instantiation recorded | 3 / 4 | explicit_non_consumer |
| [OPH-WHITNEY-CHARGED-INSTRUMENT](../paper/tex_fragments/WHITNEY_CHARGED_INSTRUMENT.tex) | Open: no same-model instantiation recorded | 4 / 3 | explicit_non_consumer |
| [OPH-WHITNEY-EPHEMERIS-CLOCK](../paper/tex_fragments/WHITNEY_EPHEMERIS_CLOCK.tex) | Open: no same-model instantiation recorded | 3 / 2 | explicit_non_consumer |
| [OPH-WHITNEY-QUANTUM-HISTORY](../paper/tex_fragments/WHITNEY_QUANTUM_HISTORY.tex) | Open: no same-model instantiation recorded | 4 / 2 | explicit_non_consumer |
| [OPH-WHITNEY-MAGNETIC-CONTINUUM](../paper/tex_fragments/WHITNEY_MAGNETIC_CONTINUUM.tex) | Open: no same-model instantiation recorded | 3 / 5 | explicit_non_consumer |
| [OPH-WHITNEY-QUANTUM-PACKET](../paper/tex_fragments/WHITNEY_QUANTUM_PACKET.tex) | Open: no same-model instantiation recorded | 3 / 5 | explicit_non_consumer |
| [OPH-WHITNEY-CHARGED-ENCLOSURE](../paper/tex_fragments/WHITNEY_CHARGED_ENCLOSURE.tex) | Open: no same-model instantiation recorded | 3 / 6 | explicit_non_consumer |
| [OPH-QFT-RESTRICTION-NATURALITY](../paper/tex_fragments/SOURCE_RESTRICTION_NATURALITY.tex) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-QFT-EFFECTIVE-QUANTUM-COMPARISON](../paper/tex_fragments/EFFECTIVE_QUANTUM_COMPARISON.tex) | Open: no same-model instantiation recorded | 8 / 4 | explicit_non_consumer |
| [OPH-SCALAR-ACTION-CLOCK](../paper/tex_fragments/SOURCE_SCALAR_CLOCK.tex) | Open: no same-model instantiation recorded | 3 / 4 | explicit_non_consumer |
| [OPH-SCALAR-SPACE-TIME-DETECTOR-LIMIT](../paper/tex_fragments/SOURCE_SCALAR_TIME_REFINEMENT.tex) | Open: no same-model instantiation recorded | 2 / 4 | explicit_non_consumer |
| [OPH-COSMO-FLRW-RECORD-DENSITY](../Lean/Geometry/SourceNetConformalRecordDensity.lean) | Open: no same-model instantiation recorded | 3 / 3 | consumed: PR-16; open: PR-15, PR-52 |
| [OPH-BH-CROSSING-READ-AREA-LAW](../Lean/Geometry/CrossingReadAreaLaw.lean) | Open: no same-model instantiation recorded | 2 / 2 | consumed: PR-16; open: PR-15, PR-52 |
| [OPH-COSMO-POSTDICTION-LEDGER](../code/cosmology/postdiction_ledger/build_cosmology_postdiction_ledger.py) | Open: no same-model instantiation recorded | 2 / 2 | explicit_non_consumer |
| [OPH-GOLDEN-SOURCE-QUADRATURE](../Lean/Geometry/GoldenSourceAssignment.lean) | Open: no same-model instantiation recorded | 4 / 5 | explicit_non_consumer |
| [OPH-SCALAR-REGIONAL-COLLAR-RECONSTRUCTION](../Lean/QFT/ScalarRegionalTimeSlice.lean) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-SEAM-PATH-TOMOGRAPHY](../Lean/Geometry/SourceSeamPathTomography.lean) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-CARTAN-SCALAR-CONTINUOUS-READOUT](../paper/tex_fragments/CARTAN_SCALAR_READOUT.tex) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-SOURCE-FEEDBACK-VERSION-TRANSPORT](../Lean/Geometry/SourceFeedbackTransport.lean) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-SCALAR-SEQUENTIAL-INSTRUMENT](../paper/tex_fragments/SOURCE_SCALAR_SEQUENTIAL_INSTRUMENT.tex) | Open: no same-model instantiation recorded | 6 / 5 | explicit_non_consumer |
| [OPH-FINITE-WORK-RESPONSE-NONIDENTIFIABILITY](../paper/tex_fragments/FINITE_EQUATION_OF_STATE.tex) | Open: no same-model instantiation recorded | 1 / 2 | explicit_non_consumer |
| [OPH-GOLDEN-SCALAR-METRIC-PRESSURE](../paper/tex_fragments/FINITE_EQUATION_OF_STATE.tex) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-WHITNEY-MAXWELL-METRIC-PRESSURE](../paper/tex_fragments/FINITE_EQUATION_OF_STATE.tex) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-FERMION-HYPERCHARGE-SOURCE-CURRENT](../paper/tex_fragments/FERMION_SOURCE_CURRENT.tex) | Open: no same-model instantiation recorded | 1 / 2 | open: PR-47, PR-52, PR-54; boundary: PR-15 |
| [OPH-QUANTUM-TRIPLE-OPERATOR-JOIN](../paper/observers_are_all_you_need.tex) | Open: no same-model instantiation recorded | 2 / 3 | consumed: PR-32, PR-34; open: PR-52, PR-58 |
| [OPH-SOURCE-FULL-FAMILY-READ-ROUTING](../Lean/Geometry/SourceReadRouting.lean) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-SOURCE-ROUTING-STORAGE-AND-RAW-COUNT](../Lean/Geometry/SourceRoutingHierarchy.lean) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-NATIVE-LOCAL-SINGLETON-PUBLICATION](../Lean/Geometry/SourceReadAcceptance.lean) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-NATIVE-REUSABLE-ACCUMULATION](../Lean/Geometry/SourceNativeAccumulator.lean) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-NATIVE-CONDITIONAL-STORED-PROGRAM](../paper/tex_fragments/SOURCE_NATIVE_STORED_PROGRAMS.tex) | Open: no same-model instantiation recorded | 1 / 6 | explicit_non_consumer |
| [OPH-FEDERATION-LOG-GLUING-NONIDENTIFICATION](../evidence/federation_recovery/recover.py) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-TEMPORAL-NATIVE-RECOVERY](../paper/tex_fragments/SOURCE_TEMPORAL_RECOVERY.tex) | Open: no same-model instantiation recorded | 1 / 5 | explicit_non_consumer |
| [OPH-NATIVE-CHECKPOINT-SELECTION](../paper/tex_fragments/SOURCE_CHECKPOINT_SELECTION.tex) | Open: no same-model instantiation recorded | 1 / 6 | explicit_non_consumer |
| [OPH-NATIVE-EVENTUAL-PUBLICATION](../paper/tex_fragments/SOURCE_EVENTUAL_PUBLICATION.tex) | Open: no same-model instantiation recorded | 1 / 6 | explicit_non_consumer |
| [OPH-FINITE-A3-SUPPORT](../paper/tex_fragments/SOURCE_EVENTUAL_PUBLICATION.tex) | Open: no same-model instantiation recorded | 1 / 4 | explicit_non_consumer |
| [OPH-REGISTERED-SOURCE-READ-ALGEBRA](../paper/tex_fragments/SOURCE_OPERATION_READS.tex) | Open: no same-model instantiation recorded | 1 / 11 | explicit_non_consumer |
| [OPH-SPIN-EXCHANGE-SELECTION-BOUNDARY](../paper/tex_fragments/SPIN_EXCHANGE_SELECTION_BOUNDARY.tex) | Open: no same-model instantiation recorded | 1 / 3 | boundary: PR-47, PR-54 |
| [OPH-CONDITIONAL-PAULI-SCALAR-Q-SELECTION](../paper/tex_fragments/PAULI_STABILITY_SELECTION.tex) | Open: no same-model instantiation recorded | 1 / 4 | boundary: PR-47, PR-54 |
| [OPH-CONDITIONAL-PAULI-OPERATION-SELECTION](../paper/tex_fragments/PAULI_SOURCE_SELECTION_BOUNDARY.tex) | Open: no same-model instantiation recorded | 1 / 4 | boundary: PR-47, PR-54 |
| [OPH-SCALAR-SEAM-SOURCE-NONENTAILMENT](../paper/tex_fragments/SOURCE_RECORD_GLUING.tex) | Open: no same-model instantiation recorded | 1 / 4 | explicit_non_consumer |
| [OPH-RECORD-GLUING-M1-SUFFICIENCY](../paper/tex_fragments/SOURCE_RECORD_GLUING.tex) | Open: no same-model instantiation recorded | 1 / 5 | explicit_non_consumer |
| [OPH-RG-LOCAL-CLOCK-REDUCTION](../paper/tex_fragments/RECORD_GLUING_PRINCIPLE.tex) | Open: no same-model instantiation recorded | 1 / 4 | boundary: PR-52 |
| [OPH-RG-INFORMATION-RESOURCE-BOUND](../paper/tex_fragments/RECORD_GLUING_PRINCIPLE.tex) | Open: no same-model instantiation recorded | 1 / 2 | explicit_non_consumer |
| [OPH-M1-SPARSE-CAUSAL-COUNT-REPLACEMENT](../code/m1_necessity/DERIVATION.md) | Open: no same-model instantiation recorded | 1 / 4 | explicit_non_consumer |
| [OPH-M1-CUT-OBSERVABLE-NONUNIVERSALITY](../code/m1_necessity/DERIVATION.md) | Open: no same-model instantiation recorded | 2 / 2 | explicit_non_consumer |
| [OPH-M1-SPARSE-ACTION-AND-TICK-OBSTRUCTION](../code/m1_necessity/DERIVATION.md) | Open: no same-model instantiation recorded | 2 / 4 | explicit_non_consumer |
| [OPH-M1-ENTROPY-PRESERVING-SPARSIFICATION](../code/m1_necessity/BALANCED_REPLACEMENT.md) | Open: no same-model instantiation recorded | 6 / 9 | explicit_non_consumer |
| [OPH-M1-MICROSCOPIC-CUT-OBSTRUCTION](../paper/tex_fragments/M1_INTERFACES.tex) | Open: no same-model instantiation recorded | 1 / 2 | explicit_non_consumer |
| [OPH-M1-VARIATIONAL-AREA-REPAIR](../paper/tex_fragments/M1_INTERFACES.tex) | Open: no same-model instantiation recorded | 2 / 4 | explicit_non_consumer |
| [OPH-M1-UNIFORM-ACTION-SPECTRAL-THRESHOLD](../paper/tex_fragments/M1_INTERFACES.tex) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-M1-GLOBAL-VACUUM-FIDELITY](../paper/tex_fragments/M1_VACUUM_FIDELITY.tex) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-M1-SPARSE-QUANTUM-TICK-REPAIR](../paper/tex_fragments/M1_VACUUM_FIDELITY.tex) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-M1-POLYNOMIAL-VACUUM-OBSTRUCTION](../paper/tex_fragments/M1_VACUUM_FIDELITY.tex) | Open: no same-model instantiation recorded | 1 / 2 | explicit_non_consumer |
| [OPH-M1-SPLIT-SIGNAL-CLOCK-OBSTRUCTION](../paper/tex_fragments/M1_VACUUM_FIDELITY.tex) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-M1-QUANTUM-VELOCITY-BUDGET](../paper/tex_fragments/M1_QUANTUM_TRANSPORT.tex) | Open: no same-model instantiation recorded | 1 / 2 | explicit_non_consumer |
| [OPH-M1-CHARGED-WEYL-SPECTRAL-LIMIT](../paper/tex_fragments/M1_QUANTUM_TRANSPORT.tex) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-M1-FINITE-HAMILTONIAN-DELAY-OBSTRUCTION](../paper/tex_fragments/M1_QUANTUM_TRANSPORT.tex) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-M1-COHERENT-VELOCITY-SEPARATION](../paper/tex_fragments/M1_QUANTUM_TRANSPORT.tex) | Open: no same-model instantiation recorded | 1 / 2 | explicit_non_consumer |
| [OPH-M1-FINITE-STENCIL-CLOCK-BOUND](../paper/tex_fragments/M1_QUANTUM_TRANSPORT.tex) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-M1-MINIMAL-TETRAHEDRAL-TRANSPORT](../paper/tex_fragments/M1_QUANTUM_TRANSPORT.tex) | Open: no same-model instantiation recorded | 2 / 5 | explicit_non_consumer |
| [OPH-COSMO-VACUUM-ENSEMBLE-CLASSIFICATION](../cosmology/oph_cosmological_vacuum_and_structure_formation.tex) | Open: no same-model instantiation recorded | 5 / 5 | open: PR-15, PR-52 |
| [OPH-COSMO-BOLTZMANN-TRANSPORT-REDUCTION](../cosmology/oph_boltzmann_transport_derivation.tex) | Open: no same-model instantiation recorded | 8 / 5 | consumed: PR-15, PR-17; open: PR-39, PR-52, PR-53 |
| [OPH-SM-BARYON-DIMENSION-SIX-CENSUS](../Lean/Screen/BaryonDimensionSix.lean) | Open: no same-model instantiation recorded | 4 / 3 | consumed: PR-55, PR-56; open: PR-57 |
| [OPH-NATIVE-REPAIR-HISTORY-SOURCE](../paper/tex_fragments/NATIVE_REPAIR_RECORD_SPECTRUM.tex) | Open: no same-model instantiation recorded | 2 / 1 | explicit_non_consumer |
| [OPH-SETTLEMENT-SCALE-TRANSFER](../paper/tex_fragments/NATIVE_REPAIR_FINITE_EVIDENCE.tex) | Open: no same-model instantiation recorded | 1 / 1 | explicit_non_consumer |
| [OPH-NATIVE-PORT-TRANSPORT-CONTROLS](../paper/tex_fragments/NATIVE_REPAIR_FINITE_EVIDENCE.tex) | Open: no same-model instantiation recorded | 1 / 1 | explicit_non_consumer |
| [OPH-REFINEMENT-ABSOLUTE-INNOVATION](../paper/tex_fragments/NATIVE_REPAIR_FINITE_EVIDENCE.tex) | Open: no same-model instantiation recorded | 1 / 1 | explicit_non_consumer |
| [OPH-LAYERED-OBSERVER-CLOCK-CONTROLS](../paper/tex_fragments/CAUSAL_MANIFOLD_OBSERVATIONS.tex) | Open: no same-model instantiation recorded | 1 / 1 | explicit_non_consumer |
| [OPH-GOLDEN-FIXED-MODE-FOURIER](../Lean/Geometry/GoldenSourceFourier.lean) | Open: no same-model instantiation recorded | 1 / 1 | explicit_non_consumer |
| [OPH-COMMON-BULK-LOG-COVARIANCE](../paper/tex_fragments/RADIAL_LIFT_THEOREMS_330.tex) | Open: no same-model instantiation recorded | 1 / 1 | explicit_non_consumer |
| [OPH-COMMON-HISTORY-FEEDBACK-PACKET](../paper/tex_fragments/COMMON_HISTORY_FEEDBACK_PACKET.tex) | Open: no same-model instantiation recorded | 2 / 2 | explicit_non_consumer |
| [OPH-SOURCE-ACTION-MEASURE-IDENTIFICATION](../paper/tex_fragments/SOURCE_ACTION_DENSITY_IDENTIFICATION.tex) | Open: no same-model instantiation recorded | 3 / 2 | explicit_non_consumer |
| [OPH-M1-MASSIVE-FLIGHT-SPECTRUM](../docs/research/MASSIVE_OPERATIONAL_CLOCKS.md) | Open: no same-model instantiation recorded | 2 / 4 | explicit_non_consumer |
| [OPH-M1-RESOLVED-LOCAL-READ-CONE](../docs/research/MASSIVE_OPERATIONAL_CLOCKS.md) | Open: no same-model instantiation recorded | 2 / 4 | explicit_non_consumer |
| [OPH-M1-OBSERVABLE-MASS-CLOCK](../docs/research/MASSIVE_OPERATIONAL_CLOCKS.md) | Open: no same-model instantiation recorded | 2 / 4 | explicit_non_consumer |
| [OPH-SOURCE-COHERENT-CODE-COMPLETION](../docs/research/COHERENT_SOURCE_CLOCKS.md) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-SOURCE-TWO-BASIS-A3-SELECTION](../docs/research/COHERENT_SOURCE_CLOCKS.md) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-SOURCE-COMPILED-MASS-CLOCK](../docs/research/COHERENT_SOURCE_CLOCKS.md) | Open: no same-model instantiation recorded | 4 / 3 | explicit_non_consumer |
| [OPH-SOURCE-LOCAL-FERMIONIC-REALIZATION](../docs/research/FERMIONIC_SOURCE_CLOCKS.md) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-SOURCE-ENCODED-FERMION-CLOCK](../docs/research/FERMIONIC_SOURCE_CLOCKS.md) | Open: no same-model instantiation recorded | 6 / 3 | explicit_non_consumer |
| [OPH-SOURCE-FERMIONIC-VACUUM-NOISE](../docs/research/FERMIONIC_SOURCE_CLOCKS.md) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-SOURCE-FIXED-RATE-PROTECTED-READS](../code/m1_fixed_noise/CONTRACT.md) | Open: no same-model instantiation recorded | 3 / 4 | explicit_non_consumer |
| [OPH-SOURCE-RAW-DECODED-NOISE-SEPARATION](../code/m1_fixed_noise/CONTRACT.md) | Open: no same-model instantiation recorded | 1 / 3 | explicit_non_consumer |
| [OPH-SOURCE-NOISY-CONTROL-HISTORIES](../code/m1_noisy_records/CONTRACT.md) | Open: no same-model instantiation recorded | 3 / 4 | explicit_non_consumer |
| [OPH-SOURCE-LIVE-RECORD-PROTECTION](../code/m1_noisy_records/CONTRACT.md) | Open: no same-model instantiation recorded | 3 / 4 | explicit_non_consumer |
| [OPH-SOURCE-NATIVE-LEAKAGE-REDUCTION](../code/m1_leakage/__init__.py) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-SOURCE-GENERAL-NOISE-REFINEMENT](../code/m1_leakage/__init__.py) | Open: no same-model instantiation recorded | 3 / 5 | explicit_non_consumer |
| [OPH-SOURCE-EXPANDER-LIVE-ARCHIVE](../code/m1_expander_archive/__init__.py) | Open: no same-model instantiation recorded | 3 / 6 | explicit_non_consumer |
| [OPH-SOURCE-FIXED-STRENGTH-PUBLIC-HISTORY](../code/m1_expander_archive/__init__.py) | Open: no same-model instantiation recorded | 3 / 6 | explicit_non_consumer |
| [OPH-SOURCE-SEA-GENERATOR](../code/m1_sea_energy/__init__.py) | Open: no same-model instantiation recorded | 3 / 4 | explicit_non_consumer |
| [OPH-SOURCE-LOCAL-SEA-LIMIT](../code/m1_sea_energy/__init__.py) | Open: no same-model instantiation recorded | 3 / 4 | explicit_non_consumer |
| [OPH-SOURCE-SEA-CLOCK](../code/m1_sea_energy/__init__.py) | Open: no same-model instantiation recorded | 4 / 7 | explicit_non_consumer |
| [OPH-QUARK-CONDITIONAL-FIXED-P-REPLAY](../code/particles/flavor/conditional_quark_mass_spec.json) | Open: no same-model instantiation recorded | 4 / 3 | explicit_non_consumer |
| [OPH-TOP-CONDITIONAL-FIXED-P-CRITICALITY](../code/particles/calibration/conditional_top_fixed_p_spec.json) | Open: no same-model instantiation recorded | 4 / 4 | explicit_non_consumer |
| [OPH-QUARK-SOURCE-W5-RESPONSE-CONSTRAINTS](../code/particles/flavor/source_w5_response_constraints_spec.json) | Open: no same-model instantiation recorded | 4 / 4 | explicit_non_consumer |
| [OPH-QUARK-NATIVE-REPAIR-FLAVOR-CONSTRAINTS](../code/particles/flavor/native_repair_flavor_constraints_spec.json) | Open: no same-model instantiation recorded | 4 / 4 | explicit_non_consumer |
| [OPH-TOP-POSITIVE-CAPACITY-BAND](../code/particles/calibration/top_positive_capacity_band_spec.json) | Open: no same-model instantiation recorded | 4 / 4 | explicit_non_consumer |
| [OPH-QUARK-CONDITIONAL-RECORD-RENEWAL](../code/particles/flavor/native_record_renewal_spec.json) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-QUARK-EXTERIOR-HIGGS-RESPONSE](../code/particles/flavor/exterior_higgs_response_spec.json) | Open: no same-model instantiation recorded | 2 / 4 | explicit_non_consumer |
| [OPH-EW-CONDITIONAL-SOURCE-VEV-MATCHING](../code/particles/calibration/source_ew_vev_matching_spec.json) | Open: no same-model instantiation recorded | 4 / 4 | explicit_non_consumer |
| [OPH-SOURCE-PAIRED-GRAVITY-NONIDENTIFICATION](../code/paired_gravity/CONTRACT.md) | Open: no same-model instantiation recorded | 1 / 6 | explicit_non_consumer |
| [OPH-MCKAY-INDEPENDENT-GOLDEN-FIELD-CONTROL](../code/a5_closure/mckay_golden_field_certificate.py) | Open: no same-model instantiation recorded | 1 / 4 | explicit_non_consumer |
| [OPH-DARK-SOURCE-DYNAMICS-LENSING](../docs/research/DARK_SOURCE_DYNAMICS_LENSING.md) | Open: no same-model instantiation recorded | 1 / 5 | explicit_non_consumer |
| [OPH-SOURCE-ADMISSIBLE-CLOCK-SPATIAL-NONSELECTION](../docs/research/SOURCE_GRAVITY_ADMISSIBILITY.md) | Open: no same-model instantiation recorded | 2 / 4 | explicit_non_consumer |
| [OPH-INDEPENDENT-ALPHA-TAU-REPRODUCTION](../docs/research/INDEPENDENT_POSTDICTION_COMPARISON.md) | Open: no same-model instantiation recorded | 2 / 4 | consumed: PR-13, PR-17, PR-37; boundary: PR-14, PR-38, PR-39, PR-50 |
| [OPH-FIXED-DISPERSION-FEASIBILITY](../docs/research/DISPERSION_FEASIBILITY.md) | Open: no same-model instantiation recorded | 2 / 3 | explicit_non_consumer |
| [OPH-MODEL-NATIVE-OBSERVER-EXECUTION](../model/observer/native.py) | Open: no same-model instantiation recorded | 9 / 9 | explicit_non_consumer |
| [OPH-MODEL-OBSERVER-RECORD-CAPACITY](../model/observer/capacity.py) | Open: no same-model instantiation recorded | 3 / 3 | explicit_non_consumer |
| [OPH-MODEL-P-RESPONSE-EXECUTION](../model/observer/p_response/experiment.py) | Open: no same-model instantiation recorded | 6 / 4 | consumed: PR-17, PR-37; open: PR-38, PR-39 |
| [OPH-MODEL-N4-LOCAL-TICK-CERTIFICATE](../model/local_tick/experiment.py) | Open: no same-model instantiation recorded | 3 / 4 | explicit_non_consumer |
| [OPH-MODEL-SCALAR-RESPONSE-EXECUTION](../model/observer/scalar_response/README.md) | Open: no same-model instantiation recorded | 7 / 8 | explicit_non_consumer |
| [OPH-MODEL-COMMON-SOURCE-EXECUTION](../model/common_source/README.md) | Open: no same-model instantiation recorded | 14 / 24 | explicit_non_consumer |
| [OPH-MODEL-SPATIAL-INTERACTION-EXECUTION](../model/spatial_interaction/README.md) | Open: no same-model instantiation recorded | 2 / 4 | explicit_non_consumer |
| [OPH-MODEL-CANONICAL-SCALAR-EXECUTION](../model/relativistic_scalar/README.md) | Open: no same-model instantiation recorded | 2 / 4 | explicit_non_consumer |
| [OPH-MODEL-SHARED-SCALAR-OBSERVATIONS](../model/scalar_observations/README.md) | supplied-first-band-control | 3 / 5 | explicit_non_consumer |
| [OPH-MODEL-PUBLISHED-OSCILLATOR-READOUT](../model/scalar_observations/experimental_readout/README.md) | Open: no same-model instantiation recorded | 3 / 7 | explicit_non_consumer |
| [OPH-MODEL-COUPLED-INTERACTION-RECOVERY](../model/coupled_reconstruction/README.md) | Open: no same-model instantiation recorded | 2 / 5 | explicit_non_consumer |
| [OPH-MODEL-ANHARMONIC-TRANSFER-COMPARISON](../model/anharmonic_transfer/README.md) | Open: no same-model instantiation recorded | 2 / 5 | explicit_non_consumer |
| [OPH-MODEL-NATIVE-PROCESS-RECOVERY](../model/OBSERVER_RELATIVE_RECOVERY.md) | native-coherent-instrument, complete-current-identification, compiled-incident-history, calibrated-complete-incident-law, native-occupation-fluctuation-field, ideal-resource-clock-attachment | 6 / 13 | explicit_non_consumer |

### Every Lean source module

The declaration count uses the existing theorem-count gate's public declaration
grammar with comments removed. Evidence membership is file-level navigation,
not an assertion that a registry claim consumes every declaration in that file.
Modules without a registry evidence edge remain visible. Every declaration in
a module retains its original quantified hypotheses and proof in the linked source.

| Lean source | Public theorem/lemma declarations | Registered evidence users |
| --- | ---: | --- |
| [Lean/Computation/FixedFederationComplexity.lean](../Lean/Computation/FixedFederationComplexity.lean) | 26 | OPH-COMPUTATION-FIXED-FEDERATION-EXECUTION-CLASSIFICATION |
| [Lean/Computation/FixedFederationCounterexamples.lean](../Lean/Computation/FixedFederationCounterexamples.lean) | 27 | OPH-COMPUTATION-FIXED-FEDERATION-PROGRESS |
| [Lean/Computation/FixedFederationExecution.lean](../Lean/Computation/FixedFederationExecution.lean) | 20 | OPH-COMPUTATION-FIXED-FEDERATION-EXECUTION-CLASSIFICATION |
| [Lean/Computation/FixedFederationExecutionExamples.lean](../Lean/Computation/FixedFederationExecutionExamples.lean) | 54 | OPH-COMPUTATION-FIXED-FEDERATION-EXECUTION-CLASSIFICATION |
| [Lean/Computation/FixedFederationFanoutControl.lean](../Lean/Computation/FixedFederationFanoutControl.lean) | 29 | OPH-COMPUTATION-FIXED-FEDERATION-EXECUTION-CLASSIFICATION |
| [Lean/Computation/FixedFederationProgress.lean](../Lean/Computation/FixedFederationProgress.lean) | 36 | OPH-COMPUTATION-FIXED-FEDERATION-PROGRESS, OPH-COMPUTATION-FIXED-FEDERATION-ENDPOINT |
| [Lean/Computation/RepairUniversality.lean](../Lean/Computation/RepairUniversality.lean) | 43 | OPH-COMPUTATION-REPAIR-PROGRAM-BOUNDARY |
| [Lean/Dynamics/CenterSpectral.lean](../Lean/Dynamics/CenterSpectral.lean) | 38 | OPH-PUBLIC-PRIVATE-DYNAMICS-FINITE |
| [Lean/Dynamics/CentralBlocks.lean](../Lean/Dynamics/CentralBlocks.lean) | 24 | OPH-PUBLIC-PRIVATE-DYNAMICS-FINITE |
| [Lean/Dynamics/ChoiCPTP.lean](../Lean/Dynamics/ChoiCPTP.lean) | 38 | OPH-QUANTUM-INSTRUMENT-SELECTION-BY-CERTAIN-STATES |
| [Lean/Dynamics/ConditionalExpectationGenerator.lean](../Lean/Dynamics/ConditionalExpectationGenerator.lean) | 6 | OPH-PUBLICIZATION-DYNAMICS-FINITE |
| [Lean/Dynamics/PrivateInner.lean](../Lean/Dynamics/PrivateInner.lean) | 6 | OPH-PUBLIC-PRIVATE-DYNAMICS-FINITE |
| [Lean/Dynamics/ProtectedCharge.lean](../Lean/Dynamics/ProtectedCharge.lean) | 8 | OPH-FINITE-CONSERVATION-WARD-PRECURSOR, OPH-MODEL-COMMON-SOURCE-EXECUTION |
| [Lean/Dynamics/PublicAutomorphism.lean](../Lean/Dynamics/PublicAutomorphism.lean) | 12 | OPH-PUBLIC-PRIVATE-DYNAMICS-FINITE |
| [Lean/Dynamics/PublicMarkov.lean](../Lean/Dynamics/PublicMarkov.lean) | 18 | OPH-PUBLIC-PRIVATE-DYNAMICS-FINITE |
| [Lean/Dynamics/StarWedderburn.lean](../Lean/Dynamics/StarWedderburn.lean) | 54 | OPH-PUBLIC-PRIVATE-DYNAMICS-FINITE |
| [Lean/Dynamics/StoneConverse.lean](../Lean/Dynamics/StoneConverse.lean) | 11 | OPH-PUBLIC-PRIVATE-DYNAMICS-FINITE |
| [Lean/Dynamics/WardLimitManifest.lean](../Lean/Dynamics/WardLimitManifest.lean) | 1 | OPH-FINITE-CONSERVATION-WARD-PRECURSOR |
| [Lean/Dynamics.lean](../Lean/Dynamics.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/Basic.lean](../Lean/EventAlgebra/Basic.lean) | 25 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/ExpectationBound.lean](../Lean/EventAlgebra/ExpectationBound.lean) | 4 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/FiniteBornFrame.lean](../Lean/EventAlgebra/FiniteBornFrame.lean) | 18 | OPH-FINITE-BORN-FRAME-RANK-GAP |
| [Lean/EventAlgebra/FiniteBuschGleason.lean](../Lean/EventAlgebra/FiniteBuschGleason.lean) | 69 | OPH-FINITE-BORN-FRAME-RANK-GAP |
| [Lean/EventAlgebra/FiniteEffectClosureBoundary.lean](../Lean/EventAlgebra/FiniteEffectClosureBoundary.lean) | 15 | OPH-FINITE-BORN-FRAME-RANK-GAP |
| [Lean/EventAlgebra/FiniteWebBornNoGo.lean](../Lean/EventAlgebra/FiniteWebBornNoGo.lean) | 34 | OPH-FINITE-BORN-FRAME-RANK-GAP |
| [Lean/EventAlgebra/FrequencyConcentration.lean](../Lean/EventAlgebra/FrequencyConcentration.lean) | 38 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean](../Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean) | 44 | OPH-QUANTUM-INSTRUMENT-SELECTION-BY-CERTAIN-STATES |
| [Lean/EventAlgebra/InterlockingContexts.lean](../Lean/EventAlgebra/InterlockingContexts.lean) | 29 | OPH-FINITE-BORN-FRAME-RANK-GAP |
| [Lean/EventAlgebra/Lueders.lean](../Lean/EventAlgebra/Lueders.lean) | 13 | OPH-QUANTUM-INSTRUMENT-SELECTION-BY-CERTAIN-STATES |
| [Lean/EventAlgebra/LuedersPhaseInstrument.lean](../Lean/EventAlgebra/LuedersPhaseInstrument.lean) | 68 | OPH-QUANTUM-LUEDERS-PHASE-INSTRUMENT, OPH-EVENTALGEBRA-SOURCE-PHASE-INSTRUMENT-OUTCOME-BRIDGE, OPH-QUANTUM-SOURCE-BOUND-INSTRUMENT-INTERFACE, OPH-QUANTUM-INSTRUMENT-SELECTION-BY-CERTAIN-STATES |
| [Lean/EventAlgebra/NoBroadcastingAdapter.lean](../Lean/EventAlgebra/NoBroadcastingAdapter.lean) | 7 | OPH-PUBLIC-RECORD-ALGEBRA-NOCLONING |
| [Lean/EventAlgebra/OperationalAdditivityBoundary.lean](../Lean/EventAlgebra/OperationalAdditivityBoundary.lean) | 41 | OPH-QUANTUM-OPERATIONAL-ADDITIVITY-BOUNDARY, OPH-QUANTUM-SOURCE-BOUND-INSTRUMENT-INTERFACE |
| [Lean/EventAlgebra/OperationalPhaseAttainment.lean](../Lean/EventAlgebra/OperationalPhaseAttainment.lean) | 19 | OPH-QUANTUM-PHASE-OPERATION-RECORDED-DECISION-ATTAINMENT |
| [Lean/EventAlgebra/OperationalPhaseInstrument.lean](../Lean/EventAlgebra/OperationalPhaseInstrument.lean) | 43 | OPH-QUANTUM-OPERATIONAL-PHASE-INTERFACE |
| [Lean/EventAlgebra/PartitionAverage.lean](../Lean/EventAlgebra/PartitionAverage.lean) | 30 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/PartitionAverageCP.lean](../Lean/EventAlgebra/PartitionAverageCP.lean) | 10 | OPH-PUBLICIZATION-DYNAMICS-FINITE |
| [Lean/EventAlgebra/PartitionPinching.lean](../Lean/EventAlgebra/PartitionPinching.lean) | 27 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/PartitionPinchingCP.lean](../Lean/EventAlgebra/PartitionPinchingCP.lean) | 2 | OPH-THERMO-FOUR-LAW-PACKAGE |
| [Lean/EventAlgebra/PauliMinimalCoupling.lean](../Lean/EventAlgebra/PauliMinimalCoupling.lean) | 14 | OPH-EVENTALGEBRA-PAULI-MINIMAL-COUPLING-SQUARE |
| [Lean/EventAlgebra/PhaseInstrumentDetermination.lean](../Lean/EventAlgebra/PhaseInstrumentDetermination.lean) | 24 | OPH-QUANTUM-PHASE-MODEL-DETERMINATION-BOUNDARY |
| [Lean/EventAlgebra/ProductSplitSeparability.lean](../Lean/EventAlgebra/ProductSplitSeparability.lean) | 61 | OPH-EVENTALGEBRA-SOURCE-REACHABILITY-DIAGONAL-CHSH, OPH-EVENTALGEBRA-PRODUCT-SPLIT-SLOT-LOCAL-CHSH |
| [Lean/EventAlgebra/PublicRecordAlgebra.lean](../Lean/EventAlgebra/PublicRecordAlgebra.lean) | 19 | OPH-PUBLIC-RECORD-ALGEBRA-NOCLONING |
| [Lean/EventAlgebra/QuantumAdequacySurface.lean](../Lean/EventAlgebra/QuantumAdequacySurface.lean) | 21 | OPH-EM-CARRIER-EVOLUTION-FLOW |
| [Lean/EventAlgebra/RecordMajorization.lean](../Lean/EventAlgebra/RecordMajorization.lean) | 21 | OPH-FINITE-QUANTUM-LIMITATION-SUITE |
| [Lean/EventAlgebra/Robertson.lean](../Lean/EventAlgebra/Robertson.lean) | 24 | OPH-FINITE-QUANTUM-LIMITATION-SUITE |
| [Lean/EventAlgebra/SchroedingerFrameFlow.lean](../Lean/EventAlgebra/SchroedingerFrameFlow.lean) | 5 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/SlotLocalitySurface.lean](../Lean/EventAlgebra/SlotLocalitySurface.lean) | 2 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/SourceBoundInstrumentInterface.lean](../Lean/EventAlgebra/SourceBoundInstrumentInterface.lean) | 33 | OPH-QUANTUM-SOURCE-BOUND-INSTRUMENT-INTERFACE, OPH-QUANTUM-INS03-EXPORT-VALIDATOR |
| [Lean/EventAlgebra/SourcePhaseBornWeightBoundary.lean](../Lean/EventAlgebra/SourcePhaseBornWeightBoundary.lean) | 20 | OPH-EVENTALGEBRA-SOURCE-PHASE-BORN-WEIGHT-BOUNDARY |
| [Lean/EventAlgebra/SourcePhaseCommonPreparationHull.lean](../Lean/EventAlgebra/SourcePhaseCommonPreparationHull.lean) | 13 | OPH-EVENTALGEBRA-SOURCE-PHASE-COMMON-PREPARATION-HULL |
| [Lean/EventAlgebra/SourcePhaseInstrumentOutcomeBridge.lean](../Lean/EventAlgebra/SourcePhaseInstrumentOutcomeBridge.lean) | 8 | OPH-EVENTALGEBRA-SOURCE-PHASE-INSTRUMENT-OUTCOME-BRIDGE, OPH-EVENTALGEBRA-SOURCE-PHASE-COMMON-PREPARATION-HULL |
| [Lean/EventAlgebra/SourcePhaseSelection.lean](../Lean/EventAlgebra/SourcePhaseSelection.lean) | 17 | OPH-EVENTALGEBRA-SOURCE-PHASE-EFFECT-SELECTION, OPH-EVENTALGEBRA-SOURCE-PHASE-INSTRUMENT-OUTCOME-BRIDGE, OPH-EVENTALGEBRA-SOURCE-PHASE-COMMON-PREPARATION-HULL, OPH-EVENTALGEBRA-SOURCE-PHASE-BORN-WEIGHT-BOUNDARY |
| [Lean/EventAlgebra/SourcePhaseSelectionInput.lean](../Lean/EventAlgebra/SourcePhaseSelectionInput.lean) | 3 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/SourcePhaseSelectionRelation.lean](../Lean/EventAlgebra/SourcePhaseSelectionRelation.lean) | 7 | OPH-EVENTALGEBRA-SOURCE-PHASE-EFFECT-SELECTION |
| [Lean/EventAlgebra/SourcePhaseSelectionSemantics.lean](../Lean/EventAlgebra/SourcePhaseSelectionSemantics.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/SourceReachabilityDelimitation.lean](../Lean/EventAlgebra/SourceReachabilityDelimitation.lean) | 34 | OPH-EVENTALGEBRA-SOURCE-REACHABILITY-DIAGONAL-CHSH, OPH-EVENTALGEBRA-SOURCE-PHASE-BORN-WEIGHT-BOUNDARY, OPH-EVENTALGEBRA-PRODUCT-SPLIT-SLOT-LOCAL-CHSH |
| [Lean/EventAlgebra/SpectralEntropyBoundary.lean](../Lean/EventAlgebra/SpectralEntropyBoundary.lean) | 11 | OPH-FINITE-QUANTUM-LIMITATION-SUITE |
| [Lean/EventAlgebra/StateExpectation.lean](../Lean/EventAlgebra/StateExpectation.lean) | 4 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/Superselection.lean](../Lean/EventAlgebra/Superselection.lean) | 9 | OPH-FINITE-QUANTUM-LIMITATION-SUITE |
| [Lean/EventAlgebra/Tsirelson.lean](../Lean/EventAlgebra/Tsirelson.lean) | 8 | No claim evidence edge; model mapping open |
| [Lean/EventAlgebra/TsirelsonSaturation.lean](../Lean/EventAlgebra/TsirelsonSaturation.lean) | 39 | OPH-QUANTUM-TSIRELSON-SATURATION, OPH-EVENTALGEBRA-SOURCE-REACHABILITY-DIAGONAL-CHSH, OPH-EVENTALGEBRA-PRODUCT-SPLIT-SLOT-LOCAL-CHSH |
| [Lean/EventAlgebra/TwoScalePublicRepair.lean](../Lean/EventAlgebra/TwoScalePublicRepair.lean) | 12 | OPH-PUBLICIZATION-DYNAMICS-FINITE |
| [Lean/EventAlgebra.lean](../Lean/EventAlgebra.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Geometry/CanonicalLorentzModule.lean](../Lean/Geometry/CanonicalLorentzModule.lean) | 18 | OPH-CANONICAL-LORENTZ-MODULE, OPH-QFT-CAUSAL-ORDER-COMPOSITION |
| [Lean/Geometry/CarrierDynamicsCompatibility.lean](../Lean/Geometry/CarrierDynamicsCompatibility.lean) | 75 | OPH-GEOMETRY-CARRIER-DYNAMICS-COMPATIBILITY, OPH-GEOMETRY-CHARGE-FIXED-INTERACTION, OPH-EM-CARRIER-MODE-EQUIVARIANCE |
| [Lean/Geometry/CausalOrderComposition.lean](../Lean/Geometry/CausalOrderComposition.lean) | 37 | OPH-QFT-CAUSAL-ORDER-COMPOSITION, OPH-GEOMETRY-EVENT-POPULATION-CHART-INTERFACE |
| [Lean/Geometry/CelestialNullCone.lean](../Lean/Geometry/CelestialNullCone.lean) | 8 | OPH-CANONICAL-LORENTZ-MODULE |
| [Lean/Geometry/CelestialSoldering.lean](../Lean/Geometry/CelestialSoldering.lean) | 7 | OPH-ALGEBRAIC-EVENT-FRAME-SOLDERING |
| [Lean/Geometry/ChargeFixedInteraction.lean](../Lean/Geometry/ChargeFixedInteraction.lean) | 69 | OPH-GEOMETRY-CHARGE-FIXED-INTERACTION, OPH-GEOMETRY-PORT-CHARGE-MINIMAL-COUPLING, OPH-GEOMETRY-WORLDLINE-HOP-TRANSPORT |
| [Lean/Geometry/CommonWorldInstrumentJoin.lean](../Lean/Geometry/CommonWorldInstrumentJoin.lean) | 56 | OPH-GEOMETRY-COMMON-WORLD-INSTRUMENT-JOIN, OPH-GEOMETRY-COMMON-WORLD-MAXWELL-CLOCK-JOIN |
| [Lean/Geometry/CommonWorldIslandBridge.lean](../Lean/Geometry/CommonWorldIslandBridge.lean) | 18 | OPH-GEOMETRY-COMMON-WORLD-ISLAND-BRIDGE |
| [Lean/Geometry/CommonWorldJointAction.lean](../Lean/Geometry/CommonWorldJointAction.lean) | 45 | OPH-GEOMETRY-COMMON-WORLD-JOINT-ACTION, OPH-GEOMETRY-SCREEN-CARRIER-MAP-CANDIDATE, OPH-GEOMETRY-CHARGE-FIXED-INTERACTION |
| [Lean/Geometry/CommonWorldKinematicsWitness.lean](../Lean/Geometry/CommonWorldKinematicsWitness.lean) | 21 | OPH-GEOMETRY-COMMON-WORLD-TWO-ISLAND-WITNESS |
| [Lean/Geometry/CommonWorldMaxwellClockJoin.lean](../Lean/Geometry/CommonWorldMaxwellClockJoin.lean) | 31 | OPH-GEOMETRY-COMMON-WORLD-MAXWELL-CLOCK-JOIN, OPH-GEOMETRY-COMMON-WORLD-JOINT-ACTION |
| [Lean/Geometry/CompositeMomentumCovariance.lean](../Lean/Geometry/CompositeMomentumCovariance.lean) | 44 | OPH-GEOMETRY-COMPOSITE-MOMENTUM-COVARIANCE |
| [Lean/Geometry/CrossingReadAreaLaw.lean](../Lean/Geometry/CrossingReadAreaLaw.lean) | 8 | OPH-BH-CROSSING-READ-AREA-LAW |
| [Lean/Geometry/EinsteinBranchInhabitant.lean](../Lean/Geometry/EinsteinBranchInhabitant.lean) | 37 | OPH-GR-EINSTEIN-BRANCH-INHABITANT |
| [Lean/Geometry/EinsteinTensorBridge.lean](../Lean/Geometry/EinsteinTensorBridge.lean) | 4 | OPH-GR-SOURCE-DIRECTION-EINSTEIN-COMPOSITION, OPH-CANONICAL-LORENTZ-MODULE |
| [Lean/Geometry/EventFrameSoldering.lean](../Lean/Geometry/EventFrameSoldering.lean) | 4 | OPH-ALGEBRAIC-EVENT-FRAME-SOLDERING |
| [Lean/Geometry/EventGermDisplacement.lean](../Lean/Geometry/EventGermDisplacement.lean) | 10 | OPH-ALGEBRAIC-EVENT-FRAME-SOLDERING |
| [Lean/Geometry/EventPopulationChartInterface.lean](../Lean/Geometry/EventPopulationChartInterface.lean) | 9 | OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE, OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM, OPH-GEOMETRY-EVENT-POPULATION-CHART-INTERFACE |
| [Lean/Geometry/FiniteLayerQuadrature.lean](../Lean/Geometry/FiniteLayerQuadrature.lean) | 6 | OPH-LAYERED-OBSERVER-CLOCK-CONTROLS |
| [Lean/Geometry/FlatDiamondError.lean](../Lean/Geometry/FlatDiamondError.lean) | 6 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/FlatDiamondNormalization.lean](../Lean/Geometry/FlatDiamondNormalization.lean) | 15 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM, OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/FlatDiamondPairIntegral.lean](../Lean/Geometry/FlatDiamondPairIntegral.lean) | 12 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/FlatDiamondVolume.lean](../Lean/Geometry/FlatDiamondVolume.lean) | 18 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/FlatLorentzVolume.lean](../Lean/Geometry/FlatLorentzVolume.lean) | 23 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/FreeEvolutionPersistence.lean](../Lean/Geometry/FreeEvolutionPersistence.lean) | 27 | OPH-GEOMETRY-FREE-EVOLUTION-PERSISTENCE |
| [Lean/Geometry/GoldenSourceAssignment.lean](../Lean/Geometry/GoldenSourceAssignment.lean) | 27 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM, OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/GoldenSourceCausalLimit.lean](../Lean/Geometry/GoldenSourceCausalLimit.lean) | 13 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/GoldenSourceCountLimit.lean](../Lean/Geometry/GoldenSourceCountLimit.lean) | 29 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/GoldenSourceFourier.lean](../Lean/Geometry/GoldenSourceFourier.lean) | 6 | OPH-GOLDEN-FIXED-MODE-FOURIER |
| [Lean/Geometry/GoldenSourcePairLimit.lean](../Lean/Geometry/GoldenSourcePairLimit.lean) | 9 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/GoldenSourceVolumeLimit.lean](../Lean/Geometry/GoldenSourceVolumeLimit.lean) | 6 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM, OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/IntegerKCombInvariance.lean](../Lean/Geometry/IntegerKCombInvariance.lean) | 30 | OPH-GW-INTEGER-K-COMB-TEMPLATE |
| [Lean/Geometry/InternalClockRestFrequency.lean](../Lean/Geometry/InternalClockRestFrequency.lean) | 30 | OPH-GEOMETRY-INTERNAL-CLOCK-REST-FREQUENCY, OPH-GEOMETRY-PROPER-TIME-INTERNAL-ACTION |
| [Lean/Geometry/InternalEnergyInertia.lean](../Lean/Geometry/InternalEnergyInertia.lean) | 43 | OPH-GEOMETRY-INTERNAL-ENERGY-INERTIA-PRECURSOR, OPH-GEOMETRY-COMPOSITE-MOMENTUM-COVARIANCE, OPH-GEOMETRY-PROPER-TIME-INTERNAL-ACTION, OPH-THERMO-MODULAR-ENERGY-ADDITIVITY |
| [Lean/Geometry/InverseSquareShellLaw.lean](../Lean/Geometry/InverseSquareShellLaw.lean) | 8 | OPH-GEOMETRY-INVERSE-SQUARE-SHELL, OPH-GEOMETRY-REST-FIBER-SHELL-TRANSPORT |
| [Lean/Geometry/LorentzOverlapCocycle.lean](../Lean/Geometry/LorentzOverlapCocycle.lean) | 9 | OPH-ALGEBRAIC-EVENT-FRAME-SOLDERING |
| [Lean/Geometry/M1Interfaces.lean](../Lean/Geometry/M1Interfaces.lean) | 8 | OPH-M1-MICROSCOPIC-CUT-OBSTRUCTION, OPH-M1-VARIATIONAL-AREA-REPAIR, OPH-M1-UNIFORM-ACTION-SPECTRAL-THRESHOLD |
| [Lean/Geometry/M1InterfacesAxiomAudit.lean](../Lean/Geometry/M1InterfacesAxiomAudit.lean) | 0 | OPH-M1-MICROSCOPIC-CUT-OBSTRUCTION, OPH-M1-VARIATIONAL-AREA-REPAIR, OPH-M1-UNIFORM-ACTION-SPECTRAL-THRESHOLD |
| [Lean/Geometry/M1Necessity.lean](../Lean/Geometry/M1Necessity.lean) | 17 | OPH-M1-SPARSE-CAUSAL-COUNT-REPLACEMENT, OPH-M1-CUT-OBSERVABLE-NONUNIVERSALITY, OPH-M1-SPARSE-ACTION-AND-TICK-OBSTRUCTION, OPH-M1-ENTROPY-PRESERVING-SPARSIFICATION |
| [Lean/Geometry/M1NecessityAxiomAudit.lean](../Lean/Geometry/M1NecessityAxiomAudit.lean) | 0 | OPH-M1-SPARSE-CAUSAL-COUNT-REPLACEMENT, OPH-M1-CUT-OBSERVABLE-NONUNIVERSALITY, OPH-M1-SPARSE-ACTION-AND-TICK-OBSTRUCTION, OPH-M1-ENTROPY-PRESERVING-SPARSIFICATION |
| [Lean/Geometry/M1QuantumTransport.lean](../Lean/Geometry/M1QuantumTransport.lean) | 13 | OPH-M1-QUANTUM-VELOCITY-BUDGET, OPH-M1-CHARGED-WEYL-SPECTRAL-LIMIT, OPH-M1-FINITE-HAMILTONIAN-DELAY-OBSTRUCTION, OPH-M1-COHERENT-VELOCITY-SEPARATION, OPH-M1-FINITE-STENCIL-CLOCK-BOUND, OPH-M1-MINIMAL-TETRAHEDRAL-TRANSPORT |
| [Lean/Geometry/M1QuantumTransportAxiomAudit.lean](../Lean/Geometry/M1QuantumTransportAxiomAudit.lean) | 0 | OPH-M1-QUANTUM-VELOCITY-BUDGET, OPH-M1-CHARGED-WEYL-SPECTRAL-LIMIT, OPH-M1-FINITE-HAMILTONIAN-DELAY-OBSTRUCTION, OPH-M1-COHERENT-VELOCITY-SEPARATION, OPH-M1-FINITE-STENCIL-CLOCK-BOUND, OPH-M1-MINIMAL-TETRAHEDRAL-TRANSPORT |
| [Lean/Geometry/M1VacuumFidelity.lean](../Lean/Geometry/M1VacuumFidelity.lean) | 10 | OPH-M1-GLOBAL-VACUUM-FIDELITY, OPH-M1-SPARSE-QUANTUM-TICK-REPAIR, OPH-M1-POLYNOMIAL-VACUUM-OBSTRUCTION, OPH-M1-SPLIT-SIGNAL-CLOCK-OBSTRUCTION |
| [Lean/Geometry/M1VacuumFidelityAxiomAudit.lean](../Lean/Geometry/M1VacuumFidelityAxiomAudit.lean) | 0 | OPH-M1-GLOBAL-VACUUM-FIDELITY, OPH-M1-SPARSE-QUANTUM-TICK-REPAIR, OPH-M1-POLYNOMIAL-VACUUM-OBSTRUCTION, OPH-M1-SPLIT-SIGNAL-CLOCK-OBSTRUCTION |
| [Lean/Geometry/MassShellKinematics.lean](../Lean/Geometry/MassShellKinematics.lean) | 15 | OPH-GEOMETRY-KINEMATIC-MASS-SHELL, OPH-GEOMETRY-INTERNAL-CLOCK-REST-FREQUENCY, OPH-GEOMETRY-INTERNAL-ENERGY-INERTIA-PRECURSOR, OPH-GEOMETRY-COMPOSITE-MOMENTUM-COVARIANCE |
| [Lean/Geometry/MetricKernelEnergy.lean](../Lean/Geometry/MetricKernelEnergy.lean) | 5 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM |
| [Lean/Geometry/NativeGeometricSourceIdentification.lean](../Lean/Geometry/NativeGeometricSourceIdentification.lean) | 12 | OPH-NATIVE-REPAIR-HISTORY-SOURCE |
| [Lean/Geometry/ObserverFrameHyperboloid.lean](../Lean/Geometry/ObserverFrameHyperboloid.lean) | 2 | OPH-CANONICAL-LORENTZ-MODULE |
| [Lean/Geometry/ObserverRestSpace.lean](../Lean/Geometry/ObserverRestSpace.lean) | 10 | OPH-CANONICAL-LORENTZ-MODULE |
| [Lean/Geometry/OrderingFractionFourDimensional.lean](../Lean/Geometry/OrderingFractionFourDimensional.lean) | 25 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM |
| [Lean/Geometry/PortChargeMinimalCoupling.lean](../Lean/Geometry/PortChargeMinimalCoupling.lean) | 42 | OPH-GEOMETRY-PORT-CHARGE-MINIMAL-COUPLING, OPH-GEOMETRY-WORLDLINE-HOP-TRANSPORT, OPH-GEOMETRY-TRANSPORTED-CHARGE-FORCE-LAW, OPH-GEOMETRY-TIMELIKE-CLASS-FORCE-LAW, OPH-EM-SOURCED-FIELD-CONSISTENCY |
| [Lean/Geometry/ProperLengthClockedChain.lean](../Lean/Geometry/ProperLengthClockedChain.lean) | 55 | OPH-GEOMETRY-PROPER-LENGTH-CLOCKED-CHAIN |
| [Lean/Geometry/ProperTimeInternalAction.lean](../Lean/Geometry/ProperTimeInternalAction.lean) | 49 | OPH-GEOMETRY-PROPER-TIME-INTERNAL-ACTION, OPH-GEOMETRY-SOURCE-CLOCK-RATE-ALONG-WORLDLINES |
| [Lean/Geometry/RecordGluingPrinciple.lean](../Lean/Geometry/RecordGluingPrinciple.lean) | 12 | OPH-RG-LOCAL-CLOCK-REDUCTION, OPH-RG-INFORMATION-RESOURCE-BOUND |
| [Lean/Geometry/RecordGluingPrincipleAxiomAudit.lean](../Lean/Geometry/RecordGluingPrincipleAxiomAudit.lean) | 0 | OPH-RG-LOCAL-CLOCK-REDUCTION, OPH-RG-INFORMATION-RESOURCE-BOUND |
| [Lean/Geometry/RefiningLatticeCausalCone.lean](../Lean/Geometry/RefiningLatticeCausalCone.lean) | 16 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM |
| [Lean/Geometry/RestFiberShellTransport.lean](../Lean/Geometry/RestFiberShellTransport.lean) | 27 | OPH-GEOMETRY-REST-FIBER-SHELL-TRANSPORT |
| [Lean/Geometry/ScreenCarrierMapCandidate.lean](../Lean/Geometry/ScreenCarrierMapCandidate.lean) | 77 | OPH-GEOMETRY-SCREEN-CARRIER-MAP-CANDIDATE, OPH-GEOMETRY-CARRIER-DYNAMICS-COMPATIBILITY, OPH-GEOMETRY-WORLDLINE-HOP-TRANSPORT |
| [Lean/Geometry/SeamStepSpeedLimit.lean](../Lean/Geometry/SeamStepSpeedLimit.lean) | 44 | OPH-GEOMETRY-SEAM-STEP-SPEED-LIMIT, OPH-GEOMETRY-TIMELIKE-CLASS-FORCE-LAW |
| [Lean/Geometry/SourceAccumulatorAxiomAudit.lean](../Lean/Geometry/SourceAccumulatorAxiomAudit.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceAccumulatorDecoder.lean](../Lean/Geometry/SourceAccumulatorDecoder.lean) | 7 | OPH-NATIVE-REUSABLE-ACCUMULATION |
| [Lean/Geometry/SourceAccumulatorProgram.lean](../Lean/Geometry/SourceAccumulatorProgram.lean) | 23 | OPH-NATIVE-REUSABLE-ACCUMULATION |
| [Lean/Geometry/SourceActionMeasure.lean](../Lean/Geometry/SourceActionMeasure.lean) | 15 | OPH-SOURCE-ACTION-MEASURE-IDENTIFICATION |
| [Lean/Geometry/SourceBankBusWitness.lean](../Lean/Geometry/SourceBankBusWitness.lean) | 8 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceBankCompiler.lean](../Lean/Geometry/SourceBankCompiler.lean) | 8 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceBankControl.lean](../Lean/Geometry/SourceBankControl.lean) | 9 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceBankExecution.lean](../Lean/Geometry/SourceBankExecution.lean) | 1 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceBankInvariant.lean](../Lean/Geometry/SourceBankInvariant.lean) | 23 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceBankLowering.lean](../Lean/Geometry/SourceBankLowering.lean) | 16 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceBankMachine.lean](../Lean/Geometry/SourceBankMachine.lean) | 22 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceBankPrecision.lean](../Lean/Geometry/SourceBankPrecision.lean) | 4 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceBankTrace.lean](../Lean/Geometry/SourceBankTrace.lean) | 2 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceBusScaling.lean](../Lean/Geometry/SourceBusScaling.lean) | 27 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceCausalBoundary.lean](../Lean/Geometry/SourceCausalBoundary.lean) | 14 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/SourceCheckpointAxiomAudit.lean](../Lean/Geometry/SourceCheckpointAxiomAudit.lean) | 0 | OPH-NATIVE-CHECKPOINT-SELECTION |
| [Lean/Geometry/SourceCheckpointEntropy.lean](../Lean/Geometry/SourceCheckpointEntropy.lean) | 6 | OPH-NATIVE-CHECKPOINT-SELECTION |
| [Lean/Geometry/SourceCheckpointMeaning.lean](../Lean/Geometry/SourceCheckpointMeaning.lean) | 10 | OPH-NATIVE-CHECKPOINT-SELECTION |
| [Lean/Geometry/SourceCheckpointNative.lean](../Lean/Geometry/SourceCheckpointNative.lean) | 9 | OPH-NATIVE-CHECKPOINT-SELECTION |
| [Lean/Geometry/SourceCheckpointPipeline.lean](../Lean/Geometry/SourceCheckpointPipeline.lean) | 5 | OPH-NATIVE-CHECKPOINT-SELECTION |
| [Lean/Geometry/SourceCheckpointPolicy.lean](../Lean/Geometry/SourceCheckpointPolicy.lean) | 12 | OPH-NATIVE-CHECKPOINT-SELECTION |
| [Lean/Geometry/SourceClockRateAlongWorldlines.lean](../Lean/Geometry/SourceClockRateAlongWorldlines.lean) | 53 | OPH-GEOMETRY-SOURCE-CLOCK-RATE-ALONG-WORLDLINES, OPH-GEOMETRY-PROPER-LENGTH-CLOCKED-CHAIN |
| [Lean/Geometry/SourceConstrainedRead.lean](../Lean/Geometry/SourceConstrainedRead.lean) | 2 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceConstrainedSelection.lean](../Lean/Geometry/SourceConstrainedSelection.lean) | 17 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceCountLimitAxiomAudit.lean](../Lean/Geometry/SourceCountLimitAxiomAudit.lean) | 0 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM, OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/SourceCountTransport.lean](../Lean/Geometry/SourceCountTransport.lean) | 13 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/SourceDerivedSpacetimeCarrier.lean](../Lean/Geometry/SourceDerivedSpacetimeCarrier.lean) | 49 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM, OPH-GEOMETRY-EVENT-POPULATION-CHART-INTERFACE |
| [Lean/Geometry/SourceEncodedMemory.lean](../Lean/Geometry/SourceEncodedMemory.lean) | 25 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceEncodedMemoryAxiomAudit.lean](../Lean/Geometry/SourceEncodedMemoryAxiomAudit.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceFeedbackTransport.lean](../Lean/Geometry/SourceFeedbackTransport.lean) | 14 | OPH-SOURCE-FEEDBACK-VERSION-TRANSPORT |
| [Lean/Geometry/SourceNativeAccumulator.lean](../Lean/Geometry/SourceNativeAccumulator.lean) | 19 | OPH-NATIVE-REUSABLE-ACCUMULATION |
| [Lean/Geometry/SourceNativeCore.lean](../Lean/Geometry/SourceNativeCore.lean) | 8 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceNativeProgramBudget.lean](../Lean/Geometry/SourceNativeProgramBudget.lean) | 3 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceNativeProgramError.lean](../Lean/Geometry/SourceNativeProgramError.lean) | 7 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceNativeProgramsAxiomAudit.lean](../Lean/Geometry/SourceNativeProgramsAxiomAudit.lean) | 0 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceNativeRecords.lean](../Lean/Geometry/SourceNativeRecords.lean) | 12 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceNativeShuttle.lean](../Lean/Geometry/SourceNativeShuttle.lean) | 19 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceNativeStoredProgram.lean](../Lean/Geometry/SourceNativeStoredProgram.lean) | 4 | OPH-NATIVE-CONDITIONAL-STORED-PROGRAM |
| [Lean/Geometry/SourceNativeUpdatesAxiomAudit.lean](../Lean/Geometry/SourceNativeUpdatesAxiomAudit.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceNetCausalCone.lean](../Lean/Geometry/SourceNetCausalCone.lean) | 13 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM |
| [Lean/Geometry/SourceNetClockSeparation.lean](../Lean/Geometry/SourceNetClockSeparation.lean) | 5 | OPH-LAYERED-OBSERVER-CLOCK-CONTROLS |
| [Lean/Geometry/SourceNetConformalRecordDensity.lean](../Lean/Geometry/SourceNetConformalRecordDensity.lean) | 22 | OPH-COSMO-FLRW-RECORD-DENSITY |
| [Lean/Geometry/SourceNetLayeredOrder.lean](../Lean/Geometry/SourceNetLayeredOrder.lean) | 43 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM |
| [Lean/Geometry/SourceNetOrderLimit.lean](../Lean/Geometry/SourceNetOrderLimit.lean) | 4 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/SourceNetVolumeError.lean](../Lean/Geometry/SourceNetVolumeError.lean) | 9 | OPH-GOLDEN-SOURCE-QUADRATURE |
| [Lean/Geometry/SourceNonlinearRecord.lean](../Lean/Geometry/SourceNonlinearRecord.lean) | 11 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceOperationReads.lean](../Lean/Geometry/SourceOperationReads.lean) | 15 | OPH-REGISTERED-SOURCE-READ-ALGEBRA |
| [Lean/Geometry/SourceOperationReadsAxiomAudit.lean](../Lean/Geometry/SourceOperationReadsAxiomAudit.lean) | 0 | OPH-REGISTERED-SOURCE-READ-ALGEBRA |
| [Lean/Geometry/SourceOrderEinsteinComposition.lean](../Lean/Geometry/SourceOrderEinsteinComposition.lean) | 14 | OPH-GR-E2E-BRANCH-ENTRY, OPH-GR-SOURCE-DIRECTION-EINSTEIN-COMPOSITION |
| [Lean/Geometry/SourceOrderFrameCompatibilityPacket.lean](../Lean/Geometry/SourceOrderFrameCompatibilityPacket.lean) | 8 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM, OPH-GR-SOURCE-DIRECTION-EINSTEIN-COMPOSITION, OPH-ALGEBRAIC-EVENT-FRAME-SOLDERING, OPH-GEOMETRY-EVENT-POPULATION-CHART-INTERFACE |
| [Lean/Geometry/SourcePassiveMemoryAxiomAudit.lean](../Lean/Geometry/SourcePassiveMemoryAxiomAudit.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourcePassiveMemoryBudget.lean](../Lean/Geometry/SourcePassiveMemoryBudget.lean) | 15 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourcePassiveReset.lean](../Lean/Geometry/SourcePassiveReset.lean) | 14 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourcePopulationQuadrature.lean](../Lean/Geometry/SourcePopulationQuadrature.lean) | 17 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM |
| [Lean/Geometry/SourcePublicationAxiomAudit.lean](../Lean/Geometry/SourcePublicationAxiomAudit.lean) | 0 | OPH-NATIVE-EVENTUAL-PUBLICATION, OPH-FINITE-A3-SUPPORT |
| [Lean/Geometry/SourcePublicationLaw.lean](../Lean/Geometry/SourcePublicationLaw.lean) | 11 | OPH-NATIVE-EVENTUAL-PUBLICATION |
| [Lean/Geometry/SourcePublicationMass.lean](../Lean/Geometry/SourcePublicationMass.lean) | 15 | OPH-NATIVE-EVENTUAL-PUBLICATION |
| [Lean/Geometry/SourcePublicationMenu.lean](../Lean/Geometry/SourcePublicationMenu.lean) | 9 | OPH-NATIVE-EVENTUAL-PUBLICATION |
| [Lean/Geometry/SourcePublicationNative.lean](../Lean/Geometry/SourcePublicationNative.lean) | 8 | OPH-NATIVE-EVENTUAL-PUBLICATION |
| [Lean/Geometry/SourcePublicationPath.lean](../Lean/Geometry/SourcePublicationPath.lean) | 11 | OPH-NATIVE-EVENTUAL-PUBLICATION |
| [Lean/Geometry/SourcePublicationPathInvariant.lean](../Lean/Geometry/SourcePublicationPathInvariant.lean) | 8 | OPH-NATIVE-EVENTUAL-PUBLICATION |
| [Lean/Geometry/SourcePublicationTail.lean](../Lean/Geometry/SourcePublicationTail.lean) | 5 | OPH-NATIVE-EVENTUAL-PUBLICATION |
| [Lean/Geometry/SourceRadiusRemoval.lean](../Lean/Geometry/SourceRadiusRemoval.lean) | 7 | OPH-NATIVE-EVENTUAL-PUBLICATION |
| [Lean/Geometry/SourceRadiusSelection.lean](../Lean/Geometry/SourceRadiusSelection.lean) | 6 | OPH-NATIVE-EVENTUAL-PUBLICATION |
| [Lean/Geometry/SourceReadAcceptance.lean](../Lean/Geometry/SourceReadAcceptance.lean) | 48 | OPH-NATIVE-LOCAL-SINGLETON-PUBLICATION |
| [Lean/Geometry/SourceReadAcceptanceAxiomAudit.lean](../Lean/Geometry/SourceReadAcceptanceAxiomAudit.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceReadAcceptanceSchedule.lean](../Lean/Geometry/SourceReadAcceptanceSchedule.lean) | 13 | OPH-NATIVE-LOCAL-SINGLETON-PUBLICATION |
| [Lean/Geometry/SourceReadRouting.lean](../Lean/Geometry/SourceReadRouting.lean) | 14 | OPH-SOURCE-FULL-FAMILY-READ-ROUTING |
| [Lean/Geometry/SourceReadSelection.lean](../Lean/Geometry/SourceReadSelection.lean) | 14 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceReadSelectionAxiomAudit.lean](../Lean/Geometry/SourceReadSelectionAxiomAudit.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceRecordGluing.lean](../Lean/Geometry/SourceRecordGluing.lean) | 20 | OPH-RECORD-GLUING-M1-SUFFICIENCY |
| [Lean/Geometry/SourceRecordGluingAxiomAudit.lean](../Lean/Geometry/SourceRecordGluingAxiomAudit.lean) | 0 | OPH-RECORD-GLUING-M1-SUFFICIENCY |
| [Lean/Geometry/SourceRecordProtection.lean](../Lean/Geometry/SourceRecordProtection.lean) | 11 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM |
| [Lean/Geometry/SourceRepairInstrument.lean](../Lean/Geometry/SourceRepairInstrument.lean) | 17 | OPH-REGISTERED-SOURCE-READ-ALGEBRA |
| [Lean/Geometry/SourceReusableBus.lean](../Lean/Geometry/SourceReusableBus.lean) | 25 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceReusableBusAxiomAudit.lean](../Lean/Geometry/SourceReusableBusAxiomAudit.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceRoutingBudget.lean](../Lean/Geometry/SourceRoutingBudget.lean) | 7 | OPH-SOURCE-ROUTING-STORAGE-AND-RAW-COUNT |
| [Lean/Geometry/SourceRoutingHierarchy.lean](../Lean/Geometry/SourceRoutingHierarchy.lean) | 12 | OPH-SOURCE-ROUTING-STORAGE-AND-RAW-COUNT |
| [Lean/Geometry/SourceRoutingRefinementAxiomAudit.lean](../Lean/Geometry/SourceRoutingRefinementAxiomAudit.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceRoutingStorage.lean](../Lean/Geometry/SourceRoutingStorage.lean) | 7 | OPH-SOURCE-ROUTING-STORAGE-AND-RAW-COUNT |
| [Lean/Geometry/SourceSeamPathTomography.lean](../Lean/Geometry/SourceSeamPathTomography.lean) | 9 | OPH-SEAM-PATH-TOMOGRAPHY |
| [Lean/Geometry/SourceSelectionControls.lean](../Lean/Geometry/SourceSelectionControls.lean) | 15 | No claim evidence edge; model mapping open |
| [Lean/Geometry/SourceSelectionLocality.lean](../Lean/Geometry/SourceSelectionLocality.lean) | 15 | OPH-SCALAR-SEAM-SOURCE-NONENTAILMENT |
| [Lean/Geometry/SourceSelectionLocalityAxiomAudit.lean](../Lean/Geometry/SourceSelectionLocalityAxiomAudit.lean) | 0 | OPH-SCALAR-SEAM-SOURCE-NONENTAILMENT |
| [Lean/Geometry/SourceStateSelection.lean](../Lean/Geometry/SourceStateSelection.lean) | 9 | OPH-FINITE-A3-SUPPORT |
| [Lean/Geometry/SourceTemporalAcceptanceAxiomAudit.lean](../Lean/Geometry/SourceTemporalAcceptanceAxiomAudit.lean) | 0 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SourceTemporalAttempts.lean](../Lean/Geometry/SourceTemporalAttempts.lean) | 7 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SourceTemporalConnected.lean](../Lean/Geometry/SourceTemporalConnected.lean) | 5 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SourceTemporalErasure.lean](../Lean/Geometry/SourceTemporalErasure.lean) | 12 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SourceTemporalEssential.lean](../Lean/Geometry/SourceTemporalEssential.lean) | 5 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SourceTemporalGuard.lean](../Lean/Geometry/SourceTemporalGuard.lean) | 18 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SourceTemporalMenu.lean](../Lean/Geometry/SourceTemporalMenu.lean) | 8 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SourceTemporalNative.lean](../Lean/Geometry/SourceTemporalNative.lean) | 8 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SourceTemporalObservation.lean](../Lean/Geometry/SourceTemporalObservation.lean) | 10 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SourceTemporalSelection.lean](../Lean/Geometry/SourceTemporalSelection.lean) | 5 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SourceTemporalTomography.lean](../Lean/Geometry/SourceTemporalTomography.lean) | 20 | OPH-TEMPORAL-NATIVE-RECOVERY |
| [Lean/Geometry/SpatialReadbackSoldering.lean](../Lean/Geometry/SpatialReadbackSoldering.lean) | 9 | OPH-ALGEBRAIC-EVENT-FRAME-SOLDERING, OPH-GEOMETRY-REST-FIBER-SHELL-TRANSPORT |
| [Lean/Geometry/TimelikeClassForceLaw.lean](../Lean/Geometry/TimelikeClassForceLaw.lean) | 46 | OPH-GEOMETRY-TIMELIKE-CLASS-FORCE-LAW |
| [Lean/Geometry/TransportedChargeForceLaw.lean](../Lean/Geometry/TransportedChargeForceLaw.lean) | 42 | OPH-GEOMETRY-TRANSPORTED-CHARGE-FORCE-LAW, OPH-GEOMETRY-TIMELIKE-CLASS-FORCE-LAW |
| [Lean/Geometry/WorldlineHopTransport.lean](../Lean/Geometry/WorldlineHopTransport.lean) | 58 | OPH-GEOMETRY-WORLDLINE-HOP-TRANSPORT, OPH-GEOMETRY-TRANSPORTED-CHARGE-FORCE-LAW, OPH-GEOMETRY-SEAM-STEP-SPEED-LIMIT, OPH-GEOMETRY-SOURCE-CLOCK-RATE-ALONG-WORLDLINES |
| [Lean/Geometry.lean](../Lean/Geometry.lean) | 0 | OPH-QFT-CAUSAL-ORDER-COMPOSITION |
| [Lean/InformationProjection/GlobalObjective.lean](../Lean/InformationProjection/GlobalObjective.lean) | 39 | OPH-THERMO-FOUR-LAW-PACKAGE |
| [Lean/InformationProjection/HistoryLaw.lean](../Lean/InformationProjection/HistoryLaw.lean) | 48 | No claim evidence edge; model mapping open |
| [Lean/InformationProjection/LogTransitionAction.lean](../Lean/InformationProjection/LogTransitionAction.lean) | 60 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/InformationProjection/PathGibbs.lean](../Lean/InformationProjection/PathGibbs.lean) | 10 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/InformationProjection/ReferenceNormalForm.lean](../Lean/InformationProjection/ReferenceNormalForm.lean) | 15 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/InformationProjection/SourceHistoryPacket.lean](../Lean/InformationProjection/SourceHistoryPacket.lean) | 101 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS, OPH-QFT-SOURCE-HISTORY-GNS-HAMILTONIAN-ATTACHMENT, OPH-QFT-SOURCE-CLOCK-CANDIDATE |
| [Lean/InformationProjection/SourceReferenceSelection.lean](../Lean/InformationProjection/SourceReferenceSelection.lean) | 14 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/InformationProjection.lean](../Lean/InformationProjection.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Locality.lean](../Lean/Locality.lean) | 0 | OPH-FINITE-LOCALITY-NOSIGNALLING |
| [Lean/Main.lean](../Lean/Main.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms/AxiomAudit.lean](../Lean/ObservableNormalForms/ObservableNormalForms/AxiomAudit.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms/ConditionalResampling.lean](../Lean/ObservableNormalForms/ObservableNormalForms/ConditionalResampling.lean) | 34 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms/Exact.lean](../Lean/ObservableNormalForms/ObservableNormalForms/Exact.lean) | 15 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms/Examples/AmdA2A5Conditional.lean](../Lean/ObservableNormalForms/ObservableNormalForms/Examples/AmdA2A5Conditional.lean) | 52 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms/Examples/ProtectedObstructions.lean](../Lean/ObservableNormalForms/ObservableNormalForms/Examples/ProtectedObstructions.lean) | 114 | OPH-CONS-PROTECTED-BEHAVIOR-PROFILE |
| [Lean/ObservableNormalForms/ObservableNormalForms/Examples/Rule90.lean](../Lean/ObservableNormalForms/ObservableNormalForms/Examples/Rule90.lean) | 7 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms/Functional.lean](../Lean/ObservableNormalForms/ObservableNormalForms/Functional.lean) | 5 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms/MechanismVariants.lean](../Lean/ObservableNormalForms/ObservableNormalForms/MechanismVariants.lean) | 8 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms/Moduli.lean](../Lean/ObservableNormalForms/ObservableNormalForms/Moduli.lean) | 56 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms/ObserverConfluence.lean](../Lean/ObservableNormalForms/ObservableNormalForms/ObserverConfluence.lean) | 13 | OPH-CONS-CROSS-SOURCE-ENDPOINT |
| [Lean/ObservableNormalForms/ObservableNormalForms/ProtectedObstructions.lean](../Lean/ObservableNormalForms/ObservableNormalForms/ProtectedObstructions.lean) | 126 | OPH-CONS-PROTECTED-BEHAVIOR-PROFILE |
| [Lean/ObservableNormalForms/ObservableNormalForms/Refinement.lean](../Lean/ObservableNormalForms/ObservableNormalForms/Refinement.lean) | 8 | OPH-CONS-REFINEMENT-MODULI |
| [Lean/ObservableNormalForms/ObservableNormalForms/Repair.lean](../Lean/ObservableNormalForms/ObservableNormalForms/Repair.lean) | 7 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms/SchedulerClassObstructions.lean](../Lean/ObservableNormalForms/ObservableNormalForms/SchedulerClassObstructions.lean) | 40 | OPH-CONS-SCHEDULER-CLASS-CAPACITY-SPLIT |
| [Lean/ObservableNormalForms/ObservableNormalForms/Stability.lean](../Lean/ObservableNormalForms/ObservableNormalForms/Stability.lean) | 4 | OPH-CONS-REFINEMENT-MODULI |
| [Lean/ObservableNormalForms/ObservableNormalForms/Stochastic.lean](../Lean/ObservableNormalForms/ObservableNormalForms/Stochastic.lean) | 38 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/ObservableNormalForms.lean](../Lean/ObservableNormalForms/ObservableNormalForms.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/ObservableNormalForms/lakefile.lean](../Lean/ObservableNormalForms/lakefile.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/A2EndpointCommutator.lean](../Lean/ObserverPatchHolography/A2EndpointCommutator.lean) | 19 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/AbstractRewriting.lean](../Lean/ObserverPatchHolography/AbstractRewriting.lean) | 8 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/BoundaryFiber.lean](../Lean/ObserverPatchHolography/BoundaryFiber.lean) | 10 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/BridgeBoundaries.lean](../Lean/ObserverPatchHolography/BridgeBoundaries.lean) | 12 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/BridgeEquivalence.lean](../Lean/ObserverPatchHolography/BridgeEquivalence.lean) | 5 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/Bridges/ObservableNormalForms.lean](../Lean/ObserverPatchHolography/Bridges/ObservableNormalForms.lean) | 3 | OPH-CONS-CROSS-SOURCE-ENDPOINT |
| [Lean/ObserverPatchHolography/CapacityClosurePrinciple.lean](../Lean/ObserverPatchHolography/CapacityClosurePrinciple.lean) | 7 | OPH-Q-N-RESERVE-CANDIDATES |
| [Lean/ObserverPatchHolography/CapacityFixedPoint.lean](../Lean/ObserverPatchHolography/CapacityFixedPoint.lean) | 4 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/CapacityNonidentifiability.lean](../Lean/ObserverPatchHolography/CapacityNonidentifiability.lean) | 15 | OPH-GR-D6-CAPACITY, OPH-GR-D6-READBACK-NOGO, OPH-Q-PN-CLOSURE, OPH-Q-N-RESERVE-CANDIDATES |
| [Lean/ObserverPatchHolography/ClebschRatio.lean](../Lean/ObserverPatchHolography/ClebschRatio.lean) | 6 | OPH-QUARK-REGISTER-CLEBSCH |
| [Lean/ObserverPatchHolography/ClosurePreflight.lean](../Lean/ObserverPatchHolography/ClosurePreflight.lean) | 4 | OPH-Q-PN-CLOSURE, OPH-Q-N-RESERVE-CANDIDATES, OPH-ALPHA-PIXEL |
| [Lean/ObserverPatchHolography/CollarClause.lean](../Lean/ObserverPatchHolography/CollarClause.lean) | 15 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/CollarLayer.lean](../Lean/ObserverPatchHolography/CollarLayer.lean) | 41 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/CollarModularT2.lean](../Lean/ObserverPatchHolography/CollarModularT2.lean) | 32 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/CollarStates.lean](../Lean/ObserverPatchHolography/CollarStates.lean) | 91 | OPH-THERMO-COLLAR-TEMPERATURE-READING |
| [Lean/ObserverPatchHolography/CollarStatesBridge.lean](../Lean/ObserverPatchHolography/CollarStatesBridge.lean) | 8 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/CollarStatesT1.lean](../Lean/ObserverPatchHolography/CollarStatesT1.lean) | 28 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/CoreAxioms.lean](../Lean/ObserverPatchHolography/CoreAxioms.lean) | 7 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/DeSitterCapacityShock.lean](../Lean/ObserverPatchHolography/DeSitterCapacityShock.lean) | 9 | OPH-GR-DS-MU2-IDENTITY, OPH-GR-DS-CAPACITY-TRANSFER |
| [Lean/ObserverPatchHolography/DirectedSeamRepair.lean](../Lean/ObserverPatchHolography/DirectedSeamRepair.lean) | 14 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/DirectedSeamRepairProgress.lean](../Lean/ObserverPatchHolography/DirectedSeamRepairProgress.lean) | 13 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/EinsteinBranch/AxiomAudit.lean](../Lean/ObserverPatchHolography/EinsteinBranch/AxiomAudit.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/EinsteinBranch/CollarPremiseDerivation.lean](../Lean/ObserverPatchHolography/EinsteinBranch/CollarPremiseDerivation.lean) | 37 | OPH-DM-COLLAR-PREMISE-DERIVATION |
| [Lean/ObserverPatchHolography/EinsteinBranch/CollarScaleSaturation.lean](../Lean/ObserverPatchHolography/EinsteinBranch/CollarScaleSaturation.lean) | 16 | OPH-DM-CONT, OPH-DM-COLLAR-SATURATION |
| [Lean/ObserverPatchHolography/EinsteinBranch/Composition.lean](../Lean/ObserverPatchHolography/EinsteinBranch/Composition.lean) | 3 | OPH-GR-E2E-BRANCH-ENTRY, OPH-GR-SOURCE-DIRECTION-EINSTEIN-COMPOSITION |
| [Lean/ObserverPatchHolography/EinsteinBranch/Consensus.lean](../Lean/ObserverPatchHolography/EinsteinBranch/Consensus.lean) | 8 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/EinsteinBranch/CosmologyLedgerBrackets.lean](../Lean/ObserverPatchHolography/EinsteinBranch/CosmologyLedgerBrackets.lean) | 36 | OPH-COSMO-POSTDICTION-LEDGER |
| [Lean/ObserverPatchHolography/EinsteinBranch/DarkSector.lean](../Lean/ObserverPatchHolography/EinsteinBranch/DarkSector.lean) | 18 | OPH-DM-CONT, OPH-DM-COLLAR-SATURATION, OPH-DM-DEEP-PROFILE-CHARACTERIZATION, OPH-DM-SPARSE-DEFECT-MECHANISM, OPH-DM-FINITE-CAP-GENERATOR-SPLIT |
| [Lean/ObserverPatchHolography/EinsteinBranch/DeepProfileClosure.lean](../Lean/ObserverPatchHolography/EinsteinBranch/DeepProfileClosure.lean) | 24 | OPH-DM-CONT, OPH-DM-DEEP-PROFILE-CHARACTERIZATION, OPH-DM-PER-CUT-COMPOSITION |
| [Lean/ObserverPatchHolography/EinsteinBranch/EdgeCenterTiltCocycle.lean](../Lean/ObserverPatchHolography/EinsteinBranch/EdgeCenterTiltCocycle.lean) | 16 | OPH-COSMO-POSTDICTION-LEDGER |
| [Lean/ObserverPatchHolography/EinsteinBranch/EdgeCenterTiltFreedom.lean](../Lean/ObserverPatchHolography/EinsteinBranch/EdgeCenterTiltFreedom.lean) | 3 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean](../Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/EinsteinBranch/FiniteCapGeneratorSplit.lean](../Lean/ObserverPatchHolography/EinsteinBranch/FiniteCapGeneratorSplit.lean) | 29 | OPH-DM-FINITE-CAP-GENERATOR-SPLIT |
| [Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean](../Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean) | 15 | OPH-DE-FIXED-CAPACITY-WLAW |
| [Lean/ObserverPatchHolography/EinsteinBranch/Manifest.lean](../Lean/ObserverPatchHolography/EinsteinBranch/Manifest.lean) | 4 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/EinsteinBranch/PerCutCollarComposition.lean](../Lean/ObserverPatchHolography/EinsteinBranch/PerCutCollarComposition.lean) | 10 | OPH-DM-PER-CUT-COMPOSITION, OPH-DM-COLLAR-PREMISE-DERIVATION |
| [Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean](../Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean) | 6 | OPH-GR-EINSTEIN-BRANCH-INHABITANT |
| [Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean](../Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean) | 5 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/EinsteinBranch/SparseRecordDefectWitness.lean](../Lean/ObserverPatchHolography/EinsteinBranch/SparseRecordDefectWitness.lean) | 34 | OPH-DM-SPARSE-DEFECT-MECHANISM |
| [Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean](../Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean) | 37 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/EinsteinBranch.lean](../Lean/ObserverPatchHolography/EinsteinBranch.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/EqualSeamSelection.lean](../Lean/ObserverPatchHolography/EqualSeamSelection.lean) | 11 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/Execution/AdaptiveRunCounterexamples.lean](../Lean/ObserverPatchHolography/Execution/AdaptiveRunCounterexamples.lean) | 8 | OPH-CONS-ADAPTIVE-REPAIR-STRATIFICATION |
| [Lean/ObserverPatchHolography/Execution/AdaptiveRunStratification.lean](../Lean/ObserverPatchHolography/Execution/AdaptiveRunStratification.lean) | 16 | OPH-CONS-ADAPTIVE-REPAIR-STRATIFICATION, OPH-CONS-ADAPTIVE-PUBLIC-ENDPOINT |
| [Lean/ObserverPatchHolography/Execution/CumulativeAttemptCapacity.lean](../Lean/ObserverPatchHolography/Execution/CumulativeAttemptCapacity.lean) | 11 | OPH-CONS-CUMULATIVE-ATTEMPT-CAPACITY, OPH-CONS-CAPACITY-BOUNDED-PUBLIC-ENDPOINT |
| [Lean/ObserverPatchHolography/Execution/CumulativeAttemptCapacityExamples.lean](../Lean/ObserverPatchHolography/Execution/CumulativeAttemptCapacityExamples.lean) | 26 | OPH-CONS-CUMULATIVE-ATTEMPT-CAPACITY |
| [Lean/ObserverPatchHolography/Execution/RankedAttemptCapacity.lean](../Lean/ObserverPatchHolography/Execution/RankedAttemptCapacity.lean) | 8 | OPH-COMPUTATION-FIXED-FEDERATION-EXECUTION-CLASSIFICATION |
| [Lean/ObserverPatchHolography/GlobalCapacityAttachment.lean](../Lean/ObserverPatchHolography/GlobalCapacityAttachment.lean) | 15 | OPH-Q-N-RESERVE-CANDIDATES |
| [Lean/ObserverPatchHolography/IcosahedralAntibridge.lean](../Lean/ObserverPatchHolography/IcosahedralAntibridge.lean) | 8 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/IcosahedralAxisNoGo.lean](../Lean/ObserverPatchHolography/IcosahedralAxisNoGo.lean) | 8 | OPH-A5-REAL-AXIS-CABIBBO-NOGO |
| [Lean/ObserverPatchHolography/IcosahedralOrbitStabilizer.lean](../Lean/ObserverPatchHolography/IcosahedralOrbitStabilizer.lean) | 4 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/KoideCirculant.lean](../Lean/ObserverPatchHolography/KoideCirculant.lean) | 8 | OPH-KOIDE-CIRCULANT-IDENTITY |
| [Lean/ObserverPatchHolography/Locality/AdaptiveScheduler.lean](../Lean/ObserverPatchHolography/Locality/AdaptiveScheduler.lean) | 28 | OPH-ADAPTIVE-SCHEDULER-LOCALITY-HELPER |
| [Lean/ObserverPatchHolography/Locality/DependencyCone.lean](../Lean/ObserverPatchHolography/Locality/DependencyCone.lean) | 16 | OPH-FINITE-LOCALITY-NOSIGNALLING |
| [Lean/ObserverPatchHolography/Locality/NoSignalling.lean](../Lean/ObserverPatchHolography/Locality/NoSignalling.lean) | 5 | OPH-FINITE-LOCALITY-NOSIGNALLING |
| [Lean/ObserverPatchHolography/Primitives.lean](../Lean/ObserverPatchHolography/Primitives.lean) | 58 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/Provenance/CausalInterval.lean](../Lean/ObserverPatchHolography/Provenance/CausalInterval.lean) | 17 | OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE |
| [Lean/ObserverPatchHolography/Provenance/FiniteCausetCompiler.lean](../Lean/ObserverPatchHolography/Provenance/FiniteCausetCompiler.lean) | 9 | OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE, OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM |
| [Lean/ObserverPatchHolography/Provenance/FiniteCausetCoverCompiler.lean](../Lean/ObserverPatchHolography/Provenance/FiniteCausetCoverCompiler.lean) | 6 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean](../Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean) | 58 | OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE |
| [Lean/ObserverPatchHolography/Provenance/MismatchProvenance.lean](../Lean/ObserverPatchHolography/Provenance/MismatchProvenance.lean) | 8 | OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE |
| [Lean/ObserverPatchHolography/Provenance/QuotientInvariance.lean](../Lean/ObserverPatchHolography/Provenance/QuotientInvariance.lean) | 10 | OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE |
| [Lean/ObserverPatchHolography/Provenance/RefinementNaturality.lean](../Lean/ObserverPatchHolography/Provenance/RefinementNaturality.lean) | 7 | OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE |
| [Lean/ObserverPatchHolography/Provenance/SeamDeltaAggregation.lean](../Lean/ObserverPatchHolography/Provenance/SeamDeltaAggregation.lean) | 7 | OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE |
| [Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean](../Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean) | 25 | OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE, OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM, OPH-QFT-RESTRICTION-NATURALITY |
| [Lean/ObserverPatchHolography/RateBridgeObstruction.lean](../Lean/ObserverPatchHolography/RateBridgeObstruction.lean) | 19 | OPH-QFT-SOURCE-RECURRENCE-CLOCK-PRECURSOR |
| [Lean/ObserverPatchHolography/RateNonidentifiability.lean](../Lean/ObserverPatchHolography/RateNonidentifiability.lean) | 39 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/RepairGapChain.lean](../Lean/ObserverPatchHolography/RepairGapChain.lean) | 3 | OPH-YM-GAP |
| [Lean/ObserverPatchHolography/RepairWordSchedule.lean](../Lean/ObserverPatchHolography/RepairWordSchedule.lean) | 12 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/Rule90.lean](../Lean/ObserverPatchHolography/Rule90.lean) | 5 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/ScalarSeamRepair.lean](../Lean/ObserverPatchHolography/ScalarSeamRepair.lean) | 14 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/SeedPi.lean](../Lean/ObserverPatchHolography/SeedPi.lean) | 4 | No claim evidence edge; model mapping open |
| [Lean/ObserverPatchHolography/W5Stabilizer.lean](../Lean/ObserverPatchHolography/W5Stabilizer.lean) | 10 | OPH-W5-STABILISER-POTENTIAL-BOUNDARY |
| [Lean/ObserverPatchHolography/YangMillsGap.lean](../Lean/ObserverPatchHolography/YangMillsGap.lean) | 5 | OPH-YM-GAP |
| [Lean/ObserverPatchHolography/YangMillsGapWitness.lean](../Lean/ObserverPatchHolography/YangMillsGapWitness.lean) | 15 | OPH-YM-GAP |
| [Lean/ObserverPatchHolography/YangMillsLemma72.lean](../Lean/ObserverPatchHolography/YangMillsLemma72.lean) | 5 | OPH-YM-GAP |
| [Lean/ObserverPatchHolography/YangMillsProp81.lean](../Lean/ObserverPatchHolography/YangMillsProp81.lean) | 6 | OPH-YM-GAP |
| [Lean/ObserverPatchHolography.lean](../Lean/ObserverPatchHolography.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Checker.lean](../Lean/OphGap/Checker.lean) | 18 | No claim evidence edge; model mapping open |
| [Lean/OphGap/OphData.lean](../Lean/OphGap/OphData.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/OphGap/OphOrigControls.lean](../Lean/OphGap/OphOrigControls.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/OphGap/OphOrigData.lean](../Lean/OphGap/OphOrigData.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/OphGap/OphOrigWitness.lean](../Lean/OphGap/OphOrigWitness.lean) | 18 | OPH-SOURCE-DIMENSIONLESS-GAP |
| [Lean/OphGap/OphWitness.lean](../Lean/OphGap/OphWitness.lean) | 8 | No claim evidence edge; model mapping open |
| [Lean/OphGap/SignedGap.lean](../Lean/OphGap/SignedGap.lean) | 25 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/E1.lean](../Lean/OphGap/Slices/E1.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/E2.lean](../Lean/OphGap/Slices/E2.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/E3.lean](../Lean/OphGap/Slices/E3.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/N1.lean](../Lean/OphGap/Slices/N1.lean) | 6 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/N2.lean](../Lean/OphGap/Slices/N2.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/N3.lean](../Lean/OphGap/Slices/N3.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/P1.lean](../Lean/OphGap/Slices/P1.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/P2.lean](../Lean/OphGap/Slices/P2.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/P3.lean](../Lean/OphGap/Slices/P3.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/P4.lean](../Lean/OphGap/Slices/P4.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/T1.lean](../Lean/OphGap/Slices/T1.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/T2.lean](../Lean/OphGap/Slices/T2.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Slices/T3.lean](../Lean/OphGap/Slices/T3.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/OphGap/Transport.lean](../Lean/OphGap/Transport.lean) | 37 | No claim evidence edge; model mapping open |
| [Lean/OphGap.lean](../Lean/OphGap.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/QFT/CPRestrictionNet.lean](../Lean/QFT/CPRestrictionNet.lean) | 8 | OPH-FINITE-LOCALITY-NOSIGNALLING, OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT/CarrierJoinTransport.lean](../Lean/QFT/CarrierJoinTransport.lean) | 78 | OPH-QUANTUM-HINGE-FIBERED-JOIN-AND-SLICE-REGION-DICHOTOMY, OPH-QFT-JOIN-RECORD-NET-MORPHISM |
| [Lean/QFT/ColimitGNSRepresentation.lean](../Lean/QFT/ColimitGNSRepresentation.lean) | 10 | OPH-QFT-COMPLETED-COLIMIT-GNS-INTERFACE, OPH-QFT-SELECTED-STATE-GNS-ATTACHMENT, OPH-QFT-SELECTED-GNS-ISING-HAMILTONIAN-COMPOSITION |
| [Lean/QFT/ColimitNormedCompletion.lean](../Lean/QFT/ColimitNormedCompletion.lean) | 80 | OPH-QFT-COLIMIT-NORMED-COMPLETION, OPH-QFT-COMPLETED-COLIMIT-GNS-INTERFACE, OPH-QFT-SELECTED-STATE-GNS-ATTACHMENT |
| [Lean/QFT/ColimitSelectedStateGNS.lean](../Lean/QFT/ColimitSelectedStateGNS.lean) | 9 | OPH-QFT-SELECTED-STATE-GNS-ATTACHMENT, OPH-QFT-SELECTED-GNS-ISING-HAMILTONIAN-COMPOSITION, OPH-QFT-SOURCE-HISTORY-GNS-HAMILTONIAN-ATTACHMENT |
| [Lean/QFT/ConjugationGauge.lean](../Lean/QFT/ConjugationGauge.lean) | 22 | OPH-FINITE-BORN-FRAME-RANK-GAP |
| [Lean/QFT/ConvexAffinityBridge.lean](../Lean/QFT/ConvexAffinityBridge.lean) | 35 | No claim evidence edge; model mapping open |
| [Lean/QFT/CoverageReceipt.lean](../Lean/QFT/CoverageReceipt.lean) | 33 | OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT/EffectiveQuantumComparison.lean](../Lean/QFT/EffectiveQuantumComparison.lean) | 13 | OPH-QFT-EFFECTIVE-QUANTUM-COMPARISON |
| [Lean/QFT/EnlargedTargetRecordEmbedding.lean](../Lean/QFT/EnlargedTargetRecordEmbedding.lean) | 53 | OPH-QFT-ENLARGED-TARGET-RECORD-EMBEDDING |
| [Lean/QFT/FiniteCausalObserverNet.lean](../Lean/QFT/FiniteCausalObserverNet.lean) | 12 | OPH-FINITE-CAUSAL-OBSERVER-NET, OPH-QFT-CAUSAL-ORDER-COMPOSITION |
| [Lean/QFT/FiniteTwoSiteIsingColimitGNS.lean](../Lean/QFT/FiniteTwoSiteIsingColimitGNS.lean) | 44 | OPH-QFT-SELECTED-GNS-ISING-HAMILTONIAN-COMPOSITION |
| [Lean/QFT/FiniteTwoSiteIsingField.lean](../Lean/QFT/FiniteTwoSiteIsingField.lean) | 40 | OPH-QFT-FINITE-TWO-SITE-ISING-HILBERT-DYNAMICS, OPH-QFT-SELECTED-GNS-ISING-HAMILTONIAN-COMPOSITION, OPH-QFT-SOURCE-HISTORY-THREE-SLOT-LOCAL-GNS-DYNAMICS |
| [Lean/QFT/FiniteUnitaryScatteringNoGo.lean](../Lean/QFT/FiniteUnitaryScatteringNoGo.lean) | 6 | OPH-FINITE-UNITARY-SCATTERING-LIMIT-NO-GO |
| [Lean/QFT/GaugeIrreducibleBorn.lean](../Lean/QFT/GaugeIrreducibleBorn.lean) | 71 | No claim evidence edge; model mapping open |
| [Lean/QFT/GlueCountermodel.lean](../Lean/QFT/GlueCountermodel.lean) | 25 | No claim evidence edge; model mapping open |
| [Lean/QFT/HistoryActionBinding.lean](../Lean/QFT/HistoryActionBinding.lean) | 38 | OPH-GAUGE-KINETIC-HISTORY-BINDING |
| [Lean/QFT/InhabitedTowerRow.lean](../Lean/QFT/InhabitedTowerRow.lean) | 21 | No claim evidence edge; model mapping open |
| [Lean/QFT/InheritanceMatrix.lean](../Lean/QFT/InheritanceMatrix.lean) | 13 | OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX |
| [Lean/QFT/InheritanceMatrixGuards.lean](../Lean/QFT/InheritanceMatrixGuards.lean) | 5 | No claim evidence edge; model mapping open |
| [Lean/QFT/JoinNetMorphism.lean](../Lean/QFT/JoinNetMorphism.lean) | 71 | OPH-QFT-JOIN-RECORD-NET-MORPHISM, OPH-QFT-OPERATOR-LAYER-MORPHISM-BOUNDARY |
| [Lean/QFT/JointInstance.lean](../Lean/QFT/JointInstance.lean) | 12 | No claim evidence edge; model mapping open |
| [Lean/QFT/JointSlotFactorisation.lean](../Lean/QFT/JointSlotFactorisation.lean) | 22 | OPH-FINITE-LOCALITY-NOSIGNALLING, OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT/LimitRestriction.lean](../Lean/QFT/LimitRestriction.lean) | 11 | OPH-QFT-RESTRICTION-NATURALITY |
| [Lean/QFT/LocallyCovariantLimit.lean](../Lean/QFT/LocallyCovariantLimit.lean) | 75 | OPH-QFT-COLIMIT-NORMED-COMPLETION, OPH-QFT-SELECTED-STATE-GNS-ATTACHMENT, OPH-QFT-RESTRICTION-NATURALITY |
| [Lean/QFT/NoncommutativeWitness.lean](../Lean/QFT/NoncommutativeWitness.lean) | 38 | No claim evidence edge; model mapping open |
| [Lean/QFT/ObserverAccessCut.lean](../Lean/QFT/ObserverAccessCut.lean) | 44 | No claim evidence edge; model mapping open |
| [Lean/QFT/ObserverEventTower.lean](../Lean/QFT/ObserverEventTower.lean) | 18 | No claim evidence edge; model mapping open |
| [Lean/QFT/ObserverEventWorld.lean](../Lean/QFT/ObserverEventWorld.lean) | 11 | No claim evidence edge; model mapping open |
| [Lean/QFT/ObserverNetDescent.lean](../Lean/QFT/ObserverNetDescent.lean) | 7 | OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT/OperationalOverlapEvidence.lean](../Lean/QFT/OperationalOverlapEvidence.lean) | 26 | OPH-OBS-PATCH-DEF |
| [Lean/QFT/OperatorLayerMorphismBoundary.lean](../Lean/QFT/OperatorLayerMorphismBoundary.lean) | 17 | OPH-QFT-OPERATOR-LAYER-MORPHISM-BOUNDARY |
| [Lean/QFT/PathTimeSliceInterface.lean](../Lean/QFT/PathTimeSliceInterface.lean) | 28 | OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX, OPH-QUANTUM-TIME-SLICE-INTERFACE-DELIMITATION |
| [Lean/QFT/PublicEndpointBridge.lean](../Lean/QFT/PublicEndpointBridge.lean) | 41 | No claim evidence edge; model mapping open |
| [Lean/QFT/RichFibreRegionalNet.lean](../Lean/QFT/RichFibreRegionalNet.lean) | 65 | OPH-FINITE-LOCALITY-NOSIGNALLING, OPH-FINITE-CAUSAL-OBSERVER-NET, OPH-QFT-CAUSAL-ORDER-COMPOSITION |
| [Lean/QFT/RichFibreWitness.lean](../Lean/QFT/RichFibreWitness.lean) | 22 | OPH-FINITE-LOCALITY-NOSIGNALLING, OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT/ScalarRegionalTimeSlice.lean](../Lean/QFT/ScalarRegionalTimeSlice.lean) | 5 | OPH-SCALAR-REGIONAL-COLLAR-RECONSTRUCTION |
| [Lean/QFT/SectorInheritance.lean](../Lean/QFT/SectorInheritance.lean) | 21 | OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX |
| [Lean/QFT/SimEarnedWitness.lean](../Lean/QFT/SimEarnedWitness.lean) | 70 | No claim evidence edge; model mapping open |
| [Lean/QFT/SourceClockCandidate.lean](../Lean/QFT/SourceClockCandidate.lean) | 34 | OPH-QFT-SOURCE-CLOCK-CANDIDATE |
| [Lean/QFT/SourceContextTomographyNoGo.lean](../Lean/QFT/SourceContextTomographyNoGo.lean) | 25 | OPH-FINITE-BORN-FRAME-RANK-GAP |
| [Lean/QFT/SourceContextWeb.lean](../Lean/QFT/SourceContextWeb.lean) | 41 | No claim evidence edge; model mapping open |
| [Lean/QFT/SourceConvexFamily.lean](../Lean/QFT/SourceConvexFamily.lean) | 24 | No claim evidence edge; model mapping open |
| [Lean/QFT/SourceCorrelationCapstone.lean](../Lean/QFT/SourceCorrelationCapstone.lean) | 26 | OPH-FINITE-LOCALITY-NOSIGNALLING, OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT/SourceDerivedEventPrecedence.lean](../Lean/QFT/SourceDerivedEventPrecedence.lean) | 7 | OPH-CONS-SOURCE-DERIVED-EVENT-PRECEDENCE |
| [Lean/QFT/SourceHistoryExpectationNet.lean](../Lean/QFT/SourceHistoryExpectationNet.lean) | 98 | OPH-QFT-SOURCE-HISTORY-STATE-PRESERVING-RECORD-NET, OPH-QFT-ENLARGED-TARGET-RECORD-EMBEDDING |
| [Lean/QFT/SourceHistoryGNSDynamics.lean](../Lean/QFT/SourceHistoryGNSDynamics.lean) | 34 | OPH-QFT-SOURCE-HISTORY-GNS-HAMILTONIAN-ATTACHMENT, OPH-QFT-SOURCE-HISTORY-THREE-SLOT-LOCAL-GNS-DYNAMICS, OPH-QFT-SOURCE-HISTORY-STATE-PRESERVING-RECORD-NET |
| [Lean/QFT/SourceHistoryThreeSlotLocalGNS.lean](../Lean/QFT/SourceHistoryThreeSlotLocalGNS.lean) | 60 | OPH-QFT-SOURCE-HISTORY-THREE-SLOT-LOCAL-GNS-DYNAMICS, OPH-QFT-SOURCE-HISTORY-STATE-PRESERVING-RECORD-NET |
| [Lean/QFT/SourceLinkSquare.lean](../Lean/QFT/SourceLinkSquare.lean) | 26 | No claim evidence edge; model mapping open |
| [Lean/QFT/SourceOperatorGeneration.lean](../Lean/QFT/SourceOperatorGeneration.lean) | 14 | OPH-FINITE-LOCALITY-NOSIGNALLING, OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT/SourceOrientedCompletion.lean](../Lean/QFT/SourceOrientedCompletion.lean) | 16 | OPH-FINITE-BORN-FRAME-RANK-GAP |
| [Lean/QFT/SourcePhaseLiftBridge.lean](../Lean/QFT/SourcePhaseLiftBridge.lean) | 16 | OPH-FINITE-BORN-FRAME-RANK-GAP |
| [Lean/QFT/SourceRecurrenceClock.lean](../Lean/QFT/SourceRecurrenceClock.lean) | 40 | OPH-QFT-SOURCE-RECURRENCE-CLOCK-PRECURSOR, OPH-EM-CARRIER-MODE-OSCILLATORS, OPH-GEOMETRY-SOURCE-CLOCK-RATE-ALONG-WORLDLINES, OPH-GEOMETRY-PROPER-LENGTH-CLOCKED-CHAIN |
| [Lean/QFT/SourceRegionalNet.lean](../Lean/QFT/SourceRegionalNet.lean) | 69 | OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT/StructuralInheritance.lean](../Lean/QFT/StructuralInheritance.lean) | 28 | OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX |
| [Lean/QFT/StructuralNetAdequacySurface.lean](../Lean/QFT/StructuralNetAdequacySurface.lean) | 22 | OPH-FINITE-LOCALITY-NOSIGNALLING, OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT/TowerAnchoredCorrelation.lean](../Lean/QFT/TowerAnchoredCorrelation.lean) | 13 | OPH-QUANTUM-TIME-SLICE-INTERFACE-DELIMITATION |
| [Lean/QFT/TowerAnchoredDiamond.lean](../Lean/QFT/TowerAnchoredDiamond.lean) | 9 | OPH-FINITE-LOCALITY-NOSIGNALLING, OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT/TripleCarrierJoin.lean](../Lean/QFT/TripleCarrierJoin.lean) | 44 | OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX, OPH-QUANTUM-TRIPLE-CARRIER-JOIN |
| [Lean/QFT/TripleCarrierOperatorJoin.lean](../Lean/QFT/TripleCarrierOperatorJoin.lean) | 23 | OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX, OPH-QUANTUM-TRIPLE-OPERATOR-JOIN |
| [Lean/QFT/TwoFactorHistoryBinding.lean](../Lean/QFT/TwoFactorHistoryBinding.lean) | 26 | OPH-GAUGE-KINETIC-HISTORY-BINDING |
| [Lean/QFT/TwoSlotCPNetWitness.lean](../Lean/QFT/TwoSlotCPNetWitness.lean) | 34 | OPH-FINITE-LOCALITY-NOSIGNALLING, OPH-FINITE-CAUSAL-OBSERVER-NET |
| [Lean/QFT.lean](../Lean/QFT.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Screen/A2HolonomyBridge.lean](../Lean/Screen/A2HolonomyBridge.lean) | 10 | OPH-SCREEN-PORT-CURRENT-INNER, OPH-SCREEN-A5-SM-RECOGNITION |
| [Lean/Screen/A5AngularBands.lean](../Lean/Screen/A5AngularBands.lean) | 6 | No claim evidence edge; model mapping open |
| [Lean/Screen/A5AngularKernels.lean](../Lean/Screen/A5AngularKernels.lean) | 22 | No claim evidence edge; model mapping open |
| [Lean/Screen/A5AngularMultiplets.lean](../Lean/Screen/A5AngularMultiplets.lean) | 6 | No claim evidence edge; model mapping open |
| [Lean/Screen/A5CarrierClassBand.lean](../Lean/Screen/A5CarrierClassBand.lean) | 27 | OPH-A5-CARRIER-CLASS-DISPERSION-BAND |
| [Lean/Screen/A5CharacterField.lean](../Lean/Screen/A5CharacterField.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/Screen/A5Commutant.lean](../Lean/Screen/A5Commutant.lean) | 13 | No claim evidence edge; model mapping open |
| [Lean/Screen/A5CouplingSymmetry.lean](../Lean/Screen/A5CouplingSymmetry.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/Screen/A5FamilyBand.lean](../Lean/Screen/A5FamilyBand.lean) | 21 | OPH-SCREEN-A5-SM-RECOGNITION |
| [Lean/Screen/A5IncidenceResponse.lean](../Lean/Screen/A5IncidenceResponse.lean) | 2 | OPH-SCREEN-PORT-CURRENT-INNER |
| [Lean/Screen/A5OPH.lean](../Lean/Screen/A5OPH.lean) | 65 | No claim evidence edge; model mapping open |
| [Lean/Screen/A5OrbitRaySeparation.lean](../Lean/Screen/A5OrbitRaySeparation.lean) | 12 | OPH-A5-PRIMITIVE-PORT-SPIN6, OPH-A5-SEAM-CURRENT-EDGE30 |
| [Lean/Screen/A5PortAction.lean](../Lean/Screen/A5PortAction.lean) | 9 | OPH-EM-CARRIER-MODE-EQUIVARIANCE |
| [Lean/Screen/A5PortModule.lean](../Lean/Screen/A5PortModule.lean) | 7 | No claim evidence edge; model mapping open |
| [Lean/Screen/A5PortSixAxesBridge.lean](../Lean/Screen/A5PortSixAxesBridge.lean) | 7 | OPH-EM-GOLDEN-SECTOR-CHARACTERS |
| [Lean/Screen/A5PrimitivePortPrediction.lean](../Lean/Screen/A5PrimitivePortPrediction.lean) | 5 | OPH-A5-PRIMITIVE-PORT-SPIN6 |
| [Lean/Screen/A5ResponseWordAlgebra.lean](../Lean/Screen/A5ResponseWordAlgebra.lean) | 67 | No claim evidence edge; model mapping open |
| [Lean/Screen/A5SixAxes.lean](../Lean/Screen/A5SixAxes.lean) | 28 | No claim evidence edge; model mapping open |
| [Lean/Screen/ActionTimeGram.lean](../Lean/Screen/ActionTimeGram.lean) | 3 | OPH-QFT-EFFECTIVE-QUANTUM-COMPARISON, OPH-SCALAR-ACTION-CLOCK |
| [Lean/Screen/AssembledActionComposition.lean](../Lean/Screen/AssembledActionComposition.lean) | 24 | OPH-SCREEN-ASSEMBLED-STRUCTURAL-ACTION, OPH-SCREEN-FERMION-SECTOR-ASSEMBLY, OPH-SCREEN-GAUGE-SECTOR-BRACKET-COMPLETION |
| [Lean/Screen/B10EdgeCenterAction.lean](../Lean/Screen/B10EdgeCenterAction.lean) | 29 | OPH-FINITE-QUANTUM-LIMITATION-SUITE, OPH-FINITE-EXTERIOR-COMPONENT-BRIDGE |
| [Lean/Screen/B10QuantumLimitations.lean](../Lean/Screen/B10QuantumLimitations.lean) | 0 | OPH-FINITE-QUANTUM-LIMITATION-SUITE, OPH-FINITE-EXTERIOR-COMPONENT-BRIDGE |
| [Lean/Screen/BaryonDimensionSix.lean](../Lean/Screen/BaryonDimensionSix.lean) | 20 | OPH-SM-BARYON-DIMENSION-SIX-CENSUS |
| [Lean/Screen/BipoSHFrameInvariant.lean](../Lean/Screen/BipoSHFrameInvariant.lean) | 6 | No claim evidence edge; model mapping open |
| [Lean/Screen/BipoSHInverseBoundary.lean](../Lean/Screen/BipoSHInverseBoundary.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/Screen/BipoSHTransferInvariant.lean](../Lean/Screen/BipoSHTransferInvariant.lean) | 2 | No claim evidence edge; model mapping open |
| [Lean/Screen/CarrierEvolutionFlow.lean](../Lean/Screen/CarrierEvolutionFlow.lean) | 78 | OPH-EM-CARRIER-EVOLUTION-FLOW, OPH-EM-FIELD-SECTOR-ENERGY-INNER-PRODUCT, OPH-EM-CURL-SECTOR-EIGENBASIS |
| [Lean/Screen/CarrierFrequencySpeed.lean](../Lean/Screen/CarrierFrequencySpeed.lean) | 13 | OPH-POSITIVE-COSINE-FREQUENCY-CONTRACTION |
| [Lean/Screen/CarrierModeEquivariance.lean](../Lean/Screen/CarrierModeEquivariance.lean) | 81 | OPH-EM-CARRIER-MODE-EQUIVARIANCE, OPH-EM-GOLDEN-SECTOR-CHARACTERS, OPH-EM-GOLDEN-SECTOR-IRREDUCIBILITY |
| [Lean/Screen/CarrierModeOscillators.lean](../Lean/Screen/CarrierModeOscillators.lean) | 81 | OPH-EM-CARRIER-MODE-OSCILLATORS, OPH-EM-CARRIER-MODE-EQUIVARIANCE, OPH-EM-CARRIER-EVOLUTION-FLOW |
| [Lean/Screen/CarrierUniqueness.lean](../Lean/Screen/CarrierUniqueness.lean) | 47 | No claim evidence edge; model mapping open |
| [Lean/Screen/CartanScalarReduction.lean](../Lean/Screen/CartanScalarReduction.lean) | 5 | OPH-SM-Q1-LOCAL-G6 |
| [Lean/Screen/CertifiedScaledStepInstrument.lean](../Lean/Screen/CertifiedScaledStepInstrument.lean) | 15 | OPH-EM-CERTIFIED-SCALED-STEP-INSTRUMENT, OPH-GEOMETRY-COMMON-WORLD-MAXWELL-CLOCK-JOIN, OPH-GEOMETRY-CARRIER-DYNAMICS-COMPATIBILITY |
| [Lean/Screen/ChargedCurrentDictionary.lean](../Lean/Screen/ChargedCurrentDictionary.lean) | 30 | OPH-SCREEN-CHARGED-CURRENT-DICTIONARY |
| [Lean/Screen/CommonEWOrderUnit.lean](../Lean/Screen/CommonEWOrderUnit.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/Screen/CommonHistoryFeedback.lean](../Lean/Screen/CommonHistoryFeedback.lean) | 4 | OPH-COMMON-HISTORY-FEEDBACK-PACKET |
| [Lean/Screen/Compact12.lean](../Lean/Screen/Compact12.lean) | 30 | No claim evidence edge; model mapping open |
| [Lean/Screen/ConeCochainBridge.lean](../Lean/Screen/ConeCochainBridge.lean) | 23 | OPH-CONE-WHITNEY-BRIDGE, OPH-WHITNEY-RADIATIVE-QUANTUM, OPH-WHITNEY-CHARGED-MATTER |
| [Lean/Screen/CurlSectorEigenbasis.lean](../Lean/Screen/CurlSectorEigenbasis.lean) | 60 | OPH-EM-CARRIER-EVOLUTION-FLOW, OPH-EM-FIELD-SECTOR-ENERGY-INNER-PRODUCT, OPH-EM-CURL-SECTOR-EIGENBASIS |
| [Lean/Screen/CurlStoneClockBridge.lean](../Lean/Screen/CurlStoneClockBridge.lean) | 19 | OPH-EM-CURL-SECTOR-EIGENBASIS |
| [Lean/Screen/DiscreteCoulombGreen.lean](../Lean/Screen/DiscreteCoulombGreen.lean) | 79 | OPH-EM-POSITION-SPACE-MAXWELL-COULOMB |
| [Lean/Screen/DiscreteGauss.lean](../Lean/Screen/DiscreteGauss.lean) | 7 | OPH-FINITE-CONSERVATION-WARD-PRECURSOR |
| [Lean/Screen/DiscreteRefinement.lean](../Lean/Screen/DiscreteRefinement.lean) | 18 | No claim evidence edge; model mapping open |
| [Lean/Screen/DispersionArmingInterface.lean](../Lean/Screen/DispersionArmingInterface.lean) | 21 | OPH-EM-DISPERSION-ARMING-INTERFACE |
| [Lean/Screen/ElectroweakBreakingComposition.lean](../Lean/Screen/ElectroweakBreakingComposition.lean) | 33 | OPH-SCREEN-ELECTROWEAK-BREAKING-COMPOSITION, OPH-SCREEN-ASSEMBLED-STRUCTURAL-ACTION, OPH-SCREEN-FERMION-SECTOR-ASSEMBLY |
| [Lean/Screen/ElectroweakPoleScaleQuotient.lean](../Lean/Screen/ElectroweakPoleScaleQuotient.lean) | 5 | OPH-WZ-STRICT-1L-POLE-MAP |
| [Lean/Screen/EqualStateWeights.lean](../Lean/Screen/EqualStateWeights.lean) | 4 | No claim evidence edge; model mapping open |
| [Lean/Screen/ExteriorComponentBridge.lean](../Lean/Screen/ExteriorComponentBridge.lean) | 11 | OPH-FINITE-EXTERIOR-COMPONENT-BRIDGE |
| [Lean/Screen/ExteriorSelection.lean](../Lean/Screen/ExteriorSelection.lean) | 10 | OPH-SCREEN-SUPER-TANNAKIAN-MATTER-LIFT, OPH-GAUGE-D7-D9, OPH-SCREEN-A5-SM-RECOGNITION |
| [Lean/Screen/FermionSectorAssembly.lean](../Lean/Screen/FermionSectorAssembly.lean) | 41 | OPH-SCREEN-FERMION-SECTOR-ASSEMBLY |
| [Lean/Screen/FieldSectorEnergyInnerProduct.lean](../Lean/Screen/FieldSectorEnergyInnerProduct.lean) | 113 | OPH-EM-CARRIER-EVOLUTION-FLOW, OPH-EM-FIELD-SECTOR-ENERGY-INNER-PRODUCT |
| [Lean/Screen/GaugeKineticInvariantForms.lean](../Lean/Screen/GaugeKineticInvariantForms.lean) | 6 | OPH-GAUGE-KINETIC-HISTORY-BINDING, OPH-SCREEN-ASSEMBLED-STRUCTURAL-ACTION, OPH-SCREEN-GAUGE-SECTOR-BRACKET-COMPLETION |
| [Lean/Screen/GaugeOrbitQuotientGap.lean](../Lean/Screen/GaugeOrbitQuotientGap.lean) | 57 | OPH-YM-GAUGE-ORBIT-QUOTIENT-DESCENT |
| [Lean/Screen/GaugeSectorBracketCompletion.lean](../Lean/Screen/GaugeSectorBracketCompletion.lean) | 39 | OPH-SCREEN-GAUGE-SECTOR-BRACKET-COMPLETION |
| [Lean/Screen/GlobalFormCharacterDescent.lean](../Lean/Screen/GlobalFormCharacterDescent.lean) | 15 | OPH-SCREEN-GLOBAL-FORM-CHARACTER-DESCENT |
| [Lean/Screen/GoldenSectorCharacters.lean](../Lean/Screen/GoldenSectorCharacters.lean) | 51 | OPH-EM-GOLDEN-SECTOR-CHARACTERS, OPH-EM-GOLDEN-SECTOR-IRREDUCIBILITY, OPH-EM-GOLDEN-SECTOR-COMPLEX-IRREDUCIBILITY |
| [Lean/Screen/GoldenSectorComplexIrreducibility.lean](../Lean/Screen/GoldenSectorComplexIrreducibility.lean) | 71 | OPH-EM-GOLDEN-SECTOR-CHARACTERS, OPH-EM-GOLDEN-SECTOR-IRREDUCIBILITY, OPH-EM-GOLDEN-SECTOR-COMPLEX-IRREDUCIBILITY |
| [Lean/Screen/GoldenSectorIrreducibility.lean](../Lean/Screen/GoldenSectorIrreducibility.lean) | 55 | OPH-EM-GOLDEN-SECTOR-CHARACTERS, OPH-EM-GOLDEN-SECTOR-IRREDUCIBILITY, OPH-EM-GOLDEN-SECTOR-COMPLEX-IRREDUCIBILITY |
| [Lean/Screen/GoldenSectorPSL2F5Representations.lean](../Lean/Screen/GoldenSectorPSL2F5Representations.lean) | 12 | OPH-EM-GOLDEN-SECTOR-CHARACTERS |
| [Lean/Screen/HolonomyInterference.lean](../Lean/Screen/HolonomyInterference.lean) | 26 | OPH-FINITE-HOLONOMY-CHARACTER-PHASE |
| [Lean/Screen/KineticFamilyCancellation.lean](../Lean/Screen/KineticFamilyCancellation.lean) | 7 | No claim evidence edge; model mapping open |
| [Lean/Screen/KineticFormDichotomy.lean](../Lean/Screen/KineticFormDichotomy.lean) | 6 | No claim evidence edge; model mapping open |
| [Lean/Screen/KogutSusskindFiberRateComparison.lean](../Lean/Screen/KogutSusskindFiberRateComparison.lean) | 78 | OPH-YM-FIBER-RATE-COMPARISON, OPH-YM-GAUGE-ORBIT-QUOTIENT-DESCENT |
| [Lean/Screen/LabeledEventReadout.lean](../Lean/Screen/LabeledEventReadout.lean) | 6 | OPH-GEOMETRY-SCREEN-CARRIER-MAP-CANDIDATE |
| [Lean/Screen/LayeredDiscreteGauss.lean](../Lean/Screen/LayeredDiscreteGauss.lean) | 15 | OPH-GEOMETRY-LAYERED-DISCRETE-GAUSS |
| [Lean/Screen/LightSignalAdequacySurface.lean](../Lean/Screen/LightSignalAdequacySurface.lean) | 25 | No claim evidence edge; model mapping open |
| [Lean/Screen/LightSignalMaxwellComposition.lean](../Lean/Screen/LightSignalMaxwellComposition.lean) | 14 | No claim evidence edge; model mapping open |
| [Lean/Screen/LocalEnergyBalance.lean](../Lean/Screen/LocalEnergyBalance.lean) | 19 | OPH-EM-LOCAL-ENERGY-BALANCE |
| [Lean/Screen/LocalFaceMaxwellAction.lean](../Lean/Screen/LocalFaceMaxwellAction.lean) | 47 | OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX, OPH-EM-LOCAL-FACE-MAXWELL-STATIC-COMPOSITION |
| [Lean/Screen/LocalGaugeJetAction.lean](../Lean/Screen/LocalGaugeJetAction.lean) | 5 | OPH-SM-Q1-LOCAL-G6 |
| [Lean/Screen/MagneticMatterStability.lean](../Lean/Screen/MagneticMatterStability.lean) | 3 | OPH-WHITNEY-MAGNETIC-CONTINUUM |
| [Lean/Screen/MatterGrammarIndexBridge.lean](../Lean/Screen/MatterGrammarIndexBridge.lean) | 8 | OPH-SM-BARYON-DIMENSION-SIX-CENSUS |
| [Lean/Screen/ModalMaxwellFactorizationBoundary.lean](../Lean/Screen/ModalMaxwellFactorizationBoundary.lean) | 17 | OPH-MODAL-MAXWELL-FACTORIZATION-BOUNDARY |
| [Lean/Screen/NeutralCurrentDictionary.lean](../Lean/Screen/NeutralCurrentDictionary.lean) | 33 | OPH-SCREEN-NEUTRAL-CURRENT-DICTIONARY |
| [Lean/Screen/NeutralPairCoupledAction.lean](../Lean/Screen/NeutralPairCoupledAction.lean) | 28 | OPH-EM-SOURCED-FIELD-CONSISTENCY |
| [Lean/Screen/NeutralPairJointStationaryWitness.lean](../Lean/Screen/NeutralPairJointStationaryWitness.lean) | 50 | OPH-EM-SOURCED-FIELD-CONSISTENCY, OPH-SERIAL-MAXWELL-READOUT |
| [Lean/Screen/OPHScreen.lean](../Lean/Screen/OPHScreen.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Screen/OrientedFaceBracketSelector.lean](../Lean/Screen/OrientedFaceBracketSelector.lean) | 18 | OPH-SCREEN-PORT-CURRENT-INNER |
| [Lean/Screen/OrientedFaceInvariantMetric.lean](../Lean/Screen/OrientedFaceInvariantMetric.lean) | 12 | OPH-SCREEN-PORT-CURRENT-INNER |
| [Lean/Screen/PSL2F5SixAxesBridge.lean](../Lean/Screen/PSL2F5SixAxesBridge.lean) | 46 | OPH-EM-GOLDEN-SECTOR-CHARACTERS |
| [Lean/Screen/PhysicalA5ForcingNoGo.lean](../Lean/Screen/PhysicalA5ForcingNoGo.lean) | 6 | OPH-SCREEN-A5-SM-RECOGNITION |
| [Lean/Screen/PortDualMetricSelection.lean](../Lean/Screen/PortDualMetricSelection.lean) | 6 | OPH-SCREEN-PORT-CURRENT-INNER |
| [Lean/Screen/PortFrameGram.lean](../Lean/Screen/PortFrameGram.lean) | 11 | No claim evidence edge; model mapping open |
| [Lean/Screen/PortGramA5Isometry.lean](../Lean/Screen/PortGramA5Isometry.lean) | 18 | OPH-GR-PORT-RESPONSE-COMPLETION |
| [Lean/Screen/PortGramRepairBand.lean](../Lean/Screen/PortGramRepairBand.lean) | 20 | OPH-GR-PORT-RESPONSE-COMPLETION, OPH-SCREEN-GALOIS-FRAME-BOUNDARY |
| [Lean/Screen/PortGramRepairCovariance.lean](../Lean/Screen/PortGramRepairCovariance.lean) | 27 | OPH-GR-PORT-RESPONSE-COMPLETION |
| [Lean/Screen/PositionSpaceMaxwellAction.lean](../Lean/Screen/PositionSpaceMaxwellAction.lean) | 38 | OPH-QFT-STRUCTURAL-INHERITANCE-MATRIX, OPH-EM-POSITION-SPACE-MAXWELL-COULOMB, OPH-EM-SOURCED-FIELD-CONSISTENCY |
| [Lean/Screen/PrimitiveHopSelection.lean](../Lean/Screen/PrimitiveHopSelection.lean) | 4 | No claim evidence edge; model mapping open |
| [Lean/Screen/PrimitivePortDualMeasure.lean](../Lean/Screen/PrimitivePortDualMeasure.lean) | 4 | No claim evidence edge; model mapping open |
| [Lean/Screen/PrimitivePortFrameQuotient.lean](../Lean/Screen/PrimitivePortFrameQuotient.lean) | 37 | OPH-GR-PORT-RESPONSE-COMPLETION, OPH-GEOMETRY-INVERSE-SQUARE-SHELL |
| [Lean/Screen/PrimitivePortMetricAttachment.lean](../Lean/Screen/PrimitivePortMetricAttachment.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/Screen/PrimitivePortOperatorSelectionBoundary.lean](../Lean/Screen/PrimitivePortOperatorSelectionBoundary.lean) | 13 | No claim evidence edge; model mapping open |
| [Lean/Screen/PrimitivePortScaleBoundary.lean](../Lean/Screen/PrimitivePortScaleBoundary.lean) | 12 | No claim evidence edge; model mapping open |
| [Lean/Screen/PrimitivePortTranslationBridge.lean](../Lean/Screen/PrimitivePortTranslationBridge.lean) | 21 | No claim evidence edge; model mapping open |
| [Lean/Screen/QuantumMatterIntegration.lean](../Lean/Screen/QuantumMatterIntegration.lean) | 11 | OPH-FINITE-EXTERIOR-COMPONENT-BRIDGE |
| [Lean/Screen/QuarkMeanDomain.lean](../Lean/Screen/QuarkMeanDomain.lean) | 3 | No claim evidence edge; model mapping open |
| [Lean/Screen/RGRepresentationFrontier.lean](../Lean/Screen/RGRepresentationFrontier.lean) | 5 | No claim evidence edge; model mapping open |
| [Lean/Screen/RegionalContinuity.lean](../Lean/Screen/RegionalContinuity.lean) | 7 | OPH-FINITE-CONSERVATION-WARD-PRECURSOR |
| [Lean/Screen/RepairWordCarrierReadout.lean](../Lean/Screen/RepairWordCarrierReadout.lean) | 20 | OPH-GR-PORT-RESPONSE-COMPLETION |
| [Lean/Screen/S2DesignSignature.lean](../Lean/Screen/S2DesignSignature.lean) | 3 | No claim evidence edge; model mapping open |
| [Lean/Screen/SL2F5PortCover.lean](../Lean/Screen/SL2F5PortCover.lean) | 7 | No claim evidence edge; model mapping open |
| [Lean/Screen/SMStructureAdequacySurface.lean](../Lean/Screen/SMStructureAdequacySurface.lean) | 28 | OPH-SCREEN-GLOBAL-FORM-CHARACTER-DESCENT |
| [Lean/Screen/SMStructureComposition.lean](../Lean/Screen/SMStructureComposition.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/Screen/ScaledMaxwellStability.lean](../Lean/Screen/ScaledMaxwellStability.lean) | 145 | OPH-EM-SCALED-MAXWELL-STABILITY-ACTION, OPH-EM-CERTIFIED-SCALED-STEP-INSTRUMENT, OPH-GEOMETRY-COMMON-WORLD-JOINT-ACTION, OPH-EM-LOCAL-ENERGY-BALANCE, OPH-GEOMETRY-PORT-CHARGE-MINIMAL-COUPLING, OPH-EM-CARRIER-MODE-OSCILLATORS, OPH-EM-CARRIER-MODE-EQUIVARIANCE, OPH-GEOMETRY-SEAM-STEP-SPEED-LIMIT, OPH-EM-GOLDEN-SECTOR-CHARACTERS, OPH-EM-CARRIER-EVOLUTION-FLOW, OPH-EM-SOURCED-FIELD-CONSISTENCY, OPH-EM-FIELD-SECTOR-ENERGY-INNER-PRODUCT, OPH-EM-CURL-SECTOR-EIGENBASIS |
| [Lean/Screen/SeamChargeContinuity.lean](../Lean/Screen/SeamChargeContinuity.lean) | 47 | OPH-EM-SOURCED-FIELD-CONSISTENCY |
| [Lean/Screen/SeamCurrentAuxiliaryOscillatorLift.lean](../Lean/Screen/SeamCurrentAuxiliaryOscillatorLift.lean) | 7 | OPH-A5-SEAM-CURRENT-EDGE30 |
| [Lean/Screen/SeamCurrentCarrierQuotient.lean](../Lean/Screen/SeamCurrentCarrierQuotient.lean) | 30 | OPH-GR-PORT-RESPONSE-COMPLETION, OPH-A5-SEAM-CURRENT-EDGE30 |
| [Lean/Screen/SeamCurrentDirichletGenerator.lean](../Lean/Screen/SeamCurrentDirichletGenerator.lean) | 34 | OPH-A5-SEAM-CURRENT-EDGE30 |
| [Lean/Screen/SeamCurrentEdge30Moment.lean](../Lean/Screen/SeamCurrentEdge30Moment.lean) | 18 | OPH-A5-SEAM-CURRENT-EDGE30 |
| [Lean/Screen/SeamCurrentEdge30Remainder.lean](../Lean/Screen/SeamCurrentEdge30Remainder.lean) | 1 | OPH-A5-SEAM-CURRENT-EDGE30 |
| [Lean/Screen/SeamCurrentFreePhotonLift.lean](../Lean/Screen/SeamCurrentFreePhotonLift.lean) | 21 | OPH-A5-SEAM-CURRENT-EDGE30, OPH-MODAL-MAXWELL-FACTORIZATION-BOUNDARY |
| [Lean/Screen/SeamCurrentHomogeneousAction.lean](../Lean/Screen/SeamCurrentHomogeneousAction.lean) | 37 | OPH-A5-SEAM-CURRENT-EDGE30 |
| [Lean/Screen/SeamCurrentPhotonLeptonThreshold.lean](../Lean/Screen/SeamCurrentPhotonLeptonThreshold.lean) | 57 | OPH-A5-SEAM-CURRENT-EDGE30 |
| [Lean/Screen/SeamCurrentPhysicalMetricAttachment.lean](../Lean/Screen/SeamCurrentPhysicalMetricAttachment.lean) | 13 | No claim evidence edge; model mapping open |
| [Lean/Screen/SeamMaxwellContinuum.lean](../Lean/Screen/SeamMaxwellContinuum.lean) | 6 | OPH-SEAM-MAXWELL-CONTINUUM |
| [Lean/Screen/SeamU1HolonomyClassification.lean](../Lean/Screen/SeamU1HolonomyClassification.lean) | 31 | No claim evidence edge; model mapping open |
| [Lean/Screen/SerialMaxwellReadout.lean](../Lean/Screen/SerialMaxwellReadout.lean) | 11 | OPH-SERIAL-MAXWELL-READOUT |
| [Lean/Screen/SourceActionTime.lean](../Lean/Screen/SourceActionTime.lean) | 5 | OPH-QFT-EFFECTIVE-QUANTUM-COMPARISON, OPH-SCALAR-ACTION-CLOCK |
| [Lean/Screen/TemporalMaxwellEvolution.lean](../Lean/Screen/TemporalMaxwellEvolution.lean) | 30 | OPH-EM-TEMPORAL-MAXWELL-EVOLUTION-CONTINUITY, OPH-EM-SCALED-MAXWELL-STABILITY-ACTION, OPH-EM-LOCAL-ENERGY-BALANCE |
| [Lean/Screen/TopThreeKernelFix.lean](../Lean/Screen/TopThreeKernelFix.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/Screen/TraceBalancedKernel.lean](../Lean/Screen/TraceBalancedKernel.lean) | 10 | OPH-SCREEN-GLOBAL-FORM-CHARACTER-DESCENT |
| [Lean/Screen/TrichotomyCases.lean](../Lean/Screen/TrichotomyCases.lean) | 6 | No claim evidence edge; model mapping open |
| [Lean/Screen/UnitSplit12.lean](../Lean/Screen/UnitSplit12.lean) | 4 | No claim evidence edge; model mapping open |
| [Lean/Screen/VolumeReadoutBridge.lean](../Lean/Screen/VolumeReadoutBridge.lean) | 15 | No claim evidence edge; model mapping open |
| [Lean/Screen/WeylYukawaConventions.lean](../Lean/Screen/WeylYukawaConventions.lean) | 16 | OPH-GAUGE-EXTERIOR-SM-WITNESS |
| [Lean/Screen/WhitneyActualCertificates.lean](../Lean/Screen/WhitneyActualCertificates.lean) | 10 | No claim evidence edge; model mapping open |
| [Lean/Screen/WhitneyAlgebraicLDL.lean](../Lean/Screen/WhitneyAlgebraicLDL.lean) | 20 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyCertificateCLI.lean](../Lean/Screen/WhitneyCertificateCLI.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Screen/WhitneyCertificatePipeline.lean](../Lean/Screen/WhitneyCertificatePipeline.lean) | 12 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyCertifiedConsumers.lean](../Lean/Screen/WhitneyCertifiedConsumers.lean) | 17 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyChargedMatter.lean](../Lean/Screen/WhitneyChargedMatter.lean) | 10 | OPH-WHITNEY-CHARGED-MATTER |
| [Lean/Screen/WhitneyConeMass.lean](../Lean/Screen/WhitneyConeMass.lean) | 9 | OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyConeMassNaturality.lean](../Lean/Screen/WhitneyConeMassNaturality.lean) | 8 | OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyConeModes.lean](../Lean/Screen/WhitneyConeModes.lean) | 18 | OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyConeNaturality.lean](../Lean/Screen/WhitneyConeNaturality.lean) | 10 | No claim evidence edge; model mapping open |
| [Lean/Screen/WhitneyCurlRow0.lean](../Lean/Screen/WhitneyCurlRow0.lean) | 3 | No claim evidence edge; model mapping open |
| [Lean/Screen/WhitneyCurlRow1.lean](../Lean/Screen/WhitneyCurlRow1.lean) | 3 | No claim evidence edge; model mapping open |
| [Lean/Screen/WhitneyCurlRow2.lean](../Lean/Screen/WhitneyCurlRow2.lean) | 3 | No claim evidence edge; model mapping open |
| [Lean/Screen/WhitneyCurlRow3.lean](../Lean/Screen/WhitneyCurlRow3.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/Screen/WhitneyExactSourceProblem.lean](../Lean/Screen/WhitneyExactSourceProblem.lean) | 3 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyFiniteCertificate.lean](../Lean/Screen/WhitneyFiniteCertificate.lean) | 14 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyFrameMorphism.lean](../Lean/Screen/WhitneyFrameMorphism.lean) | 9 | OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyFrameNaturality.lean](../Lean/Screen/WhitneyFrameNaturality.lean) | 32 | OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyGeneratedCertificate.lean](../Lean/Screen/WhitneyGeneratedCertificate.lean) | 12 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyMassPremiseInstance.lean](../Lean/Screen/WhitneyMassPremiseInstance.lean) | 5 | OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyMaxwellDynamics.lean](../Lean/Screen/WhitneyMaxwellDynamics.lean) | 31 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyNormalModeConstruction.lean](../Lean/Screen/WhitneyNormalModeConstruction.lean) | 9 | OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyOmittedCellCounterexample.lean](../Lean/Screen/WhitneyOmittedCellCounterexample.lean) | 6 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyQuantumBridge.lean](../Lean/Screen/WhitneyQuantumBridge.lean) | 24 | OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneyRotationWitness.lean](../Lean/Screen/WhitneyRotationWitness.lean) | 6 | OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneySourceAssembly.lean](../Lean/Screen/WhitneySourceAssembly.lean) | 21 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneySourceGeometry.lean](../Lean/Screen/WhitneySourceGeometry.lean) | 20 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneySourceNaturality.lean](../Lean/Screen/WhitneySourceNaturality.lean) | 30 | OPH-WHITNEY-MAXWELL-DYNAMICS, OPH-WHITNEY-RADIATIVE-QUANTUM |
| [Lean/Screen/WhitneySourceProblemCLI.lean](../Lean/Screen/WhitneySourceProblemCLI.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Screen/WhitneySpatialConsistency.lean](../Lean/Screen/WhitneySpatialConsistency.lean) | 5 | OPH-WHITNEY-SPATIAL-CONSISTENCY |
| [Lean/Screen/WhitneyTimeBridge.lean](../Lean/Screen/WhitneyTimeBridge.lean) | 15 | OPH-CONE-WHITNEY-BRIDGE |
| [Lean/Screen/Z6Descent.lean](../Lean/Screen/Z6Descent.lean) | 15 | OPH-GAUGE-D7-D9, OPH-SCREEN-A5-SM-RECOGNITION, OPH-SCREEN-GLOBAL-FORM-CHARACTER-DESCENT |
| [Lean/Screen/Z6Exact.lean](../Lean/Screen/Z6Exact.lean) | 5 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/ContinuumConnection.lean](../Lean/SourceLaw/ContinuumConnection.lean) | 7 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/ContinuumTimeLimit.lean](../Lean/SourceLaw/ContinuumTimeLimit.lean) | 11 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/CounterfamilyLocality.lean](../Lean/SourceLaw/CounterfamilyLocality.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/CounterfamilyReader.lean](../Lean/SourceLaw/CounterfamilyReader.lean) | 14 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/DynamicalResponse.lean](../Lean/SourceLaw/DynamicalResponse.lean) | 14 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/GoldenCarrier.lean](../Lean/SourceLaw/GoldenCarrier.lean) | 61 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/GoldenFamily.lean](../Lean/SourceLaw/GoldenFamily.lean) | 18 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/HistoryRegister.lean](../Lean/SourceLaw/HistoryRegister.lean) | 7 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/IdentificationLadder.lean](../Lean/SourceLaw/IdentificationLadder.lean) | 6 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/LangevinStationarity.lean](../Lean/SourceLaw/LangevinStationarity.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/MatrixLimits.lean](../Lean/SourceLaw/MatrixLimits.lean) | 31 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/MatrixSeams.lean](../Lean/SourceLaw/MatrixSeams.lean) | 5 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/OpenWindow.lean](../Lean/SourceLaw/OpenWindow.lean) | 24 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/ReflectionPositivity.lean](../Lean/SourceLaw/ReflectionPositivity.lean) | 7 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/RepairAdmissibility.lean](../Lean/SourceLaw/RepairAdmissibility.lean) | 21 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/RepairConvergence.lean](../Lean/SourceLaw/RepairConvergence.lean) | 10 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/ScalarPacket.lean](../Lean/SourceLaw/ScalarPacket.lean) | 33 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/StationaryWindow.lean](../Lean/SourceLaw/StationaryWindow.lean) | 25 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/TransferMatrixLift.lean](../Lean/SourceLaw/TransferMatrixLift.lean) | 42 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/TransferStructure.lean](../Lean/SourceLaw/TransferStructure.lean) | 19 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/TwoPointFlux.lean](../Lean/SourceLaw/TwoPointFlux.lean) | 3 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/VacuumIdentity.lean](../Lean/SourceLaw/VacuumIdentity.lean) | 18 | No claim evidence edge; model mapping open |
| [Lean/SourceLaw/WeightedSeams.lean](../Lean/SourceLaw/WeightedSeams.lean) | 8 | No claim evidence edge; model mapping open |
| [Lean/Thermodynamics/CapFirstLaw.lean](../Lean/Thermodynamics/CapFirstLaw.lean) | 9 | OPH-THERMO-FOUR-LAW-PACKAGE, OPH-GEOMETRY-INTERNAL-ENERGY-INERTIA-PRECURSOR, OPH-THERMO-GIBBS-REFERENCE-ENERGY-IDENTIFICATION |
| [Lean/Thermodynamics/CofinalSpectralTailFamily.lean](../Lean/Thermodynamics/CofinalSpectralTailFamily.lean) | 26 | OPH-THERMO-COFINAL-SPECTRAL-TAIL-COMPOSITION |
| [Lean/Thermodynamics/CoherentRefinementFamily.lean](../Lean/Thermodynamics/CoherentRefinementFamily.lean) | 26 | OPH-THERMO-COHERENT-REFINEMENT-ENVELOPE |
| [Lean/Thermodynamics/CollarTemperatureReading.lean](../Lean/Thermodynamics/CollarTemperatureReading.lean) | 23 | OPH-THERMO-COLLAR-TEMPERATURE-READING |
| [Lean/Thermodynamics/CommonObjectBinding.lean](../Lean/Thermodynamics/CommonObjectBinding.lean) | 46 | OPH-THERMO-FOUR-LAW-PACKAGE |
| [Lean/Thermodynamics/CommonReferenceObstruction.lean](../Lean/Thermodynamics/CommonReferenceObstruction.lean) | 10 | OPH-THERMO-FOUR-LAW-PACKAGE, OPH-THERMO-SOURCE-REALIZATION-AUDIT |
| [Lean/Thermodynamics/EinsteinPremiseLink.lean](../Lean/Thermodynamics/EinsteinPremiseLink.lean) | 7 | OPH-THERMO-FOUR-LAW-PACKAGE, OPH-GR-EINSTEIN-BRANCH-INHABITANT |
| [Lean/Thermodynamics/EquationOfState.lean](../Lean/Thermodynamics/EquationOfState.lean) | 17 | OPH-FINITE-WORK-RESPONSE-NONIDENTIFIABILITY, OPH-GOLDEN-SCALAR-METRIC-PRESSURE, OPH-WHITNEY-MAXWELL-METRIC-PRESSURE |
| [Lean/Thermodynamics/FiniteConditionalRepair.lean](../Lean/Thermodynamics/FiniteConditionalRepair.lean) | 35 | OPH-THERMO-FOUR-LAW-PACKAGE, OPH-THERMO-GIBBS-REFERENCE-ENERGY-IDENTIFICATION |
| [Lean/Thermodynamics/FirstLawIdentity.lean](../Lean/Thermodynamics/FirstLawIdentity.lean) | 1 | OPH-THERMO-FOUR-LAW-PACKAGE |
| [Lean/Thermodynamics/FluctuationTheorems.lean](../Lean/Thermodynamics/FluctuationTheorems.lean) | 9 | OPH-THERMO-FOUR-LAW-PACKAGE |
| [Lean/Thermodynamics/FourLawAdequacySurface.lean](../Lean/Thermodynamics/FourLawAdequacySurface.lean) | 28 | OPH-THERMO-FOUR-LAW-PACKAGE, OPH-THERMO-GIBBS-REFERENCE-ENERGY-IDENTIFICATION |
| [Lean/Thermodynamics/GibbsEntropyLimit.lean](../Lean/Thermodynamics/GibbsEntropyLimit.lean) | 9 | No claim evidence edge; model mapping open |
| [Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean](../Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean) | 36 | OPH-THERMO-GIBBS-REFERENCE-ENERGY-IDENTIFICATION, OPH-THERMO-MODULAR-ENERGY-ADDITIVITY, OPH-THERMO-COLLAR-TEMPERATURE-READING |
| [Lean/Thermodynamics/GraphDiffusion.lean](../Lean/Thermodynamics/GraphDiffusion.lean) | 42 | OPH-FINITE-TRANSPORT-GREEN-KUBO-DIFFUSION |
| [Lean/Thermodynamics/GreenKubo.lean](../Lean/Thermodynamics/GreenKubo.lean) | 48 | OPH-FINITE-TRANSPORT-GREEN-KUBO-DIFFUSION |
| [Lean/Thermodynamics/HorizonThermalitySurface.lean](../Lean/Thermodynamics/HorizonThermalitySurface.lean) | 10 | OPH-THERMO-FOUR-LAW-PACKAGE, OPH-GR-EINSTEIN-BRANCH-INHABITANT |
| [Lean/Thermodynamics/KolmogorovCriterion.lean](../Lean/Thermodynamics/KolmogorovCriterion.lean) | 18 | OPH-SCREEN-MODULAR-GEARING |
| [Lean/Thermodynamics/LowTemperatureControl.lean](../Lean/Thermodynamics/LowTemperatureControl.lean) | 36 | OPH-THERMO-FOUR-LAW-PACKAGE |
| [Lean/Thermodynamics/MixingChainRealization.lean](../Lean/Thermodynamics/MixingChainRealization.lean) | 16 | OPH-THERMO-FOUR-LAW-PACKAGE, OPH-QFT-SOURCE-RECURRENCE-CLOCK-PRECURSOR |
| [Lean/Thermodynamics/ModularEnergyAdditivity.lean](../Lean/Thermodynamics/ModularEnergyAdditivity.lean) | 33 | OPH-THERMO-MODULAR-ENERGY-ADDITIVITY |
| [Lean/Thermodynamics/NativeRepairFourCycle.lean](../Lean/Thermodynamics/NativeRepairFourCycle.lean) | 11 | OPH-NATIVE-REPAIR-HISTORY-SOURCE |
| [Lean/Thermodynamics/PhysicalCalibrationImport.lean](../Lean/Thermodynamics/PhysicalCalibrationImport.lean) | 41 | OPH-CAL-SI-ANCHOR-IMPORT, OPH-QFT-SOURCE-CLOCK-CANDIDATE |
| [Lean/Thermodynamics/PoissonizedRepair.lean](../Lean/Thermodynamics/PoissonizedRepair.lean) | 10 | OPH-PUBLICIZATION-DYNAMICS-FINITE |
| [Lean/Thermodynamics/PoissonizedRepairOperatorExp.lean](../Lean/Thermodynamics/PoissonizedRepairOperatorExp.lean) | 8 | OPH-PUBLICIZATION-DYNAMICS-FINITE |
| [Lean/Thermodynamics/RepairCurrentOrientation.lean](../Lean/Thermodynamics/RepairCurrentOrientation.lean) | 20 | OPH-FINITE-BORN-FRAME-RANK-GAP |
| [Lean/Thermodynamics/StationaryRealization.lean](../Lean/Thermodynamics/StationaryRealization.lean) | 8 | OPH-THERMO-FOUR-LAW-PACKAGE, OPH-THERMO-SOURCE-REALIZATION-AUDIT |
| [Lean/Thermodynamics/TemporalFilterMonotonicity.lean](../Lean/Thermodynamics/TemporalFilterMonotonicity.lean) | 1 | No claim evidence edge; model mapping open |
| [Lean/Time/ClockComparison.lean](../Lean/Time/ClockComparison.lean) | 11 | OPH-BOUNDED-OBSERVER-TIME-CALIBRATION |
| [Lean/Time/ClockReadout.lean](../Lean/Time/ClockReadout.lean) | 13 | OPH-BOUNDED-OBSERVER-TIME-CALIBRATION |
| [Lean/Time/ObserverHistory.lean](../Lean/Time/ObserverHistory.lean) | 5 | OPH-BOUNDED-OBSERVER-TIME-CALIBRATION |
| [Lean/Time/ProperTimeCalibration.lean](../Lean/Time/ProperTimeCalibration.lean) | 5 | OPH-BOUNDED-OBSERVER-TIME-CALIBRATION |
| [Lean/Time/SourceCountClock.lean](../Lean/Time/SourceCountClock.lean) | 5 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM |
| [Lean/Time/SourceCountClockEnclosure.lean](../Lean/Time/SourceCountClockEnclosure.lean) | 5 | OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM |
| [Lean/Time/TimeOrderLedger.lean](../Lean/Time/TimeOrderLedger.lean) | 4 | OPH-TIME-ORDER-TYPE-LEDGER |
| [Lean/Time/WorldlineRealization.lean](../Lean/Time/WorldlineRealization.lean) | 8 | OPH-BOUNDED-OBSERVER-TIME-CALIBRATION |
| [Lean/Tower/AdaptiveFixedPointEndpoint.lean](../Lean/Tower/AdaptiveFixedPointEndpoint.lean) | 3 | OPH-CONS-ADAPTIVE-PUBLIC-ENDPOINT |
| [Lean/Tower/ConsensusTower.lean](../Lean/Tower/ConsensusTower.lean) | 7 | OPH-CONSENSUS-TOWER-INTERFACE |
| [Lean/Tower/CumulativeCapacityEndpoint.lean](../Lean/Tower/CumulativeCapacityEndpoint.lean) | 2 | OPH-CONS-CAPACITY-BOUNDED-PUBLIC-ENDPOINT |
| [Lean/Tower/EventGeometryReadout.lean](../Lean/Tower/EventGeometryReadout.lean) | 18 | No claim evidence edge; model mapping open |
| [Lean/Tower/FixedFederationEndpoint.lean](../Lean/Tower/FixedFederationEndpoint.lean) | 2 | OPH-COMPUTATION-FIXED-FEDERATION-ENDPOINT |
| [Lean/Tower/FixedFederationExecutionAudit.lean](../Lean/Tower/FixedFederationExecutionAudit.lean) | 0 | OPH-COMPUTATION-FIXED-FEDERATION-EXECUTION-CLASSIFICATION |
| [Lean/Tower/FixedFederationExecutionEndpoint.lean](../Lean/Tower/FixedFederationExecutionEndpoint.lean) | 5 | OPH-COMPUTATION-FIXED-FEDERATION-EXECUTION-CLASSIFICATION |
| [Lean/Tower/FixedFederationFanoutEndpoint.lean](../Lean/Tower/FixedFederationFanoutEndpoint.lean) | 3 | OPH-COMPUTATION-FIXED-FEDERATION-EXECUTION-CLASSIFICATION |
| [Lean/Tower/FixedPointEndpoint.lean](../Lean/Tower/FixedPointEndpoint.lean) | 62 | OPH-CONS-D1, OPH-PUBLIC-WORLD-ENDPOINT-FINITE |
| [Lean/Tower/FullPremiseInhabitant.lean](../Lean/Tower/FullPremiseInhabitant.lean) | 16 | No claim evidence edge; model mapping open |
| [Lean/Tower/OperationalObserver.lean](../Lean/Tower/OperationalObserver.lean) | 26 | OPH-OBS-PATCH-DEF |
| [Lean/Tower/PublicWorldQuotient.lean](../Lean/Tower/PublicWorldQuotient.lean) | 5 | OPH-CONS-D1, OPH-PUBLIC-WORLD-ENDPOINT-FINITE |
| [Lean/Tower.lean](../Lean/Tower.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/Variational/DiscreteEulerLagrange.lean](../Lean/Variational/DiscreteEulerLagrange.lean) | 3 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/Variational/DiscreteNoether.lean](../Lean/Variational/DiscreteNoether.lean) | 6 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/Variational/EnrichmentCharacterization.lean](../Lean/Variational/EnrichmentCharacterization.lean) | 49 | OPH-MECH-ENRICHMENT-CHARACTERIZATION |
| [Lean/Variational/FiniteHistoryBridge.lean](../Lean/Variational/FiniteHistoryBridge.lean) | 3 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/Variational/FiniteRealTransfer.lean](../Lean/Variational/FiniteRealTransfer.lean) | 11 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/Variational/LegendreBridge.lean](../Lean/Variational/LegendreBridge.lean) | 44 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/Variational/MechanicsAdequacySurface.lean](../Lean/Variational/MechanicsAdequacySurface.lean) | 23 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/Variational/ModeExtremalEnrichment.lean](../Lean/Variational/ModeExtremalEnrichment.lean) | 22 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/Variational/PolynomialEnrichmentSelection.lean](../Lean/Variational/PolynomialEnrichmentSelection.lean) | 42 | OPH-MECH-POLYNOMIAL-BOOST-SELECTION |
| [Lean/Variational/RealizedHistoryLegendreNoGo.lean](../Lean/Variational/RealizedHistoryLegendreNoGo.lean) | 21 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS, OPH-GEOMETRY-INTERNAL-ENERGY-INERTIA-PRECURSOR |
| [Lean/Variational/SourceToHamiltonianComposed.lean](../Lean/Variational/SourceToHamiltonianComposed.lean) | 13 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/Variational/StationarySaddleCoverage.lean](../Lean/Variational/StationarySaddleCoverage.lean) | 3 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/Variational/TranslationInvariantComposedInstance.lean](../Lean/Variational/TranslationInvariantComposedInstance.lean) | 17 | OPH-FINITE-HISTORY-VARIATIONAL-HELPERS |
| [Lean/Variational.lean](../Lean/Variational.lean) | 0 | No claim evidence edge; model mapping open |
| [Lean/lakefile.lean](../Lean/lakefile.lean) | 0 | No claim evidence edge; model mapping open |
