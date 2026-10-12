# Observer Patch Holography

**Physics from the records observers can share.**

[Runnable model](model/README.md) · [Physics preprint](https://philpapers.org/rec/MUEFOC) · [Pragma Research](https://floatingpragma.io/research/) · [Français](README_FR.md)

Observer Patch Holography (OPH) investigates whether familiar physics can be
reconstructed from bounded observers that reach agreement. An observer patch
has local state, a boundary with ports, readback, records, and feedback for
repairing disagreement. Within the model, facts become public when records
survive comparison across overlapping patches.

The proposed universe is one structure viewed both as a timeless consistent
history and as the internal local computation of that history. Its records
and responses to interventions supply the operation; no external host or
preferred global clock is required. Space, physical duration and matter are
intended emergent readouts. Current constructions supply some geometry,
states, actions and clock maps as premises. The
[interpretation and working terms](docs/INTERPRETATION_AND_TERMS.md) distinguish
this physical proposal from the conditional results achieved in these constructions.

This repository contains the scientific work behind that idea: papers,
machine-checked proofs and executable models. It connects finite agreement
with quantum probability, spatial geometry and the symmetry pattern of the
known forces, with a conditional route to Einstein dynamics. Each construction
states the assumptions it consumes. Establishing a common physical
realization requires measurements that distinguish the proposed models.

## Run the observer model

**[model/](model/README.md) is the live home of observer and repair-law
construction.** The new [common interacting process](model/common_source/README.md)
couples quantum electric flux, charged fermions and a scalar on one evolving
register. A matter observation controls feedback; a later scalar observation
controls native repair. All 1,056 complete native/quantum histories retain
their records and terminal states. An exact finite certificate proves that
the feedback changes the scalar response, including ideal temporal splitting
error. Its [observer readouts](model/common_source/OBSERVER_READOUTS.md)
construct a GNS Hilbert representation from the actual state, execute
a readable finite modular-clock instrument, and identify interactions
between its supplied regions. The action, access, services and pulse
schedule remain declared inputs. Physical clock calibration, an interacting
spatial continuum and empirical comparison remain open.

The [history repair formulation](model/common_source/HISTORY_REPAIR.md)
recovers the same completed observations under every fair local update
schedule for these fixed inputs. Observer event order and modular readouts
remain distinct from the auxiliary solver counter.
The [native geometric readout](model/common_source/SLOW_GEOMETRY.md)
connects the inherited rank-three response carrier to recorded currents;
an additional binary controller executes displacements spanning three
directions without identifying them with physical position or momentum.

A separate [finite-resolution detector](model/common_source/SMOOTH_REFINEMENT.md)
has a quantitative error bound through its complete interacting observation
history and an explicit measurement-energy ledger. Its mathematical
occupation limit retains a fixed two-site apparatus; spatial refinement and
agreement between its clock, distance and energy readouts require further construction.

The [fixed-reference clock bridge](model/common_source/THERMAL_CLOCK.md)
connects the same Hamiltonian to modular evolution by choosing a faithful
equilibrium reference distinct from the changing observed state. It does
not establish that repair selects that reference or calibrate physical time.
The [native propagation tests](model/common_source/CAUSALITY.md) check that
information follows local seam interactions under fixed schedules. Units
may be chosen by observers; consistent relations among readouts are the target.

Its finite observer executes local repair, retains an
observation and uses that record to choose its next interaction. The
P-dependent detector computes the electromagnetic coupling and lets a
retained binary outcome control the next repair.
The [scalar detector](model/observer/scalar_response/README.md) uses the native
current to prepare a scalar field, reads one binary outcome and feeds it back
into native repair, with a controlled comparison across the existing field
refinements. Its field law, vacuum, apparatus and model clock are explicit
inputs. The local evolution extension executes the field's neighboring-site
updates at q=13 and q=34, includes the original vacuum's squeezing in the
detector probability, and passes its sampled bit to native feedback. A
separate temporal certificate extends the inherited continuum error budget;
the two finest grids remain analytic comparisons for this local extension.
The new [source-compiled field apparatus](model/observer/source_scalar/README.md)
executes a separate two-site design from blank registers through entangled
preparation, quantum evolution, pointer readout and native feedback. All
field gates are emitted from the existing twelve-port controls. Exact
certificates bound its complete ideal occupation history and finite temporal
splitting. The executed sixteen-level design has an ideal detector probability
between 0.51495 and 0.51789, excluding a fair response under the declared
apparatus. Spatial refinement, numerical execution and physical service
errors require separate bounds for that same construction.

The [specified prediction tests](model/predictions/README.md) distinguish
the current-record and contrast-record quantum instruments with exact phase
predictions, a fixed sample plan and a validator for complete independent
counts. The experimental verdict is `no_data`. A separate retrospective
audit retains the P detector's incomplete low-energy coupling contract and
its approximately 300 ppm inverse-alpha discrepancy; no fitted correction
or empirical confirmation is introduced.

Run from the repository root with Python 3.11 or newer:

```sh
python3 -m model.observer run --trials 10000 --seed 42
python3 -m model.observer.p_response --check model/observer/p_response/results.json
```

The [model README](model/README.md) gives the exact checks, reproducible
receipts and construction guide. These are completed finite software
protocols with explicit P/N premises and apparatus choices. Broader physical
identification remains part of the research below.

## Research status

We have **runnable finite observers, conditional mathematical results and
concrete repair candidates**. The native observer and P-dependent binary
detector execute complete observation, record and feedback protocols.
The common interacting construction places scalar, charged matter,
electric evolution and record-controlled interventions on one quantum work
space, with operator Gauss constraints and a certified finite response.
Its theorem attachments distinguish record creation from persistence and
exact hardware refinement from approximate occupation cutoffs. It does not
yet instantiate every compatible finite theorem or establish the required
continuum and gravitational descriptions of observer statistics.
The scalar adapter connects those interfaces to a supplied free-field
continuum family and an independently checkable conditional detector
prediction. Its local evolution supplies a checked numerical path from
neighboring field updates to the detector-controlled repair, with a temporal
refinement control. It does not derive that Hamiltonian from integer averaging
or establish a new empirical prediction of OPH.
The source-compiled two-site experiment implements preparation,
evolution and detection through the same source-control interface, including
both terminal quantum branches. Its coherent transport services, field
action and clock remain explicit inputs. This closes a finite implementation
gap without establishing a common continuum realization of observed physics.
The candidates reproduce and connect some required structures, and we can
distinguish several of them in simulations. **We do not have a
common model established as physically correct.** The
[candidate register](model/CANDIDATES.md) and
[model search](model/CONSTRUCTION.md) separate constructed results, supplied
assumptions and missing physical connections.

The [central-port construction](code/native_response/README.md#a-central-port-carrier-that-retains-the-quantum-repair)
extends the native/N4 source with twelve classical addresses while preserving
its decoded quantum channel. Its abstract source pair holds the full quantum
reference and normalized port geometry fixed while exhibiting different
finite kinetics. Occupied-port readout has a separate quantum information
cost. The [four audit controls](model/CONSTRUCTION.md#retained-audit-controls)
retain the precise carrier, marked-action, feedback and inertia boundaries.

**OPH's physical correctness is a hypothesis under investigation.** Logical
consistency, conditional theorems and successful simulations do not establish
that nature follows OPH or any particular candidate repair law. We retain
alternative laws and revise or reject physical attachments when tests fail.

The finite constructions supply model comparisons. A physical test additionally
needs an operational identification of states, records, time and parameters;
a distinguishing observable; an account of calibration and uncertainty; and
a prediction and rejection rule fixed before independent comparison data are
examined. Compare against established physics and alternative models, including
cases where they cannot be distinguished. The
[physical-test programme](https://github.com/FloatingPragma/observer-patch-holography/issues/1025)
and [frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md) track that work.
Reproducing a value used to choose the model is a fit, not an independent test.

The parallel [programme-wide obstruction search](https://github.com/FloatingPragma/observer-patch-holography/issues/1078)
seeks contradictions, joint-realization impossibility, unavoidable independent
physical inputs and empirical exclusions. Each result states which fixed
formulation or model class it rejects. A decisive negative result is valuable;
failing to find one does not confirm OPH. The issue also distinguishes scientific
refutation from a practical decision to stop pursuing the programme.

## Research focus: one repair process

**Find the simplest shared repair rule that satisfies the known constraints.**

**The completion target is one common construction with a controlled
continuum limit reproducing established physics and quantitative predictions
that survive independent tests.** Its finite observer implementation is a
completed component. Construction choices need no uniqueness theorem;
their physical consequences need compatible readouts, explicit assumptions
and discriminating evidence.

More precisely, the target is an executable observer algebra whose compatible
measurement statistics admit nontrivial 3+1-dimensional quantum fields,
their observed interactions, classical physics and gravity as effective
descriptions in the appropriate regimes. Fundamental fields or a fundamental
spacetime are not prerequisites; their observational content is the target.
Shared operational meaning defines the proposed objective structure; its
agreement with nature remains a testable requirement. The
[completion criteria](model/CONSTRUCTION.md#completion-target-and-acceptance-criteria)
and [objectivity contract](docs/INTERPRETATION_AND_TERMS.md#observer-measurements-and-operational-objectivity)
state this target, including gapless sectors and controlled errors, separately
from the current finite achievements.

The research priority is an explicit microscopic repair law whose evolution
supplies the inputs of the existing theorems on one common model. Separate
conditional results can be correct without their physical premises being
realized together. We start with the simplest compatible construction, reuse
proved results, and test each added hypothesis against alternatives. The
carrier geometry, reversible gauge response and local operational clocks must be
compatible readouts of that same process. Add a mechanism only for a
demonstrated missing connection, and count every added physical assumption.
A rejected scalar-only candidate does not establish that a simple common
source law is impossible.

Use the [established agreement mathematics](docs/AGREEMENT_RESEARCH_FOUNDATIONS.md)
as the research foundation: specify recovery and scheduling assumptions,
map linear overlap models to sheaf diffusion and boundary extension, and reuse
the resulting theorems where their premises hold. Each candidate must identify
its remaining physical obligation and a test that can distinguish its law.

Compare the [concrete candidate laws](model/CANDIDATES.md), their tested
properties and remaining choices. Follow the [model search](model/CONSTRUCTION.md), its
[model–theorem premise map](model/PREMISES.md), and the
[research issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics).
Progress means a tested connection or a rejected candidate; theorem counts
do not measure physical completion.

Paper releases and broad prose restructuring are on hold until a coherent
model milestone supports them. Audits, correctness fixes and reproducibility
checks continue. The [research policy](model/CONSTRUCTION.md#research-programme-and-paper-release-hold)
states the milestone and the correction exception.

Main-branch PDFs are current-source previews and may retain the existing
release identifier after corrections. Tagged or deposited release artifacts
retain their original bytes; a preview is not a new release.

## Start here

| To explore… | Start with… |
| --- | --- |
| Run the observer and continue its construction | [Runnable model](model/README.md) |
| The construction contracts and open physical connections | [Construction guide](model/CONSTRUCTION.md) |
| The active repair laws, assumptions and distinguishing tests | [Candidate register](model/CANDIDATES.md) |
| The existing physics argument | [*Finite Observer Consensus as a Reconstruction Principle*](https://philpapers.org/rec/MUEFOC) |
| The proofs and reproducible calculations | [Lean library](Lean/) and [reproduction guide](REPRODUCE.md) |
| The related work on learning machines | [Cadence](https://floatingpragma.io/cadence/) and [interactive demos](https://floatingpragma.io/demos/) |
| Explanations for a general audience | [Pragma Research blog](https://blog.floatingpragma.io/) |

## Mathematical foundations

The research program connects these parts of physics through one observer
architecture, with separate assumptions stated for each result.

- **Agreement between observers.** Under the stated termination and
  consistency conditions, local repairs reach a protected public record
  independent of their order. The finite theorems describe agreement,
  stability and refinement.
- **Quantum probability.** Finite event algebras support the Born probability
  rule, Lüders conditioning and exact Tsirelson results, including a Bell-state
  construction that attains 2√2. These are mathematical results on the stated
  event-algebra branch.
- **Geometry and the known forces.** A twelve-port response construction
  supplies a three-dimensional spatial readback. Complete reversible response
  and internal observer transport force the Standard Model gauge Lie type
  under the stated hypotheses. A supplied matter representation gives an
  anomaly-free fifteen-state generation; its physical realization and the
  global gauge group require additional structure.
- **Spacetime and gravity.** Declared source and reading laws support
  controlled causal-geometry constructions. The Einstein-equation implication
  assumes a common physical realization with the stated stress, entropy,
  continuum and scale data.
- **Classical and quantum matter.** A supplied charged-scalar/Maxwell action
  supports controlled nonlinear continuum motion in its real sector and an
  interacting quantum state space on a fixed mesh. Selecting the action from
  observer histories and identifying it with physical matter are separate
  requirements.
- **Constants as fixed-point problems.** Koide's relation holds exactly under
  a stated balance premise. The fine-structure calculation certifies a root
  of a declared closure map and carries diagnostic status. The capacity
  program asks whether the public capacity assigned to the universe agrees
  with the capacity reconstructed from within it. Connecting these
  constructions to measured constants requires physical identification.

The [main preprint](https://philpapers.org/rec/MUEFOC) develops these connections
and their precise scope. The [paper index](paper/) provides the specialist
accounts, including thermodynamics, field geometry and fixed-point
constructions for the constants.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## Evidence you can inspect

The [Lean library](Lean/) contains more than 12300 public theorems and lemmas
with no admitted proofs. Exact certificates and reproducible simulations
accompany the mathematical arguments. The companion
[physics simulator](https://github.com/muellerberndt/oph-physics-sim) provides
executable observer dynamics and retained evidence.

The [axiom reference](docs/AXIOM_REFERENCE.md) states the three core axioms;
the [premise register](docs/registers/PREMISE_REGISTER_V3.md) and
[model–theorem premise map](model/PREMISES.md) record additional premises and
result-specific assumptions, including physical identifications and empirical
inputs.
The [postdiction ledger](docs/POSTDICTION_LEDGER.md) records comparisons with
measured values and their input ancestry; the
[frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md) records tests whose
conditions must be fixed before the comparison data are examined.
The [falsification program](docs/OPH_FALSIFICATION_PROGRAM.md) gives the
observations that would refute particular claims.

### Reproduce the finite core

After setting up the dependencies in [REPRODUCE.md](REPRODUCE.md), check the
claim graph and selected finite algebra, record-capacity and consensus results:

```bash
python3 tools/check_claim_registry.py
python3 -m pytest -q \
  code/a5_closure/test_audit.py \
  code/capacity_readback/test_correctable_public_record_capacity.py \
  code/capacity_readback/test_reversible_public_checkpoint_packet.py \
  code/consensus/test_reference_architecture_benchmark_suite.py \
  code/consensus/test_verified_tree_packet_net.py
```

The reproduction guide describes the wider checks and the scope of each
evidence family.

## Reconstruction map

<p align="center">
  <a href="assets/prediction-chain.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg" alt="OPH reconstruction chain" width="92%">
  </a>
</p>

<p align="center"><sub>The OPH reconstruction map connects observer records, three-dimensional source geometry, causal order, clocks, fields and quantum states. Each arrow names a mathematical connection; the papers state the assumptions that join them into an effective physical description.</sub></p>

## Repository guide

| Path | Contents |
| --- | --- |
| [`flagship/`](flagship/) | Main standalone physics paper, TeX source and PDF |
| [`paper/`](paper/) | Core papers and publication index |
| [`Lean/`](Lean/) | Machine-checked mathematical development |
| [`code/`](code/) and [`evidence/`](evidence/) | Executable models, certificates and reproduction evidence |
| [`extra/`](extra/) and [`cosmology/`](cosmology/) | Focused mathematical and physical research papers |
| [`book/`](book/) | *Reverse Engineering Reality*, source and downloadable book |
| [`docs/`](docs/) | Canonical references: axioms, frozen predictions, postdictions and falsifiers |
| [`docs/registers/`](docs/registers/) | Generated premise, observation, instrument and selection registers |
| [`docs/research/`](docs/research/) | Markdown research notes, each verified by a code lane |
| [`docs/policies/`](docs/policies/) | Data, hardware-evidence and writing policies |

## Contribute

Prefer a simpler shared rule over an accumulation of separate mechanisms.
Contributions should repair an existing defect or resolve a named uncertainty
in the common model and its tests. Reuse proved results; add a proof when a
specific model obligation needs it. Standalone finite-theorem expansion and
theorem-count targets are outside this programme. Start with the
[research issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
and [reproduction guide](REPRODUCE.md). Keep exploratory logs local and report
concise decisions with their assumptions and evidence. Verified research notes
belong in [`docs/research/`](docs/research/); paper publication follows the hold above.

## From observer patches to learning machines

[Cadence](https://github.com/muellerberndt/cadence) makes the settling idea
executable as a learning architecture. Its bounded software patches carry
local state and readback, with feedback that repairs prediction errors.
Its default brain, System 1, learns from experience as animal brains do and
not by backpropagation. Optional System 2 observers add recursive feedback
inside the same settlement.

The [Cadence preprint](https://philpapers.org/rec/MUECAP-2) describes the
architecture and experiments. Its computational case rests on learning
behavior and measured resource use. [Pragma Research](https://floatingpragma.io/)
connects this work to embodied AI.

## TEO memecoin

OPH has an official memecoin on Base, TEO (Theory Of Everything). Its token
address is `0xd2a7055317D0c7b316319cDadaE592D0644A0b07`
([chart on Dexscreener](https://dexscreener.com/base/0xba55cdae026b9dd63fcb50e82094d61472948e44527ca4791033497469228049)).
Any other token using the OPH or TEO name is unrelated to this project.

The token was set up so that Bernhard Mueller can claim all the trading fees it
generates, and those fees have been a very helpful source of funding for OPH.
The token has no role in the OPH research.

Nothing here is financial advice.

## License

The repository uses split licensing. All software, including the Lean library, [`code/`](code), [`model/`](model), [`tools/`](tools), and the schemas in [`schemas/`](schemas), is licensed under [Apache-2.0](code/LICENSE). Papers, the book, documentation, figures, data, the generated ledgers in [`tracking/`](tracking), and the packaged particle data in [`pdg_data/`](pdg_data) are licensed under [CC BY-NC-SA 4.0](LICENSE); the tabulated values in `pdg_data/` carry their upstream Particle Data Group terms. The [LICENSE](LICENSE) file gives the per-directory map.
