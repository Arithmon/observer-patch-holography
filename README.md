# Observer Patch Holography

**Physics from the records observers can share.**

[Main paper](https://philpapers.org/rec/MUEFOC) · [Book](https://oph-book.floatingpragma.io/) · [Lean library](Lean/) · [Reproduce](REPRODUCE.md) · [Pragma Research](https://floatingpragma.io/) · [Français](README_FR.md)

[![Mandatory Suite](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/mandatory-suite.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/mandatory-suite.yml)
[![Lean CI](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean-ci.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean-ci.yml)
[![Runnable observer model](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/observer-model.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/observer-model.yml)

Observer Patch Holography (OPH) proposes that physics is what bounded
observers can establish together. Finite observers compare records on shared
boundaries, repair their disagreements, and treat as real whatever survives
every comparison. Three problems that resist a classical account motivate the
programme: quantum measurement, correlations that no local classical model
reproduces, and the origin of an apparent spacetime. OPH seeks one explanation
for all three in the consistency of observers' records and operations, and it
pursues that explanation as an executable observer, machine-checked
mathematics and predictions fixed in advance.

**The Universe reads itself into existence.** OPH has one hypothesis, the
closure hypothesis: the Universe is a timeless, self-referential structure
whose observers recover the specification of their own world from their
records and build hardware that realizes that same history, including its
builders. Nothing lies outside this loop, so the recovered specification owes
no further justification. It has to run and reproduce what observers measure.
Space, time and matter are effective descriptions of the records, and "our
Universe" is the public record, what survives every comparison between
observers. The [basis](docs/BASIS.md) states the hypothesis, the no-outside
rule, the agreement condition and the executable specification; the
[working terms](docs/INTERPRETATION_AND_TERMS.md) explain the timeless view
and the clocks inside it.

## Highlights

- **An observer that runs.** The smallest observer has twelve ports and two
  events. It writes a record at each event, keeps it, and lets only that
  record choose the next repair, and the record changes what is measured: the
  final occupation moves from 1/6 to 35/264. One finite quantum process joins
  electric flux, charged matter and a scalar field on a single register, and
  an exact calculation certifies that its recorded feedback changes the
  measured response. Its continued smooth protocol executes four reads on
  one disturbed state, compiles two quantum record laws of the native repair, and preserves
  an encoded qubit through record-selected recovery under declared services.
  Two observers whose regions overlap read their fields and clocks from one
  continuing history and agree on one complete law. Joined carriers form
  screens of up to 245,760 ports, and on a refining screen one repair history
  defines a continuum field with Gaussian fluctuations and waves of finite
  speed.
- **Relativistic field theory as the observers' effective description.** On
  a declared radial product over that screen, the native wave action has a
  Lorentz-signature density, and for every finite set of the observers'
  clock, field and response readouts one smooth static (1+3)-dimensional
  Lorentzian field theory reproduces them all to any chosen tolerance. The
  [integrated verdict](model/CONSTRUCTION.md#integrated-effective-description)
  states the limits: the exact native metric is not smooth in its present
  atlas, and the radial design is supplied. The carrier's own current
  records do derive a radius: their three-dimensional counter has the
  standard radial diffusion law from every initial state. Complete
  measurement histories cannot tell radial designs apart, so tying that
  radius to the field's radial direction and to physical distance is open.
- **Standard physics from finite agreement.** Finite theorems take the
  observer's readouts to familiar structure: public records that do not depend
  on the order of repairs; the Born rule, which fixes quantum probabilities,
  and the Tsirelson bound on completed records; a carrier that is
  three-dimensional in its slow response; the Standard Model gauge Lie algebra with one
  anomaly-free generation; a (1+3)-dimensional causal order in which proper
  time comes from record counts; and the Einstein equation from a realized
  null balance.
- **Mathematics you can check.** The [Lean library](Lean/) contains more than
  12600 public theorems and lemmas with no admitted proofs, and every retained
  computational receipt is replayed under continuous integration.
- **Tests fixed in advance.** Predictions and their rejection rules enter the
  [frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md) before any
  comparison data, with timestamps attested on the Bitcoin blockchain; the
  dispersion prediction for the primitive port orbit is one of them. The
  [postdiction ledger](docs/POSTDICTION_LEDGER.md) and the
  [falsification programme](docs/OPH_FALSIFICATION_PROGRAM.md) keep comparisons
  with measurements apart from tests of the mathematics.

**Agreement with the full range of observed measurement statistics has not
been established.** The theorems hold under stated premises, and the
executable observer is a candidate judged by declared acceptance criteria. The
sections below say which premises the construction supplies and which remain
inputs.

## Run the observer

From the repository root, with Python 3.11 or newer:

```sh
python3 -m model.observer run --trials 10000 --seed 42
```

The run prints one line: the Brier risk of three forecasts of the next record
and the hash of the result.

```text
{"brier": {"action_known_prior": 0.0728, "native_record": 0.0563, "within_action_shuffle": 0.0884}, "reference_matches": 10000, ...}
```

Lower is better. The forecast that uses the retained record beats the one
that knows only the chosen action and the one with records shuffled within the
same action; the exact values are 395/7392, 140/2013 and 12875/150304. The
[model directory](model/README.md) holds the implementation, its tests and its
retained results. Its joined experiments include the
[coherent screen observer](model/coherent_screen_observer/README.md), which
carries a logical qubit and an untouched reference through native clock ticks
and record-controlled feedback on screens of 15,360 and 245,760 ports, and the
[common interacting process](model/common_source/README.md), in which quantum
electric flux, charged fermions and a scalar field share one finite register
with complete instruments, retained records and native feedback; the field
action and the preparation of that process are supplied inputs. The
[repair-law guide](model/CANDIDATES.md) shows which parts of the repair law
are forced and which are free: the classical rule follows from four
requirements, its quantum record laws form one family with five coherence
multipliers, a single interference experiment reads the coherence that
separates the executed law C from its comparison law E, and nine compiled
measurements identify the whole table. The
[prediction package](model/predictions/README.md) fixes that experiment and its
count-based rejection rule before any independent counts, and the
[construction guide](model/CONSTRUCTION.md) states the acceptance criteria and
the connections that remain to be built.

One twelve-port carrier connects spatial and clock observations: current
records recover local geometry, and exact temporal correlation laws reconstruct
a detector's cyclic Hilbert space and its visible transfer modes. Repair fixes a
current-fluctuation action and a probability law relating incident hits to
moving ticks. Retained executions test that conversion against a host clock.
On larger source screens, field and feedback predictions also work at fixed
conserved occupation, without independently prepared Bernoulli occupations.
Two observers sharing a port of one carrier read their fields, local clocks
and feedback from one continuing history, and the direct native execution and
an effective description assign every joint record the same probability; a
version built from compiled source pulses reproduces that history on all 512
code addresses. On the refining spherical screens one repair history gives a
[continuum field operator](model/physical_space_time/OCCUPATION_FIELD.md)
whose Gaussian fields, decay and waves live on the two-dimensional screen
itself. The
[public-data audit](docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md) compares
the conditional scalar tilt with the ACT DR6 posterior chains and an ACT
likelihood slice, and tests a
direct attachment of native clock counts to published optical-clock ratios.
The [recovery assessment](model/OBSERVER_RELATIVE_RECOVERY.md) separates these
finite results from extended relativistic spacetime and a calibrated laboratory
identification. A [public Bell-data replay](model/predictions/BELL_RECORDS.md)
also tests the limit of interpreting the diagonal branch as a complete local
response model; it does not test the full quantum construction. The
[empirical status](docs/EMPIRICAL_STATUS.md) distinguishes those comparisons
from selection of the native repair law. Local history repair reproduces the
supplied oscillator readout; the retained data do not determine the native
record law.

## Mathematical results

- **Agreement and stable records.** Under explicit termination and consistency
  conditions, local repairs reach one shared public result whatever their
  order. Further results control which records survive refinement.
  [Consensus paper](paper/reality_as_consensus_protocol.pdf).
- **Quantum observations.** Within the stated finite quantum algebra, verified
  results cover repeatable measurements, classical records that later
  operations cannot disturb, uncertainty and the Bell correlation bound, with
  a state that attains it.
  [Finite quantum event algebras](extra/machine_checked_finite_event_algebras.pdf).
- **Three spatial directions.** The fixed twelve-port carrier and its mean
  repair rule give an exact three-dimensional slow response with Euclidean
  geometry, which an observer's current records recover.
  [Carrier construction](code/a5_closure/README.md).
- **Gauge symmetry.** Under stated premises on complete reversible response
  and observer transport, the construction forces the Standard Model gauge
  Lie algebra, and a specified matter representation gives an anomaly-free
  generation.
  [Gauge derivation](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf).
- **Thermodynamics.** One finite repair construction gives the four laws,
  fluctuation relations and the Landauer erasure bound under stated reference,
  energy-calibration and low-temperature premises.
  [Thermodynamics account](paper/observers_are_all_you_need.pdf).
- **Information capacity.** Exact results count which records stay
  distinguishable after the allowed operations, bound that count by the
  available carrier, and follow it through refinement.
  [Record capacity](code/capacity_readback/README.md).
- **Matter and fields.** A supplied charged-scalar and Maxwell action with its
  geometry and quantization give an interacting Hamiltonian after gauge
  reduction, and controlled scalar refinement connects classical detector
  readouts to a nonlinear field equation.
  [Matter continuum](paper/observers_are_all_you_need.pdf).

These are conditional mathematical results. Making their premises hold for
the observations of one common model is the central construction problem.
The [premise map](model/PREMISES.md) records which premises the executable
specification supplies and which it takes as inputs.

## Papers

- **[Observers Who Read Their Universe into Existence](flagship/from_observer_consensus_to_standard_physics.pdf)**
  is the main synthesis: the theory in one account, selected results in
  detail, the scoreboard of the construction and the closure hypothesis. Cite
  its released version, the [PhilPapers record](https://philpapers.org/rec/MUEFOC).
- **[An Executable Observer for Observer Patch Holography](paper/executable_observer_construction.pdf)**
  states the specification, its acceptance criteria, its certified components
  and the premise map.
- **The theorem library** holds the long proofs, from
  [consensus](paper/reality_as_consensus_protocol.pdf) and
  [observers and records](paper/observers_are_all_you_need.pdf) to
  [spacetime and Einstein dynamics](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf).
  The [paper index](paper/README.md) and the [paper set](docs/PAPER_SET.md)
  list every paper and its role.
- **Companion papers** carry focused results. A lepton-mass relation, the
  Koide relation, holds exactly under a stated balance premise. The capacity
  programme asks whether the public capacity assigned to the Universe agrees
  with the one reconstructed from inside. These constant calculations,
  including the fine-structure root, have diagnostic status.

## Where things are

| Path | What it holds |
| --- | --- |
| [model/](model/README.md) | The executable observer: implementation, tests, retained receipts and protocol records |
| [flagship/](flagship/) | The synthesis paper |
| [paper/](paper/README.md) | The construction paper and the theorem library, with archived papers under [paper/archive/](paper/archive/README.md) |
| [extra/](extra/README.md), [cosmology/](cosmology/) | Companion papers |
| [Lean/](Lean/) | The machine-checked library |
| [code/](code/README.md) | Exact calculations, lane verifiers and receipts behind the papers |
| [evidence/](evidence/README.md) | Retained evidence packets and custody records |
| [claims/](claims/README.md) | Established claims and their dependency views |
| [tracking/](tracking/) | Premise, observation, prediction and instrument registers |
| [docs/](docs/README.md) | The basis, the working terms, policies and registers |
| [book/](book/), [essays/](essays/) | The book and the essays |

## Reproduce and verify

The [reproduction guide](REPRODUCE.md) runs the verification suite from a
clean clone. The model verifier replays every retained receipt:

```sh
python -m pip install -r model/requirements.txt
python -m model.verify
```

The [Lean library](Lean/) builds with Lake, the [documentation index](docs/README.md)
lists the registers and policies, and the [book](book/) gives a longer road
in. The textbooks at [learn.floatingpragma.io](https://learn.floatingpragma.io/)
and the interactive simulation at
[simulation.floatingpragma.io](https://simulation.floatingpragma.io/) are
the gentlest introductions.

## Contributing

The task is to make one observation-and-repair process supply the effective
physics: connect the observer's clock, spatial interpretation and field
measurements on the same histories, then compare outcome probabilities,
correlations and responses to interventions with independent observations.
A supplied field equation remains an input even when its simulation works.
Contributions should resolve a named connection or expose a failure in the
shared model; the [construction guide](model/CONSTRUCTION.md), the
[repair-law guide](model/CANDIDATES.md) and the
[research issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
set out that work. Paper releases follow the
[release policy](model/CONSTRUCTION.md#research-programme-and-paper-release-hold);
main-branch PDFs are source previews, and deposited releases retain their
original bytes.

## License

Software, including Lean and the model code, uses [Apache-2.0](code/LICENSE).
Papers, documentation, figures, data and ledgers use [CC BY-NC-SA 4.0](LICENSE),
with upstream Particle Data Group terms for the packaged particle values.
The [license map](LICENSE) gives the directory-specific terms.

<details>
<summary>Related projects and funding</summary>

[Cadence](https://github.com/muellerberndt/cadence) is the related learning-machine
project. Its capability evidence is separate from OPH's physics evidence.
[Pragma Research](https://floatingpragma.io/) brings the projects together.

OPH's official TEO (Theory Of Everything) memecoin was deployed by the Base
community and lives on Base at
`0xd2a7055317D0c7b316319cDadaE592D0644A0b07`
([chart](https://dexscreener.com/base/0xba55cdae026b9dd63fcb50e82094d61472948e44527ca4791033497469228049)).
Its trading fees are routed to Bernhard Mueller to fund the project.
The token has no role in the theory. Other tokens using the OPH or TEO name
are unrelated. Nothing here is financial advice.

</details>

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->
