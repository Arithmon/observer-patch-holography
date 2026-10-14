# Observer Patch Holography

**Physics from the records observers can share.**

[Main paper](https://philpapers.org/rec/MUEFOC) · [Executable observer](model/README.md) · [Lean library](Lean/) · [Reproduce](REPRODUCE.md) · [Pragma Research](https://floatingpragma.io/research/) · [Français](README_FR.md)

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
  one disturbed state, compiles both quantum-native repairs, and preserves
  an encoded qubit through record-selected recovery under declared services.
- **Standard physics from finite agreement.** Finite theorems take the
  observer's readouts to familiar structure: public records that do not depend
  on the order of repairs; the Born rule, which fixes quantum probabilities,
  and the Tsirelson bound on completed records; a carrier that is
  three-dimensional on its own; the Standard Model gauge Lie algebra with one
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
retained results. Its strongest joined experiment is the
[common interacting process](model/common_source/README.md), in which quantum
electric flux, charged fermions and a scalar field share one finite register
with complete instruments, retained records and native feedback; the field
action and the preparation are supplied inputs. The
[prediction package](model/predictions/README.md) fixes two competing detector
instruments and their count-based rejection rule before any independent
counts, and the [construction guide](model/CONSTRUCTION.md) states the
acceptance criteria and the connections that remain to be built.

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
  geometry, the geometric building block of an observer's space.
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
| [claims/](claims/README.md) | The claim registry, the premise registry and the frozen-prediction register |
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
[candidate register](model/CANDIDATES.md) and the
[research issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
set out that work. Paper releases follow the
[release policy](model/CONSTRUCTION.md#research-programme-and-paper-release-hold);
main-branch PDFs are source previews, and deposited releases retain their
original bytes.

<details>
<summary>Related projects and funding</summary>

[Cadence](https://github.com/muellerberndt/cadence) is the related learning-machine
project. Its capability evidence is separate from OPH's physics evidence.
[Pragma Research](https://floatingpragma.io/) brings the projects together.

OPH's official TEO (Theory Of Everything) memecoin is on Base at
`0xd2a7055317D0c7b316319cDadaE592D0644A0b07`
([chart](https://dexscreener.com/base/0xba55cdae026b9dd63fcb50e82094d61472948e44527ca4791033497469228049)).
Bernhard Mueller can claim its trading fees, which help fund the research.
The token has no role in the theory. Other tokens using the OPH or TEO name
are unrelated. Nothing here is financial advice.

</details>

## License

Software, including Lean and the model code, uses [Apache-2.0](code/LICENSE).
Papers, documentation, figures, data and ledgers use [CC BY-NC-SA 4.0](LICENSE),
with upstream Particle Data Group terms for the packaged particle values.
The [license map](LICENSE) gives the directory-specific terms.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->
