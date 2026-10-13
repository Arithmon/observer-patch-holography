# Observer Patch Holography

**Physics from the records observers can share.**

[Main paper](https://philpapers.org/rec/MUEFOC) · [Executable observer](model/README.md) · [Lean library](Lean/) · [Pragma Research](https://floatingpragma.io/research/) · [Français](README_FR.md)

Observer Patch Holography (OPH) proposes that physics is what bounded
observers can establish together. Finite observers compare records on shared
boundaries, repair their disagreements, and treat as real whatever survives
every comparison. The project builds that picture as an executable observer,
proves what follows from it in machine-checked mathematics, and tests the
results against predictions fixed in advance.

**The Universe reads itself into existence.** OPH has one hypothesis: the
Universe is a timeless, self-referential structure whose observers recover
the specification of their own world from their records and build hardware
that realizes that same history, including its builders. Nothing lies outside
this loop, so the recovered specification owes no further justification. It
has to run and reproduce what observers measure. Space, time and matter are
effective descriptions of the records, and "our Universe" is the public
record, what survives every comparison between observers. The
[basis](docs/BASIS.md) states the hypothesis, its rule and the executable
specification; the [working terms](docs/INTERPRETATION_AND_TERMS.md) explain
the timeless view and the clocks inside it.

## Highlights

- **An observer that runs.** A twelve-port carrier with a native repair law
  writes a record at every event, keeps it, and lets only that record control
  the next repair. One finite quantum process joins electric flux, charged
  matter and a scalar field on a single register, and an exact calculation
  certifies that its recorded feedback changes the measured response.
- **Standard physics from finite agreement.** Finite theorems take the
  observer's readouts to familiar structure: public records that do not
  depend on the order of repairs, the Born rule and the Tsirelson bound on
  completed records, an intrinsically three-dimensional carrier, the Standard
  Model gauge Lie algebra with an anomaly-free generation, a (1+3)-dimensional
  causal precursor with proper time from record counts, and the Einstein
  equation from a realized null balance.
- **Mathematics you can check.** The [Lean library](Lean/) contains more than
  12300 public theorems and lemmas with no admitted proofs, and the
  computational results replay from retained receipts.
- **Tests fixed in advance.** Predictions and their rejection rules are
  registered before comparison data. The
  [postdiction ledger](docs/POSTDICTION_LEDGER.md) and the
  [falsification programme](docs/OPH_FALSIFICATION_PROGRAM.md) keep
  comparisons with measurements apart from tests of the mathematics.

**Agreement with the full range of observed measurement statistics has not
been established.** The theorems hold under stated premises, and the
executable observer is a candidate judged by declared acceptance criteria.
The sections below say which premises the construction supplies and which
remain inputs.

## The executable observer

The [model directory](model/README.md) holds the implementation, its tests and
its retained results. The observer itself is a short program: twelve ports,
a repair law that writes a coin and a current at every event, retained
records, and feedback that only those records control. Run it from the
repository root with Python 3.11 or newer:

```sh
python3 -m model.observer run --trials 10000 --seed 42
```

Its strongest joined experiment is the
[common interacting process](model/common_source/README.md): quantum electric
flux, charged fermions and a scalar field share one finite register, with
complete instruments, retained records and native feedback. An exact bound
certifies that the feedback changes the measured response, including in the
infinite-occupation limit, and the same histories supply observer Hilbert
spaces, modular clocks and geometric readouts in three directions. The field
action and the preparation are supplied inputs. The
[construction guide](model/CONSTRUCTION.md) states the acceptance criteria and
the connections that remain to be built, and the
[prediction package](model/predictions/README.md) fixes two competing
detector predictions and their rejection test before any independent counts.

## Mathematical results

- **Agreement and stable records.** Under explicit termination and consistency
  conditions, local repairs reach one shared public result whatever their
  order. Further results control which records survive refinement.
  [Consensus paper](paper/reality_as_consensus_protocol.pdf).
- **Quantum observations.** Within the stated finite quantum algebra, verified
  results cover repeatable measurements, protected classical records,
  uncertainty and the Bell correlation bound, with a state that attains it.
  [Finite quantum event algebras](extra/machine_checked_finite_event_algebras.pdf).
- **Three spatial directions.** The fixed twelve-port carrier and its mean
  repair rule give an exact three-dimensional slow response with Euclidean
  geometry, the geometric building block of an observer's space.
  [Carrier construction](code/a5_closure/README.md).
- **Gauge symmetry.** With complete reversible response and observer-transport
  premises, the construction forces the Standard Model gauge Lie algebra, and
  a specified matter representation gives an anomaly-free generation.
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
  detail, the scoreboard of the construction and the closure hypothesis. Its
  released version is on [PhilPapers](https://philpapers.org/rec/MUEFOC).
- **[An Executable Observer](paper/executable_observer_construction.pdf)**
  states the specification, its acceptance criteria, its certified components
  and the premise map.
- **The theorem library** holds the long proofs, from
  [consensus](paper/reality_as_consensus_protocol.pdf) and
  [observers and records](paper/observers_are_all_you_need.pdf) to
  [spacetime and Einstein dynamics](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf).
  The [paper index](paper/README.md) and the [paper set](docs/PAPER_SET.md)
  list every paper and its role.
- **Companion papers** carry focused results: the Koide relation holds
  exactly under a stated balance premise, and the capacity programme asks
  whether the public capacity assigned to the Universe agrees with the one
  reconstructed from inside. These constant calculations, including the
  fine-structure root, have diagnostic status.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## Explore and reproduce

Read the [main paper](https://philpapers.org/rec/MUEFOC) for the whole
account, the [book](book/) for a longer road in, and the
[documentation index](docs/README.md) for the registers and policies. The
[reproduction guide](REPRODUCE.md) replays every retained receipt and runs the
verification suite from a clean clone.

<details>
<summary>Related work, funding and licensing</summary>

[Cadence](https://github.com/muellerberndt/cadence) is the related learning-machine
project. Its capability evidence is separate from OPH's physics evidence.
[Pragma Research](https://floatingpragma.io/) brings the projects together.

OPH's official TEO (Theory Of Everything) memecoin is on Base at
`0xd2a7055317D0c7b316319cDadaE592D0644A0b07`
([chart](https://dexscreener.com/base/0xba55cdae026b9dd63fcb50e82094d61472948e44527ca4791033497469228049)).
Bernhard Mueller can claim its trading fees, which help fund the research.
The token has no role in the theory. Other tokens using the OPH or TEO name
are unrelated. Nothing here is financial advice.

Software, including Lean and the model code, uses [Apache-2.0](code/LICENSE).
Papers, documentation, figures, data and ledgers use [CC BY-NC-SA 4.0](LICENSE),
with upstream Particle Data Group terms for the packaged particle values.
The [license map](LICENSE) gives the directory-specific terms.

</details>

## Research focus

**Make one observation-and-repair process supply the effective physics.**
The task is to connect the observer's clock, spatial interpretation and field
measurements on the same histories, then compare outcome probabilities,
correlations and responses to interventions with independent observations.
A supplied field equation remains an input even when its simulation works,
and construction choices need no uniqueness theorem, only consequences that
survive the tests. Contributions should resolve a named connection or expose
a failure in the shared model. The [construction guide](model/CONSTRUCTION.md),
the [candidate register](model/CANDIDATES.md) and the
[research issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
set out that work.

Paper releases follow the [release policy](model/CONSTRUCTION.md#research-programme-and-paper-release-hold).
Main-branch PDFs are source previews; deposited releases retain their original bytes.
