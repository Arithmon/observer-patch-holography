# Observer Patch Holography

**Physics from the records observers can share.**

[Physics preprint](https://philpapers.org/rec/MUEFOC) · [Pragma Research](https://floatingpragma.io/research/) · [Français](README_FR.md)

Observer Patch Holography (OPH) asks whether bounded observers comparing
records and repairing disagreements can reproduce our observed experiences.
The goal is an executable observer whose measurements and shared
interpretations agree with ours.

**In OPH's proposed picture, reality is fundamentally computational and
observation is primary.** Reality is a timeless, self-referential structure
of observations, records and their relations. Each observer assigns meaning
to its observations, and shared meaning must survive comparisons with other
observers. There is no
fundamental, observer-independent spacetime continuum containing those
experiences. Space, time, matter and classical objects are effective
interpretations of them. The apparent continuum belongs to that description;
"objective reality" means the stable operational structure on which observers
agree.

Observers experience change from within this complete history. Its timeless
self-consistency and its internal local observation-and-repair process are
two descriptions of the same proposed structure. An observer's clock is a
readout within it, with no external time in which the whole history must run.
The [working terms](docs/INTERPRETATION_AND_TERMS.md) explain this view.

**The Universe reads itself into existence.** OPH starts from one
hypothesis: its observers recover the specification of their own world from
their records and build hardware whose execution realizes that same
history, including its builders. Because nothing lies outside this
structure, whatever observers recover and implement needs no further
justification. The carrier, the observer algebra and the repair law are
justified by being recovered, by running and by reproducing what observers
measure. Nothing unrecovered may be assumed, and the public record, what
survives every comparison between observers, is what "our Universe" means.
The executable specification is where the account enters this loop; the
loop itself has no first step. The [basis](docs/BASIS.md) states this
precisely.

This radical observers-first approach addresses foundational problems beyond
the classical picture: quantum measurement, correlations incompatible with
a local classical account, and the origin of an apparent spacetime. OPH
seeks a common explanation in the consistency of observers' records and
operations. Explaining these phenomena through one construction is the
programme's ambition; the results below establish specific parts of it.

The project has checked finite mathematical results and executable observers
that connect several of them in one process. **Agreement with the full range
of observed measurement statistics has not been established.** The work below gives
concrete pieces of that construction and tests where they fit together.

## Executable observer

The [model directory](model/README.md) contains the implementation, tests and
reproducible results. Its strongest joined experiment puts a scalar field,
charged fermions and electric flux on one finite quantum register. An
observation produces a record, that record controls feedback, and a later
observation controls local repair. Every outcome is retained, and an exact
error bound certifies that the feedback changes the measured response.
For the same declared field action, an exact infinite-occupation calculation
bounds the feedback contrast away from zero with the coordinate half-line
detector; the action and preparation remain supplied inputs.

The same histories supply observer Hilbert spaces and readable modular-clock
probes. Native currents give exact conservation and transport relations,
with geometric readouts tied to the twelve-port carrier. A smooth detector
has controlled errors as the two-site occupation cutoff increases. These are
executed and proved consequences of the specified apparatus; its field
Hamiltonian, preparation and some clock choices are supplied inputs.
See the [common process](model/common_source/README.md) for the results and
boundaries, including the distinction between repair relaxation and waves
in an observer's time.

Try the finite observer from the repository root with Python 3.11 or newer:

```sh
python3 -m model.observer run --trials 10000 --seed 42
```

The [prediction package](model/predictions/README.md) gives competing detector
predictions and a rejection test; independent experimental counts are absent.
The model also includes a finite source test that distinguishes two repair
instruments through their retained record laws using specified coherent
preparation and readout.
Simulation agreement with a supplied field law does not establish that OPH
derives the corresponding observation statistics.

## Mathematical results

- **Agreement and stable records.** Under explicit termination and consistency
  conditions, local repairs reach a shared public result independent of their
  order. Further results control record preservation and changes under
  refinement. [Consensus paper](paper/reality_as_consensus_protocol.pdf).
- **Quantum observations.** Within the stated finite quantum algebra, verified
  results cover repeatable measurements, protected classical records,
  uncertainty and the Bell correlation bound, including a state that attains
  it. [Finite quantum event algebras](extra/machine_checked_finite_event_algebras.pdf).
- **Three spatial directions.** The fixed twelve-port carrier and its mean
  repair rule give an exact three-dimensional slow response and its Euclidean
  geometry. This supplies a geometric building block for an observer's spatial
  interpretation. [Carrier construction](code/a5_closure/README.md).
- **Gauge symmetry.** With complete reversible response and observer-transport
  assumptions, the construction forces the Standard Model gauge Lie algebra.
  A specified matter representation gives an anomaly-free generation.
  [Gauge derivation](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf).
- **Thermodynamics.** A common finite repair construction gives the four laws,
  fluctuation relations and the Landauer erasure bound under stated reference,
  energy-calibration and low-temperature assumptions.
  [Thermodynamics account](paper/observers_are_all_you_need.pdf).
- **Information capacity.** Exact results count which records remain
  distinguishable after the allowed operations, bound that count by the
  available carrier, and track it through refinement.
  [Record capacity](code/capacity_readback/README.md).
- **Matter and fields.** A supplied charged-scalar/Maxwell action, geometry and
  quantization give an interacting Hamiltonian after gauge reduction. Controlled
  scalar refinement connects classical detector readouts to a nonlinear field
  equation. [Matter continuum](paper/observers_are_all_you_need.pdf).

These are conditional mathematical results. Making their premises hold for
observations in one common model is the central construction problem. The
[basis](docs/BASIS.md) explains how the theorem premises relate to the
executable specification, and the [model–theorem map](model/PREMISES.md)
records which premises the construction supplies. 

## Papers

- **[Observers Who Read Their Universe into Existence](flagship/from_observer_consensus_to_standard_physics.pdf)**
  is the main synthesis: the theory in one account, selected results in
  detail, the scoreboard of the construction and the closure hypothesis. Its released version is on [PhilPapers](https://philpapers.org/rec/MUEFOC).
- **[An Executable Observer](paper/executable_observer_construction.pdf)**
  states the specification, its acceptance criteria, its certified components
  and the premise map.
- **The theorem library** holds the long proofs, from
  [consensus](paper/reality_as_consensus_protocol.pdf) and
  [observers and records](paper/observers_are_all_you_need.pdf) to
  [spacetime and Einstein dynamics](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf).
  The [paper index](paper/README.md) and the [paper set](docs/PAPER_SET.md)
  list every paper and its role. The Koide relation holds exactly under a
  stated balance premise, and the capacity programme asks whether the public
  capacity assigned to the Universe agrees with the one reconstructed from
  inside. These constant calculations, including the fine-structure root,
  have diagnostic status.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## Explore and reproduce

Start with the [main preprint](https://philpapers.org/rec/MUEFOC),
[paper index](paper/), [book](book/) or [reproduction guide](REPRODUCE.md).
The [Lean library](Lean/) contains more than 12300 public theorems and lemmas
with no admitted proofs. The [postdiction ledger](docs/POSTDICTION_LEDGER.md)
and [falsification programme](docs/OPH_FALSIFICATION_PROGRAM.md) distinguish
comparisons with measurements from tests of the mathematical constructions.

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
A controlled continuum limit means that these observations admit an
increasingly accurate continuous description. It does not assert that a
spacetime continuum exists underneath them. The conditional route to Einstein
dynamics has the same target: observers' measurements agreeing with the
effective gravitational description.

Use the existing carrier and rules to discharge inputs by proof or execution.
A supplied field equation remains an input even when its simulation works.
Construction choices need no uniqueness theorem, but their consequences must
survive the tests. Contributions should resolve a named connection or expose
a failure in the shared model. The [construction guide](model/CONSTRUCTION.md),
[candidate register](model/CANDIDATES.md) and
[research issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
set out that work.

Paper releases follow the [model milestone policy](model/CONSTRUCTION.md#research-programme-and-paper-release-hold).
Main-branch PDFs are source previews; deposited releases retain their original bytes.
