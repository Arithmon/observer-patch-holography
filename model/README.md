# The runnable OPH model

**This is the live home of OPH's observer and repair-law construction.**
The implementation, tests, derivations and next construction steps belong
here. Start with the finite observer, run its retained results, then extend
its actual interfaces while preserving the demonstrated behavior.

The newest step is the [common interacting source](common_source/README.md):
quantum electric flux, charged CAR matter and a scalar field share one
evolving register, with complete measurements, retained records and native
feedback. For zero native preparation current, its scalar response has an
[exact nonzero infinite-occupation interval](common_source/INFINITE_HISTORY.md), including the actual
coordinate-sign measurement and the preparation, propagation and readout
errors of that calculation.
A separate moment bound controls finite occupation through retained matter
feedback for arbitrary matter/reference inputs. This
binds several theorem interfaces to one process; the interacting
spatial continuum, gravity and independent physical comparison remain open.

The [observer readout bridge](common_source/OBSERVER_READOUTS.md) extracts
Hilbert representations and modular flow from that actual observation state,
continues its records with a source-executed binary modular clock, and
computes the interaction graph of its supplied quantum regions. Its
[occupation analysis](common_source/OCCUPATION_REFINEMENT.md) establishes
a fixed two-site interacting limit, with a finite error bound before the
last sharp readout and strong convergence of that readout. A finite-cutoff
error rate for the sharp instrument remains open. These results connect
existing components; they do not supply physical time or a spatial continuum.

A separate [finite-resolution instrument](common_source/SMOOTH_REFINEMENT.md)
has a quantitative error bound for its complete interacting history,
including the last measurement. Its retained energy ledger compares actual
finite-cutoff backreaction with the infinite-model unconditional `1/8`
energy increment. The supplied energy units and two-site action remain
construction assumptions; this is not an empirical energy calibration.
The [readout checks](common_source/CALIBRATION.md) test geometric distances,
modular-parameter inference and oscillator energy convergence against their
mathematical targets, and make the remaining unit-calibration freedom explicit.

The [fixed-reference modular bridge](common_source/THERMAL_CLOCK.md) connects
continuous segments of the supplied full Hamiltonian evolution to the
modular flow of a chosen faithful Gibbs reference, with `t=-beta s`, while
the actual observed state changes. This reference and full-work algebra
differ from the earlier local-reference clock. Its clock probes continue
the finite sharp-readout histories; neither clock is included in the
infinite sharp-signal certificate or the smooth-detector occupation bound.
Reference selection by repair, physical duration and the relation to
spatial propagation remain open.

The [local history-repair formulation](common_source/HISTORY_REPAIR.md)
recovers that same completed observation history under every fair update
schedule, for the fixed source apparatus and causal graph. Observer readouts
are functions of the completed history; the solver counter is auxiliary.
The native quantum realization of this global constraint-repair algorithm
remains a separate construction requirement.

The [rank-three native carrier](common_source/SLOW_GEOMETRY.md) is connected
to the same joint observation histories. A binary controller extension
executes displacements spanning three directions while preserving the
original protocol. This supplies a finite geometric readout, with physical
translation, momentum and spatial field refinement remaining open.
The [native propagation checks](common_source/CAUSALITY.md) verify finite
dependency cones and delayed influence under fixed local seam schedules.
They distinguish causal interaction depth from an auxiliary solver counter.
The [native field relations](common_source/NATIVE_FIELDS.md) attach exact
continuity and conditional current laws to those same recorded transfers.
The existing uniform-seam calibration gives a closed graph diffusion law
for the mean load. Its decaying modes do not supply an undamped wave law;
coherent field dynamics remain a separate construction obligation.

The [completion target](CONSTRUCTION.md#completion-target-and-acceptance-criteria)
is an executable observer algebra whose compatible measurement statistics
admit the observed 3+1-dimensional quantum-field, classical and gravitational
descriptions in their respective regimes. Fundamental fields or spacetime
need not be reconstructed first. The
[operational meaning and objectivity contract](../docs/INTERPRETATION_AND_TERMS.md#observer-measurements-and-operational-objectivity)
states what inter-observer agreement must preserve. These are research
success criteria; the finite implementations below are components of that goal.

One process executes a local repair, writes an observation, keeps its
classical record and uses that record to choose the next interaction. The
P-dependent version computes the electromagnetic coupling, produces a
binary detector result and lets only that bit control the next repair.
These finite constructions are implemented and runnable under their stated
premises; they do not require a theorem that their design is unique.
The [scalar detector](observer/scalar_response/README.md) connects that
observer to the existing free-field construction: its first current prepares
a signed scalar intervention, a binary field read controls the next repair,
and the inherited refinement bounds constrain the resulting history law.
The scalar Hamiltonian, vacuum, readout and model-time convention are supplied.
Its local evolution extension executes neighboring-site scalar updates at
q=13 and q=34 and uses their detector probabilities for native feedback.
Original-vacuum squeezing and a q=13 step-halving control are included;
the inherited temporal and spatial certificates remain distinct from
floating numerical error and physical apparatus error.
Its separate [continued version](observer/scalar_response/CONTINUATION.md)
uses continuous spectral evolution and keeps the disturbed field for a
second read and a third native event, with explicit record memory and a
bound on the complete history. Combining this retained two-read process
with the local update engine remains further work.
The [source-compiled apparatus](observer/source_scalar/README.md) adds a
separate two-site quantum field: its entangled preparation from blanks,
field evolution and binary pointer all execute emitted twelve-port control
words. Four-, eight- and sixteen-level site codes retain both quantum branches,
exact complete occupation-history and finite splitting bounds, and the
native feedback. The occupation theorem covers this two-oscillator action;
it does not establish spatial refinement.
This source program performs one field read; it has not been joined to the
spectral two-read continuation.

The [prediction package](predictions/README.md) specifies two competing
quantum observer instruments, computes their exact phase-dependent
probabilities and fixes a count-based rejection protocol before independent
data. Its evaluator retains every record and outcome, including aborts.
Labeled synthetic controls exercise rejection; no independent experimental
counts are supplied. A separate retrospective coupling audit preserves the
inherited P-response gap and its physical-transport boundary.

The [source-executed C/E experiment](observer/source_discrimination/README.md)
implements both complete 78-state repair instruments and their fixed
coherent preparation/readout test with twelve-port control words. An exact
channel identity preserves arbitrary reference coherence, and the executed
terminal states feed the unchanged count protocol. Physical encoding and
service calibration remain separate from this finite construction.
E is C with one extra seam-parity record, implemented by three additional
source pulses; discarding that label does not reverse the lost interference.

## Run the observer

From the repository root, with Python 3.11 or newer:

```sh
# Exact native observation, record preservation, prediction and feedback.
python3 -m model.observer native --check model/observer/results/native.json

# A new 10,000-trial sequential run; prints the result summary.
python3 -m model.observer run --trials 10000 --seed 42

# P → electromagnetic coupling → detector bit → native feedback.
python3 -m model.observer.p_response --check model/observer/p_response/results.json

# Continue the same scalar field through two observations (NumPy/SciPy required).
python3 -m model.observer.scalar_response --continue-field \
  --check model/observer/scalar_response/continuation_results.json
```

The native and P-detector commands use the Python standard library. They run software
models, with a fresh preparation for each two-event trial. Event order is
the time parameter; the program does not assign a physical duration to a
Python invocation. The P detector evaluates the existing source-only
structured-running coupling branch, with its original endpoint status.

For the complete checks, including the local-tick certificate:

```sh
python3 -m venv .venv-model
. .venv-model/bin/activate
python3 -m pip install -r model/requirements.txt
python3 -m model.verify
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model
```

The [observer README](observer/README.md) gives exact capacity and seeded-run
replays, new output options and the full controller contract. Result writes
refuse to overwrite an existing file; use `--write PATH` to retain the full
result and example traces. The [local-tick README](local_tick/README.md)
gives its certificate replay and the incident-attempt stopping rule.
`model.verify` replays the retained results and checks their artifact hashes.
The scalar local-evolution and source-compiled receipts can also be replayed directly:

```sh
python3 -m model.observer.scalar_response.local_evolution --check model/observer/scalar_response/local_results.json
python3 -m model.observer.source_scalar --check model/observer/source_scalar/results.json
python3 -m model.predictions --check model/predictions/results.json
python3 -m model.observer.source_discrimination --check model/observer/source_discrimination/results.json
python3 -m model.common_source --check model/common_source/results.json
```

## What is implemented

| Component | Result and entry point |
| --- | --- |
| Common interacting source | [Gauge/matter/scalar observer](common_source/README.md): actual source words, operator Gauss, two committed bits, all 1,056 native/quantum paths and a certified finite feedback response; declared action and services |
| Observer Hilbert space, modular clock and regions | [Same-history readouts](common_source/OBSERVER_READOUTS.md): GNS pairings, nontrivial local modular flow, complete source-compiled clock branches and CP regional maps; access, schedule and localization declared |
| Timeless history and native rank-three carrier | [Local history repair](common_source/HISTORY_REPAIR.md) and [geometric readouts](common_source/SLOW_GEOMETRY.md): schedule-independent completed branching law and record-readable displacements, including an executed three-direction controller |
| Native observer | [Native C, retained tape and record-only controller](observer/README.md), including all 264 benchmark paths and exact feedback controls |
| Binary electromagnetic observation | [P-dependent detector](observer/p_response/README.md): evaluates α(P), writes one stable bit and uses it for the next native repair; all 528 weighted paths |
| Scalar field observation and refinement | [Scalar detector](observer/scalar_response/README.md): native-current preparation, binary field instrument and bit-only feedback; local field updates at q=13,34, with temporal and spatial comparison bounds and a [prediction protocol](observer/scalar_response/PREDICTION.md) |
| Retained scalar continuation | [Two reads of one field](observer/scalar_response/CONTINUATION.md): continuous spectral evolution, three native repairs, all 2,112 labelled histories and 209 reached record prefixes, with a separate complete-history bound |
| Source-compiled quantum field apparatus | [Two-site source program](observer/source_scalar/README.md): preparation from blanks, field evolution and quantum pointer from twelve-port gates; both terminal field branches, finite error certificates and native feedback |
| Specified alternatives and rejection test | [C/E phase experiment](predictions/README.md): exact complete outcome laws, fixed error and sample budgets, strict independent-count evaluation and an explicit no-data verdict; separate retrospective coupling audit |
| Finite capacity | [Correctable record capacity and compact controller](observer/CAPACITY.md): 11 terminal native words, 16 reachable prefixes and a four-bit classical encoding under the binary preparation |
| Direct execution | [20,000-trial native receipt](observer/results/run.json) and [200,000-trial detector receipt](observer/p_response/results.json), alongside exact references |
| Local event clock | [N4 local-tick calculation](local_tick/README.md): corrected native null and certified positive memory bounds while remote events continue |
| Central ports and occupation readout | [Native/N4 address adapter](../code/native_response/README.md#a-central-port-carrier-that-retains-the-quantum-repair): tested finite channel and readout operations, with an analytic abstract source extension and explicit record costs |

The controller's binary choice is an interface. The native record itself
contains coin and current, and the full detector audit retains more than
one bit. “Two-event” describes the finite protocol length, not the number
of possible source states or histories. Feedback preserves the committed
classical record while it can disturb quantum coherence between record values.

For the electromagnetic design, P and N are accepted closure-derived inputs.
The charge assignment, probe, calibrated exposure, record hardware and
protocol are construction choices. P enters the detector probability and
the feedback response; N bounds capacity on a declared physical carrier.
This does not claim fabricated hardware, a complete cosmological simulation,
an online learner or physical validation of the entire model.
The scalar adapter executes a supplied free-field model; native integer
averaging does not by itself derive that Hamiltonian or its quantum vacuum.
Its detector prediction is conditional on the stated apparatus and field
law, with independent software controls and a separate physical test boundary.
The two-site source program implements field preparation, evolution and
detection under its declared code services. A Gaussian-history tail bound
controls the complete ideal occupation truncation through impulse, evolution
and pointer readout, including terminal field states. Adding finite splitting
gives a common ideal two-oscillator error budget. The d=8 bound remains too
conservative to certify the small detector bias. At the executed d=16 cutoff
and 128 steps, an exact rational interval places the ideal plus probability
between 0.51495 and 0.51789, excluding a fair response. Spatial refinement,
numerical and physical service errors, apparatus identification and independent
comparison require their own evidence.

## Continue the construction

- [Construction guide](CONSTRUCTION.md): current contracts, P/N connection,
  completed results, remaining connections and research policy.
- [Finite-to-physical connections](CONSTRUCTION.md#connecting-finite-objects-to-the-physical-theorem-routes):
  the actual state, record and refinement maps required when reusing the
  existing scalar, gauge and gravity results.
- [Candidate register](CANDIDATES.md): native C/E, response extensions,
  reciprocal N4 and the separate scalar benchmark.
- [Audit controls](CONSTRUCTION.md#retained-audit-controls): retained
  counterexamples and the precise carrier, response and kinetic attachments.
- [Model–theorem premise map](PREMISES.md): generated navigation to the
  existing hypotheses and theorem attachments.
- [Mathematical foundations](../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md) and
  [working terms](../docs/INTERPRETATION_AND_TERMS.md): compatibility,
  scheduling, observer access and the distinction between history and state.

The broader [executable evidence](../code/README.md) and
[reproduction guide](../REPRODUCE.md) remain the home of supporting source,
field, geometry and parameter calculations. The new observer construction
continues here; its results do not silently change those existing claims.
The [scientific registers](../tracking/README.md) keep finite construction
evidence separate from physical observation status and parameter ancestry.

Software is licensed under [Apache-2.0](../code/LICENSE). Documentation and
data follow the repository's [licensing map](../LICENSE).
