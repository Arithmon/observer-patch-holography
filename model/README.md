# The runnable OPH model

**This is the live home of OPH's observer and repair-law construction.**
The implementation, tests, derivations and next construction steps belong
here. Start with the finite observer, run its retained results, then extend
its actual interfaces while preserving the demonstrated behavior.

One process executes a local repair, writes an observation, keeps its
classical record and uses that record to choose the next interaction. The
P-dependent version computes the electromagnetic coupling, produces a
binary detector result and lets only that bit control the next repair.
These finite constructions are implemented and runnable under their stated
premises; they do not require a theorem that their design is unique.

## Run the observer

From the repository root, with Python 3.11 or newer:

```sh
# Exact native observation, record preservation, prediction and feedback.
python3 -m model.observer native --check model/observer/results/native.json

# A new 10,000-trial sequential run; prints the result summary.
python3 -m model.observer run --trials 10000 --seed 42

# P → electromagnetic coupling → detector bit → native feedback.
python3 -m model.observer.p_response --check model/observer/p_response/results.json
```

These runtime commands use the Python standard library. They run software
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
`model.verify` regenerates the five retained receipts and checks their hashes.

## What is implemented

| Component | Result and entry point |
| --- | --- |
| Native observer | [Native C, retained tape and record-only controller](observer/README.md), including all 264 benchmark paths and exact feedback controls |
| Binary electromagnetic observation | [P-dependent detector](observer/p_response/README.md): evaluates α(P), writes one stable bit and uses it for the next native repair; all 528 weighted paths |
| Finite capacity | [Correctable record capacity and compact controller](observer/CAPACITY.md): 11 terminal native words, 16 reachable prefixes and a four-bit classical encoding under the binary preparation |
| Direct execution | [20,000-trial native receipt](observer/results/run.json) and [200,000-trial detector receipt](observer/p_response/results.json), alongside exact references |
| Local event clock | [N4 local-tick calculation](local_tick/README.md): corrected native null and certified positive memory bounds while remote events continue |

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

## Continue the construction

- [Construction guide](CONSTRUCTION.md): current contracts, P/N connection,
  completed results, remaining connections and research policy.
- [Candidate register](CANDIDATES.md): native C/E, response extensions,
  reciprocal N4 and the separate scalar benchmark.
- [Model–theorem premise map](PREMISES.md): generated navigation to the
  existing hypotheses and theorem attachments.
- [Mathematical foundations](../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md) and
  [working terms](../docs/INTERPRETATION_AND_TERMS.md): compatibility,
  scheduling, observer access and the distinction between history and state.

The broader [executable evidence](../code/README.md) and
[reproduction guide](../REPRODUCE.md) remain the home of supporting source,
field, geometry and parameter calculations. The new observer construction
continues here; its results do not silently change those existing claims.

Software is licensed under [Apache-2.0](../code/LICENSE). Documentation and
data follow the repository's [licensing map](../LICENSE).
