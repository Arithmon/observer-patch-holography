# Runnable finite observers

This directory implements a bounded observation–record–feedback process from
OPH's native local repair. The same declared process writes an observation,
retains its classical record, uses that record to choose the next interaction,
and halts when its two event slots are full. The software runs both exact
finite calculations and sampled sequential trials.

The [P-dependent binary detector](p_response/README.md) adds a concrete
electromagnetic response under accepted P/N closure inputs. Its actuator reads
one detector bit. The native current recorder described here retains a richer
coin/current letter and makes a binary seam decision. **Two events** counts
sequential interactions; it does not mean that the whole recorder has only two
possible records.

The [scalar detector](scalar_response/README.md) uses the first native current
to prepare a signed intervention in the existing scalar field, executes a
binary quantum readout and lets its retained bit select the second native
repair. Its supplied action and vacuum have an inherited continuum comparison;
the [prediction protocol](scalar_response/PREDICTION.md) states the feedback
contrast, null and reversal controls. Install `model/requirements.txt` for
this adapter; the native and P-dependent commands below remain standard-library
executions.

The [source-compiled two-site apparatus](source_scalar/README.md) executes a
separate finite field design through emitted twelve-port control words,
including entangled preparation from blanks and both pointer branches.
It retains exact initial-tail and finite temporal error certificates with
the native feedback. Coherent code services and the scalar action remain
supplied; the earlier field family's continuum bounds do not transfer to
this two-site implementation.

## Run from the repository root

The observer code uses the Python standard library. The capacity calculation
imports the existing [public-record evaluator](../../code/capacity_readback/).

```sh
# Exhaustive exact construction: 78 source states, 264 benchmark histories.
python3 -m model.observer native --check model/observer/results/native.json

# Exact record capacity and lossless classical controller encoding.
python3 -m model.observer capacity --check model/observer/results/capacity.json

# Execute 20,000 fresh, sequential two-event trials and replay the receipt.
python3 -m model.observer run --trials 20000 --seed 20261012 --check model/observer/results/run.json

# A new run, written only if the destination does not exist.
python3 -m model.observer run --trials 10000 --seed 42 --write /tmp/observer-run.json

# Independent exact checks, capacity tests, and seeded replay.
python3 -m unittest model.observer.test_native_observer model.observer.test_observer_capacity model.observer.test_observer_simulation -v
```

`--check` compares canonical bytes, and `--write` refuses to overwrite a file.
Each fresh trial has two native events, followed by an explicit halt. Running
many trials is repeated preparation, not extending one observer's memory.
The clock parameter here is native operation count, without a physical-duration
calibration. Host runtime is merely the cost of running the Python program.

## The construction

The working configurations are all nonnegative integer loads on twelve ports
with total load two. Their carrier has 78 basis states. The scored preparation
is uniform on the 66 configurations with two distinct occupied ports.
`native.states()` generates the full carrier directly; no external state
fixture is needed.

For a seam `(a,b)`, fair coin `c`, and contrast `d=x[a]-x[b]`, the local repair is

```text
j = (d - (2*c - 1)*(abs(d) mod 2)) / 2
x'[a] = x[a] - j
x'[b] = x[b] + j
record = (c,j)
```

Compatibility on a seam is specified independently: its endpoint loads differ
by at most one. The update averages them by floor and ceiling, preserving
total load and all spectator ports. Applying another seam may disturb an
earlier seam's balance; this two-event observer is not a global convergence
claim. At event 1 the seam is `E=(0,1)`. At event 2 it is `F=(0,5)` exactly when
the retained first current is positive; otherwise it is E. The program then
halts. The supplied scheduler, preparation, tape, policy and forecast table
are explicit parts of the construction.

The [derivation](DERIVATION.md) defines one fixed completely positive,
trace-preserving map on the 3354-dimensional source/tape carrier. Its coherent
writes are isometries on every nonfull sector, including inputs outside the
benchmark support. Committed classical records survive feedback. Coherence
between different record values need not survive: an exact reachable matrix
entry changes from `1/4` to `1/8`. Reading a forecast is a separate classical
readout and may dephase the record; the complete scored history probabilities
are unchanged by that readout.

## Measured behavior and controls

The stored first record supplies a forecast of positive current at event 2,
issued before that event. This table is a supplied model, with no online
learning claim. The code enumerates all 264 benchmark histories with exact
weight `1/264`; no branch is discarded.

| First record `(c,j)` | Probability | Forecast `Pr(j2 > 0 \| c,j)` |
| --- | --- | --- |
| `(0,0)` | `14/33` | `0` |
| `(1,0)` | `14/33` | `5/56` |
| `(0,1)` | `5/66` | `0` |
| `(1,-1)` | `5/66` | `1/2` |

The other two letters occur in the full carrier but not this preparation;
the forecast interface abstains on those letters.

| Exact quantity | Result |
| --- | --- |
| Positive second-current probability | `5/66` |
| Brier risk using the first record | `395/7392` |
| Brier risk using only the chosen seam | `140/2013` |
| Brier risk with independent records shuffled within the same seam | `12875/150304` |
| Mean final port-0 load with actual record feedback | `35/264` |
| Mean final port-0 load with disconnected feedback | `1/6` |
| Mean final port-0 load with independent donor-record feedback | `1/6` |

The donor selector has the same seam frequencies as the actual selector, so
the difference `-3/88` detects the record's causal effect beyond a changed
action distribution. Final source loads are external verification outputs;
the observer's policy and forecast do not read them.

The [capacity calculation](CAPACITY.md) applies the existing OPH checkpoint
confusability construction to these records. Exact terminal rereads retain
11 reachable record words on the binary preparation, or 17 on the full
78-state input family. The binary experiment has 16 reached record prefixes;
four bits encode its full classical controller history. That count excludes
the working source, program tables and Python workspace and is not a quantum
compression statement or a cosmological horizon capacity.

## Files and scope

| File | Role |
| --- | --- |
| [native.py](native.py) | Native repair, fixed controller, exact quantum and history checks |
| [capacity.py](capacity.py) | Existing capacity-evaluator adapter and compact controller |
| [simulate.py](simulate.py) | Seeded sequential execution, forecasts and paired controls |
| [results/native.json](results/native.json) | All exact benchmark histories and 18 calculation checks |
| [results/capacity.json](results/capacity.json) | Reachable records, carrier certificates and code tables |
| [results/run.json](results/run.json) | A reproducible 20,000-trial run with visible traces |
| [p_response/](p_response/) | P-dependent binary electromagnetic observation and feedback |
| [scalar_response/](scalar_response/) | Supplied scalar field, binary instrument, native continuation and controlled refinement |

These are implemented finite observer functions with designed apparatus and
protocol choices. Uniqueness of that design is not a completion requirement.
The result is a runnable construction, without a claim that a whole universe,
a continuing unbounded observer, a learned law, or a physical computer has
been implemented. The [project's interpretation](../../docs/INTERPRETATION_AND_TERMS.md)
distinguishes those broader objects from this finite execution.
