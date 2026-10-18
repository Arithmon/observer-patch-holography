# Native clock predictions before laboratory access

This protocol develops the [single-carrier audit](../../model/PHYSICAL_AUDIT.md)
under explicit apparatus assumptions. No experiment has been performed for
this protocol. Expected values below are conditional predictions, not
observed values or evidence that the apparatus assumptions hold.

The useful working hypothesis is: an identified apparatus realizes a specified
native instrument, preparation, scheduler and record interface within a
declared error; an independent clock supplies its duration correspondence.
We can derive consequences and rejection tests under that hypothesis now.
Assuming the expected observations does not supply independent confirmation
of the hypothesis or close the physical spacetime and action obligations.

## What would need to be set up

One carrier is enough. Keep three operations distinct:

1. **Identify the instrument.** Specify the load encoding, actual selective
   maps and records, preparations/resets, allowed coherent controls and
   detector errors. For an engineered emulator, these are implemented
   operations. For a proposed naturally occurring realization, they need an
   independently justified physical mapping. A detection click is not
   automatically a native attempt or a moving current.
2. **Calibrate timing separately.** Record attempt boundaries, including
   zero-current attempts, against an independent reference counter. Retain
   missed/extra-event and timestamp uncertainty. Use a fixed calibration
   subset, then hold out later endpoint readings. The existing
   [clock attachment](../../model/laboratory_clock/README.md) performs the
   bounded-error affine test. Constant duration between every attempt is an
   additional hypothesis; a few endpoint readings alone do not establish it.
3. **Predict unused observations.** Freeze the preparation, schedule,
   spatial dictionary, instrument candidate and decision rule before the
   validation trials. Predict complete arrival/readout laws, including
   timeouts and aborts. Do not estimate the clock rate from the arrival
   measurements it will subsequently be used to predict.

For a narrow classical passage demonstration, twelve occupation bits and
an independent clock can emulate the native binary dynamics. Testing the
coherent C/E distinction additionally needs a quantum encoding, record
ancillas and coherent preparation/readout. A small restricted encoding can
demonstrate that distinction without realizing all 78 source configurations.
Neither demonstration requires a federation.

The minimum retained event table contains trial ID, preparation and reset
status, instrument version, attempt number, addressed seam, coin/current
records including zeros, incident-clock count, independent timestamp and its
uncertainty, intervention/readout setting, outcome, and any abort or timeout.
The receiver's prediction uses only its declared local records and known
preparation; an audit log does not grant it hidden configuration access.

An emulator programmed to implement C is expected to reproduce C. That would
test implementation and timing. An OPH-specific physical test additionally
requires an independently constrained mapping to the physical system, with
predictions that were not installed as its programmed dynamics.

## Concrete expected passage results

Assume the unchanged native rule, independent uniform selection among the
thirty seams, fair placement coins, one initially loaded sender, and an empty
receiver at port zero. Assume a constant `tau=1 ms` per native attempt and a
known launch at an attempt boundary. This time scale is illustrative and is
neither derived from repair nor taken from the published experiments below.

The [exact passage calculation](../../model/physical_space_time/README.md)
then predicts:

| Sender / graph distance to receiver | Mean receiver incident ticks | Mean arrival duration | Arrival by 10 ms | Arrival in (10,30] ms | Arrival in (30,60] ms | Still waiting at 60 ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 / distance 1 | 22 | 132 ms | 0.13546788 | 0.15637791 | 0.13960926 | 0.56854495 |
| 2 / distance 2 | 28 | 168 ms | 0.01936592 | 0.09288837 | 0.14983358 | 0.73791213 |
| 3 / distance 3 | 30 | 180 ms | 0.00403940 | 0.05765454 | 0.14229901 | 0.79600704 |

Each row is one complete four-bin probability law, rounded to eight decimal
places. The mean is for the uncapped arrival distribution; a 60 ms experiment
does not observe that mean directly and must retain its right-censored trials.
These three settings test the shell predictions under the declared symmetry.
Testing the entire passage metric needs all 132 ordered distinct pairs, or
additional held-out symmetry checks with an explicitly restricted claim.

If the assumed/calibrated cadence is instead `[0.99,1.01] ms`, the mean-duration
ranges are `[130.68,133.32]`, `[166.32,169.68]` and `[178.2,181.8] ms. Deadline
probabilities must also change. For a deadline `T`, let `F(n)` be the native
arrival CDF after `n` attempts. With elapsed timestamp error at most `delta`,
sound cumulative bounds are

```text
F(max(0,floor((T-delta)/tau_hi)))
    <= P(recorded arrival duration <= T)
    <= F(floor((T+delta)/tau_lo)).
```

This assumes the stated error bounds and attempt cadence throughout the
trial. An independently established complete-trial total-variation error
`epsilon` widens either probability endpoint by `epsilon`, clipped to `[0,1]`.
Unknown apparatus error stays unknown. The same `tau` applies at every
deadline; choosing a different rate for each histogram cell is not permitted.
The [clock package](../../model/laboratory_clock/README.md) also converts the
native first-current waiting law under this type of conditional calibration.

The table can be reproduced with the existing model, without new data:

```python
from fractions import Fraction as F
from model.physical_space_time.construction import kernel, arrival_law

P = kernel(0)["attempt"]
for sender in (1, 2, 3):
    law = arrival_law(P, sender, 0, 60)
    p = list(map(F, law["first_arrival_probabilities"]))
    bins = [sum(p[:10]), sum(p[10:30]), sum(p[30:60]),
            F(law["timeout_probability"])]
    assert sum(bins) == 1
    print(sender, [f"{float(value):.8f}" for value in bins])
```

## Separate coherent prediction

Use the existing **total-two** coherent preparation and one-seam readout,
not the total-one passage preparation. The
[fixed prediction protocol](../../model/predictions/README.md) already
retains the complete coin/current/readout distribution and all aborts.
For 8,192 fresh trials per setting its ideal marginal expectations are:

| Preparation setting | C: plus probability / expected plus count | E: plus probability / expected plus count |
| --- | ---: | ---: |
| Relative phase 0 | `5/8` / 5,120 | `3/8` / 3,072 |
| Relative phase pi | `1/8` / 1,024 | `3/8` / 3,072 |
| Relative phase pi/2 | `3/8` / 3,072 | `3/8` / 3,072 |
| Equal incoherent mixture | `3/8` / 3,072 | `3/8` / 3,072 |

These are expected counts, not promised outcomes. Use that protocol's full
joint bins, statistical assumptions and apparatus allowance; the aggregate
table is not a replacement acceptance test. The separate
[continuing-state protocol](../../model/closed_prediction_protocol/README.md)
already supplies a finite preparation/error analysis when exact fresh
independence is unavailable. Its different trial count and resources must
not be mixed with this table.

The phase-dependent action family has further coherent predictions in the
[action audit](../../model/repair_action_audit/README.md). Arrival counts alone
cannot distinguish these laws. A successful classical timing test therefore
does not select the action's coherent component.

## Published evidence available without operating a laboratory

The following are complementary opportunities, not measurements of the
native twelve-port model. Their sources were checked on 18 October 2026.

| Primary source and public data | Useful work now | Boundary |
| --- | --- | --- |
| Spiecker et al., [*Solomon equations for qubit and two-level systems: Insights into non-Poissonian quantum jumps*](https://arxiv.org/html/2307.06900v3), [Zenodo 8114088](https://zenodo.org/records/8114088) | Reproduce measured fluxonium jump/dwell statistics and test record-analysis methods against non-Poissonian dynamics. The deposit includes `Quantum_jumps.npz`. | The device couples a qubit to a TLS environment. Its physical detection records need their own observation map; they are not native seam/coin/current records. Histogram bin width is not automatically an attempt duration or accuracy bound. |
| Dalmasso et al., [*Quantum trajectory simulation of two-dimensional non-equilibrium steady states with a trapped ion quantum processor*](https://arxiv.org/html/2605.08350v1), [Zenodo 20075236](https://zenodo.org/records/20075236) | Use hardware-produced mid-circuit records and public analysis/code to investigate complete quantum-trajectory instruments and implementation errors. | The experiment supplies coherent hopping/interactions plus source/drain resets on a square lattice. Its SSEP comparison is classical numerics. Simulated Trotter time is not elapsed laboratory time, and its instrument is not C or E. |
| Minev et al., [*To catch and reverse a quantum jump mid-flight*](https://www.nature.com/articles/s41586-019-1287-z) | Use a demonstrated monitored-jump/feedback setup as a design precedent for timing and conditional intervention. | The paper's data-availability statement offers data on reasonable request, not an immediately downloadable complete record dataset. Its three-level driven system has a different law. |

The Spiecker file was actually downloaded and inspected without executing
external analysis code: 26,376,881 bytes, repository MD5
`fbf4a2c2989c4d4d552d36676692c20e`, one Boolean array of shape
`(25,5,2,2,100,20000)`. It contains processed sampled state traces, rather
than only plotted histograms. The related
[original data/notebook](https://zenodo.org/records/7817552) documents a
2 microsecond measurement grid; the jump-data selection starts after
5 ms. The NPZ itself has no independent timestamps or clock-accuracy fields.
Thresholding, exclusions and censoring require the original processing
contract. Its uncompressed array is about 1 GB, so analysis should account
for that resource cost. No external trace was scored against OPH in this audit.

The Quantinuum archive was also downloaded and its saving code inspected.
It retains per-shot source/drain trajectories, coin outcomes and circuit
parameters, with no timestamp, accuracy or traceability fields in the saved
schema. These records are useful for instrument analysis, but their
simulation step cannot be read as an independent laboratory-clock sample.

The first two provide a practical route for work with existing data. Begin
by reproducing the publication's original observable, with its own units,
selection rules and detector model. Then derive a source-to-apparatus map
for the particular comparison being attempted. Fit any allowed clock or
detector parameters on a declared subset and test unused records or
observables, accounting for temporal dependence. This is a retrospective
benchmark when the published aggregate outcomes informed the design; it
must not be relabelled as a new prospective experiment.

If no faithful mapping exists, the dataset can still assess methods and
rule out an attempted identification. A similar waiting-time histogram
alone does not identify the native instrument. Existing ordinary quantum
models remain essential comparators.

## Work that can proceed now

Use the conditional native predictions above to develop the model, preserve
the alternative actions, and test which combinations of assumptions are
compatible. Use public measurements to constrain apparatus models where an
explicit mapping supports that use. Keep assumed calibration, synthetic
outcomes, measured external records and held-out comparisons separate in
every result. The absence of a local laboratory need not stop theory work;
the unmeasured correspondence remains a premise rather than a completed
physical derivation.
