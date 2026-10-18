# A falsifiable attachment to an independent clock

One elementary carrier is sufficient to derive and test local clock relations.
A federation is not needed for this task. The native rule already supplies
event ordering and complete current-count laws; it does not choose how many
seconds elapse between attempts. This packet adds an executable **conditional
clock calibration and held-out test**, with event-detection and timestamp
errors kept explicit. No new laboratory measurements are claimed.

The [existing clock bridge](../native_clock/CLOCK_BRIDGE.md) relates incident
hits, moving currents and fields on the native source. The
[host experiment](../native_clock/HOST_CLOCK.md) measures software execution
against a monotonic counter. This packet addresses a different missing step:
given paired event counts and independently measured timestamps, what clock
correspondence is compatible with their stated errors, and which unused
measurements can falsify it? The numerical solver never authenticates a
measurement's provenance or declares an empirical identification with OPH.

## What can be derived without a laboratory

Let an admissible native history have finite record word `r` and probability
`P(r)`. Let the source operations, scheduler and feedback depend only on
native states and records. Assign its attempt boundaries any strictly
increasing timestamps `f(0), f(1), ...`. This assignment changes none of the
conditional native operations, hence none of the complete record laws,
including terminal states and untouched references. In particular both
`f(n)=b+tau*n` and a nonlinear increasing `f` preserve these laws. This is
proved directly by induction on the identical conditional operations.

Consequently native record data alone identify neither a positive scale
`tau` in seconds per attempt nor constant cadence. This is a missing physical
correspondence, beyond the conventional choice of the unit second. External
time-triggered interventions would add a timing premise and are outside this
invariance. The argument requires no preferred universal clock and no
federation. Adding carriers with only the same record-defined scheduling
does not remove it.

The receipt checks the actual 66-configuration, twelve-port native source's
marked reduction again. The same four-bit word retains the same probability
at illustrative rates `1/1000`, `1/100` and `1` seconds per attempt, while its
stationary mean movement frequency becomes `5000/33`, `500/33` and `5/33`
per second. These are alternative supplied clock attachments, not measured
frequencies. The construction extends the existing
[duration-scale audit](../native_clock/README.md#why-the-model-does-not-predict-seconds)
to offset, nonconstant retiming and a complete error-aware attachment test.

## Exact bounded-error clock theorem

Independently declare the apparatus hypothesis

```text
t_i = b + tau*n_i,       tau > 0,
n_i in [n_i^-, n_i^+],   t_i in [t_i^-, t_i^+].
```

Here `n_i` counts all native attempts, including zero currents; `t_i` is an
independent clock reading in seconds; `b` is its arbitrary origin. The
intervals must enclose simultaneous detection and timestamp errors. A clock's
reported resolution or one standard uncertainty alone is not such a bound.
The constant-cadence equation is an apparatus hypothesis, not a consequence
of native repair. Movement-only counters have different stochastic waiting
laws and cannot be substituted into this interface.

For any `tau>0`, endpoint `i` permits exactly the real-coordinate offsets

```text
b in [t_i^- - tau*n_i^+, t_i^+ - tau*n_i^-].
```

A common offset exists precisely when every lower endpoint is below every
upper endpoint. Thus all pairwise inequalities

```text
t_i^- - t_j^+ <= tau*(n_i^+ - n_j^-)
```

give the exact feasible rate interval by rational arithmetic. A positive
coefficient gives a lower bound, a negative coefficient an upper bound, and
a zero coefficient a consistency condition. Zero rate is always excluded.
For uncertain integer counts this is a continuous outer relaxation:
infeasibility still rejects, while a feasible real count may have no integer
realization. Overlapping uncertain count intervals also receive no joint
true-count ordering constraint. Machine-readable outputs mark this outer
feasibility and never certify an integer realization in that case. With exact
counters, no count relaxation occurs.

The calibration endpoints alone determine a fixed set of `(tau,b)`. For a
new count interval `[n^-,n^+]`, the sharp outer prediction is

```text
min_tau max_i [t_i^- + tau*(n^- - n_i^+)],
max_tau min_i [t_i^+ + tau*(n^+ - n_i^-)],
```

where `tau` runs over the calibration interval. These are one-dimensional
piecewise affine extrema; the code evaluates interval ends and all relevant
line intersections exactly. This retains the correlation between origin and
rate. It then tests whether **one common calibrated pair** fits all validation
endpoints jointly. Individual validation intervals can each overlap while
the joint test fails; a retained negative control demonstrates this case.
Validation never supplies a replacement fitted rate.

If the measurement bounds hold and a common affine clock exists, this test
cannot reject it. Rejection falsifies their conjunction. Non-rejection says
only that some calibration remains possible. No statistical confidence level
is inferred from the file; one would require independently justified
simultaneous coverage for its bounds. Endpoints cannot establish constant
cadence at unobserved intervening attempts or beyond the tested range.

The proof is finite elimination of linear inequalities, not a new theorem of
clock physics. The nearest existing model results are the native count
additivity/monotonicity attachment and exact marked waiting laws in
[native_clock](../native_clock/README.md). Their source and scheduling premises
are reused unchanged. The addition is this parameter/error elimination,
held-out joint test and the following complete waiting-law conversion.

## From a calibrated cadence to a native laboratory prediction

Under the additional hypothesis of constant cadence **at every attempt**,
take elapsed duration from a known attempt boundary. If `W` is the attempt
count of the first moving current and `T` a nonnegative deadline in seconds,

```text
P(tau*W > T) = P(W > floor(T/tau)).
```

For a calibration `0 < tau_lo <= tau <= tau_hi`, monotonicity gives the exact
outer probability interval obtained at these two endpoints. The producer
reconstructs `P(W>m)` from the source marked kernel for both stationary and
actual-current-conditioned (Palm) preparations. It keeps the no-current tail;
it does not assign the mean wait to every realization or assume an exponential
clock. The named preparation and ideal scheduler remain required inputs.

The [synthetic input](endpoints_example.json) assigns endpoint errors of one
microsecond and calibrates on attempts 0 and 1,000. It yields
`tau in [499999/500000000, 500001/500000000]` seconds per attempt. Its unused
endpoints at 2,000 and 3,000 attempts pass; moving their timestamps by 0.1 s
rejects with the identical calibration. At a 10 ms deadline, the native
waiting-law conversion permits 9–10 completed attempts, with a nonzero
no-current probability in either case. All these assigned timestamps and
error bounds are software fixtures, not measurements.

A calibrated attempt scale would give dimensional **rates** for laws already
parameterized by attempts. It does not turn a Euclidean transfer generator
into real-time unitary evolution, select an energy/action normalization, or
identify the modular parameter with native attempts. An incident-clock action
also needs its joint clock law; its stochastic counter cannot be replaced
pathwise by a constant multiple of this attempt counter.

## What the actual retained measurements establish

The replay checks both original host executions: each has 264,532 native
attempts and the same complete seam/coin/current tape. The existing data
have no traceable timestamp-error bound. Their input therefore retains
`seconds_error: null`, and the new bounded-error numerical test remains
unresolved. The reported host-clock resolution is not promoted to accuracy.

As an explicitly counterfactual zero-error check, calibrating the first
baseline interval gives `2941500/5977` nanoseconds per attempt. Later baseline
endpoints already reject an **exact** affine cadence, as does the same-record
workload control. This does not contradict the existing 25% engineering
tolerance for movement-clock predictions; that was a different, deliberately
approximate contract. Its calibration, observations and verdicts are preserved.
The present retrospective attempt-cadence test was not part of the original
capture policy and is not described as a prospective experiment.

## The empirical input still required

Supply one identified apparatus with the native preparation, seam/coin
selection, recording and intervention semantics; paired independent
timestamps and attempt records including zeros; justified simultaneous
missed/extra-event and timestamp bounds; and a calibration/validation split
fixed before inspecting validation outcomes. To claim SI traceability, the
measurement also needs its documented calibration chain and associated
uncertainties, as specified by
[NIST's metrological traceability policy](https://www.nist.gov/metrology/metrological-traceability).
An SI second is the time unit; realizing and comparing it is a further
measurement operation, as distinguished in the
[BIPM SI guidance](https://www.bipm.org/en/measurement-units/si-base-units).

No such complete native-apparatus laboratory dataset is supplied here.
Existing oscillator fits use supplied apparatus evolution, and host execution
timings measure an engineered implementation. Neither silently supplies the
missing native physical identification. A traceable implementation benchmark
would still establish that apparatus's timing contract, not that its native
events are physical proper time in the Universe.

## Reproduce or test endpoint data

```sh
python -m model.laboratory_clock --check model/laboratory_clock/results.json
python -m pytest -q -W error model/laboratory_clock
python -m model.laboratory_clock --data model/laboratory_clock/endpoints_example.json
```

`--write PATH` refuses to overwrite a receipt. The input has three to 64
endpoints, at least two for calibration and one for validation, with strictly
increasing reported native counts. Time values and errors are rational
strings. Unknown counter or clock errors remain `null` and prevent the
numerical test. Nonempty provenance references are retained but never treated
as authentication. The full receipt binds source files and raw host records;
replay neither acquires nor reproduces physical timestamps.
