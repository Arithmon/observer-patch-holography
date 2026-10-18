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

## Compiled operations require separate timing coefficients

The [compiled clock history](../compiled_clock_history/README.md) executes
source preparation, a coherent intervention and terminal readout around the
same capped native repairs. Native repairs and source-pulse invocations are
different resources. The source does not equate their durations. A zero-angle
control also counts as a requested invocation in that compiler.

[`resource_clock.py`](resource_clock.py) adds the conditional timing model

```text
t_i = b + tau_native*N_i + tau_pulse*P_i,
tau_native > 0, tau_pulse > 0.
```

Each calibration endpoint identifies both exact operation counters and a
bounded independent reading. The same origin and duration coefficients must
serve all calibration and validation endpoints. Endpoints may describe
different calibration workloads under a common zeroing convention; they
are not required to form one chronological counter tape. The prediction
uses elapsed differences from a declared protocol start, so the common
origin cancels. Timestamp error must bound that difference, including both
readings. The coefficients are supplied timing hypotheses, not physical
quantities selected by repair or by compiling a matrix.

Per-pulse load/export latency can be included in the pulse coefficient if
it has the same fixed multiplicity and timing law. Other initialization,
complement-check, record-delivery and detector latency needs an explicit
bound or additional resource. It cannot disappear from the physical
comparison because a software array contracts it instantaneously.

Eliminating the common origin gives the exact closed inequalities

```text
(N_j-N_i)*tau_native + (P_j-P_i)*tau_pulse <= t_j^+ - t_i^-.
```

The implementation uses rational Fourier--Motzkin elimination to project
this two-dimensional region onto any elapsed duration
`N*tau_native + P*tau_pulse`. It preserves the correlation between the
coefficients. Independent marginal intervals generally give a wider bound.
The closed region allows zero coefficients as boundary limits; the output
separately tests whether strictly positive coefficients exist and whether
both are bounded away from zero. An unbounded pulse duration is reported as
unbounded. No value is assigned to an unobserved resource.

If every calibration workload has `P=2N`, observations determine only
`tau_native+2*tau_pulse`. Distinct positive coefficient pairs reproduce
every such observation and disagree on a subsequent native-only workload.
This is an identifiability obstruction, independent of measurement
precision. The same issue occurs when no calibration endpoint varies the
pulse count. Validation tests one common calibration against all held-out
endpoints jointly; individually compatible endpoints can fail together.
No replacement duration pair is fitted using validation data.

### Complete cost and detector law

The native endpoint replay retains the number of actual scheduled attempts,
including hidden zero-current attempts and the feedback repair. This is an
additional controller counter for the calibration protocol; it is not an
assertion that an observer can read hidden microscopic events in nature.
At cap two its maximum needs three counter bits. Forgetting this new count
recovers every full conditional matrix of the established observer history,
for both coherent settings and for each executed cap from one through three.
This full-matrix identity concerns the ideal native history. Source-pulse
cost metadata is imported from the compiler; a resource-refined execution
of its raw binary arrays is not claimed. A calibration bound for coarser
current records does not by itself certify this additional count record.
Applying an instrument-error bound here requires a proof before that
coarse-graining, with the counter service explicitly admitted.

Source pulse costs are taken from the compiler's resource interface. For
the cap-two, non-reference preparation, the complete resource law is

| Native repairs `N` | Source pulses `P` | Probability |
| ---: | ---: | ---: |
| 2 | 27 | `9/25` |
| 3 | 51 | `4/25` |
| 4 | 51 | `42/125` |
| 5 | 51 | `18/125` |

Every timeout and both detector outcomes contribute. The means are
`E[N]=408/125` and `E[P]=1059/25`; the two coherent settings have the same
resource marginal. Their complete resource/clock/detector laws differ by
total variation `1848/15625`. No duration coefficient has been used to
derive these dimensionless quantities.

For a deadline, each resource pair determines a linear duration on the
calibration region. An event is forced before the deadline if its maximum
duration, including reading error, is below the deadline; it is possible
if its minimum duration allows it. Summing the corresponding exact masses
gives valid outer probability bounds. For this componentwise ordered cost
support and a bounded calibration with strictly positive coefficient
bounds, they are sharp: the possible deadline sets are nested prefixes of
the cost ordering. The proof applies to the event jointly with either
detector outcome. It does not replace the calibration region by independent
fits for different records.

The synthetic fixture calibrates on `(N,P)=(0,0),(100,0),(100,100)`, with
assigned readings `0`, `0.1 s +/-1 us`, `0.102 s +/-1 us`. It gives

```text
tau_native in [99999/100000000, 100001/100000000] seconds,
tau_pulse  in [999/50000000, 1001/50000000] seconds.
```

An unused `(200,50)` workload at `0.201 s +/-1 us` passes. Shifting that
unused reading by `10 us` rejects with the original calibration unchanged.
The sharp expected-duration interval is
`[0.00411038544,0.00411201456] s`. All readings and errors here are assigned
software fixtures, not measured data.

At a `4.5 ms` elapsed deadline, with a further `1 us` reading-error bound,
the predicted probability of finishing is `13/25` at either setting.
The joint probability of finishing with detector one is `2451/10000` for
`gamma=0` and `543/2000` for `gamma=1`. These are unconditional masses;
timeouts contribute and are not discarded. At `4 ms` the exact nominal
cadence gives finish probability `9/25`, whereas substituting the mean
duration gives zero. A mean cannot substitute for the stopped law.

The evaluator accepts normalized ideal probability laws. The raw compiled
arrays can have a trace defect and are not silently normalized for this
use. For an authenticated normalized implementation law within total
variation `eta` of the identified ideal law, widen an event interval
`[L,U]` to `[max(0,L-eta),min(1,U+eta)]`. A complete raw-array half-nuclear
error `epsilon` alone instead permits event-weight error `2*epsilon`,
or `epsilon+abs(trace_defect)/2` if the defect is retained. Those numerical
array bounds are not physical service-error certificates. Simultaneous
coverage of calibration and timestamp bounds must also be justified before
assigning a statistical rejection level.

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
