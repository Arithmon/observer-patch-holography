# A native current counter and its measurement laws

The same executable binary observer now predicts its current-count and
waiting-record distributions. No independent random clock or Hamiltonian is
added. Its 66 native configurations reduce exactly to a three-state process
for **every finite absolute-current record word**. This reduction makes longer
record predictions inexpensive and identifies which scheduler assumptions can
be removed.

The counter is `N(n) = sum(abs(j_k), k=1..n)`, where `j_k` is the already
committed native current. Its readings are observable integer event counts.
Here `n` always means **native attempt count**, including the recorded zero
currents. A physical duration, proper time or relation to an experimental
clock has not been identified.

The [single-carrier clock comparison](COMPARISON.md) additionally derives
relative rates and explicit finite uncertainty for all twelve local current
clocks, without stationary preparation or a federation. Their asymptotic
count covariances recover all three native port-distance shells. These
unit-independent relations are testable after an operational apparatus
identification; choosing seconds is a convention, not the missing result.

The [measured host-clock attachment](HOST_CLOCK.md) now records actual native
execution against a monotonic host counter, calibrates one local clock on
one block and checks a different clock on three untouched blocks. Their
duration errors are 6.29%, 2.29% and 7.76% against a declared 25% engineering
tolerance. A host-workload control preserves every native record and rejects
that calibration in all four blocks. This is an executable apparatus check;
it supplies no traceable lab calibration or identification with proper time.

The [joint clock bridge](CLOCK_BRIDGE.md) derives how incident interactions,
including zero currents, convert to local movement counts. Their stationary
ratio is `33/5`, but their exact correlations and finite uncertainty rule out
using that mean as a pathwise conversion. The bridge retains field reads and
extra feedback hits, and replays the existing host capture with its unchanged
calibration: all three baseline blocks remain inside the propagated interval,
while the corresponding workload controls reject it. This uses the same
total-two carrier as the capture, not the full binary preparation of the
field-refinement protocol.

## Why the model does not predict seconds

The source contains integer attempt counts and dimensionless probabilities,
but no quantity with units of seconds. The executable
`duration_scale_identifiability` audit assigns illustrative values of `1/100`,
`1` and `100` seconds per attempt. They leave the exact one-attempt record law
`(28/33, 5/33)` unchanged, while implying different event rates per second
(`500/33`, `5/33` and `1/660`). Thus the native law fixes records per attempt;
it does not fix an attempt duration or a laboratory event rate.

No laboratory apparatus or public dataset has been operationally matched to
this current record. Under a separately validated constant attempt-rate map,
the minimal scale calibration is one nonzero native-attempt difference paired
with its independently measured lab interval. One pair alone does not justify
that the same rate holds over other intervals. Without such an apparatus map,
fitting generic lab event data would not test this model.

## What runs and what is measured

The program prepares the declared stationary mixture `I66/66` of the native
binary total-two carrier. At each attempt it selects one of the 30 source
seams uniformly, uses the native fair placement coin, and executes the existing
rule C. The readout retains `bit = abs(j)`: one for a moving current, zero for
a no-current outcome. Seam, coin and current direction are probability-summed
inside the complete instrument. The current counter is a function of those
retained bits; the observer never needs to inspect the hidden configuration or
its token-separation class.

The [producer](clock.py) executes the actual 120 native Kraus branches through
the shared [primitive observer runtime](../primitive_observer/runtime.py).
One five-attempt history retains all 32 bit words and their complete
unnormalized 66-dimensional terminal density matrices. Every word probability
and terminal density is compared with the exact native transition calculation.
A saved two-attempt checkpoint is serialized, independently replayed and
continued to the same complete output.

First-current laws are extracted from the executed prefixes: `1`, `01`,
`001`, `0001`, and the remaining `0000` tail. Their probabilities sum to one.
The analogous next-current law after a retained first `1` is explicitly
conditional on that first event, whose probability is `5/33`. The full parent
history, including first `0`, remains present. A physical apparatus that halts
its native execution at the first current is not implemented in this packet.
The derived stopping law describes those retained prefixes; it does not treat
later fixed-horizon terminal states as the stopped states.

## Exact source reduction

The occupied ports have three separation classes on the supplied icosahedral
carrier: adjacent, distance two, and antipodal, with 30, 30 and 6 configurations.
For **every** configuration in a class, summing each marked native transition
over each destination class gives the same row. The producer checks this
strong lumping property for all 3,960 primitive branch addresses. With rows
acting on class probability vectors, the exact marked matrices are

```text
        [13/15   0     0  ]          [1/15  1/15   0  ]
T0  =   [  0    5/6    0  ]    T1 =  [1/15  1/15  1/30]
        [  0     0    5/6 ]          [  0   1/6    0  ]
```

`T0 + T1` is stochastic. Its stationary class law is
`pi = (5/11, 5/11, 1/11)`. The full word law is
`P(b1...bn) = pi T_b1 ... T_bn 1`. This is an exact relation to the original
native kernel, not a fit to measured words. It also holds for any initial
class distribution; tests execute nonstationary basis preparations separately.
The reduction concerns this classical record family. It does not reduce the
full quantum observer algebra or choose between coherent laws C and E.

Counts follow from the coefficients of `pi (T0 + z T1)^n 1`. The retained
[receipt](results.json) includes **all** count probabilities at 1, 2, 8 and
32 attempts, including zero, rather than a conditional histogram of successes.
It also retains first-current probabilities through 16 attempts and the exact
remaining no-current tail.

## New measurement predictions

| Observable and preparation | Exact prediction |
| --- | --- |
| Moving bit at a stationary attempt | `5/33` |
| First next current, sampled at an arbitrary stationary attempt | Mean waiting count `147/22` attempts |
| Interval after an actual current, in the stationary event-conditioned law | Mean `33/5` attempts; variance `951/25` attempts squared |
| Two adjacent event-conditioned waiting intervals | Covariance `9/100` attempts squared |
| Next bit is one, given last bit one | `23/150` |
| Next bit is one, given last two bits `11` | `53/345` |

The distinction between the two waiting preparations is measurable. A
stationary sample encounters the more slowly moving states more often than a
sample taken immediately after a current. The latter, conventionally called
the *Palm* preparation, has class law
`nu = pi T1 / (5/33) = (2/5, 1/2, 1/10)`. It is obtained here by conditioning
an actual retained source event, not by assuming a fresh independent clock.

No-current outcomes leave the classical configuration unchanged. The next
current therefore has a geometric waiting law conditional on that configuration.
For integer `m >= 0`, the survival laws are exactly

```text
P(W > m | stationary sample) = (5/11)(13/15)^m + (6/11)(5/6)^m
P(W > m | after a current)   = (2/5) (13/15)^m + (3/5) (5/6)^m.
```

These formulas include the no-current tail at every finite horizon. Their
geometric sums yield the means and second moments above. They are not
independent waiting intervals: after moving, the class evolves according to
the derived jump matrix

```text
J = [1/2  1/2   0 ]
    [2/5  2/5  1/5]
    [ 0    1    0 ].
```

Writing `m = (15/2, 6, 6)` for the conditional mean waits, stationarity
`nu J = nu` gives the adjacent product moment `sum(nu_i m_i J_ij m_j)`.
Subtracting `(33/5)^2` gives the nonzero covariance `9/100`. Replacing this
process with independent equal-rate ticks loses that observation. Likewise,
the two different conditional next-bit probabilities in the table rule out
a first-order Markov description on the current bit alone. The three hidden
classes retain predictive memory without adding a clock module.

## Which scheduler premises disappear

The **stationary mean rate `5/33` does not require uniform seam choice**.
For each fixed seam exactly 20 of the 66 configurations have unequal endpoint
occupancy, and the fair coin produces a move with probability one half.
Every seam thus has the same stationary moving probability `20/(66*2)`.
Its native kernel is symmetric, so every fixed normalized mixture of these
seam kernels preserves the uniform stationary law. This proves the mean
identity for arbitrary fixed, state-independent seam probabilities. It does
not cover adaptive or state-dependent scheduling.

If all seam probabilities are positive, every configuration has a positive
escape probability `r_i`. Symmetry gives the event-conditioned configuration
law `nu_i = r_i/(66*p)`, where `p=5/33`. Consequently
`sum(nu_i/r_i) = 1/p = 33/5`. The observer's mean event interval is then fixed
without choosing uniform seam weights. Positivity matters: selecting only
one seam leaves 46 configurations forever inactive. That stationary mixture
still has mean move rate `5/33`, but its event-conditioned mean interval is
`2`, not `33/5`. This exact inactive-component negative control is retained.

The full measurement law remains scheduler-sensitive. Doubling one seam's
relative weight, with weights `(2/31,1/31,...,1/31)`, leaves the stationary
mean and mean event interval unchanged but changes the observable two-current
probability from `23/990` to `1477/63426`. The stopped configuration transition
law changes too. Mean clock normalization therefore does not eliminate the
whole scheduling premise.

An additional declared idle choice with probability `1-q`, independent at
each attempt, gives `T0' = (1-q)I + q T0`, `T1' = q T1`. The stopped moving
kernel is unchanged because `T1'(i,j)/r_i' = T1(i,j)/r_i`; its stationary
event-conditioned law also stays fixed. This removes only global idling from
the stopped transition sequence. Zero records, waiting counts and attempt
costs still change. The producer checks `q=1/2` on the full 66-state kernel;
the cancellation establishes the identity for any `0<q<=1`.

## Existing theorem attachment and limits

The existing
[`SourceRecurrenceClock.returnCount_block_additive`, `returnCount_mono` and `returnCount_le`](../../Lean/QFT/SourceRecurrenceClock.lean)
apply directly to the current record: take `i=1` and the two-valued trajectory
`omega(k)=abs(j_k)`. These count theorems accept any such record sequence; they
require no Markov assumption. They give adjacent-block additivity, monotonicity
and the bound `N(n)<=n` for this instantiated observer counter.

The nearest waiting-law result in that file, `meanReturn_hasSum`, belongs to
a different committed two-state chain. Its numerical return law is not
imported here. The advance is the checked native-kernel mapping, exact marked
reduction and newly calculated waiting/count distributions compared with full
executed CQ output. These are exact finite arithmetic checks and direct
mathematical derivations; no new Lean proof is claimed.

The [count-clock enclosure](../../Lean/Time/SourceCountClockEnclosure.lean)
assumes independently justified count-volume errors and a fourth-power
proper-time/volume relation. The present current counter supplies neither.
It therefore establishes no spatial dimension, Lorentzian geometry, proper
time calibration or physical frequency. The supplied carrier, native rule,
stationary preparation and chosen current observable remain explicit.

The related [`relative_modular_dynamics`](../relative_modular_dynamics/README.md)
packet derives a periodic Y-readout response to the Palm state's modular
parameter from an existing coherent preparation. That gives a local modular
phase readout under the declared ideal phase target. Its exact boundary test
also shows why the current native rule cannot supply that phase evolution:
real coherent preparation and real repair branches leave the Y readout
unbiased at every integer native-attempt count. Native attempt counts and the
modular parameter remain distinct model readouts, with no derived conversion
between them, and neither specifies physical duration or a lab clock.

The new laws are predictions of the executable source model. No laboratory
measurement has yet been operationally identified with this native current
record, so this packet makes no new empirical agreement claim and changes no
frozen anharmonic result. Its useful next target is a justified interpretation
of a measured event counter, with the joint count and waiting laws checked
together. Matching the mean alone would miss the distinguishing correlations.

## Reproduction

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.native_clock --check model/native_clock/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/native_clock
```

`--write PATH` creates a new receipt exclusively and refuses to overwrite an
existing one. Exact probabilities are rational strings. Floating CQ residuals
are replay diagnostics with tolerance `2e-11`, not outward-rounded numerical
certificates. Public source and theorem hashes bind the retained computation.
