# Space and time from one carrier: exact passage records and the physical gate

**Completed:** a finite signal protocol on the existing twelve-port native
carrier, exact first-arrival laws in a receiver's own incident-attempt clock,
and a metric inferred from those laws. No observer federation is required.
**Not established:** physical three-dimensional space, Lorentzian spacetime,
a spatial continuum, or a relation to a laboratory clock. The receipt sets
these acceptance fields to false even when every mathematical check passes.

This packet makes a concrete connection between geometry, propagation and a
local clock instead of identifying their parameters by name. It also shows
why the previously established rank-three response Gram and native causal
order do not yet establish physical spacetime. A local clock that ignores
remote attempts does not inherit the one-edge-per-attempt propagation bound.

```sh
python -m model.physical_space_time --check model/physical_space_time/results.json
python -m pytest -q -W error model/physical_space_time
```

## Contract and relation to existing results

The carrier's seam support is reconstructed by the existing
[current-record calibration](../relative_observer_geometry/README.md).
That calibration supplies labelled known preparations, seam addressing and
faithful current records. Removing seam-exposure weights gives its unweighted
Laplacian `L`. Neither coordinates nor a metric array enter the estimator.
The carrier incidence is still an input to the source, disclosed operationally
through those seam records; reconstruction does not derive its architecture
from self-reference.

The new passage experiment prepares **one** native unit at sender `u` and
zero elsewhere. This is an explicitly changed preparation of the same native
integer-load rule, whose total-one sector is invariant. It is not the retained
66-state binary total-two experiment. With two indistinguishable units a
receiver's first arrival does not identify which preparation site supplied
it. The total-one calibration avoids inventing a tag or giving the receiver
access to hidden source state. Its mathematical validity does not establish
a physical service capable of loading that preparation.

Select each of the thirty seams with probability `1/30`, independently of the
past, and draw each placement coin fairly conditional on the selected seam
and past. The receiver at `v` knows its initially zero occupation and reads
only its incident seam labels and signed currents. Every incident attempt,
including zero current, advances its clock. The field readback is

`n_v <- n_v + j (1[v=b] - 1[v=a])`.

The first zero-to-one change is the arrival record. Remote repair attempts
continue while the receiver waits and are not clock ticks. Known preparation,
repeated resets, independent scheduling/coins, and faithful local records
are supplied apparatus services. There is no end-to-end physical construction
or synchronization certificate for them in this packet. Exact ensemble laws
are computed; finite sampled data do not become exact expectations.

The compatibility contract is independent of a solver: retained arrival bins
and timeout must have total probability one; decoded occupations must agree
with the local current equation; all clock laws must come from the same
native branch histories. The code checks these before accepting a receipt.

## Derivation and exact result

With a single unit the unchanged native rule either leaves the unit in place
or moves it across the selected seam, each with probability one half. Thus
the **native-attempt** transition is exactly

`P = I - L/60`.

This identity is independently checked by executing all 720 native branches,
then comparing with the Laplacian reconstructed from calibration records.
For receiver `v`, split `P=A_v+B_v` into remote and incident attempts. The
column sums are `5/6` and `1/6`, respectively, independently of the state.
Consequently the sum over every unrecorded remote word converges absolutely:

`K_v = sum_(r>=0) B_v A_v^r = B_v (I-A_v)^(-1)`.

This is the completed transition for one receiver incident tick. It is the
same stopped-clock construction already used for the
[two-port total-two field law](../observer_time_fields/README.md), applied
here to a different, disclosed signal preparation and one receiver. No new
general Markov-chain theorem is claimed.

Set the arrival mean to zero at `v`. On the remaining ports the exact
first-step equations are `h=1+P^T h` for native attempts and `g=1+K_v^T g`
for receiver incident ticks. The equivalent first-step reward for counting
incident attempts is `g=(1/6)+P^T g`; uniqueness gives `g=h/6`. This is an
expectation identity, not equality of realized durations. Exact solution
for all 132 ordered distinct pairs yields:

| Native graph distance | Unordered pairs | Mean native attempts | Mean receiver incident ticks | Probability of arrival at first receiver tick |
| --- | ---: | ---: | ---: | ---: |
| 1 | 30 | 132 | 22 | `19/222` |
| 2 | 30 | 168 | 28 | `1/74` |
| 3 | 6 | 180 | 30 | `1/222` |

The completed first-arrival law is computed with an absorbing receiver,
not from unrestricted endpoint probabilities. The receipt keeps all eight
arrival bins and the full timeout probability for representatives of every
shell. No branch is postselected or renormalized. From every unhit state,
the next receiver tick reaches the receiver with probability at least
`1/222`. Hence `Pr(arrival tick > n) <= (221/222)^n`, proving almost-sure
arrival and finite mean. In native-attempt units, a shortest path uses at
most three favorable seam/coin branches, giving
`Pr(arrival attempt > 3n) <= (1-1/216000)^n`.

The matrix `H_uv=g_v(u)` is symmetric, positive off its diagonal and obeys
every triangle inequality. More precisely,

`H_uv = 60 (L^+_uu + L^+_vv - 2 L^+_uv)`.

Here `L^+` is the inverse of `L` on the constant-free subspace. The centered
Gram of the squared-distance matrix `H` is `60 L^+`, positive of rank eleven.
Thus `H` is a resistance metric and `sqrt(H)` is a Euclidean distance with
embedding dimension eleven. This is the classical hitting/commute-time and
resistance relation; see Aldous and Fill,
[Chapter 3, Corollary 3.11](https://www.stat.berkeley.edu/~aldous/RWG/Book_Ralph/Ch3.S3.html).
Here all coefficients and hypotheses are independently checked on the
actual source and its receiver-clock stopping rule.

The contribution beyond the existing inverse-response matrix is its
explicit native **arrival-record protocol**, local-clock interpretation,
complete stopped laws and stopping-error control. Its rank eleven and the
existing slow-band rank three refer to different operational dictionaries
on the same carrier. Neither rank alone selects the dimension of physical
space. A three-dimensional effective description could still be valid for
an appropriately justified restricted observable regime.

## A clock-dependent causal obstruction

In native-attempt count the exact support satisfies
`(P^n)_vu > 0` iff `dist(u,v) <= n`, with waiting allowed. The underlying
fixed-word dependence has the stronger time-ordered-path bound already
proved in [native causal cones](../common_source/CAUSALITY.md).

In receiver incident ticks, the arrival probability from the opposite port
at distance three is already `1/222` at tick one. There is no violation of
native locality: remote path edges occurred before the first recorded
receiver event. They are absent from this clock's count.

This distinction also has an exact witness **within the established binary
total-two sector**. Compare occupied ports `{1,2}` and `{1,3}`, whose initial
difference support is distance at least two from receiver zero. Apply seams
`(2,10),(0,10)` with coins `0,1`. The respective receiver traces in native
attempt count are `(0,0,1)` and `(0,0,0)`. Each receiver has counted just one
incident attempt. Reversing these two seam events eliminates the signal.

Therefore identifying one receiver incident tick with one native causal
layer is false. This is a scoped counterexample, not an impossibility theorem
for other clocks or physical spacetime. A fixed finite carrier admits trivial
diameter-sized bounds. On a growing family with arbitrarily long paths of
remote edges before a receiver edge, the same construction rules out a
size-independent distance bound per receiver tick unless additional clock or
schedule constraints are supplied.

Counting only receiver **movements** is still a different clock: before
first arrival the receiver is empty, so first arrival is always movement
tick one for every distinct sender. It loses all passage-distance information
in this experiment. Conversely, independent identity holding slots leave
`K_v` and all receiver-clock laws unchanged while multiplying expected slot
counts by the reciprocal acceptance probability. The receipt checks
acceptance `2/7` and multiplier `7/2`. These facts preclude silently
interchanging attempt, incident and movement parameters.

## Acceptance and the federation question

The finite source now supplies a defensible positive connection: native
propagation, current-defined arrival records, a local clock, and an inferred
finite metric can be derived together. This does not require federation.
Single-carrier action and clock work can likewise be pursued without one.

What remains for physical space and time is an operational description of
the relevant *joint* measurements, with a spatial extension or other
controlled observable limit; its dimensional and propagation regime;
compatible clock comparisons; and independently tested physical readouts.
A federation is one possible extension architecture, not a theorem premise
forced by the finite results. The present finite carrier does not by itself
supply the required continuum, and mere addition of carriers does not prove
it. Rank-three geometry plus a positive transfer logarithm is insufficient
to discharge those obligations.

No absolute choice of metre or second is needed for that work: changing
units is harmless. The missing assertion concerns physical relationships
among signal records, clocks and other measurements, including their
uncertainties. Choosing a unit conversion cannot establish those relations.
The current packet establishes conditional model laws, without claiming
that native loads are laboratory particles or that native clock ticks are
seconds. The supplied preparation and schedule must remain visible in any
later physical identification.
