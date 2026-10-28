# An adjacent pair is a temporary native excitation

The two occupied ports of the binary total-two carrier define an adjacent
pair when they share one of the thirty native seams. This definition uses
the supplied carrier graph, independently of the repair law or any target
mass. It selects thirty of the sixty-six configurations. It gives a finite
relational localization test; it supplies neither physical distance nor a
topological charge.

The exact result is a temporary, mobile population cluster. Under uniform
seam selection it first leaves adjacency after **15 native repair attempts
on average**, and returns after another **18 attempts on average**. A
record-controlled selector can prolong its residence, but a positive
exploration floor bounds that lifetime. The candidate also fails the
reference-like exterior test. It therefore does not establish CP-03.

## The complete native process and first exit

[`localized.py`](localized.py) uses every raw C/E partial-injection Kraus
map from the existing native instrument. It retains both fair coins and
zero maps: 240 raw maps for C and 480 for E, of which 120 and 300 are zero
on binary66. The same population kernel results after summing the records;
the quantum instruments remain different. Every kernel row is an exact
integer row divided by 60. No selected successful branch is renormalized.

Write `A` for the adjacent-pair set. In one uniform attempt from any member:

| Outcome | Probability |
| --- | --- |
| Configuration unchanged | `13/15` |
| Pair moves but remains adjacent | `1/15` |
| Pair first separates | `1/15` |

The four within-cluster moves each replace one occupied endpoint; the pair
keeps its other endpoint. The retained thirty-state killed kernel `Q`
contains all these transitions, with `Q 1=(14/15)1`. Thus the probability
of no exit during the first `n` attempts is `(14/15)^n`, and first exit at
attempt `n>=1` has mass `(14/15)^(n-1)/15`. Its mean is 15. The mean number
of moves that retain adjacency before the first exit is one. Conditioning
on survival gives a per-attempt transport probability `1/14`; that
conditional diagnostic is not the unconditioned observer output.

The first-exit flag is **offline population-path bookkeeping**. It is not a
new apparatus measurement, quantum memory or input to the controller.
Repeated projection onto `A` would add a service and can disturb quantum
outputs; no such monitored quantum survival experiment is asserted here.

There is a separate quantum one-step statement. Since every raw map is
injective, `K* P_out K` is diagonal. Summing uniform native branches gives
`P_A Gamma*(P_out) P_A = P_A/15`. Hence every state supported on `A`, including
one entangled with an untouched reference, has one-step exit probability
`1/15`. This does not identify the killed population process with a
nondisturbing quantum history measurement.

## Return and continued population motion

The pair's graph distance has three values, with class sizes `(30,30,6)`.
Every configuration in the same distance class has the following transition
counts, so this is an exact projection of the full kernel:

```text
             next distance
              1   2   3
distance 1:  56   4   0
distance 2:   4  54   2       all divided by 60
distance 3:   0  10  50
```

An exit reaches distance two. The mean hitting times of adjacency from
distances two and three are 18 and 24 attempts, verified by the exact
hitting-time equations. The stationary adjacency mass is `5/11`, agreeing
with the mean residence/excursion ratio `15/(15+18)`. A returned pair is a
new entry; its mass cannot be counted as uninterrupted survival.

## Classical feedback and the exploration floor

At every complete committed past, let each seam have conditional
probability at least `epsilon/30`. The law may be C or E, and may depend
on the selected seam and history. Coins remain fresh and fair. Each
adjacent configuration has four exiting seam/coin addresses, giving
conditional exit hazard at least `epsilon/15`. The selector's remaining
weight contributes a nonnegative hazard. Averaging over any conditional
population state preserves the bound, including correlations with the
controller. No independence between successive selector choices is needed.

If every branch executes the stated number of repairs, with no intervening
source changes except native repair and holds:

```text
P(no exit by n actual repairs) <= (1-epsilon/15)^n,
E[first exit repair count]    <= 15/epsilon.
```

At `epsilon=1/16`, the latter bound is **240 repairs**. If each continuing
attempt repairs with conditional probability at least `r`, the same proof
gives hazard `r epsilon/15` and mean at most `15/(r epsilon)`, or **3840
attempts** for both floors `1/16`. An explicit halt does not count as further
repair or continuing attempts. Its retained mass must remain separate.
The mean first-exit bound assumes continuation until exit, with no earlier
halt. Selecting only paths that continue after an informative stopping
decision does not inherit the unconditioned survival bound.

The executed two-state controller prefers seam `(0,5)` when the previous
current was positive and `(0,1)` otherwise, with uniform exploration weight
`1/16`. It starts from the occupied core `(0,1)`, memory zero, and uses only
committed current records. C and E give identical exact population results
under this decoded-current policy. At 64 attempts the never-exited mass is
approximately `0.740465`, whereas current adjacency has mass `0.859587`;
the difference is retained returned mass. The receipt stores exact rational
values at all declared horizons, without sampling paths.

As a negative control, selecting only the occupied seam `(0,1)` freezes
that starting configuration. It violates positive exploration and records
no movement. Halting also freezes it. Neither is evidence of protection
under the declared continuing law.

## Quantum coherence and the reference exterior

The new population calculation does not replace the
[adaptive complete-record coherence theorem](../observer/RECORD_RECOVERY.md#persistent-exploration-with-adaptive-classical-control).
For every persistently exploring record-only schedule without intervening
quantum controls, that theorem covers binary66, arbitrary untouched
references and arbitrary final recovery using the complete tape. At the
same exploration floor and 512 actual repairs, it bounds encoded-qubit
Bell-state recovery error from below by more than `0.3474` for C or
arbitrary C/E switching, and more than `0.4979` for E. These are trace-distance
lower bounds. They do not assert exact finite-time entanglement breaking.
The count must hold on every compared path; an expected repair count is
insufficient. Intervening coherent clocks or active quantum correction
require a separate argument.

For the fixed occupied core `(0,1)`, the exterior is the other ten ports.
The candidate's exterior is certainly empty. The uniform attempt-stationary
reference `I66/66` has that exterior vacuum with probability only `1/66`,
so the exterior trace distance is `65/66`. The shared ideal three-charge
preparation reference has vacuum weight `2/45`, giving distance `43/45`.
Its finite preparation error remains a separate bound. Even mixing uniformly
over all adjacent cores only matches the uniform reference's single-port
means; its ten-port exterior differs by trace distance `16/33`.

This is therefore an explicit finite cluster with quantified loss and
transport, not a localized deformation with a reference-like exterior.
The graph locality is supplied; a native adjacency detector, protecting
boundary invariant, physical lifetime and defect mass are not derived.
History-based sectors, larger carriers and active recovery remain separate
contracts. The [shared diagnostic contract](campaign.json) fixes this scope.

The nearest existing result is the
[native movement census](../event_circulation/README.md), which supplies
activity `2/15` for adjacent pairs and `1/6` for separated pairs. This packet
adds the killed localization kernel, exact return law, adaptive lifetime
bound and reference-exterior test. It neither repeats the source-charge
commutant calculation nor interprets a relaxation rate as rest energy.

```sh
python -m model.clock_defect_energy.localized --check
python -m pytest -q -W error model/clock_defect_energy/test_localized.py
```
