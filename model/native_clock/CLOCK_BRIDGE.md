# Joint incident-hit and movement clocks

An attempted interaction and a moving current are different observations.
This packet derives their **joint probability law** from the native rule,
including the zero-current interactions between movements. It also applies
the resulting uncertainty to the [retained host-clock records](HOST_CLOCK.md),
without capturing new timestamps or changing their calibration or tolerance.

The source is one 12-port carrier, 30 uniformly selected seams, the independent
fair placement coin, and the binary total-two sector with 66 configurations.
For one port `v`, define:

- `I_v`: attempted native interactions incident to that port, including zeros;
- `K_v`: its actual occupation changes, recorded as moving currents;
- `n`: all completed movements anywhere on the carrier.

These are record counters, not external durations. The supplied source and
scheduler remain inputs. The finite uncertainty below requires neither a
stationary initial configuration nor an additional clock process.

## Full clock conversion, with field reads and feedback

Partition the existing native transition matrix into `A+B0+B1`: remote
attempts, incident attempts with zero current, and incident attempts with a
moving current. Every row of `A` sums to `5/6`. With row distributions acting
on the left, the next-incident-hit generating kernel is exactly

```text
Q(z) = (I-A)^(-1) (B0+z B1).
```

This sums all remote attempts before the next local hit. At a finite remote
cutoff `r`, the unprocessed probability is exactly `(5/6)^r`; it is not
discarded or renormalized. The producer evaluates the resolvent in rational
arithmetic and retains joint probabilities of current count and two local
occupation reads. For the record-controlled protocol, a first field bit
`c=1` causes one extra repair on an existing incident seam. Its fair coin,
current and extra hit are all included: two scheduled hits mean `2+c` actual
incident hits, while both field reads remain in the record.

At stationary incident sampling, the chance of a moving local tick is
`5/33`. Two successive local hits both tick with probability `325/7326`,
whereas independent ticks would give `25/1089`. Equivalently, the next hit
ticks with probability `65/222` after a tick, rather than `5/33`. A constant
conversion factor therefore misses measurable correlations.

There is also an exact generating kernel at global movement boundaries.
In configuration `i`, let `p_ij` be the per-attempt moving transition
probability, `q0_i` the remote zero-current probability, and `q1_i` the
incident zero-current probability. A movement `i -> j` has local tick bit
`k_ij`. The complete joint kernel is

```text
M(z,y)[i,j] = p_ij (zy)^k_ij / (1-q0_i-q1_i z).
E[z^I_v(n) y^K_v(n)] = initial_row M(z,y)^n 1.
```

The denominator includes the infinite waiting tail. The executable evaluates
this formula and checks it against the unchanged native movement kernel at
`z=y=1`. Tests compare both stopped descriptions with independently propagated
native branches and retained omitted mass.

## A finite bound from any initial configuration

Before the next global movement, the configuration is fixed. Set
`p_i=sum_j p_ij`. The number `G` of incident zero-current hits before that
movement is geometric on the nonnegative integers:

```text
P(G=g | i) = p_i/(p_i+q1_i) [q1_i/(p_i+q1_i)]^g,
E(G | i) = q1_i/p_i,
Var(G | i) = q1_i(p_i+q1_i)/p_i^2.
```

Conditional on the starting configuration, `G` is independent of the next
moving destination. Thus an increment of `Z=I_v-(33/5)K_v` is
`G-(28/5)k`. It remains random after choosing units.

The source-derived movement kernel `J` has Palm stationary law proportional
to its escape probability. Its mean increment of `Z` is zero. Solving
`(I-J)h=E(Delta Z | configuration)` with Palm mean zero gives an exact
martingale `Z_n+h(X_n)-h(X_0)`. All 66 conditional drift and variance
identities are checked. Their constants are

```text
R = range(h) = 1254/125,
V = maximum conditional martingale variance = 29473097353/1176125000.

P(|I_v(n)-(33/5)K_v(n)| >= d) <= min(1, n V/(d-R)^2),  d>R.
```

The bound holds for arbitrary initial distributions, including the conditional
configuration at a preceding movement stopping boundary. It does not assume
that each interval starts afresh in Palm equilibrium. The finite stationary
mean of `Z` is zero; a nonstationary finite mean can differ from zero by at
most `R`. Bounded conditional second moments and the martingale strong law
give `Z_n/n -> 0` almost surely. Combined with the existing `K_v/n -> 1/6`,
this gives `I_v/K_v -> 33/5`. It does not give pathwise equality at finite `n`.

For every initial distribution, the asymptotic variance
`lim Var(Z_n)/n = 43736681/3637500` follows from the finite irreducible movement
chain, the Palm average of the martingale variance and the bounded corrector.
No Gaussian finite-sample approximation is used for the displayed tail bound.
Inserted feedback or changed scheduler weights require their own bound; the
feedback count law above is not silently substituted into this certificate.

## Reanalysis of measured host records

The existing capture stops each block at 10,000 global movements. Its first
baseline block calibrates port 2. The other three blocks remain excluded from
that fit and compare port 3. Replaying every actual native record gives:

| Block | Incident hits | Zero-current hits | Moving local ticks | `I-(33/5)K` |
| --- | ---: | ---: | ---: | ---: |
| 1 | 11168 | 9564 | 1604 | 581.6 |
| 2 | 11060 | 9366 | 1694 | -120.4 |
| 3 | 11160 | 9568 | 1592 | 652.8 |

An integer radius `d=3889`, computed from the bound at `n=10000`, gives a
three-block union bound of `1/20` under the ideal independent scheduler.
This conservative bound was obtained algebraically; it was not fitted to the
counts. The retained pseudorandom stream has no independently certified IID
entropy, so these measured comparisons are diagnostics, not a 95% physical
confidence statement. This incident-count analysis is retrospective.

Let `tau` be the unchanged first-block nanoseconds per port-2 movement and
`delta=1/4` its existing host-duration tolerance. The native count interval
and that declared engineering tolerance give the interval

```text
[(1-delta) tau max(0,(I-d)/a), (1+delta) tau (I+d)/a], a=33/5.
```

All three baseline host durations lie inside; all three corresponding
host-workload control durations lie outside. The control has exactly the
same incident and current records, and is not an independent native sample.
The interval combines an ideal-law count bound and an engineering tolerance;
it is not a calibrated measurement-accuracy interval. Native clock agreement
cannot force an implementation's host speed to remain constant.

The [field refinement protocol](../observer_time_fields/REFINEMENT.md) also
counts incident hits, preserving its feedback hit. This supplies the explicit
marked-clock construction needed to compare such readouts, but uses a
different region and preparation: here a single-port star in the total-two
sector, there a 12-port root region with product or canonical-sector
preparation. The
number `33/5`, retained host calibration and error certificate cannot be
imported into that protocol. Its positive effective Euclidean transfer
Hamiltonian has not thereby acquired a physical or modular clock.

No laboratory apparatus has been identified with these native records. A
laboratory attachment still needs paired incident/current records under its
declared preparation, detection and scheduler, plus independent timestamps
and quantified missed-event and clock errors. This packet derives the joint
conversion law and checks an existing software apparatus; it does not derive
proper time, modular phase motion or empirical OPH physics.

## Reproduction

```bash
python3 -m model.native_clock.clock_bridge --check model/native_clock/clock_bridge_results.json
python3 -m pytest -q -W error model/native_clock/test_clock_bridge.py
```

`--write PATH` creates a receipt without overwriting one. Replay checks the
retained measurements and their analysis; it does not reproduce timestamps.
