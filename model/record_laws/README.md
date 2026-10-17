# Admissible record laws for native repair

The native seam rule says where load goes when two neighbouring ports repair
their disagreement. It does not say what the event writes down. That second
choice, the record law, decides which quantum coherences survive a repair.
This package classifies every record law that keeps the classical rule, acts
on one seam and respects local load phases, and it computes what each law
predicts for the experiments that can tell them apart.

In short, the admissible laws form one convex family described by
five coherence multipliers. The executed instruments C and E are two points
of it. Two more points, C0 and T, complete a ladder of laws that record
progressively more. Exactly two laws in the whole family, E and T, are
idempotent, and every law that keeps an empty seam coherent with a full seam
completes to E when the same seam is repaired repeatedly. Classical record
data cannot distinguish any two members. A single coherent interference
experiment reads one coherence multiplier directly.

## Setup

Use the [native carrier](../observer/README.md): twelve ports, thirty seams
and the 78 integer load configurations of total two. For the seam `e=(a,b)`
write the seam key `κ(x)=(x_a,x_b)` and let `D` move one load from `a` to
`b`. A fair coin `c` breaks ties:

```text
d = x_a - x_b,     j(x,c) = [d - (2c-1)(|d| mod 2)]/2,     x -> x + j D.
```

The current `j` is the transported load. On the total-two carrier the seam
keys and their current probabilities `P(j|κ)` are:

| Seam key | Meaning | `P(-1)` | `P(0)` | `P(+1)` |
| --- | --- | ---: | ---: | ---: |
| `(0,0)` | empty seam | | 1 | |
| `(1,1)` | full seam | | 1 | |
| `(1,0)`, `(0,1)` | one load, which waits or moves | `0` or `1/2` | `1/2` | `1/2` or `0` |
| `(2,0)`, `(0,2)` | collision, which splits | `0` or `1` | | `1` or `0` |

## Theorem 1: every admissible law is a coherence table

Let `Φ` be a channel on the carrier for one seam. Require:

1. **Seam locality.** `Φ` acts on the endpoint loads and as the identity on
   the other ten ports.
2. **Exact classical law.** `Φ(|x><x|)` is the native repair distribution.
3. **Local number-phase covariance.** `Φ` commutes with
   `U_t = exp(i Σ_p t_p N_p)` for every real vector `t`.
4. **Endpoint-swap covariance.** The carrier rotation that exchanges `a`
   and `b` maps the law to itself with the current reversed.

Then there are matrices `G_j`, one for each current, such that

```text
Φ(|x><x'|) = Σ_j G_j(κ(x),κ(x')) |x + jD><x' + jD|,
```

each `G_j` is positive semidefinite on the seam keys that can carry current
`j`, and its diagonal is `G_j(κ,κ) = P(j|κ)`. Conversely every such family of
matrices defines an admissible channel.

*Proof.* Covariance gives
`Φ(|x><x'|) e^{i t·(x-x')} = U_t Φ(|x><x'|) U_t†`, so `Φ(|x><x'|)` is a
combination of units `|y><y'|` with `y-y' = x-x'`. Complete positivity puts
every nonzero Choi entry between pairs whose diagonal entries are nonzero,
and the classical law makes those pairs `(x, x+jD)`. Then `y-y'=x-x'` forces
the same current on both sides. The Choi matrix is therefore the direct sum
of the `G_j`, which is positive exactly when every `G_j` is. Trace
preservation and the classical law fix the diagonal. Locality makes `G_j`
a function of the seam keys alone. ∎

In words: a record law is a table that says, for each pair of seam
configurations that move the same way, how much of their quantum coherence
survives. Requirement 3 is the quantum form of the conservation the classical
rule obeys. A law that records the pre-repair load at one endpoint,
`(c, x_a)`, keeps the classical rule exactly and fails it: one of its Kraus
branches carries two different currents, and the receipt retains an explicit
matrix unit whose load difference changes.

## The five multipliers

Endpoint-swap covariance leaves five free multipliers on the total-two
carrier:

| Multiplier | Coherence it keeps | Bound |
| --- | --- | --- |
| `k = G_0((0,0),(1,1))` | empty seam with full seam | `\|k\| ≤ 1` |
| `g = G_0((0,0),(0,1))` | empty seam with a waiting single load | `\|g\| ≤ 1/√2` |
| `h = G_0((1,1),(0,1))` | full seam with a waiting single load | `\|h\| ≤ 1/√2` |
| `m = G_0((1,0),(0,1))` | the two waiting orientations, real | `\|m\| ≤ 1/2` |
| `n = G_-1((0,1),(0,2))` | a moving single load with a moving collision | `\|n\| ≤ 1/√2` |

The listed bounds are the two-by-two conditions. The full condition is that
the four-by-four matrix of the current-zero block is positive. One
consequence is used below: if `k = 1`, then `g = h`, because the vector that
subtracts the full-seam column from the empty-seam column lies in the kernel.

A law that also records the coin is a sum over the two coins of positive
blocks, each carrying weight `1/2`. For such laws `|g|`, `|h|` and `|n|` are
at most `1/2` and `m = 0`, since no coin lets both waiting orientations
carry the current zero.

## Four named laws and the record ladder

| Law | Each event records | Record values per seam | `k` | `g = h` | `m` | `n` |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| C0 | the current | 3 | 1 | `1/√2` | `1/2` | `1/√2` |
| C | coin and current | 6 | 1 | `1/2` | 0 | `1/2` |
| E | coin and pre-repair contrast | 10 | 1 | 0 | 0 | 0 |
| T | coin and both endpoint loads, that is coin, current and seam total | 12 | 0 | 0 | 0 | 0 |

C and E coincide branch by branch with the instruments executed by the
[native observer](../observer/README.md) and the
[prediction package](../predictions/README.md). C0 uses one Kraus map per
current with amplitude `√P(j|κ)`: it records the transported load and
nothing else.

The four laws form a ladder. Each step records one more classical quantity
after the repair, through a diagonal readout of the output:

```text
C = C0 followed by a readout of the coin,
E = C followed by a readout of the seam parity,
T = C followed by a readout of the seam total.
```

The receipt checks all three as exact operator identities on the 78-state
carrier. The coin readout has weight `1/2` per coin on a balanced output and
a sharp label on an odd one. Each step removes coherence that the previous
law kept.

Phases are part of the family. The law that gives C's waiting branches a
minus sign has the same classical law and records as C, with `g = -1/2`.
A quarter-turn variant of C0 has `g = -i/√2`.

Among the ten coin-recording laws whose records refine coin and current and
respect the swap, all ten channels are distinct, four share C's
interference values and six share E's. Only E and T are idempotent.

## Theorem 2: repetition and completion

Repairing the same seam twice composes tables entrywise:

```text
(k,g,h,m,n) followed by (k',g',h',m',n') = (k k', g g', h h', m m', n conj(h')).
```

Every coherence keeps its seam key under one more repair, so each multiplier
is applied again. The receipt checks the law for 64 ordered pairs of tables
on every pair of seam keys and for two pairs on all 6,084 matrix units.

*Consequences.*

- A law is idempotent exactly when `k,g,h,m ∈ {0,1}` and `n = n conj(h)`.
  The bounds then force `g = h = m = n = 0`, leaving `k = 1`, which is E,
  or `k = 0`, which is T.
- After `r` repairs of one seam the multipliers are `k^r, g^r, h^r, m^r` and
  `n conj(h)^(r-1)`. Every law with `k = 1` converges to E. Every law with
  `|k| < 1` converges to T. A law with `|k| = 1` and `k ≠ 1` keeps a rotating
  phase and has no limit.
- For C this is the identity `C^r = E + 2^(1-r)(C - E)`. For C0 the waiting
  coherence decays by `1/√2` per repair instead of `1/2`.

In words: E is the coarsest law for which a single repair finishes the
job, and it is where every law that keeps the empty and full seam coherent
ends up when a seam is repaired often enough.

## Theorem 3: what observations read

**Classical records.** Every member has the same action on diagonal states,
so every protocol that prepares load configurations, reads currents and
feeds them back assigns the same probabilities under every law. The
repository's [observability results](../identifiability_checks/README.md)
extend this to record-controlled quantum work for C and E.

**One coherent event.** Prepare `(|X> + e^{iθ}|Z>)/√2`, where `X` has loads
at ports 2 and 3 and `Z` at ports 1 and 2, repair seam `(0,1)` once and
measure the projector onto `(|X> + |Z>)/√2`. Summed over all records,

```text
Pr(plus) = 3/8 + Re(g e^{-iθ})/2.
```

The empty seam `X` stays put; `Z` waits with probability `1/2` and keeps
coherence `g` with `X` while it waits. The phase settings `0`, `π` and `π/2`
read `Re g` and `Im g`:

| Law | `θ = 0` | `θ = π` | `θ = π/2` |
| --- | ---: | ---: | ---: |
| C0 | `3/8 + √2/4 ≈ 0.729` | `3/8 - √2/4 ≈ 0.021` | `3/8` |
| C | `5/8` | `1/8` | `3/8` |
| E, T | `3/8` | `3/8` | `3/8` |
| any admissible law | `3/8 + Re(g)/2` | `3/8 - Re(g)/2` | `3/8 + Im(g)/2` |

This is the experiment of the [fixed prediction protocol](../predictions/README.md),
whose locked decision rule compares C with E. Read as an estimator, the same
counts measure the coherence `g` of whichever law nature uses.

**Empty against full.** With `W` the full seam, ports 0 and 1, the same
construction gives `Pr(plus) = 1/2 + Re(k e^{-iθ})/2`. C0, C and E predict
certainty at `θ = 0`; T predicts `1/2`. This second test separates the two
idempotent laws.

## Scope

The classification concerns one seam of the total-two carrier and assumes the
same law on every seam, which carrier covariance requires. Network coherence
rates beyond C and E, the [coherent continuation](../observer/COHERENT_CONTINUATION.md)
of unread records and the response extensions in the
[candidate guide](../CANDIDATES.md) are separate questions. No laboratory
data enter.

## Run

From the repository root:

```sh
python -m model.record_laws --check model/record_laws/results.json
python -m pytest -q -W error model/record_laws
```

The replay takes about half a minute with exact arithmetic in
`Q(√2)(i)`. `--write <path>` writes a fresh receipt and refuses to overwrite.
