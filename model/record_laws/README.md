# Admissible record laws for native repair

The native seam rule says where load goes when two neighbouring ports repair
their disagreement. It does not say what the event writes down. That second
choice, the record law, decides which quantum coherences survive a repair.
This package classifies the unread channels that keep the classical rule,
act on one seam and respect local load phases and endpoint symmetry. It also
specifies current-record instruments and computes their distinguishing
experimental predictions.

In short, the admissible laws form one convex family described by
five coherence multipliers. The executed instruments C and E are two points
of it. Two more points, C0 and T, complete a ladder of laws that record
progressively more. Exactly two laws in the whole family, E and T, are
idempotent, and every law that keeps an empty seam coherent with a full seam
completes to E when the same seam is repaired repeatedly. Classical current
data cannot distinguish any two members under the same current-record
interface. A single coherent interference
experiment reads one coherence multiplier directly.

The four requirements below classify unread channels. A complete recorded
instrument also specifies how the environment is read: channel covariance
guarantees a definite-current Kraus representation exists, but does not force
every possible recorded branch to have a definite current. The
[instrument control](INSTRUMENT_CONTRACT.md) checks this distinction exactly
without changing the executed C/E instruments or their predictions.

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
survives. Requirement 3 is an additional phase-covariance condition; conservation
of the classical seam total alone does not imply it. A law that records the pre-repair load at one endpoint,
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

**Classical current records.** Equip each channel with its definite-current
instrument, whose branch diagonals are the specified `P(j|κ)`. Every protocol
that prepares load configurations, reads these currents and feeds them back
then assigns the same probabilities under every law. Unread diagonal-channel
equality alone does not establish equivalence for arbitrary record bases or
additional labels; see the [instrument contract](INSTRUMENT_CONTRACT.md). The
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

## Preparation by finite native warm-up

[burn_in.py](burn_in.py) extends the existing
[C/E mixing theorem](../observer/QUANTUM_MIXING.md) to **every admissible
coherence table**, including complex phases, different laws on different
seams, and the restricted law feedback described below. The uniform binary
density `tau=P_binary/66` is the unique stationary source state of each fixed
uniformly scheduled law. Their quantum spectral gap is the same exact
`(5-sqrt(5))/60`. Equilibrium alone therefore selects none of the nine real
record-law coordinates.

More usefully, an arbitrary initial total-two density, even entangled with
an untouched reference, approaches `tau` with an explicit error. A finite
warm-up replaces an exact stationary-background preparation with that error
budget. It does not prepare the coherent probes, blank records, programme
states or source apparatus required by other experiments.

**Domination on the actual source.** Positivity of every current block gives
`|G_j(a,b)| <= sqrt(P(j|a)P(j|b))`. The right side is C0's coefficient.
Consequently C0 bounds the absolute value of every admissible channel entry.
For a binary matrix unit with signed difference `z=x-y`, let `c1` and `c2`
count seams with endpoint differences of magnitude one and two. C0's block is

```text
B_z = [1-((1-1/sqrt(2))*c1+c2/2)/30] I - L_excl(G[z=0])/60.
```

Each local block is a positive self-adjoint contraction. Its row and column
sums coincide. The exact native enumeration verifies this reduction on all
4,356 binary matrix units, symmetry, and the largest nonzero-sector sum

```text
q0 = 43/60 + 2 sqrt(2)/15 < 57/60 < r=(55+sqrt(5))/60.
```

The bound is attained for the 60 oriented adjacent-site differences. The
diagonal sector is exactly the existing exclusion chain. Its gap is the
source graph gap divided by 60, using the weighted exclusion theorem of
[Caputo, Liggett and Richthammer, §4.1.1, (4.1)](https://arxiv.org/pdf/0906.1238).
Thus every binary channel contracts the nonconstant diagonal sector by `r`
and every off-diagonal sector by at most `q0`, regardless of its phases.
Entrywise domination implies the operator 2-norm bound because the dominating
matrix has both row and column sums bounded by `q0`; normality of the admitted
channel is unnecessary. The surviving classical gap eigenmode makes the
exact binary traceless Hilbert--Schmidt contraction factor `r`. Here the
spectral gap means one minus the largest modulus of a nonstationary
eigenvalue. On all 6,084 full-carrier units, the
collision-to-collision block is diagonal and bounded in modulus by `5/6`.
The triangular binary/collision decomposition therefore gives the same
unique stationary density and spectral gap on dimension 78; it does not
give a Hilbert--Schmidt contraction bound for that nonunital full channel.

**Schedules and restricted feedback.** Uniform scheduling is sufficient,
but a weaker condition gives a quantitative result. At each attempt the seam
probabilities may be an exogenous, time-varying distribution with
`p_e >= alpha/30`, for a common `0<alpha<=1`. Decompose that distribution
into an `alpha` uniform part and a remaining distribution. Each binary
single-seam C0 map is a self-adjoint positive contraction; each classical
single-seam map fixes the uniform density. Their nonconstant contraction
factors are bounded by

```text
r_alpha = 1-alpha*(1-r),      mu_alpha = 1-alpha/6.
```

The seam choice must be independent of the source and preceding repair
outcomes. The law table may depend on the seam, attempt and preceding
**coarse current records**, with no intervening operation on the source.
For each complete current path, the absolute product of coefficients is
bounded by the corresponding product of C0 coefficients. Summing over all
paths bounds the resulting off-diagonal matrix entrywise by the product of
the exogenous C0 matrices. All diagonal path probabilities are independent
of the chosen tables, giving the same product of classical kernels.
This proves the contraction also for this particular law-only feedback.
It does not follow merely from each unread channel fixing `tau`. A changed
record basis, coherent access to the environment, other source operations,
or outcome-dependent seam scheduling requires a separate argument.

**Untouched references and full outputs.** For any initial joint density
`rho_SR` supported on the binary source, the resulting source channel `Psi_n`
obeys

```text
dist((Psi_n tensor id_R)(rho_SR), tau tensor rho_R)
    <= min(1, 33*r_alpha^n),                 dist=half trace norm.
```

The replacement map is the orthogonal Hilbert--Schmidt projection onto the
identity source operator. The preceding sector bounds give
`||Psi_n-Replacement||_(2->2) <= r_alpha^n`. A pure input has Schmidt rank
at most 66, so its output difference has dimension at most `66^2` and
`||Y||_1 <= 66||Y||_2`. Tensoring with the reference leaves the 2-norm bound
unchanged. The input Hilbert--Schmidt norm is one. Convexity proves the stated
bound for arbitrary mixed inputs and reference dimension.

For all 78 source states, choose an even integer `0<=s<=n`. After `s`
attempts the collision population is at most `mu_alpha^s`. In the proof,
replace the collision block by `tau` times its reference marginal, and
remove binary/collision cross terms. Positivity bounds those cross terms
by the square root of the collision probability; the replaced population
costs at most that probability. Keep any classical controller prefix during
this mathematical comparison. For prefix probabilities `w_h` and conditional
collision masses `q_h`, Cauchy--Schwarz gives
`sum_h w_h sqrt(q_h) <= sqrt(sum_h w_h q_h) <= mu_alpha^(s/2)`.
Conditional on each prefix the remaining
law-only protocol satisfies the same binary bound, and the untouched
reference marginals sum to the original marginal. Contractivity then gives

```text
dist(output_SR, tau tensor rho_R)
 <= min(1, mu_alpha^(s/2)+mu_alpha^s+33*r_alpha^(n-s)).
```

This split introduces **no executed measurement or reset**. The evaluator
uses the strict rational upper bound `r < 19079/20000`, verified by integer
squaring, and optimizes only over even proof splits. `coverage_of(weights)`
computes `alpha=30 min_e p_e`; it does not infer independence or certify
sampling error for estimated weights. For uniform scheduling, the retained
certificates bound binary warm-up error by `1/1000` after 221 attempts and
full-carrier error by `1/1000` after 318 attempts, with proof split `s=88`.
These are sufficient
counts, not sharp mixing times or durations in seconds.

The warm-up tape and its controller are traced out at the comparison
boundary and are unavailable to the subsequent protocol. They may remain
elsewhere as inaccessible records; this requires no physical erasure. They
need not become independent of the source: even starting
from `tau`, after repairing `(0,1)`, current `+1` has probability `5/66`
and excludes an occupied endpoint zero. A false source/tape product would
assign that joint event probability `5/396`. An exact second control changes
the seam according to the first current and moves the source away from
`tau`. A one-seam schedule also retains an exact empty/full coherence.
These controls expose residual premises without changing the native law.

Any **common complete subsequent CPTP instrument** inherits the same error,
including its joint clock, field, quantum-reference and timeout outputs,
independently of its finite horizon. For full-carrier warm-up it must be
defined on all 78 states, for example by retaining an explicit out-of-domain
outcome; a binary-only map is insufficient. Keeping or conditioning on the
warm-up tape, or conditioning on a rare future success, is outside this
unconditional bound. The two-attempt executable witness retains all records
and complex reference coherences during law-only feedback; an independent
floor/ceiling Kraus oracle checks every output before the tape is discarded.
The result derives a preparation error budget, not a physical record-law
selection, apparatus implementation, Lorentzian dynamics or interacting QFT.

## Scope

The classification concerns one seam of the total-two carrier. Requiring
carrier covariance also fixes the same law on every seam; the warm-up bound
allows seam-dependent laws within the same local family. The
[coherent continuation](../observer/COHERENT_CONTINUATION.md) of unread
environments and the response extensions in the
[candidate guide](../CANDIDATES.md) have separate contracts. No laboratory
data enter this classification or warm-up calculation.

## Run

From the repository root:

```sh
python -m model.record_laws --check model/record_laws/results.json
python -m model.record_laws.burn_in --check model/record_laws/burn_in_results.json
python -m pytest -q -W error model/record_laws
```

The replay takes about half a minute with exact arithmetic in
`Q(√2)(i)`. `--write <path>` writes a fresh receipt and refuses to overwrite.
