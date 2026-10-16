# Joint time and field laws on matched source screens

The positive incident-clock transfer in [the field construction](README.md)
and the [stopped-clock screen comparison](../native_action/INCIDENT_REFINEMENT.md)
apply to one common finite-density preparation and one common field protocol.
The resulting comparison concerns a complete four-bit record distribution,
including record-controlled native feedback. It does not infer a joint law
from a matching gap or covariance.

The preparation is an independently supplied Bernoulli product state at
density `rho=1/3`. This extends the positive-transfer dictionary beyond the
earlier fixed-total-two executions. On the two full screens it gives exactly
matching restrictions to the common ball and nonzero local fluctuations.
No physical duration per repair, fitted Hamiltonian, or independent-particle
approximation is used.

## One executable observer domain

This is a diagonal recorded branch of the
[canonical primitive observer](../primitive_observer/README.md). C and E
have the same diagonal native transition law. Their coherent continuations
are not identified by this experiment.

| Interface | Instantiation |
| --- | --- |
| State | Binary occupations on the actual full source screen. The certified common-ball calculation retains the exact exclusion action on field polynomials. |
| Preparation | Independent Bernoulli occupations of density `rho`, with occupations at root ports zero and one initially disclosed as `a,b`. Both disclosures are retained, including every label. |
| Native operation | Uniformly select a source seam and execute the existing fair-coin native repair. The unread diagonal map on a seam `e` is `Q_e=(I+SWAP_e)/2`. |
| Observer region and clock | The original twelve-port carrier is the clock region. Every attempted seam meeting it is a tick, including attempts with zero current. |
| Field access | Two counters begin at the disclosed occupations. They update using the authentic signed currents on incident seams. Subsequent field bits are processing of these records, with no new occupation measurement at each tick. |
| Read protocol | Retain `c=n_0` after the first scheduled incident hit and `d=n_1` after the `k`th scheduled incident hit, with `k>=2`. The success record is `(a,b,c,d)`. |
| Feedback | If enabled and `c=1`, apply one additional native repair on the internal root seam `(0,1)` immediately after the first read. The policy uses the committed read bit. |
| Clock accounting | Feedback is itself one additional incident tick. The final read is at total clock value `k+c` with feedback and `k` without it. The first read precedes the intervention. |
| Continuation | Remote evolution, field-conditioned state and the result of feedback all continue into the next interval. There is no reset or fresh preparation at a read. |
| Finite approximation | The common-ball calculation caps each scheduled waiting interval at `h` attempts and retains every timeout prefix. Positive integer arithmetic supplies enclosing masses; discarded numerical mass is charged. |
| Output | The four-bit success law and all timeout prefixes are a declared coarsening of the richer current/controller history. A complete law here means all outcomes of this declared output, not every source mark or terminal quantum state. |

The source identification is the one already checked by the stopped-clock
comparison:

| Quantity | Value |
| --- | ---: |
| Full screen levels | 3 and 5 |
| Full screen ports | 15,360 and 245,760 |
| Full screen seams `M` | 40,320 and 645,120 |
| Root clock region | 12 ports |
| Incident clock seams `m` | 33 |
| Support of incident marks `S=B_1(root)` | 15 ports |
| Common induced ball `Q=B_9(root)` | 222 ports, 561 internal seams |
| Outgoing seams of the common ball | 45 |
| Common-ball clock fraction `p` | `33/561=1/17` |

Retaining incident marks requires the fifteen-port support `S`. The
positive-transfer dictionary below has the narrower condition that its
field projectors and intervention commute with the remote transfer. The
two root-port fields and the internal repair `(0,1)` satisfy both
requirements. Arbitrary interventions elsewhere in `S` need not do so.

## The same positive transfer for a finite-density preparation

On the complete binary configuration space of any one of these finite
graphs, use column distributions and set

```text
A = (1/M) sum_(remote e) Q_e,
B = (1/M) sum_(incident e) Q_e,
D = I-A,                  p = m/M,
K = B D^(-1),
T = D^(-1/2) B D^(-1/2).
```

Each `Q_e` is a nonnegative symmetric orthogonal projection. Consequently
`0<=A< I`, `p I<=D<=I`, `0<=T<=I`, and `T` is symmetric stochastic and
entrywise nonnegative. The completed native clock has transfer `K`, with
the remote repairs before the terminating incident repair. The exact
similarity is `K=D^(1/2) T D^(-1/2)`.

Let `pi_rho(x)=rho^N(x)(1-rho)^(V-N(x))`, where `N(x)` is the total
occupation, and let `Omega_rho(x)=sqrt(pi_rho(x))`. Every native repair
preserves `N`, including each marked branch. Therefore its matrix commutes
with the diagonal matrix `W=diag(Omega_rho)`. Field projectors commute with
`W` as well. For every complete branch product `C`,

```text
1^T C pi_rho = Omega_rho^T C Omega_rho.
```

This follows by writing `pi_rho=W Omega_rho`, moving `W` through the
number-preserving product, and using `1^T W=Omega_rho^T`. Moreover,
`A Omega_rho=(1-p)Omega_rho`, `B Omega_rho=p Omega_rho` and
`D^s Omega_rho=p^s Omega_rho`.

Remote seams have no endpoint in the root. Root occupation projectors and
the internal repair `(0,1)` thus commute with `D` and its real powers.
Insert the similarity for each `K` in a field cylinder. Adjacent factors
cancel through these projectors and repairs; the remaining boundary
factors cancel against `Omega_rho`. This proves equality of every complete
local field word under the native `K` process and the effective `T`
process, with the same projectors and internal intervention dictionary.
The proof is branchwise and so permits the stated record-dependent choice
to apply or omit the repair. Initial conditioning is already represented
by its retained projectors; it does not require a new stationary state
after a read.

If more incident marks are retained, the existing marked construction uses
`T_m=D^(-1/2) B_m D^(-1/2)`, with `sum_m T_m=T`. Its nonnegative branch
matrices have the same complete cylinder probabilities after the same
cancellation. This identifies a marked effective instrument, not the
native process's unobserved attempt times or hidden trajectories.

The transfer `T` is determined by the graph and clock. The density changes
its boundary preparation, not the transfer itself. The argument also
follows by mixing the existing uniform fixed-number-sector identities
with their Bernoulli weights, without restricting to the two-particle
sector.

Strict positivity holds on the entire binary configuration space. For
each triangular face of the root carrier,

```text
sum_(e in triangle) Q_e >= (3/2) I.
```

Fix the spectator bits. On a triangle with zero or three occupied ports,
all three maps are identity. On its one- or two-occupation sector the
projection sum has eigenvalues `3` and `3/2`. Summing the carrier's twenty
triangles, with every one of its thirty edges counted twice, gives
`sum_(root-carrier e) Q_e >= 15 I`. All these seams are incident, so

```text
B >= (15/M) I,        T >= (15/M) I,
H = -log T.
```

The producer checks this premise against its actual root seam menu. It
enumerates the twenty source triangles and verifies the double cover of
all thirty incident carrier seams. On all eight binary triangle states it
also calls the original native endpoint/coin rule, checks identity-plus-swap,
and verifies the exact symmetric count-matrix identity
`(C-3I)(C-6I)=0` for `C=2 sum_e Q_e`. Thus the strict lower floor is not
inferred from a floating spectrum or assumed from the carrier's name.

The certified lower floors are `1/2688`, `1/43008` and `5/187` for levels
three, five and the common ball respectively. Each finite `H` is therefore
well defined, and the existing logarithm series gives, for a floor `alpha`,

```text
0 <= H - sum_(j=1)^r (I-T)^j/j
  <= (1-alpha)^(r+1)/((r+1)alpha) I.
```

The floors shrink with the full screen's seam count. This does not supply
a uniform global operator-norm estimate, a limiting infinite-screen
Hamiltonian, or a coherent native implementation of `exp(-i s H)`.
The comparison established below is of specified finite field laws.

## Exact exclusion dual for the entire four-bit law

[The backward operator](dual.py) works with squarefree monomials
`n_X=product_(i in X) n_i`. A native seam satisfies exactly

```text
Q_e n_X = (n_X + n_(swap_e X))/2.
```

It preserves the cardinality of `X`. The degree-two action is the
two-label exclusion lift: when both endpoints belong to `X`, swapping
does not create independent walks or a repeated port. This is the exact
action on binary configurations.

The observable dual of the completed distribution kernel is

```text
L = K^T = (I-A)^(-1) B.
```

For a cap of `h` attempts in each scheduled interval, it is

```text
L_h = sum_(j=0)^(h-1) A^j B,
r = 1-p,                  q = 1-r^h.
```

The incident action comes first when propagating backwards. Reversing the
order is not a harmless implementation choice, even though a stationary
commuting field dictionary can conceal that reversal in some cylinder
probabilities. Every column and row sum of `L_h` is `q` on each monomial
degree.

For `k` scheduled intervals define the following polynomials, evaluated
on the initial configuration:

```text
v    = L_h^(k-1) n_1,
F    = L_h n_0,
L0   = L_h v,
J0   = L_h (n_0 v),
J1   = L_h (n_0 Q_(0,1) v)     with feedback,
J1   = J0                      without feedback.
```

Here multiplication by `n_0` occurs after propagating through the
intervention backwards. Since `n_0^2=n_0`, the products `n_0 v` and
`n_0 Q_(0,1)v` have degree at most two. Their coefficient masses before
the first interval are `q^(k-1)`; those of `L0,J0,J1` after it are `q^k`.
The coefficient mass of `F` is `q`.

Condition on the initial retained bits `(a,b)`. The expectation of a
squarefree monomial is the product of `a` or `b` for its disclosed ports
and one factor `rho` for each other distinct port. Write this conditional
expectation as `E_ab`. The four successful conditional masses are exactly

```text
Pr(c=1,d=1, success | a,b) = E_ab[J1],
Pr(c=1,d=0, success | a,b) = q^(k-1) E_ab[F] - E_ab[J1],
Pr(c=0,d=1, success | a,b) = E_ab[L0] - E_ab[J0],
Pr(c=0,d=0, success | a,b) = q^k - q^(k-1) E_ab[F]
                                      - E_ab[L0] + E_ab[J0].
```

Multiplying by `pi_ab=rho^(a+b)(1-rho)^(2-a-b)` gives all sixteen joint
success probabilities. Feedback appears only on the `c=1` branches;
replacing both `J0` and `J1` by a single averaged controlled expression
would implement a different observer.

Thus the entire four-bit joint law, including its fourth connected
cumulant, is obtained from exact degree-one and degree-two source actions.
Conditioning accounts for the two initial bits before taking expectations.
This is not a first/second-moment closure of a general multi-time process.
The common-ball calculation uses 222 singleton and 24,531 pair monomials
instead of enumerating its `2^222` configurations.

## Timeouts are output branches

Let `j` count successful scheduled intervals before a timeout, with
`0<=j<k`. Before the first hit the only retained labels are `(a,b)`:

```text
Pr(a,b,timeout after 0) = pi_ab (1-q).
```

For `j>=1`, the already observed `c` remains in the timeout record:

```text
Pr(a,b,c=1,timeout after j)
    = pi_ab q^(j-1)(1-q) E_ab[F],
Pr(a,b,c=0,timeout after j)
    = pi_ab q^(j-1)(1-q) (q-E_ab[F]).
```

Subsequent success and timeout probabilities depend only on the scheduler,
not the state or first bit. The native intervention is trace preserving,
so these formulas apply with or without feedback. Individual prefix
probabilities retain their first-read dependence; replacing them merely
by their aggregate `q^j(1-q)` would lose output information.

There are `16+4+8(k-1)` records, hence 44 at `k=4`. Their total mass is one,
the success mass is `q^k`, and the total timeout mass is exactly
`tau=1-q^k`. Later timeout branches retain the additional native feedback
tick when `c=1`. No successful record distribution is renormalized.

## Directed arithmetic and complete error accounting

[The field-law implementation](refinement.py) computes coefficient lower
bounds on a grid `G=2^bits`. Each sparse subset action has integer entries
and row/column sum `2*number_of_subset_seams` before division by `2M`.
Division is always downward. Positive maps preserve componentwise
domination, including through successive capped intervals. The feedback
average of the two degree-one coefficients is also rounded down.

All supplied coefficient masses are at most one. For the retained ball
and `bits=52`, `2M G=1122*2^52<2^63` bounds every integer multiplication
and positive sparse accumulation. Partial completed sums remain at most
`G`. The implementation checks the corresponding bound for its declared
graph before using signed 64-bit arithmetic. No floating normalization,
eigensolver tolerance or unmeasured matrix residual enters this proof.

If a nonnegative exact coefficient vector has known mass `s` and its
computed lower vector is `u/G`, then

```text
delta = s - sum_i u_i/G,
E_ab[u/G] <= E_ab[exact] <= E_ab[u/G] + delta,
```

because each conditioned monomial lies in `[0,1]`. Here `s` is `q` or
`q^k`, including any error inherited by a previous lower vector. Exact
rational interval addition and subtraction then enclose every displayed
success and timeout formula. Negative terms use the opposite endpoint;
known probability bounds may be intersected with the interval.

Let `[ell_z,u_z]` be the intervals for every record of the complete capped
law `P_h`. The construction proves `ell_z<=P_h(z)<=u_z` and

```text
delta_cap = 1 - sum_(all z) ell_z
          = ||P_h - ell||_1.
```

This is a one-sided subprobability error, not a claim that half this L1
quantity bounds every event against an unequal-mass vector. Its exact
timeout total is stored separately from its arithmetic interval widths.

For the completed sixteen-record law `P`, every successful capped lower
mass remains a lower bound. Set

```text
mu = 1 - sum_(successful z) ell_z.
```

Then `mu` includes both unresolved waiting and arithmetic losses. For any
subset of the four retained bits, summing `ell_z` over records where those
bits are all one gives a lower indicator moment `l`; its completed moment
lies in `[l,min(1,l+mu)]`. This charges missing mass once and does not
condition on success. Increasing the cap and grid precision supplies the
corresponding convergence procedure.

Receipt endpoints may be exported by exact outward rounding to a fixed
dyadic grid. The complete deficit must then be recomputed from the
exported lower masses, or explicitly enlarged by the sum of their downward
rounding losses. Likewise the completed missing mass uses the sixteen
exported success lower masses when it certifies that exported vector.
Rounding one previously computed deficit cannot in general cover all
separate row-rounding errors.

## Composing with the actual-screen locality theorem

The supplied Bernoulli preparations have identical restrictions to `Q`.
Both initial field disclosures, every counter update and the internal
record-controlled repair are permitted local operations for the existing
complete-output theorem. The extra repair is counted in the observer's
clock and stays in the local support; it is not an additional random
scheduled hit in the stopped-word estimate.

Let `epsilon_Q(k)` be the existing one-screen first-exit bound. It is
computed directly from actual integer walk counts in the common ball,
using `m=33` and summation through degree 48 with its rational remainder.
For the completed field law on either full screen `G`, postprocessing
contractivity gives

```text
TV(P_G, P_Q) <= epsilon_Q(k),
TV(P_level3, P_level5) <= 2 epsilon_Q(k).
```

These are normalized classical laws. They are a consequence of the
stronger complete quantum-output comparison, not a new diagonal-only
locality assumption. In particular the spatial allowance is added after
the field event is extracted; it is not multiplied by a density or guessed
amplitude factor.

The retained outward pair allowances from that theorem are:

| Scheduled incident hits | Complete pair allowance |
| ---: | ---: |
| 2 | `2.731e-10` |
| 4 | `6.656e-9` |
| 8 | `4.296e-7` |

Thus a common-ball completed indicator interval `[l,u]` becomes
`[max(0,l-epsilon_Q), min(1,u+epsilon_Q)]` for either full screen.
The positive-transfer dictionary identifies these same completed field
laws with their own source-derived finite effective transfers and boundary
vectors. It does not assert that the three global transfer matrices are
equal, or close in operator norm.

The fourth joint cumulant is evaluated from its fifteen set partitions:

```text
kappa(X0,X1,X2,X3)
  = sum_(partitions pi) (-1)^(|pi|-1)(|pi|-1)!
                         product_(C in pi) E[product_(i in C) Xi].
```

All moment intervals lie in `[0,1]`; multiplication is monotone there and
negative partition coefficients reverse the endpoints. The calculation
uses all sixteen success labels and the completed missing-mass budget.
Spatial error may be added to the individual moments before this interval
calculation. Alternatively, normalized laws at TV distance `epsilon`
have fourth-cumulant difference at most `75 epsilon`: each block moment
changes by at most `epsilon`, each product changes by at most
`|pi| epsilon`, and the coefficient sum is `1+14+36+24=75`.

## Executions, resources and remaining obligations

The full-source execution reuses the actual-screen
[current-history runtime](../native_action/clock_histories.py), with both
initial field disclosures made before evolution. It retains the same
first-bit-controlled native repair and reconstructs the second field from
its own incident currents. The probability certificate is deterministic
and does not interpret a small set of seeded histories as a statistical
confidence interval.

The retained executions use seeds zero through seven with feedback both
disabled and enabled, for sixteen paired full-screen executions. Seeds five
and six exercise a feedback repair with signed current `+1` on `(0,1)`,
moving occupation from root port zero to root port one. Every execution
also resumes a serialized intermediate checkpoint and compares the full
states, records and random-generator continuation. These examples check
actual field access and continued work; the analytic law and error bounds
do not depend on their observed frequencies. A separate one-slot run
retains an externally forced timeout.

The serialized observer checkpoint wraps the existing full-source
checkpoint with both initial disclosures, both current counters, and the
first-read memory for each screen. Loading checks these original values
against the source prefix, including exact types and record ancestry.
This source replay is a validation operation, not another observation
available to the controller. After restoration, the observer advances the
saved counters and first-read memory using only the newly produced record
suffix. It neither resets the counters from terminal occupations nor
reprocesses the old prefix as part of this operational continuation.
Independent full-history reconstruction and comparisons with terminal
occupations are audit checks of the result, not hidden-state inputs to
the feedback policy.

The expected number of scheduled attempts for the completed target is
`k M/33` on each full screen and `17k` on the common ball. With feedback,
the expected additional native work is `rho`: the first field has marginal
`rho` before any controlled repair. Thus the expected complete native
work is `k M/33+rho` or `17k+rho`. These are operation counts, not physical
durations. A cap of `h` per common-ball interval permits at most `kh+1`
native attempts including feedback; the computational cost of propagating
24,753 dual coefficients is a separate resource.

The common-ball cap approximates completed-clock laws. It is not a common
cap on the two full screens' global attempts: identity-slot erasure changes
those counts. An externally capped full-source simulation retains its
own timeout and must charge its separately certified cap allowance when
compared with a completed target.

For four scheduled hits, cap 512 and a `2^52` coefficient grid, the
reproducible producer gives the following outward summaries. All quantities
are unconditioned on success.

| Quantity | Baseline | First-bit-controlled repair |
| --- | ---: | ---: |
| Complete capped record count | 44 | 44 |
| Complete capped lower-mass L1 deficit | `<1.209e-9` | `<1.212e-9` |
| Total interval-cap timeout | `<1.324e-13` | `<1.324e-13` |
| Fourth connected cumulant on either full screen | `[0.0008638252,0.0008638837]` | `[-0.0049753855,-0.0049753193]` |

The cumulant intervals include the one-screen complete-output allowance
`epsilon_Q(4)<3.328e-9`, applied at the indicator moments. In particular,
their opposite signs survive the complete numerical and spatial budgets.
The final field mean increases under the native feedback by an amount in

```text
0.0836900722 <= E[n_1(final)]_feedback - E[n_1(final)]_baseline
            <= 0.0836900880
```

on either full screen. This is a difference of the two independently
enclosed protocols; it does not assume feedback preserves the stationary
preparation. The [retained receipt](refinement_results.json) contains all
44 row intervals per protocol, lower-mass deficits recomputed after dyadic
export, completed missing-mass and timeout allowances, joint-four and
connected-four intervals, and the actual full-source execution ledgers.

The independent [small-source oracle](test_refinement.py) propagates all
binary configurations through the original native endpoint/coin rule.
It checks every success and timeout record against the dual intervals at
multiple caps and densities, including coarse grids. Its complete terminal
states distinguish reversed waiting/incident order even where a stationary
local dictionary hides that mistake. Separate controls distinguish removed
remote motion and omitted record-controlled repair. The fourth-cumulant
oracle uses an explicit partition formula independent of the producer's
partition generator.

Reproduction uses the same producer for generation and checked replay:

```text
python -m model.observer_time_fields.refinement --check model/observer_time_fields/refinement_results.json
python -m pytest model/observer_time_fields/test_dual.py model/observer_time_fields/test_refinement.py
```

The proof supplies a finite source-relative clock and effective-field
comparison for two actual screens with a compatible, non-diluting
preparation. It leaves an all-level refinement limit, physical preparation
and detector error, a native coherent realization of the effective
Hamiltonian, and dictionaries not commuting with the remote transfer as
explicit further obligations. A movement-only clock retains its negative
eigenvalue obstruction; changing the clock is not a free application of
the positive-transfer theorem.
