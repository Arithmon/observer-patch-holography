# Comparing actual screens in an incident observer clock

The complete local-output comparison in [native action](README.md) can be
stopped at a prescribed number of incident repair attempts. The resulting
bound uses the original uniform seam scheduler directly. No per-edge Poisson
parameter, duration per operation, or identification with modular time is
required. This is a finite comparison of declared source instances, not an
all-level refinement limit.

This is an attachment of the [canonical primitive observer
contract](../primitive_observer/README.md). Its four-port quantum control
uses the existing [executable runtime](../primitive_observer/runtime.py).
The full-screen histories execute the diagonal restriction of the same
native source law; they do not simulate a 245,760-qubit density matrix.

| Interface | Instantiation |
| --- | --- |
| State and preparation | The theorem admits matching local quantum restrictions and references. Full-screen examples supply Bernoulli occupations at density `1/3`, coupled on the matched ball, and disclose the initial field at root port zero. |
| Operation | Uniform full-source seam selection and the existing fair native repair. C and E coincide on the executed diagonal histories. |
| Record and field | The observer receives authentic root-incident seam, coin and signed-current marks. Its disclosed field counter integrates those currents; each later field value requires no additional occupation measurement. |
| Feedback and clock | First-field bit one requests a native repair on root seam `(0,1)`. That repair contributes its own record and incident tick. Requested scheduled-hit counts and total incident ticks remain distinct. |
| Continuation | Serialized full source states, records, counters and random-generator state continue the same seeded execution. Neither remote evolution nor the conditional state is reset. |
| Comparison and resources | Actual matched source balls, complete-output spatial bounds, directed arithmetic bounds and external cutoff tails are separate obligations. Global operation counts are an external resource ledger, not controller inputs. |

## Contract and matched source instance

On a finite graph with `M` seams, choose each successive seam independently
and uniformly and apply the same native C or E instrument. Fix a root region
`A` and let `I` be its `m>0` incident seams. Every attempted seam in `I`
increments the clock, including a zero-current attempt. Stop just after its
`k`th tick. This stopping rule depends on the seam schedule, not on whether a
repair changes occupation or on a selected quantum outcome.

To retain the incident seam, coin and current marks, the observation support
must contain both endpoints of every seam in `I`. Take `S=B_1(A)`. The
observer may retain those marks, local field readouts, its local controller
memory, and its terminal quantum state, jointly with an arbitrary untouched
reference. A feedback instrument may depend causally on those local records
and act within `S`. The controller cannot consult exterior records or the
global attempt counter. Its intervention rule uses the retained local event
order; operations scheduled by a global slot index, elapsed physical time,
or an auxiliary external clock are outside this comparison. Its operations
are finite or locally finite, with no
infinite control sequence before the stopping point. Any transported mark
labels, including signed-current orientations and retained endpoint coin
labels, must use the actual source instrument correspondence.

The [screen identification](screen_compatibility.py) at levels three and five
supplies the following finite instance:

| Quantity | Value |
| --- | ---: |
| Root region `A` | One original twelve-port carrier |
| Incident clock seams `m` | 33 |
| Observation support `S=B_1(A)` | 15 ports |
| Common induced ball `Q=B_9(A)` | 222 ports and 561 internal seams |
| Maximum full-graph degree `Delta` | 6 |
| Minimum path length from `S` to outside `Q` | 9 |

All clock seams lie inside `Q`. Its two source-address maps identify the
induced native instruments, not merely an abstract graph spectrum. The
global graphs have different seam counts. Their unobserved global waiting
counts need not agree.

## An exact stopped-word identity

Let `N_w` count ordered subsequence occurrences of a specified seam word
`w=(e_1,...,e_l)` in the schedule through the `k`th incident tick, including
that final attempt. Repeated letters are allowed. Let

```text
s(w) = number of incident letters among e_1,...,e_(l-1).
```

Then the uniform scheduler gives the exact identity

```text
E N_w = binom(k+l-1-s(w), l) / m^l
      <= binom(k+l-1, l) / m^l.                       (1)
```

A binomial coefficient with upper argument below `l` is zero. The global
seam count `M` cancels. The order of the word matters through the omission
of its final letter from `s(w)`.

**Proof.** Split the schedule into `k` independent blocks, each consisting
of remote attempts followed by its first incident attempt. Introduce a
noncommuting formal letter `X_e` for each seam. In the product of factors
`(1+X_e)` along a schedule, the coefficient of `X_(e_1)...X_(e_l)` counts
exactly the ordered occurrences of `w`. Summing the geometric remote word
in one block gives its expected subsequence series

```text
(1-R/m)^(-1) (1+C/m),
R = sum_(e not in I) X_e,       C = sum_(e in I) X_e.
```

Indeed the original remote-word weights are `M^(-n-1)`; summing the
unselected remote letters supplies the geometric normalization `m/M`.
The displayed identity is coefficientwise, so it does not require numerical
evaluation of an infinite formal series.

In the product of `k` such block series, each selected letter belongs to a
weakly increasing block index in `{1,...,k}`. After an incident selected
letter, any subsequent selected letter must belong to a strictly later
block. These are exactly the `s(w)` mandatory strict gaps. Subtracting the
number of previous strict gaps converts the assignments to weakly increasing
indices in `{1,...,k-s(w)}`. Their count is the binomial coefficient in (1),
and each assignment has coefficient `m^(-l)`. This also proves the formula
for repeated seams. Taking expectations of coefficients is justified by
nonnegativity; the geometric waiting blocks have finite moments of every
fixed order.

For example, with one incident letter `i` and one remote letter `r`,
`E N_(i,r)=binom(k,2)/m^2`, whereas
`E N_(r,i)=binom(k+1,2)/m^2`. At the first tick the former is zero and the
latter is positive. Treating the stopping time as independent of its own
incident seam events would lose this distinction.

## Complete quantum output at the stopped clock

Compare the full graph with its induced simulation on `Q`, using the same
initial reduced state on `Q` and the untouched reference. Outside operations
are omitted; equivalently, keep their original attempt slots as identity
holds. Retain all output branches. Since every clock seam is internal, this
does not alter the clock events in a coupled schedule.

For a fixed schedule, propagate the support of the complete local instrument
backwards from `S`, as in the existing complete observation theorem. Every
retained record and permitted feedback operation stays in `S`. If the
dependency cone does not cross `Q`, the full and restricted complete quantum
outputs coincide, including their reference and controller memories. A
crossing implies a temporally ordered path from `S` to outside `Q`. Reverse
that spatial path to obtain its chronological seam word, and apply (1).
On a bad schedule the trace distance of normalized complete outputs is at
most one. Averaging therefore bounds the unconditional complete output,
without a Hilbert-space dimension factor or a restriction to diagonal states.
This applies separately to C and E; it does not identify their predictions.

If every exit path has at least `d` edges, counting at most
`|S| Delta^d` possible length-`d` prefixes gives

```text
epsilon_graph(k) <= min(1,
    |S| (Delta/m)^d binom(k+d-1,d)).                 (2)
```

The stopped process has mass one: its `k`th incident hit occurs almost surely
because `m/M>0`. No conditioning on successful movement or rare field
outcomes is used. Conditional, postselected comparisons would need additional
probability denominators and are not supplied by this theorem.

## First-exit walks and a finite exact computation

The actual source neighborhood gives a stronger bound than degree counting.
Let `J` be the integer adjacency matrix of the induced graph on `Q` and let
`b_v` count full-graph seams from `v in Q` to outside `Q`. Define

```text
c_l = 1_S^T J^(l-1) b,       l >= 1.
```

This counts walks whose first `l-1` edges remain in `Q` and whose final edge
exits. Walks may repeat vertices or edges; overcounting paths is harmless.
The counts are produced by integer adjacency multiplication, including the
actual source boundary degrees. Equation (1) and a union bound give

```text
epsilon_graph(k) <= min(1,
    sum_(l>=1) c_l binom(k+l-1,l) / m^l).           (3)
```

In particular a source-free replacement of the remote transfer by a scalar
waiting factor cannot justify these coefficients. Actual remote paths are
what the count retains.

For a finite computation through degree `L`, use the exact rational sum in
(3) and bound the remainder by `c_l <= |S| Delta^l`. Put

```text
a_l = |S| (Delta/m)^l binom(k+l-1,l),
q   = (Delta/m) (k+L+1)/(L+2).
```

Whenever `q<1`, the tail after degree `L` is at most

```text
sum_(l>L) c_l binom(k+l-1,l)/m^l <= a_(L+1)/(1-q).  (4)
```

For `k>=1`, the successive ratio
`a_(l+1)/a_l=(Delta/m)(k+l)/(l+1)` is nonincreasing in `l`.
Thus every ratio after the first omitted term is at most `q`, proving (4).
All inputs and arithmetic in this certificate are integers or rationals.
The present source has `Delta/m=2/11<1`, so increasing `L` makes the
remainder vanish for every fixed finite `k`.

## Comparing the two screens

Delete operations outside the common `Q` in each graph. After forgetting
the unobserved identity slots, each restricted stopped schedule has the same
law: it is uniform on the common internal seams and stops at the same
incident seam count. The original global seam count changes waiting costs,
but not that compressed law. Matching native instruments and local controls
therefore give one common restricted experiment.

If the two initial restrictions on `Q` and the reference have trace distance
`delta`, contractivity and the two truncation bounds give

```text
D_trace(output_left, output_right)
    <= min(1, delta + epsilon_left(k) + epsilon_right(k)).  (5)
```

The factor two when the graph bounds agree is a conservative triangle
inequality. No independence of the two bad-cone events is assumed. A stronger
coupling of their boundary event schedules would require its own argument.

Even the coarse bound (2), before using actual first-exit counts, supplies
the following pair allowances at matching initial restrictions:

| Incident ticks `k` | Exact pair allowance | Decimal upper summary |
| ---: | --- | ---: |
| 1 | `15360/2357947691` | `0.000006515` |
| 2 | `153600/2357947691` | `0.000065142` |
| 4 | `307200/214358881` | `0.001433111` |

The [integer first-exit computation](incident_refinement.py) sharpens these
figures using the actual graph. Both screens have 45 outgoing boundary
seams. Their first nonzero counts, at lengths nine through fourteen, are
`(384,4992,54384,477840,3842208,28554816)`. Summation through length 48,
including the rational remainder (4), gives the following outward upper
summaries of the complete pair error:

| Incident ticks | Pair error upper bound |
| ---: | ---: |
| 1 | `2.598e-11` |
| 2 | `2.731e-10` |
| 4 | `6.656e-9` |
| 8 | `4.296e-7` |

The [retained receipt](incident_results.json) stores the full rational
bounds. At four ticks the one-screen walk-series remainder after length 48
is below `2.167e-31`. These allowances concern complete quantum histories,
not only the field statistic computed next.

## Executed stopped-clock field response

For the uncontrolled baseline, the native mean-density action of one seam
is exactly `(I+swap_e)/2`, for both candidates. Write its contributions on
the 222-port ball as `A_mean` for remote seams and `B_mean` for incident
seams, each divided by the ball's 561-seam schedule count. The completed
one-tick density transfer is

```text
K_mean = B_mean (I-A_mean)^(-1).
```

This identity retains remote evolution before every tick. It is a
first-moment reduction of the source, not a model of its complete quantum
records. Both contributions have constant row and column sums, respectively
`1-p` and `p`, so `K_mean` is doubly stochastic. Equivalently, the supplied
Bernoulli product state is preserved by each unread native seam and hence
by every schedule mixture stopped solely by seam arrivals.

On the product input at density `rho`, an initial occupation read and the
same-port read after `k` ticks satisfy

```text
Pr(n_0(0)=1, n_0(k)=1) = rho^2 + rho(1-rho) [K_mean^k]_(0,0),
Cov(n_0(0), n_0(k))    = rho(1-rho) [K_mean^k]_(0,0).       (6)
```

Indeed `E[n_i(0)n_j(0)]=rho^2+rho(1-rho)delta_(i,j)` and the final
conditional mean is linear with transfer `K_mean^k`. This calculation uses
the exact source mean action and initial product moments. It introduces no
independent-particle approximation to the complete field law. Formula (6)
does not assert stationarity under record-selected feedback; feedback is
covered separately by the complete-output comparison.

The executable computes a lower bound on the entire column
`K_mean^k e_0`. It uses a grid of `2^52`, floors every positive integer
matrix action downwards, and retains at most 512 native attempts per
scheduled interval. Each subset action has integer row and column sum
`2*number_of_subset_seams` before division by `1122`. All active and
completed integer masses together remain at most `2^52`. Consequently
`1122*2^52 < 2^63` bounds every multiply and sum, including intermediate
positive accumulations. There is no floating normalization.

Let `x` be the resulting unnormalized lower column and
`d=1-sum_i x_i`. The exact capped column has mass `1-tau`, where

```text
tau = 1-(1-(16/17)^512)^k.
```

Positivity proves componentwise domination of the capped and completed
columns. Thus `d-tau` is the finite-cap arithmetic L1 discrepancy and `d`
is the completed-column L1 discrepancy. In particular
`[K_mean^k]_(0,0)` lies in `[x_0,min(1,x_0+d)]`. This computational cap
discards unresolved column mass only inside an explicitly charged
approximation; it is not postselection on a successful field experiment.

Apply (6) to this interval first. Then add the one-screen complete-history
allowance from (3) to each side of the joint-11 probability interval. The
spatial error is not multiplied by `rho(1-rho)`: a generic complete-output
error bounds the final event probability itself. Exact baseline stationarity
sets both marginals to `rho`, so subtracting the common `rho^2` gives the
same spatial allowance for the connected field.

At `rho=1/3`, the four-tick column has total L1 allowance below `1.851e-11`,
including timeout below `1.324e-13`. Combining its directed arithmetic
interval with spatial truncation gives, for either actual full screen,

```text
0.15332581621 <= Cov(n_0(0),n_0(4)) <= 0.15332582288.
```

The decimal endpoints are outward comparisons against the exact rational
receipt. The completed target uses an expected 68 retained-ball attempts
for four ticks; full-source expectations are `4M/33` with each screen's
actual `M`. These operation counts are not physical durations. The grid
propagation and its 512-attempt cap are reported separately from those
completed-target resource expectations.

## Full-source histories and quantum controls

The [full-screen executor](clock_histories.py) couples the two *entire* seam
menus, not just the matched region. Its union contains 684,879 seam slots;
561 internal seams are shared, and each other seam belongs to one screen.
A uniform union slot applies the original repair on every screen containing
it and identity on the other. Erasing a screen's identity slots leaves
independent uniform choices among all its own source seams. The common
33 incident seams therefore give one shared scheduled clock. Endpoint-order
reversals transport both the current sign and the placement coin.

Seeds zero through three are executed without selection at two and four
scheduled hits. Every example retains its full observed history and terminal
local state. All eight restore an intermediate serialized checkpoint to
exactly the same final states, records and random-generator position. Seed
one executes the feedback repair; remote occupations change in the same
histories. These are reproducible executions, not frequency-based confidence
estimates or claims that all possible paired records are identical.

Ideal product preparation, one initial local bit disclosure and faithful
current/coin records are declared services. Terminal output reads fifteen
local bits. The stored source occupies 261,120 binary-state bytes; the receipt
reports actual attempts, records and serialized checkpoint sizes. Physical
preparation, storage and detector error are not measured by this simulation.

The external simulator cap is 500,000 union slots, with at most one extra
native feedback repair per screen. Its timeout is an explicit outcome;
a separate one-slot run forces that outcome. For `p=33/684879`, the probability
of fewer than `k` scheduled ticks before cap `N` is binomial. The exact rational
upper bound uses

```text
Pr(timeout) = (1-p)^N sum_(j<k) binom(N,j) [p/(1-p)]^j,
(1-p)^N <= exp(-Np) <= 1 / sum_(j=0..64) (Np)^j/j!.
```

It is below `8.642e-10` for two ticks and `9.114e-8` for four ticks.
This cap allowance is separate from the completed-clock spatial bound;
finite capped output cannot inherit the smaller completed-law bound without
charging its timeout. The density-column computation above has its own
512-attempt-per-interval cap and separately certified arithmetic allowance.

The [quantum control](stopped_control.py) executes C and E on a declared
four-port, three-seam path. It retains all thirteen observer outcomes when
stopping at two incident hits or an external cap of three attempts, including
exact timeout mass `20/27`. All 256 input matrix units agree with an independent
raw endpoint/coin oracle, so the exact instrument extends to arbitrary
untouched references. A Bell-reference probe runs in the canonical runtime
with a rational entrywise trace-distance allowance `10^-12`, and its actual
timeout state continues after serialization. Its current-controlled phase
changes the local root output for C; E's incident read removes that local
coherence. This small controlled path is not represented as a full-screen
quantum execution. The phase is a supplied control service, not a native
repair or an additional native incident tick.

The composition is thus explicit: source endpoint maps give the instrument;
the existing dependency-cone identity applies to its complete local output;
the stopped-word identity bounds that cone at an incident count; the checked
source map and matching product restrictions give the common finite experiment.
The density identity then extracts an observable from that experiment.
Movement-only clocks, global-counter-triggered controls, mismatched initial
restrictions and an all-level graph limit require their own hypotheses.

## Preparation, feedback and remaining obligations

The declared Bernoulli product preparation at density `rho=1/3` has exactly
matching local restrictions under the source maps and nonzero local
fluctuations. The existing native-instrument calculation proves its unread
stationarity, although record-selected feedback need not preserve it.
More generally, any compatible local preparation and untouched reference
is admitted by (5). The density is a preparation input, not a value selected
by the clock or repair law.

Uniform global total-two preparations are a different domain. Their
222-port restrictions on these two screens differ by the already certified
trace distance `831558460593/30921625444352`, approximately `0.02689`.
That `delta` must be added in (5); it cannot be discarded by identifying
their local clocks. These preparations dilute towards the vacuum as the
screen grows, unlike the stated fixed-density product domain.

An explicitly commanded native repair on an incident seam is itself an
incident tick. If the protocol stops at the `k`th total tick, including such
feedback, its uniform-scheduler prefix ends no later than the `k`th scheduled
incident hit. The bound above remains conservative, provided the intervening
controls stay in `S` and the accounting retains their extra ticks. A protocol
defined instead by `k` scheduled hits must report its additional feedback
ticks separately. Neither convention permits moving feedback outside the
observed support or consulting an exterior record.

A finite global-attempt cap is additional to this completed-clock theorem.
Its timeout outcomes must be retained and their probabilities charged. The
negative-binomial law for the global count depends on `m/M`, so equal caps
on different screens are not silently identified. Numerical propagation
errors likewise require their own complete-output or observable bounds.
The exact geometric certificate is not a floating execution certificate.

This comparison closes the auxiliary-clock requirement for the stated pair
of actual source neighborhoods and incident-clock observations. It does not
prove compatible neighborhoods at every level, a spatial continuum, a
physical clock calibration, a preferred C/E continuation, or a new empirical
prediction without the declared preparation and observation services.
