# Observations, publication and paid waiting

## The acceptance theorem

Let `possible(v,o)` mean that the source model, intervention `v`, permitted
history and declared uncertainty can produce receiver observation `o`.
Construct this relation without a publication predicate. A partial policy
`publish(o)` is sound if every published value equals every payload compatible
with `o`. For a feasible observation, publication of `v` is possible exactly
when the compatible set is `{v}`. The canonical policy returns that value and
abstains on a nonsingleton or empty set. Any sound policy agrees with it
wherever the sound policy publishes on a feasible observation.

`identifies_unique`, `sound_publication_identifies`, `canonical_some_iff`,
`canonical_sound` and `canonical_maximal` prove these statements. Feasibility
matters: a policy can satisfy the soundness implication vacuously on an
impossible observation. The canonical rule rejects that observation.

`identification_refines` proves persistence of the identified value when
additional observations shrink the compatible set without removing the actual
payload. `cut_transcript_abstains` composes the sound-policy theorem with the
native transcript locality theorem. If neither source intervention crosses
the separating cut, even the full local prefix transcript cannot justify
publication. Fixed private random coins do not remove the ambiguity.

These results derive a greatest sound partial function for a specified
observational model. They do not force publication, define a physical
instrument or prove that all observers choose the same accepted domain.

## A computable interval realization

Prepare five balanced cells `(2+a_i,2-a_i)`, with source amplitude
`a_0 = v A`, `v` in `{-1,+1}`, `A > 0`, and four blank cells. The two
interventions have identical total ideal scalar load. The source bit is not
encoded in a schedule, address, delay or receiver preparation.

A paired copy on adjacent cells executes two scalar means. For arbitrary
input amplitudes it applies the scalar mean to the amplitudes. Thus an
arbitrary word produces receiver amplitude `v A g`, where `0 <= g <= 1`
depends on the word and is independent of `v`. No gain, assigned reset or
payload refresh occurs. Raw signal attenuates; logical publication is separate
from the scalar update.

The receiver samples its own two rails and forms their contrast `z`. Suppose
the total contrast error is bounded by `E >= 0`. The conservative local
uncertainty model allows every gain in `[0,1]`, giving the possible intervals

```
v = +1: [-E, A+E]
v = -1: [-A-E, E].
```

The canonical policy for these intervals publishes +1 for `E < z <= A+E`,
publishes -1 for `-A-E <= z < -E`, and otherwise abstains. In particular it
rejects both the shared interval `[-E,E]` and impossible out-of-range values.
`bounded_is_canonical` proves exact equality with the abstract singleton rule.
The interval model overapproximates the discrete gains of finite words; no
claim of maximality over a richer instrument or finer noise model follows.

The policy's only runtime arguments are its two local samples and checkpoint
number. Amplitude and error allowances are public, payload-independent
parameters. Observation of schedule letters or future history is unnecessary.
Before influence arrives, the ideal receiver transcript is identical under
both interventions. No sound local policy can publish on that shared
transcript. Under the sufficient margin below, the interval policy publishes
at the first influence checkpoint, for every allowed error realization.

## Arbitrary native words and precision

`run_affine`, `mean_retains_positive`, `run_retains_positive` and `run_gap`
prove affine transport, persistence of positive influence, and the dyadic gap:
starting from amplitudes either zero or at least `A`, after `t` paired copies
every nonzero amplitude is at least `A/2^t`. `paired_native` identifies each
amplitude operation with two physical means and `paired_cost` counts them.

With coordinate preparation error `E0`, arbitrary signed disturbance at most
`delta` after each scalar mean (including idle coordinates), and error at most
`rho` in each receiver sample, `paired_readout_error` gives

```
E_t = E0 + 2 t delta + rho.
```

The threshold theorems prove safety for every native word. Whenever influence
is nonzero, `2 E_t < A/2^t` guarantees publication. If influence is zero,
`|z| <= E_t` forces abstention. The theorem includes the exact ambiguous
boundary. Preparation and disturbance bounds are hypotheses, not measured
hardware properties.

`paired_publication_contract` composes the signed native execution, arbitrary
implementation errors, both local sample bounds and the actual guarded
publication rule for either payload in one theorem. The output is correct
whenever published and exists when influence arrives with the stated margin.
The guard cannot discard an observation satisfying those bounds.

The captured execution uses `Q=2^14`, `A=(Q-3)/Q`, `E0=rho=1/(4Q)`, and
`delta=1/(2Q)+1/(8Q)`. Actual operations round each scalar mean to the nearest
grid value, with even ties. Odd sums occur, so rounding is exercised.
The first term in `delta` covers rounding; the second is a separately allowed
disturbance. Histories execute rounding; the Lean estimate covers additional
adversarial disturbances. Local readout endpoints are also replayed. The
strict separation inequality holds at every tested checkpoint `0 <= t <= 8`.

The positive/negative rails are

```
(0,9), (1,5), (14,40), (23,38), (45,39).
```

Every scalar seam and the injective placement are independently checked
against `code/source_routing/support_w12_l3.json`. The source is on carrier 0
and the receiver on carrier 3, which is a neighbouring carrier. General path
theorems assume the required paired embedding; this capture is not a
long-distance demonstration. Isolation and paired completion checkpoints
are supplied control restrictions on the captured scalar alphabet.

## Scheduling, availability and unconditional work

On a path of `d >= 1` paired edges, source influence occupies an initial
segment. `prefix_step`, `path_frontier` and `source_path_frontier` prove that
only the edge at the segment's frontier advances it. `advancing_letter_count`
and `uniform_frontier_sum` count exactly one advancing event identity among
the `d` letters until completion; `completed_frontier_absorbs` handles the
terminal state.

Assume the full length-H word simplex and uniform counting reference on that
alphabet, with no further temporal constraints. The existing
`RepairWordSchedule` KL theorem selects the uniform product law under exactly
these hypotheses. This is a complete-history entropy cover and reference
specified for the finite experiment. Uniform one-step marginals alone do not
imply the product law: choosing one edge uniformly and repeating it forever
has uniform marginals and never delivers for `d > 1`.

Until arrival, the chance of advancing at the next checkpoint is `1/d`,
independently of the previous progress. If `T` is first arrival, the analytic
counting argument gives

```
Pr(T > H) = d^-H sum_{j=0}^{min(d-1,H)} binom(H,j) (d-1)^(H-j)
E[min(T,H)] = sum_{t=0}^{H-1} Pr(T > t)
E[scalar means through publication or abort] = 2 E[min(T,H)].
```

For `d=1`, zero powers use the counting convention `0^0=1`: delivery occurs
after one copy. The horizon-zero history aborts with probability one and
performs zero means. In ideal exact arithmetic the d geometric waiting
stages each have mean d, so `E[T]=d^2` and the mean scalar cost is `2d^2`.
The finite tails tend to zero. `SourceReadAcceptanceSchedule` constructs all
finite histories, proves their count is `d^H`, and identifies the influence and
native cost of their stopped physical prefixes. `full_history_reward` proves
the exact transition recurrence for every reward, including indicators of
arrival or abort. `full_history_cost` proves the unconditional paid-work
recurrence over every complete history. The binomial closed form and infinite
ideal expectation argument are analytic. Independent binomial and integer
recurrence computations cross-check the finite formulas.

The fixed experimental grid is certified only through eight checkpoints.
For a separately declared finite horizon H and unit starting amplitude, take
`Q_H = 2^(H + bit_length(10H+4) + 2)` with the same Q-scaled error allowances.
Then `2 E_H < 2^-H`, so every arrival through H is distinguishable. This is
linear-in-H precision with logarithmic overhead, conditional on correspondingly
smaller physical disturbance. It proves no fixed-noise unbounded service.

For the captured four-edge example, the eight-checkpoint law publishes with
probability `7459/65536` and aborts with probability `58077/65536`; its
unconditional mean scalar cost is `32245/2048`. Both interventions have these
same counts and no wrong publications. The separate H=32 and H=64 path rows
have greater availability and stronger precision assumptions; they are exact
counting/bound certificates, not exhaustive native executions at those
horizons.

All `4^H` words remain in the finite experiment. A first publication at t
stands for all `4^(H-t)` suffix completions, which carry their original mass
but are not physically executed after stopping. Aborted histories pay the
whole `2H` means. There is no conditioning on success and no restart or free
source replenishment. Scalar reads/writes, cross-carrier means, preparation
and receiver samples are reported separately. The logical publication
controller and comparator are supplied instruments; their physical costs
are not assigned scalar-mean equivalents.

## What this establishes for the source contract

The [canonical axiom reference](../../docs/AXIOM_REFERENCE.md) assigns the
record domain and continuation interface to A1, checks accepted-data meaning
under A2, and requires a source-generated complete temporal grammar, exact
reference, observer cover and output map for an A3 policy claim. None of
these axioms permits adding successful routing as an unexplained constraint.

This construction derives the local partial output map from native
observational indistinguishability and proves its conditional service law.
Its acceptance predicate consumes observations produced by all permitted
histories, including failures. It does not remove the external choice of
paired alphabet, preparation, isolation, checkpoint interface or comparator.
No available source certificate in this construction supplies their full
A1--A3 selection or the radius-L/sqrt(q) read menu. The receipt therefore
sets both `m1_derived` and `full_axiom_source_contract` to false.

M1 additionally needs a general record/update service, source selection of
its metric menu and publication/control structure, resource bounds for whole
native workloads, and intervention-preserving q=13/q=21 compilation. Physical
clock and count-volume attachment are broader common-world requirements.
