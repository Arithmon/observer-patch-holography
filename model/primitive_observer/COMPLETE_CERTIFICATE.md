# Complete native and work output certificate

The certificate closes the numerical connection between the shared
[`runtime`](runtime.py), its [native attachment](native.py), and the exact
[prepared work histories](../common_source/SMOOTH_HISTORY.md). It compares
the actual retained native maps and work density matrices with complete
ideal outputs, including the native input's untouched reference. It covers
the prescribed, independently prepared work state; it does not extend its
vacuum theorem to an arbitrary native/work input.

[`complete_certificate.py`](complete_certificate.py) executes both native
laws C and E with protected and passive logical memory. It uses occupation
dimensions 2 and 4 and 512 product steps per Hamiltonian interval. For each
native current it executes both the empty environment and the existing
four-gap phase environment: variance `1/64`, order eight, scalar site zero,
after operation boundaries 1, 4, 6 and 9. Each environment invocation has a
fresh hidden draw. The receipt records outward guarantees, source bindings
and the resource ledger; replay regenerates the source operations and exact
work proof objects.

## Reproduce and inspect

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.primitive_observer.complete_certificate --check model/primitive_observer/complete_certificate_results.json
python -m pytest -q -W error model/primitive_observer
```

The [receipt](complete_certificate_results.json) contains twelve complete
work comparisons and sixteen native/work observer executions. The maxima
over both C/E laws and both recovery policies are:

| Occupation dimension | Executed environment | Distance to finite noiseless target | Distance to infinite noiseless target |
| --- | --- | --- | --- |
| 2 | Empty | 0.000001 | 0.247684 |
| 2 | Four phase channels | 0.021675 | 0.269357 |
| 4 | Empty | 0.000002 | 0.045012 |
| 4 | Four phase channels | 0.081372 | 0.126382 |

These are directed upper bounds on complete output half nuclear norm,
uniform over logical/reference inputs. Each numerical native attachment
contributes at most `10^-8`; the absolute complex total-mass defect is at
most `10^-7` for every such input. The finite and infinite totals above use
the finer internally calculated native bound before outward rounding.
The native contribution includes all encoder, repair, routing and recovery
contractions. No sampled set of native input states stands in for this bound.

There are **770** complete raw outcomes for C and **2,566** for E in each
execution. At dimension four, the noiseless protected E history has a
maximum of **4,535,334** native pulses on a compatible realized path; the
four phase channels add twelve pulses. The passive E history uses
**4,437,160** pulses without those channels. Exact componentwise counts and
a pulse-maximizing record are retained for every execution.

The noiseless dimension-four active-minus-passive phase-distinguishability
advantage is at least **0.259069** using the infinite target, and at least
**0.438957** using the finite target. The latter comparison can be tighter
because its common ideal reference has the same occupation cutoff as the
computed outputs. Both comparisons concern those same computed outputs;
they are separate certificates, not two different memory experiments.
With phase noise, the finite comparison still certifies an advantage of
at least **0.349223** at dimension two and **0.113479** at dimension four.
The receipt also keeps the vacuous infinite
comparisons instead of reporting their negative lower bounds as successes.

## One experiment and two comparison targets

The operation grammar, preparation, feedback, model-time coordinate and
occupation refinement remain the [canonical program](../common_source/COMPOSITION.md).
The shared interpreter compiles that program without changing its records:

```text
native logical encoding and first repair -> native current j
prepare work once, apply the j-controlled kick
evolve(1), gauge read -> r0, feedback(r0), evolve(1/2)
smooth read -> r1, evolve(1/2), smooth read -> r2, feedback(r2)
evolve(1/2), gauge read -> r3
native repair on the seam selected by r2, with declared recovery policy
```

The native/reference state may be any density operator on
`C^2 tensor H_reference`, for any finite reference dimension. Work is
prepared independently. Every raw native outcome is retained, including
reserved codes, aborts and padding output rows. Each non-aborted history
retains all sixteen work records and its disturbed work state. The four
readouts are successive measurements of one prepared work state; no
checkpoint introduces another preparation.

The finite target uses the exact finite continuous Hamiltonian at the
executed occupation dimension, exact instruments and exact native maps.
The infinite target uses the corresponding two-oscillator history and the
same native instrument, record dictionary and feedback. The finite output
is included in the infinite occupation space by the specified Fock-code
inclusion. Neither target contains the optional phase noise. Consequently
the noisy run's discrepancy includes both its numerical errors and its
departure from that noiseless target. It is not a certificate of a noisy
infinite model.

The final gauge bit is the existing model-time probe after the earlier
reads and feedback. Hamiltonian intervals total `5/2`; pointer parameters,
native attempts, source pulse counts and host computation time are separate
quantities. No physical clock calibration follows from this comparison.

## Exact data semantics

Let `a` denote the first native coin/raw record, `w=(r0,r1,r2,r3)` the work
record and `b` the final native coin/raw record. The complete stored output
is interpreted as the direct sum of

```text
(A_(a,w,b) tensor I) tau (A_(a,w,b)* tensor I) tensor R_(j(a),w).
```

Here `A` and `R` are the original binary array entries, interpreted as exact
complex rational numbers. Native maps have all 512 output rows. The first
abort rows use their retained prepared work matrix and have no later work
record. The ideal output has the same direct-sum labels, with native maps
`B` and work blocks `S`; ideal reserved/abort/padding contributions are zero.

No array is renormalized, projected onto positive matrices or replaced by a
nearby numerical eigendecomposition. The bound is a half nuclear-norm
matrix discrepancy even if a retained binary work matrix has a minute
non-Hermitian or negative component. It therefore does not assert that the
binary arrays define an exactly trace-preserving quantum channel. In
particular, the numerical bound is not capped at one on that assumption.

This interpretation certifies the factorized maps represented by the
retained arrays, uniformly in an exact normalized native/reference input.
A later floating multiplication of those maps into a supplied input or a
dense tensor product introduces new arithmetic. Such subsequent arithmetic
is not included unless its returned arrays are separately compared. The
certificate does include the source contractions that produced `A` and
every shared-runtime operation that produced `R`.

## Exact native targets

The independent logical code has columns supported on the occupation pairs

```text
logical 0: (0,2), (1,3)
logical 1: (0,3), (1,2),
```

with amplitude `1/sqrt(2)` on each support. Call its isometry `V`. The target
is derived from these integer configurations and the floor/ceiling native
endpoint rule, not by interpreting a floating target matrix as exact.

For either used seam `(0,1)` or `(2,3)`, every nonzero first-outcome effect on
this logical code is `I2/4`. For C the admitted `(coin,raw)` pairs are
`(0,1),(1,0),(0,0),(1,2)`; for E they are
`(0,5),(1,4),(0,4),(1,6)`. This gives the exact first-current weights

```text
p(-1)=1/4, p(0)=1/2, p(1)=1/4,
```

independently of the logical input and its reference. Protected recovery
restores the code after each native event, so each admitted complete native
map is `B=V/4`. Its four nonzero entries are `sqrt(2)/8`. Passive maps instead
follow both actual endpoint updates from each logical support; each path
contributes `sqrt(2)/4` to its endpoint entry. These are different target
instruments, with the same complete record contract.

All ideal entries therefore have the form `c*sqrt(2)` for rational `c`.
For a binary complex entry `x+i*y`,

```text
|x+i*y-c*sqrt(2)|^2 = x^2+y^2+2c^2-2*x*c*sqrt(2).
```

Summing this expression first and enclosing the remaining radical by exact
rational arithmetic bounds each Frobenius residual. This includes every
nominally zero entry and does not lose small errors through a floating
subtraction. Directed square roots give the norm bounds used below.

## Complete composition theorem

Write `D(X,Y)=||X-Y||_1/2`, with nuclear norm on all classical/quantum output
blocks. Suppose the independently recomputed work proof gives

```text
sum_w D(R_(j,w), S_(j,w)) <= epsilon_j,
```

for each current. Suppose also that exact directed arithmetic gives

```text
delta_r >= ||A_r-B_r||_F,
beta_r  >= ||B_r||_F,
nu_r    >= ||R_r||_1.
```

One available work bound is
`nu_r = sqrt_upper(D_work * sum_kl |(R_r)_kl|^2)`.
Define

```text
eta_native = sum_r delta_r*(delta_r+2*beta_r)*nu_r/2,
eta_work   = epsilon_-1/4 + epsilon_0/2 + epsilon_1/4,
eta        = eta_native + eta_work.
```

Then the complete output discrepancy is at most `eta`, uniformly over all
normalized mixed logical/reference inputs. First-abort rows contribute to
`eta_native` with `B=0`; they cannot be discarded because their computed
amplitudes are merely small.

To prove the native part, insert the intermediate output using ideal
native maps and the computed work blocks. For each record, the nuclear
norm inequality for products gives

```text
||A tau A* - B tau B*||_1
  <= ||A-B||_op (||A||_op+||B||_op)
  <= delta*(delta+2*beta).
```

The same inequality holds after tensoring the maps with any reference
identity; it uses `||tau||_1=1` and the operator norm of the native map, not
the Frobenius norm of its reference extension. Multiplicativity of nuclear
norm on tensor products and additivity on direct sums give `eta_native`.
This argument does not require positivity of the computed work matrices.

For the work part, the ideal native factors are positive. Holding the
first native record and selected seam fixed, completeness of the final
native instrument gives

```text
sum_b B_(a,r2,b)* B_(a,r2,b) = F_a* F_a = I2/4
```

for every admitted first record `a`, independently of `r2`. Summing final
native outcomes therefore weights each work-block discrepancy by `1/4`.
Two first records have current zero and one each has current minus or plus
one. This proves `eta_work`. There is no extra factor equal to the number
of work records or native outcomes. The theorem covers the complete
feedback connection precisely because both sides use the same selector
and the final instrument is complete for either selected seam.

The general prepared-work theorem for an arbitrary native carrier instead
uses its first-current probabilities and is bounded by `max_j epsilon_j`.
The encoded weights above must not be silently applied to a general
78-dimensional native input. Initially correlated native/work inputs need
an all-input work-channel comparison, which this prepared-output theorem
does not supply.

## Error and resource ledger

The work comparison reuses
[`exact_output.compare_densities`](../common_source/exact_output.py).
It converts each original binary matrix entry to a rational number and
compares the complete matrix with the exact polynomial history, including
the polynomial's entire occupation support and its proved remainder. The
finite-to-infinite comparison reuses the existing exact history enclosures.
These are prepared-output bounds, not sampled matrix residuals.

The triangle through the exact finite history accounts for the native
residual, actual work computation and occupation difference once each.
Because the actual work arrays are directly compared with the continuous
target, their source splitting and coefficient errors are already present.
The ideal temporal estimate must not be added again. In the noisy branch,
the same direct comparison also includes the specified environment's
change to the noiseless output; it does not separate noise from arithmetic
without another theorem.

Source resources are counted on compatible realized paths by
[`native.resources`](native.py), including the logical encoder, both native
repairs, branch-dependent recovery, work preparation, every repeated source
operation, pointers and environment invocations. Componentwise path maxima
need not occur on one common path. Expanded simulation-tree resources are
separate from the cost of one realized observer history. Reusing a computed
operator does not remove any repeated physical invocation from the count.

Blank initialization, coherent code transfer, record storage, readout and
physical operation fidelity remain supplied services. Missing physical
stage errors remain explicit obligations with a null complete physical
bound; they are not replaced by zero. Given valid preparation
trace-distance and complete-stage half-diamond service bounds, their
sum can be added by the existing telescoping and contractivity argument,
counting each executed invocation and taking a uniform bound over feedback
choices. This conditional arithmetic does not establish those premises or
infer physical durations from source counts.

## Memory and finite-theorem attachments

For ideal code inputs with relative phase plus and minus, the protected
complete output distinguishability is one. In the passive branch it is
`1-qbar`, where

```text
qbar = (q_-1+2*q_0+q_1)/4,
q_j  = probability of the second smooth record r2=0 at current j.
```

These are the existing native-memory identities with the work probabilities
from the same continued history. If both computed phase outputs are each
within `eta` of their ideal target, the difference between computed and
ideal distinguishability is at most `2*eta` by two applications of the
triangle inequality. An exact interval for `qbar` propagates through the
passive formula. A loose occupation bound may make a resulting lower
bound vacuous; the identity alone does not then certify a resolved
numerical separation. No destructive logical read or decoder is inserted
by stating this metric.

The finite [Kraus composition result](../../Lean/Dynamics/ChoiCPTP.lean)
supplies the ideal instruments' complete positivity and continuation; the
work proof supplies the prescribed preparation's quantitative output
comparison. The proof above connects them through actual record-controlled
native continuation. The
[public-record algebra](../../Lean/EventAlgebra/PublicRecordAlgebra.lean)
acts on the very same retained direct-sum labels, so classical feedback
does not erase or rewrite an earlier record. These are compatible finite
attachments. A growing direct-sum record algebra requires an explicit
fixed-tape embedding before using a theorem stated only for fixed square
matrices.

C and E remain alternative quantum instruments. Protected and passive
memory remain different recovery policies. Empty and noisy environments
remain different executed channels. Agreement of their classical native
load laws does not identify their quantum outputs. The uniform-scheduler
[transfer reconstruction](DERIVATION_TRANSFER.md) concerns another declared
native schedule and does not supply the work action or calibrate this
observer's clock.

This certificate completes an error-controlled finite native/work observer
output, with an explicit occupation comparison. It does not instantiate
the full [`OperationalObserver`](../../Lean/Tower/OperationalObserver.lean)
consensus tower, exact multiplicative spatial refinement, a physical
clock, a source-selected interacting action or an empirical identification.
Those obligations remain distinct from the numerical connection proved
here.
