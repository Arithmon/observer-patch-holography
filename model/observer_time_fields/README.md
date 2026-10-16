# Fields in the observer's own incident clock

The native repair source determines a positive effective transfer that
preserves **complete local field histories**, including a declared family of
record-controlled local interventions, when time is counted by incident
repair attempts. No duration per repair, oscillator action, or fitted
Hamiltonian enters this construction. This connects the
[incident-clock resolvent](../local_tick/README.md) and
[source exchange field](../native_effective_fields/README.md) at the level of
joint record laws, rather than identifying their separate two-point spectra.

The [joint clock process](../joint_clock_process/README.md)
preserves complete quantum histories through explicit memory on a controlled
four-port experiment. The [native action](../native_action/README.md) derives
a local Liouville generator in global attempt or auxiliary Poisson units.
Here the additional connection is a **positive transfer and its logarithm
in the incident observer clock**, preserving its local field dictionary.
The full-federation execution uses the uniform schedule on all 630 seams.

The executable instances are the original twelve-port, two-particle source
and the actual twenty-carrier federation: 240 ports and all 28,680 binary
two-particle configurations. The latter retains remote repairs, exclusion,
cross-carrier feedback, every initial readout label and finite-cap timeout.
It does not use a one-particle or mean-density closure for higher moments.

## The common observer contract

This is a classical recorded branch of the
[primitive observer specification](../primitive_observer/README.md), with
the same native seam operation and the existing source federation. C and E
have the same diagonal transition law here; this packet does not select their
coherent continuation.

| Interface | Instantiation |
| --- | --- |
| State | Binary occupations of the complete source graph at fixed total two. The continuing distribution includes the entire graph, including remote ports. |
| Preparation | The uniform distribution on that sector, with the two read-port occupations initially disclosed. This preparation is an explicit input. Conditioning on its disclosed labels retains all four outcomes. |
| Repair | Select a source seam uniformly and apply the existing fair-coin native operation. Its unread diagonal channel is `(I + SWAP_ab)/2`. |
| Clock | A scheduled repair ticks iff its seam meets the declared observer region, even when its current is zero. Remote repairs evolve the same state between ticks. |
| Field | Occupation of a read port, reconstructed from its disclosed initial value and every actual signed incident current. Endpoint projectors in the probability calculation represent this record processing. |
| Read sequence | Alternate the two read ports for four scheduled incident hits. Initial labels, successful bit words and every timeout prefix are retained. |
| Feedback | After the first field read, bit one requests one additional fair native repair on the seam joining the read ports; bit zero requests none. This uses a committed bit, never an inaccessible state. |
| Clock accounting | The feedback attempt is itself an incident tick. Read times are `(1,2,3,4)` without feedback and `(1,2+b,3+b,4+b)` with first bit `b`. The scheduled gap index is not silently called the complete clock. |
| Continuation | Every subsequent repair acts on the conditional distribution left by the preceding read and intervention. There is no reset at a tick. |
| Finite execution | Each scheduled interval has a cap on global scheduler attempts, enforced by a supplied counter/timeout service. A timeout terminates that branch with its preceding labels and fields intact. Arithmetic deficits are bounded and never normalized away. |
| Refinement | Increasing the attempt cap and integer precision converges to the complete stopped law with explicit error. Spatial carrier refinement is a separate obligation; equal clock units alone do not identify different graphs' laws. |

For the elementary source the clock region is the two read ports. For the
federation it is the two complete carriers containing the first retained
inter-carrier seam. The region includes 24 of the **actual 240** ports; it is
not an invented two-carrier universe. The two instances share the theorem
and controller, not an asserted equality of their predictions.

## A complete field-law reconstruction

All matrices below act on column distributions on the finite configuration
space. Write the uniformly scheduled unread native transfer as

```
P = A + B,
A = remote-seam contribution,
B = incident-seam contribution,
p = number of incident seams / number of all seams,
D = I - A.
```

Each seam contribution is a nonnegative symmetric positive semidefinite
matrix, since `(I+SWAP)/2` is an orthogonal projection. Consequently
`0 <= A <= (1-p) I`, `B >= 0`, and `A+B <= I`. Both the column and row sums
of `A` are `1-p`; those of `B` are `p`. Thus `p I <= D <= I` and the
completed incident-tick transfer is

```
K = sum_(n>=0) B A^n = B D^-1.
```

The remote contribution is an actual transfer, not a scalar waiting factor.
In general `K` is not symmetric. Define instead

```
T = D^-1/2 B D^-1/2,
K = D^1/2 T D^-1/2.
```

Then `0 <= T <= I`. Moreover `T 1 = 1`: `D 1=p 1` and `B 1=p 1`.
The convergent binomial series for `(I-A)^-1/2` has nonnegative coefficients,
so `T` is entrywise nonnegative as well. It is a symmetric stochastic
transfer derived entirely from the source and the chosen incident clock.

Let `F` be any multiplication operator determined by occupations in the
observer region. A remote repair fixes those occupations, hence `[F,A]=0`
and `[F,D^s]=0` for every real `s`. The same commutation holds for any native
repair supported entirely inside the region, since its ports are disjoint
from all remote seams. In particular it holds for the feedback repair in
the executed protocol.

Set `Omega=1/sqrt(number of configurations)`. In every stationary field
cylinder, replace each `K` by `D^1/2 T D^-1/2`. The adjacent similarities
cancel through every local field projector and every allowed intervention.
The remaining boundary factors cancel because
`D^s Omega=p^s Omega`. Therefore **each complete word probability is
unchanged**, for every word length. This includes the initial disclosed
labels, all zero-probability words and any causal record-dependent choice
from the commuting intervention family. The proof applies branch by branch;
it does not replace feedback by an averaged uncontrolled channel.

If incident seam, coin or current marks are retained, write their native
contributions as `B_m`. The transformed insertion is
`T_m=D^-1/2 B_m D^-1/2`; then `sum_m T_m=T` and the same cancellation proves
the corresponding marked word identities. Individual marked insertions need
not be symmetric or be field projectors. This extension does not identify
the effective hidden trajectory with the source's individual attempt times.

The identity is exact for the stated stationary preparation and the
record-induced conditional states, which are already represented by the
initial projectors. It is not a claim that an arbitrary other preparation
has the same untransformed boundary vectors. Likewise, boundary-straddling
interventions or fields outside the clock region do not automatically
commute with `D`; the retained negative controls test that distinction.

## A Hamiltonian with a source-derived time unit

Strict positivity is established separately; positivity alone would not
justify a finite logarithm. On the full federation sector each triangular
face of a carrier obeys

```
sum_(e in triangle) (I+SWAP_e)/2 >= (3/2) I.
```

To see this, fix the spectator occupations. With zero or three occupied
triangle ports all swaps act as identity. With one or two, their permutation
representation splits into the constant line and the two-dimensional
standard representation, with projection-sum eigenvalues `3` and `3/2`.
Each of the thirty edges of an icosahedral carrier belongs to two of its
twenty triangular faces. Summing gives a lower bound `15 I` per complete
carrier. The two clock carriers therefore give
`B >= (30/630) I = I/21`; additional incident seams are positive.
Since `D <= I`, also `T >= I/21`. The elementary instance supplies its own
exact finite lower-bound check.

It follows that

```
H = -log T,
H_m = sum_(j=1)^m (I-T)^j/j,
0 <= H-H_m <= (1-alpha)^(m+1) / ((m+1) alpha) I
```

when `T >= alpha I > 0`. This is a norm-controlled source reconstruction
in incident-clock units. `exp(-n H)=T^n`, with the same local read and
intervention insertions, reproduces the complete declared Euclidean field
histories, not just an energy gap. The Hamiltonian will
generally be nonlocal in the original graph coordinates because it includes
the eliminated remote waiting histories.

The mathematical continuation `exp(-i s H)` exists, but the native repairs
are stochastic operations, not that unitary. This packet neither supplies
an experimental coherent implementation of the new Hamiltonian nor
identifies its real continuation with physical native waiting. The
[compiled transfer probe](../transfer_modular_probe/README.md) illustrates
how such a separate coherent apparatus contract can be specified.

Independent identity holding also cancels from the reconstructed time and
field law. If a scheduler inserts a no-operation slot with probability
`1-a`, where `0<a<=1`, then `A_a=(1-a)I+aA`, `B_a=aB` and `D_a=aD`.
Both `K_a=K` and `T_a=T` follow exactly, so the effective `H` is unchanged.
The inserted slots are not incident repair attempts. Their physical cost,
attempt timestamps and finite-cap timeout probabilities are not invariant.
Changing relative seam rates, or doing other work during a holding slot,
is outside this identity.

## Why the clock choice matters

The existing [local movement clock](../relational_event_clock/README.md)
ticks only when one port changes occupation. Its successive completed ticks
alternate that port's binary value, so `2 n_port-1` is an eigenfunction with
eigenvalue `-1`. A positive Euclidean transfer cannot have that spectrum.
This is an explicit obstruction to applying this logarithm construction to
that clock. The positive incident clock counts an admitted repair attempt,
including its zero-current outcome; it does not discard these attempts to
make the clocks appear interchangeable.

## Arithmetic, stopping and resources

The small source uses exact rational matrices and complete field laws.
The federation implementation keeps the full two-particle sector with
sparse integer source operators and directed integer-grid propagation.
Every rounded component is a lower bound on the corresponding exact
unnormalized probability. Lost mass is accounted for across the entire
branch tree, including initialization and feedback. The final deficit thus
bounds the complete classical output's L1 error, not merely a fitted moment.
No floating eigensolve or floating normalization contributes to this bound.

For cap `h` on each of four intervals, the ideal complete timeout probability
is `1-(1-(1-p)^h)^4`, since the seam menu is state independent. The separate
arithmetic deficit bounds the computed timeout and success masses. The
completed law differs from the capped experiment by at most the timeout
probability in total variation, when a timeout is replaced by its eventual
continuation. A moment or nonlinear cumulant must additionally propagate
the appropriate complete-law uncertainty; a nonzero floating number alone
is not a certificate.

The cap counts remote as well as incident scheduler attempts. It is a supplied
execution/resource cutoff, not a duration read from the incident clock.
Implementing the capped protocol requires that global attempt counter and
timeout service; the completed incident-clock law requires neither.

The capped protocol uses at most `4h+1` native attempts, four scheduled
field reads, two disclosed initial bits, and at most one feedback attempt.
Its incident clock uses at most five ticks; remote attempts are still
charged as source operations. Processing a signed incident current is a
declared record-access service. The implementation also reports its finite
configuration and terminal-record storage. These are software/resource
contracts, not measured hardware service-error estimates.

At the retained cap 256 and grid precision 51 bits, each protocol has 124
terminal records and at most 1,025 native operations. The baseline complete
L1 arithmetic error is below `1.59e-8`; the feedback error is below `1.69e-8`.
The exact timeout probability is below `3.13e-12`. The eventual baseline
connected four-field cumulant lies inside `[-0.000027711,-0.000027693]`; the feedback
increase of the second field mean lies inside `[0.00392579,0.00392584]`.
These decimal intervals are outward summaries of the rational receipt.
Non-Gaussianity here is a finite classical field-law prediction; by itself
it establishes neither quantum dynamics nor an empirical interaction.

## What this completes, and what remains

The common source has an exact connection from its repair law to
observer-relative incident time and complete effective local field laws,
with compatible feedback and executable finite approximations. The full
federation demonstration tests exclusion and cross-carrier records together.
This completes a finite clock/field connection missing from the separate
clock and global-attempt field packets.

The architecture, seam scheduler, initial particle sector and record access
remain the stated source-instance inputs. Relating different observer
regions or spatial refinements still requires a shared record map and a
joint-law comparison; the existing
[coarse memory result](../federated_observer/MEMORY.md) supplies mean dynamics,
not those complete laws. No Lorentzian continuum, universal synchronization,
physical time calibration, C/E selection or new empirical observation is
claimed. The construction changes none of the retained laboratory fits.

Reproduce the complete packet with:

```sh
python -m model.observer_time_fields --check model/observer_time_fields/results.json
python -m pytest -q -W error model/observer_time_fields
```
