# Source histories with noisy records and control

The [fixed-rate construction](FIXED_RATE_SOURCE_READS.md) protects quantum
interfaces but explicitly leaves central records reliable. This note removes
that runtime assumption **inside the same constructed source**. Decisions run
on protected registers; a separate, explicitly maintained noisy classical
archive exposes their results. The joint named history, including rejected
outcomes, converges. Source selection and an autonomous physical clock are
not consequences of this engineering result.

The new ingredients are an exact instrument translation, a protected-to-public
interface, an arbitrary-size archive fault proof, and a lifetime resource
ledger. The executable checks are in [m1_noisy_records](../code/m1_noisy_records/README.md).

## 1. What changes, and which standard theorem is imported

Retain the proper-code M6 processors, coherent transfers, finite classical
Boolean services, resets and cofinal finite inventory/control capability of
[the source construction](COHERENT_SOURCE_CLOCKS.md). Extend its independent
Poisson dephasing by a fixed finite local jump rate mu, allowing bit errors
in stored records and faults during Boolean operations. On quantum carriers
the additional jumps preserve the native coordinate codes; a bit flip is an
example. A jump acts on a stored carrier or the operands of one native service,
and may corrupt all that service's operands. Disjoint carrier-time/service
groups have independent Poisson drivers. Thus a hit belongs to one counted
circuit location, including its bounded native helpers. Bounded spatial range
alone would not imply this condition: a common burst hitting two simultaneous
locations has probability O(p), not O(p^2), and is outside this independent
rate model. Rates lambda,mu do not decrease with q.

Use **Aharonov and Ben-Or**, [quant-ph/9906129](https://arxiv.org/abs/quant-ph/9906129),
Theorem 4 and sections 7.1--7.6: their CSS computation-code construction
requires no noiseless intermediate measurement or classical operation.
Choose a fixed computation code correcting at least its gadget spread.
Its finite circuits, blank supply, correction and universal gates are the
imported construction. The sparse-fault argument gives a bound
`p_k <= A^(2^k-1) p^(2^k)` when at least one fault per rectangle is tolerated.
Here A is a finite library constant, enlarged for the native implementation.
We do not use their numerical threshold, whose quoted optimized estimate
assumes reliable classical processing. Their G1 gates have exact finite
native Clifford+T implementations; the Toffoli decomposition is checked
including its phases. Requested variable rotations are synthesized into the
logical G1 family on bounded-size registers, using the credited
[constructive Solovay--Kitaev synthesis](https://arxiv.org/abs/quant-ph/0505030).
The exact native implementations then introduce no constant approximation
floor at the bottom of concatenation. A native T pulse is used inside those
implementations, not assumed to have a transversal logical G1 gadget.

The seven-qubit recovery census in the parent is **not** a certificate for
these different measurement-free gadgets. In particular, its one-error code
cannot simply be reused with a larger-spread coherent recovery. The imported
code below satisfies the published spread condition with a conservative
composition margin. Its constant size and complete construction, rather than
an invented native threshold, are needed.

All runtime decisions, addresses that depend on data, scratch registers and
counter arithmetic are compiled into protected gates. The fixed experiment
and the fixed order of pulses are spatially distributed in the parent's
finite preparation prehistory. They specify the circuit, not a mutable
uncoded broadcast instruction register. An independently derived or randomly
drifting global clock is not supplied here, just as it was not in the parent.

### An explicit computation code, with a distance proof

Use the standard binary Reed--Muller evaluation code RM(5,11): evaluate every
multilinear polynomial of degree at most five on all 2^11 Boolean points.
Its monomial generators are explicitly constructed in two independent ways
by the executable evidence. This is a conservative mathematical code choice,
not a new code or a practical hardware proposal.

The monomials are independent by Boolean Mobius inversion. There are
`sum_{j=0}^5 binom(11,j)=1024` of them. Two generators overlap on
`2^(11-|S union T|)` points, an even number since `|S union T|<=10`.
Thus the 2048-bit code is self-orthogonal and, by its half dimension, self-dual.
Each generator has weight divisible by four; the even pair intersections and
`wt(a XOR b)=wt(a)+wt(b)-2wt(a AND b)` make the entire code doubly even.

For completeness its distance follows without assuming a database entry.
Write a nonzero degree-r polynomial on m variables as `f=p+x_m t`. If t=0,
its weight is twice that of p. If t is nonzero, the two halves satisfy
`wt(p)+wt(p+t)>=wt(t)`. Induction on m therefore bounds every nonzero weight
below by `2^(m-r)`; a degree-r monomial attains it. The base cases r=0 and
r=m have distances 2^m and 1. Here the distance is exactly 64.

Puncture the all-ones evaluation point. The 2047-bit code C still has dimension
1024 and distance 63, witnessed by a degree-five monomial. Its shortened
subcode S consists of the words zero at that point, has dimension 1023,
and equals C-perp. The CSS stabilizers X(S),Z(S) encode
`2047-2*1023=1` qubit. All nontrivial logical Pauli supports lie in C outside
S, have weight at least 63, and the punctured monomial attains 63 outside S.
This gives the explicit **[[2047,1,63]] code correcting 31 errors**. The two
logical cosets have weights 0 and 3 modulo four, as required by the cited
punctured doubly-even construction.

For safety, combine the cited per-procedure spread bound four with a separate
leading correction by the conservative product bound `4*4=16`: a fault in
the correction can leave four errors, each spreading to at most four positions
in the subsequent ideal gate procedure. A fault in that gate has spread at
most four. Separate leading corrections use separate ancillas. Thus 31
correctable errors exceed the combined bound 16; tolerate one bad subrectangle
and classify two as potentially bad. This directly justifies the quadratic
recurrence used in (7), without transferring the parent's seven-qubit claim.

The verifier checks the full 1024-generator rank, all generator inner products,
punctured/shortened ranks, phase classes and a minimum-weight witness. Small
instances include a complete 65,536-word distance check. The all-size distance
statement is the induction above; no claim is made to enumerate 2^1024 words
or to execute the entire universal gadget library. The chosen constant block
size is intentionally large, and all five-level block and gadget factors
remain charged constants in the existence bounds.

## 2. Exact causal history compiler

For a finite primitive instrument, a public outcome y may combine several
private Kraus branches alpha. Write
`Phi_y(rho)=sum_alpha K_(y,alpha) rho K_(y,alpha)*`, with
`sum_(y,alpha) K_(y,alpha)* K_(y,alpha)=I`. Its channel dilation is

```
V psi = sum_(y,alpha) K_(y,alpha) psi
          tensor |y>_R tensor |y,alpha>_E.                  (1)
```

Tracing E leaves a classical public R, possibly correlated with the application
and an arbitrary reference. The private alpha is not a new public record or
feedback input. Coherently adding the K_(y,alpha) for one y would generally
change the channel and is forbidden.

On the accessible proper-code logical registers, implement the finite
Stinespring completion on blanks and discard E through the charged service.
For computational-basis reads this is two CNOTs and a discard. The dilation
of a native full-M6 extension is a channel representation, not permission to
perform a reversible operation on its inaccessible complement. That complement
retains the original complete CPTP transfer/reset service and its public failure
outcome. Prepared source inputs and the specified code-preserving jumps keep
the protected computation in its native coordinate codes. This construction
does not promise protected encoding of arbitrary full-M6 input states. The
complete extension is nevertheless checked, rather than silently dropping its
failure branches or renormalizing its successful subchannel.

A later classical Boolean instruction is made reversible on fresh workspace
using NOT, CNOT and Toffoli; copy its answer before discarding scratch. A
conditional quantum instruction is controlled by its logical answer register.
Each source primitive has bounded arity. Variable identifiers occupy the
already counted O(log q) bits and use bit circuits, not a single huge oracle.
Bounded source retries retain each attempt and use a padded selection circuit.
Lower-level protection syndromes are local diagnostic outputs, not global
abort commands. The fixed finite imported gadgets have no unbounded retry.

**History identity.** At every completed source-instrument cut, tracing the
discarded registers gives exactly

```
sum_h |h><h|_R tensor Phi_h(rho),                           (2)
```

where Phi_h is the original ordered composition of the outcome CP maps,
including the interventions selected by the preceding public labels. Proof: (1)
gives the first cut. A Boolean permutation maps diagonal register blocks to
diagonal blocks; its retained copies prevent the merging of distinct histories.
A controlled channel acts on each block by its corresponding branch. Applying
(1) appends the next outcome and sums its private Kraus branches incoherently.
These facts prove the induction, also after tensoring with any reference identity.
They prove equality of causal instruments, not just terminal probabilities.
An intervention at a cut may be any declared logical channel controlled by
the public label; a finite intervention program is compiled and charged by
the same rule. An unprotected external controller has no free noise guarantee.

Keep the logical record until its last feedback use; never uncompute it across
an observed cut. All control transfers follow the parent's charged neighbor
messages or local helper routes. Compilation creates no backward dependency
or instantaneous remote decision. Discards supply the usual entropy sink;
this is not a closed unitary apparatus with a one-time blank supply.

The finite control implements two noncommuting reads, a Boolean AND abort,
conditional H, X and Z, and a rotation applied only to accepted branches.
All four histories, including abort, have positive nonzero CP weight. Forward
isometries and independently formed projector instruments agree on their
entire Choi matrices, including the abort register and all input coherences.
Deleting the environmental copy, changing feedback, or removing a Toffoli
phase fails these tests. The analytic induction, not these three examples,
establishes (2) at arbitrary size.

The additional native checks use the parent's actual coherent code transfers
for both proper-code dimensions d=2,4. Their success operator is the identity
on the d-dimensional code; the 6-d complementary basis inputs all reset to
one failure state with **one shared public failure label**. The independent
oracle checks all 36 input matrix units, including code/complement coherences.
It rejects missing failure branches, coherently merged private branches,
residual coherence between public success and failure, and a TP map with the
wrong failure output. A grouped qubit-reset instrument with subsequent public
feedback is checked on every input matrix unit too. Unitary changes of the
private Kraus basis leave the channel test unchanged, as they should.

## 3. A noisy public archive with a proved refresh circuit

For integer r>=1 let n=4r+1. A public record is n ordinary noisy central bits
whose majority is its label. At the start of a refresh it has at most r wrong
bits. Allocate n private n-bit voter lanes and n blank output bits.

1. Reset the private lanes and outputs. Copy each input once into each lane,
   in n conflict-free rounds. A copy touches only one shared input and one
   private destination; it never touches two shared inputs.
2. In each lane sort its binary word with the fixed bubble sorting network.
   Each two-bit comparator maps `(a,b)` to `(a AND b,a OR b)`. Native finite
   Boolean operations and fresh scratch implement it; every such bounded
   implementation is one enlarged fault location. All gates stay in that lane.
3. Copy each lane's middle bit to its own output. At a fixed scheduled cut,
   designate the output word as the current archive word. Reuse old workspace
   with charged resets. Count idles on every live input, workspace and output.

For an explicit comparator, start two local blank bits c,d, apply
`Toffoli(a,b;c)`, then `CNOT(a;d), CNOT(b;d), CNOT(c;d)`. They hold AND and OR.
Expose c,d as the next two lane values; old a,b are local garbage and can be
reset with their outcomes retained. The finite test executes these gates on
all four inputs using the independently checked native Toffoli decomposition.
Helpers and resets remain in the lane and enlarge native counts by a fixed
factor. They never touch a shared archive input.

There are `W=n(n+2)` live service bits, up to this fixed native helper factor, and
`D=2+n+n(n-1)/2` layers. The number of service locations is at most
`V=W D=O(n^4)`. The implemented n=5 circuit has exactly 515 locations.
Pad the start of each cycle to its coarse period: the first-layer input idles
include the entire wait, and each private reset location includes its preceding
wait and reset. Their no-hit probability is bounded using that whole interval.
This accounts for storage between refreshes without an uncounted pause. The
two input/output banks alternate by a fixed schedule; the designation is not
a mutable unprotected selection bit.

**Refresh lemma.** If at most r service locations are faulty, at the end at
most r output bits are wrong. Throughout the cycle the old input word still
has the correct majority.

Proof: each faulty copy or input idle damages at most one shared input
position. There are therefore at most 2r wrong shared positions, including
the r old errors. Mark a voter damaged if a fault touches its private lane or
output. At most r voters are marked; a faulty copy can mark a lane as well
as damage one shared input, which the bound already allows. Every unmarked
voter received at most 2r wrong inputs, even if an input changed between
copies. Its majority of 4r+1 is correct and its output is correct. Input
majority is also correct because its error count never exceeds 2r. This
argument allows arbitrary operand corruption at faulty locations, not only
independent output-bit flips. Bubble sorting supplies exact majority at every
n; private-lane locality, not a favorable fault sample, proves the lemma.

The exhaustive n=5 replay covers both labels, all zero/one-error inputs and
every nonidentity output corruption of each reset, copy, comparator and idle:
**8,100 cases**, with 6,365 clean outputs and 1,735 one-error outputs. Both
implementations also track the maximum error count of the old public word
after every location: 1,130 cases reach zero, 6,090 reach one, and 880 reach
two wrong input bits; none loses its live majority. Two initial input-idle
faults added to one old error give a concrete failure with three wrong live
inputs and all five output bits wrong.
There is no assertion that two faults must always fail. Larger-instance tests
also apply multiple faults; the arbitrary-r statement is the proof above.

Under local stochastic noise of strength p, losing the invariant in a refresh
requires at least r+1 faulty locations, so

```
Pr(bad refresh) <= binom(V,r+1) p^(r+1) <= (V p)^(r+1).      (3)
```

No independence between successive refresh successes is needed: union-bound
their bad-location events, including all storage intervals. The old word
remains public until the fixed handoff. Its majority is thus available live,
not manufactured by a final reconstruction of the experiment's history.

## 4. Export from a protected quantum record without a fragile bottleneck

One must not decode R to a single bit and then fan that bit out. A single
fault there can corrupt every copy; the finite countercontrol exhibits all
seven target bits flipped by one corrupt broadcast bit.

Instead, while R is protected, prepare n independently protected blank probe
blocks and use protected logical CNOTs from R to each probe. R is diagonal
as in (2), so this copies classical information and preserves the exact joint
application/history state. Include each copy and its recovery in the protected
circuit fault bound. Then decode each probe **in its own disjoint island**:
a finite complete correction/decoder for the fixed five-level code followed
by a native classical read. This final island is allowed to be faulty. Its
only output is one archive bit and it never touches R or another island.

On a sparse protected fault path, a noiseless island recovers the right label,
even if the probe has a correctable residual error. All classical read/storage
faults before the first refresh belong to that island or its counted idle.
There are only constantly many locations per island, so wrong output bits
have local stochastic strength `p_export <= C_export p`. Failure of an
export to start within r errors requires r+1 bad islands. Thus

```
Pr(bad export) <= (n C_export p)^(r+1)                     (4)
```

apart from the already counted non-sparse protected fault event. Do not
condition the noise distribution on protected success: bound that bad event
separately, and use the unconditional disjoint-island fault events for (4).
Protected-copy faults that could corrupt the shared source are therefore not
silently treated as independent output-bit faults.

The archive is made of actual central classical bits. Private protected R
drives feedback; no physical majority decoder of the archive is trusted by
the quantum computation. Public statistical decoding names the stored record,
as with any coded communication format. A supplied finite reader can itself
be compiled into the protected computation, but its final bare bit has the
unavoidable boundary in section 7. Observing a public majority is a classical
function of the available word and cannot disturb the application state.

## 5. Paid source schedule and lifetime ledger

Fix physical L=qa, horizon, masses and read kernels as in the parents. Set
`r(q)=ceil(log2(2q))`, `n(q)=4r(q)+1`. The parent has O(q^4) primitive source
events; gate synthesis, identifiers and bookkeeping have the conservative
envelope O(q^4 log^6(2q)), depth O(q log^6(2q)). The latter follows from its
neighbor-prefix preparation, fixed local walk layers and diameter-time
report routes, with no serial traversal of all cells.

Reversible Boolean compilation has bounded overhead per bit operation. Each
public bit export adds O(log q) protected copies and independent constant-size
islands. Each live classical bit record adds O(log^2 q) private workspace.
All archive lanes refresh in parallel, locally, with depth O(log^2 q).
Keeping every named record and conservative padding therefore gives

```
live active wires W_q <= C_W q^4 log^8(2q),
scheduled depth D_q <= C_D q log^8(2q),
active spacetime locations M_q <= C_M q^5 log^16(2q).        (5)
```

The extra q in volume is essential: old records persist for O(q) coarse
periods. The quantum-only O(q^4) ledger cannot be reused. Five levels of the
fixed universal protection library add finite constant factors; the logarithmic
classical archive, all resets and all idle intervals are already in (5).
Bounds intentionally overcount workspace and idles. Fresh-blank/reset events
are bounded by the same active volume; no entropy-erasure energy is inferred.

There is a second lifetime cost. Native protection and refresh services emit
their own diagnostic/version records, in addition to the named source history.
Store them in separate passive central carriers, with at most O(log q) bits
per service for its identifiers. Their total inventory is
`R_diag=O(M_q log q)=O(q^5 log^17 q)`. Keeping every one to the final cut and
counting its idles gives **total physical storage-time volume
O(q^6 log^25 q)**. These diagnostic bits are noisy; their distribution is not
claimed equal to noiseless diagnostics. They are never consulted by control
or by the named archive's refresh. Local TP noise on an unused factor cannot
change the active reduced state, even if the factors are correlated. Therefore
they enter the physical resource ledger but need not enter the operational
bad-event union in (7). This is a proved no-feedback distinction, not omission
of their exposure. Do not protect every protection diagnostic recursively:
that would be a different, potentially unbounded requirement. Idle storage
noise does not itself emit a new explicit program record at every instant.

Allocate separate processors for simultaneous lanes and islands. Let B_q be
the actual longest local compiled depth, including native decompositions and
endpoint recovery; `B_q=O(log^8(2q))`. Use

```
eta_q = a/(1024 c B_q),   Omega_q = pi/eta_q,
local cluster radius <= c eta_q/128.                       (6)
```

Local shared processors are serialized; separate processors operate in
parallel. Finite per-owner inventories include (5) and the passive diagnostic
carriers just counted; no bounded spatial
density of hardware is being assumed. Parallel constituents of an encoded
flight cross the original link once, with recovery at its endpoints. Each
flight still takes at least its length/c. No archive service introduces a
remote edge or additional serial long-distance traversal. Fixed padding of
the original source slots accommodates (6), so wall time retains a finite
q-independent rescaling. All classical records receive at least one complete
refresh per coarse a/c interval, including during long block flights.

The Poisson no-hit decomposition used in the parent now gives
`p <= C_native (lambda+mu) a/c` for every gate, read, reset, flight and idle
location. Internal slots are shorter; using a/c only enlarges the bound.
For every specified finite set S of locations the disjoint-group construction
gives `Pr(all S faulty)<=p^|S|`. Later propagation of a fault through an ideal
gate does not mark that later gate faulty; the spread and archive proofs count
that propagation explicitly. This is the local stochastic hypothesis needed
in both (3) and the imported quantum theorem.
Classical comparators have fixed native decompositions; their internal scratch
and faults are included in C_native. Code-preserving jumps avoid an unhandled
quantum leakage channel. Complete complementary effects remain present for
general native inputs, as before.

## 6. Joint operational refinement theorem

For fixed finite lambda,mu and positive L,c, p=O(q^-1). At five levels the
probability of any non-sparse protected rectangle is

```
epsilon_quantum <= M_q A^31 p^32 = O(q^-27 log^16(2q)).     (7)
```

The constant-size published sparse-fault proof works with reference-entangled
inputs: exact correction holds on all matrix units of the code, and tensoring
with an untouched reference preserves that identity. The remaining bad-path
weight bounds trace distance, using the convention one half of trace norm.
This supplies the decoded-state/instrument bound needed here, in addition
to the classical-output formulation of the cited theorem. Known noisy blank
preparation is part of the construction; arbitrary unprotected unknown-input
encoding is not free.

There are O(q^4 log^6 q) record exports and at most O(q) coarse refreshes per
record. Combining (3)--(4) gives, with fixed constants,

```
epsilon_archive <= C q^5 log^8(2q)
                       [C' log^4(2q)/q]^(ceil(log2(2q))+1).
                                                            (8)
```

This is smaller than every inverse power of q on a cofinal tail: its logarithm
is `-(log q)^2/log 2 + O(log q log log q)`. It covers the **joint live history**,
not merely the error probability of one bit. The invariant holds throughout
every successful refresh, so taking any scheduled causal cut does not incur
an uncounted inference from a later record. Exports are padded to the
corresponding source cut; earlier private work is not advertised as a completed
public read. All associated latency is in (6).

Synthesize each requested source rotation to q^-24. There are at most
O(q^4 log^8 q) requests, giving total error O(q^-20 log^8 q); the logarithmic
synthesis length is already in B_q. Exact Boolean, Clifford and Toffoli
gates introduce no approximation floor. Contractivity and telescoping apply
to the entire instrument with retained histories. Hence

```
epsilon_joint = O(q^-20 log^16(2q)),
additional accounting error = O(q^-16 log^16(2q)),          (9)
```

using the parent's positive accounting range E_max=O(q^4). The full decoder
pullback is positive and unital; native-carrier complements retain their
failure charge. These estimates compare decoded quantum outputs **and the
joint majority-decoded classical archive**, with all named aborts. Protection
diagnostics remain available but are not required to have the noiseless
distribution. No noisy finite member is declared exactly lossless.

Equations (1)--(9) establish existence of a complete noisy-control realization
of the same source experiments. They do not need a reliable run-time classical
service. Constants may be large; faster local control and growing archive
inventory are paid cofinal resources, not fixed-device predictions.

## 7. Two boundaries that cannot be hidden in an ideal decoder

**Unmaintained whole-history failure.** For a central bit exposed to symmetric
flips at rate mu during h, `pi(h)=(1-exp(-2 mu h))/2`. A fixed odd n repetition
word left alone for a positive time has a strictly positive majority error
`sum_{j>n/2} binom(n,j) pi^j (1-pi)^(n-j)`. If a growing independent collection
of records is stored that long, the probability that its entire decoded
history is correct tends to zero. Even fixed redundancy on perfectly protected
quantum outputs is insufficient. For n=5 the exact error polynomial is
`10 pi^3 - 15 pi^4 + 6 pi^5`; the verifier independently enumerates all masks.
Growing redundancy and/or active maintenance must enter the model. This
construction explicitly supplies both, including their costs.

**Last bare-bit obstruction.** Suppose a final physical bit must remain exposed
for h after its last corrective operation. For the two possible intended
labels, any preceding preparation followed by this binary symmetric channel
has worst-case error at least pi(h). Its output distinguishability contracts
by `1-2pi(h)`; data processing prevents a preceding encoder from increasing
it beyond one. Thus (9) cannot be transferred to a bare bit without bounding
its final exposure. With h=a/c its
floor is Theta(q^-1), although the coded public record has (9). This is why
majority is an explicitly specified read statistic and is never smuggled in
as a perfect physical feedback gate. No architecture can promise a perfectly
reliable final naked bit under this noise model.
Faster terminal service can reduce the floor if its extra control cost is
charged; this is not an absolute prohibition on asymptotically reliable bare
outputs under a different terminal schedule.

The operational gain for OPH is therefore precise: the exhibited microscopic
process can supply the same causal reads, fermionic dynamics and clock outputs
with noisy decisions **and live noisy public records**, without adding a new
reliable-control axiom. Its existing source/CCG and clock-capability assumptions
remain visible rather than being relabeled as conclusions of A1--A3.
