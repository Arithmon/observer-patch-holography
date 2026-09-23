# One additional source law: response-natural classical record transport

## Target and quantifiers

The target here is the operational causal/count limit. It is not uniqueness
of Fibonacci cutoffs, a microscopic preferred radius, spontaneous execution
of one observer's measurement program, or identification of all auxiliary
repair events with spacetime volume. M1 is compared as a regulator of this
operational limit, with its finite prescriptions kept explicit.

The scalar-seam construction rules out deriving that limit from the present
three axioms alone: a source satisfying the stated data-access clauses can
have no inter-carrier record channel. This file investigates one additional
**dynamical source law**, RG. It is a proposal, not an adopted fourth core
axiom and not a property already established for the simulator.

RG is deliberately stated at the primitive process level. It mentions no
Fibonacci number, cubical population, mesoscopic radius, Lorentz cone,
count-volume limit, successful M1 history, or entropy-conditioned routing.
Its content is stronger than saying that one seam sometimes transmits a
record. The quantifiers over independent preparations, simultaneous finite
tasks, all dependencies and refinement are part of the proposed law, not
unstated extra assumptions to be discharged in another issue.

## The source object already available

The complete-response/endogenous-transport theorem gives the compact algebra
su(3) + su(2) + u(1) at a complete reference carrier. Let E be its unique
three-dimensional simple ideal,
with the positive trace form inherited from the faithful source response.
The adjoint action of its connected response group on E is SO(E). Thus E
is an intrinsic three-dimensional Euclidean vector space, and the response
acts transitively on its unit sphere. This is a mathematical response
space; identifying translation in E with record transport is the new law.
It must not be silently inherited from the gauge Lie-type theorem.
Addresses form the affine space modeled on E; choosing its origin is a
presentation choice. RG's translation and scaling action supplies a flat
homogeneous operational arena for the observer experiment. A1--A2 alone
do not identify all carrier ideals with a single global flat space, or
identify this space with the spherical screen.

Classical record registers use finite products of the already available
central port alphabet. Any finite digital payload can be encoded in these
registers. There is no claim that unknown quantum states can be copied.

## Proposed RG law

The source's central-record process admits a complete, intervention-faithful,
refinement-compatible realization of the following **generated process
theory**, with a source-observable nonnegative duration grade:

* Objects are finite independent classical registers at addresses in E.
  Finite tensor products, local preparation, finite deterministic classical
  maps, nondestructive read/copy, discard and retained versions are realized
  at some finite source cutoff. Register size is finite for each experiment;
  it is not bounded uniformly over all cutoffs or experiments.
  Local operations have zero ideal displacement and zero ideal duration;
  their actual finite source implementations and all of their events are
  retained. They may instead have a cofinally vanishing positive overhead,
  handled by the margin construction below.
* There is one nontrivial retained-record flight of displacement e and
  duration tau, with ||e||=1 and 0<tau<infinity. Translating, applying the
  connected source response, and positively scaling this primitive gives
  flights (lambda R e, lambda tau), lambda>0, R in SO(E). These are actual
  inter-carrier record processes, not rechartings of a spectator or
  changes to coordinate labels alone. Waiting preserves the record.
* Processes compose sequentially and as finite tensor products. Durations
  add on a dependency path; independent branches use their common elapsed
  time, not a serial software replay clock. The source realization retains
  every primitive/intermediate event and its resource cost. The complete
  dependency graph, including controller, timing and metadata channels,
  factors through these local operations, flights and waits. Conversely
  every finite generated process has an admitted source realization,
  preserving all independent record interventions and immutable versions.
  These realizations commute with the declared source refinements.

Here **intervention-faithful** includes insertion of a local intervention
at any designated intermediate record output, before publication. It is
not limited to changing initial preparations. The native intervention acts
only on that record at its source; all other existing record registers and
cached copies are held fixed. Lowering commutes with these inserted local
operations and subsequent reads. An immutable version is committed once
after its intervention, rather than silently overwritten later. Native
ancestry uses the complete dependency graph of the actual execution,
including every control and metadata channel.

This is one complete transport law, not a claim that its individual clauses
follow from A1--A3. In particular, it supplies cofinal classical workspace,
parallel composability, record protection and the missing link from the
internal response to record motion. Dropping any of them is a weaker
proposal and must not be called RG. The definition is also not the assertion
that M1 holds: it constrains an arbitrary local experiment before any
population, coarse sampling or read menu is chosen. It is not asserted to
be atomic, uniquely necessary, or the weakest possible extension. Counting
these substantive clauses as one law does not reduce their physical burden.
Its irreducible new content includes homogeneity and isotropic, scalable
*operational* transport, not merely the existence of a lossless channel.

For the approximate implementation version, every finite generated process
and every epsilon>0 have a finite native realization with exact classical
values and interventions, identical logical endpoint addresses, no faster
native dependency, and each logical event's time in [t,t+epsilon], where t
is its ideal scheduled time. This is a uniform statement about the whole
finite experiment, not a fixed error per primitive or per layer. One
sufficient implementation certificate is nonnegative primitive overheads
whose sum on every complete dependency path is at most epsilon: earliest
execution of a finite acyclic schedule then has the stated bound, by
induction using addition and maximum. Include waits from the initial
clock to each ideal event so that its ideal earliest time is t. The compiler
must provide the actual lowering and resource ledger. An unimplemented
promise or an annotated simulator time is not evidence that a source obeys
RG. This compiler property is explicitly part of the additional law.

Translations and scales act on the family of finite record experiments,
possibly changing cutoff and routing. They are not extra infinitesimal
reversible port generators on one finite carrier. Classical copy and
discard are instruments on finite products of central registers. An
implementation that adds a continuous spectator response or omits part of
its public tangent would fail A1, even if it realizes the reference graph.
Existence of a generated reference process theory is not, by itself, a
construction of a complete A1--A3 native model satisfying RG.

## Deriving the communication geometry

Put c=1/tau. A rotated/scaled flight has displacement lambda R e and
duration lambda tau, so its displacement norm equals c times its duration.
A wait or local operation has zero displacement. On every source dependency
path, the triangle inequality therefore gives

    ||x_final-x_initial|| <= c (t_final-t_initial).

Completeness of the primitive grammar makes this an upper bound for every
source record influence, not only for a chosen compiler. It includes
payload-dependent control and cannot be defeated by an ungraded scheduler.

Conversely, let T>=0 and ||v||<=cT. If v=0, wait. Otherwise response
transitivity supplies R e=v/||v||. Scale the one flight by lambda=||v||;
it delivers an exact retained copy in time tau||v||<=T. Wait for the
remaining duration. Thus the operational record-access relation is exactly

    (t,x) <= (s,y) iff s-t >= 0 and ||y-x|| <= c(s-t).

The Lorentz cone is a conclusion of this primitive law; no operational
boost covariance has been assumed. The law does not prove Lorentz
covariance of matter dynamics. With vanishing positive implementation
overhead, all strictly timelike accesses are attained and all spacelike
ones are excluded; the null boundary remains a limit statement.

Finite tensor composability now realizes every finite list of such reads
simultaneously, with separate protected version registers. This is the step
that a theorem about a single point-to-point channel would miss.

## Producing finite populations and actual read histories

Choose an orthonormal calculation frame in E and a bounded observation
window [0,L)^3. These are observer presentation and regulator choices, not
preferred physical axes or a new population law. For each integer q>0,
prepare one independent record at each address

    x_b=(L/q)b, b in {0,...,q-1}^3.

RG's finite tensor/preparation rule constructs these q^3 source registers
at a finite cutoff, rather than presupposing that a particular native
carrier census equals q^3. Their half-open cells partition the window,
have equal volume (L/q)^3 and assignment error at most sqrt(3)L/q.
Refinement q -> m q maps a child cell to its parent by integer division;
its register can retain or copy the designated parent version. All other
operations remain visible, with their own source costs.

Take delta_q>0 with delta_q->0 and q delta_q->infinity. At each layer run
the exhaustive finite observation of the preceding population that fits
within duration delta_q. The derived access test, not an input metric menu,
is ||x_b-x_a||<=c delta_q. Each requested version is copied along the
constructed flight and any necessary wait. At the deadline the receiver
commits the tuple of source versions it actually consumed. Independent
interventions survive every copy. Self reads are local retained versions.
This produces the complete layered read graph by an explicit finite program.
Choosing exhaustive tomography is the definition of this observer probe;
it is not a claim that all possible native executions perform that probe.

For positive lowering overhead, fix the entire finite observation horizon
first and lower its complete finite program with uniform event-time error
epsilon_q->0, for example epsilon_q<=delta_q/q. Native commits are allowed
to occur in [j delta_q,j delta_q+epsilon_q]; perfect synchronization is not
asserted. Values, writer identities and requested dependencies remain exact.
There is no sum of unaccounted per-layer errors: RG is quantified over the
whole program. Alternatively inner menus a_q=(1-eta_q)c delta_q, eta_q->0,
provide explicit flight slack, but do not by themselves prove exact commit
deadlines. Both versions have h_q/a_q->0. All lowering events are retained;
no histories are filtered for favorable M1 behavior.

The number of elementary process operations per layer is finite: one
retained-version fork/flight/wait chain per nonself read and one commit per
receiver. Native work is the sum of the compiler's actual generator bills,
including preparation, controller and storage events; no uniform physical
work bound follows merely from the duration grade. Auxiliary events are
not assigned the logical record's volume weight.

## Causal and count limits, and comparison with M1

Let h_q=sqrt(3)L/q. The existing covering-path theorem applies to the
derived menu with a_q=c delta_q (or its inner overhead margin):

    distance <= k(a_q-2h_q)  => k-layer reachability
       => distance <= k a_q.

Hence the actual generated read order agrees eventually with every fixed
strict timelike/spacelike comparison in the interior observation domain.
The limiting null boundary has zero four-dimensional Lebesgue measure;
finite shells generally have positive measure, tending to zero with their
width on a bounded window.
The half-open spatial cells and time cells give, for every bounded
continuity set A in the observation window,

    (L/q)^3 delta_q * number_of_distinguished_records_in_A
        -> integral_A dt d^3x.

This follows from cell assignment and dominated convergence, not from an
assumed Poisson law. The clipped endpoint error is at most L^3 delta_q.
Here is the source-to-native-order join, including auxiliary dependencies.
Let P_q be the order generated by the consumed writers and N_q the complete
native causal ancestry restricted to the distinguished commits. To prove
P_q subset N_q, consider a consumed writer u and the corresponding tuple
coordinate at the receiver v. Intervene at u's output boundary before publication,
holding every other record fixed. That tuple coordinate changes. If there
were no native u-to-v dependency path, the full set of ancestors of v
would be unchanged initially and closed under all native dependencies.
Induction on native updates would leave its complete local state, and
therefore v's read, unchanged, a contradiction. Compose this argument along
writer paths. This uses intermediate intervention faithfulness, not merely
equality of final outputs for every initial preparation.

For example, let two distinct versions a and b both initially store x+1.
A compiler that substitutes a read of b for a read of a agrees for every
initial x, but fails an intervention on a alone. Initial-input probing and
the value-only `primitive_lowering_extends` lemma cannot exclude this
substitution. RG's intermediate-intervention clause does. The Lean
`network_closed_region_eq` and `intervention_forces_native_ancestry` theorems
formalize the required noninterference reduction; the complete dependency
factorization and localized interventions are explicit hypotheses.

The path bound gives N_q subset C_c at the
actual native event positions and times. For any fixed strictly timelike
pair, the covering construction eventually puts it in P_q; for any fixed
spacelike or backwards pair, the upper cone bound eventually excludes it
from N_q. Thus both orders have the same almost-everywhere limit. This
does not assume that the native graph has no additional causal edges.

If both event positions move by at most H and both times by at most epsilon,
their causal margin changes by at most 2H+2c epsilon, by the reverse
triangle inequality. With H=O(1/q) and epsilon=epsilon_q+delta_q this tends
to zero. On a finite window the cone's null boundary and the endpoint
diamond boundaries have zero measure. Pull interval indicators back to the
equal-volume cells. They are bounded by one and converge almost everywhere
by the order sandwich; dominated convergence proves the weighted interval
limit for either order. Apply the same argument on the product of two
windows to the two endpoint-interval indicators and the strict precedence
indicator. The time diagonal is also null. This proves strict-pair counts,
not just weak convergence of positions. For an interior timelike diamond
of duration T at rest, the volume is pi*c^3*T^4/24 and the oriented-pair
measure is V^2/20. Therefore 2*strictPairs/[N(N-1)] tends to 1/10.

The repository's `SourceNetOrderLimit`, `SourceNetVolumeError`,
`SourceCausalBoundary`, `GoldenSourceCountLimit.partition_count_tendsto`
and `FlatDiamondNormalization` supply the covering, partition, null-set
and normalization ingredients. The native-order sandwich above is the
analytic join to RG; those pre-existing Lean modules do not themselves
assert that any native source obeys RG.

M1 is an alternative regulator of the same RG source. Use the already
proved golden assignments with q=F_n and h_q<=2 sqrt(3)L/q, and choose
delta_q=L/(c sqrt(q)). RG constructs the finite register preparation and
all required transmissions in the exact implementation version. Its
radius L/sqrt(q) follows from that regulator duration and the derived
access law. Uniformly small positive lowering overhead gives the same
logical graph at nearby native commit times and the same native-order
limits by the sandwich above. Inner menus are an optional alternative.

There is also exact finite clearance for the usual golden menu when q>1.
Every normalized squared golden displacement is in Z[phi]. If it equalled
1/q, irrationality of phi would make its phi coefficient zero, so 1/q
would be an integer, a contradiction. Hence all finitely many admitted
nonself flights have positive slack. This permits overhead smaller than
their minimum slack when only delivery before the ideal observation time
is at issue; it does not erase commit/controller costs. In the rational
grid reference below the clock satisfies delta_q^2=2/(2q+1) (L=c=1).
Boundary equality would give (2q+1)m=2q^2 for an integer m, impossible
because gcd(2q+1,q)=1 and 2q+1>2. Both clocks obey the required scaling.

This comparison is operational and asymptotic. Both finite programs act on
the same source record process, preserve independent interventions, and
have retained native histories and cost ledgers under RG. Their finite
read graphs need not be isomorphic. They agree on limiting causal/count
observables away from null boundaries, with the existing quantitative shell
bounds. Arbitrary unrelated field dynamics or quantum preparations do not
inherit equivalence from causal/count convergence alone.

More quantitatively, for a K-Lipschitz test function in the sum metric
|t-s|+||x-y||, assignment errors H in space and epsilon in time change its
weighted integral by at most K(H+epsilon) times the window volume, plus
the clipped time-cell term. Comparing two regulators adds their two
bounds. Indicators use the vanishing boundary tubes instead of a Lipschitz
bound. This specifies which observable comparisons are proved and avoids
calling two nonisomorphic finite graphs operationally identical.

**Sufficiency theorem.** Under A1--A2's complete-response theorem and RG,
there are cofinal, independently preparable finite central-record
experiments whose complete native ancestry on their distinguished records
converges to the round 3+1-dimensional causal relation, with the stated
volume, interval and strict-pair limits. M1's golden population and radius
are one regulator with these same limits. All preparations, retained
versions, interventions, clocks, parallel operations and auxiliary resource
bills are supplied by RG's single complete source law. No further source
selection premise is suppressed in this conclusion. A3 is compatible with
selecting states within a source that obeys RG; it is not used to derive RG
or to force a source to execute this particular observer experiment.

This is the closure of the stated operational-regulator target under one
explicit proposed law, not an unconditional derivation from A1--A3 or a
proof that a particular native simulator implements the law. If M1 instead
means a unique spontaneous microscopic census or a compulsory measurement
policy, that different statement is neither the target nor a consequence.

## Executable reference and proof custody

`code/source_record_gluing` executes preparations, immutable forks, flights,
waits, local additions, commits and checkpoints for q=2,3,5. It derives its
reads from primitive flight times and the observation clock, with no input
radius or read menu. Two layers, the baseline, every single-site initial +1
intervention and every first-layer commit +1 intervention replay 5,344,106
primitive events in 323 executions. Intermediate interventions occur before
the version's sole commit. The independent checker
reconstructs the allowed transactions, exact durations, consumed writers,
values and operation census. Every intervention is compared with the path
multiplicities of the graph extracted from those consumed writers. This is
a finite reference execution of RG's process grammar. It is not a native
implementation, a check of all real rotations/scales, or evidence that the
registered simulator satisfies the new physical law.
The reference commits an additive diagnostic, not the proof's lossless
tuple payload. Its exact input and intermediate response tests check that
each consumed value affects that diagnostic. The theorem's tuple-retaining
experiment is constructed analytically from RG's finite classical maps;
the additive diagnostic is not claimed to certify arbitrary payload maps
or native implementations.

`SourceRecordGluing.lean` proves the primitive cone equivalence, its join to
the covering-path theorem, value preservation under primitive lowering,
the localized-intervention ancestry reduction, event-count retention, q^3
cardinality and causal-margin stability.
The positive sufficiency theorem as a whole, especially axiom interpretation
and native realization, is the analytic argument above. The source law is
an ordinary stated hypothesis, not a Lean axiom disguised as a proof.

The source clock grade and one volume normalization define mathematical
units. A physical clock calibration, matter dynamics, curved backgrounds
and identification of every microscopic event with physical volume are
outside M1 and are not claimed by RG's sufficiency theorem.

## Non-circularity and deletion tests

RG makes a new physical claim about source operations. It does not follow
from internal rotation covariance; the scalar-seam source is the deletion
control. With only the six axis flights, path cost is L1 and the limiting
cone is not round. With a faster ungraded metadata channel the upper cone
bound fails. With destructive transfer, immutable-version rereads need not
exist. Without finite tensor composability, individual channels need not
realize a complete layer. Without cofinal workspace, growing populations
are not supplied. Without complete event accounting, the construction may
be merely a logical emulator with hidden causal dependencies.

These are explicit parts of the one proposed source law. Its sufficiency
is a mathematical theorem; whether the physical simulator or any proposed
microphysics satisfies it remains the test of the additional assumption.
No entropy optimizer is used to conceal these operational requirements.
