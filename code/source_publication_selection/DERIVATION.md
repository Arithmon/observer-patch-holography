# Eventual publication, joint demand selection, population and radius

## 1. Constructing the continuation potential

Fix a finite move alphabet A, a deterministic native response checkpoint
transition T(s,a), a normalized reference p(s,a)>0, and the indicator f(s)
that the retained record determines the declared public meaning. Publication
persists because records are retained. Checkpoints are maps over all
admissible preparations, not the realized unknown payload.

Let P act by (Pg)(s)=sum_a p(s,a)g(T(s,a)). Define H_n=P^n f and

    h(s) = sup_n H_n(s).

Each H_n is an exhaustive word sum, proved in the checkpoint package.
Persistence gives f<=Pf, so 0<=H_n<=H_(n+1)<=1. Finite summation commutes
with the increasing limit, hence Ph=h. If g>=f and Pg=g, induction gives
H_n<=g and h<=g. Thus h is the least harmonic majorant of f, constructed
from native continuations rather than supplied as a success oracle. Further,
h(s)>0 iff some H_n(s)>0, equivalently some finite word publishes. At a
published state h=1; at a state with no publishing continuation h=0.

On connected finite native support, with the receiver sampled initially,
the preceding meaning theorem identifies the latter feasibility condition:

    h(s)>0 iff equal retained records and equal current states imply
                   equal required public meanings.

The `HasRoot` condition is retained: equal records determine the current
receiver scalar. Every native advance preserves it. The alphabet must cover
the graph edges; permission of every alphabet move is needed separately to
conclude that selected histories remain on that support.

`SourcePublicationMass` and `SourcePublicationNative` kernel-check these
claims. They do not decide publication for arbitrary classical meanings.

## 2. A consistent law without a supplied deadline

For h(s)>0 define the finite cylinder and transition laws

    Q_s[w] = R_s[w] h(T_w s)/h(s),
    K(s,a) = p(s,a) h(T(s,a))/h(s).

Harmonicity proves normalization and sum_a Q_s[wa]=Q_s[w], including
zero-mass prefixes. Their product telescopes to Q_s[w]. Positive cylinders
are exactly the prefixes from which the required meaning remains viable.
This is the Doob h-transform, not a newly invented probability construction;
see Levin--Peres--Wilmer, *Markov Chains and Mixing Times*, second edition,
[section 17.6](https://pages.uoregon.edu/dlevin/MARKOV/markovmixing.pdf).
The source-specific step is constructing h from the native publication
predicate and proving its viability characterization.

For a fixed prefix w of length n, the actual marginal of the finite deadline
n+k information projection is

    Q_(n+k)[w] = R_s[w] H_k(T_w s)/H_(n+k)(s).

This follows by summing every length-k suffix, not by renaming a proposed
kernel. Its denominator is positive for all sufficiently large k. Monotone
convergence gives Q_(n+k)[w] -> Q_s[w]. By contrast, an unpublished but
recoverable prefix has zero weight at its own deadline and positive limiting
weight. Independent finite-deadline optimization need not be consistent.
The cylinder and actual-marginal statements are kernel checked in
`SourcePublicationLaw`.

Here is the analytic infinite-path argument. The finite alphabet word tree
with its consistent reference cylinders defines a probability on A^N. The
event E of publication at some finite time is the increasing union of finite
cylinder events E_n, so R_s(E)=h(s). No measurability of a hypothetical
uncountable checkpoint space is needed: only the countable reachable word
tree is used. Conditioning on E gives exactly the displayed Q cylinders.
Consequently Q(E)=1. Among all path measures V with V(E)=1, Q is the unique
relative-entropy minimizer against R. For V absolutely continuous with
respect to R,

    D(V || R) = D(V || Q) - log h(s).

The identity holds with infinite entropy as well. Singular V has infinite
cost; Q has finite cost -log h(s). Nonnegativity and the equality case of
relative entropy prove uniqueness. This extends the finite checkpoint
variational problem on an explicitly declared infinite history algebra;
it is not an identification with canonical A3's ontic-state problem.

## 3. Quantitative native completion and effective approximation

Assume N>=2 ports and every reference move probability is at least eta>0,
uniformly over reachable checkpoints. Necessarily eta<=1. The connected
calibration theorem gives a publishing word from every viable checkpoint
with at most B=N(N-1)/2 means. Pad a shorter word to B with permitted moves;
retention preserves publication. Thus H_B(s)>=delta=eta^B whenever h(s)>0.
This is an explicit, possibly very poor, bound. When N=1, a viable sampled
checkpoint is published at time zero and needs no block argument.

At a published state H_B=h=f=1. At an unpublished positive-potential state,
H_B>=delta>=delta*h. At a dead state H_B=h=f=0. Therefore

    h-H_B <= (1-delta)(h-f).

Apply the positive linear operator P^n and use P^n h=h. Iteration gives

    0 <= h(s)-H_(kB)(s) <= (1-delta)^k (h(s)-f(s)).

`SourcePublicationTail` checks the block-to-geometric inference and the
published/unpublished/dead case split. The native derivation of its uniform
block hypothesis above is analytic and uses the existing kernel-checked
connected completion theorem; it is not asserted as a new fully composed
Lean theorem.

For a viable initial checkpoint the selected first-publication time obeys
Q(T>kB)=1-H_(kB)/h <= (1-delta)^k and E_Q T <= B/delta. Publication is almost
sure, not bounded on every infinite word. A selected process can make an
arbitrarily long sequence of harmless moves. Time-prefix consistency does
not imply compatibility under changes of regulator, population or meaning.

For rational linear preparations and meanings, finite H_n is computable by
exact response elimination and finite word enumeration. For k>=1, writing
theta=(1-delta)^k gives a certified interval

    H_(kB)(s) <= h(s) <= min(1, H_(kB)(s)/(1-theta)).

These bounds shrink geometrically. Kernel intervals follow by bounding the
child numerator and positive root denominator separately. An exact sampler
can refine these rational intervals and independent dyadic random bits
until a uniform draw lies strictly between two cumulative-probability
boundaries. It terminates almost surely because there are finitely many
boundaries and the intervals converge. This is an analytic computability
argument, not a claim of efficient planning, a physical random source, or
a bounded-work sampler. The controls have a closed-form potential and need
no such expensive infinite-state approximation.

## 4. Exact four-port native law

Use the isolated path 0--1--2--3, with independent real deviations x,y at
ports 0,1 and known zero deviations elsewhere. Retain the initial scalar
at receiver 3 and every subsequent receiver sample. Edge i averages i,i+1.
Let its constant probability be p_i>0, with sum one. Require both x and y.

Before edge 0 or 1 occurs, edge 2 only averages zeros. If edge 0 occurs first,
x-y is erased with no record and publication is impossible. If edge 1 occurs
first, publication remains viable after every continuation. To prove the
second assertion, apply the response to the difference preparation (1,-1).
Immediately after edge 1 it is

    left(a)=(-2a,a,a,0), a=-1/2.

Before the first edge 2, edge 0 maps left(a) to
right(a)=(-a/2,-a/2,a,0), and edge 1 maps right(a) to left(a/4).
The other two operations are idempotent. The parameter never vanishes.
The first edge 2 records a/2, a nonzero response to (1,-1). If that record
row is (u,v), then u!=v. Every native mean preserves the total x+y; this
total and ux+vy jointly distinguish x and y. Thus the retained row plus
the current state remain viable forever. Before that sample the nonzero
difference and the conserved total give the same joint-state conclusion.
These native algebra identities and the contrast/total injectivity are
kernel checked in `SourcePublicationPath`; induction over all continuations
and the resulting probability classification are analytic.

On this safe class, uniform native completion implies reference publication
with probability one: in each block there is probability at least delta
to publish, and no trajectory can enter the dead class. The initial potential
therefore solves h=p_1+p_2 h and equals p_1/(p_0+p_1). The derived initial
kernel is

    K=(0, p_0+p_1, p_2).

After entry into the safe class K=p; while edge 2 repeats, the displayed
initial kernel continues to apply. It is not local rejection and
renormalization, whose initial probabilities are
(0,p_1/(p_1+p_2),p_2/(p_1+p_2)). For uniform p, these are respectively
(0,2/3,1/3) and (0,1/2,1/2).

Requiring only x+y on the same preparations, records and native transitions
gives h=1 everywhere reachable: conservation makes this meaning viable
forever, and the uniform block argument again gives publication. Its
selected law is the original reference. Neither experiment instantiates
the faithful full A1 response and endogenous A2 holonomy.

## 5. Permanent native rank on arbitrary paths

The preceding law extends far beyond four ports. Write the two preparation
response columns as x_i,y_i in path order and put
Delta_ij=x_i*y_j-x_j*y_i for i<j. Initially all ordered minors are
nonnegative. An adjacent mean at e,e+1 sends a minor containing both rows
to zero, leaves a minor containing neither unchanged, and averages two
ordered minors when exactly one row is included. No order reversal occurs
because the averaged indices are adjacent. All ordered minors therefore
remain nonnegative.

For a fixed nonadjacent pair i,j, one adjacent mean cannot include both
rows. Its new minor is at least half its old minor. On the path 0,...,d,
d>=2, the first safe edge 1 produces Delta_02=1/2. After any T further
native means, Delta_02>=2^(-(T+1))>0. The current state alone therefore
retains both preparation directions forever; no safety test is needed
after entry. `SourcePublicationPathInvariant` kernel-checks preservation
of all ordered two-row minors and this bound for every finite word of
adjacent means. Restricting the natural-number indices to edges below d
is exactly the finite native path, extended by zeros off that path.

Consequently the exact eventual-publication potential on **every** such
path is p_1/(p_0+p_1), independently of its length. Before entry, the first
edge among 0,1 decides safety; all other edges only average zeros. After
entry, the uniform block argument gives reference publication with
probability one. The selected first probabilities are K_0=0,
K_1=p_0+p_1 and K_i=p_i for i>=2; after entry the reference is unchanged.
For d=1 the receiver's initial sample is y and h=1, so the formula does
not apply. This is a retained-sample edge case, not a rank exception to hide.

There is a stronger analytic extension to an arbitrary number k of source
coordinates. Place them at an independent set I of path vertices (no two
adjacent), with calibrated zeros at the remaining vertices. The initial
N-by-k response matrix has every maximal ordered minor nonnegative and
the minor on I equal to one. Each adjacent pair-mean matrix is totally
nonnegative: it is a block diagonal matrix with identity blocks and the
2-by-2 all-one-half block, whose minors are nonnegative. Cauchy--Binet
preserves nonnegativity of all maximal response minors. The coefficient
of the previous I minor in its update is one if neither endpoint lies in
I and one half if exactly one does; both cannot lie in I. Other terms are
nonnegative. Hence, after **any** T native means,

    det(response on I) >= 2^(-T) > 0.

Thus the original k-coordinate preparation remains injectively encoded
in the current native state for every finite word. The topology-derived
completion is available from every prefix. Under any reference with a
uniform positive move bound, eventual recovery of all k coordinates has
probability one. The publication constraint is redundant in this class:
h=1 and the unique constrained path law is the original reference. This
removes the guard and any successful-word restriction for the stated
prepared path model. It does not require an IID selected law; a uniformly
positive state-dependent reference works as well.

For two initially adjacent interior sources at i,i+1 and a receiver outside
those two sites, the first active edge among i-1,i,i+1 decides the outcome.
The middle edge erases the difference; either existing outer edge creates
a positive nonadjacent minor, which then persists. With constant positive
edge probabilities, h is the sum of the existing outer-edge probabilities
divided by the sum over those active edges. Missing boundary edges are
omitted. For initially nonadjacent sources h=1 by the fixed-minor argument.

The independent path controls include two through five source directions,
3 through 9 ports, separated and adjacent interior placements, and every
word through each specified horizon. The producer evolves exterior-power
coefficients directly. The verifier executes independent integer-scaled
scalar preparations and recomputes every maximal minor by fraction-free
elimination. The general k-coordinate Cauchy--Binet proof is analytic, not
a Lean claim or a consequence inferred solely from these finite tests.

These statements assume exact arithmetic. The positive determinant can
decay exponentially, so they supply no fixed physical noise tolerance or
record lifetime. They also depend on the path's complete adjacent alphabet:
allowing a direct mean between two separated unknown sites can erase their
difference immediately. A complete captured support with extra seams does
not inherit this invariant without a further argument. Nor does inevitable
eventual recovery select a finite metric neighbourhood: in this model it
recovers the full declared preparation, regardless of which local read menu
one hoped to infer. This identifies the exact scope of the positive closure.

## 6. Central selection on a common demand algebra

One cannot infer a choice of public demand by comparing A3 optimizers on
unrelated spaces. Instead, explicitly include a finite demand label m in
one joint history algebra. Let pi_m>0 be its normalized prior and E_m the
event that its declared meaning is eventually published. Use the joint
reference pi_m R and impose the single constraint E={(m,w):w in E_m}.
When Z=sum_m pi_m h_m>0, the preceding entropy identity proves a unique
joint optimizer. Its demand marginal is

    Q(m)=pi_m h_m/Z.

Every recoverable demand has positive mass. Consequently this faithful
common-space completion cannot select one demand exclusively while another
recoverable demand remains admissible. For equal priors on the two-record
and aggregate demands in the uniform four-port example, the probabilities
are 1/3 and 2/3, not exclusive selection of the stronger meaning.
`SourcePublicationMenu` checks the finite joint marginal, its limiting
formula, positivity and the exact example. The infinite joint entropy
application is analytic. The prior, label algebra and publication constraint
are explicit extra structure; the theorem is not a non-entailment result
for every possible full A1--A3 realization.

This pinpoints what source selection must prove: the complete A1-generated
observable/constraint grammar must identify the required demand or exclude
the competing demands on the actual common algebra. Merely adding a label
and applying strict convexity or unique information projection does not do
that. The earlier KL-support theorem supplies the same warning for a fixed
finite model; the present construction includes native continuation masses
and a compatible infinite-time limit.

In the separated-source path class the conclusion is sharper. Every public
meaning that is a function of the full preparation is recovered almost
surely. For any finite menu of such demands, all h_m=1, so Q(m)=pi_m
exactly. This includes competing metric-neighbourhood subsets of one
declared set of source records: publication changes none of their prior
weights. `automatic_demands_unchanged` checks the algebraic conclusion.
The radius is not merely nonunique in the continuum limit; this particular
native publication mechanism exerts no selection among these demands at
any finite population. It remains a statement about this explicit common
model, not an impossibility theorem for all full-axiom source extractions.

## 7. Radius universality and alternative populations

For q=F_n, retain the actual golden sites
L*fract(phi*b), b in {0,...,q-1}^3. For any fixed 0<alpha<1 set

    a_q=L*q^(-alpha), Delta_q=a_q/c, h_q=2*sqrt(3)*L/q.

Then a_q->0 and h_q/a_q=2*sqrt(3)*q^(alpha-1)->0. Applying the proved
source-net cone theorem gives eventual exact agreement of the generated
all-neighbour read order with every strictly timelike/spacelike limiting
pair. The inner/outer cone bounds, and their shrinking shell, give the
existing interval-volume and strict-pair limits under their stated
boundary-null/contained-diamond hypotheses. For actual golden event counts,
`SourceRadiusSelection.golden_count_family` instantiates the existing
general event-count theorem, including the clipped endpoint correction,
with this duration. It does not rebrand the older arbitrary-duration
theorem as a new theorem. The positive radius, vanishing radius and mesh
ratio identities are also kernel checked. Null-cone boundary agreement
and physical clock selection are not conclusions.

Population choice is likewise not fixed by these limit tests. On the same
cube take u_q(b)=L*b/q, with cells product_i[L*b_i/q,L*(b_i+1)/q).
They partition [0,L)^3, each has volume L^3/q^3, and every point of its
cell is within sqrt(3)*L/q of its representative. Thus the proved general
partition counting/quadrature theorem, with boundary-null sets, applies
directly. Tensor with the same time cells; their endpoint discrepancy is
at most L^3*Delta_q. The same covering cone proof applies to this grid.
These are ordinary cube distances, not a periodic torus. Taking any
integer q->infinity for the grid shows that Fibonacci indexing is not
necessary for these limit properties. This grid instantiation is analytic.

The finite menus nevertheless differ. Exact ordered read counts, including
self, for exponents (1/4,1/2,3/4) are:

| q | population | 1/4 | 1/2 | 3/4 |
|---|---|---:|---:|---:|
| 5 | golden | 7919 | 2785 | 1153 |
| 5 | grid | 8671 | 4087 | 1685 |
| 13 | golden | 1473795 | 294191 | 54669 |
| 13 | grid | 1516723 | 321839 | 50653 |
| 21 | golden | 20376393 | 2911845 | 327151 |
| 21 | grid | 20141913 | 3118389 | 277255 |

For D=distance^2/L^2, the exact golden comparisons are q*D^2<=1,
q*D<=1, and q^3*D^2<=1, respectively. The producer uses integer arithmetic
in Q(sqrt(5)); the verifier uses rational enclosures with certified signs
and exact rational boundary handling. Both count complete populations
using one-dimensional multiplicities. No expanded q^6 pair tapes are
stored, and a separate direct small-population enumeration tests the tensor
counting reduction. These alternative geometric realizations are not full
A1--A3 countermodels.

## 8. Where the square root comes from, and what it does not select

The source history is explicit about introducing this parameter. Commit
`f80d0f90e` introduced the instruction to apply the complete-neighbour law
with a_q=L/sqrt(q) in `SOURCE_NET_CAUSAL_LIMIT.tex`; commit `f61ef0c0`
subsequently defined that radius in `GoldenSourceCausalLimit`. Those lines
choose the scale. They are not an extraction from A1--A3. The following is
a mathematical explanation of the scale, not a claim about the author's
unrecorded motivation.

There is a stronger explanation than an arbitrarily proposed balancing
objective. `SourceNetVolumeError.weighted_alexandrov_error` already bounds
the error of the **actual generated interval**, including the temporal
quadrature error. Normalize c=1 and a containing observation window to
T<=1. If both the covering and assignment error are bounded by h<=1,
0<a<=1, and the theorem's coverage, partition and buffered-domain
hypotheses hold (including 2h<a), its bound is

    E <= 4*pi*(T+a)*(T/2+h)^2*(h+T*h/a) + pi*T^3*a/2
      <= 18*pi*h + pi*(a/2+18*h/a).

`SourceRadiusRemoval.normalized_volume_bound` checks the latter inequality
for every admissible T,a,h. The first inequality is the existing generated
volume theorem; it is not supplied as new evidence. The a term is the time
quadrature cost, while h/a is the accumulated spatial approximation cost.
The radius-dependent part of this **certified upper bound** has the unique
minimizer a=6*sqrt(h), since

    a/2+18*h/a - 6*sqrt(h) = (a-6*sqrt(h))^2/(2*a).

`certified_balance` checks the minimum and its equality case. For
0<h<1/36 this minimizer also satisfies 2h<a<1, so the balancing regime is
actually admissible. With h=C/q it has exponent one half. Among fixed
power choices a=q^(-alpha), this bound decays at rate
q^(-min(alpha,1-alpha)); the largest certified exponent is 1/2. A different
window, sharper bound or weighting changes the coefficient. This unit-window
bound gives 6*sqrt(C)/sqrt(q), not 1/sqrt(q). Nor does optimizing an upper
bound prove that the true error is minimized there. Restoring a reference
length gives the dimensional geometric-mean form sqrt(L*h).

Thus a square-root scale is justified as a useful regulator balancing two
proved error mechanisms. It is neither selected by the cone/count limit
nor shown to be the entropy optimizer of the canonical A3 problem. It is
not established as a physical length attached to the source.

For dimensionless mesh epsilon=h/L and radius r=a/L, the deliberately
supplied minimax error objective max(r,epsilon/r) has the unique minimizer
r=sqrt(epsilon). Indeed max(r,epsilon/r)<=sqrt(epsilon) implies both
r<=sqrt(epsilon) and epsilon<=r*sqrt(epsilon), hence equality. The converse
is immediate. `SourceRadiusSelection.balance_unique` checks this.
For h=C*L/q this selects a=sqrt(C)*L/sqrt(q), not automatically the M1
coefficient. Weighting the two terms also changes the coefficient. Thus
a rate-balancing explanation can recover an exponent only after choosing
an objective and its normalization. Canonical A3 minimizes a specified
relative-entropy functional, not this error functional. No equivalence
between them is proved, so the balance criterion is not consumed as an axiom.

## 9. Can the radius be removed?

**The exact formula can be removed from the continuum target; spatial
causal control cannot simply be deleted from this construction.** Section 7
already proves that all fixed exponents 0<alpha<1 give the same stated
limits. The radius tends to zero; it is a regulator in those theorems,
not a surviving fundamental length. That does not imply all finite read
menus are interchangeable: their counts, timing and records differ.

Here is a direct test of deleting the radius while keeping the same layered
read prescription and vanishing layer duration. Permit every site in the
previous layer as a read. In the actual generated read order, e precedes f
exactly when its layer is smaller or e=f. `complete_walk` and
`complete_order` prove this for the existing generic `Walk`/`LayerPrec`.
Every fixed positive spatial separation can then be crossed in one layer.
As delta->0, for every fixed finite c one eventually has c*delta<distance;
`complete_read_exceeds_speed` checks this. It contradicts finite-speed
propagation without requiring any count interpretation.

The causal interval also has the wrong ordering fraction. Take two singleton
tips enclosing K full interior layers, each with N sites. All its KN+2
events lie between the tips. The only incomparable pairs are the
K*N*(N-1)/2 distinct same-layer pairs. Therefore its fraction of unordered
comparable pairs is exactly

    F = 1-K*N*(N-1)/((K*N+2)*(K*N+1)),
    1-1/K <= F <= 1.

It tends to one as K->infinity uniformly in N>=1, rather than the 1/10
four-dimensional Minkowski **interval** value. The comparison uses an
actual causal interval, not an arbitrary rectangular observation window.
The finite combinatorial count is analytic and independently reconstructed
from DAG transitive closure in the controls; the formula's uniform bound
and limit are kernel checked. Endpoint tips and N=1,K=1 edge cases are
included. This is a countercontrol to removing the locality restriction,
not a full A1--A3 countermodel.

Conversely, retaining only a fixed microscopic nearest-neighbour stencil
does not automatically fix the problem. In a unit cube grid of spacing
1/q, take radius a=1/q and duration a (c=1). The only one-layer moves are
self and the six axial neighbours. Reaching (1/2,1/2,0) from zero takes
q moves and hence time one. At time 3/4 the displacement is already strictly
Euclidean-timelike, since sqrt(1/2)<3/4, but remains unreachable. For every
q divisible by four, the discrepancy persists under refinement: the
generated cone is the L1 cone. If a=1/q^2, only self reads remain. The
controls reconstruct the Euclidean stencils by exact rational inequalities
and perform grid breadth-first search, independently of the producer's
closed formulas. They include both boundary exponents and exact equality
at radius 1/q. No inference from a finite sample to the general claim is
needed; the integer displacement argument proves it for all such q.

These failures explain why the successful family uses separated scales,
h/a->0 and a->0: enough microscopic directions are available per shrinking
time step to recover the round cone. They do **not** prove these conditions
necessary for every conceivable read architecture, or that all radius-free
architectures fail. The earlier generic `SourceNetLayeredOrder` already
accepts arbitrary read relations. A source-derived relation could replace
the ball rule, provided its actual propagation and counting limits are
proved. Local transport topology, eventual recoverability and A5 symmetry
alone have not supplied that join; the scalar action continuum and the
record-dependency DAG must not be silently identified.

The practical consequence is to separate two targets. Literal finite M1
still asks for a particular population and read menu. The continuum causal
target needs a source-derived propagation law in the appropriate limiting
class, not a proof that precisely L/sqrt(q) is fundamental. The results
remove that exact formula as a necessary intermediate claim for the latter
target. They do not establish membership of the complete A1--A3 source in
that class, and do not silently revise the original M1 exit criterion.

## 10. Full-basis exit and remaining derivation

The compatible temporal law and its native feasibility are derived for a
stated meaning; their arbitrary deadline and completion-probability oracle
have been removed. The joint-demand theorem and explicit radius/population
families show precisely why two proposed routes to central selection do
not finish it. They do not establish that the full canonical axioms cannot
select M1 by some other mechanism.

A derivation of literal finite M1 still requires an extraction from the complete A1 source
algebras and faithful response, compatible with endogenous A2 transport,
that fixes the golden q^3 population, required previous-version record map
and radius L/sqrt(q). It must identify the actual A3 common feasible space,
reference and complete temporal constraint grammar, and prove consistency
under regulator refinement and immutable-version interpretation. Availability
of a checkpoint is insufficient as a proof of compulsory execution; no
full-model stationary-execution counterexample is asserted here. The
physical layer clock is a separate M2 attachment. PR-52 is unconsumed.

The evidence therefore supports a substantial temporal closure and a
specific selection-obstruction theorem for the proposed common-space
mechanism. It does not meet the original full-M1 exit criterion.
