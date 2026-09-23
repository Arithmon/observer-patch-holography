# Does the full three-axiom basis force a record-transmitting seam?

This construction tests the central M1 selection question, rather than
proposing another successful routing policy. The distinction to test is
between the geometric overlap nerve and the information carried by its seam
algebras. A degree-one map of nerves does not assert that either endpoint's
record algebra embeds in a seam.

## 1. Model and proof obligations

At level n let S_n be the normalized edge-midpoint subdivision of the
icosahedral sphere. Put one carrier at each vertex. Its open-star cover has
nerve S_n. The bridge b_n is the identity and the oriented fundamental cycle
maps with degree one. The number of carriers is 10 * 4^n + 2. Each carrier
has its own twelve-port icosahedral boundary; that local boundary is not the
global triangulation.

At carrier i set

    A_i = direct_sum_(p=1..12) M_6(C),
    R_i = Z(A_i) = C^12,
    e_(i,p) = identity in block p and zero in the other blocks.

These e_(i,p) are primitive central projections, not the noncentral rank-one
coordinate projectors of the previous amplitude simulator. A carrier record
is the outcome p of its nondemolition central measurement. The internal
six-dimensional response is observable through local matrix effects.
The port interface is I_(i,p)=M_6(C), with the unital readout
pi_(i,p)(a)=a_p. A mismatch score is the trace distance of a populated
block's normalized internal state from I_6/6, and zero on an empty block.
The requested-block reset repairs this mismatch in one step.

For distinct carriers declare every nonempty pair and triple seam algebra
to be C. The inclusion of a seam into a carrier is lambda -> lambda * 1.
The unital data restriction is the normalized faithful trace

    E_i(a) = (1/72) sum_p Tr(a_p).

Triple restrictions are identity maps of C. Thus the restriction cocycles
commute. The common record domain on a seam is its scalar algebra; a local
port outcome is not declared to be a shared record of another carrier.
Local state restriction along the scalar inclusion is normalization, and
does not impose equality of the twelve private port probabilities.

The patch category contains these carrier patches, their nonempty
intersections, and their presentation charts. Closure under unions is not
an A1 requirement. In particular, a fictitious global observer with direct
access to all private records is not inserted. Local algebras of distinct
carriers embed in different tensor slots for the purpose of checking the
declared disjoint-region locality law. This algebraic embedding does not
assume a global state in defining the feasible local family.

The construction checks, separately:

* the complete response and same-response overlap holonomy;
* operational readback, record, feedback, prediction and checkpoint types;
* the full refinement squares, including record and response transport;
* a complete A2 constraint grammar and the exact A3 state problem;
* noncommunication for every program admitted by the actual source grammar.

The algebraic response checks do not by themselves prove that this is a
full-axiom model. The remaining construction and its interpretation must be
checked below. Adding a commuting record toy to an unrelated response
representation would not suffice.

## 2. A source response with twelve directions

Choose the usual exact oriented icosahedral vertex frame v_p over Q(sqrt(5))
as a calculation chart; let sigma be its Galois conjugation and r^2 its
common squared norm. The choice of chart is quotiented by proper incidence
relabeling. Define the following two 3-by-3 skew-adjoint blocks:

    D_p = diag(hat(v_p/2) + i Q_p, hat(sigma(v_p)/2)),
    Q_p = v_p v_p^T / 2 - r^2 I / 6 + I / 12.

Here hat(v) w = v cross w. The primitive response operations are generated
by these matrices on the same accessible six-dimensional internal state.
They are not an independent spectator representation. Ordered second-order
response coefficients reconstruct [D_p,D_q].

The even antipodal fields span all six imaginary symmetric matrices in the
first block. The odd fields span the two real skew blocks independently.
Consequently D has real rank twelve and its image is exactly
u(3) direct_sum so(3). It is commutator closed. Skew adjointness and
injectivity imply that -Tr(D(v)D(v)) is strictly positive for v != 0.
The accessible internal algebra is the full M_6, including effects between
the two coordinate triples. Process probes [D,E_jj]_(ij) and
[D,E_i6]_(i6) reconstruct the off-diagonal entries and all diagonal
differences. Matrix units are linear combinations of positive preparation
and effect probes, so this is ordinary linear process tomography. These
probes have response rank twelve. Restricting observations to M_3 direct_sum
M_3 would lower the rank to eleven: the relative block phase would become
invisible. The executable deletion control checks that loss. The response
therefore has no unobservable scalar direction on its stated accessible
algebra, and its center is not an operationally silent spectator.
Every accepted infinitesimal reversible response of the declared primitive
grammar is in this image: matrix multiplication and inversion generate its
connected matrix group, and no additional reversible generators are admitted.

For a proper incidence automorphism a, reconstruct the real rotation R_a
from R_a v_p = v_(a p). Then

    U_a = diag(R_a, sigma(R_a)),
    U_a D_p U_a* = D_(a p).

Both real skew blocks belong to im D. Every U_a is therefore in the same
connected response group, without using a centralizer spectator. A finite
source implementation uses Cayley factors. Rotations of order three have
no eigenvalue -1, so K=(R-I)(R+I)^(-1) is skew and
(I+K)(I-K)^(-1)=R. Every proper icosahedral rotation is a product of at most
two order-three rotations, or the identity. The source solves each pair of
Cayley generators in the original twelve D_p and executes those factors.
The verifier checks the resulting matrix, its port action and all group
products. This supplies response operations, not just a claimed list of
implementers. Recharting permutes the central port blocks and conjugates
their internal matrices by this same U_a.

The calculation uses an explicitly constructed source, not the registered
physical simulator. Selecting a mathematical witness is legitimate for an
existence or non-entailment proof; it is not evidence that the physical
simulator or the axioms select this particular witness.

## 3. Overlap groupoid and the operational grammar

Distinguish two types of overlap. Distinct carrier patches share only their
declared scalar seam. Presentation charts of one carrier describe the same
accessible algebra A_i. Their reversible overlap groupoid is the action
groupoid of the proper carrier automorphisms on that unlabeled carrier.
It has the complete carrier chart as object and the sixty proper actions as
closed arrows. The port functor is the identity on this group. Its response
functor sends a to the actually executed U_a from section 2. Thus it is
surjective on port actions, composes even before projectivization, and has
the required same-response endogenous implementers.

This is a local chart holonomy, not a channel carrying an unknown private
state from carrier i to a different carrier j. The distinction matters:
neither the A1 statement about unital seam restrictions nor the A2 formulas
for Pi and U require a faithful record-transmitting representation of a
distinct-carrier seam. If that additional requirement is intended, it must
be stated and tested; the present construction would then be excluded.
The local chart overlaps are not additional carriers or vertices of the
federation nerve. They supply the recharting morphisms in the patch category.

On a chart change the algebra map is

    (alpha_a x)_(a p) = U_a x_p U_a*.

It acts on the central ports as required. E_i alpha_a = E_i by cyclicity of
trace and port permutation. Thus its scalar overlap square commutes.
Composition gives all longer chart paths and all higher-overlap squares.
This does not ask an inner automorphism of A_i to permute its center:
the central permutation is a change of chart, whereas the response
implementer in PU(H_i) is supplied by the same internal source D_i.
Those are exactly the two distinct maps in the A2 statement.

The primitive operational grammar has the following types.

1. Reversible response: source-independent port controls in im D, their
   ordered compositions and inverses, and the executed chart paths above.
2. Record: nondemolition measurement of the twelve central projections.
   A central outcome p remains the same under internal responses.
3. Internal readback: local finite matrix effects, with their Born
   probabilities. A binary effect E can be read by the destructive
   instrument T_y(rho)=Tr(E_y rho) e_y/6, where E_0=E, E_1=1-E and
   two requested central ports store the result. This overwrites the old
   central record; it does not add an unaccounted classical register.
   Calibrated preparations and independent trials supply process tomography.
   Transport state, effect and output-port request on a chart change.
4. Repair/feedback: a requested port q and the local record p trigger an
   internal trace-reset when p=q, otherwise no change. The reset replaces
   the internal state by I_6/6. This is a bounded, unital CPTP operation.
   Its dependence on equality of ports is presentation invariant. It can
   repair a declared internal preparation mismatch; no settling requirement
   is asserted for arbitrary mismatch functions.
5. Prediction: from a supplied local preparation description and retained
   local conditioning data, compute the next instrument's probabilities.
   This is not an oracle measuring an arbitrary unknown density exactly.
   Checkpoint continuation retains the physical
   density in place and the central outcome; it does not clone an unknown
   quantum state into a classical archive. The carrier has no independent
   transcript memory. The verifier retains an execution log; the stronger
   information bound below also allows arbitrary receiver-local memory.
6. Between distinct carriers the shared instrument is identity on C.
   There is no record-to-seam encoder. A source-independent tick schedule
   may interleave local operations, but an operation at one carrier neither
   reads nor modifies another carrier's state or controller.

Reversible probe controls and irreversible feedback controls are separate
types. A record-dependent family of arbitrary blockwise unitaries is not
silently added to the reversible grammar: it would introduce additional
response directions. The stipulated feedback operation has a nontrivial
kernel on the internal matrices of its requested central block. Subsequent
linear channels cannot undo that loss for all input states. Internal
readback uses the specified measure-and-prepare instruments, also with
nontrivial kernels, not arbitrary unitary backactions. Hence including
these local read/prepare instruments adds no omitted
infinitesimal reversible response. The chart permutations normalize im D.
This proves response completeness for the whole declared grammar, rather
than checking twelve selected directions of a larger admitted tangent.

Each single carrier is a connected subfederation with the complete
readback, record, feedback, prediction and continuation interfaces. The
canonical axioms do not require every connected subfederation to have an
arbitrary remote-record service. Availability of the listed local
operations is nontrivial even though inter-carrier record transport is absent.

The source can be extended by arbitrary local controllers, local retained
transcripts and independent random tapes without restoring communication.
We use that larger class in the information bound. We do not assert that
every such extension preserves the stipulated twelve-dimensional response
grammar; the bound is stronger than needed for the actual source.

## 4. Refinement, including the information maps

The standard midpoint subdivision retains old vertices. Map each new edge
midpoint to one endpoint, using a total order in a calculation presentation.
Transport the resulting map with any relabeling; it is declared coarsening
data, not a canonical endpoint selector from an unlabeled edge. For an
ordered triangle a<b<c, only the child at a has a nondegenerate image, with
the original orientation. The other three children collapse. Consequently
the fundamental two-chain pushes forward exactly, at every level.
Compose one-step maps to define all longer coarsenings. With N_n=S_n and
b_n=id, the bridge square is an identity.

There is no hidden assertion of a Euclidean read radius. For the actual
spherical realization, let alpha_0=acos(1/sqrt(5)) be the seed edge angle.
If a spherical triangle has angular diameter at most alpha_0, normalization
of endpoint sums bounds each midpoint-to-midpoint chord by the opposite
chord divided by 2 cos(alpha_0/2). The same bound holds for corner half
edges. Indeed, for nonzero vectors x,y,

    ||x/||x|| - y/||y|||| <= ||x-y|| / min(||x||,||y||),

as follows by expanding ||x-y||^2. The endpoint sums have norm at least
2 cos(alpha_0/2). The contraction constant is less than one. Spherical
triangles lie in an open hemisphere and radial interpolation gives the
orientation-preserving homeomorphism onto the sphere. Chord meshes tend
to zero, hence so do angular meshes. A point and its coarsened image lie
in the same parent triangle; their distance is bounded by that mesh.
This supplies the controlled support estimate, not merely finite topology.

On local state families the one-step coarse state at i is the equal mixture
of the fine states over the nonempty fiber F_i of the vertex coarsening:

    rho_i = (1/|F_i|) sum_(j in F_i) rho_j.

Its Heisenberg map sends a to the constant family (a)_(j in F_i) in the
direct sum of descendant algebras. It is a unital star homomorphism.
The normalized mixture is its state pullback; it is not asserted to be a
star homomorphism in the reverse direction or a quantum copying channel.
This direct sum describes the tagged local experiment, not the joint
physical algebra of all descendants. On that joint tensor algebra the
dual of random local extraction is a -> sum_j a_j/|F_i|, a unital
completely positive map which is generally not multiplicative. Neither
map transports a private state across a fine seam.
The corresponding record interface forgets the tag in the disjoint union
of descendant records (j,p) -> p. A tagged refinement observer chooses its
descendant once and retains that tag in its continuation data. It does not
combine private records from different descendants into a new local record.

All local response matrices, port incidences and instrument formulas are
identical in the descendant charts. For every local linear instrument T_y,

    (1/|F_i|) sum_j T_y(rho_j) = T_y((1/|F_i|) sum_j rho_j).

This equality is for the unnormalized outcome states; dividing by the same
outcome probability gives the conditional state when that probability is
positive. Thus read probabilities, updates, feedback and continued words
commute, not just their unconditional averages. The direct-sum response
intertwines the constant-family inclusion. Scalar restrictions also commute
with averaging. Coarse local words lift as the same program independently
run on every descendant. Coarsening selects one tag and keeps it through
checkpoints. All descendant operations remain charged. Updating one child
and then averaging it with unchanged siblings would fail this square and
is not the lift used here. Each child branches on its own local outcomes;
the selected coarse outcome is not broadcast to its siblings.
A local word cannot acquire a new read from
another fiber. These supply the routing, record and checkpoint squares.

The fibers of a one-step map are connected stars: every member is either
the retained parent vertex or a midpoint incident to it. A coarse edge
(i,k) lifts to i--midpoint(i,k)--k. For any chosen descendant tags j in F_i
and l in F_k, the path j--i--midpoint(i,k)--k--l therefore has at most four
fine seam traversals after deleting repeated endpoints. Each traversal
carries the unique scalar datum, so their composition is the coarse scalar
interface. All traversals are retained operations, not hidden nonlocal
copies. Composing these paths gives every further refinement. The finite
tower checker tests all choices of descendant tags for every coarse edge.

At several levels, the mixture weights are the products of the one-step
fiber weights, not an incorrectly presumed uniform weight on the entire
composite fiber. All weights are positive rational functions of the declared
coarsening data and are transported under presentation change. Constant
tracial states push forward to constant tracial states at every level.
No optimizer-pushforward theorem for general A3 specifications is assumed.

## 5. Complete constraint grammar and actual A3 selection

Specify the compatible local family first. At the base chart of each
carrier it consists of twelve positive semidefinite 6-by-6 blocks with
total trace one. Other chart states are their alpha-transforms. Seam and
triple states are the unique normalized state on C. These are all the
constraints; there is no assumed global density in this definition.

For each carrier, take the full real Hermitian matrix basis in each of
the twelve blocks: six diagonal coordinates and two coordinates for each
off-diagonal pair. There are 12*36=432 real observable coordinates.
They determine every local state. Every local effect expectation factors
linearly through this basis. Every specified instrument acts linearly on
these coordinates; conditional expectations are obtained by division by
the explicitly computed outcome probability. Every finite history is
therefore obtained by composition of these same coordinate maps. Naturality
for any accepted finite history follows by induction from the generator
instrument squares. On distinct-carrier overlaps, the entire accepted
shared algebra is C, so the only state condition is normalization.
Pair, triple, chart and refinement morphisms generate the data category.
Their checked squares consequently cover every A2-visible constraint.

This is a completeness proof for an explicit full model. It is not a claim
that every A1 model has scalar seams, this instrument grammar, or these
state spaces. In particular, it does not exclude additional constraints
in the actual physical source.

After normalization there are 431 independent state coordinates: 71
diagonals, omitting one calculation-chart anchor, and the real and imaginary
parts of the 180 off-diagonal entries within blocks. The last diagonal is
one minus the other diagonal sum. Omitting any one of these 431 coordinates
hides a perturbation I_72/72 +/- T/288. For a missing diagonal take
T=E_aa-E_anchor; for a missing real or imaginary coordinate take its
Hermitian matrix unit. Each T is traceless with operator norm one, so both
states are positive. `grammar.py` checks all these omission witnesses.

The A3 cover is the set of all carrier base patches, taken modulo their
presentation charts. Set w_i=1/N and

    tau_i = direct_sum_(p=1..12) I_6/72.

This rule is exact, faithful, quotient-visible and presentation natural.
It is specified by this source's uniform carrier sampling and calibrated
full matrix trace; it is not inferred uniquely from invariance alone.
Other models may have other exact natural reference rules. There is one
fixed rule here, not an unresolved menu inside the variational problem.
The cover determines the entire compatible family: all other states are
its specified restrictions or chart transforms. We use ordinary matrix
trace in the displayed 72-dimensional faithful representation. For any
local density with eigenvalues lambda_k,

    D(rho_i || tau_i) = sum_k lambda_k log(72 lambda_k) >= 0,

with equality exactly when all lambda_k=1/72. Convexity of x log x proves
this inequality and its equality case. A density with these equal
eigenvalues is exactly I_72/72. Thus the unique minimizer of the actual
weighted local Umegaki objective is rho_i=tau_i for every carrier.
This proves existence, uniqueness, faithfulness and the selected port
weights 1/12. The weights normalize to one and their behavior on refinement
is explicit; the optimizer commutes here because every mixture of these
identical local references equals the same reference.

For an executable source realization, initialize the carriers independently
in their tracial states and run the declared local instruments. This is
one global preparation extending the already defined compatible family.
It is not invoked to establish A2 agreement or incorrectly claimed to be
selected among all global extensions by a local-cover A3 objective.
Independence is an admissible source preparation in this constructed model;
the axioms need only admit one such model for a non-entailment result.
Nonselective record measurements, internal responses, and the feedback
trace-reset preserve tau_i. Conditioning on a local central outcome p
adds its actually observed record constraint; the same entropy calculation
gives I_6/6 in that block. No remote-success constraint is inserted.

## 6. All-program read obstruction and population

Consider two source preparations differing only in the central port at
carrier s. Couple every exterior initial state and randomness tape. Each
primitive exterior operation has the same inputs in the two executions:
its own state, its own retained records and its own random tape. Induction
therefore gives identical exterior local transcripts at every finite
prefix. Equality for every prefix gives equality of the whole infinite
transcript. Any decoder or stopping rule that is a function of that
transcript gives the same answer in both worlds. It cannot read both
different source values correctly. This proves an obstruction for every
admitted routing, including unbounded local stopping, in this particular
source grammar. Ten of the fifteen Lean declarations prove the functional induction,
decoder/stopping consequences, a cross-copy deletion control and the
population arithmetic. The quantum/algebra-to-local-kernel identification
and the full axiom-model construction are analytic, not formalized in Lean.

Under the actual tracial product preparation the source central label is
uniform on twelve values and independent of the receiver's entire local
transcript. Hence any receiver decoder has success probability at most
1/12 for that label; abstention counts as failure. A binary source
intervention with equal priors gives the corresponding bound 1/2. No
global state hash, source-dependent scheduler, transport receipt, public
measurement outcome or timing side channel is inserted into the receiver's
data. Adding any such carrier of the source value changes the source.

The maximal source-record read graph consists of the self edges. Local
records are readable and persistent; no off-diagonal initial record is
readable. The graph stays so under the stated refinement. It cannot be the
M1 family, which has a q^3 record population and its nontrivial complete
metric-neighbor read menu. The population obstruction is independently
exact: 10*4^n+2 is never an integer cube. At n=0 it is 12; at every positive
level it is 2 modulo 4, whereas a cube is 0,1 or 3 modulo 4. Increasing
the support resolution does not turn this native family into M1.
The population count alone is not an obstruction to encoding q^3 logical
sites on a differently sized host. The information theorem supplies the
stronger test: an operational projection using only receiver-accessible
data cannot restore a remote independent record intervention. A global
decoder pooling all private logs would add the missing communication
interface and is not such a projection.

The conclusion is non-entailment by the stated full A1--A3 clauses, not
impossibility of M1 in every model and not impossibility on the distinct
W12 scalar-mean source. Earlier routing constructions remain valid under
their supplied communication and controller interfaces. The repair-law
RFC's A1-R/A2-R are not canonical axioms and are not silently included here.

The missing structural element is a connection between carrier records
and inter-carrier operational channels, with a rule selecting which
records must be transmitted. Degree-one spherical coverage and internal
response holonomy do not supply that connection. A3 selects the maximally
random local state in this model and cannot turn a scalar seam into a
record channel. A positive M1 derivation therefore needs a substantive
source theorem or additional structural requirement excluding this model;
calling the existing support conditions a derivation would be incorrect.

## 7. Finite evidence and what it checks

The retained response packet contains the twelve generators, 420 complex
process-tomography coordinates, 66 ordered mixed responses and their
coefficients in the reconstructed source basis, twenty executed Cayley
factors and sixty closed paths. Ordered responses are produced by composing
four primitive response jets in the bivariate truncated polynomial ring;
the verifier separately computes the commutator and checks its expansion.
The source defines a geometric primitive law, with no supplied Lie-type,
particle, coupling or empirical target. The independent checker uses a
closed per-port formula, rather than the producer's antipodal-band sum.
It checks 720 covariance identities and all 3,600 path products.

Four support stages have 12,42,162,642 carriers. Exact chain pushforwards,
vertex links, connectedness, all lifted scalar paths, nonuniform composite
fiber weights and a selective-instrument posterior control are checked.
The actual 144 two-carrier central-record intervention executions retain
6,912 local events. They use the internal preparation

    sigma_delta = (I_6 + delta (E_11-E_22))/6, delta=1/4,

so feedback has an observable effect. Reset changes delta to zero. This
positive nontracial preparation is a declared local intervention; the
unconstrained A3 state is the tracial state proved above. The record
statistics, all consumed local deviation values, feedback, predictions
and scalar seam outputs are independently replayed. The remote transcript
is unchanged for all twelve source labels, while all twelve source-local
transcripts differ. This is exact software execution, not a laboratory
quantum experiment or a sample of a purported physical source.

The entropy minimum, full grammar factorization, all-level geometric mesh
argument and axiom-model interpretation are analytic proofs in this file.
The verifier does not certify those proofs by comparing stored True flags.
The Lean audit certifies the stated all-program reductions and arithmetic,
not an encoding of the entire canonical axiom registry.

## 8. Clause audit and the missing transport requirement

The construction uses the formal clauses at the pinned axiom revision.
In particular it does not identify the local recharting groupoid with paths
that transport unknown records between distinct carrier algebras. That
identification would strengthen the current clauses. Naturality transports
state and effect together; it does not require every feasible state to be
fixed by all rechartings. The unique unconstrained A3 optimizer is fixed.

| Clause | Construction | Deletion or failure control |
| --- | --- | --- |
| A1.1--2: finite category and isotone net | Carrier open stars, scalar intersections, local presentation groupoids; scalar inclusions and disjoint tensor slots | A missing intersection loses its common normalization interface; identifying distinct slots loses the declared disjoint commutation |
| A1.3--4: records, interfaces, repair, continuation | Central C^12, block readouts, destructive internal instruments, requested-block reset, in-place checkpoint | Internal rank-one projectors are not central ports; deleting reset loses the specified mismatch repair |
| A1.5: oriented carrier boundary | Icosahedral template and complete proper actions | Missing incidence breaks directed-edge or proper-action coverage |
| A1.6: seams, triples, operational observer | C seams, identity triple restrictions, complete single-carrier observer | Deleting a seam loses that typed handshake and its declared nerve faces, not an unclaimed remote read |
| A1.7--8: support and refinement | All-level midpoint tower, fundamental-chain pushforward, mesh contraction, tagged instrument pullback and scalar routes | A constant bridge fails degree one; omitted fine traversals fail replay; unchanged prior weights fail selective conditioning |
| A1.9: presentation quotient | Transport every diagram and the declared endpoint selector | Holding numeric endpoint labels fixed under relabeling is not an admissible equivalence |
| A1.10: full response | Twelve observable skew generators, closure and actual same-response path factors | Block-diagonal observables give rank eleven; adding arbitrary controlled blockwise unitaries changes the grammar |
| A2: all accepted diagrams and holonomy | Instrument squares extend to words; the same D implements the complete local chart groupoid | An independent spectator cannot supply these covariance identities; external relabeling alone has no executed response path |
| A3: complete grammar and objective | Full state coordinates, 431 omission witnesses, all-carrier injective cover, fixed trace and counting rule | An omitted coordinate or carrier hides a feasible state perturbation; tilted block weights change the specified reference |

The geometry performs support-address and refinement operations. Its scalar
seams are not required by the current clauses to transmit a private record.
An independent calculation chart for D is allowed: only the intertwined
port response, probabilities, incidence, ranks and counts are protected
outputs. No measured target, particle label or fitted coupling is input.

A precise additional **transport requirement** is that an admitted
source-to-receiver channel on the central record simplex have a stochastic
left inverse, with its output actually accessible to the receiver. For a
finite classical channel P(y|p), this is equivalent to pairwise disjoint
supports of its output laws. If two distinct p values can both produce y,
no decoder of y is certainly correct on both interventions. Conversely,
disjoint supports define a deterministic decoder and its composition with
P is the identity. With twelve inputs and twelve outputs this forces P
to be a permutation matrix; more outputs allow redundant disjoint supports.
The Lean support theorem proves the decoder criterion. `channels.py`
checks exact controls, including an invertible noisy matrix whose inverse
has negative entries and is not a stochastic decoder.

This condition is stronger than an algebra inclusion, nonzero influence or
matrix rank. It supplies transmission, not source retention, compulsory
execution, population or demand selection. Composing retained-version
channels along a finite path supplies the transport interface of the
existing compiler when its other storage and control premises hold.

It would therefore be incorrect to advertise record transmission alone as
the single remaining sufficient assumption. A1--A3 leave the communication
service and realized demand/production law unspecified. The companion
`RECORD_GLUING.md` states one stronger proposed dynamical law, RG, at the
primitive process level and proves its sufficiency for the operational
causal/count limit. It does not put a Fibonacci population, read radius or
limiting count into the assumption. It does explicitly add scalable,
response-natural transport, cofinal classical workspace, graded complete
composition and native implementation. Those are substantive requirements,
not consequences of the decoder test. No unique spontaneous microscopic
population or compulsory observer program is claimed. The negative model
alone is not the positive exit, and RG has not been adopted as a core axiom.
