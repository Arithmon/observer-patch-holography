# What a recurrent repair instrument must change or supply

This extends [the source-to-read derivation](DERIVATION.md). It distinguishes
three objects: the registered numerical operations, quantum-channel
extensions of their *ideal* diagonal means, and the canonical A1--A3 source.
The first is executed. The second is constructed and compared. Neither is
silently identified with the third.

The companion exporter calls the same `_visible_pair_mean` and
`_terminal_complex_lift` functions as the registered driver. Extracting
these two functions preserves its threshold, arithmetic and zero-phase
convention. All seven previously captured primitive/observer cases and
earlier quantum controls were reproduced unchanged. The additional
preparations below enter at the post-unitary operation boundary; they are
interventions, not purported Gaussian-seed outputs.

## 1. Two exact obstructions to promoting the current rule

Every quantum channel is linear on density matrices, hence preserves
convex mixtures. This is a necessary condition before complete positivity
or any particular physical interpretation is considered.

**Native no-op band.** Let epsilon be the exact binary64 value of `1e-15`
and delta=2^-51. The two-slot diagonal states

    p_j = (1/2+j*delta, 1/2-j*delta), j=0,1,2

are positive, exactly representable and normalized. The native primitive
leaves p_0 and p_1 unchanged because 2*delta<=epsilon, but sends p_2 to p_0
because 4*delta>epsilon. Since p_1=(p_0+p_2)/2, the first output coordinate
violates affinity by exactly delta. There is no quantum channel agreeing
with this thresholded diagonal rule on all three states. This statement
already holds for the real-valued threshold rule; it is not an objection
based only on arbitrary floating-point rounding. The exported native
results are checked bit-for-bit. The general inequality is kernel checked.

The unthresholded ideal pair mean does admit quantum extensions, proved
below. A numerical tolerance may approximate such a map; it must not be
mistaken for its exact channel identity.

**Phase-preserving terminal lift.** A separate obstruction survives removal
of the no-op band. Interpret normalized carrier amplitudes through the
explicit block density encoding

    rho(psi) = (1/N) direct_sum_c |psi_c><psi_c|.

Other encodings are not asserted equivalent. Keep every carrier except
zero in the positive uniform amplitude state. At carrier zero compare the
two equal-weight ensembles

    {|0>, |1>}   and   {(|0>+|1>)/sqrt(2), (|0>-|1>)/sqrt(2)}.

Their average input density matrices agree exactly. After a complete
repair sweep every partner of carrier zero has initial intensity 1/12.
For the first two inputs its repaired intensities at ports 0,1 are
(13/24,1/24) and (1/24,13/24). The registered lift uses phase zero at a
zero amplitude, so both output coherences rho_(0,1) are
sqrt(13)/(24N). For the second ensemble the repaired intensities are
(7/24,7/24), and the coherences are respectively +7/(24N) and -7/(24N).
Its average output coherence is zero. Equal input ensembles therefore give
unequal average outputs, with gap sqrt(13)/(24N)>1/(8N).

This proves that the native phase-preserving projection, interpreted on
these density states, has no affine quantum-channel extension. It cannot
be repaired just by declaring a missing CPTP receipt true. A simulator
retaining a classical preparation label may implement the nonlinear
operation on that larger input; then that label is extra state and the map
does not descend to the stated density matrix. The native control checks
the four selected coherences against the analytic values with an explicit
numerical enclosure. Lean checks the ensemble-linearity and positive-gap
reductions; the density-ensemble identification is analytic.

There is also a normalization issue if each carrier vector is treated as
an independently normalized state. One seam update changes the local sum
from 1 to 1+(x_b-x_a)/2, while preserving the sum of the two carriers.
The global density encoding above accommodates changing carrier weights;
it does not prove preservation of each carrier's unit trace. Normalizing
each terminal carrier separately would introduce a different nonlinear
rule. The local trace identity is kernel checked.

## 2. Valid extensions exist, but the scalar repairs do not select one

Use the explicit global mode space H=C^(12N), with trace-one density
matrices. This is a mathematical extension space, not the A1 accessible
algebra. For a seam e={a,b}, let S_e exchange those two modes and fix the
others. Define

    Phi_e(rho) = (rho + S_e rho S_e*)/2.

The Kraus operators I/sqrt(2), S_e/sqrt(2) prove complete positivity and
sum K* K=I proves trace preservation. Lean proves these properties using
amplified positivity of permutation conjugation and convex closure of CPTP
maps. On the diagonal, Phi_e is precisely the ideal pair mean. Since
S_e^2=I, Phi_e is idempotent on the full matrix state. Distinct seam
transpositions commute, so their twirls commute as well. Their complete
composition E has exactly the original ideal terminal diagonal map R.
Neither the numerical threshold nor the phase-preserving pure lift is used.

Let D(rho)=sum_i P_i rho P_i be coordinate dephasing. Its coordinate Kraus
operators P_i make D CPTP; this is also kernel checked. Then

    E and D composed with E

are two valid channels with identical diagonals for *every* density input.
Their rules transform covariantly when the mode labels and seam matching
are relabeled together. They differ on coherences. Every convex mixture
lambda*E+(1-lambda)*D E, 0<=lambda<=1, is another CPTP extension with those
same immediate diagonal reads. Thus complete scalar repair data alone do
not determine a unique quantum instrument, even within this explicit
covariant family. No A3 objective has been inserted to select lambda.

**A subsequent native read separates the extensions.** In carrier zero
take the phase state (|0>+i|b>)/sqrt(2), where b is adjacent to zero;
let every other carrier block be I_12/12, with block weights 1/N. This is
a normalized positive density state. E distributes its single local
coherence among four ordered endpoint choices; the original local
coherence is retained with factor 1/4. D E erases it. Apply the actual
block-local U(t) next and read port zero. The difference is one eighth of
the opposite-phase probability gap, additionally divided by N:

    |p_E-p_DE| >= (2t-400t^2)/(8N).

At t=2^-17 and N=4 the lower bound is 16359/34359738368. The control
executes both channels and the native unitary. A separate checker
reconstructs the full sparse density outputs using exact Gaussian
rationals, checks every entry and the next probability, and invokes the
previous rational unitary enclosure. Choi tests also check positivity and
partial-trace normalization, with transpose as a negative control.

## 3. The native port effects cannot themselves be A1's central records

This obstruction is stronger than a mismatch of terminology. Interpret
the twelve intensity probes as P_i=|i><i| in C^12 and require the
accessible complex algebra to contain them and their native reversible
time translates. It contains Q=U* P_0 U. The already proved nonvanishing
of the twelve entries U_(0,i) gives

    P_i Q P_j = conjugate(U_(0,i))*U_(0,j) |i><j|,

with nonzero coefficient for every i,j. Closure under products and complex
scalars therefore yields all 144 matrix units: the generated algebra is
M_12(C). Its center consists only of scalar multiples of I, so none of
these twelve nontrivial coordinate projections is central. Lean checks
the coordinate-sandwich and commutator identities; the generation and
center argument is analytic.

Consequently one cannot repair the A1 central-record clause just by adding
the missing phases to this same twelve-dimensional coordinate algebra.
A separate classical record register can coexist with quantum state:
the outcome map rho -> direct_sum_i P_i rho P_i is a trace-preserving
instrument into a classical/quantum algebra, whose block labels are
central. These labels are stored outcomes, not the coherent instantaneous
coordinate projections under a unitary that mixes them. The repository's
Lueders/Kraus results support that construction. It still needs an actual
source binding, response architecture and archive/continuation semantics.
This is a necessary interface distinction, not a full A1 no-go theorem.

## 4. A derived recurrence cone, without promoting influence to a read

To test whether a valid repair extension can overcome native saturation,
compose the explicit measured cycle D E Ad_U repeatedly. After its first
cycle the state is diagonal, so subsequent cycles act on populations by

    T = R B_local,       B_(p,q)=|U_(p,q)|^2.

For this diagonal recurrence B is positive in every entry and doubly
stochastic; R averages the matching partners. The exporter executes this
transition on all basis populations for four cycles, without a target
geometry or acceptance filter. This is a declared instrument control;
the registered driver still performs its unitary stage before its repairs.

Let G be the actual carrier graph and let m(c,p) be the carrier partnered
to slot (c,p). Exact positivity gives, for k>=1,

    (T^k)_((c,p),(d,q)) > 0
      iff min(dist_G(c,d), dist_G(m(c,p),d)) <= k-1.               (1)

For k=1 the two permitted input carriers are c and m(c,p); every q is
available because B is strictly positive. For the induction step, summing
over the twelve slots of each permitted carrier takes the union of its
own and its neighbors' previous dependency sets. Nonnegative terms cannot
cancel, proving (1) at k+1. Thus the union of a carrier's twelve outputs
has exactly the radius-k *graph influence ball*. The independent verifier
compares Boolean composition, this distance formula and the support of
the executed matrix powers, retaining every slot and every population.

The tested carrier-graph diameters are 1,2,2 at N=4,8,16. Every carrier
can influence every carrier within at most two measured cycles; every
slot can influence every slot within three. This radius is a number of
specified cycles on a supplied graph, not an emergent Euclidean distance
or M1's fixed read radius. Nor is nonzero influence exact decoding.

## 5. Exact normalized information loss survives that propagation

Let V be the 11N-dimensional tangent space of per-carrier normalized
preparations. For the matching mean R, an element of ker R has opposite
values at the two endpoints of each seam. Write their oriented values as
w_e. The normalization constraints are exactly partial*w=0, where partial
is the N by 6N oriented carrier incidence matrix, including parallel seams.
The graph is connected, so rank(partial)=N-1: row sum gives the upper bound,
and a spanning tree gives N-1 independent columns by deleting leaves.
Consequently

    dim(ker R intersect V)=5N+1,   rank(R restricted to V)=6N-1.  (2)

This also closes the normalized terminal-rank calculation for the original
ideal native repairs. It is distinct from their ambient rank 6N.

For the recurrence, put eta=20t. The norm bound
||U-I||<=exp(10t)-1<=20t gives B_ii>=(1-eta)^2. As B is doubly stochastic,
||B-I||_infinity<=4*eta<1. The Neumann series makes B invertible with
||B^-1||_infinity<2. It preserves each carrier-sum constraint, so B_local
is an automorphism of V. Therefore (2) holds for T=R B_local as well.
All 5N+1 dimensions of this kernel are erased at the first cycle and
cannot be recovered from *any* later global population trajectory.

There is a constructive positive witness. The registered port-zero and
port-one matchings form a Hamiltonian carrier cycle. Orient a circulation
around it: put +1 at each outgoing slot, -1 at its paired incoming slot,
and zero elsewhere. The vector z has zero sum at every carrier and Rz=0.
Set v=B_local^-1 z. Then

    x_plus=1/12+v/48,     x_minus=1/12-v/48

are distinct, per-carrier normalized and at least 1/24 in every component,
yet T x_plus=T x_minus. Every later population record is identical. For
trace-one global densities divide both ledgers by N. The exporter solves
the inverse and executes both transitions. The independent checker proves
the exact incidence rank and circulation and checks the numerical inverse,
normalization and complete first-step outputs to residual 2^-40. The exact
erasure statement follows analytically from invertibility; a small floating
residual is not represented as exact binary64 equality.

Retaining initial records or suitable instrument outcomes may prevent this
loss. The theorem concerns unconditional population state and its later
records. It does not rule out an enlarged instrument with an archive,
nor imply that the repair's unrecorded Kraus branch is already such an archive.

## 6. Carrier counts alone do not give a closed coarse dynamics

Let C sum the twelve slot populations at each carrier. A Markov matrix K
with C T=K C exists exactly when the total transition into each output
carrier is independent of the input port within an input carrier. Necessity
follows by applying the identity to basis populations; the same columns
construct K for sufficiency. This is strong lumpability, not a fitted
coarse trajectory.

The source violates it. Choose d=m(0,0) and another port q with m(0,q)!=d;
the alternating matchings supply one. A port-zero input sends at least
(1-eta)^2/2 into carrier d. A port-q input sends at most 11*eta^2/2 there,
because the diagonal q term cannot reach d. Their coarse probabilities
differ by at least

    (1-2*eta-10*eta^2)/2 = 536706947/1073741824 > 0.

The executed transition probabilities and exact incidence choice are
independently checked at every retained population. Thus discarding port
state does not give a closed carrier Markov model of this extension. This
is an obstruction to that specific coarsening, not to every possible A1
refinement map.

## 7. Constructive retained-record completions

There are two different ways to close the information-loss gap, and they
must not be confused with each other or with existing native records.

### Exact deterministic population recovery with a minimal linear archive

Keep the ideal deterministic mean y=Rx. Orient each seam e=(a,b), and
write d_e=(x_a-x_b)/2. Then x_a=y_a+d_e and x_b=y_b-d_e. For per-carrier
normalized inputs the differences obey

    partial*d = 1-C*y,

where C sums slots at each carrier. Choose any spanning tree of the
actual carrier graph. Retain only the differences on non-tree seams.
There are 6N-(N-1)=5N+1 such coordinates. Their contribution can be
subtracted from the right side; the remaining tree differences are
uniquely recovered by deleting leaves. The last carrier equation follows
from conservation of global mass. This reconstructs every normalized
input exactly from the deterministic output and the retained differences.
It uses no radius, embedding, success-conditioned history or inverse fit.

This coordinate count is minimal for a **linear** archive. Any linear
record map H that makes (Rx,Hx) injective must be injective on ker R
intersect V, whose dimension is 5N+1. Fewer real output coordinates cannot
suffice by rank-nullity. The argument makes no claim about arbitrary
discontinuous encodings into a single real, finite bit precision, or the
quantum capacity of this classical ledger. Global access to the means and
chord records is required; this is not a local M1 decoder. Nor is the
spanning-tree choice a canonical source-selected gauge.

The producer executes rational recovery at all three populations, with
21,41,81 retained differences. The independent verifier uses Gaussian
elimination rather than the producer's leaf decoder and checks every
input, mean, chord and recovered-state digest. Additional tests recover
arbitrary positive rational preparations. The all-input result is the
incidence argument, not the finite test sample. For T=RB_local the same
archive applied to B_local*x, followed by B_local^-1, reconstructs the
initial normalized input. It preserves the exact ideal deterministic
means. The native threshold, floating rounding and terminal phase lift
remain separate operations; these controls do not silently replace them.

### Quantum recovery with retained random-unitary branches

For the coherent twirl E, let s range over all 2^(6N) binary seam words,
let S_s be the product of the indicated transpositions, and set p_s=2^(-6N).
The flagged channel with a classical register is

    J(rho) = direct_sum_s p_s S_s rho S_s*.

Its rectangular Kraus maps are sqrt(p_s)|s> tensor S_s; their adjoint
products sum to I, so J is CPTP. Tracing out the flags gives E, since
expanding the product of the individual two-branch mixtures gives exactly
this sum. Conditional inverse followed by forgetting the flag gives

    sum_s S_s* [p_s S_s rho S_s*] S_s = rho

for every density matrix, including all coherences. Lean proves the
inverse identity for every permutation and its weighted finite sum. The
verifier reconstructs executed branch permutations independently for the
zero, all-one, every single-seam and alternating words. A complete
matrix-unit test also checks all words for a two-seam system. These finite
controls do not replace the arbitrary-word proof.

This enlarged instrument retains 6N uniform independent bits per sweep,
plus its quantum state. Its branch diagonals are the original or swapped
values; **only their expectation** equals the deterministic pair means.
It is therefore not a lossless implementation of deterministic means
on each branch. Interleaving known unitary propagation still allows
recovery by inverting the entire recorded word in reverse chronological
order; discarding or dephasing the quantum state generally does not.
The flag marginal is uniform and independent of the input. Flags alone
therefore reveal no preparation data, even though flag plus quantum state
admits recovery. They cannot be renamed the required public record.

These constructions close the algebraic question of whether suitable
additional records can prevent erasure, and quantify two concrete costs.
They neither show that the current source keeps those records nor select
which enlarged interface A1--A3 requires.

## Consequence for the M1 path

These results resolve specific uncertainties rather than postulating the
desired read law: the current threshold and lift are not exact channels
on the stated density inputs; valid ideal-repair extensions exist but are
operationally distinguishable; the twelve coherent probes cannot also be
the central records; a recurrent extension has a derived graph influence
cone but irreversibly erases a known normalized information family; and
carrier-only state fails an exact closure test. A minimal linear difference
archive recovers every normalized deterministic input, while retained
random-unitary flags recover the entire quantum input through a different
channel with the same mean repair.

A source capable of M1 must specify the quantum/classical record interface,
the instrument and retained outcomes, and the actual admissible/refinement
family. The complete A1 response and endogenous A2 data still have to be
realized before A3 can select within that family. Nothing here supplies
those premises by renaming the extension controls, and no number of these
simulation cycles establishes the missing canonical selection.
