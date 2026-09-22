# Observation-generated temporal constraints

Let x belong to a real vector space of independently variable source
preparations. The preparation map P, primitive seam means A_e and local
sample functionals r_i are fixed before requesting any output. A sample
after the first k means has the response

```
L[k,i] = r_i A[e_k] ... A[e_1] P.
```

Every sample through the declared prefix is retained. The common scalar
baseline is preserved by each mean and can be subtracted locally. The
native response theorem uses precisely the existing pairAverage and run
definitions. It does not replace the word with a desired transfer matrix.

## Exact acceptance and complete linear factorization

Two preparations give the same receiver record exactly when Lx=Ly. For a
linear public meaning f, agreement of interpreted values on this fiber is
equivalent to each of the following conditions:

```
ker L is contained in ker f;
f belongs to the span of the observed response rows;
there exists c with f = c L.
```

The proof applies to every linear f, including forms absent from the finite
controls. It is a factorization theorem for this specified observer data
domain. The canonical partial policy from source_read_acceptance is applied
to the independently derived relation

```
Possible(value,samples) iff there exists x with Lx=samples and f(x)=value.
```

At every feasible sample, it publishes f(x) exactly when the factorization
holds. If it fails, there is z with Lz=0 and f(z) nonzero. The preparations
x and x+z give identical samples and different meanings. Every sound
policy therefore abstains. This necessity is stronger than crossing a cut:
a crossing in the wrong temporal order or a rank-deficient mixture can
leave the record unidentified.

The executable certificate has two alternatives. A decoder provides c and
the verifier checks cL=f. An ambiguity provides z with Lz=0 and f(z)=1.
Checking these identities proves either alternative without trusting the
producer's elimination algorithm. The finite Python certificates are
independently checked rational identities; no general Python-to-Lean
implementation equivalence is claimed.

Retaining another sample can only enlarge the available span. The earliest
prefix admitting f therefore gives a payload-independent stopping rule.
No sound rule can publish it at an earlier prefix on the unrestricted
linear domain. The positive execution witness and physical implementation
of that stopping rule are distinct obligations. For componentwise sample
errors bounded by epsilon_i, a decoder certificate gives

```
|c samples - f(x)| <= sum_i |c_i| epsilon_i.
```

This is a noise-amplification bound. An exact certificate alone does not
establish finite-precision integer identification or hardware feasibility.
Public sample rows, scale metadata, retained sample storage and decoder
arithmetic are explicit resources, external to the counted scalar means.

## The A3 consequence

Let K be a convex family of nonnegative laws on a finite history space,
with a faithful reference, positive KL weights and an attained minimum p*.
The constrained support theorem from source_read_selection proves that
every history having positive mass under any q in K has positive mass
under p*. Hence a specified f is publishable with certainty under p*
exactly when it belongs to every observation span on that feasible support.

This derives the full temporal exclusion test for the finite linear branch.
Any admitted history outside that test has an explicit invisible source
intervention. Its positive probability cannot be removed by selecting a
different faithful reference on the same feasible family. The actual A1
history grammar must exclude it, or the service must retain abstention.
A state-determining entropy cover alone is insufficient: the same-reference
cover counterexample in source_read_selection applies without modification.

The toy receipt uses all words of a three-edge chain through horizon eight.
It records first-publication counts, failed histories, expected stopped
means and samples against the original 3^H denominator. Sampling includes
the initial receiver value and every executed prefix. A failed attempt
executes the complete horizon. Read/write access for each mean is twice
the mean count. Each attempt prepares four scalar registers. These are
mathematical attempt-law costs; the external coefficient computation and
physical sampling instrument are not silently counted as scalar means.

## Native completion from connected topology

The earlier source path-tomography result supplies a recursive inverse of
one sweep with known relay values. `SourceTemporalTomography.native_endpoint`
proves that this abstract sweep is exactly the coordinate readout of
the actual pair-mean word on any simple path. Its inverse receives the
destination sample and the previously calibrated relay values; the remote
source value is reconstructed, not passed to the decoder.

Start with the receiver's own current value. Suppose a connected set K of
ports is calibrated from its retained samples. Choose an uncalibrated
vertex u adjacent to K and a simple path from u through K to the receiver.
Sweep the canonical means along this path. If its relay pre-values are
b1,...,bd and its final receiver value is y, the decoder reverses

```
y <- 2*y - bd; ...; y <- 2*y - b1.
```

The result is u's pre-value. Every relay's changed value is then computable
from the same observations, so the calibrated region grows by u. No
unintroduced vertex was previously touched, so its pre-value is its
initial value. The Lean proof works directly with observation-fiber
agreement and establishes both original-state identification and continued
calibration. It also proves that a phase cannot erase an initial distinction
from the combination of the old record and the current complete state.

`SourceTemporalConnected` derives a complete plan from finite graph
connectivity, without supplying the plan as a hypothesis. A proper rooted
calibrated region has an edge to its complement: otherwise it is closed
under every graph walk, contradicting connectivity. Attach that edge to an
existing calibrated path and recurse. The complement cardinal strictly
decreases. Exactly N-1 phases finish coverage. Each simple phase has at most
|K| means, giving at most 1+...+(N-1)=N(N-1)/2 means and N receiver samples.
Native lowering, exact sample count and the triangular work inequality are
kernel checked for every typed plan, not only the chosen spanning tree.

The captured graph has 15,360 ports and 46,050 edges. The producer's
breadth-first tree and an independent layer-based reconstruction agree on
every parent and complete path commitment. The verifier checks the actual
plan conditions, including all supported hops and complete coverage. The
plan has 15,359 phases, 425,979 means, 851,958 primitive reads and the same
number of primitive writes. Initial preparation costs 15,360 writes and
there are 15,360 receiver samples. The maximum phase length is 51. The
corresponding endpoint-noise gain 2^51 assumes exact relay calibration;
it is not a bound for the recursively propagated errors of the full plan.
Neither physical precision nor a full payload execution of this large plan
is asserted. Its paths and tree are regenerated, not stored as raw tapes.

Eight five-port histories separately verify execution, including all five
basis states and signed nonbasis preparations. The independent decoder is
given only the paths and receiver samples. It reconstructs every original
port value and the final evolved state. The finite Python topology checks
instantiate the premises of the general Lean construction analytically;
there is no claim of a generated Lean certificate for all captured paths.

## The maximal continuation constraint

Let L be a retained prefix record and M its current complete physical state
response. A meaning remains recoverable exactly when equal (L,M) data imply
equal meanings. Necessity follows because no common deterministic future
can distinguish equal current states with equal past records. Sufficiency
has a native construction: execute the connected completion protocol,
recover Mx from receiver samples, and interpret f from (Lx,Mx).

For linear preparations and meanings the exact condition is

```
ker [L; M] <= ker f,
equivalently f belongs to the row span of [L; M].
```

Thus a temporal rule preserving the ability to finish a required read must
exclude precisely the prefixes failing this test. This is a maximal
continuation constraint, not an assumed set of successful words. It is
derived from the actual prefix maps and retained record. It tests all
possible payloads without reading their realized target values. Connected
topology supplies a bounded native completion of every admitted prefix.
Applying the constrained KL support theorem to these joint rows gives the
kernel-checked `selected_native_completion_iff`: the constructed protocol
recovers f with certainty under the selected law exactly when every history
admitted by any feasible law passes the joint-span test.

The word census separates current publication, pending but recoverable
records, and permanent loss. After one chain event, one of the three words
has already erased either individual record, while the other two preserve
them for completion; none yet delivers either to the receiver. Their sum
remains recoverable after all three words. Every category uses the original
attempt denominator and includes the work of failed attempts.

Preservation does not imply spontaneous progress. Repeating an idle or
already-agreeing seam can preserve a record indefinitely without delivering
it. The construction supplies a terminating service protocol and its full
cost, not a theorem that canonical A3 selects or executes that controller.
No measure conditioned on completed words is substituted for the original
attempt law. The accepted record domain and actual source checkpoint
contract must determine which meanings are required to remain recoverable.

## A preservation guard and a progress guarantee

The continuation criterion can also define an admission rule. A checkpoint
contains the current state map and its retained receiver observations as
functions of the independently variable source records. Call it protected
when their combination is injective on that record domain. For each proposed
mean, compute the candidate state map and admit the mean exactly when the
candidate checkpoint remains protected. A rejected proposal changes no
scalar state. The test is a property of the response maps, not of a realized
payload. `linear_protected_iff` proves its joint-kernel implementation;
the proposed receiver sample is already a row of the candidate state.

The connected native completing word works from every protected checkpoint
whose retained record determines its current receiver scalar (`HasRoot`).
The initial sample establishes this additional hypothesis; each admitted
mean records the resulting receiver value, and rejection changes neither.
Thus it holds after arbitrary earlier admitted or rejected proposals. Moreover,
every prefix of that word is admitted by the guard. Otherwise a rejected
prefix would erase two source possibilities with equal retained records
and equal current states; its common remaining suffix could not separate
them, contradicting the completing-word theorem. The Lean proof establishes
equality of guarded and unguarded execution for every completing word and
then proves completion for any occurrence among arbitrary proposal noise.

This supplies a progress result for a precisely specified proposal law.
Let the source-labelled finite alphabet contain M event identities and the
completing word have B letters, with an explicit alphabet-to-seam encoding.
Under the uniform law on all proposal words, each block of B proposals has
M^B possibilities. A run still incomplete after k disjoint blocks must
avoid the completing word in every block. The number of such full proposal
words is bounded by `(M^B-1)^k`, against denominator `M^(B*k)`. Both this
count and its composition with actual guarded native executions are kernel
checked in `SourceTemporalAttempts`. The failure fraction is at most
`(1-M^(-B))^k` and tends to zero. Under the consistent product proposal law,
summing that geometric tail gives the analytic expectation bound B*M^B.
The empty-word case is immediate completion and is handled separately.

This is a deterministic transformation of the complete proposal law, not
uniform sampling from the accepted words. Rejection multiplicities matter.
The receiver retains the initial sample and each executed mean's sample;
rejected proposals incur proposal/guard work but no scalar mean or new
scalar sample. The external rational coefficient computation, rank test,
proposal metadata and guard instrument are not counted as native means.
They are explicit interface resources. The native construction does not
prove a physical implementation of that guard or derive an IID law from
canonical A3. The complete-alphabet/full-simplex/counting-reference
conditions of `RepairWordSchedule` specify the proposal-law specialization.

On the three-edge chain, the universal block is `[2,1,2,0,1,2]`: its six
means calibrate successively back from the receiver. A first proposal of
edge 0 is vetoed because it would irreversibly erase the independent source
difference. The proposal word `[0,1,2,0,1,2]` instead executes five means,
retains the veto in its six-proposal cost, and reconstructs both records.
Repeating proposal 0 forever remains protected but does not finish; its
probability is zero under the specified product law, rather than the word
being silently excluded. Exhaustive guarded replay covers all 9,841 words
through horizon eight and checks first completion, all rejection decisions,
receiver response rows, work and unconditional normalization independently.

## Specializing to M1

An arbitrary-continuation theorem identifies an additional temporal
obligation. If one event maps two preparations to exactly the same complete
scalar state, and the receiver's preceding samples agree, every subsequent
receiver transcript agrees under a common continuation. The theorem allows
an infinite word and a decoder of the entire infinite transcript. Thus
waiting cannot recover the erased distinction. Any sound partial policy
must abstain on these two possible initial values.

For the captured preparation, the initial difference is +1 at port 0 and
-1 at port 1. Their supported mean annihilates that difference on all
15,360 scalar ports. The receiver at port 45 had equal initial samples.
The individual records and their difference are lost; their sum is preserved.
Under a separately supplied uniform first-event law on the 46,050 captured
seams, this cylinder has probability 1/46050 independently of the waiting
horizon. Sound individual publication therefore has an abort probability
at least 1/46050 on each of the two coupled preparations. This probability
specialization follows from the exact seam census and the declared law;
the infinite-transcript indistinguishability is kernel checked. No canonical
A3 selection of an infinite schedule law is claimed.

This argument is specific to the unprotected preparation and read interface.
A code that retains the distinction elsewhere, an earlier discriminating
sample, or an independent record/controller channel changes the hypotheses.
The protected encoded banks are not refuted. The theorem locates a temporal
constraint they must satisfy: discriminate or preserve a required record
before an allowed operation erases all of its observable distinctions.

The source-derived available individual records are

```
Available(L) = { j : every z in ker L has z_j=0 }.
```

For direct observations of a subset S, a meaning sum_j a_j x_j is
identifiable exactly when the support of a is contained in S. Each omitted
nonzero coefficient has a one-coordinate deletion witness. A fixed added
unit changes no equality test. Thus the required direct reads follow once
the source supplies the particular affine record meaning.

This minimality extends beyond linear meanings. On a finite independent
product domain, call coordinate j essential if two inputs differing only
at j can have different public meanings. Direct access to S determines an
arbitrary meaning f exactly when S contains every essential coordinate.
Necessity is the isolated-coordinate intervention. Sufficiency replaces
unobserved coordinates one at a time: each replacement is inessential and
leaves f unchanged. `SourceTemporalEssential` proves both directions and
the unique minimal set. It does not assume that a finite list of examples
exhausts the meanings, and it does not apply unchanged to correlated
admissible preparations. The source's interpretation map supplies f.

Aggregate access has a different implication. Observing x_0+x_1 permits
that sum to be published while neither x_0 nor x_1 is identifiable. Native
captured words realize this distinction: beginning with the mean of source
ports 0 and 1 erases their difference, while transporting port 1 before
that merge permits both initial records to be recovered. Reversing the
path's temporal order supplies neither record to the receiver.

The bank compiler in source_native_programs derives finite execution for
supplied request lists. Its decoded recurrence alone does not choose the
individual-read interpretation of the history. M1 selection requires an
A1-produced meaning/readout interface whose available and required record
sets coincide with the claimed metric neighbours, plus the actual A3
cover/reference and complete temporal constraints. Neither the support
topology nor the existence of a correct compiler supplies that equality.

## What the radius investigation establishes

Three candidate arguments for the exact radius were checked against the
existing source results.

* Native accessibility does not single out a proper metric ball: the
  connected completion theorem recovers every scalar port at a finite
  regulator. A radius from a deadline would require a source-derived
  deadline, concurrent controller and relation between native work and the
  source-record metric. None is supplied by connectivity alone.
* The covering/causal-limit theorem in `SourceNetOrderLimit` consumes a
  positive read radius a and the condition h/a -> 0, with its layer duration
  a/c. With h proportional to L/q, every a=L*q^(-alpha), 0<alpha<1, shrinks
  to zero and has h/a -> 0. The square-root choice is a valid member of this
  family. The limiting-cone conclusion cannot distinguish that exponent.
  This is an inspection of the exact theorem hypotheses, not a full-axiom
  countermodel or an assertion that all these clocks are physical.
* Balancing radius against covering error can motivate a square-root scale
  only after choosing the quantity to minimize and its constants. The
  canonical A3 object-type gate does not supply such a cost objective or
  optimize between unrelated menu laws. Adding that objective here would
  assume the missing selection instead of deriving it.

The canonical A1 architecture supplies the response object and spherical
support, while its readout, record domain and continuation interfaces must
be instantiated. A2 naturality preserves an actual interpretation; it does
not choose an interpretation from every function of the records. A3 needs
the actual observable grammar, cover and faithful reference before a policy
can be read from its optimizer. The finite native-recovery and temporal-feasibility theorems do not derive
these source-selection inputs from the existence of a compiler. The full A1
response/endogeneity receipts are themselves explicitly distinguished from
the released architectural fixtures in `docs/AXIOM_REFERENCE.md`.

## Axiom obligation map

| Canonical obligation | Result on this finite interface | Retained source input |
| --- | --- | --- |
| A1 records, instruments, continuation | Observation responses derived from preparation and native prefixes | The interface's selection, full carrier response, accepted domain and refinement tower |
| A2 agreement and naturality | Complete linear fiber factorization and invariance under invertible recharting | Nonlinear/quantum data, overlap implementation and endogeneity |
| A3 complete constraints and optimizer output | Exact support-wide read criterion for an attained classical full-history KL minimum | Actual source grammar, entropy cover/reference/weights and physical publication |
| M1 metric reads | Exact required-direct-read theorem for a specified affine meaning | Metric-radius selection, mandatory publication and individual read-from identification |

These are conditional theorems on source observations. They do not constitute
a full three-axiom model or a countermodel to entailment of M1 by that basis.
