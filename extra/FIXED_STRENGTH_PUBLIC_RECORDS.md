# Fixed-strength quantum computation and live public records

The [noisy-record construction](NOISY_SOURCE_RECORDS.md) maintains actual
central bits, but its voter circuit has O(n^4) fault locations per refresh.
Its sufficient error bound requires the elementary error to shrink with q.
Increasing quantum protection alone does not repair that public interface.
This construction replaces the voter graph, and supplies a decoder whose
export error stays bounded as quantum protection grows. Together these give
convergent quantum outputs and live public histories at a fixed nonzero
elementary error strength.

This note is based on merged main and is self-contained with respect to the
native reduction and fault-norm estimates. It also composes with the native
leakage analysis of [PR #1007](https://github.com/FloatingPragma/observer-patch-holography/pull/1007).
The construction does not require that separate branch to be present when
its evidence is replayed. Source selection, the external schedule and the
logical accounting reference retain their declared status.

## 1. An explicit bounded-degree voter graph

Use vertices (x,y) in Z_m squared, N=m^2, and these eight labelled neighbours:

```
(x +/- 2y,y), (x +/- (2y+1),y),
(x,y +/- 2x), (x,y +/- (2x+1)),                  modulo m.       (1)
```

Keep multiple edges and loops. Each label is a permutation; its inverse is
another label. The symmetric adjacency A has row sum 8. We import the
Margulis--Gabber--Galil bound on the absolute nonconstant spectrum,
`||A/8 - J/N|| <= 5 sqrt(2)/8`. The exact construction and bound are
[formalized in the Archive of Formal Proofs](https://isa-afp.org/browser_info/current/AFP/Expander_Graphs/Expander_Graphs_MGG.html),
including `mgg_numerical_radius`; the cited development credits the original
Gabber--Galil and Jimbo--Maruoka proofs. This repository does not claim to
formalize that imported theorem itself.

Take length-12 labelled walks. Their adjacency B=A^12 has degree D=8^12 and
is a sum of D permutation matrices. With P=B/D,

```
||P-J/N||^2 <= (25/32)^12 < 1/16.                            (2)
```

This is an explicit family at every square size, not a search conditioned on
successful storage. D=68,719,476,736 is a very large constant. The purpose is
an existence theorem with a proved positive noise interval, not an optimized
hardware threshold. Neither D nor the depth below increases with N or q.

For a set E of bad input positions, e=|E|<=N/8, let T be voters receiving
at least half their D inputs from E, with multiplicity. Put f=1_E-e/N. For
v in T, `(P f)(v)>=1/2-e/N>=3/8`. Equation (2) implies

```
|T| (3/8)^2 <= ||P f||_2^2
              <= (1/16) e(1-e/N) <= e/16,
|T| <= 4e/9.                                                (3)
```

An even number of votes is harmless: choose the lower middle sorted bit
(ties give zero). Counting every tie as potentially wrong proves the same
bound for both logical labels. No random graph or independence of voter
outputs is used.

## 2. A physical refresh and its live invariant

Keep N shared input bits, D private vote bits per voter and N output bits.
All are noisy central carriers. One cycle is:

1. Reset every private and output bit; retain the shared inputs.
2. For each length-12 word, copy the corresponding permutation matching
   into that voter's private position. One round touches each shared input
   once. A gate never touches two shared input positions.
3. Sort each private lane with the fixed adjacent bubble network and copy
   its lower middle bit to its output. Voters execute in parallel.
4. At a fixed scheduled handoff designate the output bank as the public word.
   Until then the old input bank remains public.

COPY uses a reset destination and the existing CNOT service. An AND/OR
comparator is the finite Toffoli/CNOT construction in the parent, with fresh
scratch confined to its voter lane. No majority oracle is a physical gate.
Every layer includes idles on all live input, vote and output bits. A reset
includes the preceding wait; the first input idle includes the wait since
the previous handoff. A bounded native implementation enlarges constants.

Before that enlargement, width is W=N(D+2), depth is
`H=2+D+D(D-1)/2`, and locations V<=W H=C_A N. Hence C_A is fixed, although
large. Addresses and the permutation rounds are preloaded in charged
prehistory. The design uses a static bank schedule, not a reliable selector
bit. Active protection of a public word continues until its last named cut.

**Live refresh theorem.** Start with at most floor(N/16) wrong shared bits.
If at most floor(N/64) service locations are faulty during the complete
cycle, the new word again has at most floor(N/16) wrong bits. At every
intermediate time, the old word has the correct majority.

Proof. Let F count faulty locations. The union E of initially wrong shared
positions and shared positions ever touched by a fault has size
`e<=N/16+F<=5N/64<N/8`. A faulty gate touches at most one shared position.
This union bound covers a source that changes between its different copies.
At most F private voters are damaged. Every other voter receives the correct
bit from every source outside E and executes exact private Boolean logic.
By (3), wrong outputs number at most

```
(4/9)(N/16+F)+F <= 29N/576 < N/16.                         (4)
```

The old word has at most 5N/64 wrong bits throughout the cycle, including
inside a faulty native service. All affected operands of that service have
the same shared-position/private-voter ownership. This proves the live
claim, not merely correctness after a final decoder. Taking integer floors
only strengthens it. Consecutive cycles need no independent success events.

The finite evidence checks the actual eight matching maps and their powers,
exact walk counts and mixing certificates. It executes smaller unpowered
native schedules, including every idle and every single-location corruption,
with two independent Boolean engines. Those circuits test the implementation
and fault ownership; they are not labelled executions of the D-degree
apparatus. Compressed powered-graph controls evaluate worst-case damaged
sources and private voters using integer multiplicities. The all-size
guarantee is (1)--(4), not an extrapolation from those examples.

## 3. A fixed error threshold for the entire maintained history

First, for local stochastic faults of strength p, a bad cycle requires
`M=floor(N/64)+1` faulty locations. Thus

```
Pr(bad cycle) <= binom(C_A N,M) p^M
              <= (64 e C_A p)^(N/64)                      (5)
```

when the base is at most one. Unlike the parent's O(N^4) voter circuit, this
has a fixed positive sufficient threshold independent of N.

We also need a norm version for quantum export and general channel faults.
Use a joint local fault expansion on the complete scheduled services. For
any specified set S, the sum F(S) of paths faulty on S, with all other
locations unrestricted, satisfies `||F(S)||<=eta^|S|`. Fresh-environment
Markov channels uniformly close to the ideal complete services supply this
condition: [Stinespring continuity](https://arxiv.org/abs/quant-ph/0605009)
gives eta<=C sqrt(delta) for diamond error
delta, including the public outcome register. The mathematical environment
alignment grants no physical environment control. Classical stochastic
channels have such a dilation too. Marginal error rates of an unspecified
correlated bath alone are not this hypothesis.

For V locations, the sum of paths with at least M faults is exactly

```
sum_(|S|>=M) (-1)^(|S|-M) binom(|S|-1,M-1) F(S).            (6)
```

A path with k faults has coefficient zero for k<M and one otherwise. The
binomial expansion of `M binom(k,M) integral_0^1
x^(M-1)(1-x)^(k-M) dx=1` proves the latter identity. Consequently its norm is

```
beta(V,M) <= exp(V eta) (e V eta/M)^M.                     (7)
```

For the archive choose `eta<=1/(1024 e C_A)`. Since M>=N/64, (7) gives

```
beta_archive <= exp(-N/32).                                (8)
```

Constants refer to complete native services, with all helper operations,
idles, waits and intrinsic read/reset faults included. Noise at an event
does not disappear by declaring the event instantaneous.

## 4. Export without a growing unprotected decoding island

The quantum record R is classical in the complete compiled instrument.
Copy it, while protected, into N independently protected probe blocks. Each
copy and its recovery belongs to the protected computation. Never decode R
to one raw bit and broadcast that bit; one terminal fault would corrupt the
whole archive with a probability bounded away from zero.

Use the parent's actual [[2047,1,63]] measurement-free computation code and
its distance/spread argument. The credited
[Aharonov--Ben-Or construction](https://arxiv.org/abs/quant-ph/9906129),
sections 4.9 and 7--8, supplies the finite decoding and sparse-path gadgets.
Its classical-record decoder may retain label-dependent private garbage:
R has been dephased, and no coherence between different public labels is
promised. It preserves the conditional quantum state and arbitrary spectators.
Here the complete classical decoder is executable on all 2047 bits. It uses
[Reed's majority-logic method](https://doi.org/10.1109/TIT.1954.1057465), with
an explicit puncture treatment. For an extended RM(5,11) truth table, recover
monomial coefficients in descending degree. For a degree-j monomial, take
the parity on each of the 2^(11-j) disjoint coordinate j-cubes and majority
those parities. After higher-degree terms have been removed, each clean cube
returns that coefficient. At most 31 errors damage at most 31 cubes, while
there are at least 64 cubes; majority is therefore exact. Subtract the
recovered monomial and continue. The residual errors stay on the same sites.

Try both values of the missing all-ones evaluation point. One fill has at
most 31 errors. Re-encode each result and keep it only if its punctured
distance from the received word is at most 31. Two distinct accepted words
would contradict distance 63. The value of the recovered polynomial at the
missing point is the CSS logical label. If neither fill passes, retain a
failure diagnostic and output the fixed zero fallback. This is a complete
decoder, not rejection of the experiment. Its parity, majority and distance
calculations have finite reversible Boolean circuits with charged scratch.
The proof covers every error set through weight 31. Two independent engines
execute the full-size code; a 32-error wrong-label example shows why the
radius cannot be silently extended.

The classical decoder also has the required quantum-record meaning. Let
C_a be the punctured codewords whose missing-point value is a, and let H_a
span the computational basis words within distance 31 of C_a. Distance 63
makes H_0 and H_1 orthogonal. A reversible implementation computes the
decoder label in a fresh register while retaining its input and workspace.
On H_a it therefore has the form `|psi> -> |a> W_a|psi>`, with W_a an
isometry. This identity extends linearly with an arbitrary spectator. Every
operator on at most 31 sites sends an encoded |a> into H_a: expand it in
Paulis, whose X supports change at most those sites and whose Z supports
only change amplitudes. Discarding the private workspace thus preserves the
conditional spectator state and returns a with certainty. This argument
does not assume that damaged inputs are classical strings or that their
errors are independent. Dephasing the logical record before export is
essential; the procedure is not a coherent decoder of an unknown qubit.

Let S bound a full leading-correction/gate rectangle. Use disjoint rectangles
of the computation-code construction, not overlapping seven-qubit recoveries.
Its sparse-path induction tolerates one faulty child; two are conservatively
declared bad. The marked-location calculation (7), also on unions of disjoint
rectangles, supplies a constant A>=e S^2 with

```
eta_j <= A^-1 (A eta)^(2^j),       theta=A eta<=1/2.         (9)
```

Indeed `eta_(j+1)<=exp(S eta_j)(S eta_j)^2<=A eta_j^2`;
the bound S eta_j<=1 is preserved. Good paths have the exact decoded
instrument supplied by the computation-code proof. This is not a threshold
inferred from the finite seven-qubit Clifford census.

Decode a level-k probe **one level at a time**. At stage j, execute a fixed
complete one-level record decoder with level-(j-1) protected operations.
For j>1, first run a separate level-(j-1) protected identity/correction on
each of its 2047 child blocks, before any decoder gate couples children.
It outputs a level-(j-1) record; discard its private syndrome/workspace only
after the stage. Proceed through j=k,...,1, then perform the bare native read.
The stage includes its storage, correction and transfer. It never contacts
another export island or R. Standard preparations and the finite logical
decoder use the parent's exact universal library, so there is no fixed
approximation error inserted at every stage.

This entry correction uses the computation code's total correction property
(the cited construction's section 4.8 and the second assertion of Lemma 8):
even an arbitrary child input becomes sparsely close to some encoded qubit
when its correction path is good. It acts separately on each child. The
incoming block's at most 31 problematic children therefore become arbitrary
logical inputs at those same positions, while its other children keep their
correct logical state. The parent's spread bound of 16 fits inside that
radius. Representing all children by ideal inner decoders now puts the input
in H_a, so the preceding support argument applies. The remaining protected
decoder circuit operates on sparsely encoded children and leaves a sparsely
encoded output, supplying the next stage's invariant. For j=1 there are no
inner blocks to correct. All entry corrections are counted in S_D. Omitting
them would incorrectly apply the protected-gate simulation theorem to
possibly non-sparse child inputs. Existing correctable input damage is not
counted as a new fault or required to be absent.

Let S_D bound the number of lower-level rectangles in one decoding stage,
including its idles. The sum of their bad strengths is bounded uniformly:

```
t <= S_D sum_(j=0)^(k-1) eta_j + eta
   <= S_D eta/(1-theta)+eta <= (2 S_D+1) eta.               (10)
```

The last eta includes the final read and its wait to archive activation,
with bounded native factors absorbed into S_D. Since `2^j-1>=j`, the
geometric series proves (10) for every k. Inclusion-exclusion gives the
whole island bad norm at most exp(t)-1, hence

```
b_export <= C_E eta,       C_E=2(2 S_D+1),                 (11)
```

provided t<=1. This closes the growing-island problem: counting every native
location of an unprotected full decoder would instead give a bound growing
as s^k eta. The noisy decoder is actually constructed in protected stages.

Disjoint islands have multiplicative marked-bad bounds `b_export^r` for any
specified r islands. This follows by multiplying their marked-set expansions
on disjoint native locations; it does not condition on successful copying.
The latter bad event is counted separately. Export failure to start within
floor(N/16) errors therefore has norm

```
beta_export <= exp(N b_export) (16 e b_export)^(N/16)
            <= exp(-N/8),   b_export<=1/(256 e).             (12)
```

The output bits are physical noisy central records. Statistical majority is
the declared public interpretation, never a trusted controller operation.
Feedback uses the private protected R. Application aborts are ordinary
retained labels; protection diagnostics do not cause a global abort.

## 5. Full-carrier faults and one joint process

No code-preserving noise assumption is needed for the construction. For
d=2,4 the existing complete transfer has successful operator V_d* and
complementary maps to a failure symbol. Continue every failure with the
existing reset to |0>_d, retaining its public diagnostic flag. This implements

```
R_d(rho)=V_d* rho V_d + Tr((I-V_d V_d*)rho)|0><0|.           (13)
```

The individual complementary basis indices remain private Kraus indices.
Their coherent sum would not be this channel. Reduce each external operand
before a native gate; the pair processor's four-dimensional code sits inside
one M6 carrier and remains an internal part of that two-operand gate box.
Pack, pulse, unpack, flag and reset faults belong to the box. There is no
extra reversible full-M6 command or third data operand. The complete flagged
channels for both codes are checked on every M6 matrix unit in the evidence.
The native proper-code source in main supplies the required transfer/reset.
The noisy finite apparatus approximates that reference. Generic leakage or
coherent drift is not promoted to a new controllable source generator, and
the perturbed device is not asserted to obey the reference's exact response
identities.

The noisy leakage reduction theorem of
[Aliferis--Terhal](https://arxiv.org/abs/quant-ph/0511065), Lemma 2,
reduces these bounded noisy boxes to the local circuit-noise condition,
enlarging the library constants. Their faults are included in eta; there
is no ideal leakage detector. A correctable leaked carrier is reduced and
corrected locally, not turned into an application-wide rejection.

For completeness, unions of coherent bad regions need an intersection
argument. Assign every physical location to one protected rectangle, one
export stage or one archive cycle. The terminal read and any wait before
archive activation belong to export. Waiting after activation belongs to
the first refresh. Ancillas
and different archive words have separate lanes. Within a region use (6);
between disjoint regions the joint marked bound is the product. For their
bad majorants beta_i, inclusion-exclusion gives

```
||B_any|| <= product_i(1+beta_i)-1
          <= exp(sum_i beta_i)-1.                          (14)
```

This remains valid for successive cycles on the same physical bank because
the locations are disjoint time intervals. It uses the joint noise condition,
not a probability interpretation of fault amplitudes. Every remaining path
preserves both the protected application and the live archive invariant.
Normalization converts the small bad norm into a decoded channel error with
a fixed factor. Reference correlations and every named abort are retained.

Choose a single fixed positive eta small enough for (9), (11), (12), and (8):

```
eta <= min(1/(2A), 1/(1024 e C_A), 1/(256 e C_E)).           (15)
```

All constants are finite native-library constants, independent of q. This
is a nonempty interval, not an assumed successful M1 output or a fitted
threshold. It covers fixed nonzero full-interface Markov errors and fixed
nonzero local classical faults with sufficiently small corresponding
strengths. The source and its available controls are the same declared
capabilities used in the parents.

There is a concrete nonzero-noise inhabitant, rather than only an implication
about a possibly empty noise class. After each native quantum component use
the fresh-environment full-M6 channel
`E_p(rho)=(1-p)rho+p Tr(rho)|5><5|`, and at each classical component use an
independent bit flip with probability p. Apply the same rule to counted
storage locations. The quantum channel leaks every proper-code input with
positive probability. Its Kraus maps are `sqrt(1-p) I` and
`sqrt(p)|5><j|`, j=0,...,5. Against the ideal dilation their operator distance
is `sqrt(2-2 sqrt(1-p))<=sqrt(2p)`; a bit flip has the same bound. If r bounds
the number of components grouped into a service, choosing
`0<p<min(1/2,(eta_star/(2r))^2)` for any positive eta_star allowed by (15)
puts every complete box below eta_star. All q use this same positive p.
The complete M6 channel and dilation bound are independently checked.
This gives an explicit noisy process in the admitted source capability;
it is not a claim that physical nature selects that noise law.

## 6. Refinement, propagation and costs

For a completely integer schedule, enlarge A if necessary so that
`A>=max(3 S^2,768 C_A,192 C_E)` and restrict to the nonempty fixed interval
`0<eta<=1/(4A)`. Then theta<=1/4 and (15) holds, using e<3.
For q>=2 take

```
ell(q)=ceil(log2(2q)),  m(q)=ceil(sqrt(2048 ell(q))), N=m(q)^2,
k(q)=ceil(log2(16 ell(q))).                                (16)
```

Here theta is fixed, positive and at most 1/4. Thus N=O(log q),
`eta_k<=A^-1 (2q)^-32`, and k=O(log log q). Equations (8) and (12) make each
archive/export contribution at most (2q)^-64. The complete number of
scheduled records, refreshes and gates is polynomial in q with polylogarithmic
factors. Their joint bad norm therefore vanishes by (14).

More explicitly, let s>=2 bound one complete quantum protection level's
width/depth overhead, including native leakage reductions, and set
`kappa=ceil(log2 s)`. Staged decoding has geometric total size at most a
constant times s^k, not k times the largest stage. The parent's conservative
ledger becomes

```
active inventory                O(q^4 log^(8+kappa)(2q)),
scheduled depth                 O(q log^(8+kappa)(2q)),
active native location volume   O(q^5 log^(16+2kappa)(2q)),
retained diagnostic bits        O(q^5 log^(17+2kappa)(2q)),
complete storage-time volume    O(q^6 log^(25+3kappa)(2q)).   (17)
```

This overcounts the new constant-depth archive, and includes every old public
record's lifetime, inactive protected registers, resets, output waits and
the diagnostic record identifiers. Stored diagnostics are not reused by the
active computation. Their independent local Markov noise therefore cannot
feed back into its reduced state. The theorem does not infer a common-bath
decoupling claim from that Markov model.

The synthesis requests still number O(q^4 log^8 q): added Boolean, decoding
and protection gates use the exact fixed library. Synthesize each requested
rotation to q^-24. Noise contributes at most
`O(q^-27 log^(16+2kappa) q)` and synthesis at most `O(q^-20 log^8 q)`. Hence

```
joint decoded quantum / live named-history error
                         = O(q^-20 log^(16+2kappa)(2q)),
additional positive accounting error
                         = O(q^-16 log^(16+2kappa)(2q)).     (18)
```

For accounting use the complete local reduction (13), followed by the
computation-code decoder and public majority interpretation. Pull back the
parent's positive logical observable, assigning its named abort the declared
E_max=O(q^4). Complete decoding is TP, so the pullback remains positive and
bounded by E_max. This does not identify logical accounting with drive work.

Let B_q bound the longest local compiled depth, including the D-degree
archive's fixed constant and every decoding stage. It is
O(log^(8+kappa)(2q)). Take local service time a/(1024 c B_q) and cluster
radius at most c times that service time divided by 128. Serialize shared
processors, run disjoint lanes in parallel, and transport all constituents
of a protected block concurrently over each original link. Every constituent
still takes at least length/c; all correction and decoding occurs at
endpoints. The local work fits into fixed padded coarse slots. A common
q-independent wall-time dilation and the parent's clock/velocity comparison
survive. Faster local drives, parallel interfaces and growing inventory are
paid cofinal resources, not a claim about bounded-density fixed hardware.

## 7. Countercontrols and interpretation

A regular graph by itself does not suffice. With the identity voter graph,
N/16 old errors persist, and one damaged output outside that set leaves more
than N/16 errors. Its spectral contraction is absent. Omitting input idles
also invalidates the live claim, even if final outputs are clean.

A raw decoded broadcast has an error floor equal to its common-bit error.
An unprotected decoder of growing size need not satisfy (11). A numerical
control contrasts that growing serial exposure with the bounded recursive
series, and the complete export instrument includes all wrong labels rather
than renormalizing on majority success. No proof here grants a perfect naked
final bit after an additional noisy observation interval.

The substantive change is that vanishing elementary errors and reliable
runtime records are no longer required together: one specified apparatus
family supplies protected quantum behaviour, faulty decisions and a live
public archive at fixed subthreshold error strength. The additional expander
theorem is mathematical input, not a new OPH axiom. Source/CCG capability,
external scheduling and physical energy calibration retain the boundaries
stated in the [contract](../code/m1_expander_archive/CONTRACT.md).
