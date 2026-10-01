# Full-interface leakage, dissipation and source refinement

The [fixed-rate source construction](FIXED_RATE_SOURCE_READS.md) protects
against dephasing which preserves the native coordinate codes. A generic
environment need not preserve them. Merely extending the final accounting
observable to a failure symbol does not correct leakage: a single leaked
carrier can then reject an otherwise correct whole experiment.

Here the existing complete source interfaces implement local leakage
reduction. This gives a source realization with full-interface dissipative
or coherent noise, and a second refinement regime in which elementary
errors need not vanish at all. These are constructions in the same declared
source, not source selection or an autonomous clock from A1--A3.

## 1. Native leakage reduction, including its private branches

Let V_d embed the first d coordinates in C^6, d=2 or 4, and P_d=V_d V_d*.
The source's complete transfer already distinguishes P_d from its complement
and retains a success/failure record. Follow its failure output by the fixed
blank reset, while passing a successful code output unchanged. The resulting
instrument, with logical output in C^d, has Kraus operators

```
A_ok = V_d*,
A_j  = |0>_d <j|_6,       j=d,...,5.                       (1)
```

All A_j share one public failure label; j is a private Kraus index. With the
flag ignored only for defining the application marginal, the channel is

```
R_d(rho) = V_d* rho V_d + Tr((I-P_d) rho) |0><0|,
sum A* A = I_6,       R_d(V_d x V_d*) = x.                 (2)
```

The source transfer produces an orthogonal failure symbol. The continuation
from that symbol to a blank is the existing finite conditional reset service,
not conditioning the experiment on success. Both branches finish in a code
carrier; the public failure flag remains in an auxiliary record. Its
classification as a protection diagnostic is fixed before execution.
There is no attempt to restore unknown information which already leaked.
An error-correcting block restores logical information from other carriers.

These operations are physical finite CPTP services with counted durations.
Their faults, including a wrong flag or a leaked reset output, are ordinary
faults of the whole reduction service. They are never assumed reliable.
Sections 1--6 retain the fixed-rate parent's classical-control convention.
Section 7 composes this construction with the now merged
[noisy-record compiler](NOISY_SOURCE_RECORDS.md) and removes that runtime
convention for fixed bounded quantum generators and classical fault rates.

A mathematical isometry for (1) is

```
W psi = A_ok psi |ok,0>_E
             + sum_(j=d)^5 A_j psi |fail,j>_E.
```

It preserves the entire six-dimensional input when E is kept. Its use in a
proof does not make E accessible or add a reversible full-M6 control. Actual
reduction remains the rank-deficient complete source instrument. The protocol
adds no reversible full-carrier command to the twelve-port reference grammar.
The noisy finite apparatus is a perturbation of that reference, not an exact
realization of all its algebraic identities. In particular a generic coherent
drift is not an extra controllable source generator, and no unchanged exact
response-rank theorem is asserted for an arbitrary perturbed device.

For a two-qubit source gate, reduce each input carrier separately, pack the
two code qubits by the fixed coordinate map, execute the existing native
pulse sequence, and unpack into separate carriers. The packed processor is
C^6 with a four-dimensional code, not two freely available six-level systems.
It is an internal component of a bounded gate box with two external qubit
carriers. Leakage in its pulse, pack, unpack or reset can damage both gate
operands, and is charged to that one two-qubit box. No interaction with a
third data operand is hidden in its helpers. A faulty output is reduced again
before its next gate. Long block flights have reductions at their endpoints;
no mid-flight correction or additional long traversal is required.

The finite evidence derives (1) from the parent's transfer Kraus maps, then
checks (2) and the complete flagged instrument independently on all 36 input
matrix units for both d. A two-input native CZ box is checked on all 1,296
matrix units of C^6 tensor C^6, with the two public flags retained. Private
Kraus labels never become controller data. These channel identities include
arbitrary input references by linearity.

## 2. Why a local leaked carrier is correctable

Let J be the complete seven-qubit code isometry in the parent and D_s its
64 syndrome-decoder Kraus maps. For any channel N from one code qubit to its
full six-dimensional carrier, R_2 N is a CPTP one-qubit channel. Each of its
Kraus maps lies in the span of the four Pauli operators. Consequently

```
sum_s D_s [ (R_2 N)_i tensor identity_others ]
             (J rho J*) D_s* = rho                         (3)
```

for every site i, rho and spectator. More explicitly, the operators
`|a><b|`, 0<=a<6, 0<=b<2, span every possible code-to-carrier fault.
For every such operator F, every reduction branch A and every syndrome s,
`D_s (A F)_i J` is a scalar multiple of the logical identity. The finite
checker verifies these branch identities for all seven sites. This is a
linear-basis certificate, not a sample of leakage probabilities or input
states. It also checks the trace-preserving sum for finite amplitude-damping
and coherent-leakage channels, retaining the failure branch.

The physically meaningful full-carrier decoder is therefore

```
D_complete = D_QEC after tensor_i R_2.                     (4)
```

It is TP on the entire carrier space. Its implementation in the apparatus
uses the faulty, charged reduction/recovery services; the ideal version (4)
defines the operational comparison. It does not throw away a single leaked
carrier by testing the much larger product P_carrier as one acceptance event.

This difference is observable. A full erasure/reset of one physical code
qubit is corrected by (3), even though the old product-carrier test can reject
that branch with probability one. For N independent carriers exposed to a
fixed positive leakage hazard during a time of order a/c, the probability of
no leak can be `(1-Theta(1/q))^N`; for N=Theta(q^3) it vanishes. A demand for
zero microscopic leakage is thus stronger than the read requirement.
Correction does not require a leakage-free global final state.

For the parent's positive logical accounting E, extend application aborts
by E_max and define `E_complete=D_complete*(E direct-sum E_max)`. Positivity
and unitality give `0<=E_complete<=E_max I`. If decoded state/history error
is epsilon in half trace norm, accounting error is at most E_max epsilon.
It remains the same declared logical accounting, with E_max=O(q^4), not a
physical drive Hamiltonian or an equality of raw implementation counts.

## 3. Continuous dissipation cannot be treated as a success coin

On each counted native service let the actual generator be

```
L_t(rho) = -i[H_t+B_t,rho]
         + sum_alpha (L_alpha rho L_alpha*
                      - {L_alpha* L_alpha,rho}/2).          (5)
```

H_t is the specified native drive. B and L may mix the code and complement,
relax levels, or return population from the complement. Put
`Gamma(t)=2||B_t|| + 2 sum_alpha ||L_alpha(t)||^2`.
The completely bounded trace norm of the extra generator is at most Gamma:
left/right multiplication has bound given by the operator norms, also with
any reference tensor factor. Duhamel's formula, using contractivity of the
actual CPTP propagator and the ideal propagator, gives

```
||T_[0,h] - U_[0,h]||_diamond <= integral_0^h Gamma(t) dt.  (6)
```

The large native drive norm cancels from this estimate. It does not require
commuting the noise past a pulse. Piecewise ideal CPTP instruments and blank
preparations can be inserted: telescope their noisy intervals, including
every idle, and use contractivity again. The finite native decomposition
has a fixed number of components, which enlarges the constant. For fixed
bounded Gamma and the charged service/flight horizon h<=C a/c, the resulting
channel error delta is O(q^-1).

Unlike dephasing, amplitude damping generally has **no positive ideal-channel
weight**. For damping probability gamma in (0,1), its two qubit Kraus maps are
`K0=diag(1,sqrt(1-gamma))` and `K1=sqrt(gamma)|0><1|`.
Write the unnormalized Choi vectors in input-output order. The vector
`w=(sqrt(1-gamma),0,0,-1)` is orthogonal to both Kraus vectors but not to
`vec(I)`. For every p>0,

```
w* [J(damping)-p J(identity)] w
       = -p (1-sqrt(1-gamma))^2 < 0.                       (7)
```

Thus subtracting a purported ideal no-jump branch leaves a non-CP remainder.
The real no-jump evolution attenuates amplitudes and is not the ideal gate.
Our verifier explicitly tests this obstruction. It never calls a coherent
fault norm a failure probability.

For Markov channels use the dimension-independent Stinespring continuity
bound of [Kretschmann, Schlingemann and Werner, Theorem 1](https://arxiv.org/abs/quant-ph/0605009):
dilations can be aligned with distance at most sqrt(delta). Ancillary
unitary completions change only a fixed constant. This supplies a conservative
fault-amplitude strength

```
eta_M <= C sqrt(delta) = O(q^-1/2).                        (8)
```

The alignment is mathematical; no environment control is performed. Each
Markov service has a fresh environment. Keeping the environment in the proof
handles arbitrary reference-entangled input and avoids an unjustified
classical mixture over faults. This square-root route is sufficient; it is
not asserted optimal.

The finite checks compare driven full-M6 propagators built by a Kronecker
generator and an independent matrix-unit generator with a different matrix
exponential implementation. Complete positivity, trace preservation, code
leakage, the reset continuation and all input coherences are checked. Thermal
return branches and noncommuting drives are included, not projected out.

## 4. Memory-bearing coherent noise and the leakage reduction theorem

A second admitted noise class is a local Hamiltonian coupling to a common
bath. In each circuit location ell the perturbation K_ell may act on that
location's carrier(s), their leakage spaces and a shared bath. The bath may
retain memory and interact internally without a reset. Direct noise on
several unrelated simultaneous data locations is not covered merely because
it has bounded spatial range. Require a uniform bound on the coupling per
active location, and include each carrier's waits and transports.

In the interaction picture the integrated bound supplies a local fault
amplitude `eta_H<=C max_ell integral ||K_ell(t)||dt`. A fixed bounded coupling
therefore gives eta_H=O(q^-1). The bath dynamics and ideal native drive do not
increase the interaction-picture operator norm. The usual marked-location
Dyson expansion bounds the sum of paths faulty on any specified set of k
locations by eta_H^k, without independence of the reduced faults.

We import [Aliferis and Terhal, Lemma 2 and Theorem 1](https://arxiv.org/abs/quant-ph/0511065)
for noisy leakage-reduction units and
[Aliferis, Gottesman and Preskill, section 11](https://arxiv.org/abs/quant-ph/0504218)
for the universal concatenation argument. The source's contribution is that
(1) implements the needed reduction using admitted complete services. The
reduction is repeated before each elementary gate, including correction and
non-Clifford gadgets. A faulty reduction belongs to its gate box; it is not
excluded from the noise budget. The finite box size changes library constants,
not the fault-power recurrence. The pair processor is an internal component
of the two-carrier box described in section 1, so its non-tensor C^6 space is
not mistaken for two persistent leakage subsystems.

For fixed finite library constants A>=1 and C, on a sufficiently small
cofinal tail these results give the conservative effective strength

```
eta_k <= A^(-1) (A eta)^(2^k),
epsilon_noise <= C M0 eta_k,       M0=O(q^4 log^6(2q)).     (9)
```

A may absorb the usual exponential combinatorial factors and the bounded
native reductions. It is not the pair count of the parent's Clifford-only
census and is not a measured native threshold. Overlapping recovery circuits
are handled by the credited theorem. The parent supplies its complete
universal gadget hypotheses, preparations and non-Clifford construction;
we do not substitute a passive memory polynomial for those gadgets.

For completeness the norm comparison is operational: the small total bad
amplitude bounds the difference of decoded channels in half trace norm, with
a fixed conversion factor included in C. Correctable paths give the entire
named instrument, not just one input distribution. Application outcomes and
aborts are retained. Leakage flags created by protection are diagnostics,
not global abort commands; forcing every such flag to abort would defeat (9).

The shared-bath control executes two interactions with the same finite bath
and compares their full channel with resetting the bath between them. They
differ even on native code inputs followed by the permitted flagged readback.
It demonstrates why independent channel multiplication is insufficient;
it does not claim to simulate every bath covered by the analytic bound.

## 5. Two quantitative refinement regimes

### Fixed physical generator or coupling

For fixed bounded Markov Gamma, (8) and **five** levels give
`q^4 (q^-1/2)^32 log^6 q = O(q^-12 log^6 q)`. Four levels would leave the
accounting majorant of order log^6 q and do not establish its convergence by
this conservative route. For bounded coherent bath coupling, eta_H=O(q^-1)
and **four** levels give the same bound. This distinction between a channel
error and a coherent fault amplitude is essential.

Synthesize requested native rotations to q^-16 as in the parent. Contractivity
and telescoping over O(q^4 log^6 q) gates give the same O(q^-12 log^6 q) term.
Thus, in both cases,

```
joint decoded state/named-history error = O(q^-12 log^6(2q)),
additional accounting error            = O(q^-8 log^6(2q)). (10)
```

The rates and coupling norms do not have to decrease. The finite code member
still has noise and leakage; (10) is an operational refinement bound.

### Fixed nonzero error per service

Reducing service time does not automatically reduce a fractional pulse error.
If B=epsilon H during a pulse of fixed angle, its integrated norm is of
constant order epsilon, even as the drive grows. This is a different regime.
Suppose the **complete native boxes**, including reduction, preparation,
measurement, resets and flights, have a common strength eta with
`theta=A eta<1`. For fixed nonzero theta choose

```
k(q)=max(1,ceil(log2(16 log(2q)/log(1/theta)))).             (11)
```

Then theta^(2^k)<=(2q)^-16. Equations (9)--(10) follow again. No elementary
fault strength tends to zero in this regime: the growing redundancy supplies
the improvement. The threshold is an explicit hypothesis on a given noisy
device family, not a fitted value or a new axiom of OPH. Fixed bounded physical
rates automatically enter the first regime on a cofinal tail.

This result also permits a fixed small calibration contribution together
with fixed-rate dissipation, provided their combined amplitude bound lies
below the same threshold. A local diamond bound by itself is sufficient for
Markov channels via (8); single-time marginals alone do not bound an arbitrary
correlated bath. For that case the local coupling/fault-amplitude condition
of section 4 is required.

## 6. Causal schedule and the full storage ledger

Let s>=2 bound the native size, width and depth multiplier of a complete
protection level, including leakage reduction and all idle slots; put
`kappa=ceil(log2 s)`. It is a symbolic fixed library constant. From (11),
`s^k(q)=O(log^kappa(2q))`. A conservative envelope is

```
active quantum inventory   O(q^3 log^kappa(2q)),
scheduled depth            O(q log^(6+kappa)(2q)),
active location volume     O(q^4 log^(6+2kappa)(2q)),
retained diagnostic bits   O(q^4 log^(7+2kappa)(2q)),
total storage-time volume  O(q^5 log^(13+3kappa)(2q)).       (12)
```

The last line is a conservative inventory-times-scheduled-slots bound. It
counts the complete lifetime of the diagnostic records, even
though they are not the named output marginal. For a fixed number of levels
the extra powers kappa can be replaced by zero. Native resets and new blank
supplies are included in the active volume. No entropy-erasure energy or
bounded spatial hardware density is inferred from these counts.

Let B_q be the longest compiled local depth, bounded by
`O(log^(6+kappa)(2q))`. Set the positive local service time and drive to

```
h_local = a/(1024 c B_q),     Omega = pi/h_local,
cluster radius <= c h_local/128.                           (13)
```

Shared processors serialize their own operations; separate processors work
in parallel. Transport all constituents of an encoded block concurrently,
one counted interface per constituent. They cross the original geometric
link once, and its duration is at least length/c. Reduction and recovery
take place at endpoints. Local depth fits into the same fixed padded source
slot by (13), even for k(q) growing as log log q. A fixed q-independent
wall-time rescaling therefore survives. The source's geometry, causal reads
and the parent clock/velocity dilation ratio are preserved, while inventory
and local drive costs grow as stated. This is a cofinal apparatus family,
not a fixed-strength or bounded-inventory laboratory device.

Only the named logical marginal obeys (10). The extra leakage/repair records
can distinguish the noisy apparatus from an ideal one. They remain in
separate central carriers and never instruct an unprotected global abort.
Unknown raw-state encoding, physical source/parameter selection, autonomous
clock generation and calibration of physical energy remain outside the exit.

## 7. One construction with leakage, noisy control and a live public archive

The [merged noisy-record result](NOISY_SOURCE_RECORDS.md) uses independent
code-preserving jumps. We now replace its quantum-noise restriction by the
full-interface Markov generators of section 3, while retaining its fixed
finite local classical Poisson fault rate. The quantum services have fresh
Markov environments, including exports and resets; the classical services
have the parent's disjoint local fault drivers. The static schedule remains
specified. No runtime classical decision or archive bit is reliable.

Use its explicit [[2047,1,63]] computation code and measurement-free universal
library, not the seven-qubit code of section 2. Its 31-error capacity exceeds
the conservative correction-and-gate spread 16. Insert native reductions
before elementary quantum gates, including those implementing Boolean
decisions and feedback. Reduction/reset feed-forward, if used, is confined
to its own bounded box and is itself faulty; it is never a trusted global
controller. The leakage-reduction theorem turns these into ordinary local
fault boxes. The credited
[Aharonov--Ben-Or construction, sections 7--8](https://arxiv.org/abs/quant-ph/9906129)
supplies the measurement-free gadgets and norm fault-path argument. Sparse
paths preserve encoded matrix units with a spectator, so the same argument
controls the joint instrument, not only a final classical answer. The finite
code/spread evidence and archive circuit are recursively replayed from that
parent; we do not claim to execute its entire enormous protected machine.

There is a substantive change at the public interface. An export still uses
n disjoint decoding islands, each of fixed size because we choose a fixed
number of levels below. A fault in the protected copies belongs to the
protected-computation error, not to an assumed independent bare output.
An island with no faults decodes its correctable input; at most r faulty
islands therefore leave at most r wrong output bits. The parent's archive
refresh lemma also remains a statement about arbitrary faulty locations:
with at most r faults, a word of n=4r+1 maintains its live majority, and the
new word has at most r errors. That statement does not need fault
probabilities. On central bits, expand the actual stochastic errors in
classical basis operations; on quantum islands keep a Stinespring dilation.
Linearity then applies the same correction statement to coherent fault
paths. It is not a claim that an unmaintained bare bit is reliable.

Here is the required norm replacement for the parent's stochastic tail.
For a region of V locations write the evolution as its ideal/fault expansion,
and let F(S) be the sum of paths faulty at every location in S, with all
other locations unrestricted. The local-noise estimate gives
`||F(S)|| <= eta^|S|`. The sum of paths with at least m faults is exactly

```
sum_(|S|>=m) (-1)^(|S|-m) binom(|S|-1,m-1) F(S).
```

For a path with k faults its coefficient is zero if k<m and one otherwise,
by the binomial identity checked independently in the evidence. For k>=m,
write its sum as `m binom(k,m) integral_0^1 x^(m-1)(1-x)^(k-m) dx=1`;
the finite binomial expansion and the beta integral prove the identity at
every size. Thus the norm is at most

```
sum_(j=m)^V binom(V,j) binom(j-1,m-1) eta^j
 <= binom(V,m) eta^m (1+eta)^(V-m)
 <= exp(V eta) (V eta)^m.                                (14)
```

This also bounds bad-island paths: r+1 faulty islands require r+1 elementary
faults among O(n) locations. Use (14) on those unconditional marked-location
sums and bound the protected bad part separately. Conditioning the noise
on protected success would invalidate the argument. Small total bad norm
bounds the decoded channel error by a constant times that norm, including
arbitrary references; no division by an acceptance probability occurs.
Sequential regions can be grouped at the first violation of the archive
invariant. Equivalently, inclusion-exclusion gives a factor exponential in
the sum of their bad-norm bounds; it is absorbed into a fixed constant on
the cofinal tail where that sum tends to zero. No probabilistic union bound
on hypothetical independent coherent faults is being used.

Take r=ceil(log2(2q)), as in the parent. Classical fault probability O(q^-1)
has a dilation fault amplitude O(q^-1/2), and (8) gives the same bound for
the general quantum noise. A refresh has V=O(log^4(2q)) locations; a fixed
size island times n has fewer. Put u_q=C log^4(2q)/sqrt(q). The joint export
and live-history error is bounded by

```
epsilon_archive <= C' q^5 log^8(2q) exp(u_q) u_q^(r+1).
log(epsilon_archive)
 <= -(log q)^2/(2 log 2) + O(log q log log q).              (15)
```

This again beats every inverse power on a cofinal tail. The coefficient of
the squared logarithm is halved relative to the stochastic proof; silently
reusing its old bound would be wrong. With **six fixed protection levels**,
the active volume M_q=O(q^5 log^16(2q)) gives
`M_q O((q^-1/2)^64)=O(q^-27 log^16(2q))`. Five levels would still converge
but would not retain the parent's stronger joint-error rate. Synthesize
requested rotations to q^-24 as in the noisy-record parent. We obtain

```
joint decoded output / live named-history error = O(q^-20 log^16(2q)),
additional positive accounting error            = O(q^-16 log^16(2q)). (16)
```

Use the local-reduction accounting extension (4), now followed by the
complete computation-code decoder and majority interpretation of the public
archive. The seven-qubit finite certificate illustrates local correction;
the distance/spread proof of the actual computation code justifies this
different decoder. Every named application abort is still an output.

Six levels and the native reductions have constant overhead. The parent's
active inventory O(q^4 log^8 q), depth O(q log^8 q), active location volume
O(q^5 log^16 q) and complete diagnostic storage-time O(q^6 log^25 q) remain
valid with larger constants. All diagnostic lifetimes are counted. Their
independent local noise cannot feed back through a discarded subsystem.
The archive retains its actual noisy central bits and live majority, and
its dynamic feedback uses protected registers. Its faster local schedule,
parallel block flights and q-independent wall-time dilation also survive.

This composition is for the specified fixed-rate Markov/classical model.
It does not transfer a fixed per-operation error to the parent's unencoded
archive: then V grows while its elementary error stays constant, and neither
(14) nor the old clean-island estimate supplies convergence. The separate
fixed-error theorem in sections 5--6 has its stated classical-control
convention. Nor do we assert a common-bath live-archive theorem from marginal
noise rates alone. These distinctions prevent exchanging incompatible noise
hypotheses between the two constructions.

The combined result removes both native-code preservation and reliable
runtime records from a single constructed read/clock process. It retains
the source/CCG capability, external schedule and declared logical accounting;
none is relabeled as a consequence of A1--A3. The executable
[package](../code/m1_leakage/README.md) separates exact identities, numerical
controls, new composition estimates and credited all-size results.
