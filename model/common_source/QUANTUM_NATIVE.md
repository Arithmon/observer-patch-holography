# Coherent native repairs around one continuing quantum history

This construction joins the source-compiled C/E native instruments to the
[canonical four-read program](COMPOSITION.md). The first native record controls
the gauge kick; the second smooth record selects the seam of the last native
repair. The native output is a quantum register throughout. Its surviving
coherences and its correlations with the continuing work register are retained
inside every committed record block.

The resulting finite channel acts on the full 78-dimensional native carrier,
the work space `W_d` of dimension `4d^2`, and any untouched reference. This
extends the earlier classical native-input benchmark. It does not construct
physical load identification, unknown-state storage, or a full
`OperationalObserver`, and it does not turn approximate cutoff inclusions into
exact consensus-tower morphisms. The separate occupation comparison has the
declared scalar-vacuum premise described below.

[quantum_native.py](quantum_native.py) executes the joined instrument. Its
`output_on_joint_input` interface streams every unnormalized quantum branch
for a supplied native/work/reference input. `readout_outputs` additionally
streams the complete terminal code measurement for native/reference inputs
with the declared work preparation. These are different input contracts.

Restricting C to the uniform diagonal ensemble of the 66 binary native inputs
recovers the earlier classical benchmark, including its final port-zero load,
second-smooth plus probability and clock-bit marginal. The regression computes
these marginals from squared native-map columns and the original unnormalized
work vectors; it does not renormalize the surviving branches.

## Native channel and executed seam routing

Let `X_2` be all twelve-port nonnegative integer configurations of total load
two. Then `|X_2|=78`. For oriented seam `e=(u,v)`, coin `c`, and input `x`, put

\[
 m=x_u+x_v,\quad l=\lfloor m/2\rfloor,\quad h=\lceil m/2\rceil.
\]

The endpoint repair `F_c^e(x)` replaces the two endpoints by `(l,h)` for
`c=0`, and by `(h,l)` for `c=1`, leaving all other ports fixed. The transferred
current is `j=x_u-F_c^e(x)_u`; the original contrast is `delta=x_u-x_v`.
C records `(c,j)`, whereas E records `(c,delta)`. In either case the current
is available to the controller: for E,

\[
 p=|\delta|\bmod2,\qquad j=(\delta-(2c-1)p)/2.
\]

For candidate `A` and label `h=(c,ell)`, its native Kraus operator is

\[
 L^{A,e}_h=\frac1{\sqrt2}
   \sum_{x:\,\ell_A(x,c)=\ell}|F_c^e(x)\rangle\langle x|.
\]

The factor `1/sqrt(2)` already includes the fresh fair coin. For a fixed C
label the current specifies the inverse translation; for a fixed E label the
contrast and endpoint total determine the original endpoints. Thus each
summand family is injective, `sum_h L_h^dagger L_h=I_78`, and the operators
preserve quantum coherence precisely within their declared record groups.
These are the complete instruments from the
[source discriminator](../observer/source_discrimination/DERIVATION.md).

The existing nine-bit encoding `J:C^78 -> C^512` is unchanged. Let `pi` swap
physical port labels 1 and 5. Its source-register extension is

\[
 U_\pi J|x\rangle=J|\pi x\rangle,
 \qquad U_\pi|a\rangle=|a\rangle
   \quad(a\text{ outside the native code}).
\]

It is an involution with phase `+1` on every moved basis vector. Consequently
the source construction for seam `(0,5)` is forward routing, the existing
`(0,1)` repair, and inverse routing:

\[
 \widetilde L^{A,(0,5)}_h
       =U_\pi^\dagger\widetilde L^{A,(0,1)}_h U_\pi,
 \qquad \widetilde L^{A,e}_h J=J L^{A,e}_h.
\]

The orientation is preserved: port zero remains the first endpoint, so the
current decoder needs no sign change. [native_routing.py](native_routing.py)
implements the port permutation by disjoint encoded transpositions. A Gray
path realizes each transposition by controlled bit flips, undoing its
intermediate moves. The controlled-flip source word retains its scalar phase.
One full nine-bit macro matrix is obtained by executing its source word on all
512 columns; fixed wire addressing and actually executed source-X conjugations
reuse this matrix for the other control patterns. Inverse source words are
executed independently. These contractions avoid repeating host work; they do
not remove any repeated pulse from the source schedule.

The numerical calculation keeps all 512 native amplitudes and every raw
record code, including abort labels. An invalid first record has no decoded
current and ends in an explicit abort block. It is not assigned current zero
or conditionally normalized. Valid records can still contain amplitudes in
the 434-dimensional code complement; those amplitudes remain in the state.
The ideal native intertwining concerns code inputs only. It does not assign a
native load interpretation to arbitrary padding inputs. Whole-register
unitarity guards and residuals are numerical diagnostics, not directed
roundoff bounds or physical error-correction services.

## Complete joint Kraus law

Write `r=(r0,r1,r2,r3)` for the four common-program records. Bit `r2` is the
second smooth read; zero denotes its plus effect. The last native seam is

\[
 e(r_2)=\begin{cases}(0,5),&r_2=0,\\(0,1),&r_2=1.\end{cases}
\]

Equivalently the native actuator receives `1-r2`. The final gauge bit `r3`
is retained as the model-clock readout and does not select the native seam.
For first decoded current `j`, let `K_r^j` be the complete work Kraus operator
from `protocol.Program`, including its current-dependent initial gauge kick.
Preparation supplies the input to these operators; it is not repeated inside
them. The ideal joint branch operator is

\[
 \boxed{M_{h,r,s}=
   (L^{A,e(r_2)}_s L^{A,(0,1)}_h)\otimes K_r^{j(h)}}.
\]

Both native fair coins are already present in their `L` factors. No additional
coin probabilities multiply this expression. The output is the full CQ state

\[
 \Phi(\rho)=\bigoplus_{h,r,s}
       M_{h,r,s}\rho M_{h,r,s}^\dagger.
\]

Every branch is unnormalized, including zero branches. In particular it is
generally incorrect to replace its native component by an integer transition,
to discard the work output after reading `r`, or to tensor independently
computed native and work *density matrices* for an initially correlated input.
Tensor products of the branch operators are valid on that correlated input.

Completeness follows in the causal order of the continuation. Summing the last
native record first gives `sum_s L_s^dagger L_s=I_78` for each `r2`. Summing
`r` next gives `sum_r (K_r^j)^dagger K_r^j=I_W` for each already fixed first
current. Finally `sum_h L_h^dagger L_h=I_78`. Therefore

\[
 \sum_{h,r,s}M_{h,r,s}^\dagger M_{h,r,s}=I_{78}\otimes I_W.
\]

This proves a CPTP recorded channel for arbitrary joint native/work inputs,
including entanglement. Applying the same calculation with `I_R` proves it
with an arbitrary finite reference; equivalently expand an input into native
and work matrix units with operator-valued reference coefficients. Mixed
inputs follow by linearity or purification. The receipt's coherent native and
declared vacuum preparations form a useful subfamily, not the proof of this all-input
identity. Binary64 evaluations of the exact operators are checked numerically
against independent endpoint and work-matrix constructions.

For a committed prefix `q` with operator `T_q`, every complete continuation
has branch operators `D_(q,t) T_q` with
`sum_t D_(q,t)^dagger D_(q,t)=I`. Hence

\[
 \sum_t T_q^\dagger D_{q,t}^\dagger D_{q,t}T_q=T_q^\dagger T_q.
\]

The old-prefix probabilities and the reference reduced state are preserved
after summing continuations. This does not assert preservation of the old
conditional quantum state: later measurements and feedback generally change
it. Nor does it restore coherences between different classical record labels.
The earlier [record recovery](../observer/RECORD_RECOVERY.md) result restores
history-projected inputs in its stated scope; its retained classical tape
does not supply the cross-record coherence lost by C or E.

## A surviving-coherence test for this continuation

Take `X={2,3}`, `Z={1,2}`, and `W={0,2}`, with source addresses 0, 1 and 2.
Prepare native state `(X+exp(i theta) Z)/sqrt(2)` and a common work state
`sigma`, independent of it. The default `sigma` is the declared prepared
work vacuum. Define

\[
 q=\sum_{r:r_2=0}\operatorname{Tr}
        (K_r^0\sigma(K_r^0)^\dagger),\qquad 0\le q\le1.
\]

Consider the final native rank-one effect
`Q_plus=|(X+Z)/sqrt(2)><(X+Z)/sqrt(2)|`, forgetting the committed labels and
the work output only for this probability. The executable readout applies the
source Y rotation taking this vector to code zero, then measures all nine code
bits and retains every one of the 512 outcomes. Code-zero probability is
exactly this effect. This full basis measurement disturbs the complementary
subspace more than the older two-outcome nondemolition pointer; no equality
of their postmeasurement instruments is asserted. Raw-record aborts halt
before this read and remain explicit outcomes. Their ideal native-code weight
is zero; numerical abort weight is retained separately.
The Y basis change can produce measured code addresses outside the native
encoding even for ideal native inputs. These are ordinary outcomes of the
declared basis read, not an extra load interpretation, an abort, or evidence
of leakage before readout.

The first C coin-zero/current-zero branch preserves the `X,Z` off-diagonal,
with factor `1/2`; coin one separates it. E records different contrasts for
`X` and `Z` under both coins, so both E branches erase that off-diagonal.
On last seam `(0,5)`, both C coins fix `X,Z` with current zero. On last seam
`(0,1)`, only C coin zero keeps them in one current group; coin one takes
`Z` to `W` and separates their records. Thus, for the same candidate at both
native events,

\[
 p_C(\theta)=B+\frac{1+q}{8}\cos\theta,
 \qquad p_E(\theta)=B,
\]

where `B` is their common diagonal-input probability. The equality of the
diagonal baseline follows because C and E have the same basis transition and
the work controller uses the decoded current in either case. In particular

\[
 p_C(0)-p_C(\pi)=\frac{1+q}{4},\quad
 p_E(0)-p_E(\pi)=0,
 \quad p_C(0)-p_E(0)=\frac{1+q}{8}.
\]

The gap reverses at `theta=pi`; phase `pi/2` and the incoherent equal mixture
are controls with zero C/E gap. Quantum record frequencies alone do not reveal
this phase: the retained native interference effect is essential. This
distinguishes the two specified standard quantum instruments; it is not an
empirical observation or a disagreement with standard quantum mechanics.
The source execution retains its complete code-outcome table. Summing it to
the single plus probability is a reporting marginal, not postselection.

## Temporal and occupation comparison on the joint history

At fixed `d`, a work evolution replacement with complete-channel
half-diamond bound `epsilon_i` has that same bound after adjoining the native
register and a reference. The already committed native record selects one of
three currents. A controlled direct sum of those work channels has error at
most `max_j epsilon_i(j)`, not the sum over first-native outcomes. All other
ideal native and work stages are CPTP. Contractivity and telescoping therefore
give the canonical temporal bound

\[
 \epsilon_{\rm time}\le
  \min\left(1,\sum_{i:\,\mathrm{evolve}}\max_j\epsilon_i(j)\right).
\]

The existing palindrome estimate is uniform over these currents, so the
[continuation budget](continuation_bounds.py) already supplies this sum.
Exact ideal native repair and routing add no product-formula error. Their
floating execution discrepancies and physical implementation errors are
separate obligations; they have not been proved zero by this observation.

For the occupation comparison, assume the declared scalar vacuum family and
coupling `1/4`. Allow any normalized native/gauge/reference state tensor that
scalar vacuum; in particular the default blank gauge preparation and arbitrary
native/reference input satisfy this premise. Native repairs, routing, record
commitment and current-controlled gauge kicks commute with oscillator number
`N` and its Galerkin projector `Pi_d`. On the complete record dilation they
are isometries. They therefore preserve the summed fourth-number moment and
the projected error norm used in the [stagewise proof](COMPOSITION.md).
That proof bounds gauge operators by their operator norms and is unchanged
when native and reference factors are present. Conditional gauge controls
after either smooth read obey the same commutation identities.

Consequently the existing stagewise recurrence applies to the joint history:
the first native event is moment-neutral; both smooth reads raise the carried
moment before later Hamiltonian intervals; the final native event is again
moment-neutral. The preparation tail is included once. The projected Duhamel
error and final outside tail are combined once at the end, before the
complete-output trace-distance bound is taken. There is no division by a
native or work branch probability and no multiplication of the bound by the
number of records. This argument uses the common-domain and Hamiltonian
premises of the earlier occupation theorem; it does not prove them anew.

With `epsilon_occ` from that recurrence, the complete ideal joint recorded
output obeys

\[
 \tfrac12\|\Phi_{d,n}(\rho_d)-\Phi_\infty(\rho_\infty)\|_1
 \le\min(1,\epsilon_{\rm occ}+\epsilon_{\rm time}).
\]

The states on the two sides are the matched vacuum preparations, with the
finite state embedded in the common oscillator space and identical finite
native/gauge/reference data. Any bounded final effect, including `Q_plus`,
inherits this probability bound by contractivity. For each fixed finite
program the analytic family with `n(d)=d^2` has error `O(1/d)`; actual source
execution remains at `d=2,4`. This is a state-family occupation estimate, not
an arbitrary-scalar-input diamond bound. A supplied scalar state generally
needs a new preparation-tail and moment certificate.

## Supplied services and retained scope

The construction supplies finite source words, a complete recorded matrix
channel, and its ideal temporal/occupation comparison. Its premise ledger
still includes the native encoding, coherent proper-code load/export,
source pulse precision, blank coin/record/pointer preparation, classical
record commitment and persistence, conditional routing, and quantum work
storage. An input supplied as amplitudes is a simulator input, not a method
for reading an unknown quantum state. The existing work-history checkpoint
represents retained simulator state data; it supplies neither physical memory
nor a cloning operation, and is not a serialization API for this joined
native/work history.

All repeated source invocations belong in the realized-path resource counts,
including both directions of a selected route and the fresh native coins.
Expanded-tree enumeration and cached host matrix contractions are separate
resources. Invalid-label abort blocks retain their probability. A complete
physical error bound requires bounds for every service in the joined
construction, including the surrounding native stages; absent bounds remain
absent. The Hamiltonian durations remain model time. Neither the new routing
nor the last gauge bit supplies physical clock calibration.

## Executed results and reproduction

The receipt executes both native laws at `d=2,4`, with 512 temporal steps
per work interval. It keeps 770 pre-read branches for C and 2,566 for E;
these counts include first-record halts and all last-record aborts. The full
terminal basis read has 295,106 and 820,166 outcomes respectively, including
the unmeasured abort blocks. Factored matrices and streamed outputs avoid
storing a redundant dense history tape.

| Onsite dimension | C phase contrast | E phase contrast | Maximum C native pulses per attempted history |
| --- | --- | --- | --- |
| 2 | 0.349153278 | numerically below 2e-14 | 2,555,735 |
| 4 | 0.351351028 | numerically below 2e-14 | 5,961,279 |

The exact ideal E contrast is zero; the displayed nonzero residual is
floating execution error. C's exact contrast lies in `[1/4,1/2]` by the
formula above, independently of the particular value of `q`. E requires six
additional pulses per complete history. Each pulse has one proper-code load
and one export; the processor, work and native storage, blank supplies,
readout registers and record widths are itemized separately in the receipt.
Abort paths halt earlier, so the listed maxima cover them without charging
their unexecuted continuations as actual events.

```sh
python -m model.common_source.quantum_native --check model/common_source/quantum_native_results.json
python -m pytest -q -W error model/common_source/test_native_routing.py model/common_source/test_quantum_native.py model/common_source/test_quantum_native_theorem.py model/common_source/test_quantum_native_adversarial.py
```

The verifier compares every factored branch map with independent endpoint
averaging and product work maps, including prepared work outputs separately.
Complex native/work/reference tests check full output amplitudes, reference
marginals and old-prefix effects. The uniform binary-input restriction
recovers all five retained classical benchmark statistics. Negative controls
reject a reversed raw-current decoder, missing routing, a native phase error
invisible to basis populations, omitted aborts, changed work outputs,
discarded padding and nonfinite accumulated probabilities. Shared decoder
logic and an overflow-to-NaN normalization escape were corrected during this
audit; paired regressions reproduce the defective checks and their repairs.

A further independent audit found the same NaN comparison failure at both
input boundaries: finite but grossly unnormalized complex amplitudes could
overflow their squared norm. Input validation explicitly requires a finite
unit norm before either native/reference or joint-input execution; valid
complex states remain admitted. Resource maxima are checked against complete
record paths with their selected final seams, including the charged
zero-angle feedback requests. The receipt records complete quantum-output
retention and explicitly excludes preservation of cross-record coherence.
