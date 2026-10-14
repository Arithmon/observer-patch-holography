# Active encoded memory in the continuing observer

The [one-event native code](../observer/RECORD_RECOVERY.md) already satisfies
the Knill--Laflamme condition for every native seam, for both recorded laws C
and E. This construction supplies its missing source encoder and conditional
recovery, and inserts that recovery into the continuing native/work observer.
An unknown logical qubit stays on the same native register, including its
entanglement with a reference. The classical repair records are retained.

This is an active correction protocol. It adds specified, record-controlled
source operations between native events. It neither reverses the loss of an
arbitrary unencoded native state nor asserts that passive C/E dynamics protect
this code indefinitely. The earlier inverse repair restored history-projected
inputs; that operation alone is not the logical recovery constructed here.

## Code, branch probabilities and explicit recovery

Use the existing 78-dimensional total-load-two carrier and its nine-qubit
source encoding. Put

\[
 V|0\rangle=(|\{0,2\}\rangle+|\{1,3\}\rangle)/\sqrt2,
 \qquad
 V|1\rangle=(|\{0,3\}\rangle+|\{1,2\}\rangle)/\sqrt2.
\]

The four source addresses in this order are `2,73,74,1`. We use `V` also for
the isometry into the full 512-dimensional source register. Let
`E|l>=|l>` at source addresses zero and one, with the other eight code bits
blank. This is an input interface for an unknown logical qubit, not a
preparation of that unknown state from classical data.

For oriented seam `e`, fair coin `c` and raw repair record `a`, let
`L_h`, `h=(e,c,a)`, be the native Kraus operator, including its coin amplitude
`1/sqrt(2)`. As in the existing theorem, C groups by current and E groups by
contrast. Apply the independent endpoint floor/ceiling rule to each of the
two configurations in each logical codeword. For an admitted record, let
`S_(h,l)` be the resulting source addresses for logical input `l`.
The existing code has these properties for every one of the 30 native edges:

1. The two sets have equal cardinality `k_h` in `{0,1,2}`.
2. Each set has distinct elements and the two sets are disjoint.
3. All nonzero amplitudes have the same positive phase.

The four original occupied pairs are nonedges. For a binary native event its
effect is a linear combination of `I,N_u,N_v,N_u N_v`. Both codewords have
equal one-port occupation means and zero occupation product on a native edge;
their disjoint supports make the off-diagonal matrix element zero. Equivalently
the finite endpoint enumeration above proves

\[
 V^\dagger L_h^\dagger L_h V=p_h I_2,
 \qquad p_h=k_h/4,\qquad \sum_{c,a}p_{(e,c,a)}=1.
\]

The denominator four includes both the codeword's two-term normalization and
the fair coin. It does not include a seam-selection probability. If a policy
randomizes seams, its separately declared probability multiplies `p_h`.
This is the existing flagged Knill--Laflamme condition, with an explicit
decoder below; it is not a new general error-correction criterion.

Write `B_1=I` and let `B_2` be the source `Y(-pi/4)` pulse on code bit one.
The compiler convention is `exp(+i angle Pauli)`, so
`B_2 E|l>=(|l>+|l+2>)/sqrt(2)`. For `k_h>0`, choose a permutation `P_h` of
all 512 source addresses such that

\[
 P_h|l+2s\rangle=|x_{h,l,s}\rangle,
 \qquad x_{h,l,s}\in S_{h,l},\quad 0\le s<k_h.
\]

Complete this partial assignment to a permutation. One constructive procedure
starts with identity and fixes each prescribed image by swapping its current
image with the desired address. Distinct desired images ensure earlier fixed
assignments stay fixed. No target state is substituted into execution: these
swaps are compiled to Gray-path controlled flips, using the source macro and
phase convention of [native routing](native_routing.py).
The [current Walsh/Gray synthesis](QUANTUM_NATIVE.md#exact-gray-parity-synthesis)
reduces each nine-bit flip from 11,782 to 2,044 source pulses. It preserves
the entire chosen 512-dimensional permutation, including its complement
extension; it does not replace that extension by a different encoder that
only agrees on the two logical input columns.

Define `U_h=P_h B_(k_h)`. Likewise define `U_V=P_V B_2`, where
`P_V` maps `0,2,1,3` to `2,73,74,1`, respectively. Then

\[
 U_V E=V,\qquad L_h V=\sqrt{p_h}\,U_h E,
 \qquad R_h=U_V U_h^\dagger,
 \qquad \boxed{R_h L_h V=\sqrt{p_h}\,V}.
\]

Chronologically, recovery executes inverse `P_h`, inverse `B_(k_h)`, `B_2`,
then `P_V`. It is a full source-register unitary selected by an already
committed classical record. It uses neither the unknown logical amplitudes
nor a state-dependent normalization. No quantum scratch register or reset of
the native state is required for this operation. Valid records with `p_h=0`
receive a specified full unitary, such as identity; they remain represented.
Invalid raw records retain explicit abort blocks rather than invoking a
decoder with an invented current.

[encoded_memory.py](encoded_memory.py) executes the encoder, its inverse and
the recoveries through source words. The executed seam menu is `(0,1)`,
`(0,5)` and `(2,3)`; the all-30-edge algebraic code theorem does not assert
that every other seam's source routing has also been implemented. Source
matrices retain all 512 rows, including numerical code-complement amplitudes.
The completed recovery unitaries define a channel off the code, but exact
logical correction is asserted on `V(C^2)` only.

## Continuing correction and retained records

Consider any fixed finite number of native events. At each stage a classical
policy may choose an implemented seam from earlier committed records and
independent declared randomness. Immediately after a non-abort native event,
apply its corresponding `R_h`. The policy never accesses the unknown logical
state. Repeated use of the boxed identity gives

\[
 \left(\prod_i\sqrt{\pi_i(e_i\mid h_1,\ldots,h_{i-1})}\right)
 R_{h_m}L_{h_m}\cdots R_{h_1}L_{h_1}V
       =\sqrt{w_{\boldsymbol h}}V,
 \qquad w_{\boldsymbol h}=\prod_i
       \pi_i(e_i\mid h_1,\ldots,h_{i-1})p_{h_i}.
\]

Here `pi_i` is the normalized classical distribution for the next seam,
and equals one for a deterministic selected seam. Its square root is part
of the complete branch Kraus operator; it is not already included in the
native operator `L_h`. The selected seams, coins and raw records are all
part of the history label.

Thus, for any logical/reference density operator `rho_LR`, the complete ideal
recorded output is

\[
 \bigoplus_{\boldsymbol h}w_{\boldsymbol h}
       (V\otimes I_R)\rho_{LR}(V^\dagger\otimes I_R),
 \qquad\sum_{\boldsymbol h}w_{\boldsymbol h}=1.
\]

The tape probabilities are independent of the logical input. Recovery fixes
the old classical record projectors because it acts blockwise without changing
their labels. The logical information is restored in each ideal branch, and
tracing out the complete tape preserves the unknown qubit and its reference.
The statement is for every fixed finite protocol, with its explicit record
storage and operation count. It supplies no finite memory for an unbounded
tape or autonomous operation without fresh coin/record services.

## An exact passive failure control

Use seam `(0,1)` followed by seam `(2,3)`. Both are native edges. In the code,
one particle occupies one of ports `0,1` and the other occupies one of ports
`2,3`. The logical basis is their equal-versus-opposite pairing. The first
event's coin and record identify the first original choice, and the second
event's coin and record identify the second. This holds for both C and E.

Without correction between the two events, each of the 16 admitted histories
is compatible with exactly one of the two logical basis states. Hence the
complete passive recorded channel `Gamma_pass` satisfies

\[
 \Gamma_{\rm pass}(V|0\rangle\langle1|V^\dagger)=0,
 \qquad
 \Gamma_{\rm pass}(V|+\rangle\langle+|V^\dagger)
 =\Gamma_{\rm pass}(V|-\rangle\langle-|V^\dagger).
\]

Their input trace distance is one, so every later CPTP decoder has trace
distance error at least `1/2` on one of these phase inputs. The same channel
is entanglement breaking on this logical code; its best entanglement fidelity
for half of a Bell pair is `1/2`. This value is attainable by decoding the
classical logical basis value, so it is an optimum, not only an obstruction
for one chosen inverse circuit.

With correction after each event, each event instead has four nonzero records
of probability `1/4`. All 16 complete histories have operator `V/4` on the
logical input. The corrected channel is the identity on the logical qubit,
with Bell entanglement fidelity one. An executed inverse encoder followed by
the source logical-plus basis read gives probabilities one and zero on the
two corrected phase inputs. All other code/syndrome outcomes remain present.

The distinction is substantive: the older `(0,1)` then `(0,5)` schedule does
not expose this code's logical pairing, and remains passively correctable on
the code. It is not used as evidence that intermediate correction is needed.
The fixed `(0,1),(2,3)` challenge is separate from a work-controlled choice of
the second seam, which generally only partially erases the passive phase.

## The protected native/work observer

The protected mode first repairs seam `(0,1)` and applies its recovery. Its
decoded current controls the initial gauge kick on the same prepared work
register. The [canonical four-read program](COMPOSITION.md) then runs with
both smooth reads, record feedback and its final gauge clock probe. Smooth
bit `r2=0` chooses last seam `(2,3)`; `r2=1` chooses `(0,1)`. The last native
event is again followed by its recorded recovery. Neither the logical/native
register nor the quantum work register is prepared again during that history.

Let `K_r^j` be the complete work Kraus map, including its current kick, and
let `h,s` be the first and last native records. For an arbitrary logical/work/
reference input, the corrected branch satisfies

\[
 M_{h,r,s}(V\otimes I_W)
   =\sqrt{p_h p_s^{e(r_2)}}\,(V\otimes K_r^{j(h)}).
\]

This is an operator identity, not only equality for a classical mixture or a
few prepared phases. It proves that the native memory adds no logical noise
to the intended work instrument, including when logical and work inputs are
entangled. A work measurement may conditionally steer an entangled logical
state; the identity does not assert that every such conditional logical state
is fixed. After summing all complete outputs, the intended local work channel
preserves the logical/reference marginal. The canonical executed preparation
uses the declared scalar vacuum; supplied arbitrary work inputs retain the
finite intertwiner but not that vacuum's occupation certificate.

There is a stronger separation for the memory marginal. Replace the work
apparatus by any complete CPTP instrument acting only on the work register,
with classical controls still selecting an implemented seam and its correct
native recovery. Each complete branch remains identity on the logical input
tensor a work Kraus map. Summing those maps proves that tracing out work and
records leaves the logical/reference marginal **exactly unchanged**, even for
an initially entangled logical/work/reference input. Accurate Hamiltonian
evolution, a particular work preparation, and an occupation approximation are
unnecessary for this memory identity. They remain necessary for the claimed
full observer output. Incorrect native records, the wrong recovery, or faults
coupling the apparatus to logical information are not corrected by this claim.

For an independently prepared logical phase input and work vacuum, every
corrected branch has the same logical phase. A final source inverse encoder
and logical-plus read therefore retain phase contrast one. Reading only the
raw native source bit without decoding would be a different observable.
The full 512-outcome read retains every syndrome and code-complement outcome;
it does not project onto an assumed blank syndrome.
This optional terminal measurement consumes the logical coherence being
tested. The memory identity concerns the state before that destructive read;
it is not a nondestructive measurement of an unknown logical phase.

For the corresponding passive mode, define `q_j` as the second smooth plus
probability of the work program at first current `j`. The code's first-event
current probabilities are `1/4,1/2,1/4` for `j=-1,0,1`, so put

\[
 \bar q=(q_{-1}+2q_0+q_1)/4.
\]

For independent initial logical phase states and the same work preparation,
the complete passive outputs for logical plus/minus have trace distance
`1-qbar`. Indeed the last `(2,3)` branch records the logical pairing and
erases its phase; the last `(0,1)` branch preserves the phase in orthogonal
native outputs. All records are retained, so their weighted trace distances
add. Thus every later decoder has error at least `qbar/2` on one phase input.
The corrected outputs instead have distance one. This full-output comparison
does not assume that the passive terminal inverse-encoder read is the optimal
decoder or an optimal phase measurement.

There are also exact, strictly positive finite-cutoff obstructions independent
of the numerical values of `q_j`. The second smooth effect is
`E_plus=(I+sin(Q0,d))/2`, where `Q0,d` is the usual truncated
`(a+a^dagger)/sqrt(2)` coordinate tensored with the other work factors.

At `d=2` its site eigenvalues are `+/-1/sqrt(2)`. The elementary inequality
`|sin x|<=|x|` and `1/2<9/16` give `E_plus >= I/8`. Thus every normalized
incoming work state has `qbar>=1/8`, and the passive worst-case decoder error
is at least `1/16`.

At `d=4` the positive site eigenvalues are
`lambda_+/-=sqrt((3+/-sqrt(6))/2)`, from the characteristic polynomial
`lambda^4-3lambda^2+3/4`. The rational inequalities

\[
 (12/5)^2<6,\quad 27/10>(41/25)^2,\quad
 \lambda_+>41/25,\quad \pi/2<11/7,
 \quad 12/175>1/16
\]

put `lambda_+` more than `1/16` above `pi/2`. Since `sqrt(6)>1`,
`lambda_-<1`, more than `1/2` below `pi/2>3/2`. Both eigenvalues are
between zero and two, hence below `pi`. Consequently their sine magnitudes
are at most `cos(1/16)`. Alternating cosine Taylor terms give, for `x=1/16`,

\[
 \cos x\le1-x^2/2+x^4/24\le1-x^2/3=1-1/768.
\]

The second inequality follows from `x^2<=4`. Therefore `E_plus>=I/1536`,
and every decoder has error at least `1/3072` on one passive phase input.
These constants use only the finite detector and complete trace-preserving
preceding work operations. They are not inferred from a rounded eigenvalue
or the displayed phase contrast. The infinite coordinate spectrum does not
have this gap from `pi/2`; no state-independent infinite-occupation version of
these lower bounds is asserted.

## Error, resource and physical scope

Ideal source encoding and correction have zero error on this declared code.
Their pulse realization is conditional on the existing coherent proper-code
services and exact ideal source controls. Floating operator comparisons are
diagnostics, not directed roundoff certificates. They must test full complex
branch maps, inverse phases and complete reference outputs, rather than only
population agreement or unitarity. All zero, abort and numerical padding
outputs remain in the complete instrument.

The native encoder, routing and recoveries act on finite factors and commute
with scalar oscillator number and the Galerkin projector. The existing
stagewise occupation comparison therefore applies unchanged to this fixed
protected work history under its scalar-vacuum and coupling premises. Both
smooth reads still raise the carried moment before later evolution; the
preparation tail and final outside component are each counted once. The
complete temporal bound remains the sum of work interval bounds, uniformly
over earlier currents. No factor equal to the number of native records is
introduced. The resulting bound concerns the complete logical/native/work/
record output and its reference, on the matched prepared scalar family.
For the logical/reference memory marginal alone, these work temporal and
occupation errors need not be added: the work-independent identity above
already compares that marginal exactly. Native encoding/correction/decoding
faults retain their own distinct error obligations.

For an additional conditional service budget, each complete implemented stage
must have a specified input/reference-uniform bound on its relevant domain.
If a corrected stage is only certified on the code, a valid telescoping proof
uses ideal code-valued prefixes and complete CPTP implemented suffixes: each
term then evaluates that stage's error on a code input before contractive
later processing. It must not assume that earlier imperfect corrections leave
the actual state exactly in the code. Missing service or numerical bounds
remain missing; gate counts alone do not supply them.
The separate `native_logical_storage_through_work` service covers the encoded
native register during the full intervening work program; the inherited work
and archive storage key does not silently certify that quantum memory.

Count the initial encoder once, every attempted native repair and its actual
record-selected recovery, both directions of each selected route, all work
operations, and any terminal decoder/read. Cached macro matrices reduce host
work, not the number of pulse, load or export requests on a realized path.
Separate the maximum realized-path cost from the expanded simulation tree.
Native storage comprises nine code qubits; one unknown logical input and
the other initial blanks are distinct resources. Classical archives grow with
the finite number of events. Recovery itself needs no additional logical
copy, no unknown-state reset, and no erasure of the existing record archive.
The receipt's three known logical phase settings additionally charge three
source preparation pulses. The supplied-unknown-input interface omits that
preparation; applying it to the unknown qubit would change the channel.

The current compiled schedule has the following exact source counts. The
previous column uses the original separately expanded Walsh factors; each
replacement preserves its full endpoint operator and retains all records.

| Operation | Previous native pulses | Current native pulses |
| --- | --- | --- |
| Initial encoder or terminal decoder | 164,967 | 28,635 |
| Maximum first `(0,1)` recovery | 259,231 | 44,995 |
| Maximum `(2,3)` recovery | 306,367 | 53,179 |
| Protected C history, `d=2`, 512 temporal steps | 4,936,767 | 1,158,423 |
| Protected C history, `d=4`, 512 temporal steps | 8,342,311 | 4,563,967 |

The encoder pays for fourteen macros, eighteen conjugating X pulses and
one mixing pulse: `14*2044+18+1=28635`. The full protected-path maxima
include source preparation, both repairs and recoveries, the continuing
work history, the selected route and terminal decode/read. E adds six pulses
to each displayed complete-history maximum. Every pulse still requires one
proper-code load and one export. Fewer requests introduce no additional
quantum blanks or physical services and establish no physical error bound.
Endpoint equality does not identify the intermediate states or control
histories of the two source words.

This construction implements a finite active memory service inside the
observer model. It does not establish fault-tolerant physical memory, a
threshold against unspecified noise, calibrated duration, the passive native
geometry/clock statistics under the added controls, or a full
`OperationalObserver`/consensus-tower realization. The old unencoded observer
remains a separate mode with its existing law and benchmark.

## Reproduction and independent controls

The mode is selected explicitly with
`quantum_native.execute(native_mode='protected_qubit')`; `passive_qubit`
uses the same encoding, work program and second-seam selector without active
correction. The default `full_carrier` mode retains the previous unencoded
78-dimensional input contract. The encoded modes accept a two-dimensional
logical input, including logical/work/reference arrays through the existing
streamed joint-output API.

[encoded_observer.py](encoded_observer.py) produces the complete channel
audits, source phase reads, memory comparisons, exact passive-obstruction
constants and realized-path resource schedules. Its full-output work budget
and unknown physical-service obligations remain separate from the exact
ideal logical-memory identity.

```sh
python -m model.common_source.encoded_observer --check model/common_source/encoded_observer_results.json
python -m pytest -q -W error model/common_source/test_encoded_memory.py model/common_source/test_encoded_memory_theorem.py
```

The independent theorem tests enumerate native endpoint effects on all 30
edges, compare executed corrected branches with their exact logical target,
retain complex reference amplitudes through two events, and distinguish the
crossed-seam passive failure from the older passively correctable schedule.
The integrated tests compare complete output amplitudes on entangled logical/
work/reference inputs and check both laws' phase-distance formula. A separate
matrix-instrument control replaces the work program by a complete
measure-and-erase channel: work changes, while the logical/reference marginal
still agrees with the memory theorem. This control is not a claim that work
reset was added to the canonical source-program grammar.
