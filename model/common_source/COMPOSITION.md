# The finite executable observer specification

The common source supports one fixed finite sequence of Hamiltonian evolution,
gauge reads, smooth Ramsey reads and feedback controlled by earlier committed
bits. Each operation acts on the quantum outputs of the preceding operation.
The complete branch ensemble is retained throughout. A later gauge read after
another Hamiltonian interval is therefore a model-time readout on that same
history, including the disturbance from its earlier smooth measurements.

This is the canonical finite contract for `protocol.Program`: the same engine
executes the retained sharp and smooth benchmarks as well as the continuation
below. The construction guide defines the broader physical completion target.
The [existing premise map](../PREMISES.md) remains the complete theorem-family
navigation; the instantiated connections here have the narrower scopes below.

The [numerical output certificate](OUTPUT_CERTIFICATE.md) connects the actual
returned source amplitudes to the exact prepared continuous and infinite
histories. It treats their binary components as exact data, retains all
records and bounds their complete quantum outputs without renormalization.
It also certifies a fixed feedback intervention's final clock probability.
This is a prepared-output comparison, not an all-input source-channel or
physical service certificate; the operation grammar is unchanged.

The [mixed-state extension](OPEN_HISTORY.md) uses this same program grammar
when unobserved finite environment operations are inserted at declared
boundaries. It retains one density block per visible record, composes the
existing noise and memory results, and preserves the preparation, feedback,
clock and continuation contracts below. It introduces no second operation
language or physical noise calibration.

The new quantitative occupation result covers this entire finite sequence.
It follows a fourth-number moment through each smooth pointer and uses that
larger moment in every later evolution estimate. It does not restart a fresh
experiment at a checkpoint. The canonical sequence is

```text
prepare once
evolve(1), gauge_read -> r0, feedback(r0), evolve(1/2)
smooth_read -> r1, evolve(1/2), smooth_read -> r2, feedback(r2)
evolve(1/2), gauge_read -> r3
```

Thus its Hamiltonian model time is `5/2`, and it uses two smooth pointers and
four binary records. The pointer interaction parameter is separately `s=1`
for each smooth read; these two units are not added to Hamiltonian model time.
The protocol also accepts other fixed finite sequences in the same grammar.
Every feedback reference must name an already committed bit. Preparation is
performed once; resets and branch postselection are outside this grammar.
The runtime has explicit budgets of 64 operations and eight reads; the
canonical program retains 16 branches. The mathematical estimate permits
any fixed finite sequence and does not assert execution beyond those budgets.

The [runtime](protocol.py) supplies the finite source execution. The directed
certificate is [continuation_bounds.py](continuation_bounds.py). The older
sharp detector can remain an alternative executable protocol, but a sharp read
is rejected by this occupation certificate. The moment premise needed for
continued evolution after that read has not been proved by replacing it with
a different smooth instrument.

The separate [prepared-history certificate](SMOOTH_HISTORY.md) compares this
canonical 16-output history at d=2,4 with its infinite counterpart, for all
three native currents. It retains later feedback, evolution and the gauge
probe, with the blank-gauge and scalar-vacuum preparation fixed. This
state-specific comparison supplements the uniform moment bound, which is
one at these cutoffs; it does not cover arbitrary initial work states or
certify floating execution and physical services.

## Execute the complete observer

The [coherent native attachment](QUANTUM_NATIVE.md) surrounds this same
four-read program with source-compiled C or E repairs on the full 78-state
native carrier. It retains native/reference coherence within each committed
history, uses the first decoded current for the gauge kick, and routes the
last repair from the second smooth bit. Its complete source-code outputs
include padding and abort branches. The construction and its separate
[receipt](quantum_native_results.json) extend the native interface; they do
not replace the classical benchmark below.

### Retained classical native benchmark

```sh
python -m model.common_source.composition --check model/common_source/composition_results.json
```

The first native E-seam event produces current `j`; that current selects the
gauge preparation kick on the same work register. The quantum program retains
all sixteen records `(a,b,c,z)`. Its final native actuator reads only `1-c`:
smooth bit zero is the plus effect, whereas the native actuator uses one for
plus. The last gauge bit `z` is the clock readout. Enumerating all 66 native
binary inputs, both first coins, sixteen quantum records and both last coins
gives 4,224 outcomes per executed cutoff, including zero-weight outcomes.
Every quantum terminal state remains attached to its native outcomes during
enumeration. The receipt stores summaries; replay reconstructs all outcomes
from the sixteen retained quantum branches for each current.

With `P` the gauge pointer, `L` the smooth pointer and `F` the committed-bit
feedback, the full input operator for the quantum record is

\[
 K_{abcz}=P_z U_\tau F_c L_c U_{1/2} L_b U_{1/2} F_a P_a U_1 K_j.
\]

Here `K_j` is the current kick, and each `U` is the declared source palindrome
in execution or the same continuous Hamiltonian in the comparison. No operator
prepares the vacuum again between these reads. An independent reset control
does precisely that after the first smooth read while retaining old record
probabilities. Its terminal-state distance from the actual continuous history
exceeds 0.89 for all six executed cutoff/current cases; this is a numerical
model comparison, not an empirical effect or an outward distance enclosure.

At 512 product steps per interval, the rational ideal temporal bounds for the
complete recorded output are below 0.000999668 at d=2 and 0.005752145 at d=4.
Source executions are separately checked against independent product maps and
continuous histories. Physical services and floating errors are not folded
into those ideal bounds. The 4096-cutoff analytic family with `n=d^2` has a
combined temporal/occupation bound below 0.056279; that large apparatus is not
executed. The retained receipt stores each cutoff's shared bound once.

There is one quantum work preparation and four fresh pointer blanks per
realized history. Source word counts are 352,353 native pulses at d=2 and
3,757,897 at d=4 for this 512-step protocol, plus its two declared native
repair events. The public tape fits ten fixed-width bits: three for each
native coin/current letter and four quantum bits. This excludes work-state
storage, program tables, host workspace and physical initialization costs.
The expanded tree uses fifteen pointer invocations; that simulation total
is different from the four invocations in a single realized path.

`History` checkpoints contain the retained work/reference amplitudes and
prefix operators as simulator data. JSON restoration checks the program,
prefix, all outcomes and costs by reexecution before continuing. It supplies
neither a physical quantum-state reader nor a method to reconstruct unknown
work coherence from classical records. The complete prepared work state is
unchanged by a pause; the controller does not receive its hidden amplitudes.
Direct in-memory continuation trusts runtime-produced `History` objects and
checks their integrity and consistency; their public digest is not source
authentication. Untrusted saved data must pass the JSON replay boundary.
The program's `preparation` field describes its declared default.
An explicit `execute(initial=...)` supplies a different work/reference input;
the history records `preparation_mode='supplied_input'` and that input's
amplitudes. This override retains finite channel composition but does not
inherit the vacuum occupation bound.

## Existing finite results connected by this specification

| Result family | Actual instantiation and retained boundary |
| --- | --- |
| CAR, Gauss constraints and common action | The existing `gauge_matter` four-state invariant sector and `scalar_coupling` action supply the same `W_d` for every operation. Species, action and coupling are supplied. |
| Complete instruments and composition | Each permitted step supplies its whole Kraus family. The induction below instantiates the finite Kraus/CPTP composition identities; tests compare every branch map on an entangled reference. |
| Public records and controlled repair | The matrix prefix identity below is the quantum analogue of the finite real-state expectation-dual criterion in `Lean/Dynamics/ProtectedCharge.lean`; it preserves the classical old-prefix marginal. `records.py` retains the distinct public and quantum-block conditional expectations. A first write and postcommit preservation are different maps. |
| Preparation and source compilation | The existing proper-code compiler executes the projected vacuum, native-current kick, every evolution, complete pointer and feedback pulse. Preparation is counted once; coherent services remain explicit. |
| Observer restrictions and finite GNS | Every resulting CQ state defines the positive functional on the declared record/work algebra. The GNS construction in `OBSERVER_READOUTS.md` applies to this functional; its old four-record numerical spectrum is not reused for sixteen records. |
| Continued observations and model clocks | Two smooth reads disturb the retained state before the last Hamiltonian interval and gauge probe. The supplied full-work Gibbs identity describes each continuous interval at `t=-beta*s`; physical calibration remains open. |
| Occupation and temporal approximation | The stagewise theorem below covers that entire same-state continuation. Approximate occupation inclusions are not exact consensus-tower morphisms. |
| Classical native observer feedback | The retained 4,224-outcome benchmark uses the existing integer primitive and the specified smooth-bit decoder. Its classical native input is distinct from an arbitrary quantum M78 input. |
| Coherent native observer interface | [QUANTUM_NATIVE.md](QUANTUM_NATIVE.md) composes source-compiled C/E repairs with this same work instrument on arbitrary native/work/reference inputs. Source-executed exchange of port labels 1 and 5 routes seam `(0,1)` to `(0,5)`; the final selector is `1-r2`. Complete raw outputs retain padding and aborts. Default executed preparations use the declared field vacuum; arbitrary scalar inputs do not inherit its occupation bound. Encoding, source services and physical readout remain supplied. |
| Active encoded native memory | [ENCODED_MEMORY.md](ENCODED_MEMORY.md) attaches source-compiled record-conditioned recovery to the existing four-configuration qubit. The `protected_qubit` mode corrects after seam `(0,1)` and after the last seam, `(2,3)` for `r2=0` or `(0,1)` for `r2=1`. The same work instrument continues between them. Its branch identity holds for arbitrary logical/work/reference inputs before terminal readout; every raw outcome remains. The matched `passive_qubit` mode omits recovery. This does not restore arbitrary full-carrier coherence or inherit the passive geometry/clock calibration. |
| Alternative instruments and stronger observer contracts | Sharp readout remains executable with its own certificates; its post-read moment bound is not inherited. Regional quantum partial traces are CP but generally not multiplicative. Exact tower naturality, complete federation access and calibrated physical readouts retain their separate premise obligations. |

## The finite record and instrument contract

At a checkpoint with `k` bits, let `R_k={0,1}^k` and let `W_d` be the common
working space of dimension `4d^2`. The accessible algebra after classical
record commitment is

\[
 \mathcal A_{k,d}=\bigoplus_{r\in R_k} B(W_d),\qquad
 \mathcal C_k=\{(f(r)I)_{r\in R_k}:f:R_k\to\mathbb C\}.
\]

The public record algebra is the central subalgebra `C_k`. A state is the
collection of positive, **unnormalized** blocks `rho_r`, with
`sum_r Tr(rho_r)=1`. For a work-only branch vector the block is
`rho_r=|v_r><v_r|`. For a work-by-reference amplitude matrix `V_r`, the work
block is `rho_r=V_r V_r^dagger`; the complete referenced block is
`|vec(V_r)><vec(V_r)|` on `W_d` tensored with the untouched reference. Thus
mixed work inputs are covered by a purification, and the channel definitions
also apply directly to mixed blocks. Zero branches remain represented.
Dividing each block by its trace would change the state and invalidate both
the instrument law and the error recurrence.

For a read, let `K_(r,b)` be its work-space Kraus operators, satisfying

\[
 \sum_b K_{r,b}^\dagger K_{r,b}=I\quad\text{for every }r.
\]

The complete Schrödinger map sends
`rho_r` to the block `K_(r,b) rho_r K_(r,b)^dagger` labelled `(r,b)`.
Its Heisenberg dual is

\[
 (\Phi^*A)_r=\sum_b K_{r,b}^\dagger A_{r,b}K_{r,b}.
\]

Kraus form proves complete positivity; the displayed completeness identity
proves trace preservation. Composition preserves these properties. For any
earlier prefix `p`, define `Q_p(r)=I` if `r` extends `p`, and zero otherwise.
Completeness then proves the exact prefix identity

\[
 \Phi^*(Q_p^{\rm after})=Q_p^{\rm before},\qquad
 \sum_{r\succeq p}\operatorname{Tr}\rho_r^{\rm after}
 =\sum_{r\succeq p}\operatorname{Tr}\rho_r^{\rm before}.
\]

An evolution or committed-bit gauge feedback is the same contract with one
outcome and a record-controlled unitary. Induction therefore preserves every
old prefix marginal through every continuation. In a variable tape space,
this is an equality between observables in the appropriate before/after
algebras, not a claim that matrices of different dimensions are identical.

The block pinching `X -> sum_r Q_r X Q_r` retains within-record quantum
information. The public readback
`X -> sum_r Tr(Q_r X Q_r) Q_r/dim(W_d)` keeps only public scalars. Both are
idempotent; they are different channels. Durable classical records do not
require the continuation to preserve the whole coherent matrix algebra on
the record register. That stronger requirement would obstruct meaningful
record-conditioned control, as discussed in the construction paper's
record-algebra obstruction. Nor does knowledge of the effects alone select
the measurement's post-read quantum output: the smooth Kraus maps, including
their ancillary construction, are part of the declared protocol.

These are concrete instances of the finite mathematical interfaces in
`Lean/Dynamics/ChoiCPTP.lean` (`isCompletelyPositive_of_kraus`,
`IsCPTP.comp`, `partitionPinching_isCPTP`, `partitionAverage_isCPTP`) and the
public-record algebra in `Lean/EventAlgebra/PublicRecordAlgebra.lean`.
The Python output is not thereby a Lean proof of its floating matrices.
The exact argument above, the specified instrument matrices and numerical
execution checks have different evidence roles. In particular,
`Tower.OperationalObserver` additionally asks for a consensus-tower state,
generator, refinement, prediction and access packet. A continued finite
instrument alone does not provide that entire packet.

## The specified occupation family

All occupation estimates here use the same supplied two-site Hamiltonian,
preparation and rectangular cutoff as
[OCCUPATION_REFINEMENT.md](OCCUPATION_REFINEMENT.md):

\[
 H=H_g+N+1+\frac\kappa2(Q_0-Q_1)^2
       +g(a_0+a_0^\dagger)Z_0,\quad
 \kappa=\frac9{32},\quad g=\frac14,\quad N=n_0+n_1.
\]

Here `N` is oscillator occupation, not the cosmological capacity parameter.
`Pi_d` keeps `0<=n0,n1<d`, `H_d=Pi_d H Pi_d`, and even `d>=2` is allowed
analytically. Products such as `Q_j^2` in `H` are formed before compression.
The finite preparation is the normalized projection of the specified free
coupled vacuum to total occupation `N<=d-2`, embedded in that rectangle.
The gauge preparation and declared native-current gauge kick commute with
`N` and `Pi_d`. The theorem is uniform over native currents `-1,0,1` and the
specified on/off gauge feedback. It does not certify other initial scalar
states or other couplings merely because their finite matrices can execute.

Self-adjointness on the common finite-number core and convergence of the
rectangular Galerkin Hamiltonians are inherited from that occupation proof.
The number-sector calculation there establishes, on finite and infinite
evolution, for `w=(N+2)^4`,

\[
 M=\langle w\rangle,\qquad |\dot M|\le CM,\qquad C=13/4,
 \qquad M_0=143748821/8192000.
\]

This bound is valid for vector-valued number sectors with an untouched
ancilla or record environment. Hence it applies to a dilation retaining all
records and to their summed, unnormalized output. It is not a uniform bound
on separately normalized, arbitrarily rare branches.

Each gauge-read Kraus family commutes with `N` and `Pi_d` and is complete.
Therefore `sum_a P_a^dagger w P_a=w`. Each permitted conditional gauge
unitary also commutes with `w` and `Pi_d`. These operations preserve the
summed moment and introduce no occupation-comparison error. Conditioning the
choice of a gauge unitary on an earlier smooth bit changes neither argument:
the unitary acts blockwise on the existing record, and each block obeys the
same commutation and completeness identities. An arbitrary scalar feedback
operation would need its own moment and leakage inequalities.

## Repeated smooth reads and projected-error induction

Every smooth read uses the Ramsey instrument from
[SMOOTH_REFINEMENT.md](SMOOTH_REFINEMENT.md). Its pointer has generator
`G=Z_anc Q0/2`, run for auxiliary parameter `0<=s<=1`, with fixed preparation,
basis changes and complete binary readout. Its effects are
`(I +/- sin(Q0))/2`, with bit zero the plus effect. The finite pointer uses
`Q0,d=Pi_d Q0 Pi_d`; its generator is precisely the compression of `G`.
Pointer preparations, basis changes and recording act trivially on scalar
occupation. The previously proved number-sector estimates give

\[
 |dM/ds|\le2M,\qquad
 \|\Pi_dG(1-\Pi_d)\Psi_s\|
 \le\frac{\sqrt{M(s)}}{2(d+2)^{3/2}}.
\]

These inequalities hold for any incoming finite fourth moment, including the
entangled and disturbed output of earlier smooth measurements. Thus every
finite number of such reads still has a finite summed fourth moment. No
claim about post-sharp-read fourth moments is needed or transferred.

Write `Psi` for the complete infinite dilation vector, with all pointer
environments retained, and `Psi_d` for its finite counterpart embedded into
the same space. During the induction carry two numbers:

\[
 \sqrt{\langle\Psi,w\Psi\rangle}\le A,\qquad
 \|\Psi_d-\Pi_d\Psi\|\le e.
\]

For `r=d/2`, define

\[
 p_d=\frac{81}{80}\frac{\binom{2r}{r}}{4^r81^r},\qquad
 A_0=\sqrt{M_0},\qquad e_0=\sqrt{2p_d}.
\]

The preparation error is at most `e0`; using that full vector bound also
bounds its projection. This is a conservative initialization, consistent
with the actual normalized total-occupation preparation.

Let

\[
 L_H(d)=\frac{9}{32(d+2)}+\frac{1}{4(d+2)^{3/2}},\quad
 L_S(d)=\frac{1}{2(d+2)^{3/2}},\quad a=13/8.
\]

The stage recurrences use the **incoming** value of `A`:

| Stage | New projected error | New moment amplitude |
|---|---|---|
| Hamiltonian interval `t>=0` | `e + A (exp(a t)-1) L_H(d)/a` | `A exp(a t)` |
| Smooth Ramsey read | `e + A (exp(1)-1) L_S(d)` | `A exp(1)` |
| Gauge read or permitted feedback | `e` | `A` |

To prove the Hamiltonian row, project the infinite Schrödinger equation.
Its source term is `Pi_d H (1-Pi_d) Psi`; the existing leakage estimate bounds
it by `A exp(a u) L_H(d)`. Variation of constants and finite-unitary norm
preservation bound the new projected error by the old error plus the
integral of this source term. The pointer row is exactly the same argument
with growth exponent one and leakage coefficient `L_S`. Scalar-neutral
recording and gauge control dilations commute with `Pi_d` and are isometries,
so they cannot increase the projected error. This establishes the induction
through the whole history, including feedback and later evolution.

Only at the final checkpoint add the infinite outside component. It is
orthogonal to `Psi_d-Pi_d Psi`, and its norm is at most `A/(d+2)^2`. Therefore

\[
 E_d=\sqrt{e^2+\frac{A^2}{(d+2)^4}},\qquad
 \|\Psi_d-\Psi\|\le E_d.
\]

The code adds no outside tail at intermediate checkpoints. Adding each
intermediate tail would be unnecessarily loose; dropping the final tail
would leave part of the output uncontrolled. The final fourth moment bound
is `A^2`. Before rounding it is `M0 exp((13/4)T+2m)` for total Hamiltonian
model time `T` and `m` smooth reads. The projected residual also depends on
their order, because earlier smooth reads increase the cost of later
Hamiltonian segments.

Tracing out the dilation environment or dephasing its records is a common
channel. Consequently `min(1,E_d)` bounds the trace distance of the complete
classical/quantum output, every bounded-effect probability error and the
total variation distance of its full classical history. Trace distance uses
`||rho-sigma||_1/2`. There is no factor equal to the number of branches.
The estimate also survives any common record-controlled native continuation
by channel contractivity. A finite mixture over the declared native currents
preserves it by convexity. Those observations require the same controller
and bit dictionary on both sides of the comparison.

All certificate arithmetic is rational. Square roots and exponentials are
rounded upward using the existing directed routines. Root resolution is
`h_d=1/[10^30(d+2)^6]`; the initial probability tail is first rounded upward
with resolution `h_d^2`, avoiding exponentially long serialized fractions.
After every stage, the outgoing moment amplitude and the Duhamel increment
are rounded upward to the least multiple of `h_d`. The projected error adds
that reported increment exactly. Each rounding adds less than `h_d`;
neutral and zero-time stages are already on the grid and unchanged.
For a fixed finite program with growth factors `G_i`, the amplitude excess
obeys `Delta A_(i+1) <= G_i Delta A_i + h_d`, so it is `O(h_d)`. Its effect
on later residuals, plus finitely many increment roundings, is also
`O(h_d)`. This preserves the directed bound and the exact
preparation-plus-increments ledger while preventing denominator growth.
The Hamiltonian residual is `O(1/d)`, pointer residuals are
`O(d^-3/2)`, and the final outside term is `O(d^-2)`. Thus the complete
occupation error is `O(1/d)` at fixed protocol, with no fixed rounding floor.
No assertion uniform in growing protocol length or growing total time is
made. Extremely large moments can still exceed the evaluator's integer
serialization limit; those requests raise an error and emit no certificate.

## Temporal and supplied-service budgets

The occupation theorem compares exact finite continuous evolution with the
specified infinite history. To compare the ideal source palindrome with
that history, replace each finite Hamiltonian interval in turn. The existing
all-input channel half-diamond estimate on interval `t_j` with `n` splits is

\[
 \delta_j\le\min\left(1,\frac{S_d^3 t_j^3}{3n^2}\right),\qquad
 S_d\le3+(1+2\kappa)(2d-1)+|g|\sqrt{4d-2}.
\]

The underlying operator-norm estimate has upper bound
`min(2,S_d^3 t_j^3/(3n^2))`; the channel distance is capped at one.
Each later complete
instrument is contractive, including record-dependent control. Hence

\[
 D(\rho_{d,n}^{\rm ideal\ product},\rho_\infty)
 \le\min(1,E_d+\sum_j\delta_j).
\]

`complete_error_budget` reads each `evolve` operation's exact rational
duration; it does not infer time from an aggregate field. At the executed
dimensions `d=2,4`, it calls the existing
`scalar_coupling.time_error_bound` for every such interval. Larger even
cutoffs use the same analytic generator estimate without claiming source
execution. For fixed finite protocol, `n(d)=d^2` makes temporal error
`O(1/d)` as well. The small executed cutoffs can still have a vacuous bound
of one; a controlled asymptotic theorem does not establish small-cutoff
accuracy.

Floating source contractions, angle approximations, proper-code load/export,
preparation, pointer, feedback, record storage and other physical services
require their own error bounds. A measured floating matrix residual is not
an outward numerical certificate or a physical service bound.
`service_error_budget` therefore requires explicit obligation keys and
nonnegative exact rational error bounds. `required_service_keys(operations)`
enumerates preparation, **every** operation and storage. A larger native
composition must extend these keys for its surrounding services. Unknown
keys are rejected; missing keys are listed, and the complete bound is `None`
until the supplied mapping covers every required key. Empty mappings do not
imply zero error.
Here the preparation obligation includes the native-current gauge kick as
well as the initial scalar/gauge preparation; alternatively the caller can
give that kick a separate required key. Its cost is counted by the runtime
even though it occurs before the operation list.

When all keys are supplied, the returned value is conditional on those
inputs being valid preparation trace-distance or complete-stage
half-diamond bounds, uniform over branch controls and reference systems:

\[
 D_{\rm implemented,\infty}
 \le\min(1,E_d+\sum_j\delta_j+\epsilon_{\rm prep}
                            +\sum_\ell\epsilon_\ell).
\]

The arithmetic does not establish those error premises. Each key covers a
whole complete stage, not one selected outcome. Repeated uses must be
accounted for; alternative mutually exclusive outcomes do not justify a
spurious branch-count multiplier. Service durations and service error bounds
also do not follow from the dimensionless Hamiltonian time.

## What the later clock means

After its last smooth read, the canonical process evolves its actual
retained branch states for another model interval and reads a gauge bit.
That whole instrument is covered by the new occupation and temporal bounds.
Comparisons of settings that change the final interval are separate complete
programs with the same mathematical preceding prefix, which they reexecute.
A checkpoint binds the full program identity, so it cannot be restored into
a program with a different later clock setting. These comparisons must not
be presented as successive observations of one sampled clock record.

On the finite working algebra a supplied faithful Gibbs reference
`sigma_beta=exp(-beta H_d)/Z` gives the exact modular identity
`sigma_s(A)=exp(i t H_d) A exp(-i t H_d)` with `t=-beta s`. Attaching this
identity to a Hamiltonian segment does not substitute that reference for the
actual observation state. Its central record weights cancel blockwise.
The identity does not turn measurement or feedback into modular
automorphisms, select the reference or temperature, establish monotonic
single-shot clock readings, or calibrate physical time. The occupation
comparison proved here controls the actual evolving instrument; it is not a
separate convergence theorem for Gibbs modular data or state-dependent
modular generators.

The two sites, action, scalar vacuum, quantum probability and instruments,
feedback dictionary, source services and clock interpretation remain
supplied construction data. The result establishes composable finite
execution and a controlled occupation comparison for that specified
process. It supplies neither spatial refinement nor new empirical physics.

## Verification

`test_continuation_bounds.py` compares the rational recurrence against
independent high-precision analytic integrals, checks the whole finite
pointer-carrier moment inequality, checks later post-read growth and the
single final outside tail, and rejects sharp reads, floats, forward record
references and unsupported couplings. It recovers the prior one-pointer
bound and reuses the exact source temporal certificates at `d=2,4`.
Negative service-coverage tests ensure that absent bounds remain absent.
These finite tests corroborate the calculation; the all-cutoff statement
rests on the induction and domain premises above.
