# Concrete repair-law candidates

**Find the simplest shared repair rule that satisfies the known constraints.**

Candidate register: **12 October 2026**. Earlier law-comparison baseline:
[`39366529`](https://github.com/FloatingPragma/observer-patch-holography/commit/393665297e551f2023fb07a2853188a589ad7c33).
This page owns the concise comparison of active laws. The
[model search](CONSTRUCTION.md) supplies derivations and test summaries; the
[premise map](PREMISES.md) tracks their connections to existing
theorems. These are research candidates, with no physical selection claimed.

The [runnable model](README.md) includes two finite observer protocols on
N1: the [native record/controller](observer/README.md) and its
[P-dependent binary detector](observer/p_response/README.md). These are
explicit implementations under declared premises; construction completion
does not require proving that their design is the unique possible law.

Each candidate uses the
[established-mathematics contract](../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md): define
its agreement target independently, state its execution and fault assumptions,
and map any imported convergence or extension theorem to its actual objects.
N1/N1-E and N2/N3/N4 retain their finite channel proofs; S1 uses the declared
linear overlap structure where applicable. A sheaf heat-flow theorem does not
turn their quantum instruments or conservative evolution into that heat flow.

Clock obligations concern observer-local readings and compatible comparisons
on the same process. They impose no preferred universal clock; supplied shared
benchmark times retain their stated scope. See the
[clock convention](CONSTRUCTION.md#geometry-response-and-local-event-clocks).

## Candidates and priorities

| ID | Concrete law | Role | Strongest established property | Principal missing connection |
| --- | --- | --- | --- | --- |
| **N1** | Nearest-integer seam repair with current record, channel C | Microscopic baseline and executable observer source | One recorded process supplies conservative load repair, a finite slow geometry readout and dimensionless clock relations. Finite record-only feedback and a designed P-dependent binary response run on this instrument. | Physical realization, clock and refinement; the quantum instrument, charge assignment and detector apparatus are declared hypotheses/design choices. |
| **N1-E** | Same load repair with contrast record, channel E | Active completion comparison | Idempotent completion of C on one seam; a conditional expectation on the binary shell. | Source selection of the retained algebra and record instrument. |
| **N2** | N1 current drives local six-level response at both seam endpoints | One-way preservation benchmark | Preserves N1's load instrument and locality; inherits completion and finite binary-network convergence while protecting local block populations. | Full source grammar and selection against N3; all-input marginal preservation excludes feedback. |
| **N3** | N1 current drives the local flat endpoint response | Simpler active alternative | Preserves the same load instrument and locality; one global change of coordinates reduces it to N1 with unchanged internal variables. | A reason to select or reject its trivial closed-loop response. |
| **N4** | One coherent load–response exchange conjugates binary C/E | Reciprocity experiment with abstract central-port attachment | Private preparations change bare load probabilities. Minimal central addresses preserve its quantum channel in a complete scalar-seam source grammar. It retains normalized port geometry but exhibits memory beyond a rescaled native heat law at the same full quantum reference. | Physical source selection, preparation/readout, retained-record service and clock comparison; j is not a bare-load current. |
| **S1** | Minimal canonical inertial completion of weighted scalar agreement | Separate field/component benchmark | Explicit Hamiltonian, conditional free-field continuum and controlled preparation/observation construction. | Derivation of inertia, quantum state, controls and physical time from the microscopic process. |

N2 and N3 retain the native repair marginal. N4 tests the smallest reciprocal
extension examined after a one-term obstruction, using the same state factors.
It is an experiment, not a preferred law by default. A fixed coefficient ratio
is an assumption even without an extra fitted parameter. S1 supplies useful
field benchmarks; its dynamics are not identified with these microscopic laws.
Common refinement, physical time, interactions, gravitational premises and
observational discrimination must be joined on one process before a complete
physical model can be claimed.

## Runnable observer protocols on N1

These protocols retain N1's event law and add an explicit finite record tape,
controller and stopping rule. Their initial domain is the uniform 66-state
binary total-two source; the native implementation also checks its full
78-state extension. Two events means two successive repair invocations, not
two possible complete histories or a two-dimensional source.

The [native observer](observer/README.md) records `(c1,j1)`, issues a supplied
forecast using that record, chooses seam `(0,5)` when `j1>0` and `(0,1)`
otherwise, records the second native outcome and halts. All 264 binary-source
paths are enumerated. Its action choice is binary; its record alphabet is
richer. The [capacity calculation](observer/CAPACITY.md) gives 11 terminal
record words and 16 reachable prefixes, allowing a four-bit classical prefix
encoding. That encoding does not include the complete program or apparatus.

The [electromagnetic detector protocol](observer/p_response/README.md) accepts
P and N as closure-derived premises and computes the existing source-only
structured-running `α(P)`. Its chosen charge assignment is
`Q_v=g(P)N_v`, with `g(P)=sqrt(4π α(P))`. A prepared path qubit and calibrated
loop exposure χ produce a committed bit
`Pr(y|j1,P,χ)=[1+y sin(g(P)j1χ)]/2`. Only y reaches the next-seam controller.
All 528 weighted histories are retained. The response
`E[g(P)j2]=(5g(P)/264)sin(g(P)χ)` vanishes under disconnected and independently
shuffled detector controls. The software design runs; charge realization,
flux source, probe preparation and readout are specified hardware choices.

Both constructions preserve committed classical record values through
feedback. That condition permits disturbance of coherences between record
values. Neither finite protocol is an online learner or a complete physical
universe. The detector reads the N1 current; it is not the bare-occupation
instrument required by the N4 test below.

## N1 and N1-E: the common native baseline

Use the twelve ports and thirty seams of the source icosahedron. The tested
load space is the 78 nonnegative integer configurations with total load two.
Its 66 binary configurations form an invariant shell. For oriented seam
e=(a,b), a fair coin c assigns the ceiling to a when c=1:

```text
d = x_a − x_b,
j = [d − (2c−1)(|d| mod 2)]/2,
y = x + j(unit_b − unit_a).
```

Thus the two endpoints become the floor and ceiling of their mean. All other
loads are unchanged. The operation conserves total load and never increases
sum_p x_p². A difference of one permits a wait or swap: the binary shell has
continuing motion, rather than pathwise equal loads. Uniform independent seam
attempts are the baseline schedule; prescribed seams are diagnostic controls.
The total-two sector, schedule and initial disturbance are supplied. The
general integer rule leaves an exactly all-equal configuration fixed; it does
not explain the origin of a departure from that state.

N1 keeps the record r=(e,c,j). On the declared quantum load space, its Kraus
operator and forgotten-record channel are

```text
L_r = (1/sqrt(2)) Σ_{x:j(x,c)=j} |y(x,c)><x|,
C_e(ρ) = Σ_r L_r ρ L_r†.
```

Each branch is an injective translation, making C a completely positive,
trace-preserving channel. N1-E instead records (e,c,d), where d is the
pre-repair contrast. It has the same classical transitions and different
quantum coherences. C is not obtained by classical coarse-graining of E:
combining E's labels does not restore C's coherent Kraus sums. The exact
single-seam relation is

```text
E_e² = E_e,   E_e C_e = C_e E_e = E_e,
C_e^n = E_e + 2^(1−n)(C_e−E_e),   n ≥ 1.
```

The classical repair law alone does not select either quantum instrument.
On the binary shell, E is a trace-preserving conditional expectation. Full
uniform seam mixtures of C and E converge to the uniform binary state on the
load space. Protected external registers require explicit separate factors.
This ensemble statement does not establish confluence of individual histories.

On classical diagonal states this is the lazy symmetric simple exclusion
process with two particles. Its generator per seam attempt is
`(1/60) Σ_e(S_e−I)`. The exclusion-process corollary of
[Caputo, Liggett and Richthammer (2010), section 4.1.1](https://arxiv.org/abs/0906.1238)
gives the exact classical gap `(5−sqrt(5))/60 ≈ 0.046066` on this carrier,
equal to the one-particle gap. This is an application of existing theory.
It does not determine the full quantum-channel gap or the enlarged N2/N4
gaps: off-diagonal operators and additional response factors require their
own analysis.

Under uniform binary equilibrium and uniform seam attempts, the same records
recover the carrier Laplacian, its rank-three slow Gram and a dimensionless
rate relation. Those readouts do not identify physical space or proper time.
See the [native update](CONSTRUCTION.md#native-repair-and-its-observer-contract),
[completion](CONSTRUCTION.md#native-repair-and-its-observer-contract)
and [geometry/clock readouts](CONSTRUCTION.md#geometry-response-and-local-event-clocks).

**Selection test:** use a specified coherent preparation and retained-record
measurement to distinguish C from E. Mean loads and classical histories
cannot settle this choice. Quantum preparation, readout and finite record
capacity must be accounted for.

## N2: repair current drives endpoint response

Supply one private six-level factor at each port, together with the existing
skew-Hermitian response matrices D_p from the
[source construction](../paper/tex_fragments/SOURCE_RECORD_GLUING.tex).
The response space is (C⁶)^⊗12. With a supplied real coupling θ, replace each
N1 branch by

```text
U_(e,j) = exp[θ j(D_b−D_a)],
K_r = L_r ⊗ U_(e,j)^(a) ⊗ U_(e,j)^(b),
```

with identity on the other factors. Both endpoints receive the same unitary.
The existing matrices, their incidence-preserving port identification, the
twelve factors and the coupling are declared inputs.

| Property | Established scope |
| --- | --- |
| Load and records | For every correlated input, tracing out response preserves C and each classical record block exactly. Coherences between distinct retained records can change. |
| Locality and symmetry | Disjoint seams commute; simultaneous proper port transport and response transport give finite carrier covariance. |
| Completion | A seam-dependent controlled unitary conjugates the update to C_e tensor identity. The same conjugation transfers N1-E and the exact single-seam completion formula. |
| Network convergence | On the finite binary shell, fixed positive seam mixtures converge to the conditional expectation onto their common fixed algebra. If γ_C and γ_E are their gaps off that algebra, γ_E/2 ≤ γ_C ≤ γ_E. Its full algebra and absolute gap are uncomputed. |
| Persistent information | Each local projector diag(1,1,1,0,0,0) is protected. Arbitrary private coherences are not claimed to survive. |
| Response histories | Closed load loops can leave internal memory. For an input that factors into load and product private states, the record-conditioned response factors across ports and creates no private-to-private entanglement. |
| Backreaction | None: changing internal response cannot change unconditional native load probabilities. This follows from exact marginal preservation for all inputs. |

The convergence result concerns forgotten-record channels. It does not prove
convergence of a growing coherent history or provide unlimited record storage.
Conserved current differences also omit the constant central response
direction: attaching D_b−D_a does not derive the full gauge structure or a
physical gauge field. Loop memory alone does not establish non-Abelian
curvature; a commuting-matrix control produces loop memory too.

**Evidence:** 268 fixed checks and one identical isolated replay cover native
source binding, covariance, locality, completion and a full three-event
comparison with N3. Numerical matrix checks supplement the analytic
identities; there is no native hardware execution. The
[detailed join](CONSTRUCTION.md#geometry-response-and-local-event-clocks)
states preparations and controls.

**Next decisive work:** characterize the protected/common algebra and test
the smallest reciprocal extension compatible with local repair. A reciprocal
extension must relax exact preservation of C on specified internal inputs;
requiring it for every input rules out that feedback. Check the enlarged
source's complete event grammar before importing the gauge theorem.

## N3: the local flat alternative

Keep N1 and the same response factors. Define J=Σ_p D_p and Dbar_p=D_p−J/12.
Use opposite endpoint increments:

```text
K_flat,r = L_r ⊗ exp(−θ j Dbar_a)^(a) ⊗ exp(+θ j Dbar_b)^(b).
```

For W(x)=tensor_p exp(θ x_p Dbar_p), the same load-controlled unitary W
conjugates every seam to C_e tensor identity. Consequently every closed load
loop has identity response, and the full encoded internal algebra survives.
Load preservation, classical records, locality and proper covariance hold as
for N2. There is no response-to-load feedback.

**Simplest tested discriminator:** one local six-level state, one seam attempt
and one local two-outcome effect distinguish N2 and N3 with probability
difference **0.049095 at θ=1/8**, including both fair coin outcomes. No entangled
reference or retained history is required. All 67 fixed checks and an
identical isolated replay passed. Preparation, effect and coupling are
supplied; see the [local readout](CONSTRUCTION.md#geometry-response-and-local-event-clocks).

The earlier declared reference-assisted input, retaining all
eight coin outcomes of the three-event history gives N2 versus N3 trace
distance about **0.137390 at θ=1/8**, about **0.069358 at θ=1/16**, and zero at
θ=0. No outcome is postselected. These are finite model comparisons using
supplied quantum references and seams, not measured physical predictions.
See the [alternative and fixed tests](CONSTRUCTION.md#geometry-response-and-local-event-clocks).

**Next decisive work:** realize that local preparation/readout in the admitted
source grammar and connect it to a physical comparison. Preservation and
convergence alone cannot select N2 over this simpler alternative.

## N4: a reciprocal coherent exchange

On the binary shell let S_e swap the endpoint occupations. Use
`X=S_e(N_a−N_b)²` and `Y=iS_e(N_a−N_b)`, with the same endpoint private
factors as N2. Define

```text
g=D_b−D_a, h=D_a+D_b−J/6, J=Σ_p D_p,
G=g^(a)+g^(b), H=h^(a)+h^(b),
W_e=exp[θ(Y⊗G+X⊗H)],     K_r=W_e(L_r⊗I)W_e†.
```

The exponent also equals `θ(s_+⊗B−s_−⊗B†)`, where
`s_+=|10><01|` and `B=H−iG`. It is one coherent exchange, with no extra
state variable beyond N2. θ, the equal relative weighting, exponential law
and quantum interpretation are supplied choices.

Both response directions omit the constant central generator; this exchange
does not restore the missing central response.

**Why test it:** any single-product conjugation `exp(θ A⊗G)` with local
number-preserving A leaves C's load marginal unchanged. Its channel symmetry
blocks feedback for that whole family. N4 evades this scoped obstruction;
neither summand alone has the tested feedback. This does not prove universal
minimality among other possible laws.

**Properties:** local, total-occupation conserving, proper-carrier covariant,
with protected 3+3 block populations. Same-seam C/E completion and finite
binary-network convergence transfer by conjugacy. Its full common algebra
and absolute gap are uncomputed. The candidate is restricted to the binary
shell; full integer-pile relaxation is not established.

**Evidence:** with fixed bare load I_2/2, supplied private preparations change
endpoint occupation probability from **0.482020673 to 0.517979327 at θ=1/8**.
The one-term controls give 1/2. All 187 fixed checks and an identical isolated
replay passed. The
[derivation and tests](CONSTRUCTION.md#geometry-response-and-local-event-clocks)
give independent source and small-coupling controls without claiming an
interval bound at finite θ.

**Main cost:** the old record label j refers to the conjugated occupation,
not to a measured change of bare load. N1's current-derived geometry and
clock identities cannot simply be imported. The
[shared bare-occupation protocol](CONSTRUCTION.md#geometry-response-and-local-event-clocks)
reconstructs `K=αL_ico/396` at the full tracial reference, with `α≥1/9`
for every coupling. It retains the finite rank-three port geometry and,
at θ=1/8 and 1/16, certifies nonzero two-attempt memory for every
nonconstant carrier mode. A change of clock units cannot remove it. No
additional state factor or interaction is needed for this test. On one seam
the channel is unitarily equivalent to C; conjugating the observers too
would change the fixed comparison. The
[central-address construction](../code/native_response/README.md#a-central-port-carrier-that-retains-the-quantum-repair)
preserves the native quantum channel in a minimal twelve-block algebra and
supplies a full abstract scalar-seam grammar. Its native/N4 counterpair shares
the full quantum reference and complete collective response, so these do not
select the finite kinetics. Sampling an occupied port has a proved quantum
information cost; central control addresses and occupation records serve
different operations. Physical preparation/readout, retained-record services,
gauge dynamics and an interacting continuum require their own realization.

The [local-tick extension](local_tick/README.md) derives the separate null
when the observer counts incident attempts and all remote events remain.
It preserves a strictly positive certified N4 distinction at both couplings.
This closes the stopping-rule calculation; occupation preparation/readout
and physical calibration remain separate from the implemented N1 detector.

## S1: the canonical scalar benchmark

Weighted seam agreement supplies a scalar stiffness. Given positive repair
capacities M, grounded stiffness K_D, mass m and conjugate momentum Π, the
tested conservative completion is

```text
H = Πᵀ M^−1 Π/2 + fᵀ(K_D+m²M)f/2,
f_dot = M^−1 Π,       Π_dot = −(K_D+m²M)f.
```

Using the same capacity as inertia is a hypothesis. Together with the stated
linearity, velocity, energy and time-reversal premises it fixes this
completion. It oscillates rather than asymptotically settling. The
[construction](CONSTRUCTION.md#retained-scalar-and-entropy-constructions)
reuses the existing free scalar continuum and quantum-instrument results
under supplied canonical quantum representation, state, controls and clock.

**Evidence:** an eight-site construction composes preparation and retained
observations on one finite code with mathematical trace-distance error
**<0.001237**, before physical service errors. A direct Gaussian simulator
evaluates the same quadratic benchmark; a finite compiler simplifies its
implementation. Neither is an additional physical law or an executed complete
source history. See the [shared error and resource account](CONSTRUCTION.md#retained-scalar-and-entropy-constructions).

**Next decisive work:** connect the inertia and Hamiltonian time to the same
microscopic events as N1/N2, and execute a complete preparation/readout history
within one service budget. Equal equilibrium configurations do not select
the kinetic law; the retained [kinetic alternatives](CONSTRUCTION.md#retained-scalar-and-entropy-constructions)
provide dynamical discriminators.

## Construction principle and rejected alternatives

The [canonical repair RFC](../docs/CANONICAL_REPAIR_LAW_RFC.md) is a proposed basis
strengthening, not an adopted microscopic law. Given a complete grammar,
faithful reference, modular-invariant agreement algebras and schedule, it
constructs a mixture of conditional expectations. Those premises constrain
the candidates; the RFC does not supply the missing inputs by itself. In
particular, an odd-total integer seam reaches a nearest-balanced shell,
not an exact pathwise equality retraction.

The N2 variant with **one shared six-level processor** preserves the load
channel but fails disjoint-seam commutation. It is rejected under the declared
locality requirement. The local endpoint construction is retained instead.
Other superseded routes and their reasons are collected in
[retired approaches](CONSTRUCTION.md#remaining-construction-work).

## Keeping this register useful

Follow the [general search findings](CONSTRUCTION.md#remaining-construction-work):
use discriminating requirements, test deletion controls, count independent
assumptions and try the smallest diagnostic first. A shorter rewritten
formula does not by itself reduce the law's physical inputs.

Update the relevant card when the law, state space, parameters, schedule,
record policy or physical interpretation changes. Preserve the candidate ID
and record the changed version in Git history; create a separate ID when a
scientifically distinct alternative must be compared. Promote results only
with their scope, added assumptions, discriminator and evidence reference.
Keep exploratory logs and private contributions outside this register.

For each claimed recovery or agreement result, record the legitimate set,
preserved data, perturbation domain, schedule class, nearest primary theorem,
checked premise mapping and remaining physical discriminator. Distinguish
support extension, probability extension, harmonic residual and uniqueness.
The [bounded work products](../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md#7-work-products-and-scientific-contribution)
fit the existing owners below; they do not retire successful component evidence.

[#1076](https://github.com/FloatingPragma/observer-patch-holography/issues/1076)
owns shared-law selection;
[#1077](https://github.com/FloatingPragma/observer-patch-holography/issues/1077)
owns complete source execution. The [premise map](PREMISES.md) and
claim registry retain their existing roles. No candidate card closes a
physical premise by renaming it. Completion requires the simplest shared rule
to satisfy the required constraints on one common model, with discriminating
evidence for its physical identification.
