# Hilbert space, modular clocks and operational regions

This construction extracts observer representations and spatial readouts
from the **same executed matter/link/scalar history**. Its extra assumptions
are explicit: a finite quantum algebra with the Born state pairing, the
declared source action and preparation, and specified subsystem access.
The construction does not infer a quantum algebra from two classical bits.

Fields and spacetime are required as effective descriptions of observations:
the descriptions must give the right joint statistics, intervention responses
and compatible limits. They need not be primitive objects of the implementation.
Calling them effective does not remove these quantitative requirements.

The [local history-repair formulation](HISTORY_REPAIR.md) gives this same
experiment a completed, timeless description. Its compatibility equations
are fixed independently of the solver, and every fair local update schedule
finds the same assignment. The Hilbert/state/modular extractions below
therefore agree across those schedules. An observer's internal event ordering
and declared meaning-dependent feedback are features of the solved history;
the auxiliary repair counter is not its time readout.

## An observer Hilbert space from the actual observation state

Let `v_ab` be the four unnormalized terminal vectors already produced by
[experiment.py](experiment.py), with `p_ab = ||v_ab||²`. Retaining all records
gives the state

```text
rho_history = direct_sum_(a,b) |v_ab><v_ab|,
A_record+work = direct_sum_(a,b) M_(4d²).
```

For any positive state functional `omega`, the GNS construction gives

```text
< [A], [B] > = omega(A* B),
pi(X)[A] = [XA],             Omega = [I].
```

The quotient removes the zero-norm left ideal. In a matrix block with
`rho = V diag(lambda) V*` of rank `r`, the explicit map is
`[A] -> A V sqrt(diag(lambda))`, a `d` by `r` matrix with the Hilbert–Schmidt
inner product. This proves the state pairing and left representation
directly, and gives GNS dimension `d r`. A pure state's full matrix-algebra
GNS space is not its one-dimensional support corner. Its representation
remains faithful; its cyclic vector is not separating.

These are exact finite-dimensional formulas. The executable spectral
factorization uses a declared `1e-12` rank threshold and floating arithmetic;
its spectral values are numerical diagnostics, not outward enclosures.

All four history branches have positive weight and pure working states.
Consequently the record/work block algebra has GNS dimension `4(4d²)`.
Restricting access to just the four committed record projectors gives a
four-dimensional classical GNS space. These are different observer access
choices on the same history.

[operational.py](operational.py) also implements the observer that accesses
scalar site 0 and forgets both classical records:

```text
rho_0 = sum_(a,b) Tr_(site1,gauge/link)(|v_ab><v_ab|).
```

This is a marginal of the complete history, without postselection or an
added noise floor. Its sign-coordinate measurement reproduces the original
second-record probability. At both executed cutoffs its numerical spectrum
is strictly positive, so its full `M_d` GNS space has dimension `d²`.

| Scalar cutoff | Working Hilbert dimension | Record/work GNS dimension | Site-0 GNS dimension |
| --- | --- | --- | --- |
| 2 | 16 | 64 | 4 |
| 4 | 64 | 256 | 16 |

The regional representations are joined explicitly: the map
`J_R[A]=[embed_R(A)]` embeds each regional GNS space isometrically into the
same record/work GNS space, takes its cyclic vector to the complete history
vector, and intertwines the regional observable action. The implementation
checks all three identities for both scalar observers. Their expectations
therefore come from one representation and state, rather than separate fits.

The matrix quantum structure and access choice are premises; obtaining this
representation is a conditional reconstruction, not a derivation of the
Born rule or these premises from native integer repair.

## Nontrivial modular flow and a readable finite clock

For a faithful finite state, put `K = -log(rho_0)` and define

```text
sigma_s(A) = rho_0^(is) A rho_0^(-is)
           = exp(-is K) A exp(is K).
```

Here `s` is the dimensionless **modular parameter**. The implementation checks
the group law, multiplication, adjoints, state preservation and the finite
KMS identity, with this sign convention. Number observables move
nontrivially. The terminal sign-coordinate effect is fixed because the
unconditional state has already undergone that QND readout.

This is the finite algebraic mechanism behind the
[thermal-time proposal](https://arxiv.org/abs/gr-qc/9406019); identifying it
with an operational or physical clock requires a further instrument and
calibration. In particular, `Tr(rho_0 sigma_s(A))` is constant: the reference
state alone does not tick under its own flow.

There are two important access boundaries. The commutative record algebra
has identity modular flow. Conditioning on the sharp terminal scalar bit
also leaves a nonfaithful state on the full scalar algebra (rank at most
`d/2`). The implementation rejects a full-algebra matrix logarithm for
those states; it does not fill their kernels with invented noise. A support
corner is a possible different algebra and must be declared as such.

Discarding classical memory is not necessary for every modular observer.
The same implementation also keeps **both records** and restricts quantum
access to scalar site 1 in each record block. Each of these four actual
reduced states is numerically full rank at `d=2,4`. The algebra
`direct_sum_(a,b) M_d` then has a nontrivial blockwise modular flow and GNS
dimension `4d²`, while its central committed record projectors remain fixed.
This is a concrete finite record-aware modular observer on the same history.
The simpler clock pulse program below uses the unconditional site-0 state.

At `d=2` the actual reduced state has the particularly simple form

```text
rho_0 = (I - delta X)/2,
delta = P(b=0) - P(b=1),
omega = log(P(b=0)/P(b=1)),
K = constant I + (omega/2) X.
```

[modular_clock.py](modular_clock.py) continues every actual terminal branch
with a source `Z(pi/4)` intervention, source `X(s omega/2)` evolution and a
fresh number-bit pointer. The compiler convention is `exp(+i angle Pauli)`;
the second pulse implements the Schrödinger dual of `sigma_s`. The complete
eight-branch extension preserves the original two-bit record marginal and
has the readout law

```text
P(number=1 | s) = 1/2 + (delta/2) sin(omega s).
```

Both frequency and contrast come from the actual observation state.
The undisturbed reference state gives a constant pointer law, and full-cycle
separation gives the same law. Thus ensemble frequencies can read a modular
interval on a declared monotone phase window; one random bit cannot reveal
an exact duration or an unbounded cycle count. The pulse schedule, fresh
pointer and coherent services are supplied. This is a programmed finite
modular clock instrument, not proof that the original interacting Hamiltonian
already implements modular evolution or that an autonomous clock has emerged.

The existing OPH clock contract makes the same distinction:
[Definition 4.2a](../../paper/tex_fragments/PAPER.tex) requires a readable
transition, event correspondence and calibrated clock instrument beyond
BW-branch modular ordering. The
[finite clock certificate](../../paper/tex_fragments/OBSERVERS_TECHNICAL_SUPPLEMENT_PORT.tex)
assumes extraction/chart/clock certificates and proves invariance under
allowed execution changes. The
[time ledger](../../Lean/Time/TimeOrderLedger.lean) and
[proper-time calibration](../../Lean/Time/ProperTimeCalibration.lean)
likewise separate these interfaces. This packet does not instantiate their
geometric/BW, repair-confluence or physical calibration hypotheses.

For the present state, the local density derivative under the original
interacting `H` is nonzero, whereas its modular state is stationary. Hence
those evolutions cannot be identified on all local expectations by a
nonzero rescaling. A Gibbs assumption `rho proportional to exp(-beta H)`
would give `sigma_s = alpha_(-beta s)`; that assumption has not been shown
for this interacting observation state.

This obstruction does not invalidate the programmed clock above: its probe
first changes the actual state relative to the modular reference. The
[fixed-reference bridge](THERMAL_CLOCK.md) instead chooses
`sigma_beta = exp(-beta H)/Z` on the full working algebra, separately from
the actual observed density. Its modular flow agrees with the existing
Hamiltonian flow at `t=-beta s`. This is an explicit reference/access choice,
not a proof that the local post-measurement density was Gibbs or that repair
selects the new reference. The original negative test remains valid.

## Spatial meaning on the same carrier

[spatial.py](spatial.py) uses the supplied factorization
`H_site1 tensor H_site0 tensor H_gauge/link`. Each neutral scalar site has
a commuting regional observable algebra. The charged species and electric
link form **one shared Gauss-constrained sector**; the two encoded charge
bits do not constitute independent spatial regions.

State restriction is physical partial trace, retaining every branch weight.
The observable conditional expectation is normalized partial trace followed
by identity embedding. Its Kraus form proves complete positivity; it is
unital, idempotent and a bimodule map. It is generally not multiplicative.

The interaction graph is read from the actual Hamiltonian's tensor-centered
components, with double-commutator witnesses and deletion controls:

```text
scalar site 1 ---- scalar site 0 ---- shared charged-matter/electric-link sector
```

The scalar bond and scalar/charge interaction supply the two edges.
Removing the corresponding interaction removes its edge. The graph is
invariant under independent regional basis changes. The region assignment
is supplied; the measured interaction support is calculated. No distance,
metric, three-dimensional embedding or relativistic causal geometry follows
from these three nodes.

The native source has a separate, inherited
[rank-three slow response carrier](SLOW_GEOMETRY.md). Its exact Gram
quotient is attached to the actual recorded native currents in this same
joint experiment. The original controller spans two displacement directions;
an additional binary continuation executes three independent directions.
The full all-seam slow-band calibration and each controlled protocol have
distinct schedules. A rank-three response quotient is not by itself a
three-dimensional quantum field or a physical position/momentum pair.

The distinction matters for theorem adapters. The legacy
[FiniteCausalObserverNet](../../Lean/QFT/FiniteCausalObserverNet.lean)
requires multiplicative restriction retractions. Ordinary quantum partial
trace cannot fill those fields: a unital homomorphism `M_(nm) -> M_n` for
`m>1` would be injective by simplicity, contradicting dimension. This
regional construction supplies an explicit CP contract instead; the legacy
contract remains applicable only where its stronger hypotheses actually hold.

## What the physical tower must connect

The tower is the compatibility data that make different resolutions and
readouts describe **one common experiment**. Several distinct contracts are
currently named towers; they cannot be interchanged.

| Contract | What it connects | Status for this construction |
| --- | --- | --- |
| Exact hardware inventory | The same finite observables, states, records and generators under unital embeddings | Extra independent hardware may preserve a fixed experiment; adding interactions generally breaks exact generator naturality |
| Occupation and temporal approximation | Complete quantum instruments and their observations as cutoffs and splitting are refined | Temporal bound already certified; [interacting occupation analysis](OCCUPATION_REFINEMENT.md) controls pre-readout evolution and proves strong convergence of the final sharp readout, without a finite-cutoff readout rate |
| Bare repair tower | Normal forms, protected boundaries and coarse-graining of terminating repairs | A strictly decreasing repair mismatch is not the continuing unitary evolution |
| Physical readout tower | Common-source event, modular, stress, entropy, geometry and scale readouts, commuting maps and limiting physical identities | Regional algebra/modular/clock attachments now exist at finite cutoffs; the full physical bundle remains uninhabited by this source |

The exact matrix tower is [ConsensusTower](../../Lean/Tower/ConsensusTower.lean).
For its embedding `X -> X tensor I`, generator naturality requires
`H_fine - H_coarse tensor I` to commute with all old observables. A genuine
old/new interaction usually fails this requirement. Likewise the actual
occupation inclusion `V: d=2 -> d=4` obeys `V* H4 V = H2`, but
`||(H4 V - V H2)|blank>|| = 9/64`; two outside-cutoff amplitudes are
`9 sqrt(2)/128`. Compression equality is not exact dynamical refinement.

The [physical readout contracts](../../Lean/ObserverPatchHolography/EinsteinBranch/Consensus.lean)
and [Einstein register](../../Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean)
require actual compatible readouts plus the relevant physical premises.
Populating a record with unrelated modular, stress and geometry objects
does not prove its commutation squares or its analytic hypotheses.

The remaining major bridge is a growing **spatial** source family carrying
these same charges, interactions and instruments, with controlled local
observation limits. Two oscillators at increasing occupation cutoff remain
two oscillators. Responsive geometry additionally needs stress and geometric
readouts from that same family, the balance/Ward identities and controlled
tails, and a calibrated scale. Agreement with nature then requires
independent measured quantities. None of these requirements demands that
fields or spacetime be primitive ontology.
