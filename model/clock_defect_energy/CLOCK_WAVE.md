# Existing clock bridges and a reference-substitution check

The model has positive finite clock/field constructions. In particular,
native return histories reconstruct a reference for the **same coherent wave
generator**, and complete native clock records give a bounded-window inverse
for an observer's enacted parameter. These results must be used before
testing a different reference as a replacement for the wave reference.

An earlier version of this packet called equality of clock and field
Hamiltonians a common-clock test. That framing was too strong and is
withdrawn. Different clocks can time the same field through an operational
calibration; their Hamiltonians and spectra need not coincide. Equality of
the entire inner automorphism groups is a separate, stronger hypothesis.

The secondary substitution check compares the four references used by pending
[PR #1149](https://github.com/FloatingPragma/observer-patch-holography/pull/1149)
at commit `3d282b3d59763fa1a5171ad9a4fbc09fc0abf35a` with the existing native
wave kernel. The PR is an unmerged input at this campaign boundary, not a
claim that its full implementation has been independently audited here.

The exact nonidentities below reject that particular substitution. They do
not obstruct a common operational parameter or supersede the positive
clock/field results. Assigning SI c supplies no missing coupling or readout
map and changes none of these distinctions.

## Reuse the positive wave reference

The [native-return construction](../observer_time_fields/README.md#a-native-return-reference-for-the-coherent-wave-generator)
already uses the same 66-state kernel `K` as this packet. Its positive length
law `f(z)=exp(arcsin(z)-pi/2)` reconstructs

```text
sigma_wave = f(K)/Tr f(K) = exp(-Theta)/Tr exp(-Theta),
Theta = arccos(K).
```

Its modular unitary has positive wave phase magnitudes, and its cosine
amplitudes equal the existing arc-walk compression. This is a constructive
source-history reference, not merely an arbitrary fitted Gibbs matrix.
The admitted length sampler and full endpoint comparison are extra services.
The ideal law has infinite expected length. At finite cap h, the original
producer reconstructs `sigma_h=P_h(K)/Tr P_h(K)` and proves

```text
Theta_h=-log(P_h(K)/P_h(1)),
0 <= Theta-Theta_h <= -log(1-delta_h) I.
```

This packet reuses that producer at caps two and four, verifies exact
commutation and preservation of the constant-plus-density space, and keeps
the complete return, different-endpoint and cap probabilities. It also
retains the existing cap-512 reference and generator bounds. The conditional
working source is diagonal and differs from the reconstructed reference;
their exact squared difference is reported. Neither the old construction nor
this replay prepares that coherent reference or implements every vertex
modular intervention. This positive bridge is **existing evidence**, not a
new theorem credited to this packet.

## Readout calibration does not need equal generators

The [oriented phase inverse](../native_modular_reference/README.md#infer-an-oriented-clock-parameter-from-complete-native-records)
uses complete cosine and sine score records from two overlapping caps. Its
learned enacted gap, source-error budget, declared window and conditional
concentration theorem produce a parameter interval `[s_-,s_+]`. The clock
gap is a difference of logarithms of reference weights, not an absolute
state-vector energy level.

Here the inverse is replayed from the retained counts of
`phase_clock_results.json`, without consuming its validation-only programmed
parameter. Shared training uncertainty is counted once. Both original laws'
joint confidence failure bounds are `1/128`, under the original source and
sampling hypotheses; the seeded demonstration does not certify randomness.
No fresh source-array execution or large sampling run is performed here.

Suppose an operational coupling or calibration separately establishes
`t=kappa s` for a field generator with `||H_field||<=Omega` on its declared
band. Let `s_hat=(s_-+s_+)/2` and `delta=(s_+-s_-)/2`. Duhamel's formula and
spectral calculus give, on the clock's confidence event,

```text
||exp(-i kappa s H_field)-exp(-i kappa s_hat H_field)||
    <= min(2, |kappa| Omega delta),
||cos(kappa s H_field)-cos(kappa s_hat H_field)||
    <= min(2, |kappa| Omega delta).
```

The same bound, capped at one, bounds the associated quantum-channel
half-diamond distance and any detector-probability error. Coupling, field
implementation and preparation errors require their additional allowances.
No commutation or spectral equality with the clock Hamiltonian enters.

The finite control uses the actual elementary constant-plus-density wave
band, whose frequency norm is certified below `1/2`, and explicitly declares
`kappa=1`. That declaration is a conditional comparison, not an inferred
physical relation or a newly executed joint instrument. Extending the inverse
through the pending transported-reference history, and executing the chosen
clock/field coupling with common resources, remains a concrete join to test.

## Relevant existing theorem and construction map

| Result | Actual predicate and achieved conclusion | Boundary for this packet |
| --- | --- | --- |
| [Common-source thermal clock](../common_source/THERMAL_CLOCK.md) | A faithful Gibbs reference of the supplied interacting H on full `M_(4d^2)` gives `sigma_s=alpha_(-beta s)`; continued detector records verify changing expectations in the actual state. | The actual state need not be the reference. Intervening instruments remain separate maps; no source selection of H, beta or physical seconds. |
| [Incident-clock field dictionary](../observer_time_fields/README.md#a-complete-field-law-reconstruction) | Complete local field histories and permitted feedback are intertwined with a positive transfer and `H=-log T`, including waiting and failures. | This is an observable/history correspondence, not equality of passive repair with a unitary flow. The [degree-one/two dual](../observer_time_fields/dual.py) is an exact exclusion observable map, not independent particles. |
| [Quantum observer correspondence](../quantum_observer_correspondence/README.md) and [CP half-wait extension](../observer_correspondence/QUANTUM_CLOCK.md) | Native C/E clock marks, coherent detector responses, adaptive controls, timeouts and untouched references are preserved on the declared matrix cones; the half-wait construction supplies an ordinary CP effective instrument for its stated output interface. | The H66 all-operator correspondence and the total-one executed half-wait histories have different domains. Neither needs a clock Hamiltonian equal to a field Hamiltonian. |
| [Joint local clocks](../joint_clock_process/README.md) | A controlled four-port one-particle experiment preserves full conditional quantum/reference outputs and derives an instrument-valued union-to-local clock change, including feedback and remote evolution. | The exact random clock conversion is not a deterministic scalar-rate assumption or a full-federation uniform schedule. The [clocked common observer](../clocked_common_observer/README.md) separately executes a stopped native clock coupled to continuing work evolution. |
| [Visible event geometry and field checks](../laboratory_clock/README.md) | Ordered incident-count comparisons recover source adjacency and predict held-out joint local occupation readouts from the same native source law. Conditional timestamp fits have explicit unused-observation checks. | Event delivery/alignment, scheduler and preparation are stated. This supplies useful dimensionless clock/field checks without laboratory-unit identification or matched Hamiltonians. |
| [Closed-history Gibbs and wave references](../observer_time_fields/README.md#a-gibbs-reference-from-native-closed-histories) | Incident closed cylinders reconstruct `T^L/Tr T^L`; the separate variable-length global-attempt law reconstructs `exp(-arccos K)/Tr`. | Full endpoint closure and finite conditional error budgets are explicit. Reconstructed references are not automatically prepared source densities. |
| [Encoded faithful field clock](../observer_time_fields/REFINEMENT.md#faithful-state-clock-compatibility) | On the specified mode algebra, matching axis and `log((1+r)/(1-r))=b Omega` gives full Pauli modular action, ordered two-point/KMS law and signed response; controlled mode/copy limits yield the thermal field. | Selected reference, encoding, preparation and common b remain services. This positive join already uses the field's Bohr frequency. |
| [Product field, state and clock](../observer_time_fields/REFINEMENT.md#one-product-field-state-and-active-clock) | The same `B=-d_r^2+(A direct_sum 0)`, thermal field, retarded response and active incident-clock conversion `t=(2d/sqrt(630))tau` are joined conditionally. | Radial domain, coherent controls and growing-mode access are declared; laboratory and spatial-patch modular identifications are separate. |
| [Native full/regional/common reference clocks](../native_modular_reference/README.md#regional-algebras-under-one-native-modular-action) | Faithful references, explicit regional algebras, modular-invariant completions, common core-containing caps and full source words supply valid local inner actions and continuing reads. | Different prepared references have different flows. Ambient restriction requires the proved invariant algebra; no equality with the wave generator is assumed. |
| [Operational cap and learned clock](../native_modular_reference/README.md#bare-cap-clock-operations-on-arbitrary-spectators) | A monotone mean response has a 200-Lipschitz inverse on its declared interval; complete quadrature scores give oriented finite-confidence parameter inference. | A random individual record is not a strictly monotone clock. Learned enacted gaps and unknown ideal-reference gaps remain distinct. |
| [Event-driven transfer probe](../transfer_modular_probe/README.md#a-native-event-drives-the-modular-step) | Each zero-inclusive incident event triggers one fixed coherent clock word; full histories and finite reading precision are controlled. | Its three-class reference is distinct from the field wave. Protected registers, paused services and a declared coupling are required. |
| [Regional scalar vacuum modular response](../native_equilibration/ball_modular.py) and [fermion ball](../common_source/fermion_ball.py) | Restricted covariances determine regional modular charges of specified Gaussian/CAR vacua; same-action response is compared within fixed regulator/region conventions. | These are not clock-readout or global-to-patch flow identifications for the H66 reference. Their supplied vacuum/geometry domains remain explicit. |
| [Finite inner dynamics](../../Lean/Dynamics/PrivateInner.lean) and [finite Gibbs/KMS](../../Lean/QFT/StructuralInheritance.lean) | Automorphisms of a finite full matrix block are inner; exponential flows have commutator derivatives; a specified finite Gibbs functional satisfies the algebraic KMS identity. | They do not choose a generator, reference, clock readout or physical units. The scalar-commutant necessity used below is an application of this existing framework. |
| [ClockReadout](../../Lean/Time/ClockReadout.lean), [shared-event comparison](../../Lean/Time/ClockComparison.lean) and [proper-time calibration](../../Lean/Time/ProperTimeCalibration.lean) | Supplied monotone readings admit regraduation; two ordered shared events fix an affine map and a third tests it. A supplied affine unit-timelike history realizes proper duration. | These are operational calibration predicates, not equality of clock and field Hamiltonians. Statistical mean inversion does not by itself inhabit the individual-record interface. |
| [Curl/Stone bridge](../../Lean/Screen/CurlStoneClockBridge.lean) | The coefficient flow and its projective matrix image are intertwined; a selected clock-mode relation follows under explicit duration/rate hypotheses. | The [Maxwell/clock product](../../Lean/Geometry/CommonWorldMaxwellClockJoin.lean) is only same-index packaging unless the additional relation is supplied. [Source-count enclosures](../../Lean/Time/SourceCountClockEnclosure.lean) and [worldline rate comparisons](../../Lean/Geometry/SourceClockRateAlongWorldlines.lean) retain their separate volume/accrual hypotheses. |

## The stronger reference-substitution comparison

`spatial_dynamics.native_density()` enumerates all 3,960 endpoint/coin
branches. Its symmetric configuration count matrix `C` has column sums 60;
`K=C/60` is the compression of the admitted coherent arc walk. The
configuration wave generator is `H_w=arccos(K)`: its compressed integer-step
amplitudes are `cos(n H_w)`. This spectral generator is not the unread
repair process `K^n`, and its preparation, coherent record access, reflection
and readout remain services of the existing wave construction.

The pending complete-record clock has diagonal native reference weights

| Occupation in the transported two-port core | Per-configuration weight | Multiplicity |
| --- | --- | --- |
| 0 | `2/225` | 45 |
| 1 | `1/36` | 20 |
| 2 | `2/45` | 1 |

These are the prepared ideal target reference and its full-shell extensions,
not the unknown finitely prepared state `sigma_n`; the original preparation
and logarithm-stability bounds retain that distinction.

The producer enumerates the actual C and E two-attempt histories and derives
their current-tape port permutations. Both laws reach cores `(0,1)`, `(0,8)`,
`(1,5)` and `(5,8)`, giving precisely the four existing clock words. Their
full-shell references are permutations of the initial reference. A branch's
actual posterior reference is its supported compression, with its coin and
record probability retained by the original experiment. This diagnostic does
not substitute full-shell weights for those posterior probabilities.

Write `H_sigma=-log(sigma)`. Equality of the two continuous inner actions on
the same full matrix algebra, under any nonzero scalar change of parameter,
would imply `H_sigma=alpha H_w+beta I`. The scalar shift is harmless; it is
included rather than silently fitting absolute energy zeros.

## Exact invariants against that substitution

First, such an identity requires `[sigma,K]=0`. Actual native references
give the following squared Frobenius norms, with an explicit matrix-unit
witness for every row retained in the JSON receipt:

| Core | Nonzero commutator entries | `||[sigma,K]||_F^2` |
| --- | --- | --- |
| `(0,1)` | 160 | `157/10125000` |
| `(0,8)` | 200 | `157/8100000` |
| `(1,5)` | 160 | `157/10125000` |
| `(5,8)` | 200 | `157/8100000` |

The tracial control commutes, demonstrating the test's necessary-condition
boundary: commutation alone does not establish matching generators.

Second, four Krylov vectors `e_0,C e_0,C^2 e_0,C^3 e_0` have exact Gram
determinant **1,556,992**. Thus the symmetric `K`, and its injective spectral
function `arccos(K)`, have at least four distinct eigenvalues. Every
transported `H_sigma` has only three. Even an arbitrary full-shell unitary
change of encoding cannot make the generators affine copies. This is an
exact integer certificate, not a numerical eigenvalue or rank threshold.

Third, let `Z[p,x]=x_p-1/6` be the established native density dictionary.
The constant-plus-density subspace has orthogonal projection
`P=Z^T Z/10+11^T/66`, of rank 12. For every reached core,

```text
||(I-P) sigma Z^T||_F^2 = 2119/371250 > 0.
```

This subspace is not invariant under the transported reference, so it cannot
reduce its modular generator either. Restricting the comparison to these
field amplitudes does not silently supply a closed reference clock.

There is also a separate three-level **generator-level** check, even if an abstract
encoding selects one vector from each reference eigenspace and one from
each nonconstant density band. Those native Laplacian bands are
`5-sqrt(5), 6, 5+sqrt(5)`, and their wave frequencies are
`arccos(1-lambda/60)`. For sorted levels define

```text
R = (highest-lowest)/(middle-lowest).
```

This is not a comparison of observable Bohr frequencies: modular observables
rotate at differences of `-log(sigma)` levels, and a clock need not reproduce
all field frequencies. The modular generator-level value is
`log(5)/log(8/5)` for positive affine orientation and
`log(5)/log(25/8)` after reversal. Rational square-root, positive arcsine
series and logarithm series enclose all three ratios in disjoint intervals.
The retained receipt records their separation. This test includes arbitrary
energy offsets and either time orientation. It does not assert that the
selected density vectors already form a subspace of the reference clock.

## An additional two-mode algebra control

The native port vector `(1,-1,-1,1,0,0,0,0,0,0,0,0)` has Laplacian
eigenvalue 6. Its density lift and the constant configuration vector give
an actual orthogonal two-mode code with `K=diag(1,9/10)`. On this code,
`H_w=diag(0,theta)`, where `theta=arccos(9/10)`.

Supply instead `sigma_control=diag(2/3,1/3)`. Then exactly

```text
b = log(2)/theta,
-log(sigma_control) = b H_w + log(3/2) I.
```

With the pending clock's predual convention `U_sigma(s)=sigma^(-is)`,
`s=-t/b` matches `U_wave(t)=exp(-it H_w)` up to global phase. A balanced
two-mode probe has plus probability `19/20` after one wave step. This is
a mathematical Gibbs compatibility control on actual source modes; its
preparation and encoding are different supplied services, not an attained
repair-selected reference or a new compiled execution.

## What assigning c can test

The existing continuum product wave has unit optimal front speed in its
intrinsic energy distance. Assigning one such distance unit `ell` metres and
one wave-parameter unit `T` seconds imposes `ell/T=c`. Its older
`1/sqrt(2)` speed in a supplied flat comparison geometry is an upper bound,
not an equality available for calibration. See the
[intrinsic-front theorem](../relative_observer_geometry/README.md#intrinsic-distance-is-the-product-waves-first-response).

The existing active incident-clock conversions, including degrees five and
six, already give a positive operational relation for their stated field
protocol. They need not remove the commutator, spectral multiplicities or
gap-shape differences above. Those differences only prohibit calling the
pending reference's full modular group the same wave group by substitution.
The positive references and operational dictionaries above remain available.
Proper duration, signal identification and physical rod/clock calibration
retain their stated obligations.

The packet does not derive electromagnetic light, matter energies, an
observed Lorentz metric or a new empirical comparison. It leaves finite
compiled clock errors in their original receipt and adds no gate acceptance.

```sh
python -m model.clock_defect_energy.clock_wave --check model/clock_defect_energy/clock_wave_results.json
python -m pytest -q -W error model/clock_defect_energy/test_clock_wave.py
```
