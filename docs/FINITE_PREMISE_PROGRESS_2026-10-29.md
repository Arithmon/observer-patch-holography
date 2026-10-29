# Finite theorem attachments after the particle push

29 October 2026. CP-13 / [#1137](https://github.com/FloatingPragma/observer-patch-holography/issues/1137),
with the selected construction issues [#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139),
[#1140](https://github.com/FloatingPragma/observer-patch-holography/issues/1140) and
[#1142](https://github.com/FloatingPragma/observer-patch-holography/issues/1142).

The work supplies twelve realized finite attachments across the declared
Weyl-contact action, native coherence and continuing observer interfaces.
It discharges no canonical physical premise. Gate acceptance is recorded
separately in the [completion registry](COMPLETION_PROGRAM.md). The [source-bound map](../model/completion/finite_attachment_map.json)
assigns an owner and a relevance decision to every one of the 84 canonical
premise rows. Its seventeen inspected theorem/model attachments specify the
state, algebra, operation, time parameter and refinement domain separately.

This review starts from the existing
[theorem guide](../model/completion/THEOREM_MAP.md),
[premise map](../model/PREMISES.md) and
[canonical register](../tracking/premise_register.json). It preserves the
existing navigation coverage of 346 registered claims and 662 tracked Lean
modules. Inventory coverage is not a fresh proof audit of those modules,
and neither inventory size nor the number of mapped rows counts unproved
mathematical statements. Most canonical rows are physical attachments,
declarations or alternative-branch hypotheses.

## What is supplied on the actual finite model

The action is the existing supplied three-dimensional lattice CAR Hamiltonian
`H=dGamma[-(1/8) sigma.sin(k)]+sum_x n_x0 n_x1`. Open boxes support the
regional and stress calculations; periodic boxes support momentum fibers.
They use the same local terms and different finite boundary conditions.
No exact equivalence of their finite observation laws is assumed.

| Finite obligation | Actual attachment and independent control | Boundary |
| --- | --- | --- |
| Regional state restriction | Signed full-Fock CAR split and partial trace, retaining local N=0,1,2 sectors, all detector branches and arbitrary references. Nested restrictions compose. | Selecting the field carrier and physical regions remains a premise. |
| Regional algebra and CP expectations | Number-preserving local CAR images, explicit Kraus conditional expectations, locality and tower law; every two-mode region generates the full N=2 algebra jointly. | A mathematical attachment to the existing interface, not a compiled Lean inhabitant or a physical causal net. |
| Finite-volume observation map | Exterior-square coordinate embeddings preserve the compressed Hamiltonian; complete history errors contract under the same regional trace. | The full dynamics do not intertwine exactly. Growing boxes at fixed spacing are not a spatial continuum. |
| Native clock for an interacting energy pair | The same N=2 eta pair and dark reference give `H_S=diag(0,U)`; a supplied native-corner isometry identifies `rho_S=diag(25,8)/33` and `beta U=log(25/8)`. Complete source records retain every stop and leakage branch. | Faithfulness only on `M2`; no full-field source encoding, Gibbs preparation or physical time unit. |
| Interacting local energy conservation | Exact Gaussian-integer identities `sum h_x=H`, `J_xy=i[h_x,h_y]` and `i[H,h_x]+sum_y J_xy=0`, with complete boundary flux and interaction-deletion controls. | A finite Hamiltonian continuity identity, not a symmetric Lorentz stress tensor or Ward limit. |
| Quantum first law on that action | On the full finite N=2 Gibbs algebra, `DS_tau(delta)=beta Tr(delta H)` and `Delta S=beta Delta E-D(rho||tau)`; a noncommuting pair variation is tested. | The reference is declared. Geometric bulk/central/edge and fixed-volume area identifications are absent. |
| Finite thermal positivity and KMS | A spectral Gram factorization and direct trace check on the complete physical N=2,K=0 Gibbs fiber, with all pockets and a non-Hermitian observable. | No local continuum OS reconstruction, physical vacuum or native full-field thermal preparation. |
| One pair's energy, curvature and number current | Complete two-fermion fibers, exact CAR charge-two/even-parity identities, number continuity and rational bounds on kinetic curvature. | The supplied pair is not an identified Standard Model species or native topological defect; its same-c relativistic mass relation fails. |

These are component obligations with explicit hypotheses. Calling them
realized does not turn a supplied quantum carrier, action, source service or
isometry into a consequence of the bare repair law.

## Which existing theorems actually apply

The finite CAR regional route uses
[`CPRegionalNet`](../Lean/QFT/CPRestrictionNet.lean).
Take regions to be mode subsets, inclusion and intersection as the poset
operations, and declared disjointness to mean disjoint **nonempty** subsets;
the last qualification supplies the interface's irreflexivity field. Use
the images of number-preserving local observables in the global N=2 algebra.
The expectation onto each image is the trace-orthogonal conditional
expectation. Explicit block matrix-unit Kraus operators supply complete
positivity, including untouched-reference extensions. The Lean interface
itself demands ordinary positivity and trace preservation; its name alone
does not certify amplified positivity. The source package proves the
stronger property at the witness level.

The old
[`FiniteCausalObserverNet`](../Lean/QFT/FiniteCausalObserverNet.lean)
also demands contravariant star-homomorphic retractions. The existing
[`no_matrix_factor_over_scalar_restriction`](../Lean/QFT/NoncommutativeWitness.lean)
excludes that interface above a scalar bottom and a full matrix factor.
It remains a valid conditional theorem and an alternative historical route.
The new full-matrix CAR construction uses the attainable CP restriction
instead. Its explicit witness `E(X)=0`, `E(X^2)!=0` prevents accidental
promotion to multiplicativity. No theorem or frozen receipt was deleted.

The source-clock attachment checks the **actual operator restriction** of
the pair action, rather than matching two spectra by name. Both levels
have particle number two and the same total momentum. The off-site dark
reference is not a vacuum. The finite Gibbs/logarithm identity underlying
[`GibbsReferenceEnergyIdentification`](../Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean)
is applied on this stated `M2` algebra, with Heisenberg model coordinate
`t=-beta s` for native modular parameter `s`. A state faithful on this corner
is not faithful on the entire spatial field.

The nearest existing finite entropy result,
[`thermoFirstLawData_passes`](../Lean/Thermodynamics/EinsteinPremiseLink.lean),
acts on the classical simplex tangent space and assumes a bulk/central split.
The new matrix-algebra proof uses spectral calculus and arbitrary Hermitian
trace-zero variations, with a noncommuting test direction. It supplies a
concrete noncommutative field attachment; it is not a new Lean proof or a
direct invocation of the older scalar theorem on differently typed data.

## Failed identifications retained as controls

Three tempting transfers fail on this actual action:

1. `J* H_large J=H_small` does not imply `H_large J=J H_small`.
   Boundary hopping gives a nonzero intertwining defect. Finite observation
   comparison uses its proven error bound instead.
2. Tracing a full-sector Gibbs state does not make the regional modular
   generator equal to the bare regional Hamiltonian. The two-site state
   is faithful on its eleven-dimensional N<=2 support, not on the full
   sixteen-dimensional Fock space. Its modular generator has a nonzero
   commutator with the bare regional energy and fails even an affine fit
   allowing `I,H_A,N_A,N_A^2`. This is a finite numerical counterexample to
   that identification, not a universal modular no-go.
3. A bound-pair rest energy difference and its kinetic curvature do not
   automatically obey the same-c relativistic relation. The exact rational
   enclosure gives `U*kappa/c^2>17/16` for the declared action. Common time
   and length rescaling leaves this ratio unchanged. Different physical
   identifications or controlled limiting actions require new evidence.

Finite identities and these scoped failures replace unsupported inferences
on the active path. They do not change any historical theorem statement,
canonical disposition or original failure domain.

## What prevents the physical finite-theorem applications

The exact energy currents are a useful input to
[`WardLimitManifest`](../Lean/Dynamics/WardLimitManifest.lean).
That structure additionally requires an inhabited directed regulator family,
a genuine shrinking scale, compatible charge/current transport, one chart
realizing all components of a common stress tensor, separating tests, weak
stress convergence, a physically defined nontrivial residual and convergence
of its finite zero residual. Enlarging a fixed-spacing cube supplies none of
the missing shrinking-scale or Lorentz stress identifications. A formal
zero residual would not meet the intended physical map.

The actual
[`EinsteinRegisterData`](../Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean)
lacks, on this construction:

- geometric cap modular flow and its diamond-kernel stress identification;
- the common central/edge split and physical null translations;
- symmetric geometry/stress tensors, null-cone matching and the universal
  source shared by entropy and tomography;
- fixed-proper-volume area variation and generalized-entropy stationarity;
- Bianchi geometry and one shrinking-ball family controlling all remainders;
- the common physical tower, vacuum residue and independently justified
  gravitational normalization.

The [same-action taste packet](../model/field_tastes/README.md) further proves
`||[H,Q]||^2=U^2 w(1-w)` for each exchange-symmetric momentum pocket
with contact weight `w`. A proper pocket with `0<w<1` is not a reducing
subspace. Complete evolution retains all eight pockets and their outcomes;
this exact failed removal map does not exclude other regulators or controlled
limiting decoupling. Its finite Gibbs reflection Gram and KMS identities are
positive finite witnesses, not a complete OS field reconstruction.

The interacting field also retains its naive Weyl nodes and supplied action,
CAR dictionary and spatial geometry. No physical chiral vacuum, complete
gauge/Higgs/current attachment or controlled chiral continuum is obtained
by assigning the isolated two-particle band a name. The native moving-memory
and field-band constructions use different state/law contracts; their
protection mechanisms are not silently combined.

## Ownership, cleanup and issue closure

Every canonical premise has explicit current completion owners in the
machine-readable map. For this finite target, six rows have realized
restricted-domain instances, 25 remain required but conditional, 49 belong
to alternative branches, and four are outside this target's empirical or
baryon scope. These are relevance categories, not 84 independent proof gaps.
All canonical discharge statuses and dispositions are preserved. In
particular, a physical premise can remain unfulfilled while a restricted
instance of its interface is executable.

The twelve finite obligations are recorded separately. **Zero canonical
physical premises are marked discharged.** The reviewed map and observer
acceptance have a separate scope from those physical obligations.
No imported source, frozen receipt or historical evidence is deleted.
The active route explicitly drops the failed homomorphic, exact-intertwining,
bare-regional-Gibbs and finite-to-physical-Ward inferences.

CP-03 requires a native boundary invariant and protected defect sector.
The constructed exterior-matched coherence has a fixed core and zero energy
excess, so it does not complete that requirement. CP-04 requires the same
native excitation's maintenance action, physical clock and stress/refinement
map. CP-06 requires the interacting chiral field and its complete physical
matter/current/vacuum attachments. The CP-13 map has a different acceptance
contract: reviewed theorem interfaces and cleanup joined to accepted CP-00
and CP-01 predecessors. Their [versioned class and observer review](EXECUTABLE_PROGRAM_CLASS_V2_PROPOSAL.md)
provides the concrete predecessor specification; the completion registry
records the independent acceptance decision.

The V2 map explicitly separates the declared finite trace valuation and
counted coherent service from global source-selected Born and preparation
provenance claims. PR-03 and PR-65 retain their canonical statements and
discharge statuses. Their global derivation routes are alternative branches
for this finite observer; they remain obligations before any later claim
invokes those physical conclusions. FA-17 proves the actual finite instance.
The predecessor review is retained at commit `c35b09c099e5950f736befcf126c654542d78e15`.

```sh
python -m model.completion.finite_attachment_map --check
python -m pytest -q -W error model/completion/test_finite_attachment_map.py
python tools/build_premise_register.py --check
python tools/build_model_theorem_map.py --check
```

The validator rejects missing ownership, changed premise contracts, stale
source hashes or theorem anchors, omitted clock/refinement interfaces,
missing controls and unsupported issue or premise closure. It verifies
custody and domain declarations; the mathematical proofs and model tests
remain the cited evidence.

## Further native and complete-band attachments

FA-14 constructs a faithful native coherence with the complete exterior and
native tape equal to the same product reference. Its exact integrated contrast
lifetimes are about 6.73355 C-attempts, 10/3 E-attempts and4.45401 mixed attempts.
It has an immobile coherent core and zero modular-energy excess. The general
exact-code argument gives exterior distance at least 96/395 and a scalar
compressed inherited clock with quantified leakage. The [native package](../model/native_localization/README.md)
retains its complete source preparation and 128 raw continuation histories.

FA-15 drives the full 8,128-state interacting field using native record-derived
rate estimates, counted local source words and all four pointer outcomes.
Its ideal-word average halfdiamond bound is below 0.038565 for |s|<=1, uniform
in inputs and untouched references. It adds a supplied global field Gibbs
reference; it does not inherit that reference from native repairs.
FA-16 proves why every finite independent native product fails that exact
inheritance: the eta-calibrated integer modular gaps conflict with the
same field's pair energy 13/12<E/U<3/2. The approximate generator residual
exceeds 1/24. The [full-band service](../model/band_clock/README.md) and
[exact exclusion](COMPLETION_ACCEPTANCE_AUDIT_2026-10-29.md) keep these domains distinct.

FA-17 supplies the finite continuing observer's literal grammar membership,
complete quantum instruments, record-only controls and independent resource
accounting. It includes an actual binary-federation intercell transfer and
retains larger-request resource stops. Its declared quantum representation,
coherent services and physical apparatus boundary are explicit.
