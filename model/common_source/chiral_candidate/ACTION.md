# A concrete all-left chiral Higgs candidate

[action.py](action.py) constructs a local Weyl/Higgs mass action from the
existing exact public internal representation and interaction tensors. It
uses the computed native family anisotropy, rather than adding three
independent fitted family masses. The map from that anisotropy to Yukawa
couplings and the scalar background are declared architecture choices.
No measured mass or coupling enters the default candidate.

## Existing representation, with its actual basis

The code reads, checks retained artifact hashes and reuses:

- [the exact representation operators](../../observer_particles/results/representation_identification_operators.json)
  and [their receipt](../../observer_particles/results/representation_identification.json);
- [the exact one-Higgs tensors](../../observer_particles/results/unified_yukawa_operators.json)
  and [their receipt](../../observer_particles/results/unified_yukawa_action.json);
- [the native family feature](../particle_identification/family_attachment.py),
  which binds the native reference, source carrier and existing family band.

These producers are read-only inputs; this package never reruns them to
replace retained receipts. The corresponding public theorem interfaces are
[WeylYukawaConventions.lean](../../../Lean/Screen/WeylYukawaConventions.lean)
and [FermionSectorAssembly.lean](../../../Lean/Screen/FermionSectorAssembly.lean).
The full distinction between an internal match and physical identification
is in [SM_MATCH.md](../particle_identification/SM_MATCH.md).

The fifteen-dimensional internal space is
`M15=Lambda²(C⁵)+Lambda⁴(C⁵)`. Every component uses a left Weyl index:
`Q,u^c,d^c,L,e^c`. In particular `u^c,d^c,e^c` carry the charges of
left-handed charge-conjugated singlets; they are not right-handed fields
with those same signs. The source charge operators use `q=6Y` and
`Q_em=T3+Y`. The public retained exterior subsets fix the component order.

The producer independently forms all twelve current matrices by
differentiating the exterior action in that retained basis. It exposes
the primitive five-dimensional currents, fifteen-dimensional currents and
their three-family tensor lifts. The charge, Cartan and Casimir operators
share the same order. `ordering()` exports each component's subset, internal
type and charges, together with the three neutral-neutrino indices.

The internal/family index is `k=3*internal+family`; the spin index is
outermost, `spin*45+k`. Family axes are the existing native coordinate
axes, not electron/muon/tau assignments. The family basis metric before
normalization is a scalar multiple of the identity; choosing canonical
kinetic normalization leaves the native feature's matrix unchanged.

## Explicit trial coupling map

The actual protected-core preparation has exact reference weights
`(w0,w1,w2)=(2/225,1/36,2/45)` per configuration. The candidate reads those
weights through the public preparation producer and forms

```text
a = log(w1/w0) = log(25/8),
b = log(w2/w1) = log(8/5),
J = log(w1²/(w0 w2)) = log(125/64),
Omega = log(w2/w0) = log(5).
```

Let `F` be the native covariant family feature, in its existing orthogonal
source frame:

```text
F = diag(41/225 - 7 sqrt(5)/750,
         41/225 + 7 sqrt(5)/750,
         61/450).
```

The **chosen** trial matrices are `Y_u=aF`, `Y_d=bF`, `Y_e=JF`. They are
real and symmetric. Their family anisotropy follows from the existing
protected-core reference; selecting these three multipliers and assigning
the feature to Yukawa channels do not follow from compression alone.

The same declared Higgs doublet supplies all three exact invariant lines:
`Q u^c H`, `Q d^c H*`, `L e^c H*`. For `H=(0,v)`, `v` is the neutral
component norm, conventionally `v_SM/sqrt(2)`. The default architecture
chooses `v=Omega` in the coherent parameter's inverse units. A different
real nonnegative `v` is an explicit supplied background, rather than a
derived physical vacuum. The `v=0` control turns all these mass terms off.

Writing the retained symmetric internal component tensors as `C_f,k`,
the resulting all-left internal/family coefficient is

```text
K45(v) = v sum_f C_f,neutral tensor Y_f,
K45^T = K45,
Q45^T K45 + K45 Q45 = 0.
```

This is not an identity mass across unrelated gauge representations. The
neutral tensor pairs the oppositely charged up, down and charged-lepton
slots and leaves the neutrino slots empty. An identity mass on this all-left
table would violate electromagnetic charge conservation. The code rechecks
72 one-generation and 72 native-family current/Higgs covariance identities.
They include the Higgs variation; the frozen neutral vacuum does not
commute with every broken electroweak generator.

## Weyl, CAR and Nambu conventions

Let `epsilon=i sigma2=[[0,1],[-1,0]]`. The physical pairing tensor is
`epsilon tensor K45`, not the symmetric internal matrix by itself. It is
antisymmetric, as required for a nonzero CAR pairing of two left fields.
The declared left-Weyl rotation and boost matrices obey
`j^T epsilon+epsilon j=0`. Thus the local bilinear is a Lorentz scalar under
that declared Weyl action.

Use the sign convention

```text
L_pair = +(1/2) psi^T (epsilon tensor K45) psi + h.c.,
H_pair = +(1/2) sum_ij Delta_ij c_i† c_j† + h.c.,
Delta = epsilon tensor conjugate(K45).
```

The last identity uses `H_mass=-L_pair` and the order reversal of CAR
operators under conjugation. For the executed real neutral candidate,
`Delta=epsilon tensor K45`; this is what `pairing_kernel()` returns.
For a future complex Higgs/action kernel, conjugation must be retained.
The ordinary Weyl kinetic block and Nambu hole block belong to the separate
kinetic construction. The Nambu order `(c(p),c†(-p))` encodes particle-hole
redundancy; it does not add physical right-handed mirror multiplets.

`matrices()` returns immutable exact SymPy sparse matrices. Public helpers
include `generators`, `identification_operators`, `electric_charge`,
`mass_kernel`, `pairing_kernel`, `family_feature`, and `yukawa_matrices`.
`spin_electric_charge()` gives `I2 tensor Q45`, and `mode_labels()` labels
all ninety annihilator modes before Nambu doubling. In the retained order,
the three spin-up negative-charge lepton modes are indices `33,34,35`;
their nonzero pairing partners are `72,73,74`. The three other lepton
spin blocks, all thirty-six quark spin blocks and all six neutrino modes
remain in the full candidate.
`report()` exports the order, identities, source hashes and explicitly
conditional tree parameters. No compiled spatial chiral regulator or
nonperturbative quantum measure is supplied by these local tensor identities.

## Exact consistency and a decisive hierarchy limitation

The code checks the internal mixed/gravitational hypercharge traces,
hypercharge cubic trace, color Cartan cubic polynomial and even weak-doublet
count. These are the representation anomaly checks of the declared chiral
table, not a construction of a gauge-invariant interacting regulator measure.
There are four weak doublets per generation and twelve in the three-family
candidate. The neutral kernel consists of disjoint charged pairing blocks,
has rank 42 for positive `v`, and leaves three neutrino components massless.

Under the declared canonical kinetic interpretation, each charged sector's
tree mass parameters are `v times coefficient times eigenvalue(F)`. At the
default `v=Omega`, their dimensionless readout targets are:

| Channel | Tree mass parameters divided by `Omega`, in increasing order |
| --- | --- |
| Up type | `0.1544566473, 0.1838502904, 0.2314102039` |
| Down type | `0.0637116031, 0.0758361452, 0.0954540664` |
| Charged lepton | `0.0907450442, 0.1080141452, 0.1359561376` |

These are computed targets of this trial action, not measured Standard Model
pole masses. Every family ratio is the same, with maximum/minimum strictly
below `3/2`. Moreover `a/b<5/2`, because the equivalent rational inequality
`(25/8)²<(8/5)⁵` holds. The largest/smallest of **all** charged tree parameters
is consequently below `15/4` for every common positive `v`. No scale change
or relabeling can turn this narrow trial spectrum into the observed charged
fermion hierarchy. No fit is attempted.

All three Yukawa matrices are real multiples of the same `F`. Their left
family diagonalizations coincide; a canonical tree flavor interpretation
therefore gives the identity CKM matrix up to phases and permutations, with
no flavor CP phase. Neutrinos remain massless in this minimal attachment;
their degenerate sector supplies no physical PMNS prediction. These are
specific limitations of this coupling architecture, not a prohibition on
other source-derived interactions or on the broader OPH construction.

The candidate makes the previously separate representation, symmetry breaker
and allowed mass interaction into one exact local action. Propagation and
observer measurements can now test that action. Its local consistency does
not select the scalar vacuum, establish the physical regulator/refinement,
or make its tree parameters observed masses.

```sh
python -m model.common_source.chiral_candidate.action
python -m pytest -q model/common_source/chiral_candidate/test_action.py
```
