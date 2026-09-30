# Independent executable audit of the #628 proposal-law boundary

Date: 2026-09-30
Baseline: `FloatingPragma/observer-patch-holography` `upstream/main` at
`0c914ccd1c5ce16c7ef3058808a36c56ae9405bd`.
Branch: `arithmon/source-repair-gram-selection-0`.

## Contribution and finding

This packet does not claim the conditional uniformity theorem as a new
discovery. Upstream already has `EqualSeamSelection.lean`,
`SeamCurrentHomogeneousAction.lean`, and `docs/CANONICAL_REPAIR_LAW_RFC.md`.
They establish the conditional A2/A3 route from a supplied complete move
simplex, natural objective and unique minimizer to uniform source counting;
the RFC states that its A1-R/A2-R strengthening is proposed and not adopted.
The contribution here is an independent executable audit of the exact
#628-to-proposal boundary, the positive biased countermodel for the #628
relation alone, and receipt-checked scope for the connection to the frozen
Galois pair.

The #628 repair relation fixes admissible local conservative unit transfers,
but does not assign probabilities to the thirty seams. Distinct positive
seam laws therefore remain possible for that relation. The conditional
mathematics is exact: on a single A5 edge orbit, an A5-invariant symmetric
local conservative generator has one rate parameter and is
`L = c(5I-A)`, `c >= 0`; if `c > 0`, the certified rank-three `P3` band is
uniquely slowest and the frozen positive Gram table is `G_plus = 4 P3`.

The existing, narrower scheduler result is checked from its inputs at
runtime. #614 selects `1/30` per seam by A3 relative-entropy projection
**from an explicitly uniform move
reference** on its named reference carrier. The #628/#614 directed-seam
receipt composes the resulting uniform directed orbit with opposite balanced
placements and proves `I - Delta/60` on the total-load-one sector. It does not
prove a path law or bind that S=1 channel to the complete nonlinear integer
repair chain as a linear generator. Thus it is evidence for the conditional
selector and a named sector, not a derivation of the complete missing bridge
from the #628 relation alone.

## Exact graph and rate classification

The producer derives edges from the twenty oriented triangles in
`echosahedral_federation_reference.json`; the serialized edge list is used as
a consistency control only. It enumerates orientation-preserving
automorphisms, reconstructs vertex and unoriented-edge orbits, and solves the
edge-weight invariance equations over `Q` by exact row reduction. The separate
verifier uses a different minimum-remaining-values automorphism search and
independent rational elimination.

| Quantity | Recomputed result |
|---|---:|
| Vertices / seams / degree | 12 / 30 / 5 |
| Connected | yes |
| Proper oriented automorphism group | 60 |
| Vertex orbit / unoriented seam orbit | 12 / 30 |
| Dimension of A5-invariant edge weights | 1 |
| Invariant basis | all-ones vector on the 30 seams |

For any symmetric local proposal weights `w_ij >= 0` supported on seams,
`(L_w f)(i) = sum_{j~i} w_ij(f(i)-f(j))` and `L_w=D_w-W`. Its row sums
vanish, constants lie in its kernel, off-diagonal support is exactly the
seam graph when all weights are positive, and
`f^T L_w f = sum_{ij in E} w_ij(f_i-f_j)^2 >= 0`. The exact countermodel sets
one lexicographically first seam to rate 2 and the other 29 to rate 1. It
preserves the same admissible repair relation and conservation, remains
strictly positive, and is not a scalar multiple of `5I-A`.

## Source receipts and custody

The existing uniform scheduler packet is
`code/a5_closure/manifests/a3_scheduler_kernel_reference.json` (declared
manifest digest `sha256:89a9b7302e1ed9a5bad5542d7eedb68dda76d75de3064fdca22c966703c1a479`;
raw file SHA-256 `85435a4c36bccace0e81fb0bb9476444641a63750067557bb79a006cec8879eb`).
Its producer is `code/a5_closure/a3_scheduler_kernel_certificate.py` (SHA-256
`ffadea01ecfa8f35016c07c077f3eeef18035d99ed90ae1a4ec8d877dab27510`).
The directed repair receipt
`code/a5_closure/manifests/directed_seam_repair_reference.json` pins this
scheduler by canonical digest
`sha256:4afffa0760ff59649d1f2548d8f755f1efe4c5be14594bbde80f22d7ed3d536a`,
and pins #628 by canonical digest
`sha256:efad022da8fdb52b58e6e8356ef945227832467c162cb512d0261812f286a655`.
Its declared scope is a one-step expected scalar channel plus a local integer
macro bridge; its `uniform_s1_channel` is explicitly for sector `S=1`.

`EqualSeamSelection.lean` proves uniform weights from a supplied
`NaturalUniqueMoveProjection` (unique natural minimizer and transitive move
action); its module commentary states the canonical A1--A3 basis does not
provide the complete A1-R/A2-R fields or a source-emitted unit-counting move
reference. The #614 packet does provide a named uniform reference and A3
selection for its scheduler. These are distinct source scopes.

The existing spectral bridge is conditional and already formalized:
`PortGramRepairCovariance.normalizedKernel_tendsto_portGram` uses the declared
step `T = I-(5I-A)/60`, and `PortGramRepairBand` certifies `G_plus=4P3`,
`G_minus=4P3'`, with costs `5-sqrt(5), 6, 5+sqrt(5)`. This audit does not
identify #628 settling cost, its nonlinear transition chain, the S=1 receipt,
the A3 Hessian, or physical time with that declared operator. It makes no
A1 support-attachment claim; the registered Galois boundary still records
`PORT-GRAM-SUPPORT-ATTACHMENT` as absent and PR-53 open.

## Verdict and limits

**Verdict:** `A5_INVARIANT_REPAIR_LAW_FORCES_LOW_BAND` conditionally, together
with `REPAIR_RELATION_DOES_NOT_SELECT_PROPOSAL_LAW` for the #628 relation
alone. The #614 A3 scheduler selects the uniform law only relative to its
explicit uniform move reference. The source-level missing premise for a
full-state native dynamics claim is a source-bound identification of that
proposal/linearized generator with the complete #628 repair process; A5
invariance by itself classifies the weights but does not provide this binding.

```text
#628 fixes the admissible repair relation.
#614 fixes the uniform one-step law relative to its declared uniform reference.
The directed receipt binds the two on the S=1 expected channel.
No source theorem currently identifies that channel with the complete nonlinear
#628 integer repair process.
```

If that binding is supplied, the existing spectral chain selects `P3`, hence
`G_plus` inside the frozen Galois pair.

| Link | Status |
|---|---|
| Face incidence to 12/30/degree-five graph and proper A5 edge transitivity | PROVED (independent executable reconstruction) |
| A5-invariant symmetric local weights are one-dimensional | PROVED (exact constraint rank 29) |
| #628 relation itself selects proposal probabilities | NOT PRESENT; exact positive biased counterlaw |
| #614 A3 scheduler law on its declared uniform reference | PROVED by pinned existing receipt |
| Full integer-state process equals the declared linear generator | NOT PRESENT |
| Uniform generator has unique slow `P3`; `G_plus=4P3` in frozen pair | CONDITIONAL; existing Lean theorem chain |
| General-rate covariance limit and independent Galois negative replay | NOT ATTEMPTED |
| A1 support attachment / PR-53 | NOT PRESENT / open |

The standalone tests exercise incidence, independent orbit/rank and generator
reconstruction, receipt digests and parent pins, and actual mutations of
faces, serialized seams, the #614 probability and directed-receipt scope
flags. They reject deleted or duplicated faces/edges, a reversed oriented
face, a distance-two seam, and forged receipt claims. The full twenty-control
hostile suite, all-mode normalized-kernel replay, and a Lean formalization
were not attempted. No registries, papers, Lean files, generated surfaces or
frozen prediction files were changed.

## Custody audit (candidates were read, none edited)

| Candidate | Classification | Evidence |
|---|---|---|
| `PortFrameGram.lean` | BYTE-PINNED-PARENT | audited snapshot SHA-256 `d1cebe56450e7586eed730b52068753ebca3a6563da1453679e87fa7c7b653e3`; snapshot checked by `verify_source_current_order_sensitive_inventory.py` |
| `A5FamilyBand.lean` | BYTE-PINNED-PARENT | audited snapshot SHA-256 `63b425c890b49f49a53858b5a480b993cfee81d92978fa3ed5abf665e48645a6`; also referenced by `sm_fermion_current` receipts |
| `PortGramRepairBand.lean` | BYTE-PINNED-PARENT | audited snapshot SHA-256 `75286414b7c33492b40225dd48ca9321cf3a09ecf96b65e254af9bd02421cf72`; read only, unchanged |
| `PortGramRepairCovariance.lean` | BYTE-PINNED-PARENT | audited snapshot SHA-256 `8d96c8a01361a9cfc2f53901d73032041480c93bb50190ab12ac939a7e9a6ad8`; read only, unchanged |
| `A5PortAction.lean` | ACTIVE-INVENTORY-ONLY | present in `claims/active_surface_inventory.json`; no byte pin found in the inspected source-current snapshot |
| `record_counting_mechanism_reference.json` | BYTE-PINNED-PARENT | canonical hash consumed by the directed seam receipt; raw hash also appears in the source-order inventory |
| `record_counting_mechanism_certificate.py` | GENERATED-SURFACE | producer for the pinned #628 manifest; imported by the directed seam producer |
| claim registries, selection ledger, and papers | FROZEN/PREDICTION for this discovery scope | explicitly excluded from edits by runbook |

`BYTE-PINNED-PARENT_CHANGED: NO`.

## Files added and checks

- `source_repair_generator_certificate.py`
- `verify_source_repair_generator_independent.py`
- `tests/test_source_repair_generator_certificate.py`
- `source_repair_generator_certificate.md`

Commands run successfully:

```text
python3 -B code/a5_closure/source_repair_generator_certificate.py
python3 -B code/a5_closure/verify_source_repair_generator_independent.py
python3 -B code/a5_closure/tests/test_source_repair_generator_certificate.py
```

Five test methods passed, including the mutation controls listed above. No
Lean build or CI run was performed.
