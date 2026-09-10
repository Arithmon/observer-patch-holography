# Coset carrier certificate

## Result

The A1 twelve-port carrier is the coset geometry of the central quotient of a
supplied binary icosahedral source, and an exhaustive scan of subgroup
placements identifies exactly which relative placements reconstruct it.

**Theorem (conditional on the supplied group).** Let `S = SL(2, F_5)`, taken as
the exact finite representative of the binary icosahedral group `2I`, and let
`G = S/Z(S)`. For cyclic subgroups `C5`, `C3`, `C2` of `G`, form the coset
geometry with vertices `G/C5`, faces `G/C3`, edges `G/C2`, and incidence
`gH ~ kK` if and only if `gH ∩ kK` is nonempty.

(a) All `6 · 10 · 15 = 900` placements have twelve vertices, twenty faces, and
thirty edges, with five faces and five edges at each vertex, three vertices and
three edges on each face, and two vertices and two faces on each edge.

(b) Exactly 120 placements give a closed, coherently oriented triangulated
surface. They are the placements admitting generators `x, y, z` of orders
2, 3, 5 with `xyz = 1`; for each of them that generator choice is unique, and
they correspond bijectively to the solutions of `a^2 = b^3 = c^5 = abc` in `S`.

(c) The 120 compatible placements form two orbits of 60 under inner
automorphisms of `G` and a single orbit under `Aut(G)`, which has order 120.

(d) Every compatible placement relabels explicitly onto the committed packet:
the coset edge graph maps onto `PortFrameGram.neighbors`, the coset orientation
maps onto `CoreAxioms.orientedFaces`, the distance-three partner maps onto the
antipode `i -> 11 - i`, and the transported action of `G` is exactly the set of
sixty rows of `A5PortAction.perms`.

## Construction

The source manifest `manifests/coset_carrier_reference.json` declares the
family and prime of the source group and names the three committed Lean
declarations used for comparison. A source firewall rejects any source field
that names target structure. The verifier builds `S` and checks its order, its
centre `{+I, -I}`, its unique involution, its perfectness, and the binary
icosahedral element-order profile `1, 1, 20, 30, 24, 20, 24` on the orders
`1, 2, 3, 4, 5, 6, 10`. It forms `G`, checks the class sizes
`1, 12, 12, 15, 20` and simplicity, and checks that the preimages of the
cyclic subgroups of orders 5, 3, 2 are cyclic of orders 10, 6, 4.

Nonempty coset intersection is computed through a shared representative: the
incident pairs of each relation are exactly `{(xH, xK) : x in G}`. The coset
orientation is read from the chambers `(xC5, xC2, xC3)` that share one
representative. The committed adjacency, port rows, and oriented faces are
parsed directly from `Lean/Screen/PortFrameGram.lean`,
`Lean/Screen/A5PortAction.lean`, and
`Lean/ObserverPatchHolography/CoreAxioms.lean`; parsing fails closed, and the
receipt records a hash of each parsed table.

## Placement classification

Each placement is checked in a fixed order, and the first failure names its
class.

| Class | Placements | Per `(C5, C3)` pair | Characterization |
|---|---|---|---|
| `CARRIER` | 120 | 2 | generators with `xyz = 1` exist |
| `FACE_BOUNDARY_MISMATCH` | 180 | 3 | edge graph and face vertex sets agree with the committed packet; the coset edge-face incidence disagrees |
| `FACES_NOT_TRIANGLES` | 300 | 5 | no coset face is a triangle of the coset edge graph |
| `COLLAPSED_EDGES` | 300 | 5 | `C2` lies in the normalizer of `C5`; the thirty edges fall on six vertex pairs, five each |

Every one of the sixty `(C5, C3)` pairs splits its fifteen involutions the
same way, so the placement of `C2` carries the whole selection.

The second row is the reason the certificate checks all three incidences. For
those 180 placements, vertex-level data (the edge graph and the vertex sets of
the faces) coincide with the committed carrier, and only the edge-face
relation of the coset geometry differs. In port labels, the coset edge
`{0, 6}` is coset-incident to the faces `{0, 1, 2}` and `{6, 9, 10}`, while the
faces containing it are `{0, 3, 6}` and `{0, 4, 6}`. A check restricted to
vertex-level data accepts 300 placements; the full check accepts 120.

## Choice accounting

The orders five, three, and two fix the coset sizes and every local incidence
count, identically for all 900 placements, and they do not determine the global
incidence.

The relative placement is unique up to automorphisms of the supplied group:
the 120 compatible placements form one orbit under `Aut(G)`. Declaring the
source by the presentation `a^2 = b^3 = c^5 = abc` supplies a compatible
placement with no further choice.

Twelve ports select the stabilizer of order five, since `G/C5` has twelve
points, `G/C3` twenty, and `G/C2` thirty.

The construction outputs a coherent orientation. For every compatible
placement, 60 of the 120 graph relabellings carry it onto
`CoreAxioms.orientedFaces`, and the other 60 are exactly those relabellings
composed with the antipode, which carry it onto the fully reversed packet.

The supplied group is the residual premise. The certificate takes `S` as given
and makes no selection among candidate source groups.

## Negative controls

`negative_controls/coset_carrier_negative_controls.json` records fifteen typed
controls, each required to fail with its stated code:

| Control | Code |
|---|---|
| declared prime 3, declared prime 7 | `SOURCE_ORDER` |
| declared family `GL2` | `SOURCE_GROUP` |
| source field naming target structure | `SOURCE_FIREWALL` |
| extra source field | `SCHEMA_FIELDS` |
| comparison declaration swapped | `FIXTURE_DECLARATION` |
| representative placement of each failing class | `COLLAPSED_EDGES`, `FACES_NOT_TRIANGLES`, `FACE_BOUNDARY_MISMATCH` |
| carrier with one edge moved to a face that does not contain it | `FACE_BOUNDARY_MISMATCH` |
| carrier with one edge removed from one face | `FACE_SIDES` |
| carrier with one chamber moved to another face | `ORIENTATION_INCOHERENT` |
| committed adjacency with a degree-preserving double edge swap | `RELABEL_ADJACENCY` |
| committed packet with one face reversed | `RELABEL_ORIENTATION` |
| committed port rows with one row perturbed | `PORT_ACTION_MISMATCH` |

The bundle also records witnesses for the three failing classes and the
positive control in which every committed face is reversed.

## Reproduction

```bash
python3 code/a5_closure/coset_carrier_certificate.py all
python3 code/a5_closure/coset_carrier_certificate.py verify \
  --manifest code/a5_closure/manifests/coset_carrier_reference.json \
  --receipt code/a5_closure/receipts/coset_carrier_reference.receipt.json
python3 -m pytest -q code/a5_closure/tests/test_coset_carrier_certificate.py
```

The verifier uses the Python standard library and exact integer arithmetic.
The receipt and the control bundle are byte-identical across runs.

## Claim boundary

The result is conditional on the supplied group `SL(2, F_5)`, taken as the
exact finite representative of `2I`; that identification is classical and is
cited without formalization. The reconstruction factors through the central
quotient, so it cannot distinguish the source from its quotient. It does not
select the source group, identify ports with physical objects, or alter the
declared status of the A1 boundary packet. Selection-ledger row 4 keeps its
class, menu, and compression accounting: the declaration moves from the boundary
complex to the supplied group. The binary icosahedral double cover derived
downstream from the port frame is independent of this certificate, and this
certificate does not turn that derivation into a premise.
