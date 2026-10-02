# OPH Port-Gram support attachment

This package compares the local oriented boundary of an observer-like self-reading patch with the spherical support of its federation. Each patch has local state, twelve ports, readback records and feedback or repair operations; the public evidence consists of an exact certificate and independent replay. Its scope is finite source attachment, without a physical orientation, scale or metric identification.

## Verdict

Primary: **`SUPPORT_ATTACHMENT_REMAINS_UNSELECTED`**.

The local oriented boundary and the global degree-one support both exist. Their seed presentations are literally identical as ordered oriented face lists. The source packet nevertheless declares the local twelve-port boundary of each carrier independently of the federation nerve. It supplies no source map from port labels to support vertices. Treating the shared integer presentation as that map produces 60 valid seed attachments, not a selected one.

Secondary: **`SEED_ATTACHMENT_FAMILY_CLASSIFIED`** and **`SEED_ATTACHMENT_SET_IS_A5_TORSOR`**. These classify the raw oriented seed maps only. Every seed candidate pushes the local fundamental cycle to the designated support cycle as an exact integer chain. No source-provided common local refinement tower or local-to-global oriented homotopy was found in the pinned audit corpus.

## Reproduction and source scope

From the repository root:

```bash
python3 code/a5_closure/port_gram_support_attachment_certificate.py --json
python3 -O code/a5_closure/verify_port_gram_support_attachment_independent.py
python3 -m pytest -q code/a5_closure/tests/test_port_gram_support_attachment.py
```

The producer emits the complete 60-member seed family and its exact chain checks. The independent verifier reconstructs the finite result without importing the producer or reading that candidate list, and invokes the separately pinned independent Galois verifier for the radial degrees. It verifies the mathematical result rather than accepting an externally supplied attachment report.

Both implementations pin the twelve source files listed in their `PINNED_SOURCE_SHA256` constants. Those bytes include the local boundary, source geometry and derivation, source-realization topology and receipt, carrier-map boundary, and Galois evidence. The attachment-absence statement applies only to this explicitly listed corpus. Changed bytes fail verification and require a renewed source analysis; a matching hash fixes the analyzed text but does not itself prove that no attachment is declared.

## Reconstructed source data

### Local boundary

`CoreAxioms.orientedFaces` has 12 vertices, 20 ordered faces and 30 edges. Its 20 oriented triangles define the fundamental integer 2-cycle; each undirected edge occurs twice with opposite induced directions. The Lean declarations verify the face count, edge count and coherent edge orientation. The exact ordered triples are parsed by both new checkers.

### Global support and source tower

The three presentations agree literally, entry by entry and in order:

```text
CoreAxioms.orientedFaces
= source_selection_model.geometry.seed()
= level-0 support faces in m1_source_realization/receipt.json
```

Result: **`SUPPORT_SEED_LITERAL_ORIENTED_MATCH`**. The source realization receipt contains the four actual stages checked here:

| Level | Support vertices | Oriented support faces | Nerve carriers | `b_* z_r` |
|---:|---:|---:|---:|---|
| 0 | 12 | 20 | 13 | `[S_0]` exactly |
| 1 | 42 | 80 | 44 | `[S_1]` exactly |
| 2 | 162 | 320 | 165 | `[S_2]` exactly |
| 3 | 642 | 1,280 | 646 | `[S_3]` exactly |

At level zero, the fibre sizes are 2 over support vertex 0 and 1 over each other support vertex. The designated section selects slot zero. The source construction gives `z_r=s_*[S_r]`, and the bridge sends it to `[S_r]`; the new producer and independent verifier recompute the section/coarsening squares and integer chains across all four stages.

| Object | Current status |
|---|---|
| `K_i`, oriented boundary | DECLARED in Lean; finite incidence and orientation are formally proved there. |
| `N_r`, `S_r`, support subdivision | SOURCE-CONSTRUCTED in the source realization; EXECUTABLE-ONLY receipt/checker, not Lean-formalized. |
| `b_r`, section and `z_r` | SOURCE-CONSTRUCTED as the nerve projection, selected section and its pushed cycle; exact finite chain checks replay. |
| `iota_r` and oriented support sphere | Declared/source-level support realization; no local-ray coordinate comparison is supplied by the receipt. |
| Local radial maps and degrees | Consumed from the byte-pinned Galois reference; the existing Galois producer and separate independent verifier are replayed in independent paths. |
| Local barycentric refinement | Lean-formalized same-parent coordinate scaling within a committed face. |
| Global refinement/coarsening | SOURCE-CONSTRUCTED midpoint subdivision and fibre coarsening; executable chain checks at levels 0–3. |

The local Lean file states that it has no frequency-`n` incidence complex and does not prove compatibility with a committed refined complex. The support tower uses midpoint subdivisions and fibre maps. Consequently the seed candidate maps do not extend through the current data as a natural family: **`ATTACHMENT_ONLY_AT_SEED_LEVEL`**; extension requires a new choice or theorem.

### Source binding and provenance

Status: **`NO_ATTACHMENT_FOUND_IN_PINNED_AUDITED_CORPUS`**. The audited corpus is explicitly byte-pinned by the certificate.

Evidence: the source-clock completion constructs a nerve over the support, chooses a section, and transports the selected program host under presentation changes. It also says each primitive carrier has its own twelve-port boundary, independently of the federation nerve. No declaration or theorem in the pinned corpus maps a local port label to a support vertex. The identical integers `0..11` and identical triples establish presentation equality; they do not establish object identity or a natural transformation.

No branch-specific Galois datum, shortest-chord ranking, trace, repair selection, dynamic repair-law selector, or empirical target was used to generate or rank the attachments. The full 60-member family is emitted by `port_gram_support_attachment_certificate.py --json`.

## Exact seed attachment family

Restrict candidate identifications to bijections preserving the oriented seed and the designated section image. Independent permutation searches give the raw seed family:

```text
Aut+(K_i)                         = 60
admissible seed identifications  = 60
stabilizer                       = 1
orbits                           = 1
orbit size                       = 60
```

Thus the seed attachment set is an **`A5` torsor**. For every one of its 60 elements `h_0`, the integer chain pushforward is exactly

```text
(h_0)_*[K_i] = z_0
```

It is not merely equality in homology. The proper icosahedral relabelings act freely and transitively on this raw seed family. The fixed level-0 support row carries extra decoration: fibre sizes `[2,1,1,1,1,1,1,1,1,1,1,1]` and process link `[[0,1]]` mark support vertex 0. Its orientation-preserving automorphism subgroup has order 5 (`C5`), so quotienting the 60 raw maps by that fixed-mark subgroup gives 12 orbits; the marked host has an orbit of size 12 under `A5`. Presentation changes transport the mark, but the corpus does not declare an action identifying that transport with local port maps. Therefore this does not classify source-relative attachment equivalence. The verdict is intentionally scoped to the raw seed family. Reversing orientation leaves the unoriented graph data intact but fails the oriented chain condition.

## Geometric comparison and branch test

Local degrees consumed from the pinned Galois reference (and independently replayed by the existing Galois verifier):

```text
deg(radial_plus)  = 1
deg(radial_minus) = 7
global support chain degree = 1
```

The degree-one global support is independently present, but no source-bound map identifies its realized sphere with the local radial target. Therefore the oriented homotopy comparison is **not established** for any selected `h_r`.

Conditional statement only: if an attachment is chosen and both maps are verified continuous maps from the same oriented `S²` to the same oriented `S²` with the shared orientation convention, degree classifies their homotopy classes. Under those additional premises the global degree-one map and `radial_plus` have the same degree, while `radial_minus` has degree 7 and cannot be homotopic to it. The source has not supplied the attachment and common realization premises, so this is **`CONDITIONAL_SUPPORT_SELECTOR`**, not a source theorem and not `G_PLUS_SELECTED_BY_SOURCE_SUPPORT_DEGREE`.

The result is not an incompatibility obstruction: no seed candidate fails its cycle test. The missing datum is a source-selected port-to-support map, together with its refinement extension and realized oriented homotopy.

## Hostile controls

The tests cover the following controls. Local/support/both proper `A5` relabelings, orientation reversals, label swaps, and cycle mutations enter the pure injectable seed-classification pipeline. Remaining mutations fail their corresponding chain, pin, tower, branch-coverage, or policy gate; they are scoped mutation checks and do not all invoke the full report builder.

1. Reverse one local face orientation.
2. Reverse one support face orientation.
3. Relabel only the local frame by a proper `A5` element.
4. Relabel only support by a proper `A5` element.
5. Relabel both coherently.
6. Apply an orientation-reversing icosahedral symmetry.
7. Swap two support labels without updating faces.
8. Replace `z_r` by `-z_r`.
9. Forge support degree 1 as 7.
10. Delete a support face.
11. Alter a coarsening image.
12. Alter a section choice.
13. Break a refinement square.
14. Omit the minus branch.
15. Import `G_plus` as attachment target data.
16. Import the dynamic repair-law selector.
17. Promote integer label coincidence to source provenance.
18. Substitute graph isomorphism for oriented chain equality.
19. Claim homotopy from degree before checking the sphere/orientation premises.
20. Claim PR-53 discharge.

## Interpretation boundaries

| Claim | Result | Additional requirement |
|---|---|---|
| Local and global seeds are literally the same oriented presentation. | ESTABLISHED | Exact list equality and edge/cycle census. |
| Pinned audited corpus chooses a local-port to support map. | NO_ATTACHMENT_FOUND_IN_PINNED_AUDITED_CORPUS | Add a source theorem or source datum that selects the map without downstream branch data. |
| 60 seed maps push the local cycle to `z_0`. | ESTABLISHED as candidate mathematics | All maps explicitly enumerated; exact chain equality checked. |
| Candidate map extends naturally through current local and global refinements. | NOT PRESENT | Supply a common local incidence tower and prove the commuting squares. |
| Local radial map is homotopic to the realized global support map. | NOT PRESENT | First bind `h_r`, then supply the same oriented-sphere realization and continuity/homotopy proof. |
| Degree one selects `G_plus` in the frozen pair. | CONDITIONAL only | Requires the missing source-bound attachment premises. |
| PR-53 is discharged. | NO; unchanged | No physical direction, metric or scale follows from this finite support audit. |
