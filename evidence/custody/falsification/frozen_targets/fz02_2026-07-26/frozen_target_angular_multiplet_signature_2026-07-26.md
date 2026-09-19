# FZ-02: the angular multiplet signature, frozen target

Status: FROZEN at registration. Registered before any comparison data is
examined, under the corpus kill-condition protocol. Milestone M1 deliverable
(c) of the completion plan.

## Target statement

<!-- FZ02-TARGET-BEGIN -->
For any physical system that realizes the OPH twelve-port screen with an
A5-invariant dynamics, the rotation multiplets branch under A5 exactly as:

- l = 0: 1
- l = 1: 3
- l = 2: 5 (irreducible: the quintet does not split)
- l = 3: 3' + 4
- l = 4: 4 + 5
- l = 5: 3 + 3' + 5
- l = 6: 1 + 3 + 4 + 5

By Schur's lemma an A5-invariant self-adjoint perturbation is scalar on each
block, so the degeneracy pattern above is forced while block eigenvalues are
dynamical. Since l = 2 is irreducible, every A5-invariant ensemble mean at
l = 2 vanishes and every invariant ensemble covariance at l = 2 is scalar; a
single realization may be anisotropic or aligned. The first nonconstant
A5 invariant sits at l = 6, and its frame agrees with the l = 3 block frame.

Companion face-phase statement: the multiplicity of a nontrivial face phase
on the irreducibles (1, 3, 3', 4, 5) is (0, 1, 1, 1, 2), so the
dimension-minimal irreducible extension of a nontrivial face phase is
three-dimensional (3 or 3'). If the quotient-visible family fiber is that
extension, the generation count is N_g = 3 with one rephasing-invariant CKM
phase slot.

Scope: this registration binds the integer-l (rotational, SO(3)-restricted)
branching only. Half-integer branching through the binary icosahedral cover
is outside this registration.

Machine checks frozen with this target:
- executable: a5_angular_multiplet_reference.receipt.json (schema
  oph.a5_angular_multiplet_receipt.v1, receipt_sha256
  d95afcdef548e51f18825053ea74bc23d124672481194cff7282aee91b306b82)
- Lean: A5AngularMultiplets_pinned.lean (theorems branching_table,
  character_table_orthonormal, first_nonconstant_invariant_at_six,
  face_phase_multiplicities; standard axioms only)
<!-- FZ02-TARGET-END -->

## Decision policy (frozen with the target)

Rows and triggers. Verdict vocabulary: COMPATIBLE, FAIL, INCONCLUSIVE.
Every evaluation requires a dataset declaration, estimator, and tolerance
registered before the data is examined; a kill is published with the same
prominence as a confirmation after the audit / re-audit / repair ladder.

- FZ02-R01 (degeneracy pattern, READY): in a declared A5-invariant
  realization resolved to at least the l = 6 band, the multiplet dimension
  pattern per level must equal the frozen table. A resolved splitting of
  the l = 2 quintet into sub-multiplets under A5-preserving conditions, at
  or above the registered significance threshold of the evaluation, is a
  FAIL for the signature.
- FZ02-R02 (invariant statistics, READY): over a declared A5-invariant
  ensemble, the l = 2 covariance must be scalar. A nonscalar invariant
  covariance at three or more registered sigma is a FAIL.
- FZ02-R03 (first invariant and frame lock, READY): a nonconstant
  A5-invariant detected below l = 6 at three or more registered sigma is a
  FAIL; a measured l = 6 invariant whose frame misaligns with the l = 3
  block frame beyond the registered tolerance is a FAIL.
- FZ02-R04 (face-phase family fiber, THEOREM-TARGET): this row binds the
  family-attachment lane (issue #569), not the framework. If that lane
  produces a quotient-visible family fiber that is not the minimal
  three-dimensional extension (3 or 3'), the N_g = 3 route through the
  face phase is killed and the kill is published.
- No verdict is issued at insufficient resolution: an evaluation that
  cannot resolve the relevant band returns INCONCLUSIVE and does not count
  as a confirmation.

## Integrity

The registered content is the fenced block between FZ02-TARGET-BEGIN and
FZ02-TARGET-END. Independent verification: recompute sha256 over the block
bytes and compare with the registration manifest; verify the .ots
attestations with `ots verify` against the manifest entries; recompute the
executable receipt with
`python3 code/a5_closure/a5_harmonic_decomposition.py --receipt <out>` at
the pinned source commit and compare receipt_sha256.
