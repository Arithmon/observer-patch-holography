# FZ-02 custody erratum

This is an append-only provenance correction for the FZ-02 custody bundle.
It does not alter the frozen target, its decision policy, its kill bands, the
executable receipt, the pinned Lean module, the original registration
manifest, or any existing OpenTimestamps proof.

## Corrections

1. The `frozen_utc` value `2026-07-26T09:30:00Z` embedded in
   `registration_manifest_2026-07-26.json`,
   `fz02_freeze_registry_2026-07-26.json`, and the downstream generated
   register is a timestamp metadata error. The corrected custody record is
   the oph-meta commit
   `1e7d7c73dadeef9aa10ec60061a85cee8426c5b1`, committed at
   `2026-07-26T06:41:53Z`. This commit binds the complete artifact set below.
2. The source commit containing the frozen executable receipt and
   `Lean/Screen/A5AngularMultiplets.lean` is
   `091658ce585c107a260e7b980352be904d2419b2`, committed at
   `2026-07-26T06:41:50Z`. The `36e44fcc` value in the original custody
   metadata is its parent and predates both of those source artifacts.
3. The original registration manifest hashes the complete frozen-target
   Markdown file. Its integrity paragraph incorrectly instructs readers to
   compare the hash of only the fenced target block with that whole-file
   manifest entry. The whole-file and fenced-block hashes are distinct and
   are recorded explicitly below.

## Bound original artifacts

The original registration manifest has SHA-256
`278dedcabe9b242b21a1f0ea0d0ceb8899dd279acf9b579eb0c5bae8af858ddf`.
Its artifact entries have been independently recomputed:

| Artifact | SHA-256 scope | SHA-256 |
|---|---|---|
| `frozen_target_angular_multiplet_signature_2026-07-26.md` | complete file | `abff4289459a53be63808210b37fd4025af2abf73bc82a2b2a41594953742a57` |
| `fz02_freeze_registry_2026-07-26.json` | complete file | `bf755a072641eeff3ac000f7abd15f4eb0360175d8fa738a73c31495f98447c4` |
| `a5_angular_multiplet_reference.receipt.json` | complete file | `93d72fb361d3522145aab7bbf21edabcc4a2d960f87e66629ef831c4662a9e89` |
| `A5AngularMultiplets_pinned.lean` | complete file | `e18f8874ea1fdf3a0a74289951a5f6d8bce5763eefaba0e3c9360e0dc659d726` |

For the frozen-target Markdown file specifically:

| Content selection | SHA-256 |
|---|---|
| complete file, matching the original manifest | `abff4289459a53be63808210b37fd4025af2abf73bc82a2b2a41594953742a57` |
| fenced block including the `FZ02-TARGET-BEGIN` and `FZ02-TARGET-END` marker lines | `cdfd9ff814b19b3f70d2594299cd0047b2fdcd94a9fbcd90e8d42a46a1dee10d` |
| payload between the two marker lines | `a38b63c52e21625ce961ece55357865ee2148e625d416e6ff7885f4723ee48d8` |

Independent inspection with `ots info` confirms that each of the five
original `.ots` files binds the current whole-file SHA-256 of its paired
artifact. At the time of this erratum, each proof contains four pending
calendar attestations and no Bitcoin block attestation. The Bitcoin upgrade
therefore remains pending and must not be described as complete.

## Verification rule after this erratum

Verify the original artifacts against the whole-file hashes in
`registration_manifest_2026-07-26.json`. Verify the fenced target block
against the applicable hash above. Recompute the executable receipt at source
commit `091658ce585c107a260e7b980352be904d2419b2`. Retain the original artifacts
and `.ots` files byte-for-byte so their existing commitments remain valid.
