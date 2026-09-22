# Custody note (written 2026-07-22, after anchoring; extended 2026-09-22)

The sixteen files in this directory (seven artifacts plus one detached `.ots`
proof each, and this note) are the FZ-01 forward-test registration of
2026-07-17T07:18:00Z.

Chain of custody:

1. Anchored in `FloatingPragma/observer-patch-holography` at commit
   `cb36e4d8bcaaff8ebfe301bcad1fdff1d23de94c` ("Register and anchor the
   2026-07-17 forward-test targets"), path `falsification/frozen_targets/`.
2. Removed from that repository's working tree the next day by commit
   `0c986ccbb81bd35e5c6857a645fe2b87382d4581` ("Repo hygiene (WIP)"). The
   blobs remain reachable at the anchoring commit.
3. Extracted byte-identical from the anchoring commit into this directory on
   2026-07-22. Every artifact re-verified against the sha256 values inside
   `registration_manifest_2026-07-17.json` (all match).
4. On 2026-07-22 the eight `.ots` proofs were upgraded from pending calendar
   attestations to complete Bitcoin attestations (`ots upgrade`). Upgrading a
   proof adds the calendar's Bitcoin path; it does not alter the attested
   digest or time.

Verdict weight begins at the 2026-07-17 anchoring per the registration policy.
An earlier unanchored draft of the registry (anchored=false, no `frozen_utc`)
exists outside this repository; the anchored files here are canonical.

To verify independently: `sha256sum <artifact>` against the manifest, and
`ots verify <artifact>.ots` (requires a Bitcoin node) or `ots info` to inspect
the attestation path.

## Bitcoin attestations

`ots info` on each of the eight upgraded proofs shows the same four calendar
paths, each ending in a Bitcoin block header attestation. The block heights
are 958383, 958384, 958388 and 958864. A third party with a Bitcoin node
verifies each proof against those block headers; the merkle roots printed by
`ots info` identify the blocks without the pinned client.

## Pinned source path in the registry

`fz01_freeze_registry_2026-07-17.json` records, under `pinned_source.file`,
the absolute path of the simulator tracker on the author's machine at
freeze time. That string is part of the attested bytes and is not rewritten.
It carries no verification weight: the pinned bytes are the vendored copy
`tracker_section8_pinned_text.md` in this directory, whose sha256
`54b216b52afe25895552a04b450c3cb33e303802e99618f24361b2bf2b9e0737` the
registry records as `sha256_section8_text`, and the registry's `repo_state`
field names the simulator commit `cce09b822d764e802850cdcb9652513df181e1ca`
whose working tree the section was cut from. The receipt-portability gate
`tools/check_receipt_portability.py` scans `evidence/` and lists exactly this
field as an attested allowance.

## Status headers in the frozen target files

The four `frozen_target_*.md` files open with "Status: DRAFT. Not
registered". Those headers are anchored bytes from the drafting stage and
stay as written. The registration manifest and `claims/frozen_prediction_register.json`
attest the files as registered, and `SCIENTIFIC_STATUS_ERRATUM_2026-07-29.md`
records the later scientific-status correction. The manifest, the register
and this note supersede the headers.

## Comparator provenance for rows R04 and R05

Row FZ01-R04 freezes the comparator `alpha_s(M_Z) = 0.1179 +- 0.0009` as
transcribed from the tracker. Row FZ01-R05 freezes `Lambda_QCD(3) = 338 +- 12
MeV` as a FLAG-class average without naming the FLAG edition. Under the
registry's scoring rule the frozen comparators are not retro-edited; each
row's decision trigger (the next PDG or FLAG release) names the edition used
at comparison time, and the comparison record states both the frozen and the
then-current comparator.
