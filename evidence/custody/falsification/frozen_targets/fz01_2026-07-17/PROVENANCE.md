# Custody note (written 2026-07-22, after anchoring)

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
An earlier draft of the registry (anchored=false, no `frozen_utc`) exists at
`oph-meta/proof/epic_wins/fz01_registry/`; the anchored files here supersede it.

To verify independently: `sha256sum <artifact>` against the manifest, and
`ots verify <artifact>.ots` (requires a Bitcoin node) or `ots info` to inspect
the attestation path.
