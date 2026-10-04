# Frozen-prediction custody

`falsification/frozen_targets/` holds the custody packages named in
`claims/frozen_prediction_register.json`: each frozen
target statement, its registration manifest, decision rules, errata, the
pinned Lean and receipt copies, and a detached OpenTimestamps proof (`.ots`)
for every stamped file. Historical packages were vendored byte for byte;
later registrations are appended in their own dated directories.
`.ots.bak` files retain proofs before their Bitcoin upgrade.

The append-only `bitcoin_upgrades_2026-10-04/` directory holds upgraded
copies of all 23 archived FZ-02, FZ-10, FZ-11 and FZ-12 proofs. Its
`upgrade_audit.json` binds each copy to the unchanged original proof and
artifact and records six Bitcoin block headers. Run
`python3 tools/verify_bitcoin_upgrade_receipt.py` to check those bindings,
OpenTimestamps paths, header hashes, Merkle roots and proof of work offline.
Two public explorers reported these blocks on the best chain when checked;
this is explorer-assisted verification, not a Bitcoin full-node check.
With a configured Bitcoin Core node, independently verify an individual
copy using `ots verify -f <original-artifact> <upgraded-proof.ots>`.

`python3 tools/build_fz_registry.py --check` recomputes every manifest and
artifact hash from these copies, parses each proof with the official client
(`ots info`, local only) and checks its digest and attestation class. To check
a proof against the Bitcoin chain, run `ots verify <file>.ots` next to the
file. The timestamps anchor the bytes independently of any Git history; the
original custody commits are recorded in the register as provenance.

The FZ-13 and FZ-15 packages dated 2026-10-04 contain a human-readable
`target.md`, a complete normative `frozen_prediction.json`, and a manifest
binding both files. Their fixed-capacity and scalar-tilt statements are
conditional branch commitments. The stamped protocol fixes the physical
premises, data selection, posterior decision rule, numerical checks and
interpretation of every outcome. It does not supply an empirical score.

An OpenTimestamps calendar receives a randomized hash commitment, rather
than the document. Public calendars batch commitments into Bitcoin without
requiring a wallet or a payment from the registrant. The initial `.ots`
receipt can contain pending calendar attestations. Once inclusion is
available, upgrade a copy with the official client. The registry currently
pins the original calendar-pending `.ots` files, so an in-place upgrade in
this directory would fail its custody check:

```sh
proof_check_dir=$(mktemp -d)
cp registration_manifest_2026-10-04.json "$proof_check_dir/"
cp registration_manifest_2026-10-04.json.ots "$proof_check_dir/"
ots upgrade "$proof_check_dir/registration_manifest_2026-10-04.json.ots"
ots verify "$proof_check_dir/registration_manifest_2026-10-04.json.ots"
```

Run these from the relevant package directory. If the calendar reports
"pending confirmation," retain the original proof and try upgrading again
later; there is no need to stamp the files a second time. Once independently
verified, publish the upgraded proof and its verification record as an
append-only custody update. The original target, protocol, manifest and
calendar proof stay fixed. Upgrading changes only the copied proof. The
Python client's independent chain verification requires a Bitcoin Core
node (a pruned node is sufficient); verification using public block
explorers must be identified as explorer-assisted. `ots info` alone parses
the proof and does not verify its block against the chain.

The timestamp establishes that the exact bytes existed by the attesting
block. It does not prove authorship, scientific validity or absence of
earlier data exposure. A document's author-declared `frozen_utc` and its
later Bitcoin block timestamp are distinct facts. The new protocols require
publication and verified Bitcoin anchoring before a comparison is eligible.
See the [OpenTimestamps explanation](https://opentimestamps.org/) and
[official client documentation](https://github.com/opentimestamps/opentimestamps-client).
