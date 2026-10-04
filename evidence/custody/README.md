# Frozen-prediction custody

`falsification/frozen_targets/` holds the custody packages named in
`claims/frozen_prediction_register.json`: each frozen
target statement, its registration manifest, decision rules, errata, the
pinned Lean and receipt copies, and a detached OpenTimestamps proof (`.ots`)
for every stamped file. Historical packages were vendored byte for byte;
later registrations are appended in their own dated directories.
`.ots.bak` files retain proofs before their Bitcoin upgrade.

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
available, upgrade the proof with the official client:

```sh
ots upgrade registration_manifest_2026-10-04.json.ots
ots verify registration_manifest_2026-10-04.json.ots
```

Run these beside the unchanged manifest and retain the earlier proof.
Upgrading changes the proof, not the target, protocol or manifest. The
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
