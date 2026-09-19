# Frozen-prediction custody

`falsification/frozen_targets/` holds byte-exact copies of the custody
directories named in `claims/frozen_prediction_register.json`: each frozen
target statement, its registration manifest, decision rules, errata, the
pinned Lean and receipt copies, and a detached OpenTimestamps proof (`.ots`)
for every stamped file. `.ots.bak` files are the proofs as first anchored,
before their Bitcoin upgrade.

`python3 tools/build_fz_registry.py --check` recomputes every manifest and
artifact hash from these copies, parses each proof with the official client
(`ots info`, local only) and checks its digest and attestation class. To check
a proof against the Bitcoin chain, run `ots verify <file>.ots` next to the
file. The timestamps anchor the bytes independently of any Git history; the
original custody commits are recorded in the register as provenance.
