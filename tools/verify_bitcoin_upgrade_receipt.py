"""Verify the append-only FZ-02/10/11/12 Bitcoin proof upgrade receipt.

This is an offline check of the archived file hashes, OpenTimestamps paths,
and Bitcoin block headers. The two public explorer observations were made when
the receipt was written; this checker is not a Bitcoin full-node verifier.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CUSTODY = ROOT / "evidence" / "custody"
RECEIPT = CUSTODY / "bitcoin_upgrades_2026-10-04" / "upgrade_audit.json"
RECEIPT_SHA256 = "fb0ea4a2614a7ea33a5aa63f0c0b5419ec2ecc19c5f110d299985d90cb3ec8b0"
SOURCE_COMMIT = "0f7aa44259f53e4f5a1f4cd4904758039975755e"
UPGRADED_ROWS = {"FZ-02", "FZ-10", "FZ-11", "FZ-12"}
EXPLORERS = ["https://blockstream.info/api", "https://mempool.space/api"]
INFO_ROOT = re.compile(
    r"verify BitcoinBlockHeaderAttestation\((\d+)\)\s*\n"
    r"\s*# Bitcoin block merkle root ([0-9a-f]{64})"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(f"Bitcoin upgrade receipt: {message}")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def ots_client() -> str:
    candidates = [
        Path.home() / ".pyenv" / "versions" / "sherlock2" / "bin" / "ots",
        Path(shutil.which("ots") or ""),
    ]
    for candidate in candidates:
        if candidate.is_file():
            probe = subprocess.run(
                [str(candidate), "--version"], capture_output=True, text=True
            )
            if probe.returncode == 0:
                return str(candidate)
    raise ValueError("Bitcoin upgrade receipt: official OpenTimestamps client unavailable")


def historical_proof(relative: str) -> bytes:
    result = subprocess.run(
        ["git", "show", f"{SOURCE_COMMIT}:evidence/custody/{relative}"],
        cwd=ROOT,
        capture_output=True,
    )
    require(result.returncode == 0, f"historical proof unavailable: {relative}")
    return result.stdout


def verify_block(height: str, block: dict) -> None:
    require(
        set(block)
        == {
            "block_hash", "header_hex", "merkle_root", "block_time_utc",
            "best_chain_reported_by", "header_pow_valid",
        },
        f"block {height} keys differ",
    )
    require(block["best_chain_reported_by"] == EXPLORERS, f"block {height} explorer scope drifted")
    require(block["header_pow_valid"] is True, f"block {height} PoW claim absent")
    try:
        header = bytes.fromhex(block["header_hex"])
    except (TypeError, ValueError) as error:
        raise ValueError(f"Bitcoin upgrade receipt: invalid block {height} header") from error
    require(len(header) == 80, f"block {height} header must have 80 bytes")
    hashed = hashlib.sha256(hashlib.sha256(header).digest()).digest()
    require(hashed[::-1].hex() == block["block_hash"], f"block {height} hash mismatch")
    exponent = header[75]
    coefficient = int.from_bytes(header[72:75], "little")
    require(3 <= exponent <= 32 and 0 < coefficient < 0x800000, f"block {height} target invalid")
    target = coefficient << (8 * (exponent - 3))
    require(int.from_bytes(hashed, "little") <= target, f"block {height} PoW invalid")
    require(
        header[36:68][::-1].hex() == block["merkle_root"],
        f"block {height} Merkle root mismatch",
    )
    block_time = datetime.fromtimestamp(
        int.from_bytes(header[68:72], "little"), timezone.utc
    ).isoformat().replace("+00:00", "Z")
    require(block_time == block["block_time_utc"], f"block {height} time mismatch")


def verify(receipt_path: Path = RECEIPT) -> dict:
    payload = receipt_path.read_bytes()
    require(digest(payload) == RECEIPT_SHA256, "receipt bytes differ from fixed pin")
    audit = json.loads(payload)
    require(
        set(audit)
        == {
            "schema", "checked_at_utc", "original_custody_commit", "ots_client",
            "verification_method", "full_chain_independently_verified", "explorers",
            "packages", "blocks", "proofs",
        },
        "top-level keys differ",
    )
    require(audit["schema"] == "oph.legacy_bitcoin_proof_upgrade.v1", "schema drifted")
    require(audit["original_custody_commit"] == SOURCE_COMMIT, "source commit drifted")
    require(audit["ots_client"] == "opentimestamps-client v0.7.2", "client version drifted")
    require(audit["full_chain_independently_verified"] is False, "full-node claim is unsupported")
    require(audit["explorers"] == EXPLORERS, "explorer provenance drifted")
    require("explorer-assisted" in audit["verification_method"], "verification scope missing")
    require(set(audit["packages"]) == UPGRADED_ROWS, "upgraded package set differs")

    for height, block in audit["blocks"].items():
        require(height.isdecimal(), "non-numeric block height")
        verify_block(height, block)

    base = CUSTODY / "falsification" / "frozen_targets"
    # Only the four historical package directories, never a new freeze, can
    # appear in this dated upgrade. The original `.ots` files remain in place.
    expected = set()
    for row in UPGRADED_ROWS:
        directories = list(base.glob(f"{row.lower().replace('-', '')}_*"))
        require(len(directories) == 1, f"{row} original package missing or ambiguous")
        expected.update(
            f"falsification/frozen_targets/{proof.relative_to(base)}"
            for proof in directories[0].glob("*.ots")
        )
    actual = [item.get("original_proof") for item in audit["proofs"]]
    require(len(actual) == len(set(actual)) and set(actual) == expected, "proof set incomplete or changed")

    client = ots_client()
    seen_heights: set[str] = set()
    for item in audit["proofs"]:
        require(
            set(item)
            == {
                "original_proof", "upgraded_proof", "artifact_sha256",
                "original_proof_sha256", "upgraded_proof_sha256", "height", "merkle_root",
            },
            "proof entry keys differ",
        )
        original_rel = item["original_proof"]
        old = CUSTODY / original_rel
        artifact = CUSTODY / original_rel.removesuffix(".ots")
        relative_inside = Path(original_rel).relative_to("falsification/frozen_targets")
        new_rel = f"bitcoin_upgrades_2026-10-04/{relative_inside}"
        require(item["upgraded_proof"] == new_rel, f"upgraded path drifted: {original_rel}")
        new = CUSTODY / new_rel
        require(old.is_file() and new.is_file() and artifact.is_file(), f"missing proof or artifact: {original_rel}")
        original_bytes = old.read_bytes()
        upgraded_bytes = new.read_bytes()
        artifact_hash = digest(artifact.read_bytes())
        require(digest(original_bytes) == item["original_proof_sha256"], f"original proof hash mismatch: {original_rel}")
        require(original_bytes == historical_proof(original_rel), f"original proof was changed: {original_rel}")
        require(digest(upgraded_bytes) == item["upgraded_proof_sha256"], f"upgraded proof hash mismatch: {original_rel}")
        require(artifact_hash == item["artifact_sha256"], f"artifact hash mismatch: {original_rel}")
        info = subprocess.run([client, "info", str(new)], capture_output=True, text=True)
        require(info.returncode == 0, f"official client rejected upgraded proof: {original_rel}")
        require(
            f"File sha256 hash: {artifact_hash}" in info.stdout,
            f"proof does not bind artifact: {original_rel}",
        )
        roots = [(int(height), root) for height, root in INFO_ROOT.findall(info.stdout)]
        require(roots, f"Bitcoin attestation missing: {original_rel}")
        first_height = min(height for height, _ in roots)
        root = next(root for height, root in roots if height == first_height)
        require(
            (item["height"], item["merkle_root"]) == (first_height, root),
            f"attested height or Merkle root mismatch: {original_rel}",
        )
        block = audit["blocks"].get(str(first_height))
        require(block is not None and block["merkle_root"] == root, f"block witness missing: {original_rel}")
        seen_heights.add(str(first_height))
    require(seen_heights == set(audit["blocks"]), "unused or missing block witness")

    for row, summary in audit["packages"].items():
        package_proofs = [p for p in audit["proofs"] if f"/{row.lower().replace('-', '')}_" in p["original_proof"]]
        manifest = [p for p in package_proofs if Path(p["original_proof"]).name.startswith("registration_manifest_")]
        require(len(manifest) == 1, f"{row} registration manifest proof missing")
        require(
            summary == {"proof_count": len(package_proofs), "manifest_anchor_height": manifest[0]["height"]},
            f"{row} summary mismatch",
        )
    return audit


if __name__ == "__main__":
    result = verify()
    print(f"BITCOIN_UPGRADE_RECEIPT_VALID: {len(result['proofs'])} proofs, {len(result['blocks'])} block headers")
