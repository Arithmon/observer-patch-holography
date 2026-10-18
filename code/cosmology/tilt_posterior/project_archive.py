#!/usr/bin/env python3
"""Project the pinned ACT DR6 P-ACT-LB CAMB archive onto ordered weight/ns.

Only the two original ASCII tokens are retained from each numbered chain.
This reads tar members as streams; it never extracts files or executes the
archive's configuration. No burn-in, thinning or numerical rounding occurs.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import tarfile
import tempfile


ARCHIVE_SHA256 = "8f2801a69663ef1142bb05660e0118ef7b1e0bdf3d05afba434a0569023b98ab"
ARCHIVE_BYTES = 163494761
ARCHIVE_URL = "https://lambda.gsfc.nasa.gov/data/act/chains/lcdm/p-actbase-l-b_lcdm_camb.tar.gz"
MIRROR_URL = "https://portal.nersc.gov/project/act/dr6.02/chains/compressed/lcdm/p-actbase-l-b_lcdm_camb.tar.gz"
BURN_IN_SOURCE_URL = "https://portal.nersc.gov/project/act/dr6.02/chains/README_chains.html"
MEMBER_ROOT = "p-actbase-l-b_lcdm_camb/p-actbase-l-b_lcdm_camb"
EXPECTED_ROWS = {1: 202988, 2: 158360, 3: 236480, 4: 160522}


def sha256_file(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def parse_checkpoint(data: bytes) -> dict:
    """Read only the four primitive fields in the released checkpoint.

    A restricted parser avoids loading or evaluating YAML configuration.
    The entire source archive has already been verified against its digest.
    """
    lines = data.decode("ascii").splitlines()
    if lines[:2] != ["sampler:", "  mcmc:"]:
        raise ValueError("unexpected checkpoint layout")
    fields = {}
    for line in lines[2:]:
        if not line.startswith("    ") or ": " not in line:
            raise ValueError("unexpected checkpoint field")
        key, value = line[4:].split(": ", 1)
        if key in fields:
            raise ValueError("duplicate checkpoint field")
        if key == "converged" and value in {"true", "false"}:
            fields[key] = value == "true"
        elif key == "Rminus1_last":
            fields[key] = float(value)
            if not math.isfinite(fields[key]):
                raise ValueError("nonfinite checkpoint diagnostic")
        elif key in {"burn_in", "mpi_size"}:
            fields[key] = int(value)
        else:
            raise ValueError("unexpected checkpoint field or value")
    if set(fields) != {"converged", "Rminus1_last", "burn_in", "mpi_size"}:
        raise ValueError("incomplete checkpoint")
    return fields


def project_archive(archive: Path) -> dict:
    """Verify the release archive before reading its ordered chain members."""
    archive = Path(archive)
    archive_bytes = archive.stat().st_size
    archive_sha256 = sha256_file(archive)
    if archive_bytes != ARCHIVE_BYTES or archive_sha256 != ARCHIVE_SHA256:
        raise ValueError("source archive differs from the pinned ACT DR6 CAMB release")
    payload = {
        "schema": "oph.act_dr6.ns_projection.v1",
        "archive_sha256": archive_sha256,
        "archive_bytes": archive_bytes,
        "archive_url": ARCHIVE_URL,
        "mirror_url": MIRROR_URL,
        "burn_in_fraction": "0.5",
        "burn_in_policy": "Drop floor(0.5 * number_of_rows) rows from the beginning of each ordered chain before applying multiplicity weights.",
        "burn_in_source_url": BURN_IN_SOURCE_URL,
        "columns": ["weight", "ns"],
        "row_token_policy": "Original whitespace-delimited ASCII tokens preserved as strings; no thinning or burn-in applied to this projection.",
        "chains": [],
    }
    numbered_members = {f"{MEMBER_ROOT}.{i}.txt": i for i in EXPECTED_ROWS}
    seen = set()
    with tarfile.open(archive, "r|gz") as tar:
        for member in tar:
            if member.name not in numbered_members and member.name != f"{MEMBER_ROOT}.checkpoint":
                continue
            if member.name in seen or not member.isfile():
                raise ValueError("duplicate or nonregular source member")
            seen.add(member.name)
            stream = tar.extractfile(member)
            if stream is None:
                raise ValueError("source member cannot be read")
            with stream:
                if member.name in numbered_members:
                    source_hash = hashlib.sha256()
                    header = stream.readline()
                    source_hash.update(header)
                    columns = header.decode("ascii").lstrip("#").split()
                    if columns.count("weight") != 1 or columns.count("ns") != 1:
                        raise ValueError("ambiguous weight/ns columns")
                    weight_column, ns_column = columns.index("weight"), columns.index("ns")
                    if (weight_column, ns_column) != (0, 6):
                        raise ValueError("unexpected source column layout")
                    rows = []
                    for line in stream:
                        source_hash.update(line)
                        if not line.strip() or line.lstrip().startswith(b"#"):
                            continue
                        tokens = line.split()
                        if len(tokens) != len(columns):
                            raise ValueError("source row does not match its header")
                        rows.append([tokens[weight_column].decode("ascii"),
                                     tokens[ns_column].decode("ascii")])
                    if len(rows) != EXPECTED_ROWS[numbered_members[member.name]]:
                        raise ValueError("unexpected numbered-chain row count")
                    payload["chains"].append({
                        "member": member.name,
                        "source_bytes": member.size,
                        "source_sha256": source_hash.hexdigest(),
                        "source_column_indices_zero_based": [weight_column, ns_column],
                        "rows": rows,
                    })
                else:
                    data = stream.read()
                    payload["checkpoint"] = parse_checkpoint(data)
                    payload["checkpoint_source"] = {
                        "member": member.name,
                        "source_sha256": hashlib.sha256(data).hexdigest(),
                        "text": data.decode("ascii"),
                    }
    if seen != set(numbered_members) | {f"{MEMBER_ROOT}.checkpoint"}:
        raise ValueError("missing numbered chains or convergence checkpoint")
    payload["chains"].sort(key=lambda chain: numbered_members[chain["member"]])
    return payload


def write_projection(payload: dict, output: Path) -> None:
    """Write canonical JSON and deterministic gzip metadata atomically."""
    output = Path(output)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=output.parent, prefix=".act-projection-",
                                         delete=False) as raw:
            temporary = Path(raw.name)
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw,
                               mtime=0, compresslevel=9) as compressed:
                encoder = json.JSONEncoder(sort_keys=True, separators=(",", ":"),
                                           ensure_ascii=True)
                for chunk in encoder.iterencode(payload):
                    compressed.write(chunk.encode("ascii"))
                compressed.write(b"\n")
        os.replace(temporary, output)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="downloaded official .tar.gz archive")
    parser.add_argument("output", type=Path, help="destination projection .json.gz")
    args = parser.parse_args()
    if args.archive.resolve() == args.output.resolve():
        parser.error("projection output must differ from the source archive")
    try:
        payload = project_archive(args.archive)
        write_projection(payload, args.output)
    except (OSError, ValueError, tarfile.TarError) as exc:
        parser.exit(1, f"projection failed: {exc}\n")
    print(json.dumps({"output": str(args.output), "bytes": args.output.stat().st_size,
                      "sha256": sha256_file(args.output),
                      "rows": sum(len(chain["rows"]) for chain in payload["chains"])}))


if __name__ == "__main__":
    main()
