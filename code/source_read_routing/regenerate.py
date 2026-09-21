"""Regenerate complete event streams and independently authenticate every row.

The reference receipts and segment digests are read-only. Native production
and native verification execute separately, with exact input reconstruction,
per-segment and whole-stream checksums, resource accounting and logical replay.
Production streams use bounded pipes and need no retained multi-gigabyte file.
"""
from __future__ import annotations

import argparse
from contextlib import closing, redirect_stdout
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from prepare import prepare
from run import command, native_path
from tape import hash_file
from verify import load, require, verify

HERE = Path(__file__).resolve().parent


def phase_commitments(receipt_path, packet):
    """Authenticate the unchanged retained phase commitments before execution."""
    q, variant = packet['inputs']['q'], packet['variant']
    n, c, rounds = (packet['inputs'][key] for key in ('sites', 'carriers', 'rounds'))
    require(all(type(value) is int and value > 0 for value in (n, c, rounds)), 'phase dimensions')
    prelude = 13*c + 8*n
    per_layer, remainder = divmod(packet['costs']['events'] - prelude, rounds)
    require(remainder == 0, 'nonintegral phase census')
    segments = packet['tape']['segments']
    require(type(segments) is list and len(segments) == rounds + 1, 'complete phase commitments')
    result = []
    for phase, segment in enumerate(segments):
        name = f'q{q}_{variant}_phase{phase}.segment.json'
        require(segment['path'] == name, 'canonical phase commitment path')
        path = Path(receipt_path).parent / name
        require(hash_file(path) == segment['sha256'], 'phase manifest custody')
        manifest = load(path)
        require(manifest['schema'] == 'oph.source_read_routing.segment.v1', 'phase schema')
        expected = prelude if phase == 0 else per_layer
        require(type(manifest['events']) is int and manifest['events'] == expected, 'phase event census')
        result.append((expected * 64, manifest['decoded_sha256']))
    return result


def generated_rows(receipt_path, packet, producer, input_path, work):
    phases = phase_commitments(receipt_path, packet)
    require(packet['producer_sha256'] == hash_file(HERE / 'produce.cpp'), 'producer source pin')
    diagnostics = work / f"q{packet['inputs']['q']}_{packet['variant']}.producer.stderr"
    with diagnostics.open('wb') as errors:
        process = subprocess.Popen(command(producer, native_path(input_path), packet['variant']),
                                   stdout=subprocess.PIPE, stderr=errors)
        try:
            for expected_size, expected_hash in phases:
                digest, remaining = hashlib.sha256(), expected_size
                while remaining:
                    data = process.stdout.read(min(4 * 1024 * 1024, remaining))
                    require(bool(data) and len(data) % 64 == 0, 'truncated generated phase')
                    digest.update(data)
                    remaining -= len(data)
                    yield data
                require(digest.hexdigest() == expected_hash, 'generated phase checksum')
            require(not process.stdout.read(1), 'generated event suffix')
            require(process.wait() == 0, 'producer failed: ' + diagnostics.read_text(encoding='utf-8'))
            lines = diagnostics.read_text(encoding='utf-8').splitlines()
            require(bool(lines), 'missing producer census')
            reported = json.loads(lines[-1])
            for key, value in packet['costs'].items():
                require(type(reported.get(key)) is int and reported[key] == value, 'producer census: ' + key)
        finally:
            if process.poll() is None:
                process.kill()
            process.wait()
            process.stdout.close()


def regenerate(receipt_path, producer, verifier, work):
    receipt_path, producer, verifier, work = map(Path, (receipt_path, producer, verifier, work))
    packet = load(receipt_path)
    work.mkdir(parents=True, exist_ok=True)
    input_path = work / f"q{packet['inputs']['q']}.input"
    with redirect_stdout(sys.stderr):
        info = prepare(packet['inputs']['q'], input_path)
    require(info == packet['inputs'], 'regenerated input census differs')
    # The independent verifier reconstructs the input and read menu itself,
    # checks every event and consumed writer, and authenticates the complete
    # decoded stream. A valid prefix or a matching summary cannot pass.
    with closing(generated_rows(receipt_path, packet, producer, input_path, work)) as rows:
        return verify(receipt_path, verifier, work, input_path=input_path, raw_rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', type=Path)
    parser.add_argument('--producer', type=Path, required=True)
    parser.add_argument('--verifier', type=Path, required=True)
    parser.add_argument('--work', type=Path)
    args = parser.parse_args()
    try:
        if args.work:
            regenerate(args.receipt, args.producer, args.verifier, args.work)
        else:
            with tempfile.TemporaryDirectory(prefix='oph-routing-regenerate-') as directory:
                regenerate(args.receipt, args.producer, args.verifier, Path(directory))
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.exit(1, f'regeneration failed: {error}\n')
