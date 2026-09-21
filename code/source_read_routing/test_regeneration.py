"""Reject incomplete or forged regenerated streams without changing references."""
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import struct

import pytest

import regenerate
from pack import segment_rows
from verify import load

CONTROLS = Path(__file__).resolve().parent / 'controls'


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    receipt = tmp_path / 'q3_baseline.json'
    shutil.copyfile(CONTROLS / receipt.name, receipt)
    packet = load(receipt)
    raw = b''.join(block for segment in packet['tape']['segments']
                   for block in segment_rows(CONTROLS, segment['path']))
    for segment in packet['tape']['segments']:
        shutil.copyfile(CONTROLS / segment['path'], tmp_path / segment['path'])
    stream = tmp_path / 'events.bin'
    stream.write_bytes(raw)
    program = tmp_path / 'producer.py'
    program.write_text(
        'import pathlib, sys\n'
        f'sys.stdout.buffer.write(pathlib.Path({str(stream)!r}).read_bytes())\n'
        f'print({json.dumps(packet["costs"])!r}, file=sys.stderr)\n', encoding='utf-8')
    monkeypatch.setattr(regenerate, 'command', lambda *args: [sys.executable, str(program)])
    return receipt, packet, raw, stream, program


def collect(fixture):
    receipt, packet, _, _, _ = fixture
    return b''.join(regenerate.generated_rows(receipt, packet, Path('unused'),
                                              CONTROLS / 'q3.input', receipt.parent))


def test_complete_stream_is_identical_to_retained_control(fixture):
    assert collect(fixture) == fixture[2]


@pytest.mark.parametrize('mutation', ['empty', 'partial_row', 'missing_row', 'suffix',
                                    'value', 'nonzero_exit', 'missing_census', 'wrong_census'])
def test_generated_stream_fails_closed(fixture, mutation):
    _, packet, raw, stream, program = fixture
    if mutation == 'empty':
        stream.write_bytes(b'')
    elif mutation == 'partial_row':
        stream.write_bytes(raw[:-1])
    elif mutation == 'missing_row':
        stream.write_bytes(raw[:-64])
    elif mutation == 'suffix':
        stream.write_bytes(raw + bytes(64))
    elif mutation == 'value':
        changed = bytearray(raw)
        changed[63] ^= 1
        stream.write_bytes(changed)
    elif mutation == 'nonzero_exit':
        program.write_text(program.read_text() + 'sys.exit(7)\n', encoding='utf-8')
    elif mutation == 'missing_census':
        program.write_text('\n'.join(program.read_text().splitlines()[:2]), encoding='utf-8')
    else:
        changed = deepcopy(packet['costs'])
        changed['events'] -= 1
        program.write_text(program.read_text() + f'print({json.dumps(changed)!r}, file=sys.stderr)\n',
                           encoding='utf-8')
    with pytest.raises(ValueError):
        collect(fixture)


@pytest.mark.parametrize('mutation', ['missing_phase', 'duplicate_phase', 'path_escape',
                                    'manifest_bytes', 'phase_count', 'phase_checksum', 'producer_pin'])
def test_commitments_fail_closed_even_with_resealed_manifest(fixture, mutation):
    receipt, packet, _, _, _ = fixture
    segment = packet['tape']['segments'][0]
    path = receipt.parent / segment['path']
    if mutation == 'missing_phase':
        packet['tape']['segments'].pop()
    elif mutation == 'duplicate_phase':
        packet['tape']['segments'][1] = deepcopy(segment)
    elif mutation == 'path_escape':
        segment['path'] = '../' + segment['path']
    elif mutation == 'manifest_bytes':
        path.write_bytes(path.read_bytes() + b' ')
    elif mutation == 'producer_pin':
        packet['producer_sha256'] = '0' * 64
    else:
        manifest = load(path)
        if mutation == 'phase_count':
            manifest['events'] -= 1
        else:
            manifest['decoded_sha256'] = '0' * 64
        path.write_text(json.dumps(manifest), encoding='utf-8')
        segment['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError):
        collect(fixture)


@pytest.mark.skipif(not os.environ.get('ROUTING_VERIFIER'), reason='requires compiled native verifier')
def test_semantic_forgery_is_rejected_with_all_stream_hashes_resealed(fixture):
    receipt, packet, raw, stream, _ = fixture
    changed = bytearray(raw)
    mean_offset = next(offset for offset in range(0, len(raw), 64)
                       if struct.unpack_from('<Q', raw, offset)[0] == 5)
    changed[mean_offset + 56] ^= 1
    stream.write_bytes(changed)
    offset = 0
    for segment in packet['tape']['segments']:
        path = receipt.parent / segment['path']
        manifest = load(path)
        end = offset + manifest['events'] * 64
        manifest['decoded_sha256'] = hashlib.sha256(changed[offset:end]).hexdigest()
        path.write_text(json.dumps(manifest), encoding='utf-8')
        segment['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        offset = end
    packet['tape']['decoded_sha256'] = hashlib.sha256(changed).hexdigest()
    receipt.write_text(json.dumps(packet), encoding='utf-8')
    # The producer stream passes all its size and hash commitments. Rejection
    # must come from independently evaluating the changed scalar operation.
    assert collect(fixture) == changed
    with pytest.raises(ValueError, match='incorrect scalar operation'):
        regenerate.regenerate(receipt, Path('unused'), Path(os.environ['ROUTING_VERIFIER']),
                              receipt.parent / 'replay')
