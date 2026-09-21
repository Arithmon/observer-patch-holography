"""Outer archive/native-execution adapter; never called by blind recovery."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
ROUTING = ROOT / 'code/source_read_routing'
CONTROL = ROUTING / 'controls'
VARIANTS = ('baseline', 'source', 'branch', 'scratch')
REMOVE = ((120, 5, 123, 4), (132, 9, 135, 11))
ADD = ((120, 5, 135, 11), (123, 4, 132, 9))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def retained_bytes(variant='baseline'):
    if variant not in VARIANTS:
        raise ValueError('unknown retained variant')
    def load(name, path):
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    old_tape = sys.modules.get('tape')
    try:
        sys.modules['tape'] = load('_federation_archive_tape', ROUTING / 'tape.py')
        segment_rows = load('_federation_archive_pack', ROUTING / 'pack.py').segment_rows
        receipt = json.loads((CONTROL / f'q3_{variant}.json').read_text(encoding='utf-8'))
        pieces = []
        for segment in receipt['tape']['segments']:
            name = segment['path']
            if Path(name).name != name or sha((CONTROL / name).read_bytes()) != segment['sha256']:
                raise ValueError('segment manifest custody')
            pieces.extend(segment_rows(CONTROL, name))
        data = b''.join(pieces)
        tape = receipt['tape']
        if len(data) != tape['decoded_bytes'] or len(data) != 64 * tape['events'] or sha(data) != tape['decoded_sha256']:
            raise ValueError('whole tape custody')
        return data
    finally:
        if old_tape is None:
            sys.modules.pop('tape', None)
        else:
            sys.modules['tape'] = old_tape


def retained_rows(variant='baseline'):
    return list(struct.iter_unpack('<8Q', retained_bytes(variant)))


def read_input():
    data = (CONTROL / 'q3.input').read_bytes()
    expected = json.loads((CONTROL / 'q3_baseline.json').read_text(encoding='utf-8'))['inputs']['input_sha256']
    if sha(data) != expected:
        raise ValueError('original input custody')
    header = struct.unpack_from('<8Q', data)
    _, n, _, _, _, count, _, _ = header
    hosts = list(struct.unpack_from(f'<{n}I', data, 64))
    edges = list(struct.iter_unpack('<4I', data[64 + 4*n:64 + 4*n + 16*count]))
    return header, hosts, edges


def alternative_input():
    header, _, edges = read_input()
    changed = [tuple(row) for row in edges]
    for old, new in zip(REMOVE, ADD):
        if changed.count(old) != 1:
            raise ValueError('witness seam absent or duplicated')
        changed[changed.index(old)] = new
    data = bytearray((CONTROL / 'q3.input').read_bytes())
    start = 64 + 4*header[1]
    data[start:start + 16*len(changed)] = b''.join(struct.pack('<4I', *row) for row in changed)
    return bytes(data)


def native_replays():
    """Actually compile and run both completions; retain only small digests."""
    compiler = shutil.which('c++') or shutil.which('g++') or shutil.which('clang++')
    if compiler is None:
        raise RuntimeError('a C++17 compiler is required for native reproduction')
    result = {}
    with tempfile.TemporaryDirectory(prefix='oph-federation-recovery-') as name:
        work = Path(name)
        binary, alternate = work / 'produce', work / 'alternative.input'
        subprocess.run([compiler, '-std=c++17', '-O2', str(ROUTING / 'produce.cpp'), '-o', str(binary)], check=True, capture_output=True)
        alternate.write_bytes(alternative_input())
        for variant in VARIANTS:
            original = subprocess.run([str(binary), str(CONTROL / 'q3.input'), variant], check=True, capture_output=True).stdout
            other = subprocess.run([str(binary), str(alternate), variant], check=True, capture_output=True).stdout
            retained = retained_bytes(variant)
            if original != other or original != retained:
                raise ValueError(f'{variant}: native completion/retained mismatch')
            result[variant] = {'events': len(original)//64, 'decoded_sha256': sha(original),
                               'original_equals_alternative': True, 'original_matches_retained': True}
    return result


def source_hashes():
    paths = [p for p in Path(__file__).parent.glob('*.py')]
    paths += [ROUTING / p for p in ('produce.cpp', 'pack.py', 'tape.py')]
    paths += [CONTROL / 'q3.input', ROOT / 'code/source_routing/support_w12_l3.json']
    paths += list(CONTROL.glob('q3_*.json')) + list(CONTROL.glob('q3_*.tape'))
    return {p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in sorted(paths)}


if __name__ == '__main__':
    print(json.dumps(native_replays(), indent=2, sort_keys=True))
