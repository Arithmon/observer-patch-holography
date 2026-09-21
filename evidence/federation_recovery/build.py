"""Build the bounded complete-gluing nonidentification certificate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import struct

from custody import alternative_input, native_replays, read_input, retained_bytes, retained_rows, sha, source_hashes
from recover import instantiate, recover
from topology import analyze_completions, make_alternative

ASSUMPTIONS = [
    'authenticated_complete_retained_event_rows',
    'declared_unsigned64_rows_signed_half_unit_values',
    'declared_M1_opcode_and_register_ownership_semantics',
    'candidate_class_connected_finite_twelve_port_gluings',
    'level_three_provenance_supplied_outside_recovery',
]
NONCLAIMS = [
    'no_canonical_geometric_alternative_asserted',
    'no_source_selection_of_feedback_law',
    'no_unseen_population_or_routing_law_recovery',
    'no_complete_gluing_fixed_point',
    'no_physical_clock_or_quantum_instrument',
    'no_universe_level_existence_or_uniqueness',
]


def build():
    original_bytes = retained_bytes()
    recovered = recover(retained_rows())
    current = recovered
    rounds = []
    for iteration in (1, 2):
        rows = instantiate(current)
        raw = b''.join(struct.pack('<8Q', *row) for row in rows)
        again = recover(rows)
        if raw != original_bytes or again != recovered:
            raise ValueError('finite observed program roundtrip differs')
        rounds.append({'iteration': iteration, 'decoded_sha256': sha(raw), 'recovered_specification_equal': True})
        current = again
    _, hosts, edges = read_input()
    completions = analyze_completions(edges, make_alternative(edges), hosts)
    if not completions['witness_valid']:
        raise ValueError('topology witness failed')
    return {
        'schema': 'oph.federation_recovery.v1',
        'verdict': 'complete_gluing_not_identified',
        'first_failed_invariant': 'complete_gluing_triangle_count',
        'fixed_point_attained': False,
        'assumptions': ASSUMPTIONS, 'nonclaims': NONCLAIMS,
        'recovered': recovered, 'completions': completions,
        'native_replays': native_replays(),
        'finite_program_roundtrip': {'scope': 'observed_instruction_program_only', 'iterations': rounds},
        'alternative_input_sha256': sha(alternative_input()),
        'source_hashes': source_hashes(),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    path = Path(__file__).with_name('receipt.json')
    result = build()
    text = json.dumps(result, sort_keys=True, indent=2) + '\n'
    if args.check:
        if path.read_text(encoding='utf-8') != text:
            raise SystemExit('federation recovery receipt differs')
        print('Federation recovery receipt reproduces exactly')
    else:
        path.write_text(text, encoding='utf-8', newline='\n')
        print(f'Wrote {path.name}: {len(text.encode())} bytes')
