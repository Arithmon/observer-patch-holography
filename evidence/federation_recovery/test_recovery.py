"""Adversarial checks of the observed-word and bounded ambiguity certificate."""
from __future__ import annotations

import builtins
import copy
import importlib.util
import io
import json
from pathlib import Path
import subprocess

import pytest


HERE = Path(__file__).resolve().parent


def _load(name):
    spec = importlib.util.spec_from_file_location('_federation_test_' + name, HERE / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


verify = _load('verify')
recover = _load('recover')
custody = _load('custody')
topology = _load('topology')


@pytest.fixture(scope='module')
def baseline():
    rows = custody.retained_rows()
    recovered = recover.recover(rows)
    header, hosts, edges = custody.read_input()
    return {'rows': rows, 'recovered': recovered, 'header': header, 'hosts': hosts,
            'edges': edges, 'input_bytes': (custody.CONTROL / 'q3.input').read_bytes()}


@pytest.fixture(scope='module')
def certificate(baseline):
    # Source pins may change while this development suite is edited. The
    # public --check/verify commands additionally require the retained receipt
    # to carry the current pins. All scientific fields come from that receipt.
    packet = verify.load_json(HERE / 'receipt.json')
    packet['source_hashes'] = custody.source_hashes()
    arguments = {key: value for key, value in baseline.items() if key != 'recovered'}
    arguments.update(native_runs=custody.native_replays(), source_pins=custody.source_hashes())
    return packet, arguments


def changed_row(rows, event, field, value):
    result = list(rows)
    row = list(result[event])
    row[field] = value
    result[event] = tuple(row)
    return result


def test_exact_independent_reconstruction(baseline):
    recovered = verify.check_rows(baseline['rows'])
    assert verify.canonical(recovered) == verify.canonical(baseline['recovered'])
    assert recovered['census']['events'] == 22818
    assert recovered['invariants'] == {
        'owners': 1280, 'ports': 15360, 'population_sites': 27, 'rounds': 2,
        'observed_glued_seams': 86, 'logical_nodes': 81, 'logical_reads': 622,
        'logical_strict_relations': 1331, 'logical_height': 3}
    assert recovered['poset']['reflexive_relations'] == 1412
    assert recovered['poset']['layer_sizes'] == [27, 27, 27]


def test_independent_instruction_replay_and_second_recovery(baseline):
    actual = verify.execute_program(baseline['recovered']['program'])
    assert actual == baseline['rows'] == recover.instantiate(baseline['recovered'])
    result = verify.check_roundtrip(baseline['recovered'], baseline['rows'])
    assert result['scope'] == 'observed_instruction_program_only'
    assert [x['iteration'] for x in result['iterations']] == [1, 2]
    assert all(x['recovered_specification_equal'] is True for x in result['iterations'])


def test_signed_values_have_same_meaning(baseline):
    rows = [tuple(row[:7]) + (row[7] - (1 << 64) if row[7] >= 1 << 63 else row[7],)
            for row in baseline['rows']]
    assert verify.canonical(verify.check_rows(rows)) == verify.canonical(baseline['recovered'])
    assert recover.recover(rows) == baseline['recovered']


def test_hidden_input_firewall(baseline, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('blind recovery attempted external I/O')
    for obj, method in ((builtins, 'open'), (io, 'open'), (Path, 'read_bytes'),
                        (Path, 'read_text'), (subprocess, 'run'), (subprocess, 'Popen')):
        monkeypatch.setattr(obj, method, forbidden)
    assert recover.recover(baseline['rows']) == baseline['recovered']
    assert verify.check_rows(baseline['rows']) == baseline['recovered']
    assert verify.execute_program(baseline['recovered']['program']) == baseline['rows']


@pytest.mark.parametrize('variant', ['source', 'branch', 'scratch'])
def test_foreign_cohort_cannot_be_relabeled_baseline(variant):
    rows = custody.retained_rows(variant)
    for decoder in (recover.recover, verify.check_rows):
        with pytest.raises(ValueError):
            decoder(rows)


@pytest.mark.parametrize('decoder', [recover.recover, verify.check_rows])
def test_nonzero_mean_uses_arithmetic_not_opcode_label(baseline, decoder):
    eid = next(i for i, row in enumerate(baseline['rows']) if row[0] == 5 and row[7])
    rows = changed_row(baseline['rows'], eid, 7, baseline['rows'][eid][7] + 2)
    with pytest.raises(ValueError, match='mean|result|arithmetic'):
        decoder(rows)


@pytest.mark.parametrize('decoder', [recover.recover, verify.check_rows])
def test_stale_consumed_writer_is_rejected(baseline, decoder):
    eid = next(i for i, row in enumerate(baseline['rows']) if row[0] == 3)
    rows = changed_row(baseline['rows'], eid, 5, baseline['rows'][eid][5] - 1)
    with pytest.raises(ValueError, match='writer'):
        decoder(rows)


@pytest.mark.parametrize('field', [0, 3, 4, 5, 7])
def test_boolean_event_fields_are_not_integer_evidence(baseline, field):
    for decoder in (recover.recover, verify.check_rows):
        with pytest.raises(ValueError):
            decoder(changed_row(baseline['rows'], 0, field, False))


def rewired_program(baseline, both_ends):
    recovered = copy.deepcopy(baseline['recovered'])
    program = recovered['program']
    start = next(i for i, row in enumerate(program) if row[0] == 3)
    assert start == 16883 and program[start][3] == 1
    assert program[start+2][2] == 14
    for offset, field in ((1, 3), (2, 2), (3, 1), (5, 3)):
        program[start+offset][field] = 15
    if both_ends:
        for offset, field in ((0, 3), (2, 1), (2, 3), (4, 3)):
            program[start+offset][field] = 2
    # These rows have fully recomputed, correct arithmetic and writer IDs.
    # The defect is an inconsistent fixed port gluing, not corrupt values.
    return recover.instantiate(recovered)


@pytest.mark.parametrize('decoder', [recover.recover, verify.check_rows])
def test_arithmetically_resealed_second_port_partner(baseline, decoder):
    rows = rewired_program(baseline, False)
    with pytest.raises(ValueError, match='partner|gluing'):
        decoder(rows)


@pytest.mark.parametrize('decoder', [recover.recover, verify.check_rows])
def test_arithmetically_resealed_second_carrier_pair_seam(baseline, decoder):
    rows = rewired_program(baseline, True)
    with pytest.raises(ValueError, match='carrier|seam|partner'):
        decoder(rows)


def test_computed_instruction_cannot_supply_its_answer(baseline):
    recovered = copy.deepcopy(baseline['recovered'])
    eid = next(i for i, row in enumerate(recovered['program']) if row[0] == 5)
    recovered['program'][eid][5] = baseline['rows'][eid][7]
    for interpreter in (lambda r: verify.execute_program(r['program']), recover.instantiate):
        with pytest.raises(ValueError, match='literal'):
            interpreter(recovered)


def test_independent_topology_witness(baseline):
    h, hosts, edges = (baseline[k] for k in ('header', 'hosts', 'edges'))
    actual = verify.check_completions(h, hosts, edges)
    expected = topology.analyze_completions(edges, topology.make_alternative(edges), hosts)
    assert actual == expected
    assert (actual['original']['triangles'], actual['alternate']['triangles']) == (14000, 13990)
    assert actual['original_BFS_sha256'] == actual['alternate_BFS_sha256']
    assert actual['original']['canonical_edges_sha256'] != actual['alternate']['canonical_edges_sha256']


def test_missing_or_altered_topology_switch_is_rejected(baseline):
    altered = verify.alternative_edges(baseline['edges'])
    altered[0][1] = (altered[0][1] + 1) % 12
    with pytest.raises(ValueError, match='switch'):
        verify.check_completions(baseline['header'], baseline['hosts'], baseline['edges'], altered)


def test_duplicate_complete_graph_port_is_rejected(baseline):
    edges = [list(row) for row in baseline['edges']]
    first = edges[0]
    second = next(row for row in edges[1:] if row[0] == first[0])
    second[1] = first[1]
    with pytest.raises(ValueError, match='port'):
        verify.check_completions(baseline['header'], baseline['hosts'], edges)


def test_complete_certificate_has_only_bounded_verdict(certificate):
    packet, arguments = certificate
    result = verify.verify_packet(packet, **arguments)
    assert result == {'verified': True, 'verdict': 'complete_gluing_not_identified',
                      'fixed_point_attained': False, 'events': 22818,
                      'triangle_counts': [14000, 13990], 'physical_identification': False}


@pytest.mark.parametrize('path,value', [
    (('fixed_point_attained',), True),
    (('fixed_point_attained',), 0),
    (('recovered', 'scope', 'full_port_gluing_recovered'), True),
    (('recovered', 'scope', 'physical_identification'), True),
    (('recovered', 'scope', 'event_only'), 1),
    (('completions', 'alternate', 'triangles'), 14000),
    (('completions', 'comparison', 'all_host_BFS_parent_maps_equal'), 1),
    (('native_replays', 'baseline', 'original_equals_alternative'), 1),
    (('finite_program_roundtrip', 'scope'), 'complete_federation_fixed_point'),
    (('finite_program_roundtrip', 'iterations', 0, 'decoded_sha256'), '0'*64),
    (('finite_program_roundtrip', 'iterations', 1, 'recovered_specification_equal'), 1),
    (('alternative_input_sha256',), '0'*64),
])
def test_receipt_claim_mutants_fail_closed(certificate, path, value):
    original, arguments = certificate
    packet = copy.deepcopy(original)
    target = packet
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        verify.verify_packet(packet, **arguments)


def test_extra_completion_field_is_rejected(certificate):
    original, arguments = certificate
    packet = dict(original, complete_federation_recovered=True)
    with pytest.raises(ValueError, match='inventory'):
        verify.verify_packet(packet, **arguments)


def test_missing_roundtrip_field_is_rejected(certificate):
    original, arguments = certificate
    packet = dict(original)
    del packet['finite_program_roundtrip']
    with pytest.raises(ValueError, match='inventory'):
        verify.verify_packet(packet, **arguments)


def test_changed_used_seam_in_receipt_is_rejected(certificate):
    original, arguments = certificate
    packet = copy.deepcopy(original)
    packet['recovered']['observed_gluing'][0][1] += 1
    with pytest.raises(ValueError, match='event recovery'):
        verify.verify_packet(packet, **arguments)


def test_altered_input_tail_cannot_hide_under_same_graph(certificate):
    original, arguments = certificate
    altered = dict(arguments)
    raw = arguments['input_bytes']
    altered['input_bytes'] = raw[:-1] + bytes([raw[-1] ^ 1])
    with pytest.raises(ValueError, match='input bytes custody'):
        verify.verify_packet(original, **altered)


def test_source_receipt_pin_cannot_authenticate_itself(certificate):
    original, arguments = certificate
    packet = copy.deepcopy(original)
    packet['source_hashes']['evidence/federation_recovery/verify.py'] = '0'*64
    with pytest.raises(ValueError, match='source custody'):
        verify.verify_packet(packet, **arguments)


def test_changed_preparation_cannot_borrow_original_native_identity(certificate, baseline):
    original, arguments = certificate
    altered = copy.deepcopy(baseline['recovered'])
    program = altered['program']
    site = altered['population'][0]
    address = site['address_registers'][0]
    new_coefficient = site['address_coefficients'][0] + 1
    program[address][5] = 2 * new_coefficient
    ports = altered['carrier_ports'][site['owner']]
    program[ports[0]][5] = 2 * max(new_coefficient, 0)
    program[ports[3]][5] = 2 * max(-new_coefficient, 0)
    rows = recover.instantiate(altered)
    # Preparation/address split remains valid and the observed menu/order
    # is unchanged. Borrowing the original native digest must still fail.
    recovered = verify.check_rows(rows)
    assert recovered['invariants'] == baseline['recovered']['invariants']
    packet = copy.deepcopy(original)
    packet['recovered'] = recovered
    packet['finite_program_roundtrip'] = verify.check_roundtrip(recovered, rows)
    arguments = dict(arguments, rows=rows)
    with pytest.raises(ValueError, match='authenticated native baseline'):
        verify.verify_packet(packet, **arguments)


@pytest.mark.parametrize('text', ['{"a":1,"a":2}', '{"a":1.0}', '{"a":NaN}'])
def test_json_ambiguities_are_rejected(tmp_path, text):
    path = tmp_path / 'bad.json'
    path.write_text(text, encoding='utf-8')
    with pytest.raises(ValueError):
        verify.load_json(path)


def test_json_utf8_is_explicit(tmp_path):
    path = tmp_path / 'unicode.json'
    path.write_text(json.dumps({'description': 'bounded observer λ'}, ensure_ascii=False), encoding='utf-8')
    assert verify.load_json(path) == {'description': 'bounded observer λ'}
