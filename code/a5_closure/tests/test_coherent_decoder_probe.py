"""Hostile controls for the restricted clean-register field obstruction."""
from copy import deepcopy
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import coherent_decoder_probe as producer
import verify_coherent_decoder_independent as verifier


@pytest.fixture(scope='module')
def packet():
    return verifier.m.load(producer.REFERENCE)


def test_independent_spectral_reconstruction(packet):
    assert producer.build() == packet
    assert verifier.verify(packet)['actions_checked'] == 60
    assert packet['projector']['witness'] == ['1/20', '1/100']


@pytest.mark.parametrize('path,value', [
    (('marking','index'), 1),
    (('marking','actions_checked'), 1),
    (('projector','rank'), 2),
    (('projector','witness'), ['1/20','0']),
    (('projector','gram_scale'), ['360','0']),
    (('grammar','free_angles'), True),
    (('grammar','maximum_blocks'), 32),
    (('grammar','alphabet_size'), 1),
    (('target','arbitrary_residual_state'), True),
    (('target','postselection'), True),
    (('target','other_bits'), 'discarded without constraint'),
    (('target','receiver_codes'), [0,1]),
    (('matrix_units_required',), 3),
    (('candidate_decoders',), 1),
    (('reference_extension_required',), False),
    (('approximation_lower_bound_claimed',), True),
    (('apparatus_service_constructed',), True),
    (('physical_golden_embedding_selected',), True),
    (('general_attachment_verdict',), 'SOURCE_DERIVED_COHERENT_GV_DECODER'),
])
def test_scope_and_receipt_mutations(packet, path, value):
    changed = deepcopy(packet)
    parent = changed
    for key in path[:-1]:
        parent = parent[key]
    parent[path[-1]] = value
    with pytest.raises(ValueError, match='complete decoder obstruction receipt'):
        verifier.verify(changed)


@pytest.mark.parametrize('surface,before,after', [
    ('model/common_source/slow_geometry.py', '(minus, phi, zero)', '(one, phi, zero)'),
    ('model/observer/source_scalar/compiler.py', '"J", -math.pi,', '"J", -math.pi/2,'),
    ('model/native_interface/tomography.py', 'DIMENSION = 66', 'DIMENSION = 65'),
    ('model/compiled_native_interface/construction.py', "'Y', math.pi/8", "'Y', math.pi/7"),
])
def test_hostile_source_is_rejected_before_arithmetic(packet, tmp_path, surface, before, after):
    for path in verifier.PATHS:
        target = tmp_path/path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((verifier.ROOT/path).read_bytes())
    target = tmp_path/surface
    text = target.read_text()
    assert before in text
    target.write_text(text.replace(before, after, 1))
    with pytest.raises(ValueError, match='source custody'):
        verifier.verify(packet, root=tmp_path)


def test_marking_action_cannot_be_replaced_by_equal_trace(packet):
    changed = deepcopy(packet)
    action = changed['marking']['action']
    changed['marking']['action'] = [action.index(i) for i in range(12)]
    assert changed != packet
    with pytest.raises(ValueError, match='complete decoder obstruction receipt'):
        verifier.verify(changed)


def test_rational_projector_positive_boundary():
    # Identity is an admitted empty word and transports the rational target.
    # Nonrationality, not rank three or the mere demand for a qutrit, obstructs.
    from fractions import Fraction
    target = [[Fraction(int(i == j and i < 3)) for j in range(4)] for i in range(4)]
    assert sum(target[i][i] for i in range(4)) == 3
    assert [[sum(target[i][k]*target[k][j] for k in range(4))
             for j in range(4)] for i in range(4)] == target
