"""Controls for actual-source geometry, dependency records and cone scope."""
from itertools import product

import pytest
import sympy as sp

from . import causal_placement as c


def test_gram_is_the_existing_rank_three_response_not_an_added_mesh():
    g = c.gram()
    assert g.trace() == 12
    assert (g*g-4*g).applyfunc(sp.expand) == sp.zeros(12)
    assert set(g.diagonal()) == {1}


def test_complete_tape_parents_include_zero_current_events():
    rows = c.dependencies(((0, 1), (2, 3), (0, 1), (0, 1)))
    assert [r['parents'] for r in rows] == [(), (), (0,), (2,)]
    assert [r['height'] for r in rows] == [0, 0, 1, 2]
    assert rows[-1]['ancestors'] == (0, 2)
    assert c.interval_squared(rows[1], rows[3]) == 0
    assert 1 not in rows[3]['ancestors']


def test_nontrivial_two_endpoint_parent_merge_and_forward_cones():
    seams, _ = c.passage.recovered_carrier()
    first, last = seams[0], seams[-1]
    tape = (*seams, first, last, *reversed(seams))
    rows = c.dependencies(tape)
    assert any(len(r['parents']) == 2 for r in rows)
    for row in rows:
        for ancestor in row['ancestors']:
            assert row['height'] > rows[ancestor]['height']
            assert c.interval_squared(rows[ancestor], row).is_nonnegative


def test_all_native_converse_branches_are_retained_and_counts_are_exact():
    result = c.run()
    control = result['converse_control']
    assert [r['coins'] for r in control['all_coin_branches']] == list(map(list, product((0, 1), repeat=4)))
    assert any(record[1] == 0 for r in control['all_coin_branches'] for record in r['records'])
    assert result['finite_current_readout'] == {
        'total_two_states': 78, 'distinct_rank_three_images': 73,
        'binary_states': 66, 'distinct_binary_images': 61,
        'native_branch_telescoping_checks': 4680,
        'unit_centroid_normalization': 'raw_load_readout / (2*sqrt((5+sqrt(5))/2))',
        'normalized_load_centroid_in_convex_hull_of_unit_ports': True,
        'waiting_on_fixed_total_carrier_realizes_full_signed_control_completion': False}
    pair = result['forward_cone']['sharp_pair']
    assert c.distance_squared(*pair) == c.SPEED_SQUARED
    assert not result['acceptance']['cone_implies_ancestry']


@pytest.mark.parametrize('seam', [(True, 1), (1, 0), (0, 0), (0, 2), (0, 12), '01'])
def test_reject_foreign_seams(seam):
    with pytest.raises(ValueError):
        c.dependencies([seam])


@pytest.mark.parametrize('height', [True, -1, 0.5, 1.0, '1'])
def test_reject_inexact_or_foreign_height(height):
    with pytest.raises(ValueError):
        c.interval_squared({'seam': (0, 1), 'height': height},
                           {'seam': (0, 1), 'height': 2})
