"""Actual source-to-theorem mapping and the limits of finite matter selection."""
from copy import deepcopy

import pytest

from . import native_matter_theorems as m


@pytest.fixture(scope='module')
def report():
    return m.run()


def test_all_source_components_and_complements_are_retained(report):
    occupations = [i for row in report['components'] for i in row['occupation_masks']]
    assert sorted(occupations+report['excluded_menu_singlet_masks_retained']) == list(range(32))
    assert report['source_current_component_invariance_checks'] == 120
    assert not report['singlets_removed_by_anomaly_or_chirality_test']


def test_actual_current_selection_matches_both_theorem_sectors(report):
    assert [x['component_mask'] for x in report['admissible_sectors']] == [227, 796]
    odd, even = [set(row['occupation_masks']) for row in report['admissible_sectors']]
    assert {31 ^ mask for mask in odd} == even
    assert odd.isdisjoint(even) and len(odd) == len(even) == 15
    assert all(mask.bit_count() % 2 for mask in odd)
    assert all(not mask.bit_count() % 2 for mask in even)


def test_anomaly_cancellation_cannot_select_family_count(report):
    for count in (1, 2, 3, 4, 5, 17):
        for sector in report['admissible_sectors']:
            assert [count*a for a in sector['four_anomaly_forms']] == [0]*4
            assert count*sector['pure_color_cubic'] == 0
            assert count*sector['weak_doublet_count'] % 2 == 0
    assert not report['three_families_derived']


def test_selection_does_not_ignore_a_changed_anomaly(report):
    rows = deepcopy(report['components'])
    rows[4]['four_anomaly_forms'][0] += 1
    assert m.admissible_masks(rows) == [227]
    with pytest.raises(ValueError): m.admissible_masks(rows[:-1])


def test_source_columns_are_derived_before_theorem_comparison(monkeypatch):
    original = m._lean_table
    def changed(source, name, kind=int):
        row = original(source, name, kind)
        return (row[0]+1,)+row[1:] if name == 'charge' else row
    monkeypatch.setattr(m, '_lean_table', changed)
    with pytest.raises(ArithmeticError, match='actual source differs'):
        m.component_map()
