from fractions import Fraction as F

import numpy as np
import pytest

from . import construction as c


@pytest.fixture(scope='module')
def laws():
    return c.joint_law(), c.joint_law(effective=True)


def test_complete_joint_law_and_all_terminal_states(laws):
    receipt = c.compare_laws(*laws)
    assert receipt['record_menu_size'] == 5768
    assert receipt['positive_records'] == 689
    assert receipt['zero_records_retained'] == 5079
    assert receipt['timeout_probability'] == '369/625'
    assert receipt['total_probability'] == '1'
    assert receipt['local_restrictions_and_overlap_equal_at_every_positive_prefix']


def test_another_cap_and_passive_protocol_preserve_their_own_full_laws():
    receipt = c.compare_laws(c.joint_law(cap=1, feedback=False),
                            c.joint_law(cap=1, feedback=False, effective=True))
    assert receipt['timeout_probability'] == '21/25'


def test_timeout_postselection_is_rejected(laws):
    altered = {key: value for key, value in laws[1].items() if key[2] == 'complete'}
    with pytest.raises(ArithmeticError, match='complete record menus'):
        c.compare_laws(laws[0], altered)


def test_matching_probabilities_do_not_hide_terminal_state_corruption(laws):
    altered = dict(laws[1])
    key = next(key for key, (state, _) in altered.items() if sum(state) > 0)
    state, denominator = altered[key]
    state = state.copy()
    state[0] += 1
    state[1] -= 1
    altered[key] = state, denominator
    assert sum(state) == sum(laws[1][key][0])
    with pytest.raises(ArithmeticError, match='terminal population'):
        c.compare_laws(laws[0], altered)


def test_finite_cap_uses_transported_positive_cone(laws):
    maps, _, _, _ = c.effective_maps(2)
    assert any(value < 0 for matrix in maps.values() for value in matrix.flat)
    d = c.source()['D']
    for state, _ in laws[1].values():
        assert all(value >= 0 for value in d@state)


def test_ordinary_basis_preparation_needs_transformed_boundary():
    data = c.source()
    basis = c._zero()
    basis[0] = 1
    assert not np.array_equal(data['D']@basis, 24*basis)
    # A different preparation cannot use the uniform-state boundary unchanged.
    direct, dd, _, _ = c.direct_interval(basis, 2)
    maps, ed, _, _ = c.effective_maps(2)
    assert any(F(sum(direct[mark]), dd) != F(sum(matrix@basis), ed)
               for mark, matrix in maps.items())


def test_two_equal_marginals_do_not_determine_joint_clocks():
    control = c.clock_control()
    assert control['both_local_marginals'] == '3/4'
    assert control['joint_total_variation_from_independent_clocks'] == '1/8'
    assert control['first_interval_clock_output_TV_including_timeout'] == '2/25'


def test_zero_current_is_still_a_local_incident_tick():
    initial = (0, 0, 0)
    mark = ((1, 0), (0, 0, 0))
    word = (('scheduled', mark, initial),)
    assert c.decode_local(c.observer_view(initial, word, 0)) == ((0, 0), 1, None)
    assert c.decode_local(c.observer_view(initial, word, 1)) == ((0, 0), 0, None)


def test_false_nonincident_label_is_rejected():
    word = (('scheduled', ((0, 1), (1, 0, 0)), (1, 0, 0)),)
    with pytest.raises(ArithmeticError, match='nonincident'):
        c.observer_view((0, 0, 0), word, 0)


def test_entry_records_retain_initial_occupation_and_censoring(laws):
    masses = c.local_entry_records(laws[0])
    assert sum(map(F, masses.values())) == 1
    assert masses['initially_occupied'] == '1/6'
    assert F(masses['no_entry_timeout']) > 0
    assert any(key.startswith('first_entry_at_B_incident_tick_') for key in masses)


def test_same_source_geometry_and_action():
    result = c.geometry_and_action()
    assert result['source_C_action_branches_checked'] == 3960
    assert result['all_actual_seam_count_matrices_equal_I_plus_swap']
    assert result['dilation_dimension'] == 20
    assert result['effective_spectral_lower_bound'] == '1/20'
    assert not result['native_repairs_equal_unitary_continuation_of_H']


def test_boundary_straddling_feedback_is_outside_commuting_contract():
    data = c.source()
    index = {state: i for i, state in enumerate(data['states'])}
    violations = 0
    for seam in data['seams']:
        if not set(seam).intersection(c.UNION) or set(seam).issubset(c.UNION):
            continue
        matrix = np.zeros((66, 66), dtype=object)
        for column, state in enumerate(data['states']):
            for coin in (0, 1):
                output, _ = c.native.native(state, seam, coin)
                matrix[index[output], column] += 1
        violations += not np.array_equal(data['D']@matrix, matrix@data['D'])
    assert violations > 0


@pytest.mark.parametrize('value', [True, 0, 4, F(2), '2'])
def test_invalid_cap_is_rejected(value):
    with pytest.raises(ValueError):
        c.joint_law(cap=value)
