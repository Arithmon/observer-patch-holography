"""Independent endpoint, quantum-reference and changed-premise checks."""
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np
import pytest

from . import local_bias as m


@pytest.fixture(scope='module')
def report():
    return m.report()


def biased():
    return m.Parameters(port_weights=(25, 25)+(8,)*10)


def test_exact_native_coins_realize_independently_defined_conditional_reference():
    p = biased(); pi = m.reference(p); states = m.basis()
    for edge in m.SEAMS:
        kernel = m.seam_kernel(p, edge)
        outside = tuple(i for i in range(12) if i not in edge)
        for i, x in enumerate(states):
            fibre = [j for j, y in enumerate(states) if all(x[k] == y[k] for k in outside)]
            mass = sum(pi[j] for j in fibre)
            assert all(kernel[i][j] == (pi[j]/mass if j in fibre else 0) for j in range(66))
    assert m.coin_probability(p, (0, 5), 1) == F(25, 33)
    assert m.coin_probability(p, (0, 5), 0) == F(8, 33)


@pytest.mark.parametrize('law', ('C', 'E'))
@pytest.mark.parametrize('binary', (False, True))
def test_all_raw_branches_preserve_arbitrary_reference_mass(law, binary):
    p = biased(); states = m.basis(binary)
    rng = np.random.default_rng(1135)
    v = rng.normal(size=(len(states), 3))+1j*rng.normal(size=(len(states), 3))
    v /= np.linalg.norm(v)
    before = v.T@v.conj()
    for edge in m.SEAMS:
        after = np.zeros((3, 3), complex)
        for coin, value, pairs in m.raw_maps(law, edge, binary):
            out = np.zeros_like(v)
            for i, j in pairs:
                out[j] = np.sqrt(float(m.coin_probability(p, edge, coin)))*v[i]
            after += out.T@out.conj()
        assert np.linalg.norm(after-before) < 2e-15


@pytest.mark.parametrize('law', ('C', 'E'))
def test_modular_intertwiner_on_complex_coherent_input(law):
    p = biased(); pi = m.reference(p)
    phases = np.exp(1.4j*np.log(np.array(list(map(float, pi)))))
    rng = np.random.default_rng(72)
    v = rng.normal(size=(66, 2))+1j*rng.normal(size=(66, 2))
    for edge in m.SEAMS:
        for coin, value, pairs in m.raw_maps(law, edge):
            current = m.records.record_current(law, (*edge, coin, value))
            ratio = F(p.port_weights[edge[1]], p.port_weights[edge[0]])**current
            expected = np.exp(1.4j*np.log(float(ratio)))
            # sigma^it K sigma^-it = ratio^it K on every native column.
            for i, j in pairs:
                assert np.linalg.norm((phases[j]/phases[i]-expected)*v[i]) < 1e-14


def test_nonuniform_common_reference_and_law_endpoints(report):
    row = report['biased_control']
    pi = tuple(map(F, row['reference']))
    assert len(set(pi)) == 3 and min(pi) > 0 and sum(pi) == 1
    assert list(map(F, report['biased_reference_sector_masses'])) == [F(576, 1501), F(800, 1501), F(125, 1501)]
    assert row['classical_stationary_law_unique']
    assert not row['full_quantum_stationary_uniqueness_claimed']
    assert row['regional_factorization']['all_bare_two_port_algebras_invariant']
    for law in ('C', 'E'):
        assert row['quantum_instruments'][law]['binary_columns'] == 3960
        assert row['quantum_instruments'][law]['H78_and_H66_complete']
    assert report['quantum_distinction']['unread_output_matrix_unit_coefficient'] == {'C': '1/2', 'E': '0'}


def test_holds_mixtures_and_fixed_nonuniform_selection_preserve_reference(report):
    assert report['biased_control']['reference'] == report['held_nonuniform_scheduler_control']['reference']
    p = biased(); p_c = m.Parameters(port_weights=p.port_weights, c_weight=F(1))
    p_e = m.Parameters(port_weights=p.port_weights, c_weight=F(0))
    assert m.kernel(p_c) == m.kernel(p_e)
    assert set(map(F, report['uniform_control']['reference'])) == {F(1, 66)}


def test_unchanged_fair_coin_is_a_real_negative_control():
    pi = m.reference(biased()); fair = m.kernel(m.Parameters())
    after = tuple(sum(pi[i]*fair[i][j] for i in range(66)) for j in range(66))
    assert sum(after) == 1 and after != pi


def test_seam_reversal_changes_coin_convention_but_not_population_kernel():
    p = biased()
    assert m.seam_kernel(p, (0, 5)) == m.seam_kernel(p, (5, 0))


def test_two_actual_spectators_reject_fixed_old_reference_local_detailed_balance(report):
    rows = report['old_reference_boundary']['all_66_core_witnesses']
    assert len(rows) == 66
    for row in rows:
        core, (a, b) = set(row['core']), row['directed_seam']
        assert tuple(sorted((a, b))) in m.SEAMS
        ratios = []
        for spectator in row['spectators']:
            x = tuple(int(i in (a, spectator)) for i in range(12))
            y = tuple(int(i in (b, spectator)) for i in range(12))
            # The same local endpoint input10 and output01, differing only outside.
            assert (x[a], x[b], y[a], y[b]) == (1, 0, 0, 1)
            old_x = m.OLD_WEIGHTS[sum(x[i] for i in core)]
            old_y = m.OLD_WEIGHTS[sum(y[i] for i in core)]
            ratios.append(old_y/old_x)
        assert ratios == [F(5, 8), F(8, 25)]
    assert report['old_reference_boundary']['nonproduct_ratio_w1_squared_over_w0_w2'] == '125/64'


@pytest.mark.parametrize('arguments', [
    {'port_weights': (True,)+(1,)*11}, {'port_weights': (1,)*11},
    {'port_weights': (2,)*12}, {'port_weights': (5,)+(1,)*11},
    {'port_weights': (33,)+(9,)*11}, {'seam_weights': (0,)+(1,)*29},
    {'repair_probability': F(1, 17)}, {'repair_probability': True},
    {'c_weight': .5}, {'c_weight': F(1, 257)},
])
def test_frozen_local_family_rejects_outside_parameters(arguments):
    with pytest.raises(ValueError):
        m.Parameters(**arguments)


def test_contract_does_not_claim_full_issue_closure(tmp_path):
    contract = m.contract()
    assert len(contract['full_issue_acceptance_missing']) == 4
    for field in ('CP00_accepted', 'whole_program_frozen'):
        wrong = dict(contract); wrong[field] = True
        path = tmp_path/'wrong.json'; path.write_text(json.dumps(wrong))
        with pytest.raises(ValueError):
            m.contract(path)
    wrong = dict(contract); wrong['full_issue_acceptance_missing'] = []
    path.write_text(json.dumps(wrong))
    with pytest.raises(ValueError):
        m.contract(path)
    wrong = json.loads(json.dumps(contract))
    wrong['parameter_boundary']['port_weights']['integer_min'] = True
    path.write_text(json.dumps(wrong))
    with pytest.raises(ValueError):
        m.contract(path)


def test_theorem_names_match_bound_sources(report):
    row = report['theorem_reuse'][1]
    text = (m.ROOT/row['path']).read_text()
    assert all('theorem '+name in text for name in row['theorems'])


def test_strict_retained_report(report):
    assert m.strict_json.load(m.OUTPUT) == report
    assert report['scope']['CP00_accepted'] is False
