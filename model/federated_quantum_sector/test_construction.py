"""Independent finite controls for quantum reduction, records and locality."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import construction as c
from ..federated_observer import source


def test_complete_local_certificate_and_two_candidate_difference():
    report = c.local_certificate()
    assert report['all_outcome_matrix_units_checked'] == 96
    # Coherence between occupation on one endpoint and an untouched spectator:
    # the current record retains half; the imbalance record erases all of it.
    assert c.native_unit('C', 1, 4, (0, 1)) == {(1, 4): F(1, 2)}
    assert c.native_unit('E', 1, 4, (0, 1)) == {}


@pytest.mark.parametrize('law', ('C', 'E'))
def test_tensor_local_remote_trace_on_every_matrix_unit(law):
    # Bits0,1 form the operation's work factor; bit2 is an arbitrary reference.
    for left, right in product(range(8), repeat=2):
        actual = {}
        for (a, b), value in c.native_unit(law, left, right, (0, 1)).items():
            if a % 4 == b % 4:
                key = a//4, b//4
                actual[key] = actual.get(key, F(0))+value
        expected = {(left//4, right//4): F(1)} if left % 4 == right % 4 else {}
        assert actual == expected


@pytest.mark.parametrize('law', ('C', 'E'))
def test_disjoint_instruments_commute_outcome_by_outcome(law):
    def compose(left, right, first, second, j, k):
        result = {}
        for (a, b), v in c.native_unit(law, left, right, first, j).items():
            for key, w in c.native_unit(law, a, b, second, k).items():
                result[key] = result.get(key, F(0))+v*w
        return result
    for left, right in product(range(16), repeat=2):
        for j, k in product(c.CURRENT, repeat=2):
            assert compose(left, right, (0, 1), (2, 3), j, k) == compose(
                left, right, (2, 3), (0, 1), k, j)


def test_conditional_record_can_steer_remote_state_without_signalling():
    # Half |100><100| and half |001><001|, written as low-first work bits
    # plus remote bit2: a current+1 outcome can only come from the latter.
    initial = {(4, 4): F(1, 2), (1, 1): F(1, 2)}
    remote = {}
    for (i, j), value in initial.items():
        for (a, b), weight in c.native_unit('C', i, j, (0, 1), 1).items():
            if a % 4 == b % 4:
                remote[a//4, b//4] = remote.get((a//4, b//4), F(0))+value*weight
    assert remote == {(0, 0): F(1, 4)}
    assert remote != {(0, 0): F(1, 8), (1, 1): F(1, 8)}


@pytest.mark.parametrize('law', ('C', 'E'))
def test_general_source_block_with_common_particles_and_sector_changes(law):
    graph = source.federation(0)
    configurations = (0, 1, 2, 3, (1 << 1)|(1 << 82), (1 << 4)|(1 << 30)|(1 << 82))
    for left, right in zip(configurations, reversed(configurations), strict=True):
        assert c.source_column(graph, left, right, law) == c.sector_column(graph, left, right, law)
    # The full binary algebra has a genuine global fixed coherence; fixed-q
    # uniqueness does not extend to arbitrary coherent number-sector mixtures.
    full = (1 << graph.ports)-1
    assert c.source_column(graph, 0, full, law) == {(0, full): F(1)}
    assert c.block_killing(graph, (1,)*graph.ports, law) == 0


@pytest.mark.parametrize('law', ('C', 'E'))
def test_native_only_records_ignore_phase_even_with_continuing_interventions(law):
    rows = [c.execute_experiment(law, sign=sign, mode=mode)
            for sign, mode in ((1, 'identity'), (-1, 'identity'), (1, 'phase'), (1, 'dephase'))]
    assert len({r['current_only_law_sha256'] for r in rows}) == 1
    assert all(r['complete_record_words'] == 162 and r['exact_total_probability'] == '1' for r in rows)
    expected = ('5/8', '3/8', '3/8', '1/2') if law == 'C' else ('1/2',)*4
    assert tuple(r['exact_plus_probability'] for r in rows) == expected


@pytest.mark.parametrize('law', ('C', 'E'))
def test_complete_reference_outputs_and_serialized_continuation(law):
    result = c.execute_experiment(law, reference=True)
    assert result['untouched_reference_dimension'] == 2
    assert result['exact_plus_probability'] == '1/2'
    assert result['complete_state_and_checkpoint_residual_diagnostic'] < 1e-12
    program, exact, initial = c.experiment(law, reference=True)
    assert len(initial) == 8 and len(exact) == 162
    # Retained reference has maximally mixed marginal across every terminal
    # output together, not a fresh uncorrelated replacement on each branch.
    total = np.asarray(sum((np.asarray(m, dtype=object) for m in exact.values())), dtype=object)
    remote = [[sum(total[2*i+r, 2*i+s] for i in range(4)) for s in range(2)] for r in range(2)]
    assert remote == [[F(1, 2), F(0)], [F(0), F(1, 2)]]


def test_uniform_source_rates_and_exact_outward_probability_rounding():
    assert c.exchange_rate(0, 1, 82, 'C') == F(104, 105)
    assert c.exchange_rate(0, 1, 82, 'E') == F(619, 630)
    for rate, count in product((F(0), F(1), F(4, 7), F(104, 105)), (0, 1, 5, 631)):
        bounds = c.phase_interval(rate, count)
        truth = (1+rate**count)/2
        assert F(bounds['lower']) <= truth <= F(bounds['upper'])
        assert F(bounds['upper'])-F(bounds['lower']) <= F(1, 10**12)


@pytest.mark.parametrize('bad', (True, -1, 1.0, '1', 1 << 4096))
def test_mask_boundary(bad):
    with pytest.raises(ValueError):
        c.native_image(bad, (0, 1), 0)


@pytest.mark.parametrize('seam', ((1, 0), (0, 0), (True, 1), (0, 4096), [0, 1]))
def test_seam_boundary(seam):
    with pytest.raises(ValueError):
        c.native_image(1, seam, 0)


def test_strict_candidate_record_and_reporting_boundaries():
    for law in ('c', None, True):
        with pytest.raises(ValueError):
            c.native_unit(law, 0, 0, (0, 1))
    for value in (True, 2, 0.0):
        with pytest.raises(ValueError):
            c.native_unit('C', 0, 0, (0, 1), value)
    for rate, count in ((.5, 1), (F(1, 2), True), (F(2), 1), (F(1), 20001)):
        with pytest.raises(ValueError):
            c.phase_interval(rate, count)
    with pytest.raises(ValueError):
        c.block_killing(source.federation(0), (True,)*240, 'C')


def test_readout_keeps_complement_instead_of_postselection():
    vertices = (0, 1, 2, 3)
    initial = c._zero(4); initial[2][2] = F(1)
    outputs = c._exact_apply(c.analyzer(vertices, (0, 1)), initial, 1)
    assert outputs[0][2][2] == outputs[1][2][2] == F(1, 2)
