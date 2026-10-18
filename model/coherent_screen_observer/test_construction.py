"""Independent native and hostile-certificate controls for the screen dual."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from model.observer import native
from . import construction as c, oracle


PORTS = tuple(range(5))
EDGES = ((0, 1), (0, 2), (1, 2), (2, 3), (3, 4))
# Clock region (0,1,4): all four observed-pair masks coexist with a remote
# seam that can transport a disturbance left by a previous incident repair.
INCIDENT = (EDGES[0], EDGES[1], EDGES[2], EDGES[4])


def _output(**kwargs):
    return c.output(PORTS, EDGES, INCIDENT, cap=kwargs.pop('cap', 2),
                    bits=kwargs.pop('bits', 12), **kwargs)


def _outside(result, exact):
    return [(record, r, state) for record in result.records
            for r, state in product(range(2), range(4))
            if not result.diagonal_bounds[record][r][state][0]
            <= 2*exact[record][2*state+r, 2*state+r]
            <= result.diagonal_bounds[record][r][state][1]]


def test_independent_oracle_endpoint_and_private_record_match_original_native():
    for mask, edge, coin in product(range(32), EDGES, (0, 1)):
        loads = tuple((mask >> i) & 1 for i in range(5))+(0,)*7
        outgoing, record = native.native(loads, edge, coin)
        image, current, contrast = oracle._image(mask, edge, coin)
        assert image == sum(value << i for i, value in enumerate(outgoing))
        assert current == record[1]
        assert contrast == loads[edge[0]]-loads[edge[1]]


@pytest.mark.parametrize('law,feedback,cap,bits,rho', (
    ('C', True, 1, 8, F(1, 3)), ('C', False, 2, 12, F(1, 5)),
    ('C', True, 3, 12, F(1, 3)), ('E', False, 1, 8, F(1, 5)),
    ('E', True, 2, 12, F(1, 3)), ('E', True, 3, 12, F(1, 5))))
def test_all_twenty_one_quantum_reference_outputs_enclose_native_histories(law, feedback, cap, bits, rho):
    result = _output(law=law, feedback=feedback, cap=cap, bits=bits, rho=rho)
    exact = oracle.exact_history(5, EDGES, incident=INCIDENT, law=law,
                                 cap=cap, feedback=feedback, rho=rho)
    assert set(result.records) == set(exact) and len(exact) == 21
    assert not _outside(result, exact)
    midpoint_error = F(0)
    for record in result.records:
        matrix = exact[record]
        assert matrix[2, 5] == matrix[5, 2] == result.coherences[record]/2
        assert all(matrix[a, b] == 0 for a, b in product(range(8), repeat=2)
                   if a != b and (a, b) not in ((2, 5), (5, 2)))
        for r, state in product(range(2), range(4)):
            lo, hi = result.diagonal_bounds[record][r][state]
            assert 0 <= lo <= hi <= 1
            midpoint_error += abs(matrix[2*state+r, 2*state+r]-(lo+hi)/4)/2
    assert midpoint_error <= result.arithmetic_error
    assert sum((sum(matrix[i, i] for i in range(8)) for matrix in exact.values()), F(0)) == 1
    assert sum((sum(matrix[i, i] for i in range(8)) for record, matrix in exact.items()
                if record[0] == 'timeout'), F(0)) == result.timeout
    for r in range(2):
        assert sum((sum(matrix[2*s+r, 2*s+r] for s in range(4))
                    for matrix in exact.values()), F(0)) == F(1, 2)


def test_feedback_changes_only_the_requested_coherence_and_retains_timeout_state():
    plain, active = _output(feedback=False), _output(feedback=True)
    assert plain.columns == active.columns
    assert plain.coherences['timeout', ()] == active.coherences['timeout', ()] > 0
    for record in plain.records:
        sign = -1 if record[1] and record[1][0] == 1 else 1
        assert active.coherences[record] == sign*plain.coherences[record]
    assert active.coherences['timeout', (1,)] < 0
    assert active.coherences['complete', (1, 0)] < 0
    assert active.coherences['complete', (2, 0)] > 0


def test_candidate_and_density_are_separate_nontrivial_source_inputs():
    first, other_law = _output(law='C'), _output(law='E')
    other_density = _output(rho=F(1, 5))
    assert first.columns == other_law.columns
    assert first.coherences['complete', (1, 1)] != other_law.coherences['complete', (1, 1)]
    assert first.columns != other_density.columns
    assert first.coherence == other_density.coherence


def test_missing_mask_categories_remain_zero_outputs_and_no_remote_means_no_timeout():
    result = c.output((0, 1), ((0, 1),), ((0, 1),), cap=2, bits=8)
    assert len(result.records) == 21 and result.timeout == 0
    for record, columns in result.diagonal_bounds.items():
        if record != ('complete', (3, 3)):
            assert all(lo == hi == 0 for row in columns for lo, hi in row)
            assert result.coherences[record] == 0


class IntAlias(int):
    pass


class StrAlias(str):
    pass


class FractionAlias(F):
    pass


@pytest.mark.parametrize('change', (
    {'cap': True}, {'cap': 0}, {'cap': 1025}, {'bits': True}, {'bits': 7}, {'bits': 53},
    {'rho': True}, {'rho': 1/3}, {'rho': float('nan')}, {'rho': F(0)}, {'rho': F(1)},
    {'rho': FractionAlias(1, 3)}, {'feedback': 1}, {'law': StrAlias('C')}, {'law': 'D'}))
def test_protocol_parameters_reject_original_type_aliases_and_out_of_domain_values(change):
    with pytest.raises(ValueError):
        _output(**change)


@pytest.mark.parametrize('numerator,denominator', ((True, 3), (1, True), (1.0, 3),
    (1, 3.0), (1, 0), (-1, -3), (2, 6), (IntAlias(1), 3)))
def test_density_rejects_noncanonical_fraction_internals(numerator, denominator):
    density = F(1, 3)
    density._numerator, density._denominator = numerator, denominator
    with pytest.raises(ValueError):
        _output(rho=density)


@pytest.mark.parametrize('mutation', ('bool_port', 'bool_endpoint', 'duplicate_edge',
    'missing_pair_clock', 'foreign_incident', 'list_menu'))
def test_complete_native_graph_and_clock_are_validated_before_caches(mutation):
    ports, edges, incident = PORTS, EDGES, INCIDENT
    if mutation == 'bool_port':
        ports = (False, *PORTS[1:])
    elif mutation == 'bool_endpoint':
        edges = ((False, 1), *EDGES[1:])
    elif mutation == 'duplicate_edge':
        edges = (EDGES[0], *EDGES)
    elif mutation == 'missing_pair_clock':
        incident = INCIDENT[1:]
    elif mutation == 'foreign_incident':
        incident = (*INCIDENT, (0, 4))
    elif mutation == 'list_menu':
        incident = list(INCIDENT)
    with pytest.raises(ValueError):
        c.output(ports, edges, incident, cap=2, bits=12)


@pytest.mark.parametrize('mutation', ('record_bool', 'record_status_alias', 'drop_timeout',
    'coherence_sign', 'coherence_zero', 'changed_population', 'zero_budget', 'fraction_slots'))
def test_public_quantum_consumer_rejects_forged_complete_certificates(mutation):
    original = _output()
    if mutation in ('record_bool', 'record_status_alias', 'drop_timeout'):
        rows = list(original.records)
        if mutation == 'record_bool':
            rows[4] = ('complete', (True, 0))
        elif mutation == 'record_status_alias':
            rows[0] = (StrAlias('complete'), rows[0][1])
        else:
            rows.pop()
        changed = original._replace(records=tuple(rows))
    elif mutation in ('coherence_sign', 'coherence_zero'):
        rows = list(original.coherence)
        index = original.records.index(('complete', (1, 0)))
        rows[index] = -rows[index] if mutation == 'coherence_sign' else F(0)
        changed = original._replace(coherence=tuple(rows))
    elif mutation == 'changed_population':
        rows = list(original.columns)
        row = list(rows[0][0])
        row[0] = (F(0), F(1))
        rows[0] = (tuple(row), rows[0][1])
        changed = original._replace(columns=tuple(rows))
    elif mutation == 'zero_budget':
        assert original.arithmetic_error > 0
        changed = original._replace(arithmetic_error=F(0))
    else:
        changed = deepcopy(original)
        changed.coherence[0]._denominator = True
    with pytest.raises(ValueError):
        c.bell_blocks(changed)


def test_returned_density_cannot_poison_a_same_hash_population_cache_key():
    # Rational hashes collide modulo Python's fixed Mersenne modulus. A public
    # mutable Fraction must never also be a private cached density parameter.
    prime = (1 << 61)-1
    density = F(prime+1, prime+3)
    assert density != F(1, 3) and hash(density) == hash(F(1, 3))
    original = _output(rho=F(1, 3), bits=13)
    original.rho._numerator = density.numerator
    original.rho._denominator = density.denominator
    changed = _output(rho=density, bits=13)
    exact = oracle.exact_history(5, EDGES, incident=INCIDENT, cap=2, rho=density)
    assert not _outside(changed, exact)
    assert changed.columns != original.columns
    with pytest.raises(ValueError):
        c.validate(original)


@pytest.mark.parametrize('mutation', ('phase', 'hidden_offdiagonal', 'population'))
def test_receipt_checks_quantum_disturbance_beyond_trace_and_record_probabilities(monkeypatch, mutation):
    from . import receipt
    original = oracle.exact_history

    def corrupted(*args, **kwargs):
        rows = original(*args, **kwargs)
        matrix = rows['complete', (1, 0)]
        trace = sum(matrix[i, i] for i in range(8))
        if mutation == 'phase':
            assert matrix[2, 5] != 0
            matrix[2, 5] *= -1
            matrix[5, 2] *= -1
        elif mutation == 'hidden_offdiagonal':
            matrix[0, 7] = matrix[7, 0] = F(1, 1000)
        else:
            donor = max(range(8), key=lambda i: matrix[i, i])
            recipient = (donor+1) % 8
            amount = matrix[donor, donor]/2
            assert amount > 0
            matrix[donor, donor] -= amount
            matrix[recipient, recipient] += amount
        assert sum(matrix[i, i] for i in range(8)) == trace
        return rows

    monkeypatch.setattr(oracle, 'exact_history', corrupted)
    with pytest.raises(ArithmeticError):
        receipt.small_controls()


def test_wrong_positive_integer_action_with_the_same_mass_fails_native_oracle(monkeypatch):
    from . import receipt
    from model.observer_time_fields import dual
    original = c._step
    permutations = {}
    for degree in (1, 2):
        labels = dual.states(dual.build(PORTS, EDGES, INCIDENT, degree))
        permutations[len(labels)] = tuple(labels.index(tuple(sorted(
            1 if p == 0 else 0 if p == 1 else p for p in label))) for label in labels)

    def exchanged(*args):
        genuine = original(*args)
        # Insert an undeclared root SWAP, consistently in both exclusion
        # degrees. This is positive, integral and preserves every column mass.
        changed = genuine[list(permutations[len(genuine)])]
        assert np.all(changed >= 0) and int(changed.sum()) == int(genuine.sum())
        return changed

    c._population_data.cache_clear()
    try:
        with monkeypatch.context() as patch:
            patch.setattr(c, '_step', exchanged)
            changed = _output(bits=16)
            exact = oracle.exact_history(5, EDGES, incident=INCIDENT, cap=2)
            assert _outside(changed, exact)
            with pytest.raises(ArithmeticError):
                receipt.small_controls()
    finally:
        c._population_data.cache_clear()


@pytest.mark.parametrize('field,value', (
    ('records', None), ('records', (None,)), ('records', (('complete', [0, 0]),)),
    ('columns', {}), ('columns', (None,)), ('columns', ((None, None),)),
    ('coherence', []), ('coherence', (1,)), ('arithmetic_error', 0), ('timeout', False)))
def test_malformed_certificate_containers_raise_value_error(field, value):
    with pytest.raises(ValueError):
        c.bell_blocks(_output()._replace(**{field: value}))


@pytest.mark.parametrize('field', ('rho', 'arithmetic_error', 'timeout', 'coherence'))
def test_missing_fraction_slots_are_rejected_as_invalid_exact_input(field):
    broken = F(1, 3)
    del broken._numerator
    with pytest.raises(ValueError):
        if field == 'rho':
            _output(rho=broken)
        elif field == 'coherence':
            original = _output()
            c.bell_blocks(original._replace(coherence=(broken, *original.coherence[1:])))
        else:
            c.bell_blocks(_output()._replace(**{field: broken}))
