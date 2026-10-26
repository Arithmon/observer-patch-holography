"""Independent controls for the arbitrary-code recovery obstruction."""
from fractions import Fraction as F

import numpy as np
import pytest

from ..observer import native
from . import record_capacity as capacity


def test_rank15_is_separated_from_perfect_recovery():
    report = capacity.certificate()
    c, e = report["C"], report["E"]
    assert c == capacity.recovery_bound("C", 15, binary=True)
    assert e == capacity.recovery_bound("E", 15, binary=True)
    assert F(c["entanglement_fidelity_upper"]) == F(43, 45)
    assert F(e["entanglement_fidelity_upper"]) == F(41, 45)
    assert F(c["reference_trace_distance_lower"]) == F(2, 45)
    assert F(e["reference_trace_distance_lower"]) == F(4, 45)
    assert c["branch_rank_counts"] == {"0": 60, "10": 60, "56": 60}
    assert e["branch_rank_counts"] == {"10": 120, "46": 60}
    c, e = report["full_total_two"]["C"], report["full_total_two"]["E"]
    assert c["carrier_dimension"] == e["carrier_dimension"] == 78
    assert F(c["entanglement_fidelity_upper"]) == F(217, 225)
    assert F(e["entanglement_fidelity_upper"]) == F(41, 45)
    assert F(c["reference_trace_distance_lower"]) == F(8, 225)
    assert F(e["reference_trace_distance_lower"]) == F(4, 45)
    assert c["branch_rank_counts"] == {"1": 60, "11": 60, "66": 60}
    assert e["branch_rank_counts"] == {"1": 120, "10": 120, "56": 60}


@pytest.mark.parametrize("law", ["C", "E"])
@pytest.mark.parametrize("binary", [True, False])
def test_complete_partition_and_independent_source_effects(law, binary):
    states = native.states(binary=binary)
    partitions = {}
    for edge, coin, label, indices in capacity.branch_partition(law, binary=binary):
        a, b = edge
        expected = []
        endpoints = []
        for index, x in enumerate(states):
            # Direct floor/ceiling endpoint arithmetic, independent of native().
            record = x[a]-(x[a]+x[b]+coin)//2 if law == "C" else x[a]-x[b]
            if record == label:
                expected.append(index)
                end = list(x)
                end[a] = (x[a]+x[b]+coin)//2
                end[b] = x[a]+x[b]-end[a]
                endpoints.append(tuple(end))
        assert indices == tuple(expected)
        assert len(set(endpoints)) == len(indices)
        assert all(end in states for end in endpoints)
        partitions.setdefault((edge, coin), {})[label] = indices
    labels = {-2, -1, 0, 1, 2} if law == "E" and not binary else {-1, 0, 1}
    assert len(partitions) == 60
    for rows in partitions.values():
        assert set(rows) == labels
        assert sorted(index for indices in rows.values() for index in indices) == list(range(len(states)))
        if not binary:
            assert all(rows.values())


@pytest.mark.parametrize("law", ["C", "E"])
@pytest.mark.parametrize("binary,dimension", [(True, k) for k in (1, 6, 11, 15, 47, 66)]
                         + [(False, 15), (False, 78)])
def test_complex_code_optimum_obeys_universal_ceiling(law, binary, dimension):
    carrier_dimension = len(native.states(binary=binary))
    rng = np.random.default_rng(881+dimension)
    raw = (rng.normal(size=(carrier_dimension, dimension))
           + 1j*rng.normal(size=(carrier_dimension, dimension)))
    code, _ = np.linalg.qr(raw)
    fidelity = 0.0
    for _, _, _, indices in capacity.branch_partition(law, binary=binary):
        if indices:
            # Singular values of P V avoid squaring a small singular value.
            singular = np.linalg.svd(code[list(indices)], compute_uv=False)
            fidelity += float(np.sum(singular)**2)/(60*dimension**2)
    bound = float(F(capacity.recovery_bound(law, dimension, binary=binary)["entanglement_fidelity_upper"]))
    assert 0 <= fidelity <= bound+2e-13
    if dimension == 1:
        assert abs(fidelity-1) < 2e-13  # A known one-state code is recoverable.
    if dimension == carrier_dimension:
        assert abs(fidelity-bound) < 2e-13  # Any complete orthonormal basis attains it.


@pytest.mark.parametrize("law", ["C", "E"])
@pytest.mark.parametrize("binary", [True, False])
def test_full_carrier_matches_independent_optimal_fidelity(law, binary):
    dimension = len(native.states(binary=binary))
    # V=I: each branch contributes rank(P)^2/(60*dimension^2), attained by
    # the endpoint inverse. Rank truncation is tight for this complete input.
    optimum = sum(F(len(indices)**2, 60*dimension**2)
                  for _, _, _, indices in capacity.branch_partition(law, binary=binary))
    bound = F(capacity.recovery_bound(law, dimension, binary=binary)["entanglement_fidelity_upper"])
    assert optimum == bound
    assert optimum < F(4, 5)


@pytest.mark.parametrize("binary", [True, False])
@pytest.mark.parametrize("dimension", [True, False, 0, 79, 1.0, F(6), "15", None, np.int64(15)])
def test_invalid_dimension_refused(dimension, binary):
    with pytest.raises(ValueError):
        capacity.recovery_bound("C", dimension, binary=binary)


def test_domain_specific_dimensions():
    with pytest.raises(ValueError):
        capacity.recovery_bound("C", 67)
    assert capacity.recovery_bound("C", 67, binary=False)["logical_dimension"] == 67


@pytest.mark.parametrize("law", [None, True, "D", "", ["C"]])
def test_invalid_law_refused(law):
    with pytest.raises(ValueError):
        capacity.branch_partition(law)
    with pytest.raises(ValueError):
        capacity.recovery_bound(law, 15, binary=False)


@pytest.mark.parametrize("binary", [0, 1, None, "False", [], np.bool_(False)])
def test_invalid_binary_flags_refused(binary):
    with pytest.raises(ValueError):
        capacity.branch_partition("C", binary=binary)
    with pytest.raises(ValueError):
        capacity.recovery_bound("E", 15, binary=binary)


def test_binary_flag_is_keyword_only():
    with pytest.raises(TypeError):
        capacity.branch_partition("C", False)
    with pytest.raises(TypeError):
        capacity.recovery_bound("E", 15, False)
