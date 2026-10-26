"""Independent controls for the arbitrary-code recovery obstruction."""
from fractions import Fraction as F

import numpy as np
import pytest

from ..observer import native
from . import record_capacity as capacity


def test_rank15_is_separated_from_perfect_recovery():
    c, e = capacity.certificate().values()
    assert F(c["entanglement_fidelity_upper"]) == F(43, 45)
    assert F(e["entanglement_fidelity_upper"]) == F(41, 45)
    assert F(c["reference_trace_distance_lower"]) == F(2, 45)
    assert F(e["reference_trace_distance_lower"]) == F(4, 45)
    assert c["branch_rank_counts"] == {"0": 60, "10": 60, "56": 60}
    assert e["branch_rank_counts"] == {"10": 120, "46": 60}


@pytest.mark.parametrize("law", ["C", "E"])
def test_complete_partition_and_independent_source_effects(law):
    states = native.states(binary=True)
    for edge, coin, label, indices in capacity.branch_partition(law):
        a, b = edge
        expected = []
        for index, x in enumerate(states):
            # Direct floor/ceiling endpoint arithmetic, independent of native().
            record = x[a]-(x[a]+x[b]+coin)//2 if law == "C" else x[a]-x[b]
            if record == label:
                expected.append(index)
        assert indices == tuple(expected)


@pytest.mark.parametrize("law", ["C", "E"])
@pytest.mark.parametrize("dimension", [1, 6, 11, 15, 47, 66])
def test_complex_code_optimum_obeys_universal_ceiling(law, dimension):
    rng = np.random.default_rng(881+dimension)
    raw = rng.normal(size=(66, dimension))+1j*rng.normal(size=(66, dimension))
    code, _ = np.linalg.qr(raw)
    fidelity = 0.0
    for _, _, _, indices in capacity.branch_partition(law):
        if indices:
            # Singular values of P V avoid squaring a small singular value.
            singular = np.linalg.svd(code[list(indices)], compute_uv=False)
            fidelity += float(np.sum(singular)**2)/(60*dimension**2)
    bound = float(F(capacity.recovery_bound(law, dimension)["entanglement_fidelity_upper"]))
    assert 0 <= fidelity <= bound+2e-13
    if dimension == 1:
        assert abs(fidelity-1) < 2e-13  # A known one-state code is recoverable.


@pytest.mark.parametrize("law", ["C", "E"])
def test_full_carrier_matches_independent_optimal_fidelity(law):
    # V=I: each branch contributes rank(P)^2/(60*66^2), attained by
    # the endpoint inverse. Rank truncation is tight for this complete input.
    optimum = sum(F(len(indices)**2, 60*66**2)
                  for _, _, _, indices in capacity.branch_partition(law))
    bound = F(capacity.recovery_bound(law, 66)["entanglement_fidelity_upper"])
    assert optimum == bound
    assert optimum < F(4, 5)


@pytest.mark.parametrize("dimension", [True, False, 0, 67, 1.0, F(6), "15", None])
def test_invalid_dimension_refused(dimension):
    with pytest.raises(ValueError):
        capacity.recovery_bound("C", dimension)


@pytest.mark.parametrize("law", [None, True, "D", "", ["C"]])
def test_invalid_law_refused(law):
    with pytest.raises(ValueError):
        capacity.branch_partition(law)
