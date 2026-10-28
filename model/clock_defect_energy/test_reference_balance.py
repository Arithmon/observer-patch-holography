"""Independent complete-instrument and finite-action checks."""
from fractions import Fraction as F

from ..completion import observer
from . import reference_balance as candidate
from .joins_contract import check_receipt


def test_complete_h78_kraus_maps_independently_reconstruct_native_kernel():
    index = {state: i for i, state in enumerate(candidate.continuation.BASIS)}
    for law in ("C", "E"):
        counts = [[0]*66 for _ in range(66)]
        for seam in observer.SEAMS:
            for _, _, _, pairs in observer._raw_maps(law, seam):
                for i, j in pairs:
                    state = observer.STATES[i]
                    if state in index:
                        counts[index[state]][index[observer.STATES[j]]] += 1
        assert candidate.kernel(law) == tuple(tuple(F(x, 60) for x in row) for row in counts)


def test_fixed_clock_reference_drift_and_positive_form():
    report = candidate.report()
    assert report["native_states"] == 66
    assert report["native_movement_edges"] == 300
    assert report["native_uniform_reference_unique_by_connected_reciprocity"]
    one, other = report["results"]
    assert one["core"] == [0, 5] and other["core"] == [0, 8]
    assert [F(x["one_attempt_total_variation_from_fixed_reference"])
            for x in (one, other)] == [F(17, 750), F(17, 600)]
    assert [x["native_edges_violating_fixed_reference_detailed_balance"]
            for x in (one, other)] == [80, 100]
    assert [F(x["mode_quadratic_change_over_variance"])
            for x in (one, other)] == [F(710, 6983), F(8515, 83796)]
    assert all(F(x["mode_quadratic_change_at_fixed_reference"]) > 0
               and not x["fixed_reference_stationary"] for x in (one, other))
    assert report["scope"]["branchwise_transported_reference_covariance_preserved"]
    assert not report["scope"]["integrated_gate_accepted"]


def test_stored_exact_report_replays():
    check_receipt(candidate.OUTPUT, candidate.report())
