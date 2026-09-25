"""Independent moment and RNG checks with adversarial false-green cases."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest


SPEC = importlib.util.spec_from_file_location("source_innovation_extract", Path(__file__).with_name("extract.py"))
E = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(E)


def receipt_row(x, f, y):
    # Fit using least-squares QR/SVD, independent of the inverse normal-equation
    # formulas in the producer. Include arbitrary nonzero means.
    design = np.column_stack([np.ones(len(x)), x, f])
    coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
    simple = np.linalg.lstsq(design[:, :2], y, rcond=None)[0]
    residual = y - design @ coefficients
    simple_residual = y - design[:, :2] @ simple
    return {
        "response_b": simple[1], "joint_response_b": coefficients[1],
        "joint_fresh_c": coefficients[2], "total_power_ratio": y.var() / x.var(),
        "fresh_power_over_inherited": f.var() / x.var(),
        "innovation_power_ratio": np.mean(simple_residual**2) / x.var(),
        "joint_innovation_power_ratio": np.mean(residual**2) / x.var(),
        "squared_coherence": np.corrcoef(x, y)[0, 1]**2,
        "cells_per_group": 4, "groups": len(x),
    }


@pytest.mark.parametrize("mode", ["chain", "shuffled"])
@pytest.mark.parametrize("chunk", [64, 256, 1024])
def test_stream_replay_matches_full_original_rng_call(mode, chunk):
    carriers, seed, sizes = 1024, 129, [4, 16, 64]
    rng = np.random.default_rng(seed)
    if mode == "shuffled":
        rng.permutation(carriers // 4)
    full = rng.integers(0, 5, size=carriers * 12)
    moments, digest = E.replay_fresh(carriers, seed, mode, sizes, chunk)
    import hashlib
    assert digest == hashlib.sha256(full.tobytes()).hexdigest()
    for size in sizes:
        # Original producer does cell means then group means, not integer sums.
        grouped = full.reshape(-1, 12).mean(axis=1).reshape(-1, size).mean(axis=1)
        assert moments[size]["fresh_variance"] == pytest.approx(grouped.var(), rel=1e-13)
        assert moments[size]["fresh_mean"] == pytest.approx(grouped.mean(), abs=1e-14)


def test_covariance_recovered_from_correlated_shifted_observations():
    rng = np.random.default_rng(4371)
    x = 40 + 7 * rng.normal(size=90)
    f = -12 + .35 * x + rng.normal(size=90)
    y = 300 + .2 * x - .7 * f + rng.normal(size=90)
    row = receipt_row(x, f, y)
    result = E.recover(row, f.var())
    direct = np.cov([x, f, y], bias=True)
    assert np.allclose(result["covariance_native_load_squared"], direct, atol=1e-10)
    b = np.cov(x, y, bias=True)[0, 1] / x.var()
    residual = (y - y.mean()) - b * (x - x.mean())
    assert result["inherited_linear_residual_variance"] == pytest.approx(residual.var())
    # A slope is neither residual variance nor retained covariance power.
    assert result["inherited_linear_residual_variance"] != pytest.approx(b * b)


def test_common_amplitude_is_unidentifiable_from_ratios_without_fresh_anchor():
    rng = np.random.default_rng(55)
    x, f, noise = rng.normal(size=(3, 80))
    y = .8*x + .4*f + noise
    original = receipt_row(x, f, y)
    rescaled = receipt_row(17*x, 17*f, 17*y)
    for key in original:
        assert original[key] == pytest.approx(rescaled[key])
    a = E.recover(original, f.var())
    b = E.recover(rescaled, (17*f).var())
    assert np.allclose(np.array(b["covariance_native_load_squared"]),
                       17**2 * np.array(a["covariance_native_load_squared"]))


def test_inherited_only_control_has_zero_linear_innovation():
    rng = np.random.default_rng(101)
    x, f = rng.normal(size=(2, 100))
    y = x.copy()
    row = receipt_row(x, f, y)
    # C=0 makes Cov(x,f) not inferable from b-B. Reject rather than divide
    # by a near-zero quantity or fabricate a covariance.
    with pytest.raises(ValueError, match="degenerate"):
        E.recover(row, f.var())
    assert row["innovation_power_ratio"] < 1e-26


def test_corrupted_innovation_or_coherence_receipt_cannot_pass():
    rng = np.random.default_rng(18)
    x, f, noise = rng.normal(size=(3, 50))
    y = .5 * x + .7 * f + noise
    original = receipt_row(x, f, y)
    for field in ["innovation_power_ratio", "joint_innovation_power_ratio", "squared_coherence"]:
        bad = dict(original)
        bad[field] *= 2
        with pytest.raises(ValueError, match="inconsistency"):
            E.recover(bad, f.var())


def test_parent_replicates_are_clustered_without_iid_standard_errors():
    metrics = ["inherited_variance", "settled_variance", "inherited_linear_residual_variance",
               "joint_linear_residual_variance", "fresh_orthogonal_transferred_variance",
               "joint_response_b", "joint_fresh_c"]
    chains = []
    for parent, value in [(10, 1.), (10, 3.), (11, 20.)]:
        chains.append({"parent_level": 8, "parent_seed": parent, "mode": "chain",
                       "by_scale": [{"cells_per_group": 4, **{key: value for key in metrics}}]})
    result = E.aggregate(chains)[0]
    assert result["parent_schedule_clusters"] == 2
    assert result["independent_initial_load_realizations"] == 1
    assert result["by_scale"][0]["settled_variance"]["mean"] == 11.
    assert "se" not in result["by_scale"][0]["settled_variance"]


def test_shuffled_stream_is_not_unshuffled_stream():
    a, _ = E.replay_fresh(1024, 88, "chain", [4, 16, 64], 64)
    b, _ = E.replay_fresh(1024, 88, "shuffled", [4, 16, 64], 64)
    assert a[4]["fresh_variance"] != b[4]["fresh_variance"]
