"""Check the nested covariance identity using direct full-array projections."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location("source_nested", Path(__file__).with_name("nested.py"))
N = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(N)


@pytest.mark.parametrize("predictors", [1, 2])
def test_nested_covariance_against_direct_lifted_ols_fields(predictors):
    rng = np.random.default_rng(723)
    count = 4096
    x, f, noise = rng.normal(size=(3, count))
    y = .5*x + .2*f + .7*x*x + noise
    sizes = [4, 16, 64, 256]
    lifted = []
    for size in sizes:
        averaged = [z.reshape(-1, size).mean(axis=1) for z in [x, f, y]]
        a, b, response = averaged
        design = np.column_stack([np.ones(len(a)), a, b])[:, :predictors+1]
        fitted = design @ np.linalg.lstsq(design, response, rcond=None)[0]
        lifted.append(np.repeat(response-fitted, size))
    direct = np.cov(lifted, bias=True)
    covariance, increments = N.covariance_from_nested_residuals(np.diag(direct))
    assert np.allclose(covariance, direct, atol=1e-14)
    assert increments.sum() == pytest.approx(direct[0, 0])
    # Increments are orthogonal, even with scale-dependent regression slopes.
    fields = np.asarray(lifted)
    differences = np.vstack([fields[:-1]-fields[1:], fields[-1]])
    assert np.allclose(np.cov(differences, bias=True), np.diag(increments), atol=1e-14)


def test_nonmonotonic_variance_does_not_fake_a_valid_covariance():
    with pytest.raises(ValueError, match="decrease"):
        N.covariance_from_nested_residuals([.4, .2, .3])


def test_check_is_read_only_and_rejects_changed_receipt(tmp_path):
    # A tiny complete transform fixture makes this independent of the archive.
    source = {"chains": [{"source_path": "fixture", "parent_level": 8,
                           "mode": "chain", "parent_seed": 1, "fresh_seed": 2,
                           "by_scale": [
                               {"cells_per_group": 4,
                                "inherited_linear_residual_variance": .4,
                                "joint_linear_residual_variance": .1},
                               {"cells_per_group": 16,
                                "inherited_linear_residual_variance": .2,
                                "joint_linear_residual_variance": .03}]}]}
    (tmp_path / "receipt.json").write_text(json.dumps(source))
    (tmp_path / "nested_spec.json").write_text('{"fixture": true}')
    N.main([], tmp_path)
    path = tmp_path / "nested_receipt.json"
    before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in tmp_path.iterdir()}
    N.main(["--check"], tmp_path)
    after = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in tmp_path.iterdir()}
    assert before == after
    bad = json.loads(path.read_text())
    bad["chains"][0]["inherited_linear"]["same_location_nested_covariance"][0][0] = 999.
    path.write_text(json.dumps(bad))
    poisoned = path.read_bytes(), path.stat().st_mtime_ns
    with pytest.raises(SystemExit, match="differs"):
        N.main(["--check"], tmp_path)
    assert (path.read_bytes(), path.stat().st_mtime_ns) == poisoned


def test_unknown_verification_option_is_rejected_before_writes(tmp_path):
    with pytest.raises(SystemExit):
        N.main(["--verify"], tmp_path)
    assert not list(tmp_path.iterdir())
