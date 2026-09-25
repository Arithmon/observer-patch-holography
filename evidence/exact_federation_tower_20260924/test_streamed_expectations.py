"""The bounded verifier must exactly preserve the original seeded load stream."""
import hashlib
import importlib.util
from pathlib import Path
from fractions import Fraction
import numpy as np
import pytest

_spec = importlib.util.spec_from_file_location("tower_streamed", Path(__file__).with_name("verify_tower.py"))
tower = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tower)


@pytest.mark.parametrize("level,ports", [(0, 240), (2, 3840), (7, 24001)])
@pytest.mark.parametrize("chunk", [1, 3, 128, 1025, 100000])
def test_streamed_loads_equal_one_original_call(level, ports, chunk):
    loads = np.random.default_rng(tower.LOAD_SEED_BASE + level).integers(0, tower.LOAD_MAX+1, size=ports, dtype=np.int64)
    total = int(loads.sum()); q, r = divmod(total, ports)
    counts = [[q, ports-r]] + ([[q+1, r]] if r else [])
    old = {"V_initial": int(loads @ loads), "balanced_minimum": (ports-r)*q*q+r*(q+1)**2,
           "expected_hash": tower.sha256_of({"canonicalizer":"component_multiset", "components":[[ports,counts]]}),
           "mean_minimum": Fraction(total*total, ports), "total": total,
           "loads_int8_sha256": hashlib.sha256(loads.astype("<i1").tobytes()).hexdigest()}
    assert tower.expected(level, ports, chunk_ports=chunk) == old


def test_bad_chunk_cannot_silently_skip_loads():
    with pytest.raises(ValueError): tower.expected(0, 240, chunk_ports=0)
