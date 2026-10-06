"""A positive binary64 Gibbs mass need not have resolved relative precision."""

import numpy as np
import pytest

from quantum_information import gibbs_sectors
from maxent.information_projection import gibbs_state


@pytest.mark.parametrize("gap", [0., 20., 700., 708.])
def test_resolved_populations_are_still_constructed(gap):
    state, _ = gibbs_state([np.diag([0., gap])], [1.])
    expected = np.exp(-gap)/(1+np.exp(-gap))
    assert state[1, 1].real == pytest.approx(expected, rel=2e-14, abs=0)
    sectors = gibbs_sectors([np.zeros((1, 1))]*2, [0., gap])
    assert sectors[1][0] == pytest.approx(expected, rel=2e-14, abs=0)


@pytest.mark.parametrize("gap", [710., 720., 744., 745.])
@pytest.mark.parametrize("kind", ["state", "conditional", "sector", "joint"])
def test_state_producers_refuse_subnormal_population_precision(gap, kind):
    with pytest.raises(ValueError, match="population|precision|underflow"):
        if kind == "state":
            gibbs_state([np.diag([0., gap])], [1.])
        elif kind == "conditional":
            gibbs_sectors([np.diag([0., gap])], [0.])
        elif kind == "sector":
            gibbs_sectors([np.zeros((1, 1))]*2, [0., gap])
        else:
            # Each factor is normal; only the full-state product is subnormal.
            gibbs_sectors([np.zeros((1, 1)), np.diag([0., gap/2])], [0., gap/2])
