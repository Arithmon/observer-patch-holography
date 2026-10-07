"""Partition-function controls for the population boundary in #1048."""

import mpmath
import numpy as np
import pytest

from quantum_information.entropy_response import gibbs_entropy_response


@pytest.mark.parametrize("gap", [700., 708., 710., 720., 744., 745., 746.])
@pytest.mark.parametrize("scale", [1., 1e150])
@pytest.mark.parametrize("sign", [-1., 1.])
def test_response_uses_resolved_populations_or_refuses(gap, scale, sign):
    beta = sign*gap/scale
    ctx = mpmath.mp.clone()
    ctx.dps = 500
    # Independent scalar partition function, exact binary64 inputs. No
    # thermal constructor, matrix eigensolve or derivative ratio is reused.
    e, b = ctx.mpf(scale), ctx.mpf(beta)
    weight = ctx.exp(-abs(b*e))
    rare = weight/(1+weight)
    variance = e*e*rare*(1-rare)
    ds = -b*variance
    if rare < np.finfo(float).tiny:
        # Conservative domain: every subnormal population is refused even
        # if its scaled derivatives would be normal, e.g. gap=745, E=1e150.
        with pytest.raises(ValueError, match="population|underflow|precision"):
            gibbs_entropy_response(np.diag([0., scale]), beta)
    else:
        result = gibbs_entropy_response(np.diag([0., scale]), beta)
        assert result["energy_variance"] == pytest.approx(float(variance), rel=5e-13, abs=0)
        assert result["dt_dlambda"] == pytest.approx(float(-variance), rel=5e-13, abs=0)
        assert result["ds_dlambda"] == pytest.approx(float(ds), rel=5e-13, abs=0)


def test_maintainer_rescaled_example_has_normal_outputs():
    ctx = mpmath.mp.clone()
    ctx.dps = 500
    e, b = ctx.mpf(1e150), ctx.mpf(7.45e-148)
    p = 1/(1+ctx.exp(b*e))
    variance = e*e*p*(1-p)
    assert float(variance) == pytest.approx(2.822350730472012e-24, rel=2e-15, abs=0)
    assert float(-b*variance) == pytest.approx(-2.102651294201649e-171, rel=2e-15, abs=0)
    with pytest.raises(ValueError, match="population|precision"):
        gibbs_entropy_response(np.diag([0., 1e150]), 7.45e-148)
