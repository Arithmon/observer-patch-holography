"""Representable source spans must survive avoidable chart underflow."""

import numpy as np
import pytest

from quantum_information.algebras import (
    orthonormal_operator_basis, separation_modulus,
)


@pytest.mark.parametrize("size", (2, 4))
@pytest.mark.parametrize("phase", (1., 1j))
@pytest.mark.parametrize("power", (537, 538, 600))
@pytest.mark.parametrize("tiny_leads", (False, True))
def test_exact_shear_can_use_a_representable_chart(size, phase, power, tiny_leads):
    epsilon = 2.**-power
    u = np.zeros((size, size), complex)
    v = np.zeros_like(u)
    tiny, first, second = (0, 1, size) if tiny_leads else (size, 0, 1)
    u.flat[first], u.flat[tiny] = 1., epsilon
    v.flat[second], v.flat[tiny] = phase, epsilon
    delta = 2.**-10
    assert np.array_equal((u+delta*v)-u, delta*v)

    # U-V and U+V are exactly HS-orthogonal because ||U||=||V||.
    # Their entries 0, +/-1, +/-phase, 2epsilon give a representable chart. In
    # contrast, Gram-Schmidt on U,V creates an avoidable epsilon**2 tail.
    presentations = ([u, v], [v, u], [phase*u, u+delta*v], [u-v, u+v])
    omega = np.eye(size)[1]
    # Source coordinates (first,second,tiny) span the plane orthogonal to
    # (-epsilon,-phase*epsilon,1). Thus P_(tiny,first)=epsilon/(1+2epsilon**2),
    # which rounds to epsilon at these exponents. Other entries tend to
    # diag(1,1,0); their omitted epsilon**2 terms are below this normwise
    # allowance, so check the representable first-order component separately.
    leading = np.zeros((size*size, size*size), complex)
    leading[first, first] = leading[second, second] = 1.
    for basis in presentations:
        normalized = orthonormal_operator_basis(basis)
        columns = np.column_stack([b.reshape(-1) for b in normalized])
        projector = columns@columns.conj().T
        np.testing.assert_allclose(projector, leading, atol=3e-14, rtol=0)
        assert projector[tiny, first] == pytest.approx(epsilon, rel=5e-14, abs=0)
        assert projector[tiny, second] == pytest.approx(np.conj(phase)*epsilon, rel=5e-14, abs=0)
        restored = (projector@u.reshape(-1)).reshape(u.shape)
        assert restored.flat[tiny] == pytest.approx(epsilon, rel=5e-14, abs=0)
        # U (or V when tiny leads) is a nonzero source annihilator of e1.
        assert separation_modulus(basis, omega) == 0.
