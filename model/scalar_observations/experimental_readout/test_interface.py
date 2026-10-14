"""Connect the measured characteristic convention to the actual source readout."""
import math

import numpy as np
import pytest
from scipy.linalg import expm

from ..source import instrument


@pytest.mark.parametrize('beta', [.4j, .3+.4j, -.6+.2j])
def test_compiled_binary_readout_matches_both_characteristic_components(beta):
    dimension = 8
    a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    # An independently formed mixed state, with coherences and a high-level
    # component, prevents a Gaussian-only or blank-only identity from passing.
    v = np.arange(1, dimension+1)+1j*np.arange(dimension)[::-1]
    v /= np.linalg.norm(v)
    rho = .7*np.outer(v, v.conj())+.3*np.eye(dimension)/dimension
    characteristic = np.trace(rho@expm(beta*a.T-beta.conjugate()*a))
    phi = np.angle(beta)-math.pi/2
    rotation = np.diag(np.exp(1j*phi*np.arange(dimension)))
    rotated = rotation.conj().T@rho@rotation
    for theta, expected in ((0., characteristic.real), (-math.pi/2, characteristic.imag)):
        kraus, metadata = instrument(dimension, eta=2*abs(beta), theta=theta)
        probabilities = [np.trace(k@rotated@k.conj().T).real for k in kraus]
        assert abs(sum(probabilities)-1) < 2e-12
        assert abs(probabilities[0]-probabilities[1]-expected) < 2e-12
        assert metadata['complete_pointer_outcomes_retained']
    # This compares finite functions of the compressed coordinate; it is not
    # an occupation-error certificate for the experimental state.
