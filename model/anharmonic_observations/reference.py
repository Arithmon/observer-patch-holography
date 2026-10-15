"""Independent finite-matrix Lindblad reference, never a source prediction.

This module imports no observer compiler or source operations. It constructs
the declared finite operators directly and integrates the unsplit master
equation during each complete drive pulse. Its purpose is to check the native
unitary execution and the separately supplied diffusion splitting. Agreement
is a numerical finite-dimensional check, not an occupation or apparatus bound.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm


def operators(dimension):
    """Projected oscillator operators with spin as the outer tensor factor."""
    if type(dimension) is not int or dimension < 2:
        raise ValueError('occupation dimension must be an integer at least two')
    annihilation = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    number = np.diag(np.arange(dimension, dtype=float))
    spin_x = np.array([[0., 1.], [1., 0.]])
    return annihilation, number, spin_x


def preparation(dimension, alpha_real, thermal):
    """Direct finite displacement of the declared capped geometric mixture."""
    a, _, _ = operators(dimension)
    if not math.isfinite(alpha_real) or not math.isfinite(thermal) or thermal < 0:
        raise ValueError('finite displacement and nonnegative thermal occupation required')
    ratio = thermal / (1 + thermal)
    weights = (1-ratio)*ratio**np.arange(dimension, dtype=float)
    weights[-1] = ratio**(dimension-1)
    displacement = expm(alpha_real*(a.T-a))
    field = (displacement*weights) @ displacement.conj().T
    joint = np.zeros((2*dimension, 2*dimension), complex)
    joint[:dimension, :dimension] = field
    return joint


def lindblad_pulse(rho, *, delta, omega, duration, sign, dephasing_rate,
                   rtol=2e-11, atol=2e-13):
    """Integrate a complete constant forced-oscillator Lindblad generator.

    The dissipator is explicitly L rho L† - {L†L,rho}/2, with
    L=sqrt(dephasing_rate)*(I_spin tensor N). No split diffusion map is used.
    """
    rho = np.asarray(rho, dtype=complex)
    if rho.ndim != 2 or rho.shape[0] != rho.shape[1] or rho.shape[0] % 2:
        raise ValueError('joint spin-field density must be square with even dimension')
    dimension = rho.shape[0]//2
    a, number, spin_x = operators(dimension)
    if (any(not math.isfinite(value) for value in
            (delta, omega, duration, dephasing_rate)) or duration < 0
            or dephasing_rate < 0 or type(sign) is not int or sign not in (-1, 1)):
        raise ValueError('finite pulse parameters, nonnegative duration/rate and sign +/-1 required')
    if duration == 0:
        return rho.copy()
    joint_number = np.kron(np.eye(2), number)
    hamiltonian = delta*joint_number-sign*omega/2*np.kron(spin_x, a+a.T)
    collapse = math.sqrt(dephasing_rate)*joint_number
    loss = collapse.conj().T @ collapse

    def rhs(_time, vector):
        state = vector.reshape(rho.shape)
        commutator = -1j*(hamiltonian @ state-state @ hamiltonian)
        dissipator = collapse @ state @ collapse.conj().T-(loss @ state+state @ loss)/2
        return (commutator+dissipator).ravel()

    result = solve_ivp(rhs, (0., duration), rho.ravel(), method='DOP853',
                       rtol=rtol, atol=atol)
    if not result.success:
        raise ArithmeticError(f'Independent Lindblad integration failed: {result.message}')
    # Deliberately do not symmetrize, normalize, or clip integration residuals.
    return result.y[:, -1].reshape(rho.shape)


def step(rho, *, delta=2*math.pi*500, omega=math.pi/150e-6,
         pulse_duration=50e-6, theta=.8, phi=0., dephasing_rate=18.):
    """One complete pulse step with no spin measurement or reset."""
    dimension = len(rho)//2
    _, _, spin_x = operators(dimension)
    spin_z = np.diag([1., -1.])
    rotation_z = expm(-.5j*theta*spin_z)
    rotation_x = expm(-.5j*phi*spin_x)
    first = np.kron(rotation_x @ rotation_z @ rotation_x.conj().T, np.eye(dimension))
    second = np.kron(rotation_x.conj().T @ rotation_z @ rotation_x, np.eye(dimension))
    common = dict(delta=delta, omega=omega, dephasing_rate=dephasing_rate)
    rho = lindblad_pulse(rho, duration=pulse_duration, sign=-1, **common)
    rho = first @ rho @ first.conj().T
    rho = lindblad_pulse(rho, duration=2*pulse_duration, sign=1, **common)
    rho = second @ rho @ second.conj().T
    return lindblad_pulse(rho, duration=pulse_duration, sign=-1, **common)


def readout(rho, settings):
    """Direct characteristic-function traces and both final spin masses."""
    dimension = len(rho)//2
    a, _, _ = operators(dimension)
    blocks = [rho[index*dimension:(index+1)*dimension,
                  index*dimension:(index+1)*dimension] for index in (0, 1)]
    masses = np.array([np.trace(block).real for block in blocks])
    if masses[0] <= 0:
        raise ArithmeticError('zero selected spin probability')
    reads = []
    for setting in settings:
        beta = complex(setting['u'], setting['v'])
        displacement = expm(beta*a.T-beta.conjugate()*a)
        signals = np.array([np.real(np.exp(1j*setting['theta'])*np.trace(block @ displacement))
                            for block in blocks])
        joint = np.column_stack(((masses+signals)/2, (masses-signals)/2))
        conditional = joint[0]/masses[0]
        reads.append({**setting, 'spin_detector_joint_probabilities': joint.tolist(),
                      'selected_detector_probabilities': conditional.tolist(),
                      'selected_mean': float(conditional[0]-conditional[1])})
    q = (a+a.T)/math.sqrt(2)
    return {'spin_probabilities': masses.tolist(),
            'total_probability_mass': float(np.trace(rho).real),
            'selected_position_mean': float(np.trace(q @ blocks[0]).real/masses[0]),
            'reads': reads}


def run(dimension, alpha_real, thermal, step_counts, settings, **parameters):
    """Independent direct master-equation records at specified endpoints."""
    if (not step_counts or any(type(value) is not int or value < 0 for value in step_counts)
            or list(step_counts) != sorted(set(step_counts))):
        raise ValueError('sorted unique nonnegative step counts required')
    rho = preparation(dimension, alpha_real, thermal)
    records = []
    for count in range(step_counts[-1]+1):
        if count in step_counts:
            record = readout(rho, settings)
            record['steps'] = count
            records.append(record)
        if count < step_counts[-1]:
            rho = step(rho, **parameters)
    return records
