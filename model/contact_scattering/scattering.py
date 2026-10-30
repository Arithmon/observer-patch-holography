"""Rank-one two-fermion response and fixed-broadening volume control.

The local action is exactly H=dGamma(h)+U sum_x n_x,up n_x,down,
h(k)=-c sigma.sin(k). At a fixed total momentum its contact is rank one.
A resolvent or T matrix is a response operator, not a quantum channel.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
import math

import numpy as np

from ..bound_pair.fiber import SINGLET, SWAP
from ..interacting_field.lattice import SIGMA, real


def _side(side, maximum=256):
    if type(side) is not int or not 4 <= side <= maximum:
        raise ValueError(f'integer side 4..{maximum} required')
    return side


def _momentum(momentum):
    if np.ma.isMaskedArray(momentum):
        raise ValueError('masked momentum forbidden')
    value = np.asarray(momentum)
    if value.shape != (3,) or value.dtype.kind not in 'fiu' or not np.isfinite(value).all():
        raise ValueError('three finite real total momenta required')
    return value.astype(float)


def _frequency(frequency):
    if isinstance(frequency, (bool, np.bool_)) or not isinstance(frequency, (int, float, complex)):
        raise ValueError('finite upper-half-plane frequency required')
    z = complex(frequency)
    if not math.isfinite(z.real) or not math.isfinite(z.imag) or z.imag <= 0:
        raise ValueError('finite upper-half-plane frequency required')
    return z


def _finite(value):
    if not np.isfinite(value).all():
        raise ArithmeticError('nonfinite scattering response')
    return value


def free_blocks(k, momentum, speed):
    """All four spin states of the existing two-body symbol, at every k."""
    h1 = -speed*np.einsum('ka,aij->kij', np.sin(k), np.array(SIGMA))
    h2 = -speed*np.einsum('ka,aij->kij', np.sin(momentum-k), np.array(SIGMA))
    return (np.einsum('kab,cd->kacbd', h1, np.eye(2))+
            np.einsum('ab,kcd->kacbd', np.eye(2), h2)).reshape(-1, 4, 4)


def free_green(momentum, frequency, *, side=16, speed=.125):
    """Complete singlet resolvent quadrature; no selected momentum pocket.

    Stream complete x-planes so the host need not allocate a global fiber.
    The mathematical finite-torus interpretation requires reciprocal-grid K;
    a general K instead denotes the explicitly declared periodic quadrature.
    """
    side = _side(side)
    total, z = _momentum(momentum), _frequency(frequency)
    c = real(speed, 'coordinate speed', nonnegative=True)
    yz = np.array(list(product(range(side), repeat=2)), dtype=float)*2*np.pi/side
    answer = 0j
    for x in range(side):
        k = np.column_stack((np.full(len(yz), 2*np.pi*x/side), yz))
        blocks = free_blocks(k, total, c)
        vectors = np.linalg.solve(z*np.eye(4)-blocks,
                                  np.broadcast_to(SINGLET, (len(k), 4))[..., None])[..., 0]
        answer += np.einsum('i,ki->', SINGLET.conj(), vectors)
    return complex(_finite(answer/side**3))


def response(momentum, frequency, *, side=16, speed=.125, interaction=1.):
    """Exact scalar response formula, evaluated without an on-shell limit."""
    total, z = _momentum(momentum), _frequency(frequency)
    c = real(speed, 'coordinate speed', nonnegative=True)
    u = real(interaction, 'contact repulsion', nonnegative=True)
    g0 = free_green(total, z, side=side, speed=c)
    denominator = 1-u*g0
    tau = u/denominator
    interacting = g0/denominator
    optical = -tau.imag-abs(tau)**2*(-g0.imag)
    if -g0.imag <= 0 or -interacting.imag <= 0:
        raise ArithmeticError('upper-half-plane spectral positivity failed')
    _finite(np.array([tau, interacting, optical]))
    return {'free_green':g0, 'tau':tau, 'contact_green':interacting,
            'optical_residual':float(abs(optical)),
            'free_lorentz_density':float(-g0.imag/np.pi),
            'contact_lorentz_density':float(-interacting.imag/np.pi),
            'finite_torus_momentum':bool(np.max(abs(total*side/(2*np.pi)-
                                                   np.rint(total*side/(2*np.pi)))) < 1e-12)}


def _rational(value, name, *, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and value == 0):
        raise ValueError(f'exact {"positive" if positive else "nonnegative"} rational {name} required')
    return F(value)


def volume_enclosure(side, eta, *, speed=F(1, 8), interaction=F(1)):
    """Rational exponential quadrature bound at FIXED eta>0, uniform in E,K.

    Complex strip width a=log(1+eta/(12c)) bounds the free resolvent by2/eta.
    The Fourier alias sum is [(1+q)/(1-q)]^3-1 with q=exp(-aL), computed
    rationally. This is a mathematical truncation bound; floating errors
    of the evaluated quadrature are a separate diagnostic.
    """
    side = _side(side, 4096)
    eta = _rational(eta, 'broadening', positive=True)
    c, u = _rational(speed, 'speed'), _rational(interaction, 'interaction')
    q = (12*c/(12*c+eta))**side
    green_error = 2/eta*(((1+q)/(1-q))**3-1)
    # tau=U+U²<d,(z-H)^-1 d>, hence |tau|<=U+U²/eta.
    tau_bound = u+u*u/eta
    return {'side':side, 'eta':eta, 'alias_ratio':q,
            'free_green_error':green_error,
            'tau_error':tau_bound*tau_bound*green_error,
            'contact_green_error':(1+u/eta)**2*green_error,
            'tau_absolute_bound':tau_bound,
            'uniform_in_real_frequency_and_total_momentum':True,
            'on_shell_limit_claimed':False}


@dataclass(frozen=True)
class ContactFiber:
    """Complete four-spin fiber and its physical exchange subspace."""
    momentum: tuple[float, float, float] = (0., 0., 0.)
    side: int = 4
    speed: float = .125
    interaction: float = 1.

    def __post_init__(self):
        _side(self.side, 24)
        object.__setattr__(self, 'momentum', tuple(_momentum(self.momentum)))
        object.__setattr__(self, 'speed', real(self.speed, 'coordinate speed', nonnegative=True))
        object.__setattr__(self, 'interaction', real(self.interaction, 'repulsion', nonnegative=True))

    @property
    def volume(self):
        return self.side**3

    @property
    def dimension(self):
        return 4*self.volume

    @property
    def k(self):
        return np.array(list(product(range(self.side), repeat=3)), dtype=float)*2*np.pi/self.side

    @property
    def d(self):
        return np.tile(SINGLET, self.volume)/math.sqrt(self.volume)

    @property
    def blocks(self):
        return free_blocks(self.k, np.array(self.momentum), self.speed)

    def state(self, value):
        if np.ma.isMaskedArray(value):
            raise ValueError('masked state forbidden')
        x = np.asarray(value)
        if (x.ndim not in (1, 2) or x.shape[0] != self.dimension or
                x.dtype.kind not in 'fciu' or not np.isfinite(x).all() or
                (x.ndim == 2 and x.shape[1] == 0) or x.size > 2**20):
            raise ValueError('finite bounded complete fiber/reference columns required')
        return x.astype(complex)

    def free_resolvent(self, state, frequency):
        x, z = self.state(state), _frequency(frequency)
        shape = x.shape
        rhs = x.reshape(self.volume, 4, -1)
        return _finite(np.linalg.solve(z*np.eye(4)-self.blocks, rhs).reshape(shape))

    def resolvent(self, state, frequency):
        x, z = self.state(state), _frequency(frequency)
        r0x, r0d = self.free_resolvent(x, z), self.free_resolvent(self.d, z)
        g0 = np.vdot(self.d, r0d)
        tau = self.interaction/(1-self.interaction*g0)
        coefficient = self.d.conj()@r0x
        correction = r0d*coefficient if x.ndim == 1 else r0d[:, None]*coefficient
        return _finite(r0x+tau*correction)

    def t_apply(self, state, frequency):
        x = self.state(state)
        tau = response(self.momentum, frequency, side=self.side,
                       speed=self.speed, interaction=self.interaction)['tau']
        coefficient = self.d.conj()@x
        return self.d*tau*coefficient if x.ndim == 1 else self.d[:, None]*tau*coefficient

    def physical_project(self, state):
        """Project onto psi(K-k)=-SWAP psi(k), retaining every spin channel."""
        total = np.array(self.momentum)*self.side/(2*np.pi)
        if np.max(abs(total-np.rint(total))) >= 1e-12:
            raise ValueError('physical finite exchange requires reciprocal-grid total momentum')
        indices = np.array(list(product(range(self.side), repeat=3)))
        partner = (np.rint(total).astype(int)-indices)%self.side
        partner = (partner[:,0]*self.side+partner[:,1])*self.side+partner[:,2]
        x = self.state(state)
        swapped = np.einsum('ij,kj...->ki...', SWAP, x.reshape(self.volume, 4, -1)[partner])
        return (x-swapped.reshape(x.shape))/2
