"""Total-momentum fibers of the SAME action as interacting_field.

The periodic boxes below change the open boundary condition, not the local
couplings.  A continuous total momentum denotes a Brillouin-zone quadrature;
only reciprocal-grid K labels a fiber of the stated finite untwisted torus.
First-quantized vectors retain all four spin components at every relative
momentum, with the physical fermion exchange constraint checked explicitly.
"""
from __future__ import annotations
from itertools import product
import math
import numpy as np
from scipy.optimize import brentq
from ..interacting_field.lattice import SIGMA, real

SINGLET = np.array([0, 1, -1, 0], complex)/math.sqrt(2)
SWAP = np.eye(4)[[0, 2, 1, 3]]


def grid(side):
    if type(side) is not int or not 4 <= side <= 96:
        raise ValueError('integer side 4..96 required')
    return np.array(list(product(range(side), repeat=3)), float)*2*np.pi/side


def parameters(speed, interaction):
    c = real(speed, 'coordinate speed', nonnegative=True)
    u = real(interaction, 'repulsion', nonnegative=True)
    if not u > 2*math.sqrt(3)*c:
        raise ValueError('isolated-pair regime U > 2 sqrt(3) c required')
    return c, u


class PairFiber:
    def __init__(self, momentum, *, side=12, speed=.125, interaction=1.):
        self.speed, self.interaction = parameters(speed, interaction)
        self.k = grid(side)
        self.side, self.volume = side, side**3
        self.momentum = np.asarray(momentum, dtype=float)
        if self.momentum.shape != (3,) or not np.isfinite(self.momentum).all():
            raise ValueError('three finite coordinate momenta required')
        self.dimension = 4*self.volume
        self.d = np.tile(SINGLET, self.volume)/math.sqrt(self.volume)
        self.sin1, self.sin2 = np.sin(self.k), np.sin(self.momentum-self.k)
        h1 = -self.speed*np.einsum('ka,aij->kij', self.sin1, np.array(SIGMA))
        h2 = -self.speed*np.einsum('ka,aij->kij', self.sin2, np.array(SIGMA))
        self.blocks = (np.einsum('kab,cd->kacbd', h1, np.eye(2))+
                       np.einsum('ab,kcd->kacbd', np.eye(2), h2)).reshape(-1,4,4)
        n1, n2 = np.linalg.norm(self.sin1,axis=1), np.linalg.norm(self.sin2,axis=1)
        dot = np.sum(self.sin1*self.sin2,axis=1)
        cosine = np.divide(dot,n1*n2,out=np.zeros_like(dot),where=n1*n2>1e-24)
        # Spectral weights of the spin singlet, including degenerate points.
        self.free_energies = np.array([self.speed*(s*n1+t*n2)
                                      for s,t in product((-1,1),repeat=2)]).T
        self.weights = np.array([(1-s*t*cosine)/4
                                 for s,t in product((-1,1),repeat=2)]).T
        self.weights = np.maximum(self.weights, 0.)

    @property
    def is_finite_torus_momentum(self):
        x = self.momentum*self.side/(2*np.pi)
        return bool(np.max(np.abs(x-np.rint(x))) < 1e-12)

    def free_apply(self, state):
        a = np.asarray(state,complex)
        if a.shape != (self.dimension,) or not np.isfinite(a).all():
            raise ValueError('one finite complete fiber state required')
        return np.einsum('kab,kb->ka',self.blocks,a.reshape(-1,4)).reshape(-1)

    def apply(self, state):
        return self.free_apply(state)+self.interaction*self.d*np.vdot(self.d,state)

    def exchange_residual(self, state):
        if not self.is_finite_torus_momentum:
            raise ValueError('exchange permutation requires reciprocal-grid K')
        km = np.rint(self.momentum*self.side/(2*np.pi)).astype(int)
        indices = np.array(list(product(range(self.side),repeat=3)))
        partner = (km-indices)%self.side
        partner = (partner[:,0]*self.side+partner[:,1])*self.side+partner[:,2]
        a = np.asarray(state).reshape(-1,4)
        return float(np.linalg.norm(a+a[partner]@SWAP.T))

    def resolvent_moment(self, energy, power=1):
        if type(power) is not int or power not in (1,2):
            raise ValueError('resolvent power one or two required')
        e = real(energy,'energy')
        if e <= float(np.max(self.free_energies)):
            raise ValueError('energy must exceed the complete free spectrum')
        return float(np.sum(self.weights/(e-self.free_energies)**power)/self.volume)

    def pole(self):
        u = self.interaction
        b = 2*math.sqrt(3)*self.speed
        at_u = u*self.resolvent_moment(u)-1
        energy = u if abs(at_u)<2e-14 else brentq(
            lambda e: u*self.resolvent_moment(e)-1,u,u+b,xtol=5e-15)
        residue = 1/(u*u*self.resolvent_moment(energy,2))
        return {'energy':float(energy),'contact_residue':float(residue),
                'secular_residual':abs(u*self.resolvent_moment(energy)-1),
                'finite_torus_momentum':self.is_finite_torus_momentum}

    def eigenstate(self):
        p = self.pole()
        resolved = np.linalg.solve(p['energy']*np.eye(4)-self.blocks,
                                   np.broadcast_to(SINGLET,(self.volume,4))[...,None])[...,0]
        return (self.interaction*math.sqrt(p['contact_residue']/self.volume)*resolved).reshape(-1)


def clock_subspace(side=4, speed=.125, interaction=1.):
    if type(side) is not int or side%2:
        raise ValueError('even torus side required')
    f = PairFiber([np.pi]*3,side=side,speed=speed,interaction=interaction)
    d = f.d.copy()
    r = (math.sqrt(2)*np.sin(f.k[:,0])[:,None]*d.reshape(-1,4)).reshape(-1)
    basis = np.column_stack((r,d))
    applied = np.column_stack((f.apply(r),f.apply(d)))
    restriction = np.diag([0.,f.interaction])
    return {'basis':basis,'reference':r,'pair':d,'restriction':restriction,
            'momentum':[np.pi]*3,'side':side,'speed':f.speed,'interaction':f.interaction,
            'residuals':{'orthonormal':float(np.linalg.norm(basis.conj().T@basis-np.eye(2))),
                         'invariant':float(np.linalg.norm(applied-basis@restriction)),
                         'reference_exchange':f.exchange_residual(r),
                         'pair_exchange':f.exchange_residual(d)}}
