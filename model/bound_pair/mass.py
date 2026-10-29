"""Exact rational bounds for the infinite-lattice pair curvature at K=pi.

These are kinetic coordinate masses in the supplied lattice action.  No
vacuum, physical clock calibration or SI mass is selected by this code.
"""
from fractions import Fraction as F
from math import comb, factorial
import numpy as np
from .fiber import PairFiber, grid, parameters


def moment(n):
    return F(comb(2*n,n),4**n)


def curvature_ratio_interval(speed=F(1,8), interaction=F(1), terms=16):
    if not isinstance(speed,F) or not isinstance(interaction,F):
        raise ValueError('exact Fraction couplings required')
    if type(terms) is not int or not 0 <= terms <= 64:
        raise ValueError('0..64 geometric expansion terms required')
    if speed <= 0 or interaction <= 0 or 12*speed**2 >= interaction**2:
        raise ValueError('0 < 12 c^2 < U^2 required')
    a=(speed/interaction)**2
    lower=F(1)
    for n in range(terms+1):
        avg=F(0)
        for i in range(n+1):
            for j in range(n-i+1):
                k=n-i-j
                coeff=factorial(n)//(factorial(i)*factorial(j)*factorial(k))
                avg += coeff*(moment(i)-moment(i+1))*moment(j+1)*moment(k)
        lower += 16*a*(4*a)**n*avg
    tail=4*a*(12*a)**(terms+1)/(1-12*a)
    return lower, lower+tail


def curvature_quadrature(side=32,speed=.125,interaction=1.):
    c,u=parameters(speed,interaction)
    k=grid(side)
    sin2=np.sin(k)**2
    total=np.sum(sin2,axis=1)
    return float(2*c*c/u*np.mean(np.cos(k[:,0])**2*(
        1+4*c*c*(total-sin2[:,0])/(u*u-4*c*c*total))))


def certificate():
    c,u=F(1,8),F(1)
    lo,hi=curvature_ratio_interval(c,u)
    curvature=(c*c/u*lo,c*c/u*hi)
    kinetic_mass=(1/curvature[1],1/curvature[0])
    energy_parts=[]
    for label,k in (('eta',[np.pi]*3),('one_zero',[0,np.pi,np.pi]),('origin',[0,0,0])):
        fiber=PairFiber(k,side=8)
        pole=fiber.pole()
        v=fiber.eigenstate()
        interaction=float(fiber.interaction*abs(np.vdot(fiber.d,v))**2)
        kinetic=float(np.vdot(v,fiber.free_apply(v)).real)
        energy_parts.append({'momentum_label':label,'torus_side':8,
            'energy':pole['energy'],'interaction_energy':interaction,
            'kinetic_energy':kinetic,'d_energy_d_U':pole['contact_residue'],
            'split_residual':abs(interaction+kinetic-pole['energy'])})
    return {'coordinate_c':str(c),'coordinate_U':str(u),
            'curvature_ratio_U_kappa_over_c_squared':[str(lo),str(hi)],
            'curvature_interval':[str(x) for x in curvature],
            'kinetic_mass_interval':[str(x) for x in kinetic_mass],
            'relativistic_same_c_mass_value':str(u/(c*c)),
            'elementary_ratio_interval':['17/16','14/13'],
            'ratio_invariant_under_common_clock_and_length_rescaling':True,
            'physical_mass_calibrated':False,
            'same_state_energy_parts':energy_parts,
            'quadrature_curvature_diagnostic':curvature_quadrature()}
