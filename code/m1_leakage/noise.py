"""Driven full-M6 channels and a finite common bath, without postselection."""

import numpy as np
from scipy.linalg import expm
from .format import pack
from .interfaces import reduction


KINDS = ('relax', 'leak', 'driven_return', 'pair_leak', 'coherent_leak', 'thermal_return')


def unit(i, j, d=6):
    a = np.zeros((d, d), complex)
    a[i, j] = 1
    return a


def specification(kind):
    x, y = unit(0, 1)+unit(1, 0), -1j*unit(0, 1)+1j*unit(1, 0)
    h, b, jumps = np.zeros((6, 6), complex), np.zeros((6, 6), complex), []
    d = 2
    if kind == 'relax':
        duration, transitions = .7, [(0, 1, .4)]
    elif kind == 'leak':
        duration, transitions = .9, [(4, 1, .3)]
    elif kind == 'driven_return':
        h = .7*x+.2*(unit(0, 0)-unit(1, 1))
        duration, transitions = .6, [(3, 1, .35), (0, 3, .11)]
    elif kind == 'pair_leak':
        d, duration = 4, .8
        h = .8*np.diag([1, 1, 1, 0, 0, 0])
        transitions = [(5, 3, .4), (1, 5, .09)]
    elif kind == 'coherent_leak':
        h, b = .6*x, .25*(unit(1, 4)+unit(4, 1))+.17*unit(4, 4)
        duration, transitions = .65, []
    elif kind == 'thermal_return':
        h, duration = .4*y, .7
        transitions = [(0, 1, .2), (1, 0, .1), (4, 1, .17), (1, 4, .07)]
    else:
        raise ValueError('noise catalog')
    jumps = [np.sqrt(rate)*unit(a, b) for a, b, rate in transitions]
    return d, duration, h, b, jumps


def generator(h, jumps):
    d = len(h)
    identity = np.eye(d)
    value = -1j*(np.kron(identity, h)-np.kron(h.T, identity))
    for a in jumps:
        aa = a.conj().T@a
        value += np.kron(a.conj(), a)-.5*(np.kron(identity, aa)+np.kron(aa.T, identity))
    return value


def choi(channel, input_size, output_size=None):
    output_size = output_size or input_size
    result = np.zeros((input_size*output_size,)*2, complex)
    for i in range(input_size):
        for j in range(input_size):
            result[i*output_size:(i+1)*output_size, j*output_size:(j+1)*output_size] = (
                channel[:, i+input_size*j].reshape(output_size, output_size, order='F')/input_size)
    return result


def flagged(channel, d):
    groups = reduction(d)
    result = np.zeros(((2*d)**2, d*d), complex)
    for i in range(d):
        for j in range(d):
            state = channel[:, i+6*j].reshape(6, 6, order='F')
            out = np.zeros((2*d, 2*d), complex)
            for flag, ks in enumerate(groups):
                out[flag*d:(flag+1)*d, flag*d:(flag+1)*d] = sum(k@state@k.conj().T for k in ks)
            result[:, i+d*j] = out.reshape(-1, order='F')
    return choi(result, d, 2*d)


def channels():
    rows = []
    for kind in KINDS:
        d, t, h, b, jumps = specification(kind)
        actual = expm(t*generator(h+b, jumps))
        ideal = expm(t*generator(h, []))
        delta = choi(actual-ideal, 6)
        gamma = 2*np.linalg.norm(b, ord=2)+2*sum(np.linalg.norm(a, ord=2)**2 for a in jumps)
        rows.append(dict(kind=kind, code_dimension=d, duration=t,
                         channel_choi=pack(choi(actual, 6)), flagged_choi=pack(flagged(actual, d)),
                         diamond_upper=float(t*gamma),
                         bell_half_trace_distance=float(sum(abs(np.linalg.eigvalsh(delta)))/2)))
    return rows


def reduced_bath_channel(u):
    ks = [u[b::2, 0::2] for b in range(2)]
    return sum(np.kron(k.conj(), k) for k in ks)


def bath():
    x, z = np.array([[0, 1], [1, 0]]), np.diag([1., -1.])
    h = .4*np.kron(unit(0, 1)+unit(1, 0), np.eye(2))+np.kron(np.eye(6), .17*x)
    interaction = .27*np.kron(unit(1, 4)+unit(4, 1), z)+.11*np.kron(unit(0, 0)-unit(1, 1), x)
    u = expm(-.6j*(h+interaction))
    shared = reduced_bath_channel(u@u)
    first = reduced_bath_channel(u)
    reset = first@first
    difference = choi(shared-reset, 6)
    code_difference = flagged(shared, 2)-flagged(reset, 2)
    return dict(duration=.6, shared_choi=pack(choi(shared, 6)),
                reset_choi=pack(choi(reset, 6)),
                bell_half_trace_difference=float(sum(abs(np.linalg.eigvalsh(difference)))/2),
                code_flagged_half_trace_difference=float(sum(abs(np.linalg.eigvalsh(code_difference)))/2))


def obstruction():
    rows = []
    for gamma in (.2, .5, .8):
        k0 = np.diag([1., np.sqrt(1-gamma)])
        k1 = np.array([[0., np.sqrt(gamma)], [0., 0.]])
        j = sum(np.outer(k.T.reshape(-1), k.T.reshape(-1).conj()) for k in (k0, k1))
        identity = np.array([1., 0., 0., 1.])
        witness = np.array([np.sqrt(1-gamma), 0., 0., -1.])
        p = .25
        rows.append(dict(damping=gamma, subtracted_identity_weight=p,
                         witness=pack(witness.reshape(-1, 1)),
                         quadratic_form=float(witness@(j-p*np.outer(identity, identity))@witness)))
    return rows


def calibration():
    theta, error = .7, .04
    z = np.diag([1., -1.])
    plus = np.ones(2)/np.sqrt(2)
    rows = []
    for omega in (1, 2, 4, 16):
        t = theta/omega
        a = expm(-1j*t*omega*z/2)@plus
        b = expm(-1j*t*omega*(1+error)*z/2)@plus
        distance = sum(abs(np.linalg.eigvalsh(np.outer(a, a.conj())-np.outer(b, b.conj()))))/2
        rows.append(dict(drive=omega, duration=t, relative_error=error, half_trace_distance=float(distance)))
    return rows


def candidate():
    return dict(channels=channels(), bath=bath(), ideal_mixture_obstruction=obstruction(), calibration=calibration())
