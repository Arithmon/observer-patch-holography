"""Independent matrix-unit differential generators and complete channel checks."""

import math
import numpy as np
from m1_fixed_noise.algebra_check import exponential
from .format import close, integer, keys, need, number, unpack


def matrix_unit(i, j, d):
    result = np.zeros((d, d), complex)
    result[i, j] = 1
    return result


def specification(kind):
    h, b = np.zeros((6, 6), complex), np.zeros((6, 6), complex)
    parameters = {
        'relax': (2, .7, [(0, 1, .4)]),
        'leak': (2, .9, [(4, 1, .3)]),
        'driven_return': (2, .6, [(3, 1, .35), (0, 3, .11)]),
        'pair_leak': (4, .8, [(5, 3, .4), (1, 5, .09)]),
        'coherent_leak': (2, .65, []),
        'thermal_return': (2, .7, [(0, 1, .2), (1, 0, .1), (4, 1, .17), (1, 4, .07)])}
    d, duration, transitions = parameters[kind]
    if kind == 'driven_return':
        h[0, 1] = h[1, 0] = .7
        h[0, 0], h[1, 1] = .2, -.2
    if kind == 'pair_leak':
        h[0, 0] = h[1, 1] = h[2, 2] = .8
    if kind == 'coherent_leak':
        h[0, 1] = h[1, 0] = .6
        b[1, 4] = b[4, 1] = .25
        b[4, 4] = .17
    if kind == 'thermal_return':
        h[0, 1], h[1, 0] = -.4j, .4j
    jumps = []
    for target, source, rate in transitions:
        jump = np.zeros((6, 6), complex)
        jump[target, source] = math.sqrt(rate)
        jumps.append(jump)
    return d, duration, h, b, jumps


def evolution(h, jumps, t):
    d = len(h)
    generator = np.zeros((d*d, d*d), complex)
    for i in range(d):
        for j in range(d):
            e = matrix_unit(i, j, d)
            derivative = -1j*(h@e-e@h)
            for k in jumps:
                product = k.conj().T@k
                derivative += k@e@k.conj().T-(product@e+e@product)/2
            generator[:, i+d*j] = derivative.reshape(-1, order='F')
    return exponential(t*generator)


def choi(channel, d, out=None):
    out = out or d
    result = np.zeros((d*out, d*out), complex)
    for i in range(d):
        for j in range(d):
            value = channel[:, i+d*j].reshape((out, out), order='F')
            result[i*out:(i+1)*out, j*out:(j+1)*out] = value/d
    return result


def cptp(j, d, out):
    close(j, j.conj().T, 'Hermitian Choi matrix')
    need(np.linalg.eigvalsh(j).min() > -3e-10, 'complete positivity')
    partial = np.einsum('iaja->ij', j.reshape(d, out, d, out))
    close(partial, np.eye(d)/d, 'trace preservation on every input')


def flagged_oracle(channel, d):
    result = np.zeros((2*d*d,)*2, complex)
    for i in range(d):
        for j in range(d):
            rho = channel[:, i+6*j].reshape((6, 6), order='F')
            block = np.zeros((2*d, 2*d), complex)
            block[:d, :d] = rho[:d, :d]
            block[d, d] = sum(rho[k, k] for k in range(d, 6))
            result[i*2*d:(i+1)*2*d, j*2*d:(j+1)*2*d] = block/d
    return result


def verify_channels(rows):
    kinds = ('relax', 'leak', 'driven_return', 'pair_leak', 'coherent_leak', 'thermal_return')
    need(type(rows) is list and len(rows) == len(kinds), 'complete full-interface noise catalog')
    for row, kind in zip(rows, kinds):
        keys(row, 'kind code_dimension duration channel_choi flagged_choi diamond_upper bell_half_trace_distance')
        need(row['kind'] == kind, 'noise case order')
        d, t, h, b, jumps = specification(kind)
        integer(row['code_dimension'], d, d)
        close(number(row['duration']), t, 'frozen duration')
        actual = evolution(h+b, jumps, t)
        ideal = evolution(h, [], t)
        j = unpack(row['channel_choi'], (36, 36))
        close(j, choi(actual, 6), 'full driven channel')
        cptp(j, 6, 6)
        flags = unpack(row['flagged_choi'], (2*d*d, 2*d*d))
        close(flags, flagged_oracle(actual, d), 'unconditional flagged readout')
        cptp(flags, d, 2*d)
        upper = t*(2*np.linalg.norm(b, 2)+2*sum(np.linalg.norm(a, 2)**2 for a in jumps))
        distance = sum(abs(np.linalg.eigvalsh(choi(actual-ideal, 6))))/2
        close(number(row['diamond_upper']), upper, 'analytic diamond upper bound')
        close(number(row['bell_half_trace_distance']), distance, 'Bell-state trace distance')
        need(distance > 1e-4 and distance <= upper/2+3e-10, 'nonzero noise within analytic majorant')
        failure = sum(flags[i*2*d+d, i*2*d+d].real for i in range(d))
        if kind != 'relax':
            need(failure > 1e-4, 'leakage must actually occur')
        else:
            close(failure, 0., 'ordinary damping control')


def bath_channel(u):
    # Evolve every full-M6 matrix unit with bath |0><0|; trace that same bath.
    result = np.zeros((36, 36), complex)
    for i in range(6):
        for j in range(6):
            joint = np.outer(u[:, 2*i], u[:, 2*j].conj()).reshape(6, 2, 6, 2)
            out = joint[:, 0, :, 0]+joint[:, 1, :, 1]
            result[:, i+6*j] = out.reshape(-1, order='F')
    return result


def verify_bath(row):
    keys(row, 'duration shared_choi reset_choi bell_half_trace_difference code_flagged_half_trace_difference')
    close(number(row['duration']), .6, 'bath interaction duration')
    h = np.zeros((12, 12), complex)
    # Build the joint Hamiltonian entry by entry, without the producer's
    # Kronecker construction or its reduced-channel Kraus decomposition.
    for bath in (0, 1):
        h[bath, 2+bath] += .4
        h[2+bath, bath] += .4
        h[2+bath, 8+bath] += .27*(-1)**bath
        h[8+bath, 2+bath] += .27*(-1)**bath
    for system in range(6):
        coupling = .17+(.11 if system == 0 else -.11 if system == 1 else 0)
        h[2*system, 2*system+1] = h[2*system+1, 2*system] = coupling
    single = bath_channel(exponential(-.6j*h))
    shared_channel = bath_channel(exponential(-1.2j*h))
    shared = choi(shared_channel, 6)
    reset = choi(single@single, 6)
    close(unpack(row['shared_choi'], (36, 36)), shared, 'same bath retained')
    close(unpack(row['reset_choi'], (36, 36)), reset, 'bath reset countercontrol')
    cptp(shared, 6, 6)
    cptp(reset, 6, 6)
    distance = sum(abs(np.linalg.eigvalsh(shared-reset)))/2
    close(number(row['bell_half_trace_difference']), distance, 'observable bath memory')
    need(distance > 1e-3, 'independence must fail in the memory control')
    code_delta = flagged_oracle(shared_channel, 2)-flagged_oracle(single@single, 2)
    code_distance = sum(abs(np.linalg.eigvalsh(code_delta)))/2
    close(number(row['code_flagged_half_trace_difference']), code_distance, 'memory on native code inputs')
    need(code_distance > 1e-3, 'memory difference must survive allowed preparation and readback')


def verify_obstruction(rows):
    need(type(rows) is list and len(rows) == 3, 'damping obstruction catalog')
    for row, gamma in zip(rows, (.2, .5, .8)):
        keys(row, 'damping subtracted_identity_weight witness quadratic_form')
        close(number(row['damping']), gamma, 'damping value')
        close(number(row['subtracted_identity_weight']), .25, 'identity trial weight')
        w = unpack(row['witness'], (4, 1))[:, 0]
        close(w, [math.sqrt(1-gamma), 0, 0, -1], 'orthogonal support witness')
        value = -.25*(1-math.sqrt(1-gamma))**2
        close(number(row['quadratic_form']), value, 'negative CP witness')
        need(value < 0, 'no positive ideal-channel component')


def verify_calibration(rows):
    need(type(rows) is list and len(rows) == 4, 'calibration catalog')
    for row, drive in zip(rows, (1, 2, 4, 16)):
        keys(row, 'drive duration relative_error half_trace_distance')
        integer(row['drive'], drive, drive)
        close(number(row['duration']), .7/drive, 'fixed pulse angle')
        close(number(row['relative_error']), .04, 'fixed fractional drive error')
        close(number(row['half_trace_distance']), abs(math.sin(.04*.7/2)), 'speed cannot remove relative error')


def verify(row):
    keys(row, 'channels bath ideal_mixture_obstruction calibration')
    verify_channels(row['channels'])
    verify_bath(row['bath'])
    verify_obstruction(row['ideal_mixture_obstruction'])
    verify_calibration(row['calibration'])
