"""Independent stabilizer-projector and continuous-channel checks."""

import math

import numpy as np

from m1_fermionic_source.check import complex_array, close, exact, keys, need


def pauli(x, z):
    matrix = np.zeros((128, 128))
    for column in range(128):
        matrix[column ^ x, column] = (-1.)**((column & z).bit_count() % 2)
    return matrix


def check_decoder(rows):
    need(type(rows) is list and len(rows) == 64, 'all 64 decoder outcomes')
    masks = [sum(1 << (j-1) for j in range(1, 8) if j & (1 << bit)) for bit in range(3)]
    generators = [pauli(h, 0) for h in masks]+[pauli(0, h) for h in masks]
    # Construct the zero codeword by projecting a seed, not by the producer's
    # classical row-span construction.
    state = np.eye(128)[:, 0]
    for g in generators:
        state = (state+g@state)/2
    state /= np.linalg.norm(state)
    code = np.column_stack((state, pauli(127, 0)@state))
    maps = []
    for s, row in enumerate(rows):
        keys(row, 'syndrome correction kraus_rows denominator_squared')
        exact(row['syndrome'], s)
        exact(row['denominator_squared'], 8)
        x = 2**(s//8-1) if s//8 else 0
        z = 2**(s % 8-1) if s % 8 else 0
        exact(row['correction'], [x, z])
        project = np.eye(128)
        for k, g in enumerate(generators):
            # Multiplication by a signed permutation avoids dense cubic work.
            perm = np.argmax(np.abs(g), axis=1)
            signs = g[np.arange(128), perm]
            project = (project+(-1 if s >> k & 1 else 1)*signs[:, None]*project[perm])/2
        intended = code.T@pauli(x, z).T@project
        need(type(row['kraus_rows']) is list and len(row['kraus_rows']) == 2, 'logical output rows')
        emitted = np.zeros((2, 128))
        for j, sparse in enumerate(row['kraus_rows']):
            need(type(sparse) is list and len(sparse) == 8, 'all sparse Kraus entries')
            previous = -1
            for entry in sparse:
                need(type(entry) is list and len(entry) == 2, 'sparse entry')
                i, sign = entry
                need(type(i) is int and previous < i < 128 and type(sign) is int and sign in (-1, 1),
                     'unique canonical signed entries')
                previous = i
                emitted[j, i] = sign/math.sqrt(8)
        close(emitted, intended, 'complete decoder branch')
        maps.append(emitted)
    close(sum(k.T@k for k in maps), np.eye(128), 'unconditional trace preservation')
    # Knill-Laflamme channel identity, including arbitrary coherent linear
    # combinations of single-position errors and arbitrary spectators.
    errors = [np.eye(128)]+[pauli(x << j, z << j) for j in range(7) for x, z in ((1, 0), (0, 1), (1, 1))]
    for error in errors:
        for k in maps:
            out = k@error@code
            close(out, np.eye(2)*np.trace(out)/2, 'every single-error branch preserves the logical qubit')
    return code, maps


def check_interfaces(row):
    keys(row, 'h t cx')
    close(complex_array(row['h'], (2, 2)), np.array([[1, 1], [1, -1]])/math.sqrt(2), 'native H')
    close(complex_array(row['t'], (2, 2)), np.diag([1., (1+1j)/math.sqrt(2)]), 'native T')
    expected = np.zeros((4, 4))
    for bits in range(4):
        expected[bits ^ (2 if bits & 1 else 0), bits] = 1
    close(complex_array(row['cx'], (4, 4)), expected, 'native controlled X')


def check_memory(row, code, maps):
    keys(row, 'logical_error_by_weight raw_entangled_fidelity_by_weight single_phase_branches')
    bad, raw = [0]*8, [0]*8
    branches = []
    for mask in range(128):
        error = pauli(0, mask)
        channels = [k@error@code for k in maps]
        trace_identity = sum(abs(np.trace(v))**2/4 for v in channels)
        w = mask.bit_count()
        bad[w] += round(1-trace_identity)
        raw[w] += round(abs(np.trace(code.T@error@code))**2/4)
        if mask == 0 or mask.bit_count() == 1:
            active = [s for s, v in enumerate(channels) if np.linalg.norm(v) > .1]
            need(len(active) == 1, 'deterministic single-error syndrome')
            branches.append((mask, active[0], channels[active[0]]))
    exact(row['logical_error_by_weight'], bad)
    exact(row['raw_entangled_fidelity_by_weight'], raw)
    need(type(row['single_phase_branches']) is list and len(row['single_phase_branches']) == len(branches),
         'all single-phase histories')
    for emitted, (mask, s, matrix) in zip(row['single_phase_branches'], branches):
        keys(emitted, 'z syndrome logical')
        exact([emitted['z'], emitted['syndrome']], [mask, s])
        close(complex_array(emitted['logical'], (2, 2)), matrix, 'actual recovered logical channel')


def exponential(matrix):
    """Scaling/squaring Taylor replay, independent of scipy.linalg.expm."""
    scale = max(0, math.ceil(math.log2(max(1., np.linalg.norm(matrix, ord=np.inf)*8))))
    small = matrix/2**scale
    term = np.eye(len(matrix), dtype=complex)
    value = term.copy()
    for k in range(1, 45):
        term = term@small/k
        value += term
    for _ in range(scale):
        value = value@value
    return value


def check_noise(rows):
    need(type(rows) is list and len(rows) == 2, 'driven M6 noise catalogue')
    for row, (kind, rate, t) in zip(rows, [('mix', .31, .7), ('pair_phase', .23, .41)]):
        keys(row, 'kind rate duration zero_jump_probability channel_choi fault_choi')
        exact([row['kind'], row['rate'], row['duration']], [kind, rate, t])
        d = 6
        h = np.zeros((d, d), dtype=complex)
        if kind == 'mix':
            h[0, 1], h[1, 0] = -.7j, .7j
        else:
            for j in range(3):
                h[j, j] = math.pi
        # Build the generator on each matrix unit directly.
        generator = np.zeros((d*d, d*d), dtype=complex)
        for col in range(d*d):
            e = np.zeros((d, d), dtype=complex)
            e[col % d, col//d] = 1
            g = -1j*(h@e-e@h)+rate*(np.diag(np.diag(e))-e)
            generator[:, col] = g.reshape(-1, order='F')
        channel = exponential(t*generator)
        expected = np.zeros((d*d, d*d), dtype=complex)
        for i in range(d):
            for j in range(d):
                expected[i*d:(i+1)*d, j*d:(j+1)*d] = channel[:, i+d*j].reshape((d, d), order='F')/d
        emitted = complex_array(row['channel_choi'], (d*d, d*d))
        close(emitted, expected, 'continuous driven channel')
        u = exponential(-1j*t*h)
        bell = u.T.reshape(-1)/math.sqrt(d)
        ideal = np.outer(bell, bell.conj())
        zero = math.exp(-rate*t)
        exact(row['zero_jump_probability'], zero)
        fault = complex_array(row['fault_choi'], (d*d, d*d))
        close(emitted, zero*ideal+(1-zero)*fault, 'no-jump plus complete fault channel')
        for state in (emitted, fault):
            close(state, state.conj().T, 'Hermitian Choi matrix')
            need(np.linalg.eigvalsh(state).min() >= -2e-12, 'CP fault channel')
            marginal = np.array([[np.trace(state[i*d:(i+1)*d, j*d:(j+1)*d]) for j in range(d)] for i in range(d)])
            close(marginal, np.eye(d)/d, 'TP on the full six-dimensional interface')


def check_accounting(row, code, maps):
    keys(row, 'logical_effect logical_state cases abort_value before_encoding_phase_error '
              'raw_input_entangled_fidelity after_encoding_corrected_fidelity')
    effect = complex_array(row['logical_effect'], (2, 2))
    close(effect, np.array([[.35, .12+.05j], [.12-.05j, .75]]), 'fixed logical effect')
    psi = complex_array(row['logical_state'], (2,))
    close(psi, np.array([math.sqrt(.6), 1j*math.sqrt(.4)]), 'fixed logical input')
    physical = sum(k.T@effect@k for k in maps)
    need(np.linalg.eigvalsh(physical).min() > 0 and np.linalg.eigvalsh(physical).max() < 1,
         'positive bounded accounting on every syndrome')
    need(type(row['cases']) is list and len(row['cases']) == 4, 'accounting fault coverage')
    for emitted, phase in zip(row['cases'], (0, 1, 2, 4)):
        keys(emitted, 'phase decoded_expectation raw_code_penalty')
        exact(emitted['phase'], phase)
        state = pauli(0, phase)@code@psi
        exact(emitted['decoded_expectation'], float(np.vdot(state, physical@state).real))
        need(type(emitted['raw_code_penalty']) is float and math.isfinite(emitted['raw_code_penalty']),
             'finite floating-point penalty')
        need(abs(emitted['raw_code_penalty']-(0. if phase == 0 else 1.)) < 1e-12,
             'raw and decoded accounting must not be identified')
    exact(row['abort_value'], 1.)
    p = .125
    exact(row['before_encoding_phase_error'], p)
    # A phase before encoding acts as a logical phase, not a weight-one
    # physical error. Its Bell overlap is 1-p for every later ideal encoding.
    identity, logical_z = np.eye(2), np.diag([1., -1.])
    fidelity = (1-p)*abs(np.trace(identity)/2)**2+p*abs(np.trace(logical_z)/2)**2
    exact(row['raw_input_entangled_fidelity'], float(fidelity))
    repaired = sum(abs(np.trace(k@pauli(0, 1)@code)/2)**2 for k in maps)
    exact(row['after_encoding_corrected_fidelity'], float(repaired))
