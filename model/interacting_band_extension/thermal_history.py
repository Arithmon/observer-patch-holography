"""Reference-preserving occupation certificates for the interacting observer.

Finite Taylor paths are candidates, never assumed solutions. The full infinite
polynomial generator acting on each path has finite support. Outward residuals
and endpoint jumps certify its distance from the exact unitary history. All
record branches and all reference columns are retained.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial

import numpy as np
from scipy import sparse

from ..occupation_control import intervals as iv
from ..scalar_observations.bounds import sqrt_upper
from .thermal_forms import fixed_code_form, natural_lower_form
from .thermal_spectrum import certify_spectrum, exponential_interval


def _integer(value, minimum, maximum):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError('an original integer in the supported range is required')
    return value


def _fraction(value, *, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and not value):
        raise ValueError('an original exact nonnegative rational is required')
    return F(value)


def preparation_bracket(occupation, shell, beta=F(1, 4), coupling=F(1), *, bits=48):
    """Bound canonical-purification distance to the infinite Gibbs state."""
    occupation, shell = _integer(occupation, 2, 32), _integer(shell, 0, 60)
    beta, coupling = _fraction(beta, positive=True), _fraction(coupling)
    finite = certify_spectrum(fixed_code_form(occupation, coupling), beta, bits=bits)
    lower_form = certify_spectrum(natural_lower_form(shell, coupling), beta, bits=bits)
    n = shell+2
    e1, e2 = exponential_interval(-beta), exponential_interval(-2*beta)
    z0_hi = exponential_interval(-3*beta/2)[1]/((1-e1[1])*(1-e2[1]))
    outside_hi = z0_hi*(exponential_interval(-beta*n)[1]
        + exponential_interval(-2*beta*n)[1]-exponential_interval(-3*beta*n)[0])
    upper = lower_form['partition_function_interval'][1]+outside_hi
    lower = finite['partition_function_interval'][0]
    if not 0 < lower <= upper:
        raise ArithmeticError('the partition bracket is inconsistent')
    return {'occupation': occupation, 'lower_form_shell': shell,
            'inverse_temperature': beta, 'coupling': coupling,
            'finite_partition_interval': finite['partition_function_interval'],
            'infinite_partition_lower': lower, 'infinite_partition_upper': upper,
            'harmonic_complement_upper': outside_hi,
            'canonical_preparation_trace_distance_upper': sqrt_upper(1-lower/upper),
            'dyadic_bits': bits, 'eigensolver_accuracy_assumed': False}


def _state(value, rows):
    if (type(value) is not np.ndarray or value.ndim != 2 or value.shape[0] != rows
            or not value.shape[1] or value.dtype.kind not in 'fc'
            or not np.isfinite(value).all()):
        raise ValueError('a finite floating work/reference amplitude matrix is required')
    result = np.array(value, dtype=complex, copy=True)
    iv.frobenius_upper(*iv.complex_exact(result))  # Reject overflowing norms.
    return result


def _sqrt_rational(value):
    return iv.square_root(iv.rational(value))


def finite_gibbs_seed(occupation=8, beta=F(4), coupling=F(1), *, bits=48):
    """Canonical finite Gibbs amplitude with a certified numerical error.

    The dyadic candidate V has a nearby exact polar unitary U. The spectral
    certificate bounds ||H-U Lambda U*||. Duhamel for the heat semigroup,
    the certified partition interval, and exact scalar enclosures then bound
    sqrt(rho_d)-V exp(-beta Lambda/2)V*/sqrt(z_mid). An outward outer-product
    sum charges all rounding in the returned candidate matrix.
    """
    d = _integer(occupation, 2, 16)
    beta, coupling = _fraction(beta, positive=True), _fraction(coupling)
    spectrum = certify_spectrum(fixed_code_form(d, coupling), beta, bits=bits,
                                return_candidates=True)
    zlo, zhi = spectrum['partition_function_interval']
    if zlo <= 0:
        raise ArithmeticError('positive Gibbs normalization was not certified')
    zmid = (zlo+zhi)/2
    root_hi = sqrt_upper(zmid)
    root_lo = zmid/root_hi
    normalization = max(sqrt_upper(zhi/zmid)-1, 1-1/sqrt_upper(zmid/zlo))
    seed = np.zeros((d*d, d*d), complex)
    rounding, analytic = [], []
    for block in spectrum['blocks']:
        addresses = np.array(block['addresses'])
        scale = block['dyadic_scale']
        basis = np.array(block['dyadic_basis'], dtype=float)/scale
        eigenvalues = [F(value, scale) for value in block['dyadic_eigenvalues']]
        epsilon, delta = block['gram_frobenius_upper'], block['eigenvalue_error_upper']
        minimum = min(F(3, 2), min(eigenvalues))
        heat = beta/2*delta*exponential_interval(-beta*minimum/2)[1]
        polar = epsilon*(2+epsilon)*exponential_interval(-beta*min(eigenvalues)/2)[1]
        analytic.append((sqrt_upper(F(len(addresses)))*heat+polar)/root_lo)
        weights = []
        lo = np.zeros_like(basis); hi = np.zeros_like(basis)
        for column, value in enumerate(eigenvalues):
            elo, ehi = exponential_interval(-beta*value/2)
            lower, upper = elo/root_hi, ehi/root_lo
            interval = (iv.rational(lower)[0], iv.rational(upper)[1])
            weights.append(float((lower+upper)/2))
            vector = basis[:, column]
            product = iv.mul(iv.exact(vector[:, None]), iv.exact(vector[None, :]))
            lo, hi = iv.add((lo, hi), iv.mul(product, interval))
        candidate = (basis*np.array(weights))@basis.T
        rounding.append(iv.frobenius_upper(iv.sub(iv.exact(candidate), (lo, hi)),
                                           iv.exact(np.zeros_like(candidate))))
        seed[np.ix_(addresses, addresses)] = candidate
    allowance = normalization+sum(analytic, F())+sum(map(F, rounding), F())
    return seed, {'canonical_seed_vector_error_upper': allowance,
                  'normalization_error_upper': normalization,
                  'heat_and_polar_error_upper': sum(analytic, F()),
                  'arithmetic_error_upper': sum(map(F, rounding), F()),
                  'partition_function_interval': (zlo, zhi),
                  'reference_convention': 'common frequency-two Hermite basis',
                  'eigensolver_accuracy_assumed': False}


@dataclass(frozen=True)
class AmplitudeGenerator:
    occupation: int
    extended: int
    retained: np.ndarray
    terms: tuple
    matrix: object

    def apply_float(self, state):
        return self.matrix@state

    def apply_interval(self, state):
        shape = (self.extended**2, state.shape[1])
        re, im = (np.zeros(shape), np.zeros(shape)), (np.zeros(shape), np.zeros(shape))
        for rows, columns, coefficient, imaginary in self.terms:
            real = iv.mul(coefficient, iv.exact(state[columns].real))
            imag = iv.mul(coefficient, iv.exact(state[columns].imag))
            if imaginary:
                real, imag = iv.neg(imag), real
            for total, values in ((re, real), (im, imag)):
                new = iv.add((total[0][rows], total[1][rows]), values)
                total[0][rows], total[1][rows] = new
        return re, im


def generator(occupation, kind, coefficient, coupling=F(1)):
    """Exact -i*t*H or i*a*Q2, including every escaped oscillator row."""
    d = _integer(occupation, 2, 32)
    if kind not in ('hamiltonian', 'coordinate') or type(coefficient) not in (int, F):
        raise ValueError('Hamiltonian/coordinate and an exact signed coefficient required')
    g = _fraction(coupling)
    size = d+(4 if kind == 'hamiltonian' else 1)
    retained = np.array([i+size*j for j in range(d) for i in range(d)])
    index = {int(address): i for i, address in enumerate(retained)}
    if kind == 'hamiltonian':
        form = fixed_code_form(size, g)
        entries, metric = form.entries, form.metric
        factor = -F(coefficient)
    else:
        metric = tuple(factorial(i)*factorial(j) for j in range(size) for i in range(size))
        entries = []
        for j in range(d):
            for i in range(d):
                column = i+size*j
                entries.append((i+size*(j+1), column, F(1, 2)))
                if j:
                    entries.append((i+size*(j-1), column, F(j, 2)))
        factor = F(coefficient)
    # One entry per output row in each group: fancy indexing then accumulates
    # without losing duplicate rows. Group by exact oscillator displacement.
    groups = {}
    for row, column, value in entries:
        if column not in index or not value:
            continue
        displacement = (row % size-column % size, row//size-column//size)
        root = _sqrt_rational(value*value*F(metric[row], metric[column]))
        interval = iv.mul(root, iv.rational(factor*(1 if value > 0 else -1)))
        groups.setdefault(displacement, []).append((row, index[column], interval))
    terms, mr, mc, mv = [], [], [], []
    for entries in groups.values():
        rows = np.array([entry[0] for entry in entries])
        columns = np.array([entry[1] for entry in entries])
        bounds = tuple(np.array([entry[2][side] for entry in entries])[:, None] for side in (0, 1))
        if len(set(rows)) != len(rows):
            raise ArithmeticError('generator displacement has duplicate destination rows')
        terms.append((rows, columns, bounds, True))
        for row, column, (lo, hi) in entries:
            if row in index:
                mr.append(index[row]); mc.append(column); mv.append(1j*(lo+hi)/2)
    matrix = sparse.csr_matrix((mv, (mr, mc)), shape=(d*d, d*d))
    return AmplitudeGenerator(d, size, retained, tuple(terms), matrix)


def cell(state, operation, order=12):
    state = _state(state, operation.occupation**2)
    order = _integer(order, 2, 24)
    coefficients = [state]
    for k in range(1, order+1):
        coefficients.append(_state(operation.apply_float(coefficients[-1])/k, len(state)))
    finite, escaped = [], []
    exterior = np.ones(operation.extended**2, dtype=bool)
    exterior[operation.retained] = False
    endpoint_interval = iv.complex_exact(np.zeros_like(state))
    for k, coefficient in enumerate(coefficients):
        image = operation.apply_interval(coefficient)
        inside = tuple(tuple(side[operation.retained] for side in component) for component in image)
        derivative = (iv.complex_scale(iv.complex_exact(coefficients[k+1]), iv.exact(k+1))
                      if k < order else iv.complex_exact(np.zeros_like(state)))
        finite.append(float(iv.up(iv.frobenius_upper(*iv.complex_sub(derivative, inside))/(k+1))))
        outside = tuple(tuple(side[exterior] for side in component) for component in image)
        escaped.append(float(iv.up(iv.frobenius_upper(*outside)/(k+1))))
        endpoint_interval = tuple(iv.add(total, iv.exact(part)) for total, part in
                                  zip(endpoint_interval, (coefficient.real, coefficient.imag), strict=True))
    endpoint = sum(coefficients)
    jump = iv.frobenius_upper(*iv.complex_sub(iv.complex_exact(endpoint), endpoint_interval))
    return endpoint, {'finite_residual': iv.sum_up(finite),
                      'occupation_residual': iv.sum_up(escaped), 'endpoint_jump': jump}


def evolve(state, kind, coefficient, *, occupation, coupling=F(1), cells=8, order=12):
    cells = _integer(cells, 1, 4096)
    if type(coefficient) not in (int, F):
        raise ValueError('an exact signed evolution parameter is required')
    operation = generator(occupation, kind, F(coefficient)/cells, coupling)
    budgets = []
    for _ in range(cells):
        state, budget = cell(state, operation, order)
        budgets.append(budget)
    return state, {key: iv.sum_up([budget[key] for budget in budgets]) for key in budgets[0]}


def history(initial, sign=1, *, occupation=8, coupling=F(1), interval=F(1, 16), cells=8):
    """Two full reads and feedback, certified against the infinite history.

    Pointer diagonalization uses exp(+-i Q2/2) and fixed scalar coefficients.
    No finite-coordinate functional calculus is confused with a compression
    of the infinite detector. Returned records are [00,01,10,11].
    """
    d = _integer(occupation, 2, 32)
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('an original preparation sign +1 or -1 is required')
    interval = _fraction(interval)
    state = _state(initial, d*d)
    steps = []
    def advance(value, kind, coefficient, label):
        endpoint, budget = evolve(value, kind, coefficient, occupation=d,
                                  coupling=coupling, cells=cells)
        steps.append({'operation': label, **budget})
        return endpoint
    def read(value, label):
        plus = advance(value, 'coordinate', F(1, 2), label+'_plus')
        minus = advance(value, 'coordinate', -F(1, 2), label+'_minus')
        # cos(Q/2+pi/4), sin(Q/2+pi/4). Combining the two
        # controlled arms is an isometry; using their sum is conservative.
        c = _sqrt_rational(F(1, 8))
        midpoint = (c[0]+c[1])/2
        zero = ((1+1j)*plus+(1-1j)*minus)*midpoint
        one = ((1-1j)*plus+(1+1j)*minus)*midpoint
        # The multiplication by 1+-i is itself a rounded sum. Enclose
        # directly from the arm components, rather than trusting it.
        errors = []
        for output, imaginary_sign in ((zero, 1), (one, -1)):
            re = iv.sub(iv.add(iv.exact(plus.real), iv.exact(minus.real)),
                        iv.scale(iv.sub(iv.exact(plus.imag), iv.exact(minus.imag)), imaginary_sign))
            im = iv.add(iv.add(iv.exact(plus.imag), iv.exact(minus.imag)),
                        iv.scale(iv.sub(iv.exact(plus.real), iv.exact(minus.real)), imaginary_sign))
            errors.append(iv.frobenius_upper(*iv.complex_sub(iv.complex_exact(output),
                            (iv.mul(re, c), iv.mul(im, c)))))
        steps.append({'operation': label+'_pointer_rotation', 'finite_residual': 0.,
                      'occupation_residual': 0., 'endpoint_jump': iv.sum_up(errors)})
        return zero, one
    state = advance(state, 'coordinate', F(sign, 4), 'signed_preparation_kick')
    state = advance(state, 'hamiltonian', interval, 'first_interval')
    first = list(read(state, 'first_read'))
    first[1] = advance(first[1], 'coordinate', F(1, 4), 'record_one_feedback')
    joined = np.concatenate(first, axis=1)
    joined = advance(joined, 'hamiltonian', interval, 'second_interval')
    zero, one = read(joined, 'second_read')
    columns = state.shape[1]
    outputs = (zero[:, :columns], one[:, :columns], zero[:, columns:], one[:, columns:])
    budget = {key: iv.sum_up([step[key] for step in steps])
              for key in ('finite_residual', 'occupation_residual', 'endpoint_jump')}
    budget['vector_error_upper'] = iv.sum_up(list(budget.values()))
    output_norm = iv.frobenius_upper(*iv.complex_exact(np.concatenate(outputs, axis=1)))
    input_norm = iv.frobenius_upper(*iv.complex_exact(initial))
    budget['complete_output_trace_error_upper'] = float(iv.up(
        iv.up((input_norm+output_norm)/2)*budget['vector_error_upper']))
    return outputs, {'steps': steps, **budget, 'input_norm_upper': input_norm,
                     'output_norm_upper': output_norm, 'all_record_branches': 4,
                     'retained_reference_columns': columns}
