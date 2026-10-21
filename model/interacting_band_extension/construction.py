"""Native execution and matched comparison for one interacting band extension."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import factorial, prod, sqrt
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.optimize import least_squares
from scipy.sparse.linalg import expm_multiply

from . import algebra as a
from ..observer.source_scalar import compiler
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..scalar_observations import source as inherited_source
from ..scalar_observations.protocol import native_records


ROOT = Path(__file__).resolve().parents[2]
INTERVAL = F(1, 16)
STEPS = 64
SOURCE_PATHS = (
    'model/interacting_band_extension/algebra.py',
    'model/interacting_band_extension/construction.py',
    'model/observer/source_scalar/compiler.py',
    'model/scalar_observations/source.py',
    'model/scalar_observations/protocol.py',
    'model/scalar_observations/bounds.py',
    'model/relativistic_scalar/lattice.py',
    'model/observer/native.py',
    compiler.RESPONSE, compiler.CONTROL_RECEIPT)


def normalized(matrix):
    return np.array([[float(value)*sqrt(factorial(i)/factorial(j))
                      for j, value in enumerate(row)] for i, row in enumerate(matrix)])


@dataclass(frozen=True)
class Term:
    sites: tuple
    factors: tuple
    coefficient: F
    label: str

    def matrix(self):
        result = np.array([[1.]])
        for factor in reversed(self.factors):
            result = np.kron(result, normalized(factor))
        return float(self.coefficient)*result

    def norm_bound(self):
        return abs(self.coefficient)*prod(a.row_norm(matrix) for matrix in self.factors)


def terms(extended=True, matched=False, coupling=1):
    if type(extended) is not bool or type(matched) is not bool or (extended and matched):
        raise ValueError('admitted original or matched band selector required')
    g = a.rational(coupling)
    dims = a.DIMENSIONS if extended else a.LOW_DIMENSIONS
    result = []
    for mode, d in enumerate(dims):
        frequency_squared = F(1 + sum(x*x for x in a.WAVES[mode]))+(g/4 if matched else 0)
        quartic = g/(24 if mode == 0 else 16)
        p2, q2, q4 = (a.polynomial(d, 2, True), a.polynomial(d, 2), a.polynomial(d, 4))
        matrix = tuple(tuple(p2[i][j]/2+frequency_squared*q2[i][j]/2+quartic*q4[i][j]
                             for j in range(d)) for i in range(d))
        result.append(Term((mode,), (matrix,), F(1), f'onsite_{mode}'))
    coefficients = a.declared_quartic() if extended else a.quartic_coefficients(False)
    for powers, coefficient in sorted(coefficients.items()):
        selected = tuple(i for i, power in enumerate(powers) if power)
        if len(selected) == 1:
            continue  # Already included in its complete onsite generator.
        matrices = tuple(a.polynomial(dims[i], powers[i]) for i in selected)
        result.append(Term(selected, matrices, g*coefficient/24,
                           'quartic_'+'_'.join(f'{i}^{powers[i]}' for i in selected)))
    return tuple(result)


class Lattice:
    """Finite tensor code, with sparse independent full-Hamiltonian reference."""

    def __init__(self, extended=True, matched=False, coupling=1):
        self.extended, self.matched = extended, matched
        self.coupling = a.rational(coupling)
        self.dimensions = a.DIMENSIONS if extended else a.LOW_DIMENSIONS
        self.dimension = prod(self.dimensions)
        self.terms = terms(extended, matched, self.coupling)
        self.hamiltonian = sum((self.embed(term) for term in self.terms),
                               sparse.csr_matrix((self.dimension, self.dimension)))
        self.hamiltonian.eliminate_zeros()
        self.trace = float(self.hamiltonian.diagonal().sum())

    def embed(self, term):
        factors = dict(zip(term.sites, term.factors, strict=True))
        result = sparse.csr_matrix([[float(term.coefficient)]])
        for mode in reversed(range(len(self.dimensions))):
            factor = (sparse.csr_matrix(normalized(factors[mode])) if mode in factors
                      else sparse.eye(self.dimensions[mode], format='csr'))
            result = sparse.kron(result, factor, format='csr')
        return result

    def apply(self, state, sites, matrix):
        value = np.asarray(state, dtype=complex)
        if value.ndim not in (1, 2) or value.shape[0] != self.dimension:
            raise ValueError('one amplitude row per field configuration required')
        sites = tuple(sites)
        if not sites or len(set(sites)) != len(sites) or any(type(i) is not int or not 0 <= i < len(self.dimensions) for i in sites):
            raise ValueError('distinct admitted local tensor indices required')
        local_size = prod(self.dimensions[i] for i in sites)
        matrix = np.asarray(matrix)
        if matrix.shape != (local_size, local_size) or not np.isfinite(matrix).all():
            raise ValueError('a finite matrix matching the local code is required')
        tensor = value.reshape(self.dimensions+value.shape[1:], order='F')
        moved = np.moveaxis(tensor, sites, range(len(sites)))
        changed = (matrix@moved.reshape((local_size, -1), order='F')).reshape(moved.shape, order='F')
        return np.moveaxis(changed, range(len(sites)), sites).reshape(value.shape, order='F')

    def blank(self):
        result = np.zeros(self.dimension, complex)
        result[0] = 1
        return result

    def evolve(self, state, duration):
        t = float(a.rational(duration))
        return expm_multiply(-1j*t*self.hamiltonian, state, traceA=-1j*t*self.trace)

    def source_evolve(self, state, duration, steps):
        a.integer(steps)
        dt = a.rational(duration)/(2*steps)
        operations = [(term.sites, execute_term(term, dt)[0]) for term in self.terms]
        for _ in range(steps):
            for sites, matrix in operations+operations[::-1]:
                state = self.apply(state, sites, matrix)
        return state


@lru_cache(maxsize=256)
def execute_term(term, duration):
    """Contract the emitted local source pulses, never substituting a target."""
    duration = a.rational(duration)
    changes, eigenvalues, offset = [], [], 0
    for matrix in term.factors:
        values, basis = inherited_source._basis(normalized(matrix))
        bits = len(matrix).bit_length()-1
        changes.append(compiler.remap(compiler.compile_real_orthogonal(basis),
                                      list(range(offset, offset+bits))))
        eigenvalues.append(values)
        offset += bits
    values = np.array([float(term.coefficient)])
    for vector in reversed(eigenvalues):
        values = np.kron(values, vector)
    stages = [(f'inverse_factor_{i}', compiler.inverse(change)) for i, change in enumerate(changes)]
    stages.append(('product_phase', compiler.compile_diagonal(-float(duration)*values)))
    stages.extend((f'factor_{i}', change) for i, change in reversed(list(enumerate(changes))))
    target = expm(-1j*float(duration)*term.matrix())
    return inherited_source._execute(stages, offset, target)


def splitting_bound(lattice, duration, steps):
    """Uniform all-input operator bound for the local palindrome."""
    time = a.rational(duration)
    a.integer(steps)
    # Center each local term; the identical total scalar phase cancels in
    # the exact/product comparison. Norm bounds remain exact rationals.
    budget = F(0)
    for term in lattice.terms:
        if len(term.sites) == 1:
            matrix = term.factors[0]
            scalar = matrix[0][0]
            centered = tuple(tuple(value-(scalar if i == j else 0)
                                   for j, value in enumerate(row)) for i, row in enumerate(matrix))
            budget += abs(term.coefficient)*a.row_norm(centered)
        else:
            budget += term.norm_bound()
    error = budget**3*time**3/(3*steps**2)
    return {'centered_generator_norm_sum_upper': budget,
            'one_interval_operator_error_upper': min(F(2), error),
            'two_interval_complete_output_trace_distance_upper': min(F(1), 2*error)}


def history(lattice, sign=1, *, source=True, feedback=True, interval=INTERVAL,
            steps=STEPS, reset_high=False, return_states=False, initial=None):
    """Complete low two-bit law with the same disturbed field throughout."""
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('preparation sign must be plus or minus one')
    for flag in (source, feedback, reset_high, return_states):
        if type(flag) is not bool:
            raise ValueError('original boolean protocol flags required')
    d = lattice.dimensions[1]
    q = normalized(a.polynomial(d, 1))
    if source:
        kick = inherited_source.feedback(d, sign/4, 2)[0]
        kraus = inherited_source.instrument(d, 1, np.pi/2, 2)[0]
        control = inherited_source.feedback(d, 1/4, 2)[0]
        evolve = lambda v: lattice.source_evolve(v, interval, steps)
    else:
        kick, control = expm(1j*sign*q/4), expm(1j*q/4)
        values, basis = np.linalg.eigh(q)
        kraus = tuple((basis*function((values+np.pi/2)/2))@basis.T
                      for function in (np.cos, np.sin))
        evolve = lambda v: lattice.evolve(v, interval)
    prepared = lattice.apply(lattice.blank() if initial is None else initial, (1,), kick)
    first_state = evolve(prepared)
    branches = [lattice.apply(first_state, (1,), matrix) for matrix in kraus]
    # The reset control is a genuine trace-and-reprepare channel: retain
    # every high occupation as a Kraus label, rather than selecting vacuum.
    if reset_high:
        if not lattice.extended:
            raise ValueError('high reset needs the extended field')
        low_size = prod(a.LOW_DIMENSIONS)
        reset = []
        for branch in branches:
            columns = branch.reshape((low_size, -1), order='F')
            lifted = np.zeros((lattice.dimension, columns.shape[1]), complex)
            lifted[:low_size, :] = columns
            reset.append(lifted)
        branches = reset
    if feedback:
        branches[1] = lattice.apply(branches[1], (1,), control)
    evolved = [evolve(branch) for branch in branches]
    rows, states = [], []
    for first, branch in enumerate(evolved):
        for last, matrix in enumerate(kraus):
            terminal = lattice.apply(branch, (1,), matrix)
            rows.append({'bits': [first, last], 'probability': float(np.linalg.norm(terminal)**2)})
            states.append(terminal)
    return (rows, states) if return_states else rows


class ProductBand(Lattice):
    """Two radial modes of the same angular-constant product field.

    The interval, coupling, occupation cutoff and smooth-background
    interpretation are declared. This is not a closed nonlinear mode sector.
    """

    def __init__(self, coupling=1):
        self.coupling = g = a.rational(coupling)
        self.extended, self.matched = False, False
        self.dimensions, self.dimension = (4, 4), 16
        q2, q4, p2 = a.polynomial(4, 2), a.polynomial(4, 4), a.polynomial(4, 2, True)
        terms = []
        for index, frequency in enumerate((1, 2)):
            matrix = tuple(tuple(p2[i][j]/2+frequency**2*q2[i][j]/2+g*q4[i][j]/16
                                 for j in range(4)) for i in range(4))
            terms.append(Term((index,), (matrix,), F(1), f'product_onsite_{index}'))
        terms.append(Term((0, 1), (q2, q2), g/4, 'product_local_quartic'))
        self.terms = tuple(terms)
        self.hamiltonian = sum((self.embed(term) for term in self.terms),
                               sparse.csr_matrix((16, 16)))
        self.trace = float(self.hamiltonian.diagonal().sum())


def product_gibbs_preparation(lattice, beta=F(1, 4)):
    """Source preparation of this interacting Hamiltonian's purification.

    Keep all work/reference amplitudes. Numerical synthesis guards and
    inherited exact-code failure semantics are explicit; no postselection
    or state normalization is used to repair the returned amplitudes.
    """
    beta = a.rational(beta, positive=True)
    if not isinstance(lattice, ProductBand):
        raise ValueError('the declared sixteen-state product band is required')
    values, basis = inherited_source._basis(lattice.hamiltonian.toarray())
    weights = np.exp(-float(beta)*(values-values[0])); weights /= weights.sum()
    if not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ArithmeticError('faithful Gibbs weights exceed the numerical preparation range')
    word = compiler.compile_real_state(np.sqrt(weights))
    for bit in range(4):
        word += compiler.cnot(bit, bit+4)
    word += compiler.compile_real_orthogonal(basis)
    blank = np.zeros(256, complex); blank[0] = 1
    raw = compiler.execute(word, 8, blank).reshape(16, 16).T
    target = basis@np.diag(np.sqrt(weights))
    error = float(np.linalg.norm(raw-target))
    if not np.isfinite(error) or error > 1e-9:
        raise ArithmeticError('the interacting Gibbs preparation exceeded its source guard')
    return raw, target, {'inverse_temperature': beta,
        'retained_work_addresses': 16, 'retained_reference_addresses': 16,
        'source_resources': compiler.resources(word),
        'raw_norm_squared': float(np.linalg.norm(raw)**2),
        'purification_vector_error_numerical': error,
        'purification_vector_error_guard': F(1, 10**9),
        'ideal_minimum_density_eigenvalue': float(weights.min()),
        'physical_preparation_error_certified': False}


def product_interaction_control():
    """Conditional same-product interaction, Gibbs clock and full readouts."""
    lattice, free = ProductBand(), ProductBand(0)
    beta, interval, steps = F(1, 4), F(1, 16), 64
    raw, target, preparation = product_gibbs_preparation(lattice, beta)
    actual, reference, actual_states, reference_states = [], [], [], []
    for sign in (-1, 1):
        rows, states = history(lattice, sign, interval=interval, steps=steps,
                               initial=raw, return_states=True)
        ideal_rows, ideal_states = history(lattice, sign, source=False, interval=interval,
                                           initial=target, return_states=True)
        actual.append(rows); actual_states.append(states)
        reference.append(ideal_rows); reference_states.append(ideal_states)
    def distance(left, right):
        return sum(float(np.linalg.norm(np.outer(x.reshape(-1), x.reshape(-1).conj())
                                      -np.outer(y.reshape(-1), y.reshape(-1).conj()), ord='nuc'))/2
                   for x, y in zip(left, right, strict=True))
    complete_error = max(distance(x, y) for x, y in zip(actual_states, reference_states, strict=True))
    split = splitting_bound(lattice, interval, steps)
    certificate = split['two_interval_complete_output_trace_distance_upper']
    # Telescope the inherited fixed raw-target and polar-correction guards,
    # including the complete measurement dilations and both possible kicks.
    # This is a floating-replay allowance, not a physical service certificate.
    operation_count = 4*steps*len(lattice.terms)+4
    local_guard = F(251, 10**12)  # 2e-10 raw residual + 5.1e-11 correction.
    gain = operation_count*local_guard
    vector_guard = (preparation['purification_vector_error_guard']+gain)/(1-gain)
    numerical_guard = vector_guard*(2+vector_guard)/2
    if complete_error > float(certificate+numerical_guard):
        raise ArithmeticError('the full interacting work/reference/record law exceeded its bound')
    # Compare the same input and detectors after removing only the interaction.
    free_rows = [history(free, sign, source=False, interval=interval, initial=target)
                 for sign in (-1, 1)]
    no_feedback = [history(lattice, sign, source=False, feedback=False,
                           interval=interval, initial=target) for sign in (-1, 1)]
    h = lattice.hamiltonian.toarray(); rho = target@target.conj().T
    q = np.kron(normalized(a.polynomial(4, 1)), np.eye(4))
    t = float(interval)
    evolved = expm(1j*t*h)@q@expm(-1j*t*h)
    masses, vectors = np.linalg.eigh(rho)
    modular = (vectors*np.exp(-1j*t*np.log(masses)/float(beta)))@vectors.T
    modular_residual = float(np.linalg.norm(modular@q@modular.conj().T-evolved))
    raw_masses, raw_vectors = np.linalg.eigh(raw@raw.conj().T)
    raw_modular = (raw_vectors*np.exp(-1j*t*np.log(raw_masses)/float(beta)))@raw_vectors.conj().T
    raw_modular_residual = float(np.linalg.norm(raw_modular@q@raw_modular.conj().T-evolved))
    free_values, free_basis = np.linalg.eigh(free.hamiltonian.toarray())
    free_weights = np.exp(-float(beta)*free_values); free_weights /= free_weights.sum()
    free_rho = (free_basis*free_weights)@free_basis.T
    free_state_nonstationarity = float(np.linalg.norm(h@free_rho-free_rho@h))
    if modular_residual > 1e-10 or raw_modular_residual > 1e-9 or free_state_nonstationarity < 1e-5:
        raise ArithmeticError('the interacting clock/state identity or mismatched-state control failed')
    native = native_records(actual)
    resources = dict(preparation['source_resources'])
    counts = {key: resources[key] for key in inherited_source.COUNT_KEYS}
    for term in lattice.terms:
        meta = execute_term(term, interval/(2*steps))[1]
        for key in counts:
            counts[key] += 4*steps*meta['resources'][key]
    for count, meta in ((2, inherited_source.instrument(4, 1, np.pi/2, 2)[1]),
                        (2, inherited_source.feedback(4, .25, 2)[1])):
        for key in counts:
            counts[key] += count*meta['resources'][key]
    return {'scope': 'Two-mode interacting effective field on the declared radial product; finite ultraviolet and occupation cutoffs.',
        'spatial_domain': '(0,pi) times the same native angular screen',
        'field_modes': ['sqrt(2) sin(r) times angular constant', 'sqrt(2) sin(2r) times angular constant'],
        'spatial_measure': 'dr/pi times normalized native angular measure',
        'quartic_polynomial': '3/2*(Q1^4+Q2^4)+6*Q1^2*Q2^2',
        'effective_coupling': lattice.coupling, 'oscillator_basis_frequency': 2,
        'covariant_coupling_map': 'g_eff=lambda/pi integral mu_screen^2 sqrt(det a_screen) dx',
        'occupation_dimensions': lattice.dimensions, 'preparation': preparation,
        'model_time_per_interval': interval, 'steps_per_interval': steps,
        'source_histories_by_native_sign': actual,
        'independent_hamiltonian_histories_by_native_sign': reference,
        'complete_work_reference_record_trace_error_numerical': complete_error,
        'ideal_splitting_complete_output_bound': certificate,
        'complete_source_numerical_guard': numerical_guard,
        'guarded_local_operation_count': operation_count,
        'local_raw_plus_polar_operator_guard': local_guard,
        'interaction_removal_record_tv': max(tv(x, y) for x, y in zip(reference, free_rows)),
        'feedback_removal_record_tv': max(tv(x, y) for x, y in zip(reference, no_feedback)),
        'ideal_Gibbs_modular_identity_residual': modular_residual,
        'raw_prepared_marginal_modular_identity_residual': raw_modular_residual,
        'free_Gibbs_state_interacting_commutator_norm': free_state_nonstationarity,
        'modular_clock': 'sigma_s=alpha_(-beta*s) for this interacting H and its ideal Gibbs marginal',
        'clock_parameter': 'Supplied Hamiltonian model time; no event-counter or seconds calibration in this protocol.',
        'native_history_count': len(native),
        'native_history_probability_sum': sum(row['probability'] for row in native),
        'worst_branch_source_resources': counts,
        'additional_native_repairs_per_history': 2,
        'physical_service_error_certified': False,
        'conditional_smooth_background_requires_smooth_coefficient_density': True,
        'coupling_and_temperature_selected_by_repair': False,
        'finite_band_is_invariant_under_untruncated_local_interaction': False,
        'UV_limit_or_cutoff_microcausality_certified': False,
        'observed_QFT_or_laboratory_clock_identified': False}


def product_parameter_recovery(source_histories, raw_preparation):
    """Synthetic inverse check with an independent padded-oscillator oracle.

    Fit coupling, inverse temperature and radial frequency scale from the
    first interval only. A second interval is withheld. Rounded estimates
    are retained under explicit numerical guards, not statistical intervals.
    """
    annihilation = np.diag(np.sqrt(np.arange(1, 8)), 1)
    q = (annihilation+annihilation.T)/2
    p = 1j*(annihilation.T-annihilation)
    q2, q4, p2 = (matrix[:4, :4] for matrix in (q@q, np.linalg.matrix_power(q, 4), p@p))

    def oracle(parameters, interval=F(1, 16)):
        coupling, beta, scale = parameters
        h1, h2 = (p2/2+(scale*k)**2*q2/2+coupling*q4/16 for k in (1, 2))
        h = np.kron(np.eye(4), h1)+np.kron(h2, np.eye(4))+coupling*np.kron(q2, q2)/4
        lattice = ProductBand()
        lattice.hamiltonian, lattice.trace = sparse.csr_matrix(h), float(np.trace(h).real)
        values, vectors = np.linalg.eigh(h)
        weights = np.exp(-beta*(values-values[0])); weights /= weights.sum()
        initial = vectors*np.sqrt(weights)
        return np.array([[row['probability'] for row in history(
            lattice, sign, source=False, initial=initial, interval=interval)] for sign in (-1, 1)])

    target = np.array([[row['probability'] for row in rows] for rows in source_histories])
    select = lambda probabilities: probabilities[:, :3].ravel()
    starts = ([.6, .4, .8], [1.6, .15, 1.2], [2., .8, .7])
    fits = [least_squares(lambda x: select(oracle(x)-target), start,
                         bounds=([.01, .03, .5], [3., 2., 2.]),
                         xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=150)
            for start in starts]
    expected = np.array([1., .25, 1.])  # Used only for the synthetic acceptance check.
    if any(not fit.success or np.max(abs(fit.x-expected)) >= 1e-4
           or np.max(abs(fit.fun)) >= 1e-8 for fit in fits):
        raise ArithmeticError('the synthetic observer records failed parameter recovery')
    estimate = fits[0].x
    step = 1e-5
    jacobian = np.column_stack([(select(oracle(estimate+np.eye(3)[i]*step))
                                -select(oracle(estimate-np.eye(3)[i]*step)))/(2*step)
                               for i in range(3)])
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    withheld = np.array([[row['probability'] for row in history(
        ProductBand(), sign, initial=raw_preparation, interval=F(1, 8))] for sign in (-1, 1)])
    if singular_values[-1] <= 1e-5 or np.max(abs(oracle(estimate, F(1, 8))-withheld)) >= 1e-8:
        raise ArithmeticError('the local identifiability or withheld-interval check failed')
    return {'data_kind': 'synthetic source-executed complete probabilities, without sampling noise',
            'parameters': ['effective coupling', 'inverse temperature', 'radial frequency scale'],
            'estimates_rounded_to_1e_minus_4': np.round(estimate, 4).tolist(),
            'three_distinct_initial_guesses_agree_with_declared_inputs_within': F(1, 10000),
            'training_max_probability_residual_guard': F(1, 10**8),
            'training_model_interval': F(1, 16), 'withheld_model_interval': F(1, 8),
            'withheld_max_probability_error_guard': F(1, 10**8),
            'local_sensitivity_singular_values_rounded_to_1e_minus_6': np.round(singular_values, 6).tolist(),
            'global_identifiability_or_finite_sample_confidence_proved': False,
            'radial_topology_or_angular_tensor_identified': False,
            'absolute_energy_time_scale_identified': False,
            'physical_parameter_validation': False}


def tv(left, right):
    if [row['bits'] for row in left] != [row['bits'] for row in right]:
        raise ValueError('record alphabets must agree')
    return sum(abs(x['probability']-y['probability']) for x, y in zip(left, right, strict=True))/2


def reduced_output_distance(enlarged_states, low_states):
    """Full low quantum state jointly with each recorded bit pair.

    A thin QR evaluates the trace norm of a difference of rank-at-most17
    positive matrices; no branch is normalized or selected.
    """
    if len(enlarged_states) != len(low_states):
        raise ValueError('complete terminal record alphabets must agree')
    low_size = prod(a.LOW_DIMENSIONS)
    result = 0.
    for larger, smaller in zip(enlarged_states, low_states, strict=True):
        if np.shape(larger) != (prod(a.DIMENSIONS),) or np.shape(smaller) != (low_size,):
            raise ValueError('pure terminal branches in the declared field codes required')
        factors = np.column_stack((larger.reshape((low_size, -1), order='F'), smaller))
        _, triangular = np.linalg.qr(factors, mode='reduced')
        signs = np.ones(factors.shape[1])
        signs[-1] = -1
        matrix = (triangular*signs)@triangular.conj().T
        result += float(np.sum(abs(np.linalg.eigvalsh(matrix)))/2)
    return result


def exact_apply(lattice, vector, centered=False):
    """Sparse rational action on unnormalized Bargmann basis coefficients."""
    result = {}
    for occupations, value in vector.items():
        for term in lattice.terms:
            entries = [(occupations, value*term.coefficient)]
            for mode, matrix in zip(term.sites, term.factors, strict=True):
                following = []
                for state, coefficient in entries:
                    for row in range(len(matrix)):
                        factor = matrix[row][state[mode]]
                        if factor:
                            key = state[:mode]+(row,)+state[mode+1:]
                            following.append((key, coefficient*factor))
                entries = following
            for key, coefficient in entries:
                result[key] = result.get(key, F(0))+coefficient
    if centered:
        scalar = sum((term.coefficient*prod(matrix[0][0] for matrix in term.factors)
                      for term in lattice.terms), F(0))
        for key, value in vector.items():
            result[key] = result.get(key, F(0))-scalar*value
    return {key: value for key, value in result.items() if value}


def exact_norm_squared(vector):
    return sum((value*value*prod(factorial(n) for n in key) for key, value in vector.items()), F(0))


@lru_cache(maxsize=1)
def high_parity_instrument():
    """Source CNOTs write the two high occupations' parity into a blank bit."""
    parity = np.array([(index % 4+index//4) % 2 for index in range(16)])
    yes, no = np.diag(parity), np.diag(1-parity)
    target = np.block([[no, yes], [yes, no]])
    word = compiler.cnot(0, 4)+compiler.cnot(2, 4)
    matrix, metadata = inherited_source._execute([('high_parity_pointer', word)], 5, target)
    return (matrix[:16, :16], matrix[16:, :16]), metadata


def leakage_witness(lattice, steps=32, time=F(1, 1024), execute=True):
    """Positive odd-high occupation from the resonant interaction alone."""
    if not lattice.extended:
        raise ValueError('the leakage witness requires both new modes')
    time = a.rational(time, True)
    blank = {(0,)*9: F(1)}
    first = exact_apply(lattice, blank, centered=True)
    odd = {key: value for key, value in first.items() if (key[7]+key[8]) % 2}
    leading = exact_norm_squared(odd)
    second_squared = exact_norm_squared(exact_apply(lattice, first, centered=True))
    # sqrt(1/768) >=1/28, used as an exactly checked lower enclosure.
    root_lower = F(1, 28)*lattice.coupling
    if root_lower**2 > leading:
        raise ArithmeticError('invalid leading-amplitude lower bound')
    remainder = time*time*a.sqrt_upper(second_squared)/2
    split = splitting_bound(lattice, time, steps)['one_interval_operator_error_upper']
    ideal_margin = time*root_lower-remainder
    source_margin = ideal_margin-split
    result = {'model_time': time, 'steps': steps,
              'odd_high_first_derivative_norm_squared_exact': leading,
              'centered_second_derivative_norm_squared_exact': second_squared,
              'taylor_vector_remainder_upper': remainder,
              'source_splitting_vector_error_upper': split,
              'continuous_odd_high_probability_lower': max(F(0), ideal_margin)**2,
              'ideal_source_odd_high_probability_lower': max(F(0), source_margin)**2,
              'no_interaction_odd_high_probability_exact': '0',
              'high_read_is_separate_from_low_only_history_theorem': True}
    if execute:
        state = lattice.source_evolve(lattice.blank(), time, steps)
        reference = lattice.evolve(lattice.blank(), time)
        kraus, read_metadata = high_parity_instrument()
        probabilities = [float(np.linalg.norm(lattice.apply(state, (7, 8), matrix))**2)
                         for matrix in kraus]
        indices = np.arange(lattice.dimension)
        low = prod(a.LOW_DIMENSIONS)
        odd_mask = ((indices//low) % 4+(indices//(4*low)) % 4) % 2 == 1
        result.update({'source_odd_high_probability': probabilities[1],
                       'source_even_odd_probabilities': probabilities,
                       'source_parity_pointer': read_metadata,
                       'reference_odd_high_probability': float(np.linalg.norm(reference[odd_mask])**2)})
    return result


def checks_and_run():
    if a.quartic_coefficients() != a.declared_quartic():
        raise ArithmeticError('extended quartic failed independent exact Fourier integration')
    extended, matched, bare = Lattice(), Lattice(False, True), Lattice(False, False)
    low_size = matched.dimension
    block = extended.hamiltonian[:low_size, :low_size]
    constant = a.residual_bound()['matched_scalar_increment']
    residual = block-matched.hamiltonian-float(constant)*sparse.eye(low_size)
    projected_identity_residual = float(max(abs(residual.data), default=0))
    outputs, terminal_states = {}, {}
    for name, lattice in (('extended', extended), ('matched_low', matched), ('bare_low', bare)):
        reference = [history(lattice, sign, source=False) for sign in (-1, 1)]
        actual = None
        if name != 'bare_low':
            executions = [history(lattice, sign, return_states=True) for sign in (-1, 1)]
            actual = [execution[0] for execution in executions]
            terminal_states[name] = [execution[1] for execution in executions]
        outputs[name] = {'reference_histories_by_native_sign': reference,
                         'source_histories_by_native_sign': actual,
                         'dimension': lattice.dimension,
                         'splitting_certificate': splitting_bound(lattice, INTERVAL, STEPS)}
    high_reset = [history(extended, sign, source=False, reset_high=True) for sign in (-1, 1)]
    no_feedback = [history(extended, sign, source=False, feedback=False) for sign in (-1, 1)]
    reference = outputs['extended']['reference_histories_by_native_sign']
    actual = outputs['extended']['source_histories_by_native_sign']
    low = outputs['matched_low']['source_histories_by_native_sign']
    matching = a.history_bound(2*INTERVAL)
    product_extended = outputs['extended']['splitting_certificate']['two_interval_complete_output_trace_distance_upper']
    product_low = outputs['matched_low']['splitting_certificate']['two_interval_complete_output_trace_distance_upper']
    combined = min(F(1), matching['reduced_complete_output_trace_distance_upper']+product_extended+product_low)
    witness = leakage_witness(extended)
    resource_counts = {key: 0 for key in inherited_source.COUNT_KEYS}
    maximum_local_residual = 0.
    for term in extended.terms:
        metadata = execute_term(term, INTERVAL/(2*STEPS))[1]
        maximum_local_residual = max(maximum_local_residual, metadata['full_phase_target_residual'])
        for key in resource_counts:
            resource_counts[key] += 4*STEPS*metadata['resources'][key]
    for multiplicity, metadata in ((1, inherited_source.feedback(4, .25, 2)[1]),
                                   (1, inherited_source.feedback(4, .25, 2)[1]),
                                   (2, inherited_source.instrument(4, 1, np.pi/2, 2)[1])):
        maximum_local_residual = max(maximum_local_residual, metadata['full_phase_target_residual'])
        for key in resource_counts:
            resource_counts[key] += multiplicity*metadata['resources'][key]
    native = native_records(actual)
    comparisons = {
        'source_to_same_enlarged_reference_tv': max(tv(x, y) for x, y in zip(actual, reference)),
        'source_enlarged_to_matched_low_tv': max(tv(x, y) for x, y in zip(actual, low)),
        'source_enlarged_to_matched_complete_low_quantum_record_trace_distance': max(
            reduced_output_distance(x, y) for x, y in zip(terminal_states['extended'], terminal_states['matched_low'])),
        'reference_enlarged_to_matched_low_tv': max(tv(x, y) for x, y in zip(reference, outputs['matched_low']['reference_histories_by_native_sign'])),
        'reference_enlarged_to_unmatched_bare_low_tv': max(tv(x, y) for x, y in zip(reference, outputs['bare_low']['reference_histories_by_native_sign'])),
        'trace_and_reset_high_history_tv': max(tv(x, y) for x, y in zip(reference, high_reset)),
        'feedback_deletion_history_tv': max(tv(x, y) for x, y in zip(reference, no_feedback)),
        'ideal_source_enlarged_to_matched_source_trace_distance_upper': combined,
    }
    checks = {
        'exact_quartic_integral': a.quartic_coefficients() == a.declared_quartic(),
        'vacuum_block_is_matched_low_action_plus_scalar': projected_identity_residual < 1e-12,
        'all_complete_source_laws_normalized': all(abs(sum(row['probability'] for row in law)-1) < 2e-9 for law in actual+low),
        'source_error_within_analytic_split_budget': comparisons['source_to_same_enlarged_reference_tv'] < float(product_extended)+2e-9,
        'matched_history_within_analytic_combined_budget': comparisons['source_enlarged_to_matched_low_tv'] < float(combined)+2e-9,
        'matched_complete_quantum_output_within_analytic_combined_budget': comparisons['source_enlarged_to_matched_complete_low_quantum_record_trace_distance'] < float(combined)+2e-9,
        'strictly_positive_certified_source_resonance': witness['ideal_source_odd_high_probability_lower'] > 0,
        'source_resonance_above_certified_lower_bound': witness['source_odd_high_probability'] >= float(witness['ideal_source_odd_high_probability_lower'])-1e-15,
        'continuing_high_memory_is_observable': comparisons['trace_and_reset_high_history_tv'] > 1e-10,
        'finite_history_matching_bound_is_nonvacuous': combined < 1,
    }
    if not all(checks.values()):
        raise ArithmeticError(f'band-extension checks failed: {checks}')
    product = product_interaction_control()
    product['synthetic_parameter_recovery'] = product_parameter_recovery(
        product['source_histories_by_native_sign'], product_gibbs_preparation(ProductBand())[0])
    return a.serializable({
        'schema': 'oph.interacting_band_extension.v1', 'all_pass': True,
        'scope': 'One supplied interacting Fourier-band extension, finite occupation, matched low-only histories; no ultraviolet, Lorentzian or empirical completion.',
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'modes': a.MODES, 'occupation_dimensions': a.DIMENSIONS,
        'grid_without_quartic_aliasing_minimum': 13,
        'complete_quartic_monomials': [{'powers': key, 'coefficient': value} for key, value in sorted(a.quartic_coefficients().items())],
        'model_time_per_interval': INTERVAL, 'steps_per_interval': STEPS,
        'matched_projection_matrix_residual': projected_identity_residual,
        'matching_certificate': matching, 'executed_comparisons': comparisons,
        'nine_mode_spatial_grid_certificates': [a.grid_history_bound(n, 2*INTERVAL)
                                               for n in (13, 33, 65, 129)],
        'models': outputs, 'high_reset_reference_histories': high_reset,
        'no_feedback_reference_histories': no_feedback,
        'resonant_leakage_witness': witness,
        'same_product_interaction': product,
        'complete_native_history_count': len(native),
        'complete_native_history_probability_sum': sum(row['probability'] for row in native),
        'complete_native_history_law': native,
        'worst_branch_source_resources': resource_counts,
        'additional_native_repairs_per_history': 2,
        'maximum_full_phase_local_source_residual': maximum_local_residual,
        'checks': checks,
        'acceptance': {'one_controlled_interacting_band_extension': True,
                       'quadratic_bound_for_reduced_low_quantum_and_record_output': True,
                       'high_mode_reset_or_postselection_used_in_principal_protocol': False,
                       'occupation_limit_of_enlarged_model_certified': False,
                       'ultraviolet_band_limit_or_renormalized_qft_constructed': False,
                       'finite_occupation_code_translation_covariance_claimed': False,
                       'lorentz_invariance_or_empirical_clock_identified': False}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--output', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('--output requires a new path')
    receipt = checks_and_run()
    if args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(canonical(receipt))
    if args.check:
        compare_receipt(json.loads(args.check.read_text(encoding='utf-8')), receipt)
    display = dict(receipt['executed_comparisons'])
    key = 'ideal_source_enlarged_to_matched_source_trace_distance_upper'
    display[key] = float(F(display[key]))
    print(json.dumps({'all_pass': receipt['all_pass'], 'checks': receipt['checks'],
                      'comparison': display}, indent=2))


if __name__ == '__main__':
    main()
