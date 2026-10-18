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
            steps=STEPS, reset_high=False, return_states=False):
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
    prepared = lattice.apply(lattice.blank(), (1,), kick)
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
