"""Exact finite coarse density evolution with memory and initial hidden forcing.

The parameter is the number of uniformly selected native seam attempts.  All
vectors and arithmetic are sparse and rational; no Gaussian closure is used.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import cached_property
import hashlib
import json

from .source import federation
from .refinement import parent_map, nonlumpability_witness
from .refinement import SOURCE_PATHS as REFINEMENT_SOURCES


SOURCE_PATHS = (
    *REFINEMENT_SOURCES,
    'model/federated_observer/memory.py',
    'Lean/SourceLaw/IdentificationLadder.lean',
    'model/local_tick/README.md',
)


def _vector(values, size):
    if type(values) is not dict:
        raise ValueError('an ordinary sparse dictionary is required')
    result = {}
    for key, value in values.items():
        if (type(key) is not int or not 0 <= key < size
                or type(value) not in (int, F)):
            raise ValueError('integer site indices and exact rational values required')
        if value:
            result[key] = F(value)
    return result


def _add(*vectors):
    result = {}
    for vector in vectors:
        for key, value in vector.items():
            result[key] = result.get(key, F(0)) + value
    return {key: value for key, value in result.items() if value}


def _negative(vector):
    return {key: -value for key, value in vector.items()}


def _sup(vector):
    return max(map(abs, vector.values()), default=F(0))


def _rms2(vector, size):
    return sum((value * value for value in vector.values()), F(0)) / size


def _fraction(value):
    return {'numerator': value.numerator, 'denominator': value.denominator}


def _encode(vector):
    return [[key, value.numerator, value.denominator]
            for key, value in sorted(vector.items())]


@dataclass(frozen=True)
class Reduction:
    """R averages descendants, J repeats values, and E=JR projects fields.

    Density profiles and site observables share P because this native transfer
    is symmetric.  R produces block MEANS; block particle counts are
    fibre_size*R v.
    """
    coarse: int = 0
    fine: int = 1

    def __post_init__(self):
        # parent_map also rejects booleans and out-of-domain levels.
        parent_map(self.coarse, self.fine)

    @cached_property
    def graph(self):
        return federation(self.fine)

    @cached_property
    def parents(self):
        return parent_map(self.coarse, self.fine)

    @cached_property
    def coarse_size(self):
        return federation(self.coarse).ports

    @cached_property
    def fibre_size(self):
        return 4 ** (self.fine - self.coarse)

    @cached_property
    def fibres(self):
        rows = [[] for _ in range(self.coarse_size)]
        for vertex, parent in enumerate(self.parents):
            rows[parent].append(vertex)
        return tuple(tuple(row) for row in rows)

    def _p(self, values):
        denominator = 2 * len(self.graph.seams)
        result = {}
        for vertex, value in values.items():
            row = self.graph.neighbours[vertex]
            result[vertex] = (result.get(vertex, F(0))
                              + value * F(denominator-len(row), denominator))
            contribution = value / denominator
            for target in row:
                result[target] = result.get(target, F(0)) + contribution
        return {key: value for key, value in result.items() if value}

    def _r(self, values):
        result = {}
        for vertex, value in values.items():
            parent = self.parents[vertex]
            result[parent] = result.get(parent, F(0)) + value / self.fibre_size
        return {key: value for key, value in result.items() if value}

    def _j(self, values):
        return {vertex: value for parent, value in values.items()
                for vertex in self.fibres[parent]}

    def _f(self, values):
        return _add(values, _negative(self._j(self._r(values))))

    def _a(self, values):
        return self._r(self._p(self._j(values)))

    def _b(self, values):
        return self._r(self._p(values))

    def _c(self, values):
        return self._f(self._p(self._j(values)))

    def _d(self, values):
        return self._f(self._p(values))

    def step(self, values):
        """P=I-L_f/(2 M_f), one native attempt, on a fine density profile."""
        return self._p(_vector(values, self.graph.ports))

    def restrict(self, values):
        return self._r(_vector(values, self.graph.ports))

    def lift(self, values):
        return self._j(_vector(values, self.coarse_size))

    def hidden(self, values):
        return self._f(_vector(values, self.graph.ports))

    def a(self, values):
        """A=RPJ comes from the fine graph, not the coarse native transfer."""
        return self._a(_vector(values, self.coarse_size))

    def b(self, values):
        values = _vector(values, self.graph.ports)
        if self._r(values):
            raise ValueError('B accepts a hidden vector in ker R')
        return self._b(values)

    def c(self, values):
        return self._c(_vector(values, self.coarse_size))

    def d(self, values):
        values = _vector(values, self.graph.ports)
        if self._r(values):
            raise ValueError('D accepts a hidden vector in ker R')
        return self._d(values)

    def predict(self, initial, steps):
        """Use the initial profile and coarse history only, without full replay.

        The matrix-free kernel factors still use fine vectors; this is exact
        predictive closure with memory, not a cost-independent coarse model.
        """
        if type(steps) is not int or not 0 <= steps <= 16:
            raise ValueError('an integer horizon from zero through sixteen required')
        initial = _vector(initial, self.graph.ports)
        profiles = [self._r(initial)]
        force = self._f(initial)
        memories = []
        rows = []
        for n in range(steps):
            current = profiles[-1]
            markov = self._a(current)
            forcing = self._b(force)
            memory = _add(*(self._b(vector) for vector in memories))
            profiles.append(_add(markov, forcing, memory))
            rows.append({'attempt': n+1, 'markov': markov,
                         'forcing': forcing, 'memory': memory,
                         'hidden_forcing_rms_squared': _rms2(force, self.graph.ports)})
            # At step n these vectors are D^(n-1-j) C x_j, j<n.
            memories = [self._d(vector) for vector in memories]
            memories.append(self._c(current))
            force = self._d(force)
        return {'profiles': profiles, 'terms': rows}


def audit_preparation(reduction, initial, steps=8):
    """Compare the independently propagated fine profile at every coordinate."""
    initial = _vector(initial, reduction.graph.ports)
    if any(not 0 <= value <= 1 for value in initial.values()):
        raise ValueError('a physical binary occupation mean profile is required')
    predicted = reduction.predict(initial, steps)
    full = initial
    averaged_initial_forecast = reduction.lift(reduction.restrict(initial))
    memoryless = reduction.restrict(initial)
    initial_hidden_rms2 = _rms2(reduction.hidden(initial), reduction.graph.ports)
    rows = []
    for n, target in enumerate(predicted['profiles']):
        observed = reduction.restrict(full)
        residual = _sup(_add(target, _negative(observed)))
        assert residual == 0
        assert sum(full.values(), F(0)) == sum(initial.values(), F(0))
        assert all(0 <= value <= 1 for value in full.values())
        omitted_hidden_error = _rms2(_add(observed, _negative(
            reduction.restrict(averaged_initial_forecast))), reduction.coarse_size)
        assert omitted_hidden_error <= initial_hidden_rms2
        row = {
            'attempts': n,
            'complete_coarse_profile_residual': _fraction(residual),
            'memoryless_profile_sup_error': _fraction(_sup(_add(memoryless, _negative(observed)))),
            'total_expected_particles': _fraction(sum(full.values(), F(0))),
            'coarse_forecast_error_if_hidden_initial_profile_omitted_rms_squared':
                _fraction(omitted_hidden_error),
            'coarse_profile_sha256': hashlib.sha256(json.dumps(
                _encode(target), separators=(',', ':')).encode()).hexdigest(),
        }
        if n:
            terms = predicted['terms'][n-1]
            forcing_rms2 = _rms2(terms['forcing'], reduction.coarse_size)
            assert forcing_rms2 <= initial_hidden_rms2
            row.update({
                'memory_term_sup': _fraction(_sup(terms['memory'])),
                'initial_forcing_sup': _fraction(_sup(terms['forcing'])),
                'initial_forcing_rms_squared': _fraction(forcing_rms2),
            })
        rows.append(row)
        full = reduction.step(full)
        averaged_initial_forecast = reduction.step(averaged_initial_forecast)
        memoryless = reduction.a(memoryless)
    return {
        'initial_profile': _encode(initial),
        'initial_hidden_rms_squared': _fraction(initial_hidden_rms2),
        'all_coordinate_equality_is_exact': True,
        'finite_horizons': rows,
    }


def run():
    reduction = Reduction(0, 1)
    witness = nonlumpability_witness(0, 1)
    target = witness['target_parent_port']
    initial = reduction.lift({target: F(1, 2)})
    assert not reduction.hidden(initial)
    block = audit_preparation(reduction, initial, 8)
    assert block['finite_horizons'][1]['memoryless_profile_sup_error']['numerator'] == 0
    assert block['finite_horizons'][2]['memoryless_profile_sup_error']['numerator'] > 0
    hidden_preparations = [{v: F(1) for v in state}
                           for state in witness['initial_fixed_two_states']]
    first, second = hidden_preparations
    assert reduction.restrict(first) == reduction.restrict(second)
    difference = _add(reduction.restrict(reduction.step(first)),
                      _negative(reduction.restrict(reduction.step(second))))
    assert difference[target] != 0
    return {
        'parameter': 'integer count of uniformly selected native seam attempts on level 1',
        'fine_ports': reduction.graph.ports,
        'coarse_ports': reduction.coarse_size,
        'fibre_size': reduction.fibre_size,
        'native_attempt_denominator': 2 * len(reduction.graph.seams),
        'restriction_output': 'mean occupation per fine descendant; block count is four times this mean',
        'recurrence': 'x[n+1]=A*x[n]+B*D^n*y[0]+sum(j=0..n-1) B*D^(n-1-j)*C*x[j]',
        'block_constant_fixed_two_preparation': {
            'description': 'uniform mixture of the six two-particle configurations in one four-port descendant fibre',
            'coarse_parent_port': target,
            'certificate': block,
        },
        'same_coarse_initial_state_hidden_preparation_control': {
            'initial_coarse_profiles_equal': True,
            'target_parent_port': target,
            'one_attempt_target_mean_difference': _fraction(difference[target]),
            'one_attempt_target_count_difference': _fraction(reduction.fibre_size*difference[target]),
            'preparations': [audit_preparation(reduction, value, 4)
                             for value in hidden_preparations],
        },
        'all_input_structure': {
            'restriction_lift_identity': True,
            'uniform_site_inner_products_make_R_adjoint_to_J': True,
            'positive_contraction_P_and_hidden_compression_D': True,
            'P_spectral_lower_bound': _fraction(1-F(
                reduction.graph.max_degree, len(reduction.graph.seams))),
            'memory_kernel_C_adjoint_D_power_C_positive_semidefinite': True,
            'kernel_entries_are_not_asserted_nonnegative': True,
            'initial_forcing_rms_bound': '||B D^n y0||_coarse <= ||y0||_fine',
            'omitting_initial_hidden_profile_forecast_bound':
                '||R P^n (v0-J R v0)||_coarse <= ||y0||_fine',
        },
        'scope': {
            'exact_finite_first_moment_closure_with_memory': True,
            'initial_hidden_profile_or_conditional_preparation_rule_supplied': True,
            'predictor_stores_fine_vectors_for_kernel_application': True,
            'coarse_native_transfer_identified_with_A': False,
            'full_joint_law_or_quantum_instrument_recovered': False,
            'uniform_memory_decay_or_continuum_limit_proved': False,
        },
    }
