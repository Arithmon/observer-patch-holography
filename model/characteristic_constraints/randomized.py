"""Frozen Gram witnesses with random settings and history-dependent states.

All statistical thresholds are rational outward bounds. A physical conclusion
requires justified per-trial conditional service bounds; no laboratory records
are generated or inferred here.
"""
from dataclasses import dataclass
from fractions import Fraction as F

from ..calibrated_measurement.bounds import log_interval, root_interval
from . import bounds

SETTINGS = ('a_real', 'a_imag', 'b_real', 'b_imag', 'w_real', 'w_imag')


def _rational(value, name, nonnegative=False):
    if type(value) is not F or (nonnegative and value < 0):
        raise ValueError(f'exact {"nonnegative " if nonnegative else ""}Fraction {name} required')
    return value


@dataclass(frozen=True)
class Design:
    constant: F
    coefficients: tuple[F, ...]
    algebra_error: F = F(0)
    vector: tuple = ()
    phase: F = F(0)

    def __post_init__(self):
        _rational(self.constant, 'constant', True)
        _rational(self.algebra_error, 'algebra error', True)
        if type(self.coefficients) is not tuple or len(self.coefficients) != len(SETTINGS):
            raise ValueError('six exact witness coefficients required')
        for value in self.coefficients:
            _rational(value, 'coefficient')
        if not self.width:
            raise ValueError('a nonconstant witness required')
        certified = _witness_parameters(self.vector, self.phase)
        if (self.constant, self.coefficients, self.algebra_error) != certified:
            raise ValueError('design must equal the certified Gram witness and phase enclosure')

    @property
    def width(self):
        return sum(map(abs, self.coefficients), F(0))

    @property
    def probabilities(self):
        return tuple(abs(c)/self.width for c in self.coefficients)

    def score(self, setting, outcome):
        if type(setting) is not int or not 0 <= setting < 6 or self.coefficients[setting] == 0:
            raise ValueError('an active setting index required')
        if type(outcome) is not int or outcome not in (-1, 1):
            raise ValueError('a complete signed binary outcome required')
        sign = 1 if self.coefficients[setting] > 0 else -1
        return self.constant+self.width*sign*outcome

    def value(self, means):
        if type(means) is not tuple or len(means) != 6:
            raise ValueError('six exact real means required')
        for mean in means:
            _rational(mean, 'mean')
            if abs(mean) > 1:
                raise ValueError('binary means lie in [-1,1]')
        return self.constant+sum((c*m for c, m in zip(self.coefficients, means)), F(0))

    def allowance(self, service_errors):
        if type(service_errors) is not tuple or len(service_errors) != 6:
            raise ValueError('six conditional service-error bounds required')
        for error in service_errors:
            if error is not None:
                _rational(error, 'conditional service error', True)
        if any(error is None for error in service_errors):
            return None
        return self.algebra_error+sum((abs(c)*e for c, e in zip(self.coefficients, service_errors)), F(0))

    def specification(self):
        return {'constant': str(self.constant), 'coefficients': list(map(str, self.coefficients)),
                'setting_probabilities': list(map(str, self.probabilities)),
                'score_half_range': str(self.width), 'algebra_error_upper': str(self.algebra_error),
                'witness_vector': [[str(x), str(y)] for x, y in self.vector],
                'Weyl_phase': str(self.phase)}


def _witness_parameters(vector, phase):
    """v*Gv with <UV>=exp(i phase)*w; a,b,w are the measured means.

    Coefficients multiplying w are rationally enclosed. Their midpoint
    replacement error is added explicitly, rather than silently treating
    binary64 trigonometry as the exact Weyl phase.
    """
    if type(vector) is not tuple or len(vector) != 3:
        raise ValueError('three exact complex rational witness entries required')
    for item in vector:
        if type(item) is not tuple or len(item) != 2:
            raise ValueError('complex rational pairs required')
        for value in item:
            _rational(value, 'witness entry')
    _rational(phase, 'Weyl phase')
    ca, cb, cw = (bounds.mul(bounds.conj(vector[i]), vector[j]) for i, j in ((0, 1), (0, 2), (1, 2)))
    real = bounds.phase_real_interval(cw, phase)
    # Im(cw exp(i phase)) = Re((-i cw) exp(i phase)).
    imag = bounds.phase_real_interval((cw[1], -cw[0]), phase)
    coefficients = (2*ca[0], 2*ca[1], 2*cb[0], -2*cb[1], sum(real), -sum(imag))
    error = real[1]-real[0]+imag[1]-imag[0]
    return sum(map(bounds.norm2, vector), F(0)), coefficients, error


def witness_design(vector, phase=F(0)):
    """The only accepted coefficients are those of this positive Gram form."""
    return Design(*_witness_parameters(vector, phase), vector=vector, phase=phase)


def radius(design, trials, failure_probability):
    if type(design) is not Design or type(trials) is not int or not 1 <= trials <= 10**12:
        raise ValueError('a frozen design and 1..10^12 predetermined trials required')
    _rational(failure_probability, 'failure probability')
    if not 0 < failure_probability < 1:
        raise ValueError('failure probability must lie in (0,1)')
    return design.width*root_interval(2*log_interval(1/failure_probability)[1]/trials)[1]


def decision(design, score_sum, trials, *, failure_probability=F(1, 100), service_errors):
    """Fixed-horizon one-sided test; no optional stopping or discarded outcomes."""
    threshold = radius(design, trials, failure_probability)
    _rational(score_sum, 'score sum')
    if not trials*(design.constant-design.width) <= score_sum <= trials*(design.constant+design.width):
        raise ValueError('score sum outside the complete-outcome range')
    signed_total = (score_sum-trials*design.constant)/design.width
    if signed_total.denominator != 1 or (signed_total.numerator+trials) % 2:
        raise ValueError('score sum is not attainable by the declared binary outcomes')
    allowance = design.allowance(service_errors)
    return {'status': 'unbounded_services' if allowance is None else 'bounded',
            'mean_score': str(score_sum/trials), 'radius_upper': str(threshold),
            'allowance_upper': None if allowance is None else str(allowance),
            'failure_probability_upper': str(failure_probability),
            'rejects_joint_premises': allowance is not None and score_sum/trials+allowance+threshold < 0}


def crossing_budget(design, target_mean, *, failure_probability=F(1, 200), service_errors):
    """Budget for crossing if the future score reaches target, not test power."""
    _rational(target_mean, 'target mean')
    # Validate probability and design even if the service bounds are absent.
    radius(design, 1, failure_probability)
    if not design.constant-design.width <= target_mean <= design.constant+design.width:
        raise ValueError('target mean outside the complete-score range')
    allowance = design.allowance(service_errors)
    if allowance is None or target_mean+allowance >= 0:
        return None
    gap = -target_mean-allowance
    exact_threshold = 2*design.width**2*log_interval(1/failure_probability)[1]/gap**2
    trials = exact_threshold.__floor__()+1
    # Account also for outward root-grid rounding in the implemented decision.
    while radius(design, trials, failure_probability) >= gap:
        trials += 1
    return trials


def empirical_designs(scored, axes, interior):
    """Public means select future designs only; no test is applied to them."""
    rows = {}
    for extent, result in scored.items():
        worst = result['worst_witness']
        i, j, sign = (worst[key] for key in ('i', 'j', 'product_order_sign'))
        vector = tuple(tuple(F(x) for x in item) for item in worst['witness']['vector'])
        phase = F(sign*int(extent)**2*i*j, 100)
        design = witness_design(vector, phase)
        means = (*axes[i, 0], *axes[0, j], *interior[i, j])
        target = design.value(means)
        rows[extent] = {'coordinate_extent': int(extent), 'point': [i, j], 'product_order_sign': sign,
            'design': design.specification(), 'public_plugin_score': str(target),
            'public_plugin_score_float': float(target),
            'future_zero_service_crossing_trials': crossing_budget(design, target, service_errors=(F(0),)*6),
            'failure_probability_per_extent': '1/200',
            'family_failure_probability_for_two_extents': '1/100',
            'future_mean_must_reach_declared_target': True,
            'power_guaranteed': False, 'public_data_statistical_rejection': False,
            'physical_service_errors_certified': False}
    return rows
