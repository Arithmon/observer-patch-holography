"""Strict evidence tubes for the nine native-instrument quadratures.

This module compares a locked list of admissible tables with exact probability
tubes. It never projects noisy means onto the CP family or declares a physical
law selected. An empty compatible locked list does not rule out the family.

Counts schema (all nine settings, in core.PROBES order)::

    {"schema": "oph.repair_identification.counts.v1", "settings": [
        {"probe": "k_re", "attempted": 100,
         "counts": {"current_-1_plus": 0, "current_-1_minus": 0,
                    "current_0_plus": 75, "current_0_minus": 20,
                    "current_1_plus": 0, "current_1_minus": 0,
                    "abort": 5}}, ...]}

Every attempted trial appears in exactly one of the seven bins. In particular
the estimate above is 75/100, not 75/95. The omitted rows in this documentation
example must be supplied; the parser rejects incomplete menus and totals.

analyze_counts requires a separate sampling design and service-error bounds.
The fixed attempt counts and sampling radii must be declared before, and
independently of, all analyzed trial outcomes. Checking the design against
the counts establishes consistency, not that this premise was obeyed.
For independent attempted trials within each setting, conditional on this
design, Var(mean)<=1/(4m).
Chebyshev and a union bound give simultaneous failure probability at most
min(1, sum_i 1/(4 m_i radius_i**2)). Cross-setting independence is unnecessary.
The service bounds must relate the average attempted-trial response to the
ideal table probabilities. They are supplied assumptions, not calibration
established by this parser. Arbitrary correlated trial memory is not covered.
"""
from fractions import Fraction as F
from typing import NamedTuple

from . import core


SCHEMA = 'oph.repair_identification.counts.v1'
OUTCOMES = ('current_-1_plus', 'current_-1_minus',
            'current_0_plus', 'current_0_minus',
            'current_1_plus', 'current_1_minus', 'abort')
SOURCE_PATHS = (*core.SOURCE_PATHS, 'model/repair_identification/calibration.py')


class SamplingDesign(NamedTuple):
    """Immutable declared horizons and positive rational radii (integer pairs)."""

    planned_attempts: tuple
    sampling_radii: tuple


def _fields(value, expected):
    if (type(value) is not dict or any(type(key) is not str for key in value)
            or set(value) != set(expected)):
        raise ValueError('the complete original evidence object and its exact fields are required')


def _integer(value, *, positive=False):
    if type(value) is not int or value < int(positive):
        raise ValueError('original nonnegative counts and positive attempted totals are required')
    return value


def _probability(value):
    result = core._fraction(value)
    if not 0 <= result <= 1:
        raise ValueError('an exact rational probability or error in [0,1] is required')
    return result


def _nine(values, *, positive=False):
    if type(values) is not tuple or len(values) != len(core.PROBES):
        raise ValueError('an original tuple of nine exact rational values is required')
    result = tuple(_probability(value) for value in values)
    if positive and any(value == 0 for value in result):
        raise ValueError('sampling radii must be strictly positive')
    return result


def sampling_design(planned_attempts, sampling_radii):
    """Detach a declared design; this object does not prove its declaration time."""
    if type(planned_attempts) is not tuple or len(planned_attempts) != len(core.PROBES):
        raise ValueError('an original tuple of nine planned attempt counts is required')
    attempted = tuple(_integer(value, positive=True) for value in planned_attempts)
    radii = _nine(sampling_radii, positive=True)
    return SamplingDesign(attempted, tuple((value.numerator, value.denominator) for value in radii))


def _validate_design(design):
    if type(design) is not SamplingDesign:
        raise ValueError('a separate original SamplingDesign is required')
    if type(design.sampling_radii) is not tuple or len(design.sampling_radii) != len(core.PROBES):
        raise ValueError('nine original rational radius pairs are required')
    radii = []
    for pair in design.sampling_radii:
        if type(pair) is not tuple or len(pair) != 2:
            raise ValueError('original numerator/denominator pairs are required')
        numerator, denominator = (_integer(value, positive=True) for value in pair)
        radius = F(numerator, denominator)
        if (radius.numerator, radius.denominator) != pair:
            raise ValueError('canonical rational radius pairs are required')
        radii.append(radius)
    validated = sampling_design(design.planned_attempts, tuple(radii))
    return validated, tuple(radii)


def validate_counts(evidence):
    """Return detached integer rows (attempted, seven-bin tuple)."""
    _fields(evidence, ('schema', 'settings'))
    if type(evidence['schema']) is not str or evidence['schema'] != SCHEMA:
        raise ValueError('the complete repair-identification counts schema is required')
    settings = evidence['settings']
    if type(settings) is not list or len(settings) != len(core.PROBES):
        raise ValueError('all nine ordered probe settings are required')
    result = []
    for probe, row in zip(core.PROBES, settings, strict=True):
        _fields(row, ('probe', 'attempted', 'counts'))
        if type(row['probe']) is not str or row['probe'] != probe.name:
            raise ValueError('every original probe must occur once in canonical order')
        attempted = _integer(row['attempted'], positive=True)
        _fields(row['counts'], OUTCOMES)
        counts = tuple(_integer(row['counts'][label]) for label in OUTCOMES)
        if sum(counts) != attempted:
            raise ValueError('all currents, readouts and aborts must sum to attempted trials')
        result.append((attempted, counts))
    return tuple(result)


def analyze_tubes(centers, radii, candidates, call_count):
    """Compare locked candidates without deciding feasibility of the whole family.

    Between any two admissible tables in these tubes, the nine probability
    differences are at most 2 max(radii). Their common-current instrument
    distance is at most 14 max(radii), and N calls give at most 14N max(radii).
    These diameter bounds remain valid if no locked candidate lies in the
    tubes. No assertion that the tubes contain any admissible table follows.
    """
    centers, radii = _nine(centers), _nine(radii)
    core.history_bound(call_count, F(0))  # validates the original integer cap
    if (type(candidates) is not dict or not candidates
            or any(type(name) is not str or not name for name in candidates)):
        raise ValueError('a nonempty locked mapping of original candidate names to tables is required')
    tables = {name: core.validate_table(table) for name, table in candidates.items()}
    intervals = tuple((max(F(0), center-radius), min(F(1), center+radius))
                      for center, radius in zip(centers, radii, strict=True))
    rows, compatible, excluded = [], [], []
    for name, table in tables.items():
        predicted = core.predictions(table)
        violations = [probe.name for probe, value, (lower, upper) in
                      zip(core.PROBES, predicted, intervals, strict=True)
                      if (value-lower).sign() < 0 or (value-upper).sign() > 0]
        (excluded if violations else compatible).append(name)
        rows.append({'name': name, 'compatible': not violations,
                     'violating_probes': violations,
                     'predicted_joint_probabilities': [str(value) for value in predicted]})
    maximum_radius = max(radii)
    instrument_diameter = core.coordinate_error_bound(min(F(1), 2*maximum_radius))
    return {
        'schema': 'oph.repair_identification.probability_tubes.v1',
        'probe_order': [probe.name for probe in core.PROBES],
        'centers': [str(value) for value in centers],
        'radii': [str(value) for value in radii],
        'probability_intervals': [[str(lower), str(upper)] for lower, upper in intervals],
        'locked_candidate_checks': rows,
        'compatible_locked_candidates': compatible,
        'excluded_locked_candidates': excluded,
        'candidate_set_exhaustive': False,
        'full_family_feasibility': ('witnessed_by_locked_candidate' if compatible else 'not_decided'),
        'cp_projection_performed': False,
        'physical_law_selected': False,
        'call_cap': call_count,
        'maximum_probability_radius': str(maximum_radius),
        'common_current_instrument_diameter_bound': str(instrument_diameter),
        'complete_history_diameter_bound': str(core.history_bound(call_count, instrument_diameter)),
        'diameter_scope': 'any two admissible tables in the tubes; identical common-current controls, '
                          'all outcomes, arbitrary untouched references; no fine-record comparison',
    }


def analyze_counts(evidence, design, candidates, call_count, service_errors):
    """Conditional finite-sample evidence, with no successful-trial conditioning.

    All returned exact numerical quantities are rational strings; counts and
    schema objects remain ordinary JSON values. The separate design must have
    outcome-independent fixed horizons and radii. This parser checks only its
    consistency with the observed counts, not its prior declaration. Supplied
    service errors are mandatory, including an explicit zero-error assumption.
    """
    counts = validate_counts(evidence)
    design, sampling_radii = _validate_design(design)
    if tuple(attempted for attempted, _ in counts) != design.planned_attempts:
        raise ValueError('observed attempted counts must match the declared sampling design')
    service_errors = _nine(service_errors)
    centers, targets = [], []
    for probe, (attempted, bins) in zip(core.PROBES, counts, strict=True):
        target = bins[OUTCOMES.index(f'current_{probe.current}_plus')]
        centers.append(F(target, attempted))
        targets.append(target)
    total_radii = tuple(min(F(1), sampling+service) for sampling, service in
                        zip(sampling_radii, service_errors, strict=True))
    result = analyze_tubes(tuple(centers), total_radii, candidates, call_count)
    failure = min(F(1), sum((F(1, 4*attempted)/radius**2
                            for (attempted, _), radius in zip(counts, sampling_radii, strict=True)), F(0)))
    result['count_evidence'] = {
        'schema': SCHEMA,
        'complete_outcome_order': list(OUTCOMES),
        'attempted_by_probe': [attempted for attempted, _ in counts],
        'complete_counts_by_probe': [list(bins) for _, bins in counts],
        'target_joint_counts': targets,
        'abort_counts': [bins[-1] for _, bins in counts],
        'sampling_radii': [str(value) for value in sampling_radii],
        'declared_sampling_design': {
            'planned_attempts': list(design.planned_attempts),
            'sampling_radii': [str(value) for value in sampling_radii],
        },
        'sampling_design_consistent_with_counts': True,
        'sampling_design_prior_declaration_verified': False,
        'supplied_service_errors': [str(value) for value in service_errors],
        'simultaneous_sampling_failure_bound': str(failure),
        'simultaneous_sampling_coverage_lower_bound': str(1-failure),
        'sampling_method': 'exact rational Chebyshev variance bound and union bound',
        'premises': [
            'fixed horizons and sampling radii declared before, and independently of, all analyzed trial outcomes',
            'independent attempted trials within each setting conditional on that design; no cross-setting independence needed',
            'average response differs from the ideal table probability by at most its supplied service error',
            'every attempt, including aborts and other current/readout outcomes, is counted',
        ],
        'premises_established_by_counts': False,
    }
    return result
