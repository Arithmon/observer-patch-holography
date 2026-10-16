"""A fixed-trial C/E decision from complete native phase-readout records."""
from __future__ import annotations

from fractions import Fraction as F

import numpy as np

from ..calibrated_measurement.bounds import log_interval, root_interval
from ..compiled_native_interface import construction as compiled
from ..native_interface import construction as interface
from ..native_quantum_history import construction as quantum
from ..primitive_observer import runtime


def exact_law(candidate, attempts=4):
    candidate, attempts = quantum.candidate(candidate), quantum.attempts(attempts)
    rate = quantum.RATES[candidate]
    # Derive the mode multiplier from actual source rows, rather than only
    # accepting the retained rate constant. The orientation coherence mode
    # is the sum of the ten pair matrix units with equal coefficients.
    left, right = interface.phase_pairs()[0]
    difference = tuple(a-b for a, b in zip(quantum.STATES[left], quantum.STATES[right], strict=True))
    block = quantum.mixing.native_block(candidate, difference, binary=True)
    if len(block) != 10 or any(sum(row, F(0)) != rate for row in block):
        raise ArithmeticError('the actual source coherence multiplier changed')
    p = (1+rate**attempts)/2
    success = 2*interface.SUCCESS
    return (1-success, F(0), success*p/2, success*(1-p)/2,
            success*p/2, success*(1-p)/2)


def _law(values):
    if (type(values) is not tuple or len(values) != 6 or any(type(x) is not F or x < 0 for x in values)
            or sum(values) != 1):
        raise ValueError('a complete six-outcome exact probability tuple is required')
    return values


def trial_design(first, second, error=F(1, 100)):
    """Each candidate's likelihood-ratio error is at most affinity**N.

    Under P, Markov's inequality on sqrt(Q/P) bounds choosing Q, and
    conversely under Q. Product measures factor the affinity. Ties may
    deterministically choose P; the two weak-inequality bounds still hold.
    """
    first, second = _law(first), _law(second)
    if type(error) is not F or not 0 < error < 1:
        raise ValueError('an exact failure allowance in (0,1) is required')
    intervals = [root_interval(a*b) for a, b in zip(first, second, strict=True)]
    lower, upper = (sum((row[i] for row in intervals), F(0)) for i in (0, 1))
    if not 0 < upper < 1:
        raise ValueError('strictly distinct overlapping candidate laws are required')
    # Round the already outward exact logarithms out again for compact,
    # platform-independent receipts, retaining the decision inequality.
    grid = 10**18
    upper_log = log_interval(1/error)[1]*grid
    lower_log = log_interval(1/upper)[0]*grid
    numerator = F(-(-upper_log.numerator//upper_log.denominator), grid)
    denominator = F(lower_log.numerator//lower_log.denominator, grid)
    if denominator <= 0:
        raise ArithmeticError('the outward affinity did not resolve distinguishability')
    ratio = numerator/denominator
    trials = -(-ratio.numerator//ratio.denominator)
    assert trials*denominator >= numerator
    return {'affinity_lower': str(lower), 'affinity_upper': str(upper),
            'fixed_total_trials': trials, 'individual_error_upper': str(error),
            'log_inverse_affinity_lower': str(denominator),
            'log_inverse_error_upper': str(numerator),
            'bound': 'P_C(choose E) <= affinity^N and P_E(choose C) <= affinity^N',
            'all_failed_heralds_count_toward_trials': True,
            'physical_service_or_preparation_uncertainty_included': False}


def classify(counts, first, second, trials):
    """Exact future decision; never treat a zero-probability record as evidence."""
    first, second = _law(first), _law(second)
    if (type(trials) is not int or not 1 <= trials <= 100000
            or type(counts) is not tuple or len(counts) != 6
            or any(type(n) is not int or n < 0 for n in counts) or sum(counts) != trials):
        raise ValueError('all six nonnegative counts and the fixed trial total are required')
    if any(n and p == q == 0 for n, p, q in zip(counts, first, second, strict=True)):
        return 'outside_both_declared_models'
    p = q = F(1)
    for n, a, b in zip(counts, first, second, strict=True):
        if n:
            p *= a**n
            q *= b**n
    if p == q == 0:
        return 'outside_both_declared_models'
    return 'C' if p >= q else 'E'


def execute(attempts=4):
    # Execute the compact source word and check the full spectator lift;
    # never substitute a compressed raw matrix for its complete 512 map.
    local = compiled.local_matrix()
    full = compiled.full_matrix()
    if not np.array_equal(full, np.kron(np.eye(64), local)):
        raise ArithmeticError('the complete source tensor lift changed')
    embedding = compiled.embedding()
    target = np.kron(np.eye(64), compiled.ideal_local_matrix())
    if np.linalg.norm(target@embedding-embedding@compiled.ideal_native_matrix(), 'fro') > 1e-14:
        raise ArithmeticError('the full ideal source target no longer transports the native mixer')
    compiled.exact_effect_audit()
    bound = compiled.certify(local)
    if bound['two_call_half_diamond_upper'] > F(1, 10**10):
        raise ArithmeticError('the source mixer exceeded its fixed two-call envelope')
    source = {'full_dimension': 512, 'local_dimension': 8,
              'all_padding_coordinates_retained': True,
              'native_phase_pair_encoding_and_target_transport_checked': True,
              'two_mixer_half_diamond_upper': '1/10000000000',
              'resources_per_mixer': compiled.resources(),
              'word': compiled.word(),
              'scope': 'exact application of original binary mixer matrices versus exact target, with norm amplification',
              'coherent_proper_code_services_remain_supplied': True}
    rows = {}
    for candidate in ('C', 'E'):
        exact = exact_law(candidate, attempts)
        interface.exact_native_audit(candidate)
        program = interface.program(candidate, attempts)
        history = runtime.execute(program)
        expected = dict(zip(((0,0),(0,1),(1,0),(1,1),(2,0),(2,1)), exact, strict=True))
        actual = {branch.record: branch.probability for branch in history.branches}
        if actual.keys() != expected.keys():
            raise ArithmeticError('complete six-record law changed')
        residual = max(abs(actual[key]-float(value)) for key, value in expected.items())
        if residual > 1e-10:
            raise ArithmeticError('complete native history disagrees with its exact source law')
        checkpoint = runtime.execute(program, stop_after=2)
        restored = runtime.checkpoint_from_json(program, runtime.checkpoint_to_json(checkpoint))
        continued = runtime.execute(program, checkpoint=restored)
        distance = runtime.complete_distance(history, continued)
        if distance > 1e-10:
            raise ArithmeticError('the native experiment failed serialized continuation')
        rows[candidate] = {'exact_complete_law': list(map(str, exact)),
            'executed_complete_law': [actual[key] for key in expected],
            'runtime_probability_residual_diagnostic': residual,
            'checkpoint_complete_density_distance_diagnostic': distance,
            'conditional_match_probability': str((1+quantum.RATES[candidate]**attempts)/2),
            'runtime_is_ideal_66_sector_surrogate': True}
    first, second = exact_law('C', attempts), exact_law('E', attempts)
    design = trial_design(first, second)
    return {'native_attempts': attempts, 'fresh_source_mixer': source, 'candidates': rows,
            'complete_law_total_variation': str(sum((abs(a-b) for a, b in zip(first,second,strict=True)),F(0))/2),
            'trial_design': design, 'new_empirical_outcomes': False,
            'raw_source_maps_are_not_assumed_exactly_trace_preserving': True,
            'statistical_guarantee_applies_to_ideal_candidate_laws_only': True,
            'padding_extension_for_source_comparison': 'dephase code/complement; exact native maps on code and identity on complement with failure record',
            'hardware_bridge_authenticated': False}
