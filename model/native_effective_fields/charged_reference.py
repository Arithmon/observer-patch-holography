"""Exact relational readout of a native charge-one coherence.

The reference preparation and final exchange analyzer are declared operations.
Native C/E current instruments are the only intervening dynamics; their
uniform-seam schedule is a declared protocol, not an observer clock.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from itertools import product

from ..federated_quantum_sector import construction as quantum


CURRENT_VALUES = (-1, 0, 1)


def _initial(port):
    """(|one_port,0_ref>+|vacuum,1_ref>)/sqrt(2), at total number one."""
    one = 1 << port
    return {(one, 0, one, 0): F(1, 2),
            (0, 1, 0, 1): F(1, 2),
            (one, 0, 0, 1): F(1, 2),
            (0, 1, one, 0): F(1, 2)}


def _step(density, law, current, seams, *, parity_twist=False):
    """One exact current outcome; an optional charge parity follows repair."""
    out = defaultdict(F)
    for (left, left_ref, right, right_ref), amplitude in density.items():
        phase = -1 if parity_twist and (left.bit_count()-right.bit_count()) % 2 else 1
        for seam in seams:
            for (a, b), value in quantum.native_unit(law, left, right, seam, current).items():
                if a.bit_count() != left.bit_count() or b.bit_count() != right.bit_count():
                    raise ArithmeticError('native charged block changed occupation number')
                out[a, left_ref, b, right_ref] += amplitude*phase*value/len(seams)
    return {key: value for key, value in out.items() if value}


def _trace(density):
    return sum((value for (left, left_ref, right, right_ref), value in density.items()
                if left == right and left_ref == right_ref), F(0))


def _plus_effect(density, port):
    """The number-conserving (I+X_port,ref)/2 effect, fair off its pair."""
    one = 1 << port
    return (_trace(density) + density.get((one, 0, 0, 1), F(0))
            + density.get((0, 1, one, 0), F(0)))/2


def _word_outputs(law, seams, port, attempts, *, parity_twist=False):
    outputs = {(): _initial(port)}
    for _ in range(attempts):
        outputs = {word + (current,): _step(density, law, current, seams,
                     parity_twist=parity_twist)
                   for word, density in outputs.items() for current in CURRENT_VALUES}
    return outputs


def reference_readout_scope(seams):
    """Bind exact native current laws and a declared relational phase effect.

    A charge-parity phase after each step is an alternative quantum completion
    of the same classical endpoint/current law, with the same C or E record
    grouping. It is not a step in the chosen C or E quantum instrument. Its
    freedom shows that a physical oscillation frequency cannot be inferred
    from the classical endpoint and current rules alone.
    """
    seams = tuple(seams)
    if len(seams) != 30 or len(set(seams)) != 30:
        raise ValueError('the native thirty-seam carrier is required')
    ports = {vertex for seam in seams for vertex in seam}
    if ports != set(range(12)):
        raise ValueError('the native twelve-port carrier is required')

    # Check the exact local Kraus action on all matrix units and all current
    # outcomes. Charge parity changes only cross-number coherences.
    for law, left, right, current in product(('C', 'E'), range(4), range(4), CURRENT_VALUES):
        for (a, b) in quantum.native_unit(law, left, right, (0, 1), current):
            if a.bit_count() != left.bit_count() or b.bit_count() != right.bit_count():
                raise ArithmeticError('an endpoint Kraus branch changed number')

    rates = {}
    readouts = {}
    neutral_words = None
    for law in ('C', 'E'):
        one_step = _word_outputs(law, seams, 0, 1)
        coherence = sum((density.get((1, 0, 0, 1), F(0)) for density in one_step.values()), F(0))
        rate = 2*coherence
        expected = F(11, 12) if law == 'C' else F(5, 6)
        if rate != expected:
            raise ArithmeticError('the relational coherence disagrees with native Kraus counts')
        rates[law] = str(rate)
        rows = {}
        for attempts in (1, 2):
            direct = one_step if attempts == 1 else _word_outputs(law, seams, 0, attempts)
            twisted = _word_outputs(law, seams, 0, attempts, parity_twist=True)
            direct_words = {word: _trace(density) for word, density in direct.items()}
            twisted_words = {word: _trace(density) for word, density in twisted.items()}
            if direct_words != twisted_words or sum(direct_words.values(), F(0)) != 1:
                raise ArithmeticError('charge parity changed a native current word')
            if neutral_words is None:
                neutral_words = {attempts: direct_words}
            elif attempts not in neutral_words:
                neutral_words[attempts] = direct_words
            elif direct_words != neutral_words[attempts]:
                raise ArithmeticError('C/E current laws differed on the same preparation')
            plus = sum((_plus_effect(density, 0) for density in direct.values()), F(0))
            twisted_plus = sum((_plus_effect(density, 0) for density in twisted.values()), F(0))
            if plus != (1+rate**attempts)/2 or twisted_plus != (1+(-rate)**attempts)/2:
                raise ArithmeticError('native reference readout lost its exact charged rate')
            rows[str(attempts)] = {'plus': str(plus), 'parity_twisted_plus': str(twisted_plus),
                                   'current_word_count': len(direct_words)}
        readouts[law] = rows

    return {
        'preparation': 'declared total-one source/reference superposition',
        'analyzer': 'declared total-number-conserving pair-exchange effect (I+X)/2',
        'parameter': 'independent uniform native seam attempts',
        'source_charge_one_coherence_rate': rates,
        'complete_current_words_C_E_and_parity_twist_identical': True,
        'reference_plus_probabilities': readouts,
        'all_source_charge_U1_phase_twists_leave_current_words_identical': True,
        'relative_phase_plus_law': '(1+r_L^n cos(n theta))/2; U_theta acts on source only',
        'charge_parity_twist_is_classically_indistinguishable': True,
        'phase_frequency_selected_by_native_current_records': False,
        'reference_preparation_or_analyzer_implemented_by_current_observer': False,
        'physical_charged_spectral_pole_or_rest_energy_identified': False,
    }


def charged_spatial_pole_gate(seams):
    """Exact source C/E charge-one two-point block versus token transport.

    Every visible-current branch is checked separately. The local result then
    extends to any finite, possibly adaptive, sequence of native seam choices:
    a scalar multiple of |one_i><vac| never becomes |one_j><vac| for j != i.
    The uniform-seam rates and population comparison are separate consequences.
    """
    seams = tuple(seams)
    ports = {port for seam in seams for port in seam}
    if (len(seams) != 30 or len(set(seams)) != 30 or ports != set(range(12))
            or any(type(seam) is not tuple or len(seam) != 2
                   or any(type(port) is not int for port in seam) for seam in seams)):
        raise ValueError('the native thirty-seam twelve-port carrier is required')
    degrees = {port: sum(port in seam for seam in seams) for port in ports}
    if set(degrees.values()) != {5}:
        raise ArithmeticError('the native five-regular carrier changed')

    signed_difference_checks = 0
    for law, seam in product(('C', 'E'), seams):
        a, b = seam
        for left_local, right_local, current in product(range(4), range(4), CURRENT_VALUES):
            left = ((left_local & 1) << a) | (((left_local >> 1) & 1) << b)
            right = ((right_local & 1) << a) | (((right_local >> 1) & 1) << b)
            signed = tuple(((left >> port) & 1)-((right >> port) & 1) for port in seam)
            for output_left, output_right in quantum.native_unit(law, left, right, seam, current):
                if tuple(((output_left >> port) & 1)-((output_right >> port) & 1)
                         for port in seam) != signed:
                    raise ArithmeticError('a native current branch transported signed charge')
            signed_difference_checks += 1

    rates = {}
    checked = 0
    for law in ('C', 'E'):
        for port in sorted(ports):
            for seam in seams:
                incident = port in seam
                scalar = (F(1, 2) if law == 'C' else F(0)) if incident else F(1)
                for current in CURRENT_VALUES:
                    expected = {(1 << port, 0): scalar} if scalar and current == 0 else {}
                    actual = quantum.native_unit(law, 1 << port, 0, seam, current)
                    if actual != expected:
                        raise ArithmeticError('a charged current branch changed port or coefficient')
                    checked += 1
        rate = F(1)-F(5, 2*len(seams)) if law == 'C' else F(1)-F(5, len(seams))
        rates[law] = str(rate)

    # The diagonal one-token sector is a different block of the same exact
    # native instrument. It moves with probability 1/60 across each seam.
    population_edges = 0
    for source in sorted(ports):
        column = defaultdict(F)
        for seam in seams:
            for (left, right), weight in quantum.native_unit('C', 1 << source, 1 << source, seam).items():
                if left != right or left.bit_count() != 1:
                    raise ArithmeticError('one-token population left its sector')
                column[left.bit_length()-1] += weight/len(seams)
        expected = {port: F(1, 60) for port in ports if port != source and
                    any({source, port} == set(seam) for seam in seams)}
        expected[source] = F(11, 12)
        if dict(column) != expected or sum(column.values(), F(0)) != 1:
            raise ArithmeticError('native one-token population transport changed')
        population_edges += len(expected)-1

    return {
        'all_current_resolved_charge_one_columns_checked': checked,
        'all_current_resolved_endpoint_signed_difference_columns_checked': signed_difference_checks,
        'portwise_signed_difference_preserved_under_native_current_words': True,
        'any_current_word_or_adaptive_seam_schedule_preserves_charge_one_port': True,
        'charged_two_point_all_lags': 'G_ij(n)=delta_ij r_L^n; n is uniform native attempts',
        'charged_two_point_generating_function': 'G_ij(z)=delta_ij/(1-z r_L)',
        'charged_two_point_rates_C_E': rates,
        'scalar_temporal_Hankel_positive_form': 'sum_i |sum_n a_(i,n) r_L^n|^2 >= 0',
        'full_Osterwalder_Schrader_axioms_verified': False,
        'spatial_band_in_graph_laplacian_basis': 'flat r_L on all 12 modes',
        'one_token_population_directed_neighbor_entries': population_edges,
        'one_token_population_transition_per_neighbor': '1/60',
        'positive_speed_particle_dispersion_from_this_charge_one_block': False,
        'native_attempt_count_is_calibrated_physical_time': False,
        'charged_field_quantum_numbers_or_CAR_selected': False,
    }
