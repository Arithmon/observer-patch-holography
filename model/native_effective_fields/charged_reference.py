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
