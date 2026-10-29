"""Omit a recovery only when its admitted native code image is already correct.

This variant acts on the known flagged image of the declared two-column code.
It does not claim equality of the old and new full-register unitaries on an
arbitrary input. The baseline moving-defect and biased-reset receipts are
unchanged. Attempt, reset, archive and native instrument costs remain intact.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import json

from ..moving_defect import transport as t
from ..biased_reset import reset


@lru_cache(None)
def identity_eligible(center, law, seam, coin, raw):
    """Exact, state-independent flag decision; no logical amplitudes are read."""
    images, target, _ = t.branch(center, law, seam, coin, raw)
    if not images[0] or target not in t.LIVE:
        return False
    desired = tuple(tuple(t.circuit.encode(tuple(int(i in pair) for i in range(12)))
                          for pair in row) for row in t.supports(target))
    return tuple(map(frozenset, images)) == tuple(map(frozenset, desired))


def recover(center, law, seam, coin, raw, state):
    """Retain an eligible unnormalized flag image; otherwise use old recovery.

    The caller supplies the actual native branch of V_center's logical input.
    That input may be entangled with an arbitrary untouched reference.
    """
    value = t.routing._state(state)
    if identity_eligible(center, law, seam, coin, raw):
        return value
    return t.recover(center, law, seam, coin, raw, value)


def recovery_resources(center, law, seam, coin, raw):
    if identity_eligible(center, law, seam, coin, raw):
        return dict.fromkeys(t.routing.COUNT_KEYS, 0)
    return t.recovery_resources(center, law, seam, coin, raw)


def certificate(*, resources=True):
    """All flags checked exactly; optional rational stationary service savings.

    The stationary distribution belongs to the existing biased boundary-reset
    process. This changes no core transition or boundary reset, so it is also
    the stationary distribution of this service variant.
    """
    census = Counter()
    saved = F()
    maximum = 0
    stationary = (tuple(map(F, reset.renewal_certificate()['stationary_core_probabilities']))
                  if resources else None)
    for i, center in enumerate(t.LIVE):
        for law in ('C', 'E'):
            for seam in t.SEAMS:
                for coin in (0, 1):
                    for raw in range(1 << t.circuit.RECORD_BITS[law]):
                        label = center, law, seam, coin, raw
                        images, target, abort = t.branch(*label)
                        census['all_raw_slots'] += 1
                        if not images[0]:
                            census['raw_abort' if abort else 'valid_zero'] += 1
                            continue
                        census['positive_flags'] += 1
                        if not identity_eligible(*label):
                            continue
                        census['identity_recovery_flags'] += 1
                        census['identity_flags_with_core_motion'] += target != center
                        # The raw map is a positive partial injection. Both
                        # logical columns have equal, disjoint supports. Exact
                        # target support equality therefore implies equality of
                        # amplitudes with sqrt(p_flag) V_target, on every column.
                        if len(images[0]) != len(images[1]):
                            raise ArithmeticError('identity flag lost scalar probability')
                        if resources:
                            cost = t.recovery_resources(*label)['native_pulses']
                            saved += stationary[i]*t.branch_probability(*label)*cost/60
                            maximum = max(maximum, cost)
    result = {'schema': 'oph.recovery_elision.v1', 'census': dict(census),
              'identity': 'K_flag V_core = sqrt(p_flag) V_target on eligible flags',
              'unknown_logical_reference_preserved': True,
              'old_full_register_unitary_equality_claimed': False,
              'native_instrument_core_kernel_resets_and_archives_unchanged': True,
              'baseline_receipts_rewritten': False,
              'physical_clock_energy_or_hardware_error_certificate': False}
    if resources:
        result['stationary_saved_source_pulses_per_native_attempt'] = str(saved)
        result['largest_omitted_recovery_pulses'] = maximum
        result['resource_scope'] = 'supplied ideal service words; equal C/E and uniform seams in the biased reset process'
    return result


if __name__ == '__main__':
    print(json.dumps(certificate(), sort_keys=True, indent=2))
