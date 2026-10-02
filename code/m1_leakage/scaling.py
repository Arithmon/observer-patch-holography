"""Exact sufficient level choices and charged resource exponents."""

from .format import integer


def levels(q, contraction_bits):
    """theta=2^-contraction_bits; illustrative ratios, not native thresholds."""
    integer(q, 2, 2**4096)
    integer(contraction_bits, 1, 64)
    target = 16*(2*q-1).bit_length()
    level = 1
    while contraction_bits*2**level < target:
        level += 1
    return level


def candidate():
    fixed = []
    for kind, denominator, selected in (('markov', 2, 5), ('common_bath', 1, 4)):
        fixed.append(dict(kind=kind, amplitude_decay_denominator=denominator, selected_level=selected,
                          levels=[dict(level=k, fault_power=2**k, state_q_degree=4-2**k//denominator,
                                       accounting_q_degree=8-2**k//denominator) for k in range(1, 7)]))
    growing = []
    for q in (2, 5, 16, 63, 256, 65536, 2**32, 2**100):
        for bits in (2, 3, 4):
            k = levels(q, bits)
            growing.append(dict(q=q, contraction_bits=bits, levels=k,
                                suppression_bits=bits*2**k, required_bits=16*(2*q-1).bit_length()))
    return dict(logical_locations=[4, 6], accounting_q_degree=4,
                synthesis_per_gate_decay=16, state_decay=12, accounting_decay=8,
                fixed_rate=fixed, growing=growing,
                resources=dict(active_inventory=[3, 0, 1], depth=[1, 6, 1],
                               active_volume=[4, 6, 2], diagnostic_bits=[4, 7, 2],
                               physical_storage_volume=[5, 13, 3]),
                resource_convention='[q power, constant log power, coefficient of kappa]; fixed levels set kappa=0',
                contraction_convention='theta=A eta is illustrative; no numeric native threshold')
