"""Independent integer inequalities, including the nonvanishing-error regime."""

from .format import integer, keys, need


def verify(row):
    keys(row, 'logical_locations accounting_q_degree synthesis_per_gate_decay state_decay accounting_decay '
         'fixed_rate growing resources resource_convention contraction_convention')
    need(type(row['logical_locations']) is list and row['logical_locations'] == [4, 6]
         and all(type(x) is int for x in row['logical_locations']), 'proved parent program bound')
    integer(row['accounting_q_degree'], 4, 4)
    integer(row['synthesis_per_gate_decay'], 16, 16)
    integer(row['state_decay'], 16-4, 16-4)
    integer(row['accounting_decay'], 16-4-4, 16-4-4)
    need(type(row['fixed_rate']) is list and len(row['fixed_rate']) == 2, 'both fixed-rate noise classes')
    for item, kind, denominator, selected in zip(row['fixed_rate'], ('markov', 'common_bath'), (2, 1), (5, 4)):
        keys(item, 'kind amplitude_decay_denominator selected_level levels')
        need(item['kind'] == kind, 'noise class')
        integer(item['amplitude_decay_denominator'], denominator, denominator)
        integer(item['selected_level'], selected, selected)
        need(type(item['levels']) is list and len(item['levels']) == 6, 'all compared fixed levels')
        power = 1
        for k, record in enumerate(item['levels'], 1):
            power *= 2
            keys(record, 'level fault_power state_q_degree accounting_q_degree')
            for key, value in dict(level=k, fault_power=power, state_q_degree=4-power//denominator,
                                   accounting_q_degree=8-power//denominator).items():
                integer(record[key], value, value)
        need(item['levels'][selected-1]['state_q_degree'] == -12
             and item['levels'][selected-1]['accounting_q_degree'] == -8
             and item['levels'][selected-2]['accounting_q_degree'] == 0,
             'required protection margin includes accounting')
    catalog = [(q, b) for q in (2, 5, 16, 63, 256, 65536, 2**32, 2**100) for b in (2, 3, 4)]
    need(type(row['growing']) is list and len(row['growing']) == len(catalog), 'complete growing-level catalog')
    for item, (q, bits) in zip(row['growing'], catalog):
        keys(item, 'q contraction_bits levels suppression_bits required_bits')
        integer(item['q'], q, q)
        integer(item['contraction_bits'], bits, bits)
        integer(item['levels'], 1, 32)
        # Calculate ceil(log2(2q)) using doubling, independently of bit_length.
        ceiling, value = 0, 1
        while value < 2*q:
            value *= 2
            ceiling += 1
        required = 16*ceiling
        power = bits*(1 << item['levels'])
        integer(item['suppression_bits'], power, power)
        integer(item['required_bits'], required, required)
        need(power >= required and (item['levels'] == 1 or power//2 < required),
             'minimal sufficient growing level')
        # This exact integer comparison certifies theta^(2^k)<=(2q)^-16.
        need((1 << power) >= (2*q)**16, 'nonvanishing-noise operational suppression')
    expected = dict(active_inventory=[3, 0, 1], depth=[1, 6, 1], active_volume=[4, 6, 2],
                    diagnostic_bits=[4, 7, 2], physical_storage_volume=[5, 13, 3])
    keys(row['resources'], 'active_inventory depth active_volume diagnostic_bits physical_storage_volume')
    for name, expected_value in expected.items():
        actual = row['resources'][name]
        need(type(actual) is list and all(type(x) is int for x in actual) and actual == expected_value,
             'charged lifetime resource '+name)
    need([a+b for a, b in zip(expected['active_inventory'], expected['depth'])] == expected['active_volume'],
         'active inventory times lifetime')
    need([a+b for a, b in zip(expected['diagnostic_bits'], expected['depth'])] == expected['physical_storage_volume'],
         'diagnostic inventory times lifetime')
    need(row['resource_convention'] == '[q power, constant log power, coefficient of kappa]; fixed levels set kappa=0',
         'resource exponent meaning')
    need(row['contraction_convention'] == 'theta=A eta is illustrative; no numeric native threshold',
         'threshold boundary')
