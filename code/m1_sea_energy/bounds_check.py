"""Independent integer ceilings and lifetime algebra, with no fit exponents."""

from .format import keys, need, exact


def verify(row):
    keys(row, 'mode_exponents raw_rotations program inventory depth active diagnostic storage_time quantum_decay archive_decay synthesis_decay joint_error energy_range energy_error schedules routing')
    modes = [3, 3, 0]
    raw = [3*x for x in modes]
    program = [raw[0], raw[1]+7, 0]
    width = [program[0], program[1]+2, 1]
    depth = width.copy()
    active = [a+b for a, b in zip(width, depth)]
    diagnostic = [active[0], active[1]+1, active[2]]
    lifetime = [a+b for a, b in zip(diagnostic, depth)]
    energy = [modes[0]+1, modes[1], modes[2]]
    exponent = min(128-active[0], 128-active[0], 64-program[0])
    joint = [-exponent, active[1], active[2]]
    for key, value in dict(mode_exponents=modes, raw_rotations=raw, program=program,
                           inventory=width, depth=depth, active=active, diagnostic=diagnostic,
                           storage_time=lifetime, quantum_decay=128, archive_decay=128,
                           synthesis_decay=64, joint_error=joint, energy_range=energy,
                           energy_error=[a+b for a, b in zip(joint, energy)]).items():
        exact(row[key], value)
    need(type(row['schedules']) is list and len(row['schedules']) == 5, 'full protection schedule')
    for item, power in zip(row['schedules'], (1, 3, 8, 20, 40)):
        keys(item, 'q ell levels word_side word_bits quantum_power')
        need(all(type(x) is int and 0 < x <= 2**50 for x in item.values()), 'bounded schedule integers')
        exact([item['q'], item['ell']], [2**power, power+1])
        ell, k, m = item['ell'], item['levels'], item['word_side']
        need(2**(k-1) < 64*ell <= 2**k == item['quantum_power'], 'minimal quantum ceiling')
        need((m-1)**2 < 4096*ell <= m*m == item['word_bits'], 'minimal archive ceiling')
        need(2*2**k >= 128*ell and m*m >= 32*128*ell, 'both fixed-strength decays')
    need(type(row['routing']) is list and len(row['routing']) == 5, 'routing sizes')
    for item, n in zip(row['routing'], (2, 3, 4, 16, 128)):
        keys(item, 'modes phases mixers swaps_upper rotations_upper')
        pairs = sum(range(n))
        swaps = pairs*(n-2)*2
        exact(item, dict(modes=n, phases=n, mixers=pairs, swaps_upper=swaps,
                         rotations_upper=n+6*pairs+8*swaps))
        need(item['rotations_upper'] <= 8*n**3, 'cubic native routing bound')
