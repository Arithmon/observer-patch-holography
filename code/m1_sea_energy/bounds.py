"""Integer resource and protection schedules for the long sea preparation."""

import math


def candidate():
    schedules = []
    for power in (1, 3, 8, 20, 40):
        ell = power+1
        k = (64*ell-1).bit_length()
        m = math.isqrt(4096*ell-1)+1
        schedules.append(dict(q=2**power, ell=ell, levels=k, word_side=m, word_bits=m*m,
                              quantum_power=2**k))
    return dict(mode_exponents=[3, 3, 0], raw_rotations=[9, 9, 0], program=[9, 16, 0],
                inventory=[9, 18, 1], depth=[9, 18, 1], active=[18, 36, 2],
                diagnostic=[18, 37, 2], storage_time=[27, 55, 3],
                quantum_decay=128, archive_decay=128, synthesis_decay=64,
                joint_error=[-55, 36, 2], energy_range=[4, 3, 0], energy_error=[-51, 39, 2],
                schedules=schedules,
                routing=[dict(modes=n, phases=n, mixers=n*(n-1)//2,
                              swaps_upper=2*max(0, n-2)*(n*(n-1)//2),
                              rotations_upper=n+6*n*(n-1)//2+16*max(0, n-2)*n*(n-1)//2)
                         for n in (2, 3, 4, 16, 128)])
