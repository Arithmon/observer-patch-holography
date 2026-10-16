"""Run complete source-compiled Ramsey instruments for U,V,UV,VU."""
import numpy as np

from ..observer.source_scalar import compiler
from ..primitive_observer import runtime


def source_words():
    theta, eta = np.pi/7, np.pi/5
    u = compiler.pauli_rotation(['Z'], theta)
    v = compiler.pauli_rotation(['X'], eta)
    # Controlled exp(i theta P): exp(i theta P/2) exp(-i theta P Zpointer/2).
    cu = compiler.pauli_rotation(['Z', 'I'], theta/2)+compiler.pauli_rotation(['Z', 'Z'], -theta/2)
    cv = compiler.pauli_rotation(['X', 'I'], eta/2)+compiler.pauli_rotation(['X', 'Z'], -eta/2)
    # Words are chronological: V then U implements the matrix product UV.
    return {'U': (u, cu), 'V': (v, cv), 'UV': (v+u, cv+cu), 'VU': (u+v, cu+cv)}


def ramsey(controlled, phase):
    h = np.array([[1., 1.], [1., -1.]])/np.sqrt(2)
    pointer_h = compiler.remap(compiler.compile_real_orthogonal(h), [1])
    tape = pointer_h+controlled+compiler.compile_diagonal([0., 0., phase, phase])+pointer_h
    matrix = compiler.execute(tape, 2)
    # Pointer is the high bit, initially zero. Keep both outcomes; no projection/renormalization.
    return (matrix[:2, :2], matrix[2:, :2]), compiler.resources(tape)


def execute():
    words = source_words()
    matrices = {name: compiler.execute(word, 1) for name, (word, _) in words.items()}
    rho = np.array([[.6, .2+.1j], [.2-.1j, .4]], complex)
    bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
    reference_initial = np.outer(bell, bell)
    rows = []
    max_residual = max_mass = max_checkpoint = max_unitarity = 0.
    for name, (_, controlled) in words.items():
        target_mean = complex(np.trace(rho@matrices[name]))
        for label, phase, expected in [('real', 0., target_mean.real), ('imag', -np.pi/2, target_mean.imag)]:
            maps, resources = ramsey(controlled, phase)
            instrument = runtime.Instrument(kraus=((maps[0],), (maps[1],)), visible=True,
                resources=(('native_pulses', resources['native_pulses']),),
                binding={'source': compiler.source_summary(), 'word': name, 'quadrature': label})
            program = runtime.Program('source characteristic probe', 2, rho,
                (runtime.Event('spin', instruments=(instrument,)),))
            history = runtime.execute(program)
            probabilities = [b.probability for b in history.branches]
            measured = probabilities[0]-probabilities[1]
            max_residual = max(max_residual, abs(measured-expected))
            ref = runtime.execute(program, initial=reference_initial)
            prefix = runtime.execute(program, initial=reference_initial, stop_after=0)
            resumed = runtime.execute(program, checkpoint=runtime.checkpoint_from_json(
                program, runtime.checkpoint_to_json(prefix)))
            max_checkpoint = max(max_checkpoint, *(float(np.linalg.norm(a.density-b.density, 1))
                for a, b in zip(ref.branches, resumed.branches)))
            max_mass = max(max_mass, abs(history.probability_mass-1), abs(ref.probability_mass-1))
            max_unitarity = max(max_unitarity, float(np.linalg.norm(matrices[name].conj().T@matrices[name]-np.eye(2), 2)))
            rows.append({'word': name, 'quadrature': label, 'spin_probabilities': probabilities,
                         'spin_mean': measured, 'native_pulses': resources['native_pulses']})
    u, v, uv, vu = (complex(np.trace(rho@matrices[name])) for name in ('U', 'V', 'UV', 'VU'))
    radius = float(np.sqrt((1-abs(u)**2)*(1-abs(v)**2)))
    return {'instrument_runs': rows, 'radius': radius,
            'UV_disk_distance': abs(uv-u*v), 'VU_disk_distance': abs(vu-u*v),
            'noncommuting_order_difference': abs(uv-vu),
            'maximum_raw_unitarity_defect': max_unitarity,
            'maximum_complete_mass_defect': max_mass,
            'maximum_source_Ramsey_mean_residual': max_residual,
            'maximum_reference_checkpoint_residual': max_checkpoint,
            'finite_source_code_satisfies_exact_continuous_Weyl_relations': False,
            'physical_coherent_service_error_certified': False}
