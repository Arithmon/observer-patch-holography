"""Independent finite-CAR and integer-flux check of a dressed quantum link.

The infinite rotor is not replaced by a cyclic register.  A finite invariant
Gauss sector is reconstructed in memory from its fifteen occupied channels;
all 32768 branch amplitudes and exterior-algebra hopping signs are checked.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from verify_current import geometry_check

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / 'quantum_link_receipt.json'
F = Fraction
FIELDS = [('Q', 6, 1), ('u_c', 3, -4), ('d_c', 3, 2), ('L', 2, -3), ('e_c', 1, 6)]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()


def hashed(value):
    return sha256(canonical(value)).hexdigest()


def same(actual, expected, message):
    require(type(actual) is type(expected), message + ': type')
    if isinstance(expected, dict):
        require(actual.keys() == expected.keys(), message + ': keys')
        for key in expected:
            same(actual[key], expected[key], message + '/' + key)
    elif isinstance(expected, (tuple, list)):
        require(len(actual) == len(expected), message + ': length')
        for a, b in zip(actual, expected):
            same(a, b, message)
    else:
        require(actual == expected, message + ': value')


def rational(value):
    require(type(value) is str and len(value) < 2000, 'rational string')
    try:
        result = F(value)
    except (ValueError, ZeroDivisionError):
        raise ValueError('invalid rational') from None
    require(str(result) == value, 'canonical rational')
    return result


def unique(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'duplicate JSON key')
        result[key] = value
    return result


def load(path=OUTPUT):
    path = Path(path)
    require(path.stat().st_size <= 4000000, 'receipt size')
    def reject(_):
        raise ValueError('noninteger JSON number')
    return json.loads(path.read_text(), object_pairs_hook=unique,
                      parse_float=reject, parse_constant=reject)


def annihilate(occupation, mode):
    if not occupation & (1 << mode):
        return None, 0
    sign = -1 if (occupation & ((1 << mode) - 1)).bit_count() % 2 else 1
    return occupation ^ (1 << mode), sign


def create(occupation, mode):
    if occupation & (1 << mode):
        return None, 0
    sign = -1 if (occupation & ((1 << mode) - 1)).bit_count() % 2 else 1
    return occupation | (1 << mode), sign


def hop(occupation, target, source):
    intermediate, a = annihilate(occupation, source)
    if not a:
        return None, 0
    final, b = create(intermediate, target)
    return (final, a * b) if b else (None, 0)


def occupation_bits(branch, count=15):
    return sum(1 << (2 * i + ((branch >> i) & 1)) for i in range(count))


def charges():
    return [charge for _, count, charge in FIELDS for _ in range(count)]


def branch_algebra(quantum_charges=None, flux_sign=-1, modulus=None):
    """Check Gauss and actual CAR signs on the entire invariant branch basis."""
    quantum_charges = charges() if quantum_charges is None else quantum_charges
    require(len(quantum_charges) == 15 and sum(quantum_charges) == 0, 'neutral complete generation')
    maximum, minimum = -10**9, 10**9
    for branch in range(1 << 15):
        right_charge = sum(q for i, q in enumerate(quantum_charges) if branch & (1 << i))
        left_charge = sum(quantum_charges) - right_charge
        flux = flux_sign * right_charge
        if modulus is not None:
            flux %= modulus
        require(flux - left_charge == 0 and -flux - right_charge == 0, 'branch operator Gauss')
        maximum, minimum = max(maximum, flux), min(minimum, flux)
        occupied = occupation_bits(branch)
        require(occupied.bit_count() == 15, 'fermion particle number')
        for i, q in enumerate(quantum_charges):
            direction = (branch >> i) & 1
            source, target = 2 * i + direction, 2 * i + 1 - direction
            actual, sign = hop(occupied, target, source)
            changed = branch ^ (1 << i)
            require(actual == occupation_bits(changed) and sign == 1, 'CAR hopping sign')
            shifted_flux = flux + (q if direction else -q)
            new_right = right_charge + (-q if direction else q)
            require(shifted_flux == -new_right, 'unwrapped coupled rotor shift')
    require((minimum, maximum) == (-18, 18), 'invariant finite flux range')
    return {'branches': 1 << 15, 'flux_min': minimum, 'flux_max': maximum,
            'checked_CAR_transitions': 15 * (1 << 15)}


def amplitudes_and_fluxes(prefix):
    """Direct branch enumeration, independent of a generating convolution.

    Return integer magnitudes with one common denominator; (-i)^popcount
    carries the complete relative phase.  Unprocessed channels stay left.
    """
    require(type(prefix) is int and 0 <= prefix <= 15, 'factor prefix')
    qs = charges()
    denominator = 17 ** min(prefix, 14) * (5 if prefix == 15 else 1)
    amplitudes, fluxes = [], []
    histogram = defaultdict(int)
    for branch in range(1 << prefix):
        magnitude, flux = 1, 0
        for i in range(prefix):
            right = (branch >> i) & 1
            left_coefficient, right_coefficient = (3, 4) if i == 14 else (15, 8)
            magnitude *= right_coefficient if right else left_coefficient
            if right:
                flux -= qs[i]
        amplitudes.append(magnitude)
        fluxes.append(flux)
        histogram[flux] += magnitude * magnitude
    require(sum(histogram.values()) == denominator * denominator, 'exact state normalization')
    probabilities = {e: F(weight, denominator * denominator) for e, weight in sorted(histogram.items())}
    mean = sum(e * p for e, p in probabilities.items())
    variance = sum(e * e * p for e, p in probabilities.items()) - mean * mean
    return amplitudes, fluxes, denominator, probabilities, mean, variance


SCOPE = {
    "finite_CAR_with_quantized_integer_electric_link": True,
    "complete_one_generation_representation_census": True,
    "normalized_pure_entangled_fermion_link_state": True,
    "operator_Gauss_on_every_supported_branch": True,
    "commuting_local_Cayley_factors_exactly_replayed": True,
    "same_Hamiltonian_current_and_electric_transfer": True,
    "exact_continuous_quantum_edge_solution": True,
    "finite_cyclic_electric_cutoff": False,
    "same_preparation_as_uniform_semiclassical_history": False,
    "all_144_links_dynamically_evolved": False,
    "nonabelian_gauge_or_Yukawa_dynamics": False,
    "nonabelian_quantum_Gauss": False,
    "electric_energy_or_plaquette_factors": False,
    "Cayley_to_continuous_error_bound": False,
    "Cayley_word_identified_with_continuous_exponential": False,
    "chiral_continuum_anomaly_measure_constructed": False,
    "native_source_action_state_clock_or_quantization_selected": False,
    "physical_current_or_laboratory_identification": False,
}


PREMISES = [
    "The prepared q5 source addresses, first oriented edge and Pauli frame are supplied.",
    "One registered generation and the finite canonical anticommutation relations are supplied.",
    "Each internal component occupies the same normalized left-end spin orbital; other vertices are empty.",
    "The active electric link is ell2(Z), with integer electric unit, bilateral unitary shift and hbar=1.",
    "Each channel has one fermion and two active spatial-spin modes; the orthogonal spin modes are empty and invariant.",
    "Hopping coefficients are 1 except 2 on the e_c channel; each Cayley factor has step 1/2.",
    "The initial electric state is |0>; all other link momenta remain zero and their hopping terms are omitted.",
    "Exact state amplitudes and the finite factor schedule are mathematical data, not sampled outcomes or source-selected operations.",
]


QUANTUM_LINK = {
    "gauge_constraint_group": "U1_hypercharge",
    "ambient_space": "ell2(Z)",
    "electric_action": "E|ell>=ell|ell>",
    "shift_action": "U|ell>=|ell+1>",
    "commutator": "[E,U]=U",
    "finite_wrap": False,
    "active_edge": 0,
    "other_links": "sharp zero electric values; no other hopping factor acts",
    "Hamiltonian_channel": "kappa_i*(c_L_i^dagger*U^q_i*c_R_i+c_R_i^dagger*U^(-q_i)*c_L_i)",
    "current_channel": "i*kappa_i*q_i*(c_L_i^dagger*U^q_i*c_R_i-c_R_i^dagger*U^(-q_i)*c_L_i)",
    "Heisenberg_electric_rate": "i[H_i,E]=-J_i",
    "factor": "V_i=(I+i*step*H_i/2)^(-1)*(I-i*step*H_i/2)",
    "factor_current": "quadratic expectation in (state_before+state_after)/2; midpoint vector is not normalized",
    "discrete_current_integral": "step*sum_i J_mid_i = -(mean_E_final-mean_E_initial)",
    "instantaneous_current": "expectation of sum_i J_i in the normalized endpoint state",
    "factor_schedule": "channel ids 0 through 14, each once; no continuous-time exponential or physical clock identified",
}

SOURCE_PATHS = tuple('code/sm_fermion_current/' + name for name in
                    ('quantum_link.py', 'build_quantum_link.py', 'verify_quantum_link.py',
                     'test_quantum_link.py', 'current.py', 'verify_current.py', 'README.md')) + (
    'code/sm_local_action/jet_action.py', 'code/sm_abelian_reduction/pilot.py',
    'Lean/Screen/WeylYukawaConventions.lean', 'Lean/Screen/A5FamilyBand.lean')


def validate_custody(packet):
    require(type(packet) is dict and set(packet) == {
        'schema', 'scope', 'premises', 'carrier', 'quantum_link', 'channel_basis',
        'channels', 'initial_state', 'history', 'final_state', 'observables',
        'source_pins', 'continuous_solution'}, 'root schema')
    same(packet['schema'], 'oph.sm_fermion_quantum_link.v1', 'schema')
    same(packet['scope'], SCOPE, 'scope')
    same(packet['premises'], PREMISES, 'premises')
    same(packet['quantum_link'], QUANTUM_LINK, 'quantum link algebra')
    require(type(packet['source_pins']) is dict and set(packet['source_pins']) == set(SOURCE_PATHS),
            'source pin census')
    for path in SOURCE_PATHS:
        data = (ROOT / path).read_bytes()
        same(packet['source_pins'][path], {'bytes': len(data), 'sha256': sha256(data).hexdigest()},
             'source content hash')
    return {'custody_verified': True, 'quantum_state_verified': False}


def channel_check(received):
    result = []
    for name, count, charge in FIELDS:
        for component in range(count):
            i = len(result)
            kappa = 2 if name == 'e_c' else 1
            a, beta = (F(3, 5), F(4, 5)) if i == 14 else (F(15, 17), F(8, 17))
            # Solve the two scalar midpoint equations, rather than call Cayley.
            require(1 - a == kappa * beta / 4 and beta == kappa * (1 + a) / 4,
                    'channel midpoint equation')
            require(a * a + beta * beta == 1, 'channel unitary normalization')
            result.append({'id': i, 'multiplet': name, 'internal_component': component,
                           'integer_charge': charge, 'mode_indices': [2 * i, 2 * i + 1],
                           'kappa': str(kappa), 'step': '1/2',
                           'left_amplitude': [str(a), '0'], 'right_amplitude': ['0', str(-beta)]})
    same(received, result, 'all charged CAR channels')
    return result


def basis_check(received, geometry):
    require(geometry['edges'][0] == {'ends': [0, 16], 'axis': 0}, 'actual source edge')
    # sigma_x exchanges the two entries; its -i adjoint fixes the right mode.
    spin = [F(3, 5), F(4, 5)]
    require(sum(x * x for x in spin) == 1, 'occupied spin orbital normalization')
    expected = {
        'active_edge': 0, 'ends': [0, 16], 'axis': 0,
        'left_spin': [[str(x), '0'] for x in spin],
        'right_spin': [['0', str(-x)] for x in reversed(spin)],
        'spin_relation': 'right_spin=(i*sigma_axis)^dagger*left_spin',
        'CAR_mode_order': 'L0,R0,L1,R1,...,L14,R14',
        'spectators': 'all orthogonal spin modes and all other source sites are empty and invariant under these factors',
        'branch_embedding': '|mask> maps to ordered occupied fermion modes tensor |electric=-sum moved integer charges>',
        'one_particle_modes_on_carrier': 1920, 'active_fermion_modes': 30,
        'fermions': 15, 'invariant_sector_dimension': 32768,
    }
    same(received, expected, 'physical spin/CAR embedding')


def readout_from_branches(prefix):
    _, _, _, histogram, mean, variance = amplitudes_and_fluxes(prefix)
    purity = sum(p * p for p in histogram.values())
    return {'norm_squared': '1', 'mean_electric': str(mean),
            'electric_second_moment': str(variance + mean * mean),
            'electric_variance': str(variance), 'mean_charge_left': str(mean),
            'mean_charge_right': str(-mean), 'link_reduced_purity': str(purity),
            'flux_probabilities': [{'flux': e, 'probability': str(p)} for e, p in histogram.items()]}


def history_check(history):
    require(type(history) is list and len(history) == 16, 'complete factor history')
    summed_current = F(0)
    previous_mean = F(0)
    for prefix in range(16):
        result = readout_from_branches(prefix)
        factor = None
        if prefix:
            i = prefix - 1
            kappa, q = (2 if i == 14 else 1), charges()[i]
            a, beta = (F(3, 5), F(4, 5)) if i == 14 else (F(15, 17), F(8, 17))
            # Matrix derivative J=i*kappa*q*(|L><R|-|R><L|).
            m_left, m_imaginary = (1 + a) / 2, -beta / 2
            current = -2 * kappa * q * m_left * m_imaginary
            delta_electric = F(result['mean_electric']) - previous_mean
            require(delta_electric == -current / 2, 'current/quantum flux transfer')
            factor = {'channel': i, 'midpoint_norm_squared': str(m_left * m_left + m_imaginary * m_imaginary),
                      'midpoint_current': str(current), 'discrete_integrated_current': str(current / 2),
                      'mean_electric_change': str(delta_electric)}
            summed_current += current
        same(history[prefix], {'completed_factors': prefix, 'factor': factor, 'readout': result},
             'branch-derived factor history')
        previous_mean = F(result['mean_electric'])
    require(summed_current == F(36288, 7225), 'nonzero factor midpoint current')
    return result, summed_current


def branch_census_check(received):
    algebra = branch_algebra()
    magnitudes, fluxes, denominator, probabilities, _, _ = amplitudes_and_fluxes(15)
    stream = sha256()
    current_numerator = 0
    for mask, (magnitude, flux) in enumerate(zip(magnitudes, fluxes)):
        amplitude = F(magnitude, denominator)
        phase = mask.bit_count() % 4
        values = ([str(amplitude), '0'], ['0', str(-amplitude)],
                  [str(-amplitude), '0'], ['0', str(amplitude)])[phase]
        row = {'mask': mask, 'occupied_modes': [2 * i + ((mask >> i) & 1) for i in range(15)],
               'electric_flux': flux, 'charge_left': flux, 'charge_right': -flux,
               'amplitude': values, 'gauss_residual': [0, 0]}
        stream.update(canonical(row))
        # Sum actual off-diagonal CAR/shift matrix elements once per pair.
        for i, q in enumerate(charges()):
            if not mask & (1 << i):
                other = mask | (1 << i)
                require(fluxes[other] == flux - q, 'phase-current rotor matrix element')
                current_numerator += 2 * (2 if i == 14 else 1) * q * magnitude * magnitudes[other]
    current = F(current_numerator, denominator * denominator)
    require(current == F(47232, 7225), 'actual quantum current expectation')
    require(set(probabilities) == set(range(-18, 19)), 'complete nonwrapped integer support')
    same(received, {'branch_count': algebra['branches'], 'branch_sha256': stream.hexdigest(),
                    'branch_order': 'ascending integer mask; bit i selects right mode of channel i',
                    'branch_record_fields': ['mask', 'occupied_modes', 'electric_flux', 'charge_left',
                                             'charge_right', 'amplitude', 'gauss_residual'],
                    'norm_squared': '1', 'all_branch_Gauss_residuals': [0, 0],
                    'other_vertex_Gauss_residuals': 'zero, since their fermion modes and incident electric momenta are empty'},
         'complete complex branch census')
    return current, algebra


def continuous_derivation():
    """Derive the exact commuting-channel solution and its current in SymPy."""
    import sympy as s
    t = s.symbols('t', real=True)
    kappa = s.symbols('kappa', real=True)
    psi = s.Matrix([s.cos(kappa * t), -s.I * s.sin(kappa * t)])
    h = kappa * s.Matrix([[0, 1], [1, 0]])
    require(all(s.simplify(z) == 0 for z in s.I * psi.diff(t) - h * psi),
            'continuous exact Schrodinger state')
    require(s.simplify((psi.conjugate().T * psi)[0]) == 1,
            'continuous normalized state')
    require(psi.subs(t, 0) == s.Matrix([1, 0]), 'continuous initial state')
    expected_electric = sum(-q * s.sin((2 if i == 14 else 1) * t) ** 2 for i, q in enumerate(charges()))
    current = sum(q * (2 if i == 14 else 1) * s.sin(2 * (2 if i == 14 else 1) * t)
                  for i, q in enumerate(charges()))
    variance = sum(q * q * s.sin((2 if i == 14 else 1) * t) ** 2 *
                   s.cos((2 if i == 14 else 1) * t) ** 2 for i, q in enumerate(charges()))
    require(s.trigsimp(s.diff(expected_electric, t) + current) == 0, 'continuous Heisenberg current')
    require(s.expand_trig(expected_electric - 6 * (s.sin(t) ** 2 - s.sin(2 * t) ** 2)) == 0,
            'continuous electric law')
    require(s.expand_trig(current + 6 * s.sin(2 * t) - 12 * s.sin(4 * t)) == 0,
            'continuous current law')
    require(s.simplify(expected_electric.subs(t, s.pi / 6)) == -3, 'continuous flux witness')
    require(s.simplify(current.subs(t, s.pi / 6)) == 3 * s.sqrt(3), 'continuous current witness')
    require(s.simplify(variance.subs(t, s.pi / 6)) == s.Rational(45, 2), 'continuous variance witness')
    return expected_electric, current, variance


def continuous_point_readout():
    # Enumerate probabilities at t=pi/6 in Q; the ket itself lies in Q(i,sqrt3).
    histogram = defaultdict(int)
    for mask in range(1 << 15):
        numerator, electric = 1, 0
        for i, q in enumerate(charges()):
            bit = (mask >> i) & 1
            numerator *= (3 if bit else 1) if i == 14 else (1 if bit else 3)
            electric -= bit * q
        histogram[electric] += numerator
    denominator = 4 ** 15
    require(sum(histogram.values()) == denominator, 'continuous exact state normalization')
    probabilities = {e: F(n, denominator) for e, n in sorted(histogram.items())}
    mean = sum(e * p for e, p in probabilities.items())
    second = sum(e * e * p for e, p in probabilities.items())
    require(mean == -3 and second - mean * mean == F(45, 2), 'continuous point quantum moments')
    return {'norm_squared': '1', 'mean_electric': str(mean), 'electric_second_moment': str(second),
            'electric_variance': str(second - mean * mean), 'mean_charge_left': str(mean),
            'mean_charge_right': str(-mean), 'link_reduced_purity': str(sum(p * p for p in probabilities.values())),
            'flux_probabilities': [{'flux': e, 'probability': str(p)} for e, p in probabilities.items()]}


def continuous_check(received):
    continuous_derivation()
    expected = {
        'Hamiltonian': 'sum_i H_i on the same infinite-link invariant sector',
        'commuting_channels': True,
        'channel_left_amplitude': 'cos(kappa_i*t)',
        'channel_right_amplitude': '-i*sin(kappa_i*t)',
        'state_embedding': 'the same ordered fermion/flux branch embedding; no cyclic cutoff',
        'mean_electric': '6*(sin(t)^2-sin(2*t)^2)',
        'instantaneous_current': '-6*sin(2*t)+12*sin(4*t)',
        'electric_variance': '84*sin(t)^2*cos(t)^2+36*sin(2*t)^2*cos(2*t)^2',
        'electric_derivative': 'd(mean_electric)/dt=-instantaneous_current',
        'operator_Gauss': 'zero for every real t',
        'exact_point': {'model_time': 'pi/6', 'readout': continuous_point_readout(),
                        'instantaneous_current_Qsqrt3': ['0', '3']},
        'physical_time_calibrated': False,
        'Cayley_word_is_this_exact_time': False,
    }
    same(received, expected, 'exact continuous quantum solution')


def verify(packet):
    validate_custody(packet)
    geometry = geometry_check(packet['carrier'])
    channel_check(packet['channels'])
    basis_check(packet['channel_basis'], geometry)
    same(packet['initial_state'], {
        'occupied_modes': list(range(0, 30, 2)), 'electric_flux': 0, 'total_charge': 0,
        'norm_squared': '1', 'operator_Gauss': 'zero at every vertex',
        'preparation': 'localized generation at the left endpoint, distinct from the uniform semiclassical Slater state'},
         'gauge-constrained preparation')
    readout, summed_current = history_check(packet['history'])
    require(type(packet['final_state']) is dict and set(packet['final_state']) == {'branch_census', 'readout'},
            'final quantum state schema')
    same(packet['final_state']['readout'], readout, 'final state/current join')
    current, algebra = branch_census_check(packet['final_state']['branch_census'])
    purity = rational(readout['link_reduced_purity'])
    require(0 < purity < 1, 'pure-state electric/matter entanglement witness')
    same(packet['observables'], {
        'sum_factor_midpoint_currents': str(summed_current),
        'discrete_integrated_current': str(summed_current / 2),
        'final_instantaneous_current': str(current), 'initial_instantaneous_current': '0',
        'mean_electric_change': readout['mean_electric'], 'operator_Gauss_norm_squared': '0',
        'entanglement_witness': 'pure joint state with mixed electric reduction; exact reduced purity below 1',
        'finite_support': [-18, 18], 'ambient_electric_cutoff': None}, 'quantum operator observables')
    continuous_check(packet['continuous_solution'])
    return {'verified': True, 'independent_quantum_state_replay': True,
            'fermions': 15, 'carrier_sites': 64, 'active_quantum_links': 1,
            'branches': algebra['branches'], 'checked_CAR_transitions': algebra['checked_CAR_transitions'],
            'flux_bins': 37, 'operator_Gauss_norm_squared': '0',
            'mean_electric': readout['mean_electric'], 'final_current': str(current),
            'exact_continuous_solution': True, 'continuous_point_mean_electric': '-3',
            'operator_Gauss': True, 'exact_continuous_quantum_edge_solution': True,
            'gauge_constraint_group': 'U1_hypercharge', 'nonabelian_quantum_Gauss': False,
            'physical_current_attached': False,
            'continuous_point_current': '3*sqrt(3)', 'continuous_point_variance': '45/2',
            'cyclic_link_truncation': False, 'physical_clock_or_SM_dynamics_selected': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=OUTPUT)
    print(json.dumps(verify(load(parser.parse_args().path)), sort_keys=True, indent=2))
