"""Complete continuing source execution and exact finite-array interface bound."""
from fractions import Fraction as F
import numpy as np

from ..calibrated_measurement.bounds import root_interval
from ..primitive_observer import runtime
from . import bounds, randomized, source

VECTOR = ((F(-1), F(-1)), (F(1), F(0)), (F(0), F(1)))


def _matrix(array):
    return [[(F(complex(z).real), F(complex(z).imag)) for z in row] for row in array]


def _adjoint(a):
    return [[bounds.conj(a[j][i]) for j in range(len(a))] for i in range(len(a))]


def _add(a, b):
    return [[bounds.add(x, y) for x, y in zip(row, other)] for row, other in zip(a, b)]


def _scale(a, value):
    return [[bounds.mul(z, value) for z in row] for row in a]


def _product(a, b):
    result = []
    for i in range(len(a)):
        row = []
        for j in range(len(b[0])):
            value = (F(0), F(0))
            for k in range(len(b)):
                value = bounds.add(value, bounds.mul(a[i][k], b[k][j]))
            row.append(value)
        result.append(row)
    return result


def _norm_upper(a):
    return root_interval(sum((bounds.norm2(z) for row in a for z in row), F(0)))[1]


def finite_certificate(design, matrices, operators):
    """A positive square compared to the exact binary raw detector arrays.

    Every matrix operation here uses rational arithmetic. Each binary raw
    instrument is compared with its trace-preserving completion K T^-1/2;
    the execution itself retains the raw arrays and their mass defect.
    """
    if type(design) is not randomized.Design or len(operators) != 6:
        raise ValueError('a certified Gram design and six source probe instruments required')
    if any(np.shape(matrices[name]) != (2, 2) for name in ('U', 'V')) or any(
            len(maps) != 2 or any(np.shape(k) != (2, 2) for k in maps) for maps in operators):
        raise ValueError('two-dimensional source operators and complete binary probes required')
    identity = _matrix(np.eye(2))
    u, v = (_matrix(matrices[name]) for name in ('U', 'V'))
    polynomial = _add(_add(_scale(identity, design.vector[0]), _scale(_adjoint(u), design.vector[1])),
                      _scale(v, design.vector[2]))
    positive = _product(_adjoint(polynomial), polynomial)
    observable = _scale(identity, (design.constant, F(0)))
    masses, completion = [], F(0)
    for coefficient, maps in zip(design.coefficients, operators, strict=True):
        effects = [_product(_adjoint(_matrix(k)), _matrix(k)) for k in maps]
        difference = _add(effects[0], _scale(effects[1], (F(-1), F(0))))
        observable = _add(observable, _scale(difference, (coefficient, F(0))))
        mass = _norm_upper(_add(_add(*effects), _scale(identity, (F(-1), F(0)))))
        if mass >= 1:
            raise ValueError('raw instrument has no certified invertible completion')
        masses.append(mass)
        completion += abs(coefficient)*mass*(1+mass)/(1-mass)
    raw = _norm_upper(_add(observable, _scale(positive, (F(-1), F(0)))))
    return {'raw_positive_square_comparison_error_upper': str(raw),
            'setting_effect_sum_error_upper': list(map(str, masses)),
            'weighted_completion_error_upper': str(completion),
            'conditional_score_defect_upper': str(raw+completion),
            'conditional_score_defect_upper_float': float(raw+completion),
            'arithmetic': 'exact binary-array rational matrix arithmetic and outward rational square root',
            'bound_uniform_in_state_and_reference': True,
            'completion_is_mathematical_comparison_only': True,
            'physical_service_fidelity_certified': False}


def execute():
    design = randomized.witness_design(VECTOR)
    words = source.source_words()
    matrices = {name: source.compiler.execute(word, 1) for name, (word, _) in words.items()}
    operators = []
    for name in ('U', 'V', 'UV'):
        for phase in (0., -np.pi/2):
            operators.append(source.ramsey(words[name][1], phase)[0])
    certificate = finite_certificate(design, matrices, operators)
    active = [i for i, p in enumerate(design.probabilities) if p]
    roots = [np.sqrt(float(design.probabilities[i])) for i in active]
    raw_choice = [F(float(root))**2 for root in roots]
    choice_mass = sum(raw_choice, F(0))
    choice_bias = design.width*sum((abs(p/choice_mass-design.probabilities[i])
                                  for i, p in zip(active, raw_choice)), F(0))
    certificate['completed_randomizer_score_error_upper'] = str(choice_bias)
    certificate['raw_randomizer_mass_defect'] = str(abs(choice_mass-1))
    certificate['conditional_score_defect_including_completed_randomizer_upper'] = str(
        F(certificate['conditional_score_defect_upper'])+choice_bias)
    choice = runtime.Instrument(kraus=tuple((root*np.eye(2),)
        for root in roots), visible=True, binding={'exact_target_probabilities': [str(design.probabilities[i]) for i in active]})
    probes = tuple(runtime.Instrument(kraus=((operators[i][0],), (operators[i][1],)), visible=True,
        binding={'setting': randomized.SETTINGS[i], 'source_compiled': True}) for i in active)
    feedback = tuple(runtime.Instrument(kraus=((matrices[name],),),
        binding={'source_word': name}) for name in ('U', 'V'))
    program = runtime.Program('randomized continuing Gram witness', 2, np.eye(2)/2,
        (runtime.Event('setting_1', instruments=(choice,)),
         runtime.Event('probe_1', instruments=probes, control=0),
         runtime.Event('record_feedback', instruments=feedback, control=1),
         runtime.Event('setting_2', instruments=(choice,)),
         runtime.Event('probe_2', instruments=probes, control=2)))
    bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
    initial = np.outer(bell, bell)
    full = runtime.execute(program, initial=initial)
    prefix = runtime.execute(program, initial=initial, stop_after=3)
    resumed = runtime.execute(program, checkpoint=runtime.checkpoint_from_json(
        program, runtime.checkpoint_to_json(prefix)))
    residual = max(float(np.linalg.norm(a.density-b.density, 1)) for a, b in zip(full.branches, resumed.branches))
    rows, conditionals = [], {}
    for branch in full.branches:
        record = branch.record
        scores = [float(design.score(active[record[k]], 1-2*record[k+1])) for k in (0, 2)]
        rows.append({'record': list(record), 'probability': branch.probability, 'scores': scores})
        key = record[:2]
        weight, weighted = conditionals.get(key, (0., 0.))
        conditionals[key] = weight+branch.probability, weighted+branch.probability*scores[1]
    conditional_means = [{'first_record': list(key), 'mass': weight, 'second_score_conditional_mean': weighted/weight}
                         for key, (weight, weighted) in sorted(conditionals.items())]
    # Each setting separately has a valid quantum state, but selecting that
    # state after seeing the setting can destroy the shared Gram inequality.
    selected_mean = float(design.constant)
    for coefficient, maps in zip(design.coefficients, operators, strict=True):
        observable = maps[0].conj().T@maps[0]-maps[1].conj().T@maps[1]
        eigenvalues = np.linalg.eigvalsh(float(coefficient)*observable)
        selected_mean += float(eigenvalues[0])
    return {'design': design.specification(), 'finite_array_certificate': certificate,
            'complete_two_trial_law': rows, 'conditional_second_scores': conditional_means,
            'minimum_conditional_second_score': min(row['second_score_conditional_mean'] for row in conditional_means),
            'maximum_reference_checkpoint_residual': residual,
            'raw_complete_mass_defect': abs(full.probability_mass-1),
            'preparation_reset_between_trials': False, 'quantum_reference_retained': True,
            'setting_correlated_preparation_negative_control_mean': selected_mean,
            'setting_randomness_is_supplied_control_service': True,
            'physical_randomizer_independence_certified': False,
            'fixed_horizon_theorem_applied_to_laboratory_data': False}
