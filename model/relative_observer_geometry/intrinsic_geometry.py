"""Exact single-carrier spatial readouts and their scheduling boundary.

The estimator sees labelled preparations and currents. It removes positive
seam exposure without importing coordinates. Even preparations also permit
exact recovery from finite native shots without a probability-law oracle.
The result concerns finite structural geometry, not physical distances.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from pathlib import Path

import sympy as sp

from ..observer import native
from ..observer_reconstruction import native_geometry as ng
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json


ROOT = Path(__file__).resolve().parents[2]
PORTS = 12
SOURCE_PATHS = (
    'model/relative_observer_geometry/intrinsic_geometry.py',
    'model/observer/native.py',
    'model/observer_reconstruction/native_geometry.py',
    'model/common_source/slow_geometry.py',
    'evidence/native_geometric_source_20260925/spec.json',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
    'Lean/Screen/PrimitivePortFrameQuotient.lean',
    'Lean/Geometry/SourceDerivedSpacetimeCarrier.lean',
)


def recover_structure(records):
    """Remove seam exposure from complete, possibly biased observation laws.

    Each preparation is 2 e_p. Per-preparation and per-seam schedules may
    differ, and coin probabilities may be biased or deterministic. Every
    preparation must expose the same seam support with positive mass. The
    estimator checks currents against that disclosed preparation. It does
    not use the hidden final loads or a source adjacency table.
    """
    if type(records) not in (list, tuple) or not records:
        raise ValueError('complete positive observation laws required')
    exposures = [{} for _ in range(PORTS)]
    currents = [{} for _ in range(PORTS)]
    masses = [F(0) for _ in range(PORTS)]
    seen = set()
    for row in records:
        if type(row) is not dict or set(row) != ng.RECORD_FIELDS:
            raise ValueError('only declared observation record fields accepted')
        p = row['preparation_port']
        seam = row['seam']
        coin, current = row['coin'], row['current']
        if type(p) is not int or p not in range(PORTS):
            raise ValueError('preparation port must be a native port')
        if (type(seam) not in (list, tuple) or len(seam) != 2
                or any(type(i) is not int or i not in range(PORTS) for i in seam)
                or seam[0] >= seam[1]):
            raise ValueError('one increasing orientation per seam required')
        a, b = seam
        if type(coin) is not int or coin not in (0, 1):
            raise ValueError('a native coin is zero or one')
        if type(current) is not int or current not in (-1, 0, 1):
            raise ValueError('a native total-two current is integral')
        probability = ng._probability(row['probability'])
        key = p, a, b, coin
        if key in seen:
            raise ValueError('duplicate calibration outcome')
        seen.add(key)
        # In a 2e_p preparation every endpoint difference is even. The
        # native current equals 1_{p=a}-1_{p=b}, regardless of the coin.
        if current != int(p == a)-int(p == b):
            raise ValueError('current contradicts even-load native preparation')
        seam = a, b
        exposures[p][seam] = exposures[p].get(seam, F(0))+probability
        currents[p][seam] = currents[p].get(seam, F(0))+probability*current
        masses[p] += probability
    if any(mass != 1 for mass in masses):
        raise ValueError('each labelled preparation requires normalized complete laws')
    seams = tuple(sorted(exposures[0]))
    if not seams or any(set(exposure) != set(seams) for exposure in exposures):
        raise ValueError('the same positive seam support is required for every preparation')
    laplacian = [[F(0) for _ in range(PORTS)] for _ in range(PORTS)]
    for p in range(PORTS):
        for a, b in seams:
            conditional_current = currents[p][a, b]/exposures[p][a, b]
            laplacian[a][p] += conditional_current
            laplacian[b][p] -= conditional_current
    checks, projector = ng._band_certificate(laplacian)
    symmetric = all(laplacian[i][j] == laplacian[j][i]
                    for i in range(PORTS) for j in range(PORTS))
    graph = symmetric and all(laplacian[i][i] == 5 for i in range(PORTS))
    graph = graph and all(laplacian[i][j] in (0, -1)
                         for i in range(PORTS) for j in range(PORTS) if i != j)
    certified = graph and all(checks.values())
    return {
        'record_count': len(records), 'seams': [list(seam) for seam in seams],
        'laplacian_exact': [[str(x) for x in row] for row in laplacian],
        'seam_exposures_exact': [
            {f'{a},{b}': str(exposure[a, b]) for a, b in seams}
            for exposure in exposures],
        'rank_three_structural_gram_certified': certified,
        'band_certificate': checks,
        'gram_qsqrt5': ng._encode_q5([
            [ng._qscale(x, 4) for x in row] for row in projector]) if certified else None,
        'uniform_scheduler_or_fair_coins_required_for_this_readout': False,
        'full_response_slow_band_identified_for_the_actual_scheduler': False,
        'hidden_final_loads_or_source_coordinates_used': False,
    }


def execute_biased_calibration(seams):
    """Execute native endpoints with positive preparation-dependent exposures.

    Coin zero alone is used on even preparation labels; odd labels retain
    coins with probabilities 1/3 and 2/3. These are alternative declared
    protocols, not a change to the model's production scheduler.
    """
    records = []
    for p in range(PORTS):
        numerators = [1+(index+3*p) % 7 for index in range(len(seams))]
        total = sum(numerators)
        coins = ((0, F(1)),) if p % 2 == 0 else ((0, F(1, 3)), (1, F(2, 3)))
        for seam, weight in zip(seams, numerators, strict=True):
            for coin, coin_probability in coins:
                _, (_, current) = native.native(tuple(2*int(i == p) for i in range(PORTS)), seam, coin)
                records.append({'preparation_port': p, 'seam': list(seam),
                    'coin': coin, 'current': current,
                    'probability': str(F(weight, total)*coin_probability)})
    return records


def execute_single_shot_calibration(seams, coin_value=None):
    """Execute one native branch per preparation/seam, with fixed varied coins."""
    if coin_value is not None and (type(coin_value) is not int or coin_value not in (0, 1)):
        raise ValueError('fixed coin must be zero or one')
    records = []
    for p in range(PORTS):
        loads = tuple(2*int(i == p) for i in range(PORTS))
        for index, seam in enumerate(seams):
            coin = (p+index) % 2 if coin_value is None else coin_value
            _, (_, current) = native.native(loads, seam, coin)
            records.append({'preparation_port': p, 'seam': list(seam),
                            'coin': coin, 'current': current})
    return records


def recover_from_shots(records):
    """Recover exactly in the declared native model without probability inputs.

    Empirical frequencies are only bookkeeping weights. Conditional currents
    have zero variance on the known even preparations, so any observed coin
    frequencies give the exact conditional mean. Native rule, preparation and
    faithful record access remain premises; this is not a lab-noise theorem.
    """
    if type(records) not in (list, tuple) or not records:
        raise ValueError('nonempty native calibration shots required')
    counts, currents = Counter(), {}
    total = Counter()
    for row in records:
        if type(row) is not dict or set(row) != ng.RECORD_FIELDS-{'probability'}:
            raise ValueError('shots contain observations, not probability laws')
        p, seam = row['preparation_port'], row['seam']
        coin, current = row['coin'], row['current']
        if (type(p) is not int or p not in range(PORTS)
                or type(seam) not in (list, tuple) or len(seam) != 2
                or any(type(i) is not int or i not in range(PORTS) for i in seam)
                or type(coin) is not int or coin not in (0, 1)
                or type(current) is not int):
            raise ValueError('native preparation, seam, coin and current required')
        key = p, *seam, coin
        if key in currents and currents[key] != current:
            raise ValueError('inconsistent currents for a deterministic native calibration')
        currents[key] = current
        counts[key] += 1
        total[p] += 1
    weighted = [{'preparation_port': p, 'seam': [a, b], 'coin': coin,
                 'current': currents[p, a, b, coin],
                 'probability': str(F(count, total[p]))}
                for (p, a, b, coin), count in sorted(counts.items())]
    result = recover_structure(weighted)
    result['shot_count'] = len(records)
    result['probability_law_oracle_used'] = False
    result['empirical_frequencies_are_only_bookkeeping_weights'] = True
    result['conditional_current_variance_zero_under_native_model'] = True
    result['laboratory_noise_or_preparation_error_bounded'] = False
    return result


def _qinverse(x):
    denominator = x[0]*x[0]-5*x[1]*x[1]
    if denominator == 0:
        raise ZeroDivisionError('zero quadratic-field element')
    return x[0]/denominator, -x[1]/denominator


def qrank(matrix):
    """Gaussian rank over Q(sqrt(5)), with no floating threshold."""
    rows = [[tuple(x) for x in row] for row in matrix]
    if not rows or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError('nonempty rectangular matrix required')
    pivot = 0
    zero = F(0), F(0)
    for column in range(len(rows[0])):
        index = next((i for i in range(pivot, len(rows)) if rows[i][column] != zero), None)
        if index is None:
            continue
        rows[pivot], rows[index] = rows[index], rows[pivot]
        inv = _qinverse(rows[pivot][column])
        rows[pivot] = [ng._qmul(x, inv) for x in rows[pivot]]
        for i in range(pivot+1, len(rows)):
            factor = rows[i][column]
            rows[i] = [ng._qadd(x, ng._qscale(ng._qmul(factor, y), -1))
                       for x, y in zip(rows[i], rows[pivot], strict=True)]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def scheduling_control(laplacian, seam):
    """Exactly certify the slow-band split under one doubled seam exposure."""
    a, b = seam
    edge = [F(int(i == a)-int(i == b)) for i in range(PORTS)]
    weighted = [[laplacian[i][j]+edge[i]*edge[j]
                 for j in range(PORTS)] for i in range(PORTS)]
    lam = F(5), F(-1)
    shifted = [[ng._qadd((weighted[i][j], F(0)), ng._qscale(lam, -int(i == j)))
                for j in range(PORTS)] for i in range(PORTS)]
    multiplicity = PORTS-qrank(shifted)
    if multiplicity != 2:
        raise AssertionError('one doubled seam must leave exactly two slow eigenvectors')
    # L + bb^T >= L and bb^T annihilates constants. Hence there is no
    # nonconstant eigenvalue below lambda_min(L), while exact nullity two
    # above proves that this smallest positive eigenvalue is attained twice.
    return {
        'doubled_seam': list(seam),
        'scheduler': 'weight 2/31 on one seam, 1/31 on each other seam',
        'actual_mean_response': 'I-(L+b*b^T)/62; b=e_a-e_b',
        'smallest_positive_weighted_laplacian_eigenvalue_qsqrt5': ['5', '-1'],
        'smallest_positive_weighted_laplacian_multiplicity_exact': multiplicity,
        'slowest_nonconstant_actual_response_qsqrt5': ['57/62', '1/62'],
        'perturbation_positive_semidefinite_exact': True,
        'normalized_long_attempt_response_gram_rank': multiplicity,
        'structural_gram_rank_after_exposure_removal': 3,
        'time_parameter': 'uniformly indexed native attempts with the declared biased seam scheduler',
    }


def _inverse(matrix):
    size = len(matrix)
    rows = [list(row)+[F(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((i for i in range(column, size) if rows[i][column]), None)
        if pivot is None:
            raise ValueError('singular matrix')
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [x/scale for x in rows[column]]
        for i in range(size):
            if i != column:
                factor = rows[i][column]
                rows[i] = [x-factor*y for x, y in zip(rows[i], rows[column], strict=True)]
    return [row[size:] for row in rows]


def inverse_response_metric(laplacian):
    """Compute all-mode inverse-response geometry on the constant-free space.

    For graph L, G=(L+J/12)^-1-J/12 is the positive inverse on 1^perp.
    R_ab=(e_a-e_b)^T G (e_a-e_b) is effective resistance. Its square root
    is a Euclidean distance with centered Gram G of rank eleven. This is
    another mathematical readout of the same response, not a rod apparatus.
    """
    inverse = _inverse([[x+F(1, PORTS) for x in row] for row in laplacian])
    gram = [[x-F(1, PORTS) for x in row] for row in inverse]
    centered_identity = [[sum(laplacian[i][k]*gram[k][j] for k in range(PORTS))
                          for j in range(PORTS)] for i in range(PORTS)]
    expected = [[F(i == j)-F(1, PORTS) for j in range(PORTS)] for i in range(PORTS)]
    if centered_identity != expected or any(sum(row) for row in gram):
        raise AssertionError('inverse response must invert exactly off the constant mode')
    rank = qrank([[(x, F(0)) for x in row] for row in gram])
    resistance = [[gram[i][i]+gram[j][j]-2*gram[i][j]
                   for j in range(PORTS)] for i in range(PORTS)]
    if any(resistance[i][j] > resistance[i][k]+resistance[k][j]
           for i in range(PORTS) for j in range(PORTS) for k in range(PORTS)):
        raise AssertionError('effective resistance triangle inequality failed')
    shells = Counter(str(resistance[i][j]) for i in range(PORTS) for j in range(i+1, PORTS))
    if rank != 11 or dict(shells) != {'11/30': 30, '7/15': 30, '1/2': 6}:
        raise AssertionError('unexpected exact carrier inverse-response metric')
    return {
        'gram_exact': [[str(x) for x in row] for row in gram],
        'squared_euclidean_distance_exact': [[str(x) for x in row] for row in resistance],
        'squared_distance_shell_counts': dict(sorted(shells.items())),
        'centered_gram_rank_exact': rank,
        'L_times_G_equals_constant_free_identity_exact': True,
        'resistance_triangle_inequality_all_triples_exact': True,
        'source': 'inverse of the recovered unweighted response on the constant-free space',
        'physical_distance_or_dimension_selected': False,
    }


def wave_density_metric(spatial_tensor, measure_density=1, time_weight=1):
    """Recover a comparison Lorentz metric from a four-dimensional action.

    Inputs are exact rational coefficients at one point, not inferred native
    continuum coefficients. The continuum argument in README obtains a
    measurable positive spatial tensor from the already proved local form.
    A smooth tensor, radial geometry and actual wave apparatus are distinct
    obligations. The time weight names the chosen wave-clock coordinate.
    """
    def rational(value):
        if type(value) not in (int, F) and not isinstance(value, sp.Rational):
            raise ValueError('exact rational coefficients required')
        return sp.Rational(value)

    try:
        rows = [list(row) for row in spatial_tensor]
    except TypeError as exc:
        raise ValueError('a three-dimensional spatial tensor is required') from exc
    if len(rows) != 3 or any(len(row) != 3 for row in rows):
        raise ValueError('a three-dimensional spatial tensor is required')
    a = sp.Matrix([[rational(value) for value in row] for row in rows])
    mu, weight = rational(measure_density), rational(time_weight)
    if mu <= 0 or weight <= 0:
        raise ValueError('positive measure density and time weight required')
    if a != a.T or any(a[:k, :k].det() <= 0 for k in range(1, 4)):
        raise ValueError('a symmetric positive definite spatial tensor is required')
    coefficient = mu*sp.diag(-weight, a)
    volume = sp.sqrt(-coefficient.det())
    metric = sp.simplify(volume*coefficient.inv())
    inverse = sp.simplify(coefficient/volume)
    if sp.simplify(metric.det()-coefficient.det()) != 0:
        raise ArithmeticError('four-dimensional density determinant identity failed')
    if sp.simplify(volume*inverse-coefficient) != sp.zeros(4):
        raise ArithmeticError('metric did not reproduce the complete action density')
    return {'action_density': coefficient, 'metric': metric,
            'inverse_metric': inverse, 'metric_volume_density': volume}


def product_wave_metric_control():
    """Exact tensor/control checks, separate from finite carrier readbacks."""
    # This anisotropic rational tensor is a test coefficient, not a fitted
    # value of the subsequential native continuum tensor.
    a = sp.Matrix([[1, 0, 0], [0, 2, 1], [0, 1, 3]])
    data = wave_density_metric(a.tolist(), F(3, 2))
    doubled = wave_density_metric(a.tolist(), 3)
    rate = sp.Rational(3, 2)
    rescaled = wave_density_metric(a.tolist(), rate*sp.Rational(3, 2), 1/rate**2)
    jacobian = sp.diag(rate, 1, 1, 1)  # old coordinate t = rate * new coordinate
    if sp.simplify(rescaled['metric']-jacobian.T*data['metric']*jacobian) != sp.zeros(4):
        raise ArithmeticError('clock reparameterization must transform the same metric')
    if sp.simplify(doubled['metric']-2*data['metric']) != sp.zeros(4):
        raise ArithmeticError('action normalization must retain its conformal factor')
    encode = lambda matrix: [[str(sp.simplify(value)) for value in row]
                             for row in matrix.tolist()]
    return {
        'coefficient_fixture': 'anisotropic rational product tensor; not measured native coefficients',
        'spatial_tensor': encode(a),
        'measure_density': '3/2',
        'action_density': encode(data['action_density']),
        'metric': encode(data['metric']),
        'metric_volume_density': str(data['metric_volume_density']),
        'signature': [-1, 1, 1, 1],
        'complete_action_density_recovered_exactly': True,
        'clock_reparameterization_rate': str(rate),
        'clock_reparameterization_is_tensor_covariant': True,
        'doubling_measure_preserves_cone_but_doubles_metric': True,
        'continuum_input': 'regular strongly local H1 form and chosen product wave action',
        'coefficient_regularities_proved': ['measurable', 'uniformly elliptic in local charts'],
        'smooth_coefficients_or_curvature_proved': False,
        'native_radial_transport_or_lab_proper_time_proved': False,
        'observed_spacetime_or_Einstein_dynamics_proved': False,
    }


def run():
    uniform_records = ng.execute_calibration()
    uniform = ng.recover(uniform_records)
    structural = recover_structure(uniform_records)
    seams = tuple(tuple(seam) for seam in structural['seams'])
    biased_records = execute_biased_calibration(seams)
    biased = recover_structure(biased_records)
    if not biased['rank_three_structural_gram_certified']:
        raise AssertionError('biased calibration failed exact structural recovery')
    if biased['gram_qsqrt5'] != uniform['gram_qsqrt5']:
        raise AssertionError('seam exposure changed the structural Gram')
    if biased['laplacian_exact'] != uniform['laplacian_exact']:
        raise AssertionError('seam exposure changed the structural response')
    shot_records = execute_single_shot_calibration(seams)
    shots = recover_from_shots(shot_records)
    if shots['gram_qsqrt5'] != uniform['gram_qsqrt5']:
        raise AssertionError('finite native shots failed exact structural recovery')
    laplacian = [[F(x) for x in row] for row in structural['laplacian_exact']]
    missing = recover_structure(ng.execute_calibration({seams[0]: F(1)}))
    if missing['rank_three_structural_gram_certified']:
        raise AssertionError('one exposed seam cannot inherit the carrier certificate')
    doubled_weights = {seam: F(1+int(seam == seams[0]), 31) for seam in seams}
    doubled_records = ng.execute_calibration(doubled_weights)
    actual_biased = ng.recover(doubled_records)
    for i in range(PORTS):
        for j in range(PORTS):
            edge_i = int(i == seams[0][0])-int(i == seams[0][1])
            edge_j = int(j == seams[0][0])-int(j == seams[0][1])
            expected = F(i == j)-(laplacian[i][j]+edge_i*edge_j)/62
            if F(actual_biased['transition_matrix_exact'][i][j]) != expected:
                raise AssertionError('weighted source response does not match rank-one perturbation')
    return {
        'schema': 'oph.relative_observer_geometry.intrinsic.v1',
        'sources': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'single_carrier_recovery': biased,
        'same_structural_gram_as_uniform_protocol_exact': True,
        'finite_shot_recovery': {
            'native_branch_calls': len(shot_records),
            'shots_per_preparation_and_seam': 1,
            'probability_law_oracle_used': shots['probability_law_oracle_used'],
            'record_fields': sorted(ng.RECORD_FIELDS-{'probability'}),
            'same_structural_gram_exact': True,
            'conditional_current_variance_zero_under_native_model': True,
            'laboratory_noise_or_preparation_error_bounded': False},
        'biased_protocol': {'preparation_count': PORTS, 'native_branch_calls': len(biased_records),
            'preparation_dependent_positive_seam_exposure': True,
            'even_preparations_have_deterministic_coin_zero': True,
            'odd_preparations_have_coin_probabilities': ['1/3', '2/3'],
            'complete_exact_observation_laws_not_finite_shot_estimates': True},
        'actual_scheduler_slow_band_control': scheduling_control(laplacian, seams[0]),
        'full_response_metric': inverse_response_metric(laplacian),
        'product_wave_metric': product_wave_metric_control(),
        'incomplete_exposure_control': {'exposed_seams': len(missing['seams']),
            'rank_three_structural_gram_certified': missing['rank_three_structural_gram_certified']},
        'scope': {'federation_needed_for_finite_carrier_result': False,
            'uniform_schedule_needed_for_structural_gram_recovery': False,
            'coin_fairness_needed_for_even_preparation_readout': False,
            'arbitrary_positive_schedule_preserves_rank_three_native_slow_band': False,
            'actual_native_slow_band_requires_its_own_schedule_certificate': True,
            'complete_ensemble_probability_law_needed_for_even_preparation_structure': False,
            'known_preparation_faithful_records_and_all_seam_coverage_still_required': True,
            'observed_seam_labels_disclose_interacting_endpoint_incidence': True,
            'carrier_architecture_derived_from_self_reference': False,
            'slow_band_definition_selects_one_metric_not_all_observer_readouts': True,
            'rank_three_proves_three_dimensional_physical_bulk': False,
            'specific_federation_architecture_required_for_continuum_claim': False,
            'extension_or_limit_with_compatible_observable_laws_required_for_continuum_claim': True,
            'physical_length_time_or_lab_calibration_established': False},
        'all_pass': True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Single-carrier intrinsic geometry receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
