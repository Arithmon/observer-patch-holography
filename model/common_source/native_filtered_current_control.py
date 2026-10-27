"""An operator-filtered coherent current candidate on the protected pair.

The Hamiltonian control H_p=P d_p P is a declared control choice. Here d_p
comes from all native mixed exchange responses and P is the incidence class
projector. Actual oriented order-five rotations supply their Cayley
generators K_p. The twelve F_p=K_p+iH_p span the existing u(3)+so(3) current
tangent, with a computed invertible control remapping. The individual port
drives are different from native_sm_current's port drives.

This is an exact logical interface, not a source pulse implementation or a
selection result for unchanged C/E repair. Filtering a Hamiltonian is not
the family state-reset cooling channel. That channel erases the relative
central phase even when all its classical records are retained. The spin
lift, exterior register and Higgs attachment below keep their existing
conditional meaning; no physical spectrum or all-intervention completeness
claim is made.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache

import sympy as s

from . import native_gauge_object as objects
from . import native_persistent_response as response
from . import native_sm_current as declared


def _clean(matrix):
    return s.ImmutableMatrix(s.Matrix(matrix).applyfunc(lambda x: s.expand(s.radsimp(x))))


def _zero(matrix):
    return not any(_clean(matrix))


def _flat(matrix):
    return s.Matrix([s.re(x).expand() for x in matrix]+[s.im(x).expand() for x in matrix])


def _columns(matrices):
    return s.Matrix.hstack(*(_flat(x) for x in matrices))


def _rank(matrix):
    return matrix.to_DM(extension=True).rank()


def _entries(matrix):
    return [[str(x) for x in matrix.row(i)] for i in range(matrix.rows)]


def _operator(matrix, dimension=6):
    if (not isinstance(matrix, s.MatrixBase) or matrix.shape != (dimension, dimension)
            or any(x.has(s.Float) or x.free_symbols or x.is_finite is not True for x in matrix)):
        raise ValueError('a finite exact numerical operator on the declared domain is required')
    return _clean(matrix)


def _controls(values):
    if type(values) is not tuple or len(values) != 12:
        raise ValueError('a tuple of twelve exact real algebraic controls is required')
    result = []
    for value in values:
        if type(value) is int:
            value = s.Integer(value)
        elif type(value) is Fraction:
            if (type(value.numerator) is not int or type(value.denominator) is not int
                    or value.denominator <= 0):
                raise ValueError('original exact Fraction controls are required')
            value = s.Rational(value.numerator, value.denominator)
        elif not isinstance(value, s.Expr):
            raise ValueError('original exact real algebraic controls are required')
        if (value.has(s.Float) or value.free_symbols or value.is_real is not True
                or value.is_algebraic is not True or value.is_finite is not True):
            raise ValueError('finite exact real algebraic controls are required')
        result.append(value)
    return tuple(result)


def _combination(controls, basis):
    return _clean(sum((a*b for a, b in zip(controls, basis, strict=True)), s.zeros(6)))


def _rotation(row):
    record = row['rotation3']
    return s.ImmutableMatrix(s.SparseMatrix(*record['shape'], {
        (i, j): s.sympify(value) for i, j, value in record['entries']}))


@lru_cache(maxsize=1)
def exact_architecture():
    """Derive the filtered port drives before comparing the existing current.

    All matrices use the independently certified orthonormal native
    geometric frame. Columns of declared_from_filtered express each F_p
    in native_sm_current's port basis; its inverse maps the other way.
    This remaps simultaneous infinitesimal controls, not products of
    finite noncommuting port pulses.
    """
    geometry = objects.geometry()
    frames = objects.geometric_intertwiners()
    T = _clean(s.Matrix.hstack(*(frames[name]['frame']/s.sqrt(frames[name]['metric'])
                               for name in ('three', 'three_prime'))))
    coordinates = response.geometric_coordinates()
    chart = _clean(coordinates['left_inverse']*T)
    inverse = _clean(T.T*coordinates['frame'])
    if (not _zero(T.T*T-s.eye(6)) or not _zero(inverse*chart-s.eye(6))
            or not _zero(coordinates['frame']*chart-T)):
        raise ArithmeticError('the native geometric coordinate transport failed')

    rows = objects.geometric_projectors()['rotation_rows']
    rotations = tuple(_rotation(row) for row in rows)
    actions = tuple(_clean(s.diag(R, R.xreplace({s.sqrt(5): -s.sqrt(5)}))) for R in rotations)
    edges = set(geometry['seams'])
    near = tuple(i for i, row in enumerate(rows) if row['order'] == 5
                 and sum(tuple(sorted((p, row['port_action'][p]))) in edges for p in range(12)) == 10)
    if len(near) != 12:
        raise ArithmeticError('the incidence-defined near-fivefold class changed')
    C = _clean(sum((actions[i] for i in near), s.zeros(6)))
    high, low = 2+2*s.sqrt(5), 2-2*s.sqrt(5)
    P = _clean((C-low*s.eye(6))/(high-low))
    Q = _clean(s.eye(6)-P)
    if not _zero(P*P-P) or s.trace(P) != 3 or not _zero(C*P-high*P):
        raise ArithmeticError('the native positive class projector failed')

    # These are the exact native responses, not a retained desired-current
    # fixture. A compressed response alone need not preserve the full
    # object: the operator filter is an additional Hamiltonian control.
    _, reduced, _ = response.response_matrices(0, 1)
    derivatives, filtered, cayley, chosen = [], [], [], []
    for port in range(12):
        raw = sum((2*h for pair, h in zip(response.PORT_PAIRS, reduced, strict=True)
                   if port in pair), s.zeros(6))
        d = _clean(inverse*raw*chart)
        H = _clean(P*d*P)
        candidates = []
        v = geometry['frame'][port, :].T
        for index in near:
            R = rotations[index]
            if not _zero(R*v-v):
                continue
            skew = (R-R.T)/2
            axial = s.Matrix([skew[2, 1], skew[0, 2], skew[1, 0]])
            if s.expand((axial.T*v)[0]).is_positive is True:
                candidates.append(index)
        if len(candidates) != 1:
            raise ArithmeticError('the oriented native order-five port rotation is not unique')
        index = candidates[0]
        K = _clean((actions[index]-s.eye(6))*(actions[index]+s.eye(6)).inv())
        if (not _zero(d-d.H) or not _zero(H-H.H) or not _zero(K+K.H)
                or not _zero(K*H-H*K)):
            raise ArithmeticError('the filtered port drive lost its commuting skew/phase form')
        derivatives.append(d); filtered.append(H); cayley.append(K); chosen.append(index)
    generators = tuple(_clean(K+s.I*H) for K, H in zip(cayley, filtered, strict=True))

    # Admission verifies every original entry. Invertibility and the
    # reverse reconstruction establish equality of tangents, not equality
    # of port primitives or availability of a rational source pulse tape.
    forward = _clean(s.Matrix.hstack(*(s.Matrix(declared.admit_generator(g)) for g in generators)))
    if _rank(forward) != 12:
        raise ArithmeticError('the filtered twelve-port tangent is not faithful')
    backward = _clean(forward.to_DM(extension=True).to_field().inv().to_Matrix())
    old = declared.exact_architecture()
    for j, g in enumerate(generators):
        if not _zero(_combination(tuple(forward[:, j]), old['generators'])-g):
            raise ArithmeticError('the filtered-to-declared primitive remapping failed')
    for j, g in enumerate(old['generators']):
        if not _zero(_combination(tuple(backward[:, j]), generators)-g):
            raise ArithmeticError('the declared-to-filtered primitive remapping failed')
    return {'geometric_frame': T, 'port_frame': geometry['frame'], 'class_sum': C,
            'P': P, 'Q': Q, 'raw_derivatives': tuple(derivatives),
            'filtered_hamiltonians': tuple(filtered), 'cayley_generators': tuple(cayley),
            'generators': generators, 'rotation_rows': rows, 'rotation_actions': actions,
            'oriented_rotation_indices': tuple(chosen),
            'declared_from_filtered': forward, 'filtered_from_declared': backward,
            'native_to_fixture_port': old['native_to_fixture_port']}


def generator(controls):
    """Exact logical generator for simultaneous filtered-port controls."""
    return _combination(_controls(controls), exact_architecture()['generators'])


def to_declared_controls(controls):
    return tuple(_clean(exact_architecture()['declared_from_filtered']*s.Matrix(_controls(controls))))


def from_declared_controls(controls):
    return tuple(_clean(exact_architecture()['filtered_from_declared']*s.Matrix(_controls(controls))))


def declared_rational_controls(controls):
    """Bridge only when the exact remap fits the old rational-control API.

    Passing algebraic coefficients to the old rational-only source request
    interface would silently change its domain. No rounding is permitted.
    """
    values = to_declared_controls(controls)
    if any(x.is_Rational is not True for x in values):
        raise ValueError('the exact remap requires algebraic controls outside the declared rational source API')
    return tuple(Fraction(int(x.p), int(x.q)) for x in values)


@lru_cache(maxsize=1)
def exact_certificate():
    d = exact_architecture(); generators = d['generators']
    brackets = []
    for p in range(12):
        for q in range(p):
            bracket = _clean(generators[p]*generators[q]-generators[q]*generators[p])
            old = s.Matrix(declared.admit_generator(bracket))
            coefficients = _clean(d['filtered_from_declared']*old)
            if not _zero(_combination(tuple(coefficients), generators)-bracket):
                raise ArithmeticError('a filtered current bracket left the declared tangent')
            brackets.append({'ports': [p, q], 'filtered_coefficients': [str(x) for x in coefficients]})
    checks = 0
    for row, U in zip(d['rotation_rows'], d['rotation_actions'], strict=True):
        for p in range(12):
            if not _zero(U*generators[p]*U.T-generators[row['port_action'][p]]):
                raise ArithmeticError('filtered primitive covariance under native port transport failed')
            checks += 1
    central = s.I*d['P']
    visible = _columns(tuple(_clean(g-s.trace(g)*s.eye(6)/6) for g in generators))
    if _rank(visible) != 12:
        raise ArithmeticError('the relative central current is not visible on the logical M6')
    return {'logical_current_dimension': 12, 'logical_EndM6_response_rank': 12,
            'all_primitive_admission_and_reverse_remapping_checks': 24,
            'commutator_closure': brackets, 'native_port_covariance_identities': checks,
            'positive_class_projector': _entries(d['P']),
            'filtered_port_generators': [_entries(g) for g in generators],
            'oriented_order_five_port_actions': [d['rotation_rows'][i]['port_action']
                                                for i in d['oriented_rotation_indices']],
            'declared_from_filtered': _entries(d['declared_from_filtered']),
            'filtered_from_declared': _entries(d['filtered_from_declared']),
            'native_to_fixture_port': list(d['native_to_fixture_port']),
            'relative_center_filtered_controls': [str(x) for x in
                d['filtered_from_declared']*s.Matrix(declared.admit_generator(central))],
            'primitive_port_drives_equal_declared_primitives': all(_zero(a-b) for a, b in zip(
                generators, declared.exact_architecture()['generators'], strict=True)),
            'remapping_can_require_irrational_controls': any(not x.is_Rational for x in d['declared_from_filtered'])}


def state_cooling(operator):
    """Actual reset channel restricted to the geometric six-dimensional input.

    This linear operator map is a countercontrol, not the Hamiltonian filter.
    Its result has support in P and contains no continuing Q triplet.
    """
    value = _operator(operator); d = exact_architecture(); P, Q = d['P'], d['Q']
    return _clean(P*value*P+s.trace(Q*value)*P/3)


def state_cooling_adjoint(operator):
    value = _operator(operator); d = exact_architecture(); P, Q = d['P'], d['Q']
    return _clean(P*value*P+s.trace(P*value)*Q/3)


@lru_cache(maxsize=1)
def cooling_countercontrol():
    from . import native_family_selection as family
    d = exact_architecture(); f = family.geometry(); T = d['geometric_frame']
    if (not _zero(f['P']*T-T*d['P']) or not _zero(f['P']*f['Q'])
            or not _zero(f['P']*f['P']-f['P'])):
        raise ArithmeticError('the state-reset countercontrol differs from actual family cooling')
    # The preserve Kraus map has right support P. Every reset Kraus map
    # |a_port><color|Q has right support Q, by QP=0. Thus K exp(itP)
    # equals exp(it)K or K, branch by branch, for every classical outcome.
    # The full recorded instrument cannot detect the relative phase.
    units = tuple(s.eye(6)[:, i]*s.eye(6)[j, :] for i in range(6) for j in range(6))
    central = s.I*d['P']
    for value in units:
        if (not _zero(state_cooling(central*value-value*central))
                or not _zero(central*state_cooling(value)-state_cooling(value)*central)):
            raise ArithmeticError('the state reset unexpectedly retained the relative central response')
    old = declared.exact_architecture()['generators']
    cooled_response = tuple(s.Matrix.vstack(*(
        _flat(g*state_cooling(x)-state_cooling(x)*g) for x in units)) for g in old)
    rank = _rank(s.Matrix.hstack(*cooled_response))
    if rank != 8:
        raise ArithmeticError('the continuing cooled work-state response rank changed')
    return {'actual_family_cooling_restricted_to_same_geometric_pair': True,
            'state_reset_is_hamiltonian_operator_filter': False,
            'unread_map': 'P rho P + Tr(Q rho) P/3',
            'adjoint_map': 'P B P + Tr(P B) Q/3',
            'logical_preserve_and_reset_records': 1+12*36,
            'central_phase_invisible_in_each_classical_cooling_record': True,
            'all36_geometric_input_matrix_units_checked': True,
            'continuing_work_state_response_rank': rank,
            'full_recorded_instrument_rank_claimed': False,
            'quantum_coherent_record_processor_included': False}


@lru_cache(maxsize=1)
def theorem_interfaces():
    """Transport the supplied spin/exterior and Higgs interfaces exactly.

    The invertible tangent remapping transports the entire invariant-tensor
    equation system, including its kernel dimension. It does not select a
    spin cover, chirality, copy count, physical statistics or Yukawa values.
    """
    from . import native_matter_attachment as matter
    from . import native_matter_theorems as finite
    d = exact_architecture(); r = matter.representation(); C = d['declared_from_filtered']
    lifted, exterior, selected = [], [], []
    for p, g in enumerate(d['generators']):
        a = matter.spin_lift(g); b = matter.exterior_derivative(a)
        expected5 = sum((C[j, p]*r['lifted'][j] for j in range(12)), s.zeros(5))
        expected32 = sum((C[j, p]*r['full_currents'][j] for j in range(12)), s.zeros(32))
        if not matter.zero(a-expected5) or not matter.zero(b-expected32):
            raise ArithmeticError('the filtered current failed its same-tangent spin/exterior remapping')
        small = matter.clean(r['encoder'].T*b*r['encoder'])
        if not matter.zero(b*r['encoder']-r['encoder']*small):
            raise ArithmeticError('the supplied exterior matter sector is not invariant')
        lifted.append(a); exterior.append(b); selected.append(small)
    lines = matter.invariant_higgs_lines()
    checks = 0
    for line in lines:
        for g, a in zip(selected, lifted, strict=True):
            h = -a[3:, 3:].T if line['conjugate'] else a[3:, 3:]
            for k in range(2):
                tensors = line['components']
                equation = g.T*tensors[k]+tensors[k]*g+sum(
                    (h[j, k]*tensors[j] for j in range(2)), s.zeros(15))
                if not matter.zero(equation):
                    raise ArithmeticError('a supplied Higgs line failed a filtered primitive')
                checks += 1
    theorem = finite.run()
    return {'conditional_spin_lift_remapping_checks': 12,
            'conditional_exterior_lift_remapping_checks': 12,
            'all32_auxiliary_occupation_masks_retained': True,
            'selected_matter_dimension': len(r['selected']),
            'finite_admissible_component_masks': [x['component_mask'] for x in theorem['admissible_sectors']],
            'one_doublet_invariant_lines': len(lines), 'Higgs_tensor_substitution_checks': checks,
            'Higgs_kernel_dimensions_preserved_by_invertible_current_remapping': True,
            'supplied': ['spin cover and five-mode lift', 'exterior component menu and parity sector',
                         'nonempty chiral anomaly-free selection', 'one Higgs doublet attachment'],
            'physical_spin_statistics_or_family_count_derived': False,
            'group_path_spin_holonomy_selected_by_operator_filter': False,
            'Yukawa_coefficients_or_particle_masses_derived': False}


def run(*, include_theorem_interfaces=True):
    report = {'schema': 'oph.native_filtered_current_control.v1',
              'classification': 'declared operator-filtered coherent repair-control candidate',
              'exact': exact_certificate(), 'cooling_countercontrol': cooling_countercontrol(),
              'scope': {'Hamiltonian_operator_filter_is_an_added_control_choice': True,
                        'unchanged_C_E_selects_this_current': False,
                        'source_compilation_executed': False,
                        'source_pulse_or_record_history_certificate_provided': False,
                        'logical_control_domain': 'twelve simultaneous exact real algebraic port amplitudes',
                        'finite_noncommuting_pulse_products_identified_by_linear_remapping': False,
                        'complete_only_over_declared_logical_current_tangent': True,
                        'all_physically_accessible_interventions_classified': False,
                        'full_model_A1_complete_accessible_response_discharged': False,
                        'physical_gauge_redundancy_or_SM_spectrum_derived': False}}
    if include_theorem_interfaces:
        report['theorem_interfaces'] = theorem_interfaces()
    return report


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-theorem-interfaces', action='store_true')
    arguments = parser.parse_args()
    print(json.dumps(run(include_theorem_interfaces=not arguments.skip_theorem_interfaces),
                     indent=2, allow_nan=False))
