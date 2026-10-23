"""Exact native-port decoder for a declared Lorentzian tensor-charge map.

The directions come from the existing slow-response carrier; Lorentzian
lifting and physical stress/geometry measurements remain supplied interfaces.
"""
from fractions import Fraction as F

from .slow_geometry import (_add, _dot, _encode_vector, raw_frame,
                            registered_binding, source)


# Symmetric four-tensor coordinates; each off-diagonal entry occurs twice
# when contracted. The declared raw-coordinate metric is diag(-r^2,1,1,1),
# r^2=phi+2. Its normalized chart is the usual diag(-1,1,1,1).
TENSOR_COMPONENTS = ((0, 0), (0, 1), (0, 2), (0, 3), (1, 1),
                     (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))


def null_tomography_design():
    """Exact charges B(k_p,k_p) for the CONDITIONAL lift k_p=(1,raw_p).

    The carrier supplies raw_p, not a Lorentzian physical interpretation or
    a stress measurement. This map tests the finite directional requirement
    of EinsteinBranch.Tensor on that declared lift.
    """
    c = source()
    return [[c.q5_mul(c.q5(1 if i == j else 2), c.q5_mul(k[i], k[j]))
             for i, j in TENSOR_COMPONENTS]
            for axis in raw_frame() for k in [(c.q5(1),)+axis]]


def _null_tomography_decoder():
    """Closed inverse on the B00=0 slice, from exact icosahedral moments."""
    c = source()
    inverse_r2 = c.q5_inv(c.q5_add(c.GOLDEN, c.q5(2)))
    inverse_r4 = c.q5_mul(inverse_r2, inverse_r2)
    rows = []
    for i, j in TENSOR_COMPONENTS[1:]:
        if i == 0:
            rows.append([c.q5_mul(c.q5(F(1, 8)), c.q5_mul(inverse_r2, u[j-1]))
                         for u in raw_frame()])
        else:
            rows.append([c.q5_sub(
                c.q5_mul(c.q5(F(5, 8)), c.q5_mul(inverse_r4, c.q5_mul(u[i-1], u[j-1]))),
                c.q5_mul(c.q5(F(int(i == j), 8)), inverse_r2)) for u in raw_frame()])
    return rows


def decode_null_charges(charges):
    """Return the unique B00=0 representative, rejecting incompatible data.

    Twelve arbitrary exact readings need not be charges of one symmetric
    tensor. No projection or least-squares fit replaces this compatibility
    test. Entries are exact (a,b) pairs representing a+b*sqrt(5).
    """
    if (not isinstance(charges, (tuple, list)) or len(charges) != 12
            or any(not isinstance(q, (tuple, list)) or len(q) != 2
                   or any(type(v) not in (int, F) for v in q) for q in charges)):
        raise ValueError('twelve exact Q(sqrt5) charge pairs are required')
    values = tuple((F(a), F(b)) for a, b in charges)
    candidate = (source().q5(0),)+tuple(_dot(row, values) for row in _null_tomography_decoder())
    if tuple(_dot(row, candidate) for row in null_tomography_design()) != values:
        raise ValueError('the twelve readings are incompatible with a common symmetric tensor')
    return candidate


def _max_row_sum(matrix):
    """The exact induced infinity norm, also valid on real input errors."""
    c = source()
    largest = c.q5(0)
    for row in matrix:
        total = c.q5(0)
        for value in row:
            total = c.q5_add(total, value if c.q5_sign(value) >= 0 else c.q5_sub(c.q5(0), value))
        if c.q5_sign(c.q5_sub(total, largest)) > 0:
            largest = total
    return largest


def null_tomography_certificate():
    """Finite exact inverse/kernel certificate on the native twelve ports.

    All matrix identities are checked over Q(sqrt5); hence the displayed
    inverse also proves the statement for arbitrary real tensor components.
    This is an executable algebraic certificate, not an additional Lean proof.
    """
    c = source()
    zero, one = c.q5(0), c.q5(1)
    frame = raw_frame()
    r2 = c.q5_add(c.GOLDEN, c.q5(2))
    if c.q5_sign(r2) <= 0:
        raise AssertionError('the declared metric requires a positive spatial radius')
    inverse_r2 = c.q5_inv(r2)
    gram = c.intrinsic_gram_exact()
    if any(c.q5_mul(_dot(u, v), inverse_r2) != gram[i][j]
           for i, u in enumerate(frame) for j, v in enumerate(frame)):
        raise AssertionError('null frame differs from the native response Gram')
    design, decoder = null_tomography_design(), _null_tomography_decoder()
    sliced = [row[1:] for row in design]
    identity9 = [[c.q5(int(i == j)) for j in range(9)] for i in range(9)]
    if c.q5_matmul(decoder, sliced) != identity9:
        raise AssertionError('native null-charge decoder is not an exact left inverse')
    metric = [c.q5_sub(zero, r2), zero, zero, zero, one, zero, zero, one, zero, one]
    if any(_dot(row, metric) != zero for row in design):
        raise AssertionError('declared carrier lifts are not null')
    # Rank >=9 by the left inverse, <=9 by this nonzero metric kernel.
    # The full inverse verifies that there are no other invisible components.
    quotient = [[(inverse_r2 if i in (3, 6, 8) else zero)]+identity9[i]
                for i in range(9)]
    if c.q5_matmul(decoder, design) != quotient:
        raise AssertionError('full tensor decoder does not quotient by precisely the metric')
    projector = c.q5_matmul(sliced, decoder)
    partner = [next(j for j, v in enumerate(frame) if _add(u, v) == (zero,)*3)
               for u in frame]
    even = [[c.q5(F(int(i == j)+int(partner[i] == j), 2)) for j in range(12)]
            for i in range(12)]
    slow = c.slow_band_projector_exact()
    if projector != [[c.q5_add(even[i][j], slow[i][j]) for j in range(12)] for i in range(12)]:
        raise AssertionError('tensor charges do not match the native constant/slow/quadrupole bands')
    residual = [[c.q5_sub(c.q5(int(i == j)), projector[i][j]) for j in range(12)]
                for i in range(12)]
    lap = [[c.q5(int(x)) for x in row] for row in c.laplacian()]
    if c.q5_matmul(lap, residual) != [[c.q5_mul(c.q5(5, 1), x) for x in row] for row in residual]:
        raise AssertionError('incompatible charges are not precisely the native fast band')
    merged = c.q5_matmul(even, design)
    # Six independent spatial components survive antipodal averaging; the
    # three mixed time/space entries vanish, in addition to the metric kernel.
    spatial_decoder = decoder[3:]
    if any(row[j] != zero for row in merged for j in (1, 2, 3)):
        raise AssertionError('antipodal averaging should erase time/space mixed charges')
    identity6 = [[c.q5(int(i == j)) for j in range(6)] for i in range(6)]
    if c.q5_matmul(spatial_decoder, [row[4:] for row in merged]) != identity6:
        raise AssertionError('antipodal merged charges lost a spatial tensor component')
    return {'schema': 'oph.common_source.native_null_tomography.v1',
            'registered_binding': registered_binding(),
            'tensor_component_order': [list(pair) for pair in TENSOR_COMPONENTS],
            'declared_raw_metric_diagonal_Qsqrt5': _encode_vector((metric[0], one, one, one)),
            'design_Qsqrt5': [_encode_vector(row) for row in design],
            'B00_zero_decoder_Qsqrt5': [_encode_vector(row) for row in decoder],
            'exact_design_rank': 9, 'exact_kernel_dimension': 1,
            'exact_charge_compatibility_constraints': 3,
            'charge_image': 'constant + native slow rank3 + antipodal-even tracefree rank5',
            'incompatible_reading_band': 'native fast rank3, Laplacian eigenvalue 5+sqrt(5)',
            'antipodal_merged_rank': 6, 'antipodal_merged_kernel_dimension': 4,
            'canonical_tensor_max_error_per_charge_max_error_Qsqrt5': c.q5_str(_max_row_sum(decoder)),
            'compatibility_residual_max_error_per_charge_max_error_Qsqrt5': c.q5_str(_max_row_sum(residual)),
            'checks': {'decoder_left_inverse_exact': True, 'full_kernel_is_metric_exact': True,
                       'charge_image_matches_native_bands_exact': True,
                       'antipodal_loss_exact': True},
            'scope': {'conditional_lorentz_lift': True,
                      'physical_null_directions_derived': False,
                      'directional_stress_or_geometry_measurements_constructed': False,
                      'einstein_equation_derived': False}, 'all_pass': True}


