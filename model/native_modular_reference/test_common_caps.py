"""Exact common-reference characterization and native mixed-protocol checks."""
from fractions import Fraction as F
from itertools import combinations
from copy import deepcopy
from decimal import Decimal, localcontext
import hashlib
import json

import numpy as np
import pytest

from . import common_caps as c


@pytest.mark.parametrize("n", (3, 4, 5, 12))
def test_constructive_factorization_and_independent_spectator_ratios(n):
    relative = tuple(F(i+2, i+1) for i in range(n))
    weights = {(i, j): F(3, 7)*relative[i]*relative[j] for i, j in combinations(range(n), 2)}
    row = c.factor_characterization(n, weights)
    assert row["all_bare_two_port_algebras_invariant"]
    scale = F(row["factor_scale"]); factors = tuple(map(F, row["relative_factors"]))
    assert all(scale*factors[i]*factors[j] == w for (i, j), w in weights.items())
    for i, j in combinations(range(n), 2):
        ratios = {weights[tuple(sorted((i, k)))]/weights[tuple(sorted((j, k)))]
                  for k in range(n) if k not in (i, j)}
        assert ratios == {relative[i]/relative[j]}
    if n >= 4:
        changed = dict(weights); changed[0, 1] *= F(101, 100)
        assert not c.factor_characterization(n, changed)["all_bare_two_port_algebras_invariant"]


def test_arbitrary_three_port_weights_and_exact_native_countercontrols():
    row = c.factor_characterization(3, {(0, 1): F(2), (0, 2): F(7), (1, 2): F(11)})
    assert row["all_bare_two_port_algebras_invariant"]
    for sigma in (c.regional.source()[2], c.cut.completed_palm(c.CAPS[0]), c.reference()):
        row = c.factor_characterization(12, c.native_weights(sigma))
        assert not row["all_bare_two_port_algebras_invariant"]
        quad = row["noninvariance_four_port_witness"]
        assert len(set(quad["matching_products"])) > 1
    assert c.factor_characterization(12, c.native_weights((F(1, 66),)*66))["all_bare_two_port_algebras_invariant"]


def test_nonuniform_compatible_product_is_not_selected_by_connected_native_repair():
    # Mathematical compatible-clock control, explicitly not a sourced field.
    relative = (F(2), *(F(1) for _ in range(11)))
    states = c.regional.source()[0].states
    weights = tuple(relative[i]*relative[j] for x in states
                    for i, j in [tuple(p for p, bit in enumerate(x) if bit)])
    prior = np.array([float(w/sum(weights)) for w in weights])
    kernel = c.regional.source()[0]
    after = np.array(kernel.counts).T@prior/60
    assert np.linalg.norm(after-prior, 1) > .005
    for i, j in kernel.seams:
        if relative[i] != relative[j]:
            spectator = next(k for k in range(12) if k not in (i, j))
            assert relative[i]*relative[spectator] != relative[j]*relative[spectator]
            break
    else:
        raise AssertionError("nonconstant field has no native seam witness")


def test_actual_mixed_cap_source_and_all_1024_common_algebras():
    row = c.exact_audit()
    assert row["common_reference_cap_algebras_checked"] == 1024
    assert row["protected_intersection_modular_action"] == "identity"
    core, weights, schedule = c.cap_schedule()
    assert core == (0, 1) and sum(weights.values()) == 1
    assert schedule["native_seam_union"] == 20
    assert schedule["uniform_component_weight"] == "2/3"
    assert schedule["gap_lower_per_attempt"] == "1/600"
    assert all(e not in c.CORE for edge in weights for e in edge)
    assert set(c.reference()) == {F(17, 1125), F(23, 1500), F(1, 75)}


def test_mixed_native_complete_record_branch_covariance_and_cap_marginal_control():
    sigma = c.reference()
    kernel = c.regional.source()[0]
    index = {x: i for i, x in enumerate(kernel.states)}
    # Every branch commutes with one same reference, not different cap Gibbs
    # matrices. Labels C/E can split branches but cannot change this identity.
    for seam in c.cap_schedule()[1]:
        for coin in (0, 1):
            for i, x in enumerate(kernel.states):
                y, _ = c.cut.native.native(x, seam, coin)
                assert sigma[i] == sigma[index[y]]
    # The other cap's repair can change the first cap marginal. Only the
    # intersection is preserved throughout the mixed preparation.
    counts, scale = c.mixed_counts("C")
    cap_unit = None
    for unit in c.regional.basis_supports(c.CAPS[0], completed=False):
        vector = np.zeros(4356, np.int64)
        for i, j in unit:
            vector[66*i+j] = 1
        if np.any(counts.T@vector != scale*vector):
            cap_unit = unit; break
    assert cap_unit is not None
    for unit in c.regional.basis_supports(c.CORE, completed=False):
        vector = np.zeros(4356, np.int64)
        for i, j in unit:
            vector[66*i+j] = 1
        np.testing.assert_array_equal(counts.T@vector, scale*vector)


def test_unprotected_native_control_has_connected_tracial_selection():
    control = c.unprotected_control()
    assert control["schedule"]["protected_intersection"] == []
    assert control["schedule"]["native_seam_union"] == 28
    assert control["schedule"]["gap_lower_per_attempt"] == "1/1080"
    assert control["unique_fixed_binary_sector_state"] == "I66/66"


def test_finite_common_preparation_floor_and_reference_error():
    row = c.preparation_budget()
    assert row["mixed_cap_native_attempts_after_success"] == 18600
    assert row["actual_common_reference_modular_half_diamond_upper"] < F(99, 100000)
    assert row["minimum_work_eigenvalue_lower"] > F(1333, 100000)
    assert row["expected_native_preparation_attempts_upper"] < F(21042015, 1000)
    compact = c.preparation._compact(row)
    for key, value in row.items():
        if isinstance(value, F) and key.endswith("_lower"):
            assert F(compact[key]) <= value
        elif isinstance(value, F) and key.endswith("_upper"):
            assert F(compact[key]) >= value
    with pytest.raises(ValueError, match="faithful"):
        c.preparation_budget(attempts=0)


def test_one_complete_source_word_and_two_cap_Ramsey_readouts():
    s = F(16)
    bound = c.checked_bound(s)
    assert len(c.word(s)) == 2042
    assert bound["operator_error_upper"] == F(1, 10**9)
    outputs = [c.witness(i, s) for i in (0, 1)]
    assert [row["source_pulses"] for row in outputs] == [2056, 2068]
    for row in outputs:
        assert .98 < row["complete_bit_probabilities"][0] < .99
        assert abs(sum(row["complete_bit_probabilities"])-1) < 2e-11
    assert abs(outputs[0]["complete_bit_probabilities"][0]-outputs[1]["complete_bit_probabilities"][0]) < 2e-11
    rng = np.random.default_rng(447)
    probe = rng.normal(size=(512, 3))+1j*rng.normal(size=(512, 3))
    probe /= np.linalg.norm(probe)
    phase = np.array([float(w)**(-1j*16) for w in c.target_ratios()])[:, None]
    assert np.linalg.norm(c.matrix(s)@probe-phase*probe) < float(c.WORD_ERROR_GATE)
    # Full matrix certificate catches a terminal padding error, with no
    # native-shell projection or renormalization allowed.
    wrong = np.diag(phase[:, 0]); wrong[511, 511] += .125
    assert c.certify(s, wrong) >= F(1, 8)


def test_shared_clock_nontrivial_cap_action_and_trivial_intersection():
    sigma = c.reference()
    for cap in c.CAPS:
        ratios = {sigma[i]/sigma[j] for unit in c.regional.basis_supports(cap, completed=False) for i, j in unit}
        assert F(69, 68) in ratios and F(23, 20) in ratios
    for unit in c.regional.basis_supports(c.CORE, completed=False):
        assert {sigma[i]/sigma[j] for i, j in unit} == {F(1)}


@pytest.mark.parametrize("bad", ({(0, 1): F(1)}, {(0, 1): 1., (0, 2): F(1), (1, 2): F(1)},
                                   {(0, 1): F(0), (0, 2): F(1), (1, 2): F(1)}))
def test_factor_contract_rejects_missing_float_or_nonfaithful_weights(bad):
    with pytest.raises(ValueError, match="exact positive"):
        c.factor_characterization(3, bad)


@pytest.mark.parametrize("pair", ((False, True), (0., 1.), (F(0), F(1))))
def test_factor_characterization_rejects_noninteger_pair_key_aliases(pair):
    weights = {pair: F(1), (0, 2): F(1), (1, 2): F(1)}
    with pytest.raises(ValueError, match="exact positive"):
        c.factor_characterization(3, weights)


@pytest.mark.parametrize("bits,compact", (((False, 8), True), ((0., 8), True),
                                         ((np.int64(0), 8), True), ([0, 8], True),
                                         ((0, 8), 1), ((0, 8), "yes")))
def test_compact_mixer_does_not_bypass_original_bit_and_flag_types(bits, compact):
    with pytest.raises(ValueError, match="signed native pair"):
        c.mixer_word(1, bits, compact=compact)


@pytest.mark.parametrize("s,digest", (
    (F(0), "47f45544e2c22625f1e5375d4a1579b206732dd569c20daafa00184b4d5c6ffb"),
    (F(1), "b4ffa7a8322379289fc04467e8d416307283bdf5db47a362ad0cd9abb5081c51"),
    (F(16), "d134c2c35cad079bd235a783104c820fc8d2ee454106d6f2c4d6fcce687759be"),
))
def test_general_diagonal_compiler_preserves_retained_common_words(s, digest):
    # These are raw-word hashes from the pre-refactor fd068a2d4 producer.
    word = c.diagonal_word(s, c.target_ratios())
    assert hashlib.sha256(c.full.corner.canonical(word).encode()).hexdigest() == digest
    assert c.word(s) == word


@pytest.mark.parametrize("bad", (False, 0, -1, .5, float("nan"), float("inf"),
                                 complex(1), "1", np.int64(1)))
def test_diagonal_ratios_require_original_exact_positive_scalars(bad):
    values = (bad,)+(F(1),)*511
    for operation in (lambda: c.diagonal_word(F(1), values),
                      lambda: c.certify_diagonal(F(1), values, np.eye(512))):
        with pytest.raises(ValueError, match="512 original exact positive"):
            operation()


@pytest.mark.parametrize("values", ([F(1)]*512, (F(1),)*511, (F(1),)*513, None))
def test_diagonal_ratios_require_a_complete_original_tuple(values):
    with pytest.raises(ValueError, match="512 original exact positive"):
        c.diagonal_word(F(1), values)
    with pytest.raises(ValueError, match="512 original exact positive"):
        c.certify_diagonal(F(1), values, np.eye(512))


@pytest.mark.parametrize("bad", (True, 1., F(-1), F(0), "1", float("nan")))
def test_logarithm_validation_precedes_its_cache(bad):
    c.logarithm_interval(F(1))
    with pytest.raises(ValueError, match="original exact positive"):
        c.logarithm_interval(bad)


@pytest.mark.parametrize("ratio", (F(1, 1024), F(2, 225), F(1, 36), F(2, 45),
                                   F(8, 3), F(3), F(2)**100))
def test_range_reduced_logarithm_encloses_independent_high_precision_value(ratio):
    lo, hi = c.logarithm_interval(ratio)
    with localcontext() as context:
        context.prec = 100
        value = (Decimal(ratio.numerator)/Decimal(ratio.denominator)).ln()
        centre, rounding = F(value), F(1, 10**95)
    assert lo <= centre-rounding <= centre+rounding <= hi
    assert hi-lo < F(1, 10**27)


@pytest.mark.parametrize("ratio", (F(2, 225), F(1, 36), F(2, 45), F(8, 3)))
@pytest.mark.parametrize("parameter", (F(-16), F(16)))
def test_scaled_phase_enclosure_against_independent_long_decimal_series(ratio, parameter):
    real, imag, radius = c.phase_enclosure(parameter, ratio)
    with localcontext() as context:
        context.prec = 160
        theta = -int(parameter)*(Decimal(ratio.numerator)/Decimal(ratio.denominator)).ln()
        re, im, pr, pi = Decimal(1), Decimal(0), Decimal(1), Decimal(0)
        # Direct, unscaled Taylor series is independent of dyadic squaring.
        for k in range(1, 600):
            pr, pi = -theta*pi/k, theta*pr/k
            re += pr
            im += pi
        error = ((Decimal(real.numerator)/Decimal(real.denominator)-re)**2
                 +(Decimal(imag.numerator)/Decimal(imag.denominator)-im)**2).sqrt()
        assert F(error)+F(1, 10**100) < radius < F(1, 10**27)


def test_general_phase_keeps_zero_exact_and_integer_ratios_are_admitted():
    assert c.phase_enclosure(F(0), F(2, 225)) == (F(1), F(0), F(0))
    assert c.diagonal_word(F(16), (1,)*512) == c.word(F(0))


def test_general_diagonal_word_certifies_raw_local_reference_weights_and_padding():
    ratios = [F(1)]*512
    weights = (F(2, 225), F(1, 36), F(2, 45))
    for index, code in enumerate(c.full.support()[0]):
        ratios[code] = weights[index % 3]
    ratios = tuple(ratios)
    matrix = c.full.compiler.execute(c.diagonal_word(F(16), ratios), 9)
    assert c.certify_diagonal(F(16), ratios, matrix) < c.WORD_ERROR_GATE
    target = np.diag(np.exp(-16j*np.log(np.array(list(map(float, ratios))))))
    assert np.linalg.norm(matrix-target, "fro") < 1e-10
    wrong = matrix.copy()
    wrong[511, 511] += .125
    assert c.certify_diagonal(F(16), ratios, wrong) > F(12, 100)


@pytest.mark.parametrize("kind", ("masked", "boolean", "object", "nonfinite", "shape", "list"))
def test_general_diagonal_certificate_rejects_malformed_original_arrays(kind):
    array = np.eye(512)
    if kind == "masked":
        array = np.ma.array(array, mask=False)
    elif kind == "boolean":
        array = array.astype(bool)
    elif kind == "object":
        array = array.astype(object)
    elif kind == "nonfinite":
        array[0, 0] = float("nan")
    elif kind == "shape":
        array = array[:-1]
    else:
        array = array.tolist()
    with pytest.raises(ValueError):
        c.certify_diagonal(F(1), (F(1),)*512, array)


@pytest.fixture(scope="module")
def retained_common_receipt():
    return json.loads(c.OUTPUT.read_text(encoding="utf-8"))


@pytest.mark.parametrize("side", ("expected", "actual", "both"))
def test_negative_receipt_error_cannot_hide_in_portable_tolerance(retained_common_receipt, side):
    wrong = deepcopy(retained_common_receipt)
    wrong["cap_readouts"][0]["complete_record_half_l1_diagnostic"] = -1e-10
    c.full.corner.compare_receipt(wrong, retained_common_receipt)
    expected = wrong if side in ("expected", "both") else retained_common_receipt
    actual = wrong if side in ("actual", "both") else retained_common_receipt
    with pytest.raises(ValueError, match="nonnegative complete-record"):
        c.compare_receipt(expected, actual)


@pytest.mark.parametrize("mutation", ("negative_probability", "excess_mass", "nan", "boolean",
                                     "missing_outcome", "missing_readout", "negative_bound",
                                     "float_bound", "tuple_outcomes", "wrong_schema"))
def test_identical_malformed_receipts_fail_local_domains(retained_common_receipt, mutation):
    wrong = deepcopy(retained_common_receipt)
    row = wrong["cap_readouts"][0]
    if mutation == "negative_probability":
        row["complete_bit_probabilities"][1] = -1e-10
    elif mutation == "excess_mass":
        row["complete_bit_probabilities"] = [1., .001]
    elif mutation == "nan":
        row["complete_record_half_l1_diagnostic"] = float("nan")
    elif mutation == "boolean":
        row["complete_record_half_l1_diagnostic"] = False
    elif mutation == "missing_outcome":
        row["complete_bit_probabilities"].pop()
    elif mutation == "missing_readout":
        wrong["cap_readouts"].pop()
    elif mutation == "negative_bound":
        row["whole_sequence_half_diamond_upper"] = "-1/100"
    elif mutation == "float_bound":
        row["whole_sequence_half_diamond_upper"] = 1e-9
    elif mutation == "tuple_outcomes":
        row["complete_bit_probabilities"] = tuple(row["complete_bit_probabilities"])
    else:
        wrong["schema"] = "different"
    with pytest.raises(ValueError):
        c.compare_receipt(wrong, wrong)


def test_small_positive_receipt_roundoff_remains_portable_and_unmodified(retained_common_receipt):
    other = deepcopy(retained_common_receipt)
    row = other["cap_readouts"][0]
    row["complete_bit_probabilities"][0] += 1e-13
    row["complete_bit_probabilities"][1] -= 1e-13
    row["complete_record_half_l1_diagnostic"] += 1e-13
    before = deepcopy(other)
    c.compare_receipt(retained_common_receipt, other)
    assert other == before


def test_common_receipt_custody_includes_executed_source_and_certificate_modules():
    required = {"model/native_modular_reference/full_flow.py", "model/observer/source_scalar/compiler.py",
                "model/common_source/exact_output.py", "model/common_source/output_certificate.py",
                "model/observer/native.py", "model/observer/quantum_mixing.py",
                "model/record_laws/burn_in.py", "model/observer/source_discrimination/circuit.py"}
    assert required <= set(c.SOURCE_PATHS)
    assert all((c.ROOT/path).is_file() for path in c.SOURCE_PATHS)
