"""Independent full-matrix and reference-state tests of emitted source words."""
from copy import deepcopy
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import compiler as c


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)


def independent_source_generators():
    """The parent's closed geometric formula, without its response receipt."""
    phi, conjugate_phi = (1+math.sqrt(5))/2, (1-math.sqrt(5))/2
    base = [(0., s, t*phi) for s in (1, -1) for t in (1, -1)]
    conjugate_base = [(0., s, t*conjugate_phi) for s in (1, -1) for t in (1, -1)]
    vertices = base+[(y, z, x) for x, y, z in base]+[(z, x, y) for x, y, z in base]
    conjugates = conjugate_base+[(y, z, x) for x, y, z in conjugate_base]+[
        (z, x, y) for x, y, z in conjugate_base]

    def skew(v):
        x, y, z = np.asarray(v)/2
        return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])

    matrices = []
    for v, conjugate in zip(vertices, conjugates, strict=True):
        v = np.array(v)
        full = np.zeros((6, 6), complex)
        full[:3, :3] = skew(v)+1j*(np.outer(v, v)/2+(1/12-v@v/6)*np.eye(3))
        full[3:, 3:] = skew(conjugate)
        matrices.append(full)
    return np.array(matrices)


def assert_matrix(actual, expected, tolerance=2e-10):
    assert actual.shape == expected.shape
    assert np.linalg.norm(actual-expected, 2) < tolerance


def pauli_matrix(label):
    result = np.ones((1, 1), complex)
    for letter in reversed(label):
        result = np.kron(result, dict(I=I, X=X, Y=Y, Z=Z)[letter])
    return result


@pytest.mark.parametrize("kind", c.KINDS)
def test_full_m6_pulse_from_independent_source_formula(kind):
    small = dict(X=1j*X, Y=1j*Y, Z=1j*Z, I=1j*I, H=1j*(X+Z)/math.sqrt(2))
    target = np.zeros((6, 6), complex)
    if kind == "J":
        target = np.diag([1j, 1j, 1j, 0., 0., 0.])
    else:
        target[:2, :2] = small[kind]
    source = independent_source_generators()
    for angle in (0., 1/7, -.51, math.pi):
        generator = np.einsum("p,pij->ij", c.integrated_controls(kind, angle), source)
        assert_matrix(generator, angle*target, 3e-14)
        actual = c.native_pulse(kind, angle)
        assert_matrix(actual, expm(generator), 8e-14)
        assert_matrix(actual, expm(angle*target), 8e-14)
        assert_matrix(actual.conj().T@actual, np.eye(6), 8e-14)


def test_full_source_exponential_retains_a_resolved_off_code_coupling(monkeypatch):
    # Perturb the already-loaded source to distinguish full M6 exponentiation
    # from inserting its ideal logical target or projecting onto the code.
    generators, coefficients, exact = c._source()
    perturbed = generators.copy()
    controls = c.integrated_controls("X", .37)
    port = int(np.argmax(abs(controls)))
    perturbed[port, 0, 3] += 1e-7
    perturbed[port, 3, 0] -= 1e-7
    expected = expm(np.einsum("p,pij->ij", controls, perturbed))
    monkeypatch.setattr(c, "_source", lambda: (perturbed, coefficients, exact))
    c._native.cache_clear()
    try:
        actual = c.native_pulse("X", .37)
        assert abs(actual[3, 0]) > 1e-9
        assert_matrix(actual, expected, 2e-14)
    finally:
        c._native.cache_clear()


@pytest.mark.parametrize("d", (2, 4))
@pytest.mark.parametrize("coherent", (False, True))
def test_full_code_export_includes_complement_and_preserves_no_cross_code_coherence(d, coherent):
    operators = c.transfer_kraus(d, coherent)
    assert_matrix(sum(k.conj().T@k for k in operators), np.eye(6), 1e-15)
    for row in range(6):
        for column in range(6):
            e = np.zeros((6, 6), complex)
            e[row, column] = 1
            actual = sum(k@e@k.conj().T for k in operators)
            expected = np.zeros((d+1, d+1), complex)
            if row < d and column < d and (coherent or row == column):
                expected[row, column] = 1
            elif row == column and row >= d:
                expected[d, d] = 1
            assert_matrix(actual, expected, 1e-15)
    # A complemented input is retained as failure, with probability one.
    outside = np.eye(6)[:, d]
    probabilities = [np.linalg.norm(k@outside)**2 for k in operators]
    assert sum(probabilities) == 1
    assert np.allclose(sum(np.outer(k@outside, (k@outside).conj()) for k in operators),
                       np.diag([0]*d+[1]), rtol=0, atol=0)


@pytest.mark.parametrize("d", (2, 4))
def test_coherent_service_preserves_an_untouched_reference_and_dephasing_does_not(d):
    entangled = np.zeros((6, d), complex)
    entangled[:d] = np.eye(d)/math.sqrt(d)
    target = np.zeros((d+1, d), complex)
    target[:d] = np.eye(d)/math.sqrt(d)
    exact_density = np.outer(target.reshape(-1), target.reshape(-1).conj())
    for coherent in (False, True):
        vectors = [(k@entangled).reshape(-1) for k in c.transfer_kraus(d, coherent)]
        actual = sum(np.outer(v, v.conj()) for v in vectors)
        assert np.trace(actual) == pytest.approx(1, abs=2e-15)
        distance = np.abs(np.linalg.eigvalsh(actual-exact_density)).sum()/2
        assert distance == pytest.approx(0 if coherent else 1-1/d, abs=2e-15)


@pytest.mark.parametrize("control,target", ((0, 1), (1, 0), (0, 2), (2, 1)))
def test_phase_exact_cnot_respects_global_little_endian_bit_order(control, target):
    qubits = 3
    expected = np.zeros((8, 8), complex)
    for incoming in range(8):
        outgoing = incoming ^ ((1 << target) if incoming >> control & 1 else 0)
        expected[outgoing, incoming] = 1
    assert_matrix(c.execute(c.cnot(control, target), qubits), expected, 6e-15)


@pytest.mark.parametrize("label", ("I", "Y", "XZ", "YZ", "XIY", "XYZ", "ZZZ"))
def test_pauli_rotations_preserve_spectators_and_full_phase(label):
    for angle in (0., .29, -.41):
        expected = expm(1j*angle*pauli_matrix(label))
        assert_matrix(c.execute(c.pauli_rotation(label, angle), len(label)), expected, 2e-14)


@pytest.mark.parametrize("bits", (1, 2, 3, 4))
def test_walsh_diagonal_matches_complete_independent_target(bits):
    values = np.array([.21+math.sin(j*j+.3)/7 for j in range(1 << bits)])
    tape = c.compile_diagonal(values)
    # Serialization is lossless and no matrix is smuggled into the tape.
    assert json.loads(json.dumps(tape)) == tape
    assert all(set(event) == {"kind", "angle", "targets"} for event in tape)
    actual = c.execute(tape, bits)
    assert_matrix(actual, np.diag(np.exp(1j*values)), 5e-14)
    without_identity = [event for event in tape if event["kind"] != "I"]
    assert np.linalg.norm(c.execute(without_identity, bits)-actual, 2) > .1


@pytest.mark.parametrize("bits", (1, 2, 3))
def test_complete_real_orthogonal_basis_and_its_inverse(bits):
    size = 1 << bits
    matrix, _ = np.linalg.qr(np.random.default_rng(610+bits).normal(size=(size, size)))
    matrix[:, 0] *= -1  # Exercise an independent column sign/determinant.
    tape = c.compile_real_orthogonal(matrix)
    assert_matrix(c.execute(tape, bits), matrix)
    assert_matrix(c.execute(c.inverse(tape), bits), matrix.T)
    # An entangled input in arbitrary reference columns uses the same gate.
    reference = np.random.default_rng(91).normal(size=(size, 3)).astype(complex)
    reference += 1j*np.random.default_rng(92).normal(size=reference.shape)
    reference /= np.linalg.norm(reference)
    assert_matrix(c.execute(tape, bits, reference), matrix@reference)


@pytest.mark.parametrize("a,b,bits", ((0, 1, 2), (3, 2, 2), (7, 5, 3)))
def test_gray_rotation_has_correct_orientation_and_identity_on_unmatched_spectators(a, b, bits):
    theta = -.37
    expected = np.eye(1 << bits, dtype=complex)
    expected[np.ix_([a, b], [a, b])] = [[math.cos(theta), math.sin(theta)],
                                       [-math.sin(theta), math.cos(theta)]]
    assert_matrix(c.execute(c._gray_rotation(a, b, theta, bits), bits), expected, 3e-14)
    # Removing the projector rotates unrelated fixed-bit configurations.
    label = ["I"]*bits
    label[(a ^ b).bit_length()-1] = "Y"
    unprojected = expm(1j*theta*pauli_matrix(label))
    assert np.linalg.norm(unprojected-expected, 2) > .3


@pytest.mark.parametrize("bits", (1, 2, 3, 6))
def test_real_preparation_executes_rotations_instead_of_replacing_the_input(bits):
    size = 1 << bits
    target = np.cos(np.arange(size)*.71)
    target[np.arange(size) % 3 == 1] = 0
    target /= np.linalg.norm(target)
    tape = c.compile_real_state(target)
    blank = np.eye(size)[:, 0]
    assert np.linalg.norm(c.execute(tape, bits, blank)-target) < 2e-10
    # A second orthogonal input is propagated, not overwritten with target.
    second = np.eye(size)[:, 1]
    pair = c.execute(tape, bits, np.stack((blank, second), axis=1))
    assert_matrix(pair.conj().T@pair, np.eye(2), 2e-10)
    assert np.linalg.norm(pair[:, 0]-pair[:, 1]) == pytest.approx(math.sqrt(2), abs=2e-10)


@pytest.mark.parametrize("target", (np.array([-1., 0.]), np.array([0., -1.]),
                                  np.array([0., 0., 1., 0.])))
def test_state_preparation_preserves_negative_global_sign_and_zero_subtrees(target):
    bits = int(math.log2(len(target)))
    assert np.linalg.norm(c.execute(c.compile_real_state(target), bits, np.eye(len(target))[:, 0])-target) < 1e-14


@pytest.mark.parametrize("d", (2, 4, 8))
def test_finite_scalar_factors_against_direct_exponentials_on_every_column(d):
    bits = int(math.log2(d))
    ladder = np.diag(np.sqrt(np.arange(1, d)), 1)
    position = (ladder+ladder.T)/math.sqrt(2)
    momentum = (ladder-ladder.T)/(1j*math.sqrt(2))
    values, vectors = np.linalg.eigh(position)
    basis = c.compile_real_orthogonal(vectors)
    number = c.compile_diagonal(np.arange(d)*math.pi/2)
    angle = -.19
    position_tape = c.inverse(basis)+c.compile_diagonal(angle*values)+basis
    quadratic_tape = c.inverse(basis)+c.compile_diagonal(angle*values**2)+basis
    momentum_tape = c.inverse(number)+quadratic_tape+number
    assert_matrix(c.execute(position_tape, bits), expm(1j*angle*position))
    assert_matrix(c.execute(quadratic_tape, bits), expm(1j*angle*position@position))
    assert_matrix(c.execute(momentum_tape, bits), expm(1j*angle*momentum@momentum))
    both = c.remap(basis, list(range(bits)))+c.remap(basis, list(range(bits, 2*bits)))
    pair_tape = c.inverse(both)+c.compile_diagonal(angle*np.outer(values, values).reshape(-1))+both
    assert_matrix(c.execute(pair_tape, 2*bits), expm(1j*angle*np.kron(position, position)))
    pointer_tape = c.inverse(basis)+c.compile_diagonal(angle*np.kron([1., -1.], values))+basis
    assert_matrix(c.execute(pointer_tape, bits+1), expm(1j*angle*np.kron(Z, position)))
    # Finite squared quadrature is not compression of the unbounded square.
    compressed = position@position+np.diag([0.]*(d-1)+[d/2])
    assert np.linalg.norm(expm(1j*angle*compressed)-c.execute(quadratic_tape, bits), 2) > .1
    if d == 2:
        unitary = c.execute(quadratic_tape, bits)
        assert_matrix(unitary.conj().T@momentum@unitary, momentum, 2e-14)
        assert np.linalg.norm(momentum-2*angle*position-unitary.conj().T@momentum@unitary, 2) > .2


def test_zero_requests_and_full_service_resource_counts_remain_in_tape():
    tape = [c.pulse("X", 0., [0]), c.pulse("J", 0., [0, 1]), c.pulse("I", -.4, [1])]
    report = c.resources(tape)
    assert report["native_pulses"] == 3
    assert report["zero_angle_requests"] == 2
    assert report["proper_code_service_events"] == 6
    assert report["proper_code_load_events"] == report["proper_code_export_events"] == 3
    assert report["valid_code_failure_probability_in_exact_arithmetic"] == 0
    assert report["maximum_numerical_valid_code_leakage"] < 2e-14
    assert report["service_failure_outcomes_retained"] is True
    assert_matrix(c.execute(tape, 2), np.exp(-.4j)*np.eye(4), 5e-15)


def test_structural_schedule_counts_do_not_depend_on_numerical_eigenbasis_zeros():
    size = 8
    dense, _ = np.linalg.qr(np.random.default_rng(471).normal(size=(size, size)))
    sparse = np.eye(size)
    assert len(c.compile_real_orthogonal(dense)) == len(c.compile_real_orthogonal(sparse))
    state = np.arange(1., size+1)
    state /= np.linalg.norm(state)
    assert len(c.compile_real_state(state)) == len(c.compile_real_state(sparse[:, 0]))


def test_corrupted_source_control_is_rejected_before_pulse_execution(tmp_path, monkeypatch):
    for path in (c.RESPONSE, c.CONTROL_RECEIPT):
        dest = tmp_path/path
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes((c.ROOT/path).read_bytes())
    path = tmp_path/c.CONTROL_RECEIPT
    data = json.loads(path.read_text())
    data["evidence"]["controls"]["rows"][0]["coefficients"][0][0] = "31/7"
    path.write_text(json.dumps(data))
    monkeypatch.setattr(c, "ROOT", tmp_path)
    c._source.cache_clear()
    try:
        with pytest.raises(ValueError, match="exact source control mismatch"):
            c.source_summary()
    finally:
        c._source.cache_clear()


@pytest.mark.parametrize("mutate", (
    lambda tape: tape[0].update(kind="arbitrary_matrix"),
    lambda tape: tape[0].update(angle=float("nan")),
    lambda tape: tape[0].update(angle=True),
    lambda tape: tape[0].update(targets=[True]),
    lambda tape: tape[0].update(targets=[3]),
    lambda tape: tape[0].update(matrix=[[1, 0], [0, 1]]),
))
def test_invalid_native_requests_are_rejected(mutate):
    tape = [c.pulse("X", .2, [0])]
    mutate(tape)
    with pytest.raises(ValueError):
        c.execute(tape, 2)


def test_reject_invalid_compiler_domains():
    for bad in ([0., 1., 2.], [0., float("inf")], [False, True], [0j, 1j]):
        with pytest.raises(ValueError):
            c.compile_diagonal(bad)
    with pytest.raises(ValueError):
        c.compile_real_orthogonal(np.ones((4, 4)))
    with pytest.raises(ValueError):
        c.compile_real_orthogonal(np.eye(16))
    with pytest.raises(ValueError):
        c.compile_real_state([1., 1.])
    with pytest.raises(ValueError):
        c.execute([], 2, np.eye(2))
    with pytest.raises(ValueError):
        c.transfer_kraus(3)
    with pytest.raises(ValueError):
        c.pulse("J", .2, [0, 0])


def test_execution_leaves_inputs_and_tapes_unmodified_and_serializable():
    tape = c.compile_diagonal([.1, -.2, .3, -.4])
    copied = deepcopy(tape)
    state = np.arange(8).reshape(4, 2).astype(complex)
    old = state.copy()
    c.execute(tape, 2, state)
    assert tape == copied
    assert np.array_equal(state, old)
    json.dumps(c.resources(tape), allow_nan=False)
    json.dumps(c.source_summary(), allow_nan=False)


def test_guarded_polar_correction_retains_full_phase_and_is_reference_linear():
    # Known polar factor, including a nontrivial scalar phase. Tiny radial
    # drift imitates accumulation without introducing a gate target argument.
    rng = np.random.default_rng(206)
    raw_basis = rng.normal(size=(64, 64))+1j*rng.normal(size=(64, 64))
    unitary, _ = np.linalg.qr(raw_basis)
    unitary *= np.exp(.731j)
    radial = np.linspace(-3e-12, 4e-12, 64)
    raw = unitary@np.diag(1+radial)
    original = raw.copy()
    corrected, report = c.stabilize_executed_unitary(raw)
    assert_matrix(corrected, unitary, 3e-14)
    assert report["projection_change"] == pytest.approx(max(abs(radial)), abs=5e-15)
    assert report["raw_unitarity_defect"] > 7e-12
    assert report["restored_unitarity_defect"] < 2e-14
    assert report["projection_change"] <= report["maximum_projection_change"]
    assert report["raw_unitarity_defect"] <= report["maximum_raw_unitarity_defect"]
    assert report["physical_error_certificate"] is False
    assert np.array_equal(raw, original)
    # A single matrix acts on both columns of an untouched-reference state.
    entangled = np.eye(64)[:, :3]/math.sqrt(3)
    assert_matrix(corrected@entangled, unitary@entangled, 3e-14)
    assert np.linalg.norm(corrected@entangled) == pytest.approx(1., abs=2e-14)
    # No fitted global phase: replacing U by exp(-i*.731)U is observable here.
    assert np.linalg.norm(corrected-np.exp(-.731j)*unitary, 2) > .7
    json.dumps(report, allow_nan=False)


def test_guarded_polar_is_optional_and_keeps_the_raw_source_executor(monkeypatch):
    # Inject tiny backend drift into a raw native pulse. The raw executor
    # must still return it; stabilization is a separate explicit operation.
    angle = .19
    actual_native = c.native_pulse("X", angle)
    drifted = actual_native.copy()
    drifted[:2, :2] *= .99999999999999
    monkeypatch.setattr(c, "_native", lambda kind, value: drifted)
    raw = c.execute([c.pulse("X", angle, [0])], 1)
    assert np.array_equal(raw, drifted[:2, :2])
    assert np.linalg.norm(raw.conj().T@raw-np.eye(2), 2) > 1e-14
    corrected, report = c.stabilize_executed_unitary(raw)
    assert_matrix(corrected, expm(1j*angle*X), 1e-14)
    assert report["projection_change"] > 5e-15


@pytest.mark.parametrize("bad", (.99*np.eye(4), 1.01*np.eye(4), np.zeros((4, 4)),
                                 np.diag([1., 1., 1., 1-1e-8])))
def test_guarded_polar_rejects_genuinely_nonunitary_input(bad):
    with pytest.raises(ArithmeticError, match="fixed rounding guard"):
        c.stabilize_executed_unitary(bad)


@pytest.mark.parametrize("bad", (np.ones(4), np.ones((4, 2)), np.eye(3),
                                 np.full((2, 2), np.nan), np.ones((2, 2), dtype=bool),
                                 np.ma.array(np.eye(2), mask=False)))
def test_guarded_polar_rejects_invalid_matrix_domains(bad):
    with pytest.raises(ValueError):
        c.stabilize_executed_unitary(bad)


def test_guarded_polar_detects_a_failed_projection_backend(monkeypatch):
    def corrupt_svd(matrix, full_matrices):
        dimension = len(matrix)
        return .99*np.eye(dimension), np.ones(dimension), np.eye(dimension)
    monkeypatch.setattr(c.np.linalg, "svd", corrupt_svd)
    # Use the Frobenius norm for this hostile backend test, because NumPy's
    # spectral-norm implementation itself calls the patched SVD.
    actual_norm = np.linalg.norm
    monkeypatch.setattr(c.np.linalg, "norm", lambda matrix, ord=None: actual_norm(matrix, "fro"))
    with pytest.raises(ArithmeticError, match="polar correction"):
        c.stabilize_executed_unitary(np.eye(4))
