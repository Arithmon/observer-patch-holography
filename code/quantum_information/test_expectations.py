"""Independent analytic channels, exact obstructions and hostile map controls."""

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm, logm

from quantum_information import partial_trace, relative_entropy
from quantum_information.algebras import FiniteAlgebra, matrix_units, tensor_factor_algebra
from quantum_information.expectations import (
    apply_map, choi_matrix, expectation_diagnostics, gns_projection,
    repair_generator, state_preserving_expectation,
)

I = np.eye(2)
X = np.array([[0.,1.],[1.,0.]])
Y = np.array([[0.,-1j],[1j,0.]])
Z = np.diag([1.,-1.])


def linear_matrix(action, d):
    # Test-only extraction from independently declared actions, with explicit
    # index units rather than the producer's Gram projection.
    columns = []
    for i in range(d):
        for j in range(d):
            unit = np.zeros((d,d),complex)
            unit[i,j] = 1
            columns.append(action(unit).reshape(-1))
    return np.column_stack(columns)


def independent_choi(action, d):
    blocks = []
    for i in range(d):
        row = []
        for j in range(d):
            unit = np.zeros((d,d),complex)
            unit[i,j] = 1
            row.append(action(unit))
        blocks.append(row)
    return np.block(blocks)


def nontracial_factor_case():
    left = np.array([[.6,.1j],[-.1j,.4]])
    right = np.array([[.7,.1+.05j],[.1-.05j,.3]])
    reference = np.kron(left,right)
    algebra = tensor_factor_algebra([2,2],[0])

    def expectation(a):
        tensor = a.reshape(2,2,2,2)
        reduced = np.zeros((2,2),complex)
        for i in range(2):
            for j in range(2):
                reduced[i,j] = np.trace(right@tensor[i,:,j,:])
        return np.kron(reduced,I)
    return algebra, reference, right, linear_matrix(expectation,4)


def test_nontracial_expectation_matches_independent_weighted_partial_trace():
    algebra,rho,right,expected = nontracial_factor_case()
    actual = state_preserving_expectation(algebra,rho)
    assert np.allclose(actual,expected,rtol=0,atol=1e-13)
    diagnostics = expectation_diagnostics(expected,algebra,rho)
    assert diagnostics["passed"]
    # The observable expectation is unital; its dual is trace preserving.
    # Requiring ordinary trace preservation of E would reject a valid repair.
    assert diagnostics["ordinary_trace_preservation"] > .5
    rng = np.random.default_rng(8)
    a = rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    sigma = a@a.conj().T+np.eye(4)
    sigma /= np.trace(sigma)
    repaired = apply_map(actual,sigma,dual=True)
    assert np.allclose(repaired,np.kron(partial_trace(sigma,[2,2],[0]),right))
    assert np.trace(repaired) == pytest.approx(1)
    assert np.linalg.eigvalsh(repaired)[0] > 0


def test_choi_vectorization_with_complex_nonsymmetric_unitary():
    u = expm(-.3j*(X+.4*Y+.2*Z))
    action = lambda a: u.conj().T@a@u
    s = linear_matrix(action,2)
    assert np.allclose(choi_matrix(s,2),independent_choi(action,2),atol=1e-14)
    assert np.linalg.eigvalsh(choi_matrix(s,2))[0] > -1e-14
    rho = np.array([[.7,.1j],[-.1j,.3]])
    assert np.allclose(apply_map(s,rho,dual=True),u@rho@u.conj().T)


def test_exact_weighted_projection_counterexample_is_not_a_quantum_channel():
    # Exact rational proof of the failing positive input, independent of any
    # eigenvalue tolerance or floating-point channel construction.
    i, x, z = sp.eye(2),sp.Matrix([[0,1],[1,0]]),sp.diag(1,-1)
    positive = sp.eye(4)+sp.kronecker_product(z,x)
    image = sp.eye(4)+sp.Rational(1,2)*sp.kronecker_product(z*x,i)
    assert set(positive.eigenvals()) == {sp.Integer(0),sp.Integer(2)}
    assert image != image.conjugate().T
    reference = (sp.eye(4)+sp.Rational(1,2)*sp.kronecker_product(x,x))/4
    assert set(reference.eigenvals()) == {sp.Rational(1,8),sp.Rational(3,8)}
    # Verify the defining GNS orthogonality equations exactly, so the claimed
    # image is proved to be the projection without trusting the constructor.
    for b in (i,x,sp.Matrix([[0,-sp.I],[sp.I,0]]),z):
        lifted = sp.kronecker_product(b,i)
        assert sp.trace(reference*lifted.conjugate().T*(image-positive)) == 0
    rho = (np.eye(4)+.5*np.kron(X,X))/4
    algebra = tensor_factor_algebra([2,2],[0])
    candidate = gns_projection(algebra,rho)
    report = expectation_diagnostics(candidate,algebra,rho)
    for name in ("unitality","idempotence","range","fixes_algebra",
                 "reference_preservation","gns_self_adjointness"):
        assert report[name] < 1e-13
    assert not report["passed"]
    assert report["choi_hermiticity"] > 1
    assert np.allclose(apply_map(candidate,np.array(positive).astype(complex)),
                       np.array(image).astype(complex))
    with pytest.raises(ValueError,match="modular"):
        state_preserving_expectation(algebra,rho)


def test_nontracial_reference_cannot_be_repaired_by_unweighted_projection():
    algebra,rho,_,_ = nontracial_factor_case()
    wrong = linear_matrix(algebra.project,4)
    report = expectation_diagnostics(wrong,algebra,rho)
    assert report["choi_min_eigenvalue"] > -1e-13
    assert report["idempotence"] < 1e-13
    assert report["reference_preservation"] > .1
    assert not report["passed"]


def test_central_blocks_and_coherences_are_checked_on_the_full_algebra():
    # B = (M_2 tensor I_2) direct_sum C. Off-sector coherences must be erased.
    basis = []
    for unit in matrix_units(2):
        block = np.zeros((5,5),complex)
        block[:4,:4] = np.kron(unit,I)
        basis.append(block)
    scalar = np.zeros((5,5)); scalar[-1,-1] = 1
    algebra = FiniteAlgebra(basis+[scalar])
    _,block,right,local = nontracial_factor_case()
    rho = np.zeros((5,5),complex); rho[:4,:4] = .7*block; rho[-1,-1] = .3
    s = state_preserving_expectation(algebra,rho)
    off = np.zeros((5,5),complex); off[0,4] = 1
    assert np.linalg.norm(apply_map(s,off)) < 1e-14
    bad = s.copy(); bad[:,4] += off.reshape(-1)  # retain E_04 in the output
    assert not expectation_diagnostics(bad,algebra,rho)["passed"]
    assert np.allclose(apply_map(s,np.eye(5)),np.eye(5))
    assert np.allclose(apply_map(s,rho,dual=True),rho)


def test_pythagorean_identity_and_a3_minimizer_are_the_same_repaired_state():
    algebra,rho,_,s = nontracial_factor_case()
    rng = np.random.default_rng(17)
    x = rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    sigma = x@x.conj().T + np.eye(4)
    sigma /= np.trace(sigma)
    repaired = apply_map(s,sigma,dual=True)
    # Use SciPy logarithms as a separate numerical entropy implementation.
    d = lambda a,b: float(np.trace(a@(logm(a)-logm(b))).real)
    assert d(sigma,rho) == pytest.approx(d(sigma,repaired)+d(repaired,rho),abs=2e-13)
    assert d(sigma,repaired) > .1
    # Different full states with the same B data have the identical minimizer.
    # X tensor Z has zero B marginal and a small coefficient preserves PSD.
    alternative = repaired + .005*np.kron(X,Z)
    assert np.linalg.eigvalsh(alternative)[0] > 0
    assert np.allclose(apply_map(s,alternative,dual=True),repaired)
    assert d(alternative,rho) > d(repaired,rho)
    assert d(alternative,rho)-d(repaired,rho) == pytest.approx(d(alternative,repaired),abs=2e-13)
    assert relative_entropy(sigma,rho) == pytest.approx(d(sigma,rho),abs=2e-13)


def test_map_is_invariant_under_basis_change_and_unitary_recharting():
    algebra,rho,_,expected = nontracial_factor_case()
    basis = algebra.basis
    changed = FiniteAlgebra([basis[0],basis[1]+2j*basis[0],3*basis[2],basis[3]-basis[2]])
    assert np.allclose(state_preserving_expectation(changed,rho),expected,atol=1e-13)
    u = expm(.3j*np.kron(X,Y))
    recharted = FiniteAlgebra([u@b@u.conj().T for b in basis])
    s = state_preserving_expectation(recharted,u@rho@u.conj().T)
    for x in matrix_units(4):
        assert np.allclose(apply_map(s,u@x@u.conj().T),u@apply_map(expected,x)@u.conj().T,atol=2e-13)


def test_noncommuting_repairs_have_exact_intersection_and_gap():
    # Independent rational spectrum, rather than rounding a fitted eigenvalue.
    i,z,x = sp.eye(2),sp.diag(1,-1),sp.Matrix([[0,1],[1,0]])
    n = sp.Rational(3,5)*x+sp.Rational(4,5)*z
    # On row-major vectors, conjugation by real symmetric u is u tensor u.
    ez = (sp.eye(4)+sp.kronecker_product(z,z))/2
    en = (sp.eye(4)+sp.kronecker_product(n,n))/2
    exact = (ez+en)/2-sp.eye(4)
    assert ez*en != en*ez
    assert exact.eigenvals() == {sp.Integer(0):1,-sp.Rational(1,10):1,
                                -sp.Rational(9,10):1,sp.Integer(-1):1}
    a1,a2 = FiniteAlgebra([I,Z]),FiniteAlgebra([I,.6*X+.8*Z])
    s1,s2 = np.array(ez).astype(complex),np.array(en).astype(complex)
    generator,report = repair_generator([s1,s2],[a1,a2],I/2,[.5,.5])
    assert report["intersection_dimension"] == 1
    assert report["gap"] == pytest.approx(.1,abs=1e-14)
    assert np.allclose(generator,np.array(exact).astype(complex),atol=1e-14)
    for t in (0.,.2,2.,10.):
        channel = expm(t*generator).conj().T
        choi = choi_matrix(channel,2)
        assert np.linalg.eigvalsh(choi)[0] > -1e-13
        state = (I+.5*X+.2*Y+.3*Z)/2
        output = apply_map(channel,state)
        assert np.trace(output) == pytest.approx(1)
        assert relative_entropy(output,I/2) <= relative_entropy(state,I/2)+1e-13
        assert np.linalg.norm(output-I/2) <= np.exp(-.1*t)*np.linalg.norm(state-I/2)+1e-13


@pytest.mark.parametrize("clock", [1e-20, 1., 1e20])
def test_almost_parallel_repairs_cannot_invent_a_large_gap(clock):
    a1,a2 = FiniteAlgebra([I,Z]),FiniteAlgebra([I,Z+1e-8*X])
    s1,s2 = (state_preserving_expectation(a,I/2) for a in (a1,a2))
    with pytest.raises(ValueError,match="gap.*unresolved"):
        repair_generator([s1,s2],[a1,a2],I/2,[.5*clock,.5*clock])


@pytest.mark.parametrize("clock", [1e-250, 1e-20, 1., 1e20, 1e250])
def test_repair_clock_units_preserve_noncommuting_spectrum_and_evolution(clock):
    # Independent Bloch calculation: two equally weighted dephasings about
    # Z and (3X+4Z)/5 have decay rates (1/10, 9/10, 1), times clock.
    n = .6*X + .8*Z
    algebras = [FiniteAlgebra([I, axis]) for axis in (Z, n)]
    maps = [linear_matrix(lambda a, axis=axis: (a+axis@a@axis)/2, 2)
            for axis in (Z, n)]
    generator, report = repair_generator(maps, algebras, I/2, [.5*clock]*2)
    assert report["intersection_dimension"] == 1
    assert report["gap"]/clock == pytest.approx(.1, abs=1e-14)
    assert np.allclose(np.linalg.eigvalsh(-generator/clock),
                       [0., .1, .9, 1.], rtol=0, atol=1e-14)
    bloch_generator = (np.outer([0., 0., 1.], [0., 0., 1.])
                       + np.outer([.6, 0., .8], [.6, 0., .8]))/2 - np.eye(3)
    bloch = np.array([.3, .2, .4])
    state = (I + sum(v*a for v, a in zip(bloch, (X, Y, Z))))/2
    evolved_bloch = expm(.7*bloch_generator) @ bloch
    expected = (I + sum(v*a for v, a in zip(evolved_bloch, (X, Y, Z))))/2
    actual = apply_map(expm((.7/clock)*generator), state, dual=True)
    assert np.allclose(actual, expected, rtol=0, atol=1e-14)


@pytest.mark.parametrize("clock", [1e-20, 1e20])
def test_clock_units_preserve_nontracial_recharted_repair(clock):
    algebra, rho, _, channel = nontracial_factor_case()
    u = expm(.27j*(np.kron(X,Y)+.3*np.kron(Y,Z)))
    chart = linear_matrix(lambda a: u@a@u.conj().T, 4)
    changed_algebra = FiniteAlgebra([u@a@u.conj().T for a in algebra.basis])
    generator, report = repair_generator(
        [chart@channel@chart.conj().T], [changed_algebra],
        u@rho@u.conj().T, [clock])
    assert report["gap"]/clock == pytest.approx(1., abs=1e-13)
    assert np.allclose(generator/clock, chart@channel@chart.conj().T-np.eye(16),
                       rtol=0, atol=1e-13)


def test_repair_rate_range_cannot_erase_a_primitive_or_overflow_generator():
    algebras = [FiniteAlgebra([I, axis]) for axis in (Z, X)]
    maps = [linear_matrix(lambda a, axis=axis: (a+axis@a@axis)/2, 2)
            for axis in (Z, X)]
    for rates in ([1e-308, 1e308], [1.7e308, 1.7e308]):
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            with pytest.raises(ValueError):
                repair_generator(maps, algebras, I/2, rates)


def test_subnormal_clock_cannot_change_a_resolved_slow_decay_rate():
    # The normalized generator has a resolved gap, but rounding its physical
    # entries changes the decay of Z from 1.5e-323 to 2e-323. That is not an
    # admissible representation just because its error is below user tol.
    algebras = [FiniteAlgebra([I, axis]) for axis in (Z, X)]
    maps = [np.diag([1., 0., 0., 1.]), (np.eye(4)+np.kron(X, X))/2]
    with pytest.raises(ValueError, match="underflow"):
        repair_generator(maps, algebras, I/2, [1e-313, 1.5e-323])


def test_exact_single_dephasing_supports_a_representable_subnormal_clock():
    clock = np.nextafter(0., 1.)
    channel = np.diag([1., 0., 0., 1.])
    generator, report = repair_generator(
        [channel], [FiniteAlgebra([I, Z])], I/2, [clock])
    assert np.array_equal(generator.real, np.diag([0., -clock, -clock, 0.]))
    assert report["gap"] == clock


@pytest.mark.parametrize("bad", ([],[0.],[1.,0.],[-1.],[True],[1j],[np.nan],[np.inf],
                                 np.array(1.),np.array([[1.]])))
def test_repair_generator_rejects_empty_mismatched_or_unfair_rates(bad):
    algebra = FiniteAlgebra([I,Z])
    s = state_preserving_expectation(algebra,I/2)
    with pytest.raises(ValueError):
        repair_generator([s],[algebra],I/2,bad)


@pytest.mark.parametrize("tol", (0,-1,True,np.inf,np.nan,1.,1j,[1e-9]))
def test_invalid_tolerances_cannot_turn_off_map_verification(tol):
    algebra = FiniteAlgebra([I,Z])
    with pytest.raises(ValueError):
        expectation_diagnostics(np.zeros((4,4)),algebra,I/2,tol=tol)


@pytest.mark.parametrize("mutation", ("zero","scale","transpose","partial","wrong_algebra","nonfinite"))
def test_supplied_maps_are_recognized_without_trusting_their_claimed_properties(mutation):
    algebra = FiniteAlgebra([I,Z])
    s = state_preserving_expectation(algebra,I/2)
    if mutation == "zero": s *= 0
    elif mutation == "scale": s *= .9
    elif mutation == "transpose": s = linear_matrix(lambda a:a.T,2)
    elif mutation == "partial": s = .5*np.eye(4)+.5*s
    elif mutation == "wrong_algebra": s = state_preserving_expectation(FiniteAlgebra([I,X]),I/2)
    else: s[0,0] = np.nan
    if mutation == "nonfinite":
        with pytest.raises(ValueError): expectation_diagnostics(s,algebra,I/2)
    else:
        assert not expectation_diagnostics(s,algebra,I/2)["passed"]
        with pytest.raises(ValueError,match="supplied primitive"):
            repair_generator([s],[algebra],I/2,[1.])


def test_singular_references_are_outside_this_expectation_constructor():
    with pytest.raises(ValueError,match="faithful"):
        state_preserving_expectation(FiniteAlgebra([I,Z]),np.diag([1.,0.]))


def test_identity_repair_does_not_claim_a_nonexistent_positive_gap():
    full = FiniteAlgebra(matrix_units(2))
    generator,report = repair_generator([np.eye(4)],[full],I/2,[1.])
    assert np.linalg.norm(generator) == 0
    assert report["complement_dimension"] == 0
    assert report["gap"] is None


def test_existing_collar_alignment_criterion_matches_repair_existence():
    from collar_alignment.msa_characterizations import (
        modular_splitting_defect, random_aligned_blocks, random_generic_blocks, takesaki_defect,
    )
    rng = np.random.default_rng(35)
    algebra = tensor_factor_algebra([2,2],[0])
    spec = [(1,2,2,1)]
    for constructor,aligned in ((random_aligned_blocks,True),(random_generic_blocks,False)):
        blocks = constructor(rng,spec)
        rho = blocks[0][1]
        if aligned:
            assert modular_splitting_defect(blocks) < 1e-12
            assert takesaki_defect(blocks) < 1e-12
            s = state_preserving_expectation(algebra,rho)
            assert expectation_diagnostics(s,algebra,rho)["passed"]
        else:
            assert modular_splitting_defect(blocks) > .01
            assert takesaki_defect(blocks) > .01
            with pytest.raises(ValueError,match="modular"):
                state_preserving_expectation(algebra,rho)


def test_classical_fiber_resampling_is_the_commutative_specialization():
    # Independently specified transition table in the row-stochastic convention
    # of code/consensus/README.md. Fibers {0,1} and {2} have weights 1/6,1/3,1/2.
    kernel = np.array([[1/3,2/3,0],[1/3,2/3,0],[0,0,1]])
    reference = np.diag([1/6,1/3,1/2])
    algebra = FiniteAlgebra([np.diag([1.,1.,0.]),np.diag([0.,0.,1.])])
    s = linear_matrix(lambda a:np.diag(kernel@np.diag(a)),3)
    assert expectation_diagnostics(s,algebra,reference)["passed"]
    assert np.allclose(s,state_preserving_expectation(algebra,reference))
    assert np.allclose(np.diag(reference)@kernel,np.diag(reference))
    assert np.allclose(kernel@kernel,kernel)
    bad_kernel = np.array([[.5,.5,0],[.5,.5,0],[0,0,1]])
    wrong = linear_matrix(lambda a:np.diag(bad_kernel@np.diag(a)),3)
    assert not expectation_diagnostics(wrong,algebra,reference)["passed"]


def test_a3_completion_is_provably_nonlinear_for_incompatible_reference():
    # Diagonal data (1,0) and (0,1) force their unique positive completions.
    # Their mean data equal the diagonal of rho, whose unique D( . || rho)
    # minimizer is rho itself. No optimizer or fitted ansatz is used here.
    rho = np.array([[.5,.25],[.25,.5]])
    endpoint0,endpoint1 = np.diag([1.,0.]),np.diag([0.,1.])
    affine_candidate = (endpoint0+endpoint1)/2
    assert np.array_equal(np.diag(rho),np.diag(affine_candidate))
    assert relative_entropy(rho,rho) == 0
    assert relative_entropy(affine_candidate,rho) == pytest.approx(.5*np.log(4/3))
    assert np.linalg.norm(rho-affine_candidate) > .3
    with pytest.raises(ValueError,match="modular"):
        state_preserving_expectation(FiniteAlgebra([I,Z]),rho)


def test_affine_entropy_completion_on_an_operator_system_need_not_be_cp():
    # For specified X and Z expectations and tracial reference, entropy is
    # maximized by setting the unconstrained Y coefficient to zero. This is
    # (id+transpose)/2, positive and affine but NOT completely positive.
    # Its range span{I,X,Z} is an operator system, not an operator algebra.
    action = lambda a:(a+a.T)/2
    s = linear_matrix(action,2)
    choi = independent_choi(action,2)
    # J_id = |vec(I)><vec(I)|, J_transpose = swap, both unnormalized.
    bell = sp.Matrix([1,0,0,1])
    swap = sp.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    exact_choi = (bell*bell.T+swap)/2
    assert exact_choi.eigenvals() == {
        -sp.Rational(1,2):1,sp.Rational(1,2):2,sp.Rational(3,2):1,
    }
    assert np.array_equal(choi,np.array(exact_choi).astype(complex))
    assert np.allclose(np.linalg.eigvalsh(choi),[-.5,.5,.5,1.5])
    assert np.allclose(s@s,s)
    assert np.allclose(apply_map(s,I),I)
    for y in (-.4,0.,.4):
        sigma = (I+.3*X+y*Y+.2*Z)/2
        reduced = apply_map(s,sigma)
        assert np.linalg.eigvalsh(reduced)[0] > 0
        assert relative_entropy(reduced,I/2) <= relative_entropy(sigma,I/2)+1e-14
    with pytest.raises(ValueError,match="multiplication"):
        FiniteAlgebra([I,X,Z])


def test_entropy_identity_includes_singular_marginals_and_central_coherence():
    basis = []
    for unit in matrix_units(2):
        block = np.zeros((5,5),complex)
        block[:4,:4] = np.kron(unit,I)
        basis.append(block)
    scalar = np.diag([0.,0.,0.,0.,1.])
    algebra = FiniteAlgebra(basis+[scalar])
    rho = np.diag([.1,.1,.2,.2,.4])
    s = state_preserving_expectation(algebra,rho)
    psi = np.array([.6,0.,0.,0.,.8])
    sigma = np.outer(psi,psi)
    gamma = np.diag([.18,.18,0.,0.,.64])
    assert np.allclose(apply_map(s,sigma,dual=True),gamma,atol=1e-14)
    d_full = -.36*np.log(.1)-.64*np.log(.4)
    d_removed = -.36*np.log(.18)-.64*np.log(.64)
    d_retained = .36*np.log(1.8)+.64*np.log(1.6)
    assert d_full == pytest.approx(d_removed+d_retained,abs=1e-14)
    assert relative_entropy(sigma,rho) == pytest.approx(d_full,abs=1e-14)
    assert relative_entropy(sigma,gamma) == pytest.approx(d_removed,abs=1e-14)
    assert relative_entropy(gamma,rho) == pytest.approx(d_retained,abs=1e-14)


def test_huge_finite_maps_cannot_pass_after_diagnostic_overflow():
    algebra = FiniteAlgebra([I,Z])
    with np.errstate(over="ignore",invalid="ignore"):
        with pytest.raises(ValueError,match="finite numerical range"):
            expectation_diagnostics(np.eye(4)*1e200,algebra,I/2)


def test_boundary_entropy_identity_does_not_commute_support_with_reference():
    algebra,rho,right,s = nontracial_factor_case()
    pure_left = np.diag([1.,0.])
    input_right = np.diag([.4,.6])
    sigma = np.kron(pure_left,input_right)
    gamma = np.kron(pure_left,right)
    support = np.kron(pure_left,I)
    assert np.linalg.norm(support@rho-rho@support) > .05
    assert np.allclose(apply_map(s,sigma,dual=True),gamma,atol=1e-14)
    left = np.array([[.6,.1j],[-.1j,.4]])
    retained = -float(logm(left)[0,0].real)
    removed = float(np.trace(input_right@(logm(input_right)-logm(right))).real)
    # The support compression uses Q log(rho) Q, not log(Q rho Q).
    # Their difference is nonzero even though the entropy identity holds.
    assert abs(retained+np.log(left[0,0].real)) > .01
    assert relative_entropy(gamma,rho) == pytest.approx(retained,abs=1e-13)
    assert relative_entropy(sigma,gamma) == pytest.approx(removed,abs=1e-13)
    assert relative_entropy(sigma,rho) == pytest.approx(retained+removed,abs=1e-13)


def test_nontracial_repair_generator_after_complex_unitary_recharting():
    # Two noncommuting dephasings on R, with the entire left algebra retained.
    # A complex global rechart makes the GNS weight nontrivial and tests the
    # weighted similarity transform, which a tracial example cannot exercise.
    left = np.array([[.6,.1j],[-.1j,.4]])
    rho = np.kron(left,I/2)
    n = .6*X+.8*Z
    primitives = []
    algebras = []
    for axis in (Z,n):
        u = np.kron(I,axis)
        primitives.append(linear_matrix(lambda a,u=u:(a+u@a@u)/2,4))
        algebras.append(FiniteAlgebra(
            [np.kron(a,b) for a in matrix_units(2) for b in (I,axis)]))
    generator,report = repair_generator(primitives,algebras,rho,[.5,.5])
    expected = np.repeat([0.,.1,.9,1.],4)
    assert np.allclose(np.sort(np.linalg.eigvals(-generator).real),expected,atol=1e-13)
    assert report["intersection_dimension"] == 4
    assert report["gap"] == pytest.approx(.1,abs=1e-13)
    u = expm(.27j*(np.kron(X,Y)+.3*np.kron(Y,Z)))
    rechart = linear_matrix(lambda a:u@a@u.conj().T,4)
    changed = [rechart@s@rechart.conj().T for s in primitives]
    changed_algebras = [FiniteAlgebra([u@b@u.conj().T for b in a.basis])
                       for a in algebras]
    result,changed_report = repair_generator(
        changed,changed_algebras,u@rho@u.conj().T,[.5,.5])
    assert np.allclose(result,rechart@generator@rechart.conj().T,atol=2e-13)
    assert changed_report["intersection_dimension"] == 4
    assert changed_report["gap"] == pytest.approx(.1,abs=1e-13)
    t = .7
    channel = expm(t*result)
    assert np.linalg.eigvalsh(choi_matrix(channel,4))[0] > -1e-13
    assert np.allclose(apply_map(channel,u@rho@u.conj().T,dual=True),
                       u@rho@u.conj().T,atol=1e-13)


def test_validation_survives_optimized_python():
    import os
    from pathlib import Path
    import subprocess
    import sys

    program = '''
import numpy as np
from quantum_information.algebras import FiniteAlgebra, separation_modulus
from quantum_information.expectations import expectation_diagnostics, repair_generator
i = np.eye(2)
z = np.diag([1.,-1.])
a = FiniteAlgebra([i,z])
checks = [lambda: FiniteAlgebra([i,z,np.array([[0.,1.],[1.,0.]])]),
          lambda: separation_modulus([i],np.array([2**32,1],dtype=np.int64)),
          lambda: expectation_diagnostics(np.full((4,4),np.nan),a,i/2),
          lambda: expectation_diagnostics(np.zeros((4,4)),a,i/2,tol=np.inf),
          lambda: repair_generator([np.zeros((4,4))],[a],i/2,[1.])]
for check in checks:
    try:
        check()
    except ValueError:
        continue
    raise SystemExit("invalid input accepted under python -O")
if expectation_diagnostics(np.zeros((4,4)),a,i/2)["passed"]:
    raise SystemExit("zero map accepted under python -O")
'''
    env = dict(os.environ,PYTHONPATH=str(Path(__file__).resolve().parents[1]))
    result = subprocess.run([sys.executable,"-O","-c",program],env=env,
                            capture_output=True,text=True,timeout=30)
    assert result.returncode == 0, result.stdout+result.stderr
