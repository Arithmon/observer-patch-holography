"""Independent controls for observable-coordinate invariance (issue #1033)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from maxent_closure_acceptance import (
    duhamel_covariance, gibbs_state, independent_operator_count, i_projection,
    project_information, projection_diagnostics,
)

I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])


@pytest.mark.parametrize("scale", [1e-100, 1e-5, 1., 1e5, 1e100, -1e-5])
def test_diagonal_projection_does_not_depend_on_units(scale):
    # Analytic Bernoulli solution, without the Gibbs producer as the oracle.
    target = np.diag([.7, .3])
    assert independent_operator_count([scale*Z]) == 1
    lam, _ = i_projection(target, [scale*Z])
    assert scale*lam[0] == pytest.approx(-.5*np.log(7/3), abs=2e-11)
    rho, _ = gibbs_state([scale*Z], lam)
    np.testing.assert_allclose(rho, target, atol=2e-11, rtol=0)


def test_identity_shift_cannot_destroy_constraint_rank():
    ops = [Z + 1e8*I]
    assert independent_operator_count(ops) == 1
    lam, _ = i_projection(np.diag([.7, .3]), ops)
    assert lam[0] == pytest.approx(-.5*np.log(7/3), abs=2e-11)


def test_noncommuting_rechart_preserves_the_analytic_bloch_solution():
    target = (I + .3*X + .2*Y + .4*Z)/2
    ops = [Z, Z + 1e-5*X]
    assert independent_operator_count(ops) == 2
    lam, _ = i_projection(target, ops)
    rho, _ = gibbs_state(ops, lam)
    # MaxEnt removes only the unconstrained Y component.
    np.testing.assert_allclose(rho, (I + .3*X + .4*Z)/2, atol=2e-10, rtol=0)


def test_scalar_energy_origin_does_not_erase_resolved_offdiagonal_physics():
    # The literal spectrum of 1e20 I + X loses the splitting in float64.
    # That scalar has no effect on the exact normalized Gibbs state.
    rho, _ = gibbs_state([1e20*I + X], [.4])
    np.testing.assert_allclose(rho, (I - np.tanh(.4)*X)/2, atol=2e-14)
    cov = duhamel_covariance([1e20*I + X], [.4])
    np.testing.assert_allclose(cov, [[1/np.cosh(.4)**2]], atol=2e-13)


def test_small_raw_gradient_cannot_claim_convergence_in_wrong_units():
    # The old rank gate accepts this scale, then accepts lambda=0 immediately.
    # Raw error 4e-12 hides the physically unchanged Bloch error 4e-8.
    target = (I + 4e-8*Z)/2
    result = project_information(target, [1e-4*Z])
    assert result.iterations > 0
    assert result.normalized_residual < 1e-11
    np.testing.assert_allclose(result.state, target, atol=1e-12, rtol=0)


@pytest.mark.parametrize("scale", [1e-300, 1e-200, 1e200, 1e300])
def test_extreme_representable_units_do_not_square_into_a_rank_cutoff(scale):
    result = project_information((I + .4*Z)/2, [scale*Z])
    assert scale*result.multipliers[0] == pytest.approx(-np.arctanh(.4), abs=2e-11)
    assert result.normalized_residual < 1e-11


@pytest.mark.parametrize("seed", range(5))
def test_affine_and_unitary_recharts_preserve_state_and_covariant_multipliers(seed):
    rng = np.random.default_rng(seed)
    unitary, _ = np.linalg.qr(rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2)))
    transform = np.array([[2., -.7], [.3, 1.2]])
    units = np.array([1e-7, -1e5])
    transform *= units[:, None]
    shifts = np.array([3., -5.])*units
    base = [X, Z]
    changed = [unitary @ (sum(c*a for c, a in zip(row, base))+b*I) @ unitary.conj().T
               for row, b in zip(transform, shifts)]
    target = (I + .3*X + .2*Y + .4*Z)/2
    result = project_information(unitary @ target @ unitary.conj().T, changed)
    expected = unitary @ ((I+.3*X+.4*Z)/2) @ unitary.conj().T
    np.testing.assert_allclose(result.state, expected, atol=2e-11, rtol=0)
    np.testing.assert_allclose(transform.T @ result.multipliers,
                               -np.arctanh(.5)*np.array([.3, .4])/.5,
                               atol=3e-11, rtol=0)
    # Hessians transform covariantly, not by identical entry values.
    canonical_lam = transform.T @ result.multipliers
    expected_hessian = transform @ duhamel_covariance(base, canonical_lam) @ transform.T
    np.testing.assert_allclose(duhamel_covariance(changed, result.multipliers),
                               expected_hessian, rtol=2e-11, atol=1e-24)


def test_hessian_matches_noncommuting_frechet_exponential_control():
    from scipy.linalg import expm, expm_frechet
    ops = [X, Y, Z]
    lam = np.array([.31, -.27, .46])
    minus_h = -sum(c*a for c, a in zip(lam, ops))
    exponential = expm(minus_h)
    partition = np.trace(exponential).real
    oracle = np.empty((3, 3))
    for b, direction in enumerate(ops):
        derivative = expm_frechet(minus_h, -direction, compute_expm=False)
        derivative_rho = derivative/partition-exponential*np.trace(derivative)/partition**2
        for a, observable in enumerate(ops):
            oracle[a, b] = -np.trace(derivative_rho @ observable).real
    np.testing.assert_allclose(duhamel_covariance(ops, lam), oracle, atol=3e-14, rtol=0)


def test_logarithmic_mean_near_degeneracy_against_high_precision():
    import mpmath as mp
    for splitting in [0., 1e-14, 1e-9, .7, 350.]:
        with mp.workdps(80):
            expected = float(mp.tanh(splitting)/splitting) if splitting else 1.
        # The Hamiltonian is diagonal: its very small second eigenvalue is
        # explicit and does not require dense support inference.
        result = duhamel_covariance([Z, X], [splitting, 0.])
        assert result[1, 1] == pytest.approx(expected, rel=2e-13, abs=1e-16)


def test_qutrit_diagonal_projection_uses_all_independent_moments():
    # A full commuting algebra has an independent closed-form solution:
    # retain the input diagonal and discard its off-diagonal coherence.
    diagonal = np.array([.2, .3, .5])
    target = np.diag(diagonal).astype(complex)
    target[0, 1], target[1, 0] = .03j, -.03j
    ops = [np.diag([1., -1., 0.]), np.diag([0., 1., -1.])]
    result = project_information(target, ops)
    np.testing.assert_allclose(result.state, np.diag(diagonal), atol=1e-11, rtol=0)
    energies = -np.log(diagonal)
    energies -= np.mean(energies)
    np.testing.assert_allclose(sum(c*a for c, a in zip(result.multipliers, ops)),
                               np.diag(energies), atol=1e-10, rtol=0)


def test_noncommuting_qutrit_optimum_from_separate_matrix_exponential():
    from scipy.linalg import expm
    a = np.diag([1., 0., -1.])
    b = np.array([[0., 1., 0.], [1., 0., 1.], [0., 1., 0.]])
    optimum = expm(-.3*a+.2*b)
    optimum /= np.trace(optimum)
    unmeasured = np.array([[0., 0., 1j], [0., 0., 0.], [-1j, 0., 0.]])
    # The added coherence is trace-orthogonal to I,a,b, hence it does not
    # change either constrained moment. It is small enough to retain support.
    target = optimum+.02*unmeasured
    assert np.linalg.eigvalsh(target)[0] > .1
    result = project_information(target, [1e-6*a, a+1e-4*b])
    np.testing.assert_allclose(result.state, optimum, atol=2e-11, rtol=0)
    assert np.linalg.norm(target-result.state) > .02


def test_full_projection_ignores_scalar_energy_origin():
    result = project_information((I+.4*X)/2, [1e20*I+X])
    assert result.multipliers[0] == pytest.approx(-np.arctanh(.4), abs=1e-11)
    np.testing.assert_allclose(result.state, (I+.4*X)/2, atol=1e-11, rtol=0)


def test_numerical_rank_agrees_with_exact_rational_control():
    import sympy as sp
    a, b = sp.diag(1, -1, 0), sp.diag(0, 1, -1)
    ops = [sp.eye(3)*10**12 + a/sp.Integer(10)**20,
           b/sp.Integer(10)**30, (a+b)/sp.Integer(10)**30]
    # Use off-diagonal variation for the huge scalar: that survives float
    # storage. Tiny diagonal differences already rounded away are not data.
    ops[0] = sp.eye(3)*10**12 + sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    exact = sp.Matrix.hstack(*[m.reshape(9, 1) for m in [sp.eye(3), *ops]]).rank()-1
    assert independent_operator_count([np.array(m, complex) for m in ops]) == exact == 3


@pytest.mark.parametrize("ops", [[I], [X, 2*X], [Z, Z+1e-15*X], [X, Y, Z, X+Y]])
def test_redundant_or_unresolved_families_never_get_regularized(ops):
    with pytest.raises(ValueError, match="independent"):
        project_information(I/2, ops)


@pytest.mark.parametrize("ops", [[], [np.zeros((0, 0))], [np.ones((2, 3))],
                                [X, np.eye(3)], [np.array([[0., 1.], [0., 0.]])],
                                [X*1e-100+1e-100j*I], [[[True, 0], [0, 1]]],
                                [np.full((2, 2), np.inf)], [np.full((2, 2), np.nan)]])
def test_malformed_observables_rejected_before_solving(ops):
    for operation in [independent_operator_count, lambda a: project_information(I/2, a)]:
        with pytest.raises(ValueError):
            operation(ops)


@pytest.mark.parametrize("lam", [[True, .2], [1j, .2], [np.inf, .2], [.2], [[.2, .3]]])
def test_malformed_multipliers_rejected_by_both_thermal_consumers(lam):
    for operation in [gibbs_state, duhamel_covariance]:
        with pytest.raises(ValueError):
            operation([X, Z], lam)


def test_identity_offset_cannot_hide_nonhermiticity():
    malformed = 1e20*I + np.array([[0., 1.], [0., 0.]])
    with pytest.raises(ValueError, match="Hermitian"):
        gibbs_state([malformed], [.2])


def test_output_covariance_does_not_fabricate_a_zero_in_unrepresentable_units():
    with pytest.raises(ValueError, match="underflow"):
        duhamel_covariance([1e-200*Z], [0.])
    with pytest.raises(ValueError, match="numerical range"):
        duhamel_covariance([1e200*Z], [0.])


@pytest.mark.parametrize("ops, lam", [([1e-200*Z], [1e-200]),
                                    ([Z+1e-200*X], [1e-200])])
def test_hamiltonian_underflow_is_rejected_instead_of_erasing_terms(ops, lam):
    with pytest.raises(ValueError, match="underflow"):
        gibbs_state(ops, lam)


def test_fail_closed_on_iteration_budget_and_unusable_tolerance():
    for budget in [0, True, 1.5]:
        with pytest.raises(ValueError):
            project_information((I+.4*Z)/2, [Z], max_iter=budget)
    with pytest.raises(RuntimeError, match="converge"):
        project_information((I+.4*Z)/2, [Z], max_iter=1)
    for tol in [0., -1., 1., np.inf, True]:
        with pytest.raises(ValueError):
            project_information(I/2, [Z], tol=tol)


def test_diagnostics_reject_singular_targets_instead_of_extrapolating_faithful_theorem():
    with pytest.raises(ValueError, match="faithful"):
        projection_diagnostics(np.diag([1., 0.]), [Z], [1.])


@pytest.mark.parametrize("candidate", [[0., 0.], [.2, -.7], [-3., 1.], [.01, .02]])
def test_global_optimum_bound_for_nonstationary_candidates(candidate):
    from scipy.linalg import logm
    target = (I+.3*X+.2*Y+.4*Z)/2
    optimum = (I+.3*X+.4*Z)/2
    report = projection_diagnostics(target, [X, Z], candidate)
    gap = float(np.trace(optimum @ (logm(optimum)-logm(report.state))).real)
    distance = float(np.linalg.norm(np.linalg.eigvalsh(optimum-report.state), ord=1))
    assert 0 < gap <= report.optimality_gap_bound
    assert distance <= report.trace_distance_bound
    assert report.normalized_residual > 1e-3  # zero multipliers aren't a solution


def test_bound_and_pythagorean_accounting_separate_solver_error_from_closure_defect():
    from scipy.linalg import logm
    target = (I+.3*X+.2*Y+.4*Z)/2
    optimum = (I+.3*X+.4*Z)/2
    report = projection_diagnostics(target, [X, Z], [.2, -.7])
    entropy = lambda a, b: float(np.trace(a @ (logm(a)-logm(b))).real)
    defect = entropy(target, optimum)
    gap = entropy(optimum, report.state)
    assert defect > .01
    assert entropy(target, report.state) == pytest.approx(defect+gap, abs=2e-14)
    solved = project_information(target, [X, Z])
    assert solved.optimality_gap_bound < 1e-9
    assert entropy(target, solved.state) > .01  # no fabricated closure


def test_bound_and_residual_are_invariant_under_rescaling_and_mixing():
    target = (I+.3*X+.2*Y+.4*Z)/2
    transform = np.array([[2e-7, -3e-7], [1e5, 4e5]])
    base = [X, Z]
    changed = [sum(c*a for c, a in zip(row, base)) + offset*I
               for row, offset in zip(transform, [2e-7, 3e5])]
    lam = np.array([.2, -.7])
    reference = projection_diagnostics(target, base, lam)
    report = projection_diagnostics(target, changed, np.linalg.solve(transform.T, lam))
    assert report.normalized_residual == pytest.approx(reference.normalized_residual, rel=1e-12)
    assert report.optimality_gap_bound == pytest.approx(reference.optimality_gap_bound, rel=1e-12)


def test_large_raw_residual_is_informational_and_does_not_change_stopping_units():
    target = (I+.4*Z)/2
    result = project_information(target, [1e100*Z])
    _, legacy_raw = i_projection(target, [1e100*Z])
    assert legacy_raw == result.raw_moment_residual
    assert result.normalized_residual < 1e-11
    np.testing.assert_allclose(result.state, target, atol=1e-11, rtol=0)


def test_original_multiplier_replay_must_reject_an_incorrect_state(monkeypatch):
    import maxent.information_projection as module
    # Internal Newton solves successfully; sabotage only the original-chart
    # replay. An internal convergence flag must not bless this wrong state.
    monkeypatch.setattr(module, "gibbs_state", lambda *args: (I/2, np.log(2)))
    with pytest.raises(RuntimeError, match="replay"):
        module.project_information((I+.4*Z)/2, [Z])


def test_candidate_diagnostics_run_without_optimizer(monkeypatch):
    import maxent.information_projection as module
    def forbidden(*args, **kwargs):
        raise AssertionError("optimizer must be absent from candidate diagnostics")
    monkeypatch.setattr(module, "project_information", forbidden)
    report = module.projection_diagnostics((I+.4*Z)/2, [1e-5*Z], [-np.arctanh(.4)*1e5])
    assert report.normalized_residual < 1e-14
    assert report.optimality_gap_bound < 1e-13


def test_extended_precision_cannot_silently_create_infinite_or_zero_constraints():
    if np.finfo(np.longdouble).maxexp <= np.finfo(float).maxexp:
        pytest.skip("extended exponent range unavailable on this platform")
    for scale in [np.longdouble('1e400'), np.longdouble('1e-400')]:
        with pytest.raises(ValueError, match="range|conversion"):
            independent_operator_count([np.diag([scale, -scale])])


def test_validation_survives_optimized_python(tmp_path):
    import os
    import subprocess
    source = '''
import numpy as np
from maxent.information_projection import project_information, gibbs_state
bad = [lambda: project_information(np.eye(2)/2, [np.eye(2)]),
       lambda: gibbs_state([np.diag([1.,-1.]), np.eye(2)], [True, .2])]
for call in bad:
    try: call()
    except ValueError: pass
    else: raise RuntimeError('invalid input accepted with assertions disabled')
'''
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]))
    subprocess.run([sys.executable, "-O", "-c", source], cwd=tmp_path, env=env,
                   check=True, capture_output=True, text=True)
