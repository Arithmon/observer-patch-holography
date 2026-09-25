import ObserverPatchHolography.EinsteinBranch.EdgeCenterTiltCocycle

/-!
Kernel-checked algebraic cores. The finite-state stochastic limit theorem is
proved separately in REPORT.md; it is NOT formalized by these matrix identities.
The tilt witness instantiates the existing canonical SurvivalCocycle structure.
No measured scalar index, spectrum, or cosmological identification is assumed here.
-/
namespace Codex.Necessity
open Matrix

noncomputable def exponentialCocycle (theta : ℝ) :
    OPH.EinsteinBranch.EdgeCenterTilt.SurvivalCocycle where
  u := fun t => Real.exp (-(theta * t))
  pos := fun _ => Real.exp_pos _
  mul := by
    intro s t
    rw [mul_add, neg_add, Real.exp_add]
  cont := by fun_prop

theorem arbitrary_tilt (theta : ℝ) :
    (exponentialCocycle theta).tilt = theta := by
  apply OPH.EinsteinBranch.EdgeCenterTilt.SurvivalCocycle.tilt_eq_of_exp
  intro s
  rfl

theorem cocycle_does_not_select_tilt (theta : ℝ) :
    ∃ S : OPH.EinsteinBranch.EdgeCenterTilt.SurvivalCocycle,
      S.tilt ≠ theta := by
  refine ⟨exponentialCocycle (theta + 1), ?_⟩
  rw [arbitrary_tilt]
  linarith

theorem nonnegative_tilt_survival (theta t : ℝ) (htheta : 0 ≤ theta) (ht : 0 ≤ t) :
    0 < (exponentialCocycle theta).u t ∧ (exponentialCocycle theta).u t ≤ 1 := by
  constructor
  · exact Real.exp_pos _
  · change Real.exp (-(theta * t)) ≤ 1
    apply Real.exp_le_one_iff.mpr
    exact neg_nonpos.mpr (mul_nonneg htheta ht)

-- Six configurations of two raised loads on the four-cycle, ordered
-- {01,02,03,12,13,23}; each undirected edge attempts at rate one, swaps half.
def Q : Matrix (Fin 6) (Fin 6) ℚ :=
  !![-1, 1/2, 0, 0, 1/2, 0;
    1/2, -2, 1/2, 1/2, 0, 1/2;
    0, 1/2, -1, 0, 1/2, 0;
    0, 1/2, 0, -1, 1/2, 0;
    1/2, 0, 1/2, 1/2, -2, 1/2;
    0, 1/2, 0, 0, 1/2, -1]

def L : Matrix (Fin 4) (Fin 4) ℚ :=
  !![2, -1, 0, -1;
    -1, 2, -1, 0;
    0, -1, 2, -1;
    -1, 0, -1, 2]

def Z : Matrix (Fin 6) (Fin 4) ℚ :=
  !![1/2, 1/2, -1/2, -1/2;
    1/2, -1/2, 1/2, -1/2;
    1/2, -1/2, -1/2, 1/2;
    -1/2, 1/2, 1/2, -1/2;
    -1/2, 1/2, -1/2, 1/2;
    -1/2, -1/2, 1/2, 1/2]

def Pi : Matrix (Fin 4) (Fin 4) ℚ :=
  !![3/4, -1/4, -1/4, -1/4;
    -1/4, 3/4, -1/4, -1/4;
    -1/4, -1/4, 3/4, -1/4;
    -1/4, -1/4, -1/4, 3/4]

def K : Matrix (Fin 4) (Fin 4) ℚ :=
  !![5/16, -1/16, -3/16, -1/16;
    -1/16, 5/16, -1/16, -3/16;
    -3/16, -1/16, 5/16, -1/16;
    -1/16, -3/16, -1/16, 5/16]

def H : Matrix (Fin 6) (Fin 4) ℚ :=
  !![1/2, 1/2, -1/2, -1/2;
    1/4, -1/4, 1/4, -1/4;
    1/2, -1/2, -1/2, 1/2;
    -1/2, 1/2, 1/2, -1/2;
    -1/4, 1/4, -1/4, 1/4;
    -1/2, -1/2, 1/2, 1/2]

-- These proofs expand finite rational sums and are checked by the kernel.
-- K is the mean-zero inverse of L; Pi is the constant-mode complement.
theorem q_symmetric : Q.transpose = Q := by
  ext i j; fin_cases i <;> fin_cases j <;> norm_num [Q]

theorem q_rows_zero : ∀ i, ∑ j, Q i j = 0 := by
  intro i; fin_cases i <;> norm_num [Q, Fin.sum_univ_succ]

theorem q_offdiag_nonnegative : ∀ i j, i ≠ j → 0 ≤ Q i j := by
  intro i j; fin_cases i <;> fin_cases j <;> norm_num [Q]

theorem stationary_uniform : ∀ j, ∑ i, (1/6 : ℚ) * Q i j = 0 := by
  intro j; fin_cases j <;> norm_num [Q, Fin.sum_univ_succ]

theorem coordinate_drift : Q * Z = (-1/2 : ℚ) • (Z * L) := by
  ext i j; fin_cases i <;> fin_cases j <;>
    norm_num [Q, Z, L, Matrix.mul_apply, Fin.sum_univ_succ]

theorem static_covariance : (1/6 : ℚ) • (Z.transpose * Z) = (1/3 : ℚ) • Pi := by
  ext i j; fin_cases i <;> fin_cases j <;>
    norm_num [Z, Pi, Matrix.mul_apply, Fin.sum_univ_succ]

theorem inverse_on_meanzero : L * K = Pi ∧ K * L = Pi := by
  constructor <;> ext i j <;> fin_cases i <;> fin_cases j <;>
    norm_num [L, K, Pi, Matrix.mul_apply, Fin.sum_univ_succ]

theorem poisson_solution : Q * H = -Z := by
  ext i j; fin_cases i <;> fin_cases j <;>
    norm_num [Q, H, Z, Matrix.mul_apply, Fin.sum_univ_succ]

-- The Green-Kubo algebra; interpreting it as a time-integral limiting
-- covariance additionally uses the proved finite-state probability theorem.
def greenCov : Matrix (Fin 4) (Fin 4) ℚ :=
  (1/6 : ℚ) • (Z.transpose * H + H.transpose * Z)

theorem green_covariance : greenCov = (4/3 : ℚ) • K := by
  ext i j; fin_cases i <;> fin_cases j <;>
    norm_num [greenCov, Z, H, K, Matrix.mul_apply, Fin.sum_univ_succ]

theorem local_difference_covariance : L * greenCov * L.transpose = (4/3 : ℚ) • L := by
  rw [green_covariance]
  ext i j; fin_cases i <;> fin_cases j <;>
    norm_num [L, K, Matrix.mul_apply, Fin.sum_univ_succ]

-- Distinct spectral shapes, not a mere gain change: compare a slow mode
-- and the alternating mode. Unit lengths are unnecessary because both
-- readouts use the same nonzero test vectors.
def slow : Fin 4 → ℚ := ![1,0,-1,0]
def fast : Fin 4 → ℚ := ![1,-1,1,-1]
def power (C : Matrix (Fin 4) (Fin 4) ℚ) (v : Fin 4 → ℚ) : ℚ :=
  dotProduct v (C.mulVec v)

theorem record_spectral_shapes_differ :
    power greenCov fast / power greenCov slow = 1 ∧
    power (L * greenCov * L.transpose) fast /
      power (L * greenCov * L.transpose) slow = 4 := by
  rw [local_difference_covariance, green_covariance]
  norm_num [power, dotProduct, Matrix.mulVec, K, L, slow, fast, Fin.sum_univ_succ]

-- Pointwise core of the common temporal-filter no-redder-than-Green bound.
-- Integrating this inequality against |Fourier(w)|² is the analytic step
-- proved in temporal_filters/REPORT.md, not formalized here.
theorem filter_integrand_monotone (a b omega : ℝ)
    (ha : 0 < a) (hab : a ≤ b) :
    a^2 / (a^2 + omega^2) ≤ b^2 / (b^2 + omega^2) := by
  have hb : 0 < b := lt_of_lt_of_le ha hab
  have hdena : 0 < a^2 + omega^2 := by positivity
  have hdenb : 0 < b^2 + omega^2 := by positivity
  apply (div_le_div_iff₀ hdena hdenb).mpr
  have hsq : a^2 ≤ b^2 := by nlinarith
  have hm := mul_nonneg (sub_nonneg.mpr hsq) (sq_nonneg omega)
  nlinarith

#print axioms arbitrary_tilt
#print axioms cocycle_does_not_select_tilt
#print axioms nonnegative_tilt_survival
#print axioms q_symmetric
#print axioms q_rows_zero
#print axioms q_offdiag_nonnegative
#print axioms stationary_uniform
#print axioms coordinate_drift
#print axioms static_covariance
#print axioms inverse_on_meanzero
#print axioms poisson_solution
#print axioms green_covariance
#print axioms local_difference_covariance
#print axioms record_spectral_shapes_differ
#print axioms filter_integrand_monotone
end Codex.Necessity
