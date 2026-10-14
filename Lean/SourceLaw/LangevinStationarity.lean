import SourceLaw.HistoryRegister

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Stationary law of the declared Langevin law

The declared law is the linear Langevin equation `df = −B f ds + √(2ℏ𝕄⁻¹) dW` with drift
`B = 𝕄⁻¹ P`, where `P` is the register stiffness and `𝕄` the diagonal capacity. A centred
Gaussian with covariance `Σ` is stationary for it exactly when the Lyapunov equation
`B Σ + Σ Bᵀ = 2ℏ𝕄⁻¹` holds (standard; recorded in the companion note). This file proves the
matrix algebra:

The algebra takes a symmetric `P` with invertible determinant; uniqueness takes positive
definite `P` and positive capacities.

* `lyapunov_stationary`: `Σ = ℏ P⁻¹` solves the Lyapunov equation; for positive definite `P` and
  `ℏ > 0` this makes the Gaussian `exp(−S/ℏ)` with `S = ½ fᵀ P f` a stationary law.
* `detailed_balance`: `B Σ` is symmetric (the standard criterion for the stationary process to be
  reversible).
* `matched_noise_iff`: `ℏ P⁻¹` solves the Lyapunov equation with noise covariance `2D` exactly
  when `D = ℏ𝕄⁻¹`.
* `mobility_independent`: for any symmetric mobility `Γ`, `ℏ P⁻¹` solves the Lyapunov equation
  of the drift `Γ P` with noise `2ℏΓ`; the noise is valid when `Γ` is positive definite and
  `ℏ > 0`. This solution is fixed by the conductances through `P` and does not depend on the
  repair weighting `Γ`; uniqueness for a general `Γ` is not proved here.
* `lyapunov_unique`: for positive definite `P` and positive capacities every solution of the
  Lyapunov equation of the declared law is `ℏ P⁻¹`.
* `register_stationary_unique`: the Euclidean history register satisfies these premises.
-/

namespace OPH.SourceLaw

open scoped BigOperators
open Matrix

section Algebra

variable {𝕜 : Type*} [Field 𝕜] {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Drift matrix of the declared Langevin law: `B = 𝕄⁻¹ P`. -/
def driftMatrix (M : ι → 𝕜) (P : Matrix ι ι 𝕜) : Matrix ι ι 𝕜 :=
  diagonal (fun i => (M i)⁻¹) * P

/-- The stationary covariance `ℏ P⁻¹`. -/
noncomputable def stationaryCov (hbar : 𝕜) (P : Matrix ι ι 𝕜) : Matrix ι ι 𝕜 := hbar • P⁻¹

theorem driftMatrix_transpose (M : ι → 𝕜) (P : Matrix ι ι 𝕜) (hP : Pᵀ = P) :
    (driftMatrix M P)ᵀ = P * diagonal (fun i => (M i)⁻¹) := by
  rw [driftMatrix, transpose_mul, hP, diagonal_transpose]

theorem drift_mul_stationaryCov (M : ι → 𝕜) (P : Matrix ι ι 𝕜) (hbar : 𝕜)
    (hPU : IsUnit P.det) :
    driftMatrix M P * stationaryCov hbar P = hbar • diagonal (fun i => (M i)⁻¹) := by
  rw [driftMatrix, stationaryCov, Matrix.mul_smul, Matrix.mul_assoc,
    Matrix.mul_nonsing_inv _ hPU, Matrix.mul_one]

/-- `ℏ P⁻¹` solves the Lyapunov equation of the declared law. -/
theorem lyapunov_stationary (M : ι → 𝕜) (P : Matrix ι ι 𝕜) (hbar : 𝕜) (hP : Pᵀ = P)
    (hPU : IsUnit P.det) :
    driftMatrix M P * stationaryCov hbar P + stationaryCov hbar P * (driftMatrix M P)ᵀ
      = (2 * hbar) • diagonal (fun i => (M i)⁻¹) := by
  rw [drift_mul_stationaryCov M P hbar hPU, driftMatrix_transpose M P hP, stationaryCov,
    Matrix.smul_mul, ← Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hPU, Matrix.one_mul,
    ← add_smul]
  ring_nf

/-- Detailed balance: `B Σ` is symmetric. -/
theorem detailed_balance (M : ι → 𝕜) (P : Matrix ι ι 𝕜) (hbar : 𝕜) (hPU : IsUnit P.det) :
    (driftMatrix M P * stationaryCov hbar P)ᵀ = driftMatrix M P * stationaryCov hbar P := by
  rw [drift_mul_stationaryCov M P hbar hPU, transpose_smul, diagonal_transpose]

/-- `ℏ P⁻¹` solves the Lyapunov equation with noise covariance `2D` exactly when `D = ℏ𝕄⁻¹`. -/
theorem matched_noise_iff [CharZero 𝕜] (M : ι → 𝕜) (P D : Matrix ι ι 𝕜) (hbar : 𝕜)
    (hP : Pᵀ = P) (hPU : IsUnit P.det) :
    driftMatrix M P * stationaryCov hbar P + stationaryCov hbar P * (driftMatrix M P)ᵀ
        = (2 : 𝕜) • D ↔ D = hbar • diagonal (fun i => (M i)⁻¹) := by
  rw [lyapunov_stationary M P hbar hP hPU, ← smul_smul]
  constructor
  · intro h
    exact (smul_right_injective _ (two_ne_zero : (2 : 𝕜) ≠ 0) h).symm
  · intro h
    rw [h]

/-- For any symmetric mobility `Γ`, `ℏ P⁻¹` solves the Lyapunov equation of the drift `Γ P` with
noise `2ℏΓ`. The noise is valid when `Γ` is positive definite and `ℏ > 0`; uniqueness of the
solution for a general `Γ` is not proved here. -/
theorem mobility_independent (Γ P : Matrix ι ι 𝕜) (hbar : 𝕜) (hP : Pᵀ = P)
    (hPU : IsUnit P.det) (hΓ : Γᵀ = Γ) :
    (Γ * P) * stationaryCov hbar P + stationaryCov hbar P * (Γ * P)ᵀ = (2 * hbar) • Γ := by
  rw [stationaryCov, transpose_mul, hP, hΓ, Matrix.mul_smul, Matrix.smul_mul, Matrix.mul_assoc,
    Matrix.mul_nonsing_inv _ hPU, Matrix.mul_one, ← Matrix.mul_assoc,
    Matrix.nonsing_inv_mul _ hPU, Matrix.one_mul, ← add_smul]
  ring_nf

end Algebra

section Unique

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [DecidableEq ι] in
/-- If `A` is positive definite and `A Z + Z A = 0`, then `Z = 0`. -/
theorem sylvester_zero (A Z : Matrix ι ι ℝ) (hA : A.PosDef) (h : A * Z + Z * A = 0) : Z = 0 := by
  have hq : ∀ x : ι → ℝ, x ≠ 0 → 0 < x ⬝ᵥ (A *ᵥ x) := fun x hx => by
    simpa using hA.dotProduct_mulVec_pos hx
  have hq0 : ∀ x : ι → ℝ, 0 ≤ x ⬝ᵥ (A *ᵥ x) := fun x => by
    rcases eq_or_ne x 0 with rfl | hx
    · simp
    · exact (hq x hx).le
  have hA1 : ∑ i, ∑ j, Z i j * (A * Z) i j
      = ∑ j, (fun k => Z k j) ⬝ᵥ (A *ᵥ fun k => Z k j) := by
    simp only [dotProduct, Matrix.mulVec, Matrix.mul_apply, Finset.mul_sum]
    conv_lhs => rw [Finset.sum_comm]
  have hA2 : ∑ i, ∑ j, Z i j * (Z * A) i j
      = ∑ i, (fun k => Z i k) ⬝ᵥ (A *ᵥ fun k => Z i k) := by
    refine Finset.sum_congr rfl fun i _ => ?_
    simp only [dotProduct, Matrix.mulVec, Matrix.mul_apply, Finset.mul_sum]
    conv_lhs => rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun j _ => by ring
  have hS : ∑ j, (fun k => Z k j) ⬝ᵥ (A *ᵥ fun k => Z k j)
      + ∑ i, (fun k => Z i k) ⬝ᵥ (A *ᵥ fun k => Z i k) = 0 := by
    rw [← hA1, ← hA2, ← Finset.sum_add_distrib]
    have e : ∀ i j, Z i j * (A * Z) i j + Z i j * (Z * A) i j
        = Z i j * (A * Z + Z * A) i j := by
      intro i j
      rw [Matrix.add_apply, mul_add]
    simp only [← Finset.sum_add_distrib, e, h, Matrix.zero_apply, mul_zero, Finset.sum_const_zero]
  have h1 : 0 ≤ ∑ j, (fun k => Z k j) ⬝ᵥ (A *ᵥ fun k => Z k j) :=
    Finset.sum_nonneg fun j _ => hq0 _
  have h2 : 0 ≤ ∑ i, (fun k => Z i k) ⬝ᵥ (A *ᵥ fun k => Z i k) :=
    Finset.sum_nonneg fun i _ => hq0 _
  have h3 : ∑ j, (fun k => Z k j) ⬝ᵥ (A *ᵥ fun k => Z k j) = 0 := by linarith
  have hcol : ∀ j, (fun k => Z k j) ⬝ᵥ (A *ᵥ fun k => Z k j) = 0 := fun j =>
    (Finset.sum_eq_zero_iff_of_nonneg (fun j _ => hq0 _)).1 h3 j (Finset.mem_univ _)
  ext i j
  by_contra hne
  have hc : (fun k => Z k j) ≠ 0 := fun h0 => hne (by simpa using congrFun h0 i)
  exact (hq _ hc).ne' (hcol j)

/-- For positive definite `P` and positive capacities, every solution of the Lyapunov equation of
the declared law is `ℏ P⁻¹`. -/
theorem lyapunov_unique (M : ι → ℝ) (P X : Matrix ι ι ℝ) (hbar : ℝ) (hP : P.PosDef)
    (hM : ∀ i, 0 < M i)
    (hX : driftMatrix M P * X + X * (driftMatrix M P)ᵀ
      = (2 * hbar) • diagonal (fun i => (M i)⁻¹)) :
    X = stationaryCov hbar P := by
  have hPU : IsUnit P.det := (Matrix.isUnit_iff_isUnit_det _).1 hP.isUnit
  have hPs : Pᵀ = P := by
    have := hP.isHermitian
    rwa [Matrix.IsHermitian, conjTranspose_eq_transpose_of_trivial] at this
  have hst := lyapunov_stationary M P hbar hPs hPU
  -- the difference solves the homogeneous equation
  have hY : driftMatrix M P * (X - stationaryCov hbar P)
      + (X - stationaryCov hbar P) * (driftMatrix M P)ᵀ = 0 := by
    rw [Matrix.mul_sub, Matrix.sub_mul]
    rw [show driftMatrix M P * X - driftMatrix M P * stationaryCov hbar P
        + (X * (driftMatrix M P)ᵀ - stationaryCov hbar P * (driftMatrix M P)ᵀ)
        = (driftMatrix M P * X + X * (driftMatrix M P)ᵀ)
          - (driftMatrix M P * stationaryCov hbar P
            + stationaryCov hbar P * (driftMatrix M P)ᵀ) by abel]
    rw [hX, hst, sub_self]
  -- whitening by `S = diag √M`
  set S : Matrix ι ι ℝ := diagonal (fun i => Real.sqrt (M i)) with hSdef
  set Si : Matrix ι ι ℝ := diagonal (fun i => (Real.sqrt (M i))⁻¹) with hSidef
  have hsq : ∀ i, Real.sqrt (M i) ≠ 0 := fun i => (Real.sqrt_pos.2 (hM i)).ne'
  have hSSi : S * Si = 1 := by
    rw [hSdef, hSidef, diagonal_mul_diagonal, ← diagonal_one]
    congr 1
    funext i
    exact mul_inv_cancel₀ (hsq i)
  have hSiS : Si * S = 1 := by
    rw [hSdef, hSidef, diagonal_mul_diagonal, ← diagonal_one]
    congr 1
    funext i
    exact inv_mul_cancel₀ (hsq i)
  have hSiSi : Si * Si = diagonal (fun i => (M i)⁻¹) := by
    rw [hSidef, diagonal_mul_diagonal]
    congr 1
    funext i
    rw [← mul_inv, Real.mul_self_sqrt (hM i).le]
  have hSiT : Siᵀ = Si := by rw [hSidef, diagonal_transpose]
  have hL : ∀ Y : Matrix ι ι ℝ, S * (Si * Y) = Y := fun Y => by
    rw [← Matrix.mul_assoc, hSSi, Matrix.one_mul]
  have hR : ∀ Y : Matrix ι ι ℝ, Si * (S * Y) = Y := fun Y => by
    rw [← Matrix.mul_assoc, hSiS, Matrix.one_mul]
  -- the whitened stiffness is positive definite
  have hSiU : IsUnit Si := (Matrix.isUnit_iff_isUnit_det _).2
    (Matrix.isUnit_det_of_left_inverse hSSi)
  have hPhat : (Si * P * Si).PosDef := by
    have := hP.conjTranspose_mul_mul_same (B := Si)
      (Matrix.mulVec_injective_iff_isUnit.2 hSiU)
    rwa [conjTranspose_eq_transpose_of_trivial, hSiT] at this
  set Y := X - stationaryCov hbar P with hYdef
  have hsyl : (Si * P * Si) * (S * Y * S) + (S * Y * S) * (Si * P * Si) = 0 := by
    have key : (Si * P * Si) * (S * Y * S) + (S * Y * S) * (Si * P * Si)
        = S * (driftMatrix M P * Y + Y * (driftMatrix M P)ᵀ) * S := by
      rw [driftMatrix_transpose M P hPs, driftMatrix, ← hSiSi]
      simp only [Matrix.mul_assoc, Matrix.mul_add, Matrix.add_mul, hL, hR, hSiS, Matrix.mul_one]
    rw [key, hY, Matrix.mul_zero, Matrix.zero_mul]
  have hZ := sylvester_zero (Si * P * Si) (S * Y * S) hPhat hsyl
  have hY0 : Y = 0 := by
    have : Si * (S * Y * S) * Si = Y := by
      simp only [Matrix.mul_assoc, hR, hSSi, Matrix.mul_one]
    rw [← this, hZ, Matrix.mul_zero, Matrix.zero_mul]
  rw [hYdef] at hY0
  exact sub_eq_zero.1 hY0

end Unique

section Register

variable {ι E : Type*} [Fintype ι] [DecidableEq ι] [Fintype E]

/-- The Euclidean history register satisfies the premises (nonnegative conductances and time
coefficients, positive grounding and capacities, `Δτ > 0`): its stiffness is positive definite,
so every solution of the Lyapunov equation of the declared law is `ℏ P⁻¹`. -/
theorem register_stationary_unique (src tgt : E → ι) (c : E → ℝ) (g M J : ι → ℝ) (dτ hbar : ℝ)
    (T : ℕ) (hc : ∀ e, 0 ≤ c e) (hg : ∀ i, 0 < g i) (hJ : ∀ i, 0 ≤ J i) (hd : 0 < dτ)
    (hM : ∀ i, 0 < M i)
    (X : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ)
    (hX : driftMatrix (regCap M T) (stiffnessMatrix (regSrc src T) (regTgt tgt T)
          (regCond c J dτ T) (regGround g dτ T)) * X
        + X * (driftMatrix (regCap M T) (stiffnessMatrix (regSrc src T) (regTgt tgt T)
          (regCond c J dτ T) (regGround g dτ T)))ᵀ
        = (2 * hbar) • diagonal (fun p => (regCap M T p)⁻¹)) :
    X = stationaryCov hbar (stiffnessMatrix (regSrc src T) (regTgt tgt T)
      (regCond c J dτ T) (regGround g dτ T)) :=
  lyapunov_unique (regCap M T) _ X hbar (register_posDef src tgt c g J dτ T hc hg hJ hd)
    (fun p => hM p.1) hX

end Register

#print axioms lyapunov_stationary
#print axioms detailed_balance
#print axioms matched_noise_iff
#print axioms mobility_independent
#print axioms sylvester_zero
#print axioms lyapunov_unique
#print axioms register_stationary_unique

end OPH.SourceLaw
