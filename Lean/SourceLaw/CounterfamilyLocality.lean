import SourceLaw.HistoryRegister
import Mathlib.LinearAlgebra.Matrix.IsDiag

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Node locality of the time coupling and the `I + εA` kinetic family

In whitened coordinates the kinetic counterfamily is `H_ε = (pᵀ(I + εA)p + qᵀAq)/2`. In the
original coordinates its inverse inertia is `M⁻¹ + ε M⁻¹ V M⁻¹` and its inertia, the matrix
that couples consecutive Euclidean slabs, is `M (M + εV)⁻¹ M`.

Main results:

* `counterfamilyInertia_isDiag_iff`: with nonzero capacities, invertible `M + εV` and a nonzero
  off-diagonal entry of `V`, the inertia is diagonal exactly when `ε = 0`. A same-node time seam
  cannot carry any other member.
* `counterfamilyInertia_mul_inverseInertia`: the inertia is the inverse of the original-coordinate
  inverse inertia `M⁻¹ + ε M⁻¹ V M⁻¹`.
* `counterfamily_whitened`: with positive capacities, conjugating the whitened kinetic matrix
  `I + εA` by `M^{-1/2}` on both sides gives that inverse inertia.
* `counterfamily_uniform_register_iff`: under the same premises, with `Δτ ≠ 0` and with at least
  one time seam per node (`T > 0`), a register with node-local time seams of uniform rate `2/Δτ`
  realizes the inertia of a family member exactly when `ε = 0` and `J = M`. Inside the family,
  node locality alone forces `ε = 0` and then `J = M`
  (`counterfamily_register_representable_iff`, `counterfamilyInertia_zero`); the proof below uses
  the rate conjunct for `J = M`, but the forward direction does not need it. The rate conjunct is
  a conclusion of the reverse direction.
-/

namespace OPH.SourceLaw

open scoped BigOperators
open Matrix

variable {𝕜 : Type*} [Field 𝕜] {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Inertia of the kinetic counterfamily in original coordinates: `M (M + εV)⁻¹ M`. -/
noncomputable def counterfamilyInertia (M : ι → 𝕜) (V : Matrix ι ι 𝕜) (ε : 𝕜) :
    Matrix ι ι 𝕜 :=
  diagonal M * (diagonal M + ε • V)⁻¹ * diagonal M

/-- Inverse inertia of the counterfamily in original coordinates: `M⁻¹ + ε M⁻¹ V M⁻¹`. -/
def counterfamilyInverseInertia (M : ι → 𝕜) (V : Matrix ι ι 𝕜) (ε : 𝕜) : Matrix ι ι 𝕜 :=
  diagonal (fun i => (M i)⁻¹)
    + ε • (diagonal (fun i => (M i)⁻¹) * V * diagonal (fun i => (M i)⁻¹))

theorem diagonal_mul_diagonal_inv (M : ι → 𝕜) (hM : ∀ i, M i ≠ 0) :
    diagonal M * diagonal (fun i => (M i)⁻¹) = 1 := by
  rw [diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext i
  exact mul_inv_cancel₀ (hM i)

theorem diagonal_inv_mul_diagonal (M : ι → 𝕜) (hM : ∀ i, M i ≠ 0) :
    diagonal (fun i => (M i)⁻¹) * diagonal M = 1 := by
  rw [diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext i
  exact inv_mul_cancel₀ (hM i)

theorem det_diagonal_isUnit (M : ι → 𝕜) (hM : ∀ i, M i ≠ 0) : IsUnit (diagonal M).det := by
  rw [det_diagonal]
  exact IsUnit.mk0 _ (Finset.prod_ne_zero_iff.mpr fun i _ => hM i)

/-- At `ε = 0` the family's coupling `M (M + εV)⁻¹ M` reduces to `diag M` (a computation within
the declared family). -/
theorem counterfamilyInertia_zero (M : ι → 𝕜) (V : Matrix ι ι 𝕜) (hM : ∀ i, M i ≠ 0) :
    counterfamilyInertia M V 0 = diagonal M := by
  simp only [counterfamilyInertia, zero_smul, add_zero]
  rw [Matrix.mul_nonsing_inv _ (det_diagonal_isUnit M hM), Matrix.one_mul]

/-- The inertia is the inverse of the original-coordinate inverse inertia. -/
theorem counterfamilyInertia_mul_inverseInertia (M : ι → 𝕜) (V : Matrix ι ι 𝕜) (ε : 𝕜)
    (hM : ∀ i, M i ≠ 0) (hU : IsUnit (diagonal M + ε • V).det) :
    counterfamilyInertia M V ε * counterfamilyInverseInertia M V ε = 1 := by
  have h1 := diagonal_mul_diagonal_inv M hM
  have h2 := diagonal_inv_mul_diagonal M hM
  have hInv : diagonal (fun i => (M i)⁻¹) * (diagonal M + ε • V) * diagonal (fun i => (M i)⁻¹)
      = counterfamilyInverseInertia M V ε := by
    rw [Matrix.mul_add, Matrix.add_mul, h2, Matrix.one_mul, Matrix.mul_smul, Matrix.smul_mul]
    rfl
  rw [← hInv, counterfamilyInertia]
  calc diagonal M * (diagonal M + ε • V)⁻¹ * diagonal M
        * (diagonal (fun i => (M i)⁻¹) * (diagonal M + ε • V) * diagonal (fun i => (M i)⁻¹))
      = diagonal M * ((diagonal M + ε • V)⁻¹ * (diagonal M * diagonal (fun i => (M i)⁻¹))
          * (diagonal M + ε • V)) * diagonal (fun i => (M i)⁻¹) := by
        simp only [Matrix.mul_assoc]
    _ = diagonal M * ((diagonal M + ε • V)⁻¹ * (diagonal M + ε • V))
          * diagonal (fun i => (M i)⁻¹) := by rw [h1, Matrix.mul_one]
    _ = 1 := by rw [Matrix.nonsing_inv_mul _ hU, Matrix.mul_one, h1]

/-- Node locality: with nonzero capacities, invertible `M + εV` and a nonzero off-diagonal entry of
`V`, the counterfamily inertia is diagonal exactly when `ε = 0`. -/
theorem counterfamilyInertia_isDiag_iff (M : ι → 𝕜) (V : Matrix ι ι 𝕜) (ε : 𝕜)
    (hM : ∀ i, M i ≠ 0) (hU : IsUnit (diagonal M + ε • V).det)
    (hV : ∃ i j, i ≠ j ∧ V i j ≠ 0) :
    (counterfamilyInertia M V ε).IsDiag ↔ ε = 0 := by
  have h1 := diagonal_mul_diagonal_inv M hM
  have h2 := diagonal_inv_mul_diagonal M hM
  constructor
  · intro hD
    obtain ⟨i, j, hij, hVij⟩ := hV
    -- `(M + εV)⁻¹ = M⁻¹ J M⁻¹` is a product of diagonal matrices.
    have hX : (diagonal M + ε • V)⁻¹
        = diagonal (fun i => (M i)⁻¹) * counterfamilyInertia M V ε
            * diagonal (fun i => (M i)⁻¹) := by
      rw [counterfamilyInertia]
      calc (diagonal M + ε • V)⁻¹
          = (diagonal (fun i => (M i)⁻¹) * diagonal M) * (diagonal M + ε • V)⁻¹
              * (diagonal M * diagonal (fun i => (M i)⁻¹)) := by
            rw [h1, h2, Matrix.one_mul, Matrix.mul_one]
        _ = diagonal (fun i => (M i)⁻¹) * (diagonal M * (diagonal M + ε • V)⁻¹ * diagonal M)
              * diagonal (fun i => (M i)⁻¹) := by
            simp only [Matrix.mul_assoc]
    have hJ := (isDiag_iff_diagonal_diag _).1 hD
    have hXd : (diagonal M + ε • V)⁻¹
        = diagonal (fun k => (M k)⁻¹ * (counterfamilyInertia M V ε).diag k * (M k)⁻¹) := by
      rw [hX, ← hJ, diagonal_mul_diagonal, diagonal_mul_diagonal]
      simp only [diag_diagonal]
    -- From `X * (M + εV) = 1` with `X` diagonal, the off-diagonal entries of `M + εV` vanish.
    have hXN : (diagonal M + ε • V)⁻¹ * (diagonal M + ε • V) = 1 :=
      Matrix.nonsing_inv_mul _ hU
    rw [hXd] at hXN
    have hrow : ∀ a b, (M a)⁻¹ * (counterfamilyInertia M V ε).diag a * (M a)⁻¹
        * (diagonal M + ε • V) a b = (1 : Matrix ι ι 𝕜) a b := by
      intro a b
      rw [← hXN, diagonal_mul]
    have hxi : (M i)⁻¹ * (counterfamilyInertia M V ε).diag i * (M i)⁻¹ ≠ 0 := by
      intro h0
      have := hrow i i
      rw [h0, zero_mul, one_apply_eq] at this
      exact zero_ne_one this
    have hNij : (diagonal M + ε • V) i j = 0 := by
      have := hrow i j
      rw [one_apply_ne hij] at this
      exact (mul_eq_zero.mp this).resolve_left hxi
    have hNij' : (diagonal M + ε • V) i j = ε * V i j := by
      simp [diagonal_apply_ne _ hij]
    rw [hNij'] at hNij
    exact (mul_eq_zero.mp hNij).resolve_right hVij
  · rintro rfl
    rw [counterfamilyInertia_zero M V hM]
    exact isDiag_diagonal _

/-- With nonzero capacities, invertible `M + εV` and a nonzero off-diagonal entry of `V`, a
same-node (diagonal) time coupling carries a family member exactly when `ε = 0`. -/
theorem counterfamily_register_representable_iff (M : ι → 𝕜) (V : Matrix ι ι 𝕜) (ε : 𝕜)
    (hM : ∀ i, M i ≠ 0) (hU : IsUnit (diagonal M + ε • V).det)
    (hV : ∃ i j, i ≠ j ∧ V i j ≠ 0) :
    (∃ J : ι → 𝕜, counterfamilyInertia M V ε = diagonal J) ↔ ε = 0 := by
  rw [← counterfamilyInertia_isDiag_iff M V ε hM hU hV]
  constructor
  · rintro ⟨J, hJ⟩
    rw [hJ]
    exact isDiag_diagonal _
  · intro h
    exact ⟨_, ((isDiag_iff_diagonal_diag _).1 h).symm⟩

section Register

variable [CharZero 𝕜] {E : Type*}

/-- With nonzero capacities, invertible `M + εV`, a nonzero off-diagonal entry of `V`, `Δτ ≠ 0`
and `T > 0`, a register with same-node time seams of uniform rate `2/Δτ` realizes the inertia of
a counterfamily member exactly when `ε = 0` and `J = M`. Node locality alone forces
`ε = 0` and `J = M`; this proof of the forward direction uses the rate conjunct, which is not
needed there. -/
theorem counterfamily_uniform_register_iff (src tgt : E → ι) (c : E → 𝕜) (M J : ι → 𝕜)
    (V : Matrix ι ι 𝕜) (ε dτ : 𝕜) (T : ℕ) (hM : ∀ i, M i ≠ 0) (hd : dτ ≠ 0) (hT : 0 < T)
    (hU : IsUnit (diagonal M + ε • V).det) (hV : ∃ i j, i ≠ j ∧ V i j ≠ 0) :
    (counterfamilyInertia M V ε = diagonal J ∧
        ∀ (i : ι) (t : Fin T), seamRate (regCap M T) (regCond c J dτ T (.inr (i, t)))
          (regSrc src T (.inr (i, t))) (regTgt tgt T (.inr (i, t))) = 2 / dτ)
      ↔ (ε = 0 ∧ J = M) := by
  constructor
  · rintro ⟨hJ, hrate⟩
    refine ⟨(counterfamily_register_representable_iff M V ε hM hU hV).1 ⟨J, hJ⟩, ?_⟩
    funext i
    exact (timeSeam_rate_eq_iff src tgt c M J dτ T i ⟨0, hT⟩ (hM i) hd).1 (hrate i ⟨0, hT⟩)
  · rintro ⟨rfl, rfl⟩
    refine ⟨counterfamilyInertia_zero J V hM, fun i t => ?_⟩
    exact (timeSeam_rate_eq_iff src tgt c J J dτ T i t (hM i) hd).2 rfl

end Register

section Whitened

/-- Conjugating the whitened kinetic matrix `I + εA`, `A = M^{-1/2} V M^{-1/2}`, by `M^{-1/2}` on
both sides gives the original-coordinate inverse inertia `M⁻¹ + ε M⁻¹ V M⁻¹` (positive
capacities). -/
theorem counterfamily_whitened (M : ι → ℝ) (V : Matrix ι ι ℝ) (ε : ℝ) (hM : ∀ i, 0 < M i) :
    diagonal (fun i => (Real.sqrt (M i))⁻¹)
        * (1 + ε • (diagonal (fun i => (Real.sqrt (M i))⁻¹) * V
            * diagonal (fun i => (Real.sqrt (M i))⁻¹)))
        * diagonal (fun i => (Real.sqrt (M i))⁻¹)
      = counterfamilyInverseInertia M V ε := by
  have hS : diagonal (fun i => (Real.sqrt (M i))⁻¹) * diagonal (fun i => (Real.sqrt (M i))⁻¹)
      = diagonal (fun i => (M i)⁻¹) := by
    rw [diagonal_mul_diagonal]
    congr 1
    funext i
    rw [← mul_inv, Real.mul_self_sqrt (hM i).le]
  rw [Matrix.mul_add, Matrix.add_mul, Matrix.mul_one, hS, Matrix.mul_smul, Matrix.smul_mul]
  have hX : diagonal (fun i => (Real.sqrt (M i))⁻¹)
        * (diagonal (fun i => (Real.sqrt (M i))⁻¹) * V * diagonal (fun i => (Real.sqrt (M i))⁻¹))
        * diagonal (fun i => (Real.sqrt (M i))⁻¹)
      = diagonal (fun i => (M i)⁻¹) * V * diagonal (fun i => (M i)⁻¹) := by
    rw [← hS]
    simp only [Matrix.mul_assoc]
  rw [hX]
  rfl

end Whitened

#print axioms counterfamilyInertia_zero
#print axioms counterfamilyInertia_mul_inverseInertia
#print axioms counterfamilyInertia_isDiag_iff
#print axioms counterfamily_register_representable_iff
#print axioms counterfamily_uniform_register_iff
#print axioms counterfamily_whitened

end OPH.SourceLaw
