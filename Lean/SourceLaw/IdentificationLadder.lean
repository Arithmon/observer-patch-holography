import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.SchurComplement
import Mathlib.Tactic

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Which data determine the inertia

Consider `H = (pᵀ T p + qᵀ V q)/2` with inverse inertia `T` and potential `V` (units `ℏ = 1`).
A centred Gaussian `exp(−qᵀ G q / 2)` is an eigenfunction of `H` exactly when `G T G = V`
(a standard computation, not mechanized here and not written out in the companion note); its
position covariance is `C = (2G)⁻¹`. Below, "vacuum `G`" names this Gaussian; the Lean statements
concern only the relation `G T G = V` (`groundRel`) and do not assume that `G` is positive definite
or that the Gaussian is an eigenfunction. This file proves the algebra of an identification
ladder:

* `config_cov_indep`: for `β ≠ 0` and invertible `V`, the configuration marginal of the classical
  Gibbs law of `H` has covariance `(βV)⁻¹` for every invertible `T`, so configuration data do
  not determine the inertia.
* `groundRel_compatible`: every inverse inertia `T` is compatible with a given vacuum `G`
  (with `V = G T G`), so the vacuum alone does not determine the inertia.
* `inverse_inertia_unique`, `inverse_inertia_formula`: for an invertible vacuum `G`, the vacuum
  together with the potential determines the inverse inertia, `T = G⁻¹ V G⁻¹ = 4 C V C`.
* `markov_reader`: for a stationary reversible Gaussian chain with one-step transfer `R`,
  equal-time covariance `C0` and innovation `Q = C0 − R C0 Rᵀ`, the two-time data return
  `(C0 (R C0)⁻¹ C0 − R C0)⁻¹ = Q⁻¹ R`, the time coupling of the chain's precision.
-/

namespace OPH.SourceLaw

open scoped BigOperators
open Matrix

variable {𝕜 : Type*} [Field 𝕜] {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The Gaussian ground-state relation `G T G = V`. -/
def groundRel (G T V : Matrix ι ι 𝕜) : Prop := G * T * G = V

/-- Configuration data do not determine the inertia: for `β ≠ 0` and invertible `V`, the
configuration block of the inverse of the Gibbs precision `blockdiag(βV, βT)` is `(βV)⁻¹` for
every invertible `T`. -/
theorem config_cov_indep (β : 𝕜) (V T : Matrix ι ι 𝕜) (hβ : β ≠ 0) (hV : IsUnit V.det)
    (hT : IsUnit T.det) :
    ((fromBlocks (β • V) 0 0 (β • T))⁻¹).toBlocks₁₁ = (β • V)⁻¹ := by
  have hβV : IsUnit (β • V) := by
    rw [Matrix.isUnit_iff_isUnit_det, det_smul]
    exact (IsUnit.mk0 _ (pow_ne_zero _ hβ)).mul hV
  have hβT : IsUnit (β • T) := by
    rw [Matrix.isUnit_iff_isUnit_det, det_smul]
    exact (IsUnit.mk0 _ (pow_ne_zero _ hβ)).mul hT
  rw [inv_fromBlocks_zero₂₁_of_isUnit_iff (β • V) 0 (β • T) ⟨fun _ => hβT, fun _ => hβV⟩]
  rfl

omit [DecidableEq ι] in
/-- The vacuum alone does not determine the inertia: every inverse inertia `T` is compatible
with the vacuum `G` for the potential `G T G`. -/
theorem groundRel_compatible (G T : Matrix ι ι 𝕜) : groundRel G T (G * T * G) := rfl

/-- For an invertible vacuum `G`, two different inverse inertias with that vacuum need different
potentials. -/
theorem groundRel_family_injective (G : Matrix ι ι 𝕜) (hG : IsUnit G.det) :
    Function.Injective (fun T : Matrix ι ι 𝕜 => G * T * G) := by
  intro T T' h
  simp only at h
  have := congrArg (fun X => G⁻¹ * X * G⁻¹) h
  simp only [← Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hG, Matrix.one_mul] at this
  simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hG] using this

/-- For an invertible vacuum `G`, the vacuum and the potential determine the inverse inertia. -/
theorem inverse_inertia_unique (G T T' V : Matrix ι ι 𝕜) (hG : IsUnit G.det)
    (h : groundRel G T V) (h' : groundRel G T' V) : T = T' :=
  groundRel_family_injective G hG (h.trans h'.symm)

/-- Explicitly, `T = G⁻¹ V G⁻¹`, i.e. `T = 4 C V C` with `C = 2⁻¹ G⁻¹`, the inverse of `2G`
(the position covariance of the Gaussian vacuum). The value is an inverse inertia. -/
theorem inverse_inertia_formula [CharZero 𝕜] (G T V : Matrix ι ι 𝕜) (hG : IsUnit G.det)
    (h : groundRel G T V) :
    T = G⁻¹ * V * G⁻¹ ∧ ((2 : 𝕜) • G) * ((2 : 𝕜)⁻¹ • G⁻¹) = 1
      ∧ T = (4 : 𝕜) • (((2 : 𝕜)⁻¹ • G⁻¹) * V * ((2 : 𝕜)⁻¹ • G⁻¹)) := by
  have h1 : T = G⁻¹ * V * G⁻¹ := by
    rw [← h]
    simp only [← Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hG, Matrix.one_mul]
    rw [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hG, Matrix.mul_one]
  refine ⟨h1, ?_, ?_⟩
  · rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul, Matrix.mul_nonsing_inv _ hG]
    norm_num
  · rw [h1]
    simp only [Matrix.smul_mul, Matrix.mul_smul, smul_smul]
    norm_num

/-- Two-time data return the time coupling of a stationary reversible Gaussian chain. -/
theorem markov_reader (C0 R : Matrix ι ι 𝕜) (hC : IsUnit C0.det) (hR : IsUnit R.det)
    (hrev : R * C0 = C0 * Rᵀ) :
    (C0 * (R * C0)⁻¹ * C0 - R * C0)⁻¹ = (C0 - R * C0 * Rᵀ)⁻¹ * R := by
  have key : C0 * (R * C0)⁻¹ * C0 - R * C0 = R⁻¹ * (C0 - R * C0 * Rᵀ) := by
    rw [Matrix.mul_inv_rev, Matrix.mul_sub]
    rw [show R * C0 * Rᵀ = R * (C0 * Rᵀ) by rw [Matrix.mul_assoc], ← Matrix.mul_assoc R⁻¹ R,
      Matrix.nonsing_inv_mul _ hR, Matrix.one_mul, ← hrev]
    rw [← Matrix.mul_assoc C0 C0⁻¹, Matrix.mul_nonsing_inv _ hC, Matrix.one_mul]
  rw [key, Matrix.mul_inv_rev, Matrix.nonsing_inv_nonsing_inv _ hR]

#print axioms config_cov_indep
#print axioms groundRel_compatible
#print axioms groundRel_family_injective
#print axioms inverse_inertia_unique
#print axioms inverse_inertia_formula
#print axioms markov_reader

end OPH.SourceLaw
