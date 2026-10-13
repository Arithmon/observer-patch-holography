import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Algebra.BigOperators.Pi
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Tactic

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Weighted agreement seams and their mean repair generator

On a finite carrier with capacities `M`, a seam `e` joins the readings at `src e` and `tgt e`.
Weighted agreement moves both readings to their capacity-weighted mean and leaves every other
reading unchanged. A seam with conductance `c` acts at the rate `c (1/M_i + 1/M_j)`. A grounded
seam reconciles one reading with a fixed zero reference of unbounded capacity.

Main results:

* `repairGenerator_apply`: on a seam graph without self-seams, with nonzero capacities and
  nonzero capacity sums across each seam, the mean repair generator
  `Σ_e λ_e (f − E_e f) + Σ_i (g_i / M_i)(f − G_i f)` is `M⁻¹ K f`, where `K` is the conductance
  Laplacian of the internal seams plus the grounding conductances `g`.
* `self_seam_rate_iff`: on a seam whose two ends carry the same capacity `M`, the rate equals
  `λ` exactly when the conductance equals `λ M / 2`.
* `stiffness_quadratic` and `stiffnessMatrix_posDef`: `K` is the matrix of the quadratic form
  `Σ_e c_e (f_src − f_tgt)² + Σ_i g_i f_i²`, positive definite when the conductances are
  nonnegative and every node is grounded.
-/

namespace OPH.SourceLaw

open scoped BigOperators
open Matrix

variable {𝕜 : Type*} [Field 𝕜] {ι E : Type*}

section Agreement

variable [DecidableEq ι]

/-- Capacity-weighted agreement on the seam `(i, j)`: both readings move to their common
capacity-weighted mean and every other reading is unchanged. -/
def agree (M : ι → 𝕜) (i j : ι) (f : ι → 𝕜) : ι → 𝕜 :=
  fun k => if k = i ∨ k = j then (M i * f i + M j * f j) / (M i + M j) else f k

/-- Agreement of the reading at `i` with a fixed zero reference of unbounded capacity. -/
def groundAgree (i : ι) (f : ι → 𝕜) : ι → 𝕜 := Function.update f i 0

end Agreement

/-- The seam rate `λ = c (1/M_i + 1/M_j)`. -/
def seamRate (M : ι → 𝕜) (c : 𝕜) (i j : ι) : 𝕜 := c * (1 / M i + 1 / M j)

section Terms

variable [DecidableEq ι]

/-- One seam contributes `c / M_k` times the signed difference at the reading `k`. -/
theorem seam_term (M : ι → 𝕜) (c : 𝕜) (i j k : ι) (f : ι → 𝕜) (hij : i ≠ j)
    (hi : M i ≠ 0) (hj : M j ≠ 0) (hs : M i + M j ≠ 0) :
    seamRate M c i j * (f k - agree M i j f k)
      = c * ((if k = i then f i - f j else 0) + (if k = j then f j - f i else 0)) / M k := by
  unfold seamRate agree
  by_cases hk : k = i
  · subst hk
    simp only [true_or, if_true, hij, if_false, add_zero]
    field_simp
    ring
  · by_cases hk' : k = j
    · subst hk'
      simp only [hk, or_true, if_true, if_false, zero_add]
      field_simp
      ring
    · simp [hk, hk']

/-- A grounded seam at `i` contributes `(g_i / M_i) f_i` at `i` and nothing elsewhere. -/
theorem ground_term (M g f : ι → 𝕜) (i k : ι) :
    g i / M i * (f k - groundAgree i f k) = if k = i then g k * f k / M k else 0 := by
  unfold groundAgree
  by_cases h : k = i
  · subst h
    simp
    ring
  · simp [h]

end Terms

section Stiffness

variable [DecidableEq ι] [Fintype E]

/-- Conductance Laplacian of the internal seams plus grounding, applied to the readings `f`. -/
def stiffness (src tgt : E → ι) (c : E → 𝕜) (g : ι → 𝕜) (f : ι → 𝕜) (k : ι) : 𝕜 :=
  (∑ e, c e * ((if k = src e then f (src e) - f (tgt e) else 0)
      + (if k = tgt e then f (tgt e) - f (src e) else 0))) + g k * f k

variable [Fintype ι]

/-- The mean repair generator `Σ_e λ_e (f − E_e f) + Σ_i (g_i / M_i)(f − G_i f)`. -/
def repairGenerator (M : ι → 𝕜) (src tgt : E → ι) (c : E → 𝕜) (g : ι → 𝕜)
    (f : ι → 𝕜) : ι → 𝕜 :=
  (∑ e, seamRate M (c e) (src e) (tgt e) • (f - agree M (src e) (tgt e) f))
    + ∑ i, (g i / M i) • (f - groundAgree i f)

/-- The mean repair generator is `M⁻¹ K` on every finite seam graph without self-seams, for
nonzero capacities with nonzero sums across each seam. -/
theorem repairGenerator_apply (M : ι → 𝕜) (src tgt : E → ι) (c : E → 𝕜) (g : ι → 𝕜)
    (f : ι → 𝕜) (hM : ∀ i, M i ≠ 0) (hMM : ∀ e, M (src e) + M (tgt e) ≠ 0)
    (hst : ∀ e, src e ≠ tgt e) (k : ι) :
    repairGenerator M src tgt c g f k = stiffness src tgt c g f k / M k := by
  simp only [repairGenerator, Pi.add_apply, Finset.sum_apply, Pi.smul_apply, smul_eq_mul,
    Pi.sub_apply]
  rw [Finset.sum_congr rfl (fun e _ =>
        seam_term M (c e) (src e) (tgt e) k f (hst e) (hM _) (hM _) (hMM e)),
      Finset.sum_congr rfl (fun i _ => ground_term M g f i k)]
  simp only [Finset.sum_ite_eq, Finset.mem_univ, if_true, stiffness]
  rw [add_div, Finset.sum_div]

/-- The stiffness as a matrix: its column `l` is the stiffness of the unit reading at `l`. -/
def stiffnessMatrix (src tgt : E → ι) (c : E → 𝕜) (g : ι → 𝕜) : Matrix ι ι 𝕜 :=
  Matrix.of fun k l => stiffness src tgt c g (Pi.single l 1) k

/-- The stiffness matrix acts on readings exactly as the stiffness. -/
theorem stiffnessMatrix_mulVec (src tgt : E → ι) (c : E → 𝕜) (g : ι → 𝕜) (f : ι → 𝕜) :
    stiffnessMatrix src tgt c g *ᵥ f = stiffness src tgt c g f := by
  funext k
  simp only [stiffnessMatrix, Matrix.mulVec, dotProduct, Matrix.of_apply, stiffness,
    Finset.sum_mul, add_mul, Finset.sum_add_distrib]
  rw [Finset.sum_comm]
  congr 1
  · refine Finset.sum_congr rfl fun e _ => ?_
    by_cases hs : k = src e <;> by_cases ht : k = tgt e <;>
      simp [hs, ht, Pi.single_apply, mul_add, mul_sub, Finset.sum_add_distrib,
        Finset.sum_sub_distrib, mul_comm]
  · simp [Pi.single_apply, mul_comm]

omit [Fintype ι] in
/-- The stiffness matrix is symmetric. -/
theorem stiffnessMatrix_transpose (src tgt : E → ι) (c : E → 𝕜) (g : ι → 𝕜) :
    (stiffnessMatrix src tgt c g)ᵀ = stiffnessMatrix src tgt c g := by
  ext k l
  simp only [Matrix.transpose_apply, stiffnessMatrix, Matrix.of_apply, stiffness]
  congr 1
  · refine Finset.sum_congr rfl fun e _ => ?_
    by_cases hks : k = src e <;> by_cases hkt : k = tgt e <;>
      by_cases hls : l = src e <;> by_cases hlt : l = tgt e <;>
      simp_all [Pi.single_apply, eq_comm]
  · by_cases h : k = l
    · subst h; rfl
    · simp [h, Ne.symm h]

end Stiffness

section RealStiffness

variable [DecidableEq ι] [Fintype E] [Fintype ι]

/-- The stiffness is the matrix of `Σ_e c_e (f_src − f_tgt)² + Σ_i g_i f_i²`. -/
theorem stiffness_quadratic (src tgt : E → ι) (c : E → ℝ) (g f : ι → ℝ) :
    ∑ k, f k * stiffness src tgt c g f k
      = ∑ e, c e * (f (src e) - f (tgt e)) ^ 2 + ∑ k, g k * f k ^ 2 := by
  have h1 : ∀ k, f k * stiffness src tgt c g f k
      = (∑ e, f k * (c e * ((if k = src e then f (src e) - f (tgt e) else 0)
          + (if k = tgt e then f (tgt e) - f (src e) else 0)))) + g k * f k ^ 2 := by
    intro k
    rw [stiffness, mul_add, Finset.mul_sum]
    ring
  rw [Finset.sum_congr rfl (fun k _ => h1 k), Finset.sum_add_distrib, Finset.sum_comm]
  congr 1
  refine Finset.sum_congr rfl fun e _ => ?_
  have h2 : ∀ k, f k * (c e * ((if k = src e then f (src e) - f (tgt e) else 0)
          + (if k = tgt e then f (tgt e) - f (src e) else 0)))
      = (if k = src e then c e * f (src e) * (f (src e) - f (tgt e)) else 0)
        + (if k = tgt e then c e * f (tgt e) * (f (tgt e) - f (src e)) else 0) := by
    intro k
    by_cases hs : k = src e
    · by_cases ht : k = tgt e
      · have hst : src e = tgt e := hs.symm.trans ht
        rw [if_pos hs, if_pos ht, if_pos hs, if_pos ht, ← hst]
        ring
      · rw [if_pos hs, if_neg ht, if_pos hs, if_neg ht, hs]
        ring
    · by_cases ht : k = tgt e
      · rw [if_neg hs, if_pos ht, if_neg hs, if_pos ht, ht]
        ring
      · rw [if_neg hs, if_neg ht, if_neg hs, if_neg ht]
        ring
  rw [Finset.sum_congr rfl (fun k _ => h2 k), Finset.sum_add_distrib, Finset.sum_ite_eq',
    Finset.sum_ite_eq']
  simp only [Finset.mem_univ, if_true]
  ring

/-- If every node is grounded and conductances are nonnegative, the stiffness matrix is
positive definite. -/
theorem stiffnessMatrix_posDef (src tgt : E → ι) (c : E → ℝ) (g : ι → ℝ)
    (hc : ∀ e, 0 ≤ c e) (hg : ∀ k, 0 < g k) :
    (stiffnessMatrix src tgt c g).PosDef := by
  refine Matrix.PosDef.of_dotProduct_mulVec_pos ?_ fun x hx => ?_
  · ext k l
    simpa [Matrix.conjTranspose_apply] using
      (congrFun (congrFun (stiffnessMatrix_transpose src tgt c g) l) k).symm
  · rw [stiffnessMatrix_mulVec]
    have hq := stiffness_quadratic src tgt c g x
    simp only [star_trivial, dotProduct] at *
    rw [hq]
    obtain ⟨k, hk⟩ : ∃ k, x k ≠ 0 := Function.ne_iff.mp hx
    have h1 : 0 ≤ ∑ e, c e * (x (src e) - x (tgt e)) ^ 2 :=
      Finset.sum_nonneg fun e _ => mul_nonneg (hc e) (sq_nonneg _)
    have h2 : 0 < ∑ k, g k * x k ^ 2 :=
      Finset.sum_pos' (fun j _ => mul_nonneg (hg j).le (sq_nonneg _))
        ⟨k, Finset.mem_univ _, mul_pos (hg k) (by positivity)⟩
    linarith

end RealStiffness

section SelfSeam

variable [CharZero 𝕜]

/-- On a seam whose two ends carry the same capacity `M`, the rate is `λ` exactly when the
conductance is `λ M / 2`. The name refers to a same-node seam such as the time seam between
`(i, t)` and `(i, t + 1)`, whose two ends carry the same capacity; it is not a seam from a
reading to itself. -/
theorem self_seam_rate_iff (c M lam : 𝕜) (hM : M ≠ 0) :
    lam = c * (1 / M + 1 / M) ↔ c = lam * M / 2 := by
  constructor
  · intro h
    rw [h]
    field_simp
    ring
  · intro h
    rw [h]
    field_simp
    ring

end SelfSeam

#print axioms seam_term
#print axioms repairGenerator_apply
#print axioms stiffnessMatrix_mulVec
#print axioms stiffnessMatrix_transpose
#print axioms stiffness_quadratic
#print axioms stiffnessMatrix_posDef
#print axioms self_seam_rate_iff

end OPH.SourceLaw
