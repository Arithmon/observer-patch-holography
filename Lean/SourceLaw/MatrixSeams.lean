import SourceLaw.WeightedSeams

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The seam graph of a symmetric matrix

Any symmetric real matrix `V` is the stiffness matrix of a seam graph: the seams are the ordered
pairs of distinct nodes, the seam `(i, j)` carries conductance `-V i j / 2` (so the two seams
joining `i` and `j` together carry `-V i j`), and every node is grounded by its row sum
(`stiffnessMatrix_matSeams`). When the off-diagonal entries are nonpositive the conductances are
nonnegative (`matCond_nonneg`), and when the row sums are positive every node is grounded, which
is what the positivity theorems of the register need. Seams joining nodes that `V` does not
couple carry conductance zero and contribute nothing.
-/

namespace OPH.SourceLaw

open scoped BigOperators
open Matrix

/-- Seams of the matrix graph: ordered pairs of distinct nodes. -/
abbrev MatSeam (ι : Type*) := {p : ι × ι // p.1 ≠ p.2}

variable {ι : Type*}

/-- First endpoint of a matrix seam. -/
def matSrc (p : MatSeam ι) : ι := p.1.1

/-- Second endpoint of a matrix seam. -/
def matTgt (p : MatSeam ι) : ι := p.1.2

/-- Each ordered seam carries half of the off-diagonal coupling. -/
noncomputable def matCond (V : Matrix ι ι ℝ) (p : MatSeam ι) : ℝ := -V p.1.1 p.1.2 / 2

/-- Grounding by row sums. -/
noncomputable def rowSum [Fintype ι] (V : Matrix ι ι ℝ) (i : ι) : ℝ := ∑ j, V i j

theorem matSrc_ne_tgt (p : MatSeam ι) : matSrc p ≠ matTgt p := p.2

theorem matCond_nonneg (V : Matrix ι ι ℝ) (hV : ∀ i j, i ≠ j → V i j ≤ 0) (p : MatSeam ι) :
    0 ≤ matCond V p := by
  have := hV _ _ p.2
  unfold matCond
  linarith

variable [Fintype ι] [DecidableEq ι]

/-- For symmetric `V`, the stiffness of the matrix seam graph, grounded by the row sums, applies `V`. -/
theorem stiffness_matSeams (V : Matrix ι ι ℝ) (hV : Vᵀ = V) (f : ι → ℝ) (k : ι) :
    stiffness matSrc matTgt (matCond V) (rowSum V) f k = (V *ᵥ f) k := by
  have hsym : ∀ i j, V i j = V j i := fun i j => by
    conv_lhs => rw [← hV]
    rfl
  let F : ι × ι → ℝ := fun p => -V p.1 p.2 / 2 *
    ((if k = p.1 then f p.1 - f p.2 else 0) + (if k = p.2 then f p.2 - f p.1 else 0))
  have hF : ∀ p : ι × ι, (if p.1 ≠ p.2 then F p else 0) = F p := by
    intro p
    by_cases h : p.1 = p.2
    · simp only [F, h, ne_eq, not_true_eq_false, sub_self, ite_self, add_zero,
        mul_zero]
    · rw [if_pos h]
  have hsub : (∑ e : MatSeam ι, matCond V e * ((if k = matSrc e then f (matSrc e) - f (matTgt e)
      else 0) + (if k = matTgt e then f (matTgt e) - f (matSrc e) else 0)))
      = ∑ p : ι × ι, F p := by
    change ∑ e : MatSeam ι, F e.1 = _
    rw [← Finset.sum_subtype (Finset.univ.filter fun p : ι × ι => p.1 ≠ p.2) (by simp) F,
      Finset.sum_filter]
    exact Finset.sum_congr rfl fun p _ => hF p
  have h1 : ∑ i, ∑ j, (if k = i then -V i j / 2 * (f i - f j) else 0)
      = ∑ j, -V k j / 2 * (f k - f j) := by
    rw [Finset.sum_comm]
    simp
  have h2 : ∑ i, ∑ j, (if k = j then -V i j / 2 * (f j - f i) else 0)
      = ∑ i, -V k i / 2 * (f k - f i) := by
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.sum_ite_eq, if_pos (Finset.mem_univ _), hsym i k]
  have hpair : ∑ p : ι × ι, F p = ∑ j, -V k j * (f k - f j) := by
    rw [Fintype.sum_prod_type]
    simp only [F, mul_add, mul_ite, mul_zero, Finset.sum_add_distrib, h1, h2]
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [stiffness, hsub, hpair, rowSum, Finset.sum_mul, ← Finset.sum_add_distrib]
  simp only [mulVec, dotProduct]
  exact Finset.sum_congr rfl fun j _ => by ring

/-- A symmetric matrix is the stiffness matrix of its seam graph grounded by its row sums. -/
theorem stiffnessMatrix_matSeams (V : Matrix ι ι ℝ) (hV : Vᵀ = V) :
    stiffnessMatrix matSrc matTgt (matCond V) (rowSum V) = V := by
  rw [Matrix.ext_iff_mulVec]
  intro f
  funext k
  rw [stiffnessMatrix_mulVec, stiffness_matSeams V hV f k]

/-- A symmetric matrix with nonpositive off-diagonal entries and positive row sums is positive
definite. -/
theorem posDef_of_matSeams (V : Matrix ι ι ℝ) (hV : Vᵀ = V) (hoff : ∀ i j, i ≠ j → V i j ≤ 0)
    (hrow : ∀ i, 0 < rowSum V i) : V.PosDef := by
  rw [← stiffnessMatrix_matSeams V hV]
  exact stiffnessMatrix_posDef _ _ _ _ (matCond_nonneg V hoff) hrow

#print axioms stiffness_matSeams
#print axioms stiffnessMatrix_matSeams
#print axioms posDef_of_matSeams

end OPH.SourceLaw
