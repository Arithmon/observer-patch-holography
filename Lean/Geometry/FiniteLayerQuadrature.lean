import Mathlib

/-!
# Exact finite-layer diamond quadrature

The layer sum uses continuum spatial balls, hence the cubic tent profile.
Its even and odd formulas are exact finite sums. They are not exact formulas
for the number of sites in an arbitrary finite graph diamond.
-/
namespace OPH.FiniteLayerQuadrature
open scoped BigOperators

def layerSum (K : ℕ) : ℚ :=
  ∑ j ∈ Finset.range (K + 1), ((min j (K - j) : ℕ) : ℚ) ^ 3

def layerFactor (K : ℕ) : ℚ := 32 * layerSum K / (K : ℚ) ^ 4

lemma cube_sum (m : ℕ) :
    (∑ j ∈ Finset.range m, (j : ℚ) ^ 3) = (m : ℚ)^2 * (m - 1)^2 / 4 := by
  induction m with
  | zero => norm_num
  | succ m ih =>
    rw [Finset.sum_range_succ, ih]
    push_cast
    ring

lemma reflected_cube_sum (m : ℕ) :
    (∑ i ∈ Finset.range (m + 1), ((m - i : ℕ) : ℚ)^3) =
      ∑ i ∈ Finset.range (m + 1), (i : ℚ)^3 := by
  simpa using Finset.sum_range_reflect (fun i : ℕ => (i : ℚ)^3) (m + 1)

lemma even_layer_sum (m : ℕ) :
    layerSum (2 * m) = 2 * (∑ j ∈ Finset.range m, (j : ℚ)^3) + (m : ℚ)^3 := by
  unfold layerSum
  rw [show 2 * m + 1 = m + (m + 1) by omega, Finset.sum_range_add]
  have hl : (∑ j ∈ Finset.range m, ((min j (2 * m - j) : ℕ) : ℚ)^3) =
      ∑ j ∈ Finset.range m, (j : ℚ)^3 := by
    apply Finset.sum_congr rfl
    intro j hj
    have hj' := Finset.mem_range.mp hj
    rw [Nat.min_eq_left (by omega)]
  have hr : (∑ j ∈ Finset.range (m + 1), ((min (m + j) (2 * m - (m + j)) : ℕ) : ℚ)^3) =
      ∑ j ∈ Finset.range (m + 1), ((m - j : ℕ) : ℚ)^3 := by
    apply Finset.sum_congr rfl
    intro j hj
    have hj' := Finset.mem_range.mp hj
    have h : min (m + j) (2 * m - (m + j)) = m - j := by omega
    rw [h]
  rw [hl, hr, reflected_cube_sum, Finset.sum_range_succ]
  ring

lemma odd_layer_sum (m : ℕ) :
    layerSum (2 * m + 1) = 2 * (∑ j ∈ Finset.range (m + 1), (j : ℚ)^3) := by
  unfold layerSum
  rw [show 2 * m + 1 + 1 = (m + 1) + (m + 1) by omega, Finset.sum_range_add]
  have hl : (∑ j ∈ Finset.range (m + 1), ((min j (2 * m + 1 - j) : ℕ) : ℚ)^3) =
      ∑ j ∈ Finset.range (m + 1), (j : ℚ)^3 := by
    apply Finset.sum_congr rfl
    intro j hj
    have hj' := Finset.mem_range.mp hj
    rw [Nat.min_eq_left (by omega)]
  have hr : (∑ j ∈ Finset.range (m + 1), ((min (m + 1 + j) (2 * m + 1 - (m + 1 + j)) : ℕ) : ℚ)^3) =
      ∑ j ∈ Finset.range (m + 1), ((m - j : ℕ) : ℚ)^3 := by
    apply Finset.sum_congr rfl
    intro j hj
    have hj' := Finset.mem_range.mp hj
    have h : min (m + 1 + j) (2 * m + 1 - (m + 1 + j)) = m - j := by omega
    rw [h]
  rw [hl, hr, reflected_cube_sum]
  ring

/-- Even layer separation: the normalized volume overshoots by `4/K²`. -/
theorem even_layer_factor (m : ℕ) (hm : 0 < m) :
    layerFactor (2 * m) = 1 + 4 / ((2 * m : ℕ) : ℚ)^2 := by
  have hm' : (m : ℚ) ≠ 0 := by exact_mod_cast hm.ne'
  rw [layerFactor, even_layer_sum, cube_sum]
  push_cast
  field_simp
  ring

/-- Odd layer separation, including the degenerate zero-volume `K=1`. -/
theorem odd_layer_factor (m : ℕ) :
    layerFactor (2 * m + 1) = (1 - 1 / ((2 * m + 1 : ℕ) : ℚ)^2)^2 := by
  have h : (2 * (m : ℚ) + 1) ≠ 0 := by positivity
  rw [layerFactor, odd_layer_sum, cube_sum]
  push_cast
  field_simp
  ring

#print axioms cube_sum
#print axioms reflected_cube_sum
#print axioms even_layer_sum
#print axioms odd_layer_sum
#print axioms even_layer_factor
#print axioms odd_layer_factor
end OPH.FiniteLayerQuadrature
