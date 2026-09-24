import Mathlib.Tactic

set_option autoImplicit false

namespace OPH.M1QuantumTransport

/-!
Finite reductions for the native-time quantum transport theorem.
The compressed velocity resolution, nuclear-norm characterization, full
spectral limit and operator-algebra delay theorem are analytic proofs in
the accompanying derivation. These lemmas do not formalize those arguments.
-/

/-- The exact three-dimensional dot-product deficit is a sum of squares. -/
theorem dot_deficit (x y z a b c : ℝ) :
    (x ^ 2 + y ^ 2 + z ^ 2) * (a ^ 2 + b ^ 2 + c ^ 2) -
      (x * a + y * b + z * c) ^ 2 =
      (x * b - y * a) ^ 2 + (x * c - z * a) ^ 2 + (y * c - z * b) ^ 2 := by
  ring

/-- Every scalar trace witness of a positive velocity resolution obeys its budget. -/
theorem weighted_trace_budget {ι : Type*} [Fintype ι] (w x : ι → ℝ) (c : ℝ)
    (hw : ∀ i, 0 ≤ w i) (hs : ∑ i, w i = 1) (hx : ∀ i, x i ≤ c) :
    ∑ i, w i * x i ≤ c := by
  calc
    _ ≤ ∑ i, w i * c := Finset.sum_le_sum (fun i _ => mul_le_mul_of_nonneg_left (hx i) (hw i))
    _ = c := by rw [← Finset.sum_mul, hs]; ring

/-- An isotropic three-component trace certificate forces the factor three. -/
theorem isotropic_speed (v c : ℝ) (h : 3 * v ≤ c) : v ≤ c / 3 := by
  linarith

/-- A positive luminal speed cannot obey the three-component flight budget. -/
theorem shared_speed_impossible (c : ℝ) (hc : 0 < c) : ¬ 3 * c ≤ c := by
  linarith

/-- Three charged orthogonal flights attain the sharp Weyl speed. -/
theorem charged_three_flights (a c : ℝ) (ha : a ≠ 0) (hc : c ≠ 0) :
    a / (3 * a / c) = c / 3 := by
  field_simp

/-- Nonnegative onsite overhead cannot improve the native propagation speed. -/
theorem overhead_lowers_speed (a t h : ℝ) (ha : 0 ≤ a) (ht : 0 < t) (hh : 0 ≤ h) :
    a / (t + h) ≤ a / t := by
  exact div_le_div_of_nonneg_left ha ht (by linarith)

/-- The three-axis remainder is bounded by the Euclidean momentum norm. -/
theorem three_axis_remainder (x y z : ℝ) :
    (x + y + z) ^ 2 ≤ 3 * (x ^ 2 + y ^ 2 + z ^ 2) := by
  nlinarith [sq_nonneg (x - y), sq_nonneg (x - z), sq_nonneg (y - z)]

/-- The error accumulated over charged ticks has the explicit first-order envelope. -/
theorem charged_error (a c T K n : ℝ) (ha : 0 < a)
    (hn : n ≤ c * T / (3 * a)) :
    n * (2 * a ^ 2 * K ^ 2) + a * K ≤ (2 * c * T / 3) * a * K ^ 2 + a * K := by
  have h := mul_le_mul_of_nonneg_right hn (show 0 ≤ 2 * a ^ 2 * K ^ 2 by positivity)
  have he : c * T / (3 * a) * (2 * a ^ 2 * K ^ 2) = (2 * c * T / 3) * a * K ^ 2 := by
    field_simp
  rw [he] at h
  linarith

/-- A pure two-component spinor saturates the Weyl current bound. -/
theorem pauli_current (a b d e : ℝ) :
    (2 * (a * d + b * e)) ^ 2 + (2 * (a * e - b * d)) ^ 2 +
      (a ^ 2 + b ^ 2 - d ^ 2 - e ^ 2) ^ 2 =
      (a ^ 2 + b ^ 2 + d ^ 2 + e ^ 2) ^ 2 := by
  ring

/-- A zero sum of squared interaction coefficients forces each coefficient to vanish. -/
theorem zero_interaction_coefficients {ι : Type*} [Fintype ι] (h : ι → ℝ)
    (hz : 256 * ∑ i, h i ^ 2 = 0) : ∀ i, h i = 0 := by
  have hs : ∑ i, h i ^ 2 = 0 := by linarith
  intro i
  have hi : h i ^ 2 ≤ ∑ j, h j ^ 2 := Finset.single_le_sum (fun j _ => sq_nonneg (h j)) (Finset.mem_univ i)
  rw [hs] at hi
  nlinarith [sq_nonneg (h i)]

/-- Saturation forces every positive-weight velocity witness to attain its bound. -/
theorem saturated_weighted_bound {ι : Type*} [Fintype ι] (w f : ι → ℝ) (c : ℝ)
    (hw : ∀ i, 0 < w i) (hf : ∀ i, f i ≤ c) (hs : ∑ i, w i = 1)
    (he : ∑ i, w i * f i = c) : ∀ i, f i = c := by
  have hz : ∑ i, w i * (c - f i) = 0 := by
    simp_rw [mul_sub]
    rw [Finset.sum_sub_distrib, ← Finset.sum_mul, hs, he]
    ring
  intro i
  have hn : ∀ j, 0 ≤ w j * (c - f j) := fun j => mul_nonneg (le_of_lt (hw j)) (sub_nonneg.mpr (hf j))
  have hi : w i * (c - f i) ≤ ∑ j, w j * (c - f j) :=
    Finset.single_le_sum (fun j _ => hn j) (Finset.mem_univ i)
  rw [hz] at hi
  nlinarith [hw i, hf i]

/-- The spherical-cap budget leaves a uniform speed gap for a fixed finite stencil. -/
theorem finite_stencil_gap (N ratio : ℝ) (hN : 0 < N)
    (hcap : 2 ≤ N * (1 - ratio)) : ratio ≤ 1 - 2 / N := by
  have h : 2 / N ≤ 1 - ratio := (div_le_iff₀ hN).mpr (by nlinarith)
  linarith

/-- Approaching the full cone speed requires a diverging direction count. -/
theorem direction_resource_bound (N ratio error : ℝ) (hN : 0 ≤ N) (he : 0 < error)
    (hcap : 2 ≤ N * (1 - ratio)) (hr : 1 - error ≤ ratio) : 2 / error ≤ N := by
  have h := mul_le_mul_of_nonneg_left (show 1 - ratio ≤ error by linarith) hN
  exact (div_le_iff₀ he).mpr (by nlinarith)

end OPH.M1QuantumTransport
