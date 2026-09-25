import Geometry.SourceNetLayeredOrder
import Mathlib.Algebra.Order.BigOperators.Group.Finset

set_option autoImplicit false

/-!
# Relative action measure and record density

These statements concern a supplied linear response operator and its positive
diagonal symmetrizers. They do not identify an action density with physical
metric volume, select a source law, or fix an absolute action normalization.
The finite q5 certificate checks the connected support and all reciprocity
equations independently of its spanning-tree reconstruction algorithm.
-/

namespace OPH.SourceActionMeasure

open scoped BigOperators
open OPH.SourceNetLayeredOrder

theorem cross_of_reciprocity {mi mj ni nj a b : ℝ} (ha : a ≠ 0)
    (hm : mi * a = mj * b) (hn : ni * a = nj * b) :
    mi * nj = mj * ni := by
  apply mul_right_cancel₀ ha
  calc
    (mi * nj) * a = (mi * a) * nj := by ring
    _ = (mj * b) * nj := by rw [hm]
    _ = (nj * b) * mj := by ring
    _ = (ni * a) * mj := by rw [hn]
    _ = (mj * ni) * a := by ring

theorem ratio_of_reciprocity {mi mj ni nj a b : ℝ} (ha : a ≠ 0)
    (hni : ni ≠ 0) (hnj : nj ≠ 0)
    (hm : mi * a = mj * b) (hn : ni * a = nj * b) :
    mi / ni = mj / nj :=
  (div_eq_div_iff hni hnj).2 (cross_of_reciprocity ha hm hn)

theorem ratio_along_walk {S : Type*} {edge : S → S → Prop}
    (A : S → S → ℝ) (m n : S → ℝ)
    (hn : ∀ i, n i ≠ 0)
    (ha : ∀ i j, edge i j → A i j ≠ 0)
    (hmr : ∀ i j, edge i j → m i * A i j = m j * A j i)
    (hnr : ∀ i j, edge i j → n i * A i j = n j * A j i)
    {k : ℕ} {i j : S} (walk : Walk edge k i j) :
    m i / n i = m j / n j := by
  induction walk with
  | nil _ => rfl
  | @cons k i j z hij _ ih =>
    exact (ratio_of_reciprocity (ha i j hij) (hn i) (hn j)
      (hmr i j hij) (hnr i j hij)).trans ih

theorem normalized_symmetrizer_unique {S : Type*} [Fintype S]
    {edge : S → S → Prop} (A : S → S → ℝ) (m n : S → ℝ)
    (root : S) (connected : ∀ i, ∃ k, Walk edge k root i)
    (hn : ∀ i, n i ≠ 0)
    (ha : ∀ i j, edge i j → A i j ≠ 0)
    (hmr : ∀ i j, edge i j → m i * A i j = m j * A j i)
    (hnr : ∀ i j, edge i j → n i * A i j = n j * A j i)
    (hmass : ∑ i, m i = 1) (hnmass : ∑ i, n i = 1) : m = n := by
  let c := m root / n root
  have point : ∀ i, m i = c * n i := by
    intro i
    obtain ⟨k, walk⟩ := connected i
    exact (div_eq_iff (hn i)).1
      (ratio_along_walk A m n hn ha hmr hnr walk).symm
  have hc : c = 1 := by
    calc
      c = c * ∑ i, n i := by rw [hnmass, mul_one]
      _ = ∑ i, m i := by
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl (fun i _ => (point i).symm)
      _ = 1 := hmass
  funext i
  rw [point i, hc, one_mul]

theorem common_scale_preserves_reciprocity {mi mj a b c : ℝ}
    (h : mi * a = mj * b) : (c * mi) * a = (c * mj) * b := by
  calc
    (c * mi) * a = c * (mi * a) := by ring
    _ = c * (mj * b) := by rw [h]
    _ = (c * mj) * b := by ring

theorem common_action_scale_cancels {m k c : ℝ} (hm : m ≠ 0) (hc : c ≠ 0) :
    (c * k) / (c * m) = k / m := by
  field_simp

theorem production_volume_gauge {rate volume factor : ℝ} (hf : factor ≠ 0) :
    (rate / factor) * (factor * volume) = rate * volume := by
  field_simp

theorem constitutive_volume_gauge {density volume lapse factor : ℝ}
    (hf : factor ≠ 0) (hl : lapse ≠ 0) :
    (density / factor) * (factor * volume) / lapse = density * volume / lapse := by
  field_simp

theorem summed_production_discrepancy {S : Type*} (s : Finset S) (hs : s.Nonempty)
    (v n : S → ℝ) (scale : ℝ)
    (h : ∀ i ∈ s, 0 ≤ scale * v i - n i ∧ scale * v i - n i < 1) :
    0 ≤ ∑ i ∈ s, (scale * v i - n i) ∧
      (∑ i ∈ s, (scale * v i - n i)) < s.card := by
  constructor
  · exact Finset.sum_nonneg (fun i hi => (h i hi).1)
  · simpa using Finset.sum_lt_sum_of_nonempty hs (fun i hi => (h i hi).2)

theorem calibrated_discrepancy_identity {S : Type*} (s : Finset S)
    (v n : S → ℝ) {scale : ℝ} (hc : scale ≠ 0) :
    (∑ i ∈ s, v i) - (∑ i ∈ s, n i) / scale =
      (∑ i ∈ s, (scale * v i - n i)) / scale := by
  rw [Finset.sum_sub_distrib, ← Finset.mul_sum]
  field_simp

theorem calibrated_finite_error {S : Type*} (s : Finset S) (hs : s.Nonempty)
    (v n : S → ℝ) {scale : ℝ} (hc : 0 < scale)
    (h : ∀ i ∈ s, 0 ≤ scale * v i - n i ∧ scale * v i - n i < 1) :
    0 ≤ (∑ i ∈ s, v i) - (∑ i ∈ s, n i) / scale ∧
      (∑ i ∈ s, v i) - (∑ i ∈ s, n i) / scale < s.card / scale := by
  rw [calibrated_discrepancy_identity s v n (ne_of_gt hc)]
  obtain ⟨lower, upper⟩ := summed_production_discrepancy s hs v n scale h
  exact ⟨div_nonneg lower (le_of_lt hc), div_lt_div_of_pos_right upper hc⟩

theorem conformal_volume_is_squared_kinetic_factor (sigma : ℝ) :
    sigma ^ 4 = (sigma ^ 2) ^ 2 := by ring

theorem interior_row_sum_calibrates_conformal_factor {S : Type*}
    (s : Finset S) (stiffness : S → ℝ) {u mass : ℝ}
    (hu : u ≠ 0) (hm : mass ≠ 0)
    (hrow : ∑ j ∈ s, stiffness j = u ^ 2 * mass) :
    (∑ j ∈ s, stiffness j / (u * mass)) = u := by
  rw [← Finset.sum_div, hrow]
  field_simp

theorem calibrated_relative_profile {wi wr mi mr ui ur scale : ℝ}
    (hmi : mi ≠ 0) (hmr : mr ≠ 0) (hur : ur ≠ 0) (hc : scale ≠ 0)
    (hi : wi = scale * mi * ui) (hr : wr = scale * mr * ur) :
    (wi / mi) / (wr / mr) * ur = ui := by
  rw [hi, hr]
  field_simp

theorem kinetic_density_need_not_equal_volume_density :
    (1 : ℚ) / (1 + 2) ≠ 1 / (1 + 2 ^ 2) := by norm_num

end OPH.SourceActionMeasure
