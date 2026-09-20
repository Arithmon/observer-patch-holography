import Geometry.GoldenSourceCountLimit

set_option autoImplicit false

/-!
# Moving finite classifiers on a measurable partition

This is a dominated-convergence transport lemma for actual finite sums.
Unlike fixed continuity-set quadrature, the selected label set may change
with the population. Its hypothesis is eventual pointwise classification
agreement on the assigned cells, not convergence of counts or integrals.
The source-net specialization must derive that agreement from its read law.
-/

namespace OPH.SourceCountTransport

open MeasureTheory Set Filter
open scoped Topology BigOperators Classical

noncomputable section

variable {X : Type*} [MeasurableSpace X]

/-- The selected cell union represented as a finite sum of indicators. -/
def selected {ι : Type*} [Fintype ι] (C : ι → Set X) (P : ι → Prop) (x : X) : ℝ :=
  ∑ i, (C i).indicator (fun _ => if P i then (1 : ℝ) else 0) x

omit [MeasurableSpace X] in
theorem selected_eq {ι : Type*} [Fintype ι] (C : ι → Set X) (P : ι → Prop)
    (hd : Pairwise (fun i j => Disjoint (C i) (C j)))
    {x : X} {i : ι} (hi : x ∈ C i) :
    selected C P x = if P i then (1 : ℝ) else 0 := by
  classical
  unfold selected
  rw [Finset.sum_eq_single i]
  · exact indicator_of_mem hi _
  · intro j hj hji
    apply indicator_of_notMem
    intro hx
    exact disjoint_left.mp (hd hji) hx hi
  · simp

/-- Moving label counts converge when the actual cell classifiers eventually
agree almost everywhere. The limiting set is measured by the same measure. -/
theorem selected_count_tendsto {ι : ℕ → Type*} [∀ n, Fintype (ι n)]
    (μ : Measure X) [IsFiniteMeasure μ]
    (C : (n : ℕ) → ι n → Set X) (P : (n : ℕ) → ι n → Prop)
    (hC : ∀ n i, MeasurableSet (C n i))
    (hd : ∀ n, Pairwise (fun i j => Disjoint (C n i) (C n j)))
    (hc : ∀ n, ∀ᵐ x ∂μ, ∃ i, x ∈ C n i)
    (A : Set X) (hA : MeasurableSet A)
    (hag : ∀ᵐ x ∂μ, ∀ᶠ n in atTop, ∀ i, x ∈ C n i → (P n i ↔ x ∈ A)) :
    Tendsto (fun n => ∑ i, μ.real (C n i) * (if P n i then (1 : ℝ) else 0))
      atTop (𝓝 (μ.real A)) := by
  classical
  have hm (n : ℕ) : AEStronglyMeasurable (selected (C n) (P n)) μ := by
    unfold selected
    convert (Finset.aestronglyMeasurable_sum Finset.univ fun i _ =>
      (aestronglyMeasurable_const : AEStronglyMeasurable
        (fun _ : X => if P n i then (1 : ℝ) else 0) μ).indicator (hC n i)) using 1
    ext x
    simp
  have hb (n : ℕ) : ∀ᵐ x ∂μ, ‖selected (C n) (P n) x‖ ≤ (1 : ℝ) := by
    filter_upwards [hc n] with x hx
    obtain ⟨i, hi⟩ := hx
    rw [selected_eq _ _ (hd n) hi]
    split_ifs <;> norm_num
  have hl : ∀ᵐ x ∂μ, Tendsto (fun n => selected (C n) (P n) x) atTop
      (𝓝 (A.indicator (fun _ => (1 : ℝ)) x)) := by
    filter_upwards [ae_all_iff.mpr hc, hag] with x hx hagx
    apply tendsto_const_nhds.congr'
    filter_upwards [hagx] with n hn
    obtain ⟨i, hi⟩ := hx n
    rw [selected_eq _ _ (hd n) hi]
    simp only [hn i hi, indicator_apply]
  have ht := tendsto_integral_of_dominated_convergence (fun _ => (1 : ℝ))
    hm (integrable_const _) hb hl
  have heq (n : ℕ) : ∫ x, selected (C n) (P n) x ∂μ =
      ∑ i, μ.real (C n i) * (if P n i then (1 : ℝ) else 0) := by
    unfold selected
    rw [integral_finset_sum]
    · apply Finset.sum_congr rfl
      intro i hi
      rw [integral_indicator (hC n i), setIntegral_const]
      rfl
    · intro i hi
      exact (integrable_const _).indicator (hC n i)
  simpa only [heq, integral_indicator hA, setIntegral_const, smul_eq_mul, mul_one] using ht

/-- Replacing dominated cell weights by full weights costs at most the
total mass discrepancy, uniformly over every selected label set. -/
theorem selected_weight_error {ι : Type*} [Fintype ι]
    (w v : ι → ℝ) (hv : ∀ i, v i ≤ w i) (P : ι → Prop) :
    |(∑ i, w i * (if P i then (1 : ℝ) else 0)) -
      ∑ i, v i * (if P i then (1 : ℝ) else 0)| ≤ ∑ i, (w i - v i) := by
  classical
  rw [← Finset.sum_sub_distrib]
  simp_rw [← sub_mul]
  apply (Finset.abs_sum_le_sum_abs _ _).trans
  apply Finset.sum_le_sum
  intro i hi
  have hn : 0 ≤ w i - v i := sub_nonneg.mpr (hv i)
  split_ifs <;> simp [abs_of_nonneg hn, hn]

/-- Product cells retain exact pair masses and all diagonal terms. -/
theorem product_weight_discrepancy {ι : Type*} [Fintype ι]
    (w v : ι → ℝ) :
    (∑ i : ι × ι, (w i.1 * w i.2 - v i.1 * v i.2)) =
      (∑ i, w i) ^ 2 - (∑ i, v i) ^ 2 := by
  rw [Finset.sum_sub_distrib]
  simp only [Fintype.sum_prod_type, ← Finset.sum_mul, ← Finset.mul_sum, pow_two]

open OPH.GoldenSourceAssignment OPH.GoldenSourceCountLimit

theorem eventCell_ae_cover {n : ℕ} (hn : 0 < n) {T L δ : ℝ}
    (hL : 0 < L) (hδ : 0 < δ) :
    ∀ᵐ x ∂spaceTimeVolume T L, ∃ i, x ∈ eventCell n T L δ i := by
  have hh : ∀ᵐ x ∂spaceTimeVolume T L, x ∈ Ico 0 T ×ˢ cube L := by
    unfold spaceTimeVolume
    rw [Measure.prod_restrict]
    exact ae_restrict_mem (measurableSet_Ico.prod
      (MeasurableSet.univ_pi fun _ => measurableSet_Ico))
  filter_upwards [hh] with x hx
  rw [← eventCell_cover hn hL hδ] at hx
  exact mem_iUnion.mp hx

/-- The clipping bound applies to arbitrary moving label predicates. -/
theorem event_selected_weight_error {n : ℕ} (hn : 0 < n) {T L δ : ℝ}
    (hT : 0 ≤ T) (hL : 0 < L) (hδ : 0 < δ)
    (P : timeIndex T δ × (Fin 3 → Fin (Nat.fib n)) → Prop) :
    |(∑ i, (δ * (L ^ 3 / (Nat.fib n : ℝ) ^ 3)) * (if P i then (1 : ℝ) else 0)) -
      ∑ i, (spaceTimeVolume T L).real (eventCell n T L δ i) *
        (if P i then (1 : ℝ) else 0)| ≤ L ^ 3 * δ := by
  apply (selected_weight_error _ _ (eventCell_mass_le hn hL hδ) P).trans
  rw [Finset.sum_sub_distrib, sum_uniform_event_weight hn,
    sum_eventCell_mass hn hT hL hδ]
  have hf := (le_div_iff₀ hδ).mp (Nat.floor_le (div_nonneg hT hδ.le))
  have hh : ((⌊T / δ⌋₊ : ℝ) + 1) * δ - T ≤ δ := by nlinarith
  nlinarith [mul_le_mul_of_nonneg_right hh (pow_nonneg hL.le 3)]

def count {ι : Type*} [Fintype ι] (P : ι → Prop) : ℕ :=
  (Finset.univ.filter P).card

theorem constant_weight_count {ι : Type*} [Fintype ι] (w : ℝ) (P : ι → Prop) :
    (∑ i, w * (if P i then (1 : ℝ) else 0)) = w * (count P : ℝ) := by
  simp only [mul_ite, mul_one, mul_zero]
  rw [← Finset.sum_filter]
  simp [count, nsmul_eq_mul, mul_comm]

/-- Actual uniform event counts converge for a derived moving classifier.
The endpoint discrepancy is proved uniformly over the classifier. -/
theorem event_selected_count_tendsto {T L : ℝ} (hT : 0 ≤ T) (hL : 0 < L)
    (δ : ℕ → ℝ) (hδ : ∀ n, 0 < δ n) (hδ0 : Tendsto δ atTop (𝓝 0))
    (P : (n : ℕ) → timeIndex T (δ n) × (Fin 3 → Fin (Nat.fib (n + 2))) → Prop)
    (A : Set (ℝ × (Fin 3 → ℝ))) (hA : MeasurableSet A)
    (hag : ∀ᵐ x ∂spaceTimeVolume T L, ∀ᶠ n in atTop,
      ∀ i, x ∈ eventCell (n + 2) T L (δ n) i → (P n i ↔ x ∈ A)) :
    Tendsto (fun n => (δ n * (L ^ 3 / (Nat.fib (n + 2) : ℝ) ^ 3)) *
        (count (P n) : ℝ)) atTop (𝓝 ((spaceTimeVolume T L).real A)) := by
  letI : IsFiniteMeasure (spaceTimeVolume T L) := spaceTimeVolume_finite T L
  let U (n : ℕ) : ℝ := ∑ i,
    (δ n * (L ^ 3 / (Nat.fib (n + 2) : ℝ) ^ 3)) * (if P n i then (1 : ℝ) else 0)
  let W (n : ℕ) : ℝ := ∑ i, (spaceTimeVolume T L).real (eventCell (n + 2) T L (δ n) i) *
    (if P n i then (1 : ℝ) else 0)
  have hw : Tendsto W atTop (𝓝 ((spaceTimeVolume T L).real A)) :=
    selected_count_tendsto _ (fun n => eventCell (n + 2) T L (δ n)) P
      (fun n => eventCell_measurable (n + 2) T L (δ n))
      (fun n => eventCell_disjoint (by omega : 0 < n + 2) hL (hδ n))
      (fun n => eventCell_ae_cover (by omega : 0 < n + 2) hL (hδ n)) A hA hag
  have hd : Tendsto (fun n => U n - W n) atTop (𝓝 0) := by
    apply tendsto_iff_norm_sub_tendsto_zero.mpr
    simp only [sub_zero, Real.norm_eq_abs]
    apply squeeze_zero (fun _ => abs_nonneg _)
      (fun n => event_selected_weight_error (by omega : 0 < n + 2) hT hL (hδ n) (P n))
    simpa using hδ0.const_mul (L ^ 3)
  have hu : Tendsto U atTop (𝓝 ((spaceTimeVolume T L).real A)) := by
    simpa only [sub_add_cancel, zero_add] using hd.add hw
  simpa only [U, constant_weight_count] using hu

#print axioms selected_count_tendsto
#print axioms selected_weight_error
#print axioms product_weight_discrepancy

end
end OPH.SourceCountTransport
