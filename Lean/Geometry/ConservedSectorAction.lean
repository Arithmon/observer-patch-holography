import Geometry.SourceActionMeasure

set_option autoImplicit false

/-!
# Conserved-sector freedom in a finite action reference

The state space, transition kernel and sector map are supplied independently.
If every nonzero transition preserves the sector, reweighting a stationary
reference by any positive sector factor preserves stationarity. Conditioning
on one sector removes that factor, and every sector-only potential has zero
transition-gradient form. These are finite algebraic consequences, not a
selection of physical energy or a claim that a growing archive is stationary.

`SourceActionMeasure.ratio_along_walk` gives uniqueness of reciprocity weights
up to scale on each connected sector. Its globally normalized uniqueness
theorem requires a connected full support, which a nonconstant conserved
sector map does not supply. This module makes the surviving freedom explicit.
-/

namespace OPH.ConservedSectorAction

open scoped BigOperators

variable {S Q : Type*}

/-- A multiplicative change of reference that is constant within each sector. -/
def reweight (q : S → Q) (pi : S → ℝ) (factor : Q → ℝ) (x : S) : ℝ :=
  factor (q x) * pi x

/-- Sector reweighting preserves every detailed-balance equation. -/
theorem reweight_reciprocity (q : S → Q) (pi : S → ℝ)
    (factor : Q → ℝ) (K : S → S → ℝ)
    (hsector : ∀ x y, K x y ≠ 0 → q x = q y)
    (hbalance : ∀ x y, pi x * K x y = pi y * K y x) (x y : S) :
    reweight q pi factor x * K x y = reweight q pi factor y * K y x := by
  by_cases h : K x y = 0
  · have hz : pi y * K y x = 0 := by simpa [h] using (hbalance x y).symm
    simp only [reweight, h, mul_zero]
    rw [mul_assoc, hz, mul_zero]
  · rw [reweight, reweight, hsector x y h]
    exact OPH.SourceActionMeasure.common_scale_preserves_reciprocity (hbalance x y)

/-- Stationarity needs sector preservation but does not need reversibility. -/
theorem reweight_stationary [Fintype S] (q : S → Q) (pi : S → ℝ)
    (factor : Q → ℝ) (K : S → S → ℝ)
    (hsector : ∀ x y, K x y ≠ 0 → q x = q y)
    (hstationary : ∀ y, ∑ x, pi x * K x y = pi y) (y : S) :
    ∑ x, reweight q pi factor x * K x y = reweight q pi factor y := by
  calc
    ∑ x, reweight q pi factor x * K x y =
        ∑ x, factor (q y) * (pi x * K x y) := by
      apply Finset.sum_congr rfl
      intro x _
      by_cases h : K x y = 0
      · simp [reweight, h]
      · simp only [reweight, hsector x y h, mul_assoc]
    _ = factor (q y) * ∑ x, pi x * K x y := (Finset.mul_sum _ _ _).symm
    _ = reweight q pi factor y := by rw [hstationary]; rfl

/-- Positive sector factors preserve a faithful finite reference. -/
theorem reweight_pos (q : S → Q) (pi : S → ℝ) (factor : Q → ℝ)
    (hpi : ∀ x, 0 < pi x) (hfactor : ∀ a, 0 < factor a) (x : S) :
    0 < reweight q pi factor x := mul_pos (hfactor (q x)) (hpi x)

/-- A faithful reweighted reference has a positive normalizing mass. -/
theorem reweight_mass_pos [Fintype S] [Nonempty S] (q : S → Q) (pi : S → ℝ)
    (factor : Q → ℝ) (hpi : ∀ x, 0 < pi x) (hfactor : ∀ a, 0 < factor a) :
    0 < ∑ x, reweight q pi factor x := by
  exact Finset.sum_pos (fun x _ => reweight_pos q pi factor hpi hfactor x)
    Finset.univ_nonempty

/-- Dividing by the total mass gives a normalized reference. -/
theorem normalized_mass [Fintype S] (w : S → ℝ) (hmass : ∑ x, w x ≠ 0) :
    ∑ x, w x / (∑ y, w y) = 1 := by
  rw [← Finset.sum_div, div_self hmass]

/-- Normalizing a stationary reference preserves stationarity. -/
theorem normalized_stationary [Fintype S] (w : S → ℝ) (K : S → S → ℝ)
    (hstationary : ∀ y, ∑ x, w x * K x y = w y) (z : ℝ) (y : S) :
    ∑ x, (w x / z) * K x y = w y / z := by
  calc
    ∑ x, (w x / z) * K x y = ∑ x, (w x * K x y) / z := by
      apply Finset.sum_congr rfl
      intro x _
      ring
    _ = (∑ x, w x * K x y) / z := (Finset.sum_div _ _ _).symm
    _ = w y / z := by rw [hstationary]

/-- Conditioning on a specified sector cancels its common weight factor.
The finite set and sector membership are explicit premises. -/
theorem conditional_reweight (q : S → Q) (pi : S → ℝ) (factor : Q → ℝ)
    (sector : Finset S) (a : Q)
    (hsector : ∀ x ∈ sector, q x = a) (hfactor : factor a ≠ 0)
    (x : S) (hx : x ∈ sector) :
    reweight q pi factor x / (∑ y ∈ sector, reweight q pi factor y) =
      pi x / (∑ y ∈ sector, pi y) := by
  have hs : ∑ y ∈ sector, reweight q pi factor y =
      factor a * ∑ y ∈ sector, pi y := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro y hy
    simp only [reweight, hsector y hy]
  rw [hs, reweight, hsector x hx]
  exact mul_div_mul_left (pi x) (∑ y ∈ sector, pi y) hfactor

/-- No potential depending only on a preserved sector has native gradient cost. -/
theorem sector_potential_gradient_zero [Fintype S] (q : S → Q) (pi : S → ℝ)
    (K : S → S → ℝ) (potential : Q → ℝ)
    (hsector : ∀ x y, K x y ≠ 0 → q x = q y) :
    ∑ x, ∑ y, pi x * K x y * (potential (q y) - potential (q x)) ^ 2 = 0 := by
  apply Finset.sum_eq_zero
  intro x _
  apply Finset.sum_eq_zero
  intro y _
  by_cases h : K x y = 0
  · simp [h]
  · simp [hsector x y h]

#print axioms reweight_reciprocity
#print axioms reweight_stationary
#print axioms conditional_reweight
#print axioms sector_potential_gradient_zero

end OPH.ConservedSectorAction
