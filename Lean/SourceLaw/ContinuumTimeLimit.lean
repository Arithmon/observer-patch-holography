import SourceLaw.TransferStructure
import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
import Mathlib.Analysis.Calculus.Deriv.Slope

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The slab-thickness limit of the reconstructed mode

Transfer-matrix reconstruction of a stationary mode with transfer `r = exp(−Δτ Ω_Δτ)` and
equal-time variance `c0` gives the one-mode Hamiltonian `(𝒦 p² + 𝒱 q²)/2` with
`𝒦 = Δτ Ω_Δτ / sinh(Δτ Ω_Δτ)` and `𝒱 = Ω_Δτ sinh(Δτ Ω_Δτ) / Δτ` (standard; recorded in the
companion note). This file proves the finite-`Δτ` identities and the `Δτ → 0⁺` limits:

* `reconstructed_freq`, `reconstructed_var`: `√(𝒦 𝒱) = Ω_Δτ` and `√(𝒦/𝒱)/2 = c0`.
* `latticeFreq_tendsto`: `Ω_Δτ → √λ`.
* `kineticCoeff_tendsto`, `potentialCoeff_tendsto`: `𝒦 → 1` and `𝒱 → λ`, the canonical
  coefficients of `H_can` in whitened coordinates.
-/

namespace OPH.SourceLaw

open Filter Topology Real

/-- Kinetic coefficient of the reconstructed mode. -/
noncomputable def kineticCoeff (dτ lam : ℝ) : ℝ :=
  dτ * latticeFreq dτ lam / Real.sinh (dτ * latticeFreq dτ lam)

/-- Potential coefficient of the reconstructed mode. -/
noncomputable def potentialCoeff (dτ lam : ℝ) : ℝ :=
  latticeFreq dτ lam * Real.sinh (dτ * latticeFreq dτ lam) / dτ

section Identities

variable {dτ lam : ℝ}

theorem dτ_mul_latticeFreq (hd : 0 < dτ) :
    dτ * latticeFreq dτ lam = 2 * Real.arsinh (dτ * Real.sqrt lam / 2) := by
  unfold latticeFreq
  field_simp

theorem sinh_latticeFreq (hd : 0 < dτ) :
    Real.sinh (dτ * latticeFreq dτ lam)
      = 2 * (dτ * Real.sqrt lam / 2) * Real.sqrt (1 + (dτ * Real.sqrt lam / 2) ^ 2) := by
  rw [dτ_mul_latticeFreq hd, Real.sinh_two_mul, Real.sinh_arsinh, Real.cosh_arsinh]

theorem latticeFreq_pos (hd : 0 < dτ) (hl : 0 < lam) : 0 < latticeFreq dτ lam := by
  unfold latticeFreq
  have : 0 < Real.arsinh (dτ * Real.sqrt lam / 2) :=
    Real.arsinh_pos_iff.2 (by have := Real.sqrt_pos.2 hl; positivity)
  positivity

theorem sinh_latticeFreq_pos (hd : 0 < dτ) (hl : 0 < lam) :
    0 < Real.sinh (dτ * latticeFreq dτ lam) :=
  Real.sinh_pos_iff.2 (mul_pos hd (latticeFreq_pos hd hl))

/-- The reconstructed frequency is the lattice frequency. -/
theorem reconstructed_freq (hd : 0 < dτ) (hl : 0 < lam) :
    Real.sqrt (kineticCoeff dτ lam * potentialCoeff dτ lam) = latticeFreq dτ lam := by
  have hs := (sinh_latticeFreq_pos hd hl).ne'
  have hΩ := latticeFreq_pos hd hl
  have : kineticCoeff dτ lam * potentialCoeff dτ lam = latticeFreq dτ lam ^ 2 := by
    unfold kineticCoeff potentialCoeff
    field_simp
    try ring
  rw [this, Real.sqrt_sq hΩ.le]

/-- The ground-state position variance of the reconstructed mode Hamiltonian is the equal-time
variance `c0`. -/
theorem reconstructed_var (hd : 0 < dτ) (hl : 0 < lam) :
    Real.sqrt (kineticCoeff dτ lam / potentialCoeff dτ lam) / 2 = equalTimeVar dτ lam := by
  have hs := sinh_latticeFreq_pos hd hl
  have hΩ := latticeFreq_pos hd hl
  have hq : kineticCoeff dτ lam / potentialCoeff dτ lam
      = (dτ / Real.sinh (dτ * latticeFreq dτ lam)) ^ 2 := by
    unfold kineticCoeff potentialCoeff
    field_simp
    try ring
  rw [hq, Real.sqrt_sq (by positivity)]
  -- `c0 = Δτ / (r⁻¹ − r) = Δτ / (2 sinh(Δτ Ω))`
  have hr := transferEig_pos hd hl
  have hexp := exp_neg_latticeFreq hd hl
  have hinv : (transferEig dτ lam)⁻¹ = Real.exp (dτ * latticeFreq dτ lam) := by
    rw [← hexp, Real.exp_neg, inv_inv]
  have h2s : (transferEig dτ lam)⁻¹ - transferEig dτ lam
      = 2 * Real.sinh (dτ * latticeFreq dτ lam) := by
    rw [hinv, ← hexp, Real.sinh_eq]
    ring
  have hc : equalTimeVar dτ lam = dτ / ((transferEig dτ lam)⁻¹ - transferEig dτ lam) := by
    unfold equalTimeVar
    have : 0 < 1 - transferEig dτ lam ^ 2 := one_sub_transferEig_sq_pos hd hl
    field_simp
  rw [hc, h2s]
  field_simp

end Identities

/-- `arsinh u / u → 1` as `u → 0`. -/
theorem arsinh_slope_tendsto :
    Tendsto (fun u : ℝ => u⁻¹ * Real.arsinh u) (𝓝[≠] 0) (𝓝 1) := by
  have h := (Real.hasDerivAt_arsinh 0).tendsto_slope_zero
  simpa [Real.arsinh_zero] using h

theorem scaled_tendsto_punctured {s : ℝ} (hs : 0 < s) :
    Tendsto (fun dτ : ℝ => dτ * s / 2) (𝓝[>] 0) (𝓝[≠] 0) := by
  apply tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within
  · have : Tendsto (fun dτ : ℝ => dτ * s / 2) (𝓝 0) (𝓝 (0 * s / 2)) :=
      ((continuous_id.mul continuous_const).div_const 2).tendsto 0
    simpa using this.mono_left nhdsWithin_le_nhds
  · filter_upwards [self_mem_nhdsWithin] with dτ hdτ
    have h0 : (0 : ℝ) < dτ := hdτ
    exact (by positivity : dτ * s / 2 ≠ 0)

/-- The lattice frequency tends to the canonical frequency `√λ`. -/
theorem latticeFreq_tendsto {lam : ℝ} (hl : 0 < lam) :
    Tendsto (fun dτ => latticeFreq dτ lam) (𝓝[>] 0) (𝓝 (Real.sqrt lam)) := by
  have hs : 0 < Real.sqrt lam := Real.sqrt_pos.2 hl
  have hcomp := (arsinh_slope_tendsto.comp (scaled_tendsto_punctured hs)).const_mul (Real.sqrt lam)
  rw [mul_one] at hcomp
  refine hcomp.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with dτ hdτ
  have h0 : (0 : ℝ) < dτ := hdτ
  simp only [Function.comp, latticeFreq]
  field_simp

/-- The kinetic coefficient tends to one. -/
theorem kineticCoeff_tendsto {lam : ℝ} (hl : 0 < lam) :
    Tendsto (fun dτ => kineticCoeff dτ lam) (𝓝[>] 0) (𝓝 1) := by
  have hs : 0 < Real.sqrt lam := Real.sqrt_pos.2 hl
  have hu := scaled_tendsto_punctured hs
  have hslope := arsinh_slope_tendsto.comp hu
  have hu0 : Tendsto (fun dτ : ℝ => dτ * Real.sqrt lam / 2) (𝓝[>] 0) (𝓝 0) :=
    hu.mono_right nhdsWithin_le_nhds
  have hroot : Tendsto (fun dτ : ℝ => Real.sqrt (1 + (dτ * Real.sqrt lam / 2) ^ 2))
      (𝓝[>] 0) (𝓝 1) := by
    have : Tendsto (fun dτ : ℝ => Real.sqrt (1 + (dτ * Real.sqrt lam / 2) ^ 2)) (𝓝[>] 0)
        (𝓝 (Real.sqrt (1 + 0 ^ 2))) :=
      ((continuous_sqrt.tendsto _).comp ((tendsto_const_nhds.add (hu0.pow 2))))
    simpa using this
  have hq := hslope.div hroot one_ne_zero
  rw [div_one] at hq
  refine hq.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with dτ hdτ
  have h0 : (0 : ℝ) < dτ := hdτ
  have hsq : 0 < Real.sqrt (1 + (dτ * Real.sqrt lam / 2) ^ 2) := Real.sqrt_pos.2 (by positivity)
  simp only [Pi.div_apply, Function.comp, kineticCoeff]
  rw [sinh_latticeFreq h0, dτ_mul_latticeFreq h0]
  have hu : 0 < dτ * Real.sqrt lam / 2 := by positivity
  field_simp

/-- The potential coefficient tends to `λ`. -/
theorem potentialCoeff_tendsto {lam : ℝ} (hl : 0 < lam) :
    Tendsto (fun dτ => potentialCoeff dτ lam) (𝓝[>] 0) (𝓝 lam) := by
  have hs : 0 < Real.sqrt lam := Real.sqrt_pos.2 hl
  have hu0 : Tendsto (fun dτ : ℝ => dτ * Real.sqrt lam / 2) (𝓝[>] 0) (𝓝 0) :=
    (scaled_tendsto_punctured hs).mono_right nhdsWithin_le_nhds
  have hroot : Tendsto (fun dτ : ℝ => Real.sqrt (1 + (dτ * Real.sqrt lam / 2) ^ 2))
      (𝓝[>] 0) (𝓝 1) := by
    have : Tendsto (fun dτ : ℝ => Real.sqrt (1 + (dτ * Real.sqrt lam / 2) ^ 2)) (𝓝[>] 0)
        (𝓝 (Real.sqrt (1 + 0 ^ 2))) :=
      ((continuous_sqrt.tendsto _).comp ((tendsto_const_nhds.add (hu0.pow 2))))
    simpa using this
  have hprod := ((latticeFreq_tendsto hl).mul_const (Real.sqrt lam)).mul hroot
  rw [mul_one, Real.mul_self_sqrt hl.le] at hprod
  refine hprod.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with dτ hdτ
  have h0 : (0 : ℝ) < dτ := hdτ
  simp only [potentialCoeff]
  rw [sinh_latticeFreq h0]
  field_simp
  try ring

#print axioms reconstructed_freq
#print axioms reconstructed_var
#print axioms latticeFreq_tendsto
#print axioms kineticCoeff_tendsto
#print axioms potentialCoeff_tendsto

end OPH.SourceLaw
