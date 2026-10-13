import SourceLaw.MatrixLimits
import SourceLaw.GoldenFamily
import Mathlib.Analysis.SpecificLimits.Basic

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Slab-thickness limits along the golden refinement family

The register over a level-`q` carrier reconstructs, at slab thickness `Δτ`, a quadratic
Hamiltonian whose vacuum has equal-time covariance `C0(Δτ) = c0(A_q)` in capacity-whitened
coordinates. This file proves:

* `vacuum_charFun_tendsto`: for every positive definite `A` and every test vector `g`, the Gaussian
  characteristic function `exp(−½ gᵀ C0(Δτ) g)` tends to `exp(−½ gᵀ (½ A^{-1/2}) g)` as
  `Δτ → 0⁺`.
* `response_kernel_tendsto`: the response kernel `Ω_Δτ⁻¹ sin(s Ω_Δτ)` of the reconstructed
  frequency tends to `ω⁻¹ sin(s ω)` with `ω = A^{1/2}`.
* `golden_continuum_time_limit`: at every level `q = n + 2` of the golden family, the reconstructed
  frequency tends to `ω_q = A_q^{1/2}`, the kinetic and potential coefficients to `1` and `A_q`, the
  equal-time covariance to the canonical vacuum covariance `½ ω_q⁻¹`, and the two limits above hold.
* `joint_limit`: a family of positive definite `A_k` on index types that change with `k`, test
  vectors `g_k` with `g_kᵀ A_k^{1/2} g_k` bounded and `g_kᵀ (½ A_k^{-1/2}) g_k → c`, and any
  `Δτ_k → 0⁺` give `g_kᵀ C0(Δτ_k; A_k) g_k → c`.
* `joint_limit_L2`: the same conclusion when only the squared lengths `g_kᵀ g_k` are bounded, by
  the spectrum-uniform bound `equalTimeCov_err_L2`; no energy bound is used, so sampled test
  vectors whose energy `g_kᵀ A_k^{1/2} g_k` grows (detectors with jumps) are covered.
  `golden_joint_limit_L2` states it along the golden family for any sequence of levels.

Negative control: `joint_limit_bound_needed` exhibits one mode with `λ_k = (k + 1)²`, weight
`2(k + 1)` and `Δτ_k = 1/(k + 1)`, where the canonical value is `1` at every `k` but the weighted
equal-time variance is `2/√5` at every `k`; there the weights are unbounded, so both the bound on
`g_kᵀ A_k^{1/2} g_k` and the bound on `g_kᵀ g_k` fail.
-/

namespace OPH.SourceLaw

open Filter Topology Matrix

noncomputable section

section General

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {A : Matrix ι ι ℝ} (hA : A.PosDef)

/-- The Gaussian characteristic function of the reconstructed vacuum tends to that of the
canonical vacuum. -/
theorem vacuum_charFun_tendsto (g : ι → ℝ) :
    Tendsto (fun dτ => Real.exp (-(1 / 2) * (g ⬝ᵥ (equalTimeCov hA dτ *ᵥ g)))) (𝓝[>] 0)
      (𝓝 (Real.exp (-(1 / 2) * (g ⬝ᵥ (((1 / 2 : ℝ) • posDefInvSqrt hA) *ᵥ g))))) := by
  have hq : Continuous fun C : Matrix ι ι ℝ => Real.exp (-(1 / 2) * (g ⬝ᵥ (C *ᵥ g))) :=
    Real.continuous_exp.comp (continuous_const.mul
      (continuous_const.dotProduct (continuous_id.matrix_mulVec continuous_const)))
  have h := equalTimeCov_tendsto hA
  rw [half_invSqrt_eq hA] at h
  have h2 := (hq.tendsto ((1 / 2 : ℝ) • posDefInvSqrt hA)).comp h
  exact h2

/-- The response kernel `Ω_Δτ⁻¹ sin(s Ω_Δτ)` tends to `ω⁻¹ sin(s ω)`, `ω = A^{1/2}`. -/
theorem response_kernel_tendsto (s : ℝ) :
    Tendsto (fun dτ => hA.1.cfc (fun x => Real.sin (s * latticeFreq dτ x) / latticeFreq dτ x))
      (𝓝[>] 0) (𝓝 (hA.1.cfc (fun x => Real.sin (s * Real.sqrt x) / Real.sqrt x))) := by
  refine cfc_tendsto hA.1
    (f := fun dτ x => Real.sin (s * latticeFreq dτ x) / latticeFreq dτ x) fun i => ?_
  have hl := eig_pos hA i
  have hw := latticeFreq_tendsto hl
  have hs : Real.sqrt (hA.1.eigenvalues i) ≠ 0 := (Real.sqrt_pos.2 hl).ne'
  exact ((Real.continuous_sin.tendsto _).comp (hw.const_mul s)).div hw hs

end General

section Golden

/-- At every level `q = n + 2` of the golden family, as `Δτ → 0⁺`: frequency, coefficients,
equal-time covariance, vacuum characteristic function and response kernel tend to their canonical
values. -/
theorem golden_continuum_time_limit (n : ℕ) {m2 : ℝ} (hm : 0 < m2) (g : Site n → ℝ) (s : ℝ) :
    Tendsto (fun dτ => (golden_whitened_posDef n hm).1.cfc (latticeFreq dτ)) (𝓝[>] 0)
        (𝓝 (posDefSqrt (golden_whitened_posDef n hm)))
      ∧ Tendsto (fun dτ => (golden_whitened_posDef n hm).1.cfc (kineticCoeff dτ)) (𝓝[>] 0) (𝓝 1)
      ∧ Tendsto (fun dτ => (golden_whitened_posDef n hm).1.cfc (potentialCoeff dτ)) (𝓝[>] 0)
        (𝓝 (whitenedPacket (goldenGrid (n + 2)) m2))
      ∧ Tendsto (fun dτ => equalTimeCov (golden_whitened_posDef n hm) dτ) (𝓝[>] 0)
        (𝓝 ((1 / 2 : ℝ) • posDefInvSqrt (golden_whitened_posDef n hm)))
      ∧ Tendsto (fun dτ => Real.exp (-(1 / 2)
          * (g ⬝ᵥ (equalTimeCov (golden_whitened_posDef n hm) dτ *ᵥ g)))) (𝓝[>] 0)
        (𝓝 (Real.exp (-(1 / 2)
          * (g ⬝ᵥ (((1 / 2 : ℝ) • posDefInvSqrt (golden_whitened_posDef n hm)) *ᵥ g)))))
      ∧ Tendsto (fun dτ => (golden_whitened_posDef n hm).1.cfc
          (fun x => Real.sin (s * latticeFreq dτ x) / latticeFreq dτ x)) (𝓝[>] 0)
        (𝓝 ((golden_whitened_posDef n hm).1.cfc
          (fun x => Real.sin (s * Real.sqrt x) / Real.sqrt x))) := by
  have hC := equalTimeCov_tendsto (golden_whitened_posDef n hm)
  rw [half_invSqrt_eq (golden_whitened_posDef n hm)] at hC
  exact ⟨latticeFreq_cfc_tendsto (golden_whitened_posDef n hm),
    kineticCoeff_cfc_tendsto (golden_whitened_posDef n hm),
    potentialCoeff_cfc_tendsto (golden_whitened_posDef n hm), hC,
    vacuum_charFun_tendsto (golden_whitened_posDef n hm) g,
    response_kernel_tendsto (golden_whitened_posDef n hm) s⟩

end Golden

section Joint

variable {κ : ℕ → Type*} [∀ k, Fintype (κ k)] [∀ k, DecidableEq (κ k)]

/-- Simultaneous limit over a refinement family with level-dependent index types. -/
theorem joint_limit {A : ∀ k, Matrix (κ k) (κ k) ℝ} (hA : ∀ k, (A k).PosDef) (g : ∀ k, κ k → ℝ)
    {dτ : ℕ → ℝ} (hd : ∀ k, 0 < dτ k) (hdt : Tendsto dτ atTop (𝓝 0)) {B : ℝ}
    (hB : ∀ k, g k ⬝ᵥ ((hA k).1.cfc Real.sqrt *ᵥ g k) ≤ B) {c : ℝ}
    (hc : Tendsto (fun k => g k ⬝ᵥ ((hA k).1.cfc (fun x => 1 / (2 * Real.sqrt x)) *ᵥ g k)) atTop
      (𝓝 c)) :
    Tendsto (fun k => g k ⬝ᵥ (equalTimeCov (hA k) (dτ k) *ᵥ g k)) atTop (𝓝 c) := by
  have hbound : ∀ k, ‖g k ⬝ᵥ (equalTimeCov (hA k) (dτ k) *ᵥ g k)
      - g k ⬝ᵥ ((hA k).1.cfc (fun x => 1 / (2 * Real.sqrt x)) *ᵥ g k)‖ ≤ dτ k ^ 2 / 16 * B := by
    intro k
    rw [Real.norm_eq_abs]
    exact (equalTimeCov_quadratic_err (hA k) (hd k) (g k)).trans
      (mul_le_mul_of_nonneg_left (hB k) (by positivity))
  have hlim : Tendsto (fun k => dτ k ^ 2 / 16 * B) atTop (𝓝 0) := by
    have := ((hdt.pow 2).div_const 16).mul_const B
    simpa using this
  have hdiff := squeeze_zero_norm hbound hlim
  have := hdiff.add hc
  simpa using this

/-- Simultaneous limit with a bound on the squared lengths `g_kᵀ g_k` in place of the energy bound,
by the spectrum-uniform error bound `equalTimeCov_err_L2`. -/
theorem joint_limit_L2 {A : ∀ k, Matrix (κ k) (κ k) ℝ} (hA : ∀ k, (A k).PosDef) (g : ∀ k, κ k → ℝ)
    {dτ : ℕ → ℝ} (hd : ∀ k, 0 < dτ k) (hdt : Tendsto dτ atTop (𝓝 0)) {B : ℝ}
    (hB : ∀ k, g k ⬝ᵥ g k ≤ B) {c : ℝ}
    (hc : Tendsto (fun k => g k ⬝ᵥ ((hA k).1.cfc (fun x => 1 / (2 * Real.sqrt x)) *ᵥ g k)) atTop
      (𝓝 c)) :
    Tendsto (fun k => g k ⬝ᵥ (equalTimeCov (hA k) (dτ k) *ᵥ g k)) atTop (𝓝 c) := by
  have herr : Tendsto (fun k => g k ⬝ᵥ (equalTimeCov (hA k) (dτ k) *ᵥ g k)
      - g k ⬝ᵥ ((hA k).1.cfc (fun x => 1 / (2 * Real.sqrt x)) *ᵥ g k)) atTop (𝓝 0) := by
    have hbound : Tendsto (fun k => dτ k / 4 * B) atTop (𝓝 0) := by
      simpa using (hdt.div_const 4).mul_const B
    refine squeeze_zero_norm (fun k => ?_) hbound
    rw [Real.norm_eq_abs]
    calc _ ≤ dτ k / 4 * (g k ⬝ᵥ g k) := equalTimeCov_err_L2 (hA k) (hd k) (g k)
      _ ≤ dτ k / 4 * B := mul_le_mul_of_nonneg_left (hB k) (by linarith [hd k])
  simpa using herr.add hc

/-- Negative control: without a bound on `g_kᵀ A_k^{1/2} g_k` or on `g_kᵀ g_k` the simultaneous
limit fails (the weights `2(k + 1)` of this one-mode witness are unbounded). -/
theorem joint_limit_bound_needed :
    ∃ lam dτ w : ℕ → ℝ, (∀ k, 0 < lam k) ∧ (∀ k, 0 < dτ k) ∧ Tendsto dτ atTop (𝓝 0)
      ∧ (∀ k, w k * (1 / (2 * Real.sqrt (lam k))) = 1)
      ∧ (∀ k, w k * equalTimeVar (dτ k) (lam k) = 2 / Real.sqrt 5)
      ∧ 2 / Real.sqrt 5 ≠ 1 := by
  have hsq : ∀ k : ℕ, Real.sqrt (((k : ℝ) + 1) ^ 2) = (k : ℝ) + 1 := fun k =>
    Real.sqrt_sq (by positivity)
  have h5 : Real.sqrt (1 + 1 / 4) = Real.sqrt 5 / 2 := by
    rw [show (1 + 1 / 4 : ℝ) = 5 / 4 by norm_num, Real.sqrt_div' 5 (by norm_num : (0 : ℝ) ≤ 4),
      show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num : (0 : ℝ) ≤ 2)]
  have hs5 : 0 < Real.sqrt 5 := Real.sqrt_pos.2 (by norm_num)
  refine ⟨fun k => ((k : ℝ) + 1) ^ 2, fun k => 1 / ((k : ℝ) + 1), fun k => 2 * ((k : ℝ) + 1),
    fun k => by positivity, fun k => by positivity, tendsto_one_div_add_atTop_nhds_zero_nat,
    fun k => ?_, fun k => ?_, ?_⟩
  · rw [hsq]
    have hk : (0 : ℝ) < (k : ℝ) + 1 := by positivity
    field_simp
  · have hk : (0 : ℝ) < (k : ℝ) + 1 := by positivity
    rw [equalTimeVar_closed (by positivity) (by positivity), hsq]
    have e : (1 / ((k : ℝ) + 1)) ^ 2 * ((k : ℝ) + 1) ^ 2 / 4 = 1 / 4 := by
      field_simp
    rw [e, h5]
    field_simp
  · intro h
    have h2 : Real.sqrt 5 = 2 := by
      field_simp at h
      linarith
    have : (5 : ℝ) = 4 := by
      have := Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 5)
      rw [h2] at this
      linarith
    norm_num at this

end Joint

section GoldenJoint

/-- Along the golden family, for any sequence of levels `q_k = ℓ_k + 2` and any `Δτ_k → 0⁺`, test
vectors with bounded squared lengths whose canonical vacuum forms `g_kᵀ (½ ω_{q_k}⁻¹) g_k` converge
give the same limit for the reconstructed vacuum. -/
theorem golden_joint_limit_L2 {m2 : ℝ} (hm : 0 < m2) (lev : ℕ → ℕ) (g : ∀ k, Site (lev k) → ℝ)
    {dτ : ℕ → ℝ} (hd : ∀ k, 0 < dτ k) (hdt : Tendsto dτ atTop (𝓝 0)) {B : ℝ}
    (hB : ∀ k, g k ⬝ᵥ g k ≤ B) {c : ℝ}
    (hc : Tendsto (fun k => g k ⬝ᵥ
      (((1 / 2 : ℝ) • posDefInvSqrt (golden_whitened_posDef (lev k) hm)) *ᵥ g k)) atTop (𝓝 c)) :
    Tendsto (fun k => g k ⬝ᵥ (equalTimeCov (golden_whitened_posDef (lev k) hm) (dτ k) *ᵥ g k))
      atTop (𝓝 c) := by
  refine joint_limit_L2 (κ := fun k => Site (lev k))
    (A := fun k => whitenedPacket (goldenGrid (lev k + 2)) m2)
    (fun k => golden_whitened_posDef (lev k) hm) g hd hdt hB (hc.congr fun k => ?_)
  exact (congrArg (fun C => g k ⬝ᵥ (C *ᵥ g k))
    (half_invSqrt_eq (golden_whitened_posDef (lev k) hm))).symm

end GoldenJoint

#print axioms vacuum_charFun_tendsto
#print axioms response_kernel_tendsto
#print axioms golden_continuum_time_limit
#print axioms joint_limit
#print axioms joint_limit_bound_needed
#print axioms joint_limit_L2
#print axioms golden_joint_limit_L2

end

end OPH.SourceLaw
