import SourceLaw.TransferMatrixLift
import SourceLaw.ContinuumTimeLimit
import SourceLaw.CounterfamilyReader
import SourceLaw.StationaryWindow
import Mathlib.Topology.Instances.Matrix
import Mathlib.NumberTheory.Real.GoldenRatio

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Matrix-level slab-thickness limits of the reconstructed stationary chain

For a symmetric matrix `A`, a function of `A` is taken through its eigen-decomposition
(`Matrix.IsHermitian.cfc`). This file proves:

* `cfc_tendsto`: eigenvalue-wise convergence of a family of functions gives convergence of the
  matrix functions (entrywise topology). This is the matrix version, for `Matrix.IsHermitian.cfc`,
  of Mathlib's `tendsto_cfc_fun`.
* `cfc_quadratic_le`: an eigenvalue-wise bound `|f − g| ≤ h` gives the quadratic-form bound
  `|vᵀ f(A) v − vᵀ g(A) v| ≤ vᵀ h(A) v`.
* One mode, `Δτ → 0⁺`: the equal-time variance `c0` tends to `1/(2√λ)`; the momentum variance of
  the reconstructed vacuum, `p0 = sinh(Δτ Ω_Δτ)/(2Δτ) = √(𝒱/𝒦)/2`, tends to `√λ/2`; `p0 c0 = 1/4`;
  and `|c0 − 1/(2√λ)| ≤ Δτ² √λ / 16`. The constant `1/16` cannot be replaced by `1/32`
  (`quadratic_err_constant_needed`, at `λ = 1`, `Δτ = 1`), and at fixed `Δτ > 0` the variance
  `c0` tends to `+∞` as `λ → 0⁺` (`equalTimeVar_tendsto_atTop`), so no finite equal-time variance
  survives at a zero mode. Uniformly in `λ > 0`, `0 ≤ 1/(2√λ) − c0 ≤ Δτ/4`
  (`equalTimeVar_err_uniform`); the best uniform constant is `φ^{-5/2}/4`, `φ` the golden ratio
  (`equalTimeVar_err_sharp`), attained at `λ = 4φ/Δτ²` (`equalTimeVar_err_sharp_eq`).
* Matrix level, `A` positive definite, `Δτ → 0⁺`: `Ω_Δτ(A) → A^{1/2}`, `𝒦_Δτ(A) → 1`,
  `𝒱_Δτ(A) → A`, `C0(Δτ) → ½ A^{-1/2}`, `p0(A) → ½ A^{1/2}`, and the error bound
  `|vᵀ(C0(Δτ) − ½A^{-1/2})v| ≤ (Δτ²/16) vᵀA^{1/2}v`, with the spectrum-uniform companion
  `≤ (Δτ/4) vᵀv` that uses no energy (`equalTimeCov_err_L2`). The premise `A ≻ 0` is needed for the
  covariance by `equalTimeVar_tendsto_atTop`. The kinetic and potential coefficients extend
  continuously to a zero eigenvalue; there the Lean value of `𝒦_Δτ` is `0` by the convention
  `x / 0 = 0` while its continuous extension is `1`, so no statement here rests on that value.
* Original coordinates with positive capacities `M`: `M^{-1/2} C0(Δτ) M^{-1/2}` tends to
  `vacuumCov = ½ M^{-1/2} A^{-1/2} M^{-1/2}`; with `A = M^{-1/2} V M^{-1/2}` the matrix
  `vacuumPrecision = M^{1/2} A^{1/2} M^{1/2}` inverts `2 vacuumCov` and satisfies the ground
  relation `G M⁻¹ G = V` of the identification ladder; `vacuumPrecision_unique` shows it is the only
  positive definite solution (a positive definite and a positive semidefinite matrix with the same
  square are equal, `eq_of_mul_self_eq`), and `groundRel_two_solutions` shows that positivity is
  needed. That `vacuumCov` is the position covariance of the Gaussian ground state of
  `H = ½ Πᵀ M⁻¹ Π + ½ fᵀ V f` is the standard identification recorded in the companion note.
-/

namespace OPH.SourceLaw

open Filter Topology Matrix

noncomputable section

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {A : Matrix ι ι ℝ}

section Calculus

/-- Eigenvalue-wise convergence of a family of functions gives convergence of the matrix
functions. -/
theorem cfc_tendsto {α : Type*} {l : Filter α} (hA : A.IsHermitian) {f : α → ℝ → ℝ}
    {g : ℝ → ℝ} (h : ∀ i, Tendsto (fun t => f t (hA.eigenvalues i)) l (𝓝 (g (hA.eigenvalues i)))) :
    Tendsto (fun t => hA.cfc (f t)) l (𝓝 (hA.cfc g)) := by
  have hv : Tendsto (fun t => fun i => f t (hA.eigenvalues i)) l
      (𝓝 (fun i => g (hA.eigenvalues i))) := tendsto_pi_nhds.2 h
  have hc : Continuous (fun v : ι → ℝ => (hA.eigenvectorUnitary : Matrix ι ι ℝ) * diagonal v
      * star (hA.eigenvectorUnitary : Matrix ι ι ℝ)) :=
    (continuous_const.matrix_mul continuous_id.matrix_diagonal).matrix_mul continuous_const
  have key : (fun t => hA.cfc (f t)) = fun t => (hA.eigenvectorUnitary : Matrix ι ι ℝ)
      * diagonal (fun i => f t (hA.eigenvalues i))
      * star (hA.eigenvectorUnitary : Matrix ι ι ℝ) := by
    funext t
    rw [cfc_expand hA]
    rfl
  rw [key, cfc_expand hA]
  exact (hc.tendsto _).comp hv

/-- An eigenvalue-wise bound gives a quadratic-form bound. -/
theorem cfc_quadratic_le (hA : A.IsHermitian) {f g h : ℝ → ℝ}
    (hb : ∀ i, |f (hA.eigenvalues i) - g (hA.eigenvalues i)| ≤ h (hA.eigenvalues i)) (v : ι → ℝ) :
    |v ⬝ᵥ (hA.cfc f *ᵥ v) - v ⬝ᵥ (hA.cfc g *ᵥ v)| ≤ v ⬝ᵥ (hA.cfc h *ᵥ v) := by
  rw [dotProduct_cfc_mulVec hA, dotProduct_cfc_mulVec hA, dotProduct_cfc_mulVec hA,
    ← Finset.sum_sub_distrib]
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun i _ => ?_)
  set w := (v ᵥ* (hA.eigenvectorUnitary : Matrix ι ι ℝ)) i
  have hw : 0 ≤ w * w := mul_self_nonneg w
  have e : w * f (hA.eigenvalues i) * w - w * g (hA.eigenvalues i) * w
      = (w * w) * (f (hA.eigenvalues i) - g (hA.eigenvalues i)) := by ring
  rw [e, abs_mul, abs_of_nonneg hw]
  calc w * w * |f (hA.eigenvalues i) - g (hA.eigenvalues i)| ≤ w * w * h (hA.eigenvalues i) :=
        mul_le_mul_of_nonneg_left (hb i) hw
    _ = w * h (hA.eigenvalues i) * w := by ring

end Calculus

section Mode

variable {dτ lam : ℝ}

/-- The equal-time variance of a mode tends to the canonical vacuum variance `1/(2√λ)`. -/
theorem equalTimeVar_tendsto (hl : 0 < lam) :
    Tendsto (fun dτ => equalTimeVar dτ lam) (𝓝[>] 0) (𝓝 (1 / (2 * Real.sqrt lam))) := by
  have hs : 0 < Real.sqrt lam := Real.sqrt_pos.2 hl
  have h1 : Tendsto (fun dτ : ℝ => Real.sqrt (1 + dτ ^ 2 * lam / 4)) (𝓝[>] 0) (𝓝 1) := by
    have hc : Continuous (fun dτ : ℝ => Real.sqrt (1 + dτ ^ 2 * lam / 4)) := by fun_prop
    have := hc.tendsto 0
    simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, zero_mul, zero_div,
      add_zero, Real.sqrt_one] at this
    exact this.mono_left nhdsWithin_le_nhds
  have h2 : Tendsto (fun dτ : ℝ => 2 * Real.sqrt lam * Real.sqrt (1 + dτ ^ 2 * lam / 4)) (𝓝[>] 0)
      (𝓝 (2 * Real.sqrt lam * 1)) := tendsto_const_nhds.mul h1
  have h3 := (tendsto_const_nhds (x := (1 : ℝ))).div h2 (by positivity)
  rw [mul_one] at h3
  refine h3.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with dτ hdτ
  rw [equalTimeVar_closed hdτ hl]
  simp only [Pi.div_apply]

/-- Momentum variance of the reconstructed vacuum of one mode, `sinh(Δτ Ω_Δτ)/(2Δτ)`. -/
def momentumVar (dτ lam : ℝ) : ℝ := Real.sinh (dτ * latticeFreq dτ lam) / (2 * dτ)

theorem momentumVar_closed (hd : 0 < dτ) (hl : 0 < lam) :
    momentumVar dτ lam = Real.sqrt lam * Real.sqrt (1 + dτ ^ 2 * lam / 4) / 2 := by
  unfold momentumVar
  rw [sinh_latticeFreq hd]
  have hsq : (dτ * Real.sqrt lam / 2) ^ 2 = dτ ^ 2 * lam / 4 := by
    rw [div_pow, mul_pow, Real.sq_sqrt hl.le]
    ring
  rw [hsq]
  field_simp

/-- `√(𝒱/𝒦)/2` is the momentum variance of the reconstructed vacuum. -/
theorem reconstructed_momentum_var (hd : 0 < dτ) (hl : 0 < lam) :
    Real.sqrt (potentialCoeff dτ lam / kineticCoeff dτ lam) / 2 = momentumVar dτ lam := by
  have hW := latticeFreq_pos hd hl
  have hS := sinh_latticeFreq_pos hd hl
  have hratio : potentialCoeff dτ lam / kineticCoeff dτ lam
      = (Real.sinh (dτ * latticeFreq dτ lam) / dτ) ^ 2 := by
    unfold potentialCoeff kineticCoeff
    field_simp
  rw [hratio, Real.sqrt_sq (by positivity), momentumVar]
  field_simp

/-- The reconstructed vacuum of a mode has minimum uncertainty: `p0 c0 = 1/4`. -/
theorem momentumVar_mul_equalTimeVar (hd : 0 < dτ) (hl : 0 < lam) :
    momentumVar dτ lam * equalTimeVar dτ lam = 1 / 4 := by
  rw [momentumVar_closed hd hl, equalTimeVar_closed hd hl]
  have hs : 0 < Real.sqrt lam := Real.sqrt_pos.2 hl
  have hg : 0 < Real.sqrt (1 + dτ ^ 2 * lam / 4) := Real.sqrt_pos.2 (by positivity)
  field_simp
  norm_num

theorem momentumVar_tendsto (hl : 0 < lam) :
    Tendsto (fun dτ => momentumVar dτ lam) (𝓝[>] 0) (𝓝 (Real.sqrt lam / 2)) := by
  have h1 : Tendsto (fun dτ : ℝ => Real.sqrt (1 + dτ ^ 2 * lam / 4)) (𝓝[>] 0) (𝓝 1) := by
    have hc : Continuous (fun dτ : ℝ => Real.sqrt (1 + dτ ^ 2 * lam / 4)) := by fun_prop
    have := hc.tendsto 0
    simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, zero_mul, zero_div,
      add_zero, Real.sqrt_one] at this
    exact this.mono_left nhdsWithin_le_nhds
  have h2 := ((tendsto_const_nhds (x := Real.sqrt lam)).mul h1).div_const 2
  rw [mul_one] at h2
  refine h2.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with dτ hdτ
  rw [momentumVar_closed hdτ hl]

/-- One-mode error bound: `|c0(Δτ) − 1/(2√λ)| ≤ Δτ² √λ / 16`. -/
theorem equalTimeVar_err (hd : 0 < dτ) (hl : 0 < lam) :
    |equalTimeVar dτ lam - 1 / (2 * Real.sqrt lam)| ≤ dτ ^ 2 * Real.sqrt lam / 16 := by
  rw [equalTimeVar_closed hd hl]
  have hs : 0 < Real.sqrt lam := Real.sqrt_pos.2 hl
  set g := Real.sqrt (1 + dτ ^ 2 * lam / 4) with hg
  have hg1 : 1 ≤ g := by
    rw [hg, Real.one_le_sqrt]
    have : 0 ≤ dτ ^ 2 * lam / 4 := by positivity
    linarith
  have hgsq : g ^ 2 = 1 + dτ ^ 2 * lam / 4 := by
    rw [hg, Real.sq_sqrt (by positivity)]
  have hgpos : 0 < g := by linarith
  have hlam : Real.sqrt lam ^ 2 = lam := Real.sq_sqrt hl.le
  have hdiff : 1 / (2 * Real.sqrt lam * g) - 1 / (2 * Real.sqrt lam)
      = -((g - 1) / (2 * Real.sqrt lam * g)) := by
    field_simp
    ring
  rw [hdiff, abs_neg, abs_of_nonneg (div_nonneg (by linarith) (by positivity))]
  have step1 : (g - 1) / (2 * Real.sqrt lam * g) ≤ (g - 1) / (2 * Real.sqrt lam) := by
    apply div_le_div_of_nonneg_left (by linarith) (by positivity)
    nlinarith
  have step2 : (g - 1) / (2 * Real.sqrt lam) ≤ dτ ^ 2 * Real.sqrt lam / 16 := by
    rw [div_le_iff₀ (by positivity)]
    have hg' : g - 1 ≤ (g ^ 2 - 1) / 2 := by nlinarith
    calc g - 1 ≤ (g ^ 2 - 1) / 2 := hg'
      _ = dτ ^ 2 * lam / 8 := by rw [hgsq]; ring
      _ = dτ ^ 2 * Real.sqrt lam / 16 * (2 * Real.sqrt lam) := by
        rw [show dτ ^ 2 * Real.sqrt lam / 16 * (2 * Real.sqrt lam)
          = dτ ^ 2 * Real.sqrt lam ^ 2 / 8 by ring, hlam]
  exact step1.trans step2

/-- The constant `1/16` of `equalTimeVar_err` cannot be halved: at `λ = 1`, `Δτ = 1` the error is
`1/2 − 1/√5`, which exceeds `Δτ² √λ / 32`. -/
theorem quadratic_err_constant_needed :
    (1 : ℝ) ^ 2 * Real.sqrt 1 / 32 < |equalTimeVar 1 1 - 1 / (2 * Real.sqrt 1)| := by
  rw [equalTimeVar_closed one_pos one_pos, Real.sqrt_one]
  have hs : (16 / 15 : ℝ) < Real.sqrt (1 + 1 ^ 2 * 1 / 4) := by
    rw [Real.lt_sqrt (by norm_num)]
    norm_num
  have hlt : 1 / (2 * 1 * Real.sqrt (1 + 1 ^ 2 * 1 / 4)) < 15 / 32 := by
    rw [div_lt_iff₀ (by positivity)]
    linarith
  rw [abs_sub_comm, abs_of_pos (by linarith)]
  linarith

/-- At fixed `Δτ > 0` the equal-time variance tends to `+∞` as the eigenvalue tends to `0⁺`: no
finite stationary variance survives at a zero mode, so the premise `A ≻ 0` of the covariance limit
is needed. -/
theorem equalTimeVar_tendsto_atTop (hd : 0 < dτ) :
    Tendsto (fun lam => equalTimeVar dτ lam) (𝓝[>] 0) atTop := by
  have hc : Continuous fun lam : ℝ => 2 * Real.sqrt lam * Real.sqrt (1 + dτ ^ 2 * lam / 4) :=
    (continuous_const.mul Real.continuous_sqrt).mul
      (Continuous.sqrt (continuous_const.add ((continuous_const.mul continuous_id).div_const 4)))
  have hD : Tendsto (fun lam : ℝ => 2 * Real.sqrt lam * Real.sqrt (1 + dτ ^ 2 * lam / 4))
      (𝓝[>] 0) (𝓝[>] 0) := by
    refine tendsto_nhdsWithin_iff.2 ⟨?_, ?_⟩
    · have hv : 2 * Real.sqrt 0 * Real.sqrt (1 + dτ ^ 2 * 0 / 4) = (0 : ℝ) := by simp
      have h0 : Tendsto (fun lam : ℝ => 2 * Real.sqrt lam * Real.sqrt (1 + dτ ^ 2 * lam / 4)) (𝓝 0)
          (𝓝 (2 * Real.sqrt 0 * Real.sqrt (1 + dτ ^ 2 * 0 / 4))) := hc.tendsto 0
      rw [hv] at h0
      exact h0.mono_left nhdsWithin_le_nhds
    · filter_upwards [self_mem_nhdsWithin] with lam hl
      have hl' : (0 : ℝ) < lam := hl
      exact Set.mem_Ioi.2 (mul_pos (mul_pos two_pos (Real.sqrt_pos.2 hl'))
        (Real.sqrt_pos.2 (by positivity)))
  refine (tendsto_inv_nhdsGT_zero.comp hD).congr' ?_
  filter_upwards [self_mem_nhdsWithin] with lam hl
  have hl' : (0 : ℝ) < lam := hl
  show (2 * Real.sqrt lam * Real.sqrt (1 + dτ ^ 2 * lam / 4))⁻¹ = equalTimeVar dτ lam
  rw [equalTimeVar_closed hd hl', one_div]

/-- The gap `1/(2√λ) − c0` in closed form, with `s = √(1 + Δτ²λ/4)`. -/
theorem equalTimeVar_gap_eq (hd : 0 < dτ) (hl : 0 < lam) :
    1 / (2 * Real.sqrt lam) - equalTimeVar dτ lam
      = (Real.sqrt (1 + dτ ^ 2 * lam / 4) - 1)
        / (2 * Real.sqrt lam * Real.sqrt (1 + dτ ^ 2 * lam / 4)) := by
  rw [equalTimeVar_closed hd hl]
  have hs : 0 < Real.sqrt lam := Real.sqrt_pos.2 hl
  have hg : 0 < Real.sqrt (1 + dτ ^ 2 * lam / 4) := Real.sqrt_pos.2 (by positivity)
  field_simp

/-- The equal-time variance lies below `1/(2√λ)` by at most `Δτ/4`, uniformly in `λ > 0`. -/
theorem equalTimeVar_err_uniform (hd : 0 < dτ) (hl : 0 < lam) :
    0 ≤ 1 / (2 * Real.sqrt lam) - equalTimeVar dτ lam
      ∧ 1 / (2 * Real.sqrt lam) - equalTimeVar dτ lam ≤ dτ / 4 := by
  rw [equalTimeVar_gap_eq hd hl]
  have hs : 0 < Real.sqrt lam := Real.sqrt_pos.2 hl
  set g := Real.sqrt (1 + dτ ^ 2 * lam / 4) with hg
  have hg1 : 1 ≤ g := by
    rw [hg, Real.one_le_sqrt]
    have : 0 ≤ dτ ^ 2 * lam / 4 := by positivity
    linarith
  have hgsq : g ^ 2 = 1 + dτ ^ 2 * lam / 4 := by
    rw [hg, Real.sq_sqrt (by positivity)]
  have hlam : Real.sqrt lam ^ 2 = lam := Real.sq_sqrt hl.le
  have hgpos : 0 < g := by linarith
  refine ⟨div_nonneg (by linarith) (by positivity), ?_⟩
  rw [div_le_iff₀ (by positivity)]
  set u := dτ * Real.sqrt lam / 2 with hu
  have hu0 : 0 ≤ u := by positivity
  have hu2 : u ^ 2 = g ^ 2 - 1 := by
    rw [hu, hgsq]
    have : (dτ * Real.sqrt lam / 2) ^ 2 = dτ ^ 2 * Real.sqrt lam ^ 2 / 4 := by ring
    rw [this, hlam]
    ring
  have key : g - 1 ≤ u * g := by
    nlinarith [mul_nonneg (sub_nonneg.2 hg1) (sub_nonneg.2 hg1), mul_nonneg hu0 hgpos.le,
      sq_nonneg (u * g - (g - 1)), mul_nonneg (sub_nonneg.2 hg1) hgpos.le]
  have : dτ / 4 * (2 * Real.sqrt lam * g) = u * g := by rw [hu]; ring
  rw [this]
  exact key

/-- The best uniform constant: `1/(2√λ) − c0 ≤ Δτ φ^{-5/2}/4` for every `λ > 0`, `φ` the golden
ratio. -/
theorem equalTimeVar_err_sharp (hd : 0 < dτ) (hl : 0 < lam) :
    1 / (2 * Real.sqrt lam) - equalTimeVar dτ lam
      ≤ dτ / (4 * Real.goldenRatio ^ 2 * Real.sqrt Real.goldenRatio) := by
  rw [equalTimeVar_gap_eq hd hl]
  have hs : 0 < Real.sqrt lam := Real.sqrt_pos.2 hl
  set g := Real.sqrt (1 + dτ ^ 2 * lam / 4) with hg
  have hg1 : 1 ≤ g := by
    rw [hg, Real.one_le_sqrt]
    have : 0 ≤ dτ ^ 2 * lam / 4 := by positivity
    linarith
  have hgsq : g ^ 2 = 1 + dτ ^ 2 * lam / 4 := by
    rw [hg, Real.sq_sqrt (by positivity)]
  have hlam : Real.sqrt lam ^ 2 = lam := Real.sq_sqrt hl.le
  have hgpos : 0 < g := by linarith
  have hφ : 0 < Real.goldenRatio := Real.goldenRatio_pos
  have hφ2 : Real.goldenRatio ^ 2 = Real.goldenRatio + 1 := Real.goldenRatio_sq
  set r := Real.sqrt Real.goldenRatio with hr
  have hr0 : 0 < r := Real.sqrt_pos.2 hφ
  have hr2 : r ^ 2 = Real.goldenRatio := Real.sq_sqrt hφ.le
  set u := dτ * Real.sqrt lam / 2 with hu
  have hu0 : 0 ≤ u := by positivity
  have hu2 : u ^ 2 = g ^ 2 - 1 := by
    rw [hu, hgsq]
    have : (dτ * Real.sqrt lam / 2) ^ 2 = dτ ^ 2 * Real.sqrt lam ^ 2 / 4 := by ring
    rw [this, hlam]
    ring
  -- (g + 1) g² − (g − 1) φ⁵ = (g − φ)² (g + 1 + 2φ) ≥ 0, using φ² = φ + 1
  have hpoly : 0 ≤ (g + 1) * g ^ 2 - (g - 1) * Real.goldenRatio ^ 5 := by
    have e : (g + 1) * g ^ 2 - (g - 1) * Real.goldenRatio ^ 5
        = (g - Real.goldenRatio) ^ 2 * (g + 1 + 2 * Real.goldenRatio) := by
      have h5 : Real.goldenRatio ^ 5 = 5 * Real.goldenRatio + 3 := by
        have : Real.goldenRatio ^ 5 = (Real.goldenRatio ^ 2) ^ 2 * Real.goldenRatio := by ring
        rw [this, hφ2]
        nlinarith [hφ2]
      rw [h5]
      nlinarith [hφ2]
    rw [e]
    positivity
  -- squared form: ((g − 1) φ² r)² ≤ (u g)²
  have hsq : ((g - 1) * Real.goldenRatio ^ 2 * r) ^ 2 ≤ (u * g) ^ 2 := by
    have l : ((g - 1) * Real.goldenRatio ^ 2 * r) ^ 2 = (g - 1) ^ 2 * Real.goldenRatio ^ 5 := by
      have : ((g - 1) * Real.goldenRatio ^ 2 * r) ^ 2
          = (g - 1) ^ 2 * Real.goldenRatio ^ 4 * r ^ 2 := by ring
      rw [this, hr2]
      ring
    have rr : (u * g) ^ 2 = (g - 1) * ((g + 1) * g ^ 2) := by
      have : (u * g) ^ 2 = u ^ 2 * g ^ 2 := by ring
      rw [this, hu2]
      ring
    rw [l, rr]
    have hg0 : 0 ≤ g - 1 := by linarith
    nlinarith [mul_nonneg hg0 hpoly]
  have hlin : (g - 1) * Real.goldenRatio ^ 2 * r ≤ u * g := by
    have ha : 0 ≤ (g - 1) * Real.goldenRatio ^ 2 * r := by
      have hg0 : 0 ≤ g - 1 := by linarith
      positivity
    have hb : 0 ≤ u * g := by positivity
    nlinarith [sq_nonneg ((g - 1) * Real.goldenRatio ^ 2 * r - u * g), hsq, ha, hb]
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  have : dτ * (2 * Real.sqrt lam * g) = 4 * (u * g) := by rw [hu]; ring
  rw [this]
  nlinarith [hlin]

/-- The best uniform constant is attained at `λ = 4φ/Δτ²`. -/
theorem equalTimeVar_err_sharp_eq (hd : 0 < dτ) :
    1 / (2 * Real.sqrt (4 * Real.goldenRatio / dτ ^ 2))
        - equalTimeVar dτ (4 * Real.goldenRatio / dτ ^ 2)
      = dτ / (4 * Real.goldenRatio ^ 2 * Real.sqrt Real.goldenRatio) := by
  have hφ : 0 < Real.goldenRatio := Real.goldenRatio_pos
  have hφ2 : Real.goldenRatio ^ 2 = Real.goldenRatio + 1 := Real.goldenRatio_sq
  have hl : 0 < 4 * Real.goldenRatio / dτ ^ 2 := by positivity
  rw [equalTimeVar_gap_eq hd hl]
  generalize Real.goldenRatio = φ at *
  have hsl : Real.sqrt (4 * φ / dτ ^ 2) = 2 * Real.sqrt φ / dτ := by
    rw [show 4 * φ / dτ ^ 2 = (2 * Real.sqrt φ / dτ) ^ 2 by
      rw [div_pow, mul_pow, Real.sq_sqrt hφ.le]; ring]
    exact Real.sqrt_sq (by positivity)
  have h4 : dτ ^ 2 * (4 * φ / dτ ^ 2) = 4 * φ := by
    have hd2 : dτ ^ 2 ≠ 0 := pow_ne_zero 2 hd.ne'
    field_simp
  have hg : Real.sqrt (1 + dτ ^ 2 * (4 * φ / dτ ^ 2) / 4) = φ := by
    rw [h4, show 1 + 4 * φ / 4 = φ ^ 2 by rw [hφ2]; ring]
    exact Real.sqrt_sq hφ.le
  rw [hsl, hg]
  have hr0 : 0 < Real.sqrt φ := Real.sqrt_pos.2 hφ
  have hφ1 : φ - 1 = 1 / φ := by
    field_simp
    linear_combination hφ2
  rw [hφ1]
  have hφ0 : φ ≠ 0 := hφ.ne'
  have hr0' : Real.sqrt φ ≠ 0 := hr0.ne'
  have hd0 : dτ ≠ 0 := hd.ne'
  field_simp
  ring

end Mode

section Lift

variable (hA : A.PosDef)

theorem latticeFreq_cfc_tendsto :
    Tendsto (fun dτ => hA.1.cfc (latticeFreq dτ)) (𝓝[>] 0) (𝓝 (hA.1.cfc Real.sqrt)) :=
  cfc_tendsto hA.1 (f := fun dτ => latticeFreq dτ) fun i => latticeFreq_tendsto (eig_pos hA i)

theorem kineticCoeff_cfc_tendsto :
    Tendsto (fun dτ => hA.1.cfc (kineticCoeff dτ)) (𝓝[>] 0) (𝓝 1) := by
  rw [← cfc_one hA.1]
  exact cfc_tendsto hA.1 (f := fun dτ => kineticCoeff dτ) fun i =>
    kineticCoeff_tendsto (eig_pos hA i)

theorem potentialCoeff_cfc_tendsto :
    Tendsto (fun dτ => hA.1.cfc (potentialCoeff dτ)) (𝓝[>] 0) (𝓝 A) := by
  conv => arg 3; rw [← cfc_id hA.1]
  exact cfc_tendsto hA.1 (f := fun dτ => potentialCoeff dτ) fun i =>
    potentialCoeff_tendsto (eig_pos hA i)

/-- The equal-time covariance of the stationary chain tends to `½ A^{-1/2}`. -/
theorem equalTimeCov_tendsto :
    Tendsto (fun dτ => equalTimeCov hA dτ) (𝓝[>] 0)
      (𝓝 (hA.1.cfc (fun x => 1 / (2 * Real.sqrt x)))) :=
  cfc_tendsto hA.1 (f := fun dτ => equalTimeVar dτ) fun i => equalTimeVar_tendsto (eig_pos hA i)

/-- The momentum covariance of the reconstructed vacuum tends to `½ A^{1/2}`. -/
theorem momentumCov_tendsto :
    Tendsto (fun dτ => hA.1.cfc (momentumVar dτ)) (𝓝[>] 0)
      (𝓝 (hA.1.cfc (fun x => Real.sqrt x / 2))) :=
  cfc_tendsto hA.1 (f := fun dτ => momentumVar dτ) fun i => momentumVar_tendsto (eig_pos hA i)

/-- Quadratic-form error of the equal-time covariance. -/
theorem equalTimeCov_quadratic_err {dτ : ℝ} (hd : 0 < dτ) (v : ι → ℝ) :
    |v ⬝ᵥ (equalTimeCov hA dτ *ᵥ v) - v ⬝ᵥ (hA.1.cfc (fun x => 1 / (2 * Real.sqrt x)) *ᵥ v)|
      ≤ dτ ^ 2 / 16 * (v ⬝ᵥ (hA.1.cfc Real.sqrt *ᵥ v)) := by
  have hb := cfc_quadratic_le hA.1 (f := equalTimeVar dτ) (g := fun x => 1 / (2 * Real.sqrt x))
    (h := fun x => dτ ^ 2 * Real.sqrt x / 16) (fun i => equalTimeVar_err hd (eig_pos hA i)) v
  have hc : hA.1.cfc (fun x => dτ ^ 2 * Real.sqrt x / 16) = (dτ ^ 2 / 16) • hA.1.cfc Real.sqrt := by
    rw [cfc_smul hA.1]
    exact cfc_congr hA.1 fun i => by ring
  rw [hc, smul_mulVec, dotProduct_smul, smul_eq_mul] at hb
  exact hb

/-- The quadratic-form error bounded by the plain length, uniformly in the spectrum:
`|vᵀ(C0(Δτ) − ½A^{-1/2})v| ≤ (Δτ/4) vᵀv`. -/
theorem equalTimeCov_err_L2 {dτ : ℝ} (hd : 0 < dτ) (v : ι → ℝ) :
    |v ⬝ᵥ (equalTimeCov hA dτ *ᵥ v) - v ⬝ᵥ (hA.1.cfc (fun x => 1 / (2 * Real.sqrt x)) *ᵥ v)|
      ≤ dτ / 4 * (v ⬝ᵥ v) := by
  have hb := cfc_quadratic_le hA.1 (f := equalTimeVar dτ) (g := fun x => 1 / (2 * Real.sqrt x))
    (h := fun _ => dτ / 4) (fun i => by
      have h := equalTimeVar_err_uniform hd (eig_pos hA i)
      rw [abs_sub_comm, abs_of_nonneg h.1]
      exact h.2) v
  have hc : hA.1.cfc (fun _ => dτ / 4) = (dτ / 4) • (1 : Matrix ι ι ℝ) := cfc_const hA.1 _
  rw [hc, smul_mulVec, one_mulVec, dotProduct_smul, smul_eq_mul] at hb
  exact hb

end Lift

section Original

variable {M : ι → ℝ}

/-- Original-coordinate vacuum position covariance `½ M^{-1/2} A^{-1/2} M^{-1/2}`. -/
def vacuumCov (M : ι → ℝ) (hA : A.PosDef) : Matrix ι ι ℝ :=
  (1 / 2 : ℝ) • (jInvSqrt M * posDefInvSqrt hA * jInvSqrt M)

/-- `M^{1/2} A^{1/2} M^{1/2}`. -/
def vacuumPrecision (M : ι → ℝ) (hA : A.PosDef) : Matrix ι ι ℝ :=
  jSqrt M * posDefSqrt hA * jSqrt M

theorem half_invSqrt_eq (hA : A.PosDef) :
    hA.1.cfc (fun x => 1 / (2 * Real.sqrt x)) = (1 / 2 : ℝ) • posDefInvSqrt hA := by
  rw [posDefInvSqrt, cfc_smul hA.1]
  exact cfc_congr hA.1 fun i => by field_simp

/-- In original coordinates the equal-time covariance tends to `vacuumCov`. -/
theorem origCov_tendsto (hA : A.PosDef) :
    Tendsto (fun dτ => jInvSqrt M * equalTimeCov hA dτ * jInvSqrt M) (𝓝[>] 0)
      (𝓝 (vacuumCov M hA)) := by
  have h := ((tendsto_const_nhds (x := jInvSqrt M)).mul (equalTimeCov_tendsto hA)).mul
    (tendsto_const_nhds (x := jInvSqrt M))
  rw [half_invSqrt_eq hA, Matrix.mul_smul, Matrix.smul_mul] at h
  exact h

theorem vacuumCov_mul_vacuumPrecision (hM : ∀ i, 0 < M i) (hA : A.PosDef) :
    (2 : ℝ) • vacuumCov M hA * vacuumPrecision M hA = 1 := by
  have h1 : jInvSqrt M * jSqrt M = 1 := jInvSqrt_mul_jSqrt hM
  have h2 : posDefInvSqrt hA * posDefSqrt hA = 1 := posDefInvSqrt_mul_sqrt hA
  unfold vacuumCov vacuumPrecision
  rw [smul_smul, show (2 : ℝ) * (1 / 2) = 1 by norm_num, one_smul]
  calc jInvSqrt M * posDefInvSqrt hA * jInvSqrt M * (jSqrt M * posDefSqrt hA * jSqrt M)
      = jInvSqrt M * (posDefInvSqrt hA * (jInvSqrt M * jSqrt M) * posDefSqrt hA) * jSqrt M := by
        simp only [Matrix.mul_assoc]
    _ = 1 := by rw [h1, Matrix.mul_one, h2, Matrix.mul_one, h1]

theorem jInvSqrt_mul_self (hM : ∀ i, 0 < M i) :
    jInvSqrt M * jInvSqrt M = diagonal (fun i => (M i)⁻¹) := by
  rw [jInvSqrt, diagonal_mul_diagonal]
  congr 1
  funext i
  rw [← mul_inv, Real.mul_self_sqrt (hM i).le]

theorem jSqrt_mul_diagInv_mul_jSqrt (hM : ∀ i, 0 < M i) :
    jSqrt M * diagonal (fun i => (M i)⁻¹) * jSqrt M = 1 := by
  rw [jSqrt, diagonal_mul_diagonal, diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext i
  have hsq : Real.sqrt (M i) * Real.sqrt (M i) = M i := Real.mul_self_sqrt (hM i).le
  calc Real.sqrt (M i) * (M i)⁻¹ * Real.sqrt (M i)
      = (Real.sqrt (M i) * Real.sqrt (M i)) * (M i)⁻¹ := by ring
    _ = 1 := by rw [hsq, mul_inv_cancel₀ (hM i).ne']

/-- Ground relation: with `A = M^{-1/2} V M^{-1/2}`, `G M⁻¹ G = V` for `G = vacuumPrecision`. -/
theorem vacuumPrecision_groundRel {V : Matrix ι ι ℝ} (hM : ∀ i, 0 < M i)
    (hA : (jInvSqrt M * V * jInvSqrt M).PosDef) :
    vacuumPrecision M hA * diagonal (fun i => (M i)⁻¹) * vacuumPrecision M hA = V := by
  have hJ := jSqrt_mul_diagInv_mul_jSqrt hM
  have hS : posDefSqrt hA * posDefSqrt hA = jInvSqrt M * V * jInvSqrt M := posDefSqrt_mul_self hA
  have h1 : jSqrt M * jInvSqrt M = 1 := jSqrt_mul_jInvSqrt hM
  unfold vacuumPrecision
  calc jSqrt M * posDefSqrt hA * jSqrt M * diagonal (fun i => (M i)⁻¹)
        * (jSqrt M * posDefSqrt hA * jSqrt M)
      = jSqrt M * posDefSqrt hA * (jSqrt M * diagonal (fun i => (M i)⁻¹) * jSqrt M)
          * posDefSqrt hA * jSqrt M := by simp only [Matrix.mul_assoc]
    _ = jSqrt M * (posDefSqrt hA * posDefSqrt hA) * jSqrt M := by
          rw [hJ, Matrix.mul_one]; simp only [Matrix.mul_assoc]
    _ = V := by
          rw [hS]
          calc jSqrt M * (jInvSqrt M * V * jInvSqrt M) * jSqrt M
              = (jSqrt M * jInvSqrt M) * V * (jInvSqrt M * jSqrt M) := by
                simp only [Matrix.mul_assoc]
            _ = V := by rw [h1, jInvSqrt_mul_jSqrt hM, Matrix.one_mul, Matrix.mul_one]

/-- The spectral square root of a positive definite matrix is positive semidefinite. -/
theorem posDefSqrt_posSemidef {W : Matrix ι ι ℝ} (hW : W.PosDef) : (posDefSqrt hW).PosSemidef := by
  unfold posDefSqrt
  rw [cfc_expand hW.1]
  have hD : (diagonal (Real.sqrt ∘ hW.1.eigenvalues)).PosSemidef :=
    Matrix.PosSemidef.diagonal fun i => Real.sqrt_nonneg _
  have h := hD.mul_mul_conjTranspose_same (hW.1.eigenvectorUnitary : Matrix ι ι ℝ)
  rwa [← Matrix.star_eq_conjTranspose] at h

omit [DecidableEq ι] in
/-- A positive definite matrix and a positive semidefinite matrix with the same square are equal
(the difference solves a Sylvester equation with positive definite coefficient). -/
theorem eq_of_mul_self_eq {H K : Matrix ι ι ℝ} (hH : H.PosDef) (hK : K.PosSemidef)
    (h : H * H = K * K) : H = K := by
  have hS : (H + K).PosDef := hH.add_posSemidef hK
  have hz : (H + K) * (H - K) + (H - K) * (H + K) = 0 := by
    have e : (H + K) * (H - K) + (H - K) * (H + K) = (H * H - K * K) + (H * H - K * K) := by
      simp only [Matrix.add_mul, Matrix.mul_add, Matrix.sub_mul, Matrix.mul_sub]
      abel
    rw [e, h, sub_self, add_zero]
  exact sub_eq_zero.1 (sylvester_zero _ _ hS hz)

/-- `vacuumPrecision` is the only positive definite solution of the ground relation
`G M⁻¹ G = V`. -/
theorem vacuumPrecision_unique {V : Matrix ι ι ℝ} (hM : ∀ i, 0 < M i)
    (hA : (jInvSqrt M * V * jInvSqrt M).PosDef) {G : Matrix ι ι ℝ} (hG : G.PosDef)
    (h : G * diagonal (fun i => (M i)⁻¹) * G = V) : G = vacuumPrecision M hA := by
  have hH : (jInvSqrt M * G * jInvSqrt M).PosDef := jWhitened_posDef hM hG
  have hsq : (jInvSqrt M * G * jInvSqrt M) * (jInvSqrt M * G * jInvSqrt M)
      = posDefSqrt hA * posDefSqrt hA := by
    rw [posDefSqrt_mul_self hA, ← h, ← jInvSqrt_mul_self hM]
    simp only [Matrix.mul_assoc]
  have e := eq_of_mul_self_eq hH (posDefSqrt_posSemidef hA) hsq
  unfold vacuumPrecision
  rw [← e]
  calc G = (jSqrt M * jInvSqrt M) * G * (jInvSqrt M * jSqrt M) := by
        rw [jSqrt_mul_jInvSqrt hM, jInvSqrt_mul_jSqrt hM, Matrix.one_mul, Matrix.mul_one]
    _ = jSqrt M * (jInvSqrt M * G * jInvSqrt M) * jSqrt M := by simp only [Matrix.mul_assoc]

end Original

/-- The positivity premise of `vacuumPrecision_unique` is needed: with unit capacities and
`V = diag(1, 4)` both `diag(1, 2)` and `diag(1, −2)` solve the ground relation. -/
theorem groundRel_two_solutions :
    (!![1, 0; 0, 2] : Matrix (Fin 2) (Fin 2) ℝ) * diagonal (fun _ => (1 : ℝ)⁻¹) * !![1, 0; 0, 2]
        = !![1, 0; 0, 4]
      ∧ (!![1, 0; 0, -2] : Matrix (Fin 2) (Fin 2) ℝ) * diagonal (fun _ => (1 : ℝ)⁻¹)
          * !![1, 0; 0, -2]
        = !![1, 0; 0, 4]
      ∧ (!![1, 0; 0, 2] : Matrix (Fin 2) (Fin 2) ℝ) ≠ !![1, 0; 0, -2] := by
  refine ⟨?_, ?_, ?_⟩
  · ext i j
    fin_cases i <;> fin_cases j <;> simp [Matrix.mul_apply, Fin.sum_univ_two]
    all_goals norm_num
  · ext i j
    fin_cases i <;> fin_cases j <;> simp [Matrix.mul_apply, Fin.sum_univ_two]
    all_goals norm_num
  · intro h
    have := congrFun (congrFun h 1) 1
    norm_num at this

#print axioms cfc_tendsto
#print axioms cfc_quadratic_le
#print axioms equalTimeVar_tendsto
#print axioms momentumVar_closed
#print axioms reconstructed_momentum_var
#print axioms momentumVar_mul_equalTimeVar
#print axioms momentumVar_tendsto
#print axioms equalTimeVar_err
#print axioms quadratic_err_constant_needed
#print axioms equalTimeVar_tendsto_atTop
#print axioms equalTimeVar_gap_eq
#print axioms equalTimeVar_err_uniform
#print axioms equalTimeVar_err_sharp
#print axioms equalTimeVar_err_sharp_eq
#print axioms latticeFreq_cfc_tendsto
#print axioms kineticCoeff_cfc_tendsto
#print axioms potentialCoeff_cfc_tendsto
#print axioms equalTimeCov_tendsto
#print axioms momentumCov_tendsto
#print axioms equalTimeCov_quadratic_err
#print axioms equalTimeCov_err_L2
#print axioms half_invSqrt_eq
#print axioms origCov_tendsto
#print axioms vacuumCov_mul_vacuumPrecision
#print axioms jInvSqrt_mul_self
#print axioms jSqrt_mul_diagInv_mul_jSqrt
#print axioms vacuumPrecision_groundRel
#print axioms posDefSqrt_posSemidef
#print axioms eq_of_mul_self_eq
#print axioms vacuumPrecision_unique
#print axioms groundRel_two_solutions

end

end OPH.SourceLaw
