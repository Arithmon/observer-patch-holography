import SourceLaw.IdentificationLadder
import Mathlib.Analysis.SpecialFunctions.Arsinh
import Mathlib.Analysis.SpecialFunctions.Sqrt

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Transfer structure of the stationary chain of a register mode, one mode at a time

In whitened coordinates the interior rows of the register precision have slab blocks
`2/Δτ + Δτλ` and time-neighbour blocks `−1/Δτ` on each eigenmode `λ > 0` of the whitened
stiffness (the two end slabs differ). The stationary chain of that mode is the Gaussian
autoregression with one-step transfer `r` and equal-time variance `c0`. This file defines `r`,
the lattice frequency `Ω_Δτ` and `c0` in closed form and proves, for `Δτ > 0` and `λ > 0`:

* `transferEig_pos`, `transferEig_lt_one`, `transferEig_add_inv`: `0 < r < 1`,
  `r + r⁻¹ = 2 + Δτ²λ`.
* `exp_neg_latticeFreq`, `cosh_latticeFreq`: `r = exp(−Δτ Ω_Δτ)` and
  `cosh(Δτ Ω_Δτ) = 1 + Δτ²λ/2` with `Ω_Δτ = (2/Δτ) arsinh(Δτ√λ/2)`.
* `equalTimeVar_closed`: `c0 = 1 / (2 √λ √(1 + Δτ²λ/4))`.
* `chain_innovation`, `chain_precision_diag`, `chain_precision_off`: the autoregression with
  transfer `r` and innovation `Δτ r` has exactly the interior precision blocks of the register
  mode.
* `markov_reader_unit`, `contact_identity_mode`: the two-time data return the unit time
  coupling (in original coordinates the time coefficient `J`, see `history_reader_general`), and
  the lag-zero row of the interior precision holds.
-/

namespace OPH.SourceLaw

open Real

/-- The slab coefficient `a = 1 + Δτ²λ/2`. -/
noncomputable def slabCoeff (dτ lam : ℝ) : ℝ := 1 + dτ ^ 2 * lam / 2

/-- One-step transfer of the stationary mode: the root of `r + r⁻¹ = 2a` below one. -/
noncomputable def transferEig (dτ lam : ℝ) : ℝ :=
  slabCoeff dτ lam - Real.sqrt (slabCoeff dτ lam ^ 2 - 1)

/-- Lattice frequency `Ω_Δτ = (2/Δτ) arsinh(Δτ √λ / 2)`. -/
noncomputable def latticeFreq (dτ lam : ℝ) : ℝ := 2 / dτ * Real.arsinh (dτ * Real.sqrt lam / 2)

/-- Equal-time variance of the stationary mode, `c0 = Δτ r / (1 − r²)`. -/
noncomputable def equalTimeVar (dτ lam : ℝ) : ℝ :=
  dτ * transferEig dτ lam / (1 - transferEig dτ lam ^ 2)

section Mode

variable {dτ lam : ℝ}

theorem slabCoeff_gt_one (hd : 0 < dτ) (hl : 0 < lam) : 1 < slabCoeff dτ lam := by
  unfold slabCoeff
  have : 0 < dτ ^ 2 * lam := by positivity
  linarith

theorem slabCoeff_sq_sub_one_pos (hd : 0 < dτ) (hl : 0 < lam) : 0 < slabCoeff dτ lam ^ 2 - 1 := by
  have h := slabCoeff_gt_one hd hl
  nlinarith

theorem sqrt_slab_sq (hd : 0 < dτ) (hl : 0 < lam) :
    Real.sqrt (slabCoeff dτ lam ^ 2 - 1) ^ 2 = slabCoeff dτ lam ^ 2 - 1 :=
  Real.sq_sqrt (slabCoeff_sq_sub_one_pos hd hl).le

theorem transferEig_pos (hd : 0 < dτ) (hl : 0 < lam) : 0 < transferEig dτ lam := by
  unfold transferEig
  have ha := slabCoeff_gt_one hd hl
  have hs : Real.sqrt (slabCoeff dτ lam ^ 2 - 1) < slabCoeff dτ lam := by
    rw [Real.sqrt_lt' (by linarith)]
    linarith
  linarith

theorem transferEig_lt_one (hd : 0 < dτ) (hl : 0 < lam) : transferEig dτ lam < 1 := by
  unfold transferEig
  have ha := slabCoeff_gt_one hd hl
  have hs : slabCoeff dτ lam - 1 < Real.sqrt (slabCoeff dτ lam ^ 2 - 1) := by
    rw [Real.lt_sqrt (by linarith)]
    nlinarith
  linarith

/-- `(a − s)(a + s) = 1` with `s = √(a² − 1)`. -/
theorem transferEig_mul_conj (hd : 0 < dτ) (hl : 0 < lam) :
    transferEig dτ lam * (slabCoeff dτ lam + Real.sqrt (slabCoeff dτ lam ^ 2 - 1)) = 1 := by
  unfold transferEig
  have h := sqrt_slab_sq hd hl
  nlinarith

theorem transferEig_inv (hd : 0 < dτ) (hl : 0 < lam) :
    (transferEig dτ lam)⁻¹ = slabCoeff dτ lam + Real.sqrt (slabCoeff dτ lam ^ 2 - 1) :=
  inv_eq_of_mul_eq_one_right (transferEig_mul_conj hd hl)

theorem transferEig_add_inv (hd : 0 < dτ) (hl : 0 < lam) :
    transferEig dτ lam + (transferEig dτ lam)⁻¹ = 2 + dτ ^ 2 * lam := by
  rw [transferEig_inv hd hl]
  unfold transferEig slabCoeff
  ring

/-- The transfer eigenvalue is the exponential of minus the lattice frequency times `Δτ`. -/
theorem exp_neg_latticeFreq (hd : 0 < dτ) (hl : 0 < lam) :
    Real.exp (-(dτ * latticeFreq dτ lam)) = transferEig dτ lam := by
  set x := dτ * Real.sqrt lam / 2 with hx
  have hx0 : 0 ≤ x := by positivity
  have hxsq : x ^ 2 = dτ ^ 2 * lam / 4 := by
    rw [hx, div_pow, mul_pow, Real.sq_sqrt hl.le]
    ring
  have harg : -(dτ * latticeFreq dτ lam) = -(2 * Real.arsinh x) := by
    unfold latticeFreq
    rw [← hx]
    field_simp
  rw [harg, show -(2 * Real.arsinh x) = -(Real.arsinh x) + -(Real.arsinh x) by ring,
    Real.exp_add, Real.exp_neg, Real.exp_arsinh]
  -- `(x + √(1 + x²))⁻² = a − √(a² − 1)` with `a = 1 + 2x²`
  have hq : 0 < x + Real.sqrt (1 + x ^ 2) := by positivity
  have hs1 : Real.sqrt (1 + x ^ 2) ^ 2 = 1 + x ^ 2 := Real.sq_sqrt (by positivity)
  have ha : slabCoeff dτ lam = 1 + 2 * x ^ 2 := by
    unfold slabCoeff
    rw [hxsq]
    ring
  have hsa : Real.sqrt (slabCoeff dτ lam ^ 2 - 1) = 2 * x * Real.sqrt (1 + x ^ 2) := by
    rw [ha]
    rw [show (1 + 2 * x ^ 2) ^ 2 - 1 = (2 * x * Real.sqrt (1 + x ^ 2)) ^ 2 by
      rw [mul_pow, hs1]; ring]
    exact Real.sqrt_sq (by positivity)
  unfold transferEig
  rw [hsa, ha, ← mul_inv]
  apply inv_eq_of_mul_eq_one_right
  linear_combination (1 - 2 * x ^ 2 - 2 * x * Real.sqrt (1 + x ^ 2)) * hs1

/-- Dispersion: `cosh(Δτ Ω_Δτ) = 1 + Δτ²λ/2`. -/
theorem cosh_latticeFreq (hd : 0 < dτ) (hl : 0 < lam) :
    Real.cosh (dτ * latticeFreq dτ lam) = 1 + dτ ^ 2 * lam / 2 := by
  have harg : dτ * latticeFreq dτ lam = 2 * Real.arsinh (dτ * Real.sqrt lam / 2) := by
    unfold latticeFreq
    field_simp
  rw [harg, Real.cosh_two_mul, Real.cosh_sq', Real.sinh_arsinh]
  rw [div_pow, mul_pow, Real.sq_sqrt hl.le]
  ring

theorem one_sub_transferEig_sq_pos (hd : 0 < dτ) (hl : 0 < lam) :
    0 < 1 - transferEig dτ lam ^ 2 := by
  have h0 := transferEig_pos hd hl
  have h1 := transferEig_lt_one hd hl
  nlinarith

theorem equalTimeVar_pos (hd : 0 < dτ) (hl : 0 < lam) : 0 < equalTimeVar dτ lam := by
  unfold equalTimeVar
  have := transferEig_pos hd hl
  have := one_sub_transferEig_sq_pos hd hl
  positivity

/-- `1 − r² = 2 r √(a² − 1)`. -/
theorem one_sub_transferEig_sq (hd : 0 < dτ) (hl : 0 < lam) :
    1 - transferEig dτ lam ^ 2 = 2 * transferEig dτ lam * Real.sqrt (slabCoeff dτ lam ^ 2 - 1) := by
  have h := transferEig_mul_conj hd hl
  unfold transferEig at *
  nlinarith [sqrt_slab_sq hd hl]

/-- Closed form of the equal-time variance. -/
theorem equalTimeVar_closed (hd : 0 < dτ) (hl : 0 < lam) :
    equalTimeVar dτ lam = 1 / (2 * Real.sqrt lam * Real.sqrt (1 + dτ ^ 2 * lam / 4)) := by
  have hr := transferEig_pos hd hl
  have hsa : Real.sqrt (slabCoeff dτ lam ^ 2 - 1)
      = dτ * Real.sqrt lam * Real.sqrt (1 + dτ ^ 2 * lam / 4) := by
    rw [show slabCoeff dτ lam ^ 2 - 1
        = (dτ * Real.sqrt lam * Real.sqrt (1 + dτ ^ 2 * lam / 4)) ^ 2 by
      rw [mul_pow, mul_pow, Real.sq_sqrt hl.le, Real.sq_sqrt (by positivity)]
      unfold slabCoeff
      ring]
    exact Real.sqrt_sq (by positivity)
  unfold equalTimeVar
  rw [one_sub_transferEig_sq hd hl, hsa]
  have : 0 < Real.sqrt lam := Real.sqrt_pos.2 hl
  have : 0 < Real.sqrt (1 + dτ ^ 2 * lam / 4) := Real.sqrt_pos.2 (by positivity)
  field_simp

/-- The innovation variance of the autoregression is `Δτ r`. -/
theorem chain_innovation (hd : 0 < dτ) (hl : 0 < lam) :
    equalTimeVar dτ lam * (1 - transferEig dτ lam ^ 2) = dτ * transferEig dτ lam := by
  unfold equalTimeVar
  have := one_sub_transferEig_sq_pos hd hl
  field_simp

/-- Slab block of the chain precision: `q⁻¹ + r q⁻¹ r = 2/Δτ + Δτλ` with `q = Δτ r`. -/
theorem chain_precision_diag (hd : 0 < dτ) (hl : 0 < lam) :
    (dτ * transferEig dτ lam)⁻¹ + transferEig dτ lam * (dτ * transferEig dτ lam)⁻¹
        * transferEig dτ lam = 2 / dτ + dτ * lam := by
  have hr := transferEig_pos hd hl
  have hsum := transferEig_add_inv hd hl
  have hinv : (transferEig dτ lam)⁻¹ = 2 + dτ ^ 2 * lam - transferEig dτ lam := by linarith
  have key : (dτ * transferEig dτ lam)⁻¹ + transferEig dτ lam * (dτ * transferEig dτ lam)⁻¹
      * transferEig dτ lam = ((transferEig dτ lam)⁻¹ + transferEig dτ lam) / dτ := by
    field_simp
  rw [key, add_comm, hsum]
  field_simp

/-- The time-neighbour block of the chain precision is `−q⁻¹ r`; here `q⁻¹ r = 1/Δτ`. -/
theorem chain_precision_off (hd : 0 < dτ) (hl : 0 < lam) :
    (dτ * transferEig dτ lam)⁻¹ * transferEig dτ lam = 1 / dτ := by
  have hr := transferEig_pos hd hl
  field_simp

/-- The two-time data return the unit time coupling: `Δτ (c0² / (r c0) − r c0)⁻¹ = 1`. -/
theorem markov_reader_unit (hd : 0 < dτ) (hl : 0 < lam) :
    dτ * (equalTimeVar dτ lam ^ 2 / (transferEig dτ lam * equalTimeVar dτ lam)
      - transferEig dτ lam * equalTimeVar dτ lam)⁻¹ = 1 := by
  have hr := transferEig_pos hd hl
  have hc := equalTimeVar_pos hd hl
  have hin := chain_innovation hd hl
  have key : equalTimeVar dτ lam ^ 2 / (transferEig dτ lam * equalTimeVar dτ lam)
      - transferEig dτ lam * equalTimeVar dτ lam
      = equalTimeVar dτ lam * (1 - transferEig dτ lam ^ 2) / transferEig dτ lam := by
    field_simp
  rw [key, hin]
  field_simp

/-- Lag-zero row of the precision: `(2/Δτ + Δτλ) c0 − (2/Δτ) r c0 = 1`. -/
theorem contact_identity_mode (hd : 0 < dτ) (hl : 0 < lam) :
    (2 / dτ + dτ * lam) * equalTimeVar dτ lam - 2 / dτ * (transferEig dτ lam * equalTimeVar dτ lam)
      = 1 := by
  have hr := transferEig_pos hd hl
  have hsum := transferEig_add_inv hd hl
  have hin := chain_innovation hd hl
  have hinv : (transferEig dτ lam)⁻¹ = 2 + dτ ^ 2 * lam - transferEig dτ lam := by linarith
  have key : (2 / dτ + dτ * lam) * equalTimeVar dτ lam
      - 2 / dτ * (transferEig dτ lam * equalTimeVar dτ lam)
      = equalTimeVar dτ lam * (1 - transferEig dτ lam ^ 2) / (dτ * transferEig dτ lam) := by
    have e : 2 / dτ + dτ * lam = ((transferEig dτ lam)⁻¹ + transferEig dτ lam) / dτ := by
      rw [add_comm ((transferEig dτ lam)⁻¹), hsum]
      field_simp
    rw [e]
    field_simp
    ring
  rw [key, hin]
  field_simp

end Mode

#print axioms transferEig_pos
#print axioms transferEig_lt_one
#print axioms transferEig_add_inv
#print axioms exp_neg_latticeFreq
#print axioms cosh_latticeFreq
#print axioms equalTimeVar_closed
#print axioms chain_innovation
#print axioms chain_precision_diag
#print axioms chain_precision_off
#print axioms markov_reader_unit
#print axioms contact_identity_mode

end OPH.SourceLaw
