import SourceLaw.StationaryWindow
import SourceLaw.MatrixLimits
import Mathlib.Analysis.SpecificLimits.Basic

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The open register: the end effect vanishes on interior slabs

The declared register over `T + 1` slabs, with time coupling equal to the capacities and open ends,
has stiffness `P_T`; the stationary covariance of the declared law on it is `ℏ P_T⁻¹`. In
capacity-whitened coordinates `P_T` is `historyPrecision 1 A Δτ T`. This file inverts `P_T`
explicitly and proves (`Δτ > 0`, `A` positive definite, `r = r(λ)`, `c0 = c0(λ)` per mode):

* `greenMode_inv`: one mode of the open precision times the kernel
  `G(u, s) = c0 r^|u−s| + a_s r^u + b_s r^(T−u)` (`greenMode`) is the identity; `a_s`, `b_s` are the
  coefficients of the two boundary solutions `r^u`, `r^(T−u)` (`bdA`, `bdB`).
* `greenMode_diag_bound`: `0 ≤ G(t, t) − c0 ≤ c0 r ((r^t)² + (r^(T−t))²)/(1 − r)`;
  `greenMode_interior_tendsto`: `G(t_k, t_k) → c0` whenever `t_k → ∞` and `T_k − t_k → ∞`.
* `greenMode_end_tendsto`, `greenMode_end_ne`: at an end slab `G(0, 0) → c0 (1 + r)`, which differs
  from `c0`.
* `register_inv`: `P_T⁻¹ = greenMat` (slab blocks `G(u, s)(A)`) for `T ≥ 1`;
  `interior_block_tendsto`: the equal-time block of `P_{T_k}⁻¹` at slab `t_k` tends to the
  stationary-chain covariance `C0(Δτ)`; `end_block_tendsto`, `end_block_ne`: at the first slab the
  block tends to `C0(Δτ)(1 + R)`, which differs from `C0(Δτ)`.
* Original coordinates: `register_conj`, `declared_register_inv`: for the register of a seam graph
  with capacities `M` the inverse is `(M^{-1/2} ⊗ 1) greenMat (M^{-1/2} ⊗ 1)`;
  `declared_interior_block_tendsto` and `declared_end_block_tendsto` give the limits
  `M^{-1/2} C0(Δτ) M^{-1/2}` and `M^{-1/2} C0(Δτ)(1 + R) M^{-1/2}`.
-/

namespace OPH.SourceLaw

open scoped BigOperators Kronecker
open Filter Topology Matrix

section Mode

variable {dτ lam : ℝ}

/-- One mode of the open register precision (no end-slab correction). -/
noncomputable def regPrecMode (dτ x : ℝ) (T : ℕ) (t u : Fin (T + 1)) : ℝ :=
  dτ * (if t = u then 1 else 0) * x + dτ⁻¹ * pathLap T t u

/-- The boundary solution `r^u` against one row of the open precision. -/
theorem regPrec_pow (hd : 0 < dτ) (hl : 0 < lam) {T : ℕ} (hT : 1 ≤ T) (t : Fin (T + 1)) :
    ∑ u : Fin (T + 1), regPrecMode dτ lam T t u * transferEig dτ lam ^ (u : ℕ)
      = dτ⁻¹ * (1 - transferEig dτ lam)
        * ((if (t : ℕ) = 0 then (transferEig dτ lam)⁻¹ else 0)
          - (if (t : ℕ) = T then transferEig dτ lam ^ T else 0)) := by
  have hq := transfer_quadratic hd hl
  have hr := transferEig_pos hd hl
  set r := transferEig dτ lam
  have hL := pathLap_sum_apply T (fun n => r ^ n) t
  simp only [regPrecMode, add_mul, Finset.sum_add_distrib, mul_assoc, ← Finset.mul_sum]
  rw [hL]
  have hdiag : ∑ u : Fin (T + 1), (if t = u then (1 : ℝ) else 0) * (lam * r ^ (u : ℕ))
      = lam * r ^ (t : ℕ) := by
    rw [Finset.sum_eq_single t]
    · simp
    · intro b _ hb; rw [if_neg (Ne.symm hb), zero_mul]
    · intro h; exact absurd (Finset.mem_univ t) h
  rw [hdiag]
  have ht := t.isLt
  have hdτ : dτ ≠ 0 := hd.ne'
  have hr0 : r ≠ 0 := hr.ne'
  rcases Nat.eq_zero_or_pos (t : ℕ) with h0 | hpos
  · have hlt : (t : ℕ) < T := by omega
    rw [if_pos hlt, if_neg (by omega), if_pos h0, if_neg (by omega), h0]
    field_simp
    simp only [pow_one, zero_add]
    nlinarith [hq]
  · by_cases hTt : (t : ℕ) = T
    · rw [if_neg (by omega), if_pos hpos, if_neg (by omega), if_pos hTt, hTt]
      obtain ⟨k, hk⟩ : ∃ k, T = k + 1 := ⟨T - 1, by omega⟩
      rw [hk, Nat.add_sub_cancel, pow_succ]
      field_simp
      nlinarith [hq, pow_pos hr k]
    · have hlt : (t : ℕ) < T := lt_of_le_of_ne (by omega) hTt
      rw [if_pos hlt, if_pos hpos, if_neg (by omega), if_neg hTt]
      obtain ⟨k, hk⟩ : ∃ k, (t : ℕ) = k + 1 := ⟨(t : ℕ) - 1, by omega⟩
      rw [hk, Nat.add_sub_cancel]
      field_simp
      linear_combination (r ^ k) * hq

/-- The reversed boundary solution `r^(T−u)` against one row of the open precision. -/
theorem regPrec_powRev (hd : 0 < dτ) (hl : 0 < lam) {T : ℕ} (hT : 1 ≤ T) (t : Fin (T + 1)) :
    ∑ u : Fin (T + 1), regPrecMode dτ lam T t u * transferEig dτ lam ^ (T - (u : ℕ))
      = dτ⁻¹ * (1 - transferEig dτ lam)
        * ((if (t : ℕ) = T then (transferEig dτ lam)⁻¹ else 0)
          - (if (t : ℕ) = 0 then transferEig dτ lam ^ T else 0)) := by
  have hq := transfer_quadratic hd hl
  have hr := transferEig_pos hd hl
  set r := transferEig dτ lam
  have hL := pathLap_sum_apply T (fun n => r ^ (T - n)) t
  simp only [regPrecMode, add_mul, Finset.sum_add_distrib, mul_assoc, ← Finset.mul_sum]
  rw [hL]
  have hdiag : ∑ u : Fin (T + 1), (if t = u then (1 : ℝ) else 0) * (lam * r ^ (T - (u : ℕ)))
      = lam * r ^ (T - (t : ℕ)) := by
    rw [Finset.sum_eq_single t]
    · simp
    · intro b _ hb; rw [if_neg (Ne.symm hb), zero_mul]
    · intro h; exact absurd (Finset.mem_univ t) h
  rw [hdiag]
  have ht := t.isLt
  have hdτ : dτ ≠ 0 := hd.ne'
  have hr0 : r ≠ 0 := hr.ne'
  rcases Nat.eq_zero_or_pos (t : ℕ) with h0 | hpos
  · have hlt : (t : ℕ) < T := by omega
    rw [if_pos hlt, if_neg (by omega), if_neg (by omega), if_pos h0, h0]
    obtain ⟨k, hk⟩ : ∃ k, T = k + 1 := ⟨T - 1, by omega⟩
    rw [hk, Nat.sub_zero, show k + 1 - (0 + 1) = k by omega]
    field_simp
    linear_combination (r ^ k) * hq
  · by_cases hTt : (t : ℕ) = T
    · rw [if_neg (by omega), if_pos hpos, if_pos hTt, if_neg (by omega), hTt]
      rw [Nat.sub_self, show T - (T - 1) = 1 by omega]
      field_simp
      linear_combination hq
    · have hlt : (t : ℕ) < T := lt_of_le_of_ne (by omega) hTt
      rw [if_pos hlt, if_pos hpos, if_neg hTt, if_neg (by omega)]
      obtain ⟨k, hk⟩ : ∃ k, T - ((t : ℕ) + 1) = k := ⟨_, rfl⟩
      rw [show T - (t : ℕ) = k + 1 by omega, show T - ((t : ℕ) - 1) = k + 2 by omega, hk]
      field_simp
      linear_combination (r ^ k) * hq

/-- Denominator `1 − r² (r^s)² (r^(T−s))² = 1 − r^(2T+2)` of the boundary coefficients. -/
noncomputable def bdDen (r : ℝ) (T s : ℕ) : ℝ := 1 - r ^ 2 * (r ^ s) ^ 2 * (r ^ (T - s)) ^ 2

/-- Coefficient of the boundary solution `r^u` in column `s`. -/
noncomputable def bdA (r c : ℝ) (T s : ℕ) : ℝ :=
  c * r * r ^ s * (1 + r * (r ^ (T - s)) ^ 2) / bdDen r T s

/-- Coefficient of the boundary solution `r^(T−u)` in column `s`. -/
noncomputable def bdB (r c : ℝ) (T s : ℕ) : ℝ :=
  c * r * r ^ (T - s) * (1 + r * (r ^ s) ^ 2) / bdDen r T s

/-- Explicit inverse of one mode of the open register precision. -/
noncomputable def greenMode (dτ lam : ℝ) (T : ℕ) (u s : Fin (T + 1)) : ℝ :=
  equalTimeVar dτ lam * transferEig dτ lam ^ tdist u s
    + bdA (transferEig dτ lam) (equalTimeVar dτ lam) T s * transferEig dτ lam ^ (u : ℕ)
    + bdB (transferEig dτ lam) (equalTimeVar dτ lam) T s * transferEig dτ lam ^ (T - (u : ℕ))

theorem bdDen_pos (hd : 0 < dτ) (hl : 0 < lam) (T s : ℕ) :
    0 < bdDen (transferEig dτ lam) T s := by
  have hr := transferEig_pos hd hl
  have hr1 := transferEig_lt_one hd hl
  unfold bdDen
  have hx : (transferEig dτ lam ^ s) ^ 2 ≤ 1 := by
    apply pow_le_one₀ (by positivity); exact pow_le_one₀ hr.le hr1.le
  have hy : (transferEig dτ lam ^ (T - s)) ^ 2 ≤ 1 := by
    apply pow_le_one₀ (by positivity); exact pow_le_one₀ hr.le hr1.le
  have hr2 : transferEig dτ lam ^ 2 < 1 := by nlinarith
  have : transferEig dτ lam ^ 2 * (transferEig dτ lam ^ s) ^ 2 * (transferEig dτ lam ^ (T - s)) ^ 2
      < 1 := by
    calc _ ≤ transferEig dτ lam ^ 2 * 1 * 1 := by gcongr
      _ < 1 := by simpa using hr2
  linarith

/-- The explicit inverse: one row of the open precision against one column of `greenMode`. -/
theorem greenMode_inv (hd : 0 < dτ) (hl : 0 < lam) {T : ℕ} (hT : 1 ≤ T) (t s : Fin (T + 1)) :
    ∑ u : Fin (T + 1), regPrecMode dτ lam T t u * greenMode dτ lam T u s
      = if t = s then 1 else 0 := by
  have hwin := window_mode hd hl T t s
  have hp := regPrec_pow hd hl hT t
  have hpr := regPrec_powRev hd hl hT t
  have hD := bdDen_pos hd hl T (s : ℕ)
  have hr := transferEig_pos hd hl
  set r := transferEig dτ lam with hr_def
  set c := equalTimeVar dτ lam with hc_def
  have hsplit : ∑ u : Fin (T + 1), regPrecMode dτ lam T t u * greenMode dτ lam T u s
      = ∑ u : Fin (T + 1), regPrecMode dτ lam T t u * (c * r ^ tdist u s)
        + bdA r c T s * ∑ u : Fin (T + 1), regPrecMode dτ lam T t u * r ^ (u : ℕ)
        + bdB r c T s * ∑ u : Fin (T + 1), regPrecMode dτ lam T t u * r ^ (T - (u : ℕ)) := by
    simp only [greenMode, mul_add, Finset.sum_add_distrib, Finset.mul_sum]
    congr 1
    · congr 1
      exact Finset.sum_congr rfl fun u _ => by ring
    · exact Finset.sum_congr rfl fun u _ => by ring
  have hwp : ∀ u, windowPrecMode dτ lam T t u = regPrecMode dτ lam T t u
      + dτ⁻¹ * endWeight T t * (if t = u then 1 else 0) * (1 - r) := by
    intro u
    simp only [windowPrecMode, regPrecMode]
    ring
  have htoep : ∑ u : Fin (T + 1), regPrecMode dτ lam T t u * (c * r ^ tdist u s)
      = (if t = s then 1 else 0) - dτ⁻¹ * endWeight T t * (1 - r) * (c * r ^ tdist t s) := by
    rw [← hwin]
    simp only [hwp, add_mul, Finset.sum_add_distrib]
    have he : ∑ u : Fin (T + 1), dτ⁻¹ * endWeight T t * (if t = u then 1 else 0) * (1 - r)
        * (c * r ^ tdist u s) = dτ⁻¹ * endWeight T t * (1 - r) * (c * r ^ tdist t s) := by
      rw [Finset.sum_eq_single t]
      · simp only [if_true]; ring
      · intro b _ hb; rw [if_neg (Ne.symm hb)]; ring
      · intro h; exact absurd (Finset.mem_univ t) h
    rw [he]
    ring
  rw [hsplit, htoep, hp, hpr]
  simp only [bdA, bdB]
  have hsT : (s : ℕ) ≤ T := Nat.lt_succ_iff.mp s.isLt
  have hspl : r ^ T = r ^ (s : ℕ) * r ^ (T - (s : ℕ)) := (pow_mul_pow_sub r hsT).symm
  have hDne : bdDen r T s ≠ 0 := hD.ne'
  unfold bdDen at hDne ⊢
  have hr0 : r ≠ 0 := hr.ne'
  have ht := t.isLt
  rcases Nat.eq_zero_or_pos (t : ℕ) with h0 | hpos
  · have hne : (t : ℕ) ≠ T := by omega
    have hew : endWeight T t = 1 := by
      unfold endWeight; rw [if_pos h0, if_neg hne]; norm_num
    have htd : tdist (t : ℕ) (s : ℕ) = (s : ℕ) := by unfold tdist; rw [h0]; simp
    rw [hew, htd, if_pos h0, if_pos h0, if_neg hne, if_neg hne, hspl]
    field_simp
    ring
  · have hne0 : (t : ℕ) ≠ 0 := by omega
    by_cases hTt : (t : ℕ) = T
    · have hew : endWeight T t = 1 := by
        unfold endWeight; rw [if_neg hne0, if_pos hTt]; norm_num
      have htd : tdist (t : ℕ) (s : ℕ) = T - (s : ℕ) := by
        unfold tdist; rw [hTt]; split_ifs with h <;> omega
      rw [hew, htd, if_neg hne0, if_neg hne0, if_pos hTt, if_pos hTt, hspl]
      generalize r ^ (s : ℕ) = x at hDne ⊢
      generalize r ^ (T - (s : ℕ)) = y at hDne ⊢
      have hD2 : 1 - r ^ 2 * y ^ 2 * x ^ 2 ≠ 0 := by
        rwa [show 1 - r ^ 2 * y ^ 2 * x ^ 2 = 1 - r ^ 2 * x ^ 2 * y ^ 2 by ring]
      field_simp
      ring
    · have hew : endWeight T t = 0 := by
        unfold endWeight; rw [if_neg hne0, if_neg hTt]; norm_num
      rw [hew, if_neg hne0, if_neg hne0, if_neg hTt, if_neg hTt]
      ring

/-- Diagonal value of the explicit inverse. -/
theorem greenMode_diag {T : ℕ} (t : Fin (T + 1)) :
    greenMode dτ lam T t t = equalTimeVar dτ lam
      + bdA (transferEig dτ lam) (equalTimeVar dτ lam) T t * transferEig dτ lam ^ (t : ℕ)
      + bdB (transferEig dτ lam) (equalTimeVar dτ lam) T t
        * transferEig dτ lam ^ (T - (t : ℕ)) := by
  unfold greenMode
  have : tdist (t : ℕ) (t : ℕ) = 0 := by unfold tdist; simp
  rw [this, pow_zero, mul_one]

/-- The end effect at slab `t`: `0 ≤ G(t, t) − c0 ≤ c0 r ((r^t)² + (r^(T−t))²)/(1 − r)`. -/
theorem greenMode_diag_bound (hd : 0 < dτ) (hl : 0 < lam) {T : ℕ} (t : Fin (T + 1)) :
    0 ≤ greenMode dτ lam T t t - equalTimeVar dτ lam
      ∧ greenMode dτ lam T t t - equalTimeVar dτ lam
        ≤ equalTimeVar dτ lam * transferEig dτ lam
          * ((transferEig dτ lam ^ (t : ℕ)) ^ 2 + (transferEig dτ lam ^ (T - (t : ℕ))) ^ 2)
          / (1 - transferEig dτ lam) := by
  have hr := transferEig_pos hd hl
  have hr1 := transferEig_lt_one hd hl
  have hc : 0 < equalTimeVar dτ lam := by
    rw [equalTimeVar_closed hd hl]; positivity
  have hx0 : 0 < transferEig dτ lam ^ (t : ℕ) := pow_pos hr _
  have hx1 : transferEig dτ lam ^ (t : ℕ) ≤ 1 := pow_le_one₀ hr.le hr1.le
  have hy0 : 0 < transferEig dτ lam ^ (T - (t : ℕ)) := pow_pos hr _
  have hy1 : transferEig dτ lam ^ (T - (t : ℕ)) ≤ 1 := pow_le_one₀ hr.le hr1.le
  have hD := bdDen_pos hd hl T t
  unfold bdDen at hD
  rw [greenMode_diag]
  unfold bdA bdB bdDen
  generalize transferEig dτ lam ^ (t : ℕ) = x at hx0 hx1 hD ⊢
  generalize transferEig dτ lam ^ (T - (t : ℕ)) = y at hy0 hy1 hD ⊢
  generalize equalTimeVar dτ lam = c at hc ⊢
  generalize transferEig dτ lam = r at hr hr1 hD ⊢
  have hDne : 1 - r ^ 2 * x ^ 2 * y ^ 2 ≠ 0 := hD.ne'
  have hsum : c + c * r * x * (1 + r * y ^ 2) / (1 - r ^ 2 * x ^ 2 * y ^ 2) * x
        + c * r * y * (1 + r * x ^ 2) / (1 - r ^ 2 * x ^ 2 * y ^ 2) * y - c
      = c * r * (x ^ 2 * (1 + r * y ^ 2) + y ^ 2 * (1 + r * x ^ 2))
        / (1 - r ^ 2 * x ^ 2 * y ^ 2) := by
    field_simp
    ring
  rw [hsum]
  have hN : 0 ≤ c * r * (x ^ 2 * (1 + r * y ^ 2) + y ^ 2 * (1 + r * x ^ 2)) := by positivity
  refine ⟨div_nonneg hN hD.le, ?_⟩
  have h1 : 0 ≤ x ^ 2 * (1 - y ^ 2) := mul_nonneg (sq_nonneg x) (by nlinarith)
  have h2 : 0 ≤ y ^ 2 * (1 - x ^ 2) := mul_nonneg (sq_nonneg y) (by nlinarith)
  have hnum : x ^ 2 * (1 + r * y ^ 2) + y ^ 2 * (1 + r * x ^ 2) ≤ (1 + r) * (x ^ 2 + y ^ 2) := by
    nlinarith [mul_nonneg hr.le (add_nonneg h1 h2)]
  have hxy : x ^ 2 * y ^ 2 ≤ 1 := by
    have hx2 : x ^ 2 ≤ 1 := by nlinarith
    have hy2 : y ^ 2 ≤ 1 := by nlinarith
    calc x ^ 2 * y ^ 2 ≤ 1 * 1 := by gcongr
      _ = 1 := by ring
  have hden : 1 - r ^ 2 ≤ 1 - r ^ 2 * x ^ 2 * y ^ 2 := by
    have : r ^ 2 * (x ^ 2 * y ^ 2) ≤ r ^ 2 * 1 := mul_le_mul_of_nonneg_left hxy (sq_nonneg r)
    nlinarith
  have hr2 : 0 < 1 - r ^ 2 := by nlinarith
  have hr1' : 0 < 1 - r := by linarith
  calc c * r * (x ^ 2 * (1 + r * y ^ 2) + y ^ 2 * (1 + r * x ^ 2)) / (1 - r ^ 2 * x ^ 2 * y ^ 2)
      ≤ c * r * ((1 + r) * (x ^ 2 + y ^ 2)) / (1 - r ^ 2 * x ^ 2 * y ^ 2) :=
        div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left hnum (by positivity)) hD.le
    _ ≤ c * r * ((1 + r) * (x ^ 2 + y ^ 2)) / (1 - r ^ 2) :=
        div_le_div_of_nonneg_left (by positivity) hr2 hden
    _ = c * r * (x ^ 2 + y ^ 2) / (1 - r) := by
        have hr1ne : 1 - r ≠ 0 := hr1'.ne'
        have hr2ne : 1 - r ^ 2 ≠ 0 := hr2.ne'
        field_simp
        ring

/-- Interior slabs: the end effect vanishes when both distances to the ends grow. -/
theorem greenMode_interior_tendsto (hd : 0 < dτ) (hl : 0 < lam) {T t : ℕ → ℕ}
    (ht : ∀ k, t k ≤ T k) (h1 : Tendsto t atTop atTop)
    (h2 : Tendsto (fun k => T k - t k) atTop atTop) :
    Tendsto (fun k => greenMode dτ lam (T k) ⟨t k, Nat.lt_succ_of_le (ht k)⟩
      ⟨t k, Nat.lt_succ_of_le (ht k)⟩) atTop (𝓝 (equalTimeVar dτ lam)) := by
  have hr := transferEig_pos hd hl
  have hr1 := transferEig_lt_one hd hl
  have hp : Tendsto (fun n : ℕ => transferEig dτ lam ^ n) atTop (𝓝 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one hr.le hr1
  have hx := hp.comp h1
  have hy := hp.comp h2
  have hub : Tendsto (fun k => equalTimeVar dτ lam * transferEig dτ lam
      * ((transferEig dτ lam ^ t k) ^ 2 + (transferEig dτ lam ^ (T k - t k)) ^ 2)
      / (1 - transferEig dτ lam)) atTop (𝓝 0) := by
    have := ((tendsto_const_nhds (x := equalTimeVar dτ lam * transferEig dτ lam)).mul
      ((hx.pow 2).add (hy.pow 2))).div_const (1 - transferEig dτ lam)
    simpa using this
  have hdiff : Tendsto (fun k => greenMode dτ lam (T k) ⟨t k, Nat.lt_succ_of_le (ht k)⟩
      ⟨t k, Nat.lt_succ_of_le (ht k)⟩ - equalTimeVar dτ lam) atTop (𝓝 0) :=
    squeeze_zero (fun k => (greenMode_diag_bound hd hl _).1)
      (fun k => (greenMode_diag_bound hd hl _).2) hub
  have := hdiff.add_const (equalTimeVar dτ lam)
  simpa using this

/-- At the middle slab of a window of `2n + 1` slabs the equal-time variance tends to `c0`. -/
theorem greenMode_middle_tendsto (hd : 0 < dτ) (hl : 0 < lam) :
    Tendsto (fun n : ℕ => greenMode dτ lam (2 * n) ⟨n, by omega⟩ ⟨n, by omega⟩) atTop
      (𝓝 (equalTimeVar dτ lam)) := by
  have h2 : Tendsto (fun n : ℕ => 2 * n - n) atTop atTop := by
    have e : (fun n : ℕ => 2 * n - n) = id := by funext n; simp only [id]; omega
    rw [e]
    exact tendsto_id
  exact greenMode_interior_tendsto hd hl (T := fun n => 2 * n) (t := fun n => n)
    (fun n => by show n ≤ 2 * n; omega) tendsto_id h2

/-- At an end slab the equal-time variance tends to `c0 (1 + r)`. -/
theorem greenMode_end_tendsto (hd : 0 < dτ) (hl : 0 < lam) :
    Tendsto (fun T : ℕ => greenMode dτ lam T 0 0) atTop
      (𝓝 (equalTimeVar dτ lam * (1 + transferEig dτ lam))) := by
  have hr := transferEig_pos hd hl
  have hr1 := transferEig_lt_one hd hl
  set r := transferEig dτ lam
  set c := equalTimeVar dτ lam
  have hx : Tendsto (fun n : ℕ => r ^ n) atTop (𝓝 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one hr.le hr1
  have hform : ∀ T : ℕ, greenMode dτ lam T 0 0
      = c + c * r * (1 + r * (r ^ T) ^ 2) / (1 - r ^ 2 * (r ^ T) ^ 2)
        + c * r * r ^ T * (1 + r) / (1 - r ^ 2 * (r ^ T) ^ 2) * r ^ T := by
    intro T
    rw [greenMode_diag]
    unfold bdA bdB bdDen
    simp
    rfl
  simp only [hform]
  have hD : Tendsto (fun T : ℕ => 1 - r ^ 2 * (r ^ T) ^ 2) atTop (𝓝 1) := by
    have := (tendsto_const_nhds (x := r ^ 2)).mul (hx.pow 2)
    simpa using (tendsto_const_nhds (x := (1 : ℝ))).sub this
  have h1 : Tendsto (fun T : ℕ => c * r * (1 + r * (r ^ T) ^ 2)) atTop (𝓝 (c * r)) := by
    have := (tendsto_const_nhds (x := c * r)).mul
      ((tendsto_const_nhds (x := (1 : ℝ))).add ((tendsto_const_nhds (x := r)).mul (hx.pow 2)))
    simpa using this
  have h2 : Tendsto (fun T : ℕ => c * r * r ^ T * (1 + r)) atTop (𝓝 0) := by
    have := ((tendsto_const_nhds (x := c * r)).mul hx).mul (tendsto_const_nhds (x := 1 + r))
    simpa using this
  have := ((tendsto_const_nhds (x := c)).add (h1.div hD one_ne_zero)).add
    ((h2.div hD one_ne_zero).mul hx)
  convert this using 2
  ring

/-- The end-slab limit differs from the stationary value. -/
theorem greenMode_end_ne (hd : 0 < dτ) (hl : 0 < lam) :
    equalTimeVar dτ lam * (1 + transferEig dτ lam) ≠ equalTimeVar dτ lam := by
  have hr := transferEig_pos hd hl
  have hc : 0 < equalTimeVar dτ lam := by
    rw [equalTimeVar_closed hd hl]; positivity
  intro h
  nlinarith

end Mode

section Matrix

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {A : Matrix ι ι ℝ} (hA : A.PosDef) {dτ : ℝ}

/-- The equal-time block of slab `t`. -/
def slabBlock {T : ℕ} (X : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ) (t : Fin (T + 1)) :
    Matrix ι ι ℝ :=
  Matrix.of fun i k => X (i, t) (k, t)

omit [Fintype ι] [DecidableEq ι] in
theorem slabBlock_blockOf {T : ℕ} (F : Fin (T + 1) → Fin (T + 1) → Matrix ι ι ℝ)
    (t : Fin (T + 1)) : slabBlock (blockOf F) t = F t t := by
  ext i k
  rfl

theorem regPrec_blocks (T : ℕ) :
    historyPrecision 1 A dτ T = blockOf fun t u => hA.1.cfc fun x => regPrecMode dτ x T t u := by
  ext ⟨i, t⟩ ⟨k, u⟩
  have hc : hA.1.cfc (fun x => regPrecMode dτ x T t u)
      = (dτ * (if t = u then 1 else 0)) • A + (dτ⁻¹ * pathLap T t u) • (1 : Matrix ι ι ℝ) := by
    rw [smul_eq_cfc hA.1, ← cfc_one hA.1, cfc_smul hA.1, cfc_add hA.1]
    exact cfc_congr hA.1 fun i => by unfold regPrecMode; ring
  simp only [blockOf, Matrix.of_apply, hc, historyPrecision_apply, Matrix.add_apply,
    Matrix.smul_apply, smul_eq_mul]
  by_cases h : t = u
  · subst h
    rw [Matrix.one_apply_eq, if_pos rfl]
    ring
  · rw [Matrix.one_apply_ne h, if_neg h]
    ring

/-- The explicit inverse of the open register precision, by slab blocks. -/
noncomputable def greenMat (dτ : ℝ) (T : ℕ) : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ :=
  blockOf fun u s => hA.1.cfc fun x => greenMode dτ x T u s

theorem historyPrecision_mul_greenMat (hd : 0 < dτ) {T : ℕ} (hT : 1 ≤ T) :
    historyPrecision 1 A dτ T * greenMat hA dτ T = 1 := by
  rw [regPrec_blocks hA, greenMat, blockOf_mul, ← blockOf_delta]
  congr 1
  funext t s
  simp only [cfc_mul hA.1]
  rw [cfc_sum hA.1]
  have hconst : (if t = s then (1 : Matrix ι ι ℝ) else 0)
      = hA.1.cfc (fun _ => if t = s then (1 : ℝ) else 0) := by
    split_ifs
    · exact (cfc_one hA.1).symm
    · rw [cfc_const hA.1 0, zero_smul]
  rw [hconst]
  exact cfc_congr hA.1 fun i => greenMode_inv hd (eig_pos hA i) hT t s

/-- `P_T⁻¹ = greenMat` in whitened coordinates. -/
theorem register_inv (hd : 0 < dτ) {T : ℕ} (hT : 1 ≤ T) :
    (historyPrecision 1 A dτ T)⁻¹ = greenMat hA dτ T :=
  Matrix.inv_eq_right_inv (historyPrecision_mul_greenMat hA hd hT)

/-- Interior slabs: the equal-time block of `P_{T_k}⁻¹` at slab `t_k` tends to `C0(Δτ)`. -/
theorem interior_block_tendsto (hd : 0 < dτ) {T t : ℕ → ℕ} (ht : ∀ k, t k ≤ T k)
    (h1 : Tendsto t atTop atTop) (h2 : Tendsto (fun k => T k - t k) atTop atTop) :
    Tendsto (fun k => slabBlock (historyPrecision 1 A dτ (T k))⁻¹
      ⟨t k, Nat.lt_succ_of_le (ht k)⟩) atTop (𝓝 (equalTimeCov hA dτ)) := by
  have hev : ∀ᶠ k in atTop, hA.1.cfc (fun x => greenMode dτ x (T k)
      ⟨t k, Nat.lt_succ_of_le (ht k)⟩ ⟨t k, Nat.lt_succ_of_le (ht k)⟩)
      = slabBlock (historyPrecision 1 A dτ (T k))⁻¹ ⟨t k, Nat.lt_succ_of_le (ht k)⟩ := by
    filter_upwards [h2.eventually_ge_atTop 1] with k hk
    rw [register_inv hA hd (by omega), greenMat, slabBlock_blockOf]
  refine Tendsto.congr' hev ?_
  exact cfc_tendsto hA.1 (f := fun k x => greenMode dτ x (T k) ⟨t k, Nat.lt_succ_of_le (ht k)⟩
    ⟨t k, Nat.lt_succ_of_le (ht k)⟩) fun i => greenMode_interior_tendsto hd (eig_pos hA i) ht h1 h2

theorem cfc_end_value :
    hA.1.cfc (fun x => equalTimeVar dτ x * (1 + transferEig dτ x))
      = equalTimeCov hA dτ * (1 + transferMat hA dτ) := by
  rw [equalTimeCov, transferMat, ← cfc_one hA.1, cfc_add hA.1, cfc_mul hA.1]

/-- At the first slab the block tends to `C0(Δτ)(1 + R)`. -/
theorem end_block_tendsto (hd : 0 < dτ) :
    Tendsto (fun T : ℕ => slabBlock (historyPrecision 1 A dτ T)⁻¹ 0) atTop
      (𝓝 (equalTimeCov hA dτ * (1 + transferMat hA dτ))) := by
  have hev : ∀ᶠ T in atTop, hA.1.cfc (fun x => greenMode dτ x T 0 0)
      = slabBlock (historyPrecision 1 A dτ T)⁻¹ 0 := by
    filter_upwards [eventually_ge_atTop 1] with T hT
    rw [register_inv hA hd hT, greenMat, slabBlock_blockOf]
  rw [← cfc_end_value hA]
  exact Tendsto.congr' hev (cfc_tendsto hA.1 (f := fun T x => greenMode dτ x T 0 0) fun i =>
    greenMode_end_tendsto hd (eig_pos hA i))

/-- The end-slab limit is not the stationary-chain covariance. -/
theorem end_block_ne [Nonempty ι] (hd : 0 < dτ) :
    equalTimeCov hA dτ * (1 + transferMat hA dτ) ≠ equalTimeCov hA dτ := by
  intro h
  have hpos : (hA.1.cfc (fun x => equalTimeVar dτ x * transferEig dτ x)).PosDef :=
    cfc_posDef hA.1 fun i => by
      have hl := eig_pos hA i
      have hc : 0 < equalTimeVar dτ (hA.1.eigenvalues i) := by
        rw [equalTimeVar_closed hd hl]; positivity
      exact mul_pos hc (transferEig_pos hd hl)
  have h0 : hA.1.cfc (fun x => equalTimeVar dτ x * transferEig dτ x) = 0 := by
    have e1 : hA.1.cfc (fun x => equalTimeVar dτ x * transferEig dτ x)
        = equalTimeCov hA dτ * transferMat hA dτ := (cfc_mul hA.1 _ _).symm
    have e : equalTimeCov hA dτ * transferMat hA dτ
        = equalTimeCov hA dτ * (1 + transferMat hA dτ) - equalTimeCov hA dτ := by
      rw [Matrix.mul_add, Matrix.mul_one]
      abel
    rw [e1, e, h, sub_self]
  have hdet := hpos.det_pos
  rw [h0, Matrix.det_zero ‹Nonempty ι›] at hdet
  exact lt_irrefl 0 hdet

end Matrix

section Original

variable {ι E : Type*} [Fintype ι] [DecidableEq ι] [Fintype E] {dτ : ℝ}

/-- Slab blocks of a slab-wise diagonal conjugation. -/
theorem slabBlock_conj (d : ι → ℝ) {T : ℕ} (X : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ)
    (t : Fin (T + 1)) :
    slabBlock ((diagonal d ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) * X
        * (diagonal d ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))) t
      = diagonal d * slabBlock X t * diagonal d := by
  have hk : diagonal d ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)
      = diagonal (fun p : ι × Fin (T + 1) => d p.1) := by
    rw [← diagonal_one, diagonal_kronecker_diagonal]
    congr 1
    funext p
    simp
  rw [hk]
  ext i k
  simp only [slabBlock, Matrix.of_apply, Matrix.mul_diagonal, Matrix.diagonal_mul]

theorem slabBlock_conj_jInvSqrt (M : ι → ℝ) {T : ℕ}
    (X : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ) (t : Fin (T + 1)) :
    slabBlock ((jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) * X
        * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))) t
      = jInvSqrt M * slabBlock X t * jInvSqrt M :=
  slabBlock_conj _ X t

/-- In whitened coordinates the register stiffness with time coupling `M` is the open history
precision. -/
theorem register_conj (src tgt : E → ι) (c : E → ℝ) (g M : ι → ℝ) (hM : ∀ i, 0 < M i) (T : ℕ) :
    stiffnessMatrix (regSrc src T) (regTgt tgt T) (regCond c M dτ T) (regGround g dτ T)
      = (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
        * historyPrecision 1 (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M) dτ T
        * (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) := by
  have h1 : jSqrt M * jInvSqrt M = 1 := jSqrt_mul_jInvSqrt hM
  have hV : stiffnessMatrix src tgt c g
      = jSqrt M * (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M) * jSqrt M := by
    calc stiffnessMatrix src tgt c g
        = (jSqrt M * jInvSqrt M) * stiffnessMatrix src tgt c g
            * (jInvSqrt M * jSqrt M) := by
          rw [h1, jInvSqrt_mul_jSqrt hM, Matrix.one_mul, Matrix.mul_one]
      _ = jSqrt M * (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M) * jSqrt M := by
          simp only [Matrix.mul_assoc]
  rw [historyPrecision_conj, register_eq_historyPrecision, Matrix.mul_one,
    jSqrt_mul_jSqrt hM, ← hV]

variable (src tgt : E → ι) (c : E → ℝ) (g M : ι → ℝ)
  (hA : (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M).PosDef)

/-- The inverse of the declared register stiffness, for `T ≥ 1`. -/
theorem declared_register_inv (hM : ∀ i, 0 < M i) (hd : 0 < dτ) {T : ℕ} (hT : 1 ≤ T) :
    (stiffnessMatrix (regSrc src T) (regTgt tgt T) (regCond c M dτ T) (regGround g dτ T))⁻¹
      = (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) * greenMat hA dτ T
        * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) := by
  apply Matrix.inv_eq_right_inv
  rw [register_conj src tgt c g M hM T]
  have hK1 : (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
      * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) = 1 := by
    rw [← mul_kronecker_mul, jSqrt_mul_jInvSqrt hM, Matrix.mul_one, one_kronecker_one]
  have hG := historyPrecision_mul_greenMat hA hd hT
  calc (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
        * historyPrecision 1 (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M) dτ T
        * (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
        * ((jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) * greenMat hA dτ T
          * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)))
      = (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
        * (historyPrecision 1 (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M) dτ T
          * ((jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
            * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)))
          * greenMat hA dτ T)
        * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) := by
        simp only [Matrix.mul_assoc]
    _ = 1 := by rw [hK1, Matrix.mul_one, hG, Matrix.mul_one, hK1]

/-- Original coordinates, interior slabs: the block tends to `M^{-1/2} C0(Δτ) M^{-1/2}`. -/
theorem declared_interior_block_tendsto (hM : ∀ i, 0 < M i) (hd : 0 < dτ) {T t : ℕ → ℕ}
    (ht : ∀ k, t k ≤ T k)
    (h1 : Tendsto t atTop atTop) (h2 : Tendsto (fun k => T k - t k) atTop atTop) :
    Tendsto (fun k => slabBlock (stiffnessMatrix (regSrc src (T k)) (regTgt tgt (T k))
        (regCond c M dτ (T k)) (regGround g dτ (T k)))⁻¹ ⟨t k, Nat.lt_succ_of_le (ht k)⟩) atTop
      (𝓝 (jInvSqrt M * equalTimeCov hA dτ * jInvSqrt M)) := by
  have hev : ∀ᶠ k in atTop, jInvSqrt M * hA.1.cfc (fun x => greenMode dτ x (T k)
      ⟨t k, Nat.lt_succ_of_le (ht k)⟩ ⟨t k, Nat.lt_succ_of_le (ht k)⟩) * jInvSqrt M
      = slabBlock (stiffnessMatrix (regSrc src (T k)) (regTgt tgt (T k))
        (regCond c M dτ (T k)) (regGround g dτ (T k)))⁻¹ ⟨t k, Nat.lt_succ_of_le (ht k)⟩ := by
    filter_upwards [h2.eventually_ge_atTop 1] with k hk
    rw [declared_register_inv src tgt c g M hA hM hd (by omega), slabBlock_conj_jInvSqrt,
      greenMat, slabBlock_blockOf]
  refine Tendsto.congr' hev ?_
  exact ((tendsto_const_nhds (x := jInvSqrt M)).mul (cfc_tendsto hA.1
    (f := fun k x => greenMode dτ x (T k) ⟨t k, Nat.lt_succ_of_le (ht k)⟩
      ⟨t k, Nat.lt_succ_of_le (ht k)⟩) fun i =>
      greenMode_interior_tendsto hd (eig_pos hA i) ht h1 h2)).mul tendsto_const_nhds

/-- Original coordinates, first slab: the block tends to `M^{-1/2} C0(Δτ)(1 + R) M^{-1/2}`. -/
theorem declared_end_block_tendsto (hM : ∀ i, 0 < M i) (hd : 0 < dτ) :
    Tendsto (fun T : ℕ => slabBlock (stiffnessMatrix (regSrc src T) (regTgt tgt T)
        (regCond c M dτ T) (regGround g dτ T))⁻¹ 0) atTop
      (𝓝 (jInvSqrt M * (equalTimeCov hA dτ * (1 + transferMat hA dτ)) * jInvSqrt M)) := by
  have hev : ∀ᶠ T in atTop, jInvSqrt M * hA.1.cfc (fun x => greenMode dτ x T 0 0) * jInvSqrt M
      = slabBlock (stiffnessMatrix (regSrc src T) (regTgt tgt T)
        (regCond c M dτ T) (regGround g dτ T))⁻¹ 0 := by
    filter_upwards [eventually_ge_atTop 1] with T hT
    rw [declared_register_inv src tgt c g M hA hM hd hT, slabBlock_conj_jInvSqrt, greenMat,
      slabBlock_blockOf]
  rw [← cfc_end_value hA]
  refine Tendsto.congr' hev ?_
  exact ((tendsto_const_nhds (x := jInvSqrt M)).mul (cfc_tendsto hA.1
    (f := fun T x => greenMode dτ x T 0 0) fun i =>
      greenMode_end_tendsto hd (eig_pos hA i))).mul tendsto_const_nhds

end Original

#print axioms regPrec_pow
#print axioms regPrec_powRev
#print axioms bdDen_pos
#print axioms greenMode_inv
#print axioms greenMode_diag
#print axioms greenMode_diag_bound
#print axioms greenMode_interior_tendsto
#print axioms greenMode_middle_tendsto
#print axioms greenMode_end_tendsto
#print axioms greenMode_end_ne
#print axioms slabBlock_blockOf
#print axioms regPrec_blocks
#print axioms historyPrecision_mul_greenMat
#print axioms register_inv
#print axioms interior_block_tendsto
#print axioms cfc_end_value
#print axioms end_block_tendsto
#print axioms end_block_ne
#print axioms slabBlock_conj
#print axioms slabBlock_conj_jInvSqrt
#print axioms register_conj
#print axioms declared_register_inv
#print axioms declared_interior_block_tendsto
#print axioms declared_end_block_tendsto

end OPH.SourceLaw
