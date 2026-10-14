import SourceLaw.RepairAdmissibility
import SourceLaw.TransferMatrixLift
import SourceLaw.LangevinStationarity
import Mathlib.Analysis.Matrix.Order

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The declared register and the stationary chain on a finite window

The stationary chain whose interior precision is that of the register has the block Toeplitz
covariance `C0 R^|u−s|` on any window of `T + 1` consecutive slabs. The finite open register of the
declared law has the same interior precision but differs at its two end slabs. This file computes
the exact end-slab term that closes that gap (`Δτ > 0` and a positive definite whitened stiffness
`A` throughout).

Main results:

* `window_mode`: for one mode, the open-register precision with the end-slab correction
  `(1 − r)/Δτ` times the Toeplitz kernel `c0 r^|u−s|` is the identity.
* `stationary_window`: in whitened coordinates,
  `(P + E((1 − R)/Δτ)) Σ = 1`, where `P` is the register precision, `E(X)` adds `X` at each end
  slab and `Σ` is the window covariance of the stationary chain.
* `window_correction_unique`: among terms `D` added equally at both end slabs, the one that makes
  the window exact is unique, `D = (1 − R)/Δτ`.
* `window_boundary_not_diag`, `no_nodeLocal_window_correction`: when the whitened stiffness is not
  diagonal, that correction is not node-local, so no node-local term `D` added equally at both end
  slabs (`endSlabs`) makes the finite register precision the inverse of the stationary window
  covariance.
* `windowPrec_posDef`, `corrected_law_stationary`, `corrected_law_stationary_unique`: the
  corrected precision `P_c = P + E((1 − R)/Δτ)` is positive definite, `ℏ P_c⁻¹ = ℏ Σ`, and every
  solution of the Lyapunov equation of the declared Langevin law (whitened coordinates, unit
  capacities) with precision `P_c` is `ℏ Σ`.
* `declared_register_window`: in original coordinates, the stiffness of the Euclidean history
  register with time coupling `M`, plus the end-slab correction, is the inverse of the stationary
  window covariance.
-/

namespace OPH.SourceLaw

open scoped BigOperators Kronecker
open Matrix

/-- Slab distance `|u − s|`. -/
def tdist (u s : ℕ) : ℕ := if u ≤ s then s - u else u - s

/-- End-slab weight: one for each end of the window that slab `t` is. -/
def endWeight (T : ℕ) (t : Fin (T + 1)) : ℝ :=
  (if (t : ℕ) = 0 then 1 else 0) + (if (t : ℕ) = T then 1 else 0)

theorem endWeight_nonneg (T : ℕ) (t : Fin (T + 1)) : 0 ≤ endWeight T t := by
  unfold endWeight
  split_ifs <;> norm_num

theorem endWeight_zero_ne (T : ℕ) : endWeight T 0 ≠ 0 := by
  unfold endWeight
  simp only [Fin.val_zero, if_true]
  split_ifs <;> norm_num

/-- The path Laplacian applied to a slab profile. -/
theorem pathLap_sum_apply (T : ℕ) (v : ℕ → ℝ) (t : Fin (T + 1)) :
    ∑ u : Fin (T + 1), pathLap T t u * v u
      = (if (t : ℕ) < T then v t - v (t + 1) else 0)
        + (if 0 < (t : ℕ) then v t - v (t - 1) else 0) := by
  have hmv : ∑ u : Fin (T + 1), pathLap T t u * v u
      = (pathLap T *ᵥ fun u : Fin (T + 1) => v u) t := rfl
  have ht := t.isLt
  rw [hmv, pathLap, stiffnessMatrix_mulVec]
  simp only [stiffness, one_mul, zero_mul, add_zero, Finset.sum_add_distrib, Fin.val_castSucc,
    Fin.val_succ]
  congr 1
  · split_ifs with h
    · rw [Finset.sum_eq_single ⟨t, h⟩]
      · simp
      · intro w _ hw
        rw [if_neg]
        intro hc
        have := congrArg Fin.val hc
        simp only [Fin.val_castSucc] at this
        exact hw (Fin.ext (by show (w : ℕ) = (t : ℕ); omega))
      · intro hn
        exact absurd (Finset.mem_univ _) hn
    · refine Finset.sum_eq_zero fun w _ => if_neg ?_
      intro hc
      have := congrArg Fin.val hc
      simp only [Fin.val_castSucc] at this
      have := w.isLt
      omega
  · split_ifs with h
    · rw [Finset.sum_eq_single ⟨(t : ℕ) - 1, by omega⟩]
      · simp [Nat.sub_add_cancel (show 1 ≤ (t : ℕ) by omega)]
      · intro w _ hw
        rw [if_neg]
        intro hc
        have := congrArg Fin.val hc
        simp only [Fin.val_succ] at this
        exact hw (Fin.ext (by show (w : ℕ) = (t : ℕ) - 1; omega))
      · intro hn
        exact absurd (Finset.mem_univ _) hn
    · refine Finset.sum_eq_zero fun w _ => if_neg ?_
      intro hc
      have := congrArg Fin.val hc
      simp only [Fin.val_succ] at this
      omega

section Mode

variable {dτ lam : ℝ}

/-- `(2 + Δτ² λ) r = 1 + r²`. -/
theorem transfer_quadratic (hd : 0 < dτ) (hl : 0 < lam) :
    (2 + dτ ^ 2 * lam) * transferEig dτ lam = 1 + transferEig dτ lam ^ 2 := by
  have hr := (transferEig_pos hd hl).ne'
  rw [← transferEig_add_inv hd hl, add_mul, inv_mul_cancel₀ hr]
  ring

/-- One row of the corrected window precision against the Toeplitz kernel. -/
theorem window_row_mode (hd : 0 < dτ) (hl : 0 < lam) {n a b : ℕ} (ha : a ≤ n) (hb : b ≤ n) :
    dτ * lam * (equalTimeVar dτ lam * transferEig dτ lam ^ tdist a b)
      + dτ⁻¹ * ((if a < n then equalTimeVar dτ lam * transferEig dτ lam ^ tdist a b
            - equalTimeVar dτ lam * transferEig dτ lam ^ tdist (a + 1) b else 0)
          + (if 0 < a then equalTimeVar dτ lam * transferEig dτ lam ^ tdist a b
            - equalTimeVar dτ lam * transferEig dτ lam ^ tdist (a - 1) b else 0))
      + dτ⁻¹ * ((if a = 0 then 1 else 0) + (if a = n then 1 else 0)) * (1 - transferEig dτ lam)
          * (equalTimeVar dτ lam * transferEig dτ lam ^ tdist a b)
      = if a = b then 1 else 0 := by
  have key := transfer_quadratic hd hl
  have hcon := contact_identity_mode hd hl
  have hu : dτ⁻¹ * dτ = 1 := inv_mul_cancel₀ hd.ne'
  set r := transferEig dτ lam with hr_def
  set c := equalTimeVar dτ lam with hc_def
  rcases lt_trichotomy a b with hab | rfl | hab
  · obtain ⟨m, rfl⟩ : ∃ m, b = a + m + 1 := ⟨b - a - 1, by omega⟩
    have t1 : tdist a (a + m + 1) = m + 1 := by unfold tdist; split_ifs <;> omega
    have t2 : tdist (a + 1) (a + m + 1) = m := by unfold tdist; split_ifs <;> omega
    rw [t1, t2, if_pos (by omega : a < n), if_neg (by omega : ¬ a = n),
      if_neg (by omega : ¬ a = a + m + 1)]
    rcases Nat.eq_zero_or_pos a with ha0 | ha0
    · rw [if_neg (by omega : ¬ 0 < a), if_pos ha0]
      linear_combination (dτ⁻¹ * c * r ^ m) * key - (lam * c * r ^ (m + 1) * dτ) * hu
    · have t3 : tdist (a - 1) (a + m + 1) = m + 2 := by unfold tdist; split_ifs <;> omega
      rw [if_pos ha0, if_neg (by omega : ¬ a = 0), t3]
      linear_combination (dτ⁻¹ * c * r ^ m) * key - (lam * c * r ^ (m + 1) * dτ) * hu
  · have t0 : tdist a a = 0 := by unfold tdist; split_ifs <;> omega
    rw [t0, if_pos rfl, pow_zero, mul_one]
    rcases Nat.eq_zero_or_pos a with ha0 | ha0
    · rw [if_neg (by omega : ¬ 0 < a), if_pos ha0]
      rcases lt_or_eq_of_le ha with han | han
      · have t2 : tdist (a + 1) a = 1 := by unfold tdist; split_ifs <;> omega
        rw [if_pos han, if_neg (by omega : ¬ a = n), t2, pow_one]
        linear_combination hcon
      · rw [if_neg (by omega : ¬ a < n), if_pos han]
        linear_combination hcon
    · have t3 : tdist (a - 1) a = 1 := by unfold tdist; split_ifs <;> omega
      rw [if_pos ha0, if_neg (by omega : ¬ a = 0), t3, pow_one]
      rcases lt_or_eq_of_le ha with han | han
      · have t2 : tdist (a + 1) a = 1 := by unfold tdist; split_ifs <;> omega
        rw [if_pos han, if_neg (by omega : ¬ a = n), t2, pow_one]
        linear_combination hcon
      · rw [if_neg (by omega : ¬ a < n), if_pos han]
        linear_combination hcon
  · obtain ⟨m, rfl⟩ : ∃ m, a = b + m + 1 := ⟨a - b - 1, by omega⟩
    have t1 : tdist (b + m + 1) b = m + 1 := by unfold tdist; split_ifs <;> omega
    have t3 : tdist (b + m + 1 - 1) b = m := by unfold tdist; split_ifs <;> omega
    rw [t1, t3, if_pos (by omega : 0 < b + m + 1), if_neg (by omega : ¬ b + m + 1 = 0),
      if_neg (by omega : ¬ b + m + 1 = b)]
    rcases lt_or_eq_of_le ha with han | han
    · have t2 : tdist (b + m + 1 + 1) b = m + 2 := by unfold tdist; split_ifs <;> omega
      rw [if_pos han, if_neg (by omega : ¬ b + m + 1 = n), t2]
      linear_combination (dτ⁻¹ * c * r ^ m) * key - (lam * c * r ^ (m + 1) * dτ) * hu
    · rw [if_neg (by omega : ¬ b + m + 1 < n), if_pos han]
      linear_combination (dτ⁻¹ * c * r ^ m) * key - (lam * c * r ^ (m + 1) * dτ) * hu

/-- One mode of the corrected window precision. -/
noncomputable def windowPrecMode (dτ x : ℝ) (T : ℕ) (t u : Fin (T + 1)) : ℝ :=
  dτ * (if t = u then 1 else 0) * x + dτ⁻¹ * pathLap T t u
    + dτ⁻¹ * endWeight T t * (if t = u then 1 else 0) * (1 - transferEig dτ x)

/-- For one mode, the corrected window precision inverts the Toeplitz kernel `c0 r^|u−s|`. -/
theorem window_mode (hd : 0 < dτ) (hl : 0 < lam) (T : ℕ) (t s : Fin (T + 1)) :
    ∑ u : Fin (T + 1), windowPrecMode dτ lam T t u
        * (equalTimeVar dτ lam * transferEig dτ lam ^ tdist u s)
      = if t = s then 1 else 0 := by
  have hrow := window_row_mode hd hl (n := T) (a := t) (b := s) (Nat.lt_succ_iff.mp t.isLt)
    (Nat.lt_succ_iff.mp s.isLt)
  have hL := pathLap_sum_apply T (fun n => equalTimeVar dτ lam * transferEig dτ lam ^ tdist n s) t
  have e : (if t = s then (1 : ℝ) else 0) = if (t : ℕ) = (s : ℕ) then 1 else 0 := by
    by_cases h : t = s
    · rw [if_pos h, if_pos (congrArg Fin.val h)]
    · rw [if_neg h, if_neg (fun h' => h (Fin.ext h'))]
  rw [e]
  simp only [windowPrecMode, add_mul, Finset.sum_add_distrib]
  have s1 : ∑ u : Fin (T + 1), dτ * (if t = u then (1 : ℝ) else 0) * lam
      * (equalTimeVar dτ lam * transferEig dτ lam ^ tdist u s)
      = dτ * lam * (equalTimeVar dτ lam * transferEig dτ lam ^ tdist t s) := by
    rw [Finset.sum_eq_single t]
    · rw [if_pos rfl]
      ring
    · intro u _ hu
      rw [if_neg (Ne.symm hu)]
      ring
    · intro hn
      exact absurd (Finset.mem_univ t) hn
  have s2 : ∑ u : Fin (T + 1), dτ⁻¹ * pathLap T t u
      * (equalTimeVar dτ lam * transferEig dτ lam ^ tdist u s)
      = dτ⁻¹ * ∑ u : Fin (T + 1), pathLap T t u
          * (equalTimeVar dτ lam * transferEig dτ lam ^ tdist u s) := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun u _ => by ring
  have s3 : ∑ u : Fin (T + 1), dτ⁻¹ * endWeight T t * (if t = u then (1 : ℝ) else 0)
      * (1 - transferEig dτ lam) * (equalTimeVar dτ lam * transferEig dτ lam ^ tdist u s)
      = dτ⁻¹ * endWeight T t * (1 - transferEig dτ lam)
          * (equalTimeVar dτ lam * transferEig dτ lam ^ tdist t s) := by
    rw [Finset.sum_eq_single t]
    · rw [if_pos rfl]
      ring
    · intro u _ hu
      rw [if_neg (Ne.symm hu)]
      ring
    · intro hn
      exact absurd (Finset.mem_univ t) hn
  rw [s1, s2, s3, hL]
  unfold endWeight
  linear_combination hrow

end Mode

section Matrix

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {T : ℕ}

/-- A matrix over `ι × slabs` given by its slab blocks. -/
def blockOf (F : Fin (T + 1) → Fin (T + 1) → Matrix ι ι ℝ) :
    Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ :=
  Matrix.of fun p q => F p.2 q.2 p.1 q.1

omit [DecidableEq ι] in
theorem blockOf_mul (F G : Fin (T + 1) → Fin (T + 1) → Matrix ι ι ℝ) :
    blockOf F * blockOf G = blockOf fun t s => ∑ u, F t u * G u s := by
  ext ⟨i, t⟩ ⟨k, s⟩
  simp only [blockOf, Matrix.mul_apply, Matrix.of_apply, Fintype.sum_prod_type, Matrix.sum_apply]
  exact Finset.sum_comm

omit [Fintype ι] in
theorem blockOf_delta :
    blockOf (T := T) (fun t s => if t = s then (1 : Matrix ι ι ℝ) else 0) = 1 := by
  ext ⟨i, t⟩ ⟨k, s⟩
  simp only [blockOf, Matrix.of_apply]
  by_cases h : t = s
  · subst h
    rw [if_pos rfl, Matrix.one_apply, Matrix.one_apply]
    by_cases hik : i = k
    · subst hik
      simp
    · simp [hik]
  · rw [if_neg h, Matrix.zero_apply, Matrix.one_apply_ne]
    intro hc
    exact h (congrArg Prod.snd hc)

section Calculus

variable {A : Matrix ι ι ℝ} (hA : A.IsHermitian)

theorem cfc_const (c : ℝ) : hA.cfc (fun _ => c) = c • (1 : Matrix ι ι ℝ) := by
  rw [← cfc_one hA, cfc_smul]
  exact cfc_congr hA fun _ => (mul_one c).symm

theorem cfc_sum {κ : Type*} [DecidableEq κ] (s : Finset κ) (f : κ → ℝ → ℝ) :
    ∑ u ∈ s, hA.cfc (f u) = hA.cfc (fun x => ∑ u ∈ s, f u x) := by
  refine Finset.induction_on s ?_ ?_
  · simp only [Finset.sum_empty]
    rw [cfc_const hA 0, zero_smul]
  · intro a s ha ih
    rw [Finset.sum_insert ha, ih, cfc_add]
    exact cfc_congr hA fun i => by rw [Finset.sum_insert ha]

end Calculus

variable {A : Matrix ι ι ℝ} (hA : A.PosDef) {dτ : ℝ}

/-- End-slab correction: `X` added at each end slab of the window. -/
def endSlabs (T : ℕ) (X : Matrix ι ι ℝ) : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ :=
  X ⊗ₖ diagonal (endWeight T)

/-- Covariance of the stationary chain on a window of `T + 1` slabs: blocks `C0 R^|u−s|`. -/
noncomputable def windowCov (dτ : ℝ) (T : ℕ) : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ :=
  blockOf fun u s => equalTimeCov hA dτ * transferMat hA dτ ^ tdist u s

/-- Register precision with the end-slab correction `(1 − R)/Δτ`. -/
noncomputable def windowPrec (dτ : ℝ) (T : ℕ) : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ :=
  historyPrecision 1 A dτ T + endSlabs T (dτ⁻¹ • (1 - transferMat hA dτ))

theorem cfc_window_affine (α β γ : ℝ) :
    hA.1.cfc (fun x => α * x + β + γ * (1 - transferEig dτ x))
      = α • A + β • (1 : Matrix ι ι ℝ) + γ • (1 - transferMat hA dτ) := by
  rw [transferMat, smul_eq_cfc hA.1, ← cfc_one hA.1, cfc_sub, cfc_smul, cfc_smul, cfc_add,
    cfc_add]
  exact cfc_congr hA.1 fun i => by ring

theorem windowPrec_blocks (T : ℕ) :
    windowPrec hA dτ T = blockOf fun t u => hA.1.cfc fun x => windowPrecMode dτ x T t u := by
  ext ⟨i, t⟩ ⟨k, u⟩
  have hc : hA.1.cfc (fun x => windowPrecMode dτ x T t u)
      = (dτ * (if t = u then 1 else 0)) • A + (dτ⁻¹ * pathLap T t u) • (1 : Matrix ι ι ℝ)
        + (dτ⁻¹ * endWeight T t * (if t = u then 1 else 0)) • (1 - transferMat hA dτ) := by
    rw [← cfc_window_affine hA]
    exact cfc_congr hA.1 fun i => by unfold windowPrecMode; ring
  simp only [windowPrec, endSlabs, blockOf, Matrix.add_apply, historyPrecision_apply,
    kroneckerMap_apply, Matrix.of_apply, hc, Matrix.smul_apply, smul_eq_mul]
  by_cases h : t = u
  · subst h
    rw [Matrix.one_apply_eq, diagonal_apply_eq, if_pos rfl]
    ring
  · rw [Matrix.one_apply_ne h, diagonal_apply_ne _ h, if_neg h]
    ring

theorem windowCov_blocks (T : ℕ) :
    windowCov hA dτ T = blockOf fun u s => hA.1.cfc fun x =>
      equalTimeVar dτ x * transferEig dτ x ^ tdist u s := by
  unfold windowCov
  congr 1
  funext u s
  rw [equalTimeCov, transferMat, cfc_pow, cfc_mul]

/-- The corrected register precision is the inverse of the stationary window covariance. -/
theorem stationary_window (hd : 0 < dτ) (T : ℕ) : windowPrec hA dτ T * windowCov hA dτ T = 1 := by
  rw [windowPrec_blocks, windowCov_blocks, blockOf_mul, ← blockOf_delta]
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
  exact cfc_congr hA.1 fun i => window_mode hd (eig_pos hA i) T t s

theorem windowCov_mul_windowPrec (hd : 0 < dτ) (T : ℕ) :
    windowCov hA dτ T * windowPrec hA dτ T = 1 :=
  mul_eq_one_comm.1 (stationary_window hA hd T)

/-- Among terms `D` added equally at both end slabs, the one that makes the window exact is unique:
`D = (1 − R)/Δτ`. -/
theorem window_correction_unique (hd : 0 < dτ) (T : ℕ) (D : Matrix ι ι ℝ)
    (h : (historyPrecision 1 A dτ T + endSlabs T D) * windowCov hA dτ T = 1) :
    D = dτ⁻¹ • (1 - transferMat hA dτ) := by
  have hc := windowCov_mul_windowPrec hA hd T
  have heq : historyPrecision 1 A dτ T + endSlabs T D = windowPrec hA dτ T := by
    calc historyPrecision 1 A dτ T + endSlabs T D
        = (historyPrecision 1 A dτ T + endSlabs T D) * (windowCov hA dτ T * windowPrec hA dτ T) := by
          rw [hc, Matrix.mul_one]
      _ = (historyPrecision 1 A dτ T + endSlabs T D) * windowCov hA dτ T * windowPrec hA dτ T := by
          rw [Matrix.mul_assoc]
      _ = windowPrec hA dτ T := by rw [h, Matrix.one_mul]
  have hE : endSlabs T D = endSlabs T (dτ⁻¹ • (1 - transferMat hA dτ)) := by
    unfold windowPrec at heq
    exact add_left_cancel heq
  ext i k
  have h0 := congrFun (congrFun hE (i, 0)) (k, 0)
  simp only [endSlabs, kroneckerMap_apply, diagonal_apply_eq] at h0
  exact mul_right_cancel₀ (endWeight_zero_ne T) h0

/-- If the whitened stiffness is not diagonal, the end-slab correction is not node-local. -/
theorem window_boundary_not_diag (hd : 0 < dτ) (hA' : ¬ A.IsDiag) :
    ¬ (dτ⁻¹ • (1 - transferMat hA dτ)).IsDiag := by
  intro hD
  apply hA'
  have h1 : (1 - transferMat hA dτ).IsDiag := by
    have h := hD.smul dτ
    rwa [smul_smul, mul_inv_cancel₀ hd.ne', one_smul] at h
  have hR : (transferMat hA dτ).IsDiag := by
    have h := isDiag_one.sub h1
    rwa [sub_sub_cancel] at h
  have hRi : (transferMat hA dτ)⁻¹.IsDiag := by
    rw [← (isDiag_iff_diagonal_diag _).1 hR, inv_diagonal]
    exact isDiag_diagonal _
  have hsum := transferMat_add_inv hA hd
  have hA2 : (dτ ^ 2 • A).IsDiag := by
    have h3 : dτ ^ 2 • A
        = transferMat hA dτ + (transferMat hA dτ)⁻¹ - (2 : ℝ) • (1 : Matrix ι ι ℝ) := by
      rw [hsum]
      abel
    rw [h3]
    exact (hR.add hRi).sub (isDiag_one.smul _)
  have h4 := hA2.smul (dτ ^ 2)⁻¹
  rwa [smul_smul, inv_mul_cancel₀ (pow_ne_zero 2 hd.ne'), one_smul] at h4

/-- When the whitened stiffness couples distinct nodes, adding a node-local term `D` at both end
slabs never makes the finite register precision the inverse of the stationary window covariance. -/
theorem no_nodeLocal_window_correction (hd : 0 < dτ) (T : ℕ) (hA' : ¬ A.IsDiag)
    (D : Matrix ι ι ℝ) (hD : D.IsDiag) :
    (historyPrecision 1 A dτ T + endSlabs T D) * windowCov hA dτ T ≠ 1 := by
  intro h
  have hD' := window_correction_unique hA hd T D h
  rw [hD'] at hD
  exact window_boundary_not_diag hA hd hA' hD

theorem pathLap_posSemidef (T : ℕ) : (pathLap (𝕜 := ℝ) T).PosSemidef := by
  refine Matrix.PosSemidef.of_dotProduct_mulVec_nonneg ?_ fun x => ?_
  · ext a b
    simpa [Matrix.conjTranspose_apply] using
      (congrFun (congrFun (stiffnessMatrix_transpose (fun u : Fin T => u.castSucc)
        (fun u : Fin T => u.succ) (fun _ => (1 : ℝ)) (fun _ => 0)) b) a).symm
  · rw [pathLap, stiffnessMatrix_mulVec]
    have hq := stiffness_quadratic (fun u : Fin T => u.castSucc) (fun u : Fin T => u.succ)
      (fun _ => (1 : ℝ)) (fun _ => 0) x
    simp only [star_trivial, dotProduct] at *
    rw [hq]
    apply add_nonneg <;> exact Finset.sum_nonneg fun _ _ => by positivity

/-- The corrected precision is positive definite. -/
theorem windowPrec_posDef (hd : 0 < dτ) (T : ℕ) : (windowPrec hA dτ T).PosDef := by
  unfold windowPrec historyPrecision endSlabs
  have h1 : (dτ • (A ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))).PosDef :=
    (hA.kronecker Matrix.PosDef.one).smul hd
  have h2 : (dτ⁻¹ • ((1 : Matrix ι ι ℝ) ⊗ₖ pathLap (𝕜 := ℝ) T)).PosSemidef :=
    (Matrix.PosSemidef.one.kronecker (pathLap_posSemidef T)).smul (inv_nonneg.2 hd.le)
  have h3 : ((dτ⁻¹ • (1 - transferMat hA dτ)) ⊗ₖ diagonal (endWeight T)).PosSemidef :=
    (((one_sub_transferMat_posDef hA hd).posSemidef).smul (inv_nonneg.2 hd.le)).kronecker
      (Matrix.PosSemidef.diagonal fun t => endWeight_nonneg T t)
  exact (h1.add_posSemidef h2).add_posSemidef h3

/-- The stationary covariance `ℏ P_c⁻¹` of the corrected law (precision `P_c`, the register
precision plus the end-slab term) is `ℏ Σ`. -/
theorem corrected_law_stationary (hd : 0 < dτ) (T : ℕ) (hbar : ℝ) :
    stationaryCov hbar (windowPrec hA dτ T) = hbar • windowCov hA dτ T := by
  rw [stationaryCov, Matrix.inv_eq_right_inv (stationary_window hA hd T)]

/-- Every solution of the Lyapunov equation of the declared Langevin law (whitened coordinates,
unit capacities) with the corrected precision is the stationary window covariance `ℏ Σ`. -/
theorem corrected_law_stationary_unique (hd : 0 < dτ) (T : ℕ) (hbar : ℝ)
    (X : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ)
    (hX : driftMatrix (fun _ => (1 : ℝ)) (windowPrec hA dτ T) * X
        + X * (driftMatrix (fun _ => (1 : ℝ)) (windowPrec hA dτ T))ᵀ
      = (2 * hbar) • diagonal (fun _ => (1 : ℝ)⁻¹)) :
    X = hbar • windowCov hA dτ T := by
  rw [lyapunov_unique (fun _ => 1) _ X hbar (windowPrec_posDef hA hd T) (fun _ => one_pos) hX,
    corrected_law_stationary hA hd T hbar]

end Matrix

section Original

variable {ι E : Type*} [Fintype ι] [DecidableEq ι] [Fintype E] {T : ℕ}

omit [DecidableEq ι] in
/-- Conjugating the history precision slab-wise by a matrix `S` (a change of coordinates when `S`
is symmetric). -/
theorem historyPrecision_conj (S J K : Matrix ι ι ℝ) (dτ : ℝ) (T : ℕ) :
    (S ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) * historyPrecision J K dτ T
        * (S ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
      = historyPrecision (S * J * S) (S * K * S) dτ T := by
  simp only [historyPrecision, Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.smul_mul,
    ← mul_kronecker_mul, Matrix.mul_one, Matrix.one_mul]

omit [DecidableEq ι] in
theorem endSlabs_conj (S X : Matrix ι ι ℝ) (T : ℕ) :
    (S ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) * endSlabs T X
        * (S ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
      = endSlabs T (S * X * S) := by
  simp only [endSlabs, ← mul_kronecker_mul, Matrix.mul_one, Matrix.one_mul]

/-- In original coordinates: the stiffness of the Euclidean history register with time coupling
equal to the capacity, plus the end-slab correction, is the inverse of the stationary window
covariance. -/
theorem declared_register_window (src tgt : E → ι) (c : E → ℝ) (g M : ι → ℝ) (hM : ∀ i, 0 < M i)
    (hA : (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M).PosDef) {dτ : ℝ} (hd : 0 < dτ)
    (T : ℕ) :
    (stiffnessMatrix (regSrc src T) (regTgt tgt T) (regCond c M dτ T) (regGround g dτ T)
        + endSlabs T (jSqrt M * (dτ⁻¹ • (1 - transferMat hA dτ)) * jSqrt M))
      * ((jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) * windowCov hA dτ T
          * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))) = 1 := by
  have h1 : jSqrt M * jInvSqrt M = 1 := jSqrt_mul_jInvSqrt hM
  have hK1 : (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
      * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) = 1 := by
    rw [← mul_kronecker_mul, h1, Matrix.mul_one, one_kronecker_one]
  have hV : stiffnessMatrix src tgt c g
      = jSqrt M * (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M) * jSqrt M := by
    calc stiffnessMatrix src tgt c g
        = (jSqrt M * jInvSqrt M) * stiffnessMatrix src tgt c g
            * (jInvSqrt M * jSqrt M) := by
          rw [h1, jInvSqrt_mul_jSqrt hM, Matrix.one_mul, Matrix.mul_one]
      _ = jSqrt M * (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M) * jSqrt M := by
          simp only [Matrix.mul_assoc]
  have hreg : stiffnessMatrix (regSrc src T) (regTgt tgt T) (regCond c M dτ T) (regGround g dτ T)
      = (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
        * historyPrecision 1 (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M) dτ T
        * (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) := by
    rw [historyPrecision_conj, register_eq_historyPrecision, Matrix.mul_one,
      jSqrt_mul_jSqrt hM, ← hV]
  rw [hreg, ← endSlabs_conj, ← Matrix.add_mul, ← Matrix.mul_add]
  have hw := stationary_window hA hd T
  unfold windowPrec at hw
  calc (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
        * (historyPrecision 1 (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M) dτ T
          + endSlabs T (dτ⁻¹ • (1 - transferMat hA dτ)))
        * (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
        * ((jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) * windowCov hA dτ T
          * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)))
      = (jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
        * ((historyPrecision 1 (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M) dτ T
          + endSlabs T (dτ⁻¹ • (1 - transferMat hA dτ)))
          * ((jSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
            * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)))
          * windowCov hA dτ T)
        * (jInvSqrt M ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ)) := by
        simp only [Matrix.mul_assoc]
    _ = 1 := by rw [hK1, Matrix.mul_one, hw, Matrix.mul_one, hK1]

end Original

#print axioms pathLap_sum_apply
#print axioms window_row_mode
#print axioms window_mode
#print axioms blockOf_mul
#print axioms stationary_window
#print axioms window_correction_unique
#print axioms window_boundary_not_diag
#print axioms no_nodeLocal_window_correction
#print axioms windowPrec_posDef
#print axioms corrected_law_stationary
#print axioms corrected_law_stationary_unique
#print axioms historyPrecision_conj
#print axioms declared_register_window

end OPH.SourceLaw
