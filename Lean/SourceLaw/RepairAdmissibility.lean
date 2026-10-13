import SourceLaw.CounterfamilyLocality
import SourceLaw.MatrixSeams
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.Topology.Order.DenselyOrdered
import Mathlib.Topology.Order.OrderClosed
import Mathlib.Topology.Algebra.Ring.Real

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Repair admissibility of a general time coupling

The Euclidean history register couples consecutive slabs through same-node time seams. Here the
time coupling is an arbitrary matrix `J`: the history precision is `Δτ K ⊗ 1 + Δτ⁻¹ J ⊗ L`,
with `K` the spatial stiffness and `L` the free-end path Laplacian of the slab index. A precision
is a *repair stiffness* when every off-diagonal entry is nonpositive. For a symmetric matrix this
is exactly realizability as the stiffness of agreement seams of nonnegative conductance together
with grounding by the row sums, whose sign is not constrained
(`isRepairStiffness_iff_matCond_nonneg`, `repairStiffness_realized`).

Main results:

* `register_eq_historyPrecision`: for a node-local coupling `J = diag j` the general precision is
  the stiffness of the Euclidean history register.
* `historyPrecision_repair_iff`: at a fixed slab thickness `Δτ > 0`, with at least two time seams
  per node (`T ≥ 2`), the precision is a repair stiffness exactly when every entry of `J` is
  nonnegative and `2 J_ik ≤ −Δτ² K_ik` off the diagonal.
* `fixed_step_admits_nondiagonal`: with `T ≥ 2`, at every fixed `Δτ > 0` a nondiagonal coupling
  on two nodes is admissible for some repair-stiffness `K`, so the small-`Δτ` hypothesis below
  cannot be dropped.
* `eventually_repair_isDiag`: with at least one time seam per node (`T > 0`), if the precision is
  a repair stiffness for every small enough `Δτ > 0`, then `J` is diagonal. Node locality of the
  time coupling follows from admissibility at all small `Δτ`.
* `counterfamily_eventually_repair_iff`: if `T > 0`, the capacities are positive, `M + εV` is
  invertible and `V`, which is also the spatial stiffness, couples two distinct nodes with
  nonpositive off-diagonal entries, the counterfamily inertia `M (M + εV)⁻¹ M` is admissible
  for all small `Δτ` exactly when `ε = 0`.
* `eventually_repair_clock_iff`: if `T > 0`, the capacities are positive and `K` is a repair
  stiffness, admissibility for all small `Δτ` together with the uniform time seam rate `2/Δτ`
  holds exactly when `J = diag M`.
-/

namespace OPH.SourceLaw

open scoped BigOperators Kronecker
open Matrix Filter Topology

section Path

variable {𝕜 : Type*} [Field 𝕜]

/-- Free-end path Laplacian of the slab index: unit seams between `t` and `t + 1`, no grounding. -/
def pathLap (T : ℕ) : Matrix (Fin (T + 1)) (Fin (T + 1)) 𝕜 :=
  stiffnessMatrix (fun u : Fin T => u.castSucc) (fun u : Fin T => u.succ) (fun _ => (1 : 𝕜))
    (fun _ => 0)

/-- Consecutive slabs are coupled with weight `-1`. -/
theorem pathLap_adj (T : ℕ) (u : Fin T) : pathLap (𝕜 := 𝕜) T u.castSucc u.succ = -1 := by
  simp only [pathLap, stiffnessMatrix, Matrix.of_apply, stiffness, one_mul, zero_mul, add_zero]
  rw [Finset.sum_eq_single u]
  · simp [Fin.ext_iff]
  · intro v _ hv
    have hv' : (v : ℕ) ≠ u := fun h => hv (Fin.ext h)
    simp only [Pi.single_apply, Fin.ext_iff, Fin.val_castSucc, Fin.val_succ]
    split_ifs <;> first | (exfalso; omega) | simp
  · intro h
    exact absurd (Finset.mem_univ u) h

theorem sum_castSucc_indicator (T : ℕ) (t : Fin (T + 1)) :
    (∑ v : Fin T, if t = v.castSucc then (1 : 𝕜) else 0) = if (t : ℕ) < T then 1 else 0 := by
  split_ifs with h
  · rw [Finset.sum_eq_single ⟨t, h⟩]
    · simp
    · intro v _ hv
      rw [if_neg]
      intro hc
      have := congrArg Fin.val hc
      simp only [Fin.val_castSucc] at this
      exact hv (Fin.ext (by show (v : ℕ) = t; omega))
    · intro hn
      exact absurd (Finset.mem_univ _) hn
  · refine Finset.sum_eq_zero fun v _ => ?_
    rw [if_neg]
    intro hc
    have := congrArg Fin.val hc
    simp only [Fin.val_castSucc] at this
    have := v.isLt
    omega

theorem sum_succ_indicator (T : ℕ) (t : Fin (T + 1)) :
    (∑ v : Fin T, if t = v.succ then (1 : 𝕜) else 0) = if 0 < (t : ℕ) then 1 else 0 := by
  have ht := t.isLt
  split_ifs with h
  · rw [Finset.sum_eq_single ⟨(t : ℕ) - 1, by omega⟩]
    · rw [if_pos]
      apply Fin.ext
      simp only [Fin.val_succ]
      omega
    · intro v _ hv
      rw [if_neg]
      intro hc
      have := congrArg Fin.val hc
      simp only [Fin.val_succ] at this
      exact hv (Fin.ext (by show (v : ℕ) = (t : ℕ) - 1; omega))
    · intro hn
      exact absurd (Finset.mem_univ _) hn
  · refine Finset.sum_eq_zero fun v _ => ?_
    rw [if_neg]
    intro hc
    have := congrArg Fin.val hc
    simp only [Fin.val_succ] at this
    omega

/-- The diagonal of the path Laplacian counts the time neighbours of a slab. -/
theorem pathLap_diag (T : ℕ) (t : Fin (T + 1)) :
    pathLap (𝕜 := 𝕜) T t t = (if (t : ℕ) < T then 1 else 0) + (if 0 < (t : ℕ) then 1 else 0) := by
  simp only [pathLap, stiffnessMatrix, Matrix.of_apply, stiffness, one_mul, zero_mul, add_zero]
  have h1 : ∀ v : Fin T, (if t = v.castSucc then
      (Pi.single t (1 : 𝕜) : Fin (T + 1) → 𝕜) v.castSucc
        - (Pi.single t (1 : 𝕜) : Fin (T + 1) → 𝕜) v.succ else 0)
        = if t = v.castSucc then 1 else 0 := by
    intro v
    split_ifs with h
    · subst h
      simp [Fin.ext_iff]
    · rfl
  have h2 : ∀ v : Fin T, (if t = v.succ then
      (Pi.single t (1 : 𝕜) : Fin (T + 1) → 𝕜) v.succ
        - (Pi.single t (1 : 𝕜) : Fin (T + 1) → 𝕜) v.castSucc else 0)
        = if t = v.succ then 1 else 0 := by
    intro v
    split_ifs with h
    · subst h
      simp [Fin.ext_iff]
    · rfl
  rw [Finset.sum_add_distrib, Finset.sum_congr rfl fun v _ => h1 v,
    Finset.sum_congr rfl fun v _ => h2 v, sum_castSucc_indicator, sum_succ_indicator]

end Path

section General

variable {𝕜 : Type*} [Field 𝕜] {ι : Type*}

/-- History precision of a general time coupling `J`: slab blocks `Δτ K` and the time Laplacian
weighted by `J / Δτ`. -/
def historyPrecision (J K : Matrix ι ι 𝕜) (dτ : 𝕜) (T : ℕ) :
    Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) 𝕜 :=
  dτ • (K ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) 𝕜)) + dτ⁻¹ • (J ⊗ₖ pathLap T)

theorem historyPrecision_apply (J K : Matrix ι ι 𝕜) (dτ : 𝕜) (T : ℕ) (i k : ι)
    (t s : Fin (T + 1)) :
    historyPrecision J K dτ T (i, t) (k, s)
      = dτ * (K i k * (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) 𝕜) t s)
        + dτ⁻¹ * (J i k * pathLap T t s) := by
  simp only [historyPrecision, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul,
    kroneckerMap_apply]

end General

section Register

variable {𝕜 : Type*} [Field 𝕜] {ι E : Type*} [DecidableEq ι] [Fintype ι] [Fintype E]

omit [Fintype ι] in
theorem stiffness_zero (src tgt : E → ι) (c : E → 𝕜) (g : ι → 𝕜) (k : ι) :
    stiffness src tgt c g 0 k = 0 := by
  simp [stiffness]

/-- The register stiffness at `(i, t)` is `Δτ` times the spatial stiffness of slab `t` plus
`J_i / Δτ` times the path stiffness of the time series of node `i`. -/
theorem register_stiffness_slice (src tgt : E → ι) (c : E → 𝕜) (g j : ι → 𝕜) (dτ : 𝕜) (T : ℕ)
    (f : ι × Fin (T + 1) → 𝕜) (i : ι) (t : Fin (T + 1)) :
    stiffness (regSrc src T) (regTgt tgt T) (regCond c j dτ T) (regGround g dτ T) f (i, t)
      = dτ * stiffness src tgt c g (fun x => f (x, t)) i
        + j i / dτ * stiffness (fun u : Fin T => u.castSucc) (fun u : Fin T => u.succ)
            (fun _ => (1 : 𝕜)) (fun _ => 0) (fun u => f (i, u)) t := by
  simp only [stiffness, Fintype.sum_sum_type, Fintype.sum_prod_type, regSrc, regTgt, regCond,
    regGround]
  have hsp : ∀ e : E, (∑ u : Fin (T + 1), dτ * c e *
      ((if (i, t) = (src e, u) then f (src e, u) - f (tgt e, u) else 0)
        + (if (i, t) = (tgt e, u) then f (tgt e, u) - f (src e, u) else 0)))
      = dτ * (c e * ((if i = src e then f (src e, t) - f (tgt e, t) else 0)
        + (if i = tgt e then f (tgt e, t) - f (src e, t) else 0))) := by
    intro e
    rw [Finset.sum_eq_single t]
    · simp only [Prod.mk.injEq, and_true]
      ring
    · intro u _ hu
      simp [Prod.mk.injEq, Ne.symm hu]
    · intro h
      exact absurd (Finset.mem_univ t) h
  have htm : ∀ m : ι, (∑ u : Fin T, j m / dτ *
      ((if (i, t) = (m, u.castSucc) then f (m, u.castSucc) - f (m, u.succ) else 0)
        + (if (i, t) = (m, u.succ) then f (m, u.succ) - f (m, u.castSucc) else 0)))
      = if m = i then j i / dτ * ∑ u : Fin T, 1 *
          ((if t = u.castSucc then f (i, u.castSucc) - f (i, u.succ) else 0)
            + (if t = u.succ then f (i, u.succ) - f (i, u.castSucc) else 0)) else 0 := by
    intro m
    by_cases hm : m = i
    · subst hm
      rw [if_pos rfl, Finset.mul_sum]
      refine Finset.sum_congr rfl fun u _ => ?_
      simp only [Prod.mk.injEq, true_and, one_mul]
    · rw [if_neg hm]
      refine Finset.sum_eq_zero fun u _ => ?_
      simp [Prod.mk.injEq, Ne.symm hm]
  rw [Finset.sum_congr rfl fun e _ => hsp e, Finset.sum_congr rfl fun m _ => htm m,
    Finset.sum_ite_eq', if_pos (Finset.mem_univ i), ← Finset.mul_sum]
  ring

/-- For a node-local coupling `J = diag j` the general history precision is the stiffness of the
Euclidean history register. -/
theorem register_eq_historyPrecision (src tgt : E → ι) (c : E → 𝕜) (g j : ι → 𝕜) (dτ : 𝕜)
    (T : ℕ) :
    stiffnessMatrix (regSrc src T) (regTgt tgt T) (regCond c j dτ T) (regGround g dτ T)
      = historyPrecision (diagonal j) (stiffnessMatrix src tgt c g) dτ T := by
  ext ⟨i, t⟩ ⟨k, s⟩
  rw [historyPrecision_apply]
  simp only [stiffnessMatrix, Matrix.of_apply]
  rw [register_stiffness_slice]
  have h1 : (fun x => (Pi.single (k, s) (1 : 𝕜) : ι × Fin (T + 1) → 𝕜) (x, t))
      = if t = s then (Pi.single k (1 : 𝕜) : ι → 𝕜) else 0 := by
    funext x
    by_cases hts : t = s
    · subst hts
      simp [Pi.single_apply]
    · simp [Prod.ext_iff, hts]
  have h2 : (fun u => (Pi.single (k, s) (1 : 𝕜) : ι × Fin (T + 1) → 𝕜) (i, u))
      = if i = k then (Pi.single s (1 : 𝕜) : Fin (T + 1) → 𝕜) else 0 := by
    funext u
    by_cases hik : i = k
    · subst hik
      simp [Pi.single_apply]
    · simp [Prod.ext_iff, hik]
  rw [h1, h2]
  by_cases hts : t = s
  · subst hts
    rw [if_pos rfl, Matrix.one_apply_eq]
    by_cases hik : i = k
    · subst hik
      rw [if_pos rfl, diagonal_apply_eq]
      simp only [pathLap, stiffnessMatrix, Matrix.of_apply]
      ring
    · rw [if_neg hik, stiffness_zero, diagonal_apply_ne _ hik]
      ring
  · rw [if_neg hts, stiffness_zero, Matrix.one_apply_ne hts]
    by_cases hik : i = k
    · subst hik
      rw [if_pos rfl, diagonal_apply_eq]
      simp only [pathLap, stiffnessMatrix, Matrix.of_apply]
      ring
    · rw [if_neg hik, stiffness_zero, diagonal_apply_ne _ hik]
      ring

end Register

section Admissible

variable {ι : Type*}

/-- A precision is a repair stiffness when every off-diagonal entry is nonpositive. -/
def IsRepairStiffness {κ : Type*} (P : Matrix κ κ ℝ) : Prop := ∀ a b, a ≠ b → P a b ≤ 0

/-- Seams of nonnegative conductance always give a repair stiffness. -/
theorem stiffnessMatrix_isRepairStiffness {E : Type*} [DecidableEq ι] [Fintype E]
    (src tgt : E → ι) (c : E → ℝ) (g : ι → ℝ) (hc : ∀ e, 0 ≤ c e) :
    IsRepairStiffness (stiffnessMatrix src tgt c g) := by
  intro k l hkl
  simp only [stiffnessMatrix, Matrix.of_apply, stiffness, Pi.single_apply, if_neg hkl, mul_zero,
    add_zero]
  refine Finset.sum_nonpos fun e _ => mul_nonpos_of_nonneg_of_nonpos (hc e) ?_
  split_ifs <;> first | (norm_num; done) | simp_all

/-- A matrix is a repair stiffness exactly when its matrix seam graph has nonnegative
conductances. -/
theorem isRepairStiffness_iff_matCond_nonneg (P : Matrix ι ι ℝ) :
    IsRepairStiffness P ↔ ∀ p : MatSeam ι, 0 ≤ matCond P p := by
  constructor
  · intro h p
    exact matCond_nonneg P h p
  · intro h a b hab
    have := h ⟨(a, b), hab⟩
    simp only [matCond] at this
    linarith

/-- A symmetric repair stiffness is the stiffness of its matrix seam graph, grounded by its row
sums, and every conductance of that graph is nonnegative. -/
theorem repairStiffness_realized [Fintype ι] [DecidableEq ι] (P : Matrix ι ι ℝ) (hP : Pᵀ = P)
    (h : IsRepairStiffness P) :
    stiffnessMatrix matSrc matTgt (matCond P) (rowSum P) = P ∧ ∀ p, 0 ≤ matCond P p :=
  ⟨stiffnessMatrix_matSeams P hP, (isRepairStiffness_iff_matCond_nonneg P).1 h⟩

theorem pathLap_isRepairStiffness (T : ℕ) : IsRepairStiffness (pathLap (𝕜 := ℝ) T) :=
  stiffnessMatrix_isRepairStiffness _ _ _ _ fun _ => zero_le_one

variable [DecidableEq ι] {T : ℕ}

omit [DecidableEq ι] in
theorem historyPrecision_same_slab (J K : Matrix ι ι ℝ) (dτ : ℝ) (i k : ι) (t : Fin (T + 1)) :
    historyPrecision J K dτ T (i, t) (k, t) = dτ * K i k + dτ⁻¹ * (J i k * pathLap T t t) := by
  rw [historyPrecision_apply, Matrix.one_apply_eq, mul_one]

omit [DecidableEq ι] in
theorem historyPrecision_adj (J K : Matrix ι ι ℝ) (dτ : ℝ) (i k : ι) (u : Fin T) :
    historyPrecision J K dτ T (i, u.castSucc) (k, u.succ) = -(J i k / dτ) := by
  have hne : u.castSucc ≠ u.succ := Fin.castSucc_lt_succ.ne
  rw [historyPrecision_apply, Matrix.one_apply_ne hne, pathLap_adj]
  ring

omit [DecidableEq ι] in
/-- A nonnegative coupling across consecutive slabs is forced by any admissible step. -/
theorem repair_coupling_nonneg {J K : Matrix ι ι ℝ} {dτ : ℝ} (hT : 0 < T) (hd : 0 < dτ)
    (h : IsRepairStiffness (historyPrecision J K dτ T)) (i k : ι) : 0 ≤ J i k := by
  let u0 : Fin T := ⟨0, hT⟩
  have hne : ((i, u0.castSucc) : ι × Fin (T + 1)) ≠ (k, u0.succ) := by
    intro hc
    have := congrArg (fun p => (p.2 : ℕ)) hc
    simp [u0] at this
  have h1 := h _ _ hne
  rw [historyPrecision_adj] at h1
  have h2 : 0 ≤ J i k / dτ := by linarith
  have h3 := mul_nonneg h2 hd.le
  rwa [div_mul_cancel₀ _ hd.ne'] at h3

/-- At a fixed slab thickness, with at least two time seams, the history precision is a repair
stiffness exactly when every entry of `J` is nonnegative and `2 J_ik ≤ −Δτ² K_ik` off the
diagonal. The bound is attained by an interior slab. -/
theorem historyPrecision_repair_iff (hT : 2 ≤ T) {dτ : ℝ} (hd : 0 < dτ) (J K : Matrix ι ι ℝ) :
    IsRepairStiffness (historyPrecision J K dτ T)
      ↔ (∀ i, 0 ≤ J i i) ∧ ∀ i k, i ≠ k → 0 ≤ J i k ∧ 2 * J i k ≤ -(dτ ^ 2 * K i k) := by
  have hT1 : 0 < T := by omega
  let t1 : Fin (T + 1) := ⟨1, by omega⟩
  have hdeg : pathLap (𝕜 := ℝ) T t1 t1 = 2 := by
    rw [pathLap_diag]
    have h1 : (t1 : ℕ) < T := by show 1 < T; omega
    have h2 : 0 < (t1 : ℕ) := by show 0 < 1; omega
    rw [if_pos h1, if_pos h2]
    norm_num
  constructor
  · intro h
    refine ⟨fun i => repair_coupling_nonneg hT1 hd h i i, fun i k hik => ?_⟩
    refine ⟨repair_coupling_nonneg hT1 hd h i k, ?_⟩
    have hne : ((i, t1) : ι × Fin (T + 1)) ≠ (k, t1) := by
      intro hc
      exact hik (congrArg Prod.fst hc)
    have h1 := h _ _ hne
    rw [historyPrecision_same_slab, hdeg] at h1
    have key : dτ * (dτ * K i k + dτ⁻¹ * (J i k * 2)) = dτ ^ 2 * K i k + 2 * J i k := by
      rw [mul_add, ← mul_assoc dτ dτ⁻¹, mul_inv_cancel₀ hd.ne', one_mul]
      ring
    have h2 := mul_nonpos_of_nonneg_of_nonpos hd.le h1
    rw [key] at h2
    linarith
  · rintro ⟨hdiag, hoff⟩ ⟨i, t⟩ ⟨k, s⟩ hne
    have hJ : ∀ a b, 0 ≤ J a b := fun a b =>
      if hab : a = b then hab ▸ hdiag a else (hoff a b hab).1
    rw [historyPrecision_apply]
    by_cases hts : t = s
    · subst hts
      have hik : i ≠ k := fun h => hne (by rw [h])
      rw [Matrix.one_apply_eq, mul_one, pathLap_diag]
      obtain ⟨h0, h2⟩ := hoff i k hik
      have hdeg2 : ((if (t : ℕ) < T then (1 : ℝ) else 0) + (if 0 < (t : ℕ) then 1 else 0)) ≤ 2 := by
        split_ifs <;> norm_num
      have e1 : J i k * ((if (t : ℕ) < T then (1 : ℝ) else 0) + (if 0 < (t : ℕ) then 1 else 0))
          ≤ J i k * 2 := mul_le_mul_of_nonneg_left hdeg2 h0
      have e2 : dτ⁻¹ * (J i k * ((if (t : ℕ) < T then (1 : ℝ) else 0)
          + (if 0 < (t : ℕ) then 1 else 0))) ≤ dτ⁻¹ * (-(dτ ^ 2 * K i k)) := by
        apply mul_le_mul_of_nonneg_left _ (inv_nonneg.mpr hd.le)
        linarith
      have e3 : dτ⁻¹ * (-(dτ ^ 2 * K i k)) = -(dτ * K i k) := by
        field_simp
      linarith
    · rw [Matrix.one_apply_ne hts, mul_zero, mul_zero, zero_add]
      exact mul_nonpos_of_nonneg_of_nonpos (inv_nonneg.mpr hd.le)
        (mul_nonpos_of_nonneg_of_nonpos (hJ i k) (pathLap_isRepairStiffness T t s hts))

/-- A node-local nonnegative coupling over a repair stiffness is admissible at every step. -/
theorem historyPrecision_diagonal_repair {K : Matrix ι ι ℝ} {j : ι → ℝ} {dτ : ℝ}
    (hK : IsRepairStiffness K) (hj : ∀ i, 0 ≤ j i) (hd : 0 < dτ) :
    IsRepairStiffness (historyPrecision (diagonal j) K dτ T) := by
  rintro ⟨i, t⟩ ⟨k, s⟩ hne
  rw [historyPrecision_apply]
  by_cases hts : t = s
  · subst hts
    have hik : i ≠ k := fun h => hne (by rw [h])
    rw [Matrix.one_apply_eq, mul_one, diagonal_apply_ne _ hik, zero_mul, mul_zero, add_zero]
    exact mul_nonpos_of_nonneg_of_nonpos hd.le (hK i k hik)
  · rw [Matrix.one_apply_ne hts, mul_zero, mul_zero, zero_add]
    refine mul_nonpos_of_nonneg_of_nonpos (inv_nonneg.mpr hd.le) ?_
    by_cases hik : i = k
    · subst hik
      rw [diagonal_apply_eq]
      exact mul_nonpos_of_nonneg_of_nonpos (hj i) (pathLap_isRepairStiffness T t s hts)
    · simp [diagonal_apply_ne _ hik]

omit [DecidableEq ι] in
/-- With at least one time seam per node (`T > 0`), admissibility at every small enough step forces
a node-local time coupling. -/
theorem eventually_repair_isDiag (hT : 0 < T) (J K : Matrix ι ι ℝ)
    (h : ∀ᶠ dτ in 𝓝[>] (0 : ℝ), IsRepairStiffness (historyPrecision J K dτ T)) : J.IsDiag := by
  intro i k hik
  have hpos : ∀ᶠ dτ in 𝓝[>] (0 : ℝ), 0 < dτ := self_mem_nhdsWithin
  obtain ⟨dτ, hrep, hd⟩ := (h.and hpos).exists
  have h0 : 0 ≤ J i k := repair_coupling_nonneg hT hd hrep i k
  have hdeg : pathLap (𝕜 := ℝ) T 0 0 = 1 := by
    rw [pathLap_diag]
    simp [hT]
  have hup : ∀ᶠ dτ in 𝓝[>] (0 : ℝ), J i k ≤ -(dτ ^ 2 * K i k) := by
    filter_upwards [h, hpos] with dτ hrep hd
    have hne : ((i, 0) : ι × Fin (T + 1)) ≠ (k, 0) := by
      intro hc
      exact hik (congrArg Prod.fst hc)
    have h1 := hrep _ _ hne
    rw [historyPrecision_same_slab, hdeg, mul_one] at h1
    have key : dτ * (dτ * K i k + dτ⁻¹ * J i k) = dτ ^ 2 * K i k + J i k := by
      rw [mul_add, ← mul_assoc dτ dτ⁻¹, mul_inv_cancel₀ hd.ne', one_mul]
      ring
    have h2 := mul_nonpos_of_nonneg_of_nonpos hd.le h1
    rw [key] at h2
    linarith
  have hlim : Tendsto (fun dτ : ℝ => -(dτ ^ 2 * K i k)) (𝓝[>] 0) (𝓝 0) := by
    have : Tendsto (fun dτ : ℝ => -(dτ ^ 2 * K i k)) (𝓝 0) (𝓝 (-((0 : ℝ) ^ 2 * K i k))) :=
      ((continuous_pow 2).mul continuous_const).neg.tendsto 0
    simpa using this.mono_left nhdsWithin_le_nhds
  have hle : J i k ≤ 0 := ge_of_tendsto hlim hup
  linarith

/-- With at least two time seams per node (`T ≥ 2`), at every fixed step some nondiagonal coupling
on two nodes is admissible for a repair-stiffness `K`, so the small-step hypothesis of
`eventually_repair_isDiag` cannot be dropped. -/
theorem fixed_step_admits_nondiagonal (hT : 2 ≤ T) {dτ : ℝ} (hd : 0 < dτ) :
    ∃ J K : Matrix (Fin 2) (Fin 2) ℝ, ¬ J.IsDiag ∧ IsRepairStiffness K
      ∧ IsRepairStiffness (historyPrecision J K dτ T) := by
  refine ⟨!![1, dτ ^ 2 / 2; dτ ^ 2 / 2, 1], !![2, -1; -1, 2], ?_, ?_, ?_⟩
  · intro hD
    have := hD (i := 0) (j := 1) (by decide)
    simp at this
    exact absurd this (by positivity)
  · intro a b hab
    fin_cases a <;> fin_cases b <;> simp_all
  · rw [historyPrecision_repair_iff hT hd]
    refine ⟨fun i => ?_, fun i k hik => ?_⟩
    · fin_cases i <;> simp
    · fin_cases i <;> fin_cases k <;> simp_all <;> constructor <;> nlinarith [sq_nonneg dτ]

end Admissible

section Counterfamily

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {T : ℕ}

/-- If `T > 0`, the capacities are positive, `M + εV` is invertible and `V`, which is also the
spatial stiffness, couples two distinct nodes with nonpositive off-diagonal entries, the
counterfamily inertia `M (M + εV)⁻¹ M` is admissible for all small steps exactly when
`ε = 0`. -/
theorem counterfamily_eventually_repair_iff (hT : 0 < T) (M : ι → ℝ) (V : Matrix ι ι ℝ) (ε : ℝ)
    (hM : ∀ i, 0 < M i) (hU : IsUnit (diagonal M + ε • V).det)
    (hV : ∃ i j, i ≠ j ∧ V i j ≠ 0) (hVoff : IsRepairStiffness V) :
    (∀ᶠ dτ in 𝓝[>] (0 : ℝ),
        IsRepairStiffness (historyPrecision (counterfamilyInertia M V ε) V dτ T))
      ↔ ε = 0 := by
  constructor
  · intro h
    exact (counterfamilyInertia_isDiag_iff M V ε (fun i => (hM i).ne') hU hV).1
      (eventually_repair_isDiag hT _ _ h)
  · rintro rfl
    rw [counterfamilyInertia_zero M V (fun i => (hM i).ne')]
    filter_upwards [self_mem_nhdsWithin] with dτ hd
    exact historyPrecision_diagonal_repair hVoff (fun i => (hM i).le) hd

omit [Fintype ι] in
/-- For `T > 0`, positive capacities and a repair-stiffness `K`: admissibility for all small steps
together with the uniform time seam rate `2/Δτ` (the rate of the register's time seams,
`timeSeam_rate`) holds exactly when the time coupling is the capacity matrix. -/
theorem eventually_repair_clock_iff (hT : 0 < T) (M : ι → ℝ) (J K : Matrix ι ι ℝ)
    (hM : ∀ i, 0 < M i) (hK : IsRepairStiffness K) {dτ0 : ℝ} (hd0 : dτ0 ≠ 0) :
    ((∀ᶠ dτ in 𝓝[>] (0 : ℝ), IsRepairStiffness (historyPrecision J K dτ T))
        ∧ ∀ i, J i i / dτ0 * (1 / M i + 1 / M i) = 2 / dτ0)
      ↔ J = diagonal M := by
  constructor
  · rintro ⟨h, hrate⟩
    have hD := eventually_repair_isDiag hT J K h
    rw [← (isDiag_iff_diagonal_diag J).1 hD]
    congr 1
    funext i
    have hi := hrate i
    have key : J i i / dτ0 * (1 / M i + 1 / M i) = 2 / dτ0 * (J i i / M i) := by ring
    rw [key] at hi
    have h2d : (2 : ℝ) / dτ0 ≠ 0 := div_ne_zero two_ne_zero hd0
    have h1 : J i i / M i = 1 := mul_left_cancel₀ h2d (by rw [hi, mul_one])
    show J i i = M i
    rwa [div_eq_one_iff_eq (hM i).ne'] at h1
  · rintro rfl
    refine ⟨?_, fun i => ?_⟩
    · filter_upwards [self_mem_nhdsWithin] with dτ hd
      exact historyPrecision_diagonal_repair hK (fun i => (hM i).le) hd
    · rw [diagonal_apply_eq]
      have := (hM i).ne'
      field_simp
      ring

end Counterfamily

#print axioms pathLap_adj
#print axioms pathLap_diag
#print axioms register_stiffness_slice
#print axioms register_eq_historyPrecision
#print axioms stiffnessMatrix_isRepairStiffness
#print axioms isRepairStiffness_iff_matCond_nonneg
#print axioms repairStiffness_realized
#print axioms historyPrecision_repair_iff
#print axioms historyPrecision_diagonal_repair
#print axioms eventually_repair_isDiag
#print axioms fixed_step_admits_nondiagonal
#print axioms counterfamily_eventually_repair_iff
#print axioms eventually_repair_clock_iff

end OPH.SourceLaw
