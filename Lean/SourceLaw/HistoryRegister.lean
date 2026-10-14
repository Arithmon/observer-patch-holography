import SourceLaw.WeightedSeams

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The Euclidean history register

The register holds one real working reading for every node `i` and every Euclidean slab
index `t : Fin (T + 1)`. Its seams are the spatial seams of the carrier in every slab
(conductance `Δτ c_e`), the grounding in every slab (`Δτ g_i`), and the same-node time seams
between `(i, t)` and `(i, t + 1)` (conductance `J_i / Δτ`). Every time copy of node `i` carries the
capacity `M_i`. The register entries are working data, not records, and the slab index is a
declared Euclidean index, not a physical clock.

Main results:

* `timeSeam_rate_eq_iff`: a time seam has rate `2/Δτ` exactly when `J_i = M_i`.
* `uniform_time_rate_iff`: with at least one time seam per node (`T > 0`), all time seams share
  one rate `λ` exactly when `J_i = λ Δτ M_i / 2` for every node.
* `register_drift`: on a carrier without self-seams, with nonzero capacities and nonzero
  `M_i + M_j` for every pair of nodes, the mean repair generator over all register seams is
  `𝕄⁻¹ P f`, where
  `P` is the register stiffness. The declared Langevin law takes `𝕄⁻¹ P` as its drift matrix
  (`driftMatrix`), so its deterministic part is `−𝕄⁻¹ P f`; that is a definition, not part of
  the statement.
* `registerAction_eq`: half the register stiffness form is the time-difference kinetic term
  with matrix `diag J` plus the slab potentials.
* `register_posDef`: with nonnegative conductances and time coefficients, `Δτ > 0` and every node
  grounded, the register stiffness is positive definite.
-/

namespace OPH.SourceLaw

open scoped BigOperators
open Matrix

variable {𝕜 : Type*} [Field 𝕜] {ι E : Type*}

/-- Seams of the register with `T + 1` slabs: spatial seams in every slab and same-node time
seams between `(i, t)` and `(i, t + 1)`. -/
abbrev RegEdge (ι E : Type*) (T : ℕ) := (E × Fin (T + 1)) ⊕ (ι × Fin T)

/-- First endpoint of a register seam. -/
def regSrc (src : E → ι) (T : ℕ) : RegEdge ι E T → ι × Fin (T + 1)
  | .inl (e, t) => (src e, t)
  | .inr (i, t) => (i, t.castSucc)

/-- Second endpoint of a register seam. -/
def regTgt (tgt : E → ι) (T : ℕ) : RegEdge ι E T → ι × Fin (T + 1)
  | .inl (e, t) => (tgt e, t)
  | .inr (i, t) => (i, t.succ)

/-- Spatial seams carry `Δτ c_e`; time seams carry `J_i / Δτ`. -/
def regCond (c : E → 𝕜) (J : ι → 𝕜) (dτ : 𝕜) (T : ℕ) : RegEdge ι E T → 𝕜
  | .inl (e, _) => dτ * c e
  | .inr (i, _) => J i / dτ

/-- Grounding in every slab: `Δτ g_i`. -/
def regGround (g : ι → 𝕜) (dτ : 𝕜) (T : ℕ) (p : ι × Fin (T + 1)) : 𝕜 := dτ * g p.1

/-- Every time copy of node `i` carries the capacity `M_i`. -/
def regCap (M : ι → 𝕜) (T : ℕ) (p : ι × Fin (T + 1)) : 𝕜 := M p.1

theorem reg_src_ne_tgt (src tgt : E → ι) (hst : ∀ e, src e ≠ tgt e) (T : ℕ) :
    ∀ s : RegEdge ι E T, regSrc src T s ≠ regTgt tgt T s := by
  rintro (⟨e, t⟩ | ⟨i, t⟩)
  · simp [regSrc, regTgt, hst e]
  · simp only [regSrc, regTgt, ne_eq, Prod.mk.injEq, true_and]
    exact Fin.castSucc_lt_succ.ne

/-- The rate of a time seam: `(J_i / Δτ)(1/M_i + 1/M_i)`. -/
theorem timeSeam_rate (src tgt : E → ι) (c : E → 𝕜) (M J : ι → 𝕜) (dτ : 𝕜) (T : ℕ)
    (i : ι) (t : Fin T) :
    seamRate (regCap M T) (regCond c J dτ T (.inr (i, t)))
        (regSrc src T (.inr (i, t))) (regTgt tgt T (.inr (i, t)))
      = J i / dτ * (1 / M i + 1 / M i) := rfl

section Rate

variable [CharZero 𝕜]

/-- A time seam has rate `2/Δτ` exactly when its conductance coefficient equals the
capacity. -/
theorem timeSeam_rate_eq_iff (src tgt : E → ι) (c : E → 𝕜) (M J : ι → 𝕜) (dτ : 𝕜) (T : ℕ)
    (i : ι) (t : Fin T) (hM : M i ≠ 0) (hd : dτ ≠ 0) :
    seamRate (regCap M T) (regCond c J dτ T (.inr (i, t)))
        (regSrc src T (.inr (i, t))) (regTgt tgt T (.inr (i, t))) = 2 / dτ ↔ J i = M i := by
  rw [timeSeam_rate]
  have key : J i / dτ * (1 / M i + 1 / M i) = 2 / dτ * (J i / M i) := by ring
  rw [key]
  have h2d : (2 : 𝕜) / dτ ≠ 0 := div_ne_zero two_ne_zero hd
  constructor
  · intro h
    have h1 : J i / M i = 1 := mul_left_cancel₀ h2d (by rw [h, mul_one])
    rwa [div_eq_one_iff_eq hM] at h1
  · intro h
    rw [h, div_self hM, mul_one]

/-- With at least one time seam per node (`T > 0`), all time seams share the rate `λ` exactly
when `J_i = λ Δτ M_i / 2` for every node. -/
theorem uniform_time_rate_iff (src tgt : E → ι) (c : E → 𝕜) (M J : ι → 𝕜) (dτ : 𝕜) (T : ℕ)
    (hM : ∀ i, M i ≠ 0) (hd : dτ ≠ 0) (hT : 0 < T) (lam : 𝕜) :
    (∀ (i : ι) (t : Fin T), seamRate (regCap M T) (regCond c J dτ T (.inr (i, t)))
        (regSrc src T (.inr (i, t))) (regTgt tgt T (.inr (i, t))) = lam)
      ↔ ∀ i, J i = lam * dτ * M i / 2 := by
  have key : ∀ i, J i / dτ * (1 / M i + 1 / M i) = lam ↔ J i = lam * dτ * M i / 2 := by
    intro i
    have hi := hM i
    constructor
    · intro h
      rw [← h]
      field_simp
      ring
    · intro h
      rw [h]
      field_simp
      ring
  constructor
  · intro h i
    exact (key i).1 (by rw [← timeSeam_rate src tgt c M J dτ T i ⟨0, hT⟩]; exact h i ⟨0, hT⟩)
  · intro h i t
    rw [timeSeam_rate]
    exact (key i).2 (h i)

end Rate

section Drift

variable [DecidableEq ι] [Fintype ι] [Fintype E]

/-- On a carrier without self-seams (`hst`), with nonzero capacities (`hM`) and nonzero
`M_i + M_j` for every pair of nodes (`hMM`), the mean repair generator summed over every register
seam (spatial, grounded and time seams) is the register stiffness divided by the capacity,
`(P f)_k / M_k`. The declared Langevin law
takes the matrix `𝕄⁻¹ P` of this map as its drift matrix (`driftMatrix`, with
`stiffnessMatrix_mulVec`); that step is a definition, not part of this statement. -/
theorem register_drift (src tgt : E → ι) (c : E → 𝕜) (g M J : ι → 𝕜) (dτ : 𝕜) (T : ℕ)
    (f : ι × Fin (T + 1) → 𝕜) (hM : ∀ i, M i ≠ 0) (hMM : ∀ i j, M i + M j ≠ 0)
    (hst : ∀ e, src e ≠ tgt e) (k : ι × Fin (T + 1)) :
    repairGenerator (regCap M T) (regSrc src T) (regTgt tgt T) (regCond c J dτ T)
        (regGround g dτ T) f k
      = stiffness (regSrc src T) (regTgt tgt T) (regCond c J dτ T) (regGround g dτ T) f k
          / regCap M T k :=
  repairGenerator_apply _ _ _ _ _ f (fun p => hM p.1) (fun _ => hMM _ _)
    (reg_src_ne_tgt src tgt hst T) k

end Drift

section Action

variable [DecidableEq ι] [Fintype ι] [Fintype E]

/-- The register action `S f = ½ ⟨f, P f⟩` with `P` the register stiffness. -/
noncomputable def registerAction (src tgt : E → ι) (c : E → ℝ) (g J : ι → ℝ) (dτ : ℝ) (T : ℕ)
    (f : ι × Fin (T + 1) → ℝ) : ℝ :=
  (1 / 2) * ∑ p, f p * stiffness (regSrc src T) (regTgt tgt T) (regCond c J dτ T)
    (regGround g dτ T) f p

/-- Half the register stiffness form is the time-difference kinetic term with matrix `diag J`
plus the slab potentials. -/
theorem registerAction_eq (src tgt : E → ι) (c : E → ℝ) (g J : ι → ℝ) (dτ : ℝ) (T : ℕ)
    (f : ι × Fin (T + 1) → ℝ) :
    registerAction src tgt c g J dτ T f
      = ∑ t : Fin T, ∑ i, J i / dτ * (f (i, t.succ) - f (i, t.castSucc)) ^ 2 / 2
        + ∑ t : Fin (T + 1), dτ * (∑ e, c e * (f (src e, t) - f (tgt e, t)) ^ 2
            + ∑ i, g i * f (i, t) ^ 2) / 2 := by
  unfold registerAction
  rw [stiffness_quadratic, Fintype.sum_sum_type]
  simp only [regCond, regSrc, regTgt, regGround, Fintype.sum_prod_type]
  have hA : (∑ e, ∑ t : Fin (T + 1), dτ * c e * (f (src e, t) - f (tgt e, t)) ^ 2)
      = ∑ t : Fin (T + 1), dτ * ∑ e, c e * (f (src e, t) - f (tgt e, t)) ^ 2 := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun t _ => ?_
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun e _ => by ring
  have hB : (∑ i, ∑ t : Fin T, J i / dτ * (f (i, t.castSucc) - f (i, t.succ)) ^ 2)
      = ∑ t : Fin T, ∑ i, J i / dτ * (f (i, t.succ) - f (i, t.castSucc)) ^ 2 := by
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun t _ => Finset.sum_congr rfl fun i _ => by ring
  have hC : (∑ i, ∑ t : Fin (T + 1), dτ * g i * f (i, t) ^ 2)
      = ∑ t : Fin (T + 1), dτ * ∑ i, g i * f (i, t) ^ 2 := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun t _ => ?_
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  rw [hA, hB, hC]
  have hR1 : ∑ t : Fin T, ∑ i, J i / dτ * (f (i, t.succ) - f (i, t.castSucc)) ^ 2 / 2
      = 1 / 2 * ∑ t : Fin T, ∑ i, J i / dτ * (f (i, t.succ) - f (i, t.castSucc)) ^ 2 := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun t _ => ?_
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring
  have hR2 : ∑ t : Fin (T + 1), dτ * (∑ e, c e * (f (src e, t) - f (tgt e, t)) ^ 2
        + ∑ i, g i * f (i, t) ^ 2) / 2
      = 1 / 2 * ∑ t : Fin (T + 1), dτ * ∑ e, c e * (f (src e, t) - f (tgt e, t)) ^ 2
        + 1 / 2 * ∑ t : Fin (T + 1), dτ * ∑ i, g i * f (i, t) ^ 2 := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun t _ => by ring
  rw [hR1, hR2]
  ring

/-- With every node grounded, nonnegative conductances and time coefficients and positive slab
thickness, the register stiffness is positive definite. -/
theorem register_posDef (src tgt : E → ι) (c : E → ℝ) (g J : ι → ℝ) (dτ : ℝ) (T : ℕ)
    (hc : ∀ e, 0 ≤ c e) (hg : ∀ i, 0 < g i) (hJ : ∀ i, 0 ≤ J i) (hd : 0 < dτ) :
    (stiffnessMatrix (regSrc src T) (regTgt tgt T) (regCond c J dτ T)
      (regGround g dτ T)).PosDef := by
  refine stiffnessMatrix_posDef _ _ _ _ ?_ ?_
  · rintro (⟨e, t⟩ | ⟨i, t⟩)
    · exact mul_nonneg hd.le (hc e)
    · exact div_nonneg (hJ i) hd.le
  · intro p
    exact mul_pos hd (hg p.1)

end Action

#print axioms timeSeam_rate_eq_iff
#print axioms uniform_time_rate_iff
#print axioms register_drift
#print axioms registerAction_eq
#print axioms register_posDef

end OPH.SourceLaw
