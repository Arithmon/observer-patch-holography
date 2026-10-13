import SourceLaw.GoldenCarrier
import Mathlib.NumberTheory.Real.GoldenRatio

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The golden refinement family at every level

The source scalar packet uses, at its levels `q = F_k` (Fibonacci numbers; OPH's scalar adapter
accepts `q = 13, 34, 89, 233`), the fractional parts `fract(b φ)` for `0 ≤ b < q` in increasing
order, followed by the Dirichlet marker `1`. This file defines that grid for every `q` and proves:

* `goldenFract_injective`, `goldenPoints_card`: the points `fract(b φ)` are pairwise distinct
  (irrationality of `φ`), so there are exactly `q` of them. The first is the statement of OPH's
  `OPH.GoldenSourceAssignment.orbit_injective`, restated because this development imports only
  Mathlib.
* `goldenGrid_strictMono`, `goldenGrid_zero`, `goldenGrid_last`: the grid is strictly increasing,
  starts at `0` and ends at the marker `1`.
* `golden_packet_posDef`, `golden_whitened_posDef`: at every level `q = n + 2` the packet stiffness
  and its capacity-whitened form are positive definite, so the carrier premises of the register
  theorems hold along the whole family (for the levels `q = 13, 34, 89, 233` of the source family
  and every other `q ≥ 2`).
* `whitenedPacket_kroneckerSum`: on every carrier with strictly increasing coordinates the
  capacity-whitened packet is the Kronecker sum `B ⊗ 1 ⊗ 1 + 1 ⊗ B ⊗ 1 + 1 ⊗ 1 ⊗ B + m² 1` of the
  line operator `B = W^{-1/2} K W^{-1/2}` (with `m² = 1` this is the operator `A_q` of the source
  packet).
* `golden13_eq_goldenGrid`: the explicit `q = 13` carrier of the companion module `GoldenCarrier`
  is this grid at `q = 13`; `golden34_whitened_posDef` is the `q = 34` instance.
* `golden_calibrated_separation`: at every level `q ≥ 3` no `αI + βA_q` equals `A_q + εA_q²` for
  `ε ≠ 0`; between the diagonal neighbours `(0,0,0)` and `(1,1,0)` the entry of `A_q` is `0` and
  that of `A_q²` is `2 B₀₁²` (`golden_distance_two_entry`), with `B₀₁ ≠ 0` (`lineOp_zero_one_ne`).

Negative control: `rational_slope_collision` shows two coinciding points for the rational slope
`3/2`, the case the irrationality argument excludes.
-/

namespace OPH.SourceLaw

open scoped Kronecker
open Matrix

noncomputable section

/-- The points `fract(b φ)` for `0 ≤ b < q`. -/
def goldenPoints (q : ℕ) : Finset ℝ :=
  (Finset.range q).image fun b : ℕ => Int.fract ((b : ℝ) * Real.goldenRatio)

/-- The points `fract(b φ)` are pairwise distinct. This is the statement of OPH's
`OPH.GoldenSourceAssignment.orbit_injective` (`Lean/Geometry/GoldenSourceAssignment.lean`). -/
theorem goldenFract_injective :
    Function.Injective fun b : ℕ => Int.fract ((b : ℝ) * Real.goldenRatio) := by
  intro a b h
  simp only at h
  rw [Int.fract_eq_fract] at h
  obtain ⟨z, hz⟩ := h
  by_contra hab
  have hne : ((a : ℤ) - b : ℤ) ≠ 0 := by omega
  have hirr : Irrational (((a : ℤ) - b : ℤ) * Real.goldenRatio) :=
    Irrational.intCast_mul Real.goldenRatio_irrational hne
  apply hirr.ne_int z
  push_cast
  linarith

theorem goldenPoints_card (q : ℕ) : (goldenPoints q).card = q := by
  rw [goldenPoints, Finset.card_image_of_injective _ goldenFract_injective, Finset.card_range]

/-- The golden grid of level `q`: the increasing enumeration of `goldenPoints q`, then `1`. -/
def goldenGrid (q : ℕ) (i : Fin (q + 1)) : ℝ :=
  if h : (i : ℕ) < q then (goldenPoints q).orderEmbOfFin (goldenPoints_card q) ⟨i, h⟩ else 1

theorem goldenPoints_lt_one {q : ℕ} {x : ℝ} (hx : x ∈ goldenPoints q) : x < 1 := by
  rw [goldenPoints, Finset.mem_image] at hx
  obtain ⟨b, _, rfl⟩ := hx
  exact Int.fract_lt_one _

theorem goldenPoints_nonneg {q : ℕ} {x : ℝ} (hx : x ∈ goldenPoints q) : 0 ≤ x := by
  rw [goldenPoints, Finset.mem_image] at hx
  obtain ⟨b, _, rfl⟩ := hx
  exact Int.fract_nonneg _

theorem goldenGrid_strictMono (q : ℕ) : StrictMono (goldenGrid q) := by
  intro i j hij
  by_cases hjq : (j : ℕ) < q
  · have hiq : (i : ℕ) < q := lt_trans hij hjq
    simp only [goldenGrid, dif_pos hiq, dif_pos hjq]
    exact ((goldenPoints q).orderEmbOfFin (goldenPoints_card q)).strictMono
      (show (⟨i, hiq⟩ : Fin q) < ⟨j, hjq⟩ from hij)
  · have hiq : (i : ℕ) < q := by
      have := j.isLt
      omega
    simp only [goldenGrid, dif_pos hiq, dif_neg hjq]
    exact goldenPoints_lt_one (Finset.orderEmbOfFin_mem _ _ _)

theorem goldenGrid_zero {q : ℕ} (hq : 1 ≤ q) : goldenGrid q 0 = 0 := by
  have h0 : ((0 : Fin (q + 1)) : ℕ) < q := by simp; omega
  simp only [goldenGrid, dif_pos h0]
  have hmem : (0 : ℝ) ∈ goldenPoints q := by
    rw [goldenPoints, Finset.mem_image]
    exact ⟨0, Finset.mem_range.2 (by omega), by simp⟩
  have hz : ((⟨0, by omega⟩ : Fin q)) = ⟨((0 : Fin (q + 1)) : ℕ), h0⟩ := by simp
  rw [← hz, Finset.orderEmbOfFin_zero (goldenPoints_card q) (by omega)]
  apply le_antisymm
  · exact Finset.min'_le _ _ hmem
  · exact Finset.le_min' _ _ _ fun y hy => goldenPoints_nonneg hy

theorem goldenGrid_last (q : ℕ) : goldenGrid q (Fin.last q) = 1 := by
  simp [goldenGrid]

/-- At every level `q = n + 2` the packet stiffness is positive definite. -/
theorem golden_packet_posDef (n : ℕ) {m2 : ℝ} (hm : 0 < m2) :
    (packetV (goldenGrid (n + 2)) m2).PosDef :=
  packetV_posDef (goldenGrid_strictMono (n + 2)) hm

/-- At every level `q = n + 2` the capacity-whitened packet is positive definite. -/
theorem golden_whitened_posDef (n : ℕ) {m2 : ℝ} (hm : 0 < m2) :
    (whitenedPacket (goldenGrid (n + 2)) m2).PosDef :=
  whitenedPacket_posDef (goldenGrid_strictMono (n + 2)) hm

section Kronecker

variable {n : ℕ} {x : Fin (n + 3) → ℝ}

/-- `W^{-1/2}` on the interior nodes of a line. -/
def lineInvSqrt (x : Fin (n + 3) → ℝ) : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
  diagonal fun i => (√(weight x i))⁻¹

/-- The line operator `B = W^{-1/2} K W^{-1/2}`. -/
def lineOp (x : Fin (n + 3) → ℝ) : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
  lineInvSqrt x * dirichlet x * lineInvSqrt x

theorem capInvSqrt_kronecker (hx : StrictMono x) :
    diagonal (fun s : Site n => (√(capacity3 x s))⁻¹)
      = lineInvSqrt x ⊗ₖ (lineInvSqrt x ⊗ₖ lineInvSqrt x) := by
  rw [lineInvSqrt, diagonal_kronecker_diagonal, diagonal_kronecker_diagonal]
  congr 1
  funext s
  obtain ⟨a, b, c⟩ := s
  simp only [capacity3]
  rw [Real.sqrt_mul (weight_pos hx a).le, Real.sqrt_mul (weight_pos hx b).le, mul_inv, mul_inv]

theorem lumped_whiten (hx : StrictMono x) : lineInvSqrt x * lumped x * lineInvSqrt x = 1 := by
  rw [lineInvSqrt, lumped, diagonal_mul_diagonal, diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext i
  have hw := weight_pos hx i
  have hsq : √(weight x i) * √(weight x i) = weight x i := Real.mul_self_sqrt hw.le
  calc (√(weight x i))⁻¹ * weight x i * (√(weight x i))⁻¹
      = weight x i * (√(weight x i) * √(weight x i))⁻¹ := by ring
    _ = 1 := by rw [hsq, mul_inv_cancel₀ hw.ne']

/-- The whitened packet is the Kronecker sum of the line operator plus `m² 1`. -/
theorem whitenedPacket_kroneckerSum (hx : StrictMono x) (m2 : ℝ) :
    whitenedPacket x m2
      = lineOp x ⊗ₖ (1 : Matrix (Fin (n + 1) × Fin (n + 1)) (Fin (n + 1) × Fin (n + 1)) ℝ)
        + (1 : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) ⊗ₖ (lineOp x ⊗ₖ
          (1 : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ))
        + (1 : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) ⊗ₖ
          ((1 : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) ⊗ₖ lineOp x)
        + m2 • (1 : Matrix (Site n) (Site n) ℝ) := by
  have hW := lumped_whiten hx
  have key : ∀ X Y Z : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ,
      (lineInvSqrt x ⊗ₖ (lineInvSqrt x ⊗ₖ lineInvSqrt x)) * (X ⊗ₖ (Y ⊗ₖ Z))
          * (lineInvSqrt x ⊗ₖ (lineInvSqrt x ⊗ₖ lineInvSqrt x))
        = (lineInvSqrt x * X * lineInvSqrt x) ⊗ₖ
          ((lineInvSqrt x * Y * lineInvSqrt x) ⊗ₖ (lineInvSqrt x * Z * lineInvSqrt x)) := by
    intro X Y Z
    rw [← mul_kronecker_mul, ← mul_kronecker_mul, ← mul_kronecker_mul, ← mul_kronecker_mul]
  unfold whitenedPacket packetV packetStiffness packetMass lineOp
  rw [capInvSqrt_kronecker hx]
  simp only [Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.smul_mul, key, hW,
    one_kronecker_one]

end Kronecker

/-- The explicit `q = 13` carrier is the golden grid at `q = 13`. -/
theorem golden13_eq_goldenGrid : golden13 = goldenGrid 13 := by
  obtain ⟨hmono, hlast, hval, _⟩ := golden13_sorted_orbit
  have hbound : ∀ i : Fin 13, goldenMult13 i < 13 := by decide
  have hmem : ∀ i : Fin 13, golden13 i.castSucc ∈ goldenPoints 13 := by
    intro i
    rw [hval i, goldenPoints, Finset.mem_image]
    exact ⟨goldenMult13 i, Finset.mem_range.2 (hbound i), rfl⟩
  have hsm : StrictMono fun i : Fin 13 => golden13 i.castSucc :=
    hmono.comp Fin.strictMono_castSucc
  have huniq := Finset.orderEmbOfFin_unique (goldenPoints_card 13) hmem hsm
  funext i
  by_cases hi : (i : ℕ) < 13
  · have hcast : i = (⟨i, hi⟩ : Fin 13).castSucc := Fin.ext rfl
    rw [hcast]
    simp only [goldenGrid, Fin.val_castSucc, dif_pos hi]
    exact congrFun huniq ⟨i, hi⟩
  · have hlast' : i = Fin.last 13 := Fin.ext (by have := i.isLt; simp; omega)
    rw [hlast', hlast, goldenGrid_last]

/-- The `q = 34` level: the capacity-whitened packet with `m² = 1` is positive definite. -/
theorem golden34_whitened_posDef : (whitenedPacket (goldenGrid 34) 1).PosDef :=
  golden_whitened_posDef 32 one_pos

/-- Negative control: with slope `3/2` the points for `b = 0` and `b = 2` coincide. -/
theorem rational_slope_collision :
    Int.fract (((0 : ℕ) : ℝ) * (3 / 2)) = Int.fract (((2 : ℕ) : ℝ) * (3 / 2)) := by
  norm_num [Int.fract]

section Calibration

/-- The entry `B₀₁` of the line operator is nonzero on every strictly increasing grid with at least
two interior nodes. -/
theorem lineOp_zero_one_ne {n : ℕ} {x : Fin (n + 4) → ℝ} (hx : StrictMono x) :
    lineOp x 0 1 ≠ 0 := by
  unfold lineOp lineInvSqrt
  rw [Matrix.mul_diagonal, Matrix.diagonal_mul]
  have hw0 : 0 < weight x 0 := weight_pos hx 0
  have hw1 : 0 < weight x 1 := weight_pos hx 1
  have hg : 0 < gap x (1 : Fin (n + 2)).castSucc := gap_pos hx _
  have hd : dirichlet x 0 1 = -(1 / gap x (0 : Fin (n + 2)).succ) := by
    simp [dirichlet]
  rw [hd]
  have hg' : 0 < gap x (0 : Fin (n + 2)).succ := gap_pos hx _
  have h0 : 0 < (Real.sqrt (weight x 0))⁻¹ := inv_pos.2 (Real.sqrt_pos.2 hw0)
  have h1 : 0 < (Real.sqrt (weight x 1))⁻¹ := inv_pos.2 (Real.sqrt_pos.2 hw1)
  have : (Real.sqrt (weight x 0))⁻¹ * -(1 / gap x (0 : Fin (n + 2)).succ)
      * (Real.sqrt (weight x 1))⁻¹ < 0 := by
    have : 0 < 1 / gap x (0 : Fin (n + 2)).succ := by positivity
    nlinarith [mul_pos (mul_pos h0 this) h1]
  exact this.ne

/-- At every golden level `q = n + 3` the whitened packet `A_q` has entry `0` between the diagonal
neighbours `(0,0,0)` and `(1,1,0)`, and `A_q²` has entry `2 B₀₁²` there. -/
theorem golden_distance_two_entry (n : ℕ) (m2 : ℝ) :
    whitenedPacket (goldenGrid (n + 3)) m2
        ((0 : Fin (n + 2)), ((0 : Fin (n + 2)), (0 : Fin (n + 2))))
        ((1 : Fin (n + 2)), ((1 : Fin (n + 2)), (0 : Fin (n + 2)))) = 0
      ∧ (whitenedPacket (goldenGrid (n + 3)) m2 * whitenedPacket (goldenGrid (n + 3)) m2)
          ((0 : Fin (n + 2)), ((0 : Fin (n + 2)), (0 : Fin (n + 2))))
          ((1 : Fin (n + 2)), ((1 : Fin (n + 2)), (0 : Fin (n + 2))))
        = 2 * (lineOp (goldenGrid (n + 3)) 0 1) ^ 2 := by
  rw [whitenedPacket_kroneckerSum (goldenGrid_strictMono (n + 3))]
  constructor
  · simp [Matrix.kroneckerMap_apply]
  · simp only [Matrix.add_mul, Matrix.mul_add, ← Matrix.mul_kronecker_mul, Matrix.smul_mul,
      Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one]
    simp [Matrix.kroneckerMap_apply]
    ring

/-- At every golden level `q ≥ 3` no `αI + βA_q` equals `A_q + εA_q²` for `ε ≠ 0`. -/
theorem golden_calibrated_separation (n : ℕ) (m2 ε : ℝ) :
    (∃ α β : ℝ, whitenedPacket (goldenGrid (n + 3)) m2
          + ε • (whitenedPacket (goldenGrid (n + 3)) m2 * whitenedPacket (goldenGrid (n + 3)) m2)
        = α • (1 : Matrix (Site (n + 1)) (Site (n + 1)) ℝ)
          + β • whitenedPacket (goldenGrid (n + 3)) m2)
      ↔ ε = 0 := by
  constructor
  · rintro ⟨α, β, h⟩
    obtain ⟨h0, h2⟩ := golden_distance_two_entry n m2
    have hB := lineOp_zero_one_ne (goldenGrid_strictMono (n + 3))
    have := congrFun (congrFun h ((0 : Fin (n + 2)), ((0 : Fin (n + 2)), (0 : Fin (n + 2)))))
      ((1 : Fin (n + 2)), ((1 : Fin (n + 2)), (0 : Fin (n + 2))))
    rw [Matrix.add_apply, Matrix.smul_apply, Matrix.add_apply, Matrix.smul_apply, Matrix.smul_apply,
      h0, h2] at this
    simp at this
    rcases this with h | h
    · exact h
    · exact absurd h hB
  · rintro rfl
    exact ⟨0, 1, by simp⟩

end Calibration

#print axioms goldenFract_injective
#print axioms goldenPoints_card
#print axioms goldenPoints_lt_one
#print axioms goldenPoints_nonneg
#print axioms goldenGrid_strictMono
#print axioms goldenGrid_zero
#print axioms goldenGrid_last
#print axioms golden_packet_posDef
#print axioms golden_whitened_posDef
#print axioms capInvSqrt_kronecker
#print axioms lumped_whiten
#print axioms whitenedPacket_kroneckerSum
#print axioms golden13_eq_goldenGrid
#print axioms golden34_whitened_posDef
#print axioms rational_slope_collision
#print axioms lineOp_zero_one_ne
#print axioms golden_distance_two_entry
#print axioms golden_calibrated_separation

end

end OPH.SourceLaw
