import SourceLaw.MatrixSeams
import SourceLaw.TwoPointFlux
import Mathlib.LinearAlgebra.Matrix.Kronecker

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The scalar packet, recomputed from one-dimensional coordinates

The source scalar packet builds its carrier from sorted coordinates
`0 = x_0 < x_1 < ... < x_{q-1} < x_q = 1`, the last being the Dirichlet marker: gaps
`h_j = x_{j+1} - x_j`, interior lumped weights `w_i = (h_{i-1} + h_i)/2`, the Dirichlet stiffness
`K` of the interior nodes, tensor masses `M = W ⊗ W ⊗ W`,
`K_D = K ⊗ W ⊗ W + W ⊗ K ⊗ W + W ⊗ W ⊗ K` and `V = K_D + m² M`. Interior nodes are indexed by
`Fin (n + 1)` (so `q = n + 2`), coordinates by `Fin (n + 3)`, sites by `Site n`.

Main results, for strictly increasing coordinates:

* `capacity3_eq_cell`: the capacities are the cell volumes `w_a w_b w_c`.
* `packetV_offdiag`: off the diagonal, `V` is minus the two-point flux conductance `tpfa`: for
  axis neighbours the face measure (the product of the two transverse weights) over the centre
  distance (the gap between them), and zero for every other pair. `packetV_axis_x`, `_y`, `_z`
  state the neighbour values with `fluxCond`.
* `packetV_rowSum`: the row sums are the Dirichlet stubs of the boundary layers plus `m² M`.
* `stiffnessMatrix_packet`: `V` is the stiffness matrix of its seam graph (`edge3`, `ground3`).
* `edge3_nonneg`, `ground3_pos`, `capacity3_pos`, `packetV_posDef`, `edge3_src_ne_tgt`: the
  carrier-side premises of the register theorems (nonnegative conductances, positive grounding and
  capacities, positive definite `V`, distinct seam ends) hold on every such carrier with `m² > 0`;
  the time coefficient and `Δτ` are separate inputs.
-/

namespace OPH.SourceLaw

open scoped BigOperators Kronecker
open Matrix

noncomputable section

variable {n : ℕ}

section OneDim

/-- Gap `h_j = x_{j+1} - x_j`. -/
def gap (x : Fin (n + 3) → ℝ) (j : Fin (n + 2)) : ℝ := x j.succ - x j.castSucc

/-- Lumped weight of interior node `i` (coordinate `x_{i+1}`): half of its two gaps. -/
def weight (x : Fin (n + 3) → ℝ) (i : Fin (n + 1)) : ℝ := (gap x i.castSucc + gap x i.succ) / 2

/-- Dirichlet stiffness of the interior nodes. -/
def dirichlet (x : Fin (n + 3) → ℝ) : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
  Matrix.of fun i j =>
    if i = j then 1 / gap x i.castSucc + 1 / gap x i.succ
    else if (j : ℕ) = i + 1 then -(1 / gap x i.succ)
    else if (i : ℕ) = j + 1 then -(1 / gap x i.castSucc)
    else 0

/-- Reciprocal centre distance between neighbouring interior nodes, zero otherwise. -/
def chainCond (x : Fin (n + 3) → ℝ) (i j : Fin (n + 1)) : ℝ :=
  if (j : ℕ) = i + 1 then 1 / gap x i.succ
  else if (i : ℕ) = j + 1 then 1 / gap x i.castSucc
  else 0

/-- Dirichlet stubs: the two boundary nodes are tied to the zero reference across their outer
gaps. -/
def bnd (x : Fin (n + 3) → ℝ) (i : Fin (n + 1)) : ℝ :=
  (if (i : ℕ) = 0 then 1 / gap x i.castSucc else 0)
    + (if (i : ℕ) = n then 1 / gap x i.succ else 0)

/-- Diagonal lumped mass. -/
def lumped (x : Fin (n + 3) → ℝ) : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ := diagonal (weight x)

theorem succ_eq_castSucc_of_val {i j : Fin (n + 1)} (h : (i : ℕ) = j + 1) :
    j.succ = i.castSucc :=
  Fin.ext (by simp [h])

theorem dirichlet_offdiag (x : Fin (n + 3) → ℝ) {i j : Fin (n + 1)} (hij : i ≠ j) :
    dirichlet x i j = -chainCond x i j := by
  simp only [dirichlet, chainCond, of_apply, if_neg hij]
  split_ifs <;> ring

theorem dirichlet_transpose (x : Fin (n + 3) → ℝ) : (dirichlet x)ᵀ = dirichlet x := by
  ext i j
  rw [transpose_apply]
  by_cases hij : i = j
  · rw [hij]
  · rw [dirichlet_offdiag x hij, dirichlet_offdiag x (Ne.symm hij), chainCond, chainCond]
    by_cases h1 : (i : ℕ) = j + 1
    · have h2 : ¬ (j : ℕ) = i + 1 := by omega
      rw [if_pos h1, if_neg h2, if_pos h1, succ_eq_castSucc_of_val h1]
    · by_cases h2 : (j : ℕ) = i + 1
      · rw [if_neg h1, if_pos h2, if_pos h2, succ_eq_castSucc_of_val h2]
      · rw [if_neg h1, if_neg h2, if_neg h2, if_neg h1]

theorem dirichlet_apply_split (x : Fin (n + 3) → ℝ) (i j : Fin (n + 1)) :
    dirichlet x i j = (if i = j then 1 / gap x i.castSucc + 1 / gap x i.succ else 0)
      + (if (j : ℕ) = i + 1 then -(1 / gap x i.succ) else 0)
      + (if (i : ℕ) = j + 1 then -(1 / gap x i.castSucc) else 0) := by
  simp only [dirichlet, of_apply]
  by_cases hij : i = j
  · subst hij
    simp
  · rw [if_neg hij, if_neg hij, zero_add]
    by_cases h1 : (j : ℕ) = i + 1
    · have h2 : ¬ (i : ℕ) = j + 1 := by omega
      rw [if_pos h1, if_pos h1, if_neg h2, add_zero]
    · rw [if_neg h1, if_neg h1, zero_add]

theorem sum_fin_ite_val (m : ℕ) (D : ℝ) :
    ∑ j : Fin (n + 1), (if (j : ℕ) = m then D else 0) = if m < n + 1 then D else 0 := by
  calc ∑ j : Fin (n + 1), (if (j : ℕ) = m then D else 0)
      = ∑ j ∈ Finset.range (n + 1), (if j = m then D else 0) :=
        Fin.sum_univ_eq_sum_range (fun j => if j = m then D else 0) (n + 1)
    _ = if m < n + 1 then D else 0 := by
      rw [Finset.sum_ite_eq']
      simp only [Finset.mem_range]

/-- The row sums of the Dirichlet stiffness are the boundary stubs. -/
theorem dirichlet_rowSum (x : Fin (n + 3) → ℝ) (i : Fin (n + 1)) :
    ∑ j, dirichlet x i j = bnd x i := by
  have hi := i.isLt
  have h3 : ∑ j : Fin (n + 1), (if (i : ℕ) = j + 1 then -(1 / gap x i.castSucc) else 0)
      = if 0 < (i : ℕ) then -(1 / gap x i.castSucc) else 0 := by
    by_cases h0 : (i : ℕ) = 0
    · have hne : ∀ j : Fin (n + 1), ¬ (i : ℕ) = j + 1 := fun j => by omega
      simp [h0]
    · have hiff : ∀ j : Fin (n + 1), ((i : ℕ) = j + 1) ↔ ((j : ℕ) = (i : ℕ) - 1) :=
        fun j => by omega
      simp only [hiff]
      rw [sum_fin_ite_val, if_pos (by omega), if_pos (by omega)]
  simp only [dirichlet_apply_split, Finset.sum_add_distrib]
  rw [Finset.sum_ite_eq, if_pos (Finset.mem_univ _), sum_fin_ite_val, h3, bnd]
  split_ifs <;> first | ring1 | omega

variable {x : Fin (n + 3) → ℝ}

theorem gap_pos (hx : StrictMono x) (j : Fin (n + 2)) : 0 < gap x j :=
  sub_pos.2 (hx Fin.castSucc_lt_succ)

theorem weight_pos (hx : StrictMono x) (i : Fin (n + 1)) : 0 < weight x i :=
  div_pos (add_pos (gap_pos hx _) (gap_pos hx _)) two_pos

theorem bnd_nonneg (hx : StrictMono x) (i : Fin (n + 1)) : 0 ≤ bnd x i := by
  have h1 := gap_pos hx i.castSucc
  have h2 := gap_pos hx i.succ
  unfold bnd
  split_ifs <;> positivity

theorem chainCond_nonneg (hx : StrictMono x) (i j : Fin (n + 1)) : 0 ≤ chainCond x i j := by
  have h1 := gap_pos hx i.castSucc
  have h2 := gap_pos hx i.succ
  unfold chainCond
  split_ifs <;> positivity

end OneDim

/-- Sites of the three-dimensional carrier. -/
abbrev Site (n : ℕ) := Fin (n + 1) × Fin (n + 1) × Fin (n + 1)

section Packet

/-- Capacity of a site: the product of its three lumped weights. -/
def capacity3 (x : Fin (n + 3) → ℝ) (s : Site n) : ℝ :=
  weight x s.1 * (weight x s.2.1 * weight x s.2.2)

/-- Tensor mass `W ⊗ W ⊗ W`. -/
def packetMass (x : Fin (n + 3) → ℝ) : Matrix (Site n) (Site n) ℝ :=
  lumped x ⊗ₖ (lumped x ⊗ₖ lumped x)

/-- Tensor stiffness `K ⊗ W ⊗ W + W ⊗ K ⊗ W + W ⊗ W ⊗ K`. -/
def packetStiffness (x : Fin (n + 3) → ℝ) : Matrix (Site n) (Site n) ℝ :=
  dirichlet x ⊗ₖ (lumped x ⊗ₖ lumped x) + lumped x ⊗ₖ (dirichlet x ⊗ₖ lumped x)
    + lumped x ⊗ₖ (lumped x ⊗ₖ dirichlet x)

/-- The packet operator `V = K_D + m² M`. -/
def packetV (x : Fin (n + 3) → ℝ) (m2 : ℝ) : Matrix (Site n) (Site n) ℝ :=
  packetStiffness x + m2 • packetMass x

/-- Two-point flux conductance between two sites: for axis neighbours, the product of the two
transverse weights (face measure) over the gap between them (centre distance); zero otherwise. -/
def tpfa (x : Fin (n + 3) → ℝ) (k l : Site n) : ℝ :=
  (if k.2 = l.2 then chainCond x k.1 l.1 * (weight x k.2.1 * weight x k.2.2) else 0)
    + (if k.1 = l.1 ∧ k.2.2 = l.2.2 then
        weight x k.1 * (chainCond x k.2.1 l.2.1 * weight x k.2.2) else 0)
    + (if k.1 = l.1 ∧ k.2.1 = l.2.1 then
        weight x k.1 * (weight x k.2.1 * chainCond x k.2.2 l.2.2) else 0)

/-- Spatial seam conductances of the packet graph. -/
def edge3 (x : Fin (n + 3) → ℝ) (m2 : ℝ) : MatSeam (Site n) → ℝ := matCond (packetV x m2)

/-- Grounding of the packet graph: Dirichlet stubs plus `m² M`. -/
def ground3 (x : Fin (n + 3) → ℝ) (m2 : ℝ) : Site n → ℝ := rowSum (packetV x m2)

/-- The capacities are the cell volumes. -/
theorem capacity3_eq_cell (x : Fin (n + 3) → ℝ) : packetMass x = diagonal (capacity3 x) := by
  unfold packetMass lumped
  rw [diagonal_kronecker_diagonal, diagonal_kronecker_diagonal]
  rfl

theorem packetV_apply (x : Fin (n + 3) → ℝ) (m2 : ℝ) (k l : Site n) :
    packetV x m2 k l
      = dirichlet x k.1 l.1 * (lumped x k.2.1 l.2.1 * lumped x k.2.2 l.2.2)
        + lumped x k.1 l.1 * (dirichlet x k.2.1 l.2.1 * lumped x k.2.2 l.2.2)
        + lumped x k.1 l.1 * (lumped x k.2.1 l.2.1 * dirichlet x k.2.2 l.2.2)
        + m2 * (lumped x k.1 l.1 * (lumped x k.2.1 l.2.1 * lumped x k.2.2 l.2.2)) := by
  obtain ⟨a, b, c⟩ := k
  obtain ⟨a', b', c'⟩ := l
  simp only [packetV, packetStiffness, packetMass, add_apply, smul_apply, smul_eq_mul,
    kronecker_apply]

/-- Off the diagonal, the packet operator is minus the two-point flux conductance. -/
theorem packetV_offdiag (x : Fin (n + 3) → ℝ) (m2 : ℝ) (k l : Site n) (hkl : k ≠ l) :
    packetV x m2 k l = -tpfa x k l := by
  obtain ⟨a, b, c⟩ := k
  obtain ⟨a', b', c'⟩ := l
  rw [packetV_apply]
  simp only [lumped, diagonal_apply, tpfa]
  by_cases ha : a = a' <;> by_cases hb : b = b' <;> by_cases hc : c = c' <;>
    simp_all [dirichlet_offdiag]

theorem packetV_axis_x (x : Fin (n + 3) → ℝ) (m2 : ℝ) (a : Fin n) (b c : Fin (n + 1)) :
    packetV x m2 (a.castSucc, b, c) (a.succ, b, c)
      = -fluxCond (weight x b * weight x c) (gap x a.castSucc.succ) := by
  have hne : a.castSucc ≠ a.succ := Fin.castSucc_lt_succ.ne
  rw [packetV_offdiag x m2 _ _ (by simp [hne])]
  simp [tpfa, chainCond, fluxCond, hne]
  ring

theorem packetV_axis_y (x : Fin (n + 3) → ℝ) (m2 : ℝ) (a : Fin (n + 1)) (b : Fin n)
    (c : Fin (n + 1)) :
    packetV x m2 (a, b.castSucc, c) (a, b.succ, c)
      = -fluxCond (weight x a * weight x c) (gap x b.castSucc.succ) := by
  have hne : b.castSucc ≠ b.succ := Fin.castSucc_lt_succ.ne
  rw [packetV_offdiag x m2 _ _ (by simp [hne])]
  simp [tpfa, chainCond, fluxCond, hne]
  ring

theorem packetV_axis_z (x : Fin (n + 3) → ℝ) (m2 : ℝ) (a b : Fin (n + 1)) (c : Fin n) :
    packetV x m2 (a, b, c.castSucc) (a, b, c.succ)
      = -fluxCond (weight x a * weight x b) (gap x c.castSucc.succ) := by
  have hne : c.castSucc ≠ c.succ := Fin.castSucc_lt_succ.ne
  rw [packetV_offdiag x m2 _ _ (by simp [hne])]
  simp [tpfa, chainCond, fluxCond, hne]
  ring

theorem rowSum_kronecker {ι κ : Type*} [Fintype ι] [Fintype κ] (A : Matrix ι ι ℝ)
    (B : Matrix κ κ ℝ) (i : ι) (j : κ) :
    rowSum (A ⊗ₖ B) (i, j) = rowSum A i * rowSum B j := by
  simp only [rowSum, Fintype.sum_prod_type, kronecker_apply, Finset.sum_mul_sum]

theorem rowSum_add {ι : Type*} [Fintype ι] (A B : Matrix ι ι ℝ) (i : ι) :
    rowSum (A + B) i = rowSum A i + rowSum B i := by
  simp only [rowSum, add_apply, Finset.sum_add_distrib]

theorem rowSum_smul {ι : Type*} [Fintype ι] (r : ℝ) (A : Matrix ι ι ℝ) (i : ι) :
    rowSum (r • A) i = r * rowSum A i := by
  simp only [rowSum, smul_apply, smul_eq_mul, Finset.mul_sum]

theorem rowSum_diagonal {ι : Type*} [Fintype ι] [DecidableEq ι] (d : ι → ℝ) (i : ι) :
    rowSum (diagonal d) i = d i := by
  simp [rowSum, diagonal_apply]

/-- Row sums of the packet operator: Dirichlet stubs of the boundary layers plus `m² M`. -/
theorem packetV_rowSum (x : Fin (n + 3) → ℝ) (m2 : ℝ) (k : Site n) :
    rowSum (packetV x m2) k
      = bnd x k.1 * (weight x k.2.1 * weight x k.2.2)
        + weight x k.1 * (bnd x k.2.1 * weight x k.2.2)
        + weight x k.1 * (weight x k.2.1 * bnd x k.2.2) + m2 * capacity3 x k := by
  obtain ⟨a, b, c⟩ := k
  have hK : ∀ i, rowSum (dirichlet x) i = bnd x i := dirichlet_rowSum x
  have hW : ∀ i, rowSum (lumped x) i = weight x i := fun i => rowSum_diagonal _ i
  simp only [packetV, packetStiffness, packetMass, rowSum_add, rowSum_smul, rowSum_kronecker, hK,
    hW, capacity3]

theorem lumped_transpose (x : Fin (n + 3) → ℝ) : (lumped x)ᵀ = lumped x := diagonal_transpose _

theorem packetV_transpose (x : Fin (n + 3) → ℝ) (m2 : ℝ) : (packetV x m2)ᵀ = packetV x m2 := by
  simp only [packetV, packetStiffness, packetMass, transpose_add, transpose_smul,
    ← kroneckerMap_transpose, dirichlet_transpose, lumped_transpose]

/-- The packet operator is the stiffness matrix of its seam graph. -/
theorem stiffnessMatrix_packet (x : Fin (n + 3) → ℝ) (m2 : ℝ) :
    stiffnessMatrix matSrc matTgt (edge3 x m2) (ground3 x m2) = packetV x m2 :=
  stiffnessMatrix_matSeams _ (packetV_transpose x m2)

theorem edge3_src_ne_tgt (e : MatSeam (Site n)) : matSrc e ≠ matTgt e := matSrc_ne_tgt e

variable {x : Fin (n + 3) → ℝ}

theorem capacity3_pos (hx : StrictMono x) (s : Site n) : 0 < capacity3 x s :=
  mul_pos (weight_pos hx _) (mul_pos (weight_pos hx _) (weight_pos hx _))

theorem tpfa_nonneg (hx : StrictMono x) (k l : Site n) : 0 ≤ tpfa x k l := by
  have hw := weight_pos hx
  have hc := chainCond_nonneg hx
  unfold tpfa
  refine add_nonneg (add_nonneg ?_ ?_) ?_ <;> split_ifs <;>
    first
    | exact le_refl 0
    | exact mul_nonneg (hc _ _) (mul_nonneg (hw _).le (hw _).le)
    | exact mul_nonneg (hw _).le (mul_nonneg (hc _ _) (hw _).le)
    | exact mul_nonneg (hw _).le (mul_nonneg (hw _).le (hc _ _))

theorem packetV_offdiag_nonpos (hx : StrictMono x) (m2 : ℝ) (k l : Site n) (hkl : k ≠ l) :
    packetV x m2 k l ≤ 0 := by
  rw [packetV_offdiag x m2 k l hkl]
  exact neg_nonpos.2 (tpfa_nonneg hx k l)

/-- The coupling between the axis neighbours `(a, b, c)` and `(a + 1, b, c)` is nonzero. -/
theorem packetV_axis_ne_zero (hx : StrictMono x) (m2 : ℝ) (a : Fin n) (b c : Fin (n + 1)) :
    packetV x m2 (a.castSucc, b, c) (a.succ, b, c) ≠ 0 := by
  rw [packetV_axis_x]
  have h1 := gap_pos hx a.castSucc.succ
  have h2 := weight_pos hx b
  have h3 := weight_pos hx c
  unfold fluxCond
  have : 0 < weight x b * weight x c / gap x a.castSucc.succ := by positivity
  linarith

theorem rowSum_packetV_pos (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) (k : Site n) :
    0 < rowSum (packetV x m2) k := by
  rw [packetV_rowSum]
  have hb1 := bnd_nonneg hx k.1
  have hb2 := bnd_nonneg hx k.2.1
  have hb3 := bnd_nonneg hx k.2.2
  have hw1 := weight_pos hx k.1
  have hw2 := weight_pos hx k.2.1
  have hw3 := weight_pos hx k.2.2
  have hc := capacity3_pos hx k
  positivity

theorem edge3_nonneg (hx : StrictMono x) (m2 : ℝ) (e : MatSeam (Site n)) : 0 ≤ edge3 x m2 e :=
  matCond_nonneg _ (packetV_offdiag_nonpos hx m2) e

theorem ground3_pos (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) (k : Site n) :
    0 < ground3 x m2 k :=
  rowSum_packetV_pos hx hm k

theorem packetV_posDef (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) : (packetV x m2).PosDef :=
  posDef_of_matSeams _ (packetV_transpose x m2) (packetV_offdiag_nonpos hx m2)
    (rowSum_packetV_pos hx hm)

end Packet

end

#print axioms capacity3_eq_cell
#print axioms packetV_offdiag
#print axioms packetV_axis_x
#print axioms packetV_axis_y
#print axioms packetV_axis_z
#print axioms packetV_rowSum
#print axioms dirichlet_rowSum
#print axioms stiffnessMatrix_packet
#print axioms edge3_nonneg
#print axioms ground3_pos
#print axioms packetV_posDef

end OPH.SourceLaw
