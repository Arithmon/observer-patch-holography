import SourceLaw.ScalarPacket
import SourceLaw.CounterfamilyLocality
import SourceLaw.LangevinStationarity
import SourceLaw.ReflectionPositivity
import SourceLaw.ContinuumTimeLimit
import SourceLaw.TransferMatrixLift
import SourceLaw.RepairAdmissibility
import SourceLaw.CounterfamilyReader
import SourceLaw.StationaryWindow
import Mathlib.NumberTheory.Real.GoldenRatio
import Mathlib.Analysis.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
import Mathlib.Algebra.Group.Pi.Units

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The golden scalar carriers q = 3 and q = 13

The source scalar packet places its nodes on the sorted golden orbit
`{fract(b φ) : 0 ≤ b < q}` followed by the Dirichlet marker 1 (the nodes are not replaced by
uniform-grid labels). The only external input for the nodes is the sort order, the
multipliers `b_0, ..., b_{q-1}` (`goldenMult3`, `goldenMult13`), checked against the source code
`grid()` (integer pairs equal, values inside the 10⁻⁵⁰ interval enclosures). Lean proves that
these multipliers give the sorted orbit (`golden3_sorted_orbit`, `golden13_sorted_orbit`); every
other quantity is recomputed from the packet definitions in `ScalarPacket`. The instances fix
the mass parameter at `m² = 1`, the source packet's dimensionless mass, and `golden_ladder`
takes the Gibbs factor `β = 1`; the general packet theorems hold for every `m² > 0`. q = 13
(Fibonacci index 7) is the smallest declared level; q = 3 (index 4) is an extension below the
declared levels.

On these carriers the premises of the general theorems hold and their instances fire:

* bindings to the source values: the sorted orbit, the lumped weights and an axis coupling in
  closed form;
* `packet_register_drift`: the mean repair generator of the register is `(P f)_k / M_k`; the
  declared law takes this as its drift matrix by definition;
* rate lemma instances: capacity-matched inertia `J = M` gives the same time rate `2/Δτ` at every
  site, unit inertia `J = 1` gives different rates at two sites (four distinct capacities);
* counterfamily: the inertia `M (M + εV)⁻¹ M` of the `I + εA` family is not node-local at
  `ε = 1/4` and is node-local at `ε = 0`;
* `packet_register_posDef`, `packet_stationary_unique`: the register stiffness is positive
  definite and every solution of the Lyapunov equation of the declared law is `ℏ P⁻¹`;
* `packet_modes`: every whitened mode is positive, so the per-mode transfer structure, bond
  reflection positivity and the `Δτ → 0` limit apply to each mode;
* `golden_repair`: the register with capacity-matched inertia is a repair stiffness at every
  step, the counterfamily member `ε = 1/4` is not admissible for all small steps, and
  admissibility for all small steps with the uniform time rate fixes the inertia to `M`;
* `golden_counterfamily_reader`: at `ε = 1/4` the two-time reader of the counterfamily chain
  returns the counterfamily inertia, which is not node-local; at `ε = 0` the reader is node-local;
* `golden_window`: the corrected register precision inverts the stationary window covariance,
  no node-local term added equally at both end slabs does, and in original coordinates the
  declared register plus the end-slab correction inverts the stationary window covariance.
-/

namespace OPH.SourceLaw

open scoped BigOperators Kronecker
open Matrix Filter Topology

noncomputable section

section GoldenLine

/-- Coordinate `i` of a golden line: `fract(b_i φ)` for `i < q`, then the Dirichlet marker 1. -/
def goldenCoord (q : ℕ) (mult : ℕ → ℕ) (i : ℕ) : ℝ :=
  if i < q then Int.fract ((mult i : ℝ) * Real.goldenRatio) else 1

theorem goldenRatio_bounds :
    (1.6180339 : ℝ) < Real.goldenRatio ∧ Real.goldenRatio < 1.6180340 := by
  have h1 : (2.2360679 : ℝ) < √5 := by
    rw [Real.lt_sqrt (by norm_num)]
    norm_num
  have h2 : √5 < (2.2360680 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num)]
    norm_num
  have hφ : Real.goldenRatio = (1 + √5) / 2 := rfl
  constructor <;> rw [hφ] <;> linarith

theorem fract_mul_goldenRatio {b : ℕ} {z : ℤ} (h1 : (z : ℝ) ≤ b * Real.goldenRatio)
    (h2 : (b : ℝ) * Real.goldenRatio < z + 1) :
    Int.fract ((b : ℝ) * Real.goldenRatio) = b * Real.goldenRatio - z := by
  rw [← Int.self_sub_floor, Int.floor_eq_iff.2 ⟨h1, h2⟩]

/-- Sorted golden-orbit multipliers for q = 3 (Fibonacci index 4). -/
def goldenMult3 : ℕ → ℕ
  | 0 => 0
  | 1 => 2
  | 2 => 1
  | _ => 0

/-- `⌊b_i φ⌋` for q = 3; proof witnesses, checked by `golden3_closed`. -/
def goldenFloor3 : ℕ → ℤ
  | 0 => 0
  | 1 => 3
  | 2 => 1
  | _ => 0

/-- Sorted golden-orbit multipliers for q = 13 (Fibonacci index 7). -/
def goldenMult13 : ℕ → ℕ
  | 0 => 0
  | 1 => 5
  | 2 => 10
  | 3 => 2
  | 4 => 7
  | 5 => 12
  | 6 => 4
  | 7 => 9
  | 8 => 1
  | 9 => 6
  | 10 => 11
  | 11 => 3
  | 12 => 8
  | _ => 0

/-- `⌊b_i φ⌋` for q = 13; proof witnesses, checked by `golden13_closed`. -/
def goldenFloor13 : ℕ → ℤ
  | 0 => 0
  | 1 => 8
  | 2 => 16
  | 3 => 3
  | 4 => 11
  | 5 => 19
  | 6 => 6
  | 7 => 14
  | 8 => 1
  | 9 => 9
  | 10 => 17
  | 11 => 4
  | 12 => 12
  | _ => 0

theorem golden3_closed (i : ℕ) (hi : i < 3) :
    goldenCoord 3 goldenMult3 i = goldenMult3 i * Real.goldenRatio - goldenFloor3 i := by
  obtain ⟨h1, h2⟩ := goldenRatio_bounds
  rw [goldenCoord, if_pos hi]
  interval_cases i <;>
    exact fract_mul_goldenRatio (by norm_num [goldenMult3, goldenFloor3] <;> linarith)
      (by norm_num [goldenMult3, goldenFloor3] <;> linarith)

theorem golden13_closed (i : ℕ) (hi : i < 13) :
    goldenCoord 13 goldenMult13 i = goldenMult13 i * Real.goldenRatio - goldenFloor13 i := by
  obtain ⟨h1, h2⟩ := goldenRatio_bounds
  rw [goldenCoord, if_pos hi]
  interval_cases i <;>
    exact fract_mul_goldenRatio (by norm_num [goldenMult13, goldenFloor13] <;> linarith)
      (by norm_num [goldenMult13, goldenFloor13] <;> linarith)

theorem golden3_lt_succ (i : ℕ) (hi : i < 3) :
    goldenCoord 3 goldenMult3 i < goldenCoord 3 goldenMult3 (i + 1) := by
  obtain ⟨h1, h2⟩ := goldenRatio_bounds
  rw [golden3_closed i hi]
  rcases Nat.lt_or_ge (i + 1) 3 with h | h
  · rw [golden3_closed (i + 1) h]
    interval_cases i <;> norm_num [goldenMult3, goldenFloor3] <;> linarith
  · have h2' : i = 2 := by omega
    subst h2'
    rw [goldenCoord, if_neg (by norm_num)]
    norm_num [goldenMult3, goldenFloor3]
    linarith

theorem golden13_lt_succ (i : ℕ) (hi : i < 13) :
    goldenCoord 13 goldenMult13 i < goldenCoord 13 goldenMult13 (i + 1) := by
  obtain ⟨h1, h2⟩ := goldenRatio_bounds
  rw [golden13_closed i hi]
  rcases Nat.lt_or_ge (i + 1) 13 with h | h
  · rw [golden13_closed (i + 1) h]
    interval_cases i <;> norm_num [goldenMult13, goldenFloor13] <;> linarith
  · have h12 : i = 12 := by omega
    subst h12
    rw [goldenCoord, if_neg (by norm_num)]
    norm_num [goldenMult13, goldenFloor13]
    linarith

/-- The q = 3 golden line. -/
def golden3 : Fin (1 + 3) → ℝ := fun i => goldenCoord 3 goldenMult3 i

/-- The q = 13 golden line. -/
def golden13 : Fin (11 + 3) → ℝ := fun i => goldenCoord 13 goldenMult13 i

theorem golden3_strictMono : StrictMono golden3 := by
  rw [Fin.strictMono_iff_lt_succ]
  intro i
  simpa [golden3] using golden3_lt_succ i i.isLt

theorem golden13_strictMono : StrictMono golden13 := by
  rw [Fin.strictMono_iff_lt_succ]
  intro i
  simpa [golden13] using golden13_lt_succ i i.isLt

/-- The q = 3 coordinates are the golden orbit in increasing order, then the marker 1. -/
theorem golden3_sorted_orbit :
    StrictMono golden3 ∧ golden3 (Fin.last 3) = 1
      ∧ (∀ i : Fin 3, golden3 i.castSucc = Int.fract ((goldenMult3 i : ℝ) * Real.goldenRatio))
      ∧ ((List.range 3).map goldenMult3).Perm (List.range 3) := by
  refine ⟨golden3_strictMono, ?_, fun i => ?_, by decide⟩
  · simp [golden3, goldenCoord]
  · simp [golden3, goldenCoord, i.isLt]

/-- The q = 13 coordinates are the golden orbit in increasing order, then the marker 1. -/
theorem golden13_sorted_orbit :
    StrictMono golden13 ∧ golden13 (Fin.last 13) = 1
      ∧ (∀ i : Fin 13, golden13 i.castSucc = Int.fract ((goldenMult13 i : ℝ) * Real.goldenRatio))
      ∧ ((List.range 13).map goldenMult13).Perm (List.range 13) := by
  refine ⟨golden13_strictMono, ?_, fun i => ?_, by decide⟩
  · simp [golden13, goldenCoord]
  · simp [golden13, goldenCoord, i.isLt]

end GoldenLine

section Values

theorem weight_eq {n : ℕ} (x : Fin (n + 3) → ℝ) (i : Fin (n + 1)) :
    weight x i = (x i.succ.succ - x i.castSucc.castSucc) / 2 := by
  unfold weight gap
  rw [Fin.succ_castSucc]
  ring

theorem golden3_weight0 : weight golden3 0 = -1 / 4 + 1 / 4 * √5 := by
  rw [weight_eq]
  simp only [golden3, Fin.val_succ, Fin.val_castSucc, Fin.val_zero]
  rw [golden3_closed _ (by norm_num), golden3_closed _ (by norm_num)]
  norm_num [goldenMult3, goldenFloor3, Real.goldenRatio]
  ring

theorem golden3_weight1 : weight golden3 1 = 3 / 2 - 1 / 2 * √5 := by
  rw [weight_eq]
  simp only [golden3, Fin.val_succ, Fin.val_castSucc, Fin.val_one]
  rw [golden3_closed 1 (by norm_num)]
  norm_num [goldenCoord, goldenMult3, goldenFloor3, Real.goldenRatio]
  ring

theorem golden13_weight0 : weight golden13 0 = -11 / 2 + 5 / 2 * √5 := by
  rw [weight_eq]
  simp only [golden13, Fin.val_succ, Fin.val_castSucc, Fin.val_zero]
  rw [golden13_closed _ (by norm_num), golden13_closed _ (by norm_num)]
  norm_num [goldenMult13, goldenFloor13, Real.goldenRatio]
  ring

theorem golden13_weight1 : weight golden13 1 = 7 / 4 - 3 / 4 * √5 := by
  rw [weight_eq]
  simp only [golden13, Fin.val_succ, Fin.val_castSucc, Fin.val_one]
  rw [golden13_closed _ (by norm_num), golden13_closed _ (by norm_num)]
  norm_num [goldenMult13, goldenFloor13, Real.goldenRatio]
  ring

theorem golden3_weight_lt : weight golden3 0 < weight golden3 1 := by
  rw [golden3_weight0, golden3_weight1]
  have h : (2.2360679 : ℝ) < √5 := by
    rw [Real.lt_sqrt (by norm_num)]
    norm_num
  have h' : √5 < (2.2360680 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num)]
    norm_num
  linarith

theorem golden13_weight_gt : weight golden13 1 < weight golden13 0 := by
  rw [golden13_weight0, golden13_weight1]
  have h : (2.2360679 : ℝ) < √5 := by
    rw [Real.lt_sqrt (by norm_num)]
    norm_num
  have h' : √5 < (2.2360680 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num)]
    norm_num
  linarith

/-- The q = 3 axis coupling between the corner site and its neighbour is exactly `-1/4`. -/
theorem golden3_axis_coupling : packetV golden3 1 (0, 0, 0) (1, 0, 0) = -1 / 4 := by
  have h := packetV_axis_x golden3 1 (0 : Fin 1) 0 0
  simp only [Fin.castSucc_zero, Fin.succ_zero_eq_one] at h
  rw [h, golden3_weight0, fluxCond]
  have hg : gap golden3 1 = 2 - (1 + √5) / 2 := by
    simp only [gap, golden3, Fin.val_succ, Fin.val_castSucc, Fin.val_one]
    rw [golden3_closed _ (by norm_num), golden3_closed _ (by norm_num)]
    norm_num [goldenMult3, goldenFloor3, Real.goldenRatio]
    ring
  rw [hg]
  have h5 : √5 * √5 = 5 := Real.mul_self_sqrt (by norm_num)
  have hpos : (0 : ℝ) < 2 - (1 + √5) / 2 := by
    have : √5 < 3 := by
      rw [Real.sqrt_lt' (by norm_num)]
      norm_num
    linarith
  have hne : (2 : ℝ) - (1 + √5) / 2 ≠ 0 := hpos.ne'
  have hsq : (-1 / 4 + 1 / 4 * √5) * (-1 / 4 + 1 / 4 * √5) = (2 - (1 + √5) / 2) / 4 := by
    linear_combination (1 / 16) * h5
  rw [hsq, div_right_comm, div_self hne]
  norm_num

end Values

section PacketInstances

variable {n : ℕ} {x : Fin (n + 3) → ℝ}

/-- On a packet carrier the mean repair generator of the register is `(P f)_k / M_k`; the declared
law takes this as its drift matrix by definition, which is not part of this statement. -/
theorem packet_register_drift (hx : StrictMono x) (m2 : ℝ) (J : Site n → ℝ) (dτ : ℝ) (T : ℕ)
    (f : Site n × Fin (T + 1) → ℝ) (k : Site n × Fin (T + 1)) :
    repairGenerator (regCap (capacity3 x) T) (regSrc matSrc T) (regTgt matTgt T)
        (regCond (edge3 x m2) J dτ T) (regGround (ground3 x m2) dτ T) f k
      = stiffness (regSrc matSrc T) (regTgt matTgt T) (regCond (edge3 x m2) J dτ T)
          (regGround (ground3 x m2) dτ T) f k / regCap (capacity3 x) T k :=
  register_drift matSrc matTgt (edge3 x m2) (ground3 x m2) (capacity3 x) J dτ T f
    (fun i => (capacity3_pos hx i).ne')
    (fun i j => (add_pos (capacity3_pos hx i) (capacity3_pos hx j)).ne') matSrc_ne_tgt k

/-- Capacity-matched inertia gives the time rate `2/Δτ` at every site. -/
theorem packet_matched_rate (hx : StrictMono x) (m2 : ℝ) {dτ : ℝ} (hd : dτ ≠ 0) (T : ℕ)
    (s : Site n) (t : Fin T) :
    seamRate (regCap (capacity3 x) T) (regCond (edge3 x m2) (capacity3 x) dτ T (.inr (s, t)))
        (regSrc matSrc T (.inr (s, t))) (regTgt matTgt T (.inr (s, t))) = 2 / dτ :=
  (timeSeam_rate_eq_iff matSrc matTgt (edge3 x m2) (capacity3 x) (capacity3 x) dτ T s t
    (capacity3_pos hx s).ne' hd).2 rfl

/-- With unit inertia the time rate at a site is `2/(Δτ M)`. -/
theorem packet_unit_rate (hx : StrictMono x) (m2 : ℝ) {dτ : ℝ} (hd : dτ ≠ 0) (T : ℕ)
    (s : Site n) (t : Fin T) :
    seamRate (regCap (capacity3 x) T) (regCond (edge3 x m2) (fun _ => 1) dτ T (.inr (s, t)))
        (regSrc matSrc T (.inr (s, t))) (regTgt matTgt T (.inr (s, t)))
      = 2 / (dτ * capacity3 x s) := by
  rw [timeSeam_rate]
  have := (capacity3_pos hx s).ne'
  field_simp
  ring

/-- Unit inertia gives different time rates at two sites of different capacity. -/
theorem packet_unit_rate_ne (hx : StrictMono x) (m2 : ℝ) {dτ : ℝ} (hd : dτ ≠ 0) (T : ℕ)
    (t : Fin T) {s1 s2 : Site n} (h : capacity3 x s1 ≠ capacity3 x s2) :
    seamRate (regCap (capacity3 x) T) (regCond (edge3 x m2) (fun _ => 1) dτ T (.inr (s1, t)))
        (regSrc matSrc T (.inr (s1, t))) (regTgt matTgt T (.inr (s1, t)))
      ≠ seamRate (regCap (capacity3 x) T) (regCond (edge3 x m2) (fun _ => 1) dτ T (.inr (s2, t)))
        (regSrc matSrc T (.inr (s2, t))) (regTgt matTgt T (.inr (s2, t))) := by
  rw [packet_unit_rate hx m2 hd T s1 t, packet_unit_rate hx m2 hd T s2 t]
  intro heq
  apply h
  have h1 := (capacity3_pos hx s1).ne'
  have h2 := (capacity3_pos hx s2).ne'
  rw [div_eq_div_iff (mul_ne_zero hd h1) (mul_ne_zero hd h2)] at heq
  have h3 : dτ * capacity3 x s1 = dτ * capacity3 x s2 := by linarith
  exact mul_left_cancel₀ hd h3

/-- If `weight x 0 < weight x 1`, four sites carry four strictly increasing capacities. -/
theorem capacity_chain (hx : StrictMono x) (h01 : weight x 0 < weight x 1) :
    capacity3 x (0, 0, 0) < capacity3 x (1, 0, 0) ∧ capacity3 x (1, 0, 0) < capacity3 x (1, 1, 0)
      ∧ capacity3 x (1, 1, 0) < capacity3 x (1, 1, 1) := by
  have h0 := weight_pos hx 0
  have h1 := weight_pos hx 1
  simp only [capacity3]
  refine ⟨?_, ?_, ?_⟩
  · exact mul_lt_mul_of_pos_right h01 (mul_pos h0 h0)
  · exact mul_lt_mul_of_pos_left (mul_lt_mul_of_pos_right h01 h0) h1
  · exact mul_lt_mul_of_pos_left (mul_lt_mul_of_pos_left h01 h1) h1

/-- The counterfamily member `ε > 0` is not node-local on a packet carrier with a coupled pair. -/
theorem packet_counterfamily_not_diag (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) (a : Fin n)
    {ε : ℝ} (hε : 0 < ε) :
    ¬ (counterfamilyInertia (capacity3 x) (packetV x m2) ε).IsDiag := by
  have hMpos := capacity3_pos hx
  have hPD : (diagonal (capacity3 x) + ε • packetV x m2).PosDef :=
    (Matrix.PosDef.diagonal hMpos).add ((packetV_posDef hx hm).smul hε)
  have hU : IsUnit (diagonal (capacity3 x) + ε • packetV x m2).det :=
    (Matrix.isUnit_iff_isUnit_det _).1 hPD.isUnit
  have hV : ∃ i j, i ≠ j ∧ packetV x m2 i j ≠ 0 :=
    ⟨(a.castSucc, 0, 0), (a.succ, 0, 0), by simp [Fin.castSucc_lt_succ.ne],
      packetV_axis_ne_zero hx m2 a 0 0⟩
  rw [counterfamilyInertia_isDiag_iff _ _ _ (fun i => (hMpos i).ne') hU hV]
  exact hε.ne'

/-- The member `ε = 0` is node-local: its inertia is the capacity matrix. -/
theorem packet_counterfamily_zero (hx : StrictMono x) (m2 : ℝ) :
    counterfamilyInertia (capacity3 x) (packetV x m2) 0 = diagonal (capacity3 x) :=
  counterfamilyInertia_zero _ _ fun i => (capacity3_pos hx i).ne'

/-- The register over a packet carrier with capacity-matched inertia is positive definite. -/
theorem packet_register_posDef (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) {dτ : ℝ} (hd : 0 < dτ)
    (T : ℕ) :
    (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T) (regCond (edge3 x m2) (capacity3 x) dτ T)
      (regGround (ground3 x m2) dτ T)).PosDef :=
  register_posDef matSrc matTgt (edge3 x m2) (ground3 x m2) (capacity3 x) dτ T (edge3_nonneg hx m2)
    (ground3_pos hx hm) (fun i => (capacity3_pos hx i).le) hd

/-- On a packet carrier, every solution of the Lyapunov equation of the declared register law is
`ℏ P⁻¹`. -/
theorem packet_stationary_unique (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) {dτ : ℝ}
    (hd : 0 < dτ) (hbar : ℝ) (T : ℕ) (X : Matrix (Site n × Fin (T + 1)) (Site n × Fin (T + 1)) ℝ)
    (hX : driftMatrix (regCap (capacity3 x) T) (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
          (regCond (edge3 x m2) (capacity3 x) dτ T) (regGround (ground3 x m2) dτ T)) * X
        + X * (driftMatrix (regCap (capacity3 x) T) (stiffnessMatrix (regSrc matSrc T)
          (regTgt matTgt T) (regCond (edge3 x m2) (capacity3 x) dτ T)
          (regGround (ground3 x m2) dτ T)))ᵀ
        = (2 * hbar) • diagonal (fun p => (regCap (capacity3 x) T p)⁻¹)) :
    X = stationaryCov hbar (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
      (regCond (edge3 x m2) (capacity3 x) dτ T) (regGround (ground3 x m2) dτ T)) :=
  register_stationary_unique matSrc matTgt (edge3 x m2) (ground3 x m2) (capacity3 x) (capacity3 x)
    dτ hbar T (edge3_nonneg hx m2) (ground3_pos hx hm) (fun i => (capacity3_pos hx i).le) hd
    (capacity3_pos hx) X hX

/-- Whitened spatial operator `M^{-1/2} V M^{-1/2}`. -/
def whitenedPacket (x : Fin (n + 3) → ℝ) (m2 : ℝ) : Matrix (Site n) (Site n) ℝ :=
  diagonal (fun s => (√(capacity3 x s))⁻¹) * packetV x m2
    * diagonal (fun s => (√(capacity3 x s))⁻¹)

theorem whitenedPacket_posDef (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) :
    (whitenedPacket x m2).PosDef := by
  have hSu : IsUnit (diagonal (fun s => (√(capacity3 x s))⁻¹)) := by
    rw [isUnit_diagonal, Pi.isUnit_iff]
    intro s
    exact isUnit_iff_ne_zero.2 (inv_ne_zero (Real.sqrt_pos.2 (capacity3_pos hx s)).ne')
  have hSH : (diagonal (fun s => (√(capacity3 x s))⁻¹))ᴴ
      = diagonal (fun s => (√(capacity3 x s))⁻¹) := by
    rw [diagonal_conjTranspose]
    simp
  have h := (packetV_posDef hx hm).conjTranspose_mul_mul_same (mulVec_injective_of_isUnit hSu)
  rw [hSH] at h
  exact h

/-- Every whitened mode of a packet carrier is positive; the per-mode transfer structure, bond
reflection positivity and the `Δτ → 0` frequency limit hold for each mode. -/
theorem packet_modes (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) {dτ : ℝ} (hd : 0 < dτ)
    (i : Site n) :
    let lam := (whitenedPacket_posDef hx hm).1.eigenvalues i
    0 < lam ∧ 0 < transferEig dτ lam ∧ transferEig dτ lam < 1
      ∧ transferEig dτ lam + (transferEig dτ lam)⁻¹ = 2 + dτ ^ 2 * lam
      ∧ (∀ (k : ℕ) (t : Fin k → ℕ) (v : Fin k → ℝ),
          0 ≤ bondGram (equalTimeVar dτ lam) (transferEig dτ lam) t v)
      ∧ Filter.Tendsto (fun d => latticeFreq d lam) (nhdsWithin 0 (Set.Ioi 0))
          (nhds (Real.sqrt lam)) := by
  intro lam
  have hl : 0 < lam := (whitenedPacket_posDef hx hm).eigenvalues_pos i
  exact ⟨hl, transferEig_pos hd hl, transferEig_lt_one hd hl, transferEig_add_inv hd hl,
    fun k t v => chain_bond_rp hd hl t v, latticeFreq_tendsto hl⟩

end PacketInstances

section LadderAndLaw

variable {n : ℕ} {x : Fin (n + 3) → ℝ}

/-- Capacity square roots `M^{1/2}` on the sites. -/
def capSqrt (x : Fin (n + 3) → ℝ) : Matrix (Site n) (Site n) ℝ :=
  diagonal fun s => √(capacity3 x s)

/-- Inverse inertia of the canonical completion, `M⁻¹`. -/
def invInertia (x : Fin (n + 3) → ℝ) : Matrix (Site n) (Site n) ℝ :=
  diagonal fun s => (capacity3 x s)⁻¹

/-- The matrix `G = M^{1/2} A^{1/2} M^{1/2}` of the Gaussian ground state `exp(−qᵀ G q / 2)` of
`½(pᵀM⁻¹p + qᵀVq)`, with `A` the whitened packet operator; its position covariance is `(2G)⁻¹`. -/
def groundPrecision (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) : Matrix (Site n) (Site n) ℝ :=
  capSqrt x * (whitenedPacket_posDef hx hm).1.cfc Real.sqrt * capSqrt x

theorem capSqrt_inv_capSqrt (hx : StrictMono x) : capSqrt x * invInertia x * capSqrt x = 1 := by
  rw [capSqrt, invInertia, diagonal_mul_diagonal, diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext s
  have h := capacity3_pos hx s
  have hs : √(capacity3 x s) * √(capacity3 x s) = capacity3 x s := Real.mul_self_sqrt h.le
  calc √(capacity3 x s) * (capacity3 x s)⁻¹ * √(capacity3 x s)
      = (√(capacity3 x s) * √(capacity3 x s)) * (capacity3 x s)⁻¹ := by ring
    _ = 1 := by rw [hs, mul_inv_cancel₀ h.ne']

theorem capSqrt_whitening (hx : StrictMono x) :
    capSqrt x * diagonal (fun s => (√(capacity3 x s))⁻¹) = 1
      ∧ diagonal (fun s => (√(capacity3 x s))⁻¹) * capSqrt x = 1 := by
  have h : ∀ s, √(capacity3 x s) ≠ 0 := fun s => (Real.sqrt_pos.2 (capacity3_pos hx s)).ne'
  constructor
  · rw [capSqrt, diagonal_mul_diagonal, ← diagonal_one]
    congr 1
    funext s
    exact mul_inv_cancel₀ (h s)
  · rw [capSqrt, diagonal_mul_diagonal, ← diagonal_one]
    congr 1
    funext s
    exact inv_mul_cancel₀ (h s)

/-- On a packet carrier the ground-state relation `G M⁻¹ G = V` holds. -/
theorem groundPrecision_rel (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) :
    groundRel (groundPrecision hx hm) (invInertia x) (packetV x m2) := by
  set hA := whitenedPacket_posDef hx hm with hAdef
  have hsq : hA.1.cfc Real.sqrt * hA.1.cfc Real.sqrt = whitenedPacket x m2 :=
    (cfc_mul hA.1 _ _).trans ((cfc_congr hA.1 fun i => Real.mul_self_sqrt
      (hA.eigenvalues_pos i).le).trans (cfc_id hA.1))
  unfold groundRel groundPrecision
  calc capSqrt x * hA.1.cfc Real.sqrt * capSqrt x * invInertia x
        * (capSqrt x * hA.1.cfc Real.sqrt * capSqrt x)
      = capSqrt x * hA.1.cfc Real.sqrt * (capSqrt x * invInertia x * capSqrt x)
          * hA.1.cfc Real.sqrt * capSqrt x := by
        simp only [Matrix.mul_assoc]
    _ = capSqrt x * (hA.1.cfc Real.sqrt * hA.1.cfc Real.sqrt) * capSqrt x := by
        rw [capSqrt_inv_capSqrt hx, Matrix.mul_one]
        simp only [Matrix.mul_assoc]
    _ = capSqrt x * whitenedPacket x m2 * capSqrt x := by rw [hsq]
    _ = (capSqrt x * diagonal (fun s => (√(capacity3 x s))⁻¹)) * packetV x m2
          * (diagonal (fun s => (√(capacity3 x s))⁻¹) * capSqrt x) := by
        simp only [whitenedPacket, Matrix.mul_assoc]
    _ = packetV x m2 := by
        rw [(capSqrt_whitening hx).1, (capSqrt_whitening hx).2, Matrix.one_mul, Matrix.mul_one]

theorem groundPrecision_isUnit (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) :
    IsUnit (groundPrecision hx hm).det := by
  have hA := whitenedPacket_posDef hx hm
  have hQ : IsUnit (hA.1.cfc Real.sqrt) :=
    (cfc_posDef hA.1 fun i => Real.sqrt_pos.2 (hA.eigenvalues_pos i)).isUnit
  have hS : IsUnit (capSqrt x) := by
    rw [capSqrt, isUnit_diagonal, Pi.isUnit_iff]
    exact fun s => isUnit_iff_ne_zero.2 (Real.sqrt_pos.2 (capacity3_pos hx s)).ne'
  rw [← Matrix.isUnit_iff_isUnit_det]
  exact (hS.mul hQ).mul hS

/-- For the ground-state matrix `G` of the canonical completion (`groundPrecision`, built from `M`
and `V`), `G T G = V` holds exactly when `T = M⁻¹`: a consistency check of the ladder on the
carrier, not an identification of the inertia from independent vacuum data. -/
theorem packet_inverse_inertia (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2)
    (T : Matrix (Site n) (Site n) ℝ) :
    groundRel (groundPrecision hx hm) T (packetV x m2) ↔ T = invInertia x := by
  constructor
  · intro h
    exact inverse_inertia_unique _ _ _ _ (groundPrecision_isUnit hx hm) h (groundPrecision_rel hx hm)
  · rintro rfl
    exact groundPrecision_rel hx hm

theorem packet_inverse_inertia_formula (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) :
    invInertia x = (groundPrecision hx hm)⁻¹ * packetV x m2 * (groundPrecision hx hm)⁻¹ :=
  (inverse_inertia_formula _ _ _ (groundPrecision_isUnit hx hm) (groundPrecision_rel hx hm)).1

/-- The configuration marginal of the Gibbs precision does not see the kinetic block. -/
theorem packet_config_cov (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) {β : ℝ} (hβ : β ≠ 0) :
    ((fromBlocks (β • packetV x m2) 0 0 (β • invInertia x))⁻¹).toBlocks₁₁ = (β • packetV x m2)⁻¹ :=
  config_cov_indep β _ _ hβ ((Matrix.isUnit_iff_isUnit_det _).1 (packetV_posDef hx hm).isUnit)
    (by
      rw [invInertia, det_diagonal]
      exact IsUnit.mk0 _ (Finset.prod_ne_zero_iff.2 fun s _ => inv_ne_zero (capacity3_pos hx s).ne'))

/-- The two-time reader on a packet carrier: the general reader formula applies to the stationary
chain whose interior precision is that of the register (`R = r(A)`, `C0 = c0(A)`), whose
innovation covariance is invertible. -/
theorem packet_markov_reader (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) {dτ : ℝ} (hd : 0 < dτ) :
    (equalTimeCov (whitenedPacket_posDef hx hm) dτ
        * (transferMat (whitenedPacket_posDef hx hm) dτ
          * equalTimeCov (whitenedPacket_posDef hx hm) dτ)⁻¹
        * equalTimeCov (whitenedPacket_posDef hx hm) dτ
        - transferMat (whitenedPacket_posDef hx hm) dτ * equalTimeCov (whitenedPacket_posDef hx hm) dτ)⁻¹
      = (equalTimeCov (whitenedPacket_posDef hx hm) dτ
          - transferMat (whitenedPacket_posDef hx hm) dτ * equalTimeCov (whitenedPacket_posDef hx hm) dτ
            * (transferMat (whitenedPacket_posDef hx hm) dτ)ᵀ)⁻¹
        * transferMat (whitenedPacket_posDef hx hm) dτ
      ∧ IsUnit (equalTimeCov (whitenedPacket_posDef hx hm) dτ
          - transferMat (whitenedPacket_posDef hx hm) dτ * equalTimeCov (whitenedPacket_posDef hx hm) dτ
            * (transferMat (whitenedPacket_posDef hx hm) dτ)ᵀ).det := by
  refine ⟨markov_reader _ _
    ((Matrix.isUnit_iff_isUnit_det _).1 (equalTimeCov_posDef _ hd).isUnit)
    ((Matrix.isUnit_iff_isUnit_det _).1 (transferMat_posDef _ hd).isUnit) ?_, innovation_isUnit _ hd⟩
  rw [transferMat_transpose]
  exact transferMat_comm _

/-- The declared law on a packet carrier, for the register stiffness with capacity-matched
inertia: `ℏ P⁻¹` solves the Lyapunov equation, `B Σ` is symmetric, the matched-noise equivalence
holds, and `ℏ P⁻¹` also solves the Lyapunov equation for the unit mobility `Γ = 1`. -/
theorem packet_langevin (hx : StrictMono x) {m2 : ℝ} (hm : 0 < m2) {dτ : ℝ} (hd : 0 < dτ)
    (hbar : ℝ) (T : ℕ) :
    let P := stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
      (regCond (edge3 x m2) (capacity3 x) dτ T) (regGround (ground3 x m2) dτ T)
    (driftMatrix (regCap (capacity3 x) T) P * stationaryCov hbar P
        + stationaryCov hbar P * (driftMatrix (regCap (capacity3 x) T) P)ᵀ
        = (2 * hbar) • diagonal (fun p => (regCap (capacity3 x) T p)⁻¹))
      ∧ (driftMatrix (regCap (capacity3 x) T) P * stationaryCov hbar P)ᵀ
        = driftMatrix (regCap (capacity3 x) T) P * stationaryCov hbar P
      ∧ (∀ D, driftMatrix (regCap (capacity3 x) T) P * stationaryCov hbar P
          + stationaryCov hbar P * (driftMatrix (regCap (capacity3 x) T) P)ᵀ = (2 : ℝ) • D
          ↔ D = hbar • diagonal (fun p => (regCap (capacity3 x) T p)⁻¹))
      ∧ ((1 : Matrix (Site n × Fin (T + 1)) (Site n × Fin (T + 1)) ℝ) * P) * stationaryCov hbar P
          + stationaryCov hbar P * ((1 : Matrix (Site n × Fin (T + 1)) (Site n × Fin (T + 1)) ℝ) * P)ᵀ
          = (2 * hbar) • (1 : Matrix (Site n × Fin (T + 1)) (Site n × Fin (T + 1)) ℝ) := by
  intro P
  have hP : Pᵀ = P := stiffnessMatrix_transpose _ _ _ _
  have hPU : IsUnit P.det :=
    (Matrix.isUnit_iff_isUnit_det _).1 (packet_register_posDef hx hm hd T).isUnit
  exact ⟨lyapunov_stationary _ P hbar hP hPU, detailed_balance _ P hbar hPU,
    fun D => matched_noise_iff _ P D hbar hP hPU, mobility_independent 1 P hbar hP hPU transpose_one⟩

/-- On a packet carrier the seam-and-grounding sum of one slab is the quadratic form `fᵀVf` of the
packet operator (in the register action it carries the factor `Δτ/2`). -/
theorem packet_slab_quadratic (x : Fin (n + 3) → ℝ) (m2 : ℝ) (f : Site n → ℝ) :
    ∑ e, edge3 x m2 e * (f (matSrc e) - f (matTgt e)) ^ 2 + ∑ k, ground3 x m2 k * f k ^ 2
      = f ⬝ᵥ (packetV x m2 *ᵥ f) := by
  rw [← stiffness_quadratic matSrc matTgt (edge3 x m2) (ground3 x m2) f,
    ← stiffnessMatrix_packet x m2, stiffnessMatrix_mulVec]
  rfl

/-- The register action on a packet carrier: time kinetic term plus `Δτ fᵀVf / 2` per slab. -/
theorem packet_registerAction (x : Fin (n + 3) → ℝ) (m2 : ℝ) (J : Site n → ℝ) (dτ : ℝ) (T : ℕ)
    (f : Site n × Fin (T + 1) → ℝ) :
    registerAction matSrc matTgt (edge3 x m2) (ground3 x m2) J dτ T f
      = ∑ t : Fin T, ∑ i, J i / dτ * (f (i, t.succ) - f (i, t.castSucc)) ^ 2 / 2
        + ∑ t : Fin (T + 1), dτ * ((fun i => f (i, t)) ⬝ᵥ (packetV x m2 *ᵥ fun i => f (i, t))) / 2 := by
  rw [registerAction_eq]
  congr 1
  refine Finset.sum_congr rfl fun t _ => ?_
  exact congrArg (fun z => dτ * z / 2) (packet_slab_quadratic x m2 fun i => f (i, t))

end LadderAndLaw

section GoldenInstances

/-- The premises of the register theorems hold on the q = 3 and q = 13 golden carriers. -/
theorem golden_premises :
    (∀ s, 0 < capacity3 golden3 s) ∧ (∀ e, 0 ≤ edge3 golden3 1 e)
      ∧ (∀ s, 0 < ground3 golden3 1 s) ∧ (∀ e : MatSeam (Site 1), matSrc e ≠ matTgt e)
      ∧ (∀ s, 0 < capacity3 golden13 s) ∧ (∀ e, 0 ≤ edge3 golden13 1 e)
      ∧ (∀ s, 0 < ground3 golden13 1 s) ∧ (∀ e : MatSeam (Site 11), matSrc e ≠ matTgt e) :=
  ⟨capacity3_pos golden3_strictMono, edge3_nonneg golden3_strictMono 1,
    ground3_pos golden3_strictMono one_pos, edge3_src_ne_tgt,
    capacity3_pos golden13_strictMono, edge3_nonneg golden13_strictMono 1,
    ground3_pos golden13_strictMono one_pos, edge3_src_ne_tgt⟩

/-- Positive instance: capacity-matched inertia gives `2/Δτ` at every site of both carriers. -/
theorem golden_matched_rate {dτ : ℝ} (hd : dτ ≠ 0) (T : ℕ) (t : Fin T) :
    (∀ s : Site 1, seamRate (regCap (capacity3 golden3) T)
        (regCond (edge3 golden3 1) (capacity3 golden3) dτ T (.inr (s, t)))
        (regSrc matSrc T (.inr (s, t))) (regTgt matTgt T (.inr (s, t))) = 2 / dτ)
      ∧ ∀ s : Site 11, seamRate (regCap (capacity3 golden13) T)
        (regCond (edge3 golden13 1) (capacity3 golden13) dτ T (.inr (s, t)))
        (regSrc matSrc T (.inr (s, t))) (regTgt matTgt T (.inr (s, t))) = 2 / dτ :=
  ⟨fun s => packet_matched_rate golden3_strictMono 1 hd T s t,
    fun s => packet_matched_rate golden13_strictMono 1 hd T s t⟩

/-- Four distinct capacities on each carrier. -/
theorem golden_capacities :
    (capacity3 golden3 (0, 0, 0) < capacity3 golden3 (1, 0, 0)
      ∧ capacity3 golden3 (1, 0, 0) < capacity3 golden3 (1, 1, 0)
      ∧ capacity3 golden3 (1, 1, 0) < capacity3 golden3 (1, 1, 1))
    ∧ (capacity3 golden13 (1, 1, 1) < capacity3 golden13 (1, 1, 0)
      ∧ capacity3 golden13 (1, 1, 0) < capacity3 golden13 (1, 0, 0)
      ∧ capacity3 golden13 (1, 0, 0) < capacity3 golden13 (0, 0, 0)) := by
  refine ⟨capacity_chain golden3_strictMono golden3_weight_lt, ?_⟩
  have h0 := weight_pos golden13_strictMono 0
  have h1 := weight_pos golden13_strictMono 1
  have h01 := golden13_weight_gt
  simp only [capacity3]
  refine ⟨?_, ?_, ?_⟩
  · exact mul_lt_mul_of_pos_left (mul_lt_mul_of_pos_left h01 h1) h1
  · exact mul_lt_mul_of_pos_left (mul_lt_mul_of_pos_right h01 h0) h1
  · exact mul_lt_mul_of_pos_right h01 (mul_pos h0 h0)

/-- Negative instance: unit inertia gives different time rates at two sites of each carrier. -/
theorem golden_unit_rate_ne {dτ : ℝ} (hd : dτ ≠ 0) (T : ℕ) (t : Fin T) :
    seamRate (regCap (capacity3 golden3) T)
        (regCond (edge3 golden3 1) (fun _ => 1) dτ T (.inr ((0, 0, 0), t)))
        (regSrc matSrc T (.inr ((0, 0, 0), t))) (regTgt matTgt T (.inr ((0, 0, 0), t)))
      ≠ seamRate (regCap (capacity3 golden3) T)
        (regCond (edge3 golden3 1) (fun _ => 1) dτ T (.inr ((1, 0, 0), t)))
        (regSrc matSrc T (.inr ((1, 0, 0), t))) (regTgt matTgt T (.inr ((1, 0, 0), t)))
    ∧ seamRate (regCap (capacity3 golden13) T)
        (regCond (edge3 golden13 1) (fun _ => 1) dτ T (.inr ((0, 0, 0), t)))
        (regSrc matSrc T (.inr ((0, 0, 0), t))) (regTgt matTgt T (.inr ((0, 0, 0), t)))
      ≠ seamRate (regCap (capacity3 golden13) T)
        (regCond (edge3 golden13 1) (fun _ => 1) dτ T (.inr ((1, 0, 0), t)))
        (regSrc matSrc T (.inr ((1, 0, 0), t))) (regTgt matTgt T (.inr ((1, 0, 0), t))) :=
  ⟨packet_unit_rate_ne golden3_strictMono 1 hd T t golden_capacities.1.1.ne,
    packet_unit_rate_ne golden13_strictMono 1 hd T t golden_capacities.2.2.2.ne'⟩

/-- Counterfamily on both carriers: `ε = 1/4` is not node-local, `ε = 0` is. -/
theorem golden_counterfamily :
    ¬ (counterfamilyInertia (capacity3 golden3) (packetV golden3 1) (1 / 4)).IsDiag
      ∧ (counterfamilyInertia (capacity3 golden3) (packetV golden3 1) 0).IsDiag
      ∧ ¬ (counterfamilyInertia (capacity3 golden13) (packetV golden13 1) (1 / 4)).IsDiag
      ∧ (counterfamilyInertia (capacity3 golden13) (packetV golden13 1) 0).IsDiag := by
  refine ⟨packet_counterfamily_not_diag golden3_strictMono one_pos 0 (by norm_num), ?_,
    packet_counterfamily_not_diag golden13_strictMono one_pos 0 (by norm_num), ?_⟩
  · rw [packet_counterfamily_zero golden3_strictMono 1]
    exact isDiag_diagonal _
  · rw [packet_counterfamily_zero golden13_strictMono 1]
    exact isDiag_diagonal _

/-- The registers over both carriers are positive definite and their modes are positive. -/
theorem golden_register {dτ : ℝ} (hd : 0 < dτ) (T : ℕ) :
    (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
      (regCond (edge3 golden3 1) (capacity3 golden3) dτ T) (regGround (ground3 golden3 1) dτ T)).PosDef
    ∧ (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
      (regCond (edge3 golden13 1) (capacity3 golden13) dτ T)
        (regGround (ground3 golden13 1) dτ T)).PosDef
    ∧ (∀ i, 0 < (whitenedPacket_posDef golden3_strictMono one_pos).1.eigenvalues i)
    ∧ (∀ i, 0 < (whitenedPacket_posDef golden13_strictMono one_pos).1.eigenvalues i) :=
  ⟨packet_register_posDef golden3_strictMono one_pos hd T,
    packet_register_posDef golden13_strictMono one_pos hd T,
    fun i => (packet_modes golden3_strictMono one_pos hd i).1,
    fun i => (packet_modes golden13_strictMono one_pos hd i).1⟩

/-- The matrix transfer structure fires on both carriers: with `A` the whitened packet operator,
`R + R⁻¹ = 2 + Δτ² A`, `C0 ≻ 0`, `1 − R ≻ 0` and block bond reflection positivity. -/
theorem golden_transfer_lift {dτ : ℝ} (hd : 0 < dτ) :
    (transferMat (whitenedPacket_posDef golden3_strictMono one_pos) dτ
        + (transferMat (whitenedPacket_posDef golden3_strictMono one_pos) dτ)⁻¹
        = (2 : ℝ) • (1 : Matrix (Site 1) (Site 1) ℝ) + dτ ^ 2 • whitenedPacket golden3 1
      ∧ (equalTimeCov (whitenedPacket_posDef golden3_strictMono one_pos) dτ).PosDef
      ∧ (1 - transferMat (whitenedPacket_posDef golden3_strictMono one_pos) dτ).PosDef
      ∧ ∀ {m : ℕ} (t : Fin m → ℕ) (v : Fin m → Site 1 → ℝ),
          0 ≤ ∑ a, ∑ b, v a ⬝ᵥ ((equalTimeCov (whitenedPacket_posDef golden3_strictMono one_pos) dτ
            * transferMat (whitenedPacket_posDef golden3_strictMono one_pos) dτ ^ (t a + t b + 1))
              *ᵥ v b))
    ∧ (transferMat (whitenedPacket_posDef golden13_strictMono one_pos) dτ
        + (transferMat (whitenedPacket_posDef golden13_strictMono one_pos) dτ)⁻¹
        = (2 : ℝ) • (1 : Matrix (Site 11) (Site 11) ℝ) + dτ ^ 2 • whitenedPacket golden13 1
      ∧ (equalTimeCov (whitenedPacket_posDef golden13_strictMono one_pos) dτ).PosDef
      ∧ (1 - transferMat (whitenedPacket_posDef golden13_strictMono one_pos) dτ).PosDef
      ∧ ∀ {m : ℕ} (t : Fin m → ℕ) (v : Fin m → Site 11 → ℝ),
          0 ≤ ∑ a, ∑ b, v a ⬝ᵥ ((equalTimeCov (whitenedPacket_posDef golden13_strictMono one_pos) dτ
            * transferMat (whitenedPacket_posDef golden13_strictMono one_pos) dτ ^ (t a + t b + 1))
              *ᵥ v b)) :=
  ⟨⟨transferMat_add_inv _ hd, equalTimeCov_posDef _ hd, one_sub_transferMat_posDef _ hd,
      fun t v => block_bond_rp _ hd t v⟩,
    ⟨transferMat_add_inv _ hd, equalTimeCov_posDef _ hd, one_sub_transferMat_posDef _ hd,
      fun t v => block_bond_rp _ hd t v⟩⟩

theorem golden3_posDef : (whitenedPacket golden3 1).PosDef :=
  whitenedPacket_posDef golden3_strictMono one_pos

theorem golden13_posDef : (whitenedPacket golden13 1).PosDef :=
  whitenedPacket_posDef golden13_strictMono one_pos

/-- Identification ladder on both carriers: for the ground-state matrix `G` of the canonical
completion and the potential, `G K G = V` holds exactly for `K = M⁻¹` (a consistency check, as in
`packet_inverse_inertia`), and the configuration marginal does not see the kinetic block. -/
theorem golden_ladder :
    (∀ K, groundRel (groundPrecision golden3_strictMono one_pos) K (packetV golden3 1)
        ↔ K = invInertia golden3)
      ∧ (∀ K, groundRel (groundPrecision golden13_strictMono one_pos) K (packetV golden13 1)
        ↔ K = invInertia golden13)
      ∧ ((fromBlocks ((1 : ℝ) • packetV golden3 1) 0 0 ((1 : ℝ) • invInertia golden3))⁻¹).toBlocks₁₁
        = ((1 : ℝ) • packetV golden3 1)⁻¹
      ∧ ((fromBlocks ((1 : ℝ) • packetV golden13 1) 0 0 ((1 : ℝ) • invInertia golden13))⁻¹).toBlocks₁₁
        = ((1 : ℝ) • packetV golden13 1)⁻¹ :=
  ⟨packet_inverse_inertia golden3_strictMono one_pos, packet_inverse_inertia golden13_strictMono one_pos,
    packet_config_cov golden3_strictMono one_pos one_ne_zero,
    packet_config_cov golden13_strictMono one_pos one_ne_zero⟩

/-- The general two-time reader applies on both carriers. -/
theorem golden_reader {dτ : ℝ} (hd : 0 < dτ) :
    (equalTimeCov golden3_posDef dτ * (transferMat golden3_posDef dτ * equalTimeCov golden3_posDef dτ)⁻¹
        * equalTimeCov golden3_posDef dτ - transferMat golden3_posDef dτ * equalTimeCov golden3_posDef dτ)⁻¹
      = (equalTimeCov golden3_posDef dτ - transferMat golden3_posDef dτ * equalTimeCov golden3_posDef dτ
          * (transferMat golden3_posDef dτ)ᵀ)⁻¹ * transferMat golden3_posDef dτ
    ∧ (equalTimeCov golden13_posDef dτ
        * (transferMat golden13_posDef dτ * equalTimeCov golden13_posDef dτ)⁻¹
        * equalTimeCov golden13_posDef dτ - transferMat golden13_posDef dτ * equalTimeCov golden13_posDef dτ)⁻¹
      = (equalTimeCov golden13_posDef dτ - transferMat golden13_posDef dτ * equalTimeCov golden13_posDef dτ
          * (transferMat golden13_posDef dτ)ᵀ)⁻¹ * transferMat golden13_posDef dτ :=
  ⟨(packet_markov_reader golden3_strictMono one_pos hd).1,
    (packet_markov_reader golden13_strictMono one_pos hd).1⟩

theorem golden13_matched_whitened :
    (jInvSqrt (capacity3 golden13) * packetV golden13 1 * jInvSqrt (capacity3 golden13)).PosDef :=
  jWhitened_posDef (capacity3_pos golden13_strictMono) (packetV_posDef golden13_strictMono one_pos)

theorem golden13_unit_whitened :
    (jInvSqrt (fun _ : Site 11 => (1 : ℝ)) * packetV golden13 1
      * jInvSqrt (fun _ : Site 11 => (1 : ℝ))).PosDef :=
  jWhitened_posDef (fun _ => one_pos) (packetV_posDef golden13_strictMono one_pos)

/-- On the q = 13 carrier the two-time reader of the stationary chain returns the time-seam
coefficient, both for capacity-matched inertia (`J = M`) and for unit inertia (`J = 1`); since the
capacities are not all equal (`golden_capacities`), the two readings differ. -/
theorem golden13_history_reader {dτ : ℝ} (hd : 0 < dτ) :
    dτ • (lag0 golden13_matched_whitened dτ * (lag1 golden13_matched_whitened dτ)⁻¹
        * lag0 golden13_matched_whitened dτ - lag1 golden13_matched_whitened dτ)⁻¹
      = diagonal (capacity3 golden13)
    ∧ dτ • (lag0 golden13_unit_whitened dτ * (lag1 golden13_unit_whitened dτ)⁻¹
        * lag0 golden13_unit_whitened dτ - lag1 golden13_unit_whitened dτ)⁻¹
      = diagonal (fun _ : Site 11 => (1 : ℝ)) :=
  ⟨history_reader_general (capacity3_pos golden13_strictMono) _ hd,
    history_reader_general (fun _ => one_pos) _ hd⟩

/-- The declared law's stationary covariance solves the Lyapunov equation on both carriers. -/
theorem golden_law {dτ : ℝ} (hd : 0 < dτ) (hbar : ℝ) (T : ℕ) :
    driftMatrix (regCap (capacity3 golden3) T) (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
          (regCond (edge3 golden3 1) (capacity3 golden3) dτ T) (regGround (ground3 golden3 1) dτ T))
        * stationaryCov hbar (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
          (regCond (edge3 golden3 1) (capacity3 golden3) dτ T) (regGround (ground3 golden3 1) dτ T))
      + stationaryCov hbar (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
          (regCond (edge3 golden3 1) (capacity3 golden3) dτ T) (regGround (ground3 golden3 1) dτ T))
        * (driftMatrix (regCap (capacity3 golden3) T) (stiffnessMatrix (regSrc matSrc T)
          (regTgt matTgt T) (regCond (edge3 golden3 1) (capacity3 golden3) dτ T)
          (regGround (ground3 golden3 1) dτ T)))ᵀ
      = (2 * hbar) • diagonal (fun p => (regCap (capacity3 golden3) T p)⁻¹)
    ∧ driftMatrix (regCap (capacity3 golden13) T) (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
          (regCond (edge3 golden13 1) (capacity3 golden13) dτ T) (regGround (ground3 golden13 1) dτ T))
        * stationaryCov hbar (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
          (regCond (edge3 golden13 1) (capacity3 golden13) dτ T) (regGround (ground3 golden13 1) dτ T))
      + stationaryCov hbar (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
          (regCond (edge3 golden13 1) (capacity3 golden13) dτ T) (regGround (ground3 golden13 1) dτ T))
        * (driftMatrix (regCap (capacity3 golden13) T) (stiffnessMatrix (regSrc matSrc T)
          (regTgt matTgt T) (regCond (edge3 golden13 1) (capacity3 golden13) dτ T)
          (regGround (ground3 golden13 1) dτ T)))ᵀ
      = (2 * hbar) • diagonal (fun p => (regCap (capacity3 golden13) T p)⁻¹) :=
  ⟨(packet_langevin golden3_strictMono one_pos hd hbar T).1,
    (packet_langevin golden13_strictMono one_pos hd hbar T).1⟩

end GoldenInstances

section SuccessorInstances

/-- The whitened packet operator couples distinct sites, so it is not diagonal. -/
theorem whitenedPacket_not_diag {n : ℕ} {x : Fin (n + 3) → ℝ} (hx : StrictMono x) (m2 : ℝ)
    (a : Fin n) : ¬ (whitenedPacket x m2).IsDiag := by
  intro h
  have hne : ((a.castSucc, 0, 0) : Site n) ≠ (a.succ, 0, 0) := by
    simp [Fin.castSucc_lt_succ.ne]
  have h0 := h hne
  simp only [whitenedPacket, diagonal_mul, mul_diagonal] at h0
  have hs1 := (Real.sqrt_pos.2 (capacity3_pos hx (a.castSucc, 0, 0))).ne'
  have hs2 := (Real.sqrt_pos.2 (capacity3_pos hx (a.succ, 0, 0))).ne'
  rcases mul_eq_zero.1 h0 with h1 | h1
  · rcases mul_eq_zero.1 h1 with h2 | h2
    · exact inv_ne_zero hs1 h2
    · exact packetV_axis_ne_zero hx m2 a 0 0 h2
  · exact inv_ne_zero hs2 h1

/-- Repair admissibility on both carriers: the register with capacity-matched inertia is a repair
stiffness at every step, the counterfamily member `ε = 1/4` is not admissible for all small steps,
and admissibility for all small steps together with the uniform time rate `2/Δτ` holds exactly for
the capacity matrix. -/
theorem golden_repair (T : ℕ) (hT : 0 < T) {dτ0 : ℝ} (hd0 : dτ0 ≠ 0) :
    (∀ {dτ : ℝ}, 0 < dτ → IsRepairStiffness
        (historyPrecision (diagonal (capacity3 golden3)) (packetV golden3 1) dτ T))
      ∧ ¬ (∀ᶠ dτ in 𝓝[>] (0 : ℝ), IsRepairStiffness (historyPrecision
          (counterfamilyInertia (capacity3 golden3) (packetV golden3 1) (1 / 4)) (packetV golden3 1) dτ T))
      ∧ (∀ J : Matrix (Site 1) (Site 1) ℝ,
          ((∀ᶠ dτ in 𝓝[>] (0 : ℝ), IsRepairStiffness (historyPrecision J (packetV golden3 1) dτ T))
            ∧ ∀ i, J i i / dτ0 * (1 / capacity3 golden3 i + 1 / capacity3 golden3 i) = 2 / dτ0)
          ↔ J = diagonal (capacity3 golden3))
      ∧ (∀ {dτ : ℝ}, 0 < dτ → IsRepairStiffness
        (historyPrecision (diagonal (capacity3 golden13)) (packetV golden13 1) dτ T))
      ∧ ¬ (∀ᶠ dτ in 𝓝[>] (0 : ℝ), IsRepairStiffness (historyPrecision
          (counterfamilyInertia (capacity3 golden13) (packetV golden13 1) (1 / 4))
          (packetV golden13 1) dτ T))
      ∧ (∀ J : Matrix (Site 11) (Site 11) ℝ,
          ((∀ᶠ dτ in 𝓝[>] (0 : ℝ), IsRepairStiffness (historyPrecision J (packetV golden13 1) dτ T))
            ∧ ∀ i, J i i / dτ0 * (1 / capacity3 golden13 i + 1 / capacity3 golden13 i) = 2 / dτ0)
          ↔ J = diagonal (capacity3 golden13)) :=
  ⟨fun hd => historyPrecision_diagonal_repair (packetV_offdiag_nonpos golden3_strictMono 1)
      (fun i => (capacity3_pos golden3_strictMono i).le) hd,
    fun h => golden_counterfamily.1 (eventually_repair_isDiag hT _ _ h),
    fun J => eventually_repair_clock_iff hT _ J _ (capacity3_pos golden3_strictMono)
      (packetV_offdiag_nonpos golden3_strictMono 1) hd0,
    fun hd => historyPrecision_diagonal_repair (packetV_offdiag_nonpos golden13_strictMono 1)
      (fun i => (capacity3_pos golden13_strictMono i).le) hd,
    fun h => golden_counterfamily.2.2.1 (eventually_repair_isDiag hT _ _ h),
    fun J => eventually_repair_clock_iff hT _ J _ (capacity3_pos golden13_strictMono)
      (packetV_offdiag_nonpos golden13_strictMono 1) hd0⟩

/-- Two-time correlations on both carriers: at `ε = 1/4` the reader of the counterfamily chain
returns the counterfamily inertia, which is not node-local; at `ε = 0` the reader is node-local. -/
theorem golden_counterfamily_reader {dτ : ℝ} (hd : 0 < dτ) :
    cfReader (capacity3_pos golden3_strictMono) (packetV_posDef golden3_strictMono one_pos)
        (by norm_num : (0 : ℝ) ≤ 1 / 4) dτ
        = counterfamilyInertia (capacity3 golden3) (packetV golden3 1) (1 / 4)
      ∧ ¬ (cfReader (capacity3_pos golden3_strictMono) (packetV_posDef golden3_strictMono one_pos)
        (by norm_num : (0 : ℝ) ≤ 1 / 4) dτ).IsDiag
      ∧ (cfReader (capacity3_pos golden3_strictMono) (packetV_posDef golden3_strictMono one_pos)
        le_rfl dτ).IsDiag
      ∧ cfReader (capacity3_pos golden13_strictMono) (packetV_posDef golden13_strictMono one_pos)
        (by norm_num : (0 : ℝ) ≤ 1 / 4) dτ
        = counterfamilyInertia (capacity3 golden13) (packetV golden13 1) (1 / 4)
      ∧ ¬ (cfReader (capacity3_pos golden13_strictMono) (packetV_posDef golden13_strictMono one_pos)
        (by norm_num : (0 : ℝ) ≤ 1 / 4) dτ).IsDiag
      ∧ (cfReader (capacity3_pos golden13_strictMono) (packetV_posDef golden13_strictMono one_pos)
        le_rfl dτ).IsDiag := by
  refine ⟨counterfamily_history_reader _ _ _ hd, ?_, ?_, counterfamily_history_reader _ _ _ hd, ?_, ?_⟩
  · rw [counterfamily_history_reader _ _ _ hd]
    exact golden_counterfamily.1
  · rw [counterfamily_history_reader _ _ _ hd]
    exact golden_counterfamily.2.1
  · rw [counterfamily_history_reader _ _ _ hd]
    exact golden_counterfamily.2.2.1
  · rw [counterfamily_history_reader _ _ _ hd]
    exact golden_counterfamily.2.2.2

/-- The spatial stiffness of the q = 3 carrier whitened by `M^{-1/2}` is positive definite (the
premise `hA` of `declared_register_window` for capacity-matched inertia). -/
theorem golden3_register_whitened :
    (jInvSqrt (capacity3 golden3) * stiffnessMatrix matSrc matTgt (edge3 golden3 1) (ground3 golden3 1)
      * jInvSqrt (capacity3 golden3)).PosDef := by
  rw [stiffnessMatrix_packet]
  exact jWhitened_posDef (capacity3_pos golden3_strictMono) (packetV_posDef golden3_strictMono one_pos)

/-- The spatial stiffness of the q = 13 carrier whitened by `M^{-1/2}` is positive definite (the
premise `hA` of `declared_register_window` for capacity-matched inertia). -/
theorem golden13_register_whitened :
    (jInvSqrt (capacity3 golden13)
      * stiffnessMatrix matSrc matTgt (edge3 golden13 1) (ground3 golden13 1)
      * jInvSqrt (capacity3 golden13)).PosDef := by
  rw [stiffnessMatrix_packet]
  exact jWhitened_posDef (capacity3_pos golden13_strictMono) (packetV_posDef golden13_strictMono one_pos)

/-- Window on both carriers: the corrected register precision inverts the stationary window
covariance, no node-local term added equally at both end slabs does, and in original
coordinates the declared register plus the end-slab correction inverts the stationary window
covariance. -/
theorem golden_window {dτ : ℝ} (hd : 0 < dτ) (T : ℕ) :
    windowPrec golden3_posDef dτ T * windowCov golden3_posDef dτ T = 1
      ∧ (∀ D : Matrix (Site 1) (Site 1) ℝ, D.IsDiag →
          (historyPrecision 1 (whitenedPacket golden3 1) dτ T + endSlabs T D)
            * windowCov golden3_posDef dτ T ≠ 1)
      ∧ (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
            (regCond (edge3 golden3 1) (capacity3 golden3) dτ T) (regGround (ground3 golden3 1) dτ T)
          + endSlabs T (jSqrt (capacity3 golden3)
            * (dτ⁻¹ • (1 - transferMat golden3_register_whitened dτ)) * jSqrt (capacity3 golden3)))
        * ((jInvSqrt (capacity3 golden3) ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
          * windowCov golden3_register_whitened dτ T
          * (jInvSqrt (capacity3 golden3) ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))) = 1
      ∧ windowPrec golden13_posDef dτ T * windowCov golden13_posDef dτ T = 1
      ∧ (∀ D : Matrix (Site 11) (Site 11) ℝ, D.IsDiag →
          (historyPrecision 1 (whitenedPacket golden13 1) dτ T + endSlabs T D)
            * windowCov golden13_posDef dτ T ≠ 1)
      ∧ (stiffnessMatrix (regSrc matSrc T) (regTgt matTgt T)
            (regCond (edge3 golden13 1) (capacity3 golden13) dτ T)
            (regGround (ground3 golden13 1) dτ T)
          + endSlabs T (jSqrt (capacity3 golden13)
            * (dτ⁻¹ • (1 - transferMat golden13_register_whitened dτ)) * jSqrt (capacity3 golden13)))
        * ((jInvSqrt (capacity3 golden13) ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))
          * windowCov golden13_register_whitened dτ T
          * (jInvSqrt (capacity3 golden13) ⊗ₖ (1 : Matrix (Fin (T + 1)) (Fin (T + 1)) ℝ))) = 1 :=
  ⟨stationary_window _ hd T,
    fun D hD => no_nodeLocal_window_correction _ hd T (whitenedPacket_not_diag golden3_strictMono 1 0)
      D hD,
    declared_register_window matSrc matTgt _ _ _ (capacity3_pos golden3_strictMono) _ hd T,
    stationary_window _ hd T,
    fun D hD => no_nodeLocal_window_correction _ hd T (whitenedPacket_not_diag golden13_strictMono 1 0)
      D hD,
    declared_register_window matSrc matTgt _ _ _ (capacity3_pos golden13_strictMono) _ hd T⟩

end SuccessorInstances

end

#print axioms golden3_sorted_orbit
#print axioms whitenedPacket_posDef
#print axioms golden13_sorted_orbit
#print axioms golden3_weight0
#print axioms golden3_weight1
#print axioms golden13_weight0
#print axioms golden13_weight1
#print axioms golden3_axis_coupling
#print axioms packet_register_drift
#print axioms packet_unit_rate_ne
#print axioms packet_counterfamily_not_diag
#print axioms packet_stationary_unique
#print axioms packet_modes
#print axioms golden_premises
#print axioms golden_matched_rate
#print axioms golden_capacities
#print axioms golden_unit_rate_ne
#print axioms golden_counterfamily
#print axioms golden_register
#print axioms golden_transfer_lift
#print axioms groundPrecision_rel
#print axioms packet_inverse_inertia
#print axioms packet_inverse_inertia_formula
#print axioms packet_config_cov
#print axioms packet_markov_reader
#print axioms packet_langevin
#print axioms packet_slab_quadratic
#print axioms packet_registerAction
#print axioms golden_ladder
#print axioms golden_reader
#print axioms golden_law
#print axioms golden13_history_reader
#print axioms whitenedPacket_not_diag
#print axioms golden_repair
#print axioms golden_counterfamily_reader
#print axioms golden_window

end OPH.SourceLaw
