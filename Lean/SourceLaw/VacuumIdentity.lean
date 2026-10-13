import SourceLaw.RepairConvergence
import SourceLaw.OpenWindow
import SourceLaw.ContinuumConnection
import SourceLaw.MatrixLimits

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The reference state reached by the declared dynamics is the reconstructed vacuum

Three limits compose on the register of a seam graph with capacities `M` (time coupling `M`, open
ends, `T + 1` slabs). In sample time `s → ∞` the covariance flow of the declared law takes every
initial matrix to `ℏ P_T⁻¹` (`covFlow_tendsto`). As the window grows around an interior slab, the
equal-time block of `ℏ P_T⁻¹` there tends to `ℏ M^{-1/2} C0(Δτ) M^{-1/2}`
(`declared_interior_block_tendsto`). As `Δτ → 0⁺` that tends to `ℏ vacuumCov` (`origCov_tendsto`);
that `vacuumCov` is the position covariance of the Gaussian ground state of
`½ Πᵀ M⁻¹ Π + ½ fᵀ V f` at `ℏ = 1` is the standard identification recorded in the companion note.

* `iterated_reference_limit`: the three limits with their values, for every carrier with
  nonnegative conductances, positive grounding and positive capacities; its positivity
  hypothesis holds on every such carrier (`carrier_whitened_posDef`).
* `first_slab_reference_limit`: at the first slab the block of `ℏ P_T⁻¹` tends instead to
  `ℏ M^{-1/2} C0(Δτ) (1 + R) M^{-1/2}`, `R = r(A)`, which differs from the interior value for
  `ℏ ≠ 0`.
* `vacuumCov_whitened`: in capacity-whitened coordinates `vacuumCov` is `½ A^{-1/2}`.
* `golden_reference_limit`: at every level `q = n + 2` of the golden family the three limits hold
  and the final value, in whitened coordinates, is `ℏ · ½ ω_q⁻¹` with `ω_q = A_q^{1/2}`, the
  position covariance of the Fock vacuum of `dΓ(ω_q)` at `ℏ = 1` (standard identification).
* Two-site carrier (`M = 1`, one seam of conductance `κ = 9/32`, ground `1`): `twoSite_stiffness`
  gives the stiffness `[[41/32, −9/32], [−9/32, 41/32]]`, `twoSite_vacuum_cov` gives
  `vacuumCov = [[9/20, 1/20], [1/20, 9/20]]`, which by the standard identification is the position
  covariance of the coupled vacuum of `H = (P0² + P1² + Q0² + Q1²)/2 + κ (Q0 − Q1)²/2`, and
  `twoSite_reference_limit` states the three limits with that final value.
  `twoSite_calibrated_equiv`: on this carrier every family response has the form `αI + βV` with
  real `α`, `β`, namely `γ(V + εV²) = −(25γε/16) 1 + γ(1 + 41ε/16) V`; for `γ > 0` the
  coefficient of `V` is positive exactly when `ε > −16/41`.
-/

namespace OPH.SourceLaw

open Filter Topology Matrix

noncomputable section

section General

variable {ι E : Type*} [Fintype ι] [DecidableEq ι] [Fintype E]

omit [Fintype ι] [DecidableEq ι] in
theorem slabBlock_smul {T : ℕ} (a : ℝ) (X : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ)
    (t : Fin (T + 1)) : slabBlock (a • X) t = a • slabBlock X t := by
  ext i k
  simp [slabBlock]

/-- The three limits of the declared dynamics on the register of a seam graph, with their values:
sample time to `ℏ P_T⁻¹`, interior slabs of a growing window to `ℏ M^{-1/2} C0(Δτ) M^{-1/2}`, and
`Δτ → 0⁺` to `ℏ vacuumCov`. -/
def ReferenceLimits (src tgt : E → ι) (c : E → ℝ) (g M : ι → ℝ) (hc : ∀ e, 0 ≤ c e)
    (hg : ∀ i, 0 < g i) (hM : ∀ i, 0 < M i)
    (hA : (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M).PosDef) (hbar : ℝ) : Prop :=
  (∀ (dτ : ℝ) (hd : 0 < dτ) (T : ℕ) (S0 : Matrix (ι × Fin (T + 1)) (ι × Fin (T + 1)) ℝ),
      Tendsto (fun s => covFlow (M := regCap M T) (fun p => hM p.1)
          (register_posDef src tgt c g M dτ T hc hg (fun i => (hM i).le) hd) hbar S0 s) atTop
        (𝓝 (stationaryCov hbar (stiffnessMatrix (regSrc src T) (regTgt tgt T)
          (regCond c M dτ T) (regGround g dτ T)))))
    ∧ (∀ (dτ : ℝ), 0 < dτ → ∀ (T t : ℕ → ℕ) (ht : ∀ k, t k ≤ T k), Tendsto t atTop atTop →
      Tendsto (fun k => T k - t k) atTop atTop →
      Tendsto (fun k => slabBlock (stationaryCov hbar (stiffnessMatrix (regSrc src (T k))
          (regTgt tgt (T k)) (regCond c M dτ (T k)) (regGround g dτ (T k))))
          ⟨t k, Nat.lt_succ_of_le (ht k)⟩) atTop
        (𝓝 (hbar • (jInvSqrt M * equalTimeCov hA dτ * jInvSqrt M))))
    ∧ Tendsto (fun dτ => hbar • (jInvSqrt M * equalTimeCov hA dτ * jInvSqrt M)) (𝓝[>] 0)
        (𝓝 (hbar • vacuumCov M hA))

/-- The three limits hold on every carrier with nonnegative conductances, positive grounding and
positive capacities. -/
theorem iterated_reference_limit (src tgt : E → ι) (c : E → ℝ) (g M : ι → ℝ)
    (hc : ∀ e, 0 ≤ c e) (hg : ∀ i, 0 < g i) (hM : ∀ i, 0 < M i)
    (hA : (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M).PosDef) (hbar : ℝ) :
    ReferenceLimits src tgt c g M hc hg hM hA hbar := by
  refine ⟨fun dτ hd T S0 => covFlow_tendsto _ _ hbar S0, fun dτ hd T t ht h1 h2 => ?_,
    (origCov_tendsto (M := M) hA).const_smul hbar⟩
  have h := (declared_interior_block_tendsto src tgt c g M hA hM hd ht h1 h2).const_smul hbar
  refine h.congr fun k => ?_
  rw [stationaryCov, slabBlock_smul]

/-- On every carrier with nonnegative conductances, positive grounding and positive capacities the
whitened stiffness is positive definite, so the positivity hypothesis of `iterated_reference_limit`
always holds. -/
theorem carrier_whitened_posDef (src tgt : E → ι) (c : E → ℝ) (g M : ι → ℝ) (hc : ∀ e, 0 ≤ c e)
    (hg : ∀ i, 0 < g i) (hM : ∀ i, 0 < M i) :
    (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M).PosDef :=
  jWhitened_posDef hM (stiffnessMatrix_posDef src tgt c g hc hg)

/-- At the first slab the block of `ℏ P_T⁻¹` tends to `ℏ M^{-1/2} C0(Δτ) (1 + R) M^{-1/2}`, with
`R = r(A)` the transfer matrix; for `ℏ ≠ 0` on a nonempty carrier this differs from the interior
value `ℏ M^{-1/2} C0(Δτ) M^{-1/2}`. -/
theorem first_slab_reference_limit [Nonempty ι] (src tgt : E → ι) (c : E → ℝ) (g M : ι → ℝ)
    (hM : ∀ i, 0 < M i) (hA : (jInvSqrt M * stiffnessMatrix src tgt c g * jInvSqrt M).PosDef)
    (hbar : ℝ) {dτ : ℝ} (hd : 0 < dτ) :
    Tendsto (fun T : ℕ => slabBlock (stationaryCov hbar (stiffnessMatrix (regSrc src T)
        (regTgt tgt T) (regCond c M dτ T) (regGround g dτ T))) 0) atTop
        (𝓝 (hbar • (jInvSqrt M * (equalTimeCov hA dτ * (1 + transferMat hA dτ)) * jInvSqrt M)))
      ∧ (hbar ≠ 0 →
        hbar • (jInvSqrt M * (equalTimeCov hA dτ * (1 + transferMat hA dτ)) * jInvSqrt M)
          ≠ hbar • (jInvSqrt M * equalTimeCov hA dτ * jInvSqrt M)) := by
  refine ⟨?_, fun hb h => ?_⟩
  · have h := (declared_end_block_tendsto src tgt c g M hA hM hd).const_smul hbar
    refine h.congr fun T => ?_
    rw [stationaryCov, slabBlock_smul]
  · have h1 := smul_right_injective (Matrix ι ι ℝ) hb h
    have h2 : jSqrt M * (jInvSqrt M * (equalTimeCov hA dτ * (1 + transferMat hA dτ)) * jInvSqrt M)
          * jSqrt M = jSqrt M * (jInvSqrt M * equalTimeCov hA dτ * jInvSqrt M) * jSqrt M := by
      rw [h1]
    have e : ∀ X : Matrix ι ι ℝ, jSqrt M * (jInvSqrt M * X * jInvSqrt M) * jSqrt M = X := by
      intro X
      calc jSqrt M * (jInvSqrt M * X * jInvSqrt M) * jSqrt M
          = (jSqrt M * jInvSqrt M) * X * (jInvSqrt M * jSqrt M) := by
            simp only [Matrix.mul_assoc]
        _ = X := by
            rw [jSqrt_mul_jInvSqrt hM, jInvSqrt_mul_jSqrt hM, Matrix.one_mul, Matrix.mul_one]
    rw [e, e] at h2
    exact end_block_ne hA hd h2

/-- In capacity-whitened coordinates the vacuum covariance is `½ A^{-1/2}`. -/
theorem vacuumCov_whitened {M : ι → ℝ} (hM : ∀ i, 0 < M i) {A : Matrix ι ι ℝ} (hA : A.PosDef) :
    jSqrt M * vacuumCov M hA * jSqrt M = (1 / 2 : ℝ) • posDefInvSqrt hA := by
  unfold vacuumCov
  rw [Matrix.mul_smul, Matrix.smul_mul]
  congr 1
  calc jSqrt M * (jInvSqrt M * posDefInvSqrt hA * jInvSqrt M) * jSqrt M
      = (jSqrt M * jInvSqrt M) * posDefInvSqrt hA * (jInvSqrt M * jSqrt M) := by
        simp only [Matrix.mul_assoc]
    _ = posDefInvSqrt hA := by
        rw [jSqrt_mul_jInvSqrt hM, jInvSqrt_mul_jSqrt hM, Matrix.one_mul, Matrix.mul_one]

/-- Functions of equal symmetric matrices agree. -/
theorem cfc_eq_of_eq {A B : Matrix ι ι ℝ} (hA : A.IsHermitian) (hB : B.IsHermitian) (h : A = B)
    (f : ℝ → ℝ) : hA.cfc f = hB.cfc f := by
  subst h
  rfl

end General

section Golden

/-- The whitened register stiffness of a golden level is positive definite. -/
theorem golden_register_whitened (n : ℕ) {m2 : ℝ} (hm : 0 < m2) :
    (jInvSqrt (capacity3 (goldenGrid (n + 2)))
      * stiffnessMatrix matSrc matTgt (edge3 (goldenGrid (n + 2)) m2)
        (ground3 (goldenGrid (n + 2)) m2)
      * jInvSqrt (capacity3 (goldenGrid (n + 2)))).PosDef := by
  rw [stiffnessMatrix_packet]
  exact jWhitened_posDef (capacity3_pos (goldenGrid_strictMono (n + 2)))
    (packetV_posDef (goldenGrid_strictMono (n + 2)) hm)

theorem golden_register_whitened_eq (n : ℕ) (m2 : ℝ) :
    jInvSqrt (capacity3 (goldenGrid (n + 2)))
      * stiffnessMatrix matSrc matTgt (edge3 (goldenGrid (n + 2)) m2)
        (ground3 (goldenGrid (n + 2)) m2)
      * jInvSqrt (capacity3 (goldenGrid (n + 2)))
      = whitenedPacket (goldenGrid (n + 2)) m2 := by
  rw [stiffnessMatrix_packet]
  rfl

/-- At every golden level the three limits hold, and the final value in whitened coordinates is
`ℏ · ½ ω_q⁻¹`, `ω_q = A_q^{1/2}`. -/
theorem golden_reference_limit (n : ℕ) {m2 : ℝ} (hm : 0 < m2) (hbar : ℝ) :
    ReferenceLimits matSrc matTgt (edge3 (goldenGrid (n + 2)) m2) (ground3 (goldenGrid (n + 2)) m2)
        (capacity3 (goldenGrid (n + 2))) (edge3_nonneg (goldenGrid_strictMono (n + 2)) m2)
        (ground3_pos (goldenGrid_strictMono (n + 2)) hm)
        (capacity3_pos (goldenGrid_strictMono (n + 2))) (golden_register_whitened n hm) hbar
      ∧ jSqrt (capacity3 (goldenGrid (n + 2)))
          * (hbar • vacuumCov (capacity3 (goldenGrid (n + 2))) (golden_register_whitened n hm))
          * jSqrt (capacity3 (goldenGrid (n + 2)))
        = hbar • ((1 / 2 : ℝ) • posDefInvSqrt (golden_whitened_posDef n hm)) := by
  refine ⟨iterated_reference_limit _ _ _ _ _ _ _ _ _ hbar, ?_⟩
  rw [Matrix.mul_smul, Matrix.smul_mul,
    vacuumCov_whitened (capacity3_pos (goldenGrid_strictMono (n + 2)))
      (golden_register_whitened n hm)]
  unfold posDefInvSqrt
  rw [cfc_eq_of_eq (golden_register_whitened n hm).1 (golden_whitened_posDef n hm).1
    (golden_register_whitened_eq n m2)]

end Golden

section TwoSite

/-- First end of the single seam. -/
def twoSiteSrc : Fin 1 → Fin 2 := fun _ => 0

/-- Second end of the single seam. -/
def twoSiteTgt : Fin 1 → Fin 2 := fun _ => 1

/-- Conductance `κ = 9/32` of the single seam. -/
def twoSiteCond : Fin 1 → ℝ := fun _ => 9 / 32

/-- Unit grounding at both nodes. -/
def twoSiteGround : Fin 2 → ℝ := fun _ => 1

/-- The stiffness `1 + κ L₂` with `κ = 9/32`. -/
def twoSiteV : Matrix (Fin 2) (Fin 2) ℝ := !![41 / 32, -9 / 32; -9 / 32, 41 / 32]

theorem twoSite_stiffness :
    stiffnessMatrix twoSiteSrc twoSiteTgt twoSiteCond twoSiteGround = twoSiteV := by
  ext k l
  fin_cases k <;> fin_cases l <;>
    simp [stiffnessMatrix, stiffness, twoSiteSrc, twoSiteTgt, twoSiteCond, twoSiteGround, twoSiteV,
      Pi.single_apply] <;> norm_num

theorem twoSite_whitened :
    jInvSqrt (fun _ : Fin 2 => (1 : ℝ))
      * stiffnessMatrix twoSiteSrc twoSiteTgt twoSiteCond twoSiteGround
      * jInvSqrt (fun _ : Fin 2 => (1 : ℝ)) = twoSiteV := by
  have h1 : jInvSqrt (fun _ : Fin 2 => (1 : ℝ)) = 1 := by
    rw [jInvSqrt, ← diagonal_one]
    congr 1
    funext i
    simp
  rw [h1, Matrix.one_mul, Matrix.mul_one, twoSite_stiffness]

theorem twoSiteV_posDef : twoSiteV.PosDef := by
  apply posDef_of_matSeams
  · ext i j; fin_cases i <;> fin_cases j <;> simp [twoSiteV, Matrix.transpose_apply]
  · intro i j hij
    fin_cases i <;> fin_cases j <;> simp_all [twoSiteV] <;> norm_num
  · intro i
    fin_cases i <;> simp [rowSum, twoSiteV, Fin.sum_univ_two] <;> norm_num

theorem twoSite_register_whitened :
    (jInvSqrt (fun _ : Fin 2 => (1 : ℝ))
      * stiffnessMatrix twoSiteSrc twoSiteTgt twoSiteCond twoSiteGround
      * jInvSqrt (fun _ : Fin 2 => (1 : ℝ))).PosDef := by
  rw [twoSite_whitened]
  exact twoSiteV_posDef

/-- Every eigenvalue of the two-site stiffness is `1` or `25/16`. -/
theorem twoSiteV_eig (i : Fin 2) :
    (twoSiteV_posDef.1.eigenvalues i - 1) * (twoSiteV_posDef.1.eigenvalues i - 25 / 16) = 0 := by
  set hA := twoSiteV_posDef.1
  have e1 : hA.cfc (fun x => (x - 1) * (x - 25 / 16))
      = (hA.cfc id - hA.cfc (fun _ => 1)) * (hA.cfc id - (25 / 16 : ℝ) • hA.cfc (fun _ => 1)) := by
    rw [cfc_sub hA, cfc_smul hA, cfc_sub hA, cfc_mul hA]
    exact cfc_congr hA fun i => by simp only [id]; ring
  have hpoly : hA.cfc (fun x => (x - 1) * (x - 25 / 16)) = 0 := by
    rw [e1, cfc_one hA, cfc_id hA]
    ext i j; fin_cases i <;> fin_cases j <;>
      simp [twoSiteV, Matrix.mul_apply, Fin.sum_univ_two, Matrix.one_apply] <;> norm_num
  set U := (hA.eigenvectorUnitary : Matrix (Fin 2) (Fin 2) ℝ)
  have hU : star U * U = 1 := eigU_star_mul hA
  set d := (fun x : ℝ => (x - 1) * (x - 25 / 16)) ∘ hA.eigenvalues
  have hconj : star U * (U * diagonal d * star U) * U = diagonal d := by
    rw [Matrix.mul_assoc, Matrix.mul_assoc, Matrix.mul_assoc, hU, Matrix.mul_one,
      ← Matrix.mul_assoc, hU, Matrix.one_mul]
  have hd : diagonal d = 0 := by
    rw [← hconj, ← cfc_expand hA, hpoly, Matrix.mul_zero, Matrix.zero_mul]
  have := congrFun (congrFun hd i) i
  simpa [d, diagonal_apply_eq] using this

theorem twoSiteV_half_invSqrt :
    twoSiteV_posDef.1.cfc (fun x => 1 / (2 * Real.sqrt x))
      = !![9 / 20, 1 / 20; 1 / 20, 9 / 20] := by
  set hA := twoSiteV_posDef.1
  have hcongr : hA.cfc (fun x => 1 / (2 * Real.sqrt x))
      = hA.cfc (fun x => (1 / 2 : ℝ) + (-(8 / 45 : ℝ)) * (x - 1)) := by
    apply cfc_congr hA
    intro i
    rcases mul_eq_zero.1 (twoSiteV_eig i) with h | h
    · have h1 : hA.eigenvalues i = 1 := by linarith
      rw [h1]; simp
    · have h1 : hA.eigenvalues i = 25 / 16 := by linarith
      rw [h1]
      have : Real.sqrt (25 / 16 : ℝ) = 5 / 4 := by
        rw [show (25 / 16 : ℝ) = (5 / 4) ^ 2 by norm_num]
        exact Real.sqrt_sq (by norm_num)
      rw [this]; norm_num
  have hlin : hA.cfc (fun x => (1 / 2 : ℝ) + (-(8 / 45 : ℝ)) * (x - 1))
      = (1 / 2 : ℝ) • hA.cfc (fun _ => 1)
        + (-(8 / 45 : ℝ)) • (hA.cfc id - hA.cfc (fun _ => 1)) := by
    rw [cfc_sub hA, cfc_smul hA, cfc_smul hA, cfc_add hA]
    exact cfc_congr hA fun i => by simp only [id]; ring
  rw [hcongr, hlin, cfc_one hA, cfc_id hA]
  ext i j; fin_cases i <;> fin_cases j <;> simp [twoSiteV] <;> norm_num

/-- The vacuum covariance of the two-site carrier is `[[9/20, 1/20], [1/20, 9/20]]`. -/
theorem twoSite_vacuum_cov :
    vacuumCov (fun _ : Fin 2 => (1 : ℝ)) twoSite_register_whitened
      = !![9 / 20, 1 / 20; 1 / 20, 9 / 20] := by
  have h1 : jInvSqrt (fun _ : Fin 2 => (1 : ℝ)) = 1 := by
    rw [jInvSqrt, ← diagonal_one]
    congr 1
    funext i
    simp
  unfold vacuumCov
  simp only [h1, Matrix.one_mul, Matrix.mul_one]
  rw [posDefInvSqrt, cfc_eq_of_eq _ twoSiteV_posDef.1 twoSite_whitened, cfc_smul,
    ← twoSiteV_half_invSqrt]
  exact cfc_congr _ fun i => by field_simp

/-- On the two-site carrier the three limits hold and the final value is
`ℏ [[9/20, 1/20], [1/20, 9/20]]`. -/
theorem twoSite_reference_limit (hbar : ℝ) :
    ReferenceLimits twoSiteSrc twoSiteTgt twoSiteCond twoSiteGround (fun _ => 1)
        (fun _ => by norm_num [twoSiteCond]) (fun _ => one_pos) (fun _ => one_pos)
        twoSite_register_whitened hbar
      ∧ hbar • vacuumCov (fun _ : Fin 2 => (1 : ℝ)) twoSite_register_whitened
        = hbar • !![9 / 20, 1 / 20; 1 / 20, 9 / 20] := by
  refine ⟨iterated_reference_limit _ _ _ _ _ _ _ _ _ hbar, ?_⟩
  rw [twoSite_vacuum_cov]

/-- On the two-site carrier every family response has the calibrated same-capacity form `αI + βV`
with real `α`, `β`: `γ(V + εV²) = −(25γε/16) 1 + γ(1 + 41ε/16) V`, because
`V² = (41/16) V − (25/16) 1`. For `γ > 0` the coefficient of `V` is positive exactly when
`ε > −16/41`. -/
theorem twoSite_calibrated_equiv (ε γ : ℝ) :
    γ • (twoSiteV + ε • (twoSiteV * twoSiteV))
      = (-(25 * γ * ε / 16)) • (1 : Matrix (Fin 2) (Fin 2) ℝ)
        + (γ * (1 + 41 * ε / 16)) • twoSiteV := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [twoSiteV] <;> ring

end TwoSite

#print axioms slabBlock_smul
#print axioms iterated_reference_limit
#print axioms carrier_whitened_posDef
#print axioms first_slab_reference_limit
#print axioms vacuumCov_whitened
#print axioms cfc_eq_of_eq
#print axioms golden_register_whitened
#print axioms golden_register_whitened_eq
#print axioms golden_reference_limit
#print axioms twoSite_stiffness
#print axioms twoSite_whitened
#print axioms twoSiteV_posDef
#print axioms twoSite_register_whitened
#print axioms twoSiteV_eig
#print axioms twoSiteV_half_invSqrt
#print axioms twoSite_vacuum_cov
#print axioms twoSite_reference_limit
#print axioms twoSite_calibrated_equiv

end

end OPH.SourceLaw
