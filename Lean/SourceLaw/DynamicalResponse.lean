import SourceLaw.MatrixLimits
import SourceLaw.CounterfamilyLocality
import SourceLaw.MatrixSeams

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The two-time reader fixes the reconstructed dynamical response

For a positive definite time coupling `J` and stiffness `V`, the stationary chain whitened by
`J^{-1/2}` has whitened stiffness `A_J = J^{-1/2} V J^{-1/2}`. Its transfer-matrix reconstruction
has Hamiltonian `(pᵀ 𝒦 p + qᵀ 𝒱 q)/2` with `𝒦 𝒱 = Ω_Δτ²` (Heisenberg equations `q̈ = −𝒦 𝒱 q` in
model time, recorded in the companion note). This file proves:

* `reconstructed_acceleration_eq`: `𝒦(A_J) 𝒱(A_J) = Ω_Δτ(A_J)²` for `Δτ > 0`.
* `reconstructed_acceleration_tendsto`: in original coordinates the acceleration matrix
  `J^{-1/2} Ω_Δτ(A_J)² J^{1/2}` tends to `J⁻¹ V` as `Δτ → 0⁺`.
* `counterfamily_acceleration`, `counterfamily_acceleration_whitened`: on the `I + εA` family,
  `J_ε⁻¹ V = M⁻¹ V + ε M⁻¹ V M⁻¹ V`, which is `A + ε A²` in whitened coordinates.
* `counterfamily_dynamical_response`: for positive capacities, positive definite `V` and `ε ≥ 0`,
  the two-time reader of the family's stationary chain returns `J_ε`, and the reconstructed
  acceleration matrix tends to `M⁻¹ V + ε M⁻¹ V M⁻¹ V`.
* `counterfamilyInertia_posDef_of_shift`, `counterfamily_dynamical_response_shift`: the same holds
  on every member with `M + εV` positive definite, negative `ε` included; on the path
  `path3_shift_posDef_iff` gives that domain exactly as `ε > −1/4`.
* `calibrated_response_iff`: some `αI + βA` with real `α`, `β` (the calibrated same-capacity
  response, read as in OPH's control script, which fits `α` and `β` on two modes) equals the
  family response `A + εA²` exactly when `ε = 0` or `A²` lies in the span of `1` and `A`;
  `path3_calibrated_separation`: on the three-site path this happens only at `ε = 0`, because
  the distance-two entry of `αI + βA` is `0` and that of `A + εA²` is `ε`.
* `path3_stiffness_response`, `path3_dynamical_response`: on the three-site path `A = 1 + L_path3`
  with unit capacities (the matrix of the OPH control `kinetic_control.py`), the distance-two
  entry of `γ (A + ε A²)` is `γ ε`; on the family over this path the reader returns `J_ε`, and the
  distance-two entry of `A + ε A²`, which `counterfamily_dynamical_response` identifies as the
  limit of the acceleration matrix, is `ε`, the stiffness-response coefficient of that control at
  common clock factor 1.
-/

namespace OPH.SourceLaw

open Filter Topology Matrix

noncomputable section

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

theorem kineticCoeff_mul_potentialCoeff {dτ lam : ℝ} (hd : 0 < dτ) (hl : 0 < lam) :
    kineticCoeff dτ lam * potentialCoeff dτ lam = latticeFreq dτ lam ^ 2 := by
  have hs := (sinh_latticeFreq_pos hd hl).ne'
  unfold kineticCoeff potentialCoeff
  field_simp

/-- The acceleration matrix of the reconstructed Hamiltonian is `Ω_Δτ²`. -/
theorem reconstructed_acceleration_eq {A : Matrix ι ι ℝ} (hA : A.PosDef) {dτ : ℝ} (hd : 0 < dτ) :
    hA.1.cfc (kineticCoeff dτ) * hA.1.cfc (potentialCoeff dτ)
      = hA.1.cfc (fun x => latticeFreq dτ x ^ 2) := by
  rw [cfc_mul hA.1]
  exact cfc_congr hA.1 fun i => kineticCoeff_mul_potentialCoeff hd (eig_pos hA i)

/-- In original coordinates the reconstructed acceleration matrix tends to `J⁻¹ V`. -/
theorem reconstructed_acceleration_tendsto {J V : Matrix ι ι ℝ} (hJ : J.PosDef)
    (hAJ : (posDefInvSqrt hJ * V * posDefInvSqrt hJ).PosDef) :
    Tendsto (fun dτ => posDefInvSqrt hJ * hAJ.1.cfc (fun x => latticeFreq dτ x ^ 2) * posDefSqrt hJ)
      (𝓝[>] 0) (𝓝 (J⁻¹ * V)) := by
  have hlim : Tendsto (fun dτ => hAJ.1.cfc (fun x => latticeFreq dτ x ^ 2)) (𝓝[>] 0)
      (𝓝 (hAJ.1.cfc id)) := by
    refine cfc_tendsto hAJ.1 (f := fun dτ x => latticeFreq dτ x ^ 2) fun i => ?_
    have hl := hAJ.eigenvalues_pos i
    have := (latticeFreq_tendsto hl).pow 2
    simpa [Real.sq_sqrt hl.le] using this
  rw [cfc_id hAJ.1] at hlim
  have hmul := (tendsto_const_nhds (x := posDefInvSqrt hJ)).mul (hlim.mul
    (tendsto_const_nhds (x := posDefSqrt hJ)))
  have h1 : posDefInvSqrt hJ * posDefSqrt hJ = 1 := posDefInvSqrt_mul_sqrt hJ
  have h2 : posDefInvSqrt hJ * posDefInvSqrt hJ = J⁻¹ := by
    refine (Matrix.inv_eq_left_inv ?_).symm
    calc posDefInvSqrt hJ * posDefInvSqrt hJ * J
        = posDefInvSqrt hJ * posDefInvSqrt hJ * (posDefSqrt hJ * posDefSqrt hJ) := by
          rw [posDefSqrt_mul_self hJ]
      _ = posDefInvSqrt hJ * (posDefInvSqrt hJ * posDefSqrt hJ) * posDefSqrt hJ := by
          simp only [Matrix.mul_assoc]
      _ = 1 := by rw [h1, Matrix.mul_one, h1]
  have hJinv : posDefInvSqrt hJ * (posDefInvSqrt hJ * V * posDefInvSqrt hJ) * posDefSqrt hJ
      = J⁻¹ * V := by
    calc posDefInvSqrt hJ * (posDefInvSqrt hJ * V * posDefInvSqrt hJ) * posDefSqrt hJ
        = (posDefInvSqrt hJ * posDefInvSqrt hJ) * V * (posDefInvSqrt hJ * posDefSqrt hJ) := by
          simp only [Matrix.mul_assoc]
      _ = J⁻¹ * V := by rw [h2, h1, Matrix.mul_one]
  rw [← hJinv]
  simpa only [Matrix.mul_assoc] using hmul

/-- The acceleration matrix of the `I + εA` family in original coordinates. -/
theorem counterfamily_acceleration (M : ι → ℝ) (V : Matrix ι ι ℝ) (ε : ℝ) (hM : ∀ i, M i ≠ 0)
    (hU : IsUnit (diagonal M + ε • V).det) :
    (counterfamilyInertia M V ε)⁻¹ * V
      = diagonal (fun i => (M i)⁻¹) * V
        + ε • (diagonal (fun i => (M i)⁻¹) * V * diagonal (fun i => (M i)⁻¹) * V) := by
  have h := counterfamilyInertia_mul_inverseInertia M V ε hM hU
  rw [Matrix.inv_eq_right_inv h, counterfamilyInverseInertia, Matrix.add_mul, Matrix.smul_mul]

/-- The acceleration matrix of the family in whitened coordinates is `A + ε A²`. -/
theorem counterfamily_acceleration_whitened {M : ι → ℝ} (V : Matrix ι ι ℝ) (ε : ℝ)
    (hM : ∀ i, 0 < M i) (hU : IsUnit (diagonal M + ε • V).det) :
    jSqrt M * ((counterfamilyInertia M V ε)⁻¹ * V) * jInvSqrt M
      = jInvSqrt M * V * jInvSqrt M
        + ε • ((jInvSqrt M * V * jInvSqrt M) * (jInvSqrt M * V * jInvSqrt M)) := by
  rw [counterfamily_acceleration M V ε (fun i => (hM i).ne') hU, ← jInvSqrt_mul_self hM]
  have h1 : jSqrt M * jInvSqrt M = 1 := jSqrt_mul_jInvSqrt hM
  have e1 : jSqrt M * (jInvSqrt M * jInvSqrt M * V) * jInvSqrt M = jInvSqrt M * V * jInvSqrt M := by
    calc jSqrt M * (jInvSqrt M * jInvSqrt M * V) * jInvSqrt M
        = (jSqrt M * jInvSqrt M) * jInvSqrt M * V * jInvSqrt M := by simp only [Matrix.mul_assoc]
      _ = jInvSqrt M * V * jInvSqrt M := by rw [h1, Matrix.one_mul]
  have e2 : jSqrt M * (jInvSqrt M * jInvSqrt M * V * (jInvSqrt M * jInvSqrt M) * V) * jInvSqrt M
      = (jInvSqrt M * V * jInvSqrt M) * (jInvSqrt M * V * jInvSqrt M) := by
    calc jSqrt M * (jInvSqrt M * jInvSqrt M * V * (jInvSqrt M * jInvSqrt M) * V) * jInvSqrt M
        = (jSqrt M * jInvSqrt M) * jInvSqrt M * V * jInvSqrt M * jInvSqrt M * V * jInvSqrt M := by
          simp only [Matrix.mul_assoc]
      _ = (jInvSqrt M * V * jInvSqrt M) * (jInvSqrt M * V * jInvSqrt M) := by
          rw [h1, Matrix.one_mul]; simp only [Matrix.mul_assoc]
  rw [Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.smul_mul, e1, e2]

/-- The reader returns the family's inertia, and the reconstructed acceleration matrix of that
chain tends to `M⁻¹ V + ε M⁻¹ V M⁻¹ V`. -/
theorem counterfamily_dynamical_response {M : ι → ℝ} {V : Matrix ι ι ℝ} {ε : ℝ}
    (hM : ∀ i, 0 < M i) (hV : V.PosDef) (hε : 0 ≤ ε) :
    (∀ {dτ : ℝ}, 0 < dτ → cfReader hM hV hε dτ = counterfamilyInertia M V ε) ∧
    Tendsto (fun dτ => posDefInvSqrt (counterfamilyInertia_posDef hM hV hε)
        * (whitenedBy_posDef (counterfamilyInertia_posDef hM hV hε) hV).1.cfc
            (fun x => latticeFreq dτ x ^ 2)
        * posDefSqrt (counterfamilyInertia_posDef hM hV hε)) (𝓝[>] 0)
      (𝓝 (diagonal (fun i => (M i)⁻¹) * V
        + ε • (diagonal (fun i => (M i)⁻¹) * V * diagonal (fun i => (M i)⁻¹) * V))) := by
  refine ⟨fun hd => counterfamily_history_reader hM hV hε hd, ?_⟩
  have hU : IsUnit (diagonal M + ε • V).det :=
    (Matrix.isUnit_iff_isUnit_det _).1 (capacity_shift_posDef hM hV hε).isUnit
  rw [← counterfamily_acceleration M V ε (fun i => (hM i).ne') hU]
  exact reconstructed_acceleration_tendsto (counterfamilyInertia_posDef hM hV hε)
    (whitenedBy_posDef (counterfamilyInertia_posDef hM hV hε) hV)

/-- The family's inertia `M (M + εV)⁻¹ M` is positive definite whenever `M + εV` is, negative `ε`
included. -/
theorem counterfamilyInertia_posDef_of_shift {M : ι → ℝ} {V : Matrix ι ι ℝ} {ε : ℝ}
    (hM : ∀ i, 0 < M i) (hS : (diagonal M + ε • V).PosDef) :
    (counterfamilyInertia M V ε).PosDef := by
  have hInv := hS.inv
  have hDu : IsUnit (diagonal M) := by
    rw [isUnit_diagonal, Pi.isUnit_iff]
    exact fun i => isUnit_iff_ne_zero.2 (hM i).ne'
  have hDH : (diagonal M)ᴴ = diagonal M := by
    rw [diagonal_conjTranspose]
    simp
  have h := hInv.conjTranspose_mul_mul_same (mulVec_injective_of_isUnit hDu)
  rw [hDH] at h
  exact h

/-- On every member of the family with `M + εV` positive definite, negative `ε` included, the
two-time reader of the member's stationary chain returns its inertia and the reconstructed
acceleration matrix tends to `M⁻¹ V + ε M⁻¹ V M⁻¹ V`. -/
theorem counterfamily_dynamical_response_shift {M : ι → ℝ} {V : Matrix ι ι ℝ} {ε : ℝ}
    (hM : ∀ i, 0 < M i) (hV : V.PosDef) (hS : (diagonal M + ε • V).PosDef) :
    (∀ {dτ : ℝ}, 0 < dτ →
      dτ • (lagW0 (counterfamilyInertia_posDef_of_shift hM hS)
              (whitenedBy_posDef (counterfamilyInertia_posDef_of_shift hM hS) hV) dτ
            * (lagW1 (counterfamilyInertia_posDef_of_shift hM hS)
              (whitenedBy_posDef (counterfamilyInertia_posDef_of_shift hM hS) hV) dτ)⁻¹
            * lagW0 (counterfamilyInertia_posDef_of_shift hM hS)
              (whitenedBy_posDef (counterfamilyInertia_posDef_of_shift hM hS) hV) dτ
            - lagW1 (counterfamilyInertia_posDef_of_shift hM hS)
              (whitenedBy_posDef (counterfamilyInertia_posDef_of_shift hM hS) hV) dτ)⁻¹
        = counterfamilyInertia M V ε) ∧
    Tendsto (fun dτ => posDefInvSqrt (counterfamilyInertia_posDef_of_shift hM hS)
        * (whitenedBy_posDef (counterfamilyInertia_posDef_of_shift hM hS) hV).1.cfc
            (fun x => latticeFreq dτ x ^ 2)
        * posDefSqrt (counterfamilyInertia_posDef_of_shift hM hS)) (𝓝[>] 0)
      (𝓝 (diagonal (fun i => (M i)⁻¹) * V
        + ε • (diagonal (fun i => (M i)⁻¹) * V * diagonal (fun i => (M i)⁻¹) * V))) := by
  refine ⟨fun hd => history_reader_posDef _ _ hd, ?_⟩
  have hU : IsUnit (diagonal M + ε • V).det :=
    (Matrix.isUnit_iff_isUnit_det _).1 hS.isUnit
  rw [← counterfamily_acceleration M V ε (fun i => (hM i).ne') hU]
  exact reconstructed_acceleration_tendsto (counterfamilyInertia_posDef_of_shift hM hS)
    (whitenedBy_posDef (counterfamilyInertia_posDef_of_shift hM hS) hV)

/-- Calibration test: some `αI + βA` with real `α`, `β` (the calibrated same-capacity response, read
as in OPH's control script) equals the family response `A + εA²` exactly when `ε = 0` or `A²` lies
in the span of `1` and `A`. -/
theorem calibrated_response_iff (A : Matrix ι ι ℝ) (ε : ℝ) :
    (∃ α β : ℝ, A + ε • (A * A) = α • (1 : Matrix ι ι ℝ) + β • A) ↔
      (ε = 0 ∨ ∃ a b : ℝ, A * A = a • (1 : Matrix ι ι ℝ) + b • A) := by
  constructor
  · rintro ⟨α, β, h⟩
    by_cases hε : ε = 0
    · exact Or.inl hε
    · right
      refine ⟨ε⁻¹ * α, ε⁻¹ * (β - 1), ?_⟩
      have h2 : ε • (A * A) = α • (1 : Matrix ι ι ℝ) + (β - 1) • A := by
        rw [eq_sub_of_add_eq' h, sub_smul, one_smul]
        abel
      have h3 : A * A = ε⁻¹ • (ε • (A * A)) := by
        rw [smul_smul, inv_mul_cancel₀ hε, one_smul]
      rw [h3, h2, smul_add, smul_smul, smul_smul]
  · rintro (h | ⟨a, b, h⟩)
    · exact ⟨0, 1, by simp [h]⟩
    · refine ⟨ε * a, 1 + ε * b, ?_⟩
      rw [h, smul_add, smul_smul, smul_smul, add_smul, one_smul]
      abel

/-- The three-site path stiffness `1 + L_path3` of the OPH control. -/
def path3A : Matrix (Fin 3) (Fin 3) ℝ := !![2, -1, 0; -1, 3, -1; 0, -1, 2]

theorem path3A_posDef : path3A.PosDef := by
  apply posDef_of_matSeams
  · ext i j; fin_cases i <;> fin_cases j <;> simp [path3A, Matrix.transpose_apply]
  · intro i j hij
    fin_cases i <;> fin_cases j <;> simp_all [path3A]
  · intro i
    fin_cases i <;> simp [rowSum, path3A, Fin.sum_univ_three]
    all_goals norm_num

/-- On the path with unit capacities the kinetic matrix `1 + εA` of the family is positive definite
exactly when `ε > −1/4` (the eigenvalues of `A` are `1, 2, 4`); this covers the control cases
`ε = 0, ±1/16` and excludes the singular `ε = −1/4`. -/
theorem path3_shift_posDef_iff (ε : ℝ) :
    (diagonal (fun _ : Fin 3 => (1 : ℝ)) + ε • path3A).PosDef ↔ -(1 / 4) < ε := by
  constructor
  · intro h
    have hx : (![1, -2, 1] : Fin 3 → ℝ) ≠ 0 := by
      intro hx
      have := congrFun hx 0
      simp at this
    have hv := h.dotProduct_mulVec_pos hx
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_three, path3A, diagonal] at hv
    nlinarith
  · intro hε
    refine Matrix.PosDef.of_dotProduct_mulVec_pos ?_ fun x hx => ?_
    · ext i j
      fin_cases i <;> fin_cases j <;> simp [path3A]
    · have hx' : 0 < x 0 ^ 2 + x 1 ^ 2 + x 2 ^ 2 := by
        by_contra hneg
        apply hx
        funext i
        fin_cases i <;> simp <;> nlinarith [sq_nonneg (x 0), sq_nonneg (x 1), sq_nonneg (x 2)]
      simp [Matrix.mulVec, dotProduct, Fin.sum_univ_three, path3A, diagonal]
      nlinarith [sq_nonneg (x 0 + x 1 + x 2), sq_nonneg (x 0 - x 2), sq_nonneg (x 0 - x 1),
        sq_nonneg (x 1 - x 2), mul_pos (show (0 : ℝ) < 1 / 4 + ε by linarith) hx']

/-- The distance-two entry of `γ (A + ε A²)` on the path is `γ ε` (by computation). -/
theorem path3_stiffness_response (γ ε : ℝ) :
    (γ • (path3A + ε • (path3A * path3A))) 0 2 = γ * ε := by
  simp [path3A]

/-- On the path with unit capacities and `ε ≥ 0`, the reader returns the family's inertia, and the
distance-two entry of `A + ε A²` (the limit of the acceleration matrix by
`counterfamily_dynamical_response`) is `ε`. -/
theorem path3_dynamical_response {ε : ℝ} (hε : 0 ≤ ε) :
    (∀ {dτ : ℝ}, 0 < dτ →
      cfReader (M := fun _ : Fin 3 => (1 : ℝ)) (fun _ => one_pos) path3A_posDef hε dτ
      = counterfamilyInertia (fun _ => 1) path3A ε) ∧
    (diagonal (fun _ : Fin 3 => ((1 : ℝ))⁻¹) * path3A
        + ε • (diagonal (fun _ : Fin 3 => ((1 : ℝ))⁻¹) * path3A
          * diagonal (fun _ : Fin 3 => ((1 : ℝ))⁻¹) * path3A) : Matrix (Fin 3) (Fin 3) ℝ) 0 2
      = ε := by
  refine ⟨(counterfamily_dynamical_response (fun _ => one_pos) path3A_posDef hε).1, ?_⟩
  simp [path3A]

/-- On the three-site path no `αI + βA` equals `A + εA²` for `ε ≠ 0`: the distance-two entry of
`αI + βA` is `0` and that of `A + εA²` is `ε`. -/
theorem path3_calibrated_separation (ε : ℝ) :
    (∃ α β : ℝ, path3A + ε • (path3A * path3A) = α • (1 : Matrix (Fin 3) (Fin 3) ℝ) + β • path3A)
      ↔ ε = 0 := by
  constructor
  · rintro ⟨α, β, h⟩
    have := congrFun (congrFun h 0) 2
    simpa [path3A, Matrix.mul_apply, Fin.sum_univ_three] using this
  · rintro rfl
    exact ⟨0, 1, by simp⟩

#print axioms kineticCoeff_mul_potentialCoeff
#print axioms reconstructed_acceleration_eq
#print axioms reconstructed_acceleration_tendsto
#print axioms counterfamily_acceleration
#print axioms counterfamily_acceleration_whitened
#print axioms counterfamily_dynamical_response
#print axioms counterfamilyInertia_posDef_of_shift
#print axioms counterfamily_dynamical_response_shift
#print axioms path3A_posDef
#print axioms path3_shift_posDef_iff
#print axioms path3_stiffness_response
#print axioms path3_dynamical_response
#print axioms calibrated_response_iff
#print axioms path3_calibrated_separation

end

end OPH.SourceLaw
