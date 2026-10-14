import SourceLaw.TransferMatrixLift
import SourceLaw.CounterfamilyLocality

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The two-time history reader for a positive definite time coupling

The stationary chain whose precision has time coupling `W` and potential `V` is whitened by any
symmetric `S` with `S⁻¹ S⁻¹ = W`: in the whitened coordinates the coupling is the identity and
the potential is `A = S V S`, so the chain covariances are `C(0) = S c0(A) S` and
`C(1) = S r(A) c0(A) S`. The coupling `W` need not be node-local. The header and the docstrings
below take `Δτ > 0` and a positive definite whitened potential `A` throughout.

Main results:

* `history_reader_whitened`: for any `S` with `S⁻¹ S⁻¹ = W` (symmetry is not needed here) the
  two-time reader `Δτ (C(0) C(1)⁻¹ C(0) − C(1))⁻¹` returns `W`.
* `contact_identity_whitened`: for a symmetric such `S` the slab row of the precision
  `(2/Δτ) W C(0) − (1/Δτ) W (C(1) + C(1)ᵀ) + Δτ V C(0) = 1` holds.
* `history_reader_posDef`: with the spectral square root, the reader returns every positive
  definite coupling `W`.
* `counterfamily_history_reader`: for the kinetic counterfamily with positive capacities `M`,
  positive definite `V` and `ε ≥ 0`, whose coupling is the inertia `M (M + εV)⁻¹ M`, the reader
  `cfReader` of the counterfamily chain returns that inertia.
* `counterfamily_reader_isDiag_iff`: if moreover `V` couples two distinct nodes, the reader is
  node-local exactly when `ε = 0`, so these two-time correlations separate `ε = 0` from every
  `ε > 0`.

The two correlations `lagW0`, `lagW1` use unit-`ℏ` normalization. For raw correlations
`C̃j = ℏ Cj` with supplied `ℏ > 0`, divide by `ℏ` before using the reader; equivalently,
the raw-correlation reader must be multiplied by `ℏ` to return the coupling.

These correlations are defined as the whitened stationary chain's covariances
conjugated by `W^{-1/2}`; this file does not prove that they are the covariances of a history
register with coupling `W`. For a diagonal coupling `W = diag M`, `declared_register_window`
shows that the precision of the declared register with time coupling `M`, plus an end-slab term,
is the inverse of the stationary chain's window covariance conjugated by `M^{-1/2}`; identifying
those window blocks with `lagW0`, `lagW1` is not proved here. The names `history_reader_*` and
`counterfamily_history_reader` refer to these stationary-chain correlations.
-/

namespace OPH.SourceLaw

open scoped BigOperators
open Matrix

noncomputable section

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

section Whitening

variable {S Si W V : Matrix ι ι ℝ}

/-- For any whitening `S` with `S⁻¹ S⁻¹ = W`, a positive definite whitened potential `S V S` and
`Δτ > 0`, the two-time reader returns `W`. -/
theorem history_reader_whitened (h1 : Si * S = 1) (h2 : S * Si = 1) (hW : Si * Si = W)
    (hA : (S * V * S).PosDef) {dτ : ℝ} (hd : 0 < dτ) :
    dτ • ((S * equalTimeCov hA dτ * S)
        * (S * (transferMat hA dτ * equalTimeCov hA dτ) * S)⁻¹ * (S * equalTimeCov hA dτ * S)
        - S * (transferMat hA dτ * equalTimeCov hA dτ) * S)⁻¹ = W := by
  set C := equalTimeCov hA dτ
  set R := transferMat hA dτ
  have hRC : IsUnit (R * C).det := by
    rw [det_mul]
    exact ((Matrix.isUnit_iff_isUnit_det _).1 (transferMat_posDef hA hd).isUnit).mul
      ((Matrix.isUnit_iff_isUnit_det _).1 (equalTimeCov_posDef hA hd).isUnit)
  have hinv1 : (S * (R * C) * S)⁻¹ = Si * (R * C)⁻¹ * Si := conj_inv h2 hRC
  have hcore := reader_core hA hd
  have hmid : S * C * S * (Si * (R * C)⁻¹ * Si) * (S * C * S) - S * (R * C) * S
      = dτ • (S * S) := by
    calc S * C * S * (Si * (R * C)⁻¹ * Si) * (S * C * S) - S * (R * C) * S
        = S * (C * (S * Si) * (R * C)⁻¹ * (Si * S) * C - R * C) * S := by
          simp only [Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_assoc]
      _ = S * (dτ • (1 : Matrix ι ι ℝ)) * S := by
          rw [h2, h1, Matrix.mul_one, Matrix.mul_one, hcore]
      _ = dτ • (S * S) := by
          rw [Matrix.mul_smul, Matrix.smul_mul, Matrix.mul_one]
  rw [hinv1, hmid, ← hW]
  have hinv2 : (dτ • (S * S))⁻¹ = dτ⁻¹ • (Si * Si) := by
    apply Matrix.inv_eq_right_inv
    calc dτ • (S * S) * (dτ⁻¹ • (Si * Si)) = (dτ * dτ⁻¹) • (S * (S * Si) * Si) := by
          rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul]
          simp only [Matrix.mul_assoc]
      _ = 1 := by rw [mul_inv_cancel₀ hd.ne', one_smul, h2, Matrix.mul_one, h2]
  rw [hinv2, smul_smul, mul_inv_cancel₀ hd.ne', one_smul]

/-- The slab row of the precision with time coupling `W` and potential `V`, for a symmetric `S`
with `S⁻¹ S⁻¹ = W`: `(2/Δτ) W C(0) − (1/Δτ) W (C(1) + C(1)ᵀ) + Δτ V C(0) = 1`. -/
theorem contact_identity_whitened (h1 : Si * S = 1) (h2 : S * Si = 1) (hW : Si * Si = W)
    (hSt : Sᵀ = S) (hA : (S * V * S).PosDef) {dτ : ℝ} (hd : 0 < dτ) :
    (2 / dτ) • (W * (S * equalTimeCov hA dτ * S))
      - (1 / dτ) • (W * (S * (transferMat hA dτ * equalTimeCov hA dτ) * S
          + (S * (transferMat hA dτ * equalTimeCov hA dτ) * S)ᵀ))
      + dτ • (V * (S * equalTimeCov hA dτ * S)) = 1 := by
  set C := equalTimeCov hA dτ
  set R := transferMat hA dτ
  have hl : ∀ X : Matrix ι ι ℝ, Si * (S * X) = X := fun X => by
    rw [← Matrix.mul_assoc, h1, Matrix.one_mul]
  have hJS : W = Si * Si := hW.symm
  have hV : V = Si * (S * V * S) * Si := by
    calc V = (Si * S) * V * (S * Si) := by rw [h1, h2, Matrix.one_mul, Matrix.mul_one]
      _ = Si * (S * V * S) * Si := by simp only [Matrix.mul_assoc]
  have hlag1t : (S * (R * C) * S)ᵀ = S * (R * C) * S := by
    rw [transpose_mul, transpose_mul, transpose_mul, hSt, transferMat_transpose hA,
      equalTimeCov_transpose hA, ← transferMat_comm hA]
    simp only [Matrix.mul_assoc]
    rfl
  have e1 : W * (S * C * S) = Si * C * S := by
    rw [hJS]
    simp only [Matrix.mul_assoc, hl]
  have e2 : W * (S * (R * C) * S) = Si * (R * C) * S := by
    rw [hJS]
    simp only [Matrix.mul_assoc, hl]
  have e3 : V * (S * C * S) = Si * ((S * V * S) * C) * S := by
    conv_lhs => rw [hV]
    simp only [Matrix.mul_assoc, hl]
  have key := contact_identity_matrix hA hd
  have hconj : Si * (((2 / dτ) • (1 : Matrix ι ι ℝ) + dτ • (S * V * S)) * C
      - (2 / dτ) • (R * C)) * S = 1 := by
    rw [key, Matrix.mul_one, h1]
  rw [hlag1t, Matrix.mul_add, e1, e2, e3, ← hconj]
  simp only [Matrix.add_mul, Matrix.mul_add, Matrix.sub_mul, Matrix.mul_sub, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.one_mul, Matrix.mul_assoc]
  module

end Whitening

section PosDefCoupling

variable {W : Matrix ι ι ℝ} (hW : W.PosDef)

/-- Spectral `W^{-1/2}` of a positive definite coupling. -/
def posDefInvSqrt : Matrix ι ι ℝ := hW.1.cfc fun x => (Real.sqrt x)⁻¹

/-- Spectral `W^{1/2}`. -/
def posDefSqrt : Matrix ι ι ℝ := hW.1.cfc Real.sqrt

theorem posDefSqrt_mul_invSqrt : posDefSqrt hW * posDefInvSqrt hW = 1 := by
  rw [posDefSqrt, posDefInvSqrt, cfc_mul, ← cfc_one hW.1]
  exact cfc_congr hW.1 fun i => mul_inv_cancel₀ (Real.sqrt_pos.2 (hW.eigenvalues_pos i)).ne'

theorem posDefInvSqrt_mul_sqrt : posDefInvSqrt hW * posDefSqrt hW = 1 := by
  rw [posDefSqrt, posDefInvSqrt, cfc_mul, ← cfc_one hW.1]
  exact cfc_congr hW.1 fun i => inv_mul_cancel₀ (Real.sqrt_pos.2 (hW.eigenvalues_pos i)).ne'

theorem posDefSqrt_mul_self : posDefSqrt hW * posDefSqrt hW = W := by
  rw [posDefSqrt, cfc_mul]
  conv_rhs => rw [← cfc_id hW.1]
  exact cfc_congr hW.1 fun i => Real.mul_self_sqrt (hW.eigenvalues_pos i).le

theorem posDefInvSqrt_transpose : (posDefInvSqrt hW)ᵀ = posDefInvSqrt hW :=
  cfc_transpose hW.1 _

/-- Whitening by a positive definite coupling keeps a positive definite potential positive
definite. -/
theorem whitenedBy_posDef {V : Matrix ι ι ℝ} (hV : V.PosDef) :
    (posDefInvSqrt hW * V * posDefInvSqrt hW).PosDef := by
  have hdet : (posDefInvSqrt hW).det * (posDefSqrt hW).det = 1 := by
    rw [← det_mul, posDefInvSqrt_mul_sqrt, det_one]
  have hSu : IsUnit (posDefInvSqrt hW) :=
    (Matrix.isUnit_iff_isUnit_det _).2 (isUnit_iff_ne_zero.2 (left_ne_zero_of_mul_eq_one hdet))
  have hSH : (posDefInvSqrt hW)ᴴ = posDefInvSqrt hW := by
    rw [conjTranspose_eq_transpose_of_trivial, posDefInvSqrt_transpose]
  have h := hV.conjTranspose_mul_mul_same (mulVec_injective_of_isUnit hSu)
  rw [hSH] at h
  exact h

variable {V : Matrix ι ι ℝ} (hA : (posDefInvSqrt hW * V * posDefInvSqrt hW).PosDef) (dτ : ℝ)

/-- Equal-time covariance attached to the time coupling `W`: the whitened stationary chain's
equal-time covariance conjugated by `W^{-1/2}` (a definition; see the module note). -/
def lagW0 : Matrix ι ι ℝ := posDefInvSqrt hW * equalTimeCov hA dτ * posDefInvSqrt hW

/-- One-slab covariance `E[f_{t+1} f_tᵀ]` attached to `W`, defined in the same way. -/
def lagW1 : Matrix ι ι ℝ :=
  posDefInvSqrt hW * (transferMat hA dτ * equalTimeCov hA dτ) * posDefInvSqrt hW

variable {dτ}

/-- The two-time reader returns every positive definite time coupling `W`. Here the reader acts
on `lagW0`, `lagW1`, the lag covariances of the whitened stationary chain conjugated by
`W^{-1/2}` on both sides, not on a finite register (`history_reader_general` reads `lag0`,
`lag1`; `history_reader_whitened` reads the correlations conjugated by `S`). -/
theorem history_reader_posDef (hd : 0 < dτ) :
    dτ • (lagW0 hW hA dτ * (lagW1 hW hA dτ)⁻¹ * lagW0 hW hA dτ - lagW1 hW hA dτ)⁻¹ = W :=
  history_reader_whitened (posDefSqrt_mul_invSqrt hW) (posDefInvSqrt_mul_sqrt hW)
    (posDefSqrt_mul_self hW) hA hd

/-- The slab row of the precision with a positive definite time coupling. -/
theorem contact_identity_posDef (hd : 0 < dτ) :
    (2 / dτ) • (W * lagW0 hW hA dτ) - (1 / dτ) • (W * (lagW1 hW hA dτ + (lagW1 hW hA dτ)ᵀ))
      + dτ • (V * lagW0 hW hA dτ) = 1 :=
  contact_identity_whitened (posDefSqrt_mul_invSqrt hW) (posDefInvSqrt_mul_sqrt hW)
    (posDefSqrt_mul_self hW) (posDefInvSqrt_transpose hW) hA hd

end PosDefCoupling

section Counterfamily

variable {M : ι → ℝ} {V : Matrix ι ι ℝ} {ε : ℝ}

omit [Fintype ι] in
theorem capacity_shift_posDef (hM : ∀ i, 0 < M i) (hV : V.PosDef) (hε : 0 ≤ ε) :
    (diagonal M + ε • V).PosDef :=
  (Matrix.PosDef.diagonal hM).add_posSemidef (hV.posSemidef.smul hε)

/-- The counterfamily inertia `M (M + εV)⁻¹ M` is positive definite for `ε ≥ 0`. -/
theorem counterfamilyInertia_posDef (hM : ∀ i, 0 < M i) (hV : V.PosDef) (hε : 0 ≤ ε) :
    (counterfamilyInertia M V ε).PosDef := by
  have hInv := (capacity_shift_posDef hM hV hε).inv
  have hDu : IsUnit (diagonal M) := by
    rw [isUnit_diagonal, Pi.isUnit_iff]
    exact fun i => isUnit_iff_ne_zero.2 (hM i).ne'
  have hDH : (diagonal M)ᴴ = diagonal M := by
    rw [diagonal_conjTranspose]
    simp
  have h := hInv.conjTranspose_mul_mul_same (mulVec_injective_of_isUnit hDu)
  rw [hDH] at h
  exact h

/-- The two-time history reader `Δτ (C(0) C(1)⁻¹ C(0) − C(1))⁻¹` of the counterfamily chain. -/
def cfReader (hM : ∀ i, 0 < M i) (hV : V.PosDef) (hε : 0 ≤ ε) (dτ : ℝ) : Matrix ι ι ℝ :=
  dτ • (lagW0 (counterfamilyInertia_posDef hM hV hε)
          (whitenedBy_posDef (counterfamilyInertia_posDef hM hV hε) hV) dτ
        * (lagW1 (counterfamilyInertia_posDef hM hV hε)
          (whitenedBy_posDef (counterfamilyInertia_posDef hM hV hε) hV) dτ)⁻¹
        * lagW0 (counterfamilyInertia_posDef hM hV hε)
          (whitenedBy_posDef (counterfamilyInertia_posDef hM hV hε) hV) dτ
        - lagW1 (counterfamilyInertia_posDef hM hV hε)
          (whitenedBy_posDef (counterfamilyInertia_posDef hM hV hε) hV) dτ)⁻¹

/-- The two-time reader of the counterfamily chain returns its inertia `M (M + εV)⁻¹ M`
(positive capacities, positive definite `V`, `ε ≥ 0`). -/
theorem counterfamily_history_reader (hM : ∀ i, 0 < M i) (hV : V.PosDef) (hε : 0 ≤ ε)
    {dτ : ℝ} (hd : 0 < dτ) : cfReader hM hV hε dτ = counterfamilyInertia M V ε := by
  unfold cfReader
  exact history_reader_posDef _ _ hd

/-- The counterfamily chain satisfies the slab row of its own precision. -/
theorem counterfamily_contact_identity (hM : ∀ i, 0 < M i) (hV : V.PosDef) (hε : 0 ≤ ε)
    {dτ : ℝ} (hd : 0 < dτ) :
    (2 / dτ) • (counterfamilyInertia M V ε * lagW0 (counterfamilyInertia_posDef hM hV hε)
          (whitenedBy_posDef (counterfamilyInertia_posDef hM hV hε) hV) dτ)
      - (1 / dτ) • (counterfamilyInertia M V ε * (lagW1 (counterfamilyInertia_posDef hM hV hε)
          (whitenedBy_posDef (counterfamilyInertia_posDef hM hV hε) hV) dτ
        + (lagW1 (counterfamilyInertia_posDef hM hV hε)
          (whitenedBy_posDef (counterfamilyInertia_posDef hM hV hε) hV) dτ)ᵀ))
      + dτ • (V * lagW0 (counterfamilyInertia_posDef hM hV hε)
          (whitenedBy_posDef (counterfamilyInertia_posDef hM hV hε) hV) dτ) = 1 :=
  contact_identity_posDef _ _ hd

/-- Two-time correlations distinguish the counterfamily: if `V` couples two distinct nodes, the
two-time reader of the counterfamily chain is node-local exactly when `ε = 0` (for `ε ≥ 0`). -/
theorem counterfamily_reader_isDiag_iff (hM : ∀ i, 0 < M i) (hV : V.PosDef) (hε : 0 ≤ ε)
    (hVne : ∃ i j, i ≠ j ∧ V i j ≠ 0) {dτ : ℝ} (hd : 0 < dτ) :
    (cfReader hM hV hε dτ).IsDiag ↔ ε = 0 := by
  rw [counterfamily_history_reader hM hV hε hd]
  exact counterfamilyInertia_isDiag_iff M V ε (fun i => (hM i).ne')
    ((Matrix.isUnit_iff_isUnit_det _).1 (capacity_shift_posDef hM hV hε).isUnit) hVne

end Counterfamily

end

#print axioms history_reader_whitened
#print axioms contact_identity_whitened
#print axioms history_reader_posDef
#print axioms contact_identity_posDef
#print axioms counterfamilyInertia_posDef
#print axioms counterfamily_history_reader
#print axioms counterfamily_contact_identity
#print axioms counterfamily_reader_isDiag_iff

end OPH.SourceLaw
