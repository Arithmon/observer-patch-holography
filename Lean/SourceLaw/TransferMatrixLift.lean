import SourceLaw.ReflectionPositivity
import Mathlib.Analysis.Matrix.HermitianFunctionalCalculus
import Mathlib.Analysis.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
import Mathlib.Algebra.Group.Pi.Units

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The transfer structure of the stationary chain as matrices

For a symmetric positive definite whitened operator `A`, a function of `A` is taken through its
eigen-decomposition (`Matrix.IsHermitian.cfc`, which acts on bare functions because the spectrum
is finite). With the per-mode transfer eigenvalue and equal-time variance of `TransferStructure`,
`R := r(A)` and `C0 := c0(A)` satisfy, for `Δτ > 0`:

* `transferMat_add_inv`: `R + R⁻¹ = 2 + Δτ² A`;
* `equalTimeCov_mul_one_sub`, `equalTimeCov_eq`: `C0 (1 − R²) = Δτ R`, so `C0 = Δτ R (1 − R²)⁻¹`;
* `equalTimeCov_posDef`, `transferMat_posDef`, `one_sub_transferMat_posDef`, `transferMat_comm`:
  `C0 ≻ 0`, `0 ≺ R ≺ 1`, and `R`, `C0` commute;
* `markov_reader_matrix`: `Δτ (C0 (R C0)⁻¹ C0 − R C0)⁻¹ = 1`, i.e. `(C0 (R C0)⁻¹ C0 − R C0)⁻¹`
  is the time-seam conductance `Δτ⁻¹` in whitened coordinates;
* `contact_identity_matrix`: `(2/Δτ + Δτ A) C0 − (2/Δτ) R C0 = 1`;
* `block_bond_rp`: the block bond reflection Gram `Σ_{a,b} ⟨v_a, C0 R^{t_a+t_b+1} v_b⟩` is
  nonnegative, being the sum over modes of the per-mode bond Grams.
-/

namespace OPH.SourceLaw

open scoped BigOperators
open Matrix

noncomputable section

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {A : Matrix ι ι ℝ}

section Calculus

variable (hA : A.IsHermitian)

theorem conj_mul_conj {U D1 D2 : Matrix ι ι ℝ} (h : star U * U = 1) :
    U * D1 * star U * (U * D2 * star U) = U * (D1 * D2) * star U := by
  calc U * D1 * star U * (U * D2 * star U) = U * D1 * (star U * U) * D2 * star U := by
        simp only [Matrix.mul_assoc]
    _ = U * (D1 * D2) * star U := by
        rw [h, Matrix.mul_one]
        simp only [Matrix.mul_assoc]

theorem eigU_star_mul :
    star (hA.eigenvectorUnitary : Matrix ι ι ℝ) * hA.eigenvectorUnitary = 1 :=
  Unitary.star_mul_self_of_mem hA.eigenvectorUnitary.2

theorem eigU_mul_star :
    (hA.eigenvectorUnitary : Matrix ι ι ℝ) * star (hA.eigenvectorUnitary : Matrix ι ι ℝ) = 1 :=
  Unitary.mul_star_self_of_mem hA.eigenvectorUnitary.2

theorem cfc_expand (f : ℝ → ℝ) :
    hA.cfc f = (hA.eigenvectorUnitary : Matrix ι ι ℝ) * diagonal (f ∘ hA.eigenvalues)
      * star (hA.eigenvectorUnitary : Matrix ι ι ℝ) := by
  rw [IsHermitian.cfc, Unitary.conjStarAlgAut_apply, RCLike.ofReal_real_eq_id]
  rfl

theorem cfc_congr {f g : ℝ → ℝ} (h : ∀ i, f (hA.eigenvalues i) = g (hA.eigenvalues i)) :
    hA.cfc f = hA.cfc g := by
  rw [cfc_expand, cfc_expand, show f ∘ hA.eigenvalues = g ∘ hA.eigenvalues from funext h]

theorem cfc_mul (f g : ℝ → ℝ) : hA.cfc f * hA.cfc g = hA.cfc (fun x => f x * g x) := by
  rw [cfc_expand, cfc_expand, cfc_expand, conj_mul_conj (eigU_star_mul hA),
    diagonal_mul_diagonal]
  rfl

theorem cfc_add (f g : ℝ → ℝ) : hA.cfc f + hA.cfc g = hA.cfc (fun x => f x + g x) := by
  rw [cfc_expand, cfc_expand, cfc_expand, ← Matrix.add_mul, ← Matrix.mul_add, diagonal_add]
  rfl

theorem cfc_sub (f g : ℝ → ℝ) : hA.cfc f - hA.cfc g = hA.cfc (fun x => f x - g x) := by
  rw [cfc_expand, cfc_expand, cfc_expand, ← Matrix.sub_mul, ← Matrix.mul_sub, diagonal_sub]
  rfl

theorem cfc_smul (c : ℝ) (f : ℝ → ℝ) : c • hA.cfc f = hA.cfc (fun x => c * f x) := by
  rw [cfc_expand, cfc_expand, ← Matrix.smul_mul, ← Matrix.mul_smul, ← diagonal_smul]
  rfl

theorem cfc_one : hA.cfc (fun _ => 1) = 1 := by
  rw [cfc_expand]
  have h1 : diagonal ((fun _ : ℝ => (1 : ℝ)) ∘ hA.eigenvalues) = 1 := diagonal_one
  rw [h1, Matrix.mul_one, eigU_mul_star]

theorem cfc_id : hA.cfc id = A := by
  rw [IsHermitian.cfc]
  exact hA.spectral_theorem.symm

theorem smul_eq_cfc (c : ℝ) : c • A = hA.cfc (fun x => c * x) :=
  (congrArg (c • ·) (cfc_id hA)).symm.trans (cfc_smul hA c id)

theorem cfc_pow (f : ℝ → ℝ) (k : ℕ) : hA.cfc f ^ k = hA.cfc (fun x => f x ^ k) := by
  induction k with
  | zero =>
    simp only [pow_zero]
    exact (cfc_one hA).symm
  | succ k ih =>
    rw [pow_succ, ih, cfc_mul]
    simp only [pow_succ]

theorem cfc_inv {f : ℝ → ℝ} (h : ∀ i, f (hA.eigenvalues i) ≠ 0) :
    (hA.cfc f)⁻¹ = hA.cfc (fun x => (f x)⁻¹) := by
  apply Matrix.inv_eq_right_inv
  rw [cfc_mul, ← cfc_one hA]
  exact cfc_congr hA fun i => mul_inv_cancel₀ (h i)

theorem cfc_transpose (f : ℝ → ℝ) : (hA.cfc f)ᵀ = hA.cfc f := by
  rw [cfc_expand, transpose_mul, transpose_mul, diagonal_transpose, star_eq_conjTranspose,
    conjTranspose_eq_transpose_of_trivial, transpose_transpose, Matrix.mul_assoc]

theorem cfc_posDef {f : ℝ → ℝ} (h : ∀ i, 0 < f (hA.eigenvalues i)) : (hA.cfc f).PosDef := by
  rw [cfc_expand, star_eq_conjTranspose]
  have hU : IsUnit (hA.eigenvectorUnitary : Matrix ι ι ℝ) := by
    rw [Matrix.isUnit_iff_isUnit_det]
    have hdet : (hA.eigenvectorUnitary : Matrix ι ι ℝ).det
        * (star (hA.eigenvectorUnitary : Matrix ι ι ℝ)).det = 1 := by
      rw [← det_mul, eigU_mul_star, det_one]
    exact isUnit_iff_ne_zero.2 (left_ne_zero_of_mul_eq_one hdet)
  exact (Matrix.PosDef.diagonal fun i => h i).mul_mul_conjTranspose_same
    (vecMul_injective_of_isUnit hU)

/-- The quadratic form of a function of `A` is a sum over modes. -/
theorem dotProduct_cfc_mulVec (f : ℝ → ℝ) (v w : ι → ℝ) :
    v ⬝ᵥ (hA.cfc f *ᵥ w)
      = ∑ i, (v ᵥ* (hA.eigenvectorUnitary : Matrix ι ι ℝ)) i * f (hA.eigenvalues i)
          * (w ᵥ* (hA.eigenvectorUnitary : Matrix ι ι ℝ)) i := by
  rw [cfc_expand, ← mulVec_mulVec, ← mulVec_mulVec, dotProduct_mulVec]
  have hT : star (hA.eigenvectorUnitary : Matrix ι ι ℝ)
      = (hA.eigenvectorUnitary : Matrix ι ι ℝ)ᵀ := by
    rw [star_eq_conjTranspose, conjTranspose_eq_transpose_of_trivial]
  rw [hT, mulVec_transpose]
  simp only [dotProduct, mulVec_diagonal, Function.comp]
  exact Finset.sum_congr rfl fun i _ => by ring

end Calculus

omit [DecidableEq ι] in
/-- Reordering a mode-diagonal block Gram into per-mode bond Grams. -/
theorem sum_bondGram_reorder {m : ℕ} (t : Fin m → ℕ) (Y : Fin m → ι → ℝ) (C R : ι → ℝ) :
    ∑ a, ∑ b, ∑ i, Y a i * Y b i * (C i * R i ^ (t a + t b + 1))
      = ∑ i, bondGram (C i) (R i) t (fun a => Y a i) := by
  simp only [bondGram]
  conv_rhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Finset.sum_comm]

section Lift

variable (hA : A.PosDef) (dτ : ℝ)

/-- Transfer matrix `R = r(A)`. -/
def transferMat : Matrix ι ι ℝ := hA.1.cfc (transferEig dτ)

/-- Equal-time covariance `C0 = c0(A)` of the stationary chain. -/
def equalTimeCov : Matrix ι ι ℝ := hA.1.cfc (equalTimeVar dτ)

variable {dτ}

theorem eig_pos (i : ι) : 0 < hA.1.eigenvalues i := hA.eigenvalues_pos i

theorem transferMat_inv (hd : 0 < dτ) :
    (transferMat hA dτ)⁻¹ = hA.1.cfc (fun x => (transferEig dτ x)⁻¹) :=
  cfc_inv hA.1 (f := transferEig dτ) fun i => (transferEig_pos hd (eig_pos hA i)).ne'

theorem transferMat_add_inv (hd : 0 < dτ) :
    transferMat hA dτ + (transferMat hA dτ)⁻¹ = (2 : ℝ) • (1 : Matrix ι ι ℝ) + dτ ^ 2 • A := by
  rw [transferMat_inv hA hd, transferMat, cfc_add, ← cfc_one hA.1, cfc_smul, smul_eq_cfc hA.1,
    cfc_add]
  exact cfc_congr hA.1 fun i => by
    rw [transferEig_add_inv hd (eig_pos hA i)]
    ring

theorem transferMat_comm :
    transferMat hA dτ * equalTimeCov hA dτ = equalTimeCov hA dτ * transferMat hA dτ := by
  rw [transferMat, equalTimeCov, cfc_mul, cfc_mul]
  exact cfc_congr hA.1 fun i => mul_comm _ _

theorem equalTimeCov_mul_one_sub (hd : 0 < dτ) :
    equalTimeCov hA dτ * (1 - transferMat hA dτ * transferMat hA dτ) = dτ • transferMat hA dτ := by
  rw [transferMat, equalTimeCov, cfc_mul, ← cfc_one hA.1, cfc_sub, cfc_mul, cfc_smul]
  exact cfc_congr hA.1 fun i => by
    have h := chain_innovation (dτ := dτ) (lam := hA.1.eigenvalues i) hd (eig_pos hA i)
    linear_combination h

theorem one_sub_transferMat_sq_isUnit (hd : 0 < dτ) :
    IsUnit (1 - transferMat hA dτ * transferMat hA dτ) := by
  rw [transferMat, cfc_mul, ← cfc_one hA.1, cfc_sub]
  exact (cfc_posDef hA.1 (f := fun x => 1 - transferEig dτ x * transferEig dτ x) fun i => by
    have h := one_sub_transferEig_sq_pos hd (eig_pos hA i)
    show 0 < 1 - transferEig dτ (hA.1.eigenvalues i) * transferEig dτ (hA.1.eigenvalues i)
    rw [← pow_two]
    exact h).isUnit

theorem equalTimeCov_eq (hd : 0 < dτ) :
    equalTimeCov hA dτ
      = dτ • transferMat hA dτ * (1 - transferMat hA dτ * transferMat hA dτ)⁻¹ := by
  rw [← equalTimeCov_mul_one_sub hA hd, Matrix.mul_assoc,
    Matrix.mul_nonsing_inv _ ((Matrix.isUnit_iff_isUnit_det _).1
      (one_sub_transferMat_sq_isUnit hA hd)), Matrix.mul_one]

theorem equalTimeCov_posDef (hd : 0 < dτ) : (equalTimeCov hA dτ).PosDef :=
  cfc_posDef hA.1 fun i => equalTimeVar_pos hd (eig_pos hA i)

theorem transferMat_posDef (hd : 0 < dτ) : (transferMat hA dτ).PosDef :=
  cfc_posDef hA.1 fun i => transferEig_pos hd (eig_pos hA i)

theorem one_sub_transferMat_posDef (hd : 0 < dτ) : (1 - transferMat hA dτ).PosDef := by
  rw [transferMat, ← cfc_one hA.1, cfc_sub]
  exact cfc_posDef hA.1 (f := fun x => 1 - transferEig dτ x) fun i =>
    sub_pos.2 (transferEig_lt_one hd (eig_pos hA i))

theorem transferMat_transpose : (transferMat hA dτ)ᵀ = transferMat hA dτ := cfc_transpose hA.1 _

theorem equalTimeCov_transpose : (equalTimeCov hA dτ)ᵀ = equalTimeCov hA dτ := cfc_transpose hA.1 _

/-- The innovation covariance `C0 − R C0 Rᵀ` of the stationary chain is invertible (mode by mode it
is `Δτ r > 0`, `chain_innovation`). -/
theorem innovation_isUnit (hd : 0 < dτ) :
    IsUnit (equalTimeCov hA dτ - transferMat hA dτ * equalTimeCov hA dτ * (transferMat hA dτ)ᵀ).det := by
  rw [transferMat_transpose, transferMat, equalTimeCov, cfc_mul, cfc_mul, cfc_sub]
  have hp := cfc_posDef hA.1 (f := fun x => equalTimeVar dτ x
      - transferEig dτ x * equalTimeVar dτ x * transferEig dτ x) fun i => by
    have h := chain_innovation (dτ := dτ) (lam := hA.1.eigenvalues i) hd (eig_pos hA i)
    have hr := transferEig_pos hd (eig_pos hA i)
    have he : equalTimeVar dτ (hA.1.eigenvalues i)
        - transferEig dτ (hA.1.eigenvalues i) * equalTimeVar dτ (hA.1.eigenvalues i)
          * transferEig dτ (hA.1.eigenvalues i) = dτ * transferEig dτ (hA.1.eigenvalues i) := by
      linear_combination h
    show 0 < equalTimeVar dτ (hA.1.eigenvalues i)
        - transferEig dτ (hA.1.eigenvalues i) * equalTimeVar dτ (hA.1.eigenvalues i)
          * transferEig dτ (hA.1.eigenvalues i)
    rw [he]
    exact mul_pos hd hr
  exact (Matrix.isUnit_iff_isUnit_det _).1 hp.isUnit

/-- `(C0 (R C0)⁻¹ C0 − R C0)⁻¹` is the time-seam conductance `1/Δτ` (times the identity) in
whitened coordinates. -/
theorem markov_reader_matrix (hd : 0 < dτ) :
    dτ • (equalTimeCov hA dτ * (transferMat hA dτ * equalTimeCov hA dτ)⁻¹ * equalTimeCov hA dτ
      - transferMat hA dτ * equalTimeCov hA dτ)⁻¹ = 1 := by
  have hr := fun i => transferEig_pos hd (eig_pos hA i)
  have hc := fun i => equalTimeVar_pos hd (eig_pos hA i)
  rw [transferMat, equalTimeCov, cfc_mul,
    cfc_inv hA.1 (f := fun x => transferEig dτ x * equalTimeVar dτ x)
      fun i => (mul_pos (hr i) (hc i)).ne', cfc_mul, cfc_mul, cfc_sub]
  have key : ∀ i, dτ * (equalTimeVar dτ (hA.1.eigenvalues i)
      * (transferEig dτ (hA.1.eigenvalues i) * equalTimeVar dτ (hA.1.eigenvalues i))⁻¹
      * equalTimeVar dτ (hA.1.eigenvalues i)
      - transferEig dτ (hA.1.eigenvalues i) * equalTimeVar dτ (hA.1.eigenvalues i))⁻¹ = 1 := by
    intro i
    have h := markov_reader_unit (dτ := dτ) (lam := hA.1.eigenvalues i) hd (eig_pos hA i)
    have h1 := (hr i).ne'
    have h2 := (hc i).ne'
    have hne : equalTimeVar dτ (hA.1.eigenvalues i)
        * (transferEig dτ (hA.1.eigenvalues i) * equalTimeVar dτ (hA.1.eigenvalues i))⁻¹
        * equalTimeVar dτ (hA.1.eigenvalues i)
        = equalTimeVar dτ (hA.1.eigenvalues i) ^ 2
          / (transferEig dτ (hA.1.eigenvalues i) * equalTimeVar dτ (hA.1.eigenvalues i)) := by
      field_simp
    rw [hne]
    exact h
  have hnz : ∀ i, equalTimeVar dτ (hA.1.eigenvalues i)
      * (transferEig dτ (hA.1.eigenvalues i) * equalTimeVar dτ (hA.1.eigenvalues i))⁻¹
      * equalTimeVar dτ (hA.1.eigenvalues i)
      - transferEig dτ (hA.1.eigenvalues i) * equalTimeVar dτ (hA.1.eigenvalues i) ≠ 0 := by
    intro i h0
    have h := key i
    rw [h0, _root_.inv_zero, mul_zero] at h
    exact zero_ne_one h
  rw [cfc_inv hA.1 (f := fun x => equalTimeVar dτ x * (transferEig dτ x * equalTimeVar dτ x)⁻¹
      * equalTimeVar dτ x - transferEig dτ x * equalTimeVar dτ x) hnz, cfc_smul, ← cfc_one hA.1]
  exact cfc_congr hA.1 key

/-- Contact identity: `(2/Δτ + Δτ A) C0 − (2/Δτ) R C0 = 1`. -/
theorem contact_identity_matrix (hd : 0 < dτ) :
    ((2 / dτ) • (1 : Matrix ι ι ℝ) + dτ • A) * equalTimeCov hA dτ
      - (2 / dτ) • (transferMat hA dτ * equalTimeCov hA dτ) = 1 := by
  rw [transferMat, equalTimeCov, smul_eq_cfc hA.1, ← cfc_one hA.1, cfc_smul, cfc_add, cfc_mul,
    cfc_mul, cfc_smul, cfc_sub]
  exact cfc_congr hA.1 fun i => by
    have h := contact_identity_mode (dτ := dτ) (lam := hA.1.eigenvalues i) hd (eig_pos hA i)
    linear_combination h

/-- Block bond reflection positivity of the stationary chain covariance. -/
theorem block_bond_rp (hd : 0 < dτ) {m : ℕ} (t : Fin m → ℕ) (v : Fin m → ι → ℝ) :
    0 ≤ ∑ a, ∑ b, v a ⬝ᵥ ((equalTimeCov hA dτ * transferMat hA dτ ^ (t a + t b + 1)) *ᵥ v b) := by
  have key : ∑ a, ∑ b, v a ⬝ᵥ ((equalTimeCov hA dτ * transferMat hA dτ ^ (t a + t b + 1)) *ᵥ v b)
      = ∑ i, bondGram (equalTimeVar dτ (hA.1.eigenvalues i)) (transferEig dτ (hA.1.eigenvalues i)) t
          (fun a => (v a ᵥ* (hA.1.eigenvectorUnitary : Matrix ι ι ℝ)) i) := by
    rw [← sum_bondGram_reorder]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun b _ => ?_
    rw [transferMat, equalTimeCov, cfc_pow, cfc_mul, dotProduct_cfc_mulVec]
    exact Finset.sum_congr rfl fun i _ => by ring
  rw [key]
  exact Finset.sum_nonneg fun i _ => chain_bond_rp hd (eig_pos hA i) t _

/-- The core of the two-time reader: `C0 (R C0)⁻¹ C0 − R C0 = Δτ`. -/
theorem reader_core (hd : 0 < dτ) :
    equalTimeCov hA dτ * (transferMat hA dτ * equalTimeCov hA dτ)⁻¹ * equalTimeCov hA dτ
      - transferMat hA dτ * equalTimeCov hA dτ = dτ • (1 : Matrix ι ι ℝ) := by
  have hr := fun i => transferEig_pos hd (eig_pos hA i)
  have hc := fun i => equalTimeVar_pos hd (eig_pos hA i)
  rw [transferMat, equalTimeCov, cfc_mul,
    cfc_inv hA.1 (f := fun x => transferEig dτ x * equalTimeVar dτ x)
      fun i => (mul_pos (hr i) (hc i)).ne', cfc_mul, cfc_mul, cfc_sub, ← cfc_one hA.1, cfc_smul]
  refine cfc_congr hA.1 fun i => ?_
  have h := chain_innovation (dτ := dτ) (lam := hA.1.eigenvalues i) hd (eig_pos hA i)
  have h1 := (hr i).ne'
  have h2 := (hc i).ne'
  show equalTimeVar dτ (hA.1.eigenvalues i)
      * (transferEig dτ (hA.1.eigenvalues i) * equalTimeVar dτ (hA.1.eigenvalues i))⁻¹
      * equalTimeVar dτ (hA.1.eigenvalues i)
      - transferEig dτ (hA.1.eigenvalues i) * equalTimeVar dτ (hA.1.eigenvalues i) = dτ * 1
  field_simp
  linear_combination h

end Lift

section GeneralInertia

/-- `J^{-1/2}` for a positive diagonal time-seam coefficient `J`. -/
def jInvSqrt (J : ι → ℝ) : Matrix ι ι ℝ := diagonal fun i => (√(J i))⁻¹

/-- `J^{1/2}`. -/
def jSqrt (J : ι → ℝ) : Matrix ι ι ℝ := diagonal fun i => √(J i)

variable {J : ι → ℝ} {V : Matrix ι ι ℝ}

theorem jSqrt_mul_jInvSqrt (hJ : ∀ i, 0 < J i) : jSqrt J * jInvSqrt J = 1 := by
  rw [jSqrt, jInvSqrt, diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext i
  exact mul_inv_cancel₀ (Real.sqrt_pos.2 (hJ i)).ne'

theorem jInvSqrt_mul_jSqrt (hJ : ∀ i, 0 < J i) : jInvSqrt J * jSqrt J = 1 := by
  rw [jSqrt, jInvSqrt, diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext i
  exact inv_mul_cancel₀ (Real.sqrt_pos.2 (hJ i)).ne'

theorem jSqrt_mul_jSqrt (hJ : ∀ i, 0 < J i) : jSqrt J * jSqrt J = diagonal J := by
  rw [jSqrt, diagonal_mul_diagonal]
  congr 1
  funext i
  exact Real.mul_self_sqrt (hJ i).le

/-- Whitening by `J` keeps a positive definite operator positive definite. -/
theorem jWhitened_posDef (hJ : ∀ i, 0 < J i) (hV : V.PosDef) :
    (jInvSqrt J * V * jInvSqrt J).PosDef := by
  have hSu : IsUnit (jInvSqrt J) := by
    rw [jInvSqrt, isUnit_diagonal, Pi.isUnit_iff]
    exact fun i => isUnit_iff_ne_zero.2 (inv_ne_zero (Real.sqrt_pos.2 (hJ i)).ne')
  have hSH : (jInvSqrt J)ᴴ = jInvSqrt J := by
    rw [jInvSqrt, diagonal_conjTranspose]
    simp
  have h := hV.conjTranspose_mul_mul_same (mulVec_injective_of_isUnit hSu)
  rw [hSH] at h
  exact h

/-- Equal-time covariance of the stationary chain with time coefficient `J`, in original
coordinates. -/
def lag0 (hA : (jInvSqrt J * V * jInvSqrt J).PosDef) (dτ : ℝ) : Matrix ι ι ℝ :=
  jInvSqrt J * equalTimeCov hA dτ * jInvSqrt J

/-- One-slab covariance `E[f_{t+1} f_tᵀ]` of the same stationary chain, in original coordinates. -/
def lag1 (hA : (jInvSqrt J * V * jInvSqrt J).PosDef) (dτ : ℝ) : Matrix ι ι ℝ :=
  jInvSqrt J * (transferMat hA dτ * equalTimeCov hA dτ) * jInvSqrt J

theorem conj_inv {S Si X : Matrix ι ι ℝ} (h1 : S * Si = 1) (hX : IsUnit X.det) :
    (S * X * S)⁻¹ = Si * X⁻¹ * Si := by
  apply Matrix.inv_eq_right_inv
  calc S * X * S * (Si * X⁻¹ * Si) = S * X * (S * Si) * X⁻¹ * Si := by simp only [Matrix.mul_assoc]
    _ = 1 := by
        rw [h1, Matrix.mul_one, Matrix.mul_assoc S X, Matrix.mul_nonsing_inv _ hX, Matrix.mul_one, h1]

/-- Contact identity in original coordinates for any positive time coefficient `J`:
`(2/Δτ) J C(0) − (1/Δτ) J (C(1) + C(1)ᵀ) + Δτ V C(0) = 1`. -/
theorem contact_identity_general (hJ : ∀ i, 0 < J i) (hA : (jInvSqrt J * V * jInvSqrt J).PosDef)
    {dτ : ℝ} (hd : 0 < dτ) :
    (2 / dτ) • (diagonal J * lag0 hA dτ) - (1 / dτ) • (diagonal J * (lag1 hA dτ + (lag1 hA dτ)ᵀ))
      + dτ • (V * lag0 hA dτ) = 1 := by
  set S := jInvSqrt J
  set Si := jSqrt J
  set C := equalTimeCov hA dτ
  set R := transferMat hA dτ
  have h1 : Si * S = 1 := jSqrt_mul_jInvSqrt hJ
  have h2 : S * Si = 1 := jInvSqrt_mul_jSqrt hJ
  have hl : ∀ X : Matrix ι ι ℝ, Si * (S * X) = X := fun X => by
    rw [← Matrix.mul_assoc, h1, Matrix.one_mul]
  have hSt : Sᵀ = S := diagonal_transpose _
  have hJS : diagonal J = Si * Si := (jSqrt_mul_jSqrt hJ).symm
  have hV : V = Si * (S * V * S) * Si := by
    calc V = (Si * S) * V * (S * Si) := by rw [h1, h2, Matrix.one_mul, Matrix.mul_one]
      _ = Si * (S * V * S) * Si := by simp only [Matrix.mul_assoc]
  have hlag1t : (S * (R * C) * S)ᵀ = S * (R * C) * S := by
    rw [transpose_mul, transpose_mul, transpose_mul, hSt, transferMat_transpose hA,
      equalTimeCov_transpose hA, ← transferMat_comm hA]
    simp only [Matrix.mul_assoc]
    rfl
  have e1 : diagonal J * (S * C * S) = Si * C * S := by
    rw [hJS]; simp only [Matrix.mul_assoc, hl]
  have e2 : diagonal J * (S * (R * C) * S) = Si * (R * C) * S := by
    rw [hJS]; simp only [Matrix.mul_assoc, hl]
  have e3 : V * (S * C * S) = Si * ((S * V * S) * C) * S := by
    conv_lhs => rw [hV]
    simp only [Matrix.mul_assoc, hl]
  have key := contact_identity_matrix hA hd
  have hconj : Si * (((2 / dτ) • (1 : Matrix ι ι ℝ) + dτ • (S * V * S)) * C
      - (2 / dτ) • (R * C)) * S = 1 := by
    rw [key, Matrix.mul_one, h1]
  change (2 / dτ) • (diagonal J * (S * C * S))
      - (1 / dτ) • (diagonal J * (S * (R * C) * S + (S * (R * C) * S)ᵀ))
      + dτ • (V * (S * C * S)) = 1
  rw [hlag1t, Matrix.mul_add, e1, e2, e3, ← hconj]
  simp only [Matrix.add_mul, Matrix.mul_add, Matrix.sub_mul, Matrix.mul_sub, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.one_mul, Matrix.mul_assoc]
  module

/-- The two-time reader of the stationary chain returns the time coefficient for any positive `J`
(the time-seam conductance is `J/Δτ`): `Δτ (C(0) C(1)⁻¹ C(0) − C(1))⁻¹ = J`, where `C(0)`,
`C(1)` are the chain's lag covariances `lag0`, `lag1`, not covariances of the finite register. -/
theorem history_reader_general (hJ : ∀ i, 0 < J i) (hA : (jInvSqrt J * V * jInvSqrt J).PosDef)
    {dτ : ℝ} (hd : 0 < dτ) :
    dτ • (lag0 hA dτ * (lag1 hA dτ)⁻¹ * lag0 hA dτ - lag1 hA dτ)⁻¹ = diagonal J := by
  set S := jInvSqrt J
  set Si := jSqrt J
  set C := equalTimeCov hA dτ
  set R := transferMat hA dτ
  have h1 : Si * S = 1 := jSqrt_mul_jInvSqrt hJ
  have h2 : S * Si = 1 := jInvSqrt_mul_jSqrt hJ
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
  change dτ • ((S * C * S) * (S * (R * C) * S)⁻¹ * (S * C * S) - S * (R * C) * S)⁻¹ = diagonal J
  rw [hinv1, hmid, ← jSqrt_mul_jSqrt hJ]
  have hinv2 : (dτ • (S * S))⁻¹ = dτ⁻¹ • (Si * Si) := by
    apply Matrix.inv_eq_right_inv
    calc dτ • (S * S) * (dτ⁻¹ • (Si * Si)) = (dτ * dτ⁻¹) • (S * (S * Si) * Si) := by
          rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul]
          simp only [Matrix.mul_assoc]
      _ = 1 := by rw [mul_inv_cancel₀ hd.ne', one_smul, h2, Matrix.mul_one, h2]
  rw [hinv2, smul_smul, mul_inv_cancel₀ hd.ne', one_smul]

end GeneralInertia

end

#print axioms transferMat_add_inv
#print axioms equalTimeCov_mul_one_sub
#print axioms equalTimeCov_eq
#print axioms equalTimeCov_posDef
#print axioms transferMat_posDef
#print axioms one_sub_transferMat_posDef
#print axioms transferMat_comm
#print axioms markov_reader_matrix
#print axioms contact_identity_matrix
#print axioms block_bond_rp
#print axioms innovation_isUnit
#print axioms transferMat_transpose
#print axioms reader_core
#print axioms contact_identity_general
#print axioms history_reader_general

end OPH.SourceLaw
