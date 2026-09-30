import WhitneyFiniteCertificate
import Mathlib.LinearAlgebra.Matrix.Block

set_option autoImplicit false
set_option maxHeartbeats 4000000
set_option maxRecDepth 8192

open scoped BigOperators Matrix
open OrderDual

namespace OPH.WhitneyAlgebraicLDL

open OPH.WhitneyFiniteCertificate

noncomputable section

/-- The standard coordinate vector, kept explicit so the algebraic
Gram--Schmidt construction does not import a real inner-product structure. -/
def standard {n : Nat} (i : Fin n) : Fin n → Q5 :=
  fun j => if j = i then 1 else 0

/-- Algebraic Gram--Schmidt in the original coordinate order.  For a positive
definite real embedding, every denominator is nonzero; the construction itself
stays entirely in `Q5`. -/
noncomputable def ortho {n : Nat}
    (B : LinearMap.BilinForm Q5 (Fin n → Q5)) : Fin n → (Fin n → Q5) :=
  wellFounded_lt.fix fun i previous =>
    standard i - ∑ j : Finset.Iio i,
      (B (previous j (Finset.mem_Iio.mp j.property)) (standard i) /
        B (previous j (Finset.mem_Iio.mp j.property))
          (previous j (Finset.mem_Iio.mp j.property))) •
        previous j (Finset.mem_Iio.mp j.property)

theorem ortho_eq {n : Nat} (B : LinearMap.BilinForm Q5 (Fin n → Q5))
    (i : Fin n) :
    ortho B i = standard i - ∑ j : Finset.Iio i,
      (B (ortho B j) (standard i) / B (ortho B j) (ortho B j)) • ortho B j := by
  rw [ortho, wellFounded_lt.fix_eq]

theorem ortho_above {n : Nat} (B : LinearMap.BilinForm Q5 (Fin n → Q5)) :
    ∀ (i j : Fin n), i < j → ortho B i j = 0
  | i, j, hij => by
    rw [ortho_eq]
    simp only [Pi.sub_apply, Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
    have hstd : standard i j = 0 := by simp [standard, ne_of_gt hij]
    rw [hstd, zero_sub]
    apply neg_eq_zero.mpr
    apply Finset.sum_eq_zero
    intro k _
    have hki : (k : Fin n) < i := Finset.mem_Iio.mp k.property
    rw [ortho_above B (k : Fin n) j (hki.trans hij)]
    simp
termination_by i => i
decreasing_by exact hki

theorem ortho_diag {n : Nat} (B : LinearMap.BilinForm Q5 (Fin n → Q5))
    (i : Fin n) : ortho B i i = 1 := by
  rw [ortho_eq]
  simp only [Pi.sub_apply, Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
  have hstd : standard i i = 1 := by simp [standard]
  rw [hstd]
  have hsum :
      (∑ j : Finset.Iio i,
        (B (ortho B (j : Fin n)) (standard i) /
          B (ortho B (j : Fin n)) (ortho B (j : Fin n))) *
          ortho B (j : Fin n) i) = 0 := by
    apply Finset.sum_eq_zero
    intro j _
    rw [ortho_above B (j : Fin n) i (Finset.mem_Iio.mp j.property)]
    simp
  rw [hsum, sub_zero]

theorem ortho_ne_zero {n : Nat} (B : LinearMap.BilinForm Q5 (Fin n → Q5))
    (i : Fin n) : ortho B i ≠ 0 := by
  intro h
  have hi := congrFun h i
  rw [ortho_diag B i] at hi
  exact one_ne_zero hi

theorem ortho_pair {n : Nat} (B : LinearMap.BilinForm Q5 (Fin n → Q5))
    (hB : B.IsSymm) (hself : ∀ i, B (ortho B i) (ortho B i) ≠ 0) :
    ∀ (a b : Fin n), a < b → B (ortho B a) (ortho B b) = 0 := by
  intro a b hab
  let wf : WellFounded ((· < ·) : Fin n → Fin n → Prop) := wellFounded_lt
  induction b using wf.induction generalizing a with
  | h b ih =>
    rw [ortho_eq B b, LinearMap.BilinForm.sub_right]
    simp only [map_sum, map_smul, smul_eq_mul]
    let a' : Finset.Iio b := ⟨a, Finset.mem_Iio.mpr hab⟩
    have hsum :
        (∑ j : Finset.Iio b,
          (B (ortho B j) (standard b) / B (ortho B j) (ortho B j)) *
            B (ortho B a) (ortho B j)) = B (ortho B a) (standard b) := by
      rw [Fintype.sum_eq_single a']
      · dsimp [a']
        field_simp [hself a]
      · intro j hja
        have hne : (j : Fin n) ≠ a := by
          intro h
          apply hja
          apply Subtype.ext
          exact h
        have hzero : B (ortho B a) (ortho B (j : Fin n)) = 0 := by
          rcases lt_or_gt_of_ne hne with hja' | haj'
          · rw [hB.eq]
            exact ih a hab (j : Fin n) hja'
          · exact ih (j : Fin n) (Finset.mem_Iio.mp j.property) a haj'
        rw [hzero, mul_zero]
    rw [hsum, sub_self]

def columnMatrix {n : Nat} (B : LinearMap.BilinForm Q5 (Fin n → Q5)) :
    Matrix (Fin n) (Fin n) Q5 := fun i j => ortho B j i

theorem columnMatrix_upper {n : Nat}
    (B : LinearMap.BilinForm Q5 (Fin n → Q5)) :
    (columnMatrix B).BlockTriangular id := by
  intro i j hij
  exact ortho_above B j i hij

theorem columnMatrix_diag {n : Nat}
    (B : LinearMap.BilinForm Q5 (Fin n → Q5)) (i : Fin n) :
    columnMatrix B i i = 1 := ortho_diag B i

theorem columnMatrix_det {n : Nat}
    (B : LinearMap.BilinForm Q5 (Fin n → Q5)) :
    Matrix.det (columnMatrix B) = 1 := by
  rw [Matrix.det_of_upperTriangular (columnMatrix_upper B)]
  simp [columnMatrix_diag]

theorem columnMatrix_isUnit {n : Nat}
    (B : LinearMap.BilinForm Q5 (Fin n → Q5)) :
    IsUnit (columnMatrix B) := by
  rw [Matrix.isUnit_iff_isUnit_det, columnMatrix_det]
  exact isUnit_one

def matrixForm {n : Nat} (A : Matrix (Fin n) (Fin n) Q5) :
    LinearMap.BilinForm Q5 (Fin n → Q5) := Matrix.toBilin' A

theorem matrixForm_symm {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (hA : (evalMatrix A).PosDef) : (matrixForm A).IsSymm := by
  have hAsymm : A.IsSymm := by
    apply Matrix.IsSymm.ext
    intro i j
    apply eval_injective
    have h := congrFun (congrFun hA.1.eq i) j
    simpa [evalMatrix, Matrix.conjTranspose_apply] using h
  constructor
  intro x y
  simp only [matrixForm, Matrix.toBilin'_apply]
  calc
    (∑ i, ∑ j, x i * A i j * y j) =
        ∑ j, ∑ i, y j * A j i * x i := by
          rw [Finset.sum_comm]
          apply Finset.sum_congr rfl
          intro i _
          apply Finset.sum_congr rfl
          intro j _
          rw [hAsymm.apply]
          ring
    _ = ∑ i, ∑ j, y i * A i j * x j := rfl

theorem ortho_self_exactPositive {n : Nat}
    (A : Matrix (Fin n) (Fin n) Q5) (hA : (evalMatrix A).PosDef) (i : Fin n) :
    ExactPositive (matrixForm A (ortho (matrixForm A) i) (ortho (matrixForm A) i)) := by
  apply (exactPositive_iff_eval_pos _).mpr
  change 0 < eval (Matrix.toBilin' A (ortho (matrixForm A) i) (ortho (matrixForm A) i))
  rw [eval_toBilin']
  apply hA.dotProduct_mulVec_pos
  intro hz
  have hv0 : ortho (matrixForm A) i = 0 := by
    funext j
    have heval := congrFun hz j
    change eval (ortho (matrixForm A) i j) = 0 at heval
    exact (map_eq_zero_iff eval eval_injective).mp heval
  exact (ortho_ne_zero (matrixForm A) i) hv0

theorem ortho_self_ne_zero {n : Nat}
    (A : Matrix (Fin n) (Fin n) Q5) (hA : (evalMatrix A).PosDef) (i : Fin n) :
    matrixForm A (ortho (matrixForm A) i) (ortho (matrixForm A) i) ≠ 0 := by
  intro hz
  have hp := (exactPositive_iff_eval_pos _).mp (ortho_self_exactPositive A hA i)
  rw [hz, map_zero] at hp
  exact lt_irrefl 0 hp

theorem ortho_pair_matrixForm {n : Nat}
    (A : Matrix (Fin n) (Fin n) Q5) (hA : (evalMatrix A).PosDef)
    (i j : Fin n) (hij : i ≠ j) :
    matrixForm A (ortho (matrixForm A) i) (ortho (matrixForm A) j) = 0 := by
  rcases lt_or_gt_of_ne hij with hij' | hji'
  · exact ortho_pair (matrixForm A) (matrixForm_symm A hA)
      (ortho_self_ne_zero A hA) i j hij'
  · rw [(matrixForm_symm A hA).eq]
    exact ortho_pair (matrixForm A) (matrixForm_symm A hA)
      (ortho_self_ne_zero A hA) j i hji'

def orthoDiagonal {n : Nat} (A : Matrix (Fin n) (Fin n) Q5) : Fin n → Q5 :=
  fun i => matrixForm A (ortho (matrixForm A) i) (ortho (matrixForm A) i)

theorem column_congruence {n : Nat} (A : Matrix (Fin n) (Fin n) Q5) :
    (columnMatrix (matrixForm A)).transpose * A * columnMatrix (matrixForm A) =
      fun i j => matrixForm A (ortho (matrixForm A) i) (ortho (matrixForm A) j) := by
  funext i j
  simp only [Matrix.mul_apply, Matrix.transpose_apply, columnMatrix, matrixForm,
    Matrix.toBilin'_apply]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro l _
  rw [Finset.sum_mul]

theorem column_congruence_diagonal {n : Nat}
    (A : Matrix (Fin n) (Fin n) Q5) (hA : (evalMatrix A).PosDef) :
    (columnMatrix (matrixForm A)).transpose * A * columnMatrix (matrixForm A) =
      Matrix.diagonal (orthoDiagonal A) := by
  rw [column_congruence]
  funext i j
  by_cases hij : i = j
  · subst j
    simp [orthoDiagonal]
  · rw [Matrix.diagonal_apply_ne _ hij]
    exact ortho_pair_matrixForm A hA i j hij

theorem upper_inverse_diag_one {n : Nat} (M : Matrix (Fin n) (Fin n) Q5)
    [Invertible M] (hupper : M.BlockTriangular id) (hdiag : ∀ i, M i i = 1)
    (i : Fin n) : M⁻¹ i i = 1 := by
  have hupperInv : M⁻¹.BlockTriangular id :=
    Matrix.blockTriangular_inv_of_blockTriangular hupper
  have hprod := congrFun (congrFun (Matrix.inv_mul_of_invertible M) i) i
  simp only [Matrix.mul_apply, Matrix.one_apply, ↓reduceIte] at hprod
  rw [Fintype.sum_eq_single i] at hprod
  · simpa [hdiag i] using hprod
  · intro k hki
    rcases lt_or_gt_of_ne hki with hki' | hik'
    · rw [hupperInv hki', zero_mul]
    · rw [hupper hik', mul_zero]

/-- The external receipt schema: a unit-lower factor, an exact diagonal,
positive pivots, and the exact `A = L D Lᵀ` identity. -/
def UnitLower {n : Nat} (L : Matrix (Fin n) (Fin n) Q5) : Prop :=
  L.BlockTriangular toDual ∧ ∀ i, L i i = 1

instance {n : Nat} (L : Matrix (Fin n) (Fin n) Q5) : Decidable (UnitLower L) := by
  unfold UnitLower Matrix.BlockTriangular
  infer_instance

def LDLValid {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (c : Certificate n) : Prop := Valid A c ∧ UnitLower c.factor

instance {n : Nat} (A : Matrix (Fin n) (Fin n) Q5) (c : Certificate n) :
    Decidable (LDLValid A c) := by
  unfold LDLValid
  infer_instance

def ldlCheck {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (c : Certificate n) : Bool := decide (LDLValid A c)

theorem ldlCheck_iff_valid {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (c : Certificate n) : ldlCheck A c = true ↔ LDLValid A c := by
  simp [ldlCheck]

theorem ldl_checked_sound {n : Nat} {A : Matrix (Fin n) (Fin n) Q5}
    {c : Certificate n} (hc : ldlCheck A c = true) : (evalMatrix A).PosDef :=
  valid_sound ((ldlCheck_iff_valid A c).mp hc).1

theorem posDef_ldl_complete {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (hA : (evalMatrix A).PosDef) : ∃ c : Certificate n, ldlCheck A c = true := by
  let B := matrixForm A
  let U := columnMatrix B
  let d := orthoDiagonal A
  letI : Invertible U := (columnMatrix_isUnit B).invertible
  let V : Matrix (Fin n) (Fin n) Q5 := U⁻¹
  let L : Matrix (Fin n) (Fin n) Q5 := V.transpose
  have hUupper : U.BlockTriangular id := columnMatrix_upper B
  have hVupper : V.BlockTriangular id := by
    exact Matrix.blockTriangular_inv_of_blockTriangular hUupper
  have hVdiag : ∀ i, V i i = 1 := by
    intro i
    exact upper_inverse_diag_one U hUupper (columnMatrix_diag B) i
  have hLunit : UnitLower L := by
    constructor
    · intro i j hij
      exact hVupper hij
    · intro i
      exact hVdiag i
  have hdiag : U.transpose * A * U = Matrix.diagonal d := by
    exact column_congruence_diagonal A hA
  have hfac : A = L * Matrix.diagonal d * L.transpose := by
    symm
    calc
      L * Matrix.diagonal d * L.transpose =
          U⁻¹.transpose * (U.transpose * A * U) * (U⁻¹.transpose).transpose := by
            simp only [L, V, hdiag]
      _ = (U⁻¹.transpose * U.transpose) * A * (U * U⁻¹) := by
        noncomm_ring
      _ = (U * U⁻¹).transpose * A * (U * U⁻¹) := by
        rw [Matrix.transpose_mul]
      _ = A := by rw [Matrix.mul_inv_of_invertible]; simp
  have hdetL : Matrix.det L ≠ 0 := by
    rw [Matrix.det_of_lowerTriangular L hLunit.1]
    simp [hLunit.2]
  have hd : ∀ i, ExactPositive (d i) := by
    intro i
    exact ortho_self_exactPositive A hA i
  refine ⟨⟨L, d⟩, (ldlCheck_iff_valid A _).2 ?_⟩
  exact ⟨⟨hfac, hdetL, hd⟩, hLunit⟩

theorem posDef_iff_exists_ldl_certificate {n : Nat}
    (A : Matrix (Fin n) (Fin n) Q5) :
    (evalMatrix A).PosDef ↔ ∃ c : Certificate n, ldlCheck A c = true := by
  constructor
  · exact posDef_ldl_complete A
  · rintro ⟨c, hc⟩
    exact ldl_checked_sound hc

#print axioms posDef_ldl_complete
#print axioms posDef_iff_exists_ldl_certificate

end

end OPH.WhitneyAlgebraicLDL
