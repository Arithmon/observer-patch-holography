import Mathlib.Algebra.QuadraticAlgebra.Basic
import Mathlib.Data.Real.StarOrdered
import Mathlib.LinearAlgebra.BilinearForm.Orthogonal
import Mathlib.LinearAlgebra.Matrix.BilinearForm
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.NumberTheory.Real.Irrational
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

set_option autoImplicit false
set_option maxHeartbeats 4000000
set_option maxRecDepth 8192

open scoped BigOperators Matrix
open Module

namespace OPH.WhitneyFiniteCertificate

/-!
# Exact finite certificates over Q(sqrt(5))

`Q5` is the exact coefficient field used by the Whitney mass and stability
receipts.  The checker recomputes a congruence factorization and decides the
sign of every diagonal entry by rational arithmetic.  It never imports a
producer-supplied positivity flag.

The certificate factor is not required to be lower triangular.  Existing
unit-lower LDL receipts are accepted as a special case.  Allowing an arbitrary
invertible congruence factor makes certificate existence invariant under lawful
coordinate changes and lets completeness use the algebraic orthogonal-basis
theorem over the same exact field.
-/

private theorem five_not_rat_square :
    ∀ r : ℚ, r ^ 2 ≠ (5 : ℚ) + 0 * r := by
  intro r hr
  have hs : IsSquare (5 : ℚ) := by
    refine ⟨r, ?_⟩
    simpa [pow_two] using hr.symm
  have hsNat : IsSquare (5 : ℕ) := Rat.isSquare_natCast_iff.mp hs
  exact (show Nat.Prime 5 by decide).not_isSquare hsNat

instance q5Irreducible : Fact (∀ r : ℚ, r ^ 2 ≠ (5 : ℚ) + 0 * r) :=
  ⟨five_not_rat_square⟩

/-- Exact coefficients `a + b sqrt(5)`. -/
abbrev Q5 := QuadraticAlgebra ℚ 5 0

/-- The real embedding selecting the positive square root of five. -/
noncomputable def eval : Q5 →+* ℝ where
  toFun z := (z.re : ℝ) + (z.im : ℝ) * Real.sqrt 5
  map_zero' := by norm_num
  map_one' := by
    change ((1 : ℚ) : ℝ) + ((0 : ℚ) : ℝ) * Real.sqrt 5 = 1
    norm_num
  map_add' x y := by
    simp only [QuadraticAlgebra.re_add, QuadraticAlgebra.im_add, Rat.cast_add]
    ring
  map_mul' x y := by
    let s : ℝ := Real.sqrt 5
    change
      (((x.re * y.re + 5 * x.im * y.im : ℚ) : ℚ) : ℝ) +
          (((x.re * y.im + x.im * y.re + 0 * x.im * y.im : ℚ) : ℚ) : ℝ) * s =
        ((x.re : ℝ) + (x.im : ℝ) * s) * ((y.re : ℝ) + (y.im : ℝ) * s)
    push_cast
    have hs : s ^ 2 = (5 : ℝ) := by
      dsimp [s]
      norm_num [Real.sq_sqrt]
    calc
      (x.re : ℝ) * y.re + 5 * (x.im : ℝ) * y.im +
          ((x.re : ℝ) * y.im + (x.im : ℝ) * y.re + 0 * x.im * y.im) * s =
        (x.re : ℝ) * y.re + ((x.re : ℝ) * y.im + x.im * y.re) * s +
          (x.im : ℝ) * y.im * 5 := by ring
      _ = (x.re : ℝ) * y.re + ((x.re : ℝ) * y.im + x.im * y.re) * s +
          (x.im : ℝ) * y.im * s ^ 2 := by rw [hs]
      _ = ((x.re : ℝ) + x.im * s) * ((y.re : ℝ) + y.im * s) := by ring

@[simp] theorem eval_mk (a b : ℚ) :
    eval (QuadraticAlgebra.mk a b : Q5) = (a : ℝ) + (b : ℝ) * Real.sqrt 5 := by
  rfl

theorem eval_injective : Function.Injective eval := by
  exact (eval : Q5 →+* ℝ).injective

instance : DecidableEq Q5 := fun x y =>
  if hre : x.re = y.re then
    if him : x.im = y.im then isTrue (QuadraticAlgebra.ext hre him)
    else isFalse (fun h => him (congrArg QuadraticAlgebra.im h))
  else isFalse (fun h => hre (congrArg QuadraticAlgebra.re h))

/-- Rational sign predicate for `a + b sqrt(5)`. -/
def ExactPositive (x : Q5) : Prop :=
  (x.im = 0 ∧ 0 < x.re) ∨
  (0 < x.im ∧ (0 ≤ x.re ∨ x.re ^ 2 < 5 * x.im ^ 2)) ∨
  (x.im < 0 ∧ 0 < x.re ∧ 5 * x.im ^ 2 < x.re ^ 2)

instance (x : Q5) : Decidable (ExactPositive x) := by
  unfold ExactPositive
  infer_instance

def positive (x : Q5) : Bool := decide (ExactPositive x)

theorem positive_iff_exact (x : Q5) : positive x = true ↔ ExactPositive x := by
  simp [positive]

private theorem sqrt_five_pos : (0 : ℝ) < Real.sqrt 5 := Real.sqrt_pos.2 (by norm_num)

private theorem sqrt_five_sq : (Real.sqrt 5) ^ 2 = (5 : ℝ) := by
  norm_num [Real.sq_sqrt]

/-- The rational sign test is both sound and complete for the real embedding. -/
theorem exactPositive_iff_eval_pos (x : Q5) : ExactPositive x ↔ 0 < eval x := by
  rcases x with ⟨a, b⟩
  simp only [QuadraticAlgebra.im, QuadraticAlgebra.re, eval_mk]
  have hs := sqrt_five_pos
  have hs2 := sqrt_five_sq
  constructor
  · rintro (⟨rfl, ha⟩ | ⟨hb, ha | hsq⟩ | ⟨hb, ha, hsq⟩)
    · norm_num at ha ⊢
      exact_mod_cast ha
    · have hbR : (0 : ℝ) < (b : ℝ) := by exact_mod_cast hb
      have haR : (0 : ℝ) ≤ (a : ℝ) := by exact_mod_cast ha
      positivity
    · have hbR : (0 : ℝ) < (b : ℝ) := by exact_mod_cast hb
      by_cases ha0 : 0 ≤ a
      · have haR : (0 : ℝ) ≤ (a : ℝ) := by exact_mod_cast ha0
        positivity
      · have haR : (a : ℝ) < 0 := by exact_mod_cast (lt_of_not_ge ha0)
        have hsqR : (a : ℝ) ^ 2 < 5 * (b : ℝ) ^ 2 := by exact_mod_cast hsq
        have hlt : -(a : ℝ) < (b : ℝ) * Real.sqrt 5 := by
          rw [← sq_lt_sq₀ (by linarith)
            (mul_nonneg (le_of_lt hbR) (le_of_lt hs))]
          nlinarith
        linarith
    · have hbR : (b : ℝ) < 0 := by exact_mod_cast hb
      have haR : (0 : ℝ) < (a : ℝ) := by exact_mod_cast ha
      have hsqR : 5 * (b : ℝ) ^ 2 < (a : ℝ) ^ 2 := by exact_mod_cast hsq
      have hlt : (-(b : ℝ)) * Real.sqrt 5 < (a : ℝ) := by
        rw [← sq_lt_sq₀
          (mul_nonneg (by linarith) (le_of_lt hs)) (le_of_lt haR)]
        nlinarith
      linarith
  · intro h
    by_cases hb0 : b = 0
    · left
      refine ⟨hb0, ?_⟩
      subst b
      norm_num at h
      exact_mod_cast h
    rcases lt_or_gt_of_ne hb0 with hb | hb
    · right; right
      have hbR : (b : ℝ) < 0 := by exact_mod_cast hb
      have haR : (0 : ℝ) < (a : ℝ) := by nlinarith
      refine ⟨hb, by exact_mod_cast haR, ?_⟩
      have hlt : (-(b : ℝ)) * Real.sqrt 5 < (a : ℝ) := by nlinarith
      have hsqR := (sq_lt_sq₀
        (mul_nonneg (by linarith) (le_of_lt hs)) (le_of_lt haR)).2 hlt
      have : 5 * (b : ℝ) ^ 2 < (a : ℝ) ^ 2 := by nlinarith
      exact_mod_cast this
    · right; left
      refine ⟨hb, ?_⟩
      by_cases ha0 : 0 ≤ a
      · exact Or.inl ha0
      · right
        have hbR : (0 : ℝ) < (b : ℝ) := by exact_mod_cast hb
        have haR : (a : ℝ) < 0 := by exact_mod_cast (lt_of_not_ge ha0)
        have hlt : -(a : ℝ) < (b : ℝ) * Real.sqrt 5 := by nlinarith
        have hsqR := (sq_lt_sq₀ (by linarith : (0 : ℝ) ≤ -(a : ℝ))
          (mul_nonneg (le_of_lt hbR) (le_of_lt hs))).2 hlt
        have : (a : ℝ) ^ 2 < 5 * (b : ℝ) ^ 2 := by nlinarith
        exact_mod_cast this

theorem positive_iff_eval_pos (x : Q5) : positive x = true ↔ 0 < eval x :=
  (positive_iff_exact x).trans (exactPositive_iff_eval_pos x)

noncomputable def evalMatrix {n : Nat} (A : Matrix (Fin n) (Fin n) Q5) :
    Matrix (Fin n) (Fin n) ℝ := A.map eval

structure Certificate (n : Nat) where
  factor : Matrix (Fin n) (Fin n) Q5
  diagonal : Fin n → Q5

def Valid {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (c : Certificate n) : Prop :=
  A = c.factor * Matrix.diagonal c.diagonal * c.factor.transpose ∧
  Matrix.det c.factor ≠ 0 ∧
  ∀ i, ExactPositive (c.diagonal i)

instance {n : Nat} (A : Matrix (Fin n) (Fin n) Q5) (c : Certificate n) :
    Decidable (Valid A c) := by
  unfold Valid
  infer_instance

def check {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (c : Certificate n) : Bool := decide (Valid A c)

theorem check_iff_valid {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (c : Certificate n) : check A c = true ↔ Valid A c := by
  simp [check]

theorem evalMatrix_diagonal {n : Nat} (d : Fin n → Q5) :
    evalMatrix (Matrix.diagonal d) = Matrix.diagonal (fun i => eval (d i)) := by
  ext i j
  by_cases h : i = j <;> simp [evalMatrix, h]

theorem evalMatrix_mul {n : Nat} (A B : Matrix (Fin n) (Fin n) Q5) :
    evalMatrix (A * B) = evalMatrix A * evalMatrix B := by
  simpa [evalMatrix] using
    (Matrix.map_mul (L := A) (M := B) (f := eval))

theorem evalMatrix_transpose {n : Nat} (A : Matrix (Fin n) (Fin n) Q5) :
    evalMatrix A.transpose = (evalMatrix A).transpose := by
  rfl

/-- A checked exact congruence certificate proves Mathlib positive definiteness. -/
theorem valid_sound {n : Nat} {A : Matrix (Fin n) (Fin n) Q5}
    {c : Certificate n} (hc : Valid A c) : (evalMatrix A).PosDef := by
  rcases hc with ⟨hfac, hdet, hpos⟩
  let L : Matrix (Fin n) (Fin n) ℝ := evalMatrix c.factor
  let D : Matrix (Fin n) (Fin n) ℝ := Matrix.diagonal (fun i => eval (c.diagonal i))
  have hD : D.PosDef := by
    rw [Matrix.posDef_diagonal_iff]
    intro i
    exact (exactPositive_iff_eval_pos _).mp (hpos i)
  have hdetL : Matrix.det L ≠ 0 := by
    intro hz
    apply hdet
    apply eval_injective
    rw [eval.map_det]
    simpa [L, evalMatrix] using hz
  have hunitL : IsUnit L := by
    rw [Matrix.isUnit_iff_isUnit_det]
    exact isUnit_iff_ne_zero.mpr hdetL
  have hcong : (L * D * Lᴴ).PosDef :=
    hD.mul_mul_conjTranspose_same (Matrix.vecMul_injective_of_isUnit hunitL)
  rw [hfac]
  simpa [L, D, evalMatrix_mul, evalMatrix_diagonal, evalMatrix_transpose,
    Matrix.conjTranspose_eq_transpose_of_trivial] using hcong

theorem checked_sound {n : Nat} {A : Matrix (Fin n) (Fin n) Q5}
    {c : Certificate n} (hc : check A c = true) : (evalMatrix A).PosDef :=
  valid_sound ((check_iff_valid A c).mp hc)

noncomputable def evalVector {n : Nat} (x : Fin n → Q5) : Fin n → ℝ :=
  fun i => eval (x i)

theorem eval_toBilin' {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (x y : Fin n → Q5) :
    eval (Matrix.toBilin' A x y) =
      dotProduct (evalVector x) ((evalMatrix A).mulVec (evalVector y)) := by
  simp only [Matrix.toBilin'_apply, map_sum, map_mul, dotProduct, Matrix.mulVec,
    evalVector, evalMatrix, Matrix.map_apply]
  apply Finset.sum_congr rfl
  intro i _
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j _
  ring

/-- Positive definiteness supplies an exact certificate over the same Q5
carrier.  The proof uses algebraic orthogonalization, so no real-valued factor
is imported and no fixed square-root enclosure is needed. -/
theorem posDef_complete {n : Nat} (A : Matrix (Fin n) (Fin n) Q5)
    (hA : (evalMatrix A).PosDef) : ∃ c : Certificate n, check A c = true := by
  let B : LinearMap.BilinForm Q5 (Fin n → Q5) := Matrix.toBilin' A
  have hAsymm : A.IsSymm := by
    apply Matrix.IsSymm.ext
    intro i j
    apply eval_injective
    have h := congrFun (congrFun hA.1.eq i) j
    simpa [evalMatrix, Matrix.conjTranspose_apply] using h
  have hBsymm : B.IsSymm := by
    constructor
    intro x y
    simp only [B, Matrix.toBilin'_apply]
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
  obtain ⟨v, hv⟩ : ∃ v : Basis (Fin n) Q5 (Fin n → Q5), B.IsOrthoᵢ v := by
    obtain ⟨v, hv⟩ := LinearMap.BilinForm.exists_orthogonal_basis (B := B)
      (LinearMap.BilinForm.isSymm_iff.mp hBsymm)
    have hdim : finrank Q5 (Fin n → Q5) = n := Module.finrank_fin_fun Q5
    let e : Fin (finrank Q5 (Fin n → Q5)) ≃ Fin n := finCongr hdim
    let v' : Basis (Fin n) Q5 (Fin n → Q5) := v.reindex e
    refine ⟨v', ?_⟩
    intro i j hij
    have hne : e.symm i ≠ e.symm j := fun h => hij (e.symm.injective h)
    simpa [v', Basis.reindex_apply] using hv hne
  let b : Basis (Fin n) Q5 (Fin n → Q5) := Pi.basisFun Q5 (Fin n)
  let U : Matrix (Fin n) (Fin n) Q5 := b.toMatrix v
  let V : Matrix (Fin n) (Fin n) Q5 := v.toMatrix b
  let d : Fin n → Q5 := fun i => B (v i) (v i)
  have hdiag : LinearMap.BilinForm.toMatrix v B = Matrix.diagonal d := by
    funext i j
    by_cases hij : i = j
    · subst j
      simp [d]
    · have hij0 : B (v i) (v j) = 0 := hv hij
      simp [LinearMap.BilinForm.toMatrix_apply, d, hij, hij0]
  have hbase : LinearMap.BilinForm.toMatrix b B = A := by
    simpa [b, B, LinearMap.BilinForm.toMatrix_basisFun] using
      (LinearMap.BilinForm.toMatrix'_toBilin' A)
  have hcong : U.transpose * A * U = Matrix.diagonal d := by
    calc
      U.transpose * A * U =
          U.transpose * LinearMap.BilinForm.toMatrix b B * U := by rw [hbase]
      _ = LinearMap.BilinForm.toMatrix v B := by
        exact LinearMap.BilinForm.toMatrix_mul_basis_toMatrix (b := b) v B
      _ = Matrix.diagonal d := hdiag
  have hUV : U * V = 1 := by
    simpa [U, V, b] using Basis.toMatrix_mul_toMatrix_flip b v
  have hfactor : A = V.transpose * Matrix.diagonal d * V := by
    symm
    calc
      V.transpose * Matrix.diagonal d * V =
          V.transpose * (U.transpose * A * U) * V := by rw [hcong]
      _ = (U * V).transpose * A * (U * V) := by
        simp [Matrix.transpose_mul, Matrix.mul_assoc]
      _ = A := by rw [hUV]; simp
  have hdetV : IsUnit (Matrix.det V) :=
    Matrix.isUnit_det_of_left_inverse hUV
  have hdetFactor : Matrix.det V.transpose ≠ 0 :=
    (Matrix.isUnit_det_transpose V hdetV).ne_zero
  have hd : ∀ i, ExactPositive (d i) := by
    intro i
    apply (exactPositive_iff_eval_pos _).mpr
    rw [eval_toBilin']
    apply hA.dotProduct_mulVec_pos
    intro hz
    have hv0 : v i = 0 := by
      funext j
      have heval := congrFun hz j
      change eval (v i j) = 0 at heval
      exact (map_eq_zero_iff eval eval_injective).mp heval
    exact (v.ne_zero i) hv0
  refine ⟨⟨V.transpose, d⟩, (check_iff_valid A _).2 ?_⟩
  exact ⟨hfactor, hdetFactor, hd⟩

/-- General checker semantics and certificate existence completeness on every
supported finite dimension. -/
theorem posDef_iff_exists_certificate {n : Nat}
    (A : Matrix (Fin n) (Fin n) Q5) :
    (evalMatrix A).PosDef ↔ ∃ c : Certificate n, check A c = true := by
  constructor
  · exact posDef_complete A
  · rintro ⟨c, hc⟩
    exact checked_sound hc

#print axioms exactPositive_iff_eval_pos
#print axioms valid_sound
#print axioms checked_sound
#print axioms posDef_complete
#print axioms posDef_iff_exists_certificate

end OPH.WhitneyFiniteCertificate
