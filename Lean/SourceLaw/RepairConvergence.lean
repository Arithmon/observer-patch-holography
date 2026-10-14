import SourceLaw.MatrixLimits
import SourceLaw.LangevinStationarity
import Mathlib.Analysis.SpecialFunctions.ExpDeriv

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Convergence of the moment flow of the declared law

The declared law `df = −𝕄⁻¹ P f ds + √(2ℏ 𝕄⁻¹) dW` (sample time `s`, not physical time) moves the
mean by `E(s) = exp(−s B)`, `B = 𝕄⁻¹ P`, and the covariance by the moment equation
`dΣ/ds = −BΣ − ΣBᵀ + 2ℏ 𝕄⁻¹`. For positive capacities `𝕄` and positive definite `P`, with
`Ã = 𝕄^{-1/2} P 𝕄^{-1/2}`, this file builds `E(s) = 𝕄^{-1/2} exp(−s Ã) 𝕄^{1/2}` through the
eigen-decomposition of `Ã` and proves:

* `flow_zero`, `flow_add`, `flow_hasDerivAt`: `E(0) = 1`, `E(s + t) = E(s) E(t)`, and entry by entry
  `dE/ds = −B E(s)`.
* `flow_tendsto_zero`, `mean_tendsto`: `E(s) → 0` and `E(s) m₀ → 0` as `s → ∞`.
* `covFlow_zero`, `covFlow_hasDerivAt`, `covFlow_tendsto`: for every initial matrix `Σ₀`,
  `Σ(s) = E(s)(Σ₀ − ℏP⁻¹)E(s)ᵀ + ℏP⁻¹` has `Σ(0) = Σ₀`, satisfies the moment equation entry by
  entry, and tends to `ℏP⁻¹`.

Derivatives are stated entry by entry (`entry_hasDerivAt_mul`), so no norm on matrices is used.
-/

namespace OPH.SourceLaw

open Filter Topology Matrix

noncomputable section

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

section Entry

omit [DecidableEq ι] in
/-- Entrywise product rule for matrix-valued functions of one real variable. -/
theorem entry_hasDerivAt_mul {F G : ℝ → Matrix ι ι ℝ} {F' G' : Matrix ι ι ℝ} {s : ℝ}
    (hF : ∀ i j, HasDerivAt (fun t => F t i j) (F' i j) s)
    (hG : ∀ i j, HasDerivAt (fun t => G t i j) (G' i j) s) (i j : ι) :
    HasDerivAt (fun t => (F t * G t) i j) ((F' * G s + F s * G') i j) s := by
  simp only [Matrix.mul_apply, Matrix.add_apply]
  rw [← Finset.sum_add_distrib]
  exact HasDerivAt.fun_sum fun k _ => (hF i k).mul (hG k j)

/-- Entrywise derivative of `exp(−s A)` built from the eigen-decomposition. -/
theorem cfc_exp_hasDerivAt {A : Matrix ι ι ℝ} (hA : A.IsHermitian) (s : ℝ) (i j : ι) :
    HasDerivAt (fun t => hA.cfc (fun x => Real.exp (-(t * x))) i j)
      (hA.cfc (fun x => -x * Real.exp (-(s * x))) i j) s := by
  have e : ∀ f : ℝ → ℝ, hA.cfc f i j = ∑ k, (hA.eigenvectorUnitary : Matrix ι ι ℝ) i k
      * f (hA.eigenvalues k) * star (hA.eigenvectorUnitary : Matrix ι ι ℝ) k j := by
    intro f
    rw [cfc_expand hA, Matrix.mul_apply]
    simp only [Matrix.mul_diagonal, Function.comp]
  simp only [e]
  refine HasDerivAt.fun_sum fun k _ => ?_
  have h1 : HasDerivAt (fun t => Real.exp (-(t * hA.eigenvalues k)))
      (Real.exp (-(s * hA.eigenvalues k)) * (-(1 * hA.eigenvalues k))) s :=
    ((hasDerivAt_id s).mul_const (hA.eigenvalues k)).neg.exp
  convert (h1.const_mul ((hA.eigenvectorUnitary : Matrix ι ι ℝ) i k)).mul_const
    (star (hA.eigenvectorUnitary : Matrix ι ι ℝ) k j) using 1
  ring

end Entry

section Flow

variable {M : ι → ℝ} {P : Matrix ι ι ℝ}

/-- The mean flow `E(s) = 𝕄^{-1/2} exp(−s Ã) 𝕄^{1/2}`. -/
def flow (M : ι → ℝ) (hM : ∀ i, 0 < M i) (hP : P.PosDef) (s : ℝ) : Matrix ι ι ℝ :=
  jInvSqrt M * (jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-(s * x))) * jSqrt M

variable (hM : ∀ i, 0 < M i) (hP : P.PosDef)

theorem flow_zero : flow M hM hP 0 = 1 := by
  have h : (jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-(0 * x))) = 1 := by
    rw [← cfc_one (jWhitened_posDef hM hP).1]
    exact cfc_congr _ fun i => by simp
  rw [flow, h, Matrix.mul_one, jInvSqrt_mul_jSqrt hM]

theorem flow_add (s t : ℝ) : flow M hM hP (s + t) = flow M hM hP s * flow M hM hP t := by
  have hc : (jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-((s + t) * x)))
      = (jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-(s * x)))
        * (jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-(t * x))) := by
    rw [cfc_mul (jWhitened_posDef hM hP).1]
    exact cfc_congr _ fun i => by rw [← Real.exp_add]; congr 1; ring
  have h1 : jSqrt M * jInvSqrt M = 1 := jSqrt_mul_jInvSqrt hM
  rw [flow, flow, flow, hc]
  calc jInvSqrt M * ((jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-(s * x)))
        * (jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-(t * x)))) * jSqrt M
      = jInvSqrt M * (jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-(s * x)))
          * (jSqrt M * jInvSqrt M) * (jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-(t * x)))
          * jSqrt M := by
        rw [h1, Matrix.mul_one]; simp only [Matrix.mul_assoc]
    _ = jInvSqrt M * (jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-(s * x))) * jSqrt M
          * (jInvSqrt M * (jWhitened_posDef hM hP).1.cfc (fun x => Real.exp (-(t * x)))
            * jSqrt M) := by
        simp only [Matrix.mul_assoc]

/-- Entry by entry, `dE/ds = −B E(s)` with `B = 𝕄⁻¹ P`. -/
theorem flow_hasDerivAt (s : ℝ) (i j : ι) :
    HasDerivAt (fun t => flow M hM hP t i j) ((-(driftMatrix M P) * flow M hM hP s) i j) s := by
  set hA := jWhitened_posDef hM hP
  have hent : ∀ C : Matrix ι ι ℝ, (jInvSqrt M * C * jSqrt M) i j = (√(M i))⁻¹ * C i j * √(M j) := by
    intro C
    rw [jInvSqrt, jSqrt, Matrix.mul_diagonal, Matrix.diagonal_mul]
  have hc : hA.1.cfc (fun x => -x * Real.exp (-(s * x)))
      = -((jInvSqrt M * P * jInvSqrt M) * hA.1.cfc (fun x => Real.exp (-(s * x)))) := by
    have h1 := cfc_mul hA.1 id (fun x => Real.exp (-(s * x)))
    rw [cfc_id hA.1] at h1
    have h2 : hA.1.cfc (fun x => -x * Real.exp (-(s * x)))
        = (-1 : ℝ) • hA.1.cfc (fun x => id x * Real.exp (-(s * x))) := by
      rw [cfc_smul hA.1]
      exact cfc_congr hA.1 fun i => by simp only [id]; ring
    rw [h2, ← h1, neg_one_smul]
  have hmat : jInvSqrt M * hA.1.cfc (fun x => -x * Real.exp (-(s * x))) * jSqrt M
      = -(driftMatrix M P) * flow M hM hP s := by
    rw [hc, driftMatrix, flow, ← jInvSqrt_mul_self hM]
    simp only [Matrix.mul_neg, Matrix.neg_mul, Matrix.mul_assoc]
  rw [← hmat, hent]
  have e : (fun t => flow M hM hP t i j)
      = fun t => (√(M i))⁻¹ * hA.1.cfc (fun x => Real.exp (-(t * x))) i j * √(M j) := by
    funext t
    rw [flow, hent]
  rw [e]
  exact ((cfc_exp_hasDerivAt hA.1 s i j).const_mul _).mul_const _

/-- The mean flow tends to zero. -/
theorem flow_tendsto_zero : Tendsto (fun s => flow M hM hP s) atTop (𝓝 0) := by
  set hA := jWhitened_posDef hM hP
  have h0 : hA.1.cfc (fun _ => (0 : ℝ)) = 0 := by
    have h := cfc_smul hA.1 0 (fun _ => 1)
    rw [zero_smul] at h
    rw [h]
    exact cfc_congr hA.1 fun i => by ring
  have hC : Tendsto (fun s => hA.1.cfc (fun x => Real.exp (-(s * x)))) atTop (𝓝 0) := by
    rw [← h0]
    exact cfc_tendsto hA.1 (f := fun s x => Real.exp (-(s * x))) fun i =>
      Real.tendsto_exp_neg_atTop_nhds_zero.comp (tendsto_id.atTop_mul_const (hA.eigenvalues_pos i))
  have h := ((tendsto_const_nhds (x := jInvSqrt M)).mul hC).mul (tendsto_const_nhds (x := jSqrt M))
  rw [Matrix.mul_zero, Matrix.zero_mul] at h
  exact h

/-- Every mean tends to zero under the flow. -/
theorem mean_tendsto (m0 : ι → ℝ) : Tendsto (fun s => flow M hM hP s *ᵥ m0) atTop (𝓝 0) := by
  have hc : Continuous fun X : Matrix ι ι ℝ => X *ᵥ m0 :=
    continuous_id.matrix_mulVec continuous_const
  have h := (hc.tendsto 0).comp (flow_tendsto_zero hM hP)
  rw [Matrix.zero_mulVec] at h
  exact h

/-- The covariance flow `Σ(s) = E(s)(Σ₀ − ℏP⁻¹)E(s)ᵀ + ℏP⁻¹`. -/
def covFlow (hbar : ℝ) (S0 : Matrix ι ι ℝ) (s : ℝ) : Matrix ι ι ℝ :=
  flow M hM hP s * (S0 - stationaryCov hbar P) * (flow M hM hP s)ᵀ + stationaryCov hbar P

theorem covFlow_zero (hbar : ℝ) (S0 : Matrix ι ι ℝ) : covFlow hM hP hbar S0 0 = S0 := by
  rw [covFlow, flow_zero hM hP, transpose_one, Matrix.one_mul, Matrix.mul_one, sub_add_cancel]

/-- Entry by entry, `Σ(s)` solves the moment equation `dΣ/ds = −BΣ − ΣBᵀ + 2ℏ𝕄⁻¹`. -/
theorem covFlow_hasDerivAt (hbar : ℝ) (S0 : Matrix ι ι ℝ) (s : ℝ) (i j : ι) :
    HasDerivAt (fun t => covFlow hM hP hbar S0 t i j)
      ((-(driftMatrix M P) * covFlow hM hP hbar S0 s
        - covFlow hM hP hbar S0 s * (driftMatrix M P)ᵀ
        + (2 * hbar) • diagonal (fun i => (M i)⁻¹)) i j) s := by
  have hPU : IsUnit P.det := (Matrix.isUnit_iff_isUnit_det _).1 hP.isUnit
  have hPs : Pᵀ = P := by
    have := hP.isHermitian
    rwa [Matrix.IsHermitian, conjTranspose_eq_transpose_of_trivial] at this
  have hE := flow_hasDerivAt hM hP s
  have hX : ∀ a b, HasDerivAt (fun _ : ℝ => (S0 - stationaryCov hbar P) a b)
      ((0 : Matrix ι ι ℝ) a b) s := fun a b => hasDerivAt_const s _
  have hET : ∀ a b, HasDerivAt (fun t => (flow M hM hP t)ᵀ a b)
      ((-(driftMatrix M P) * flow M hM hP s)ᵀ a b) s := fun a b => hE b a
  have h1 := entry_hasDerivAt_mul (G := fun _ => S0 - stationaryCov hbar P) hE hX
  have h2 := entry_hasDerivAt_mul
    (F := fun t => flow M hM hP t * (S0 - stationaryCov hbar P))
    (G := fun t => (flow M hM hP t)ᵀ) h1 hET i j
  have h3 := h2.add (hasDerivAt_const s (stationaryCov hbar P i j))
  have key : -(driftMatrix M P) * covFlow hM hP hbar S0 s
        - covFlow hM hP hbar S0 s * (driftMatrix M P)ᵀ
        + (2 * hbar) • diagonal (fun i => (M i)⁻¹)
      = ((-(driftMatrix M P) * flow M hM hP s) * (S0 - stationaryCov hbar P)
          + flow M hM hP s * 0) * (flow M hM hP s)ᵀ
        + flow M hM hP s * (S0 - stationaryCov hbar P)
          * (-(driftMatrix M P) * flow M hM hP s)ᵀ := by
    rw [← lyapunov_stationary M P hbar hPs hPU, covFlow]
    simp only [Matrix.mul_zero, add_zero, Matrix.transpose_mul, Matrix.transpose_neg]
    noncomm_ring
  rw [key]
  convert h3 using 1
  simp only [add_zero]

/-- For every initial matrix the covariance flow tends to `ℏ P⁻¹`. -/
theorem covFlow_tendsto (hbar : ℝ) (S0 : Matrix ι ι ℝ) :
    Tendsto (fun s => covFlow hM hP hbar S0 s) atTop (𝓝 (stationaryCov hbar P)) := by
  have hE := flow_tendsto_zero hM hP
  have hET : Tendsto (fun s => (flow M hM hP s)ᵀ) atTop (𝓝 0) := by
    have h := (continuous_id.matrix_transpose.tendsto (0 : Matrix ι ι ℝ)).comp hE
    simpa using h
  have h := ((hE.mul (tendsto_const_nhds (x := S0 - stationaryCov hbar P))).mul hET).add
    (tendsto_const_nhds (x := stationaryCov hbar P))
  rw [Matrix.zero_mul, Matrix.zero_mul, zero_add] at h
  exact h

end Flow

#print axioms entry_hasDerivAt_mul
#print axioms cfc_exp_hasDerivAt
#print axioms flow_zero
#print axioms flow_add
#print axioms flow_hasDerivAt
#print axioms flow_tendsto_zero
#print axioms mean_tendsto
#print axioms covFlow_zero
#print axioms covFlow_hasDerivAt
#print axioms covFlow_tendsto

end

end OPH.SourceLaw
