import Mathlib.Analysis.Normed.Module.Basic
import Mathlib.Data.Fintype.Card
import Mathlib.Tactic

set_option autoImplicit false

/-!
# Local-clock and classical-information reductions for record gluing

The source realization, local clock expansion and time-faithful refinement are
physical hypotheses. These theorems prove algebraic, path, limiting and resource
consequences; they do not construct a complete native A1--A3 source.
-/

namespace OPH.RecordGluingPrinciple

def quadratic (a b c d e f x y z : ℝ) : ℝ :=
  a*x*x + b*y*y + c*z*z + 2*d*x*y + 2*e*x*z + 2*f*y*z

/-- Two half-turns and a cyclic rotation suffice for quadratic isotropy. -/
theorem finite_symmetry_quadratic (a b c d e f : ℝ)
    (h₁ : ∀ x y z, quadratic a b c d e f x (-y) (-z) = quadratic a b c d e f x y z)
    (h₂ : ∀ x y z, quadratic a b c d e f (-x) y (-z) = quadratic a b c d e f x y z)
    (h₃ : ∀ x y z, quadratic a b c d e f y z x = quadratic a b c d e f x y z) :
    a = c ∧ b = c ∧ d = 0 ∧ e = 0 ∧ f = 0 := by
  have hd := h₁ 1 1 0
  have he := h₁ 1 0 1
  have hf := h₂ 0 1 1
  have ha := h₃ 1 0 0
  have hb := h₃ 0 1 0
  simp only [quadratic] at *
  constructor
  · nlinarith
  constructor
  · nlinarith
  exact ⟨by nlinarith, by nlinarith, by nlinarith⟩

variable {E : Type*} [NormedAddCommGroup E]

/-- Summing all step costs bounds displacement; auxiliary steps are retained. -/
theorem graded_path_bound (p : ℕ → E) (cost : ℕ → ℝ) (k : ℕ) (tau : ℝ)
    (ht : 0 ≤ tau) (h : ∀ j < k, tau * ‖p (j+1) - p j‖ ≤ cost j) :
    tau * ‖p k - p 0‖ ≤ ∑ j ∈ Finset.range k, cost j := by
  induction k with
  | zero => simp
  | succ k ih =>
    have hi := ih (fun j hj => h j (Nat.lt_succ_of_lt hj))
    have hs := h k (Nat.lt_succ_self k)
    have triangle : ‖p (k+1) - p 0‖ ≤ ‖p (k+1) - p k‖ + ‖p k - p 0‖ := by
      simpa using norm_add_le (p (k+1) - p k) (p k - p 0)
    have scaled := mul_le_mul_of_nonneg_left triangle ht
    rw [Finset.sum_range_succ]
    nlinarith

/-- Local lower bounds and arbitrarily faithful timing exclude coarse shortcuts. -/
theorem refinement_clock_bound (tau radius duration : ℝ) (hr : 0 ≤ radius)
    (h : ∀ epsilon : ℝ, 0 < epsilon →
      (tau - epsilon) * radius ≤ duration + epsilon) :
    tau * radius ≤ duration := by
  by_contra hn
  have gap : 0 < tau * radius - duration := by linarith
  let epsilon := (tau * radius - duration) / (2 * (radius + 1))
  have hp : 0 < epsilon := by dsimp [epsilon]; positivity
  have he := h epsilon hp
  have identity : epsilon * (2 * (radius + 1)) = tau * radius - duration := by
    dsimp [epsilon]
    field_simp
  nlinarith

noncomputable def subdividedTime (tau kappa radius count : ℝ) : ℝ :=
  count * (tau * (radius / count) + kappa * (radius / count)^2)

theorem subdivided_time_formula (tau kappa radius count : ℝ) (hn : count ≠ 0) :
    subdividedTime tau kappa radius count = tau * radius + kappa * radius^2 / count := by
  dsimp [subdividedTime]
  field_simp

theorem finite_positive_overhead (tau kappa radius count : ℝ)
    (hk : 0 < kappa) (hr : 0 < radius) (hn : 0 < count) :
    tau * radius < subdividedTime tau kappa radius count := by
  rw [subdivided_time_formula tau kappa radius count (ne_of_gt hn)]
  have : 0 < kappa * radius^2 / count := by positivity
  linarith

/-- Strict timelike deadlines can be met at a sufficiently large finite subdivision. -/
theorem exists_subdivision_before (tau kappa radius deadline : ℝ)
    (hmargin : tau * radius < deadline) :
    ∃ n : ℕ, 0 < n ∧ subdividedTime tau kappa radius n < deadline := by
  obtain ⟨n, hn⟩ := exists_nat_gt (max 0 (kappa * radius^2 / (deadline - tau * radius)))
  have hnreal : (0 : ℝ) < n := lt_of_le_of_lt (le_max_left _ _) hn
  have hnb : kappa * radius^2 / (deadline - tau * radius) < n :=
    lt_of_le_of_lt (le_max_right _ _) hn
  have hm : 0 < deadline - tau * radius := by linarith
  have hc : kappa * radius^2 < (n : ℝ) * (deadline - tau * radius) := (div_lt_iff₀ hm).mp hnb
  refine ⟨n, by exact_mod_cast hnreal, ?_⟩
  rw [subdivided_time_formula tau kappa radius n (ne_of_gt hnreal)]
  have := (div_lt_iff₀ hnreal).mpr (by nlinarith : kappa * radius^2 < (deadline - tau * radius) * n)
  linarith

/-- A decoder for every independently prepared tuple forces transcript injectivity. -/
theorem lossless_injective {A B : Type*} (encode : A → B) (decode : B → A)
    (h : ∀ x, decode (encode x) = x) : Function.Injective encode := by
  intro x y he
  calc
    x = decode (encode x) := (h x).symm
    _ = decode (encode y) := congrArg decode he
    _ = y := h y

theorem transcript_capacity {A B : Type*} [Fintype A] [Fintype B]
    (encode : A → B) (decode : B → A) (h : ∀ x, decode (encode x) = x) :
    Fintype.card A ≤ Fintype.card B :=
  Fintype.card_le_of_injective encode (lossless_injective encode decode h)

theorem independent_records_capacity (d m : ℕ) {T : Type*} [Fintype T]
    (encode : (Fin m → Fin d) → T) (decode : T → (Fin m → Fin d))
    (h : ∀ x, decode (encode x) = x) : d^m ≤ Fintype.card T := by
  simpa using transcript_capacity encode decode h

/-- Retain input and XOR the computed answer into an output register. -/
def reversibleEvaluation {A I : Type*} (f : A → I → Bool) (state : A × (I → Bool)) :
    A × (I → Bool) := (state.1, fun j => Bool.xor (state.2 j) (f state.1 j))

theorem reversible_evaluation_involutive {A I : Type*} (f : A → I → Bool) :
    Function.Involutive (reversibleEvaluation f) := by
  intro state
  rcases state with ⟨x, y⟩
  simp only [reversibleEvaluation, Prod.mk.injEq, true_and]
  funext j
  cases y j <;> cases f x j <;> rfl

theorem reversible_evaluation_blank {A I : Type*} (f : A → I → Bool) (x : A) :
    reversibleEvaluation f (x, fun _ => false) = (x, f x) := by
  simp [reversibleEvaluation]

/-- An interior cube gives a strict, boundary-independent read-count lower bound. -/
theorem dense_inner_cube (n : ℕ) (hn : 2 ≤ n) :
    3 * (n / 2)^2 < n^2 ∧ n ≤ 2 * (n / 2) + 1 ∧ n^2 ≤ 2 * (n^2 - 2 * (n / 2)) := by
  have hm : 2 * (n / 2) ≤ n := Nat.mul_div_le n 2
  have hrem : n ≤ 2 * (n / 2) + 1 := by omega
  have hn2 : n ≤ n^2 := by nlinarith
  have hsub : 2 * (n / 2) ≤ n^2 := le_trans hm hn2
  have hrestore := Nat.sub_add_cancel hsub
  exact ⟨by nlinarith [sq_nonneg (n - 2 * (n / 2) : ℤ)], hrem, by nlinarith⟩

end OPH.RecordGluingPrinciple
