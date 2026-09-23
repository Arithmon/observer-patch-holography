import Mathlib.Algebra.BigOperators.Group.List.Basic
import Mathlib.Tactic

set_option autoImplicit false

/-!
# Sparse read families and M1 necessity

These declarations check the binary route subroutine, layer obstruction,
finite cut algebra and positive-action normalization bounds. The explicit
three-dimensional stencil, boundary tube, measure convergence and spectral
interpretation are proved analytically in code/m1_necessity/DERIVATION.md.
No native source, physical clock or full-axiom model is assumed to be produced.
-/

namespace OPH.M1Necessity

theorem cube_face_count (n : ℤ) :
    (2*n+1)^3-(2*n-1)^3 = 24*n^2+2 := by ring

/-- Every remainder below the dyadic radius has a route of at most n steps. -/
theorem binary_route (n m : ℕ) (hm : m < 2^n) :
    ∃ steps : List ℕ, steps.sum = m ∧ steps.length ≤ n ∧
      ∀ a ∈ steps, ∃ b < n, a = 2^b := by
  induction n generalizing m with
  | zero =>
      have : m = 0 := by simpa using hm
      subst m
      exact ⟨[], by simp, by simp, by simp⟩
  | succ n ih =>
      by_cases hsmall : m < 2^n
      · obtain ⟨steps, hs, hl, hb⟩ := ih m hsmall
        exact ⟨steps, hs, by omega, fun a ha => by
          obtain ⟨b, hbn, rfl⟩ := hb a ha
          exact ⟨b, by omega, rfl⟩⟩
      · have hpow : 2^(n+1) = 2^n*2 := by rw [pow_succ]
        have hrem : m-2^n < 2^n := by omega
        obtain ⟨steps, hs, hl, hb⟩ := ih (m-2^n) hrem
        refine ⟨2^n :: steps, ?_, by simpa using Nat.succ_le_succ hl, ?_⟩
        · simp only [List.sum_cons, hs]
          omega
        · intro a ha
          simp only [List.mem_cons] at ha
          rcases ha with rfl | ha
          · exact ⟨n, Nat.lt_succ_self n, rfl⟩
          · obtain ⟨b, hbn, rfl⟩ := hb a ha
            exact ⟨b, by omega, rfl⟩

/-- Repeated radius flights plus binary remainders account for every step. -/
theorem axis_route (n m : ℕ) :
    ∃ steps : List ℕ, steps.sum = m ∧ steps.length ≤ m / 2^n+n ∧
      ∀ a ∈ steps, ∃ b ≤ n, a = 2^b := by
  have hp : 0 < 2^n := by positivity
  obtain ⟨small, hs, hl, hb⟩ := binary_route n (m % 2^n) (Nat.mod_lt _ hp)
  refine ⟨List.replicate (m / 2^n) (2^n) ++ small, ?_, ?_, ?_⟩
  · simp only [List.sum_append, List.sum_replicate, smul_eq_mul, hs]
    simpa [Nat.mul_comm, Nat.add_comm] using Nat.mod_add_div m (2^n)
  · simp only [List.length_append, List.length_replicate]
    omega
  · intro a ha
    rcases List.mem_append.mp ha with ha | ha
    · have he : a = 2^n := (List.mem_replicate.mp ha).2
      exact ⟨n, le_rfl, he⟩
    · obtain ⟨b, hbn, he⟩ := hb a ha
      exact ⟨b, Nat.le_of_lt hbn, he⟩

theorem route_budget (r R n k residual corrections : ℝ)
    (hR : 0 < R) (hn : 0 < n) (hk : k ≤ r/R)
    (he : residual ≤ 3*r/n+3*r/R+2*R)
    (hc : corrections ≤ residual/R+3*n) :
    k+corrections ≤ (1+3/n+3/R)*r/R+3*n+2 := by
  have hdiv := div_le_div_of_nonneg_right he (le_of_lt hR)
  have identity : r/R+(3*r/n+3*r/R+2*R)/R+3*n =
      (1+3/n+3/R)*r/R+3*n+2 := by field_simp; ring
  linarith

theorem strict_layer_increase {E : Type*} (height : E → ℕ) (edge : E → E → Prop)
    (grade : ∀ a b, edge a b → height b = height a+1)
    {a b : E} (path : Relation.TransGen edge a b) : height a < height b := by
  induction path with
  | single h => have := grade _ _ h; omega
  | tail _ h ih => have := grade _ _ h; omega

/-- Adjacent-layer comparability cannot survive deletion of its sole direct edge. -/
theorem adjacent_layers_require_edge {E : Type*} (height : E → ℕ)
    (edge : E → E → Prop) (grade : ∀ a b, edge a b → height b = height a+1)
    {a b : E} (hab : height b = height a+1) :
    Relation.TransGen edge a b ↔ edge a b := by
  constructor
  · intro h
    cases h with
    | single e => exact e
    | tail path e =>
        have := strict_layer_increase height edge grade path
        have := grade _ _ e
        omega
  · exact Relation.TransGen.single

theorem sum_lengths (r : ℕ) :
    2 * ∑ j ∈ Finset.range (r+1), (j : ℝ) = (r : ℝ)*(r+1) := by
  induction r with
  | zero => simp
  | succ r ih =>
      rw [Finset.sum_range_succ]
      push_cast
      nlinarith

theorem sum_squared_lengths (r : ℕ) :
    6 * ∑ j ∈ Finset.range (r+1), (j : ℝ)^2 = (r : ℝ)*(r+1)*(2*r+1) := by
  induction r with
  | zero => simp
  | succ r ih =>
      rw [Finset.sum_range_succ]
      push_cast
      nlinarith

theorem axial_coordinate_cut (q : ℝ) (r : ℕ) :
    ∑ j ∈ Finset.range (r+1), (j : ℝ)*q^2 = q^2*r*(r+1)/2 := by
  rw [← Finset.sum_mul]
  have := sum_lengths r
  nlinarith

theorem axial_diagonal_cut (q : ℝ) (r : ℕ) :
    ∑ j ∈ Finset.range (r+1), 2*q*(j : ℝ)*(q-j) =
      q^2*r*(r+1)-q*r*(r+1)*(2*r+1)/3 := by
  have hsum : (∑ j ∈ Finset.range (r+1), 2*q*(j : ℝ)*(q-j)) =
      2*q^2*(∑ j ∈ Finset.range (r+1), (j : ℝ))-
      2*q*(∑ j ∈ Finset.range (r+1), (j : ℝ)^2) := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    apply Finset.sum_congr rfl
    intro j hj
    ring
  rw [hsum]
  have h1 := sum_lengths r
  have h2 := sum_squared_lengths r
  have e1 : (∑ j ∈ Finset.range (r+1), (j : ℝ)) = (r : ℝ)*(r+1)/2 := by linarith
  have e2 : (∑ j ∈ Finset.range (r+1), (j : ℝ)^2) = (r : ℝ)*(r+1)*(2*r+1)/6 := by linarith
  rw [e1, e2]
  ring

theorem sparse_dense_cut_separation (sparse dense D R q : ℝ)
    (hD : 0 ≤ D) (hR : 0 < R) (hq : 0 < q)
    (hs : sparse ≤ D*R*q^2) (hd : R^4*q^2/32 ≤ dense) :
    sparse ≤ (32*D/R^3)*dense := by
  have hc : 0 ≤ 32*D/R^3 := by positivity
  have := mul_le_mul_of_nonneg_left hd hc
  have he : (32*D/R^3)*(R^4*q^2/32) = D*R*q^2 := by field_simp
  rw [he] at this
  exact le_trans hs this

/-- A common entropy conversion does not change the orientation discrepancy. -/
theorem scalar_calibration_cancels (x y ax ay scale : ℝ)
    (hx : x ≠ 0) (hax : ax ≠ 0) (hay : ay ≠ 0) (hs : scale ≠ 0) :
    ((scale*y)/ay)/((scale*x)/ax) = (y/ay)/(x/ax) := by field_simp

theorem isotropic_moment_normalization (moment diagonal h : ℝ)
    (hm : moment ≠ 0) (hh : h ≠ 0) (hd : 3*diagonal = moment) :
    (6/(h^2*moment))*h^2*diagonal = 2 := by
  field_simp
  nlinarith

/-- The exact fourth moment bounds the spatial consistency remainder. -/
theorem normalized_fourth_remainder (moment fourth h R M : ℝ)
    (hm : 0 < moment) (hh : 0 < h) (hM : 0 ≤ M)
    (hf : fourth ≤ R^2*moment) :
    (6/(h^2*moment))*h^4*fourth*M/24 ≤ (h*R)^2*M/4 := by
  have hp : 0 ≤ (6/(h^2*moment))*h^4*M/24 := by positivity
  have bound := mul_le_mul_of_nonneg_left hf hp
  have equality : (6/(h^2*moment))*h^4*M/24*(R^2*moment) = (h*R)^2*M/4 := by
    field_simp
    ring
  rw [equality] at bound
  nlinarith

/-- Positive isotropic weights cannot satisfy centered-leapfrog stability at a/c. -/
theorem full_radius_tick_instability (a c dt mass eigenvalue : ℝ)
    (ha : 0 < a) (hc : 0 < c) (ht : dt = a/c) (hm : 0 ≤ mass)
    (he : 6 ≤ a^2*eigenvalue) : 4 < dt^2*(c^2*eigenvalue+mass) := by
  subst dt
  have positive_mass : 0 ≤ (a/c)^2*mass := by positivity
  have identity : (a/c)^2*c^2*eigenvalue = a^2*eigenvalue := by field_simp
  nlinarith

/-- The explicit larger radius and spaced reads retain the dense cut coefficient. -/
theorem balanced_crossing_scale (s : ℝ) : s^4*(s^3)^4 = (s^8)^2 := by ring

/-- Removing a positive fraction of a fixed-radius stencil has a cut cost. -/
theorem removed_fraction_enclosure (removed weight eta R : ℝ)
    (hR : 0 < R) (heta : 0 < eta)
    (bound : removed ≤ (2*eta*R+1)^3+weight/(eta*R)) :
    removed/R^3 ≤ (2*eta+1/R)^3+weight/(eta*R^4) := by
  have divided := div_le_div_of_nonneg_right bound (by positivity : 0 ≤ R^3)
  have equality : ((2*eta*R+1)^3+weight/(eta*R))/R^3 =
      (2*eta+1/R)^3+weight/(eta*R^4) := by
    field_simp
  rwa [equality] at divided

end OPH.M1Necessity
