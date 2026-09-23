import Mathlib.Data.Finset.Card
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic

set_option autoImplicit false

/-!
# Finite crossing changes and the sparse interface scale

The finite graph bound controls binary volume correction and thin connecting
paths. The remaining statements check the normalization and scale algebra.
The compactness, perimeter, isoperimetric and spectral limit proofs are
analytic in code/m1_interfaces/DERIVATION.md, not formalized by these lemmas.
-/

namespace OPH.M1Interfaces

variable {V : Type*} [Fintype V] [DecidableEq V]

def crossing (edges : Finset (V × V)) (label : V → Bool) : Finset (V × V) :=
  edges.filter (fun e => label e.1 ≠ label e.2)

def changed (f g : V → Bool) : Finset V :=
  Finset.univ.filter (fun v => f v ≠ g v)

def incident (edges : Finset (V × V)) (v : V) : Finset (V × V) :=
  edges.filter (fun e => e.1 = v ∨ e.2 = v)

/-- Every newly crossing edge touches a changed site. -/
theorem crossing_cover (edges : Finset (V × V)) (f g : V → Bool) :
    crossing edges f ⊆ crossing edges g ∪ (changed f g).biUnion (incident edges) := by
  intro e he
  obtain ⟨hedge, hf⟩ := Finset.mem_filter.mp he
  by_cases hg : g e.1 ≠ g e.2
  · exact Finset.mem_union_left _ (Finset.mem_filter.mpr ⟨hedge, hg⟩)
  · have heq : g e.1 = g e.2 := not_ne_iff.mp hg
    apply Finset.mem_union_right
    by_cases hfirst : f e.1 = g e.1
    · have hsecond : f e.2 ≠ g e.2 := by
        intro h
        exact hf (hfirst.trans (heq.trans h.symm))
      exact Finset.mem_biUnion.mpr ⟨e.2, by simp [changed, hsecond],
        by simp [incident, hedge]⟩
    · exact Finset.mem_biUnion.mpr ⟨e.1, by simp [changed, hfirst],
        by simp [incident, hedge]⟩

/-- Changing k labels increases an unordered cut by at most D*k. -/
theorem cut_card_le (edges : Finset (V × V)) (f g : V → Bool) (D : ℕ)
    (degree : ∀ v, (incident edges v).card ≤ D) :
    (crossing edges f).card ≤ (crossing edges g).card + D * (changed f g).card := by
  calc
    (crossing edges f).card ≤
        (crossing edges g ∪ (changed f g).biUnion (incident edges)).card :=
      Finset.card_le_card (crossing_cover edges f g)
    _ ≤ (crossing edges g).card + ((changed f g).biUnion (incident edges)).card :=
      Finset.card_union_le _ _
    _ ≤ (crossing edges g).card + ∑ v ∈ changed f g, (incident edges v).card :=
      Nat.add_le_add_left (Finset.card_biUnion_le) _
    _ ≤ (crossing edges g).card + ∑ _v ∈ changed f g, D :=
      Nat.add_le_add_left (Finset.sum_le_sum (fun v _ => degree v)) _
    _ = (crossing edges g).card + D * (changed f g).card := by simp [Nat.mul_comm]

/-- The same bound controls the absolute change, including added bridges. -/
theorem absolute_cut_change (edges : Finset (V × V)) (f g : V → Bool) (D : ℕ)
    (degree : ∀ v, (incident edges v).card ≤ D) :
    |((crossing edges f).card : ℝ) - (crossing edges g).card| ≤
      (D : ℝ) * (changed f g).card := by
  have hs : changed g f = changed f g := by ext v; simp [changed, ne_comm]
  have h₁ : ((crossing edges f).card : ℝ) ≤
      (crossing edges g).card + (D : ℝ) * (changed f g).card := by
    exact_mod_cast cut_card_le edges f g D degree
  have hn := cut_card_le edges g f D degree
  rw [hs] at hn
  have h₂ : ((crossing edges g).card : ℝ) ≤
      (crossing edges f).card + (D : ℝ) * (changed f g).card := by exact_mod_cast hn
  exact abs_le.mpr ⟨by linarith, by linarith⟩

/-- Exact volume correction of order q² has vanishing cost when D/q² does. -/
theorem normalized_volume_correction (cost D changedCount C q : ℝ)
    (hq : 0 < q) (hD : 0 ≤ D) (hc : cost ≤ D*changedCount)
    (hn : changedCount ≤ C*q^2) : cost/q^4 ≤ C*D/q^2 := by
  have h := le_trans hc (mul_le_mul_of_nonneg_left hn hD)
  have hd := div_le_div_of_nonneg_right h (by positivity : 0 ≤ q^4)
  have he : D*(C*q^2)/q^4 = C*D/q^2 := by field_simp
  rwa [he] at hd

theorem compactness_scale_identity (q m K r : ℝ) (hq : q = m^2)
    (hK : K^4 = m^3) : r^4/(m*q) = (r/K)^4 := by
  rw [div_pow, hK, hq]
  congr 1
  ring

theorem alias_action_scale_identity (q m K r : ℝ) (hm : m ≠ 0)
    (hK : K ≠ 0) (hq : q = m^2) :
    q^2*r^5/(m^4*K^5) = (r/K)^5 := by
  rw [hq, div_pow]
  field_simp

/-- Two nonnegative symbol contributions control the sum below saturation. -/
theorem sum_of_saturations (a b c : ℝ) (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 ≤ c) :
    min (a+b) c ≤ min a c + min b c := by
  by_cases hac : a ≤ c
  · rw [min_eq_left hac]
    by_cases hbc : b ≤ c
    · rw [min_eq_left hbc]
      exact min_le_left _ _
    · rw [min_eq_right (le_of_not_ge hbc)]
      have := min_le_right (a+b) c
      linarith
  · rw [min_eq_right (le_of_not_ge hac)]
    have := min_le_right (a+b) c
    have hmin : 0 ≤ min b c := le_min hb hc
    linarith

theorem explicit_repair_factors (s : ℝ) (hs : s ≠ 0) :
    (s^8*s^16)/(s^7)^4 = 1/s^4 ∧
    (s^7)^4/(s^16)^2 = 1/s^4 ∧
    (s^7/s^6)^5 = s^5 := by
  constructor
  · field_simp
  constructor
  · field_simp
  · field_simp

end OPH.M1Interfaces
